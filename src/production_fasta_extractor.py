#!/usr/bin/env python3
"""
Production-Grade FastaAAExtractor
Advanced protein extraction with provenance tracking and CARD RGI integration

This module provides a production-ready protein extraction system that:
1. Integrates seamlessly with enhanced CARD RGI coordinate system
2. Maintains accession-based naming consistency across pipeline
3. Provides comprehensive provenance tracking and metadata
4. Generates extraction manifests for complete traceability
5. Handles batch processing with integrity checking

Author: GenomeAMRAnalyzer Pipeline
Version: 2.0 - Production Grade
"""

import os
import sys
import time
import hashlib
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Set
from dataclasses import dataclass
import json
import csv

# BioPython imports
try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False
    print("Warning: BioPython not available. Some features may be limited.")


@dataclass
class EnhancedGeneCoordinate:
    """Gene coordinates with full CARD RGI provenance"""
    accession: str
    gene_name: str
    contig_id: str
    start: int
    end: int
    strand: str
    cut_off: str
    pass_bitscore: float
    best_hit_aro: str
    model_type: str
    drug_class: str
    resistance_mechanism: str
    amr_gene_family: str
    analysis_timestamp: str
    rgi_version: str
    card_version: str
    coordinate_file_path: str = ""

    @property
    def length(self) -> int:
        """Return gene length in base pairs"""
        return abs(self.end - self.start) + 1
    
    @property
    def is_forward_strand(self) -> bool:
        """Check if gene is on forward strand"""
        return self.strand in ['+', '1']


@dataclass
class ProteinExtractionResult:
    """Result of protein extraction with complete provenance"""
    accession: str
    gene_name: str
    protein_sequence: str
    extraction_success: bool = True
    error_message: str = ""
    
    # Source information
    source_genome_file: str = ""
    coordinate_source_file: str = ""
    
    # Sequence metadata
    sequence_length_aa: int = 0
    sequence_checksum: str = ""
    start_position: int = 0
    end_position: int = 0
    strand: str = "+"
    contig_id: str = ""
    
    # Provenance metadata
    extraction_timestamp: str = ""
    extraction_method: str = "coordinate_based"
    rgi_version: str = ""
    card_version: str = ""
    
    # Quality metrics
    pass_bitscore: float = 0.0
    cut_off: str = ""


class ProductionFastaExtractor:
    """
    Production-grade protein extractor with full provenance tracking
    """
    
    def __init__(self, output_dir: str, coordinate_manifest_path: Optional[str] = None):
        """Initialize production extractor"""
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Create output subdirectories
        (self.output_dir / "proteins").mkdir(exist_ok=True)
        (self.output_dir / "logs").mkdir(exist_ok=True)
        (self.output_dir / "manifests").mkdir(exist_ok=True)
        
        self.coordinate_manifest_path = coordinate_manifest_path
        
        # Setup logging
        self._setup_logging()
        
        # Data storage
        self.extraction_results: Dict[str, List[ProteinExtractionResult]] = {}
        self.failed_extractions: Dict[str, str] = {}
        self.processed_accessions: Set[str] = set()
        
        # Statistics
        self.stats = {
            'total_accessions': 0,
            'successful_accessions': 0,
            'failed_accessions': 0,
            'total_proteins_extracted': 0,
            'total_extraction_failures': 0
        }
        
        self.logger.info("Production FastaExtractor initialized")
    
    def _setup_logging(self):
        """Setup comprehensive logging"""
        log_file = self.output_dir / "logs" / "protein_extraction.log"
        
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        
        self.logger = logging.getLogger('ProductionFastaExtractor')
        self.logger.setLevel(logging.INFO)
        
        # File handler
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(formatter)
        self.logger.addHandler(file_handler)
    
    def extract_accession_coordinates_from_rgi(self, accession: str, 
                                             coordinate_file: str) -> List[EnhancedGeneCoordinate]:
        """Load coordinates for specific accession from RGI coordinate file"""
        coordinates = []
        
        try:
            with open(coordinate_file, 'r') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    if row['accession'] == accession:
                        coord = EnhancedGeneCoordinate(
                            accession=row['accession'],
                            gene_name=row['gene_name'],
                            contig_id=row['contig_id'],
                            start=int(row['start']),
                            end=int(row['end']),
                            strand=row['strand'],
                            cut_off=row['cut_off'],
                            pass_bitscore=float(row['pass_bitscore']),
                            best_hit_aro=row['best_hit_aro'],
                            model_type=row['model_type'],
                            drug_class=row['drug_class'],
                            resistance_mechanism=row['resistance_mechanism'],
                            amr_gene_family=row['amr_gene_family'],
                            analysis_timestamp=row['analysis_timestamp'],
                            rgi_version=row['rgi_version'],
                            card_version=row['card_version'],
                            coordinate_file_path=coordinate_file
                        )
                        coordinates.append(coord)
            
            self.logger.info(f"Loaded {len(coordinates)} coordinates for {accession}")
            return coordinates
            
        except Exception as e:
            self.logger.error(f"Failed to load coordinates for {accession}: {e}")
            return []
    
    def extract_proteins_from_genome(self, accession: str, genome_file: str, 
                                   coordinates: List[EnhancedGeneCoordinate],
                                   target_genes: Optional[List[str]] = None) -> List[ProteinExtractionResult]:
        """Extract proteins from genome using CARD coordinates"""
        
        if not BIOPYTHON_AVAILABLE:
            raise RuntimeError("BioPython is required for protein extraction")
        
        extraction_results = []
        current_timestamp = time.strftime('%Y-%m-%d %H:%M:%S')
        
        # Filter coordinates by target genes if specified
        if target_genes:
            target_genes_lower = [g.lower() for g in target_genes]
            coordinates = [c for c in coordinates if c.gene_name.lower() in target_genes_lower]
        
        if not coordinates:
            self.logger.warning(f"No target coordinates found for {accession}")
            return []
        
        try:
            # Load genome file
            genome_records = {}
            for record in SeqIO.parse(genome_file, "fasta"):
                genome_records[record.id] = record
            
            self.logger.info(f"Loaded {len(genome_records)} contigs for {accession}")
            
            # Extract each protein
            for coord in coordinates:
                try:
                    # Find matching contig
                    target_record = None
                    if coord.contig_id in genome_records:
                        target_record = genome_records[coord.contig_id]
                    else:
                        # Try to find by partial match or first record
                        for record_id, record in genome_records.items():
                            if coord.contig_id in record_id or record_id in coord.contig_id:
                                target_record = record
                                break
                        
                        if not target_record and len(genome_records) == 1:
                            target_record = list(genome_records.values())[0]
                    
                    if not target_record:
                        result = ProteinExtractionResult(
                            accession=accession,
                            gene_name=coord.gene_name,
                            protein_sequence="",
                            extraction_success=False,
                            error_message=f"Contig {coord.contig_id} not found in genome",
                            source_genome_file=genome_file,
                            coordinate_source_file=coord.coordinate_file_path,
                            extraction_timestamp=current_timestamp
                        )
                        extraction_results.append(result)
                        continue
                    
                    # Extract nucleotide sequence
                    try:
                        # Adjust coordinates (1-based to 0-based)
                        start = max(0, coord.start - 1)
                        end = min(len(target_record.seq), coord.end)
                        
                        nucleotide_seq = target_record.seq[start:end]
                        
                        # Handle reverse strand
                        if not coord.is_forward_strand:
                            nucleotide_seq = nucleotide_seq.reverse_complement()
                        
                        # Translate to protein
                        protein_seq = nucleotide_seq.translate()
                        
                        # Remove stop codons
                        protein_seq_str = str(protein_seq).rstrip('*')
                        
                        # Calculate checksum
                        checksum = hashlib.sha256(protein_seq_str.encode()).hexdigest()[:16]
                        
                        result = ProteinExtractionResult(
                            accession=accession,
                            gene_name=coord.gene_name,
                            protein_sequence=protein_seq_str,
                            extraction_success=True,
                            source_genome_file=genome_file,
                            coordinate_source_file=coord.coordinate_file_path,
                            sequence_length_aa=len(protein_seq_str),
                            sequence_checksum=checksum,
                            start_position=coord.start,
                            end_position=coord.end,
                            strand=coord.strand,
                            contig_id=coord.contig_id,
                            extraction_timestamp=current_timestamp,
                            rgi_version=coord.rgi_version,
                            card_version=coord.card_version,
                            pass_bitscore=coord.pass_bitscore,
                            cut_off=coord.cut_off
                        )
                        
                        extraction_results.append(result)
                        self.logger.debug(f"Extracted {coord.gene_name} for {accession}: {len(protein_seq_str)} AA")
                        
                    except Exception as e:
                        result = ProteinExtractionResult(
                            accession=accession,
                            gene_name=coord.gene_name,
                            protein_sequence="",
                            extraction_success=False,
                            error_message=f"Translation error: {str(e)}",
                            source_genome_file=genome_file,
                            coordinate_source_file=coord.coordinate_file_path,
                            extraction_timestamp=current_timestamp
                        )
                        extraction_results.append(result)
                        
                except Exception as e:
                    self.logger.error(f"Error extracting {coord.gene_name} for {accession}: {e}")
            
            self.logger.info(f"Extracted {len([r for r in extraction_results if r.extraction_success])} proteins for {accession}")
            return extraction_results
            
        except Exception as e:
            self.logger.error(f"Failed to process genome {accession}: {e}")
            return []
    
    def save_accession_proteins(self, accession: str, 
                              extraction_results: List[ProteinExtractionResult]) -> str:
        """Save extracted proteins for specific accession to FASTA file"""
        
        output_file = self.output_dir / "proteins" / f"{accession}_proteins.fasta"
        
        try:
            with open(output_file, 'w') as f:
                for result in extraction_results:
                    if result.extraction_success and result.protein_sequence:
                        # Create comprehensive FASTA header
                        header = (f">{accession}|{result.gene_name}|{result.start_position}-{result.end_position}|"
                                f"{result.strand}|{result.contig_id}|{result.sequence_checksum}")
                        f.write(f"{header}\n")
                        f.write(f"{result.protein_sequence}\n")
            
            successful_count = len([r for r in extraction_results if r.extraction_success])
            self.logger.info(f"Saved {successful_count} proteins for {accession} to {output_file}")
            return str(output_file)
            
        except Exception as e:
            self.logger.error(f"Failed to save proteins for {accession}: {e}")
            return ""
    
    def generate_extraction_manifest(self, output_file: Optional[str] = None) -> str:
        """Generate comprehensive extraction manifest"""
        
        if output_file is None:
            output_file = str(self.output_dir / "manifests" / "extraction_manifest.json")
        
        manifest = {
            "manifest_version": "1.0",
            "generated_timestamp": time.strftime('%Y-%m-%d %H:%M:%S'),
            "extractor_version": "2.0",
            "total_accessions_processed": len(self.extraction_results),
            "successful_accessions": len([acc for acc, results in self.extraction_results.items() 
                                        if any(r.extraction_success for r in results)]),
            "failed_accessions": len(self.failed_extractions),
            "total_proteins_extracted": sum(len([r for r in results if r.extraction_success]) 
                                          for results in self.extraction_results.values()),
            "accessions": {}
        }
        
        # Add detailed data for each accession
        for accession, results in self.extraction_results.items():
            successful_results = [r for r in results if r.extraction_success]
            manifest["accessions"][accession] = {
                "protein_file": str(self.output_dir / "proteins" / f"{accession}_proteins.fasta"),
                "proteins_extracted": len(successful_results),
                "extraction_failures": len(results) - len(successful_results),
                "proteins": [
                    {
                        "gene_name": r.gene_name,
                        "sequence_length_aa": r.sequence_length_aa,
                        "start_position": r.start_position,
                        "end_position": r.end_position,
                        "strand": r.strand,
                        "contig_id": r.contig_id,
                        "sequence_checksum": r.sequence_checksum,
                        "extraction_timestamp": r.extraction_timestamp,
                        "pass_bitscore": r.pass_bitscore
                    }
                    for r in successful_results
                ]
            }
        
        # Add failed accessions
        if self.failed_extractions:
            manifest["failed_accessions"] = self.failed_extractions
        
        # Write JSON manifest
        Path(output_file).parent.mkdir(parents=True, exist_ok=True)
        with open(output_file, 'w') as f:
            json.dump(manifest, f, indent=2)
        
        self.logger.info(f"Extraction manifest saved to: {output_file}")
        return output_file
    
    def process_batch_from_coordinates(self, genome_dir: str, coordinate_dir: str, 
                                     target_genes: Optional[List[str]] = None) -> bool:
        """Process batch of genomes using coordinate directory structure"""
        
        genome_path = Path(genome_dir)
        coord_path = Path(coordinate_dir)
        
        if not genome_path.exists() or not coord_path.exists():
            self.logger.error(f"Input directories not found: {genome_path}, {coord_path}")
            return False
        
        # Find coordinate files
        coord_files = list(coord_path.glob("*_coordinates.csv"))
        
        if not coord_files:
            self.logger.error(f"No coordinate files found in {coord_path}")
            return False
        
        self.logger.info(f"Processing {len(coord_files)} accessions")
        self.stats['total_accessions'] = len(coord_files)
        
        for coord_file in coord_files:
            # Extract accession from filename
            accession = coord_file.stem.replace('_coordinates', '')
            
            if accession in self.processed_accessions:
                self.logger.warning(f"Accession {accession} already processed")
                continue
            
            # Find corresponding genome file
            genome_file = None
            for pattern in [f"{accession}.fasta", f"{accession}_*.fasta"]:
                matches = list(genome_path.glob(pattern))
                if matches:
                    genome_file = matches[0]
                    break
            
            if not genome_file:
                self.logger.error(f"Genome file not found for {accession}")
                self.failed_extractions[accession] = "Genome file not found"
                self.stats['failed_accessions'] += 1
                continue
            
            try:
                # Load coordinates
                coordinates = self.extract_accession_coordinates_from_rgi(accession, str(coord_file))
                
                if not coordinates:
                    self.logger.warning(f"No coordinates found for {accession}")
                    self.failed_extractions[accession] = "No coordinates found"
                    self.stats['failed_accessions'] += 1
                    continue
                
                # Extract proteins
                results = self.extract_proteins_from_genome(
                    accession, str(genome_file), coordinates, target_genes
                )
                
                if results:
                    # Save results
                    self.extraction_results[accession] = results
                    protein_file = self.save_accession_proteins(accession, results)
                    
                    # Update statistics
                    successful_extractions = len([r for r in results if r.extraction_success])
                    self.stats['total_proteins_extracted'] += successful_extractions
                    self.stats['total_extraction_failures'] += len(results) - successful_extractions
                    
                    if successful_extractions > 0:
                        self.stats['successful_accessions'] += 1
                        self.logger.info(f"Successfully processed {accession}: {successful_extractions} proteins")
                    else:
                        self.stats['failed_accessions'] += 1
                        self.failed_extractions[accession] = "All protein extractions failed"
                else:
                    self.stats['failed_accessions'] += 1
                    self.failed_extractions[accession] = "Extraction processing failed"
                
                self.processed_accessions.add(accession)
                
            except Exception as e:
                self.logger.error(f"Error processing {accession}: {e}")
                self.failed_extractions[accession] = str(e)
                self.stats['failed_accessions'] += 1
        
        # Generate final manifest
        self.generate_extraction_manifest()
        
        # Log final statistics
        self.logger.info(f"Batch processing complete. Processed {len(self.processed_accessions)} accessions")
        self.logger.info(f"Successful: {self.stats['successful_accessions']}, Failed: {self.stats['failed_accessions']}")
        self.logger.info(f"Total proteins extracted: {self.stats['total_proteins_extracted']}")
        
        return self.stats['successful_accessions'] > 0


def main():
    """Command line interface for production protein extraction"""
    import argparse
    
    parser = argparse.ArgumentParser(description="Production FastaAAExtractor with CARD RGI integration")
    parser.add_argument("--genome-dir", required=True, help="Directory containing genome FASTA files")
    parser.add_argument("--coordinate-dir", required=True, help="Directory containing CARD coordinate CSV files")
    parser.add_argument("--output-dir", required=True, help="Output directory for extracted proteins")
    parser.add_argument("--target-genes", nargs="+", help="Target gene names to extract")
    parser.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    
    args = parser.parse_args()
    
    # Setup logging
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    
    # Initialize extractor
    extractor = ProductionFastaExtractor(args.output_dir)
    
    # Process batch
    success = extractor.process_batch_from_coordinates(
        args.genome_dir, args.coordinate_dir, args.target_genes
    )
    
    if success:
        print("Protein extraction completed successfully")
        return 0
    else:
        print("Protein extraction failed")
        return 1


if __name__ == "__main__":
    sys.exit(main())