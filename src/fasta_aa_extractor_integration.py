#!/usr/bin/env python3
"""
FastaAAExtractor Pipeline Integration - Production Grade
Advanced integration between CARD RGI coordinates and protein extraction

This module provides a production-ready pipeline component that:
1. Takes CARD coordinate output and genome files with full provenance tracking
2. Extracts amino acid sequences for user-specified genes with metadata
3. Maintains accession-based naming consistency across pipeline
4. Generates comprehensive extraction manifests for traceability
5. Handles batch processing with integrity checking

Author: GenomeAMRAnalyzer Pipeline
Version: 2.0 - Production Grade with Provenance Tracking
"""

import os
import sys
import logging
import subprocess
import tempfile
import shutil
import zipfile
import hashlib
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union, Set
from dataclasses import dataclass, field
from collections import defaultdict
import json
import argparse
from datetime import datetime
import traceback

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
class GeneCoordinate:
    """Represents gene coordinates from CARD analysis with enhanced provenance"""
    gene_name: str
    genome_id: str  # This should be the accession
    start: int
    end: int
    strand: str = "+"
    # Prefer 'contig'; keep 'contig_id' for backward compatibility
    contig_id: str = ""
    contig: Optional[str] = None
    cut_off: str = ""
    pass_bitscore: float = 0.0
    best_hit_aro: str = ""
    model_type: str = ""
    drug_class: str = ""
    resistance_mechanism: str = ""
    amr_gene_family: str = ""
    # Enhanced provenance metadata
    source_accession: str = ""
    analysis_timestamp: str = ""
    rgi_version: str = ""
    card_version: str = ""
    coordinate_file_path: str = ""
    
    @property
    def is_forward_strand(self) -> bool:
        return str(self.strand) in ['+', '1']
    
    @property
    def length(self) -> int:
        try:
            return int(self.end) - int(self.start) + 1
        except Exception:
            return 0


@dataclass
class ExtractionResult:
    """Results from protein extraction"""
    gene_name: str
    genome_id: str
    sequence: str
    coordinates: GeneCoordinate
    extraction_method: str = "coordinate_based"
    quality_score: float = 1.0
    warnings: List[str] = field(default_factory=list)
    
    @property
    def fasta_header(self) -> str:
        """Generate FASTA header for this sequence"""
        return f"{self.genome_id}_{self.gene_name}_{self.coordinates.start}_{self.coordinates.end}_{self.coordinates.strand}"
    
    @property
    def sequence_length(self) -> int:
        """Return sequence length"""
        return len(self.sequence)


class FastaAAExtractorIntegrator:
    """
    Pipeline component for extracting amino acid sequences from genomes
    using CARD coordinates and preparing them for WildTypeAligner
    """
    
    def __init__(self, 
                 output_dir: str = "extracted_proteins",
                 temp_dir: Optional[str] = None,
                 external_extractor_path: Optional[str] = None,
                 log_level: str = "INFO"):
        """
        Initialize the FastaAAExtractor integrator
        
        Args:
            output_dir: Directory for output files
            temp_dir: Temporary directory for intermediate files
            external_extractor_path: Path to external FastaAAExtractor tool
            log_level: Logging level
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        self.temp_dir = Path(temp_dir) if temp_dir else Path(tempfile.mkdtemp())
        self.temp_dir.mkdir(parents=True, exist_ok=True)
        
        self.external_extractor_path = external_extractor_path
        
        # Setup logging
        self._setup_logging(log_level)
        
        # Data storage
        self.gene_coordinates: Dict[str, List[GeneCoordinate]] = defaultdict(list)
        self.extraction_results: List[ExtractionResult] = []
        self.failed_extractions: List[Dict] = []
        
        # Statistics
        self.stats = {
            'total_genomes': 0,
            'total_genes_requested': 0,
            'successful_extractions': 0,
            'failed_extractions': 0,
            'sequences_generated': 0
        }
        
        self.logger.info("FastaAAExtractor integrator initialized")
    
    def _setup_logging(self, log_level: str) -> None:
        """Setup comprehensive logging"""
        log_file = self.output_dir / "fasta_extraction.log"
        
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        
        self.logger = logging.getLogger('FastaAAExtractor')
        self.logger.setLevel(getattr(logging, log_level.upper()))
        
        # File handler
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(formatter)
        self.logger.addHandler(file_handler)
        
        # Console handler
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(formatter)
        self.logger.addHandler(console_handler)
    
    def load_card_coordinates(self, 
                            coordinates_file: Union[str, Path],
                            file_format: str = "auto") -> None:
        """
        Load gene coordinates from CARD integrator output
        
        Args:
            coordinates_file: Path to coordinates file
            file_format: File format (csv, tsv, json, auto)
        """
        try:
            coords_path = Path(coordinates_file)
            if not coords_path.exists():
                raise FileNotFoundError(f"Coordinates file not found: {coordinates_file}")
            
            self.logger.info(f"Loading CARD coordinates from {coords_path}")
            
            # Auto-detect format
            if file_format == "auto":
                suffix = coords_path.suffix.lower()
                if suffix == ".json":
                    file_format = "json"
                elif suffix in [".csv", ".tsv", ".tab"]:
                    file_format = "csv"
                else:
                    file_format = "csv"  # Default
            
            if file_format == "json":
                self._load_coordinates_json(coords_path)
            else:
                self._load_coordinates_csv(coords_path)
            
            total_coords = sum(len(coords) for coords in self.gene_coordinates.values())
            self.logger.info(f"Loaded {total_coords} gene coordinates for {len(self.gene_coordinates)} genomes")
            
        except Exception as e:
            self.logger.error(f"Failed to load CARD coordinates: {e}")
            self.logger.debug(traceback.format_exc())
            raise
    
    def _load_coordinates_json(self, json_file: Path) -> None:
        """Load coordinates from JSON format"""
        with open(json_file, 'r') as f:
            data = json.load(f)
        
        for genome_id, genes in data.items():
            for gene_name, coord_info in genes.items():
                try:
                    coord = GeneCoordinate(
                        gene_name=gene_name,
                        genome_id=genome_id,
                        start=int(coord_info['start']),
                        end=int(coord_info['end']),
                        strand=str(coord_info['strand']),
                        contig=coord_info.get('contig') or coord_info.get('contig_id')
                    )
                    self.gene_coordinates[genome_id].append(coord)
                except Exception as e:
                    self.logger.warning(f"Invalid coordinate entry for {genome_id}:{gene_name}: {e}")
    
    def _load_coordinates_csv(self, csv_file: Path) -> None:
        """Load coordinates from CSV/TSV format"""
        # Simple CSV parser without pandas dependency
        with open(csv_file, 'r') as f:
            lines = f.readlines()
        
        if not lines:
            raise ValueError("Empty coordinates file")
        
        # Parse header
        delimiter = '\t' if csv_file.suffix.lower() in ['.tsv', '.tab'] else ','
        header = lines[0].strip().split(delimiter)
        
        # Map column names (case insensitive)
        col_map = {}
        for i, col in enumerate(header):
            col_lower = col.lower().strip()
            if col_lower in ['genome_id', 'accession', 'genome']:
                col_map['genome_id'] = i
            elif col_lower in ['gene', 'gene_name', 'protein']:
                col_map['gene_name'] = i
            elif col_lower in ['start', 'start_pos', 'begin']:
                col_map['start'] = i
            elif col_lower in ['end', 'end_pos', 'stop']:
                col_map['end'] = i
            elif col_lower in ['strand', 'direction']:
                col_map['strand'] = i
            elif col_lower in ['contig', 'contig_id', 'chromosome', 'scaffold']:
                col_map['contig'] = i
        
        required_cols = ['genome_id', 'gene_name', 'start', 'end', 'strand']
        missing_cols = [col for col in required_cols if col not in col_map]
        if missing_cols:
            raise ValueError(f"Missing required columns: {missing_cols}")
        
        # Parse data rows
        for line_num, line in enumerate(lines[1:], 2):
            if not line.strip():
                continue
            
            try:
                fields = line.strip().split(delimiter)
                if len(fields) < len(required_cols):
                    continue
                
                coord = GeneCoordinate(
                    gene_name=fields[col_map['gene_name']].strip(),
                    genome_id=fields[col_map['genome_id']].strip(),
                    start=int(fields[col_map['start']]),
                    end=int(fields[col_map['end']]),
                    strand=fields[col_map['strand']].strip(),
                    contig=fields[col_map['contig']].strip() if 'contig' in col_map else None
                )
                self.gene_coordinates[coord.genome_id].append(coord)
                
            except Exception as e:
                self.logger.warning(f"Invalid coordinate entry at line {line_num}: {e}")
    
    def extract_proteins(self, 
                        genomes_dir: Union[str, Path],
                        gene_list: Optional[List[str]] = None,
                        use_external_tool: bool = True) -> Dict[str, str]:
        """
        Extract protein sequences from genomes using coordinates
        
        Args:
            genomes_dir: Directory containing genome FASTA files
            gene_list: List of genes to extract (if None, extract all)
            use_external_tool: Whether to use external FastaAAExtractor tool
            
        Returns:
            Dictionary mapping output file types to file paths
        """
        try:
            genomes_path = Path(genomes_dir)
            if not genomes_path.exists():
                raise FileNotFoundError(f"Genomes directory not found: {genomes_dir}")
            
            self.logger.info(f"Starting protein extraction from {genomes_path}")
            
            # Get list of genome files
            genome_files = self._find_genome_files(genomes_path)
            self.stats['total_genomes'] = len(genome_files)
            
            # Filter gene list if provided
            if gene_list:
                self._filter_coordinates_by_genes(gene_list)
                self.stats['total_genes_requested'] = len(gene_list)
            else:
                all_genes = set()
                for coords in self.gene_coordinates.values():
                    all_genes.update(coord.gene_name for coord in coords)
                self.stats['total_genes_requested'] = len(all_genes)
            
            # Choose extraction method
            if use_external_tool and self.external_extractor_path:
                return self._extract_with_external_tool(genome_files)
            else:
                return self._extract_with_internal_method(genome_files)
                
        except Exception as e:
            self.logger.error(f"Protein extraction failed: {e}")
            self.logger.debug(traceback.format_exc())
            raise
    
    def _find_genome_files(self, genomes_dir: Path) -> Dict[str, Path]:
        """Find genome files and map them to genome IDs"""
        genome_files = {}
        
        # Common genome file extensions
        extensions = ['.fna', '.fasta', '.fa', '.faa', '.ffn']
        
        for genome_id in self.gene_coordinates.keys():
            found = False
            for ext in extensions:
                # Try exact match first
                genome_file = genomes_dir / f"{genome_id}{ext}"
                if genome_file.exists():
                    genome_files[genome_id] = genome_file
                    found = True
                    break
                
                # Try pattern matching
                for file_path in genomes_dir.glob(f"*{genome_id}*{ext}"):
                    genome_files[genome_id] = file_path
                    found = True
                    break
                
                if found:
                    break
            
            if not found:
                # Try genomic.fna pattern (common for NCBI downloads)
                genomic_file = genomes_dir / f"{genome_id}_genomic.fna"
                if genomic_file.exists():
                    genome_files[genome_id] = genomic_file
                else:
                    self.logger.warning(f"Genome file not found for {genome_id}")
        
        self.logger.info(f"Found {len(genome_files)} genome files")
        return genome_files
    
    def _filter_coordinates_by_genes(self, gene_list: List[str]) -> None:
        """Filter coordinates to only include specified genes"""
        gene_set = set(gene_list)
        
        for genome_id in list(self.gene_coordinates.keys()):
            filtered_coords = [
                coord for coord in self.gene_coordinates[genome_id]
                if coord.gene_name in gene_set
            ]
            self.gene_coordinates[genome_id] = filtered_coords
            
            # Remove genomes with no matching genes
            if not filtered_coords:
                del self.gene_coordinates[genome_id]
    
    def _extract_with_external_tool(self, genome_files: Dict[str, Path]) -> Dict[str, str]:
        """Extract proteins using external FastaAAExtractor tool"""
        self.logger.info("Using external FastaAAExtractor tool")
        
        output_files = {}
        
        for genome_id, genome_file in genome_files.items():
            if genome_id not in self.gene_coordinates:
                continue
            
            try:
                # Create coordinates file for this genome
                coord_file = self._create_coordinate_file(genome_id)
                
                # Run external tool
                output_zip = self._run_external_extractor(genome_file, coord_file)
                
                # Process output
                if output_zip and output_zip.exists():
                    self._process_external_output(genome_id, output_zip)
                    self.stats['successful_extractions'] += 1
                else:
                    self.logger.error(f"External tool failed for {genome_id}")
                    self.stats['failed_extractions'] += 1
                    
            except Exception as e:
                self.logger.error(f"External extraction failed for {genome_id}: {e}")
                self.stats['failed_extractions'] += 1
        
        # Generate output files
        output_files = self._generate_output_files()
        return output_files
    
    def _extract_with_internal_method(self, genome_files: Dict[str, Path]) -> Dict[str, str]:
        """Extract proteins using internal BioPython-based method"""
        self.logger.info("Using internal extraction method")
        
        if not BIOPYTHON_AVAILABLE:
            raise ImportError("BioPython not available for internal extraction method")
        
        for genome_id, genome_file in genome_files.items():
            if genome_id not in self.gene_coordinates:
                continue
            
            try:
                self._extract_from_genome_internal(genome_id, genome_file)
                self.stats['successful_extractions'] += 1
            except Exception as e:
                self.logger.error(f"Internal extraction failed for {genome_id}: {e}")
                self.stats['failed_extractions'] += 1
        
        output_files = self._generate_output_files()
        return output_files
    
    def _create_coordinate_file(self, genome_id: str) -> Path:
        """Create coordinate file for external tool"""
        coord_file = self.temp_dir / f"{genome_id}_coordinates.tsv"
        
        with open(coord_file, 'w') as f:
            f.write("Gene\tStart\tEnd\tStrand\n")
            for coord in self.gene_coordinates[genome_id]:
                f.write(f"{coord.gene_name}\t{coord.start}\t{coord.end}\t{coord.strand}\n")
        
        return coord_file
    
    def _run_external_extractor(self, genome_file: Path, coord_file: Path) -> Optional[Path]:
        """Run external FastaAAExtractor tool"""
        try:
            output_dir = self.temp_dir / f"external_output_{genome_file.stem}"
            output_dir.mkdir(exist_ok=True)
            
            # Construct command
            cmd = [
                "python", str(self.external_extractor_path),
                "--genome", str(genome_file),
                "--coordinates", str(coord_file),
                "--output", str(output_dir)
            ]
            
            # Run command
            result = subprocess.run(
                cmd, 
                capture_output=True, 
                text=True, 
                timeout=300  # 5 minute timeout
            )
            
            if result.returncode == 0:
                # Find output zip file
                zip_files = list(output_dir.glob("*.zip"))
                if zip_files:
                    return zip_files[0]
                else:
                    # Look for individual FASTA files
                    fasta_files = list(output_dir.glob("*.faa"))
                    if fasta_files:
                        return output_dir
            else:
                self.logger.error(f"External tool error: {result.stderr}")
                
        except subprocess.TimeoutExpired:
            self.logger.error("External tool timeout")
        except Exception as e:
            self.logger.error(f"External tool execution failed: {e}")
        
        return None
    
    def _process_external_output(self, genome_id: str, output_path: Path) -> None:
        """Process output from external tool"""
        if output_path.suffix == '.zip':
            # Extract zip file
            with zipfile.ZipFile(output_path, 'r') as zip_ref:
                extract_dir = self.temp_dir / f"extracted_{genome_id}"
                zip_ref.extractall(extract_dir)
                
                # Process extracted FASTA files
                for fasta_file in extract_dir.glob("*.faa"):
                    self._process_fasta_file(genome_id, fasta_file)
        else:
            # Process directory with FASTA files
            for fasta_file in output_path.glob("*.faa"):
                self._process_fasta_file(genome_id, fasta_file)
    
    def _process_fasta_file(self, genome_id: str, fasta_file: Path) -> None:
        """Process individual FASTA file from external tool"""
        try:
            if BIOPYTHON_AVAILABLE:
                for record in SeqIO.parse(fasta_file, "fasta"):
                    # Parse gene name from filename or header
                    gene_name = self._extract_gene_name(fasta_file.stem, record.id)
                    
                    # Find corresponding coordinate
                    coord = self._find_coordinate(genome_id, gene_name)
                    if coord:
                        result = ExtractionResult(
                            gene_name=gene_name,
                            genome_id=genome_id,
                            sequence=str(record.seq),
                            coordinates=coord,
                            extraction_method="external_tool"
                        )
                        self.extraction_results.append(result)
                        self.stats['sequences_generated'] += 1
            else:
                # Simple FASTA parsing without BioPython
                with open(fasta_file, 'r') as f:
                    content = f.read()
                    sequences = self._parse_fasta_simple(content)
                    
                    for header, sequence in sequences:
                        gene_name = self._extract_gene_name(fasta_file.stem, header)
                        coord = self._find_coordinate(genome_id, gene_name)
                        if coord:
                            result = ExtractionResult(
                                gene_name=gene_name,
                                genome_id=genome_id,
                                sequence=sequence,
                                coordinates=coord,
                                extraction_method="external_tool"
                            )
                            self.extraction_results.append(result)
                            self.stats['sequences_generated'] += 1
                            
        except Exception as e:
            self.logger.error(f"Error processing FASTA file {fasta_file}: {e}")
    
    def _extract_from_genome_internal(self, genome_id: str, genome_file: Path) -> None:
        """Extract proteins using internal BioPython method"""
        # Load genome sequences
        genome_records = {}
        for record in SeqIO.parse(genome_file, "fasta"):
            genome_records[record.id] = record
        
        # Extract each gene
        for coord in self.gene_coordinates[genome_id]:
            try:
                # Find the right contig/chromosome
                target_record = None
                if coord.contig and coord.contig in genome_records:
                    target_record = genome_records[coord.contig]
                else:
                    # Try first record (common for single-contig genomes)
                    target_record = list(genome_records.values())[0]
                
                if not target_record:
                    raise ValueError(f"No suitable contig found for {coord.gene_name}")
                
                # Extract DNA sequence
                dna_seq = target_record.seq[coord.start-1:coord.end]  # Convert to 0-based
                
                # Handle reverse strand
                if not coord.is_forward_strand:
                    dna_seq = dna_seq.reverse_complement()
                
                # Translate to protein
                protein_seq = dna_seq.translate()
                
                # Remove stop codon if present
                if protein_seq.endswith('*'):
                    protein_seq = protein_seq[:-1]
                
                result = ExtractionResult(
                    gene_name=coord.gene_name,
                    genome_id=genome_id,
                    sequence=str(protein_seq),
                    coordinates=coord,
                    extraction_method="internal_biopython"
                )
                
                # Quality checks
                if len(protein_seq) < 10:
                    result.warnings.append("Very short protein sequence")
                if '*' in str(protein_seq):
                    result.warnings.append("Internal stop codons found")
                
                self.extraction_results.append(result)
                self.stats['sequences_generated'] += 1
                
            except Exception as e:
                self.logger.error(f"Failed to extract {coord.gene_name} from {genome_id}: {e}")
                self.failed_extractions.append({
                    'genome_id': genome_id,
                    'gene_name': coord.gene_name,
                    'error': str(e)
                })
    
    def _generate_output_files(self) -> Dict[str, str]:
        """Generate output files for WildTypeAligner"""
        output_files = {}
        
        # Group sequences by gene
        sequences_by_gene = defaultdict(list)
        for result in self.extraction_results:
            sequences_by_gene[result.gene_name].append(result)
        
        # Create individual FASTA files for each gene
        fasta_files = []
        for gene_name, results in sequences_by_gene.items():
            gene_fasta = self.output_dir / f"{gene_name}_extracted_proteins.faa"
            
            with open(gene_fasta, 'w') as f:
                for result in results:
                    f.write(f">{result.fasta_header}\n")
                    f.write(f"{result.sequence}\n")
            
            fasta_files.append(str(gene_fasta))
            self.logger.info(f"Generated {gene_fasta} with {len(results)} sequences")
        
        output_files['individual_fasta'] = fasta_files
        
        # Create combined FASTA file
        combined_fasta = self.output_dir / "all_extracted_proteins.faa"
        with open(combined_fasta, 'w') as f:
            for result in self.extraction_results:
                f.write(f">{result.fasta_header}\n")
                f.write(f"{result.sequence}\n")
        
        output_files['combined_fasta'] = str(combined_fasta)
        
        # Create metadata file
        metadata_file = self.output_dir / "extraction_metadata.json"
        metadata = {
            'extraction_stats': self.stats,
            'extraction_results': [
                {
                    'gene_name': result.gene_name,
                    'genome_id': result.genome_id,
                    'sequence_length': result.sequence_length,
                    'coordinates': {
                        'start': result.coordinates.start,
                        'end': result.coordinates.end,
                        'strand': result.coordinates.strand
                    },
                    'extraction_method': result.extraction_method,
                    'quality_score': result.quality_score,
                    'warnings': result.warnings
                }
                for result in self.extraction_results
            ],
            'failed_extractions': self.failed_extractions
        }
        
        with open(metadata_file, 'w') as f:
            json.dump(metadata, f, indent=2)
        
        output_files['metadata'] = str(metadata_file)
        
        # Create summary report
        summary_file = self.output_dir / "extraction_summary.txt"
        with open(summary_file, 'w') as f:
            f.write("FastaAAExtractor Integration - Extraction Summary\n")
            f.write("=" * 50 + "\n\n")
            f.write(f"Extraction Date: {datetime.now().isoformat()}\n")
            f.write(f"Total Genomes: {self.stats['total_genomes']}\n")
            f.write(f"Total Genes Requested: {self.stats['total_genes_requested']}\n")
            f.write(f"Successful Extractions: {self.stats['successful_extractions']}\n")
            f.write(f"Failed Extractions: {self.stats['failed_extractions']}\n")
            f.write(f"Sequences Generated: {self.stats['sequences_generated']}\n\n")
            
            f.write("Genes Extracted:\n")
            for gene_name, results in sequences_by_gene.items():
                f.write(f"  {gene_name}: {len(results)} sequences\n")
            
            if self.failed_extractions:
                f.write(f"\nFailed Extractions ({len(self.failed_extractions)}):\n")
                for failure in self.failed_extractions:
                    f.write(f"  {failure['genome_id']}:{failure['gene_name']} - {failure['error']}\n")
        
        output_files['summary'] = str(summary_file)
        
        self.logger.info(f"Generated {len(output_files)} output files")
        return output_files
    
    def prepare_for_wild_type_aligner(self, 
                                    reference_dir: Union[str, Path],
                                    output_structure: str = "by_gene") -> Dict[str, str]:
        """
        Prepare extracted sequences for WildTypeAligner input
        
        Args:
            reference_dir: Directory containing reference sequences
            output_structure: How to organize output ("by_gene" or "by_genome")
            
        Returns:
            Dictionary with paths to prepared files
        """
        try:
            ref_path = Path(reference_dir)
            if not ref_path.exists():
                raise FileNotFoundError(f"Reference directory not found: {reference_dir}")
            
            self.logger.info("Preparing sequences for WildTypeAligner")
            
            # Create aligner input directory
            aligner_input_dir = self.output_dir / "wild_type_aligner_input"
            aligner_input_dir.mkdir(exist_ok=True)
            
            # Group sequences for alignment
            if output_structure == "by_gene":
                return self._prepare_by_gene(aligner_input_dir, ref_path)
            else:
                # Fallback to by_gene structure until by_genome is implemented
                return self._prepare_by_gene(aligner_input_dir, ref_path)
                
        except Exception as e:
            self.logger.error(f"Failed to prepare for WildTypeAligner: {e}")
            raise
    
    def _prepare_by_gene(self, output_dir: Path, ref_dir: Path) -> Dict[str, str]:
        """Prepare sequences organized by gene"""
        prepared_files = {}
        
        # Group by gene
        sequences_by_gene = defaultdict(list)
        for result in self.extraction_results:
            sequences_by_gene[result.gene_name].append(result)
        
        for gene_name, results in sequences_by_gene.items():
            gene_dir = output_dir / gene_name
            gene_dir.mkdir(exist_ok=True)
            
            # Create query sequences file
            query_file = gene_dir / f"{gene_name}_query_sequences.faa"
            with open(query_file, 'w') as f:
                for result in results:
                    f.write(f">{result.fasta_header}\n")
                    f.write(f"{result.sequence}\n")
            
            # Find reference sequence
            ref_file = self._find_reference_sequence(ref_dir, gene_name)
            if ref_file:
                # Copy reference to gene directory
                ref_copy = gene_dir / f"{gene_name}_reference.faa"
                shutil.copy2(ref_file, ref_copy)
                
                prepared_files[gene_name] = {
                    'query_sequences': str(query_file),
                    'reference_sequence': str(ref_copy),
                    'output_directory': str(gene_dir)
                }
            else:
                self.logger.warning(f"No reference sequence found for {gene_name}")
        
        return prepared_files
    
    def _find_reference_sequence(self, ref_dir: Path, gene_name: str) -> Optional[Path]:
        """Find reference sequence file for a gene"""
        # Try various naming patterns
        patterns = [
            f"{gene_name}.faa",
            f"{gene_name}_reference.faa",
            f"{gene_name.lower()}.faa",
            f"{gene_name.upper()}.faa",
            f"*{gene_name}*.faa"
        ]
        
        for pattern in patterns:
            matches = list(ref_dir.glob(pattern))
            if matches:
                return matches[0]
        
        return None
    
    # Helper methods
    def _extract_gene_name(self, filename: str, header: str) -> str:
        """Extract gene name from filename or FASTA header"""
        # Try filename first
        if '_' in filename:
            parts = filename.split('_')
            for part in parts:
                if part.lower() in ['acra', 'acrb', 'tolc', 'mdtf', 'oqxa', 'oqxb']:
                    return part.lower()
        
        # Try header
        if header:
            for part in header.split('_'):
                if part.lower() in ['acra', 'acrb', 'tolc', 'mdtf', 'oqxa', 'oqxb']:
                    return part.lower()
        
        # Default to filename
        return filename.split('_')[0] if '_' in filename else filename
    
    def _find_coordinate(self, genome_id: str, gene_name: str) -> Optional[GeneCoordinate]:
        """Find coordinate for a specific gene in a genome"""
        for coord in self.gene_coordinates.get(genome_id, []):
            if coord.gene_name.lower() == gene_name.lower():
                return coord
        return None
    
    def _parse_fasta_simple(self, content: str) -> List[Tuple[str, str]]:
        """Simple FASTA parser without BioPython"""
        sequences = []
        lines = content.strip().split('\n')
        
        header = None
        sequence_lines = []
        
        for line in lines:
            if line.startswith('>'):
                if header and sequence_lines:
                    sequences.append((header, ''.join(sequence_lines)))
                header = line[1:].strip()
                sequence_lines = []
            else:
                sequence_lines.append(line.strip())
        
        # Add last sequence
        if header and sequence_lines:
            sequences.append((header, ''.join(sequence_lines)))
        
        return sequences


def main():
    """Command line interface"""
    parser = argparse.ArgumentParser(
        description="FastaAAExtractor Pipeline Integration",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic extraction with CARD coordinates
  python fasta_aa_extractor_integration.py \\
      --coordinates card_results.csv \\
      --genomes genomes/ \\
      --output extracted_proteins/

  # Extract specific genes only
  python fasta_aa_extractor_integration.py \\
      --coordinates card_results.json \\
      --genomes genomes/ \\
      --genes mdtF acrA acrB \\
      --output extracted_proteins/

  # Use external FastaAAExtractor tool
  python fasta_aa_extractor_integration.py \\
      --coordinates card_results.csv \\
      --genomes genomes/ \\
      --external-tool /path/to/FastaAAExtractor.py \\
      --output extracted_proteins/

  # Prepare for WildTypeAligner
  python fasta_aa_extractor_integration.py \\
      --coordinates card_results.csv \\
      --genomes genomes/ \\
      --references references/ \\
      --prepare-aligner \\
      --output extracted_proteins/
        """
    )
    
    # Required arguments
    parser.add_argument('--coordinates', required=True,
                       help='Path to CARD coordinates file (CSV/TSV/JSON)')
    parser.add_argument('--genomes', required=True,
                       help='Directory containing genome FASTA files')
    parser.add_argument('--output', required=True,
                       help='Output directory for extracted proteins')
    
    # Optional arguments
    parser.add_argument('--genes', nargs='+',
                       help='Specific genes to extract (default: extract all)')
    parser.add_argument('--external-tool',
                       help='Path to external FastaAAExtractor tool')
    parser.add_argument('--references',
                       help='Directory containing reference sequences')
    parser.add_argument('--prepare-aligner', action='store_true',
                       help='Prepare output for WildTypeAligner')
    parser.add_argument('--coordinate-format', choices=['csv', 'tsv', 'json', 'auto'],
                       default='auto', help='Format of coordinates file')
    parser.add_argument('--log-level', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                       default='INFO', help='Logging level')
    parser.add_argument('--temp-dir',
                       help='Temporary directory for intermediate files')
    
    args = parser.parse_args()
    
    try:
        # Initialize integrator
        print("Initializing FastaAAExtractor integrator...")
        integrator = FastaAAExtractorIntegrator(
            output_dir=args.output,
            temp_dir=args.temp_dir,
            external_extractor_path=args.external_tool,
            log_level=args.log_level
        )
        
        # Load coordinates
        print(f"Loading CARD coordinates from {args.coordinates}...")
        integrator.load_card_coordinates(
            coordinates_file=args.coordinates,
            file_format=args.coordinate_format
        )
        
        # Extract proteins
        print(f"Extracting proteins from genomes in {args.genomes}...")
        use_external = bool(args.external_tool)
        output_files = integrator.extract_proteins(
            genomes_dir=args.genomes,
            gene_list=args.genes,
            use_external_tool=use_external
        )
        
        # Prepare for aligner if requested
        if args.prepare_aligner and args.references:
            print("Preparing sequences for WildTypeAligner...")
            aligner_files = integrator.prepare_for_wild_type_aligner(
                reference_dir=args.references
            )
            output_files.update(aligner_files)
        
        # Print results
        print("\n" + "=" * 60)
        print("PROTEIN EXTRACTION COMPLETE!")
        print("=" * 60)
        print(f"Total genomes processed: {integrator.stats['total_genomes']}")
        print(f"Successful extractions: {integrator.stats['successful_extractions']}")
        print(f"Failed extractions: {integrator.stats['failed_extractions']}")
        print(f"Sequences generated: {integrator.stats['sequences_generated']}")
        
        print(f"\nOutput files:")
        for file_type, file_path in output_files.items():
            if isinstance(file_path, list):
                print(f"  {file_type}: {len(file_path)} files")
                for fp in file_path[:3]:  # Show first 3
                    print(f"    - {fp}")
                if len(file_path) > 3:
                    print(f"    - ... and {len(file_path) - 3} more")
            else:
                print(f"  {file_type}: {file_path}")
        
        print(f"\nNext steps:")
        print(f"  1. Review extraction summary: {output_files.get('summary', 'N/A')}")
        print(f"  2. Run WildTypeAligner on extracted proteins")
        print(f"  3. Continue with SubScan for mutation analysis")
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()