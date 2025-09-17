from .configuration_manager import config_manager
#!/usr/bin/env python3
"""
SEPI-Enhanced WildTypeAligner - Intelligent reference acquisition and pairwise alignment

This module integrates SEPI 2.0 for intelligent reference protein acquisition and
EMBOSS WATER-style pairwise alignment for mutation analysis.

Workflow:
1. Analyze extracted proteins to determine genus/species
2. Use SEPI 2.0 to fetch appropriate reference sequences  
3. Perform pairwise alignment using EMBOSS WATER
4. Generate alignment reports for downstream mutation analysis

Author: GenomeAMRAnalyzer Pipeline
Version: 1.0.0
"""

import os
import sys
import argparse
import logging
import subprocess
import json
import csv
import hashlib
from pathlib import Path
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass
import tempfile
import re
from collections import Counter
import time

# Add SEPI path for import
sepi_path = Path(__file__).parent.parent / "MetaDataHarvester" / "sepi2.0"
sys.path.append(str(sepi_path))

try:
    from Bio import SeqIO, Entrez
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


@dataclass
class SpeciesInfo:
    """Information about detected species"""
    genus: str
    species: str
    strain: str = ""
    confidence: float = 0.0
    method: str = "header_parsing"


@dataclass
class ReferenceProtein:
    """Reference protein information"""
    protein_name: str
    accession: str
    organism: str
    sequence: str
    source_strain: str = ""
    assembly_level: str = ""
    ncbi_url: str = ""


@dataclass
class AlignmentResult:
    """Pairwise alignment result"""
    query_id: str
    reference_id: str
    query_sequence: str
    reference_sequence: str
    alignment_score: float
    identity_percent: float
    similarity_percent: float
    gaps: int
    alignment_length: int
    query_coverage: float
    reference_coverage: float
    alignment_file: str = ""


@dataclass
class WildTypeAlignerConfig:
    """Configuration for SEPI-Enhanced WildTypeAligner"""
    input_dir: str
    output_dir: str
    email: str
    sepi_path: str
    target_genes: List[str]
    reference_organism: Optional[str] = None
    assembly_level: str = "complete_genome"
    max_references: int = 3
    use_water: bool = True
    water_executable: str = "water"
    gap_open: float = 10.0
    gap_extend: float = 0.5


class SEPIEnhancedWildTypeAligner:
    """
    SEPI-Enhanced WildTypeAligner with intelligent reference acquisition
    """

    def __init__(self, config: WildTypeAlignerConfig):
        """Initialize the SEPI-enhanced aligner"""
        self.config = config
        self.setup_logging()
        self.setup_directories()
        
        # Validate dependencies
        self.validate_dependencies()
        
        # Species detection cache
        self.species_cache = {}
        
        # Statistics
        self.stats = {
            'proteins_processed': 0,
            'species_detected': 0,
            'references_fetched': 0,
            'alignments_performed': 0,
            'failed_alignments': 0
        }

    def setup_logging(self):
        """Setup logging configuration"""
        log_dir = Path(self.config.output_dir) / "logs"
        log_dir.mkdir(parents=True, exist_ok=True)

        # Add file handler
        file_handler = logging.FileHandler(log_dir / "wildtype_aligner.log")
        file_handler.setFormatter(logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        ))
        
        self.logger = logging.getLogger('SEPIWildTypeAligner')
        self.logger.addHandler(file_handler)
        self.logger.setLevel(logging.INFO)

    def setup_directories(self):
        """Create necessary output directories"""
        output_path = Path(self.config.output_dir)
        
        # Create subdirectories
        subdirs = [
            "references", "alignments", "reports", "temp", "logs"
        ]
        
        for subdir in subdirs:
            (output_path / subdir).mkdir(parents=True, exist_ok=True)

        self.logger.info(f"Output directory: {output_path.absolute()}")

    def validate_dependencies(self):
        """Validate required dependencies"""
        dependencies_ok = True
        
        # Check BioPython
        if not BIOPYTHON_AVAILABLE:
            self.logger.error("BioPython not available - required for sequence processing")
            dependencies_ok = False
        
        # Check SEPI
        sepi_script = Path(self.config.sepi_path) / "sepi.py"
        if not sepi_script.exists():
            self.logger.error(f"SEPI script not found: {sepi_script}")
            dependencies_ok = False
        
        # Check WATER (optional - we'll create fallback)
        if self.config.use_water:
            try:
                result = subprocess.run(
                    [self.config.water_executable, "-help"],
                    capture_output=True, text=True, timeout=10
                )
                if result.returncode != 0:
                    self.logger.warning("EMBOSS WATER not available - will use BioPython alignment")
                    self.config.use_water = False
            except (FileNotFoundError, subprocess.TimeoutExpired):
                self.logger.warning("EMBOSS WATER not available - will use BioPython alignment")
                self.config.use_water = False
        
        if not dependencies_ok:
            raise RuntimeError("Required dependencies not available")

    def run_alignment_pipeline(self) -> bool:
        """
        Main method to run the SEPI-enhanced alignment pipeline
        
        Returns:
            bool: True if successful, False otherwise
        """
        try:
            self.logger.info("Starting SEPI-Enhanced WildType Alignment Pipeline")
            
            # Find protein files
            protein_files = self._find_protein_files()
            if not protein_files:
                self.logger.error("No protein FASTA files found")
                return False
            
            self.logger.info(f"Found {len(protein_files)} protein files to process")
            
            # Process each protein file
            all_alignments = []
            for protein_file in protein_files:
                alignments = self._process_protein_file(protein_file)
                all_alignments.extend(alignments)
            
            # Generate summary reports
            self._generate_summary_reports(all_alignments)
            
            # Print statistics
            self._print_statistics()
            
            return True

        except Exception as e:
            self.logger.error(f"Error in alignment pipeline: {e}")
            return False

    def _find_protein_files(self) -> List[Path]:
        """Find protein FASTA files to process"""
        input_path = Path(self.config.input_dir)
        
        if not input_path.exists():
            raise FileNotFoundError(f"Input directory not found: {input_path}")
        
        # Look for FASTA files
        protein_files = []
        patterns = ['*.faa', '*.fa', '*.fasta']
        
        if input_path.is_file():
            # Single file
            protein_files = [input_path]
        else:
            # Directory - search recursively
            for pattern in patterns:
                protein_files.extend(input_path.rglob(pattern))
        
        return sorted(protein_files)

    def _process_protein_file(self, protein_file: Path) -> List[AlignmentResult]:
        """Process a single protein file"""
        self.logger.info(f"Processing protein file: {protein_file.name}")
        
        try:
            # Parse protein sequences
            sequences = list(SeqIO.parse(protein_file, "fasta"))
            if not sequences:
                self.logger.warning(f"No sequences found in {protein_file}")
                return []
            
            alignments = []
            
            # Group sequences by gene name
            gene_groups = self._group_sequences_by_gene(sequences)
            
            for gene_name, gene_sequences in gene_groups.items():
                if self.config.target_genes and gene_name not in self.config.target_genes:
                    continue
                
                self.logger.info(f"Processing gene: {gene_name} ({len(gene_sequences)} sequences)")
                
                # Detect species from sequences
                species_info = self._detect_species(gene_sequences)
                
                # Fetch reference sequences using SEPI
                references = self._fetch_references_with_sepi(gene_name, species_info)
                
                if not references:
                    self.logger.warning(f"No references found for {gene_name}")
                    continue
                
                # Perform alignments
                for seq_record in gene_sequences:
                    for reference in references:
                        alignment = self._perform_alignment(seq_record, reference)
                        if alignment:
                            alignments.append(alignment)
                
            return alignments
            
        except Exception as e:
            self.logger.error(f"Error processing {protein_file}: {e}")
            return []

    def _group_sequences_by_gene(self, sequences: List[SeqRecord]) -> Dict[str, List[SeqRecord]]:
        """Group sequences by gene name"""
        gene_groups = {}
        
        for seq_record in sequences:
            gene_name = self._extract_gene_name(seq_record.id)
            
            if gene_name not in gene_groups:
                gene_groups[gene_name] = []
            gene_groups[gene_name].append(seq_record)
        
        return gene_groups

    def _extract_gene_name(self, sequence_id: str) -> str:
        """Extract gene name from sequence ID"""
        # Handle our naming convention: KWV17775.1_genome_acrA_1000_2000_+
        parts = sequence_id.split('_')
        
        # Look for gene names in parts
        for part in parts:
            part_lower = part.lower()
            if part_lower in ['acra', 'acrb', 'tolc', 'mara', 'marr', 'rama', 'ramr', 'soxs', 'rob', 'envr']:
                return part_lower
        
        # Fallback: use second-to-last part or "unknown"
        if len(parts) >= 2:
            return parts[-3] if len(parts) >= 3 else parts[0]
        
        return "unknown"

    def _detect_species(self, sequences: List[SeqRecord]) -> SpeciesInfo:
        """Detect species from sequence headers and content"""
        
        # Cache key for species detection
        seq_ids = [seq.id for seq in sequences[:3] if seq.id]  # Use first 3 sequences, filter None
        cache_key = hashlib.md5('_'.join(seq_ids).encode()).hexdigest()
        
        if cache_key in self.species_cache:
            return self.species_cache[cache_key]
        
        # Method 1: Parse from sequence headers
        organisms = []
        for seq_record in sequences[:5]:  # Check first 5 sequences
            header_organism = self._parse_organism_from_header(seq_record.id)
            if header_organism:
                organisms.append(header_organism)
        
        if organisms:
            # Find most common organism
            organism_counts = Counter(organisms)
            most_common = organism_counts.most_common(1)[0]
            organism_parts = most_common[0].split()
            
            species_info = SpeciesInfo(
                genus=organism_parts[0] if organism_parts else "Unknown",
                species=organism_parts[1] if len(organism_parts) > 1 else "unknown",
                strain=organism_parts[2] if len(organism_parts) > 2 else "",
                confidence=most_common[1] / len(organisms),
                method="header_parsing"
            )
        else:
            # Fallback: Use config organism or default
            if self.config.reference_organism:
                parts = self.config.reference_organism.split()
                species_info = SpeciesInfo(
                    genus=parts[0] if parts else "Escherichia",
                    species=parts[1] if len(parts) > 1 else "coli",
                    confidence=0.5,
                    method="config_fallback"
                )
            else:
                species_info = SpeciesInfo(
                    genus="Escherichia",
                    species="coli",
                    confidence=0.3,
                    method="default_fallback"
                )
        
        # Cache the result
        self.species_cache[cache_key] = species_info
        self.stats['species_detected'] += 1
        
        self.logger.info(f"Detected species: {species_info.genus} {species_info.species} "
                        f"(confidence: {species_info.confidence:.2f}, method: {species_info.method})")
        
        return species_info

    def _parse_organism_from_header(self, header: str) -> Optional[str]:
        """Parse organism information from sequence header"""
        header_lower = header.lower()
        
        # Look for genome names in our naming convention
        if '_genome' in header_lower:
            # Extract accession part (e.g., KWV17775.1 from KWV17775.1_genome_acrA_1000_2000_+)
            accession = header.split('_')[0]
            
            # Map accession prefixes to organisms (this could be enhanced)
            prefix_mappings = {
                'KWV': 'Escherichia coli',
                'APQ': 'Klebsiella pneumoniae', 
                'AQU': 'Enterobacter cloacae',
                # Add more mappings as needed
            }
            
            for prefix, organism in prefix_mappings.items():
                if accession.startswith(prefix):
                    return organism
        
        return None

    def _fetch_references_with_sepi(self, gene_name: str, species_info: SpeciesInfo) -> List[ReferenceProtein]:
        """Use SEPI 2.0 to fetch reference sequences"""
        try:
            self.logger.info(f"Fetching references for {gene_name} from {species_info.genus} {species_info.species}")
            
            # Prepare SEPI command
            organism = f"{species_info.genus} {species_info.species}"
            output_name = f"{gene_name}_{species_info.genus}_{species_info.species}_refs"
            output_path = Path(self.config.output_dir) / "references" / output_name
            
            sepi_cmd = [
                sys.executable, 
                str(Path(self.config.sepi_path) / "sepi.py"),
                "--organism", organism,
                "--proteins", gene_name,
                "--assembly_level", self.config.assembly_level,
                "--output", str(output_path),
                "--email", self.config.email,
                "--multi_fasta"
            ]
            
            # Run SEPI
            self.logger.debug(f"Running SEPI command: {' '.join(sepi_cmd)}")
            
            result = subprocess.run(
                sepi_cmd,
                capture_output=True,
                text=True,
                timeout=300,  # 5 minute timeout
                cwd=Path(self.config.sepi_path)
            )
            
            if result.returncode == 0:
                # Parse SEPI results
                references = self._parse_sepi_results(output_path, gene_name)
                self.stats['references_fetched'] += len(references)
                self.logger.info(f"SEPI fetched {len(references)} references for {gene_name}")
                return references[:self.config.max_references]
            else:
                self.logger.error(f"SEPI failed: {result.stderr}")
                return []
                
        except Exception as e:
            self.logger.error(f"Error running SEPI for {gene_name}: {e}")
            return []

    def _parse_sepi_results(self, output_path: Path, gene_name: str) -> List[ReferenceProtein]:
        """Parse SEPI output files to extract reference proteins"""
        references = []
        
        try:
            # Look for multi-FASTA file
            multi_fasta = output_path.with_suffix('.fasta')
            if multi_fasta.exists():
                for record in SeqIO.parse(multi_fasta, "fasta"):
                    ref_protein = ReferenceProtein(
                        protein_name=gene_name,
                        accession=record.id,
                        organism=self._extract_organism_from_header(record.description),
                        sequence=str(record.seq),
                        source_strain=self._extract_strain_from_header(record.description)
                    )
                    references.append(ref_protein)
            
            # Also check for CSV metadata
            csv_file = Path(str(output_path) + "_accessions.csv")
            if csv_file.exists():
                self._enrich_references_with_metadata(references, csv_file)
            
        except Exception as e:
            self.logger.error(f"Error parsing SEPI results: {e}")
        
        return references

    def _extract_organism_from_header(self, header: str) -> str:
        """Extract organism from FASTA header"""
        # This is a simplified extraction - could be enhanced
        if '[' in header and ']' in header:
            start = header.find('[')
            end = header.find(']', start)
            return header[start+1:end]
        return "Unknown"

    def _extract_strain_from_header(self, header: str) -> str:
        """Extract strain information from FASTA header"""
        # Look for strain indicators
        strain_patterns = [r'strain\s+([^\s\]]+)', r'str\.\s+([^\s\]]+)']
        for pattern in strain_patterns:
            match = re.search(pattern, header, re.IGNORECASE)
            if match:
                return match.group(1)
        return ""

    def _enrich_references_with_metadata(self, references: List[ReferenceProtein], csv_file: Path):
        """Enrich reference proteins with metadata from SEPI CSV"""
        try:
            with open(csv_file, 'r') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    accession = row.get('Accession Number', '')
                    for ref in references:
                        if accession in ref.accession:
                            ref.source_strain = row.get('Source Strain', ref.source_strain)
                            ref.ncbi_url = row.get('NCBI URL', '')
                            break
        except Exception as e:
            self.logger.warning(f"Could not enrich metadata: {e}")

    def _perform_alignment(self, query: SeqRecord, reference: ReferenceProtein) -> Optional[AlignmentResult]:
        """Perform pairwise alignment between query and reference"""
        try:
            if self.config.use_water:
                return self._perform_water_alignment(query, reference)
            else:
                return self._perform_biopython_alignment(query, reference)
        except Exception as e:
            self.logger.error(f"Alignment failed for {query.id} vs {reference.accession}: {e}")
            self.stats['failed_alignments'] += 1
            return None

    def _perform_water_alignment(self, query: SeqRecord, reference: ReferenceProtein) -> Optional[AlignmentResult]:
        """Perform alignment using EMBOSS WATER"""
        try:
            # Create temporary files
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as query_file:
                SeqIO.write(query, query_file, "fasta")
                query_path = query_file.name
            
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as ref_file:
                ref_record = SeqRecord(Seq(reference.sequence), id=reference.accession)
                SeqIO.write(ref_record, ref_file, "fasta")
                ref_path = ref_file.name
            
            # Output alignment file
            alignment_dir = Path(self.config.output_dir) / "alignments"
            alignment_file = alignment_dir / f"{query.id}_vs_{reference.accession}.water"
            
            # Run WATER
            water_cmd = [
                self.config.water_executable,
                "-asequence", query_path,
                "-bsequence", ref_path,
                "-gapopen", str(self.config.gap_open),
                "-gapextend", str(self.config.gap_extend),
                "-outfile", str(alignment_file),
                "-aformat", "srspair"  # SubScan-compatible format
            ]
            
            result = subprocess.run(water_cmd, capture_output=True, text=True, timeout=60)
            
            # Clean up temp files
            os.unlink(query_path)
            os.unlink(ref_path)
            
            if result.returncode == 0:
                # Parse WATER output
                alignment_result = self._parse_water_output(alignment_file, query.id, reference.accession, str(query.seq), reference.sequence)
                self.stats['alignments_performed'] += 1
                return alignment_result
            else:
                self.logger.error(f"WATER failed: {result.stderr}")
                return None
                
        except Exception as e:
            self.logger.error(f"WATER alignment error: {e}")
            return None

    def _perform_biopython_alignment(self, query: SeqRecord, reference: ReferenceProtein) -> Optional[AlignmentResult]:
        """Perform alignment using BioPython (fallback method)"""
        try:
            from Bio import pairwise2
            from Bio.pairwise2 import format_alignment
            
            # Perform global alignment
            alignments = pairwise2.align.globalxx(str(query.seq), reference.sequence)
            
            if alignments:
                best_alignment = alignments[0]
                
                # Calculate statistics
                aligned_query = best_alignment[0]
                aligned_ref = best_alignment[1]
                score = best_alignment[2]
                
                # Calculate identity and similarity
                matches = sum(1 for a, b in zip(aligned_query, aligned_ref) if a == b and a != '-')
                gaps = aligned_query.count('-') + aligned_ref.count('-')
                alignment_length = len(aligned_query)
                
                identity_percent = (matches / alignment_length) * 100 if alignment_length > 0 else 0
                query_coverage = (len(aligned_query.replace('-', '')) / len(query.seq)) * 100
                ref_coverage = (len(aligned_ref.replace('-', '')) / len(reference.sequence)) * 100
                
                # Save alignment to file
                alignment_dir = Path(self.config.output_dir) / "alignments"
                alignment_file = alignment_dir / f"{query.id}_vs_{reference.accession}.bio"
                
                with open(alignment_file, 'w') as f:
                    f.write(format_alignment(*best_alignment))
                
                alignment_result = AlignmentResult(
                    query_id=query.id,
                    reference_id=reference.accession,
                    query_sequence=str(query.seq),
                    reference_sequence=reference.sequence,
                    alignment_score=score,
                    identity_percent=identity_percent,
                    similarity_percent=identity_percent,  # Simplified
                    gaps=gaps,
                    alignment_length=alignment_length,
                    query_coverage=query_coverage,
                    reference_coverage=ref_coverage,
                    alignment_file=str(alignment_file)
                )
                
                self.stats['alignments_performed'] += 1
                return alignment_result
            
        except ImportError:
            self.logger.error("BioPython pairwise2 not available")
        except Exception as e:
            self.logger.error(f"BioPython alignment error: {e}")
        
        return None

    def _parse_water_output(self, alignment_file: Path, query_id: str, ref_id: str, query_seq: str, ref_seq: str) -> AlignmentResult:
        """Parse EMBOSS WATER output file"""
        # This is a simplified parser - could be enhanced for more detailed parsing
        try:
            with open(alignment_file, 'r') as f:
                content = f.read()
            
            # Extract basic statistics (simplified)
            score = 0.0
            identity = 0.0
            similarity = 0.0
            gaps = 0
            
            # Look for statistics in WATER output
            for line in content.split('\n'):
                if 'Score:' in line:
                    score = float(re.search(r'Score:\s*([\d.]+)', line).group(1))
                elif 'Identity:' in line:
                    match = re.search(r'Identity:\s*(\d+)/(\d+)\s*\(([\d.]+)%\)', line)
                    if match:
                        identity = float(match.group(3))
                elif 'Similarity:' in line:
                    match = re.search(r'Similarity:\s*(\d+)/(\d+)\s*\(([\d.]+)%\)', line)
                    if match:
                        similarity = float(match.group(3))
                elif 'Gaps:' in line:
                    match = re.search(r'Gaps:\s*(\d+)', line)
                    if match:
                        gaps = int(match.group(1))
            
            return AlignmentResult(
                query_id=query_id,
                reference_id=ref_id,
                query_sequence=query_seq,
                reference_sequence=ref_seq,
                alignment_score=score,
                identity_percent=identity,
                similarity_percent=similarity,
                gaps=gaps,
                alignment_length=len(query_seq),  # Simplified
                query_coverage=100.0,  # Simplified
                reference_coverage=100.0,  # Simplified
                alignment_file=str(alignment_file)
            )
            
        except Exception as e:
            self.logger.error(f"Error parsing WATER output: {e}")
            # Return basic result
            return AlignmentResult(
                query_id=query_id,
                reference_id=ref_id,
                query_sequence=query_seq,
                reference_sequence=ref_seq,
                alignment_score=0.0,
                identity_percent=0.0,
                similarity_percent=0.0,
                gaps=0,
                alignment_length=len(query_seq),
                query_coverage=100.0,
                reference_coverage=100.0,
                alignment_file=str(alignment_file)
            )

    def _generate_summary_reports(self, alignments: List[AlignmentResult]):
        """Generate summary reports"""
        try:
            reports_dir = Path(self.config.output_dir) / "reports"
            
            # Generate CSV report
            csv_file = reports_dir / "alignment_summary.csv"
            with open(csv_file, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow([
                    'Query ID', 'Reference ID', 'Alignment Score', 'Identity %', 
                    'Similarity %', 'Gaps', 'Query Coverage %', 'Reference Coverage %',
                    'Alignment File'
                ])
                
                for alignment in alignments:
                    writer.writerow([
                        alignment.query_id,
                        alignment.reference_id,
                        alignment.alignment_score,
                        f"{alignment.identity_percent:.2f}",
                        f"{alignment.similarity_percent:.2f}",
                        alignment.gaps,
                        f"{alignment.query_coverage:.2f}",
                        f"{alignment.reference_coverage:.2f}",
                        alignment.alignment_file
                    ])
            
            # Generate JSON report
            json_file = reports_dir / "alignment_summary.json"
            alignment_data = []
            for alignment in alignments:
                alignment_data.append({
                    'query_id': alignment.query_id,
                    'reference_id': alignment.reference_id,
                    'alignment_score': alignment.alignment_score,
                    'identity_percent': alignment.identity_percent,
                    'similarity_percent': alignment.similarity_percent,
                    'gaps': alignment.gaps,
                    'query_coverage': alignment.query_coverage,
                    'reference_coverage': alignment.reference_coverage,
                    'alignment_file': alignment.alignment_file
                })
            
            with open(json_file, 'w') as f:
                json.dump({
                    'summary': self.stats,
                    'alignments': alignment_data
                }, f, indent=2)
            
            self.logger.info(f"Generated summary reports: {csv_file}, {json_file}")
            
        except Exception as e:
            self.logger.error(f"Error generating reports: {e}")

    def _print_statistics(self):
        """Print alignment statistics"""
        self.logger.info("=" * 60)
        self.logger.info("SEPI-ENHANCED WILDTYPE ALIGNMENT STATISTICS")
        self.logger.info("=" * 60)
        self.logger.info(f"Proteins processed: {self.stats['proteins_processed']}")
        self.logger.info(f"Species detected: {self.stats['species_detected']}")
        self.logger.info(f"References fetched: {self.stats['references_fetched']}")
        self.logger.info(f"Alignments performed: {self.stats['alignments_performed']}")
        self.logger.info(f"Failed alignments: {self.stats['failed_alignments']}")

        success_rate = (self.stats['alignments_performed'] / max(1, self.stats['alignments_performed'] + self.stats['failed_alignments']) * 100)
        self.logger.info(f"Alignment success rate: {success_rate:.1f}%")


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description="SEPI-Enhanced WildType Aligner with intelligent reference acquisition",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic alignment with automatic species detection
  python sepi_wildtype_aligner.py --input proteins/ --output alignments/ --email your@email.com

  # Specify target organism and genes
  python sepi_wildtype_aligner.py --input proteins/ --output alignments/ --email your@email.com --organism "Escherichia coli" --genes acrA acrB

  # Use high-quality references only
  python sepi_wildtype_aligner.py --input proteins/ --output alignments/ --email your@email.com --assembly-level complete_genome
        """
    )

    parser.add_argument(
        '--input',
        required=True,
        help='Input directory containing protein FASTA files'
    )

    parser.add_argument(
        '--output',
        required=True,
        help='Output directory for alignments and reports'
    )

    parser.add_argument(
        '--email',
        required=True,
        help='Email address for NCBI API (required)'
    )

    parser.add_argument(
        '--sepi-path',
        default='MetaDataHarvester/sepi2.0',
        help='Path to SEPI 2.0 installation'
    )

    parser.add_argument(
        '--genes',
        nargs='+',
        help='Target genes to process (e.g., acrA acrB tolC)'
    )

    parser.add_argument(
        '--organism',
        help='Reference organism for SEPI (auto-detected if not specified)'
    )

    parser.add_argument(
        '--assembly-level',
        choices=['complete_genome', 'chromosome', 'scaffold', 'contig'],
        default='complete_genome',
        help='Assembly level filter for references'
    )

    parser.add_argument(
        '--max-references',
        type=int,
        default=3,
        help='Maximum number of reference sequences per gene'
    )

    parser.add_argument(
        '--no-water',
        action='store_true',
        help='Disable EMBOSS WATER and use BioPython alignment'
    )

    args = parser.parse_args()

    # Create configuration
    config = WildTypeAlignerConfig(
        input_dir=args.input,
        output_dir=args.output,
        email=args.email,
        sepi_path=args.sepi_path,
        target_genes=args.genes or [],
        reference_organism=args.organism,
        assembly_level=args.assembly_level,
        max_references=args.max_references,
        use_water=not args.no_water
    )

    # Run aligner
    aligner = SEPIEnhancedWildTypeAligner(config)
    success = aligner.run_alignment_pipeline()

    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()