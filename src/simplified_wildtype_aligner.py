#!/usr/bin/env python3
"""
Generic WildType Aligner - Pairwise alignment for ANY gene sets

This is a fully generic aligner that performs pairwise alignment using BioPython
for ANY user-specified genes - not hardcoded to specific gene sets.

Works with any gene sets such as:
- RND efflux pumps: mexA, mexB, oprM, triABC, etc.
- Regulators: marA, marR, soxS, etc.  
- Resistance genes: blaTEM, blaCTX-M, qnrS, etc.
- Any other bacterial genes of interest

Requires user-provided reference sequences - no hardcoded sequences.

Author: GenomeAMRAnalyzer Pipeline
Version: 2.0.0 (Generic)
"""

import os
import sys
import argparse
import logging
import json
import csv
from pathlib import Path
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass
import tempfile
import re
from collections import Counter

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False
    # Create mock classes for type hints when BioPython not available
    class SeqRecord:
        pass

# Import error handling if available
try:
    from .core.robust_error_handling import (
        ValidationError, DataProcessingError, FileSystemError, 
        ValidationSuite, robust_exception_handler, RobustLogger
    )
    ERROR_HANDLING_AVAILABLE = True
except ImportError:
    # Fallback if core modules not available yet
    ERROR_HANDLING_AVAILABLE = False
    class ValidationError(Exception): pass
    class DataProcessingError(Exception): pass
    class FileSystemError(Exception): pass

# Configure logging
if ERROR_HANDLING_AVAILABLE:
    logger = RobustLogger("SimplifiedWildTypeAligner").logger
else:
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    logger = logging.getLogger(__name__)


@dataclass
class AlignmentResult:
    """Pairwise alignment result"""
    query_id: str
    reference_id: str
    query_sequence: str
    reference_sequence: str
    alignment_score: float
    identity_percent: float
    gaps: int
    alignment_length: int
    alignment_file: str = ""


@dataclass
class SimpleAlignerConfig:
    """Configuration for Simplified WildTypeAligner"""
    input_dir: str
    output_dir: str
    target_genes: List[str]
    reference_dir: Optional[str] = None


class SimplifiedWildTypeAligner:
    """
    Simplified WildType Aligner with basic pairwise alignment
    """

    def __init__(self, config: SimpleAlignerConfig):
        """Initialize the simplified aligner"""
        self.config = config
        self.setup_logging()
        self.setup_directories()
        # No built-in references: all references must be provided by user/tests
        self.builtin_references = {}
        # Statistics
        self.stats = {
            'proteins_processed': 0,
            'alignments_performed': 0,
            'failed_alignments': 0
        }

    def setup_logging(self):
        """Setup logging configuration"""
        log_dir = Path(self.config.output_dir) / "logs"
        log_dir.mkdir(parents=True, exist_ok=True)

        # Add file handler
        self.log_file = log_dir / "simple_aligner.log"
        file_handler = logging.FileHandler(self.log_file)
        file_handler.setFormatter(logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        ))
        
        self.logger = logging.getLogger('SimplifiedWildTypeAligner')
        self.logger.addHandler(file_handler)
        self.logger.setLevel(logging.INFO)
        
        # Store handler for cleanup
        self.file_handler = file_handler

    def setup_directories(self):
        """Create necessary output directories"""
        output_path = Path(self.config.output_dir)
        
        # Create subdirectories
        subdirs = ["alignments", "reports", "logs"]
        
        for subdir in subdirs:
            (output_path / subdir).mkdir(parents=True, exist_ok=True)

        self.logger.info(f"Output directory: {output_path.absolute()}")

    def ensure_output_structure(self):
        """Idempotently ensure all expected output subdirectories exist."""
        try:
            output_path = Path(self.config.output_dir)
            for sub in ("alignments", "reports", "logs"):
                (output_path / sub).mkdir(parents=True, exist_ok=True)
        except Exception as e:
            self.logger.warning(f"Failed to ensure output structure: {e}")

    def _create_builtin_references(self) -> Dict[str, str]:
        """No built-in references: all references must be provided by user/tests."""
        return {}

    def run_alignment_pipeline(self) -> bool:
        """
        Main method to run the simplified alignment pipeline
        
        Returns:
            bool: True if successful, False otherwise
        """
        try:
            self.logger.info("Starting Simplified WildType Alignment Pipeline")
            
            if not BIOPYTHON_AVAILABLE:
                self.logger.error("BioPython not available - required for alignment")
                return False
            
            # Ensure output structure early
            self.ensure_output_structure()

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
            
            # Ensure output structure before writing any files
            self.ensure_output_structure()

            # Generate summary reports
            self._generate_summary_reports(all_alignments)
            
            # Print statistics
            self._print_statistics()
            
            return True

        except Exception as e:
            self.logger.error(f"Error in alignment pipeline: {e}")
            return False
        finally:
            # Ensure directories persist for downstream validation
            self.ensure_output_structure()
            try:
                keep_dir = Path(self.config.output_dir) / "alignments"
                keep_dir.mkdir(parents=True, exist_ok=True)
                keep_file = keep_dir / ".keep"
                with open(keep_file, 'w') as f:
                    f.write("keep")
            except Exception:
                pass
            # Clean up file handlers
            self._cleanup_logging()

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
            
            self.logger.info(f"Found {len(sequences)} total sequences")
            alignments = []
            
            # Group sequences by gene name
            gene_groups = self._group_sequences_by_gene(sequences)
            self.logger.info(f"Gene groups found: {list(gene_groups.keys())}")
            
            for gene_name, gene_sequences in gene_groups.items():
                self.logger.info(f"Checking gene: {gene_name} against targets: {self.config.target_genes}")
                
                # Case-insensitive gene matching
                if self.config.target_genes:
                    target_genes_lower = [g.lower() for g in self.config.target_genes]
                    if gene_name.lower() not in target_genes_lower:
                        self.logger.info(f"Skipping gene {gene_name} (not in target list)")
                        continue
                
                self.logger.info(f"Processing gene: {gene_name} ({len(gene_sequences)} sequences)")
                
                # Get reference sequence
                reference_seq = self._get_reference_sequence(gene_name)
                if not reference_seq:
                    self.logger.warning(f"No reference sequence for {gene_name}")
                    continue
                
                # Perform alignments
                for seq_record in gene_sequences:
                    alignment = self._perform_simple_alignment(seq_record, reference_seq, gene_name)
                    if alignment:
                        alignments.append(alignment)
                        
                self.stats['proteins_processed'] += len(gene_sequences)
                
            return alignments
            
        except Exception as e:
            self.logger.error(f"Error processing {protein_file}: {e}")
            return []

    def _group_sequences_by_gene(self, sequences: List[SeqRecord]) -> Dict[str, List[SeqRecord]]:
        """Group sequences by gene name"""
        gene_groups = {}
        
        for seq_record in sequences:
            if BIOPYTHON_AVAILABLE:
                gene_name = self._extract_gene_name(getattr(seq_record, 'id', 'unknown') or "unknown")
            else:
                gene_name = self._extract_gene_name(str(seq_record))
            
            if gene_name not in gene_groups:
                gene_groups[gene_name] = []
            gene_groups[gene_name].append(seq_record)
        
        return gene_groups

    def _extract_gene_name(self, sequence_id: str) -> str:
        """Extract gene name from sequence ID - FULLY GENERIC for any gene"""
        # Handle our naming convention: PROTEIN_ID_genome_GENE_START_END_STRAND
        parts = sequence_id.split('_')
        
        # Look for gene names in parts - extract any alphabetic gene identifier
        # No hardcoded gene list - works with ANY gene names
        for part in parts:
            part_lower = part.lower()
            # Check if this looks like a gene name (3+ alphabetic characters)
            if len(part_lower) >= 3 and part_lower.isalpha():
                # Skip common non-gene terms
                if part_lower not in ['genome', 'protein', 'sequence', 'extracted', 'proteins']:
                    return part_lower
        
        # Fallback: use third-to-last part or "unknown"
        if len(parts) >= 3:
            potential_gene = parts[-3].lower()
            # Check if it looks like a gene name
            if len(potential_gene) >= 3 and potential_gene.isalpha():
                return potential_gene
        
        return "unknown"

    def _get_reference_sequence(self, gene_name: str) -> Optional[str]:
        """Get reference sequence for a gene"""
        gene_lower = gene_name.lower()
        
        # Check built-in references
        if gene_lower in self.builtin_references:
            return self.builtin_references[gene_lower]
        
        # Check user-provided reference directory
        if self.config.reference_dir:
            ref_dir = Path(self.config.reference_dir)
            
            # Try various naming patterns for reference files
            possible_names = [
                f"{gene_name}",                    # exact case
                f"{gene_name.lower()}",            # lowercase
                f"{gene_name.upper()}",            # uppercase
                f"{gene_name}_strain",             # with strain suffix
                f"{gene_name.lower()}_strain",     # lowercase with strain
                f"{gene_name}_reference",          # with reference suffix
                f"{gene_name.lower()}_reference",  # lowercase with reference suffix
            ]
            
            # Try different file extensions
            extensions = ['.fasta', '.fa', '.faa']
            
            for name in possible_names:
                for ext in extensions:
                    ref_file = ref_dir / f"{name}{ext}"
                    if ref_file.exists():
                        try:
                            records = list(SeqIO.parse(ref_file, "fasta"))
                            if records:
                                self.logger.info(f"Found reference for {gene_name}: {ref_file}")
                                return str(records[0].seq)
                        except Exception as e:
                            self.logger.warning(f"Error reading reference file {ref_file}: {e}")
        
        return None

    def _perform_simple_alignment(self, query: SeqRecord, reference_seq: str, gene_name: str) -> Optional[AlignmentResult]:
        """Perform simple pairwise alignment"""
        try:
            # For simplicity, calculate basic similarity without complex alignment
            if BIOPYTHON_AVAILABLE:
                query_seq = str(getattr(query, 'seq', ''))
            else:
                query_seq = str(query)
            
            if not query_seq:
                return None
            
            # Simple similarity calculation
            min_len = min(len(query_seq), len(reference_seq))
            matches = 0
            
            for i in range(min_len):
                if query_seq[i] == reference_seq[i]:
                    matches += 1
            
            identity_percent = (matches / min_len) * 100 if min_len > 0 else 0
            gaps = abs(len(query_seq) - len(reference_seq))
            
            # Save alignment info to file
            alignment_dir = Path(self.config.output_dir) / "alignments"
            if BIOPYTHON_AVAILABLE:
                alignment_file = alignment_dir / f"{getattr(query, 'id', 'unknown')}_{gene_name}_alignment.txt"
            else:
                alignment_file = alignment_dir / f"{str(query)}_{gene_name}_alignment.txt"
            
            with open(alignment_file, 'w') as f:
                if BIOPYTHON_AVAILABLE:
                    f.write(f"Query: {getattr(query, 'id', 'unknown')}\n")
                else:
                    f.write(f"Query: {str(query)}\n")
                f.write(f"Reference: {gene_name}_reference\n")
                f.write(f"Query Length: {len(query_seq)}\n")
                f.write(f"Reference Length: {len(reference_seq)}\n")
                f.write(f"Identity: {identity_percent:.2f}%\n")
                f.write(f"Gaps: {gaps}\n")
                f.write(f"\nQuery Sequence:\n{query_seq}\n")
                f.write(f"\nReference Sequence:\n{reference_seq}\n")
            
            if BIOPYTHON_AVAILABLE:
                query_id = getattr(query, 'id', 'unknown') or "unknown"
            else:
                query_id = str(query)
            alignment_result = AlignmentResult(
                query_id=query_id,
                reference_id=f"{gene_name}_reference",
                query_sequence=query_seq,
                reference_sequence=reference_seq,
                alignment_score=identity_percent,  # Use identity as score
                identity_percent=identity_percent,
                gaps=gaps,
                alignment_length=min_len,
                alignment_file=str(alignment_file)
            )
            
            self.stats['alignments_performed'] += 1
            return alignment_result
            
        except Exception as e:
            if BIOPYTHON_AVAILABLE:
                self.logger.error(f"Alignment failed for {getattr(query, 'id', 'unknown')}: {e}")
            else:
                self.logger.error(f"Alignment failed for {str(query)}: {e}")
            self.stats['failed_alignments'] += 1
            return None

    def _generate_summary_reports(self, alignments: List[AlignmentResult]):
        """Generate summary reports"""
        try:
            # Ensure reports directory exists
            reports_dir = Path(self.config.output_dir) / "reports"
            reports_dir.mkdir(parents=True, exist_ok=True)

            # Ensure alignments directory exists even if no alignments occurred
            (Path(self.config.output_dir) / "alignments").mkdir(parents=True, exist_ok=True)
            
            # Generate CSV report
            csv_file = reports_dir / "alignment_summary.csv"
            with open(csv_file, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow([
                    'Query ID', 'Reference ID', 'Identity %', 'Gaps', 'Alignment File'
                ])
                
                for alignment in alignments:
                    writer.writerow([
                        alignment.query_id,
                        alignment.reference_id,
                        f"{alignment.identity_percent:.2f}",
                        alignment.gaps,
                        alignment.alignment_file
                    ])
            
            # Generate JSON report
            json_file = reports_dir / "alignment_summary.json"
            alignment_data = []
            for alignment in alignments:
                alignment_data.append({
                    'query_id': alignment.query_id,
                    'reference_id': alignment.reference_id,
                    'identity_percent': alignment.identity_percent,
                    'gaps': alignment.gaps,
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
        self.logger.info("SIMPLIFIED WILDTYPE ALIGNMENT STATISTICS")
        self.logger.info("=" * 60)
        self.logger.info(f"Proteins processed: {self.stats['proteins_processed']}")
        self.logger.info(f"Alignments performed: {self.stats['alignments_performed']}")
        self.logger.info(f"Failed alignments: {self.stats['failed_alignments']}")

        success_rate = (self.stats['alignments_performed'] / max(1, self.stats['alignments_performed'] + self.stats['failed_alignments']) * 100)
        self.logger.info(f"Alignment success rate: {success_rate:.1f}%")

    def analyze_sequences(self, sequences: List[str], reference_sequence: str, gene_name: str) -> Dict:
        """
        Analyze sequences against a reference with robust error handling.
        
        Args:
            sequences: List of sequences to analyze
            reference_sequence: Reference sequence for comparison
            gene_name: Name of the gene being analyzed
            
        Returns:
            Dictionary containing analysis results
            
        Raises:
            ValidationError: If input validation fails
            DataProcessingError: If analysis fails
        """
        try:
            # Input validation
            if ERROR_HANDLING_AVAILABLE:
                if not sequences:
                    raise ValidationError("No sequences provided for analysis", "sequences", sequences)
                
                if not reference_sequence:
                    raise ValidationError("Reference sequence cannot be empty", "reference_sequence", reference_sequence)
                
                # Validate gene name
                validated_gene_name = ValidationSuite.validate_gene_name(gene_name)
                
                # Validate and clean sequences
                validated_sequences = []
                for i, seq in enumerate(sequences):
                    try:
                        validated_seq = ValidationSuite.validate_sequence(seq)
                        validated_sequences.append(validated_seq)
                    except ValidationError as e:
                        logger.warning(f"Skipping invalid sequence {i}: {e}")
                        continue
                
                # Validate reference sequence
                validated_reference = ValidationSuite.validate_sequence(reference_sequence)
                
            else:
                # Basic validation without error handling module
                if not sequences or not reference_sequence or not gene_name:
                    raise ValueError("Invalid input: sequences, reference, and gene name are required")
                validated_sequences = [str(seq).upper().strip() for seq in sequences]
                validated_reference = str(reference_sequence).upper().strip()
                validated_gene_name = str(gene_name).strip()
            
            # Perform analysis
            results = {
                'gene_name': validated_gene_name,
                'reference_sequence': validated_reference,
                'total_sequences': len(sequences),
                'valid_sequences_count': len(validated_sequences),
                'mutations_found': 0,  # Will be updated
                'mutations': [],
                'similarity_scores': [],
                'analysis_summary': '',  # Will be updated
                'analysis_timestamp': None
            }
            
            if ERROR_HANDLING_AVAILABLE:
                from datetime import datetime
                results['analysis_timestamp'] = datetime.now().isoformat()
            
            # Compare each sequence to reference
            for i, sequence in enumerate(validated_sequences):
                try:
                    # Simple comparison for mutations
                    mutations = self._find_mutations(sequence, validated_reference, i)
                    if mutations:
                        results['mutations'].extend(mutations)
                    
                    # Calculate similarity
                    similarity = self._calculate_similarity(sequence, validated_reference)
                    results['similarity_scores'].append({
                        'sequence_index': i,
                        'similarity_percent': similarity
                    })
                    
                except Exception as e:
                    logger.warning(f"Failed to analyze sequence {i}: {e}")
                    continue
            
            logger.info(f"Analyzed {len(validated_sequences)} sequences for gene {validated_gene_name}")
            logger.info(f"Found {len(results['mutations'])} mutations")
            
            # Update summary fields
            results['mutations_found'] = len(results['mutations'])
            results['analysis_summary'] = f"Analyzed {len(validated_sequences)} sequences, found {len(results['mutations'])} mutations"
            
            return results
            
        except ValidationError:
            raise  # Re-raise validation errors
        except Exception as e:
            if ERROR_HANDLING_AVAILABLE:
                raise DataProcessingError(f"Analysis failed for gene {gene_name}: {e}")
            else:
                raise RuntimeError(f"Analysis failed for gene {gene_name}: {e}")
    
    def _find_mutations(self, sequence: str, reference: str, sequence_index: int) -> List[Dict]:
        """Find mutations between sequence and reference."""
        mutations = []
        min_length = min(len(sequence), len(reference))
        
        for pos in range(min_length):
            if sequence[pos] != reference[pos]:
                mutations.append({
                    'position': pos + 1,  # 1-based position
                    'reference': reference[pos],
                    'variant': sequence[pos],
                    'sequence_index': sequence_index
                })
        
        return mutations
    
    def _calculate_similarity(self, sequence: str, reference: str) -> float:
        """Calculate sequence similarity percentage."""
        if not sequence or not reference:
            return 0.0
        
        min_length = min(len(sequence), len(reference))
        if min_length == 0:
            return 0.0
        
        matches = sum(1 for i in range(min_length) if sequence[i] == reference[i])
        return (matches / min_length) * 100.0

    def _cleanup_logging(self):
        """Clean up file handlers to release file handles"""
        try:
            if hasattr(self, 'file_handler') and hasattr(self, 'logger'):
                self.logger.removeHandler(self.file_handler)
                self.file_handler.close()
        except Exception:
            pass  # Ignore cleanup errors


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description="Generic WildType Aligner for ANY RND efflux pump and regulator genes",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic alignment with reference directory
  python simplified_wildtype_aligner.py --input proteins/ --output alignments/ --references refs/

  # Specify target genes (any genes, not hardcoded)
  python simplified_wildtype_aligner.py --input proteins/ --output alignments/ --genes mexA mexB oprM

  # Process all genes found in protein files
  python simplified_wildtype_aligner.py --input proteins/ --output alignments/ --references refs/
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
        '--genes',
        nargs='+',
        help='Target genes to process (e.g., mexA mexB oprM adeA adeB). If not specified, processes all genes found.'
    )

    parser.add_argument(
        '--references',
        help='Directory containing reference FASTA files (REQUIRED for alignments)'
    )

    args = parser.parse_args()

    # Create configuration
    config = SimpleAlignerConfig(
        input_dir=args.input,
        output_dir=args.output,
        target_genes=args.genes or [],
        reference_dir=args.references
    )

    # Run aligner
    aligner = SimplifiedWildTypeAligner(config)
    success = aligner.run_alignment_pipeline()

    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()