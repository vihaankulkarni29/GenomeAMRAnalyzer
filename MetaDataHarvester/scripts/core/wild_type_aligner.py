import sys; import os; sys.path.append(os.path.join(os.path.dirname(__file__), '\\\\\\\\\src')); from configuration_manager import config_manager
#!/usr/bin/env python3
"""
WildTypeAligner - Streamlined sequence alignment for AMR analysis

This module performs pairwise sequence alignment between harvested proteins
and their corresponding wild-type reference sequences for mutation detection.

Usage:
    python wild_type_aligner.py --fasta-dir results/fasta --metadata results/metadata.csv --output-dir results/alignments

Author: MetaDataHarvester Pipeline
Version: 2.1 - Streamlined and Robust
"""

import os
import sys
import logging
from pathlib import Path
from typing import List, Dict, Optional
import pandas as pd
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio.SeqRecord import SeqRecord
import argparse


class Config:
    """Configuration for WildTypeAligner"""
    def __init__(self, reference_dir: str = "references", output_dir: str = "results/alignments",
                 gap_open: float = -10.0, gap_extend: float = -0.5):
        self.reference_dir = Path(reference_dir)
        self.output_dir = Path(output_dir)
        self.gap_open = gap_open
        self.gap_extend = gap_extend


class WildTypeAligner:
    """
    Streamlined sequence aligner for AMR analysis
    """

    def __init__(self, config: Config):
        self.config = config
        self.reference_sequences = {}
        self.setup_logging()
        self.load_reference_sequences()

    def setup_logging(self):
        """Setup logging configuration"""
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        self.logger = logging.getLogger('WildTypeAligner')

    def load_reference_sequences(self):
        """Load all reference sequences from the reference directory"""
        if not self.config.reference_dir.exists():
            self.logger.error(f"Reference directory not found: {self.config.reference_dir}")
            return

        ref_files = list(self.config.reference_dir.glob("*.fasta")) + list(self.config.reference_dir.glob("*.fa"))

        for ref_file in ref_files:
            try:
                for record in SeqIO.parse(ref_file, "fasta"):
                    protein_family = self.extract_protein_family(ref_file.name, record.description)
                    if protein_family:
                        self.reference_sequences[protein_family.lower()] = record
                        self.logger.info(f"Loaded reference for {protein_family}: {record.id}")
            except Exception as e:
                self.logger.error(f"Error loading {ref_file}: {e}")
                continue

        self.logger.info(f"Loaded {len(self.reference_sequences)} reference sequences")

    def extract_protein_family(self, filename: str, description: str) -> Optional[str]:
        """Extract protein family from filename or description"""
        # Protein families to look for
        families = [
            'acra', 'acrb', 'tolc', 'acrz', 'acrr', 'acrf', 'acre', 'acrd',
            'mdtb', 'mdtc', 'oqxa', 'oqxb', 'mara', 'marr', 'rama', 'ramr',
            'soxs', 'rob', 'envr', 'eefa', 'eefb', 'eefc'
        ]

        # Check filename first
        filename_lower = filename.lower()
        for family in families:
            if family in filename_lower:
                return family

        # Check description
        description_lower = description.lower()
        for family in families:
            if family in description_lower:
                return family

        return None

    def align_sequences(self, seq1: SeqRecord, seq2: SeqRecord) -> tuple[str, str, float]:
        """Perform pairwise global alignment"""
        try:
            aligner = PairwiseAligner()
            aligner.mode = 'global'
            aligner.gap_score = self.config.gap_open
            aligner.extend_gap_score = self.config.gap_extend
            aligner.match_score = 2
            aligner.mismatch_score = -1

            alignments = aligner.align(str(seq1.seq), str(seq2.seq))
            if not alignments:
                raise ValueError("No alignment generated")

            best = alignments[0]
            return str(best[0]), str(best[1]), best.score

        except Exception as e:
            self.logger.error(f"Alignment failed: {e}")
            raise

    def save_alignment(self, seq1: SeqRecord, seq2: SeqRecord,
                      aligned_seq1: str, aligned_seq2: str, score: float, output_file: str):
        """Save alignment in needle-like format"""
        try:
            with open(output_file, 'w') as f:
                f.write("# EMBOSS needle-like format\n")
                f.write(f"# Sequence 1: {seq1.id}\n")
                f.write(f"# Sequence 2: {seq2.id}\n")
                f.write(f"# Score: {score:.1f}\n")
                f.write("#\n")

                # Write alignment in blocks
                block_size = 60
                for i in range(0, len(aligned_seq1), block_size):
                    seq1_block = aligned_seq1[i:i + block_size]
                    seq2_block = aligned_seq2[i:i + block_size]

                    f.write(f"{seq1.id}\t{seq1_block}\n")

                    # Match line
                    match_line = ""
                    for a, b in zip(seq1_block, seq2_block):
                        if a == b and a != '-':
                            match_line += "|"
                        elif a == '-' or b == '-':
                            match_line += " "
                        else:
                            match_line += "."
                    f.write(f"\t\t{match_line}\n")
                    f.write(f"{seq2.id}\t{seq2_block}\n\n")

        except Exception as e:
            self.logger.error(f"Failed to save alignment: {e}")
            raise

    def process_fasta_directory(self, fasta_dir: str, metadata_df: pd.DataFrame, output_dir: str) -> List[str]:
        """Process all FASTA files and generate alignments"""
        fasta_path = Path(fasta_dir)
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        alignment_files = []

        if not self.reference_sequences:
            self.logger.error("No reference sequences loaded")
            return alignment_files

        fasta_files = list(fasta_path.glob("*.fasta")) + list(fasta_path.glob("*.fa"))
        self.logger.info(f"Found {len(fasta_files)} FASTA files to process")

        for fasta_file in fasta_files:
            try:
                # Parse all sequences in the FASTA file
                sequences = list(SeqIO.parse(fasta_file, "fasta"))
                self.logger.info(f"Processing {len(sequences)} sequences from {fasta_file.name}")

                for seq_record in sequences:
                    # Find protein family from metadata
                    accession = seq_record.id.split('.')[0]
                    metadata_row = metadata_df[
                        metadata_df['Accession_Number'].str.contains(accession, na=False)
                    ]

                    if metadata_row.empty:
                        self.logger.warning(f"No metadata for {seq_record.id}")
                        continue

                    protein_family = metadata_row.iloc[0]['Protein_Family']
                    if not protein_family:
                        self.logger.warning(f"No protein family for {seq_record.id}")
                        continue

                    # Handle special cases for protein family mapping
                    protein_family_lower = protein_family.lower()
                    if protein_family_lower == 'rnd_efflux':
                        # Map RND_Efflux to acra reference
                        protein_family_lower = 'acra'
                        self.logger.info(f"Mapped RND_Efflux to acra reference for {seq_record.id}")
                        protein_family = 'acra'  # Also update the original for logging

                    if protein_family_lower not in self.reference_sequences:
                        self.logger.warning(f"No reference for {protein_family} (available: {list(self.reference_sequences.keys())})")
                        continue

                    # Get reference and align
                    reference_seq = self.reference_sequences[protein_family.lower()]
                    self.logger.info(f"🔬 Processing {protein_family.upper()} sequence: {seq_record.id}")

                    aligned_seq1, aligned_seq2, score = self.align_sequences(seq_record, reference_seq)

                    # Save alignment with sanitized filename
                    safe_id = seq_record.id.replace('|', '_').replace('/', '_').replace('\\', '_')
                    safe_ref_id = reference_seq.id.replace('|', '_').replace('/', '_').replace('\\', '_')
                    output_file = output_path / f"{safe_id}_vs_{safe_ref_id}.needle"

                    self.save_alignment(seq_record, reference_seq, aligned_seq1, aligned_seq2, score, str(output_file))
                    alignment_files.append(str(output_file))

                    self.logger.info(f"Processed alignment: {output_file.name}")

            except Exception as e:
                self.logger.error(f"Error processing {fasta_file}: {e}")
                continue

        self.logger.info(f"Generated {len(alignment_files)} alignment files")
        return alignment_files

    def run_alignment(self, fasta_dir: str, metadata_file: str, output_dir: str) -> str:
        """Run the complete alignment workflow"""
        try:
            # Load metadata
            metadata_df = pd.read_csv(metadata_file)
            self.logger.info(f"Loaded metadata for {len(metadata_df)} proteins")

            # Process alignments
            alignment_files = self.process_fasta_directory(fasta_dir, metadata_df, output_dir)

            self.logger.info("Alignment workflow complete!")
            return output_dir

        except Exception as e:
            self.logger.error(f"Alignment workflow failed: {e}")
            raise


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description="Streamlined sequence alignment for AMR analysis"
    )

    parser.add_argument(
        "--fasta-dir",
        required=True,
        help="Directory containing FASTA files"
    )

    parser.add_argument(
        "--metadata",
        required=True,
        help="Metadata CSV file"
    )

    parser.add_argument(
        "--output-dir",
        required=True,
        help="Output directory for alignments"
    )

    parser.add_argument(
        "--reference-dir",
        default="references",
        help="Reference sequences directory"
    )

    parser.add_argument(
        "--gap-open",
        type=float,
        default=-10.0,
        help="Gap opening penalty"
    )

    parser.add_argument(
        "--gap-extend",
        type=float,
        default=-0.5,
        help="Gap extension penalty"
    )

    args = parser.parse_args()

    # Create configuration
    config = Config(
        reference_dir=args.reference_dir,
        output_dir=args.output_dir,
        gap_open=args.gap_open,
        gap_extend=args.gap_extend
    )

    # Run alignment
    print("Initializing WildTypeAligner...")
    aligner = WildTypeAligner(config)

    print("Starting sequence alignment...")
    output_dir = aligner.run_alignment(args.fasta_dir, args.metadata, args.output_dir)

    print("\n============================================================")
    print("ALIGNMENT COMPLETE!")
    print("============================================================")
    print(f"Output directory: {output_dir}")
    print(f"Reference sequences used: {len(aligner.reference_sequences)}")

    # Count generated files
    alignment_files = list(Path(output_dir).glob("*.needle"))
    print(f"Alignment files generated: {len(alignment_files)}")

    if alignment_files:
        print("\nSample alignment files:")
        for i, file in enumerate(alignment_files[:5]):
            print(f"   {i+1}. {file.name}")

    print("\nNext steps:")
    print("   1. Review alignment quality in .needle files")
    print("   2. Run SubScan for substitution analysis")
    print("   3. Analyze results for AMR patterns")


if __name__ == "__main__":
    main()