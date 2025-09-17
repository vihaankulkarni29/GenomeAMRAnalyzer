#!/usr/bin/env python3
"""
SubScan - Streamlined amino acid substitution analysis

This module analyzes alignment files to identify amino acid substitutions
between resistant and reference protein sequences for AMR research.

Usage:
    python sub_scan.py --alignment-dir results/alignments --metadata results/metadata.csv --output-dir results/substitutions

Author: MetaDataHarvester Pipeline
Version: 2.1 - Streamlined and Robust
"""

import os
import sys
import logging
from pathlib import Path
from typing import List, Dict, Optional
import pandas as pd
import argparse


class Config:
    """Configuration for SubScan"""
    def __init__(self, alignment_dir: str, metadata_file: str, output_dir: str):
        self.alignment_dir = Path(alignment_dir)
        self.metadata_file = Path(metadata_file)
        self.output_dir = Path(output_dir)
        self.amino_acids = set('ACDEFGHIKLMNPQRSTVWY')


class SubScan:
    """
    Streamlined substitution scanner for AMR analysis
    """

    def __init__(self, config: Config):
        self.config = config
        self.setup_logging()

    def setup_logging(self):
        """Setup logging configuration"""
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s'
        )
        self.logger = logging.getLogger('SubScan')

    def parse_needle_file(self, needle_file: str) -> tuple[str, str]:
        """Parse needle alignment file to extract sequences"""
        try:
            with open(needle_file, 'r') as f:
                content = f.read()

            # Extract sequences from alignment
            lines = content.split('\n')
            seq1_lines = []
            seq2_lines = []

            for line in lines:
                line = line.strip()
                if not line or line.startswith('#') or line.startswith('-'):
                    continue

                # Look for sequence lines (typically start with sequence name)
                if len(line.split()) >= 2:
                    parts = line.split()
                    if len(parts) >= 2 and any(char.isalpha() for char in parts[-1]):
                        # This looks like a sequence line
                        seq_name = parts[0]
                        seq_data = parts[-1]

                        if not seq1_lines:
                            seq1_lines.append(seq_data)
                        elif seq_name == seq1_lines[0] if seq1_lines else False:
                            seq1_lines.append(seq_data)
                        else:
                            seq2_lines.append(seq_data)

            resistant_seq = ''.join(seq1_lines)
            reference_seq = ''.join(seq2_lines)

            return resistant_seq, reference_seq

        except Exception as e:
            self.logger.error(f"Failed to parse {needle_file}: {e}")
            return "", ""

    def identify_substitutions(self, resistant_seq: str, reference_seq: str) -> List[Dict]:
        """Identify amino acid substitutions"""
        substitutions = []

        if len(resistant_seq) != len(reference_seq):
            self.logger.warning("Sequence lengths don't match")
            return substitutions

        # Track positions (excluding gaps)
        res_pos = 0
        ref_pos = 0

        for i in range(len(resistant_seq)):
            res_aa = resistant_seq[i]
            ref_aa = reference_seq[i]

            # Update position counters
            if res_aa != '-':
                res_pos += 1
            if ref_aa != '-':
                ref_pos += 1

            # Check for substitutions
            if (res_aa != ref_aa and
                res_aa != '-' and ref_aa != '-' and
                res_aa in self.amino_acids and ref_aa in self.amino_acids):

                substitution = {
                    'position': ref_pos,
                    'reference_aa': ref_aa,
                    'resistant_aa': res_aa,
                    'substitution': f"{ref_aa}{ref_pos}{res_aa}"
                }
                substitutions.append(substitution)

        return substitutions

    def process_alignment_file(self, alignment_file: str, metadata_df: pd.DataFrame) -> List[Dict]:
        """Process a single alignment file"""
        try:
            # Parse alignment
            resistant_seq, reference_seq = self.parse_needle_file(alignment_file)

            if not resistant_seq or not reference_seq:
                return []

            # Extract accession from filename
            filename = Path(alignment_file).name
            accession = filename.split('_')[0] if '_' in filename else filename.split('.')[0]

            # Find metadata
            metadata_row = metadata_df[
                metadata_df['Accession_Number'].str.contains(accession, na=False)
            ]

            if metadata_row.empty:
                self.logger.warning(f"No metadata for {accession}")
                return []

            protein_info = metadata_row.iloc[0]

            # Identify substitutions
            substitutions = self.identify_substitutions(resistant_seq, reference_seq)

            # Create records
            records = []
            for sub in substitutions:
                record = {
                    'Accession_Number': protein_info['Accession_Number'],
                    'Organism': protein_info['Organism'],
                    'Protein_Name': protein_info['Protein_Name'],
                    'Protein_Family': protein_info['Protein_Family'],
                    'Substitution': sub['substitution'],
                    'Resistant_Amino_Acid': sub['resistant_aa'],
                    'Sensitive_Amino_Acid': sub['reference_aa'],
                    'Residue_Position': sub['position'],
                    'Alignment_File': filename
                }
                records.append(record)

            self.logger.info(f"Found {len(substitutions)} substitutions in {filename}")
            return records

        except Exception as e:
            self.logger.error(f"Error processing {alignment_file}: {e}")
            return []

    def process_alignment_directory(self, alignment_dir: str, metadata_file: str, output_dir: str) -> str:
        """Process all alignment files in directory"""
        try:
            # Load metadata
            metadata_df = pd.read_csv(metadata_file)
            self.logger.info(f"Loaded metadata for {len(metadata_df)} proteins")

            # Find alignment files
            alignment_path = Path(alignment_dir)
            alignment_files = list(alignment_path.glob("*.needle")) + list(alignment_path.glob("*.aln"))

            self.logger.info(f"Found {len(alignment_files)} alignment files")

            # Process each file
            all_substitutions = []
            for alignment_file in alignment_files:
                records = self.process_alignment_file(str(alignment_file), metadata_df)
                all_substitutions.extend(records)

            # Save results
            output_path = Path(output_dir)
            output_path.mkdir(parents=True, exist_ok=True)

            if all_substitutions:
                df = pd.DataFrame(all_substitutions)
                output_file = output_path / "substitutions.csv"
                df.to_csv(output_file, index=False)
                self.logger.info(f"Saved {len(all_substitutions)} substitutions to {output_file}")
            else:
                self.logger.warning("No substitutions found")

            return str(output_path)

        except Exception as e:
            self.logger.error(f"Failed to process alignment directory: {e}")
            raise


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description="Streamlined amino acid substitution analysis"
    )

    parser.add_argument(
        "--alignment-dir",
        required=True,
        help="Directory containing alignment files"
    )

    parser.add_argument(
        "--metadata",
        required=True,
        help="Metadata CSV file"
    )

    parser.add_argument(
        "--output-dir",
        required=True,
        help="Output directory for results"
    )

    args = parser.parse_args()

    # Create configuration
    config = Config(
        alignment_dir=args.alignment_dir,
        metadata_file=args.metadata,
        output_dir=args.output_dir
    )

    # Run analysis
    print("Initializing SubScan...")
    scanner = SubScan(config)

    print("Starting substitution analysis...")
    output_dir = scanner.process_alignment_directory(
        args.alignment_dir,
        args.metadata,
        args.output_dir
    )

    print("\n============================================================")
    print("SUBSTITUTION ANALYSIS COMPLETE!")
    print("============================================================")
    print(f"Output directory: {output_dir}")

    # Count results
    results_file = Path(output_dir) / "substitutions.csv"
    if results_file.exists():
        df = pd.read_csv(results_file)
        print(f"Substitutions found: {len(df)}")

        if len(df) > 0:
            print("\nTop substitutions:")
            top_subs = df['Substitution'].value_counts().head(5)
            for sub, count in top_subs.items():
                print(f"   {sub}: {count} times")
    else:
        print("No substitutions found")

    print("\nNext steps:")
    print("   1. Review substitution results in CSV file")
    print("   2. Analyze patterns for AMR research")
    print("   3. Use data for ML model training")


if __name__ == "__main__":
    main()