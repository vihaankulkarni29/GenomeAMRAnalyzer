import sys; import os; sys.path.append(os.path.join(os.path.dirname(__file__), '\\\\\\src')); from configuration_manager import config_manager
#!/usr/bin/env python3
"""
SetupSpeciesReferences - Build species-specific reference database
Creates genus-specific reference sequences to eliminate false positives

Author: MetaDataHarvester Pipeline
Version: 2.0 - Species-Specific References
"""

import os
import sys
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
from Bio import SeqIO, Entrez
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIWWW, NCBIXML
import pandas as pd
import json
import time
from collections import defaultdict


class SpeciesReferenceBuilder:
    """
    Builds species-specific reference database for accurate mutation detection
    """

    def __init__(self, email: str, output_dir: str = "references"):
        self.email = email
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)

        # Setup logging
        self.logger = logging.getLogger(__name__)
        self._setup_logging()

        # Initialize NCBI
        Entrez.email = email

        # Define target genera and proteins
        self.target_genera = [
            'Escherichia', 'Pseudomonas', 'Klebsiella', 'Acinetobacter',
            'Salmonella', 'Shigella', 'Enterobacter', 'Burkholderia'
        ]

        self.target_proteins = {
            'AcrA': 'multidrug efflux RND transporter periplasmic adaptor subunit AcrA',
            'AcrB': 'multidrug efflux RND transporter permease subunit AcrB',
            'TolC': 'multidrug efflux RND transporter outer membrane lipoprotein TolC',
            'MexA': 'multidrug efflux RND transporter periplasmic adaptor subunit MexA',
            'MexB': 'multidrug efflux RND transporter permease subunit MexB',
            'OprM': 'multidrug efflux RND transporter outer membrane lipoprotein OprM'
        }

    def _setup_logging(self):
        """Setup comprehensive logging"""
        log_file = self.output_dir / "reference_builder.log"
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )

    def build_species_reference_database(self) -> None:
        """
        Build complete species-specific reference database
        """
        self.logger.info("Building species-specific reference database")
        self.logger.info(f"Target genera: {self.target_genera}")
        self.logger.info(f"Target proteins: {list(self.target_proteins.keys())}")

        total_references = 0

        for genus in self.target_genera:
            self.logger.info(f"Processing genus: {genus}")
            genus_dir = self.output_dir / genus
            genus_dir.mkdir(exist_ok=True)

            genus_references = 0

            for protein_name, protein_description in self.target_proteins.items():
                try:
                    reference_seq = self._find_genus_specific_reference(genus, protein_name, protein_description)

                    if reference_seq:
                        # Save reference sequence
                        fasta_file = genus_dir / f"{protein_name}_reference.fasta"
                        SeqIO.write([reference_seq], fasta_file, "fasta")

                        # Save metadata
                        metadata_file = genus_dir / f"{protein_name}_reference.json"
                        metadata = {
                            'genus': genus,
                            'protein_name': protein_name,
                            'accession': reference_seq.id,
                            'description': reference_seq.description,
                            'sequence_length': len(reference_seq.seq),
                            'source': 'NCBI_refseq',
                            'selection_criteria': 'representative_wild_type'
                        }

                        with open(metadata_file, 'w') as f:
                            json.dump(metadata, f, indent=2)

                        genus_references += 1
                        self.logger.info(f"  ✓ {genus} {protein_name}: {reference_seq.id}")

                    else:
                        self.logger.warning(f"  ✗ No reference found for {genus} {protein_name}")

                except Exception as e:
                    self.logger.error(f"  ✗ Error processing {genus} {protein_name}: {e}")

            self.logger.info(f"Completed {genus}: {genus_references} references created")
            total_references += genus_references

        self.logger.info(f"Reference database build complete: {total_references} total references")

    def _find_genus_specific_reference(self, genus: str, protein_name: str,
                                     protein_description: str) -> Optional[SeqRecord]:
        """
        Find representative reference sequence for specific genus and protein
        """
        # Search for reference sequences in RefSeq
        query = f'"{genus}"[Organism] AND "{protein_description}"[Title] AND refseq[Filter]'

        try:
            # Search NCBI
            handle = Entrez.esearch(db="protein", term=query, retmax=50)
            search_results = Entrez.read(handle)
            handle.close()

            accession_ids = search_results.get('IdList', [])

            if not accession_ids:
                self.logger.debug(f"No RefSeq entries found for {genus} {protein_name}")
                return None

            # Fetch sequences and select best representative
            handle = Entrez.efetch(db="protein", id=accession_ids, rettype="fasta", retmode="text")
            sequences = list(SeqIO.parse(handle, "fasta"))
            handle.close()

            if not sequences:
                return None

            # Select best representative sequence
            best_reference = self._select_representative_sequence(sequences, genus, protein_name)

            return best_reference

        except Exception as e:
            self.logger.warning(f"Error searching for {genus} {protein_name}: {e}")
            return None

    def _select_representative_sequence(self, sequences: List[SeqRecord],
                                      genus: str, protein_name: str) -> SeqRecord:
        """
        Select most representative sequence from candidates
        """
        if not sequences:
            return None

        # Score sequences based on multiple criteria
        scored_sequences = []

        for seq in sequences:
            score = 0
            reasons = []

            # Prefer complete sequences
            if len(seq.seq) > 300:  # Reasonable minimum length
                score += 2
                reasons.append("complete_sequence")

            # Prefer well-annotated sequences
            if "multidrug" in seq.description.lower():
                score += 2
                reasons.append("resistance_annotated")

            # Prefer reference strains
            if any(strain in seq.description.lower() for strain in
                   ['k-12', 'mg1655', 'pao1', 'atcc', 'nctc']):
                score += 1
                reasons.append("reference_strain")

            # Prefer sequences with proper genus in description
            if genus.lower() in seq.description.lower():
                score += 1
                reasons.append("genus_match")

            scored_sequences.append((seq, score, reasons))

        # Sort by score (highest first)
        scored_sequences.sort(key=lambda x: x[1], reverse=True)

        best_seq, best_score, reasons = scored_sequences[0]

        self.logger.debug(
            f"Selected {best_seq.id} for {genus} {protein_name} "
            f"(score: {best_score}, reasons: {reasons})"
        )

        return best_seq

    def validate_reference_database(self) -> Dict[str, Dict]:
        """
        Validate completeness and quality of reference database
        """
        validation_results = {}

        for genus_dir in self.output_dir.iterdir():
            if genus_dir.is_dir() and genus_dir.name in self.target_genera:
                genus = genus_dir.name
                validation_results[genus] = {}

                for protein_name in self.target_proteins.keys():
                    fasta_file = genus_dir / f"{protein_name}_reference.fasta"
                    metadata_file = genus_dir / f"{protein_name}_reference.json"

                    status = {
                        'sequence_exists': fasta_file.exists(),
                        'metadata_exists': metadata_file.exists(),
                        'sequence_valid': False,
                        'metadata_valid': False,
                        'warnings': []
                    }

                    # Validate sequence
                    if status['sequence_exists']:
                        try:
                            records = list(SeqIO.parse(fasta_file, 'fasta'))
                            if len(records) == 1 and len(records[0].seq) > 0:
                                status['sequence_valid'] = True
                            else:
                                status['warnings'].append("Invalid sequence format")
                        except Exception as e:
                            status['warnings'].append(f"Sequence parsing error: {e}")

                    # Validate metadata
                    if status['metadata_exists']:
                        try:
                            with open(metadata_file, 'r') as f:
                                metadata = json.load(f)
                                required_fields = ['genus', 'protein_name', 'accession']
                                if all(field in metadata for field in required_fields):
                                    status['metadata_valid'] = True
                                else:
                                    status['warnings'].append("Missing required metadata fields")
                        except Exception as e:
                            status['warnings'].append(f"Metadata parsing error: {e}")

                    validation_results[genus][protein_name] = status

        return validation_results

    def generate_reference_report(self) -> None:
        """Generate comprehensive reference database report"""
        report_file = self.output_dir / "reference_database_report.txt"

        # Validate database
        validation_results = self.validate_reference_database()

        with open(report_file, 'w') as f:
            f.write("SPECIES-SPECIFIC REFERENCE DATABASE REPORT\n")
            f.write("=" * 50 + "\n\n")

            total_expected = len(self.target_genera) * len(self.target_proteins)
            total_present = 0

            f.write("DATABASE COMPLETENESS:\n")
            f.write("-" * 25 + "\n")

            for genus in self.target_genera:
                genus_results = validation_results.get(genus, {})
                genus_present = sum(1 for status in genus_results.values() if status['sequence_valid'])

                f.write(f"{genus}: {genus_present}/{len(self.target_proteins)} proteins\n")
                total_present += genus_present

            f.write(f"Overall completeness: {total_present/total_expected:.1f}\n")
            f.write("\n")

            f.write("DETAILED VALIDATION:\n")
            f.write("-" * 25 + "\n")

            for genus in self.target_genera:
                f.write(f"\n{genus}:\n")

                for protein_name in self.target_proteins.keys():
                    status = validation_results.get(genus, {}).get(protein_name, {})

                    if status.get('sequence_valid'):
                        f.write(f"  ✓ {protein_name}: Valid\n")
                    else:
                        f.write(f"  ✗ {protein_name}: Missing or invalid\n")
                        if status.get('warnings'):
                            for warning in status['warnings']:
                                f.write(f"    - {warning}\n")

            f.write("\nREFERENCE SEQUENCES BY GENUS:\n")
            f.write("-" * 35 + "\n")

            for genus_dir in self.output_dir.iterdir():
                if genus_dir.is_dir() and genus_dir.name in self.target_genera:
                    genus = genus_dir.name
                    fasta_files = list(genus_dir.glob("*_reference.fasta"))

                    f.write(f"{genus}: {len(fasta_files)} reference sequences\n")

                    for fasta_file in fasta_files:
                        try:
                            records = list(SeqIO.parse(fasta_file, 'fasta'))
                            if records:
                                record = records[0]
                                f.write(f"  - {fasta_file.stem}: {record.id} ({len(record.seq)} aa)\n")
                        except Exception as e:
                            f.write(f"  - {fasta_file.stem}: Error reading sequence\n")

        self.logger.info(f"Generated reference database report: {report_file}")

    def create_phylogenetic_mapping(self) -> None:
        """
        Create phylogenetic relationship mapping for reference selection fallback
        """
        phylogenetic_file = self.output_dir / "phylogenetic_mapping.json"

        phylogenetic_map = {
            'Escherichia': {
                'close_relatives': ['Shigella', 'Salmonella'],
                'distant_relatives': ['Klebsiella', 'Enterobacter'],
                'phylum': 'Proteobacteria'
            },
            'Pseudomonas': {
                'close_relatives': ['Acinetobacter', 'Burkholderia'],
                'distant_relatives': ['Stenotrophomonas'],
                'phylum': 'Proteobacteria'
            },
            'Klebsiella': {
                'close_relatives': ['Escherichia', 'Salmonella', 'Enterobacter'],
                'distant_relatives': ['Shigella'],
                'phylum': 'Proteobacteria'
            },
            'Acinetobacter': {
                'close_relatives': ['Pseudomonas', 'Burkholderia'],
                'distant_relatives': ['Stenotrophomonas'],
                'phylum': 'Proteobacteria'
            }
        }

        with open(phylogenetic_file, 'w') as f:
            json.dump(phylogenetic_map, f, indent=2)

        self.logger.info(f"Created phylogenetic mapping: {phylogenetic_file}")

    def backup_existing_references(self) -> None:
        """Backup any existing reference files before rebuilding"""
        backup_dir = self.output_dir / "backup"
        backup_dir.mkdir(exist_ok=True)

        # Find all existing reference files
        existing_files = list(self.output_dir.glob("**/*_reference.fasta"))
        existing_files.extend(self.output_dir.glob("**/*_reference.json"))

        if existing_files:
            self.logger.info(f"Backing up {len(existing_files)} existing reference files")

            for file_path in existing_files:
                # Create relative path for backup
                relative_path = file_path.relative_to(self.output_dir)
                backup_path = backup_dir / relative_path
                backup_path.parent.mkdir(parents=True, exist_ok=True)

                # Copy file
                import shutil
                shutil.copy2(file_path, backup_path)

            self.logger.info(f"Backup completed in: {backup_dir}")


def main():
    """Main entry point"""
    import argparse

    parser = argparse.ArgumentParser(
        description="Build species-specific reference database for AMR analysis"
    )

    parser.add_argument(
        "--email",
        required=True,
        help="NCBI email for API access"
    )

    parser.add_argument(
        "--output-dir",
        default="references",
        help="Output directory for reference database"
    )

    parser.add_argument(
        "--backup",
        action="store_true",
        help="Backup existing references before rebuilding"
    )

    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="Only validate existing database, don't rebuild"
    )

    args = parser.parse_args()

    # Setup logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    # Initialize builder
    builder = SpeciesReferenceBuilder(args.email, args.output_dir)

    if args.backup:
        builder.backup_existing_references()

    if not args.validate_only:
        # Build reference database
        builder.build_species_reference_database()
        builder.create_phylogenetic_mapping()

    # Validate and report
    builder.generate_reference_report()

    print(f"\nSpecies-specific reference database setup complete!")
    print(f"Database location: {args.output_dir}")
    print(f"See reference_database_report.txt for details")


if __name__ == "__main__":
    main()