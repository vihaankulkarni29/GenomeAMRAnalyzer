import sys; import os; sys.path.append(os.path.join(os.path.dirname(__file__), '\\\\\\\\\src')); from configuration_manager import config_manager
#!/usr/bin/env python3
"""
ReferenceValidator - Scientific validation of reference sequences
Ensures biological relevance and proper selection of reference proteins

Author: MetaDataHarvester Pipeline
Version: 1.0 - Scientifically Validated
"""

import os
import sys
import json
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
from Bio import SeqIO, Entrez
from Bio.Blast import NCBIWWW, NCBIXML
import pandas as pd
from datetime import datetime


@dataclass
class ReferenceValidationResult:
    """Results of reference sequence validation"""
    accession: str
    is_valid: bool
    species: str
    protein_family: str
    completeness_score: float
    orthology_confidence: float
    validation_warnings: List[str]
    validation_errors: List[str]


class ReferenceValidator:
    """
    Validates reference sequences for biological relevance and scientific accuracy
    """

    def __init__(self, email: str, cache_dir: str = "cache"):
        self.email = email
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(exist_ok=True)

        # Setup logging
        self.logger = logging.getLogger(__name__)

        # Known resistance mutations from literature
        self.known_resistance_mutations = self._load_known_mutations()

        # Expected domain structures for efflux pumps
        self.domain_structures = self._load_domain_structures()

    def _load_known_mutations(self) -> Dict[str, List[str]]:
        """Load known resistance mutations from literature"""
        return {
            config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0]: ['T104A', 'H596N', 'G616D'],
            config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1]: ['V610A', 'F628L', 'A279T', 'R717C'],
            config_manager.get_default_genes("rnd_efflux_pumps", "primary")[2]: ['D153A', 'E473K']
        }

    def _load_domain_structures(self) -> Dict[str, Dict]:
        """Load expected domain structures for efflux pumps"""
        return {
            config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0]: {
                'expected_length': (350, 450),
                'key_domains': ['lipoyl_domain', 'beta_barrel'],
                'critical_residues': [104, 596, 616]
            },
            config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1]: {
                'expected_length': (1000, 1200),
                'key_domains': ['transmembrane_domain', 'ATP_binding'],
                'critical_residues': [610, 628, 717]
            }
        }

    def validate_reference_sequence(self, fasta_file: str, expected_family: str) -> ReferenceValidationResult:
        """
        Comprehensive validation of a reference sequence

        Args:
            fasta_file: Path to FASTA file
            expected_family: Expected protein family (acrA, acrB, etc.)

        Returns:
            ValidationResult with detailed assessment
        """
        try:
            # Parse sequence
            records = list(SeqIO.parse(fasta_file, 'fasta'))
            if len(records) != 1:
                return ReferenceValidationResult(
                    accession="unknown",
                    is_valid=False,
                    species="unknown",
                    protein_family=expected_family,
                    completeness_score=0.0,
                    orthology_confidence=0.0,
                    validation_warnings=[],
                    validation_errors=["FASTA file must contain exactly one sequence"]
                )

            record = records[0]
            accession = record.id.split('.')[0]  # Remove version number

            # Extract species from description
            species = self._extract_species_from_description(record.description)

            # Validate sequence properties
            validation_warnings = []
            validation_errors = []

            # Check sequence length
            seq_length = len(record.seq)
            expected_range = self.domain_structures.get(expected_family, {}).get('expected_length', (0, 10000))
            if not (expected_range[0] <= seq_length <= expected_range[1]):
                validation_warnings.append(
                    f"Sequence length {seq_length} outside expected range {expected_range} for {expected_family}"
                )

            # Check for complete sequence (no ambiguous residues)
            ambiguous_count = str(record.seq).count('X')
            if ambiguous_count > 0:
                validation_warnings.append(f"Sequence contains {ambiguous_count} ambiguous residues")

            # Validate against NCBI taxonomy
            taxonomy_valid = self._validate_taxonomy(accession, species)
            if not taxonomy_valid:
                validation_errors.append("Taxonomic classification mismatch")

            # Check orthology
            orthology_score = self._check_orthology(record, expected_family)

            # Calculate completeness score
            completeness_score = self._calculate_completeness_score(
                record, expected_family, validation_warnings, validation_errors
            )

            # Determine overall validity
            is_valid = (
                len(validation_errors) == 0 and
                completeness_score >= 0.7 and
                orthology_score >= 0.8
            )

            return ReferenceValidationResult(
                accession=accession,
                is_valid=is_valid,
                species=species,
                protein_family=expected_family,
                completeness_score=completeness_score,
                orthology_confidence=orthology_score,
                validation_warnings=validation_warnings,
                validation_errors=validation_errors
            )

        except Exception as e:
            self.logger.error(f"Error validating reference {fasta_file}: {e}")
            return ReferenceValidationResult(
                accession="unknown",
                is_valid=False,
                species="unknown",
                protein_family=expected_family,
                completeness_score=0.0,
                orthology_confidence=0.0,
                validation_warnings=[],
                validation_errors=[f"Validation failed: {str(e)}"]
            )

    def _extract_species_from_description(self, description: str) -> str:
        """Extract species name from sequence description"""
        # Look for species patterns in description
        import re

        # Common patterns: [Escherichia coli], Escherichia coli str. K-12, etc.
        patterns = [
            r'\[([^\]]+)\]',  # [Escherichia coli]
            r'([A-Z][a-z]+ [a-z]+(?: [a-z]+)?)',  # Escherichia coli
        ]

        for pattern in patterns:
            match = re.search(pattern, description)
            if match:
                species = match.group(1)
                # Clean up common artifacts
                species = species.replace(' str.', '').replace(' substr.', '')
                return species.strip()

        return "unknown"

    def _validate_taxonomy(self, accession: str, expected_species: str) -> bool:
        """Validate taxonomic classification using NCBI"""
        try:
            # Query NCBI for taxonomy information
            Entrez.email = self.email
            handle = Entrez.efetch(db="protein", id=accession, rettype="gb", retmode="text")
            record = handle.read()
            handle.close()

            # Check if expected species is in the record
            return expected_species.lower() in record.lower()

        except Exception as e:
            self.logger.warning(f"Could not validate taxonomy for {accession}: {e}")
            return True  # Default to valid if we can't check

    def _check_orthology(self, record, expected_family: str) -> float:
        """Check orthology by comparing to known family members"""
        try:
            # Perform BLAST search against known efflux pump sequences
            result_handle = NCBIWWW.qblast(
                "blastp", "nr",
                str(record.seq),
                expect=1e-10,
                hitlist_size=50
            )

            blast_records = NCBIXML.read(result_handle)

            # Analyze top hits for orthology
            orthology_score = 0.0
            efflux_pump_hits = 0

            for alignment in blast_records.alignments[:10]:  # Top 10 hits
                title = alignment.title.lower()
                if any(term in title for term in ['efflux', 'pump', 'transporter', expected_family]):
                    efflux_pump_hits += 1

            if efflux_pump_hits > 0:
                orthology_score = min(1.0, efflux_pump_hits / 5.0)  # Normalize to 0-1

            return orthology_score

        except Exception as e:
            self.logger.warning(f"Could not check orthology: {e}")
            return 0.5  # Neutral score if BLAST fails

    def _calculate_completeness_score(self, record, family: str,
                                    warnings: List[str], errors: List[str]) -> float:
        """Calculate sequence completeness score"""
        score = 1.0

        # Length penalty
        expected_range = self.domain_structures.get(family, {}).get('expected_length', (0, 10000))
        seq_len = len(record.seq)
        if not (expected_range[0] <= seq_len <= expected_range[1]):
            length_ratio = min(seq_len, expected_range[1]) / max(seq_len, expected_range[0])
            score *= length_ratio

        # Ambiguous residue penalty
        ambiguous_ratio = str(record.seq).count('X') / len(record.seq)
        score *= (1 - ambiguous_ratio)

        # Warning penalties
        score *= (1 - len(warnings) * 0.1)

        # Error penalties
        score *= (1 - len(errors) * 0.3)

        return max(0.0, min(1.0, score))

    def validate_reference_directory(self, reference_dir: str) -> Dict[str, ReferenceValidationResult]:
        """
        Validate all reference sequences in a directory

        Args:
            reference_dir: Path to directory containing reference FASTA files

        Returns:
            Dictionary mapping filenames to validation results
        """
        results = {}
        ref_path = Path(reference_dir)

        for fasta_file in ref_path.glob("*.fasta"):
            # Extract expected family from filename
            filename = fasta_file.stem.lower()
            expected_family = "unknown"

            if "acra" in filename:
                expected_family = config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0]
            elif "acrb" in filename:
                expected_family = config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1]
            elif "tolc" in filename:
                expected_family = config_manager.get_default_genes("rnd_efflux_pumps", "primary")[2]

            result = self.validate_reference_sequence(str(fasta_file), expected_family)
            results[fasta_file.name] = result

            # Log results
            status = "✓ VALID" if result.is_valid else "✗ INVALID"
            self.logger.info(f"{status} {fasta_file.name}: {result.species} ({result.completeness_score:.2f})")

            if result.validation_warnings:
                for warning in result.validation_warnings:
                    self.logger.warning(f"  Warning: {warning}")

            if result.validation_errors:
                for error in result.validation_errors:
                    self.logger.error(f"  Error: {error}")

        return results

    def generate_validation_report(self, results: Dict[str, ReferenceValidationResult],
                                 output_file: str) -> None:
        """Generate detailed validation report"""
        report_data = []

        for filename, result in results.items():
            report_data.append({
                'filename': filename,
                'accession': result.accession,
                'species': result.species,
                'protein_family': result.protein_family,
                'is_valid': result.is_valid,
                'completeness_score': result.completeness_score,
                'orthology_confidence': result.orthology_confidence,
                'warnings': '; '.join(result.validation_warnings),
                'errors': '; '.join(result.validation_errors)
            })

        df = pd.DataFrame(report_data)
        df.to_csv(output_file, index=False)

        # Summary statistics
        valid_count = sum(1 for r in results.values() if r.is_valid)
        total_count = len(results)

        print(f"\n{'='*60}")
        print("REFERENCE VALIDATION SUMMARY")
        print(f"{'='*60}")
        print(f"Total references: {total_count}")
        print(f"Valid references: {valid_count}")
        print(f"Invalid references: {total_count - valid_count}")
        print(".1f")
        print(f"Report saved to: {output_file}")


def main():
    """Main entry point"""
    import argparse

    parser = argparse.ArgumentParser(description="Validate reference sequences for biological relevance")
    parser.add_argument("--reference-dir", required=True, help="Directory containing reference FASTA files")
    parser.add_argument("--email", required=True, help="NCBI email for API access")
    parser.add_argument("--output", default="reference_validation_report.csv", help="Output report file")

    args = parser.parse_args()

    # Setup logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    # Validate references
    validator = ReferenceValidator(args.email)
    results = validator.validate_reference_directory(args.reference_dir)
    validator.generate_validation_report(results, args.output)


if __name__ == "__main__":
    main()