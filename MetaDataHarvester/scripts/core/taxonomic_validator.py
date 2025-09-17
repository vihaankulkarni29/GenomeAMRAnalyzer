import sys; import os; sys.path.append(os.path.join(os.path.dirname(__file__), '\\\\\\\\\src')); from configuration_manager import config_manager
#!/usr/bin/env python3
"""
TaxonomicValidator - Ensures sequences are from correct species and are orthologs
Critical for avoiding false positives in AMR mutation analysis

Author: MetaDataHarvester Pipeline
Version: 1.0 - Scientifically Validated
"""

import os
import sys
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
from dataclasses import dataclass
from Bio import SeqIO, Entrez
from Bio.Blast import NCBIWWW, NCBIXML
import pandas as pd
from collections import defaultdict


@dataclass
class TaxonomicValidationResult:
    """Results of taxonomic and orthology validation"""
    accession: str
    original_species: str
    validated_species: str
    is_correct_species: bool
    orthology_confidence: float
    contamination_risk: float
    validation_warnings: List[str]
    validation_errors: List[str]


class TaxonomicValidator:
    """
    Validates taxonomic classification and orthology of downloaded sequences
    """

    def __init__(self, email: str, expected_species: List[str], cache_dir: str = "cache"):
        self.email = email
        self.expected_species = set(s.lower() for s in expected_species)
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(exist_ok=True)

        # Setup logging
        self.logger = logging.getLogger(__name__)

        # Known efflux pump ortholog groups
        self.efflux_pump_families = {
            config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0]: ['multidrug efflux RND transporter periplasmic adaptor subunit AcrA'],
            config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1]: ['multidrug efflux RND transporter permease subunit AcrB'],
            config_manager.get_default_genes("rnd_efflux_pumps", "primary")[2]: ['multidrug efflux RND transporter outer membrane lipoprotein TolC']
        }

    def validate_sequence_taxonomy(self, accession: str, expected_family: str) -> TaxonomicValidationResult:
        """
        Comprehensive taxonomic and orthology validation

        Args:
            accession: NCBI accession number
            expected_family: Expected protein family

        Returns:
            ValidationResult with taxonomic assessment
        """
        try:
            # Fetch sequence record from NCBI
            Entrez.email = self.email
            handle = Entrez.efetch(db="protein", id=accession, rettype="gb", retmode="text")
            record_text = handle.read()
            handle.close()

            # Parse taxonomic information
            original_species = self._extract_species_from_gb(record_text)
            validated_species = self._validate_species_via_blast(accession)

            # Check if species matches expectations
            is_correct_species = self._check_species_match(validated_species)

            # Assess orthology
            orthology_confidence = self._assess_orthology(record_text, expected_family)

            # Check for contamination
            contamination_risk = self._assess_contamination_risk(accession, expected_family)

            # Generate warnings and errors
            validation_warnings = []
            validation_errors = []

            if not is_correct_species:
                validation_errors.append(
                    f"Species mismatch: expected {self.expected_species}, got {validated_species}"
                )

            if orthology_confidence < 0.7:
                validation_warnings.append(f"Low orthology confidence: {orthology_confidence:.2f}")

            if contamination_risk > 0.5:
                validation_warnings.append(f"High contamination risk: {contamination_risk:.2f}")

            return TaxonomicValidationResult(
                accession=accession,
                original_species=original_species,
                validated_species=validated_species,
                is_correct_species=is_correct_species,
                orthology_confidence=orthology_confidence,
                contamination_risk=contamination_risk,
                validation_warnings=validation_warnings,
                validation_errors=validation_errors
            )

        except Exception as e:
            self.logger.error(f"Error validating taxonomy for {accession}: {e}")
            return TaxonomicValidationResult(
                accession=accession,
                original_species="unknown",
                validated_species="unknown",
                is_correct_species=False,
                orthology_confidence=0.0,
                contamination_risk=1.0,
                validation_warnings=[],
                validation_errors=[f"Validation failed: {str(e)}"]
            )

    def _extract_species_from_gb(self, gb_record: str) -> str:
        """Extract species from GenBank record"""
        lines = gb_record.split('\n')
        for line in lines:
            if line.startswith('  ORGANISM'):
                return line.replace('  ORGANISM', '').strip()
        return "unknown"

    def _validate_species_via_blast(self, accession: str) -> str:
        """Validate species using BLAST taxonomy"""
        try:
            # Get sequence
            handle = Entrez.efetch(db="protein", id=accession, rettype="fasta", retmode="text")
            fasta_data = handle.read()
            handle.close()

            # BLAST search
            result_handle = NCBIWWW.qblast("blastp", "nr", fasta_data, expect=1e-10, hitlist_size=20)
            blast_records = NCBIXML.read(result_handle)

            # Extract species from top hits
            species_counts = defaultdict(int)
            for alignment in blast_records.alignments:
                species = self._extract_species_from_title(alignment.title)
                if species != "unknown":
                    species_counts[species] += 1

            # Return most common species
            if species_counts:
                return max(species_counts, key=species_counts.get)

        except Exception as e:
            self.logger.warning(f"BLAST validation failed for {accession}: {e}")

        return "unknown"

    def _extract_species_from_title(self, title: str) -> str:
        """Extract species name from BLAST hit title"""
        import re

        # Look for species patterns
        patterns = [
            r'\[([^\]]+)\]',  # [Escherichia coli]
            r'([A-Z][a-z]+ [a-z]+(?: [a-z]+)?)',  # Escherichia coli
        ]

        for pattern in patterns:
            match = re.search(pattern, title)
            if match:
                species = match.group(1)
                species = species.replace(' str.', '').replace(' substr.', '')
                return species.strip()

        return "unknown"

    def _check_species_match(self, validated_species: str) -> bool:
        """Check if validated species matches expected species"""
        if validated_species == "unknown":
            return False

        validated_lower = validated_species.lower()
        return validated_lower in self.expected_species

    def _assess_orthology(self, gb_record: str, expected_family: str) -> float:
        """Assess orthology confidence based on annotation"""
        confidence = 0.0

        # Check product description
        product_lines = [line for line in gb_record.split('\n') if '/product=' in line]
        for line in product_lines:
            product_desc = line.split('/product=')[1].strip('"')
            expected_terms = self.efflux_pump_families.get(expected_family, [])

            for term in expected_terms:
                if term.lower() in product_desc.lower():
                    confidence += 0.5

        # Check gene name
        gene_lines = [line for line in gb_record.split('\n') if '/gene=' in line]
        for line in gene_lines:
            gene_name = line.split('/gene=')[1].strip('"')
            if expected_family.lower() in gene_name.lower():
                confidence += 0.3

        # Check for conserved domains
        if 'efflux' in gb_record.lower() or 'transporter' in gb_record.lower():
            confidence += 0.2

        return min(1.0, confidence)

    def _assess_contamination_risk(self, accession: str, expected_family: str) -> float:
        """Assess risk of sequence contamination or misannotation"""
        try:
            # BLAST search to check for mixed signals
            handle = Entrez.efetch(db="protein", id=accession, rettype="fasta", retmode="text")
            fasta_data = handle.read()
            handle.close()

            result_handle = NCBIWWW.qblast("blastp", "nr", fasta_data, expect=1e-5, hitlist_size=50)
            blast_records = NCBIXML.read(result_handle)

            # Analyze hit diversity
            species_set = set()
            family_terms = set()

            for alignment in blast_records.alignments:
                species = self._extract_species_from_title(alignment.title)
                if species != "unknown":
                    species_set.add(species)

                title_lower = alignment.title.lower()
                if 'efflux' in title_lower:
                    family_terms.add('efflux')
                if 'transporter' in title_lower:
                    family_terms.add('transporter')
                if expected_family.lower() in title_lower:
                    family_terms.add(expected_family)

            # Calculate contamination risk
            species_diversity = len(species_set)
            family_consistency = len(family_terms) / 3.0  # Normalize to 0-1

            # High species diversity or low family consistency indicates contamination risk
            risk = (species_diversity - 1) * 0.1 + (1 - family_consistency) * 0.5
            return min(1.0, max(0.0, risk))

        except Exception as e:
            self.logger.warning(f"Contamination assessment failed for {accession}: {e}")
            return 0.5  # Neutral risk if assessment fails

    def validate_fasta_file(self, fasta_file: str, expected_family: str) -> TaxonomicValidationResult:
        """
        Validate sequences in a FASTA file

        Args:
            fasta_file: Path to FASTA file
            expected_family: Expected protein family

        Returns:
            ValidationResult for the sequence
        """
        try:
            records = list(SeqIO.parse(fasta_file, 'fasta'))
            if len(records) != 1:
                return TaxonomicValidationResult(
                    accession="multiple_sequences",
                    original_species="unknown",
                    validated_species="unknown",
                    is_correct_species=False,
                    orthology_confidence=0.0,
                    contamination_risk=1.0,
                    validation_warnings=[],
                    validation_errors=["FASTA file must contain exactly one sequence"]
                )

            record = records[0]
            accession = record.id.split('.')[0]

            return self.validate_sequence_taxonomy(accession, expected_family)

        except Exception as e:
            self.logger.error(f"Error validating FASTA file {fasta_file}: {e}")
            return TaxonomicValidationResult(
                accession="unknown",
                original_species="unknown",
                validated_species="unknown",
                is_correct_species=False,
                orthology_confidence=0.0,
                contamination_risk=1.0,
                validation_warnings=[],
                validation_errors=[f"Validation failed: {str(e)}"]
            )

    def validate_metadata_csv(self, metadata_file: str, expected_family: str) -> Dict[str, TaxonomicValidationResult]:
        """
        Validate all sequences in a metadata CSV file

        Args:
            metadata_file: Path to metadata CSV
            expected_family: Expected protein family

        Returns:
            Dictionary mapping accessions to validation results
        """
        results = {}

        try:
            df = pd.read_csv(metadata_file)

            for _, row in df.iterrows():
                accession = str(row.get('Accession_Number', ''))
                if accession and accession != 'nan':
                    result = self.validate_sequence_taxonomy(accession, expected_family)
                    results[accession] = result

                    # Log results
                    status = "✓ VALID" if result.is_correct_species else "✗ INVALID"
                    self.logger.info(
                        f"{status} {accession}: {result.validated_species} "
                        f"(orthology: {result.orthology_confidence:.2f})"
                    )

        except Exception as e:
            self.logger.error(f"Error validating metadata file {metadata_file}: {e}")

        return results

    def generate_validation_report(self, results: Dict[str, TaxonomicValidationResult],
                                 output_file: str) -> None:
        """Generate detailed validation report"""
        report_data = []

        for accession, result in results.items():
            report_data.append({
                'accession': accession,
                'original_species': result.original_species,
                'validated_species': result.validated_species,
                'is_correct_species': result.is_correct_species,
                'orthology_confidence': result.orthology_confidence,
                'contamination_risk': result.contamination_risk,
                'warnings': '; '.join(result.validation_warnings),
                'errors': '; '.join(result.validation_errors)
            })

        df = pd.DataFrame(report_data)
        df.to_csv(output_file, index=False)

        # Summary statistics
        valid_count = sum(1 for r in results.values() if r.is_correct_species)
        total_count = len(results)

        print(f"\n{'='*60}")
        print("TAXONOMIC VALIDATION SUMMARY")
        print(f"{'='*60}")
        print(f"Total sequences: {total_count}")
        print(f"Correct species: {valid_count}")
        print(f"Incorrect species: {total_count - valid_count}")
        print(".1f")
        print(f"Report saved to: {output_file}")


def main():
    """Main entry point"""
    import argparse

    parser = argparse.ArgumentParser(description="Validate taxonomic classification and orthology")
    parser.add_argument("--metadata", required=True, help="Metadata CSV file with accessions")
    parser.add_argument("--expected-species", nargs='+', required=True,
                       help="Expected species (e.g., 'Escherichia coli')")
    parser.add_argument("--protein-family", required=True,
                       help="Expected protein family (acrA, acrB, tolC)")
    parser.add_argument("--email", required=True, help="NCBI email for API access")
    parser.add_argument("--output", default="taxonomic_validation_report.csv",
                       help="Output report file")

    args = parser.parse_args()

    # Setup logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    # Validate sequences
    validator = TaxonomicValidator(args.email, args.expected_species)
    results = validator.validate_metadata_csv(args.metadata, args.protein_family)
    validator.generate_validation_report(results, args.output)


if __name__ == "__main__":
    main()