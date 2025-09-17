import sys; import os; sys.path.append(os.path.join(os.path.dirname(__file__), '\\\\\\\\\src')); from configuration_manager import config_manager
#!/usr/bin/env python3
"""
CARDIntegrator - Integration with CARD database for known resistance mutations
Critical for validating AMR research findings against established knowledge

Author: MetaDataHarvester Pipeline
Version: 1.0 - Scientifically Validated
"""

import os
import sys
import logging
import requests
import json
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
from dataclasses import dataclass
import pandas as pd
from Bio import SeqIO
import time
from urllib.parse import urljoin


@dataclass
class CARDMutation:
    """CARD database mutation entry"""
    aro_id: str
    gene_name: str
    mutation: str
    drug_class: str
    resistance_mechanism: str
    evidence_level: str
    publications: List[str]


@dataclass
class CARDValidationResult:
    """Results of CARD database validation"""
    mutation: str
    is_known_resistance: bool
    card_entries: List[CARDMutation]
    confidence_score: float
    clinical_significance: str
    validation_warnings: List[str]


class CARDIntegrator:
    """
    Integrates with CARD database to validate mutations against known resistance
    """

    def __init__(self, cache_dir: str = "card_cache", offline_mode: bool = False):
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(exist_ok=True)
        self.offline_mode = offline_mode

        # Setup logging
        self.logger = logging.getLogger(__name__)

        # CARD API endpoints
        self.base_url = "https://card.mcmaster.ca/api/v1/"

        # Cache for CARD data
        self.card_cache_file = self.cache_dir / "card_data.json"
        self.card_data = self._load_card_cache()

        # Resistance gene mapping
        self.resistance_genes = {
            config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0]: [config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0], 'AcrA'],
            config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1]: [config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1], 'AcrB'],
            config_manager.get_default_genes("rnd_efflux_pumps", "primary")[2]: [config_manager.get_default_genes("rnd_efflux_pumps", "primary")[2], 'TolC'],
            'marA': ['marA', 'MarA'],
            'soxS': ['soxS', 'SoxS'],
            'ramA': ['ramA', 'RamA']
        }

    def _load_card_cache(self) -> Dict:
        """Load cached CARD data"""
        if self.card_cache_file.exists():
            try:
                with open(self.card_cache_file, 'r') as f:
                    return json.load(f)
            except Exception as e:
                self.logger.warning(f"Could not load CARD cache: {e}")

        return {}

    def _save_card_cache(self) -> None:
        """Save CARD data to cache"""
        try:
            with open(self.card_cache_file, 'w') as f:
                json.dump(self.card_data, f, indent=2)
        except Exception as e:
            self.logger.warning(f"Could not save CARD cache: {e}")

    def fetch_card_data(self, gene_name: str) -> List[CARDMutation]:
        """
        Fetch CARD data for a specific gene

        Args:
            gene_name: Gene name to search for

        Returns:
            List of CARD mutations for the gene
        """
        if self.offline_mode:
            return self.card_data.get(gene_name, [])

        try:
            # Search for gene in CARD
            search_url = urljoin(self.base_url, "aro/search")
            params = {
                'q': gene_name,
                'limit': 100
            }

            response = requests.get(search_url, params=params, timeout=30)
            response.raise_for_status()

            search_results = response.json()

            mutations = []
            for result in search_results.get('results', []):
                aro_id = result.get('aro_id', '')

                # Get detailed information for this ARO
                detail_url = urljoin(self.base_url, f"aro/{aro_id}")
                detail_response = requests.get(detail_url, timeout=30)

                if detail_response.status_code == 200:
                    detail_data = detail_response.json()

                    # Extract resistance information
                    resistance_info = self._extract_resistance_info(detail_data)
                    mutations.extend(resistance_info)

            # Cache results
            self.card_data[gene_name] = mutations
            self._save_card_cache()

            return mutations

        except Exception as e:
            self.logger.error(f"Error fetching CARD data for {gene_name}: {e}")
            return self.card_data.get(gene_name, [])

    def _extract_resistance_info(self, card_data: Dict) -> List[CARDMutation]:
        """Extract resistance information from CARD entry"""
        mutations = []

        try:
            aro_id = card_data.get('aro_id', '')
            gene_name = card_data.get('name', '')

            # Look for resistance information
            resistance_data = card_data.get('resistance', {})

            for drug_class, mechanisms in resistance_data.items():
                for mechanism, details in mechanisms.items():
                    if 'mutations' in details:
                        for mutation in details['mutations']:
                            card_mutation = CARDMutation(
                                aro_id=aro_id,
                                gene_name=gene_name,
                                mutation=mutation,
                                drug_class=drug_class,
                                resistance_mechanism=mechanism,
                                evidence_level=details.get('evidence_level', 'Unknown'),
                                publications=details.get('publications', [])
                            )
                            mutations.append(card_mutation)

        except Exception as e:
            self.logger.warning(f"Error extracting resistance info: {e}")

        return mutations

    def validate_mutation(self, mutation: str, gene_name: str) -> CARDValidationResult:
        """
        Validate a mutation against CARD database

        Args:
            mutation: Mutation in format "A123T" (ref_pos_var)
            gene_name: Gene name (acrA, acrB, etc.)

        Returns:
            Validation result with CARD information
        """
        # Fetch CARD data for this gene
        card_mutations = self.fetch_card_data(gene_name)

        # Look for matching mutations
        matching_entries = []
        confidence_score = 0.0
        clinical_significance = "Unknown"

        # Parse mutation
        try:
            if len(mutation) >= 4:  # At least ref + pos + var
                ref_aa = mutation[0]
                pos = ''.join(c for c in mutation[1:-1] if c.isdigit())
                var_aa = mutation[-1]

                if pos:
                    pos_int = int(pos)

                    # Look for exact matches
                    for card_mut in card_mutations:
                        if self._mutation_matches(card_mut.mutation, ref_aa, pos_int, var_aa):
                            matching_entries.append(card_mut)
                            confidence_score = max(confidence_score, 0.9)  # High confidence for exact match

                    # Look for position-based matches (less specific)
                    if not matching_entries:
                        for card_mut in card_mutations:
                            if str(pos_int) in card_mut.mutation:
                                matching_entries.append(card_mut)
                                confidence_score = max(confidence_score, 0.6)  # Medium confidence

        except Exception as e:
            self.logger.warning(f"Error parsing mutation {mutation}: {e}")

        # Determine clinical significance
        if matching_entries:
            # Check evidence levels
            evidence_levels = [entry.evidence_level for entry in matching_entries]
            if any(level in ['Strong', 'Moderate'] for level in evidence_levels):
                clinical_significance = "Established Resistance"
            elif any(level == 'Weak' for level in evidence_levels):
                clinical_significance = "Potential Resistance"
            else:
                clinical_significance = "Reported"
        else:
            clinical_significance = "Novel"

        # Generate warnings
        validation_warnings = []
        if not matching_entries:
            validation_warnings.append("Mutation not found in CARD database")
        if confidence_score < 0.7:
            validation_warnings.append("Low confidence in CARD matching")

        return CARDValidationResult(
            mutation=mutation,
            is_known_resistance=len(matching_entries) > 0,
            card_entries=matching_entries,
            confidence_score=confidence_score,
            clinical_significance=clinical_significance,
            validation_warnings=validation_warnings
        )

    def _mutation_matches(self, card_mutation: str, ref_aa: str, pos: int, var_aa: str) -> bool:
        """Check if CARD mutation matches our mutation"""
        try:
            # Normalize CARD mutation format
            card_mut_norm = card_mutation.replace(' ', '').upper()

            # Check various formats
            patterns = [
                f"{ref_aa}{pos}{var_aa}",  # A123T
                f"{ref_aa}{pos}{var_aa}",  # a123t
                f"{ref_aa}@{pos}@{var_aa}",  # Some CARD formats
            ]

            for pattern in patterns:
                if pattern.upper() in card_mut_norm:
                    return True

            return False

        except Exception:
            return False

    def validate_mutation_list(self, mutations: List[str], gene_name: str) -> Dict[str, CARDValidationResult]:
        """
        Validate a list of mutations against CARD

        Args:
            mutations: List of mutations in format ["A123T", "G456D", ...]
            gene_name: Gene name

        Returns:
            Dictionary mapping mutations to validation results
        """
        results = {}

        for mutation in mutations:
            result = self.validate_mutation(mutation, gene_name)
            results[mutation] = result

            # Rate limiting for API calls
            time.sleep(0.1)

        return results

    def validate_csv_file(self, csv_file: str, gene_name: str,
                         mutation_column: str = "mutation") -> Dict[str, CARDValidationResult]:
        """
        Validate mutations from a CSV file

        Args:
            csv_file: Path to CSV file
            gene_name: Gene name
            mutation_column: Column name containing mutations

        Returns:
            Dictionary mapping mutations to validation results
        """
        try:
            df = pd.read_csv(csv_file)
            mutations = df[mutation_column].dropna().unique().tolist()
            return self.validate_mutation_list(mutations, gene_name)

        except Exception as e:
            self.logger.error(f"Error validating CSV file {csv_file}: {e}")
            return {}

    def generate_validation_report(self, results: Dict[str, CARDValidationResult],
                                 output_file: str) -> None:
        """Generate detailed CARD validation report"""
        report_data = []

        for mutation, result in results.items():
            # Get drug classes and mechanisms
            drug_classes = list(set(entry.drug_class for entry in result.card_entries))
            mechanisms = list(set(entry.resistance_mechanism for entry in result.card_entries))
            publications = list(set(pub for entry in result.card_entries for pub in entry.publications))

            report_data.append({
                'mutation': mutation,
                'is_known_resistance': result.is_known_resistance,
                'confidence_score': result.confidence_score,
                'clinical_significance': result.clinical_significance,
                'drug_classes': '; '.join(drug_classes),
                'resistance_mechanisms': '; '.join(mechanisms),
                'publications': '; '.join(publications[:5]),  # Limit to first 5
                'warnings': '; '.join(result.validation_warnings)
            })

        df = pd.DataFrame(report_data)
        df.to_csv(output_file, index=False)

        # Summary statistics
        known_resistance = sum(1 for r in results.values() if r.is_known_resistance)
        total_mutations = len(results)

        print(f"\n{'='*60}")
        print("CARD DATABASE VALIDATION SUMMARY")
        print(f"{'='*60}")
        print(f"Total mutations: {total_mutations}")
        print(f"Known resistance mutations: {known_resistance}")
        print(".1f")
        print(f"Report saved to: {output_file}")

    def get_resistance_summary(self, gene_name: str) -> Dict[str, List[str]]:
        """
        Get summary of resistance mutations for a gene from CARD

        Args:
            gene_name: Gene name

        Returns:
            Dictionary with resistance information
        """
        card_mutations = self.fetch_card_data(gene_name)

        summary = {
            'drug_classes': [],
            'resistance_mechanisms': [],
            'mutations': []
        }

        for mutation in card_mutations:
            if mutation.drug_class not in summary['drug_classes']:
                summary['drug_classes'].append(mutation.drug_class)
            if mutation.resistance_mechanism not in summary['resistance_mechanisms']:
                summary['resistance_mechanisms'].append(mutation.resistance_mechanism)
            if mutation.mutation not in summary['mutations']:
                summary['mutations'].append(mutation.mutation)

        return summary


def main():
    """Main entry point"""
    import argparse

    parser = argparse.ArgumentParser(description="Validate mutations against CARD database")
    parser.add_argument("--csv-file", help="CSV file containing mutations")
    parser.add_argument("--mutation-column", default="mutation", help="Column with mutations")
    parser.add_argument("--gene-name", required=True, help="Gene name (acrA, acrB, etc.)")
    parser.add_argument("--output", default="card_validation_report.csv", help="Output report file")
    parser.add_argument("--offline", action="store_true", help="Run in offline mode")

    args = parser.parse_args()

    # Setup logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    # Validate mutations
    integrator = CARDIntegrator(offline_mode=args.offline)

    if args.csv_file:
        results = integrator.validate_csv_file(args.csv_file, args.gene_name, args.mutation_column)
    else:
        # Interactive mode
        mutations = []
        print("Enter mutations (one per line, empty line to finish):")
        while True:
            mutation = input().strip()
            if not mutation:
                break
            mutations.append(mutation)

        results = integrator.validate_mutation_list(mutations, args.gene_name)

    integrator.generate_validation_report(results, args.output)


if __name__ == "__main__":
    main()