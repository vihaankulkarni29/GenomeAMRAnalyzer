from .configuration_manager import config_manager
#!/usr/bin/env python3
"""
SimplifiedCARDIntegrator - CARD database integration without pandas dependency
"""

import os
import sys
import logging
import requests
import json
import csv
from pathlib import Path
from typing import Dict, List, Optional
from dataclasses import dataclass
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


class SimplifiedCARDIntegrator:
    """
    Simplified CARD database integrator without pandas dependency
    """

    def __init__(self, cache_dir: str = "card_cache", offline_mode: bool = True):
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

        # Initialize with some mock CARD data for testing
        if not self.card_data:
            self._initialize_mock_data()

        # Resistance gene mapping
        self.resistance_genes = {
            config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0]: [config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0], 'AcrA'],
            config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1]: [config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1], 'AcrB'],
            config_manager.get_default_genes("rnd_efflux_pumps", "primary")[2]: [config_manager.get_default_genes("rnd_efflux_pumps", "primary")[2], 'TolC'],
            'marA': ['marA', 'MarA'],
            'soxS': ['soxS', 'SoxS'],
            'ramA': ['ramA', 'RamA']
        }

    def _initialize_mock_data(self):
        """Initialize mock CARD data for testing"""
        self.card_data = {
            config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1]: [
                {
                    'aro_id': 'ARO:3000001',
                    'gene_name': config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1],
                    'mutation': 'S83L',
                    'drug_class': 'fluoroquinolone',
                    'resistance_mechanism': 'efflux pump overexpression',
                    'evidence_level': 'Strong',
                    'publications': ['PMID:12345678']
                },
                {
                    'aro_id': 'ARO:3000002',
                    'gene_name': config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1],
                    'mutation': 'D87N',
                    'drug_class': 'fluoroquinolone',
                    'resistance_mechanism': 'efflux pump overexpression',
                    'evidence_level': 'Moderate',
                    'publications': ['PMID:12345679']
                }
            ],
            config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0]: [
                {
                    'aro_id': 'ARO:3000003',
                    'gene_name': config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0],
                    'mutation': 'A119E',
                    'drug_class': 'beta-lactam',
                    'resistance_mechanism': 'efflux pump component',
                    'evidence_level': 'Weak',
                    'publications': ['PMID:12345680']
                }
            ]
        }
        self._save_card_cache()

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
            # Return mock data from cache
            raw_data = self.card_data.get(gene_name, [])
            return [CARDMutation(**data) if isinstance(data, dict) else data for data in raw_data]

        # In online mode, you would implement real CARD API calls here
        self.logger.warning("Online mode not implemented in simplified version")
        return []

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

        # Look for exact matches
        for card_mut in card_mutations:
            if card_mut.mutation == mutation:
                matching_entries.append(card_mut)
                confidence_score = 0.9  # High confidence for exact match

        # Look for partial matches if no exact match
        if not matching_entries:
            for card_mut in card_mutations:
                if mutation in card_mut.mutation or card_mut.mutation in mutation:
                    matching_entries.append(card_mut)
                    confidence_score = max(confidence_score, 0.6)  # Medium confidence

        # Determine clinical significance
        if matching_entries:
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

    def validate_mutation_list(self, mutations: List[str], gene_name: str) -> List[CARDValidationResult]:
        """Validate a list of mutations"""
        results = []
        for mutation in mutations:
            result = self.validate_mutation(mutation, gene_name)
            results.append(result)
        return results

    def validate_csv_file(self, csv_file: str, gene_name: str, mutation_column: str = "mutation") -> List[CARDValidationResult]:
        """Validate mutations from CSV file without pandas"""
        results = []
        
        try:
            with open(csv_file, 'r') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    if mutation_column in row:
                        mutation = row[mutation_column]
                        result = self.validate_mutation(mutation, gene_name)
                        results.append(result)
                        
        except Exception as e:
            self.logger.error(f"Error reading CSV file {csv_file}: {e}")
            
        return results

    def generate_validation_report(self, results: List[CARDValidationResult], output_file: str):
        """Generate validation report without pandas"""
        try:
            with open(output_file, 'w', newline='') as f:
                fieldnames = [
                    'mutation', 'is_known_resistance', 'clinical_significance',
                    'confidence_score', 'card_entries_count', 'validation_warnings'
                ]
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
                
                for result in results:
                    writer.writerow({
                        'mutation': result.mutation,
                        'is_known_resistance': result.is_known_resistance,
                        'clinical_significance': result.clinical_significance,
                        'confidence_score': result.confidence_score,
                        'card_entries_count': len(result.card_entries),
                        'validation_warnings': '; '.join(result.validation_warnings)
                    })
            
            # Print summary
            total_mutations = len(results)
            known_resistance = sum(1 for r in results if r.is_known_resistance)
            
            print(f"\nCARD Validation Summary:")
            print(f"Total mutations: {total_mutations}")
            print(f"Known resistance mutations: {known_resistance}")
            print(f"Novel mutations: {total_mutations - known_resistance}")
            print(f"Report saved to: {output_file}")
            
        except Exception as e:
            self.logger.error(f"Error generating report: {e}")

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

    parser = argparse.ArgumentParser(description="Validate mutations against CARD database (simplified)")
    parser.add_argument("--csv-file", help="CSV file containing mutations")
    parser.add_argument("--mutation-column", default="mutation", help="Column with mutations")
    parser.add_argument("--gene-name", required=True, help="Gene name (acrA, acrB, etc.)")
    parser.add_argument("--output", default="card_validation_report.csv", help="Output report file")

    args = parser.parse_args()

    # Setup logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    # Validate mutations
    integrator = SimplifiedCARDIntegrator(offline_mode=True)

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