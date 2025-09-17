#!/usr/bin/env python3
"""
MICIntegrator - Minimum Inhibitory Concentration data collection and integration
Links AMR mutations with clinical MIC data for resistance phenotype correlation

Author: MetaDataHarvester Pipeline
Version: 2.0 - MIC Integration
"""

import os
import sys
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
from dataclasses import dataclass
import pandas as pd
import numpy as np
from Bio import Entrez
import requests
import json
import time
from collections import defaultdict


@dataclass
class MICData:
    """MIC data for a specific antibiotic-organism combination"""
    antibiotic: str
    organism: str
    mic_value: float
    mic_unit: str
    breakpoint_s: Optional[float]  # Susceptible breakpoint
    breakpoint_r: Optional[float]  # Resistant breakpoint
    source: str
    study_id: Optional[str]
    confidence: float

    @property
    def resistance_category(self) -> str:
        """Determine resistance category based on MIC and breakpoints"""
        if self.breakpoint_s is None or self.breakpoint_r is None:
            return "Unknown"

        if self.mic_value <= self.breakpoint_s:
            return "Susceptible"
        elif self.mic_value >= self.breakpoint_r:
            return "Resistant"
        else:
            return "Intermediate"


@dataclass
class MICCorrelation:
    """Correlation between mutation and MIC data"""
    mutation: str
    antibiotic: str
    correlation_coefficient: float
    p_value: float
    sample_size: int
    effect_size: float
    confidence_interval: Tuple[float, float]


class MICIntegrator:
    """
    Integrates MIC data from various sources to correlate with AMR mutations
    """

    def __init__(self, email: str, cache_dir: str = "mic_cache"):
        self.email = email
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(exist_ok=True)

        # Setup logging
        self.logger = logging.getLogger(__name__)
        self._setup_logging()

        # Initialize NCBI
        Entrez.email = email

        # MIC data sources
        self.data_sources = {
            'EUCAST': self._fetch_eucast_data,
            'CLSI': self._fetch_clsi_data,
            'CARD': self._fetch_card_mic_data,
            'PATRIC': self._fetch_patric_mic_data,
            'PubMed': self._fetch_pubmed_mic_data
        }

    def _setup_logging(self):
        """Setup logging configuration"""
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )

    def collect_mic_data(self, organism: str, antibiotic: str) -> List[MICData]:
        """
        Collect MIC data for specific organism-antibiotic combination

        Args:
            organism: Bacterial species name
            antibiotic: Antibiotic name

        Returns:
            List of MIC data points
        """
        self.logger.info(f"Collecting MIC data for {organism} vs {antibiotic}")

        mic_data = []

        # Try each data source
        for source_name, fetch_function in self.data_sources.items():
            try:
                self.logger.debug(f"Trying {source_name}...")
                source_data = fetch_function(organism, antibiotic)
                if source_data:
                    mic_data.extend(source_data)
                    self.logger.info(f"Found {len(source_data)} MIC values from {source_name}")
            except Exception as e:
                self.logger.warning(f"Failed to fetch from {source_name}: {e}")

        # Remove duplicates and validate
        mic_data = self._deduplicate_mic_data(mic_data)
        mic_data = self._validate_mic_data(mic_data)

        self.logger.info(f"Collected {len(mic_data)} validated MIC data points")
        return mic_data

    def _fetch_eucast_data(self, organism: str, antibiotic: str) -> List[MICData]:
        """Fetch MIC data from EUCAST database"""
        # EUCAST API integration
        # This would connect to EUCAST's breakpoint database
        return []  # Placeholder

    def _fetch_clsi_data(self, organism: str, antibiotic: str) -> List[MICData]:
        """Fetch MIC data from CLSI database"""
        # CLSI API integration
        return []  # Placeholder

    def _fetch_card_mic_data(self, organism: str, antibiotic: str) -> List[MICData]:
        """Fetch MIC data from CARD database"""
        try:
            # Query CARD for MIC data associated with resistance genes
            card_url = "https://card.mcmaster.ca/api/v1/genes"
            params = {
                'search': f'{organism} {antibiotic}',
                'limit': 100
            }

            response = requests.get(card_url, params=params, timeout=30)
            if response.status_code == 200:
                data = response.json()
                return self._parse_card_mic_data(data, organism, antibiotic)

        except Exception as e:
            self.logger.warning(f"CARD MIC fetch failed: {e}")

        return []

    def _fetch_patric_mic_data(self, organism: str, antibiotic: str) -> List[MICData]:
        """Fetch MIC data from PATRIC database"""
        try:
            # PATRIC API for antimicrobial resistance data
            patric_url = "https://www.patricbrc.org/api/genome/"
            # Implementation would query PATRIC's AMR data
            return []
        except Exception as e:
            self.logger.warning(f"PATRIC MIC fetch failed: {e}")
        return []

    def _fetch_pubmed_mic_data(self, organism: str, antibiotic: str) -> List[MICData]:
        """Fetch MIC data from PubMed studies"""
        try:
            # Search PubMed for MIC studies
            query = f'"{organism}"[Organism] AND "{antibiotic}"[Title/Abstract] AND MIC[Title/Abstract]'
            handle = Entrez.esearch(db="pubmed", term=query, retmax=50)
            search_results = Entrez.read(handle)
            handle.close()

            pmids = search_results.get('IdList', [])
            if pmids:
                return self._extract_mic_from_abstracts(pmids, organism, antibiotic)

        except Exception as e:
            self.logger.warning(f"PubMed MIC fetch failed: {e}")

        return []

    def _parse_card_mic_data(self, data: Dict, organism: str, antibiotic: str) -> List[MICData]:
        """Parse MIC data from CARD API response"""
        mic_data = []

        # Parse CARD's JSON response for MIC information
        # This would extract MIC values from CARD's resistance gene annotations

        return mic_data

    def _extract_mic_from_abstracts(self, pmids: List[str], organism: str, antibiotic: str) -> List[MICData]:
        """Extract MIC values from PubMed abstracts"""
        mic_data = []

        try:
            # Fetch abstracts
            handle = Entrez.efetch(db="pubmed", id=pmids, rettype="abstract", retmode="text")
            abstracts = handle.read()
            handle.close()

            # Parse MIC values from text
            # Use regex to find MIC values like "MIC = 8 μg/mL"
            import re
            mic_pattern = r'MIC\s*[=<>]\s*(\d+(?:\.\d+)?)\s*(μg/mL|mg/L|ug/mL)'
            matches = re.findall(mic_pattern, abstracts, re.IGNORECASE)

            for match in matches:
                try:
                    value, unit = match
                    mic_data.append(MICData(
                        antibiotic=antibiotic,
                        organism=organism,
                        mic_value=float(value),
                        mic_unit=unit,
                        breakpoint_s=None,  # Would need to lookup from standards
                        breakpoint_r=None,
                        source="PubMed",
                        study_id=pmids[0] if pmids else None,
                        confidence=0.5  # Lower confidence for text extraction
                    ))
                except ValueError:
                    continue

        except Exception as e:
            self.logger.warning(f"Abstract MIC extraction failed: {e}")

        return mic_data

    def _deduplicate_mic_data(self, mic_data: List[MICData]) -> List[MICData]:
        """Remove duplicate MIC entries"""
        seen = set()
        unique_data = []

        for mic in mic_data:
            key = (mic.antibiotic, mic.organism, mic.mic_value, mic.source)
            if key not in seen:
                seen.add(key)
                unique_data.append(mic)

        return unique_data

    def _validate_mic_data(self, mic_data: List[MICData]) -> List[MICData]:
        """Validate MIC data quality"""
        validated_data = []

        for mic in mic_data:
            if self._is_valid_mic_value(mic):
                validated_data.append(mic)

        return validated_data

    def _is_valid_mic_value(self, mic: MICData) -> bool:
        """Check if MIC value is biologically reasonable"""
        # Basic validation rules
        if mic.mic_value <= 0:
            return False

        # Check for reasonable ranges based on antibiotic
        antibiotic_ranges = {
            'ciprofloxacin': (0.001, 128),
            'tetracycline': (0.25, 256),
            'ampicillin': (0.5, 512),
            'gentamicin': (0.25, 256),
            'cefotaxime': (0.06, 512)
        }

        if mic.antibiotic.lower() in antibiotic_ranges:
            min_val, max_val = antibiotic_ranges[mic.antibiotic.lower()]
            if min_val <= mic.mic_value <= max_val:
                return True

        # Default range check
        return 0.001 <= mic.mic_value <= 1024

    def correlate_mutations_with_mic(self, mutations_df: pd.DataFrame,
                                   mic_data: List[MICData]) -> List[MICCorrelation]:
        """
        Correlate mutations with MIC data

        Args:
            mutations_df: DataFrame with mutation data
            mic_data: List of MIC data points

        Returns:
            List of mutation-MIC correlations
        """
        correlations = []

        # Group MIC data by antibiotic
        mic_by_antibiotic = defaultdict(list)
        for mic in mic_data:
            mic_by_antibiotic[mic.antibiotic].append(mic)

        # For each mutation, calculate correlation with MIC
        for _, mutation_row in mutations_df.iterrows():
            mutation_key = f"{mutation_row['reference_aa']}{mutation_row['position']}{mutation_row['variant_aa']}"

            for antibiotic, mic_list in mic_by_antibiotic.items():
                try:
                    correlation = self._calculate_mutation_mic_correlation(
                        mutation_row, mic_list, antibiotic
                    )
                    if correlation:
                        correlations.append(correlation)
                except Exception as e:
                    self.logger.warning(f"Correlation calculation failed for {mutation_key}: {e}")

        return correlations

    def _calculate_mutation_mic_correlation(self, mutation_row: pd.Series,
                                          mic_list: List[MICData],
                                          antibiotic: str) -> Optional[MICCorrelation]:
        """Calculate correlation between mutation presence and MIC values"""
        try:
            # This would implement statistical correlation analysis
            # For now, return a placeholder
            return MICCorrelation(
                mutation=f"{mutation_row['reference_aa']}{mutation_row['position']}{mutation_row['variant_aa']}",
                antibiotic=antibiotic,
                correlation_coefficient=0.0,
                p_value=1.0,
                sample_size=len(mic_list),
                effect_size=0.0,
                confidence_interval=(0.0, 0.0)
            )
        except Exception:
            return None

    def generate_mic_report(self, mic_data: List[MICData],
                          correlations: List[MICCorrelation]) -> None:
        """Generate comprehensive MIC integration report"""
        report_file = self.cache_dir / "mic_integration_report.txt"

        with open(report_file, 'w') as f:
            f.write("MIC DATA INTEGRATION REPORT\n")
            f.write("=" * 40 + "\n\n")

            # MIC data summary
            f.write("MIC DATA SUMMARY:\n")
            f.write("-" * 20 + "\n")
            f.write(f"Total MIC data points: {len(mic_data)}\n")

            if mic_data:
                antibiotics = set(m.antibiotic for m in mic_data)
                organisms = set(m.organism for m in mic_data)
                sources = set(m.source for m in mic_data)

                f.write(f"Antibiotics: {', '.join(antibiotics)}\n")
                f.write(f"Organisms: {', '.join(organisms)}\n")
                f.write(f"Data sources: {', '.join(sources)}\n\n")

                # MIC distribution
                mic_values = [m.mic_value for m in mic_data]
                f.write("MIC VALUE DISTRIBUTION:\n")
                f.write("-" * 25 + "\n")
                f.write(f"Mean MIC: {np.mean(mic_values):.2f}\n")
                f.write(f"Median MIC: {np.median(mic_values):.2f}\n")
                f.write(f"Min MIC: {np.min(mic_values):.2f}\n")
                f.write(f"Max MIC: {np.max(mic_values):.2f}\n\n")

            # Correlation summary
            f.write("MUTATION-MIC CORRELATIONS:\n")
            f.write("-" * 30 + "\n")
            f.write(f"Total correlations analyzed: {len(correlations)}\n")

            if correlations:
                significant_correlations = [c for c in correlations if c.p_value < 0.05]
                f.write(f"Significant correlations (p<0.05): {len(significant_correlations)}\n\n")

                # Top correlations
                sorted_correlations = sorted(correlations,
                                           key=lambda x: abs(x.correlation_coefficient),
                                           reverse=True)

                f.write("TOP CORRELATIONS:\n")
                f.write("-" * 20 + "\n")
                for i, corr in enumerate(sorted_correlations[:10]):
                    f.write(f"{i+1}. {corr.mutation} vs {corr.antibiotic}: ")
                    f.write(".3f")
                    f.write(f" (p={corr.p_value:.3f}, n={corr.sample_size})\n")

        self.logger.info(f"MIC integration report saved: {report_file}")


def main():
    """Main entry point"""
    import argparse

    parser = argparse.ArgumentParser(
        description="MIC data collection and integration system"
    )

    parser.add_argument(
        "--organism",
        required=True,
        help="Bacterial organism name"
    )

    parser.add_argument(
        "--antibiotic",
        required=True,
        help="Antibiotic name"
    )

    parser.add_argument(
        "--email",
        required=True,
        help="NCBI email for API access"
    )

    parser.add_argument(
        "--output-dir",
        default="mic_data",
        help="Output directory for MIC data"
    )

    args = parser.parse_args()

    # Setup logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    # Initialize MIC integrator
    integrator = MICIntegrator(args.email, args.output_dir)

    # Collect MIC data
    mic_data = integrator.collect_mic_data(args.organism, args.antibiotic)

    # Generate report
    integrator.generate_mic_report(mic_data, [])

    print(f"\nMIC data collection complete!")
    print(f"Collected {len(mic_data)} MIC data points")
    print(f"Results saved to: {args.output_dir}")


if __name__ == "__main__":
    main()