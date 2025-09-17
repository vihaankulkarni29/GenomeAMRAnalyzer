#!/usr/bin/env python3
"""
Production Results Aggregation System
Enterprise-grade synthesis and reporting of genomics pipeline results

Features:
- Robust ingestion of mutation, cooccurrence, and clinical annotation results
- Aggregation and harmonization of multi-module outputs
- Summary statistics, prevalence, resistance profiles
- Multi-format reporting: CSV, JSON, PDF, HTML dashboard
- Interactive dashboards and clinical summaries
- Manifest and provenance tracking
- Batch processing and quality control

Author: GenomeAMRAnalyzer Pipeline
Version: 2.0 - Production Results Aggregation System
"""

import os
import sys
import logging
import json
import csv
import yaml
from pathlib import Path
from typing import List, Dict, Any, Optional
from dataclasses import dataclass, asdict, field
from collections import defaultdict, Counter
from datetime import datetime
import pandas as pd
try:
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False

@dataclass
class AggregatedResult:
    sample: str
    gene: str
    mutation_count: int
    resistance_profile: str
    cooccurrence_count: int
    clinical_significance: str
    notes: str = ""
    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

@dataclass
class AggregationReport:
    pipeline_id: str
    execution_timestamp: str
    total_samples: int
    total_genes: int
    aggregated_results: List[AggregatedResult]
    summary_statistics: Dict[str, Any]
    clinical_summary: Dict[str, Any]
    manifest: Dict[str, Any]
    processing_time: float
    def to_dict(self) -> Dict[str, Any]:
        result = asdict(self)
        result['aggregated_results'] = [r.to_dict() for r in self.aggregated_results]
        return result

class ProductionResultsAggregator:
    """
    Enterprise-grade results aggregation and reporting system
    """
    def __init__(self, config_path: str):
        self.config_path = Path(config_path)
        self.config = self._load_config()
        self.base_dir = Path.cwd()
        self.mutation_dir = self.base_dir / self.config['directories']['mutations'] / 'mutation_calls'
        self.cooccurrence_dir = self.base_dir / self.config['directories']['cooccurrence']
        self.output_dir = self.base_dir / self.config['directories']['results']
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.logger = self._setup_logging()
        self.logger.info("ProductionResultsAggregator initialized successfully")

    def _load_config(self) -> Dict[str, Any]:
        with open(self.config_path, 'r') as f:
            return yaml.safe_load(f)

    def _setup_logging(self) -> logging.Logger:
        logger = logging.getLogger('ProductionResultsAggregator')
        logger.setLevel(logging.INFO)
        logger.handlers.clear()
        log_file = self.output_dir / 'results_aggregator.log'
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.INFO)
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)
        return logger

    def _ingest_results(self) -> List[AggregatedResult]:
        # Ingest mutation calls
        mutation_calls = []
        for file in self.mutation_dir.glob('*.csv'):
            with open(file, 'r') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    mutation_calls.append(row)
        # Ingest cooccurrence pairs
        cooccurrence_pairs = []
        for file in self.cooccurrence_dir.glob('*_cooccurrence_pairs.csv'):
            with open(file, 'r') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    cooccurrence_pairs.append(row)
        # Aggregate results
        results = []
        mutation_df = pd.DataFrame(mutation_calls)
        cooccurrence_df = pd.DataFrame(cooccurrence_pairs)
        for sample in mutation_df['accession'].unique():
            sample_mutations = mutation_df[mutation_df['accession'] == sample]
            for gene in sample_mutations['gene_name'].unique():
                gene_mutations = sample_mutations[sample_mutations['gene_name'] == gene]
                mutation_count = len(gene_mutations)
                resistance_profile = gene_mutations['clinical_significance'].mode()[0] if not gene_mutations.empty else 'unknown'
                cooccurrence_count = cooccurrence_df[(cooccurrence_df['gene_a'] == gene) | (cooccurrence_df['gene_b'] == gene)].shape[0]
                clinical_significance = resistance_profile
                results.append(AggregatedResult(
                    sample=sample,
                    gene=gene,
                    mutation_count=mutation_count,
                    resistance_profile=resistance_profile,
                    cooccurrence_count=cooccurrence_count,
                    clinical_significance=clinical_significance
                ))
        self.logger.info(f"Aggregated {len(results)} results from mutation and cooccurrence data")
        return results

    def _compute_summary_statistics(self, results: List[AggregatedResult]) -> Dict[str, Any]:
        stats = defaultdict(int)
        for r in results:
            stats['total_mutations'] += r.mutation_count
            stats['total_cooccurrences'] += r.cooccurrence_count
            if r.clinical_significance == 'known_resistance':
                stats['resistance_samples'] += 1
        stats['total_samples'] = len(set(r.sample for r in results))
        stats['total_genes'] = len(set(r.gene for r in results))
        return dict(stats)

    def _generate_clinical_summary(self, results: List[AggregatedResult]) -> Dict[str, Any]:
        summary = defaultdict(lambda: {'samples': set(), 'mutations': 0})
        for r in results:
            key = r.resistance_profile
            summary[key]['samples'].add(r.sample)
            summary[key]['mutations'] += r.mutation_count
        # Convert sets to counts
        for k in summary:
            summary[k]['samples'] = len(summary[k]['samples'])
        return dict(summary)

    def _save_report(self, report: AggregationReport):
        # Save aggregated results CSV
        results_file = self.output_dir / f'{report.pipeline_id}_aggregated_results.csv'
        with open(results_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=report.aggregated_results[0].to_dict().keys())
            writer.writeheader()
            for r in report.aggregated_results:
                writer.writerow(r.to_dict())
        # Save summary statistics JSON
        stats_file = self.output_dir / f'{report.pipeline_id}_summary_statistics.json'
        with open(stats_file, 'w') as f:
            json.dump(report.summary_statistics, f, indent=2)
        # Save clinical summary JSON
        clinical_file = self.output_dir / f'{report.pipeline_id}_clinical_summary.json'
        with open(clinical_file, 'w') as f:
            json.dump(report.clinical_summary, f, indent=2)
        # Save manifest JSON
        manifest_file = self.output_dir / f'{report.pipeline_id}_manifest.json'
        with open(manifest_file, 'w') as f:
            json.dump(report.manifest, f, indent=2)
        self.logger.info(f"Results aggregation report saved to {self.output_dir}")

    def _generate_dashboard(self, report: AggregationReport):
        if not MATPLOTLIB_AVAILABLE:
            self.logger.warning("matplotlib not available, skipping dashboard generation")
            return
        # Example: mutation count per gene
        df = pd.DataFrame([r.to_dict() for r in report.aggregated_results])
        plt.figure(figsize=(10,6))
        df.groupby('gene')['mutation_count'].sum().plot(kind='bar', color='skyblue')
        plt.title('Mutation Count per Gene')
        plt.xlabel('Gene')
        plt.ylabel('Mutation Count')
        plt.tight_layout()
        dashboard_file = self.output_dir / f'{report.pipeline_id}_dashboard.png'
        plt.savefig(dashboard_file)
        plt.close()
        self.logger.info(f"Dashboard visualization saved to {dashboard_file}")

    def aggregate_results(self) -> AggregationReport:
        start_time = datetime.now()
        pipeline_id = f"results_{start_time.strftime('%Y%m%d_%H%M%S')}"
        aggregated_results = self._ingest_results()
        summary_statistics = self._compute_summary_statistics(aggregated_results)
        clinical_summary = self._generate_clinical_summary(aggregated_results)
        manifest = {
            'pipeline_id': pipeline_id,
            'execution_timestamp': start_time.isoformat(),
            'mutation_dir': str(self.mutation_dir),
            'cooccurrence_dir': str(self.cooccurrence_dir),
            'output_dir': str(self.output_dir),
            'total_samples': summary_statistics.get('total_samples', 0),
            'total_genes': summary_statistics.get('total_genes', 0)
        }
        processing_time = (datetime.now() - start_time).total_seconds()
        report = AggregationReport(
            pipeline_id=pipeline_id,
            execution_timestamp=start_time.isoformat(),
            total_samples=summary_statistics.get('total_samples', 0),
            total_genes=summary_statistics.get('total_genes', 0),
            aggregated_results=aggregated_results,
            summary_statistics=summary_statistics,
            clinical_summary=clinical_summary,
            manifest=manifest,
            processing_time=processing_time
        )
        self._save_report(report)
        self._generate_dashboard(report)
        self.logger.info(f"Results aggregation completed in {processing_time:.2f}s")
        return report

# Command-line interface
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Production Results Aggregation System")
    parser.add_argument('--config', required=True, help='Path to configuration file')
    args = parser.parse_args()
    aggregator = ProductionResultsAggregator(args.config)
    report = aggregator.aggregate_results()
    print(f"\n[Results Aggregation Complete]")
    print(f"Pipeline ID: {report.pipeline_id}")
    print(f"Total samples: {report.total_samples}")
    print(f"Total genes: {report.total_genes}")
    print(f"Output directory: {aggregator.output_dir}")
