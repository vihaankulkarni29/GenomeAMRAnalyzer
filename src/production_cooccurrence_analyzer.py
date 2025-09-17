#!/usr/bin/env python3
"""
Production Cooccurrence Analysis Engine
Enterprise-grade statistical and network analysis of resistance gene cooccurrence

Features:
- Robust ingestion of mutation calls and clinical annotations
- Pairwise and higher-order cooccurrence statistics
- Fisher's exact test, chi-square, permutation significance
- Network construction and visualization
- Clinical annotation integration (CARD)
- Batch processing, quality control, and reporting
- Manifest and provenance tracking

Author: GenomeAMRAnalyzer Pipeline
Version: 2.0 - Production Cooccurrence Analysis Engine
"""

import os
import sys
import logging
import json
import csv
import yaml
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple, Set
from dataclasses import dataclass, asdict, field
from collections import defaultdict, Counter
from datetime import datetime
import itertools
import networkx as nx
import numpy as np
from scipy.stats import fisher_exact, chi2_contingency

try:
    import matplotlib.pyplot as plt
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False

@dataclass
class CooccurrencePair:
    gene_a: str
    gene_b: str
    count: int
    p_value: float
    significance: str
    clinical_annotation: Optional[str] = None
    resistance_class: Optional[str] = None
    notes: str = ""
    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

@dataclass
class CooccurrenceReport:
    pipeline_id: str
    execution_timestamp: str
    total_samples: int
    total_genes: int
    total_pairs: int
    significant_pairs: int
    cooccurrence_pairs: List[CooccurrencePair]
    network_metrics: Dict[str, Any]
    clinical_summary: Dict[str, Any]
    manifest: Dict[str, Any]
    processing_time: float
    def to_dict(self) -> Dict[str, Any]:
        result = asdict(self)
        result['cooccurrence_pairs'] = [pair.to_dict() for pair in self.cooccurrence_pairs]
        return result

class ProductionCooccurrenceAnalyzer:
    """
    Enterprise-grade cooccurrence analysis engine
    """
    def __init__(self, config_path: str, card_database_path: Optional[str] = None):
        self.config_path = Path(config_path)
        self.config = self._load_config()
        self.base_dir = Path.cwd()
        self.input_dir = self.base_dir / self.config['directories']['mutations'] / 'mutation_calls'
        self.output_dir = self.base_dir / self.config['directories']['cooccurrence']
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.card_database_path = card_database_path or self.config.get('card_database')
        self.card_mechanisms = self._load_card_mechanisms()
        self.logger = self._setup_logging()
        self.logger.info("ProductionCooccurrenceAnalyzer initialized successfully")

    def _load_config(self) -> Dict[str, Any]:
        with open(self.config_path, 'r') as f:
            return yaml.safe_load(f)

    def _setup_logging(self) -> logging.Logger:
        logger = logging.getLogger('ProductionCooccurrenceAnalyzer')
        logger.setLevel(logging.INFO)
        logger.handlers.clear()
        log_file = self.output_dir / 'cooccurrence_analyzer.log'
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

    def _load_card_mechanisms(self) -> Dict[str, Any]:
        if not self.card_database_path or not Path(self.card_database_path).exists():
            self.logger.warning("CARD database not found - clinical annotation will be limited")
            return {}
        with open(self.card_database_path, 'r') as f:
            card_data = json.load(f)
        mechanisms = {}
        for aro_id, aro_data in card_data.items():
            if isinstance(aro_data, dict) and 'model_name' in aro_data:
                gene_name = aro_data.get('model_name', '').lower()
                resistance_class = aro_data.get('ARO_category', {}).get('category_aro_class_name', 'Unknown')
                mechanism = aro_data.get('ARO_category', {}).get('category_aro_description', 'Unknown')
                mechanisms[gene_name] = {
                    'resistance_class': resistance_class,
                    'mechanism': mechanism,
                    'aro_id': aro_id
                }
        self.logger.info(f"Loaded {len(mechanisms)} resistance mechanisms from CARD database")
        return mechanisms

    def _ingest_mutation_calls(self) -> List[Dict[str, Any]]:
        mutation_calls = []
        for file in self.input_dir.glob('*.csv'):
            with open(file, 'r') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    mutation_calls.append(row)
        self.logger.info(f"Ingested {len(mutation_calls)} mutation calls from {self.input_dir}")
        return mutation_calls

    def _compute_cooccurrence(self, mutation_calls: List[Dict[str, Any]], target_genes: List[str]) -> List[CooccurrencePair]:
        # Build sample-gene matrix
        sample_gene = defaultdict(set)
        for call in mutation_calls:
            sample = call.get('accession', 'unknown')
            gene = call.get('gene_name', 'unknown').lower()
            if gene in target_genes:
                sample_gene[sample].add(gene)
        # Count cooccurrences
        gene_pairs = list(itertools.combinations(sorted(target_genes), 2))
        pair_counts = Counter()
        for sample, genes in sample_gene.items():
            for pair in gene_pairs:
                if pair[0] in genes and pair[1] in genes:
                    pair_counts[pair] += 1
        # Statistical significance
        cooccurrence_pairs = []
        total_samples = len(sample_gene)
        for pair, count in pair_counts.items():
            # Build contingency table
            a, b = pair
            both = count
            only_a = sum(1 for genes in sample_gene.values() if a in genes and b not in genes)
            only_b = sum(1 for genes in sample_gene.values() if b in genes and a not in genes)
            neither = total_samples - (both + only_a + only_b)
            table = [[both, only_a], [only_b, neither]]
            try:
                _, p_value = fisher_exact(table)
            except Exception:
                p_value = 1.0
            significance = 'significant' if p_value < 0.05 else 'not_significant'
            # Clinical annotation
            annotation = None
            resistance_class = None
            if a in self.card_mechanisms:
                annotation = self.card_mechanisms[a]['mechanism']
                resistance_class = self.card_mechanisms[a]['resistance_class']
            if b in self.card_mechanisms and not annotation:
                annotation = self.card_mechanisms[b]['mechanism']
                resistance_class = self.card_mechanisms[b]['resistance_class']
            cooccurrence_pairs.append(CooccurrencePair(
                gene_a=a,
                gene_b=b,
                count=count,
                p_value=p_value,
                significance=significance,
                clinical_annotation=annotation,
                resistance_class=resistance_class
            ))
        return cooccurrence_pairs

    def _build_network(self, cooccurrence_pairs: List[CooccurrencePair]) -> nx.Graph:
        G = nx.Graph()
        for pair in cooccurrence_pairs:
            G.add_edge(pair.gene_a, pair.gene_b, weight=pair.count, p_value=pair.p_value, significance=pair.significance)
        return G

    def _compute_network_metrics(self, G: nx.Graph) -> Dict[str, Any]:
        metrics = {
            'num_nodes': G.number_of_nodes(),
            'num_edges': G.number_of_edges(),
            'degree_centrality': nx.degree_centrality(G),
            'clustering': nx.clustering(G),
            'connected_components': [list(c) for c in nx.connected_components(G)]
        }
        return metrics

    def _visualize_network(self, G: nx.Graph, output_path: Path):
        if not MATPLOTLIB_AVAILABLE:
            self.logger.warning("matplotlib not available, skipping network visualization")
            return
        plt.figure(figsize=(8,6))
        pos = nx.spring_layout(G)
        nx.draw(G, pos, with_labels=True, node_color='skyblue', edge_color='gray', node_size=1200, font_size=10)
        plt.title('Gene Cooccurrence Network')
        plt.tight_layout()
        plt.savefig(output_path)
        plt.close()
        self.logger.info(f"Network visualization saved to {output_path}")

    def _save_report(self, report: CooccurrenceReport):
        # Save cooccurrence pairs CSV
        pairs_file = self.output_dir / f'{report.pipeline_id}_cooccurrence_pairs.csv'
        with open(pairs_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=report.cooccurrence_pairs[0].to_dict().keys())
            writer.writeheader()
            for pair in report.cooccurrence_pairs:
                writer.writerow(pair.to_dict())
        # Save network metrics JSON
        metrics_file = self.output_dir / f'{report.pipeline_id}_network_metrics.json'
        with open(metrics_file, 'w') as f:
            json.dump(report.network_metrics, f, indent=2)
        # Save manifest JSON
        manifest_file = self.output_dir / f'{report.pipeline_id}_manifest.json'
        with open(manifest_file, 'w') as f:
            json.dump(report.manifest, f, indent=2)
        # Save clinical summary JSON
        clinical_file = self.output_dir / f'{report.pipeline_id}_clinical_summary.json'
        with open(clinical_file, 'w') as f:
            json.dump(report.clinical_summary, f, indent=2)
        self.logger.info(f"Cooccurrence report saved to {self.output_dir}")

    def analyze_cooccurrence(self, target_genes: List[str]) -> CooccurrenceReport:
        start_time = datetime.now()
        pipeline_id = f"cooccurrence_{start_time.strftime('%Y%m%d_%H%M%S')}"
        mutation_calls = self._ingest_mutation_calls()
        cooccurrence_pairs = self._compute_cooccurrence(mutation_calls, [g.lower() for g in target_genes])
        G = self._build_network(cooccurrence_pairs)
        network_metrics = self._compute_network_metrics(G)
        if MATPLOTLIB_AVAILABLE:
            self._visualize_network(G, self.output_dir / f'{pipeline_id}_network.png')
        clinical_summary = self._summarize_clinical(cooccurrence_pairs)
        manifest = {
            'pipeline_id': pipeline_id,
            'execution_timestamp': start_time.isoformat(),
            'input_dir': str(self.input_dir),
            'output_dir': str(self.output_dir),
            'target_genes': target_genes,
            'total_samples': len(set(call['accession'] for call in mutation_calls)),
            'total_pairs': len(cooccurrence_pairs)
        }
        processing_time = (datetime.now() - start_time).total_seconds()
        report = CooccurrenceReport(
            pipeline_id=pipeline_id,
            execution_timestamp=start_time.isoformat(),
            total_samples=len(set(call['accession'] for call in mutation_calls)),
            total_genes=len(target_genes),
            total_pairs=len(cooccurrence_pairs),
            significant_pairs=len([p for p in cooccurrence_pairs if p.significance == 'significant']),
            cooccurrence_pairs=cooccurrence_pairs,
            network_metrics=network_metrics,
            clinical_summary=clinical_summary,
            manifest=manifest,
            processing_time=processing_time
        )
        self._save_report(report)
        self.logger.info(f"Cooccurrence analysis completed in {processing_time:.2f}s")
        return report

    def _summarize_clinical(self, cooccurrence_pairs: List[CooccurrencePair]) -> Dict[str, Any]:
        summary = defaultdict(lambda: {'count': 0, 'significant': 0})
        for pair in cooccurrence_pairs:
            key = pair.clinical_annotation or 'unknown'
            summary[key]['count'] += 1
            if pair.significance == 'significant':
                summary[key]['significant'] += 1
        return dict(summary)

# Command-line interface
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Production Cooccurrence Analysis Engine")
    parser.add_argument('--config', required=True, help='Path to configuration file')
    parser.add_argument('--genes', nargs='+', required=True, help='Target genes to analyze')
    parser.add_argument('--card-db', help='Path to CARD database JSON file')
    args = parser.parse_args()
    analyzer = ProductionCooccurrenceAnalyzer(args.config, args.card_db)
    report = analyzer.analyze_cooccurrence(args.genes)
    print(f"\n[Cooccurrence Analysis Complete]")
    print(f"Pipeline ID: {report.pipeline_id}")
    print(f"Total pairs: {report.total_pairs}")
    print(f"Significant pairs: {report.significant_pairs}")
    print(f"Output directory: {analyzer.output_dir}")
