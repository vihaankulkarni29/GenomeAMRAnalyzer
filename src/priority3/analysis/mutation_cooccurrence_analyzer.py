"""
Mutation Co-occurrence Analyzer
==============================

Enterprise-grade co-occurrence analysis for AMR genomics:
- Matrix-based mutation co-occurrence and association analysis
- Statistical significance testing (Fisher's exact, chi-square, correlation)
- Integration with MIC phenotype data
- Visualization-ready output (heatmaps, network graphs)
- Provenance and reproducibility tracking
- Performance optimization for large datasets

Engineering Principles:
- Statistical rigor and reproducibility
- Complete traceability to source data
- Modular, extensible design for new methods
- Memory-efficient for large-scale studies
"""

import os
import logging
import numpy as np
import pandas as pd
from typing import List, Dict, Any, Optional, Tuple, Set
from pathlib import Path
from dataclasses import dataclass, field
from collections import defaultdict, Counter
from scipy.stats import fisher_exact, chi2_contingency, spearmanr
import json

from ..db.repositories import GenomeRepository

@dataclass
class CooccurrenceResult:
    mutation_a: str
    mutation_b: str
    count_a: int
    count_b: int
    count_both: int
    p_value: float
    odds_ratio: Optional[float]
    method: str
    significant: bool
    notes: Optional[str] = None

class MutationCooccurrenceAnalyzer:

    def compute_multiway_cooccurrence(self, mutation_matrix: pd.DataFrame, mutation_groups: Dict[str, List[str]], max_order: int = 3) -> List[Dict[str, Any]]:
        """
        Compute co-occurrence statistics for all combinations of mutation groups (e.g., per protein).
        mutation_groups: dict mapping group name (e.g., protein) to list of mutation IDs
        max_order: maximum order of co-occurrence (e.g., 3 for triple-wise)
        Returns: List of dicts with group names, mutation sets, and co-occurrence stats
        """
        from itertools import combinations
        results = []
        group_names = list(mutation_groups.keys())
        for order in range(2, max_order+1):
            for group_combo in combinations(group_names, order):
                # Get all mutations in these groups
                mutation_set = set()
                for g in group_combo:
                    mutation_set.update(mutation_groups[g])
                if not mutation_set:
                    continue
                # Compute presence/absence for each genome
                present = mutation_matrix[list(mutation_set)].any(axis=1)
                count_present = int(present.sum())
                count_total = len(present)
                if count_present >= self.min_count:
                    results.append({
                        'group_combo': group_combo,
                        'mutation_ids': list(mutation_set),
                        'order': order,
                        'count_present': count_present,
                        'count_total': count_total,
                        'fraction_present': count_present / count_total if count_total else 0
                    })
        return results
    """
    Matrix-based co-occurrence and association analysis for AMR mutations.
    """
    def __init__(self, 
                 db_path: str = "priority3.db",
                 min_count: int = 5,
                 significance_level: float = 0.01,
                 output_dir: str = "cooccurrence_results"):
        self.repository = GenomeRepository(db_path)
        self.min_count = min_count
        self.significance_level = significance_level
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        self.logger = logging.getLogger("CooccurrenceAnalyzer")

    def load_mutation_matrix(self, 
                            genome_accessions: Optional[List[str]] = None
                            ) -> Tuple[pd.DataFrame, List[str], List[str]]:
        """
        Load binary mutation presence/absence matrix for all genomes.
        Rows: genomes, Columns: mutation IDs
        Returns: (matrix, genome_accessions, mutation_ids)
        """
        # Query all mutation artifacts from database
        mutation_data = []
        genome_ids = []
        mutation_set = set()
        
        # Get all genomes if not specified
        if genome_accessions is None:
            genome_records = self.repository.list_genomes()
            genome_accessions = [g.accession for g in genome_records]
        
        for accession in genome_accessions:
            artifacts = self.repository.list_artifacts(accession=accession, artifact_type="subscan_analysis")
            for artifact in artifacts:
                try:
                    with open(artifact.path, 'r') as f:
                        data = json.load(f)
                        muts = [m['mutation_id'] for m in data.get('mutations', [])]
                        mutation_data.append(set(muts))
                        genome_ids.append(accession)
                        mutation_set.update(muts)
                except Exception as e:
                    self.logger.warning(f"Failed to load artifact for {accession}: {e}")
        
        mutation_ids = sorted(mutation_set)
        # Build binary matrix
        matrix = np.zeros((len(genome_ids), len(mutation_ids)), dtype=int)
        for i, muts in enumerate(mutation_data):
            for j, mut in enumerate(mutation_ids):
                if mut in muts:
                    matrix[i, j] = 1
        df = pd.DataFrame(matrix, index=genome_ids, columns=mutation_ids)
        return df, genome_ids, mutation_ids

    def compute_cooccurrence(self, 
                            mutation_matrix: pd.DataFrame
                            ) -> List[CooccurrenceResult]:
        """
        Compute pairwise co-occurrence statistics for all mutation pairs.
        Returns: List of CooccurrenceResult
        """
        results = []
        mutation_ids = list(mutation_matrix.columns)
        n = len(mutation_ids)
        for i in range(n):
            for j in range(i+1, n):
                mut_a = mutation_ids[i]
                mut_b = mutation_ids[j]
                a = mutation_matrix[mut_a]
                b = mutation_matrix[mut_b]
                count_a = int(a.sum())
                count_b = int(b.sum())
                count_both = int((a & b).sum())
                # Only test if both mutations are present in enough genomes
                if count_a >= self.min_count and count_b >= self.min_count:
                    # 2x2 contingency table
                    table = np.array([
                        [count_both, count_a - count_both],
                        [count_b - count_both, len(a) - count_a - count_b + count_both]
                    ])
                    try:
                        test_result = fisher_exact(table)
                        if isinstance(test_result, tuple) and len(test_result) == 2:
                            odds_ratio, p_value = test_result
                            if isinstance(p_value, (float, int)):
                                p_value = float(p_value)
                            else:
                                p_value = float('nan')
                            if isinstance(odds_ratio, (float, int)):
                                odds_ratio = float(odds_ratio)
                            else:
                                odds_ratio = float('nan')
                        else:
                            odds_ratio, p_value = float('nan'), float('nan')
                        significant = (isinstance(p_value, float) and not isinstance(p_value, tuple) and p_value < self.significance_level)
                        results.append(CooccurrenceResult(
                            mutation_a=mut_a,
                            mutation_b=mut_b,
                            count_a=count_a,
                            count_b=count_b,
                            count_both=count_both,
                            p_value=p_value,
                            odds_ratio=odds_ratio,
                            method="fisher_exact",
                            significant=significant
                        ))
                    except Exception as e:
                        self.logger.warning(f"Fisher test failed for {mut_a}, {mut_b}: {e}")
        return results

    def correlate_with_phenotype(self, 
                                mutation_matrix: pd.DataFrame,
                                mic_data: Dict[str, Dict[str, Any]],
                                antibiotic: str
                                ) -> List[Dict[str, Any]]:
        """
        Correlate mutation presence/absence with MIC values for a given antibiotic.
        Returns: List of correlation results per mutation.
        """
        results = []
        for mut in mutation_matrix.columns:
            mic_values = []
            mut_present = []
            for genome in mutation_matrix.index:
                mic = None
                if genome in mic_data:
                    for rec in mic_data[genome].get('mic_records', []):
                        if rec['antibiotic_standardized'].lower() == antibiotic.lower():
                            mic = rec['mic_value']
                            break
                if mic is not None:
                    mic_values.append(mic)
                    mut_present.append(mutation_matrix.at[genome, mut])
            if len(mic_values) > 5:
                try:
                    corr, pval = spearmanr(mut_present, mic_values)
                    if not isinstance(pval, (float, int)):
                        pval = float('nan')
                        significant = False
                    else:
                        significant = pval < self.significance_level
                    results.append({
                        'mutation_id': mut,
                        'antibiotic': antibiotic,
                        'n': len(mic_values),
                        'correlation': corr,
                        'p_value': pval,
                        'significant': significant
                    })
                except Exception as e:
                    self.logger.warning(f"Correlation failed for {mut}: {e}")
        return results

    def save_results(self, results: List[CooccurrenceResult], filename: str):
        """
        Save co-occurrence results to CSV for downstream analysis/visualization.
        """
        df = pd.DataFrame([r.__dict__ for r in results])
        out_path = self.output_dir / filename
        df.to_csv(out_path, index=False)
        self.logger.info(f"Saved co-occurrence results to {out_path}")

    def close(self):
        self.repository.close()
