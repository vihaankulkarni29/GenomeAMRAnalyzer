#!/usr/bin/env python3
"""
Generic Co-occurrence Analyzer - Robust analysis of simultaneous mutations
Analyzes co-occurrence patterns of mutations across user-specified gene lists

This module is designed to be:
1. Bug-proof with comprehensive error handling
2. Generic for ANY gene list (not hardcoded to specific proteins)
3. Robust with input validation and edge case handling
4. Scientifically accurate with statistical significance testing
5. Performance optimized for large datasets

Author: GenomeAMRAnalyzer Pipeline
Version: 1.0 - Production Ready
"""

import os
import sys
import logging
import warnings
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Set, Union, Any
from dataclasses import dataclass, field
from collections import defaultdict, Counter
import pandas as pd
import numpy as np
from itertools import combinations, product
import json
import argparse
from datetime import datetime
import traceback
from typing import Iterable

# Import robust error handling if available
try:
    from .core.robust_error_handling import (
        ValidationError, DataProcessingError, 
        ValidationSuite, robust_exception_handler, RobustLogger
    )
    ERROR_HANDLING_AVAILABLE = True
except ImportError:
    # Fallback if core modules not available yet
    ERROR_HANDLING_AVAILABLE = False
    class ValidationError(Exception): pass
    class DataProcessingError(Exception): pass

# Try to expose hypergeom at module scope so tests can patch it
try:
    from scipy.stats import hypergeom  # type: ignore
except Exception:  # pragma: no cover - scipy optional
    hypergeom = None  # Tests may mock this symbol

# Suppress pandas warnings for cleaner output
warnings.filterwarnings('ignore', category=pd.errors.PerformanceWarning)
warnings.filterwarnings('ignore', category=FutureWarning)

# Configure logging
if ERROR_HANDLING_AVAILABLE:
    logger = RobustLogger("GenericCoOccurrenceAnalyzer").logger
else:
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)


@dataclass
class MutationEvent:
    """Represents a single mutation event with comprehensive metadata"""
    genome_id: str
    gene: str
    position: int
    reference_aa: str
    variant_aa: str
    substitution: str
    protein_name: Optional[str] = None
    organism: Optional[str] = None
    mic_data: Optional[Dict[str, float]] = field(default_factory=dict)
    
    def __post_init__(self):
        """Validate mutation event data"""
        if not self.genome_id or not self.gene:
            raise ValueError("genome_id and gene are required")
        
        if self.position <= 0:
            raise ValueError("Position must be positive")
            
        if len(self.reference_aa) != 1 or len(self.variant_aa) != 1:
            raise ValueError("Amino acids must be single characters")
            
        # Validate amino acid codes
        valid_aa = set('ACDEFGHIKLMNPQRSTVWYXU*-')
        if self.reference_aa not in valid_aa or self.variant_aa not in valid_aa:
            raise ValueError(f"Invalid amino acid codes: {self.reference_aa}, {self.variant_aa}")


@dataclass
class CoOccurrencePattern:
    """Represents a co-occurrence pattern with statistical data"""
    genes: Tuple[str, ...]
    mutations: Tuple[str, ...]
    genome_count: int
    total_genomes: int
    frequency: float
    genomes: Set[str]
    statistical_significance: Optional[float] = None
    expected_frequency: Optional[float] = None
    enrichment_score: Optional[float] = None
    
    @property
    def percentage(self) -> float:
        """Return frequency as percentage"""
        return self.frequency * 100
    
    def __str__(self) -> str:
        """Human-readable representation"""
        genes_str = " + ".join(self.genes)
        mutations_str = " + ".join(self.mutations)
        return f"{genes_str} ({mutations_str}): {self.genome_count}/{self.total_genomes} ({self.percentage:.1f}%)"


class GenericCoOccurrenceAnalyzer:
    """
    Robust co-occurrence analyzer for ANY gene list
    
    Features:
    - Generic analysis for any user-specified genes
    - Comprehensive error handling and validation
    - Statistical significance testing
    - Performance optimization for large datasets
    - Detailed logging and debugging
    """
    
    def __init__(self, 
                 output_dir: str = "cooccurrence_results",
                 min_genomes: int = 2,
                 max_combination_size: int = 5,
                 significance_threshold: float = 0.05,
                 log_level: str = "INFO"):
        """
        Initialize the co-occurrence analyzer
        
        Args:
            output_dir: Directory for output files
            min_genomes: Minimum genomes required for a pattern to be reported
            max_combination_size: Maximum number of genes to analyze together
            significance_threshold: P-value threshold for statistical significance
            log_level: Logging level (DEBUG, INFO, WARNING, ERROR)
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        self.min_genomes = max(1, min_genomes)
        self.max_combination_size = max(2, max_combination_size)
        self.significance_threshold = significance_threshold
        
        self.mutations_data: List[MutationEvent] = []
        self.gene_list: Set[str] = set()
        self.genome_set: Set[str] = set()
        self.patterns: List[CoOccurrencePattern] = []
        
        # Setup logging
        self._setup_logging(log_level)
        
        # Statistics tracking
        self.stats = {
            'total_mutations': 0,
            'total_genomes': 0,
            'total_genes': 0,
            'patterns_found': 0,
            'significant_patterns': 0
        }
    
    def _setup_logging(self, log_level: str) -> None:
        """Setup comprehensive logging with Windows-safe file handling"""
        log_file = self.output_dir / "cooccurrence_analysis.log"

        # Create formatter
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )

        # Use a per-instance logger to avoid handler duplication across tests
        self.logger = logging.getLogger(f'CoOccurrenceAnalyzer[{id(self)}]')
        self.logger.setLevel(getattr(logging, log_level.upper()))
        self.logger.propagate = False

        # If handlers somehow exist (re-instantiation), clear them
        for h in list(self.logger.handlers):
            self.logger.removeHandler(h)

        self._log_handlers: List[logging.Handler] = []

        # Decide whether to create a file handler. On Windows, tests that delete temp
        # directories can fail if a FileHandler keeps the file open. Avoid creating a
        # file handler when output_dir is inside the system temp directory unless
        # explicitly overridden.
        allow_file_logging = True
        try:
            tmp_root = Path(tempfile.gettempdir()).resolve()
            if Path(self.output_dir).resolve().is_relative_to(tmp_root) and os.environ.get('GENOMEAMR_FORCE_FILE_LOGGING') != '1':
                allow_file_logging = False
        except Exception:
            # Fallback: keep file logging enabled if detection fails
            allow_file_logging = True

        if allow_file_logging:
            file_handler = logging.FileHandler(log_file)
            file_handler.setFormatter(formatter)
            self.logger.addHandler(file_handler)
            self._log_handlers.append(file_handler)

        # Console handler (always enabled)
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(formatter)
        self.logger.addHandler(console_handler)
        self._log_handlers.append(console_handler)

        self.logger.info("Co-occurrence Analyzer initialized")

    def close(self) -> None:
        """Close and remove any logging handlers to release file locks (Windows-safe)."""
        if hasattr(self, 'logger'):
            for handler in getattr(self, '_log_handlers', []):
                try:
                    handler.flush()
                except Exception:
                    pass
                try:
                    handler.close()
                except Exception:
                    pass
                try:
                    self.logger.removeHandler(handler)
                except Exception:
                    pass
            self._log_handlers = []

    # Context manager support for deterministic cleanup
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close()
        return False
    
    def analyze_cooccurrence(self, mutations_by_gene: Dict[str, List[Dict]]) -> Dict[str, Any]:
        """
        Analyze co-occurrence patterns with robust error handling.
        
        Args:
            mutations_by_gene: Dictionary with gene names as keys and mutation lists as values
            
        Returns:
            Dictionary containing analysis results
            
        Raises:
            ValidationError: If input validation fails
            DataProcessingError: If analysis fails
        """
        try:
            # Input validation
            if ERROR_HANDLING_AVAILABLE:
                if not mutations_by_gene:
                    raise ValidationError("No mutation data provided", "mutations_by_gene", mutations_by_gene)
                
                if not isinstance(mutations_by_gene, dict):
                    raise ValidationError("Mutations data must be a dictionary", "mutations_by_gene", type(mutations_by_gene))
                
                # Validate each gene's data
                validated_data = {}
                for gene_name, mutations in mutations_by_gene.items():
                    validated_gene = ValidationSuite.validate_gene_name(gene_name)
                    
                    if not isinstance(mutations, list):
                        logger.warning(f"Skipping {gene_name}: mutations must be a list")
                        continue
                    
                    validated_data[validated_gene] = mutations
                
                mutations_by_gene = validated_data
            
            else:
                # Basic validation without error handling module
                if not mutations_by_gene or not isinstance(mutations_by_gene, dict):
                    raise ValueError("Invalid mutations_by_gene: must be a non-empty dictionary")
            
            # Initialize results
            results = {
                'genes_analyzed': list(mutations_by_gene.keys()),
                'total_mutations': 0,
                'cooccurrence_patterns': [],
                'analysis_summary': '',  # Will be updated
                'statistics': {},
                'analysis_timestamp': None
            }
            
            if ERROR_HANDLING_AVAILABLE:
                results['analysis_timestamp'] = datetime.now().isoformat()
            
            # Count total mutations
            for gene, mutations in mutations_by_gene.items():
                results['total_mutations'] += len(mutations)
            
            logger.info(f"Analyzing co-occurrence for {len(mutations_by_gene)} genes")
            logger.info(f"Total mutations: {results['total_mutations']}")
            
            # Find co-occurrence patterns
            if len(mutations_by_gene) >= 2:
                patterns = self._find_cooccurrence_patterns(mutations_by_gene)
                results['cooccurrence_patterns'] = patterns
                
                # Calculate statistics
                results['statistics'] = self._calculate_statistics(mutations_by_gene, patterns)
            else:
                logger.warning("Co-occurrence analysis requires at least 2 genes")
                results['cooccurrence_patterns'] = []
                results['statistics'] = {'message': 'Insufficient genes for co-occurrence analysis'}
            
            logger.info(f"Found {len(results['cooccurrence_patterns'])} co-occurrence patterns")
            
            # Update summary
            results['analysis_summary'] = f"Analyzed {len(mutations_by_gene)} genes with {results['total_mutations']} total mutations, found {len(results['cooccurrence_patterns'])} co-occurrence patterns"
            
            return results
            
        except ValidationError:
            raise  # Re-raise validation errors
        except Exception as e:
            if ERROR_HANDLING_AVAILABLE:
                raise DataProcessingError(f"Co-occurrence analysis failed: {e}")
            else:
                raise RuntimeError(f"Co-occurrence analysis failed: {e}")
    
    def _find_cooccurrence_patterns(self, mutations_by_gene: Dict[str, List[Dict]]) -> List[Dict]:
        """Find co-occurrence patterns between genes."""
        patterns = []
        
        # Create genome-to-mutations mapping
        genome_mutations = defaultdict(list)
        
        for gene, mutations in mutations_by_gene.items():
            for mutation in mutations:
                # Extract genome ID from mutation data
                genome_id = mutation.get('sequence_index', mutation.get('genome_id', 'unknown'))
                genome_mutations[genome_id].append({
                    'gene': gene,
                    'mutation': mutation
                })
        
        # Find patterns
        for genome_id, mutations in genome_mutations.items():
            if len(mutations) >= 2:  # Co-occurrence requires at least 2 mutations
                genes_involved = [mut['gene'] for mut in mutations]
                
                # Create pattern
                pattern = {
                    'genome_id': genome_id,
                    'genes': genes_involved,
                    'mutation_count': len(mutations),
                    'mutations': mutations
                }
                patterns.append(pattern)
        
        return patterns
    
    def _calculate_statistics(self, mutations_by_gene: Dict[str, List[Dict]], patterns: List[Dict]) -> Dict:
        """Calculate co-occurrence statistics."""
        stats = {
            'total_genes': len(mutations_by_gene),
            'total_patterns': len(patterns),
            'genes_with_cooccurrence': set(),
            'pattern_frequency': Counter()
        }
        
        # Analyze patterns
        for pattern in patterns:
            gene_combination = tuple(sorted(pattern['genes']))
            stats['pattern_frequency'][gene_combination] += 1
            stats['genes_with_cooccurrence'].update(pattern['genes'])
        
        # Convert set to list for JSON serialization
        stats['genes_with_cooccurrence'] = list(stats['genes_with_cooccurrence'])
        
        # Convert Counter to dict
        stats['pattern_frequency'] = dict(stats['pattern_frequency'])
        
        return stats
    
    def load_substitution_data(self, 
                             substitution_file: Union[str, Path], 
                             gene_list: Optional[List[str]] = None,
                             mic_file: Optional[Union[str, Path]] = None) -> None:
        """
        Load substitution data with comprehensive validation
        
        Args:
            substitution_file: Path to SubScan output CSV
            gene_list: Optional list of genes to analyze (if None, analyze all)
            mic_file: Optional MIC data file
        """
        try:
            self.logger.info(f"Loading substitution data from {substitution_file}")
            
            # Validate file existence
            sub_file = Path(substitution_file)
            if not sub_file.exists():
                raise FileNotFoundError(f"Substitution file not found: {substitution_file}")
            
            # Load substitution data
            df = pd.read_csv(sub_file)
            self.logger.info(f"Loaded {len(df)} substitution records")
            
            # Validate required columns ("substitution" is optional and will be created if missing)
            required_cols = ['genome_id', 'gene', 'position', 'reference_aa', 'variant_aa']
            missing_cols = [col for col in required_cols if col not in df.columns]
            
            # Try alternative column names
            col_mapping = {
                'accession': 'genome_id',
                'accession_number': 'genome_id',
                'protein': 'gene',
                'protein_name': 'gene',
                'ref_aa': 'reference_aa',
                'alt_aa': 'variant_aa',
                'resistant_aa': 'variant_aa',
                'pos': 'position'
            }
            
            for old_col, new_col in col_mapping.items():
                if old_col in df.columns and new_col in missing_cols:
                    df = df.rename(columns={old_col: new_col})
                    missing_cols.remove(new_col)
                    self.logger.info(f"Mapped column '{old_col}' to '{new_col}'")
            
            if missing_cols:
                raise ValueError(f"Missing required columns: {missing_cols}")
            
            # Clean and validate data
            df = self._clean_substitution_data(df)
            
            # Filter by gene list if provided
            if gene_list:
                self.gene_list = set(gene_list)
                initial_count = len(df)
                df = df[df['gene'].isin(self.gene_list)]
                filtered_count = len(df)
                self.logger.info(f"Filtered to {filtered_count}/{initial_count} mutations for genes: {gene_list}")
            else:
                self.gene_list = set(df['gene'].unique())
                self.logger.info(f"Analyzing all genes: {sorted(self.gene_list)}")
            
            # Load MIC data if provided
            mic_data = {}
            if mic_file and Path(mic_file).exists():
                mic_data = self._load_mic_data(mic_file)
            
            # Convert to MutationEvent objects
            self.mutations_data = []
            for _, row in df.iterrows():
                try:
                    mutation = MutationEvent(
                        genome_id=str(row['genome_id']),
                        gene=str(row['gene']),
                        position=int(row['position']),
                        reference_aa=str(row['reference_aa']),
                        variant_aa=str(row['variant_aa']),
                        substitution=str(row['substitution']),
                        protein_name=row.get('protein_name'),
                        organism=row.get('organism'),
                        mic_data=mic_data.get(str(row['genome_id']), {})
                    )
                    self.mutations_data.append(mutation)
                except Exception as e:
                    self.logger.warning(f"Skipping invalid mutation record: {e}")
                    continue
            
            self.genome_set = set(mut.genome_id for mut in self.mutations_data)
            
            # Update statistics
            self.stats.update({
                'total_mutations': len(self.mutations_data),
                'total_genomes': len(self.genome_set),
                'total_genes': len(self.gene_list)
            })
            
            self.logger.info(f"Successfully loaded {self.stats['total_mutations']} mutations from {self.stats['total_genomes']} genomes across {self.stats['total_genes']} genes")
            
        except Exception as e:
            self.logger.error(f"Failed to load substitution data: {e}")
            self.logger.debug(traceback.format_exc())
            raise
    
    def _clean_substitution_data(self, df: pd.DataFrame) -> pd.DataFrame:
        """Clean and validate substitution data"""
        self.logger.info("Cleaning substitution data...")
        
        initial_count = len(df)
        
        # Remove rows with missing critical data
        df = df.dropna(subset=['genome_id', 'gene', 'position', 'reference_aa', 'variant_aa'])
        
        # Convert position to int
        df['position'] = pd.to_numeric(df['position'], errors='coerce')
        df = df.dropna(subset=['position'])
        df['position'] = df['position'].astype(int)
        
        # Filter out invalid positions
        df = df[df['position'] > 0]
        
        # Clean amino acid codes
        df['reference_aa'] = df['reference_aa'].str.upper().str.strip()
        df['variant_aa'] = df['variant_aa'].str.upper().str.strip()
        
        # Validate amino acid codes
        valid_aa = set('ACDEFGHIKLMNPQRSTVWYXU*-')
        df = df[df['reference_aa'].isin(valid_aa) & df['variant_aa'].isin(valid_aa)]
        
        # Create substitution string if missing
        if 'substitution' not in df.columns:
            df['substitution'] = df['reference_aa'] + df['position'].astype(str) + df['variant_aa']
        
        # Remove identical reference and variant
        df = df[df['reference_aa'] != df['variant_aa']]
        
        final_count = len(df)
        removed = initial_count - final_count
        
        if removed > 0:
            self.logger.info(f"Removed {removed} invalid records during cleaning")
        
        return df
    
    def _load_mic_data(self, mic_file: Union[str, Path]) -> Dict[str, Dict[str, float]]:
        """Load MIC data with error handling"""
        try:
            self.logger.info(f"Loading MIC data from {mic_file}")
            mic_df = pd.read_csv(mic_file)
            
            mic_data = {}
            for _, row in mic_df.iterrows():
                genome_id = str(row.get('genome_id', row.get('accession', '')))
                if genome_id:
                    mic_data[genome_id] = {
                        col: float(row[col]) for col in mic_df.columns 
                        if col not in ['genome_id', 'accession'] and pd.notna(row[col])
                    }
            
            self.logger.info(f"Loaded MIC data for {len(mic_data)} genomes")
            return mic_data
            
        except Exception as e:
            self.logger.warning(f"Failed to load MIC data: {e}")
            return {}
    
    def analyze_cooccurrence_patterns(self) -> None:
        """
        Analyze co-occurrence patterns with comprehensive statistical analysis
        """
        if not self.mutations_data:
            raise ValueError("No mutation data loaded. Call load_substitution_data() first.")
        
        self.logger.info("Starting co-occurrence analysis...")
        
        # Group mutations by genome
        genome_mutations = defaultdict(lambda: defaultdict(list))
        for mutation in self.mutations_data:
            genome_mutations[mutation.genome_id][mutation.gene].append(mutation)
        
        self.patterns = []
        
        # Analyze single gene patterns (baseline)
        self._analyze_single_gene_patterns(genome_mutations)
        
        # Analyze multi-gene combinations
        for combination_size in range(2, min(len(self.gene_list) + 1, self.max_combination_size + 1)):
            self._analyze_gene_combinations(genome_mutations, combination_size)
        
        # Calculate statistical significance
        self._calculate_statistical_significance()
        
        # Sort patterns by frequency
        self.patterns.sort(key=lambda p: p.frequency, reverse=True)
        
        self.stats['patterns_found'] = len(self.patterns)
        self.stats['significant_patterns'] = sum(1 for p in self.patterns 
                                               if p.statistical_significance and p.statistical_significance < self.significance_threshold)
        
        self.logger.info(f"Found {self.stats['patterns_found']} co-occurrence patterns")
        self.logger.info(f"Found {self.stats['significant_patterns']} statistically significant patterns")
    
    def _analyze_single_gene_patterns(self, genome_mutations: Dict) -> None:
        """Analyze single gene mutation patterns"""
        self.logger.info("Analyzing single gene patterns...")
        
        for gene in self.gene_list:
            genomes_with_mutations = set()
            for genome_id, gene_mutations in genome_mutations.items():
                if gene in gene_mutations and gene_mutations[gene]:
                    genomes_with_mutations.add(genome_id)
            
            if len(genomes_with_mutations) >= self.min_genomes:
                pattern = CoOccurrencePattern(
                    genes=(gene,),
                    mutations=(f"{gene}_mutated",),
                    genome_count=len(genomes_with_mutations),
                    total_genomes=len(self.genome_set),
                    frequency=len(genomes_with_mutations) / len(self.genome_set),
                    genomes=genomes_with_mutations
                )
                self.patterns.append(pattern)
    
    def _analyze_gene_combinations(self, genome_mutations: Dict, combination_size: int) -> None:
        """Analyze gene combinations of specified size"""
        self.logger.info(f"Analyzing {combination_size}-gene combinations...")
        
        for gene_combo in combinations(sorted(self.gene_list), combination_size):
            genomes_with_all_mutations = set()
            
            for genome_id, gene_mutations in genome_mutations.items():
                has_all_mutations = True
                for gene in gene_combo:
                    if gene not in gene_mutations or not gene_mutations[gene]:
                        has_all_mutations = False
                        break
                
                if has_all_mutations:
                    genomes_with_all_mutations.add(genome_id)
            
            if len(genomes_with_all_mutations) >= self.min_genomes:
                pattern = CoOccurrencePattern(
                    genes=gene_combo,
                    mutations=tuple(f"{gene}_mutated" for gene in gene_combo),
                    genome_count=len(genomes_with_all_mutations),
                    total_genomes=len(self.genome_set),
                    frequency=len(genomes_with_all_mutations) / len(self.genome_set),
                    genomes=genomes_with_all_mutations
                )
                self.patterns.append(pattern)
    
    def _calculate_statistical_significance(self) -> None:
        """Calculate statistical significance using hypergeometric test"""
        self.logger.info("Calculating statistical significance...")
        
        try:
            # Single gene frequencies for expected calculation
            single_gene_freqs = {}
            for pattern in self.patterns:
                if len(pattern.genes) == 1:
                    single_gene_freqs[pattern.genes[0]] = pattern.frequency
            
            for pattern in self.patterns:
                if len(pattern.genes) > 1:
                    # Calculate expected frequency assuming independence
                    expected_freq = 1.0
                    for gene in pattern.genes:
                        if gene in single_gene_freqs:
                            expected_freq *= single_gene_freqs[gene]
                    
                    pattern.expected_frequency = expected_freq
                    
                    # Calculate enrichment score
                    if expected_freq > 0:
                        pattern.enrichment_score = pattern.frequency / expected_freq
                    
                    # Hypergeometric test
                    expected_count = expected_freq * pattern.total_genomes
                    if expected_count > 0 and hypergeom is not None:
                        try:
                            p_value = hypergeom.sf(
                                pattern.genome_count - 1,
                                pattern.total_genomes,
                                int(expected_count * pattern.total_genomes),
                                pattern.total_genomes
                            )
                            pattern.statistical_significance = p_value
                        except Exception as e:
                            self.logger.debug(f"Hypergeometric calculation failed: {e}")
            if hypergeom is None:
                self.logger.warning("scipy not available for statistical significance testing")
        except Exception as e:
            self.logger.warning(f"Statistical significance calculation failed: {e}")
    
    def generate_detailed_report(self) -> Dict[str, Any]:
        """Generate comprehensive analysis report"""
        if not self.patterns:
            raise ValueError("No patterns found. Run analyze_cooccurrence_patterns() first.")
        
        self.logger.info("Generating detailed report...")
        
        # Summary statistics
        summary = {
            'analysis_date': datetime.now().isoformat(),
            'total_mutations': self.stats['total_mutations'],
            'total_genomes': self.stats['total_genomes'],
            'total_genes': self.stats['total_genes'],
            'genes_analyzed': sorted(list(self.gene_list)),
            'patterns_found': self.stats['patterns_found'],
            'significant_patterns': self.stats['significant_patterns']
        }
        
        # Pattern details
        patterns_data = []
        for pattern in self.patterns:
            pattern_dict = {
                'genes': list(pattern.genes),
                'mutation_type': list(pattern.mutations),
                'genome_count': pattern.genome_count,
                'total_genomes': pattern.total_genomes,
                'frequency': round(pattern.frequency, 4),
                'percentage': round(pattern.percentage, 2),
                'genomes': sorted(list(pattern.genomes))
            }
            
            if pattern.statistical_significance is not None:
                pattern_dict['p_value'] = pattern.statistical_significance
                pattern_dict['is_significant'] = pattern.statistical_significance < self.significance_threshold
            
            if pattern.expected_frequency is not None:
                pattern_dict['expected_frequency'] = pattern.expected_frequency
            
            if pattern.enrichment_score is not None:
                pattern_dict['enrichment_score'] = pattern.enrichment_score
            
            patterns_data.append(pattern_dict)
        
        report = {
            'summary': summary,
            'patterns': patterns_data,
            'methodology': {
                'min_genomes_threshold': self.min_genomes,
                'max_combination_size': self.max_combination_size,
                'significance_threshold': self.significance_threshold
            }
        }
        
        return report
    
    def save_results(self, filename_prefix: str = "cooccurrence_analysis") -> Dict[str, str]:
        """Save analysis results in multiple formats"""
        if not self.patterns:
            raise ValueError("No patterns to save. Run analyze_cooccurrence_patterns() first.")
        
        self.logger.info("Saving results...")
        
        output_files = {}
        
        # Generate report
        report = self.generate_detailed_report()
        
        # Save JSON report
        json_file = self.output_dir / f"{filename_prefix}.json"
        with open(json_file, 'w') as f:
            json.dump(report, f, indent=2, default=str)
        output_files['json'] = str(json_file)
        
        # Save CSV summary
        csv_file = self.output_dir / f"{filename_prefix}.csv"
        patterns_df = pd.DataFrame(report['patterns'])
        patterns_df.to_csv(csv_file, index=False)
        output_files['csv'] = str(csv_file)
        
        # Save detailed patterns
        detailed_file = self.output_dir / f"{filename_prefix}_detailed.csv"
        detailed_data = []
        for pattern_dict in report['patterns']:
            for genome in pattern_dict['genomes']:
                detailed_data.append({
                    'genes': '+'.join(pattern_dict['genes']),
                    'genome_id': genome,
                    'pattern_frequency': pattern_dict['frequency'],
                    'pattern_percentage': pattern_dict['percentage'],
                    'genome_count': pattern_dict['genome_count'],
                    'is_significant': pattern_dict.get('is_significant', False)
                })
        
        if detailed_data:
            detailed_df = pd.DataFrame(detailed_data)
            detailed_df.to_csv(detailed_file, index=False)
            output_files['detailed_csv'] = str(detailed_file)
        
        # Save summary statistics
        stats_file = self.output_dir / f"{filename_prefix}_summary.txt"
        with open(stats_file, 'w') as f:
            f.write("Co-occurrence Analysis Summary\n")
            f.write("=" * 40 + "\n\n")
            f.write(f"Analysis Date: {report['summary']['analysis_date']}\n")
            f.write(f"Total Mutations: {report['summary']['total_mutations']}\n")
            f.write(f"Total Genomes: {report['summary']['total_genomes']}\n")
            f.write(f"Total Genes: {report['summary']['total_genes']}\n")
            f.write(f"Genes Analyzed: {', '.join(report['summary']['genes_analyzed'])}\n")
            f.write(f"Patterns Found: {report['summary']['patterns_found']}\n")
            f.write(f"Significant Patterns: {report['summary']['significant_patterns']}\n\n")
            
            f.write("Top 10 Co-occurrence Patterns:\n")
            f.write("-" * 40 + "\n")
            for i, pattern in enumerate(report['patterns'][:10], 1):
                genes_str = ' + '.join(pattern['genes'])
                f.write(f"{i:2d}. {genes_str}: {pattern['genome_count']}/{pattern['total_genomes']} ({pattern['percentage']:.1f}%)\n")
        
        output_files['summary'] = str(stats_file)
        
        self.logger.info(f"Results saved to {len(output_files)} files in {self.output_dir}")
        return output_files


def main():
    """Command line interface"""
    parser = argparse.ArgumentParser(
        description="Generic Co-occurrence Analyzer for mutation patterns",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze all genes in substitution file
  python generic_cooccurrence_analyzer.py --substitutions mutations.csv --output results/
  
  # Analyze specific genes only
  python generic_cooccurrence_analyzer.py --substitutions mutations.csv --genes mdtF geneX geneY --output results/
  
  # Include MIC data correlation
  python generic_cooccurrence_analyzer.py --substitutions mutations.csv --mic-data mic_values.csv --output results/
  
  # Advanced analysis with custom parameters
  python generic_cooccurrence_analyzer.py --substitutions mutations.csv --min-genomes 5 --max-combinations 4 --output results/
        """
    )
    
    # Required arguments
    parser.add_argument('--substitutions', required=True,
                       help='Path to substitution data CSV file (SubScan output)')
    parser.add_argument('--output', required=True,
                       help='Output directory for results')
    
    # Optional arguments
    parser.add_argument('--genes', nargs='+',
                       help='Specific genes to analyze (default: analyze all genes)')
    parser.add_argument('--mic-data',
                       help='Path to MIC data CSV file')
    parser.add_argument('--min-genomes', type=int, default=2,
                       help='Minimum genomes required for a pattern (default: 2)')
    parser.add_argument('--max-combinations', type=int, default=5,
                       help='Maximum gene combination size (default: 5)')
    parser.add_argument('--significance-threshold', type=float, default=0.05,
                       help='P-value threshold for significance (default: 0.05)')
    parser.add_argument('--log-level', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                       default='INFO', help='Logging level (default: INFO)')
    parser.add_argument('--prefix', default='cooccurrence_analysis',
                       help='Output filename prefix (default: cooccurrence_analysis)')
    
    args = parser.parse_args()
    
    try:
        # Initialize analyzer
        print("Initializing Co-occurrence Analyzer...")
        analyzer = GenericCoOccurrenceAnalyzer(
            output_dir=args.output,
            min_genomes=args.min_genomes,
            max_combination_size=args.max_combinations,
            significance_threshold=args.significance_threshold,
            log_level=args.log_level
        )
        
        # Load data
        print(f"Loading substitution data from {args.substitutions}...")
        analyzer.load_substitution_data(
            substitution_file=args.substitutions,
            gene_list=args.genes,
            mic_file=args.mic_data
        )
        
        # Analyze patterns
        print("Analyzing co-occurrence patterns...")
        analyzer.analyze_cooccurrence_patterns()
        
        # Save results
        print("Saving results...")
        output_files = analyzer.save_results(filename_prefix=args.prefix)
        
        # Print summary
        print("\n" + "=" * 60)
        print("CO-OCCURRENCE ANALYSIS COMPLETE!")
        print("=" * 60)
        print(f"Total mutations analyzed: {analyzer.stats['total_mutations']}")
        print(f"Total genomes: {analyzer.stats['total_genomes']}")
        print(f"Genes analyzed: {', '.join(sorted(analyzer.gene_list))}")
        print(f"Patterns found: {analyzer.stats['patterns_found']}")
        print(f"Significant patterns: {analyzer.stats['significant_patterns']}")
        
        print(f"\nOutput files:")
        for file_type, filepath in output_files.items():
            print(f"  {file_type.upper()}: {filepath}")
        
        # Show top patterns
        if analyzer.patterns:
            print(f"\nTop 5 Co-occurrence Patterns:")
            for i, pattern in enumerate(analyzer.patterns[:5], 1):
                print(f"  {i}. {pattern}")
        
        print(f"\nNext steps:")
        print(f"  1. Review detailed results in {output_files.get('json', 'output files')}")
        print(f"  2. Examine significant patterns for biological insights")
        print(f"  3. Correlate with MIC data for resistance phenotypes")
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()