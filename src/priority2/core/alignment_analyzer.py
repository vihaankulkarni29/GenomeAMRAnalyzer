"""
Alignment Analyzer
-----------------
Module for analyzing sequence alignment results and extracting AMR-relevant information.
Provides comprehensive reporting and statistical analysis of alignment data.
"""

import os
import json
import pandas as pd
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple
from dataclasses import dataclass, asdict
import logging

@dataclass
class AlignmentMetrics:
    """Metrics for a single alignment."""
    query_name: str
    query_length: int
    target_name: str
    target_length: int
    alignment_length: int
    identity: float
    coverage: float
    mapq: int
    strand: str

@dataclass
class SampleReport:
    """Report for a single sample."""
    sample_name: str
    total_sequences: int
    aligned_sequences: int
    alignment_rate: float
    top_hits: List[AlignmentMetrics]
    coverage_stats: Dict[str, float]

class AlignmentAnalyzer:
    """
    Analyzes alignment results and generates comprehensive reports.
    """
    
    def __init__(self, amr_database_info: Optional[Dict[str, Any]] = None, logger: Optional[logging.Logger] = None):
        self.amr_database_info = amr_database_info or {}
        self.logger = logger or self._setup_logger()
        
    def _setup_logger(self) -> logging.Logger:
        logger = logging.getLogger("AlignmentAnalyzer")
        if not logger.handlers:
            handler = logging.StreamHandler()
            formatter = logging.Formatter('%(asctime)s %(levelname)s %(name)s: %(message)s')
            handler.setFormatter(formatter)
            logger.addHandler(handler)
        logger.setLevel(logging.INFO)
        return logger

    def parse_paf_file(self, paf_file: str) -> List[AlignmentMetrics]:
        """
        Parse a PAF alignment file and extract metrics.
        """
        alignments = []
        try:
            with open(paf_file, 'r') as f:
                for line in f:
                    if line.strip():
                        metrics = self._parse_paf_line(line.strip())
                        if metrics:
                            alignments.append(metrics)
            self.logger.info(f"Parsed {len(alignments)} alignments from {paf_file}")
        except Exception as e:
            self.logger.error(f"Failed to parse PAF file {paf_file}: {e}")
        return alignments

    def _parse_paf_line(self, line: str) -> Optional[AlignmentMetrics]:
        """Parse a single PAF line."""
        try:
            fields = line.split('\t')
            if len(fields) >= 12:
                query_name = fields[0]
                query_length = int(fields[1])
                query_start = int(fields[2])
                query_end = int(fields[3])
                strand = fields[4]
                target_name = fields[5]
                target_length = int(fields[6])
                target_start = int(fields[7])
                target_end = int(fields[8])
                matches = int(fields[9])
                alignment_length = int(fields[10])
                mapq = int(fields[11])
                
                # Calculate metrics
                identity = (matches / alignment_length) * 100 if alignment_length > 0 else 0
                coverage = ((query_end - query_start) / query_length) * 100 if query_length > 0 else 0
                
                return AlignmentMetrics(
                    query_name=query_name,
                    query_length=query_length,
                    target_name=target_name,
                    target_length=target_length,
                    alignment_length=alignment_length,
                    identity=identity,
                    coverage=coverage,
                    mapq=mapq,
                    strand=strand
                )
        except (ValueError, IndexError) as e:
            self.logger.warning(f"Failed to parse PAF line: {e}")
        return None

    def analyze_sample(self, paf_file: str, min_identity: float = 80.0, min_coverage: float = 50.0) -> SampleReport:
        """
        Analyze alignments for a single sample and generate report.
        """
        sample_name = Path(paf_file).stem
        alignments = self.parse_paf_file(paf_file)
        
        if not alignments:
            return SampleReport(
                sample_name=sample_name,
                total_sequences=0,
                aligned_sequences=0,
                alignment_rate=0.0,
                top_hits=[],
                coverage_stats={}
            )
        
        # Filter high-quality alignments
        high_quality = [
            aln for aln in alignments 
            if aln.identity >= min_identity and aln.coverage >= min_coverage
        ]
        
        # Calculate statistics
        total_sequences = len(set(aln.query_name for aln in alignments))
        aligned_sequences = len(set(aln.query_name for aln in high_quality))
        alignment_rate = (aligned_sequences / total_sequences) * 100 if total_sequences > 0 else 0
        
        # Get top hits (highest identity)
        top_hits = sorted(high_quality, key=lambda x: x.identity, reverse=True)[:10]
        
        # Coverage statistics
        coverages = [aln.coverage for aln in high_quality]
        coverage_stats = {
            "mean": sum(coverages) / len(coverages) if coverages else 0,
            "min": min(coverages) if coverages else 0,
            "max": max(coverages) if coverages else 0,
            "median": sorted(coverages)[len(coverages)//2] if coverages else 0
        }
        
        return SampleReport(
            sample_name=sample_name,
            total_sequences=total_sequences,
            aligned_sequences=aligned_sequences,
            alignment_rate=alignment_rate,
            top_hits=top_hits,
            coverage_stats=coverage_stats
        )

    def analyze_batch(self, paf_files: List[str], output_dir: str) -> Dict[str, SampleReport]:
        """
        Analyze multiple PAF files and generate batch report.
        """
        os.makedirs(output_dir, exist_ok=True)
        reports = {}
        
        for paf_file in paf_files:
            try:
                report = self.analyze_sample(paf_file)
                reports[report.sample_name] = report
                self.logger.info(f"Analyzed {report.sample_name}: {report.aligned_sequences}/{report.total_sequences} sequences aligned")
            except Exception as e:
                self.logger.error(f"Failed to analyze {paf_file}: {e}")
        
        # Save batch summary
        self._save_batch_summary(reports, output_dir)
        return reports

    def _save_batch_summary(self, reports: Dict[str, SampleReport], output_dir: str):
        """Save batch analysis summary."""
        summary_data = []
        for sample_name, report in reports.items():
            summary_data.append({
                "Sample": sample_name,
                "Total_Sequences": report.total_sequences,
                "Aligned_Sequences": report.aligned_sequences,
                "Alignment_Rate_%": round(report.alignment_rate, 2),
                "Mean_Coverage_%": round(report.coverage_stats.get("mean", 0), 2),
                "Top_Hit_Target": report.top_hits[0].target_name if report.top_hits else "None",
                "Top_Hit_Identity_%": round(report.top_hits[0].identity, 2) if report.top_hits else 0
            })
        
        # Save as CSV
        df = pd.DataFrame(summary_data)
        csv_path = os.path.join(output_dir, "batch_summary.csv")
        df.to_csv(csv_path, index=False)
        
        # Save as JSON
        json_path = os.path.join(output_dir, "batch_summary.json")
        with open(json_path, 'w') as f:
            json.dump([asdict(report) for report in reports.values()], f, indent=2)
        
        self.logger.info(f"Batch summary saved to {csv_path} and {json_path}")

    def identify_amr_genes(self, alignments: List[AlignmentMetrics], min_identity: float = 95.0) -> List[str]:
        """
        Identify potential AMR genes based on high-confidence alignments.
        """
        amr_genes = []
        for aln in alignments:
            if aln.identity >= min_identity and aln.coverage >= 80.0:
                amr_genes.append(aln.target_name)
        return list(set(amr_genes))

    def generate_detailed_report(self, sample_report: SampleReport, output_path: str):
        """
        Generate a detailed HTML/text report for a single sample.
        """
        report_content = f"""
# AMR Analysis Report: {sample_report.sample_name}

## Summary
- Total Sequences: {sample_report.total_sequences}
- Aligned Sequences: {sample_report.aligned_sequences}
- Alignment Rate: {sample_report.alignment_rate:.2f}%

## Coverage Statistics
- Mean Coverage: {sample_report.coverage_stats.get('mean', 0):.2f}%
- Min Coverage: {sample_report.coverage_stats.get('min', 0):.2f}%
- Max Coverage: {sample_report.coverage_stats.get('max', 0):.2f}%

## Top Hits
"""
        for i, hit in enumerate(sample_report.top_hits[:5], 1):
            report_content += f"""
{i}. Target: {hit.target_name}
   - Identity: {hit.identity:.2f}%
   - Coverage: {hit.coverage:.2f}%
   - MAPQ: {hit.mapq}
"""
        
        with open(output_path, 'w') as f:
            f.write(report_content)
        
        self.logger.info(f"Detailed report saved to {output_path}")
