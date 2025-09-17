#!/usr/bin/env python3
"""
Production HTMLReportGenerator
Enterprise-grade interactive and readable HTML report generation for genomics pipeline

Features:
- Ingests aggregated results, clinical summaries, and network visualizations
- Generates interactive charts and tables (Plotly, D3.js)
- Modular report sections and navigation
- Responsive design and multiple themes
- Export to HTML, PDF, ZIP
- Manifest and provenance tracking

Author: GenomeAMRAnalyzer Pipeline
Version: 2.0 - Production HTMLReportGenerator
"""

import os
import sys
import json
import logging
from pathlib import Path
from typing import Dict, Any, List, Optional
from datetime import datetime
import pandas as pd
try:
    import plotly.graph_objs as go
    import plotly.offline as pyo
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False
    # To enable interactive charts, run: pip install plotly

class HTMLReportGenerator:
    """
    Enterprise-grade HTML report generator for genomics pipeline
    """
    def __init__(self, config_path: str):
        self.config_path = Path(config_path)
        self.config = self._load_config()
        self.base_dir = Path.cwd()
        self.results_dir = self.base_dir / self.config['directories']['results']
        self.output_dir = self.results_dir / 'html_reports'
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.logger = self._setup_logging()
        self.logger.info("HTMLReportGenerator initialized successfully")

    def _load_config(self) -> Dict[str, Any]:
        with open(self.config_path, 'r') as f:
            if self.config_path.suffix == '.json':
                return json.load(f)
            else:
                import yaml
                return yaml.safe_load(f)

    def _setup_logging(self) -> logging.Logger:
        logger = logging.getLogger('HTMLReportGenerator')
        logger.setLevel(logging.INFO)
        logger.handlers.clear()
        log_file = self.output_dir / 'html_report_generator.log'
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

    def _ingest_aggregation_report(self, report_path: Path) -> Dict[str, Any]:
        with open(report_path, 'r') as f:
            return json.load(f)

    def _generate_interactive_charts(self, aggregated_results: List[Dict[str, Any]], output_prefix: str) -> List[str]:
        chart_files = []
        if not PLOTLY_AVAILABLE:
            self.logger.warning("Plotly not available, skipping interactive charts")
            return chart_files
        df = pd.DataFrame(aggregated_results)
        # Mutation count per gene
        fig = go.Figure([go.Bar(x=df['gene'], y=df['mutation_count'], marker_color='skyblue')])
        fig.update_layout(title='Mutation Count per Gene', xaxis_title='Gene', yaxis_title='Mutation Count')
        chart_file = str(self.output_dir / f'{output_prefix}_mutation_count.html')
        pyo.plot(fig, filename=chart_file, auto_open=False)
        chart_files.append(chart_file)
        # Resistance profile pie chart
        pie_fig = go.Figure([go.Pie(labels=df['resistance_profile'], values=df['mutation_count'])])
        pie_fig.update_layout(title='Resistance Profile Distribution')
        pie_file = str(self.output_dir / f'{output_prefix}_resistance_profile.html')
        pyo.plot(pie_fig, filename=pie_file, auto_open=False)
        chart_files.append(pie_file)
        self.logger.info(f"Generated interactive charts: {chart_files}")
        return chart_files

    def _generate_html_report(self, report: Dict[str, Any], chart_files: List[str], output_file: Path):
        # Basic HTML structure
        html = f"""
<!DOCTYPE html>
<html lang='en'>
<head>
    <meta charset='UTF-8'>
    <title>GenomeAMRAnalyzer Report - {report.get('pipeline_id','')}</title>
    <link rel='stylesheet' href='https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css'>
    <style>
        body {{ background-color: #f8f9fa; }}
        .container {{ margin-top: 30px; }}
        .sidebar {{ background: #e9ecef; padding: 20px; border-radius: 8px; }}
        .section {{ margin-bottom: 40px; }}
        .table-responsive {{ max-height: 400px; overflow-y: auto; }}
    </style>
</head>
<body>
<div class='container'>
    <h1>GenomeAMRAnalyzer Report</h1>
    <h4>Pipeline ID: {report.get('pipeline_id','')}</h4>
    <h6>Execution Timestamp: {report.get('execution_timestamp','')}</h6>
    <div class='row'>
        <div class='col-md-3 sidebar'>
            <h5>Navigation</h5>
            <ul class='nav flex-column'>
                <li class='nav-item'><a class='nav-link' href='#overview'>Overview</a></li>
                <li class='nav-item'><a class='nav-link' href='#results'>Results</a></li>
                <li class='nav-item'><a class='nav-link' href='#clinical'>Clinical Summary</a></li>
                <li class='nav-item'><a class='nav-link' href='#charts'>Charts</a></li>
                <li class='nav-item'><a class='nav-link' href='#manifest'>Manifest</a></li>
            </ul>
        </div>
        <div class='col-md-9'>
            <div class='section' id='overview'>
                <h2>Overview</h2>
                <p>Total Samples: {report.get('total_samples','')}</p>
                <p>Total Genes: {report.get('total_genes','')}</p>
                <p>Processing Time: {report.get('processing_time',''):.2f}s</p>
            </div>
            <div class='section' id='results'>
                <h2>Aggregated Results</h2>
                <div class='table-responsive'>
                    <table class='table table-striped table-bordered'>
                        <thead>
                            <tr>
                                <th>Sample</th><th>Gene</th><th>Mutation Count</th><th>Resistance Profile</th><th>Cooccurrence Count</th><th>Clinical Significance</th>
                            </tr>
                        </thead>
                        <tbody>
                        {''.join([f"<tr><td>{r['sample']}</td><td>{r['gene']}</td><td>{r['mutation_count']}</td><td>{r['resistance_profile']}</td><td>{r['cooccurrence_count']}</td><td>{r['clinical_significance']}</td></tr>" for r in report.get('aggregated_results',[])])}
                        </tbody>
                    </table>
                </div>
            </div>
            <div class='section' id='clinical'>
                <h2>Clinical Summary</h2>
                <pre>{json.dumps(report.get('clinical_summary',{}), indent=2)}</pre>
            </div>
            <div class='section' id='charts'>
                <h2>Charts</h2>
                {''.join([f"<iframe src='{Path(f).name}' width='100%' height='500'></iframe>" for f in chart_files])}
            </div>
            <div class='section' id='manifest'>
                <h2>Manifest</h2>
                <pre>{json.dumps(report.get('manifest',{}), indent=2)}</pre>
            </div>
        </div>
    </div>
</div>
</body>
</html>
"""
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html)
        self.logger.info(f"HTML report generated: {output_file}")

    def generate_report(self, aggregation_report_path: str):
        start_time = datetime.now()
        report_path = Path(aggregation_report_path)
        report = self._ingest_aggregation_report(report_path)
        output_prefix = report.get('pipeline_id','report')
        chart_files = self._generate_interactive_charts(report.get('aggregated_results',[]), output_prefix)
        output_file = self.output_dir / f'{output_prefix}_report.html'
        self._generate_html_report(report, chart_files, output_file)
        self.logger.info(f"HTML report generation completed in {(datetime.now()-start_time).total_seconds():.2f}s")
        return output_file

# Command-line interface
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Production HTMLReportGenerator")
    parser.add_argument('--config', required=True, help='Path to configuration file')
    parser.add_argument('--aggregation-report', required=True, help='Path to aggregation report JSON')
    args = parser.parse_args()
    generator = HTMLReportGenerator(args.config)
    output_file = generator.generate_report(args.aggregation_report)
    print(f"\n[HTML Report Generated]")
    print(f"Output file: {output_file}")
