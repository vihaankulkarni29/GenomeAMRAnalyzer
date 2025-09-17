#!/usr/bin/env python3
"""
HTMLReportGenerator - Streamlined HTML report generation

This module generates interactive HTML reports from substitution analysis results
for AMR research with charts, tables, and data visualizations.

Usage:
    python html_report_generator.py --input results/substitutions.csv --output results/report.html

Author: MetaDataHarvester Pipeline
Version: 2.1 - Streamlined and Robust
"""

import os
import sys
import json
from pathlib import Path
from typing import Dict, List
import pandas as pd
from datetime import datetime
import argparse


class HTMLReportGenerator:
    """
    Streamlined HTML report generator for AMR analysis
    """

    def __init__(self, input_file: str, output_file: str):
        self.input_file = Path(input_file)
        self.output_file = Path(output_file)
        self.data = None
        self.stats = {}

        # Load and process data
        self.load_data()
        self.calculate_stats()

    def load_data(self):
        """Load substitution data from CSV"""
        try:
            self.data = pd.read_csv(self.input_file)
            print(f"Loaded {len(self.data)} substitution records")
        except Exception as e:
            print(f"Error loading data: {e}")
            sys.exit(1)

    def calculate_stats(self):
        """Calculate basic statistics"""
        if self.data is None or self.data.empty:
            self.stats = {
                'total_substitutions': 0,
                'unique_proteins': 0,
                'unique_organisms': 0,
                'protein_families': {},
                'organisms': {},
            }
            return

        self.stats = {
            'total_substitutions': len(self.data),
            'unique_proteins': self.data['Accession_Number'].nunique() if 'Accession_Number' in self.data.columns else 0,
            'unique_organisms': self.data['Organism'].nunique() if 'Organism' in self.data.columns else 0,
            'protein_families': self.data['Protein_Family'].value_counts().to_dict() if 'Protein_Family' in self.data.columns else {},
            'organisms': self.data['Organism'].value_counts().head(10).to_dict() if 'Organism' in self.data.columns else {},
        }

        # Position statistics
        if 'Residue_Position' in self.data.columns:
            positions = self.data['Residue_Position'].dropna()
            self.stats['position_stats'] = {
                'min': int(positions.min()) if not positions.empty else 0,
                'max': int(positions.max()) if not positions.empty else 0,
                'mean': round(positions.mean(), 1) if not positions.empty else 0,
            }

        # Common substitutions
        if 'Substitution' in self.data.columns:
            self.stats['common_substitutions'] = self.data['Substitution'].value_counts().head(10).to_dict()

    def generate_chart_data(self) -> Dict:
        """Generate data for charts"""
        chart_data = {
            'protein_families': self.stats.get('protein_families', {}),
            'organisms': self.stats.get('organisms', {}),
            'common_substitutions': self.stats.get('common_substitutions', {}),
        }

        # Position distribution
        if 'Residue_Position' in self.data.columns:
            positions = self.data['Residue_Position'].dropna().values
            if len(positions) > 0:
                hist_bins = list(range(0, int(max(positions)) + 50, 50))
                hist_counts = [0] * (len(hist_bins) - 1)

                for pos in positions:
                    for i in range(len(hist_bins) - 1):
                        if hist_bins[i] <= pos < hist_bins[i + 1]:
                            hist_counts[i] += 1
                            break

                chart_data['position_histogram'] = {
                    'bins': hist_bins[:-1],
                    'counts': hist_counts
                }

        return chart_data

    def generate_table_rows(self) -> str:
        """Generate HTML table rows"""
        rows = []
        for _, row in self.data.iterrows():
            rows.append(f"""
                <tr>
                    <td>{row.get('Accession_Number', 'N/A')}</td>
                    <td>{row.get('Organism', 'N/A')}</td>
                    <td><span class="badge bg-primary">{row.get('Protein_Family', 'N/A')}</span></td>
                    <td>{row.get('Substitution', 'N/A')}</td>
                    <td>{row.get('Residue_Position', 'N/A')}</td>
                    <td>{row.get('Alignment_File', 'N/A')}</td>
                </tr>
            """)
        return '\n'.join(rows)

    def generate_html(self) -> str:
        """Generate complete HTML report"""
        chart_data = self.generate_chart_data()
        table_rows = self.generate_table_rows()

        html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>AMR Substitution Analysis Report</title>
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css" rel="stylesheet">
    <link href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0/css/all.min.css" rel="stylesheet">
    <style>
        body {{ font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; background-color: #f8f9fa; }}
        .navbar {{ background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); border-bottom: 3px solid #5a67d8; }}
        .card {{ border: none; border-radius: 10px; box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1); transition: transform 0.2s; }}
        .card:hover {{ transform: translateY(-2px); }}
        .stat-card {{ background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; }}
        .chart-container {{ position: relative; height: 400px; width: 100%; }}
        .table-responsive {{ border-radius: 10px; overflow: hidden; }}
    </style>
</head>
<body>
    <!-- Navigation -->
    <nav class="navbar navbar-expand-lg navbar-dark">
        <div class="container">
            <a class="navbar-brand" href="#">
                <i class="fas fa-dna"></i> AMR Substitution Analysis Report
            </a>
            <div class="navbar-nav ms-auto">
                <span class="navbar-text">
                    Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
                </span>
            </div>
        </div>
    </nav>

    <div class="container mt-4">
        <!-- Summary Statistics -->
        <div class="row mb-4">
            <div class="col-md-3">
                <div class="card stat-card text-white">
                    <div class="card-body text-center">
                        <i class="fas fa-exchange-alt fa-2x mb-2"></i>
                        <h3>{self.stats['total_substitutions']:,}</h3>
                        <p class="mb-0">Total Substitutions</p>
                    </div>
                </div>
            </div>
            <div class="col-md-3">
                <div class="card stat-card text-white">
                    <div class="card-body text-center">
                        <i class="fas fa-virus fa-2x mb-2"></i>
                        <h3>{self.stats['unique_proteins']:,}</h3>
                        <p class="mb-0">Unique Proteins</p>
                    </div>
                </div>
            </div>
            <div class="col-md-3">
                <div class="card stat-card text-white">
                    <div class="card-body text-center">
                        <i class="fas fa-bacteria fa-2x mb-2"></i>
                        <h3>{self.stats['unique_organisms']:,}</h3>
                        <p class="mb-0">Organisms</p>
                    </div>
                </div>
            </div>
            <div class="col-md-3">
                <div class="card stat-card text-white">
                    <div class="card-body text-center">
                        <i class="fas fa-chart-line fa-2x mb-2"></i>
                        <h3>{len(self.stats.get('protein_families', {}))}</h3>
                        <p class="mb-0">Protein Families</p>
                    </div>
                </div>
            </div>
        </div>

        <!-- Charts Row -->
        <div class="row mb-4">
            <div class="col-md-6">
                <div class="card">
                    <div class="card-header">
                        <h5><i class="fas fa-chart-bar"></i> Substitution Positions</h5>
                    </div>
                    <div class="card-body">
                        <div class="chart-container">
                            <canvas id="positionChart"></canvas>
                        </div>
                    </div>
                </div>
            </div>
            <div class="col-md-6">
                <div class="card">
                    <div class="card-header">
                        <h5><i class="fas fa-dna"></i> Protein Families</h5>
                    </div>
                    <div class="card-body">
                        <div class="chart-container">
                            <canvas id="familyChart"></canvas>
                        </div>
                    </div>
                </div>
            </div>
        </div>

        <!-- Data Table -->
        <div class="card">
            <div class="card-header">
                <h5><i class="fas fa-table"></i> Substitution Data
                    <button class="btn btn-sm btn-primary float-end" onclick="exportData()">
                        <i class="fas fa-download"></i> Export CSV
                    </button>
                </h5>
            </div>
            <div class="card-body">
                <div class="table-responsive">
                    <table class="table table-striped table-hover" id="dataTable">
                        <thead class="table-dark">
                            <tr>
                                <th>Accession</th>
                                <th>Organism</th>
                                <th>Protein Family</th>
                                <th>Substitution</th>
                                <th>Position</th>
                                <th>Alignment File</th>
                            </tr>
                        </thead>
                        <tbody>
                            {table_rows}
                        </tbody>
                    </table>
                </div>
            </div>
        </div>
    </div>

    <!-- Embedded Data -->
    <script>
        const chartData = {json.dumps(chart_data)};
        const tableData = {self.data.to_json(orient='records')};
    </script>

    <!-- JavaScript -->
    <script>
        // Initialize charts
        document.addEventListener('DOMContentLoaded', function() {{
            // Position histogram
            if (chartData.position_histogram) {{
                const positionCtx = document.getElementById('positionChart').getContext('2d');
                new Chart(positionCtx, {{
                    type: 'bar',
                    data: {{
                        labels: chartData.position_histogram.bins.map(b => `${{b}}-${{b+49}}`),
                        datasets: [{{
                            label: 'Substitutions',
                            data: chartData.position_histogram.counts,
                            backgroundColor: 'rgba(102, 126, 234, 0.6)',
                            borderColor: 'rgba(102, 126, 234, 1)',
                            borderWidth: 1
                        }}]
                    }},
                    options: {{
                        responsive: true,
                        maintainAspectRatio: false,
                        scales: {{
                            y: {{
                                beginAtZero: true
                            }}
                        }}
                    }}
                }});
            }}

            // Protein families pie chart
            const familyCtx = document.getElementById('familyChart').getContext('2d');
            new Chart(familyCtx, {{
                type: 'doughnut',
                data: {{
                    labels: Object.keys(chartData.protein_families),
                    datasets: [{{
                        data: Object.values(chartData.protein_families),
                        backgroundColor: [
                            '#FF6384', '#36A2EB', '#FFCE56', '#4BC0C0',
                            '#9966FF', '#FF9F40', '#FF6384', '#C9CBCF'
                        ]
                    }}]
                }},
                options: {{
                    responsive: true,
                    maintainAspectRatio: false
                }}
            }});
        }});

        function exportData() {{
            const csvContent = [
                ['Accession', 'Organism', 'Protein_Family', 'Substitution', 'Residue_Position', 'Alignment_File'],
                ...JSON.parse(tableData).map(item => [
                    item.Accession_Number || '',
                    item.Organism || '',
                    item.Protein_Family || '',
                    item.Substitution || '',
                    item.Residue_Position || '',
                    item.Alignment_File || ''
                ])
            ].map(row => row.join(',')).join('\\n');

            const blob = new Blob([csvContent], {{ type: 'text/csv' }});
            const url = window.URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.href = url;
            a.download = 'substitutions.csv';
            a.click();
            window.URL.revokeObjectURL(url);
        }}
    </script>

</body>
</html>"""
        return html

    def generate_report(self):
        """Generate and save the HTML report"""
        try:
            html_content = self.generate_html()

            with open(self.output_file, 'w', encoding='utf-8') as f:
                f.write(html_content)

            print(f"HTML report generated successfully: {self.output_file}")
            print(f"Report includes {self.stats['total_substitutions']} substitutions")

        except Exception as e:
            print(f"Error generating report: {e}")
            raise

    def generate_enhanced_report(self, mutations, clinical_studies, mic_data, output_file):
        """Generate enhanced report with clinical and MIC data"""
        self.output_file = Path(output_file)

        # Convert mutations to DataFrame for compatibility
        mutation_data = []
        for mutation in mutations:
            mutation_data.append({
                'position': mutation.position,
                'reference_aa': mutation.reference_aa,
                'variant_aa': mutation.variant_aa,
                'mutation': mutation.mutation_key,
                'frequency': mutation.frequency,
                'confidence_score': mutation.confidence_score,
                'functional_impact': mutation.functional_impact,
                'resistance_risk': mutation.resistance_risk,
                'known_resistance': mutation.known_resistance_association,
                'clinical_evidence': len(mutation.clinical_evidence),
                'species_count': len(mutation.species_context)
            })

        self.data = pd.DataFrame(mutation_data)
        self.calculate_enhanced_stats(clinical_studies, mic_data)
        self.generate_enhanced_html(clinical_studies, mic_data)

    def calculate_enhanced_stats(self, clinical_studies, mic_data):
        """Calculate enhanced statistics including clinical and MIC data"""
        base_stats = self.calculate_stats()

        # Add clinical and MIC statistics
        self.stats.update({
            'clinical_studies': len(clinical_studies),
            'mic_data_points': len(mic_data),
            'significant_mutations': sum(1 for _, row in self.data.iterrows()
                                       if row.get('confidence_score', 0) > 0.8),
            'high_risk_mutations': sum(1 for _, row in self.data.iterrows()
                                     if str(row.get('resistance_risk', '')).lower() == 'high'),
            'known_resistance_mutations': sum(1 for _, row in self.data.iterrows()
                                            if pd.notna(row.get('known_resistance')))
        })

    def generate_enhanced_html(self, clinical_studies, mic_data):
        """Generate enhanced HTML with clinical and MIC data"""
        try:
            html_content = self.generate_enhanced_html_content(clinical_studies, mic_data)

            with open(self.output_file, 'w', encoding='utf-8') as f:
                f.write(html_content)

            print(f"Enhanced HTML report generated successfully: {self.output_file}")
            print(f"Report includes {len(self.data)} mutations, {len(clinical_studies)} studies, {len(mic_data)} MIC points")

        except Exception as e:
            print(f"Error generating enhanced report: {e}")
            raise

    def generate_enhanced_html_content(self, clinical_studies, mic_data):
        """Generate enhanced HTML content"""
        chart_data = self.generate_enhanced_chart_data(clinical_studies, mic_data)
        table_rows = self.generate_enhanced_table_rows()

        html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Enhanced AMR Scientific Analysis Report</title>
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css" rel="stylesheet">
    <link href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0/css/all.min.css" rel="stylesheet">
    <style>
        body {{ font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; background-color: #f8f9fa; }}
        .navbar {{ background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); border-bottom: 3px solid #5a67d8; }}
        .card {{ border: none; border-radius: 10px; box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1); transition: transform 0.2s; }}
        .card:hover {{ transform: translateY(-2px); }}
        .stat-card {{ background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; }}
        .high-risk {{ background: linear-gradient(135deg, #ff6b6b 0%, #ee5a24 100%); }}
        .clinical-card {{ background: linear-gradient(135deg, #4ecdc4 0%, #44a08d 100%); }}
        .mic-card {{ background: linear-gradient(135deg, #45b7d1 0%, #96c93d 100%); }}
        .chart-container {{ position: relative; height: 400px; width: 100%; }}
        .table-responsive {{ border-radius: 10px; overflow: hidden; }}
        .resistance-high {{ background-color: #ffe6e6; }}
        .resistance-medium {{ background-color: #fff3cd; }}
        .confidence-high {{ color: #28a745; font-weight: bold; }}
        .confidence-medium {{ color: #ffc107; font-weight: bold; }}
        .confidence-low {{ color: #dc3545; font-weight: bold; }}
    </style>
</head>
<body>
    <!-- Navigation -->
    <nav class="navbar navbar-expand-lg navbar-dark">
        <div class="container">
            <a class="navbar-brand" href="#">
                <i class="fas fa-microscope"></i> Enhanced AMR Scientific Analysis Report
            </a>
            <div class="navbar-nav ms-auto">
                <span class="navbar-text">
                    Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
                </span>
            </div>
        </div>
    </nav>

    <div class="container mt-4">
        <!-- Enhanced Summary Statistics -->
        <div class="row mb-4">
            <div class="col-md-3">
                <div class="card stat-card text-white">
                    <div class="card-body text-center">
                        <i class="fas fa-dna fa-2x mb-2"></i>
                        <h3>{self.stats.get('total_substitutions', 0):,}</h3>
                        <p class="mb-0">Mutations Called</p>
                    </div>
                </div>
            </div>
            <div class="col-md-3">
                <div class="card high-risk text-white">
                    <div class="card-body text-center">
                        <i class="fas fa-exclamation-triangle fa-2x mb-2"></i>
                        <h3>{self.stats.get('high_risk_mutations', 0)}</h3>
                        <p class="mb-0">High Risk</p>
                    </div>
                </div>
            </div>
            <div class="col-md-3">
                <div class="card clinical-card text-white">
                    <div class="card-body text-center">
                        <i class="fas fa-book-medical fa-2x mb-2"></i>
                        <h3>{self.stats.get('clinical_studies', 0)}</h3>
                        <p class="mb-0">Clinical Studies</p>
                    </div>
                </div>
            </div>
            <div class="col-md-3">
                <div class="card mic-card text-white">
                    <div class="card-body text-center">
                        <i class="fas fa-vial fa-2x mb-2"></i>
                        <h3>{self.stats.get('mic_data_points', 0)}</h3>
                        <p class="mb-0">MIC Data Points</p>
                    </div>
                </div>
            </div>
        </div>

        <!-- Charts Row -->
        <div class="row mb-4">
            <div class="col-md-6">
                <div class="card">
                    <div class="card-header">
                        <h5><i class="fas fa-chart-bar"></i> Mutation Confidence Distribution</h5>
                    </div>
                    <div class="card-body">
                        <div class="chart-container">
                            <canvas id="confidenceChart"></canvas>
                        </div>
                    </div>
                </div>
            </div>
            <div class="col-md-6">
                <div class="card">
                    <div class="card-header">
                        <h5><i class="fas fa-shield-alt"></i> Resistance Risk Assessment</h5>
                    </div>
                    <div class="card-body">
                        <div class="chart-container">
                            <canvas id="resistanceChart"></canvas>
                        </div>
                    </div>
                </div>
            </div>
        </div>

        <!-- Clinical Studies Section -->
        <div class="row mb-4">
            <div class="col-12">
                <div class="card">
                    <div class="card-header">
                        <h5><i class="fas fa-book-medical"></i> Relevant Clinical Studies
                            <span class="badge bg-primary ms-2">{len(clinical_studies)}</span>
                        </h5>
                    </div>
                    <div class="card-body">
                        <div class="table-responsive">
                            <table class="table table-striped table-hover">
                                <thead class="table-dark">
                                    <tr>
                                        <th>PMID</th>
                                        <th>Title</th>
                                        <th>Journal</th>
                                        <th>Study Type</th>
                                        <th>Relevance Score</th>
                                    </tr>
                                </thead>
                                <tbody>
                                    {self.generate_clinical_studies_rows(clinical_studies)}
                                </tbody>
                            </table>
                        </div>
                    </div>
                </div>
            </div>
        </div>

        <!-- Enhanced Data Table -->
        <div class="card">
            <div class="card-header">
                <h5><i class="fas fa-table"></i> Intelligent Mutation Analysis
                    <button class="btn btn-sm btn-primary float-end" onclick="exportData()">
                        <i class="fas fa-download"></i> Export CSV
                    </button>
                </h5>
            </div>
            <div class="card-body">
                <div class="table-responsive">
                    <table class="table table-striped table-hover" id="dataTable">
                        <thead class="table-dark">
                            <tr>
                                <th>Mutation</th>
                                <th>Position</th>
                                <th>Confidence</th>
                                <th>Resistance Risk</th>
                                <th>Functional Impact</th>
                                <th>Clinical Evidence</th>
                                <th>Species Count</th>
                                <th>Known Resistance</th>
                            </tr>
                        </thead>
                        <tbody>
                            {table_rows}
                        </tbody>
                    </table>
                </div>
            </div>
        </div>
    </div>

    <!-- Embedded Data -->
    <script>
        const chartData = {json.dumps(chart_data)};
        const tableData = {self.data.to_json(orient='records')};
    </script>

    <!-- JavaScript -->
    <script>
        // Initialize charts
        document.addEventListener('DOMContentLoaded', function() {{
            // Confidence distribution
            const confidenceCtx = document.getElementById('confidenceChart').getContext('2d');
            new Chart(confidenceCtx, {{
                type: 'bar',
                data: {{
                    labels: ['High (>0.8)', 'Medium (0.6-0.8)', 'Low (<0.6)'],
                    datasets: [{{
                        label: 'Mutations',
                        data: [
                            {sum(1 for _, row in self.data.iterrows() if row.get('confidence_score', 0) > 0.8)},
                            {sum(1 for _, row in self.data.iterrows() if 0.6 <= row.get('confidence_score', 0) <= 0.8)},
                            {sum(1 for _, row in self.data.iterrows() if row.get('confidence_score', 0) < 0.6)}
                        ],
                        backgroundColor: ['#28a745', '#ffc107', '#dc3545']
                    }}]
                }},
                options: {{
                    responsive: true,
                    maintainAspectRatio: false,
                    scales: {{
                        y: {{
                            beginAtZero: true
                        }}
                    }}
                }}
            }});

            // Resistance risk chart
            const resistanceCtx = document.getElementById('resistanceChart').getContext('2d');
            new Chart(resistanceCtx, {{
                type: 'doughnut',
                data: {{
                    labels: ['High', 'Medium', 'Low', 'Unknown'],
                    datasets: [{{
                        data: [
                            {sum(1 for _, row in self.data.iterrows() if str(row.get('resistance_risk', '')).lower() == 'high')},
                            {sum(1 for _, row in self.data.iterrows() if str(row.get('resistance_risk', '')).lower() == 'medium')},
                            {sum(1 for _, row in self.data.iterrows() if str(row.get('resistance_risk', '')).lower() == 'low')},
                            {sum(1 for _, row in self.data.iterrows() if str(row.get('resistance_risk', '')).lower() == 'unknown')}
                        ],
                        backgroundColor: ['#dc3545', '#ffc107', '#28a745', '#6c757d']
                    }}]
                }},
                options: {{
                    responsive: true,
                    maintainAspectRatio: false
                }}
            }});
        }});

        function exportData() {{
            const csvContent = [
                ['Mutation', 'Position', 'Confidence', 'Resistance_Risk', 'Functional_Impact', 'Clinical_Evidence', 'Species_Count', 'Known_Resistance'],
                ...JSON.parse(tableData).map(item => [
                    item.mutation || '',
                    item.position || '',
                    item.confidence_score || '',
                    item.resistance_risk || '',
                    item.functional_impact || '',
                    item.clinical_evidence || '',
                    item.species_count || '',
                    item.known_resistance || ''
                ])
            ].map(row => row.join(',')).join('\\n');

            const blob = new Blob([csvContent], {{ type: 'text/csv' }});
            const url = window.URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.href = url;
            a.download = 'enhanced_mutations.csv';
            a.click();
            window.URL.revokeObjectURL(url);
        }}
    </script>

</body>
</html>"""
        return html

    def generate_enhanced_chart_data(self, clinical_studies, mic_data):
        """Generate enhanced chart data"""
        return {
            'confidence_distribution': {
                'high': sum(1 for _, row in self.data.iterrows() if row.get('confidence_score', 0) > 0.8),
                'medium': sum(1 for _, row in self.data.iterrows() if 0.6 <= row.get('confidence_score', 0) <= 0.8),
                'low': sum(1 for _, row in self.data.iterrows() if row.get('confidence_score', 0) < 0.6)
            },
            'resistance_risks': {
                'high': sum(1 for _, row in self.data.iterrows() if str(row.get('resistance_risk', '')).lower() == 'high'),
                'medium': sum(1 for _, row in self.data.iterrows() if str(row.get('resistance_risk', '')).lower() == 'medium'),
                'low': sum(1 for _, row in self.data.iterrows() if str(row.get('resistance_risk', '')).lower() == 'low'),
                'unknown': sum(1 for _, row in self.data.iterrows() if str(row.get('resistance_risk', '')).lower() == 'unknown')
            }
        }

    def generate_enhanced_table_rows(self):
        """Generate enhanced table rows"""
        rows = []
        for _, row in self.data.iterrows():
            confidence_class = 'confidence-high' if row.get('confidence_score', 0) > 0.8 else \
                             'confidence-medium' if row.get('confidence_score', 0) > 0.6 else 'confidence-low'

            resistance_class = 'resistance-high' if str(row.get('resistance_risk', '')).lower() == 'high' else \
                             'resistance-medium' if str(row.get('resistance_risk', '')).lower() == 'medium' else ''

            rows.append(f"""
                <tr class="{resistance_class}">
                    <td><strong>{row.get('mutation', 'N/A')}</strong></td>
                    <td>{row.get('position', 'N/A')}</td>
                    <td><span class="{confidence_class}">{row.get('confidence_score', 0):.3f}</span></td>
                    <td><span class="badge bg-danger">{row.get('resistance_risk', 'Unknown')}</span></td>
                    <td>{row.get('functional_impact', 'N/A')}</td>
                    <td>{row.get('clinical_evidence', 0)}</td>
                    <td>{row.get('species_count', 0)}</td>
                    <td>{row.get('known_resistance', 'None')}</td>
                </tr>
            """)
        return '\n'.join(rows)

    def generate_clinical_studies_rows(self, clinical_studies):
        """Generate clinical studies table rows"""
        rows = []
        for study in clinical_studies[:20]:  # Show top 20
            rows.append(f"""
                <tr>
                    <td><a href="https://pubmed.ncbi.nlm.nih.gov/{study.pmid}/" target="_blank">{study.pmid}</a></td>
                    <td>{study.title[:80]}{'...' if len(study.title) > 80 else ''}</td>
                    <td>{study.journal}</td>
                    <td><span class="badge bg-info">{study.study_type}</span></td>
                    <td>{study.relevance_score:.1f}</td>
                </tr>
            """)
        return '\n'.join(rows)


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description="Generate interactive HTML reports from substitution analysis"
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Path to substitutions CSV file"
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Path for output HTML file"
    )

    args = parser.parse_args()

    # Generate report
    print("Generating HTML report...")
    generator = HTMLReportGenerator(args.input, args.output)
    generator.generate_report()

    print("\n============================================================")
    print("HTML REPORT GENERATED!")
    print("============================================================")
    print(f"Output file: {args.output}")
    print("Open the HTML file in your web browser to view the interactive report")


if __name__ == "__main__":
    main()