#!/usr/bin/env python3
"""
Example script demonstrating the new interactive Plotly charts in Enhanced HTML Reporter
"""

import sys
import os
from pathlib import Path

# Add src directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from enhanced_html_reporter import EnhancedHTMLReportGenerator

def create_sample_data():
    """Create sample data for demonstration"""
    
    # Sample genome data
    genomes = [
        {
            'accession': 'GCF_000005825.2',
            'organism': 'Escherichia coli K-12',
            'strain': 'MG1655',
            'biosample': 'SAMN02604091',
            'bioproject': 'PRJNA225',
            'download_date': '2025-09-18'
        },
        {
            'accession': 'GCF_000008865.2',
            'organism': 'Escherichia coli O157:H7',
            'strain': 'EDL933',
            'biosample': 'SAMN02603478',
            'bioproject': 'PRJNA224',
            'download_date': '2025-09-18'
        }
    ]
    
    # Sample mutation data with frequency patterns
    mutations = [
        {'genome_id': 'GCF_000005825.2', 'gene': 'gyrA', 'position': 83, 'reference_aa': 'S', 'variant_aa': 'L', 'substitution': 'S83L'},
        {'genome_id': 'GCF_000005825.2', 'gene': 'parC', 'position': 80, 'reference_aa': 'S', 'variant_aa': 'I', 'substitution': 'S80I'},
        {'genome_id': 'GCF_000008865.2', 'gene': 'gyrA', 'position': 83, 'reference_aa': 'S', 'variant_aa': 'L', 'substitution': 'S83L'},
        {'genome_id': 'GCF_000008865.2', 'gene': 'gyrA', 'position': 87, 'reference_aa': 'D', 'variant_aa': 'N', 'substitution': 'D87N'},
        {'genome_id': 'GCF_000005825.2', 'gene': 'acrR', 'position': 45, 'reference_aa': 'A', 'variant_aa': 'V', 'substitution': 'A45V'},
        {'genome_id': 'GCF_000008865.2', 'gene': 'acrR', 'position': 45, 'reference_aa': 'A', 'variant_aa': 'V', 'substitution': 'A45V'},
    ]
    
    # Sample co-occurrence data
    cooccurrence = [
        {
            'mutation_a': 'gyrA_S83L',
            'mutation_b': 'parC_S80I',
            'count_a': 15,
            'count_b': 12,
            'count_both': 8,
            'p_value': 0.001,
            'odds_ratio': 3.5,
            'significant': True
        }
    ]
    
    # Sample MIC data
    mic_data = [
        {
            'accession': 'GCF_000005825.2',
            'antibiotic': 'Ciprofloxacin',
            'mic_value': 2.0,
            'mic_unit': 'μg/mL',
            'resistance_profile': 'Intermediate',
            'quality_score': 0.95
        }
    ]
    
    # Sample statistics
    stats = {
        'total_genomes': 2,
        'total_mutations': 6,
        'unique_mutations': 4,
        'analysis_time': '2025-09-18 19:45:00'
    }
    
    return genomes, mutations, cooccurrence, mic_data, stats

def main():
    """Generate example interactive HTML report"""
    
    print("🧬 Enhanced HTML Reporter - Interactive Charts Demo")
    print("=" * 55)
    
    # Create sample data
    genomes, mutations, cooccurrence, mic_data, stats = create_sample_data()
    
    # Initialize the enhanced reporter
    output_dir = Path("example_reports")
    reporter = EnhancedHTMLReportGenerator(str(output_dir))
    
    # Generate interactive report with Plotly charts
    report_path = reporter.generate_template_based_report(
        run_id="demo_2025_09_18",
        genomes=genomes,
        mutations=mutations,
        mic_data=mic_data,
        cooccurrence=cooccurrence,
        stats=stats,
        artifact_links={}
    )
    
    print(f"✅ Interactive HTML report generated: {report_path}")
    print()
    print("📊 Features included in this report:")
    print("  • Interactive Plotly.js mutation frequency bar chart")
    print("  • Responsive DataTables for all data")
    print("  • Modern, professional styling")
    print("  • Hover tooltips and zoom functionality")
    print("  • Mobile-friendly responsive design")
    print()
    print(f"🌐 Open the report in your browser: file://{Path(report_path).absolute()}")

if __name__ == "__main__":
    main()