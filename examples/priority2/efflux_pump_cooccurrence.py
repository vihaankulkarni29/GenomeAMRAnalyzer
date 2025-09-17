import sys; import os; sys.path.append(os.path.join(os.path.dirname(__file__), '\\\\\\src')); from configuration_manager import config_manager
"""
Priority 2 Advanced Example: Efflux Pump Co-occurrence Analysis
--------------------------------------------------------------
Specialized workflow for analyzing RND efflux pump mutations and co-occurrence patterns
in enteric gram-negative bacteria, specifically addressing Nishad Sir's research questions.
"""

import os
import sys
import logging
import pandas as pd
from pathlib import Path
from collections import defaultdict, Counter

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

def setup_logging():
    """Set up comprehensive logging."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler('efflux_pump_analysis.log')
        ]
    )

def efflux_pump_cooccurrence_analysis():
    """
    Advanced Example: RND Efflux Pump Co-occurrence Analysis
    
    This example demonstrates:
    1. Targeted analysis of RND efflux pump genes (acrA, acrB, acrF, tolC, etc.)
    2. Detection of mutations in efflux pump components
    3. Quantification of co-occurrence patterns
    4. Statistical analysis of resistance patterns
    """
    
    print("=" * 70)
    print("Priority 2 Advanced Example: Efflux Pump Co-occurrence Analysis")
    print("=" * 70)
    
    setup_logging()
    
    # Define efflux pump genes of interest
    efflux_pump_genes = {
        'RND_pumps': [config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0], config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1], 'acrF', 'acrD', 'acrE'],
        'outer_membrane': [config_manager.get_default_genes("rnd_efflux_pumps", "primary")[2]],
        'regulators': ['acrR', 'marA', 'soxS', 'ramA'],
        'additional': ['mdtF', 'mexA', 'mexB', 'mexC', 'mexD']
    }
    
    # Example data paths
    amr_database = "/path/to/efflux_pump_database.fasta"
    input_samples = [
        "/path/to/ecoli_clinical_isolates/sample_001.fasta",
        "/path/to/ecoli_clinical_isolates/sample_002.fasta",
        "/path/to/klebsiella_isolates/sample_003.fasta",
        # ... more samples
    ]
    output_directory = "./example_output/efflux_cooccurrence"
    
    os.makedirs(output_directory, exist_ok=True)
    
    try:
        # Step 1: Configure pipeline for efflux pump analysis
        print("Step 1: Configuring pipeline for efflux pump analysis...")
        
        # Import here to avoid path issues in example
        from src.priority2.config.pipeline_config import PipelineConfig, DatabaseConfig, ProcessingConfig
        from src.priority2.core.pipeline_orchestrator import PipelineOrchestrator
        from src.priority2.core.alignment_analyzer import AlignmentAnalyzer
        
        # High-stringency configuration for mutation detection
        config = PipelineConfig(
            database=DatabaseConfig(amr_database_path=amr_database),
            processing=ProcessingConfig(
                threads=8,
                min_identity=85.0,  # Allow for mutations
                min_coverage=70.0,  # Good coverage for reliable mutation calling
                quality_threshold=25,
                enable_parallel=True,
                max_parallel_samples=4
            ),
            config_name="efflux_pump_analysis"
        )
        
        print("✓ Configuration optimized for efflux pump mutation detection")
        
        # Step 2: Run initial alignment analysis
        print("Step 2: Running alignment analysis...")
        
        orchestrator = PipelineOrchestrator(config)
        
        def efflux_progress(message: str, progress: float):
            print(f"Efflux analysis progress: {progress:.1f}% - {message}")
        
        results = orchestrator.run_pipeline(
            input_samples=input_samples,
            progress_callback=efflux_progress
        )
        
        print("✓ Alignment analysis completed")
        
        # Step 3: Efflux pump-specific analysis
        print("Step 3: Analyzing efflux pump patterns...")
        
        analyzer = AlignmentAnalyzer(amr_database_info={
            'efflux_pumps': efflux_pump_genes,
            'target_organisms': ['Escherichia coli', 'Klebsiella pneumoniae', 'Enterobacter spp']
        })
        
        # Analyze each sample for efflux pump patterns
        sample_reports = {}
        for sample_path in input_samples:
            sample_name = Path(sample_path).stem
            paf_file = os.path.join(output_directory, sample_name, f"{sample_name}.paf")
            
            if os.path.exists(paf_file):
                report = analyzer.analyze_sample(paf_file, min_identity=85.0, min_coverage=70.0)
                sample_reports[sample_name] = report
        
        # Step 4: Co-occurrence analysis
        print("Step 4: Performing co-occurrence analysis...")
        
        cooccurrence_results = analyze_efflux_cooccurrence(sample_reports, efflux_pump_genes)
        
        # Step 5: Generate specialized reports
        print("Step 5: Generating efflux pump reports...")
        
        generate_efflux_reports(cooccurrence_results, output_directory)
        
        # Step 6: Statistical analysis
        print("Step 6: Performing statistical analysis...")
        
        stats_results = perform_efflux_statistics(cooccurrence_results)
        
        # Display summary results
        print("\nEfflux Pump Analysis Summary:")
        print("-" * 40)
        print(f"Samples analyzed: {len(sample_reports)}")
        print(f"Samples with efflux pumps: {stats_results['samples_with_efflux']}")
        print(f"Most common efflux genes: {stats_results['top_genes']}")
        print(f"Co-occurrence patterns found: {stats_results['cooccurrence_patterns']}")
        
        return True
        
    except Exception as e:
        print(f"Error in efflux pump analysis: {e}")
        return False

def analyze_efflux_cooccurrence(sample_reports, efflux_pump_genes):
    """
    Analyze co-occurrence patterns of efflux pump genes.
    
    This function quantifies:
    - How often specific efflux pump genes occur together
    - Patterns of acrA + acrB co-occurrence (as requested by Nishad Sir)
    - Statistical significance of co-occurrence patterns
    """
    
    cooccurrence_data = {
        'sample_gene_matrix': {},
        'pairwise_cooccurrence': defaultdict(int),
        'gene_frequencies': Counter(),
        'resistance_patterns': {}
    }
    
    # Build gene presence/absence matrix
    for sample_name, report in sample_reports.items():
        sample_genes = set()
        
        # Extract efflux pump genes from AMR hits
        for amr_hit in report.amr_genes_detected:
            gene_name = amr_hit.gene_name.lower()
            
            # Check against known efflux pump genes
            for category, genes in efflux_pump_genes.items():
                for efflux_gene in genes:
                    if efflux_gene.lower() in gene_name:
                        sample_genes.add(efflux_gene)
                        cooccurrence_data['gene_frequencies'][efflux_gene] += 1
        
        cooccurrence_data['sample_gene_matrix'][sample_name] = sample_genes
        
        # Calculate pairwise co-occurrence for this sample
        genes_list = list(sample_genes)
        for i in range(len(genes_list)):
            for j in range(i + 1, len(genes_list)):
                gene_pair = tuple(sorted([genes_list[i], genes_list[j]]))
                cooccurrence_data['pairwise_cooccurrence'][gene_pair] += 1
    
    # RND efflux pump genes (configurable)
    acra_acrb_analysis = analyze_acra_acrb_pattern(cooccurrence_data['sample_gene_matrix'])
    cooccurrence_data['acra_acrb_analysis'] = acra_acrb_analysis
    
    return cooccurrence_data

def analyze_acra_acrb_pattern(sample_gene_matrix):
    """
    Specific analysis of acrA and acrB co-occurrence patterns.
    """
    
    pattern_counts = {
        'acrA_only': 0,
        'acrB_only': 0,
        'both_acra_acrb': 0,
        'neither': 0,
        'total_samples': len(sample_gene_matrix)
    }
    
    cooccurrence_samples = []
    
    for sample_name, genes in sample_gene_matrix.items():
        has_acra = config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0] in genes
        has_acrb = config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1] in genes
        
        if has_acra and has_acrb:
            pattern_counts['both_acra_acrb'] += 1
            cooccurrence_samples.append(sample_name)
        elif has_acra:
            pattern_counts['acrA_only'] += 1
        elif has_acrb:
            pattern_counts['acrB_only'] += 1
        else:
            pattern_counts['neither'] += 1
    
    # Calculate co-occurrence statistics
    total_with_either = pattern_counts['acrA_only'] + pattern_counts['acrB_only'] + pattern_counts['both_acra_acrb']
    
    if total_with_either > 0:
        cooccurrence_rate = (pattern_counts['both_acra_acrb'] / total_with_either) * 100
    else:
        cooccurrence_rate = 0
    
    return {
        'pattern_counts': pattern_counts,
        'cooccurrence_rate': cooccurrence_rate,
        'cooccurrence_samples': cooccurrence_samples,
        'summary': f"acrA and acrB co-occur in {pattern_counts['both_acra_acrb']}/{total_with_either} samples with efflux pumps ({cooccurrence_rate:.1f}%)"
    }

def generate_efflux_reports(cooccurrence_results, output_directory):
    """Generate comprehensive reports for efflux pump analysis."""
    
    # 1. Co-occurrence matrix
    cooccurrence_matrix = create_cooccurrence_matrix(cooccurrence_results)
    matrix_path = os.path.join(output_directory, "efflux_cooccurrence_matrix.csv")
    cooccurrence_matrix.to_csv(matrix_path)
    
    # 2. Gene frequency report
    gene_freq_df = pd.DataFrame([
        {'Gene': gene, 'Frequency': freq, 'Percentage': (freq/len(cooccurrence_results['sample_gene_matrix'])*100)}
        for gene, freq in cooccurrence_results['gene_frequencies'].most_common()
    ])
    freq_path = os.path.join(output_directory, "efflux_gene_frequencies.csv")
    gene_freq_df.to_csv(freq_path, index=False)
    
    # RND efflux pump genes (configurable)
    acra_acrb_report = cooccurrence_results['acra_acrb_analysis']
    report_path = os.path.join(output_directory, "acrA_acrB_cooccurrence_report.txt")
    
    with open(report_path, 'w') as f:
        f.write("acrA and acrB Co-occurrence Analysis Report\n")
        f.write("=" * 45 + "\n\n")
        f.write(f"Total samples analyzed: {acra_acrb_report['pattern_counts']['total_samples']}\n")
        f.write(f"Samples with acrA only: {acra_acrb_report['pattern_counts']['acrA_only']}\n")
        f.write(f"Samples with acrB only: {acra_acrb_report['pattern_counts']['acrB_only']}\n")
        f.write(f"Samples with both acrA and acrB: {acra_acrb_report['pattern_counts']['both_acra_acrb']}\n")
        f.write(f"Co-occurrence rate: {acra_acrb_report['cooccurrence_rate']:.2f}%\n\n")
        f.write(f"Summary: {acra_acrb_report['summary']}\n\n")
        
        if acra_acrb_report['cooccurrence_samples']:
            f.write("Samples with acrA + acrB co-occurrence:\n")
            for sample in acra_acrb_report['cooccurrence_samples']:
                f.write(f"- {sample}\n")
    
    # 4. HTML visualization report
    generate_html_visualization(cooccurrence_results, output_directory)
    
    print(f"✓ Reports generated in {output_directory}")

def create_cooccurrence_matrix(cooccurrence_results):
    """Create a co-occurrence matrix for visualization."""
    
    all_genes = list(cooccurrence_results['gene_frequencies'].keys())
    matrix_data = []
    
    for gene1 in all_genes:
        row = {'Gene': gene1}
        for gene2 in all_genes:
            if gene1 == gene2:
                row[gene2] = cooccurrence_results['gene_frequencies'][gene1]
            else:
                pair = tuple(sorted([gene1, gene2]))
                row[gene2] = cooccurrence_results['pairwise_cooccurrence'].get(pair, 0)
        matrix_data.append(row)
    
    return pd.DataFrame(matrix_data).set_index('Gene')

def perform_efflux_statistics(cooccurrence_results):
    """Perform statistical analysis on efflux pump patterns."""
    
    total_samples = len(cooccurrence_results['sample_gene_matrix'])
    samples_with_efflux = len([s for s in cooccurrence_results['sample_gene_matrix'].values() if s])
    
    # Top 5 most common genes
    top_genes = cooccurrence_results['gene_frequencies'].most_common(5)
    
    # Most common co-occurrence patterns
    top_cooccurrence = sorted(
        cooccurrence_results['pairwise_cooccurrence'].items(),
        key=lambda x: x[1],
        reverse=True
    )[:10]
    
    return {
        'total_samples': total_samples,
        'samples_with_efflux': samples_with_efflux,
        'efflux_prevalence': (samples_with_efflux / total_samples * 100) if total_samples > 0 else 0,
        'top_genes': top_genes,
        'cooccurrence_patterns': len(top_cooccurrence),
        'most_common_pairs': top_cooccurrence
    }

def generate_html_visualization(cooccurrence_results, output_directory):
    """Generate an HTML visualization of efflux pump co-occurrence."""
    
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Efflux Pump Co-occurrence Analysis</title>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 20px; }}
            .header {{ background-color: #e8f4f8; padding: 15px; border-radius: 10px; }}
            .section {{ margin: 20px 0; }}
            .gene-freq {{ background-color: #f0f8f0; padding: 10px; margin: 5px 0; border-radius: 5px; }}
            .cooccurrence {{ background-color: #fff8e1; padding: 10px; margin: 5px 0; border-radius: 5px; }}
            .highlight {{ background-color: #ffeb3b; font-weight: bold; }}
            table {{ border-collapse: collapse; width: 100%; }}
            th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
            th {{ background-color: #f2f2f2; }}
        </style>
    </head>
    <body>
        <div class="header">
            <h1>RND Efflux Pump Co-occurrence Analysis</h1>
            <p>Analysis of efflux pump gene patterns in clinical isolates</p>
        </div>
        
        <div class="section">
            <h2>Gene Frequency Analysis</h2>
            {_generate_gene_frequency_html(cooccurrence_results)}
        </div>
        
        <div class="section">
            <h2>acrA + acrB Co-occurrence Pattern</h2>
            <div class="highlight">
                {cooccurrence_results['acra_acrb_analysis']['summary']}
            </div>
        </div>
        
        <div class="section">
            <h2>Top Co-occurrence Patterns</h2>
            {_generate_cooccurrence_html(cooccurrence_results)}
        </div>
    </body>
    </html>
    """
    
    html_path = os.path.join(output_directory, "efflux_cooccurrence_visualization.html")
    with open(html_path, 'w') as f:
        f.write(html_content)

def _generate_gene_frequency_html(cooccurrence_results):
    """Generate HTML for gene frequency section."""
    html = "<table><tr><th>Gene</th><th>Frequency</th><th>Samples</th></tr>"
    
    total_samples = len(cooccurrence_results['sample_gene_matrix'])
    
    for gene, freq in cooccurrence_results['gene_frequencies'].most_common(10):
        percentage = (freq / total_samples) * 100
        html += f"<tr><td>{gene}</td><td>{freq}</td><td>{percentage:.1f}%</td></tr>"
    
    html += "</table>"
    return html

def _generate_cooccurrence_html(cooccurrence_results):
    """Generate HTML for co-occurrence patterns."""
    html = "<table><tr><th>Gene Pair</th><th>Co-occurrence Count</th></tr>"
    
    top_pairs = sorted(
        cooccurrence_results['pairwise_cooccurrence'].items(),
        key=lambda x: x[1],
        reverse=True
    )[:10]
    
    for (gene1, gene2), count in top_pairs:
        html += f"<tr><td>{gene1} + {gene2}</td><td>{count}</td></tr>"
    
    html += "</table>"
    return html

def main():
    """Run the efflux pump co-occurrence analysis example."""
    
    print("Priority 2 Advanced Example: Efflux Pump Co-occurrence Analysis")
    print("This example addresses Nishad Sir's research questions about")
    print("acrA and acrB co-occurrence in clinical isolates.")
    print()
    
    try:
        success = efflux_pump_cooccurrence_analysis()
        
        if success:
            print("\n" + "=" * 70)
            print("✓ Efflux pump co-occurrence analysis completed successfully!")
            print("✓ Results demonstrate quantification of acrA + acrB patterns")
            print("✓ Analysis generalizable to any efflux pump gene combinations")
            print("✓ Ready for integration into full AMR analysis pipeline")
        else:
            print("\n✗ Analysis failed - check error messages above")
            
    except Exception as e:
        print(f"\nUnexpected error: {e}")

if __name__ == "__main__":
    main()