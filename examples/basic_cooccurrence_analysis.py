import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))
from configuration_manager import config_manager
#!/usr/bin/env python3
"""
Example: Basic Co-occurrence Analysis
=====================================

This example demonstrates how to use the Generic Co-occurrence Analyzer
to analyze mutation patterns across RND efflux pump genes.
"""

import sys
import os
from pathlib import Path
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))

from generic_cooccurrence_analyzer import GenericCoOccurrenceAnalyzer, MutationEvent

def create_sample_data():
    """Create sample mutation data for demonstration"""
    
    # Sample mutations across different genes and genomes
    mutations = [
        # Genome 1: mdtF and acrA mutations
        MutationEvent(
            genome_id="GCF_000005825.2",
            gene="mdtF",
            position=123,
            reference_aa="A",
            variant_aa="V",
            substitution="A123V"
        ),
        MutationEvent(
            genome_id="GCF_000005825.2",
            gene=config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0],
            position=456,
            reference_aa="R",
            variant_aa="H",
            substitution="R456H"
        ),
        
        # Genome 2: mdtF mutation only
        MutationEvent(
            genome_id="GCF_000006945.2",
            gene="mdtF",
            position=123,
            reference_aa="A",
            variant_aa="V",
            substitution="A123V"
        ),
        
        # RND efflux pump genes (configurable)
        MutationEvent(
            genome_id="GCF_000007123.1",
            gene=config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0],
            position=456,
            reference_aa="R",
            variant_aa="H",
            substitution="R456H"
        ),
        MutationEvent(
            genome_id="GCF_000007123.1",
            gene=config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1],
            position=789,
            reference_aa="D",
            variant_aa="N",
            substitution="D789N"
        ),
        
        # Genome 4: All three genes mutated (high resistance?)
        MutationEvent(
            genome_id="GCF_000008456.1",
            gene="mdtF",
            position=123,
            reference_aa="A",
            variant_aa="V",
            substitution="A123V"
        ),
        MutationEvent(
            genome_id="GCF_000008456.1",
            gene=config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0],
            position=456,
            reference_aa="R",
            variant_aa="H",
            substitution="R456H"
        ),
        MutationEvent(
            genome_id="GCF_000008456.1",
            gene=config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1],
            position=789,
            reference_aa="D",
            variant_aa="N",
            substitution="D789N"
        ),
    ]
    
    return mutations

def main():
    """Run basic co-occurrence analysis example"""
    
    print("GenomeAMRAnalyzer - Basic Co-occurrence Analysis Example")
    print("=" * 60)
    
    # Define genes of interest (RND efflux pump components)
    genes_of_interest = ["mdtF", config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0], config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1]]
    
    # Initialize analyzer
    print(f"\n1. Initializing analyzer for genes: {genes_of_interest}")
    analyzer = GenericCoOccurrenceAnalyzer(output_dir="example_results")
    
    # Create and load sample data
    print("\n2. Creating sample mutation data...")
    mutations = create_sample_data()
    
    # Add mutations to analyzer
    print(f"   Loading {len(mutations)} mutations...")
    analyzer.mutations_data = mutations
    
    # Update gene list and genome set
    for mutation in mutations:
        analyzer.gene_list.add(mutation.gene)
        analyzer.genome_set.add(mutation.genome_id)
    
    # Update stats
    analyzer.stats['total_mutations'] = len(mutations)
    analyzer.stats['total_genomes'] = len(analyzer.genome_set)
    analyzer.stats['total_genes'] = len(analyzer.gene_list)
    
    # Display mutation summary
    print(f"\n3. Mutation Summary:")
    print(f"   Total genomes analyzed: {analyzer.stats['total_genomes']}")
    print(f"   Total mutations found: {analyzer.stats['total_mutations']}")
    print(f"   Genes with mutations: {list(analyzer.gene_list)}")
    
    print(f"\n   Mutations per gene:")
    gene_counts = {}
    for mutation in mutations:
        gene_counts[mutation.gene] = gene_counts.get(mutation.gene, 0) + 1
    
    for gene, count in gene_counts.items():
        print(f"     {gene}: {count} mutations")
    
    # Analyze co-occurrence patterns
    print(f"\n4. Analyzing co-occurrence patterns...")
    analyzer.analyze_cooccurrence_patterns()
    
    # Generate detailed report
    print(f"\n5. Generating detailed analysis report...")
    results = analyzer.generate_detailed_report()
    
    # Display key results
    print(f"\n6. Analysis Results:")
    print(f"   Patterns found: {analyzer.stats.get('patterns_found', 0)}")
    print(f"   Significant patterns: {analyzer.stats.get('significant_patterns', 0)}")
    
    if analyzer.patterns:
        print(f"\n   Co-occurrence patterns detected:")
        for i, pattern in enumerate(analyzer.patterns[:3], 1):  # Show first 3 patterns
            genes = list(pattern.genes)
            print(f"     Pattern {i}: {genes}")
            print(f"       Genomes: {len(pattern.genomes)}")
            print(f"       Frequency: {pattern.frequency:.3f}")
    
    # Save results
    print(f"\n7. Saving results...")
    output_files = analyzer.save_results("example_analysis")
    
    print(f"   Results saved to:")
    for file_type, file_path in output_files.items():
        print(f"     {file_type}: {Path(file_path).name}")
    
    print(f"\n8. Analysis Complete!")
    print(f"   Results suggest potential coordinated resistance mechanisms")
    print(f"   involving {len(genes_of_interest)} RND efflux pump components.")

if __name__ == "__main__":
    main()