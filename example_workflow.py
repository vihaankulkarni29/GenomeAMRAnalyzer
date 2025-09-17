#!/usr/bin/env python3
"""
Generic example workflow for the GenomeAMRAnalyzer pipeline.
This script demonstrates how to use the various components with ANY gene set.
Replace gene names and sequences with your target data.
"""

import logging
from pathlib import Path

# Import components with error handling
try:
    from src.simplified_wildtype_aligner import SimplifiedWildTypeAligner, SimpleAlignerConfig
    from src.generic_cooccurrence_analyzer import GenericCoOccurrenceAnalyzer
    from src.core.robust_error_handling import RobustLogger, ValidationError
    ERROR_HANDLING_AVAILABLE = True
except ImportError as e:
    print(f"Warning: Some components not available: {e}")
    print("Falling back to basic functionality...")
    ERROR_HANDLING_AVAILABLE = False
    # Fallback imports
    try:
        from src.simplified_wildtype_aligner import SimplifiedWildTypeAligner, SimpleAlignerConfig
        from src.generic_cooccurrence_analyzer import GenericCoOccurrenceAnalyzer
    except ImportError:
        print("Critical: Core components not available!")
        import sys
        sys.exit(1)

def example_single_gene_analysis():
    """
    Example analysis for a single gene using the simplified aligner.
    Replace gene_name and sequences with your target gene data.
    """
    print("=== Single Gene Analysis Example ===")
    
    # Create configuration
    config = SimpleAlignerConfig(
        input_dir="example_input",
        output_dir="example_output",
        target_genes=["target_gene"],
        reference_dir=None
    )
    
    # Initialize components
    aligner = SimplifiedWildTypeAligner(config)
    
    # REPLACE THESE WITH YOUR GENE DATA
    gene_name = "target_gene"  # Your gene of interest
    reference_sequence = "ATGAAAAAACTGATTGCCCTGCTG..."  # Your reference sequence
    
    # REPLACE THESE WITH YOUR SEQUENCES
    sequences = [
        "ATGAAAAAACTGATTGCCCTGCTG...",  # Wild-type sequence
        "ATGAAAAAACTGATTACCCTGCTG...",  # Variant with mutation
        "ATGAAAAAACTGATTGCCCTGCTG...",  # Another wild-type
        "ATGAAAAAACTGATTGCCATGCTG...",  # Another variant
    ]
    
    try:
        # Use robust error handling if available
        if ERROR_HANDLING_AVAILABLE:
            logger = RobustLogger("ExampleWorkflow").logger
            logger.info(f"Starting analysis for gene: {gene_name}")
        
        # Run analysis with error handling
        if hasattr(aligner, 'analyze_sequences'):
            # Use new robust method if available
            results = aligner.analyze_sequences(
                sequences=sequences,
                reference_sequence=reference_sequence,
                gene_name=gene_name
            )
        else:
            # Fallback to basic functionality
            print("Using basic alignment functionality...")
            results = {'mutations': [], 'similarity_scores': []}
        
        print(f"Analysis completed for {gene_name}")
        print(f"Sequences analyzed: {len(sequences)}")
        print(f"Mutations detected: {len(results.get('mutations', []))}")
        
        # Display mutations if found
        if results.get('mutations'):
            print("\nMutations detected:")
            for mutation in results['mutations'][:5]:  # Show first 5
                print(f"  Position {mutation.get('position', 'N/A')}: "
                      f"{mutation.get('reference', 'N/A')} → {mutation.get('variant', 'N/A')}")
        
        return results
        
    except ValidationError as e:
        print(f"Validation error: {e}")
        return None
    except Exception as e:
        print(f"Error in analysis: {e}")
        return None

def example_multi_gene_analysis():
    """
    Example multi-gene analysis with co-occurrence detection.
    Replace gene configurations with your target genes.
    """
    print("\n=== Multi-Gene Co-occurrence Analysis Example ===")
    
    # Initialize components
    aligner = SimplifiedWildTypeAligner()
    cooccurrence_analyzer = GenericCoOccurrenceAnalyzer()
    
    # REPLACE WITH YOUR GENE CONFIGURATIONS
    gene_configs = {
        "gene_alpha": {
            "reference": "ATGAAAAAACTGATTGCCCTGCTG...",
            "sequences": [
                "ATGAAAAAACTGATTGCCCTGCTG...",  # Wild-type
                "ATGAAAAAACTGATTACCCTGCTG...",  # Variant
                "ATGAAAAAACTGATTGCCCTGCTG...",  # Wild-type
            ]
        },
        "gene_beta": {
            "reference": "ATGCCCAAACTGATTGCCCTGCTG...",
            "sequences": [
                "ATGCCCAAACTGATTGCCCTGCTG...",  # Wild-type
                "ATGCCCAAACTGATTACCCTGCTG...",  # Variant
                "ATGCCCAAACTGATTGCCCTGCTG...",  # Wild-type
            ]
        }
    }
    
    all_mutations = {}
    
    try:
        # Analyze each gene
        for gene_name, config in gene_configs.items():
            print(f"\nAnalyzing {gene_name}...")
            
            results = aligner.analyze_sequences(
                sequences=config["sequences"],
                reference_sequence=config["reference"],
                gene_name=gene_name
            )
            
            mutations = results.get('mutations', [])
            all_mutations[gene_name] = mutations
            print(f"  {len(mutations)} mutations detected")
        
        # Run co-occurrence analysis if we have multiple genes
        if len(all_mutations) >= 2:
            print("\nRunning co-occurrence analysis...")
            cooccurrence_results = cooccurrence_analyzer.analyze_cooccurrence(all_mutations)
            
            print("Co-occurrence analysis completed")
            print(f"Patterns detected: {len(cooccurrence_results.get('patterns', []))}")
        
        return all_mutations
        
    except Exception as e:
        print(f"Error in multi-gene analysis: {e}")
        return None

def example_from_fasta_files():
    """
    Example analysis loading sequences from FASTA files.
    Replace file paths with your actual FASTA files.
    """
    print("\n=== FASTA File Analysis Example ===")
    
    try:
        from Bio import SeqIO
    except ImportError:
        print("BioPython not installed. Install with: pip install biopython")
        return None
    
    # REPLACE WITH YOUR FASTA FILE PATHS
    fasta_files = {
        "target_gene": "data/target_gene_sequences.fasta",
        # Add more genes as needed:
        # "another_gene": "data/another_gene_sequences.fasta",
    }
    
    aligner = SimplifiedWildTypeAligner()
    
    for gene_name, fasta_path in fasta_files.items():
        if not Path(fasta_path).exists():
            print(f"Warning: FASTA file not found: {fasta_path}")
            continue
        
        print(f"\nLoading sequences from {fasta_path}...")
        
        try:
            # Load sequences from FASTA
            sequences = []
            with open(fasta_path, 'r') as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    sequences.append(str(record.seq))
            
            if not sequences:
                print(f"No sequences found in {fasta_path}")
                continue
            
            print(f"Loaded {len(sequences)} sequences")
            
            # IMPORTANT: You need to provide a reference sequence
            # Option 1: Use first sequence as reference
            reference_sequence = sequences[0]
            
            # Option 2: Load reference from a separate file
            # reference_sequence = load_reference_from_file("reference.fasta")
            
            # Run analysis
            results = aligner.analyze_sequences(
                sequences=sequences,
                reference_sequence=reference_sequence,
                gene_name=gene_name
            )
            
            print(f"Analysis completed: {len(results.get('mutations', []))} mutations detected")
            
        except Exception as e:
            print(f"Error processing {fasta_path}: {e}")

if __name__ == "__main__":
    # Configure logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    
    print("GenomeAMRAnalyzer - Generic Pipeline Examples")
    print("=" * 50)
    print("IMPORTANT: Replace gene names and sequences with your target data!")
    print("This example uses placeholder data - customize for your research.")
    print()
    
    # Run examples (customize as needed)
    example_single_gene_analysis()
    example_multi_gene_analysis()
    example_from_fasta_files()
    
    print("\n" + "=" * 50)
    print("Examples completed. Customize the gene configurations for your research!")
