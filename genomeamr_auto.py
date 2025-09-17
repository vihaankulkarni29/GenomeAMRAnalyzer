#!/usr/bin/env python3
"""
GenomeAMRAnalyzer: Zero-Setup AMR Analysis Tool

Single command to run complete AMR analysis:
1. Automatically installs required tools (RGI, CARD DB)
2. Resolves NCBI URLs to genome accessions  
3. Downloads genomes and runs full pipeline
4. Generates comprehensive HTML reports

Usage:
    python genomeamr_auto.py --url "NCBI_URL" --genes config/genes_erythromycin.txt --email your@email.com
    python genomeamr_auto.py --accessions my_genomes.txt --genes config/genes_default.txt --email your@email.com
"""

import sys
import argparse
from pathlib import Path

# Add src to path for imports
project_root = Path(__file__).resolve().parent
sys.path.insert(0, str(project_root / "src"))

def main():
    parser = argparse.ArgumentParser(
        description="GenomeAMRAnalyzer: Zero-Setup AMR Analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze erythromycin resistance from NCBI search
  python genomeamr_auto.py \\
    --url "https://www.ncbi.nlm.nih.gov/nuccore/(Escherichia%%20coli%%20AND%%20erythromycin%%20resistance)" \\
    --genes config/genes_erythromycin.txt \\
    --email your.name@institution.edu

  # Analyze your genome list with default genes
  python genomeamr_auto.py \\
    --accessions my_genomes.txt \\
    --genes config/genes_default.txt \\
    --email your.name@institution.edu

Pre-built gene lists:
  config/genes_default.txt      - Broad AMR analysis (35+ genes)
  config/genes_erythromycin.txt - Macrolide resistance
  config/genes_betalactam.txt   - Beta-lactamase resistance
        """
    )
    
    # Input sources
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("--url", help="NCBI nuccore URL to analyze")
    input_group.add_argument("--accessions", help="File with genome accessions (one per line)")
    
    # Required parameters
    parser.add_argument("--genes", required=True, help="Gene list file (e.g., config/genes_erythromycin.txt)")
    parser.add_argument("--email", required=True, help="Your email for NCBI access")
    
    # Optional parameters
    parser.add_argument("--output", default="results", help="Output directory (default: results)")
    parser.add_argument("--config", default="config/snakemake_config.yaml", help="Pipeline config file")
    parser.add_argument("--max-genomes", type=int, default=50, help="Maximum genomes to analyze")
    parser.add_argument("--skip-install", action="store_true", help="Skip automatic tool installation")
    
    # Parse arguments
    args = parser.parse_args()
    
    # Validate inputs
    if not Path(args.genes).exists():
        print(f"❌ Gene list file not found: {args.genes}")
        return 1
    
    if not Path(args.config).exists():
        print(f"❌ Config file not found: {args.config}")
        return 1
    
    if args.accessions and not Path(args.accessions).exists():
        print(f"❌ Accessions file not found: {args.accessions}")
        return 1
    
    # Show banner
    print("🧬 GenomeAMRAnalyzer: Zero-Setup AMR Analysis")
    print("=" * 50)
    
    # Show analysis parameters
    input_source = args.url if args.url else args.accessions
    input_type = "NCBI URL" if args.url else "Accession file"
    
    print(f"📝 Input: {input_type}")
    print(f"   {input_source}")
    print(f"🎯 Genes: {Path(args.genes).name} ({count_genes(args.genes)} genes)")
    print(f"📧 Email: {args.email}")
    print(f"📁 Output: {args.output}/")
    print()
    
    # Import and run the main pipeline
    try:
        from run_pipeline import main as run_main
        
        # Build arguments for run_pipeline.py
        run_args = [
            "--config", args.config,
            "--genes-file", args.genes,
            "--email", args.email,
        ]
        
        if not args.skip_install:
            run_args.append("--auto-install")
        
        if args.url:
            run_args.extend(["--url", args.url])
        else:
            run_args.extend(["--accessions-file", args.accessions])
        
        # Override sys.argv to pass arguments to run_pipeline
        original_argv = sys.argv
        sys.argv = ["run_pipeline.py"] + run_args
        
        try:
            result = run_main()
            return result
        finally:
            sys.argv = original_argv
            
    except ImportError as e:
        print(f"❌ Failed to import pipeline: {e}")
        return 1
    except Exception as e:
        print(f"❌ Pipeline failed: {e}")
        return 1


def count_genes(genes_file: str) -> int:
    """Count number of genes in gene list file."""
    try:
        with open(genes_file) as f:
            return len([line.strip() for line in f if line.strip() and not line.startswith("#")])
    except:
        return 0


def show_quick_examples():
    """Show quick usage examples."""
    print("""
🚀 Quick Examples:

# 1. Erythromycin resistance analysis
python genomeamr_auto.py \\
  --url "https://www.ncbi.nlm.nih.gov/nuccore/(Escherichia coli AND erythromycin resistance)" \\
  --genes config/genes_erythromycin.txt \\
  --email your@email.com

# 2. Beta-lactam resistance analysis  
python genomeamr_auto.py \\
  --url "https://www.ncbi.nlm.nih.gov/nuccore/(clinical isolates AND ESBL)" \\
  --genes config/genes_betalactam.txt \\
  --email your@email.com

# 3. Your own genome list
python genomeamr_auto.py \\
  --accessions my_genomes.txt \\
  --genes config/genes_default.txt \\
  --email your@email.com

📁 Gene Lists Available:
  config/genes_default.txt      - 35+ genes (broad AMR)
  config/genes_erythromycin.txt - Macrolide resistance
  config/genes_betalactam.txt   - Beta-lactamase resistance
    """)


if __name__ == "__main__":
    if len(sys.argv) == 1:
        show_quick_examples()
    else:
        sys.exit(main())