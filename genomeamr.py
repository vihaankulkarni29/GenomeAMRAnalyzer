#!/usr/bin/env python3
"""
GenomeAMRAnalyzer - User-Friendly CLI Interface

A beginner-friendly command-line interface for antimicrobial resistance gene analysis.
Provides guided workflows, example data, and easy setup validation.

Usage:
    python genomeamr.py --quick-test                    # Test with sample data
    python genomeamr.py --examples                      # Show usage examples
    python genomeamr.py --check-setup                   # Validate environment
    python genomeamr    # Print header
    print("🔬 GenomeAMRAnalyzer - Antimicrobial Resistance Gene Analysis")
    print("=" * 65)
    print("🧪 Production-ready pipeline for bacterial genome AMR analysis")
    print("✨ Beginner-friendly interface with guided workflows")
    print()tutorial                      # Interactive tutorial
    python genomeamr.py --accessions GCF_000005825.2    # Analyze specific genome
"""

import argparse
import os
import sys
import subprocess
import tempfile
import webbrowser
from pathlib import Path
from typing import List, Optional

# Add src directory to path for imports
sys.path.insert(0, str(Path(__file__).parent / "src"))

def check_environment() -> bool:
    """Check if the environment is properly set up."""
    try:
        print("🔍 Checking environment setup...")
    except UnicodeEncodeError:
        print("Checking environment setup...")
    
    required_packages = [
        "numpy", "pandas", "matplotlib", "seaborn", "plotly", 
        "Bio", "requests", "yaml", "networkx"
    ]
    
    package_names = {
        "Bio": "biopython",
        "yaml": "pyyaml"
    }
    
    missing_packages = []
    for package in required_packages:
        try:
            __import__(package)
            display_name = package_names.get(package, package)
            try:
                print(f"✅ {display_name}")
            except UnicodeEncodeError:
                print(f"OK {display_name}")
        except ImportError:
            display_name = package_names.get(package, package)
            missing_packages.append(display_name)
            try:
                print(f"❌ {display_name} - MISSING")
            except UnicodeEncodeError:
                print(f"MISSING {display_name}")
    
    if missing_packages:
        try:
            print(f"\n⚠️  Missing packages: {', '.join(missing_packages)}")
            print("📥 Install with: pip install -r requirements.txt")
        except UnicodeEncodeError:
            print(f"\nMissing packages: {', '.join(missing_packages)}")
            print("Install with: pip install -r requirements.txt")
        return False
    
    try:
        print("\n✅ All required packages are available!")
    except UnicodeEncodeError:
        print("\nAll required packages are available!")
    return True

def create_sample_accessions() -> str:
    """Create a temporary file with sample accession numbers."""
    sample_accessions = [
        "GCF_000005825.2",  # E. coli K-12 MG1655
        "GCF_000006945.2",  # Salmonella enterica subsp. enterica serovar Typhimurium str. LT2
    ]
    
    temp_file = tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False)
    for accession in sample_accessions:
        temp_file.write(f"{accession}\n")
    temp_file.close()
    
    return temp_file.name

def run_quick_test() -> bool:
    """Run a quick test with sample data."""
    print("🚀 Running quick test with sample genomes...")
    print("📋 Using E. coli K-12 MG1655 and Salmonella Typhimurium LT2")
    
    # Check if main pipeline exists
    pipeline_script = Path(__file__).parent / "run_pipeline.py"
    if not pipeline_script.exists():
        print("❌ Main pipeline script not found!")
        return False
    
    # Create sample accession file
    accession_file = create_sample_accessions()
    
    # Create test output directory
    output_dir = Path("quick_test_output")
    output_dir.mkdir(exist_ok=True)
    
    try:
        # Run the pipeline with sample data
        cmd = [
            sys.executable, str(pipeline_script),
            "--config", "config/snakemake_config.yaml",
            "--output-dir", str(output_dir),
            "--mock-mode"  # Use mock mode for quick testing
        ]
        
        print(f"🔧 Running: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        
        if result.returncode == 0:
            print("✅ Quick test completed successfully!")
            print(f"📊 Results saved to: {output_dir}")
            
            # Check if HTML report was generated
            html_report = output_dir / "reports" / "pipeline_report.html"
            if html_report.exists():
                print(f"🌐 HTML report: {html_report}")
                print("💡 Open the HTML file in your browser to view results")
                
                # Try to open report automatically
                try:
                    import webbrowser
                    webbrowser.open(f"file://{html_report.absolute()}")
                    print("🌐 Opening report in default browser...")
                except Exception:
                    print("💡 Manually open the HTML file to view results")
            
            return True
        else:
            print("❌ Quick test failed!")
            print(f"Error: {result.stderr}")
            return False
            
    except subprocess.TimeoutExpired:
        print("⏰ Quick test timed out - this might indicate installation issues")
        return False
    except Exception as e:
        print(f"❌ Error during quick test: {e}")
        return False
    finally:
        # Clean up temporary file
        if os.path.exists(accession_file):
            os.unlink(accession_file)

def show_examples():
    """Display usage examples for different scenarios."""
    examples = {
        "🔬 **Research Scientists**": [
            "# Publication-ready analysis of E. coli strains",
            "python genomeamr.py --accessions GCF_000005825.2,GCF_000006945.2",
            "",
            "# Batch analysis with custom gene focus",
            "echo 'GCF_000005825.2' > my_genomes.txt",
            "python run_pipeline.py --config config/snakemake_config.yaml",
            "",
            "# Advanced co-occurrence network analysis",
            "python genomeamr.py --accessions GCF_000005825.2 --network-analysis",
        ],
        
        "🎓 **Students & Educators**": [
            "# Interactive tutorial with step-by-step guidance",
            "python genomeamr.py --tutorial",
            "",
            "# Quick demonstration with sample data",
            "python genomeamr.py --quick-test",
            "",
            "# Learn about different analysis steps",
            "python genomeamr.py --explain-pipeline",
        ],
        
        "🏥 **Clinical Laboratories**": [
            "# Fast screening for clinical isolates",
            "python genomeamr.py --clinical-mode --genes acrA,acrB,tolC",
            "",
            "# Batch analysis for outbreak investigation",
            "python run_pipeline.py --config clinical_config.yaml",
            "",
            "# Resistance profiling with clinical focus",
            "python genomeamr.py --accessions YOUR_ISOLATE --clinical-report",
        ],
        
        "💻 **Bioinformaticians**": [
            "# Full pipeline with custom configuration",
            "python run_pipeline.py --config custom_config.yaml",
            "",
            "# Programmatic analysis with Python API",
            "from src.generic_cooccurrence_analyzer import CooccurrenceAnalyzer",
            "analyzer = CooccurrenceAnalyzer()",
            "",
            "# High-throughput processing",
            "python run_pipeline.py --batch-mode --parallel 8",
        ]
    }
    
    print("📚 GenomeAMRAnalyzer - Usage Examples by User Type\n")
    print("=" * 60)
    
    for category, commands in examples.items():
        print(f"\n{category}")
        print("-" * (len(category) - 4))  # Adjust for emoji
        for cmd in commands:
            if cmd.startswith("#"):
                print(f"  💡 {cmd[2:]}")
            elif cmd.strip():
                print(f"  $ {cmd}")
            else:
                print()
    
    print("\n" + "=" * 60)
    print("🚀 **Getting Started Recommendations:**")
    print("  1. New users: python genomeamr.py --tutorial")
    print("  2. Quick test: python genomeamr.py --quick-test") 
    print("  3. Check setup: python genomeamr.py --check-setup")
    print("  4. See use cases: Read USE_CASES.md and QUICKSTART.md")

def run_tutorial():
    """Run an interactive tutorial."""
    print("🎓 Welcome to GenomeAMRAnalyzer Tutorial!")
    print("=" * 50)
    print("This tutorial will guide you through:")
    print("  • Environment validation")
    print("  • Running your first analysis") 
    print("  • Understanding the results")
    print("  • Next steps for your research")
    
    steps = [
        {
            "title": "Step 1: Environment Check",
            "description": "Let's make sure everything is set up correctly",
            "action": lambda: check_environment()
        },
        {
            "title": "Step 2: Quick Test Analysis",
            "description": "Run a test analysis with sample genomes (E. coli & Salmonella)",
            "action": lambda: run_quick_test()
        },
        {
            "title": "Step 3: Understanding Results",
            "description": "Learn how to interpret the analysis output",
            "action": lambda: explain_results()
        },
        {
            "title": "Step 4: Next Steps",
            "description": "Plan your own research analysis",
            "action": lambda: suggest_next_steps()
        }
    ]
    
    for i, step in enumerate(steps, 1):
        print(f"\n{'='*20} {step['title']} {'='*20}")
        print(step['description'])
        
        response = input(f"\n🚀 Ready for step {i}? (y/n/q): ").lower().strip()
        if response in ['y', 'yes']:
            success = step['action']()
            if not success and i < len(steps):
                print("\n⚠️  Issue encountered. Continue anyway? (y/n): ", end="")
                if input().lower().strip() not in ['y', 'yes']:
                    break
        elif response in ['q', 'quit']:
            print("Tutorial paused. Run again anytime with: python genomeamr.py --tutorial")
            break
        else:
            print("Skipping this step...")
    
    print("\n🎉 Tutorial complete! You're ready to use GenomeAMRAnalyzer!")
    print("📚 Next: Check out USE_CASES.md for real-world examples")

def explain_results():
    """Explain how to interpret results."""
    print("\n📊 Understanding Your Results")
    print("=" * 30)
    
    print("📁 **Output Directory Structure:**")
    structure = [
        "results/",
        "├── card_analysis/          # 🧬 CARD resistance gene detection",
        "├── protein_extraction/     # 🔬 Extracted protein sequences", 
        "├── alignments/            # 📏 Sequence alignments with references",
        "├── mutations/             # 🧪 Statistical mutation analysis",
        "├── cooccurrence/          # 🕸️  Gene co-occurrence networks",
        "└── reports/               # 📈 Interactive HTML reports",
        "    └── pipeline_report.html  # 🌐 Main interactive dashboard"
    ]
    
    for line in structure:
        print(f"  {line}")
    
    print("\n🔍 **Key Analysis Components:**")
    explanations = [
        "🧬 **CARD Results**: Comprehensive resistance gene identification using CARD database",
        "🔬 **Protein Extraction**: High-quality protein sequences from genomic coordinates", 
        "📏 **Sequence Alignments**: EMBOSS Water-style pairwise alignments with references",
        "🧪 **Mutation Analysis**: Statistical significance testing for sequence variations",
        "🕸️  **Co-occurrence Networks**: Gene interaction patterns and resistance mechanisms",
        "📈 **Interactive Reports**: Publication-ready visualizations and statistical summaries"
    ]
    
    for explanation in explanations:
        print(f"  {explanation}")
    
    print("\n💡 **How to Use Your Results:**")
    usage_tips = [
        "📖 **Start with**: pipeline_report.html for overview and key findings",
        "🔬 **Research**: Dive into individual analysis folders for detailed data",
        "📊 **Publication**: Use generated figures and statistical tables directly",
        "🔄 **Iteration**: Modify parameters based on initial results",
        "📋 **Documentation**: Check logs/ folders for detailed processing records"
    ]
    
    for tip in usage_tips:
        print(f"  {tip}")
    
    return True

def suggest_next_steps():
    """Suggest next steps for users."""
    print("\n🚀 **Your Research Journey - Next Steps**")
    print("=" * 40)
    
    suggestions = {
        "📚 **Learning More**": [
            "• Read USE_CASES.md for real-world research examples",
            "• Check QUICKSTART.md for 5-minute setup guide", 
            "• Review docs/ folder for detailed documentation",
            "• Explore config/ for customization options"
        ],
        
        "🔬 **For Research**": [
            "• Prepare your own genome accession list",
            "• Customize config/snakemake_config.yaml for your study",
            "• Consider computational resources for large datasets",
            "• Plan validation with known positive/negative controls"
        ],
        
        "🏥 **For Clinical Use**": [
            "• Validate pipeline with characterized clinical isolates",
            "• Set up automated reporting workflows",
            "• Configure clinical-relevant gene panels",
            "• Establish quality control procedures"
        ],
        
        "🎓 **For Teaching**": [
            "• Use tutorial mode with students",
            "• Create custom example datasets",
            "• Develop exercises around different analysis steps",
            "• Compare results with literature examples"
        ]
    }
    
    for category, items in suggestions.items():
        print(f"\n{category}")
        for item in items:
            print(f"  {item}")
    
    print("\n💬 **Getting Help:**")
    help_options = [
        "🐛 **Issues**: Report bugs at GitHub Issues tracker",
        "❓ **Questions**: Check documentation or ask the community",
        "🤝 **Collaboration**: Contact for research partnerships",
        "📖 **Updates**: Watch GitHub repository for new features"
    ]
    
    for option in help_options:
        print(f"  {option}")
    
    return True

def run_analysis(accessions: List[str], output_dir: Optional[str] = None):
    """Run analysis for specified accessions."""
    print(f"🔬 Analyzing {len(accessions)} genome(s)...")
    print(f"📋 Accessions: {', '.join(accessions)}")
    
    # Validate accessions format
    for acc in accessions:
        if not acc.strip() or len(acc.strip()) < 10:
            print(f"⚠️  Warning: '{acc}' might not be a valid accession format")
    
    # Create accession file
    accession_file = tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False)
    for accession in accessions:
        accession_file.write(f"{accession.strip()}\n")
    accession_file.close()
    
    # Set output directory
    if not output_dir:
        output_dir = f"analysis_results_{len(accessions)}_genomes"
    
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    print(f"📁 Output directory: {output_path.absolute()}")
    
    try:
        # Run the main pipeline
        pipeline_script = Path(__file__).parent / "run_pipeline.py"
        cmd = [
            sys.executable, str(pipeline_script),
            "--config", "config/snakemake_config.yaml",
            "--output-dir", str(output_path)
        ]
        
        print(f"🔧 Running pipeline...")
        print("⏰ This may take several minutes depending on genome size and network speed")
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode == 0:
            print("✅ Analysis completed successfully!")
            print(f"📊 Results saved to: {output_path.absolute()}")
            
            # Look for HTML report
            html_report = output_path / "reports" / "pipeline_report.html"
            if html_report.exists():
                print(f"🌐 Interactive report: {html_report}")
                try:
                    import webbrowser
                    webbrowser.open(f"file://{html_report.absolute()}")
                    print("🌐 Opening report in default browser...")
                except Exception:
                    print("💡 Manually open the HTML file to view results")
            
            # Summarize what was found
            summarize_results(output_path)
            
        else:
            print("❌ Analysis failed!")
            print(f"Error: {result.stderr}")
            if result.stdout:
                print(f"Output: {result.stdout}")
            
    except Exception as e:
        print(f"❌ Error during analysis: {e}")
    finally:
        # Clean up
        if os.path.exists(accession_file.name):
            os.unlink(accession_file.name)

def summarize_results(output_path: Path):
    """Provide a quick summary of analysis results."""
    print("\n📋 **Quick Results Summary:**")
    
    try:
        # Check different result directories
        card_dir = output_path / "card_analysis" 
        if card_dir.exists():
            card_files = list(card_dir.glob("*.json"))
            print(f"  🧬 CARD Analysis: {len(card_files)} genome(s) processed")
        
        align_dir = output_path / "alignments"
        if align_dir.exists():
            align_files = list(align_dir.glob("*.txt"))
            print(f"  📏 Alignments: {len(align_files)} alignment file(s) generated")
        
        reports_dir = output_path / "reports"
        if reports_dir.exists():
            report_files = list(reports_dir.glob("*.html"))
            print(f"  📈 Reports: {len(report_files)} interactive report(s) created")
            
    except Exception as e:
        print(f"  ⚠️  Could not summarize results: {e}")
    
    print("  💡 Check the HTML report for detailed findings and visualizations")

def main():
    parser = argparse.ArgumentParser(
        description="GenomeAMRAnalyzer - User-Friendly AMR Analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
🚀 Quick Start Examples:
  python genomeamr.py --quick-test                     # Test with sample data (fastest)
  python genomeamr.py --tutorial                       # Interactive step-by-step guide
  python genomeamr.py --examples                       # Show usage examples by user type
  python genomeamr.py --accessions GCF_000005825.2     # Analyze E. coli K-12
  python genomeamr.py --check-setup                    # Validate your environment

📚 Documentation:
  • QUICKSTART.md - 5-minute setup guide
  • USE_CASES.md - Real-world research examples  
  • docs/ - Comprehensive documentation

🆘 Support: https://github.com/yourusername/GenomeAMRAnalyzer
        """
    )
    
    # Main action arguments (mutually exclusive)
    action_group = parser.add_mutually_exclusive_group(required=True)
    action_group.add_argument("--quick-test", action="store_true",
                            help="🚀 Run quick test with sample genomes (E. coli & Salmonella)")
    action_group.add_argument("--tutorial", action="store_true",
                            help="🎓 Interactive tutorial with step-by-step guidance")
    action_group.add_argument("--examples", action="store_true",
                            help="📚 Show usage examples for different user types")
    action_group.add_argument("--check-setup", action="store_true", 
                            help="🔍 Check if environment is properly configured")
    action_group.add_argument("--accessions", type=str,
                            help="🔬 Comma-separated genome accessions (e.g., GCF_000005825.2)")
    
    # Optional arguments
    parser.add_argument("--output-dir", type=str,
                       help="📁 Output directory for results (default: auto-generated)")
    parser.add_argument("--clinical-mode", action="store_true",
                       help="🏥 Use clinical-focused analysis parameters")
    parser.add_argument("--genes", type=str,
                       help="🎯 Focus on specific genes (e.g., acrA,acrB,tolC)")
    
    # Parse arguments
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
        
    args = parser.parse_args()
    
    # Print header with Windows-compatible characters
    try:
        print("🔬 GenomeAMRAnalyzer - Antimicrobial Resistance Gene Analysis")
    except UnicodeEncodeError:
        print("GenomeAMRAnalyzer - Antimicrobial Resistance Gene Analysis")
    
    print("=" * 65)
    
    try:
        print("🧪 Production-ready pipeline for bacterial genome AMR analysis")
        print("✨ Beginner-friendly interface with guided workflows")
    except UnicodeEncodeError:
        print("Production-ready pipeline for bacterial genome AMR analysis")
        print("Beginner-friendly interface with guided workflows")
    
    print()
    
    # Execute based on arguments
    if args.quick_test:
        if not check_environment():
            print("❌ Environment check failed. Install requirements first:")
            print("   pip install -r requirements.txt")
            sys.exit(1)
        success = run_quick_test()
        sys.exit(0 if success else 1)
        
    elif args.tutorial:
        run_tutorial()
        
    elif args.examples:
        show_examples()
        
    elif args.check_setup:
        success = check_environment()
        if success:
            try:
                print("\n🎉 Your environment is ready for GenomeAMRAnalyzer!")
                print("🚀 Next step: python genomeamr.py --quick-test")
            except UnicodeEncodeError:
                print("\nYour environment is ready for GenomeAMRAnalyzer!")
                print("Next step: python genomeamr.py --quick-test")
        sys.exit(0 if success else 1)
        
    elif args.accessions:
        if not check_environment():
            print("❌ Environment check failed. Install requirements first:")
            print("   pip install -r requirements.txt")
            sys.exit(1)
        
        accession_list = [acc.strip() for acc in args.accessions.split(",")]
        run_analysis(accession_list, args.output_dir)
    
    print("\n" + "="*65)
    print("💡 **Need help?**")
    print("   • Tutorial: python genomeamr.py --tutorial")
    print("   • Examples: python genomeamr.py --examples") 
    print("   • Setup check: python genomeamr.py --check-setup")
    print("   • Documentation: Read QUICKSTART.md and USE_CASES.md")

if __name__ == "__main__":
    main()
