# Quick Start Guide for GenomeAMRAnalyzer

## 🚀 Get Started in 5 Minutes

### Step 1: One-Command Installation
```bash
# Install with pip (recommended)
pip install genomeamranalyzer

# Or clone and install
git clone https://github.com/yourusername/GenomeAMRAnalyzer.git
cd GenomeAMRAnalyzer
pip install -e .
```

### Step 2: Quick Test Run
```bash
# Test with example data (no setup required)
genomeamr --quick-test

# Your first analysis with real data
genomeamr --accessions GCF_000005825.2,GCF_000006945.2 --genes acrA,acrB,tolC
```

### Step 3: View Results
```bash
# Results automatically open in your browser
# Or manually open: reports/pipeline_report.html
```

## 🎯 Common Use Cases

### For Researchers: Analyze Your Genomes
```bash
# Upload your accession list
echo "GCF_000005825.2" > my_genomes.txt
echo "GCF_000006945.2" >> my_genomes.txt

# Run analysis
genomeamr --accession-file my_genomes.txt --output my_results/
```

### For Students: Learn AMR Analysis
```bash
# Educational mode with detailed explanations
genomeamr --tutorial --genes acrA,acrB --explain-steps
```

### For Lab Managers: Batch Processing
```bash
# Process multiple studies
genomeamr --batch-mode --config lab_settings.yaml
```

## 📋 No Configuration Required

GenomeAMRAnalyzer works out-of-the-box with smart defaults:
- ✅ Automatic genome download
- ✅ Built-in reference databases
- ✅ Standard resistance genes pre-configured
- ✅ Interactive HTML reports
- ✅ No external tools required

## 🆘 Need Help?

### Quick Help
```bash
genomeamr --help
genomeamr --examples
```

### Common Issues
- **No internet?** Use `--offline-mode`
- **Slow download?** Use `--max-genomes 5`
- **Different genes?** Use `--genes yourGene1,yourGene2`

### Get Support
- 📖 [Full Documentation](https://github.com/yourusername/GenomeAMRAnalyzer/wiki)
- 💬 [Ask Questions](https://github.com/yourusername/GenomeAMRAnalyzer/discussions)
- 🐛 [Report Issues](https://github.com/yourusername/GenomeAMRAnalyzer/issues)

---

**Ready to analyze antimicrobial resistance? Get started with one command!**