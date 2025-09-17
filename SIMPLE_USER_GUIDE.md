# 🔬 GenomeAMRAnalyzer - Simple User Guide
*For Project Investigators & Non-Technical Users*

## 📋 What You Need (Just 3 Things!)

1. **Your Email Address** (for NCBI database access)
2. **Genes You Want to Study** (we provide common lists)
3. **Genomes to Analyze** (NCBI search URL OR list of genome IDs)

---

## 🚀 How to Run (Copy & Paste Commands)

### Option 1: Search NCBI for Genomes
```bash
python genomeamr_auto.py \
  --url "https://www.ncbi.nlm.nih.gov/nuccore/(Escherichia%20coli%20AND%20erythromycin%20resistance)" \
  --genes config/genes_erythromycin.txt \
  --email YOUR_EMAIL@institution.edu
```
*Replace `YOUR_EMAIL@institution.edu` with your actual email*

### Option 2: Use Your Own Genome List
```bash
python genomeamr_auto.py \
  --accessions my_genomes.txt \
  --genes config/genes_default.txt \
  --email YOUR_EMAIL@institution.edu
```

---

## 📁 Available Gene Sets (Pre-Made for You)

- `config/genes_default.txt` - Common resistance genes
- `config/genes_erythromycin.txt` - Erythromycin resistance
- `config/genes_betalactam.txt` - Beta-lactam resistance (if available)

---

## 🎯 Example Commands for Common Research

### Studying E. coli Antibiotic Resistance
```bash
python genomeamr_auto.py \
  --url "https://www.ncbi.nlm.nih.gov/nuccore/(Escherichia%20coli%20AND%20antibiotic%20resistance)" \
  --genes config/genes_default.txt \
  --email researcher@university.edu
```

### Studying Specific Bacterial Genomes
```bash
# First, create a file called "my_genomes.txt" with genome IDs (one per line):
# GCF_000005825.2
# GCF_000006945.2
# GCF_000008865.2

python genomeamr_auto.py \
  --accessions my_genomes.txt \
  --genes config/genes_default.txt \
  --email researcher@university.edu
```

---

## ⏱️ What Happens When You Run It?

1. **Automatic Setup** (2-5 minutes)
   - Downloads necessary tools
   - Sets up databases
   - *You don't need to do anything!*

2. **Analysis** (5-30 minutes depending on genome count)
   - Downloads genomes
   - Finds resistance genes
   - Analyzes mutations
   - Creates reports

3. **Results Ready!**
   - HTML report opens in your browser
   - All files saved in `reports/` folder

---

## 📊 Understanding Your Results

### Main Report: `reports/pipeline_report.html`
- **Summary**: How many genomes, genes found, mutations detected
- **Gene Co-occurrence**: Which resistance genes appear together
- **Detailed Tables**: All findings with genome names and coordinates

### Key Files Created:
- `reports/pipeline_report.html` - **Main results (open this first!)**
- `cooccurrence_results/` - Gene interaction analysis
- `alignments/` - Protein sequence comparisons
- `proteins/` - Extracted resistance proteins

---

## ❓ Troubleshooting

### "Command not found" or "Python not found"
- Contact your IT support to install Python
- Or use Anaconda/Miniconda (easier for beginners)

### "Permission denied" or "Access denied"
- Run as Administrator (Windows) or use `sudo` (Mac/Linux)
- Or ask IT support for help

### Tool installation fails
- The tool will still run in "simulation mode"
- Results will be generated with available tools
- Contact the developer (you!) for assistance

### Need different genes?
- Create a text file with gene names (one per line)
- Use: `--genes your_gene_file.txt`

---

## 🆘 Need Help?

1. **First**: Check if `reports/pipeline_report.html` was created
2. **Second**: Look for error messages in the terminal
3. **Third**: Contact the developer (the person who gave you this tool)

---

## 📞 Contact Information

**Developer**: [Your Name]  
**Email**: [Your Email]  
**Project**: GenomeAMRAnalyzer

---

*Last Updated: September 2025*