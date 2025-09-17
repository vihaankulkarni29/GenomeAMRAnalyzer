# PhD-Ready AMR Genomics Pipeline Setup Guide

## Overview
This pipeline is engineered for **production research use** by PhD students and researchers analyzing antimicrobial resistance (AMR) genomics data. It provides robust, error-tolerant processing of NCBI genome data with comprehensive mutation analysis.

## ✅ **What's Fixed for Production Use**

### **Robust Error Handling**
- ✅ **Per-genome processing**: One failed genome doesn't stop the entire pipeline
- ✅ **Graceful degradation**: Missing tools (like RGI) are detected and handled
- ✅ **Comprehensive logging**: Every error is logged with context and suggestions
- ✅ **Always generates reports**: Even partial results produce useful output

### **Flexible Accession Handling**
- ✅ **Smart conversion**: Automatically converts NZ_/NC_/CP_ → GCF_/GCA_ when possible
- ✅ **Fallback strategies**: If conversion fails, tries individual accession searches
- ✅ **Deduplication**: Removes duplicate accessions automatically
- ✅ **Validation**: Checks accession formats and warns about issues

### **Research-Grade Features**
- ✅ **Checkpointing**: Resume interrupted runs without re-downloading
- ✅ **Batch processing**: Handles large accession lists efficiently
- ✅ **MIC integration**: Links resistance phenotypes to genotypes
- ✅ **Co-occurrence analysis**: Identifies mutation patterns across proteins

## 🚀 **Installation (One-Time Setup)**

### **1. Install Core Dependencies**
```powershell
# Install Python packages
pip install -r requirements.txt

# Install RGI (CARD Resistance Gene Identifier) - REQUIRED for protein analysis
conda install -c bioconda rgi

# Download CARD database
rgi load --card_json
```

### **2. Set Up NCBI Access**
Edit `accession_pipeline.yaml`:
```yaml
ncbi:
  email: "your.email@university.edu"  # REQUIRED - use your real email
  api_key: ""  # Optional - get from NCBI for faster access
```

### **3. Prepare Your Data**
Create two files:
- `accession_list.txt` - One genome accession per line (any format: GCF_, GCA_, NZ_, NC_, CP_)
- `gene_list.txt` - One gene symbol per line (e.g., acrA, acrB, tolC)

## 📊 **Running the Pipeline**

### **Simple One-Click Run**
```powershell
.\run_realworld_pipeline.bat
```

### **What It Does Automatically**
1. **Validates** your input files and configuration
2. **Converts** accessions to proper assembly format when possible
3. **Downloads** genomes from NCBI with retry logic
4. **Analyzes** resistance genes using CARD/RGI
5. **Extracts** your target proteins
6. **Aligns** proteins to reference sequences
7. **Identifies** mutations and co-occurrence patterns
8. **Harvests** MIC data when available
9. **Generates** comprehensive HTML report

### **Output Files**
- `report/html_reports/amr_report_*.html` - Main research report
- `logs/pipeline.log` - Detailed execution log
- `logs/accession_conversion.log` - Conversion details
- `cooccurrence_results/` - Mutation co-occurrence data
- `genome_data/` - Downloaded genome files

## 🔬 **For PhD Research Use**

### **Handling Large Datasets**
- The pipeline is designed for **hundreds of genomes**
- Automatic **rate limiting** respects NCBI guidelines
- **Checkpointing** allows resuming multi-hour runs
- **Error isolation** means one bad genome doesn't stop analysis

### **Customizing for Your Research**
Edit `accession_pipeline.yaml`:
```yaml
input:
  accessions_file: "my_study_accessions.txt"
  gene_list_file: "efflux_pump_genes.txt"
  max_genomes: 100  # Safety limit for query-based searches

subscan:
  score_threshold: 0.7  # Adjust mutation calling sensitivity
  
cooccurrence:
  min_mutation_frequency: 0.05  # Minimum frequency for co-occurrence analysis
```

## 🛠 **Troubleshooting**

### **Common Issues and Solutions**

#### "RGI binary not found"
```powershell
conda install -c bioconda rgi
rgi load --card_json
```

#### "No genomes were downloaded"
- Check your internet connection
- Verify accessions are valid (recently submitted accessions may not be in assembly database yet)
- Check `logs/pipeline.log` for specific NCBI errors

#### "Conversion failed for all accessions"
- This is normal for very new accessions (< 6 months old)
- The pipeline will automatically try individual searches as fallback
- Check `accession_list_conversion_log.json` for details

#### "No mutations found"
- Ensure your gene list matches genes present in CARD database
- Check reference files in `reference_proteins/` directory
- Lower the `score_threshold` in config for more sensitive detection

### **Getting Help**
1. Check `logs/pipeline.log` for detailed error messages
2. Review the HTML report's "Pipeline Status" section
3. Check individual stage logs in respective output directories

## 📈 **Pipeline Performance**

### **Typical Processing Times**
- **10 genomes**: 15-30 minutes
- **50 genomes**: 1-2 hours  
- **100 genomes**: 2-4 hours

### **Resource Requirements**
- **RAM**: 4-8 GB recommended
- **Storage**: ~1 GB per 100 genomes
- **Network**: Stable internet for NCBI downloads

## 🎯 **Research Output Quality**

### **Data Integrity**
- ✅ All downloads are checksummed and validated
- ✅ Mutations are called against high-quality reference sequences
- ✅ MIC data is linked via verified BioSample connections
- ✅ Co-occurrence analysis uses statistical significance testing

### **Reproducibility**
- ✅ Complete execution logs with timestamps
- ✅ Version tracking of databases and tools
- ✅ Deterministic output (same inputs = same results)
- ✅ Full provenance tracking from accession to mutation

## 📋 **Validation Checklist**

Before publishing results, verify:
- [ ] All target genomes were successfully downloaded
- [ ] RGI analysis completed for >90% of genomes
- [ ] Target genes were found in >80% of genomes
- [ ] Reference alignments show good coverage (>70%)
- [ ] MIC data is available for sufficient subset
- [ ] Co-occurrence patterns are statistically significant

The pipeline generates detailed statistics for each of these metrics in the final report.

---

**This pipeline is production-ready for PhD research.** It handles real-world data challenges, provides comprehensive error reporting, and generates publication-quality results with full provenance tracking.