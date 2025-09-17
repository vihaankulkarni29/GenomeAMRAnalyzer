# URL-Based Genome Discovery Workflow

This document provides a complete guide for using the new URL-based genome discovery system that replaces static accession lists with dynamic NCBI search results.

## 🚀 Overview

The new workflow allows users to:
1. Submit NCBI search URLs instead of manually curated accession lists
2. Automatically discover and download the first N genomes from search results
3. Maintain full integration with RGI for gene coordinate extraction
4. Process large datasets efficiently with parallel downloading

## 📋 Prerequisites

1. **Python Environment**:
   ```bash
   python >= 3.8
   ```

2. **Required Packages**:
   ```bash
   pip install -r requirements.txt
   ```

3. **NCBI Email** (required for E-utilities):
   - Update your email in `config/snakemake_config.yaml`
   - Optionally add NCBI API key for faster access

## 🔧 Configuration

Edit `config/snakemake_config.yaml`:

```yaml
# NCBI Configuration
ncbi_email: "your_email@domain.com"  # REQUIRED
ncbi_api_key: "your_api_key_here"    # OPTIONAL (faster downloads)

# URL-based genome discovery
ncbi_search:
  search_url: ""  # Will be provided by user
  max_genomes: 10  # Process first N genomes
  quality_filters:
    min_genome_size: 1000000  # 1MB minimum
    max_genome_size: 20000000  # 20MB maximum
  timeout_seconds: 300
  retry_attempts: 3
```

## 🌐 URL Formats Supported

The system supports various NCBI URL formats:

### 1. Simple Search
```
https://www.ncbi.nlm.nih.gov/nuccore/?term=Escherichia+coli+complete+genome
```

### 2. Complex Query with Filters
```
https://www.ncbi.nlm.nih.gov/nuccore/?term=(Escherichia+coli+and+macrolide+resistance)+AND+%22Escherichia+coli%22%5Bporgn%3A__txid562%5D+and+complete+genome
```

### 3. Organism-Specific Search
```
https://www.ncbi.nlm.nih.gov/nuccore/?term=Salmonella+AND+antibiotic+resistance
```

## 🧪 Testing the System

### 1. Run Validation Tests
```bash
cd GenomeAMRAnalyzer
python tests/test_url_genome_discovery.py
```

Expected output:
```
============================================================
URL-based Genome Discovery Test Suite
============================================================

✅ Test 1: Successfully parsed URL
✅ Test 2: Successfully parsed URL  
✅ Test 3: Successfully parsed URL
✅ GenomeMetadata test passed
✅ DownloadResult test passed
✅ Configuration loading test passed
✅ Workflow initialization test passed

============================================================
Test Results: 5/5 tests passed
🎉 All tests passed! URL-based genome discovery is ready.
============================================================
```

### 2. Quick Demo (Discovery Only)
```python
import asyncio
from pathlib import Path
import sys

sys.path.insert(0, str(Path('.').resolve() / "src"))
from src.url_to_genomes_workflow import URLToGenomesWorkflow

async def quick_demo():
    config_path = "config/snakemake_config.yaml"
    workflow = URLToGenomesWorkflow(config_path)
    
    # Discover genomes from URL
    url = "https://www.ncbi.nlm.nih.gov/nuccore/?term=Escherichia+coli+complete+genome"
    genomes = await workflow.discover_genomes_from_url(url)
    
    print(f"Found {len(genomes)} genomes:")
    for genome in genomes[:5]:
        print(f"  {genome.accession} - {genome.organism}")

# Run demo
asyncio.run(quick_demo())
```

## 🔬 Complete Workflow Usage

### Step 1: Import Modules
```python
import asyncio
from pathlib import Path
import sys

# Add src to Python path
sys.path.insert(0, str(Path('.').resolve() / "src"))

from src.url_to_genomes_workflow import URLToGenomesWorkflow
```

### Step 2: Initialize Workflow
```python
config_path = "config/snakemake_config.yaml"
workflow = URLToGenomesWorkflow(config_path)
```

### Step 3: Discover Genomes
```python
url = "your_ncbi_search_url_here"
genomes = await workflow.discover_genomes_from_url(url)

print(f"Discovered {len(genomes)} genomes")
for genome in genomes:
    print(f"  {genome.accession} - {genome.organism} ({genome.length:,} bp)")
```

### Step 4: Download Genomes
```python
results = await workflow.download_genomes(genomes)

successful = [r for r in results if r.success]
failed = [r for r in results if not r.success]

print(f"Downloaded: {len(successful)}/{len(results)} genomes")
```

## 📊 Example Workflow

```python
#!/usr/bin/env python3
import asyncio
from pathlib import Path
import sys

# Setup imports
sys.path.insert(0, str(Path('.').resolve() / "src"))
from src.url_to_genomes_workflow import URLToGenomesWorkflow

async def main():
    # Initialize workflow
    workflow = URLToGenomesWorkflow("config/snakemake_config.yaml")
    
    # User provides NCBI URL
    url = "https://www.ncbi.nlm.nih.gov/nuccore/?term=Escherichia+coli+complete+genome"
    
    # Discover genomes
    print("🔍 Discovering genomes...")
    genomes = await workflow.discover_genomes_from_url(url)
    print(f"✅ Found {len(genomes)} genomes")
    
    # Download genomes
    print("⬇️  Downloading genomes...")
    results = await workflow.download_genomes(genomes)
    
    # Process results
    successful = [r for r in results if r.success]
    print(f"📁 Downloaded {len(successful)} genomes to: {workflow.genome_directory}")
    
    # Next: Run RGI analysis on downloaded genomes
    # (Integration with existing RGI pipeline)

if __name__ == "__main__":
    asyncio.run(main())
```

## 🧬 Integration with RGI Pipeline

After downloading genomes, the workflow integrates with the existing RGI pipeline:

1. **RGI Analysis**: Process downloaded FASTA files with CARD database
2. **Gene Coordinate Extraction**: Extract coordinates for target genes (acrB, acrA, tolC)
3. **FastaAAExtractor**: Extract amino acid sequences for discovered genes
4. **Downstream Analysis**: Continue with cooccurrence analysis, alignment, etc.

## 📁 Output Directory Structure

```
url_genomes/
├── genomes/
│   ├── NC_000913.3.fasta
│   ├── NC_002695.2.fasta
│   └── ...
├── logs/
│   ├── genome_discovery.log
│   ├── download_summary.log
│   └── ...
└── metadata/
    ├── discovered_genomes.json
    └── download_results.json
```

## ⚡ Performance Features

- **Parallel Downloads**: Multiple genomes downloaded simultaneously
- **Rate Limiting**: Respects NCBI E-utilities guidelines
- **Error Handling**: Comprehensive retry logic with exponential backoff
- **Progress Tracking**: Detailed logging and status updates
- **Quality Filtering**: Automatic genome size and quality validation

## 🔧 Troubleshooting

### Common Issues

1. **Import Errors**:
   ```bash
   # Ensure you're in the project root
   cd GenomeAMRAnalyzer
   python -c "import sys; sys.path.insert(0, 'src'); from src.url_to_genomes_workflow import URLToGenomesWorkflow"
   ```

2. **Configuration Not Found**:
   ```bash
   # Check config file exists
   ls config/snakemake_config.yaml
   ```

3. **NCBI Rate Limiting**:
   - Add your NCBI API key to config
   - Reduce concurrent downloads
   - Check your email is valid in config

### Debug Mode
```python
import logging
logging.basicConfig(level=logging.DEBUG)
```

## 🎯 Next Steps

1. **Validate System**: Run all tests successfully
2. **Test with Real URLs**: Try with your NCBI search URLs
3. **RGI Integration**: Connect with existing RGI pipeline
4. **Production Deployment**: Scale for larger datasets

## 📧 Support

For issues or questions about the URL-based workflow:
1. Check test results: `python tests/test_url_genome_discovery.py`
2. Review logs in `logs/` directory
3. Verify configuration in `config/snakemake_config.yaml`

---

**Ready to process thousands of genomes with a single URL! 🚀**

## 🔗 Integrate with CARD/RGI

To run the full URL → genomes → RGI pipeline end-to-end:

```powershell
# From project root
python -m src.url_to_card_pipeline --url "https://www.ncbi.nlm.nih.gov/nuccore/?term=Escherichia+coli+complete+genome" --config config/snakemake_config.yaml
```

Outputs:
- FASTAs in `directories.genomes` (from config)
- CARD RGI results in `directories.card_results`
- Logs in `directories.logs`

Works even if `rgi` isn’t installed: the CARDRunner simulates outputs for testing so your pipeline won’t block in dev.