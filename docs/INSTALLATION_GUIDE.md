## 📋 System Requirements

### Minimum Requirements
- **Python**: 3.9+ (recommended: 3.11+)
- **RAM**: 4GB minimum, 8GB recommended
- **Storage**: 10GB free space for analysis data
- **OS**: Windows 10+, macOS 10.15+, Linux (Ubuntu 18.04+)

### Recommended for Production
- **CPU**: 4+ cores for parallel processing
- **RAM**: 16GB+ for large genome datasets
- **Storage**: SSD recommended, 50GB+ for extensive analyses
- **Network**: Stable internet for NCBI downloads and CARD updates

### External Dependencies (Optional)
- **EMBOSS**: For advanced alignment algorithms (auto-downloaded if needed)
- **CARD RGI**: Integrated resistance gene analysis (fallback mode available)

## 🔧 Installation Guide

### Option 1: Standard Installation
```bash
# Clone repository
git clone https://github.com/yourusername/GenomeAMRAnalyzer.git
cd GenomeAMRAnalyzer

# Create virtual environment (recommended)
python -m venv genomeamr-env
source genomeamr-env/bin/activate  # Linux/macOS
# genomeamr-env\Scripts\activate  # Windows

# Install dependencies
pip install -r requirements.txt

# Verify installation
python -c "import src.card_runner; print('Installation successful!')"
```

### Option 2: Package Installation
```bash
# Install from source as package
pip install -e .

# Use command-line tools
genomeamr --help
genomeamr-download --help
```

### Option 3: Development Installation
```bash
# Install with development tools
pip install -e ".[dev,visualization,performance]"

# Run tests
pytest tests/

# Code formatting
black src/
flake8 src/
```

## 🚀 Usage Examples

### Example 1: Quick Analysis
```bash
# Analyze pre-downloaded genomes
python run_pipeline.py \
    --config config/snakemake_config.yaml \
    --skip-download

# Results will be in: reports/html/
```

### Example 2: Full Pipeline with Custom Genes
```bash
# Create custom accession list
echo -e "GCF_000005825.2\nGCF_000006945.2" > my_accessions.txt

# Run full analysis
python run_pipeline.py \
    --config config/snakemake_config.yaml
```

### Example 3: Individual Module Usage
```bash
# Download genomes only
python src/simple_genome_downloader.py \
    --accession-list accessions.txt \
    --output-dir genomes/ \
    --batch-size 5

# Run CARD analysis only  
python src/card_runner.py \
    --input-dir genomes/ \
    --output-dir card_results/ \
    --genes acrA acrB tolC mexA mexB oprM

# Generate report only
python src/html_report_generator.py \
    --config config/snakemake_config.yaml \
    --aggregation-report results/manifest.json
```

### Example 4: High-Performance Computing
```bash
# Configure for HPC environment
export GENOMEAMR_THREADS=16
export GENOMEAMR_MEMORY=32

# Run with high-performance settings
python run_pipeline.py \
    --config config/hpc_config.yaml
```

## ⚙️ Configuration

### Basic Configuration (`config/snakemake_config.yaml`)
```yaml
# Environment settings
environment:
  deployment_type: "local"  # local, hpc, cloud
  max_parallel_jobs: 4
  memory_limit_gb: 8

# NCBI configuration  
ncbi:
  email: "your.email@institution.edu"  # REQUIRED
  api_key: ""  # Recommended for production
  max_genomes: 50

# Analysis parameters
analysis:
  target_genes:
    - "acrA"
    - "acrB"  
    - "tolC"
  quality_control:
    min_alignment_identity: 70.0
    min_coverage: 5
```

### Production Configuration Tips
1. **Email Configuration**: Always set a valid email for NCBI requests
2. **API Keys**: Use NCBI API keys for higher rate limits
3. **Resource Limits**: Adjust memory/CPU based on system capabilities
4. **Gene Selection**: Customize target genes for your research focus

## 📊 Output Structure

```
GenomeAMRAnalyzer/
├── genome_data/           # Downloaded genomes
├── card_results/          # CARD RGI analysis
├── proteins/              # Extracted protein sequences
├── alignments/            # Sequence alignments
├── mutations/             # Mutation analysis
├── cooccurrence/          # Co-occurrence analysis
├── reports/               # Final reports
│   ├── html/             # Interactive HTML reports
│   ├── summary.json      # Analysis summary
│   └── manifest.json     # Complete provenance
└── logs/                  # Execution logs
```

### Key Output Files
- `reports/html/*_report.html`: Interactive analysis reports
- `reports/summary.json`: High-level analysis summary
- `cooccurrence/*_pairs.csv`: Gene co-occurrence statistics
- `mutations/*_calls.csv`: Detailed mutation calls
- `logs/pipeline.log`: Complete execution log

## 🔬 Scientific Workflow

### Phase 1: Data Acquisition
1. **Genome Download**: Fetch bacterial genomes from NCBI using accession numbers
2. **Quality Control**: Validate genome completeness and size thresholds
3. **Metadata Extraction**: Parse organism information and assembly statistics

### Phase 2: Resistance Analysis
1. **CARD RGI**: Identify antimicrobial resistance genes using CARD database
2. **Coordinate Extraction**: Extract genomic coordinates for target genes
3. **Protein Extraction**: Generate amino acid sequences from coordinates

### Phase 3: Comparative Analysis
1. **Reference Alignment**: Align extracted proteins to reference sequences
2. **Mutation Detection**: Statistical mutation calling with confidence scoring
3. **Clinical Annotation**: Annotate mutations with resistance mechanisms

### Phase 4: Network Analysis
1. **Co-occurrence Statistics**: Calculate pairwise gene associations
2. **Network Construction**: Build resistance gene interaction networks
3. **Statistical Testing**: Fisher's exact tests and permutation analysis

### Phase 5: Reporting
1. **Data Aggregation**: Combine results across all analysis phases
2. **Visualization**: Generate interactive charts and network diagrams
3. **HTML Reports**: Create comprehensive analysis reports

## 🧬 SEPI Integration

GenomeAMRAnalyzer integrates with [SEPI 2.0](https://github.com/vihaankulkarni29/sepi2.0) for dynamic reference protein fetching:

### Features
- **Automatic Reference Fetching**: Download species-specific reference proteins
- **Cache Management**: Local storage of fetched references for reuse
- **Fallback Support**: Graceful degradation if SEPI unavailable
- **Quality Scoring**: Reference quality assessment and validation

### Usage
```python
# Automatic integration - no manual configuration needed
# SEPI fetches references during alignment phase
```

## 🚨 Troubleshooting

### Common Issues

#### Installation Problems
```bash
# Issue: Package conflicts
Solution: Use virtual environment
python -m venv fresh-env
source fresh-env/bin/activate
pip install -r requirements.txt

# Issue: Missing BioPython
pip install biopython>=1.80
```

#### Runtime Errors
```bash
# Issue: NCBI download failures
Solution: Check internet connection and email config
# Edit config/snakemake_config.yaml:
ncbi:
  email: "your.email@domain.com"  # Valid email required

# Issue: Memory errors
Solution: Reduce batch size or increase system memory
analysis:
  batch_size: 5  # Reduce from default 10
```

#### Output Issues
```bash
# Issue: Empty results
Solution: Check target gene names and genome quality
# Verify genes exist in genomes:
grep -r "acrA" card_results/

# Issue: Missing HTML report
Solution: Check all pipeline steps completed
ls reports/html/  # Should contain .html files
```

### Performance Optimization

#### For Large Datasets (>100 genomes)
1. Increase parallel processing:
   ```yaml
   environment:
     max_parallel_jobs: 8
     memory_limit_gb: 16
   ```

2. Use batch processing:
   ```yaml
   performance:
     batch_size: 20
     parallel_downloads: 5
   ```

#### For Limited Resources
1. Reduce memory usage:
   ```yaml
   performance:
     batch_size: 5
     memory_per_job: 1
   ```

2. Enable cleanup:
   ```yaml
   performance:
     cleanup_temp: true
     compress_outputs: true
   ```

## 📚 Advanced Usage

### Custom Gene Analysis
```yaml
# Add custom resistance genes to config
analysis:
  target_genes:
    - "custom_gene_1"
    - "custom_gene_2"
    - "efflux_pump_x"
```

### Integration with Existing Workflows
```python
# Use as Python library
from src.card_runner import CARDRunner
from src.production_cooccurrence_analyzer import ProductionCooccurrenceAnalyzer

# Initialize components
card = CARDRunner("config.yaml")
cooccur = ProductionCooccurrenceAnalyzer("config.yaml")

# Run analysis
results = card.analyze_genomes("genome_dir/")
report = cooccur.analyze_cooccurrence(["acrA", "acrB"])
```

### High-Performance Computing
```bash
# SLURM job script example
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00

module load python/3.11
source genomeamr-env/bin/activate

export GENOMEAMR_THREADS=16
export GENOMEAMR_MEMORY=64

python run_pipeline.py --config config/hpc_config.yaml
```

## 🤝 Contributing

We welcome contributions! Please see our [Contributing Guidelines](CONTRIBUTING.md) for details.

### Development Setup
```bash
# Clone and setup development environment
git clone https://github.com/yourusername/GenomeAMRAnalyzer.git
cd GenomeAMRAnalyzer
pip install -e ".[dev]"

# Run tests
pytest tests/ -v

# Code quality checks
black src/
flake8 src/
mypy src/
```

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 📞 Support

- **GitHub Issues**: [Report bugs and request features](https://github.com/yourusername/GenomeAMRAnalyzer/issues)
- **Discussions**: [Community support and questions](https://github.com/yourusername/GenomeAMRAnalyzer/discussions)
- **Documentation**: [Wiki and guides](https://github.com/yourusername/GenomeAMRAnalyzer/wiki)

## 🏆 Citation

If you use GenomeAMRAnalyzer in your research, please cite:

```bibtex
@software{genomeamranalyzer2025,
  title={GenomeAMRAnalyzer: Production-grade antimicrobial resistance gene analysis pipeline},
  author={GenomeAMRAnalyzer Development Team},
  version={2.0.0},
  year={2025},
  url={https://github.com/yourusername/GenomeAMRAnalyzer}
}
```

## 🔬 Scientific Background

### Antimicrobial Resistance (AMR)
Antimicrobial resistance represents one of the most pressing challenges in modern medicine. This pipeline focuses on:

- **Efflux Pump Systems**: RND family transporters (AcrAB-TolC, MexAB-OprM)
- **Resistance Mechanisms**: Gene mutations, overexpression, horizontal transfer
- **Clinical Relevance**: Correlation with treatment outcomes and epidemiology
- **Population Genomics**: Resistance evolution and transmission patterns

### Technical Innovation
- **Statistical Rigor**: Multiple testing correction, confidence intervals
- **Network Analysis**: Graph theory applications to resistance patterns
- **Clinical Integration**: CARD database integration for mechanism annotation
- **Reproducibility**: Complete provenance tracking and version control

---

**GenomeAMRAnalyzer v2.0.0** - Production-ready antimicrobial resistance analysis for the genomics era.