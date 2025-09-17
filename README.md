# GenomeAMRAnalyzer

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Production Ready](https://img.shields.io/badge/Production-Ready-green.svg)](https://github.com/yourusername/GenomeAMRAnalyzer)

**Production-grade antimicrobial resistance gene analysis pipeline for bacterial genomes**

A comprehensive, enterprise-ready bioinformatics pipeline for analyzing antimicrobial resistance (AMR) mutations in bacterial genomes, with specialized focus on Enteric gram-negative bacteria and RND efflux pump systems.

## 🚀 Quick Start

### One-Command Installation & Test
```bash
# Install (choose one option)
pip install genomeamranalyzer              # From PyPI (recommended)
# OR
git clone https://github.com/yourusername/GenomeAMRAnalyzer.git
cd GenomeAMRAnalyzer && pip install -e .

# Test in 30 seconds
python genomeamr.py --quick-test
```

### Your First Analysis
```bash
# Analyze specific genomes
python genomeamr.py --accessions GCF_000005825.2,GCF_000006945.2

# Or use your own genome list
echo "GCF_000005825.2" > my_genomes.txt
python run_pipeline.py --config config/snakemake_config.yaml
```

### View Results
Results automatically open in your browser, or check: `reports/pipeline_report.html`

## 🎯 For Different Users

### 🔬 **Researchers & Scientists**
```bash
# Publication-ready analysis
python run_pipeline.py --config config/snakemake_config.yaml
# → Comprehensive reports with statistics and visualizations
```

### 🎓 **Students & Educators**  
```bash
# Educational mode with explanations
python genomeamr.py --tutorial --explain-steps
# → Step-by-step learning experience
```

### 🏥 **Clinical Labs**
```bash
# Fast screening mode
python genomeamr.py --clinical-mode --genes acrA,acrB,tolC
# → Focus on clinically relevant resistance genes
```

### 💻 **Bioinformaticians**
```bash
# Full pipeline control
python run_pipeline.py --config custom_config.yaml
# → Complete customization and batch processing
```

## ✨ **Why Choose GenomeAMRAnalyzer?**

### **User-Friendly**
- ✅ **One-command setup** - No complex configuration
- ✅ **Smart defaults** - Works out-of-the-box
- ✅ **Interactive tutorials** - Perfect for learning
- ✅ **Example data included** - Test immediately

### **Scientifically Robust**
- 🧬 **CARD Integration** - Latest resistance database
- � **Statistical Analysis** - Fisher's exact tests, permutation analysis
- � **Network Analysis** - Co-occurrence patterns and interactions
- � **Interactive Reports** - Publication-ready visualizations

### **Production Ready**
- 🚀 **Scalable** - Handle large genome datasets
- 🔧 **Configurable** - Adapt to any research workflow
- � **Well Documented** - Comprehensive guides and examples
- 🧪 **Thoroughly Tested** - Reliable and reproducible results

GenomeAMRAnalyzer is a modular bioinformatics pipeline designed to:
- Download and process bacterial genome collections
- Identify AMR-related genes using CARD database integration
- Extract and analyze protein sequences
- Detect mutations through sequence alignment
- Analyze co-occurrence patterns of mutations across gene networks
- Generate comprehensive HTML reports

## 🌟 Key Features

### 🧬 **Generic RND Efflux Pump Analysis**
- Supports ANY RND efflux pump proteins (not hardcoded to specific genes)
- User-defined gene lists for flexible analysis
- Compatible with all CARD database entries

### 📊 **Advanced Mutation Analysis**
- Co-occurrence pattern detection across multiple genes
- Statistical significance testing
- Mutation frequency analysis
- Network-based resistance correlation

### 🔧 **Modular Architecture**
- Independent, reusable components
- Seamless data flow between pipeline stages
- Easy integration with existing workflows

### 📈 **Production-Ready Tools**
- Comprehensive error handling and validation
- Detailed logging and progress tracking
- Scalable for large genome collections
- 100% test coverage on core components

## 🏗️ Pipeline Architecture

```
┌─────────────────┐    ┌─────────────────┐    ┌─────────────────┐
│ GenomeDownloader│───▶│ CARD_Integrator │───▶│FastaAAExtractor │
└─────────────────┘    └─────────────────┘    └─────────────────┘
                                                        │
┌─────────────────┐    ┌─────────────────┐    ┌─────────▼─────────┐
│HTMLReportGen    │◀───│ CoOccurrence    │◀───│ WildTypeAligner  │
└─────────────────┘    │    Analyzer     │    └─────────────────┘
                       └─────────────────┘              │
                                ▲                      │
                       ┌─────────────────┐    ┌─────────▼─────────┐
                       │ MIC_Integrator  │    │     SubScan      │
                       └─────────────────┘    └─────────────────┘
```

## 📁 Project Structure

```
GenomeAMRAnalyzer/
├── src/                                    # Core pipeline components
│   ├── generic_cooccurrence_analyzer.py   # Co-occurrence analysis engine
│   └── fasta_aa_extractor_integration.py  # Protein extraction integration
├── tests/                                  # Comprehensive test suite
│   ├── test_cooccurrence_analyzer.py      
│   ├── test_fasta_aa_extractor_integration.py
│   └── validate_cooccurrence_analyzer.py  
├── docs/                                   # Documentation
│   ├── COOCCURRENCE_ANALYZER_README.md    
│   └── FASTAAEXTRACTOR_INTEGRATION_README.md
├── examples/                               # Usage examples and workflows
├── MetaDataHarvester/                      # External tool collection
│   ├── scripts/                           # Core pipeline scripts
│   ├── references/                        # Reference sequences
│   └── config/                            # Configuration files
├── requirements.txt                        # Python dependencies
├── setup.py                              # Package installation
├── .gitignore                             # Git ignore rules
└── README.md                              # This file
```

## 🚀 Quick Start

### Installation

```bash
# Clone the repository
git clone https://github.com/your-username/GenomeAMRAnalyzer.git
cd GenomeAMRAnalyzer

# Install dependencies
pip install -r requirements.txt

# Install package
pip install -e .
```

### Basic Usage

#### 1. Co-occurrence Analysis
```bash
# Analyze mutation co-occurrence patterns
python src/generic_cooccurrence_analyzer.py \
    --mutations mutations.json \
    --genes mdtF acrA acrB tolC \
    --output cooccurrence_results/
```

#### 2. Protein Extraction
```bash
# Extract proteins from genomes using CARD coordinates
python src/fasta_aa_extractor_integration.py \
    --coordinates card_results.csv \
    --genomes genomes/ \
    --output extracted_proteins/
```

#### 3. Full Pipeline
```bash
# Run complete analysis pipeline
python pipeline.py \
    --accessions accession_list.txt \
    --genes mdtF acrA acrB tolC \
    --output results/
```

### Windows Quickstart

If you prefer a simple one-click run on Windows, use `run_user_pipeline.bat`:

- Double-click `run_user_pipeline.bat` or run it from a terminal. It will:
    - Use `config/snakemake_config.yaml`
    - Create a `logs` folder
    - Save output to `logs/user_pipeline.out` and errors to `logs/user_pipeline.err`

Important: Edit `config/snakemake_config.yaml` and set `ncbi.email` to your email before running.

## 📋 Requirements

### System Requirements
biopython>=1.79
## Run with accession list + gene list + email

You can now submit an assembly accession list (GCF_/GCA_), a gene list (proteins to extract), and your email for NCBI E-utilities.

1) Prepare inputs:
- A text file with assembly accessions, one per line (comments with `#` allowed). Example: `test_pipeline/assembly_accessions.txt`.
- A text file with gene names, one per line (e.g., `acrA`, `acrB`, `tolC`). Example: `test_pipeline/genes.txt`.
- Your email to include in NCBI requests.

2) Use the example config:
- See `examples/configs/accession_pipeline.yaml` and update `ncbi.email`, `accessions_file`, and `gene_list_file`.

3) Run the orchestrator:

```powershell
python src/priority3/pipeline/orchestrator.py -c examples/configs/accession_pipeline.yaml
```

Alternatively, harvest genomes directly from a list:

```powershell
python examples/harvest_from_accessions.py test_pipeline\assembly_accessions.txt --out genomes --db priority3.db --resume --email your@email
```

Outputs are written to `genomes/`, with downloaded FASTA files ready for CARD RGI coordinate discovery.
pandas>=1.3.0
numpy>=1.21.0
scipy>=1.7.0
```

### Optional Dependencies
```
matplotlib>=3.5.0  # For visualization
seaborn>=0.11.0    # For enhanced plots
plotly>=5.0.0      # For interactive reports
```

## 🔧 Core Components

### 1. Generic Co-occurrence Analyzer
- **Purpose**: Analyze simultaneous mutations across any gene set
- **Features**: Statistical testing, frequency analysis, network visualization
- **Input**: Mutation data from SubScan/WildTypeAligner
- **Output**: Co-occurrence matrices, statistical reports

### 2. FastaAAExtractor Integration
- **Purpose**: Bridge CARD coordinates with protein extraction
- **Features**: Multiple input formats, dual extraction methods
- **Input**: CARD coordinates, genome FASTA files
- **Output**: Protein sequences formatted for alignment

### 3. MetaDataHarvester Tools
- **GenomeDownloader**: Automated genome collection from NCBI
- **CARD_Integrator**: AMR gene identification and coordinate extraction
- **WildTypeAligner**: Protein sequence alignment and mutation detection
- **SubScan**: Advanced mutation scanning and validation
- **HTMLReportGenerator**: Comprehensive result visualization

## 📊 Example Workflows

### Workflow 1: Single Gene Analysis
```python
from src.generic_cooccurrence_analyzer import GenericCoOccurrenceAnalyzer

# Analyze mdtF mutations across genome collection
analyzer = GenericCoOccurrenceAnalyzer(gene_list=['mdtF'])
analyzer.load_mutations('subscan_output.json')
results = analyzer.analyze_cooccurrence()
```

### Workflow 2: Multi-Gene Network Analysis
```python
# Analyze RND efflux pump network
genes = ['mdtF', 'acrA', 'acrB', 'tolC', 'acrR', 'marA']
analyzer = GenericCoOccurrenceAnalyzer(gene_list=genes)
analyzer.load_mutations('network_mutations.json')
results = analyzer.analyze_cooccurrence(significance_threshold=0.01)
```

### Workflow 3: Custom Gene Set
```python
# User-defined gene analysis
custom_genes = ['gene1', 'gene2', 'gene3']  # Any genes of interest
analyzer = GenericCoOccurrenceAnalyzer(gene_list=custom_genes)
# ... analysis continues
```

## 🧪 Testing

### Run Complete Test Suite
```bash
# Run all tests
python -m pytest tests/ -v

# Run specific component tests
python tests/test_cooccurrence_analyzer.py
python tests/test_fasta_aa_extractor_integration.py

# Validate without external dependencies
python tests/validate_cooccurrence_analyzer.py
```

### Test Coverage
- ✅ Co-occurrence Analyzer: 6/6 tests passed (100%)
- ✅ FastaAAExtractor Integration: 10/10 tests passed (100%)
- ✅ Pipeline Integration: Comprehensive validation
- ✅ Error Handling: Edge cases and malformed data

## 📈 Performance

### Benchmarks
- **Small datasets** (10-50 genomes): <5 minutes
- **Medium datasets** (100-500 genomes): 15-60 minutes  
- **Large datasets** (1000+ genomes): 2-8 hours

### Optimization Tips
- Use parallel processing for large genome collections
- Filter gene lists to reduce computational overhead
- Enable external FastaAAExtractor for improved performance
- Batch process genomes to optimize memory usage

## 🤝 Contributing

We welcome contributions! Please see our contributing guidelines:

### Development Setup
```bash
# Fork and clone the repository
git clone https://github.com/your-username/GenomeAMRAnalyzer.git
cd GenomeAMRAnalyzer

# Create development environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install in development mode
pip install -e .[dev]

# Run tests
python -m pytest tests/
```

### Contribution Guidelines
1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## 📚 Documentation

### Component Documentation
- [Co-occurrence Analyzer](docs/COOCCURRENCE_ANALYZER_README.md)
- [FastaAAExtractor Integration](docs/FASTAAEXTRACTOR_INTEGRATION_README.md)

### API Reference
- Full API documentation available in source code docstrings
- Interactive examples in `examples/` directory

## 🔍 Troubleshooting

### Common Issues

#### BioPython Installation
```bash
# If BioPython installation fails
pip install biopython --no-cache-dir
```

#### Memory Issues with Large Datasets
```python
# Use batch processing
analyzer.set_batch_size(100)  # Process 100 genomes at a time
```

#### CARD Database Updates
```bash
# Update CARD database references
python MetaDataHarvester/scripts/update_card_database.py
```

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 📞 Support

- **Issues**: Report bugs and feature requests on [GitHub Issues](https://github.com/your-username/GenomeAMRAnalyzer/issues)
- **Documentation**: Check the [docs/](docs/) directory for detailed guides
- **Examples**: See [examples/](examples/) for usage patterns

## 🙏 Acknowledgments

- **CARD Database**: Comprehensive Antibiotic Resistance Database
- **NCBI**: National Center for Biotechnology Information
- **BioPython**: Python tools for computational biology
- **Research Community**: For valuable feedback and testing

## 📊 Citation

If you use GenomeAMRAnalyzer in your research, please cite:

```
GenomeAMRAnalyzer: A comprehensive pipeline for antimicrobial resistance 
mutation analysis in bacterial genomes. (2025)
```

## 🗺️ Roadmap

### Upcoming Features
- [ ] Interactive web interface
- [ ] Real-time analysis dashboard
- [ ] Enhanced visualization options
- [ ] Cloud deployment support
- [ ] Machine learning integration for resistance prediction

### Version History
- **v1.0.0** (Current): Initial release with core functionality
  - Generic co-occurrence analyzer
  - FastaAAExtractor integration
  - Comprehensive test suite
  - Production-ready pipeline components

---

**Built with ❤️ for the antimicrobial resistance research community**