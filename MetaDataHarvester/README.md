# MetaDataHarvester

A **scientifically validated** bioinformatics pipeline for analyzing antimicrobial resistance (AMR) mutations in bacterial efflux pump proteins, specifically focusing on AcrA/AcrB systems in Gram-negative bacteria.

## 🚀 **NEW: SEPI 2.0 Integration**

**10x Efficiency Boost!** Now integrated with [SEPI 2.0](https://github.com/vihaankulkarni29/sepi2.0) for automated protein extraction:

- ✅ **Automated RND Protein Extraction** - No more manual curation
- ✅ **Direct Genome Processing** - From accessions to results
- ✅ **Batch Processing** - Handle hundreds of genomes simultaneously
- ✅ **Quality Assurance** - SEPI's advanced validation algorithms

**See [SEPI_INTEGRATION_README.md](SEPI_INTEGRATION_README.md) for complete details!**

## 🔬 Scientific Validation Features

- **✅ Reference Sequence Validation**: Ensures biological relevance and proper selection
- **✅ Taxonomic Verification**: Confirms species classification and orthology
- **✅ Multiple Sequence Alignment**: Statistically valid mutation detection with evolutionary context
- **✅ Statistical Significance Testing**: Proper statistical analysis of mutation frequencies
- **✅ CARD Database Integration**: Cross-references with known resistance mutations
- **✅ Quality Control Metrics**: Comprehensive validation and reproducibility features

## 🧬 Species-Specific Enhancements

- **✅ Genus-Specific References**: Eliminates false positives from phylogenetic divergence
- **✅ Enhanced Metadata Collection**: Taxonomic, clinical, and functional annotations
- **✅ Phylogenetic Reference Selection**: Intelligent fallback to closely related species
- **✅ Alignment Quality Assessment**: Comprehensive validation of alignment reliability
- **✅ Clinical Correlation**: MIC data and resistance phenotype integration

## 📋 Pipeline Components

1. **NCBIProteinHarvester** - Downloads protein sequences from NCBI
2. **WildTypeAligner** - Creates sequence alignments against reference proteins
3. **SubScan** - Analyzes alignments for amino acid substitutions
4. **HTMLReportGenerator** - Creates interactive HTML reports

## 🛠️ Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/MetaDataHarvester.git
cd MetaDataHarvester

# Install dependencies
pip install -r requirements.txt
```

## 📖 Usage

### 🚀 **SEPI-Enhanced Workflow (Highly Recommended)**

```bash
# 1. Create accession list
echo "CP000034.1" > my_genomes.txt
echo "AE005174.2" >> my_genomes.txt

# 2. Run complete pipeline with SEPI automation
python scripts/Scientific_AMR_Workflow.py \
    --use-sepi \
    --accession-list my_genomes.txt \
    --protein-family acrA \
    --reference-db references \
    --email your.email@example.com \
    --output-dir sepi_results

# 3. View interactive results
# Open sepi_results/scientific_report.html in your browser
```

### Traditional Scientifically Validated Workflow

```bash
# Complete scientifically validated analysis (manual protein curation)
python scripts/Scientific_AMR_Workflow.py \
    --protein-fasta-dir protein_data/fasta \
    --metadata-file protein_data/metadata.csv \
    --protein-family acrA \
    --reference-db references \
    --email your.email@example.com
```

### Species-Specific Enhanced Workflow (Recommended)

```bash
# 1. Build species-specific reference database
python scripts/setup_species_references.py \
    --email your.email@example.com \
    --output-dir references

# 2. Enhanced metadata collection with taxonomic data
python scripts/data/enhanced_ncbi_harvester.py \
    --accession-file data/sequence.txt \
    --email your.email@example.com \
    --output-dir enhanced_protein_data

# 3. Species-specific alignment (eliminates false positives!)
python scripts/core/species_specific_aligner.py \
    --fasta-file enhanced_protein_data/sequences.fasta \
    --metadata-file enhanced_protein_data/enhanced_metadata.csv \
    --reference-db references
```

### Individual Components (Advanced Users)

```bash
# Scientifically validated components
python scripts/core/reference_validator.py --reference-dir references --email your.email@example.com
python scripts/core/taxonomic_validator.py --metadata data.csv --expected-species "E. coli" --protein-family acrA --email your.email@example.com
python scripts/core/msa_analyzer.py --fasta-dir sequences/ --reference ref.fasta --protein-family acrA
python scripts/core/card_integrator.py --csv-file mutations.csv --gene-name acrA

# Legacy components (for backward compatibility)
python scripts/data/ncbi_protein_harvester.py --accession-file data/sequence.txt --output-dir protein_data --email your.email@example.com
```

### Workflow Script

For automated execution of the complete pipeline:

```bash
python scripts/AcrAB_Mutation_Workflow.py --accession-file data/sequence.txt --email your.email@example.com
```

## 📁 Project Structure

```
MetaDataHarvester/
├── scripts/                 # Pipeline scripts
│   ├── core/               # Core analysis modules
│   ├── data/               # Data acquisition tools
│   └── AcrAB_Mutation_Workflow.py
├── config/                 # Configuration files
├── references/             # Reference protein sequences
├── data/                   # Input data files
├── docs/                   # Documentation
├── USAGE_GUIDE.md         # Detailed usage instructions
└── README.md              # This file
```

## 🔬 Research Applications

- **Antimicrobial Resistance Studies**: Identify mutations in efflux pump proteins
- **Comparative Genomics**: Analyze protein conservation across bacterial strains
- **Drug Resistance Mechanisms**: Study AcrA/AcrB system variations
- **Clinical Isolate Analysis**: Compare resistant vs. susceptible strains

## 📊 Output Formats

- **FASTA files**: Downloaded protein sequences
- **CSV files**: Mutation analysis results
- **NEEDLE files**: Sequence alignments
- **HTML reports**: Interactive dashboards with charts

## 🤝 Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Submit a pull request

## 📄 License

This project is licensed under the MIT License - see the LICENSE file for details.

## 📧 Contact

For questions or support, please open an issue on GitHub.

## 🔗 Dependencies

- Python 3.7+
- BioPython
- Pandas
- NumPy
- Matplotlib
- Chart.js (for HTML reports)

## 📚 Documentation

See `USAGE_GUIDE.md` for detailed usage instructions and examples.