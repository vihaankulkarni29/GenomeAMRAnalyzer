# Production Pipeline Usage Guide
## Senior Bioinformatician-Level Genome Processing System

### Overview
This production pipeline system handles large-scale genome processing (100-1000 genomes) with bulletproof naming conventions, complete provenance tracking, and comprehensive integrity checking. The system follows a senior bioinformatician approach with:

- **Accession-centric naming** across all pipeline stages
- **Complete provenance tracking** from URL to proteins  
- **Batch integrity checking** with comprehensive manifests
- **Production-grade error handling** and recovery
- **Comprehensive reporting** for reproducibility

### Pipeline Stages

#### Stage 1: URL Discovery & Genome Download
- Discovers genomes from NCBI search URLs
- Downloads with accession-based naming: `{accession}_{organism}.fasta`
- Generates metadata manifests with checksums and provenance
- **Output**: `genomes/` directory with FASTA files and `genome_manifest.json`

#### Stage 2: CARD RGI Analysis  
- Runs CARD RGI on downloaded genomes
- Generates coordinate files: `{accession}_coordinates.csv`
- Maintains strict accession-based naming
- **Output**: `coordinates/` directory with CSV files and `coordinate_manifest.json`

#### Stage 3: Protein Extraction
- Extracts proteins using RGI coordinates
- Creates protein files: `{accession}_proteins.fasta`
- Maintains full provenance from coordinates to proteins
- **Output**: `proteins/` directory with FASTA files and `extraction_manifest.json`

#### Stage 4: Master Manifest Generation
- Generates comprehensive pipeline manifest
- Provides complete traceability across all stages
- Includes performance metrics and integrity checksums
- **Output**: `manifests/` directory with master manifest and reports

### Quick Start

#### Basic Usage
```bash
# Run complete pipeline for efflux pump genes
python src/production_pipeline_orchestrator.py \
  --source-url "https://www.ncbi.nlm.nih.gov/nuccore/?term=Pseudomonas+aeruginosa+complete+genome" \
  --target-genes mexA mexB oprM \
  --output-dir /data/pipeline_output \
  --email your.email@institution.edu \
  --api-key your_ncbi_api_key
```

#### Advanced Usage with Configuration
```bash
# Use configuration file for production settings
python src/production_pipeline_orchestrator.py \
  --source-url "https://www.ncbi.nlm.nih.gov/nuccore/?term=Escherichia+coli+complete+genome" \
  --target-genes mexA mexB mexC mexD mexE mexF oprM oprN \
  --output-dir /data/large_batch_output \
  --config config/production_pipeline_config.yaml \
  --max-concurrent 30 \
  --rgi-threads 16
```

### Output Structure

The pipeline creates a standardized output structure:

```
pipeline_output/
├── genomes/                          # Downloaded genome FASTA files
│   ├── GCF_000001405.40.fasta       # Accession-named genomes
│   ├── GCF_000005825.2_E_coli.fasta # With organism if available
│   └── genome_manifest.json         # Complete download manifest
├── coordinates/                      # CARD RGI coordinate files  
│   ├── GCF_000001405.40_coordinates.csv
│   ├── GCF_000005825.2_coordinates.csv
│   └── coordinate_manifest.json     # RGI analysis manifest
├── proteins/                        # Extracted protein files
│   ├── proteins/
│   │   ├── GCF_000001405.40_proteins.fasta
│   │   └── GCF_000005825.2_proteins.fasta
│   └── manifests/
│       └── extraction_manifest.json # Protein extraction manifest
├── manifests/                       # Master pipeline manifests
│   └── pipeline_12345_master_manifest.json
├── logs/                           # Comprehensive logging
│   └── pipeline_12345.log
└── reports/                        # Performance and accession reports
    ├── pipeline_12345_accession_report.json
    └── pipeline_12345_performance_report.json
```

### Naming Conventions

#### Accession-Centric Naming
All files follow strict accession-based naming:

- **Genomes**: `{accession}.fasta` or `{accession}_{organism}.fasta`
- **Coordinates**: `{accession}_coordinates.csv`  
- **Proteins**: `{accession}_proteins.fasta`
- **Manifests**: Linked by accession across all stages

#### FASTA Headers
Protein FASTA headers include complete provenance:
```
>{accession}|{gene_name}|{start}-{end}|{strand}|{contig}|{checksum}
```

### Manifest System

#### Master Manifest
Complete pipeline overview with:
- Processing statistics and success rates
- File paths to all component manifests  
- Integrity checksums and provenance tracking
- Performance metrics and timing data

#### Component Manifests
- **Genome Manifest**: Download metadata, checksums, organism info
- **Coordinate Manifest**: RGI results, gene coordinates, analysis provenance
- **Extraction Manifest**: Protein sequences, extraction metadata, quality metrics

### Quality Control

#### Integrity Checking
- File checksum validation across all stages
- Cross-validation between manifests
- Accession consistency checking
- Complete provenance chain validation

#### Error Handling
- Graceful failure handling with detailed error messages
- Partial processing capability (continue with successful accessions)
- Comprehensive failure reporting and recovery guidance
- Automatic retry for transient failures

### Performance Optimization

#### Large Batch Processing (1000+ genomes)
- Chunk-based processing to manage memory
- Parallel download and processing
- Progress monitoring and estimation
- Resource usage optimization

#### Configuration Tuning
- Adjust `max_concurrent_downloads` based on network capacity
- Tune `rgi_threads` to match available CPU cores
- Configure memory limits for large genome sets
- Enable compression for storage optimization

### Reproducibility

#### Complete Provenance
Every output file includes:
- Source URL and discovery timestamp
- Download metadata and checksums
- RGI version and analysis parameters
- Extraction timestamps and methods

#### Manifest-Based Reproduction
Master manifests enable complete pipeline reproduction:
- Exact source URLs and parameters
- Component versions and configurations
- Success/failure tracking per accession
- Performance metrics for optimization

### Troubleshooting

#### Common Issues
1. **NCBI API Rate Limiting**: Use API key and reduce concurrent downloads
2. **RGI Installation**: Ensure CARD database is properly installed
3. **Memory Issues**: Reduce chunk size or increase memory limits
4. **Disk Space**: Enable compression and cleanup temporary files

#### Monitoring Progress
- Monitor log files in real-time: `tail -f logs/pipeline_*.log`
- Check intermediate manifests for stage completion
- Review performance reports for optimization opportunities

### Integration with Existing Workflows

#### Snakemake Integration
The pipeline outputs can be directly integrated with Snakemake workflows:
```python
# Use master manifest to define downstream rules
configfile: "manifests/pipeline_12345_master_manifest.json"

rule downstream_analysis:
    input: config["extraction_manifest_path"]
    output: "analysis_results.txt"
    shell: "analyze_proteins.py {input}"
```

#### HPC Cluster Usage
Configure for HPC environments:
- Use appropriate resource allocation
- Configure cluster-specific temporary directories
- Adjust thread counts for node specifications
- Use shared storage for output directories

### Support and Documentation

For additional support:
- Review component documentation in `docs/` directory
- Check test suites for usage examples
- Examine configuration templates for advanced options
- Monitor GitHub issues for known limitations

This production pipeline represents a senior bioinformatician-level approach to large-scale genome processing with complete traceability and bulletproof naming conventions.