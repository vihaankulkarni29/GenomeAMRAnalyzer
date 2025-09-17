# FastaAAExtractor Pipeline Integration

## Overview

The FastaAAExtractor Integration is a comprehensive pipeline component that bridges CARD coordinates with protein sequence extraction for downstream analysis. This tool seamlessly connects genome coordinate data with amino acid sequence extraction, preparing data for the WildTypeAligner component.

## Features

### Core Functionality
- **Multiple Input Formats**: Supports CSV, TSV, and JSON coordinate files from CARD integrator
- **Dual Extraction Methods**: Internal BioPython-based extraction or external FastaAAExtractor tool integration
- **Flexible Gene Selection**: Extract all genes or filter to specific gene lists
- **WildTypeAligner Preparation**: Automatically formats output for downstream alignment analysis
- **Comprehensive Error Handling**: Robust validation and graceful error recovery

### Key Components

#### 1. GeneCoordinate Class
```python
@dataclass
class GeneCoordinate:
    gene_name: str
    genome_id: str
    start: int
    end: int
    strand: str
    contig: Optional[str] = None
    protein_id: Optional[str] = None
    confidence: float = 1.0
```

#### 2. ExtractionResult Class
```python
@dataclass
class ExtractionResult:
    gene_name: str
    genome_id: str
    sequence: str
    coordinates: GeneCoordinate
    extraction_method: str = "coordinate_based"
    quality_score: float = 1.0
    warnings: List[str] = field(default_factory=list)
```

#### 3. FastaAAExtractorIntegrator Class
Main integration class handling:
- Coordinate loading and validation
- Genome file discovery
- Protein extraction workflows
- Output file generation
- WildTypeAligner preparation

## Installation

### Requirements
- Python 3.8+
- BioPython (optional, for internal extraction method)
- Standard library modules: os, sys, logging, subprocess, tempfile, etc.

### Setup
```bash
# Copy integration files to your workspace
cp fasta_aa_extractor_integration.py /your/pipeline/directory/
cp test_fasta_aa_extractor_integration.py /your/pipeline/directory/

# Install optional dependencies
pip install biopython  # For internal extraction method
```

## Usage

### Command Line Interface

#### Basic Usage
```bash
# Extract proteins using CARD coordinates
python fasta_aa_extractor_integration.py \
    --coordinates card_results.csv \
    --genomes genomes/ \
    --output extracted_proteins/
```

#### Advanced Usage
```bash
# Extract specific genes with external tool
python fasta_aa_extractor_integration.py \
    --coordinates card_results.json \
    --genomes genomes/ \
    --genes mdtF acrA acrB tolC \
    --external-tool /path/to/FastaAAExtractor.py \
    --references references/ \
    --prepare-aligner \
    --output extracted_proteins/
```

### Python API

#### Basic Integration
```python
from fasta_aa_extractor_integration import FastaAAExtractorIntegrator

# Initialize integrator
integrator = FastaAAExtractorIntegrator(
    output_dir="extracted_proteins",
    log_level="INFO"
)

# Load CARD coordinates
integrator.load_card_coordinates("card_results.csv")

# Extract proteins
output_files = integrator.extract_proteins(
    genomes_dir="genomes/",
    gene_list=["mdtF", "acrA", "acrB"],
    use_external_tool=False
)

# Prepare for WildTypeAligner
aligner_files = integrator.prepare_for_wild_type_aligner(
    reference_dir="references/"
)
```

## Input Formats

### CARD Coordinates (CSV/TSV)
```csv
genome_id,gene_name,start,end,strand,contig
GCF_000005825.2,mdtF,1245,2678,+,chromosome
GCF_000005825.2,acrA,2890,3825,+,chromosome
GCF_000005825.2,acrB,3840,6951,+,chromosome
```

### CARD Coordinates (JSON)
```json
{
  "GCF_000005825.2": {
    "mdtF": {
      "start": 1245,
      "end": 2678,
      "strand": "+",
      "contig": "chromosome",
      "confidence": 0.95
    },
    "acrA": {
      "start": 2890,
      "end": 3825,
      "strand": "+",
      "contig": "chromosome",
      "confidence": 0.98
    }
  }
}
```

### Genome Files
- Standard FASTA format (.fna, .fasta, .fa)
- Supports both single-contig and multi-contig genomes
- Automatic file discovery based on genome IDs

### Reference Sequences
- Amino acid FASTA files (.faa)
- One file per gene (e.g., `mdtF_reference.faa`)
- Used for WildTypeAligner preparation

## Output Files

### Generated Outputs

1. **Individual Gene FASTA Files**
   ```
   mdtF_extracted_proteins.faa
   acrA_extracted_proteins.faa
   acrB_extracted_proteins.faa
   ```

2. **Combined FASTA File**
   ```
   all_extracted_proteins.faa
   ```

3. **Metadata File**
   ```json
   {
     "extraction_stats": {
       "total_genomes": 150,
       "successful_extractions": 148,
       "failed_extractions": 2,
       "sequences_generated": 444
     },
     "extraction_results": [...]
   }
   ```

4. **Summary Report**
   ```
   FastaAAExtractor Integration - Extraction Summary
   ==================================================
   Extraction Date: 2024-01-15T10:30:00
   Total Genomes: 150
   Successful Extractions: 148
   Failed Extractions: 2
   Sequences Generated: 444
   ```

### WildTypeAligner Preparation
```
wild_type_aligner_input/
├── mdtF/
│   ├── mdtF_query_sequences.faa
│   ├── mdtF_reference.faa
│   └── output/
├── acrA/
│   ├── acrA_query_sequences.faa
│   ├── acrA_reference.faa
│   └── output/
└── acrB/
    ├── acrB_query_sequences.faa
    ├── acrB_reference.faa
    └── output/
```

## Pipeline Integration

### Data Flow
```
CARD Integrator → FastaAAExtractor Integration → WildTypeAligner
     ↓                        ↓                        ↓
Coordinates         Protein Sequences        Alignments & Mutations
```

### Integration Points

#### From CARD Integrator
- Accepts coordinate files in multiple formats
- Validates coordinate data integrity
- Maps coordinates to genome files

#### To WildTypeAligner
- Generates properly formatted FASTA files
- Organizes sequences by gene for parallel processing
- Provides reference sequences for alignment

#### To Co-occurrence Analyzer
- Metadata files track extraction quality
- Failed extractions logged for downstream filtering

## Quality Control

### Validation Features
- Coordinate range validation
- Strand direction verification
- Sequence length checks
- Internal stop codon detection
- Reference sequence availability

### Error Handling
- Graceful handling of malformed coordinates
- Missing genome file recovery
- Partial extraction success reporting
- Comprehensive logging system

## Testing

### Comprehensive Test Suite
Run the complete test suite:
```bash
python test_fasta_aa_extractor_integration.py
```

### Test Coverage
- ✅ Gene coordinate validation
- ✅ Integrator initialization
- ✅ CSV/TSV/JSON coordinate loading
- ✅ Genome file discovery
- ✅ Gene filtering functionality
- ✅ Internal extraction simulation
- ✅ Reference sequence discovery
- ✅ WildTypeAligner preparation
- ✅ Error handling and edge cases

### Test Results
```
Total Tests: 10
Passed: 10
Failed: 0
Errors: 0
Success Rate: 100.0%
```

## Performance Considerations

### Scalability
- Processes genomes in batches
- Memory-efficient sequence handling
- Parallel processing support (via external tools)
- Temporary file cleanup

### Optimization Tips
- Use external FastaAAExtractor for large datasets
- Filter gene lists to reduce processing time
- Batch genome processing for memory efficiency
- Enable parallel processing where available

## Troubleshooting

### Common Issues

#### Missing BioPython
```
Warning: BioPython not available. Some features may be limited.
Solution: pip install biopython
```

#### Genome File Not Found
```
Warning: Genome file not found for GCF_123456789.1
Solution: Check genome file naming conventions and extensions
```

#### Coordinate Validation Errors
```
Error: Invalid coordinates: start=400, end=100
Solution: Verify coordinate data integrity in CARD output
```

#### External Tool Execution
```
Error: External tool execution failed
Solution: Check external tool path and permissions
```

### Debugging
- Enable DEBUG log level for detailed output
- Check extraction summary for failed extractions
- Validate input coordinate formats
- Verify genome file accessibility

## Advanced Configuration

### Custom Extraction Methods
```python
# Custom external tool integration
integrator = FastaAAExtractorIntegrator(
    external_extractor_path="/custom/path/extractor.py",
    temp_dir="/custom/temp",
    log_level="DEBUG"
)
```

### Batch Processing
```python
# Process multiple coordinate files
for coord_file in coordinate_files:
    integrator.load_card_coordinates(coord_file)
    
output_files = integrator.extract_proteins(
    genomes_dir="genomes/",
    use_external_tool=True
)
```

### Quality Filtering
```python
# Filter results by quality score
high_quality_results = [
    result for result in integrator.extraction_results 
    if result.quality_score > 0.8 and not result.warnings
]
```

## Integration with Pipeline Components

### GenomeDownloader Output
- Automatically discovers genome files from GenomeDownloader output
- Supports standard NCBI file naming conventions
- Handles both individual and batch downloads

### CARD Integrator Output
- Direct compatibility with CARD integrator coordinate formats
- Preserves gene confidence scores and metadata
- Maintains coordinate validation integrity

### WildTypeAligner Input
- Generates EMBOSS WATER compatible FASTA files
- Organizes sequences for parallel alignment processing
- Provides reference sequences for accurate alignment

### SubScan Integration
- Extraction metadata supports mutation quality filtering
- Failed extractions excluded from downstream analysis
- Quality scores propagated to mutation analysis

## Version History

### v1.0 (Current)
- Initial production release
- Full CARD coordinate integration
- Dual extraction method support
- WildTypeAligner preparation
- Comprehensive test suite (100% pass rate)

## Support

### Documentation
- Complete API documentation in source code
- Comprehensive test examples
- Integration guides for pipeline components

### Contact
- Part of GenomeAMRAnalyzer Pipeline project
- For issues and feature requests, consult pipeline documentation

## License

This component is part of the GenomeAMRAnalyzer Pipeline project. See main project documentation for licensing information.