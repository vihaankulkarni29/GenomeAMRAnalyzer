# Single Sequence Analyzer Usage Guide

## Quick Start for Your Workflow

### Step 1: Prepare Your Input Files
You need two files for each bacterial genome:
1. **Genome FASTA file** (`.fna`) - contains the DNA sequence
2. **GFF annotation file** (`.gff`) - contains gene annotations

### Step 2: Run the Analysis
```bash
# Basic usage
python single_sequence_analyzer.py \
    --genome-file your_bacteria.fna \
    --gff-file your_bacteria.gff

# With custom output directory
python single_sequence_analyzer.py \
    --genome-file your_bacteria.fna \
    --gff-file your_bacteria.gff \
    --output-dir my_analysis_results
```

### Step 3: Check Your Results
The tool will create a directory with:
- `reports/single_sequence_analysis.html` - Open this in your browser
- `substitutions/all_substitutions.csv` - All mutations found
- `reports/analysis_summary.json` - Detailed analysis data

## What the Tool Does (Exactly Your Workflow)

### 1. Protein Extraction
- Extracts **AcrA** protein sequence from your genome
- Extracts **AcrB** protein sequence from your genome
- Uses the GFF file to find the correct gene locations

### 2. Reference Alignment
- Aligns your AcrA against the reference AcrA sequence
- Aligns your AcrB against the reference AcrB sequence
- Uses EMBOSS needle for high-quality global alignments

### 3. Mutation Detection
- Identifies all amino acid substitutions
- Specifically tracks your target mutations:
  - **T104A** in AcrA (Alanine replaces Threonine at position 104)
  - **H596N** in AcrB (Asparagine replaces Histidine at position 596)

### 4. Analysis & Reporting
- Counts how many times each target mutation occurs
- Checks if both mutations occur together in the same bacteria
- If target mutations not found, shows what mutations ARE present
- Generates HTML report with all findings

## Example Output

### Console Output:
```
2025-09-03 08:18:36 - INFO - Starting single sequence analysis
2025-09-03 08:18:36 - INFO - Step 1: Extracting AcrA and AcrB proteins
2025-09-03 08:18:36 - INFO - Extracted AcrA protein to results/proteins/bacteria_acrA.fasta
2025-09-03 08:18:36 - INFO - Extracted AcrB protein to results/proteins/bacteria_acrB.fasta
2025-09-03 08:18:36 - INFO - Step 2: Aligning proteins against references
2025-09-03 08:18:36 - INFO - AcrA alignment completed
2025-09-03 08:18:36 - INFO - AcrB alignment completed
2025-09-03 08:18:36 - INFO - Step 3: Identifying amino acid substitutions
2025-09-03 08:18:36 - INFO - Found 3 mutations in bacteria_acrA.needle
2025-09-03 08:18:36 - INFO - Found 5 mutations in bacteria_acrB.needle
2025-09-03 08:18:36 - INFO - TARGET MUTATION ANALYSIS:
2025-09-03 08:18:36 - INFO -   T104A in AcrA: 1 times
2025-09-03 08:18:36 - INFO -   H596N in AcrB: 1 times
2025-09-03 08:18:36 - INFO -   Co-occurring: YES
```

### HTML Report Shows:
- ✅ **T104A found in AcrA** (1 time at position 104)
- ✅ **H596N found in AcrB** (1 time at position 596)
- ✅ **Co-occurring: YES** - Both mutations present together
- 📊 All other mutations found in the sequence
- 📈 Statistics and recommendations

## Scaling to Multiple Sequences

### Option 1: Batch Processing Script
```bash
#!/bin/bash
# batch_analyze.sh

GENOMES_DIR="path/to/your/genomes"
OUTPUT_BASE="batch_results"

for genome_file in $GENOMES_DIR/*.fna; do
    # Get base name without extension
    base_name=$(basename "$genome_file" .fna)
    gff_file="$GENOMES_DIR/${base_name}.gff"
    output_dir="$OUTPUT_BASE/${base_name}_analysis"

    echo "Analyzing $base_name..."

    python single_sequence_analyzer.py \
        --genome-file "$genome_file" \
        --gff-file "$gff_file" \
        --output-dir "$output_dir"

    echo "Completed $base_name"
done
```

### Option 2: Python Batch Script
```python
# batch_analyzer.py
import os
import glob
from single_sequence_analyzer import SingleSequenceAnalyzer

genomes_dir = "path/to/your/genomes"
output_base = "batch_results"

# Find all genome files
genome_files = glob.glob(f"{genomes_dir}/*.fna")

for genome_file in genome_files:
    base_name = os.path.basename(genome_file).replace('.fna', '')
    gff_file = f"{genomes_dir}/{base_name}.gff"
    output_dir = f"{output_base}/{base_name}_analysis"

    print(f"Analyzing {base_name}...")

    analyzer = SingleSequenceAnalyzer(genome_file, gff_file, output_dir)
    results = analyzer.analyze_sequence()

    if results['success']:
        print(f"✓ Completed {base_name}")
    else:
        print(f"✗ Failed {base_name}: {results.get('error', 'Unknown error')}")
```

## Data Collection for ML Training

### What You Get for Each Sequence:
1. **Target Mutations Present/Absent** - Binary classification data
2. **Co-occurrence Patterns** - Which mutations appear together
3. **All Mutations Found** - Complete mutation profile
4. **Protein Sequences** - For sequence-based ML models
5. **Alignment Data** - Quality metrics and statistics

### Building Your Dataset:
After analyzing multiple sequences, you'll have:
- **Features**: Mutation presence/absence, positions, types
- **Labels**: Antimicrobial resistance profiles
- **Metadata**: Species, collection date, source

### Example Dataset Structure:
```csv
Isolate_ID,Species,AcrA_T104A,AcrB_H596N,Co_occurring,Resistance_Profile,Other_Mutations,...
ISO001,E_coli,1,1,1,CRE,T104A,H596N,L4S_acrB,...
ISO002,K_pneumoniae,0,1,0,ESBL,H596N_acrB,V8A_acrA,...
ISO003,E_coli,1,0,0,AMP,T104A_acrA,G12D_acrA,...
```

## Troubleshooting

### Common Issues:

#### 1. "No AcrA/AcrB proteins found"
- Check if genes are properly annotated in GFF file
- Genes might use different names (acrA vs AcrA)
- Genome might not contain these genes

#### 2. "Alignment failed"
- Check reference sequences exist in `references/` directory
- Verify protein sequences are valid amino acids
- Check for sequence format issues

#### 3. "No mutations found"
- This is normal! Many sequences won't have mutations
- The tool will still report what it found
- Use this data to build your "wild-type" baseline

#### 4. Memory/Performance Issues
- Large genomes: Use `--log-level WARNING` to reduce output
- Many sequences: Process in smaller batches
- Slow alignments: Check reference sequence quality

## Integration with Your Research

### For Nishad Sir's Analysis:
1. **Run on multiple resistant isolates** - Get baseline mutation patterns
2. **Compare resistant vs sensitive** - Identify discriminatory mutations
3. **Build classification models** - Predict resistance from mutations
4. **Generate publication figures** - Mutation frequency plots, co-occurrence networks

### Next Steps After Single Sequence:
1. **Collect 50-100 sequences** from different sources
2. **Label with resistance data** (MIC values, susceptibility)
3. **Train ML models** to predict resistance
4. **Validate predictions** on held-out test set
5. **Generate insights** for clinical decision-making

## Quick Commands Reference

```bash
# Analyze single sequence
python single_sequence_analyzer.py --genome-file genome.fna --gff-file genome.gff

# Custom output
python single_sequence_analyzer.py --genome-file genome.fna --gff-file genome.gff --output-dir my_results

# Quiet mode (less logging)
python single_sequence_analyzer.py --genome-file genome.fna --gff-file genome.gff --log-level WARNING

# Check help
python single_sequence_analyzer.py --help
```

This tool gives you exactly what you need: targeted analysis of your specific mutations with the ability to scale up for robust ML training data! 🎯