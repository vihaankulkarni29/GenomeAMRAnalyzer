# Docker Data Directories

## data_input/
Place your input files here before running the Docker pipeline:
- **Accession list file** (e.g., `accessions.txt`) - List of genome accession numbers
- **Gene list file** (e.g., `genes.txt`) - List of genes to analyze

## data_output/
After running the pipeline, all results will be saved here:
- Analysis reports (HTML, CSV, JSON)
- Log files
- Extracted protein sequences
- AMR gene analysis results

## Usage
1. Copy your input files to `data_input/`
2. Run `./run_docker.sh` (Linux/Mac) or `run_docker.bat` (Windows)
3. Find results in `data_output/`