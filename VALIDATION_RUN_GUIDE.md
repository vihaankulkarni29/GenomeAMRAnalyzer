# End-to-End Pipeline Validation Checklist

This guide outlines the steps to perform a full, manual validation of the GenomeAMRAnalyzer pipeline using the Docker environment.

### ✅ 1. Prepare Validation Dataset
- [ ] In the `data_input/` directory, create a new sub-directory named `validation_set`.
- [ ] Copy 2-3 known genome FASTA files (e.g., from a previous successful run) into `validation_set/`.
- [ ] Create a `validation_genes.txt` file inside `validation_set/` and add the names of 2-3 genes you expect to find in the sample genomes.

### ✅ 2. Execute the Pipeline
- [ ] Open your terminal in the project root.
- [ ] Run the command: `./run_docker.sh`
- [ ] When prompted, enter the paths to your validation files (e.g., `validation_set/genome1.fasta`, `validation_set/validation_genes.txt`).

### ✅ 3. Deeply Inspect the HTML Report
- [ ] Once the run is complete, open the `amr_report.html` file from the `data_output/` directory in a web browser.
- [ ] **Critical Review Checklist:**
    - [ ] **Run Summary:** Does the summary table show the correct number of genomes and genes?
    - [ ] **Risk Profiling:** Is the color-coded risk assessment logical based on the findings?
    - [ ] **Plotly Chart:** Is the mutation frequency chart rendered correctly and is it interactive?
    - [ ] **Genomic Context Tables (VF & Plasmid):** Are the virulence factor and plasmid tables populated correctly? Do they show "None detected" if nothing was found?
    - [ ] **Main Results Table:** Is the data scientifically plausible? Are the confidence scores calculated as expected?
    - [ ] **Data Integrity:** Cross-reference a few specific results in the report with the raw `.tsv` files in `data_output/` to ensure the data has been parsed and presented accurately.
