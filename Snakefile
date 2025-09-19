# AMR Genomics Pipeline - Snakemake Workflow
# ===========================================
# 
# Research Question: "How are antimicrobial resistance genomes developing mutations?"
# Pipeline Flow: NCBI → Genome Harvester → CARD → FastaAAExtractor → Aligner → SubScan → Co-occurrence → Report
#
# Senior Bioinformatics Design:
# - Preserves robust Python error handling
# - Adds bioinformatics tool management
# - Enables reproducible research
# - Supports scalable execution

import os
from pathlib import Path

# Configuration
configfile: "config/snakemake_config.yaml"

# Define samples from accession list
ACCESSIONS = []
if Path("accession_list.txt").exists():
    with open("accession_list.txt", 'r') as f:
        ACCESSIONS = [line.strip() for line in f if line.strip() and not line.startswith('#')]

# Define target genes from gene list
GENES = []
if Path("gene_list.txt").exists():
    with open("gene_list.txt", 'r') as f:
        GENES = [line.strip() for line in f if line.strip() and not line.startswith('#')]

# Final target: HTML report
rule all:
    input:
        "reports/amr_analysis_report.html",
        "results/pipeline_summary.json"

# Rule 1: Accession Conversion & Validation
rule convert_accessions:
    input:
        accessions="accession_list.txt"
    output:
        converted="results/converted_accessions.txt",
        log="logs/accession_conversion.log"
    conda:
        "envs/python_base.yaml"
    log:
        "logs/convert_accessions.log"
    script:
        "scripts/robust_accession_converter.py"

# Rule 2: Genome Harvesting
rule harvest_genomes:
    input:
        accessions="results/converted_accessions.txt"
    output:
        genomes=expand("genomes/{accession}.fasta", accession=ACCESSIONS),
        checkpoint="results/harvest_checkpoint.json"
    params:
        output_dir="genomes",
        email=config["ncbi_email"]
    conda:
        "envs/python_base.yaml"
    log:
        "logs/genome_harvest.log"
    script:
        "scripts/genome_harvester.py"

# Rule 3: Abricate Analysis (replaces RGI)
rule run_abricate:
    input:
        genome="genomes/{accession}.fasta"
    output:
        abricate_tsv="abricate_results/{accession}_abricate.tsv",
        coords_csv="card_results/{accession}_coordinates.csv"
    params:
        card_db="card"
    conda:
        "envs/abricate.yaml"
    log:
        "logs/abricate_{accession}.log"
    shell:
        """
        # Run Abricate with CARD database
        abricate --db {params.card_db} --nopath {input.genome} > {output.abricate_tsv} 2> {log}
        
        # Convert to coordinate CSV format
        python -m src.abricate_to_coords --in-tsv {output.abricate_tsv} --out-csv {output.coords_csv} 2>> {log}
        """

# Rule 4: Protein Extraction
rule extract_proteins:
    input:
        genome="genomes/{accession}.fasta",
        coords="card_results/{accession}_coordinates.csv"
        genome="genomes/{accession}.fasta",
        rgi_txt="card_results/{accession}_rgi.txt",
        gene_list="gene_list.txt"
    output:
        proteins=expand("proteins/{{accession}}_{gene}.faa", gene=GENES)
    conda:
        "envs/python_bio.yaml"
    log:
        "logs/extract_proteins_{accession}.log"
    script:
        "scripts/extract_proteins_snakemake.py"

# Rule 5: Sequence Alignment
rule align_sequences:
    input:
        proteins=expand("proteins/{accession}_{gene}.faa", accession=ACCESSIONS, gene=GENES),
        reference="references/{gene}_reference.faa"
    output:
        alignment="alignments/{gene}_alignment.aln"
    conda:
        "envs/alignment.yaml"
    log:
        "logs/align_{gene}.log"
    shell:
        """
        muscle -align proteins/*_{wildcards.gene}.faa \
               -ref {input.reference} \
               -output {output.alignment} 2> {log}
        """

# Rule 6: Mutation Detection (SubScan)
rule detect_mutations:
    input:
        alignment="alignments/{gene}_alignment.aln"
    output:
        mutations="mutations/{gene}_mutations.csv"
    conda:
        "envs/python_bio.yaml"
    log:
        "logs/mutations_{gene}.log"
    script:
        "scripts/subscan_snakemake.py"

# Rule 7: Co-occurrence Analysis
rule cooccurrence_analysis:
    input:
        mutations=expand("mutations/{gene}_mutations.csv", gene=GENES)
    output:
        cooccurrence="results/cooccurrence_analysis.csv",
        matrix="results/cooccurrence_matrix.json"
    conda:
        "envs/python_analysis.yaml"
    log:
        "logs/cooccurrence.log"
    script:
        "scripts/cooccurrence_snakemake.py"

# Rule 8: MIC Data Harvesting
rule harvest_mic_data:
    input:
        genomes=expand("genomes/{accession}.fasta", accession=ACCESSIONS)
    output:
        mic_data="results/mic_metadata.csv"
    conda:
        "envs/python_base.yaml"
    log:
        "logs/mic_harvest.log"
    script:
        "scripts/mic_harvester_snakemake.py"

# Rule 9: Generate Final Report
rule generate_report:
    input:
        genomes=expand("genomes/{accession}.fasta", accession=ACCESSIONS),
        mutations=expand("mutations/{gene}_mutations.csv", gene=GENES),
        cooccurrence="results/cooccurrence_analysis.csv",
        mic_data="results/mic_metadata.csv"
    output:
        report="reports/amr_analysis_report.html",
        summary="results/pipeline_summary.json"
    conda:
        "envs/python_reporting.yaml"
    log:
        "logs/generate_report.log"
    script:
        "scripts/report_generator_snakemake.py"

# Utility rule: Clean intermediate files
rule clean:
    shell:
        """
        rm -rf genomes/* card_results/* proteins/* alignments/* mutations/*
        rm -rf logs/* results/* reports/*
        """

# Utility rule: Setup conda environments
rule setup_environments:
    output:
        touch("envs/.setup_complete")
    shell:
        """
        conda env create -f envs/rgi.yaml
        conda env create -f envs/alignment.yaml  
        conda env create -f envs/python_bio.yaml
        conda env create -f envs/python_analysis.yaml
        conda env create -f envs/python_reporting.yaml
        """