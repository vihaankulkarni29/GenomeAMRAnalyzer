# SubScan: Targeted Antimicrobial Resistance Gene Mutation Analysis Pipeline
# Official Docker image for production use

# Use the official Miniconda3 image as a base
FROM continuumio/miniconda3

# Metadata
LABEL maintainer="GenomeAMRAnalyzer Team <genomeamr@example.com>"
LABEL version="1.0.0"
LABEL description="SubScan: A tool for targeted, reference-driven mutation analysis of antimicrobial resistance genes."
LABEL org.opencontainers.image.title="SubScan"
LABEL org.opencontainers.image.description="Comprehensive pipeline for AMR gene mutation analysis with interactive reporting"
LABEL org.opencontainers.image.version="1.0.0"
LABEL org.opencontainers.image.source="https://github.com/vihaankulkarni29/GenomeAMRAnalyzer"
LABEL org.opencontainers.image.documentation="https://github.com/vihaankulkarni29/GenomeAMRAnalyzer/blob/main/README.md"
LABEL org.opencontainers.image.licenses="MIT"

# Set the working directory inside the container
WORKDIR /app

# Copy the entire project context into the container's working directory
COPY . .

# Configure Conda channels, create the environment, and install Abricate (replaces RGI)
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda config --set channel_priority strict && \
    conda create -n GenomeAMRAnalyzer_env python=3.11 abricate -y && \
    conda clean -afy

# Install Python packages from requirements.txt into the created Conda environment
# Note: We use the specific pip from our new environment to ensure correct installation
RUN /opt/conda/envs/GenomeAMRAnalyzer_env/bin/pip install --no-cache-dir -r requirements.txt

# Set the entrypoint to run the main script within the Conda environment
# This makes the container act like an executable for our pipeline
ENTRYPOINT ["conda", "run", "-n", "GenomeAMRAnalyzer_env", "python", "genomeamr_auto.py"]
