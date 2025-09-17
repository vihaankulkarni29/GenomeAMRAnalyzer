# GenomeAMRAnalyzer Dockerfile
# Maintainer: Vihaan (2025) | Built for robust, reproducible AMR analysis

# 1. Use official lightweight Ubuntu as base
FROM ubuntu:22.04

# 2. Set noninteractive mode for apt
ENV DEBIAN_FRONTEND=noninteractive

# 3. Install system dependencies
RUN apt-get update && apt-get install -y \
    wget \
    git \
    build-essential \
    python3 \
    python3-pip \
    python3-venv \
    ca-certificates \
    libbz2-dev \
    liblzma-dev \
    libffi-dev \
    libssl-dev \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# 4. Create working directory
WORKDIR /app

# 5. Copy pipeline code and configs
COPY . /app


# 6. Install Miniconda (for RGI and bioinformatics tools)
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh && \
    /opt/conda/bin/conda clean -afy

# 7. Add conda to PATH
ENV PATH="/opt/conda/bin:$PATH"

# 7.5 Accept Conda Terms of Service for default channels
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r


# 8. Create and activate pipeline environment, install dependencies from bioconda and conda-forge
RUN conda create -y -n genomeamr -c conda-forge -c bioconda python=3.10 snakemake biopython rgi && \
    conda clean -afy

# 9. Activate environment by default
SHELL ["/bin/bash", "-c"]
RUN echo "conda activate genomeamr" >> ~/.bashrc

# 10. Set default command to bash (for interactive use)
CMD ["bash"]

# 11. Usage notes:
# - Mount your data with -v /host/data:/app/data
# - Run pipeline: snakemake --cores 4 --configfile config/snakemake_config.yaml
# - For RGI: rgi main ...
