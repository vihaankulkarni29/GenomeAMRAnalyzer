# Task 1: Create the Dockerfile

# Use the official Miniconda3 image as a base
FROM continuumio/miniconda3

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
