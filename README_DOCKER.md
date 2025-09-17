# GenomeAMRAnalyzer: Docker & Pipeline Usage Guide

## 1. Prerequisites
- Install Docker Desktop (https://www.docker.com/products/docker-desktop/)
- Ensure WSL2 integration is enabled (recommended for Windows)

## 2. Build the Docker Image
Open PowerShell in your project directory and run:

    docker build -t genome-amr-analyzer .

- This command reads the Dockerfile and creates a portable image with all dependencies.
- The `-t` flag names your image for easy reference.

## 3. Run the Pipeline in a Container
To execute the pipeline, mount your data and run Snakemake:

    docker run -it --rm -v ${PWD}:/app genome-amr-analyzer snakemake --cores 4 --configfile config/snakemake_config.yaml

- `-it` gives you an interactive shell.
- `--rm` cleans up the container after exit.
- `-v ${PWD}:/app` mounts your project folder inside the container.
- You can adjust `--cores` for parallel execution.

## 4. RGI and Other Tools
To run RGI or other tools interactively:

    docker run -it --rm -v ${PWD}:/app genome-amr-analyzer bash

Then inside the container:

    rgi main ...

## 5. Data & Results
- All files in your project folder are accessible inside the container at `/app`.
- Output files will appear in your local folder after the run.

## 6. Troubleshooting
- If you encounter permission issues, try running PowerShell as Administrator.
- For large datasets, increase memory allocation in Docker Desktop settings.

## 7. Reproducibility & Sharing
- Share your Dockerfile and README for publication-grade reproducibility.
- Optionally, push your image to Docker Hub for collaborators:

    docker tag genome-amr-analyzer yourusername/genome-amr-analyzer:latest
    docker push yourusername/genome-amr-analyzer:latest

## 8. Advanced: WSL2 Setup (Optional)
- WSL2 allows you to run Linux tools natively on Windows. For advanced users, you can run the pipeline directly in Ubuntu via WSL2.
- See https://docs.microsoft.com/en-us/windows/wsl/ for setup instructions.

---

**This setup ensures your pipeline runs identically on any system, supporting robust, publication-ready AMR analysis.**
