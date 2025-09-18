@echo off
REM Task 3: Interactive Docker Runner Script for GenomeAMRAnalyzer (Windows)
REM This script provides a user-friendly interface for the Docker workflow

echo 🧬 GenomeAMRAnalyzer Docker Runner
echo ==================================

REM Check for Docker
echo 🔍 Checking for Docker...
docker --version >nul 2>&1
if errorlevel 1 (
    echo ❌ ERROR: Docker is not installed or not in PATH.
    echo Please install Docker from https://www.docker.com/get-started
    pause
    exit /b 1
)

docker-compose --version >nul 2>&1
if errorlevel 1 (
    echo ❌ ERROR: docker-compose is not installed or not in PATH.
    echo Please install docker-compose.
    pause
    exit /b 1
)

echo ✅ Docker found!

REM Create I/O directories
echo 📁 Creating data directories...
if not exist "data_input" mkdir data_input
if not exist "data_output" mkdir data_output
echo ✅ Created data_input and data_output directories

REM Build the Docker image
echo 🏗️  Building Docker image (this may take a few minutes on first run)...
docker-compose build
if errorlevel 1 (
    echo ❌ ERROR: Failed to build Docker image.
    pause
    exit /b 1
)
echo ✅ Docker image built successfully!

REM Prompt for inputs
echo.
echo 📋 Input File Configuration
echo Please ensure your files are placed in data_input\ directory
echo.

REM Check if data_input has any files
dir /b data_input 2>nul | findstr . >nul
if errorlevel 1 (
    echo ⚠️  WARNING: data_input directory is empty!
    echo Please place your accession list and gene list files in data_input\
    echo Then run this script again.
    pause
    exit /b 1
)

echo Files found in data_input:
dir data_input

echo.
set /p accession_file="📄 Enter the filename of your accession list (e.g., accessions.txt): "
set /p gene_file="🧬 Enter the filename of your gene list (e.g., genes.txt): "

REM Validate inputs
if not exist "data_input\%accession_file%" (
    echo ❌ ERROR: File 'data_input\%accession_file%' not found!
    pause
    exit /b 1
)

if not exist "data_input\%gene_file%" (
    echo ❌ ERROR: File 'data_input\%gene_file%' not found!
    pause
    exit /b 1
)

echo ✅ Input files validated!

REM Prompt for email
set /p email="📧 Enter your email address for NCBI queries: "

REM Run the pipeline
echo.
echo 🚀 Starting GenomeAMRAnalyzer pipeline...
echo Input files: %accession_file%, %gene_file%
echo Output directory: data_output
echo.

REM Run the container with the user's inputs
docker-compose run --rm analyzer --accessions "/app/data_input/%accession_file%" --genes "/app/data_input/%gene_file%" --email "%email%" --output-dir "/app/data_output"

if errorlevel 1 (
    echo.
    echo ❌ Pipeline failed. Check the error messages above.
    echo Common issues:
    echo - Invalid accession format in accession list
    echo - Network connectivity problems
    echo - Invalid email address
    pause
    exit /b 1
) else (
    echo.
    echo 🎉 Pipeline completed successfully!
    echo 📊 Results are available in data_output\
    echo.
    echo Output files:
    dir data_output
    pause
)