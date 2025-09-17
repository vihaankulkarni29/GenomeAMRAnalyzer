@echo off
echo ==========================================
echo AMR PIPELINE - QUICK RGI SETUP
echo ==========================================
echo.
echo CRITICAL: Your pipeline needs RGI for research functionality
echo.
echo Current Status:
echo ✅ 11 genomes downloaded successfully
echo ❌ No RGI = No resistance gene analysis  
echo ❌ No protein extraction = No mutations
echo ❌ No research insights = Tool is redundant
echo.
echo ==========================================
echo SOLUTION: Install Bioinformatics Tools
echo ==========================================
echo.
echo 1. Install Miniconda from: https://docs.conda.io/en/latest/miniconda.html
echo 2. Run: conda create -n amr_analysis python=3.9
echo 3. Run: conda activate amr_analysis  
echo 4. Run: conda install -c bioconda rgi
echo 5. Run: rgi download -d card
echo 6. Re-run this pipeline
echo.
echo ==========================================
echo ALTERNATIVE: Use existing genomes for now
echo ==========================================
echo.
echo Would you like to:
echo [1] See detailed installation guide
echo [2] Continue anyway (limited functionality)
echo [3] Exit to install RGI first
echo.
set /p choice="Enter choice (1-3): "

if "%choice%"=="1" (
    echo Opening installation guide...
    start notepad RGI_INSTALLATION_GUIDE.md
    goto :end
)

if "%choice%"=="2" (
    echo.
    echo WARNING: Continuing without RGI provides limited research value.
    echo Your 11 downloaded genomes are available in genome_data/
    echo.
    echo Genomes successfully processed from your NZ_CP151* accessions:
    dir genome_data\*.fasta /b
    echo.
    echo To get research insights, please install RGI and re-run.
    goto :end
)

if "%choice%"=="3" (
    echo.
    echo Excellent choice! Install RGI first for full functionality.
    echo See RGI_INSTALLATION_GUIDE.md for detailed steps.
    goto :end
)

:end
pause