REM ==============================================================
REM  GenomeAMRAnalyzer - One-Click Setup for Windows
REM  For Project Investigators & Non-Technical Users
REM ==============================================================

@echo off
cls
echo.
echo ========================================
echo  GenomeAMRAnalyzer - Quick Start
echo ========================================
echo.

REM Check if Python is available
python --version >nul 2>&1
if %errorlevel% neq 0 (
    echo ERROR: Python is not installed or not in PATH
    echo.
    echo Please install Python from: https://www.python.org/downloads/
    echo Make sure to check "Add Python to PATH" during installation
    echo.
    pause
    exit /b 1
)

echo Python found! Setting up your analysis...
echo.

REM Install requirements if needed
if exist requirements.txt (
    echo Installing Python packages...
    pip install -r requirements.txt
    echo.
)

echo Ready! Choose your analysis:
echo.
echo 1. Erythromycin resistance (E. coli)
echo 2. Beta-lactam resistance (clinical isolates)  
echo 3. Custom analysis (your own inputs)
echo 4. Show all options
echo.

set /p choice="Enter your choice (1-4): "

if "%choice%"=="1" goto erythromycin
if "%choice%"=="2" goto betalactam
if "%choice%"=="3" goto custom
if "%choice%"=="4" goto help
goto invalid

:erythromycin
echo.
set /p email="Enter your email address: "
echo.
echo Running erythromycin resistance analysis...
python genomeamr_auto.py --url "https://www.ncbi.nlm.nih.gov/nuccore/(Escherichia coli AND erythromycin resistance)" --genes config/genes_erythromycin.txt --email %email%
goto end

:betalactam
echo.
set /p email="Enter your email address: "
echo.
echo Running beta-lactam resistance analysis...
python genomeamr_auto.py --url "https://www.ncbi.nlm.nih.gov/nuccore/(clinical isolates AND ESBL)" --genes config/genes_betalactam.txt --email %email%
goto end

:custom
echo.
echo For custom analysis, you need:
echo 1. Your email address
echo 2. Either a genome list file OR an NCBI search URL
echo 3. A genes file (or use config/genes_default.txt)
echo.
set /p email="Enter your email address: "
set /p genomes="Enter genome file path OR NCBI URL: "
set /p genes="Enter genes file path (or press Enter for default): "
if "%genes%"=="" set genes=config/genes_default.txt
echo.
echo Running custom analysis...

REM Check if genomes parameter is a file or URL
echo %genomes% | findstr /i "http" >nul
if %errorlevel% equ 0 (
    python genomeamr_auto.py --url "%genomes%" --genes "%genes%" --email %email%
) else (
    python genomeamr_auto.py --accessions "%genomes%" --genes "%genes%" --email %email%
)
goto end

:help
echo.
echo All available options:
python genomeamr_auto.py --help
goto end

:invalid
echo.
echo Invalid choice. Please enter 1, 2, 3, or 4.
pause
goto start

:end
echo.
echo ========================================
echo Analysis complete!
echo Check the 'reports' folder for results.
echo ========================================
echo.
pause