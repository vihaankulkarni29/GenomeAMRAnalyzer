#!/bin/bash

#==============================================================
# GenomeAMRAnalyzer - One-Click Setup for Mac/Linux
# For Project Investigators & Non-Technical Users  
#==============================================================

clear
echo
echo "========================================"
echo " GenomeAMRAnalyzer - Quick Start"
echo "========================================"
echo

# Check if Python is available
PYTHON_CMD="python"
if ! command -v python &> /dev/null; then
    if command -v python3 &> /dev/null; then
        PYTHON_CMD="python3"
        echo "Using python3..."
    else
        echo "ERROR: Python is not installed"
        echo
        echo "Please install Python:"
        echo "  Mac: brew install python"
        echo "  Ubuntu/Debian: sudo apt install python3 python3-pip"
        echo "  CentOS/RHEL: sudo yum install python3 python3-pip"
        echo "  Windows: Use Windows PowerShell and run quick_start.bat instead"
        echo
        echo "Press Enter to exit..."
        read -r
        exit 1
    fi
fi

echo "Python found! Setting up your analysis..."
echo

# Install requirements if needed
if [ -f requirements.txt ]; then
    echo "Installing Python packages..."
    
    # Try pip install with fallback options
    if pip install -r requirements.txt 2>/dev/null; then
        echo "✅ Packages installed successfully"
    elif pip install -r requirements.txt --user 2>/dev/null; then
        echo "✅ Packages installed to user directory"
    elif pip install -r requirements.txt --break-system-packages 2>/dev/null; then
        echo "✅ Packages installed (system override)"
    else
        echo "⚠️ Warning: Package installation failed, but analysis may still work"
        echo "   If you encounter import errors, please install packages manually:"
        echo "   pip install -r requirements.txt --user"
    fi
    echo
fi

# Validate required files exist
missing_files=()
if [ ! -f "genomeamr_auto.py" ]; then
    missing_files+=("genomeamr_auto.py")
fi
if [ ! -f "config/genes_default.txt" ]; then
    missing_files+=("config/genes_default.txt")
fi
if [ ! -f "config/genes_erythromycin.txt" ]; then
    missing_files+=("config/genes_erythromycin.txt")
fi

if [ ${#missing_files[@]} -gt 0 ]; then
    echo "❌ ERROR: Missing required files:"
    for file in "${missing_files[@]}"; do
        echo "  - $file"
    done
    echo
    echo "Please make sure you're running this script from the GenomeAMRAnalyzer directory."
    echo "Press Enter to exit..."
    read -r
    exit 1
fi

echo "Ready! Choose your analysis:"
echo
echo "1. Erythromycin resistance (E. coli)"
echo "2. Beta-lactam resistance (clinical isolates)"
echo "3. Custom analysis (your own inputs)"
echo "4. Show all options"
echo

read -p "Enter your choice (1-4): " choice

case $choice in
    1)
        echo
        read -p "Enter your email address: " email
        echo
        echo "Running erythromycin resistance analysis..."
        $PYTHON_CMD genomeamr_auto.py \
            --url "https://www.ncbi.nlm.nih.gov/nuccore/(Escherichia coli AND erythromycin resistance)" \
            --genes config/genes_erythromycin.txt \
            --email "$email"
        
        # Check if command succeeded
        if [ $? -eq 0 ]; then
            echo "✅ Analysis completed successfully!"
        else
            echo "❌ Analysis failed. Please check the error messages above."
        fi
        ;;
    2)
        echo
        read -p "Enter your email address: " email
        echo
        echo "Running beta-lactam resistance analysis..."
        $PYTHON_CMD genomeamr_auto.py \
            --url "https://www.ncbi.nlm.nih.gov/nuccore/(clinical isolates AND ESBL)" \
            --genes config/genes_betalactam.txt \
            --email "$email"
        
        # Check if command succeeded
        if [ $? -eq 0 ]; then
            echo "✅ Analysis completed successfully!"
        else
            echo "❌ Analysis failed. Please check the error messages above."
        fi
        ;;
    3)
        echo
        echo "For custom analysis, you need:"
        echo "1. Your email address"
        echo "2. Either a genome list file OR an NCBI search URL"
        echo "3. A genes file (or use config/genes_default.txt)"
        echo
        read -p "Enter your email address: " email
        read -p "Enter genome file path OR NCBI URL: " genomes
        read -p "Enter genes file path (or press Enter for default): " genes
        
        if [ -z "$genes" ]; then
            genes="config/genes_default.txt"
        fi
        
        echo
        echo "Running custom analysis..."
        
        # Check if genomes parameter is a URL
        if [[ $genomes == http* ]]; then
            $PYTHON_CMD genomeamr_auto.py --url "$genomes" --genes "$genes" --email "$email"
        else
            $PYTHON_CMD genomeamr_auto.py --accessions "$genomes" --genes "$genes" --email "$email"
        fi
        
        # Check if command succeeded
        if [ $? -eq 0 ]; then
            echo "✅ Analysis completed successfully!"
        else
            echo "❌ Analysis failed. Please check the error messages above."
        fi
        ;;
    4)
        echo
        echo "All available options:"
        $PYTHON_CMD genomeamr_auto.py --help
        echo
        echo "Press Enter to continue..."
        read -r
        exit 0
        ;;
    *)
        echo
        echo "Invalid choice. Please enter 1, 2, 3, or 4."
        echo "Press Enter to exit..."
        read -r
        exit 1
        ;;
esac

echo
echo "========================================"
echo "Analysis complete!"
echo "Check the 'reports' folder for results."
echo "========================================"
echo
echo "Press Enter to exit..."
read -r