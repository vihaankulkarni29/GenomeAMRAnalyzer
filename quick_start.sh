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
if ! command -v python &> /dev/null; then
    if ! command -v python3 &> /dev/null; then
        echo "ERROR: Python is not installed"
        echo
        echo "Please install Python:"
        echo "  Mac: brew install python"
        echo "  Ubuntu/Debian: sudo apt install python3 python3-pip"
        echo "  CentOS/RHEL: sudo yum install python3 python3-pip"
        echo
        exit 1
    else
        alias python=python3
    fi
fi

echo "Python found! Setting up your analysis..."
echo

# Install requirements if needed
if [ -f requirements.txt ]; then
    echo "Installing Python packages..."
    pip install -r requirements.txt
    echo
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
        python genomeamr_auto.py \
            --url "https://www.ncbi.nlm.nih.gov/nuccore/(Escherichia coli AND erythromycin resistance)" \
            --genes config/genes_erythromycin.txt \
            --email "$email"
        ;;
    2)
        echo
        read -p "Enter your email address: " email
        echo
        echo "Running beta-lactam resistance analysis..."
        python genomeamr_auto.py \
            --url "https://www.ncbi.nlm.nih.gov/nuccore/(clinical isolates AND ESBL)" \
            --genes config/genes_betalactam.txt \
            --email "$email"
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
            python genomeamr_auto.py --url "$genomes" --genes "$genes" --email "$email"
        else
            python genomeamr_auto.py --accessions "$genomes" --genes "$genes" --email "$email"
        fi
        ;;
    4)
        echo
        echo "All available options:"
        python genomeamr_auto.py --help
        ;;
    *)
        echo
        echo "Invalid choice. Please enter 1, 2, 3, or 4."
        exit 1
        ;;
esac

echo
echo "========================================"
echo "Analysis complete!"
echo "Check the 'reports' folder for results."
echo "========================================"
echo