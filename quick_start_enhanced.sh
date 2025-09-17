#!/bin/bash

#==============================================================
# GenomeAMRAnalyzer - Enhanced Quick Start for Mac/Linux
# With better error handling and user feedback
#==============================================================

# Colors for better output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

clear
echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE} GenomeAMRAnalyzer - Quick Start${NC}"
echo -e "${BLUE}========================================${NC}"
echo

# Function to pause and wait for user input
pause_for_user() {
    echo
    echo -e "${YELLOW}Press Enter to continue...${NC}"
    read -r
}

# Function to pause and exit
pause_and_exit() {
    echo
    echo -e "${YELLOW}Press Enter to exit...${NC}"
    read -r
    exit $1
}

# Detect operating system
if [[ "$OSTYPE" == "msys" || "$OSTYPE" == "win32" ]]; then
    echo -e "${YELLOW}⚠️ You appear to be on Windows.${NC}"
    echo "For the best experience on Windows, please use:"
    echo "  quick_start.bat (double-click or run in Command Prompt)"
    echo
    echo "If you prefer to continue with this bash script, press Enter..."
    read -r
fi

# Check if Python is available
PYTHON_CMD=""
if command -v python &> /dev/null; then
    PYTHON_CMD="python"
    echo -e "${GREEN}✅ Python found${NC}"
elif command -v python3 &> /dev/null; then
    PYTHON_CMD="python3"
    echo -e "${GREEN}✅ Python3 found${NC}"
else
    echo -e "${RED}❌ ERROR: Python is not installed${NC}"
    echo
    echo "Please install Python:"
    echo "  Mac: brew install python"
    echo "  Ubuntu/Debian: sudo apt install python3 python3-pip"
    echo "  CentOS/RHEL: sudo yum install python3 python3-pip"
    echo "  Windows: Use quick_start.bat instead of this script"
    pause_and_exit 1
fi

# Validate required files exist
echo "Checking required files..."
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
    echo -e "${RED}❌ ERROR: Missing required files:${NC}"
    for file in "${missing_files[@]}"; do
        echo "  - $file"
    done
    echo
    echo "Please make sure you're running this script from the GenomeAMRAnalyzer directory."
    echo "Current directory: $(pwd)"
    pause_and_exit 1
fi

echo -e "${GREEN}✅ All required files found${NC}"

# Install requirements if needed
if [ -f requirements.txt ]; then
    echo "Installing Python packages..."
    
    # Try different installation methods
    if $PYTHON_CMD -m pip install -r requirements.txt 2>/dev/null; then
        echo -e "${GREEN}✅ Packages installed successfully${NC}"
    elif $PYTHON_CMD -m pip install -r requirements.txt --user 2>/dev/null; then
        echo -e "${GREEN}✅ Packages installed to user directory${NC}"
    elif $PYTHON_CMD -m pip install -r requirements.txt --break-system-packages 2>/dev/null; then
        echo -e "${YELLOW}✅ Packages installed (system override)${NC}"
    else
        echo -e "${YELLOW}⚠️ Warning: Package installation failed${NC}"
        echo "   Analysis may still work with existing packages"
        echo "   If you encounter import errors, manually install:"
        echo "   $PYTHON_CMD -m pip install -r requirements.txt --user"
        pause_for_user
    fi
    echo
fi

echo -e "${GREEN}Ready! Choose your analysis:${NC}"
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
        if [ -z "$email" ]; then
            echo -e "${RED}❌ Email is required for NCBI access${NC}"
            pause_and_exit 1
        fi
        echo
        echo -e "${BLUE}Running erythromycin resistance analysis...${NC}"
        $PYTHON_CMD genomeamr_auto.py \
            --url "https://www.ncbi.nlm.nih.gov/nuccore/(Escherichia coli AND erythromycin resistance)" \
            --genes config/genes_erythromycin.txt \
            --email "$email"
        
        # Check if command succeeded
        if [ $? -eq 0 ]; then
            echo -e "${GREEN}✅ Analysis completed successfully!${NC}"
        else
            echo -e "${RED}❌ Analysis failed. Please check the error messages above.${NC}"
        fi
        ;;
    2)
        echo
        read -p "Enter your email address: " email
        if [ -z "$email" ]; then
            echo -e "${RED}❌ Email is required for NCBI access${NC}"
            pause_and_exit 1
        fi
        echo
        echo -e "${BLUE}Running beta-lactam resistance analysis...${NC}"
        $PYTHON_CMD genomeamr_auto.py \
            --url "https://www.ncbi.nlm.nih.gov/nuccore/(clinical isolates AND ESBL)" \
            --genes config/genes_betalactam.txt \
            --email "$email"
        
        # Check if command succeeded
        if [ $? -eq 0 ]; then
            echo -e "${GREEN}✅ Analysis completed successfully!${NC}"
        else
            echo -e "${RED}❌ Analysis failed. Please check the error messages above.${NC}"
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
        if [ -z "$email" ]; then
            echo -e "${RED}❌ Email is required for NCBI access${NC}"
            pause_and_exit 1
        fi
        
        read -p "Enter genome file path OR NCBI URL: " genomes
        if [ -z "$genomes" ]; then
            echo -e "${RED}❌ Genome source is required${NC}"
            pause_and_exit 1
        fi
        
        read -p "Enter genes file path (or press Enter for default): " genes
        if [ -z "$genes" ]; then
            genes="config/genes_default.txt"
        fi
        
        echo
        echo -e "${BLUE}Running custom analysis...${NC}"
        
        # Check if genomes parameter is a URL
        if [[ $genomes == http* ]]; then
            $PYTHON_CMD genomeamr_auto.py --url "$genomes" --genes "$genes" --email "$email"
        else
            # Check if file exists
            if [ ! -f "$genomes" ]; then
                echo -e "${RED}❌ Genome file not found: $genomes${NC}"
                pause_and_exit 1
            fi
            $PYTHON_CMD genomeamr_auto.py --accessions "$genomes" --genes "$genes" --email "$email"
        fi
        
        # Check if command succeeded
        if [ $? -eq 0 ]; then
            echo -e "${GREEN}✅ Analysis completed successfully!${NC}"
        else
            echo -e "${RED}❌ Analysis failed. Please check the error messages above.${NC}"
        fi
        ;;
    4)
        echo
        echo -e "${BLUE}All available options:${NC}"
        $PYTHON_CMD genomeamr_auto.py --help
        pause_and_exit 0
        ;;
    *)
        echo
        echo -e "${RED}❌ Invalid choice. Please enter 1, 2, 3, or 4.${NC}"
        pause_and_exit 1
        ;;
esac

echo
echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}Analysis complete!${NC}"
echo -e "${GREEN}Check the 'reports' folder for results.${NC}"
echo -e "${GREEN}========================================${NC}"
pause_and_exit 0