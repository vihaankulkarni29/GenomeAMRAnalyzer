#!/bin/bash

# Simple test script for debugging quick_start.sh
echo "Testing GenomeAMRAnalyzer files..."

# Check if we're in the right directory
if [ ! -f "genomeamr_auto.py" ]; then
    echo "❌ ERROR: genomeamr_auto.py not found"
    echo "Current directory: $(pwd)"
    echo "Files in directory:"
    ls -la
    exit 1
fi

echo "✅ genomeamr_auto.py found"

# Test the Python script
echo "Testing genomeamr_auto.py help..."
python genomeamr_auto.py --help

# Test gene files
echo "Checking gene files..."
if [ -f "config/genes_default.txt" ]; then
    echo "✅ config/genes_default.txt found"
    echo "Contents (first 5 lines):"
    head -5 config/genes_default.txt
else
    echo "❌ config/genes_default.txt not found"
fi

if [ -f "config/genes_erythromycin.txt" ]; then
    echo "✅ config/genes_erythromycin.txt found"
    echo "Contents (first 5 lines):"
    head -5 config/genes_erythromycin.txt
else
    echo "❌ config/genes_erythromycin.txt not found"
fi

echo "Test complete!"