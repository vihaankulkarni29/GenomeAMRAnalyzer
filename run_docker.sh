#!/bin/bash

# Task 3: Interactive Docker Runner Script for GenomeAMRAnalyzer
# This script provides a user-friendly interface for the Docker workflow

set -e  # Exit on any error

echo "🧬 GenomeAMRAnalyzer Docker Runner"
echo "=================================="

# Check for Docker
echo "🔍 Checking for Docker..."
if ! command -v docker &> /dev/null; then
    echo "❌ ERROR: Docker is not installed or not in PATH."
    echo "Please install Docker from https://www.docker.com/get-started"
    exit 1
fi

if ! command -v docker-compose &> /dev/null; then
    echo "❌ ERROR: docker-compose is not installed or not in PATH."
    echo "Please install docker-compose."
    exit 1
fi

echo "✅ Docker found!"

# Create I/O directories
echo "📁 Creating data directories..."
mkdir -p ./data_input
mkdir -p ./data_output
echo "✅ Created ./data_input and ./data_output directories"

# Build the Docker image
echo "🏗️  Building Docker image (this may take a few minutes on first run)..."
if ! docker-compose build; then
    echo "❌ ERROR: Failed to build Docker image."
    exit 1
fi
echo "✅ Docker image built successfully!"

# Prompt for inputs
echo ""
echo "📋 Input File Configuration"
echo "Please ensure your files are placed in ./data_input/ directory"
echo ""

# Check if data_input has any files
if [ ! "$(ls -A ./data_input)" ]; then
    echo "⚠️  WARNING: ./data_input directory is empty!"
    echo "Please place your accession list and gene list files in ./data_input/"
    echo "Then run this script again."
    exit 1
fi

echo "Files found in ./data_input:"
ls -la ./data_input/

echo ""
read -p "📄 Enter the filename of your accession list (e.g., accessions.txt): " accession_file
read -p "🧬 Enter the filename of your gene list (e.g., genes.txt): " gene_file

# Validate inputs
if [ ! -f "./data_input/$accession_file" ]; then
    echo "❌ ERROR: File './data_input/$accession_file' not found!"
    exit 1
fi

if [ ! -f "./data_input/$gene_file" ]; then
    echo "❌ ERROR: File './data_input/$gene_file' not found!"
    exit 1
fi

echo "✅ Input files validated!"

# Prompt for email
read -p "📧 Enter your email address for NCBI queries: " email

# Run the pipeline
echo ""
echo "🚀 Starting GenomeAMRAnalyzer pipeline..."
echo "Input files: $accession_file, $gene_file"
echo "Output directory: ./data_output"
echo ""

# Run the container with the user's inputs
# Note: We use paths inside the container (/app/data_input/)
if docker-compose run --rm analyzer \
    --accessions "/app/data_input/$accession_file" \
    --genes "/app/data_input/$gene_file" \
    --email "$email" \
    --output "/app/data_output"; then
    
    echo ""
    echo "🎉 Pipeline completed successfully!"
    echo "📊 Results are available in ./data_output/"
    echo ""
    echo "Output files:"
    ls -la ./data_output/
    
else
    echo ""
    echo "❌ Pipeline failed. Check the error messages above."
    echo "Common issues:"
    echo "- Invalid accession format in accession list"
    echo "- Network connectivity problems"
    echo "- Invalid email address"
    exit 1
fi