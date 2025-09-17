#!/usr/bin/env python3
"""
Input Validation Script for AMR Pipeline
=========================================

Validates accession list and gene list before pipeline execution.
"""

import sys
from pathlib import Path

def validate_accessions(file_path):
    """Validate accession list for completeness and format"""
    if not Path(file_path).exists():
        print(f'ERROR: Accession file not found: {file_path}')
        return False
        
    with open(file_path, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]
        
    if not lines:
        print('ERROR: Accession file is empty')
        return False
        
    print(f'Found {len(lines)} accessions')
    
    # Check for duplicates
    unique_lines = list(set(lines))
    if len(unique_lines) != len(lines):
        print(f'WARNING: Found {len(lines) - len(unique_lines)} duplicate accessions')
        
    # Validate format
    valid_formats = 0
    for line in lines[:5]:  # Check first 5
        if any(line.startswith(prefix) for prefix in ['GCF_', 'GCA_', 'NZ_', 'NC_', 'CP']):
            valid_formats += 1
            
    if valid_formats == 0:
        print('ERROR: No valid accession formats found')
        return False
        
    return True

def validate_genes(file_path):
    """Validate gene list"""
    if not Path(file_path).exists():
        print(f'ERROR: Gene file not found: {file_path}')
        return False
        
    with open(file_path, 'r') as f:
        genes = [line.strip() for line in f if line.strip()]
        
    if not genes:
        print('ERROR: Gene file is empty')
        return False
        
    print(f'Found {len(genes)} genes: {genes}')
    return True

def main():
    """Main validation function"""
    print('Validating input files...')
    
    acc_valid = validate_accessions('accession_list.txt')
    gene_valid = validate_genes('gene_list.txt')
    
    if not acc_valid or not gene_valid:
        print('❌ Input validation failed')
        sys.exit(1)
    else:
        print('✅ Input validation passed')

if __name__ == '__main__':
    main()