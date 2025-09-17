#!/usr/bin/env python3
"""
Pipeline Results Validation Script
==================================

Validates pipeline outputs and generates execution summary.
"""

import os
import json
from pathlib import Path

def validate_pipeline_results():
    """Comprehensive validation of pipeline outputs"""
    print('\n=== Pipeline Results Summary ===')
    
    # Check genome downloads
    genome_dir = Path('genome_data')
    if genome_dir.exists():
        genome_files = list(genome_dir.glob('*.fasta')) + list(genome_dir.glob('*.fa')) + list(genome_dir.glob('*.fna'))
        print(f'Downloaded genomes: {len(genome_files)}')
    else:
        print('Downloaded genomes: 0 (directory not found)')
        
    # Check CARD results
    card_dir = Path('card_results')
    if card_dir.exists():
        card_files = list(card_dir.glob('*_rgi.txt'))
        print(f'CARD RGI results: {len(card_files)}')
    else:
        print('CARD RGI results: 0 (directory not found)')
        
    # Check extracted proteins
    protein_dir = Path('extracted_proteins')
    if protein_dir.exists():
        protein_files = list(protein_dir.glob('*.faa'))
        print(f'Extracted proteins: {len(protein_files)}')
    else:
        print('Extracted proteins: 0 (directory not found)')
        
    # Check alignments
    align_dir = Path('alignments')
    if align_dir.exists():
        align_files = list(align_dir.glob('*.aln'))
        print(f'Alignment files: {len(align_files)}')
    else:
        print('Alignment files: 0 (directory not found)')
        
    # Check cooccurrence results
    cooccur_dir = Path('cooccurrence_results')
    if cooccur_dir.exists():
        cooccur_files = list(cooccur_dir.glob('*.csv')) + list(cooccur_dir.glob('*.json'))
        print(f'Co-occurrence results: {len(cooccur_files)}')
    else:
        print('Co-occurrence results: 0 (directory not found)')
        
    # Check reports
    report_dir = Path('reports')
    if report_dir.exists():
        report_files = list(report_dir.glob('*.html'))
        print(f'HTML reports: {len(report_files)}')
    else:
        print('HTML reports: 0 (directory not found)')
        
    print('\n=== Pipeline Execution Complete ===')
    print('Check individual stage logs for detailed results.')
    
    # Check for database
    if Path('priority3.db').exists():
        print('✅ Pipeline database created successfully')
    else:
        print('⚠️  Pipeline database not found')

def main():
    """Main validation function"""
    validate_pipeline_results()

if __name__ == '__main__':
    main()