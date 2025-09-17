"""
GenomeAMRAnalyzer: A comprehensive pipeline for antimicrobial resistance analysis

This package provides tools for analyzing antimicrobial resistance mutations
in bacterial genomes, with a focus on RND efflux pump systems.
"""

__version__ = "1.0.0"
__author__ = "GenomeAMRAnalyzer Team"

# Core components
from .generic_cooccurrence_analyzer import GenericCoOccurrenceAnalyzer, MutationEvent
from .fasta_aa_extractor_integration import FastaAAExtractorIntegrator, GeneCoordinate, ExtractionResult

__all__ = [
    'GenericCoOccurrenceAnalyzer',
    'MutationEvent', 
    'FastaAAExtractorIntegrator',
    'GeneCoordinate',
    'ExtractionResult'
]