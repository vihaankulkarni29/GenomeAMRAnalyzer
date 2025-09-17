"""
Priority 3 - Advanced Genomic Data Collection Pipeline
====================================================

Enterprise-grade genomic data harvesting with:
- NCBI GenBank integration with checkpointing
- MIC metadata enrichment 
- Advanced mutation analysis
- Intelligent co-occurrence detection
- Publication-ready HTML reports
- MCP-powered data processing acceleration

NO biological interpretation - pure data collection and processing.
"""

from typing import TYPE_CHECKING

__version__ = "0.1.0"
__author__ = "GenomeAMRAnalyzer Team"

# Avoid importing heavy modules at package import time to keep tests lightweight
if TYPE_CHECKING:
    pass

__all__ = []