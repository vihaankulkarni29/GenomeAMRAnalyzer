"""
Priority 2 Core Module
---------------------
High-performance AMR genome analysis pipeline components.
"""

__version__ = "2.0.0"
__author__ = "AMR Analysis Pipeline Team"

from .enhanced_sequence_processor import EnhancedSequenceProcessor, ProcessingStats
from .data_quality_controller import DataQualityController
from .performance_optimizer import PerformanceOptimizer
from .export_manager import ExportManager
from .alignment_analyzer import AlignmentAnalyzer
from .pipeline_orchestrator import PipelineOrchestrator, PipelineConfig, PipelineResults

__all__ = [
    'EnhancedSequenceProcessor',
    'ProcessingStats',
    'DataQualityController', 
    'PerformanceOptimizer',
    'ExportManager',
    'AlignmentAnalyzer',
    'PipelineOrchestrator',
    'PipelineConfig',
    'PipelineResults'
]