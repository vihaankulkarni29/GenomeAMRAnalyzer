"""
Priority 2 AMR Analysis Pipeline
-------------------------------
High-performance, scalable pipeline for antimicrobial resistance genome analysis
with focus on RND efflux pumps and co-occurrence pattern detection.
"""

__version__ = "2.0.0"
__author__ = "AMR Analysis Pipeline Team"

# Import main modules
from .core import *
from .pipelines import *
from .integrations import *
from .config import *

# Main pipeline interface
from .core.pipeline_orchestrator import PipelineOrchestrator
from .config.pipeline_config import PipelineConfig, ConfigManager

__all__ = [
    'PipelineOrchestrator',
    'PipelineConfig', 
    'ConfigManager'
]