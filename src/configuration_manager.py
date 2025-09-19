
"""
Configuration Management System for GenomeAMRAnalyzer
Replaces hardcoded values with configurable parameters
"""

import json
import os
from pathlib import Path
from typing import Dict, List, Any, Optional

class ConfigurationManager:
    """
    Centralized configuration management system
    """
    
    def __init__(self, config_file: Optional[str] = None):
        """Initialize configuration manager"""
        self.config_file = config_file or self._find_config_file()
        self.config = self._load_configuration()
    
    def _find_config_file(self) -> str:
        """Find configuration file in standard locations"""
        search_paths = [
            "genomeamr_config.json",
            "priority1_fixes_output/genomeamr_config.json",
            "config/genomeamr_config.json",
            os.path.expanduser("~/.genomeamr/config.json")
        ]
        
        for path in search_paths:
            if os.path.exists(path):
                return path
        
        # Return default path if none found
        return "genomeamr_config.json"
    
    def _load_configuration(self) -> Dict[str, Any]:
        """Load configuration from file"""
        try:
            with open(self.config_file, 'r', encoding='utf-8') as f:
                return json.load(f)
        except FileNotFoundError:
            # Return minimal default configuration
            return {
                "default_genes": {
                    "rnd_efflux_pumps": {
                        "primary": ["acrA", "acrB", "acrE"]
                    }
                },
                "reference_strains": {
                    "primary": {
                        "id": "MG1655", 
                        "species": "Escherichia coli"
                    }
                },
                "analysis_parameters": {
                    "sequence_processing": {
                        "max_length": 50000
                    }
                }
            }
    
    def get_default_genes(self, category: str = "rnd_efflux_pumps", level: str = "primary") -> List[str]:
        """Get default genes for specified category and level"""
        return self.config.get("default_genes", {}).get(category, {}).get(level, [])
    
    def get_reference_strain(self, strain_type: str = "primary") -> Dict[str, str]:
        """Get reference strain information"""
        if strain_type == "primary":
            return self.config.get("reference_strains", {}).get("primary", {})
        else:
            alternatives = self.config.get("reference_strains", {}).get("alternatives", [])
            for strain in alternatives:
                if strain.get("id") == strain_type:
                    return strain
            return {}
    
    def get_analysis_parameter(self, category: str, parameter: str, default=None):
        """Get analysis parameter"""
        return self.config.get("analysis_parameters", {}).get(category, {}).get(parameter, default)

# Global configuration manager instance
config_manager = ConfigurationManager()

def get_config() -> ConfigurationManager:
    """Get global configuration manager instance"""
    return config_manager

def replace_hardcoded_genes(default_replacement: Optional[List[str]] = None) -> List[str]:
    """Replace hardcoded gene references with configuration"""
    if default_replacement:
        return default_replacement
    return config_manager.get_default_genes()

def replace_hardcoded_strain(default_replacement: Optional[str] = None) -> str:
    """Replace hardcoded strain references with configuration"""
    if default_replacement:
        return default_replacement
    strain_info = config_manager.get_reference_strain()
    return strain_info.get("id", "MG1655")
