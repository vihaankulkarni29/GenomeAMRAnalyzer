#!/usr/bin/env python3
"""
SEPI 2.0 Configuration Manager - Production Grade
Robust configuration system for WildType Aligner integration

This module provides:
1. Centralized SEPI 2.0 configuration management
2. Environment-specific settings and overrides
3. Validation and error handling for all configurations
4. Dynamic configuration updates and caching
5. Integration with ProductionWildTypeAligner

Author: GenomeAMRAnalyzer Pipeline
Version: 2.0 - Production Configuration System
"""

import os
import sys
import yaml
import json
import logging
import hashlib
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Any, Union
from dataclasses import dataclass, asdict, field
from datetime import datetime
import subprocess


@dataclass
class SEPIEnvironmentConfig:
    """SEPI environment and system configuration"""
    sepi_script_path: str
    cache_directory: str
    temp_directory: str
    max_cache_size_mb: int = 1000
    cache_retention_days: int = 30
    concurrent_downloads: int = 3
    timeout_seconds: int = 300
    retry_attempts: int = 3
    
    def __post_init__(self):
        """Validate and create directories"""
        # Convert to Path objects for internal use, but keep as strings for serialization
        sepi_path = Path(self.sepi_script_path).resolve()
        cache_path = Path(self.cache_directory).resolve()
        temp_path = Path(self.temp_directory).resolve()
        
        # Update with resolved paths
        self.sepi_script_path = str(sepi_path)
        self.cache_directory = str(cache_path)
        self.temp_directory = str(temp_path)
        
        # Create directories if they don't exist
        cache_path.mkdir(parents=True, exist_ok=True)
        temp_path.mkdir(parents=True, exist_ok=True)
    
    @property
    def sepi_script_path_obj(self) -> Path:
        """Get sepi script path as Path object"""
        return Path(self.sepi_script_path)
    
    @property
    def cache_directory_obj(self) -> Path:
        """Get cache directory as Path object"""
        return Path(self.cache_directory)
    
    @property
    def temp_directory_obj(self) -> Path:
        """Get temp directory as Path object"""
        return Path(self.temp_directory)


@dataclass
class NCBIAccessConfig:
    """NCBI API access configuration"""
    email: str
    api_key: Optional[str] = None
    max_requests_per_second: float = 3.0
    batch_size: int = 200
    entrez_tool: str = "GenomeAMRAnalyzer"
    entrez_email: str = ""
    
    def __post_init__(self):
        """Set default email for Entrez"""
        if not self.entrez_email:
            self.entrez_email = self.email


@dataclass
class OrganismPreferences:
    """Organism-specific configuration preferences"""
    organism_name: str
    assembly_level: str = "complete_genome"
    strain_preferences: List[str] = field(default_factory=list)
    exclude_strains: List[str] = field(default_factory=list)
    max_assemblies: int = 5
    quality_filters: Dict[str, Any] = field(default_factory=dict)


@dataclass
class GeneTargetConfig:
    """Gene-specific target configuration"""
    gene_name: str
    organism_preferences: List[str] = field(default_factory=list)
    alternative_names: List[str] = field(default_factory=list)
    quality_threshold: float = 0.8
    min_sequence_length: int = 50
    max_sequence_length: int = 5000
    required_annotations: List[str] = field(default_factory=list)


@dataclass
class SEPIQualityConfig:
    """Quality control and validation configuration"""
    min_sequence_length: int = 50
    max_sequence_length: int = 5000
    require_annotation: bool = True
    validate_protein_sequence: bool = True
    exclude_partial_sequences: bool = True
    min_assembly_quality: str = "scaffold"
    checksum_validation: bool = True


class SEPIConfigurationManager:
    """
    Production-grade SEPI 2.0 configuration manager
    """
    
    def __init__(self, config_file: Optional[str] = None, 
                 workspace_root: Optional[str] = None):
        """Initialize configuration manager"""
        
        self.workspace_root = Path(workspace_root) if workspace_root else Path.cwd()
        self.config_file = config_file
        
        # Setup logging
        self.logger = logging.getLogger('SEPIConfigurationManager')
        
        # Configuration storage
        self.environment_config: Optional[SEPIEnvironmentConfig] = None
        self.ncbi_config: Optional[NCBIAccessConfig] = None
        self.quality_config: Optional[SEPIQualityConfig] = None
        self.organism_configs: Dict[str, OrganismPreferences] = {}
        self.gene_configs: Dict[str, GeneTargetConfig] = {}
        
        # Load configuration
        self._initialize_configuration()
        
        # Validate configuration
        self._validate_configuration()
    
    def _initialize_configuration(self):
        """Initialize configuration from multiple sources"""
        
        # Load base configuration
        if self.config_file and Path(self.config_file).exists():
            self._load_from_file(self.config_file)
        else:
            self._load_default_configuration()
        
        # Apply environment overrides
        self._apply_environment_overrides()
        
        # Load gene-specific configurations
        self._load_gene_configurations()
        
        # Load organism-specific configurations
        self._load_organism_configurations()
    
    def _load_default_configuration(self):
        """Load robust default configuration"""
        
        # Environment configuration
        sepi_path = self._find_sepi_script()
        cache_dir = self.workspace_root / ".sepi_cache"
        temp_dir = self.workspace_root / "temp" / "sepi"
        
        self.environment_config = SEPIEnvironmentConfig(
            sepi_script_path=str(sepi_path),
            cache_directory=str(cache_dir),
            temp_directory=str(temp_dir),
            max_cache_size_mb=2000,
            cache_retention_days=60,
            concurrent_downloads=5,
            timeout_seconds=600,
            retry_attempts=3
        )
        
        # NCBI configuration with workspace email lookup
        default_email = self._get_workspace_email()
        self.ncbi_config = NCBIAccessConfig(
            email=default_email,
            api_key=self._get_workspace_api_key(),
            max_requests_per_second=3.0,
            batch_size=100,
            entrez_tool="GenomeAMRAnalyzer-SEPI2.0"
        )
        
        # Quality configuration
        self.quality_config = SEPIQualityConfig(
            min_sequence_length=30,
            max_sequence_length=10000,
            require_annotation=True,
            validate_protein_sequence=True,
            exclude_partial_sequences=True,
            min_assembly_quality="complete_genome",
            checksum_validation=True
        )
        
        self.logger.info("Loaded default SEPI configuration")
    
    def _find_sepi_script(self) -> Path:
        """Intelligently locate SEPI 2.0 script"""
        
        # Search locations in priority order
        search_paths = [
            self.workspace_root / "MetaDataHarvester" / "sepi2.0" / "sepi.py",
            self.workspace_root / "sepi2.0" / "sepi.py", 
            Path("MetaDataHarvester/sepi2.0/sepi.py"),
            Path("sepi2.0/sepi.py"),
            Path("sepi.py")
        ]
        
        for path in search_paths:
            if path.exists() and path.is_file():
                self.logger.info(f"Found SEPI script at: {path}")
                return path
        
        # Search in PATH
        try:
            result = subprocess.run(["where", "sepi.py"], 
                                  capture_output=True, text=True, shell=True)
            if result.returncode == 0:
                sepi_path = Path(result.stdout.strip().split('\n')[0])
                if sepi_path.exists():
                    self.logger.info(f"Found SEPI in PATH: {sepi_path}")
                    return sepi_path
        except Exception:
            pass
        
        # Default fallback
        default_path = self.workspace_root / "MetaDataHarvester" / "sepi2.0" / "sepi.py"
        self.logger.warning(f"SEPI script not found, using default: {default_path}")
        return default_path
    
    def _get_workspace_email(self) -> str:
        """Get email from workspace configuration"""
        
        config_files = [
            self.workspace_root / "config" / "snakemake_config.yaml",
            self.workspace_root / "config.yaml",
            self.workspace_root / ".sepi_config.yaml"
        ]
        
        for config_file in config_files:
            if config_file.exists():
                try:
                    with open(config_file, 'r') as f:
                        config = yaml.safe_load(f)
                    
                    # Check multiple possible email keys
                    email_keys = ['ncbi_email', 'email', 'user_email', 'entrez_email']
                    for key in email_keys:
                        if key in config and config[key]:
                            return config[key]
                            
                except Exception as e:
                    self.logger.warning(f"Failed to read email from {config_file}: {e}")
        
        # Environment variable fallback
        env_email = os.environ.get('NCBI_EMAIL') or os.environ.get('EMAIL')
        if env_email:
            return env_email
        
        # Default placeholder
        return "researcher@institution.edu"
    
    def _get_workspace_api_key(self) -> Optional[str]:
        """Get NCBI API key from workspace configuration"""
        
        config_files = [
            self.workspace_root / "config" / "snakemake_config.yaml",
            self.workspace_root / "config.yaml"
        ]
        
        for config_file in config_files:
            if config_file.exists():
                try:
                    with open(config_file, 'r') as f:
                        config = yaml.safe_load(f)
                    
                    api_keys = ['ncbi_api_key', 'api_key', 'entrez_api_key']
                    for key in api_keys:
                        if key in config and config[key]:
                            return config[key]
                            
                except Exception as e:
                    self.logger.warning(f"Failed to read API key from {config_file}: {e}")
        
        # Environment variable fallback
        return os.environ.get('NCBI_API_KEY')
    
    def _load_gene_configurations(self):
        """Load gene-specific configurations with AMR focus"""
        
        # AMR-specific gene configurations (senior bioinformatician knowledge)
        amr_genes = {
            # E. coli AcrAB-TolC system
            'acrA': GeneTargetConfig(
                gene_name='acrA',
                organism_preferences=['Escherichia coli K-12 MG1655', 'Escherichia coli', 'Escherichia'],
                alternative_names=['acra', 'AcrA', 'membrane fusion protein AcrA'],
                quality_threshold=0.95,
                min_sequence_length=300,
                max_sequence_length=450,
                required_annotations=['membrane fusion protein', 'efflux', 'AcrA']
            ),
            'acrB': GeneTargetConfig(
                gene_name='acrB',
                organism_preferences=['Escherichia coli K-12 MG1655', 'Escherichia coli', 'Escherichia'],
                alternative_names=['acrb', 'AcrB', 'multidrug efflux pump AcrB'],
                quality_threshold=0.95,
                min_sequence_length=1000,
                max_sequence_length=1200,
                required_annotations=['multidrug efflux', 'RND', 'AcrB']
            ),
            'tolC': GeneTargetConfig(
                gene_name='tolC',
                organism_preferences=['Escherichia coli K-12 MG1655', 'Escherichia coli', 'Escherichia'],
                alternative_names=['tolc', 'TolC', 'outer membrane channel TolC'],
                quality_threshold=0.95,
                min_sequence_length=400,
                max_sequence_length=500,
                required_annotations=['outer membrane', 'channel', 'TolC']
            ),
            
            # Pseudomonas MexAB-OprM system
            'mexA': GeneTargetConfig(
                gene_name='mexA',
                organism_preferences=['Pseudomonas aeruginosa PAO1', 'Pseudomonas aeruginosa', 'Pseudomonas'],
                alternative_names=['mexa', 'MexA', 'membrane fusion protein MexA'],
                quality_threshold=0.95,
                min_sequence_length=350,
                max_sequence_length=450,
                required_annotations=['membrane fusion protein', 'efflux', 'MexA']
            ),
            'mexB': GeneTargetConfig(
                gene_name='mexB',
                organism_preferences=['Pseudomonas aeruginosa PAO1', 'Pseudomonas aeruginosa', 'Pseudomonas'],
                alternative_names=['mexb', 'MexB', 'multidrug efflux pump MexB'],
                quality_threshold=0.95,
                min_sequence_length=1000,
                max_sequence_length=1200,
                required_annotations=['multidrug efflux', 'RND', 'MexB']
            ),
            'oprM': GeneTargetConfig(
                gene_name='oprM',
                organism_preferences=['Pseudomonas aeruginosa PAO1', 'Pseudomonas aeruginosa', 'Pseudomonas'],
                alternative_names=['oprm', 'OprM', 'outer membrane protein OprM'],
                quality_threshold=0.95,
                min_sequence_length=450,
                max_sequence_length=550,
                required_annotations=['outer membrane', 'efflux', 'OprM']
            ),
            
            # Additional AMR genes
            'acrR': GeneTargetConfig(
                gene_name='acrR',
                organism_preferences=['Escherichia coli K-12 MG1655', 'Escherichia coli'],
                alternative_names=['acrr', 'AcrR', 'AcrAB operon repressor'],
                quality_threshold=0.90,
                min_sequence_length=150,
                max_sequence_length=250,
                required_annotations=['repressor', 'transcriptional regulator', 'AcrR']
            ),
            'marA': GeneTargetConfig(
                gene_name='marA',
                organism_preferences=['Escherichia coli K-12 MG1655', 'Escherichia coli'],
                alternative_names=['mara', 'MarA', 'multiple antibiotic resistance protein MarA'],
                quality_threshold=0.90,
                min_sequence_length=120,
                max_sequence_length=180,
                required_annotations=['transcriptional activator', 'mar', 'MarA']
            )
        }
        
        self.gene_configs.update(amr_genes)
        self.logger.info(f"Loaded {len(amr_genes)} AMR gene configurations")
    
    def _load_organism_configurations(self):
        """Load organism-specific preferences"""
        
        organisms = {
            'escherichia_coli': OrganismPreferences(
                organism_name='Escherichia coli',
                assembly_level='complete_genome',
                strain_preferences=['K-12 MG1655', 'K-12', 'MG1655', 'str. K-12'],
                exclude_strains=['O157:H7', 'STEC'],
                max_assemblies=3,
                quality_filters={
                    'min_genome_size': 4000000,
                    'max_genome_size': 6000000,
                    'min_n50': 1000000
                }
            ),
            'pseudomonas_aeruginosa': OrganismPreferences(
                organism_name='Pseudomonas aeruginosa',
                assembly_level='complete_genome',
                strain_preferences=['PAO1', 'PA14', 'LESB58'],
                exclude_strains=[],
                max_assemblies=3,
                quality_filters={
                    'min_genome_size': 5500000,
                    'max_genome_size': 7500000,
                    'min_n50': 2000000
                }
            ),
            'klebsiella_pneumoniae': OrganismPreferences(
                organism_name='Klebsiella pneumoniae',
                assembly_level='complete_genome',
                strain_preferences=['ATCC 700721', 'MGH 78578'],
                exclude_strains=[],
                max_assemblies=5,
                quality_filters={
                    'min_genome_size': 5000000,
                    'max_genome_size': 6500000
                }
            ),
            'acinetobacter_baumannii': OrganismPreferences(
                organism_name='Acinetobacter baumannii',
                assembly_level='complete_genome',
                strain_preferences=['ATCC 17978', 'AYE'],
                exclude_strains=[],
                max_assemblies=5,
                quality_filters={
                    'min_genome_size': 3500000,
                    'max_genome_size': 4500000
                }
            )
        }
        
        self.organism_configs.update(organisms)
        self.logger.info(f"Loaded {len(organisms)} organism configurations")
    
    def _apply_environment_overrides(self):
        """Apply environment variable overrides"""
        
        # Environment overrides for sensitive data
        env_overrides = {
            'SEPI_EMAIL': lambda: setattr(self.ncbi_config, 'email', os.environ['SEPI_EMAIL']),
            'SEPI_API_KEY': lambda: setattr(self.ncbi_config, 'api_key', os.environ['SEPI_API_KEY']),
            'SEPI_CACHE_DIR': lambda: setattr(self.environment_config, 'cache_directory', os.environ['SEPI_CACHE_DIR']),
            'SEPI_TEMP_DIR': lambda: setattr(self.environment_config, 'temp_directory', os.environ['SEPI_TEMP_DIR']),
            'SEPI_MAX_CONCURRENT': lambda: setattr(self.environment_config, 'concurrent_downloads', int(os.environ['SEPI_MAX_CONCURRENT']))
        }
        
        for env_var, override_func in env_overrides.items():
            if env_var in os.environ:
                try:
                    override_func()
                    self.logger.info(f"Applied environment override: {env_var}")
                except Exception as e:
                    self.logger.warning(f"Failed to apply environment override {env_var}: {e}")
    
    def _validate_configuration(self):
        """Comprehensive configuration validation"""
        
        errors = []
        warnings = []
        
        # Ensure all configs are initialized
        if not self.environment_config:
            errors.append("Environment configuration not initialized")
            return
        
        if not self.ncbi_config:
            errors.append("NCBI configuration not initialized")
            return
            
        if not self.quality_config:
            errors.append("Quality configuration not initialized")
            return
        
        # Validate environment configuration
        if not self.environment_config.sepi_script_path_obj.exists():
            errors.append(f"SEPI script not found: {self.environment_config.sepi_script_path}")
        
        # Validate NCBI configuration
        if not self.ncbi_config.email or '@' not in self.ncbi_config.email:
            errors.append("Valid email address required for NCBI access")
        
        # Test SEPI availability
        try:
            result = subprocess.run(
                [sys.executable, self.environment_config.sepi_script_path, "--help"],
                capture_output=True, text=True, timeout=10
            )
            if result.returncode != 0:
                errors.append("SEPI script is not functional")
        except Exception as e:
            errors.append(f"Cannot execute SEPI script: {e}")
        
        # Validate cache directory permissions
        try:
            test_file = self.environment_config.cache_directory_obj / ".test_write"
            test_file.write_text("test")
            test_file.unlink()
        except Exception:
            errors.append(f"Cache directory not writable: {self.environment_config.cache_directory}")
        
        # Report validation results
        if errors:
            error_msg = "Configuration validation failed:\n" + "\n".join(f"  - {e}" for e in errors)
            self.logger.error(error_msg)
            raise ValueError(error_msg)
        
        if warnings:
            warning_msg = "Configuration warnings:\n" + "\n".join(f"  - {w}" for w in warnings)
            self.logger.warning(warning_msg)
        
        self.logger.info("Configuration validation successful")
    
    def get_gene_config(self, gene_name: str) -> GeneTargetConfig:
        """Get configuration for specific gene"""
        gene_key = gene_name.lower()
        
        if gene_key in self.gene_configs:
            return self.gene_configs[gene_key]
        
        # Return default configuration for unknown genes
        return GeneTargetConfig(
            gene_name=gene_name,
            organism_preferences=[],
            quality_threshold=0.8,
            min_sequence_length=50,
            max_sequence_length=5000
        )
    
    def get_organism_config(self, organism_name: str) -> OrganismPreferences:
        """Get configuration for specific organism"""
        organism_key = organism_name.lower().replace(' ', '_')
        
        if organism_key in self.organism_configs:
            return self.organism_configs[organism_key]
        
        # Return default configuration for unknown organisms
        return OrganismPreferences(
            organism_name=organism_name,
            assembly_level='complete_genome',
            max_assemblies=5
        )
    
    def generate_sepi_config_file(self, gene_name: str, organism: str, 
                                 output_file: Optional[str] = None) -> str:
        """Generate SEPI-compatible YAML configuration file"""
        
        if not self.ncbi_config or not self.environment_config:
            raise ValueError("Configuration not properly initialized")
        
        gene_config = self.get_gene_config(gene_name)
        organism_config = self.get_organism_config(organism)
        
        # Build SEPI configuration
        sepi_config = {
            'organism': organism_config.organism_name,
            'proteins': [gene_name] + gene_config.alternative_names,
            'output_name': f"{gene_name}_{organism.replace(' ', '_')}",
            'settings': {
                'assembly_level': organism_config.assembly_level,
                'multi_fasta': True,
                'html_report': True,
                'max_assemblies': organism_config.max_assemblies,
                'quality_filters': organism_config.quality_filters
            },
            'user_email': self.ncbi_config.email
        }
        
        # Add API key if available
        if self.ncbi_config.api_key:
            sepi_config['api_key'] = self.ncbi_config.api_key
        
        # Generate output file path
        if output_file is None:
            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
            output_file = str(self.environment_config.temp_directory_obj / f"sepi_config_{gene_name}_{timestamp}.yaml")
        
        # Write configuration file
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_path, 'w') as f:
            yaml.dump(sepi_config, f, default_flow_style=False, indent=2)
        
        self.logger.info(f"Generated SEPI config: {output_path}")
        return str(output_path)
    
    def save_configuration(self, output_file: str):
        """Save complete configuration to file"""
        
        if not all([self.environment_config, self.ncbi_config, self.quality_config]):
            raise ValueError("Configuration not properly initialized")
        
        # Type assertions for mypy
        assert self.environment_config is not None
        assert self.ncbi_config is not None 
        assert self.quality_config is not None
        
        config_data = {
            'sepi_configuration_version': '2.0',
            'generated_timestamp': datetime.now().isoformat(),
            'environment': asdict(self.environment_config),
            'ncbi_access': asdict(self.ncbi_config),
            'quality_control': asdict(self.quality_config),
            'organisms': {k: asdict(v) for k, v in self.organism_configs.items()},
            'genes': {k: asdict(v) for k, v in self.gene_configs.items()}
        }
        
        Path(output_file).parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_file, 'w') as f:
            yaml.dump(config_data, f, default_flow_style=False, indent=2)
        
        self.logger.info(f"Configuration saved to: {output_file}")
    
    def _load_from_file(self, config_file: str):
        """Load configuration from YAML file"""
        
        try:
            with open(config_file, 'r') as f:
                config_data = yaml.safe_load(f)
            
            # Load environment config
            if 'environment' in config_data:
                self.environment_config = SEPIEnvironmentConfig(**config_data['environment'])
            
            # Load NCBI config
            if 'ncbi_access' in config_data:
                self.ncbi_config = NCBIAccessConfig(**config_data['ncbi_access'])
            
            # Load quality config
            if 'quality_control' in config_data:
                self.quality_config = SEPIQualityConfig(**config_data['quality_control'])
            
            # Load organism configs
            if 'organisms' in config_data:
                for k, v in config_data['organisms'].items():
                    self.organism_configs[k] = OrganismPreferences(**v)
            
            # Load gene configs
            if 'genes' in config_data:
                for k, v in config_data['genes'].items():
                    self.gene_configs[k] = GeneTargetConfig(**v)
            
            self.logger.info(f"Configuration loaded from: {config_file}")
            
        except Exception as e:
            self.logger.error(f"Failed to load configuration from {config_file}: {e}")
            raise


def create_default_sepi_configuration(workspace_root: str, 
                                    config_output: Optional[str] = None) -> SEPIConfigurationManager:
    """
    Create and save default SEPI configuration for workspace
    """
    
    # Initialize configuration manager
    config_manager = SEPIConfigurationManager(workspace_root=workspace_root)
    
    # Save configuration
    if config_output is None:
        config_output = str(Path(workspace_root) / "config" / "sepi_configuration.yaml")
    
    config_manager.save_configuration(config_output)
    
    return config_manager


if __name__ == "__main__":
    """Command line interface for SEPI configuration management"""
    import argparse
    
    parser = argparse.ArgumentParser(description="SEPI 2.0 Configuration Manager")
    parser.add_argument("--workspace", required=True, help="Workspace root directory")
    parser.add_argument("--config-file", help="Configuration file to load")
    parser.add_argument("--save-config", help="Save configuration to file")
    parser.add_argument("--gene", help="Generate configuration for specific gene")
    parser.add_argument("--organism", help="Organism for gene configuration")
    parser.add_argument("--validate", action="store_true", help="Validate configuration only")
    
    args = parser.parse_args()
    
    # Setup logging
    logging.basicConfig(level=logging.INFO, 
                       format='%(asctime)s - %(levelname)s - %(message)s')
    
    try:
        # Initialize configuration manager
        config_manager = SEPIConfigurationManager(
            config_file=args.config_file,
            workspace_root=args.workspace
        )
        
        if args.validate:
            print("✅ Configuration validation successful")
        
        if args.save_config:
            config_manager.save_configuration(args.save_config)
            print(f"✅ Configuration saved to: {args.save_config}")
        
        if args.gene and args.organism:
            config_file = config_manager.generate_sepi_config_file(args.gene, args.organism)
            print(f"✅ SEPI configuration generated: {config_file}")
            
    except Exception as e:
        print(f"❌ Configuration failed: {e}")
        sys.exit(1)