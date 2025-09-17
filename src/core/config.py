"""
Configuration management system for GenomeAMRAnalyzer.
Supports YAML and JSON configuration files for maximum flexibility.
"""
import yaml
import json
from pathlib import Path
from typing import Dict, Any, List, Optional, Union
from dataclasses import dataclass, asdict, field
import logging

from .error_handling import ConfigurationError, ValidationError, ValidationSuite

@dataclass
class GeneConfig:
    """Configuration for a single gene analysis."""
    name: str
    reference_sequence: str
    description: Optional[str] = None
    organism: Optional[str] = None
    uniprot_id: Optional[str] = None
    ncbi_id: Optional[str] = None
    gene_family: Optional[str] = None
    
    def __post_init__(self):
        """Validate gene configuration after initialization."""
        # Validate gene name
        self.name = ValidationSuite.validate_gene_name(self.name)
        
        # Validate reference sequence
        if self.reference_sequence:
            self.reference_sequence = ValidationSuite.validate_sequence(
                self.reference_sequence, min_length=10, max_length=50000
            )
        else:
            raise ConfigurationError("Reference sequence is required", "reference_sequence")
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        return asdict(self)

@dataclass
class AnalysisConfig:
    """Configuration for analysis parameters."""
    alignment_method: str = "muscle"
    similarity_threshold: float = 0.8
    min_sequence_length: int = 10
    max_sequence_length: int = 50000
    gap_penalty: float = -1.0
    match_score: float = 2.0
    mismatch_penalty: float = -1.0
    
    def __post_init__(self):
        """Validate analysis parameters."""
        if not 0.0 <= self.similarity_threshold <= 1.0:
            raise ConfigurationError(
                f"Similarity threshold must be between 0 and 1, got {self.similarity_threshold}",
                "similarity_threshold"
            )
        
        if self.min_sequence_length < 1:
            raise ConfigurationError(
                f"Minimum sequence length must be positive, got {self.min_sequence_length}",
                "min_sequence_length"
            )
        
        if self.max_sequence_length < self.min_sequence_length:
            raise ConfigurationError(
                "Maximum sequence length must be >= minimum sequence length",
                "max_sequence_length"
            )

@dataclass
class OutputConfig:
    """Configuration for output settings."""
    output_directory: str = "results"
    output_formats: List[str] = field(default_factory=lambda: ["fasta", "csv", "json", "html"])
    save_alignments: bool = True
    save_mutations: bool = True
    save_statistics: bool = True
    create_visualizations: bool = True
    
    def __post_init__(self):
        """Validate output configuration."""
        # Validate output formats
        valid_formats = {"fasta", "csv", "json", "html", "excel", "tsv"}
        invalid_formats = set(self.output_formats) - valid_formats
        if invalid_formats:
            raise ConfigurationError(
                f"Invalid output formats: {invalid_formats}. Valid formats: {valid_formats}",
                "output_formats"
            )

@dataclass
class PerformanceConfig:
    """Configuration for performance settings."""
    max_workers: int = 4
    batch_size: int = 1000
    use_cache: bool = True
    cache_size_mb: int = 500
    memory_limit_gb: Optional[float] = None
    enable_parallel_processing: bool = True
    
    def __post_init__(self):
        """Validate performance configuration."""
        if self.max_workers < 1:
            raise ConfigurationError("max_workers must be positive", "max_workers")
        
        if self.batch_size < 1:
            raise ConfigurationError("batch_size must be positive", "batch_size")
        
        if self.cache_size_mb < 0:
            raise ConfigurationError("cache_size_mb must be non-negative", "cache_size_mb")

@dataclass
class LoggingConfig:
    """Configuration for logging settings."""
    log_level: str = "INFO"
    log_file: Optional[str] = None
    log_format: str = "detailed"
    enable_performance_logging: bool = True
    
    def __post_init__(self):
        """Validate logging configuration."""
        valid_levels = {"DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"}
        if self.log_level.upper() not in valid_levels:
            raise ConfigurationError(
                f"Invalid log level: {self.log_level}. Valid levels: {valid_levels}",
                "log_level"
            )
        
        valid_formats = {"simple", "detailed", "json"}
        if self.log_format not in valid_formats:
            raise ConfigurationError(
                f"Invalid log format: {self.log_format}. Valid formats: {valid_formats}",
                "log_format"
            )

@dataclass
class PipelineConfig:
    """Main pipeline configuration."""
    # Core configuration sections
    genes: List[GeneConfig] = field(default_factory=list)
    analysis: AnalysisConfig = field(default_factory=AnalysisConfig)
    output: OutputConfig = field(default_factory=OutputConfig)
    performance: PerformanceConfig = field(default_factory=PerformanceConfig)
    logging: LoggingConfig = field(default_factory=LoggingConfig)
    
    # Metadata
    project_name: str = "GenomeAMRAnalysis"
    version: str = "1.0.0"
    description: str = "Genome antimicrobial resistance analysis"
    author: str = ""
    
    def __post_init__(self):
        """Validate complete pipeline configuration."""
        if not self.genes:
            raise ConfigurationError("At least one gene configuration is required", "genes")
        
        # Validate project name
        self.project_name = ValidationSuite.validate_gene_name(self.project_name)
    
    @classmethod
    def from_file(cls, config_path: Union[str, Path]) -> 'PipelineConfig':
        """
        Load configuration from YAML or JSON file.
        
        Args:
            config_path: Path to configuration file
            
        Returns:
            PipelineConfig instance
            
        Raises:
            ConfigurationError: If file cannot be loaded or parsed
        """
        config_path = Path(config_path)
        
        if not config_path.exists():
            raise ConfigurationError(f"Configuration file not found: {config_path}")
        
        try:
            with open(config_path, 'r', encoding='utf-8') as f:
                if config_path.suffix.lower() in ['.yaml', '.yml']:
                    data = yaml.safe_load(f)
                elif config_path.suffix.lower() == '.json':
                    data = json.load(f)
                else:
                    raise ConfigurationError(
                        f"Unsupported config file format: {config_path.suffix}. Use .yaml, .yml, or .json"
                    )
        except yaml.YAMLError as e:
            raise ConfigurationError(f"YAML parsing error: {e}")
        except json.JSONDecodeError as e:
            raise ConfigurationError(f"JSON parsing error: {e}")
        except Exception as e:
            raise ConfigurationError(f"Error reading config file: {e}")
        
        return cls.from_dict(data)
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'PipelineConfig':
        """
        Create configuration from dictionary.
        
        Args:
            data: Configuration dictionary
            
        Returns:
            PipelineConfig instance
        """
        try:
            # Parse gene configurations
            genes = []
            for gene_data in data.get('genes', []):
                genes.append(GeneConfig(**gene_data))
            
            # Parse other sections with defaults
            analysis = AnalysisConfig(**data.get('analysis', {}))
            output = OutputConfig(**data.get('output', {}))
            performance = PerformanceConfig(**data.get('performance', {}))
            logging_config = LoggingConfig(**data.get('logging', {}))
            
            # Create main config
            config_data = {
                'genes': genes,
                'analysis': analysis,
                'output': output,
                'performance': performance,
                'logging': logging_config,
                'project_name': data.get('project_name', 'GenomeAMRAnalysis'),
                'version': data.get('version', '1.0.0'),
                'description': data.get('description', 'Genome antimicrobial resistance analysis'),
                'author': data.get('author', '')
            }
            
            return cls(**config_data)
            
        except TypeError as e:
            raise ConfigurationError(f"Invalid configuration structure: {e}")
        except Exception as e:
            raise ConfigurationError(f"Error creating configuration: {e}")
    
    def to_file(self, config_path: Union[str, Path], format_type: str = "yaml") -> None:
        """
        Save configuration to file.
        
        Args:
            config_path: Path where to save configuration
            format_type: File format ('yaml' or 'json')
        """
        config_path = Path(config_path)
        config_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Convert to dictionary
        data = self.to_dict()
        
        try:
            with open(config_path, 'w', encoding='utf-8') as f:
                if format_type.lower() == 'yaml':
                    yaml.dump(data, f, default_flow_style=False, sort_keys=False, indent=2)
                elif format_type.lower() == 'json':
                    json.dump(data, f, indent=2, ensure_ascii=False)
                else:
                    raise ConfigurationError(f"Unsupported format: {format_type}")
        except Exception as e:
            raise ConfigurationError(f"Error saving configuration: {e}")
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        return {
            'project_name': self.project_name,
            'version': self.version,
            'description': self.description,
            'author': self.author,
            'genes': [gene.to_dict() for gene in self.genes],
            'analysis': asdict(self.analysis),
            'output': asdict(self.output),
            'performance': asdict(self.performance),
            'logging': asdict(self.logging)
        }
    
    def validate(self) -> List[str]:
        """
        Perform comprehensive validation and return any warnings.
        
        Returns:
            List of warning messages
        """
        warnings = []
        
        # Check for duplicate gene names
        gene_names = [gene.name for gene in self.genes]
        if len(gene_names) != len(set(gene_names)):
            warnings.append("Duplicate gene names detected")
        
        # Check output directory permissions
        try:
            output_path = Path(self.output.output_directory)
            output_path.mkdir(parents=True, exist_ok=True)
            if not output_path.exists() or not os.access(output_path, os.W_OK):
                warnings.append(f"Output directory not writable: {output_path}")
        except Exception:
            warnings.append("Cannot create output directory")
        
        # Performance warnings
        if self.performance.max_workers > 16:
            warnings.append("High worker count may cause memory issues")
        
        if self.performance.batch_size > 10000:
            warnings.append("Large batch size may cause memory issues")
        
        return warnings

def create_example_config() -> PipelineConfig:
    """
    Create an example configuration for demonstration.
    
    Returns:
        Example PipelineConfig
    """
    genes = [
        GeneConfig(
            name="example_gene_1",
            reference_sequence="ATGAAAAAACTGATTGCCCTGCTGAACGTGCTGCCCAAAGATAACGCCCGCGCGAACCTGTGA",
            description="Example gene for resistance analysis",
            organism="Escherichia coli",
            gene_family="resistance"
        ),
        GeneConfig(
            name="example_gene_2", 
            reference_sequence="ATGCCCAAACTGATTGCCCTGCTGAACGTGCTGCCCAAAGATAACGCCCGCGCGAACCTGTGA",
            description="Another example gene",
            organism="Escherichia coli",
            gene_family="resistance"
        )
    ]
    
    return PipelineConfig(
        genes=genes,
        project_name="Example_AMR_Analysis",
        description="Example antimicrobial resistance analysis project",
        author="Research Team"
    )

def save_example_config(config_path: Union[str, Path] = "config/example_config.yaml") -> None:
    """
    Save an example configuration file.
    
    Args:
        config_path: Where to save the example configuration
    """
    config = create_example_config()
    config.to_file(config_path, "yaml")
    print(f"Example configuration saved to: {config_path}")

if __name__ == "__main__":
    # Create example configuration
    save_example_config()