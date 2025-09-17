"""
Pipeline Configuration Management
--------------------------------
Centralized configuration system for AMR genome analysis pipeline.
Handles database paths, processing parameters, and workflow settings with validation.
"""

import os
import json
import yaml
from pathlib import Path
from typing import Dict, Any, Optional, List, Union
from dataclasses import dataclass, field, asdict
from enum import Enum
import logging

class AlignmentPreset(Enum):
    """Minimap2 alignment presets."""
    MAP_ONT = "map-ont"  # Oxford Nanopore
    MAP_PB = "map-pb"    # PacBio
    MAP_HiFi = "map-hifi" # PacBio HiFi
    SR = "sr"            # Short reads
    SPLICE = "splice"    # Spliced alignment

class OutputFormat(Enum):
    """Supported output formats."""
    PAF = "paf"
    SAM = "sam"
    JSON = "json"
    CSV = "csv"

@dataclass
class DatabaseConfig:
    """Configuration for reference databases."""
    amr_database_path: str
    metadata_path: Optional[str] = None
    index_path: Optional[str] = None
    database_version: Optional[str] = None
    last_updated: Optional[str] = None
    
    def __post_init__(self):
        if not os.path.exists(self.amr_database_path):
            raise ValueError(f"AMR database not found: {self.amr_database_path}")

@dataclass
class ProcessingConfig:
    """Configuration for sequence processing."""
    threads: int = 4
    preset: AlignmentPreset = AlignmentPreset.MAP_ONT
    min_identity: float = 80.0
    min_coverage: float = 50.0
    quality_threshold: int = 20
    max_memory_gb: int = 8
    chunk_size: int = 1000
    enable_parallel: bool = True
    max_parallel_samples: int = 2
    
    def __post_init__(self):
        if self.threads < 1:
            self.threads = 1
        if not 0 <= self.min_identity <= 100:
            raise ValueError("min_identity must be between 0 and 100")
        if not 0 <= self.min_coverage <= 100:
            raise ValueError("min_coverage must be between 0 and 100")

@dataclass
class OutputConfig:
    """Configuration for output generation."""
    base_output_dir: str
    formats: List[OutputFormat] = field(default_factory=lambda: [OutputFormat.PAF, OutputFormat.JSON])
    save_intermediate: bool = True
    compress_outputs: bool = False
    generate_reports: bool = True
    report_format: str = "html"
    
    def __post_init__(self):
        os.makedirs(self.base_output_dir, exist_ok=True)

@dataclass
class LoggingConfig:
    """Configuration for logging."""
    level: str = "INFO"
    log_file: Optional[str] = None
    console_output: bool = True
    detailed_logs: bool = False
    
    def get_log_level(self) -> int:
        """Convert string level to logging constant."""
        return getattr(logging, self.level.upper(), logging.INFO)

@dataclass
class PipelineConfig:
    """Main pipeline configuration."""
    database: DatabaseConfig
    processing: ProcessingConfig = field(default_factory=ProcessingConfig)
    output: OutputConfig = field(default_factory=lambda: OutputConfig("./output"))
    logging: LoggingConfig = field(default_factory=LoggingConfig)
    pipeline_version: str = "2.0.0"
    config_name: str = "default"
    
    @classmethod
    def from_file(cls, config_path: str) -> 'PipelineConfig':
        """Load configuration from YAML or JSON file."""
        path_obj = Path(config_path)
        
        if not path_obj.exists():
            raise FileNotFoundError(f"Configuration file not found: {config_path}")
        
        with open(path_obj, 'r') as f:
            if path_obj.suffix.lower() in ['.yml', '.yaml']:
                data = yaml.safe_load(f)
            elif path_obj.suffix.lower() == '.json':
                data = json.load(f)
            else:
                raise ValueError(f"Unsupported config format: {path_obj.suffix}")
        
        return cls.from_dict(data)
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'PipelineConfig':
        """Create configuration from dictionary."""
        # Extract nested configurations
        db_config = DatabaseConfig(**data.get('database', {}))
        
        proc_data = data.get('processing', {})
        if 'preset' in proc_data:
            proc_data['preset'] = AlignmentPreset(proc_data['preset'])
        proc_config = ProcessingConfig(**proc_data)
        
        out_data = data.get('output', {})
        if 'formats' in out_data:
            out_data['formats'] = [OutputFormat(fmt) for fmt in out_data['formats']]
        out_config = OutputConfig(**out_data)
        
        log_config = LoggingConfig(**data.get('logging', {}))
        
        return cls(
            database=db_config,
            processing=proc_config,
            output=out_config,
            logging=log_config,
            pipeline_version=data.get('pipeline_version', '2.0.0'),
            config_name=data.get('config_name', 'default')
        )
    
    def to_file(self, config_path: str, format: str = 'yaml'):
        """Save configuration to file."""
        path_obj = Path(config_path)
        data = self.to_dict()
        
        with open(path_obj, 'w') as f:
            if format.lower() in ['yml', 'yaml']:
                yaml.dump(data, f, default_flow_style=False, indent=2)
            elif format.lower() == 'json':
                json.dump(data, f, indent=2)
            else:
                raise ValueError(f"Unsupported format: {format}")
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert configuration to dictionary."""
        data = asdict(self)
        
        # Convert enums to strings
        data['processing']['preset'] = self.processing.preset.value
        data['output']['formats'] = [fmt.value for fmt in self.output.formats]
        
        return data
    
    def validate(self) -> List[str]:
        """Validate configuration and return list of issues."""
        issues = []
        
        # Check database accessibility
        if not os.path.exists(self.database.amr_database_path):
            issues.append(f"AMR database not accessible: {self.database.amr_database_path}")
        
        # Check output directory writability
        try:
            test_file = Path(self.output.base_output_dir) / "test_write.tmp"
            test_file.touch()
            test_file.unlink()
        except Exception:
            issues.append(f"Output directory not writable: {self.output.base_output_dir}")
        
        # Check resource limits
        cpu = os.cpu_count() or 1
        if self.processing.threads > cpu:
            issues.append(f"Thread count ({self.processing.threads}) exceeds CPU count ({cpu})")
        
        return issues

class ConfigManager:
    """Manages pipeline configurations with validation and presets."""
    
    def __init__(self, config_dir: Optional[str] = None):
        self.config_dir = Path(config_dir) if config_dir else Path.home() / ".amr_analyzer"
        self.config_dir.mkdir(exist_ok=True)
        self.logger = self._setup_logger()
    
    def _setup_logger(self) -> logging.Logger:
        logger = logging.getLogger("ConfigManager")
        if not logger.handlers:
            handler = logging.StreamHandler()
            formatter = logging.Formatter('%(asctime)s %(levelname)s %(name)s: %(message)s')
            handler.setFormatter(formatter)
            logger.addHandler(handler)
        logger.setLevel(logging.INFO)
        return logger
    
    def create_default_config(self, amr_db_path: str, output_dir: str) -> PipelineConfig:
        """Create a default configuration."""
        return PipelineConfig(
            database=DatabaseConfig(amr_database_path=amr_db_path),
            output=OutputConfig(base_output_dir=output_dir),
            config_name="default"
        )
    
    def create_high_throughput_config(self, amr_db_path: str, output_dir: str) -> PipelineConfig:
        """Create configuration optimized for high-throughput processing."""
        return PipelineConfig(
            database=DatabaseConfig(amr_database_path=amr_db_path),
            processing=ProcessingConfig(
                threads=min(16, (os.cpu_count() or 1)),
                enable_parallel=True,
                max_parallel_samples=4,
                chunk_size=5000
            ),
            output=OutputConfig(
                base_output_dir=output_dir,
                save_intermediate=False,
                compress_outputs=True
            ),
            config_name="high_throughput"
        )
    
    def create_sensitive_config(self, amr_db_path: str, output_dir: str) -> PipelineConfig:
        """Create configuration for sensitive detection."""
        return PipelineConfig(
            database=DatabaseConfig(amr_database_path=amr_db_path),
            processing=ProcessingConfig(
                min_identity=70.0,
                min_coverage=40.0,
                quality_threshold=10
            ),
            output=OutputConfig(
                base_output_dir=output_dir,
                generate_reports=True,
                report_format="html"
            ),
            config_name="sensitive"
        )
    
    def save_config(self, config: PipelineConfig, name: Optional[str] = None) -> str:
        """Save configuration to config directory."""
        name = name or config.config_name
        config_path = self.config_dir / f"{name}.yaml"
        config.to_file(str(config_path))
        self.logger.info(f"Configuration saved to: {config_path}")
        return str(config_path)
    
    def load_config(self, name: str) -> PipelineConfig:
        """Load configuration by name."""
        config_path = self.config_dir / f"{name}.yaml"
        if not config_path.exists():
            config_path = self.config_dir / f"{name}.json"
        
        if not config_path.exists():
            raise FileNotFoundError(f"Configuration '{name}' not found")
        
        return PipelineConfig.from_file(str(config_path))
    
    def list_configs(self) -> List[str]:
        """List available configurations."""
        configs = []
        for ext in ['*.yaml', '*.yml', '*.json']:
            configs.extend([p.stem for p in self.config_dir.glob(ext)])
        return sorted(set(configs))
    
    def validate_config(self, config: PipelineConfig) -> bool:
        """Validate configuration and log issues."""
        issues = config.validate()
        if issues:
            for issue in issues:
                self.logger.warning(f"Configuration issue: {issue}")
            return False
        self.logger.info("Configuration validation passed")
        return True
