"""
Robust error handling and validation system for GenomeAMRAnalyzer.
Provides comprehensive exception management, logging, and validation capabilities.
"""
import logging
import traceback
import functools
import sys
import os
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Union, Tuple
from datetime import datetime
import json


class PipelineError(Exception):
    """Base exception for all pipeline errors."""
    def __init__(self, message: str, error_code: str = None, context: Dict = None):
        super().__init__(message)
        self.message = message
        self.error_code = error_code or "GENERIC_ERROR"
        self.context = context or {}
        self.timestamp = datetime.now().isoformat()


class ValidationError(PipelineError):
    """Exception for data validation errors."""
    def __init__(self, message: str, field: str = None, value: Any = None):
        super().__init__(message, "VALIDATION_ERROR")
        self.field = field
        self.value = value


class ConfigurationError(PipelineError):
    """Exception for configuration errors."""
    def __init__(self, message: str, config_key: str = None):
        super().__init__(message, "CONFIG_ERROR")
        self.config_key = config_key


class DataProcessingError(PipelineError):
    """Exception for data processing errors."""
    def __init__(self, message: str, stage: str = None):
        super().__init__(message, "DATA_PROCESSING_ERROR")
        self.stage = stage


class FileSystemError(PipelineError):
    """Exception for file system related errors."""
    def __init__(self, message: str, file_path: str = None):
        super().__init__(message, "FILESYSTEM_ERROR")
        self.file_path = file_path


class RobustLogger:
    """Production-grade logging system with multiple fallbacks."""
    
    def __init__(self, name: str = "GenomeAMRAnalyzer", 
                 log_level: str = "INFO",
                 log_file: Optional[str] = None):
        
        self.name = name
        self.logger = logging.getLogger(name)
        self.logger.setLevel(getattr(logging, log_level.upper()))
        
        # Clear existing handlers to avoid duplicates
        self.logger.handlers.clear()
        
        # Set up formatters
        detailed_formatter = logging.Formatter(
            '%(asctime)s | %(name)s | %(levelname)s | %(filename)s:%(lineno)d | %(message)s'
        )
        simple_formatter = logging.Formatter(
            '%(asctime)s | %(levelname)s | %(message)s'
        )
        
        # Console handler with fallback
        try:
            console_handler = logging.StreamHandler(sys.stdout)
            console_handler.setFormatter(simple_formatter)
            console_handler.setLevel(logging.INFO)
            self.logger.addHandler(console_handler)
        except Exception:
            # Fallback to stderr if stdout fails
            try:
                console_handler = logging.StreamHandler(sys.stderr)
                console_handler.setFormatter(simple_formatter)
                self.logger.addHandler(console_handler)
            except Exception:
                pass  # Give up on console logging
        
        # File handler if specified
        if log_file:
            try:
                # Ensure log directory exists
                log_path = Path(log_file)
                log_path.parent.mkdir(parents=True, exist_ok=True)
                
                file_handler = logging.FileHandler(log_file)
                file_handler.setFormatter(detailed_formatter)
                file_handler.setLevel(logging.DEBUG)
                self.logger.addHandler(file_handler)
                
            except Exception as e:
                self.logger.warning(f"Could not set up file logging: {e}")
    
    def log_exception(self, exception: Exception, context: Dict = None):
        """Log exception with full context and traceback."""
        error_info = {
            'exception_type': type(exception).__name__,
            'exception_message': str(exception),
            'traceback': traceback.format_exc(),
            'context': context or {},
            'timestamp': datetime.now().isoformat()
        }
        
        self.logger.error(f"Exception occurred: {json.dumps(error_info, indent=2)}")


class ValidationSuite:
    """Comprehensive validation suite for all data types."""
    
    @staticmethod
    def validate_sequence(sequence: str, min_length: int = 1, max_length: int = 100000) -> str:
        """
        Validate and clean a biological sequence.
        
        Args:
            sequence: Input sequence string
            min_length: Minimum allowed length
            max_length: Maximum allowed length
            
        Returns:
            Cleaned sequence string
            
        Raises:
            ValidationError: If sequence is invalid
        """
        if not isinstance(sequence, str):
            raise ValidationError(f"Sequence must be string, got {type(sequence)}", "sequence", sequence)
        
        if not sequence or not sequence.strip():
            raise ValidationError("Sequence cannot be empty", "sequence", sequence)
        
        # Clean sequence
        cleaned_sequence = sequence.strip().upper()
        
        # Remove any whitespace
        cleaned_sequence = ''.join(cleaned_sequence.split())
        
        # Validate length
        if len(cleaned_sequence) < min_length:
            raise ValidationError(
                f"Sequence too short: {len(cleaned_sequence)} < {min_length}",
                "sequence_length", len(cleaned_sequence)
            )
        
        if len(cleaned_sequence) > max_length:
            raise ValidationError(
                f"Sequence too long: {len(cleaned_sequence)} > {max_length}",
                "sequence_length", len(cleaned_sequence)
            )
        
        # Validate characters (allow standard nucleotide codes + gaps)
        valid_chars = set('ATCGRYSWKMBDHVN-')
        invalid_chars = set(cleaned_sequence) - valid_chars
        if invalid_chars:
            raise ValidationError(
                f"Invalid characters in sequence: {invalid_chars}",
                "sequence_characters", invalid_chars
            )
        
        return cleaned_sequence
    
    @staticmethod
    def validate_gene_name(gene_name: str) -> str:
        """
        Validate and sanitize gene name.
        
        Args:
            gene_name: Input gene name
            
        Returns:
            Sanitized gene name
            
        Raises:
            ValidationError: If gene name is invalid
        """
        if not isinstance(gene_name, str):
            raise ValidationError(f"Gene name must be string, got {type(gene_name)}", "gene_name", gene_name)
        
        if not gene_name or not gene_name.strip():
            raise ValidationError("Gene name cannot be empty", "gene_name", gene_name)
        
        # Clean and sanitize
        cleaned_name = gene_name.strip()
        
        # Replace spaces with underscores
        cleaned_name = cleaned_name.replace(' ', '_')
        
        # Remove dangerous characters for file systems
        import re
        cleaned_name = re.sub(r'[<>:"/\\|?*]', '_', cleaned_name)
        
        # Ensure reasonable length
        if len(cleaned_name) > 100:
            raise ValidationError(f"Gene name too long: {len(cleaned_name)} > 100", "gene_name", cleaned_name)
        
        return cleaned_name
    
    @staticmethod
    def validate_file_path(file_path: Union[str, Path], must_exist: bool = True, 
                          must_be_file: bool = True) -> Path:
        """
        Validate file path with comprehensive checks.
        
        Args:
            file_path: Path to validate
            must_exist: Whether file must exist
            must_be_file: Whether path must be a file (vs directory)
            
        Returns:
            Validated Path object
            
        Raises:
            ValidationError: If path is invalid
        """
        if not file_path:
            raise ValidationError("File path cannot be empty", "file_path", file_path)
        
        try:
            path = Path(file_path)
        except Exception as e:
            raise ValidationError(f"Invalid path format: {e}", "file_path", file_path)
        
        if must_exist:
            if not path.exists():
                raise ValidationError(f"Path does not exist: {path}", "file_path", str(path))
            
            if must_be_file and not path.is_file():
                raise ValidationError(f"Path is not a file: {path}", "file_path", str(path))
            
            if not must_be_file and not path.is_dir():
                raise ValidationError(f"Path is not a directory: {path}", "file_path", str(path))
        
        # Check permissions if exists
        if must_exist and not os.access(path, os.R_OK):
            raise ValidationError(f"No read permission for: {path}", "file_path", str(path))
        
        return path
    
    @staticmethod
    def validate_numeric_parameter(value: Any, param_name: str, 
                                 min_value: float = None, max_value: float = None,
                                 allow_none: bool = False) -> float:
        """
        Validate numeric parameters with range checking.
        
        Args:
            value: Value to validate
            param_name: Parameter name for error messages
            min_value: Minimum allowed value
            max_value: Maximum allowed value
            allow_none: Whether None values are acceptable
            
        Returns:
            Validated numeric value
            
        Raises:
            ValidationError: If value is invalid
        """
        if value is None:
            if allow_none:
                return None
            else:
                raise ValidationError(f"{param_name} cannot be None", param_name, value)
        
        try:
            numeric_value = float(value)
        except (ValueError, TypeError):
            raise ValidationError(f"{param_name} must be numeric, got {type(value)}", param_name, value)
        
        if min_value is not None and numeric_value < min_value:
            raise ValidationError(
                f"{param_name} must be >= {min_value}, got {numeric_value}", 
                param_name, numeric_value
            )
        
        if max_value is not None and numeric_value > max_value:
            raise ValidationError(
                f"{param_name} must be <= {max_value}, got {numeric_value}", 
                param_name, numeric_value
            )
        
        return numeric_value


def robust_exception_handler(logger: RobustLogger = None):
    """
    Decorator for robust exception handling with automatic logging and recovery.
    
    Args:
        logger: Logger instance to use
        
    Returns:
        Decorator function
    """
    def decorator(func: Callable) -> Callable:
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            nonlocal logger
            if logger is None:
                logger = RobustLogger()
            
            try:
                return func(*args, **kwargs)
                
            except PipelineError:
                # Re-raise pipeline errors as-is
                raise
                
            except Exception as e:
                # Convert unexpected errors to pipeline errors
                context = {
                    'function': func.__name__,
                    'args': str(args)[:200] + '...' if len(str(args)) > 200 else str(args),
                    'kwargs': str(kwargs)[:200] + '...' if len(str(kwargs)) > 200 else str(kwargs)
                }
                
                logger.log_exception(e, context)
                
                # Determine appropriate pipeline error type
                if "file" in str(e).lower() or "path" in str(e).lower():
                    raise FileSystemError(f"File system error in {func.__name__}: {e}")
                elif "memory" in str(e).lower():
                    raise DataProcessingError(f"Memory error in {func.__name__}: {e}")
                else:
                    raise DataProcessingError(f"Unexpected error in {func.__name__}: {e}")
        
        return wrapper
    return decorator


def setup_global_error_handling():
    """Set up global exception handling for the entire pipeline."""
    def handle_exception(exc_type, exc_value, exc_traceback):
        if issubclass(exc_type, KeyboardInterrupt):
            # Allow keyboard interrupts to pass through
            sys.__excepthook__(exc_type, exc_value, exc_traceback)
            return
        
        logger = RobustLogger()
        
        error_report = {
            'error_type': exc_type.__name__,
            'error_message': str(exc_value),
            'traceback': traceback.format_tb(exc_traceback),
            'timestamp': datetime.now().isoformat()
        }
        
        logger.logger.critical(f"Unhandled exception: {json.dumps(error_report, indent=2)}")
        
        # Save error report to file
        try:
            error_file = Path("logs") / f"crash_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
            error_file.parent.mkdir(parents=True, exist_ok=True)
            
            with open(error_file, 'w') as f:
                json.dump(error_report, f, indent=2, default=str)
                
            print(f"\nCrash report saved to: {error_file}")
            
        except Exception:
            print("\nCould not save crash report")
        
        # Call original exception handler
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
    
    sys.excepthook = handle_exception


# Initialize global error handling
setup_global_error_handling()
