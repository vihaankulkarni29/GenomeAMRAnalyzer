"""
Comprehensive dependency management and validation for GenomeAMRAnalyzer.
Provides automatic detection, validation, and installation of required packages.
"""
import sys
import importlib
import subprocess
from typing import List, Dict, Optional, Tuple
import logging
from pathlib import Path

# Required packages with version constraints
REQUIRED_PACKAGES = {
    'pandas': '>=1.5.0',
    'numpy': '>=1.24.0', 
    'biopython': '>=1.80',
    'requests': '>=2.28.0',
    'pyyaml': '>=6.0',
    'matplotlib': '>=3.6.0',
    'seaborn': '>=0.12.0',
    'scikit-learn': '>=1.2.0',
    'jinja2': '>=3.1.0',
}

# Optional packages for enhanced functionality
OPTIONAL_PACKAGES = {
    'plotly': '>=5.0.0',
    'dash': '>=2.0.0',
    'jupyter': '>=1.0.0',
    'psutil': '>=5.9.0'
}

class DependencyError(Exception):
    """Exception raised for dependency-related errors."""
    pass

class DependencyChecker:
    """Comprehensive dependency validation and management."""
    
    def __init__(self, logger: Optional[logging.Logger] = None):
        self.logger = logger or logging.getLogger(__name__)
        self.missing_packages = []
        self.version_conflicts = []
        self.import_errors = {}
    
    def check_all_dependencies(self, include_optional: bool = False) -> bool:
        """
        Check all required dependencies and optionally optional ones.
        
        Args:
            include_optional: Whether to check optional packages
            
        Returns:
            True if all required dependencies are satisfied
            
        Raises:
            DependencyError: If critical dependencies are missing
        """
        self.logger.info("Checking dependencies...")
        
        # Check required packages
        packages_to_check = REQUIRED_PACKAGES.copy()
        if include_optional:
            packages_to_check.update(OPTIONAL_PACKAGES)
        
        all_satisfied = True
        
        for package_name, version_req in packages_to_check.items():
            is_optional = package_name in OPTIONAL_PACKAGES
            
            try:
                is_satisfied = self._check_single_package(package_name, version_req)
                if not is_satisfied and not is_optional:
                    all_satisfied = False
                    
            except Exception as e:
                self.import_errors[package_name] = str(e)
                if not is_optional:
                    all_satisfied = False
                    self.logger.error(f"Failed to check {package_name}: {e}")
                else:
                    self.logger.warning(f"Optional package {package_name} not available: {e}")
        
        # Report results
        self._report_dependency_status()
        
        if not all_satisfied:
            raise DependencyError(
                f"Missing required dependencies: {self.missing_packages}. "
                f"Run: pip install -r requirements.txt"
            )
        
        return all_satisfied
    
    def _check_single_package(self, package_name: str, version_req: str) -> bool:
        """Check a single package and its version."""
        try:
            # Import the package
            if package_name == 'biopython':
                import Bio
                module = Bio
                version = Bio.__version__
            elif package_name == 'scikit-learn':
                import sklearn
                module = sklearn
                version = sklearn.__version__
            elif package_name == 'pyyaml':
                import yaml
                module = yaml
                version = getattr(yaml, '__version__', 'unknown')
            else:
                module = importlib.import_module(package_name)
                version = getattr(module, '__version__', 'unknown')
            
            # Check version if specified
            if version != 'unknown' and version_req.startswith('>='):
                required_version = version_req[2:]
                if not self._version_satisfies(version, required_version):
                    self.version_conflicts.append(
                        f"{package_name}: found {version}, requires {version_req}"
                    )
                    return False
            
            self.logger.debug(f"✓ {package_name} {version}")
            return True
            
        except ImportError:
            self.missing_packages.append(package_name)
            self.logger.warning(f"✗ {package_name}: Not installed")
            return False
        except Exception as e:
            self.import_errors[package_name] = str(e)
            self.logger.error(f"✗ {package_name}: Unexpected error - {e}")
            return False
    
    def _version_satisfies(self, current_version: str, required_version: str) -> bool:
        """Check if current version satisfies requirement."""
        try:
            # Simple version comparison for basic cases
            current_parts = [int(x) for x in current_version.split('.')]
            required_parts = [int(x) for x in required_version.split('.')]
            
            # Pad with zeros to make same length
            max_len = max(len(current_parts), len(required_parts))
            current_parts.extend([0] * (max_len - len(current_parts)))
            required_parts.extend([0] * (max_len - len(required_parts)))
            
            return current_parts >= required_parts
        except Exception:
            # If version parsing fails, assume it's satisfied
            return True
    
    def _report_dependency_status(self):
        """Report comprehensive dependency status."""
        if not self.missing_packages and not self.version_conflicts and not self.import_errors:
            self.logger.info("✓ All dependencies satisfied")
            return
        
        if self.missing_packages:
            self.logger.error(f"Missing packages: {', '.join(self.missing_packages)}")
        
        if self.version_conflicts:
            self.logger.warning("Version conflicts:")
            for conflict in self.version_conflicts:
                self.logger.warning(f"  {conflict}")
        
        if self.import_errors:
            self.logger.error("Import errors:")
            for package, error in self.import_errors.items():
                self.logger.error(f"  {package}: {error}")
    
    def install_missing_packages(self, auto_install: bool = False) -> bool:
        """
        Install missing packages using pip.
        
        Args:
            auto_install: If True, automatically install without prompting
            
        Returns:
            True if installation successful
        """
        if not self.missing_packages:
            self.logger.info("No packages to install")
            return True
        
        if not auto_install:
            try:
                response = input(f"Install missing packages {self.missing_packages}? (y/n): ")
                if response.lower() != 'y':
                    return False
            except KeyboardInterrupt:
                return False
        
        try:
            for package in self.missing_packages:
                self.logger.info(f"Installing {package}...")
                subprocess.run([
                    sys.executable, "-m", "pip", "install", package
                ], check=True, capture_output=True, text=True)
            
            self.logger.info("All packages installed successfully")
            
            # Re-check dependencies
            self.missing_packages = []
            self.version_conflicts = []
            self.import_errors = {}
            return self.check_all_dependencies()
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Failed to install packages: {e}")
            return False
        except Exception as e:
            self.logger.error(f"Unexpected error during installation: {e}")
            return False

def safe_import(package_name: str, alternative_name: Optional[str] = None):
    """
    Safely import a package with fallback options.
    
    Args:
        package_name: Primary package to import
        alternative_name: Alternative package name if primary fails
        
    Returns:
        Imported module or None if import fails
    """
    try:
        return importlib.import_module(package_name)
    except ImportError:
        if alternative_name:
            try:
                return importlib.import_module(alternative_name)
            except ImportError:
                pass
        
        logging.warning(f"Could not import {package_name}")
        return None

def validate_runtime_environment() -> Dict[str, str]:
    """
    Validate the runtime environment and return system information.
    
    Returns:
        Dictionary with system information
    """
    import platform
    
    env_info = {
        'python_version': sys.version,
        'platform': platform.platform(),
        'architecture': platform.architecture()[0],
        'processor': platform.processor(),
        'python_executable': sys.executable,
    }
    
    # Check memory availability if psutil is available
    try:
        import psutil
        env_info['total_memory_gb'] = round(psutil.virtual_memory().total / (1024**3), 2)
        env_info['available_memory_gb'] = round(psutil.virtual_memory().available / (1024**3), 2)
    except ImportError:
        env_info['memory_info'] = 'psutil not available'
    
    return env_info

# Auto-check dependencies on import (but only in main execution)
def auto_dependency_check():
    """Automatically check dependencies when module is imported."""
    try:
        checker = DependencyChecker()
        checker.check_all_dependencies(include_optional=False)
    except DependencyError as e:
        logging.error(f"Dependency check failed: {e}")
        logging.info("To install missing dependencies, run: pip install -r requirements.txt")
    except Exception as e:
        logging.warning(f"Could not perform dependency check: {e}")

# Only run auto-check if not being imported
if __name__ != "__main__":
    auto_dependency_check()
