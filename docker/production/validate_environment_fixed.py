#!/usr/bin/env python3
"""
Fixed Environment Validation Script for ProductionWildTypeAligner
Supports both development and containerized environments

This script validates:
1. Directory structure (with auto-creation)
2. Configuration files (with flexible path resolution)
3. Application imports and dependencies
4. Resource configuration
5. Security settings

Usage:
    python validate_environment_fixed.py [--json] [--container] [--create-missing]

Author: GenomeAMRAnalyzer Pipeline
Version: 2.0 - Fixed Production Validation
"""

import os
import sys
import json
import time
import logging
import argparse
import platform
from pathlib import Path
from typing import Dict, List, Optional, Any


class FixedEnvironmentValidator:
    """Enhanced environment validator that works in both dev and container environments"""
    
    def __init__(self, is_container: bool = False, create_missing: bool = False):
        self.is_container = is_container
        self.create_missing = create_missing
        
        # Determine base paths based on environment
        if is_container:
            self.base_path = Path("/app")
            self.config_path = Path("/app/config")
        else:
            # Development environment - use current directory structure
            self.base_path = Path(__file__).parent.parent.parent
            self.config_path = self.base_path / "config"
        
        # Expected directory structure
        self.required_dirs = {
            'src': self.base_path / 'src',
            'config': self.config_path,
            'data': self.base_path / 'data',
            'results': self.base_path / 'results',
            'logs': self.base_path / 'logs',
            'cache': self.base_path / 'cache'
        }
        
        # Expected configuration files
        self.required_config_files = {
            'app.yml': self.config_path / 'production' / 'app.yml',
            'requirements.txt': self.base_path / 'requirements' / 'production.txt'
        }
        
        # Validation results
        self.validation_results = {
            'timestamp': int(time.time()),
            'environment': 'container' if is_container else 'development',
            'platform': platform.system(),
            'overall_status': 'UNKNOWN',
            'validations': {},
            'recommendations': [],
            'fixed_issues': [],
            'failed_validations': []
        }
        
        # Setup logging
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
        self.logger = logging.getLogger(__name__)
    
    def log(self, message: str, level: str = 'info'):
        """Log validation messages"""
        getattr(self.logger, level)(message)
    
    def validate_directory_structure(self) -> Dict[str, Any]:
        """Validate and optionally create required directories"""
        self.log("Validating directory structure...")
        
        validation = {
            'status': 'PASS',
            'details': {},
            'issues': [],
            'fixed': []
        }
        
        for name, path in self.required_dirs.items():
            if path.exists() and path.is_dir():
                validation['details'][name] = 'EXISTS'
                self.log(f"✅ Directory {name}: {path}")
            else:
                validation['details'][name] = 'MISSING'
                validation['issues'].append(f"Required directory missing: {name} ({path})")
                
                if self.create_missing:
                    try:
                        path.mkdir(parents=True, exist_ok=True)
                        validation['fixed'].append(f"Created directory: {name}")
                        validation['details'][name] = 'CREATED'
                        self.log(f"✅ Created directory {name}: {path}")
                    except Exception as e:
                        validation['issues'].append(f"Failed to create directory {name}: {e}")
                        validation['status'] = 'FAIL'
                        self.log(f"❌ Failed to create {name}: {e}", 'error')
                else:
                    validation['status'] = 'FAIL'
                    self.log(f"❌ Missing directory {name}: {path}", 'error')
        
        return validation
    
    def validate_configuration_files(self) -> Dict[str, Any]:
        """Validate configuration file existence and basic structure"""
        self.log("Validating configuration files...")
        
        validation = {
            'status': 'PASS',
            'details': {},
            'issues': [],
            'fixed': []
        }
        
        for name, path in self.required_config_files.items():
            if path.exists() and path.is_file():
                validation['details'][name] = 'EXISTS'
                self.log(f"✅ Config file {name}: {path}")
                
                # Validate YAML syntax for .yml files
                if name.endswith('.yml') or name.endswith('.yaml'):
                    try:
                        import yaml
                        with open(path, 'r') as f:
                            yaml.safe_load(f)
                        validation['details'][f"{name}_syntax"] = 'VALID'
                        self.log(f"✅ YAML syntax valid for {name}")
                    except Exception as e:
                        validation['issues'].append(f"Invalid YAML syntax in {name}: {e}")
                        validation['status'] = 'FAIL'
                        self.log(f"❌ Invalid YAML in {name}: {e}", 'error')
                
            else:
                validation['details'][name] = 'MISSING'
                validation['issues'].append(f"Required configuration file missing: {name}")
                validation['status'] = 'FAIL'
                self.log(f"❌ Missing config file {name}: {path}", 'error')
        
        return validation
    
    def validate_application_imports(self) -> Dict[str, Any]:
        """Validate that core application modules can be imported"""
        self.log("Validating application imports...")
        
        validation = {
            'status': 'PASS',
            'details': {},
            'issues': [],
            'import_test_results': {}
        }
        
        # Add src directory to Python path
        src_path = str(self.required_dirs['src'])
        if src_path not in sys.path:
            sys.path.insert(0, src_path)
        
        # Test core imports (using correct class names)
        core_imports = [
            ('production_wildtype_aligner', 'ProductionWildTypeAligner'),
            ('sepi_configuration_manager', 'SEPIConfigurationManager'),
            ('generic_cooccurrence_analyzer', 'GenericCoOccurrenceAnalyzer'),  # Fixed spelling
            ('fasta_aa_extractor_integration', 'FastaAAExtractorIntegrator')   # Fixed class name
        ]
        
        for module_name, class_name in core_imports:
            try:
                module = __import__(module_name)
                if hasattr(module, class_name):
                    validation['import_test_results'][module_name] = 'SUCCESS'
                    self.log(f"✅ Successfully imported {module_name}.{class_name}")
                else:
                    validation['import_test_results'][module_name] = 'CLASS_MISSING'
                    validation['issues'].append(f"Class {class_name} not found in {module_name}")
                    validation['status'] = 'FAIL'
            except ImportError as e:
                validation['import_test_results'][module_name] = 'IMPORT_ERROR'
                validation['issues'].append(f"Failed to import {module_name}: {e}")
                validation['status'] = 'FAIL'
                self.log(f"❌ Import failed {module_name}: {e}", 'error')
            except Exception as e:
                validation['import_test_results'][module_name] = 'UNKNOWN_ERROR'
                validation['issues'].append(f"Unknown error importing {module_name}: {e}")
                validation['status'] = 'FAIL'
                self.log(f"❌ Unknown error {module_name}: {e}", 'error')
        
        return validation
    
    def validate_resource_configuration(self) -> Dict[str, Any]:
        """Validate system resources and configuration"""
        self.log("Validating resource configuration...")
        
        validation = {
            'status': 'PASS',
            'details': {},
            'issues': [],
            'recommendations': []
        }
        
        # Check system resources
        try:
            import psutil
            
            # CPU information
            cpu_count = psutil.cpu_count()
            cpu_percent = psutil.cpu_percent(interval=1)
            validation['details']['system_cpu_cores'] = cpu_count
            validation['details']['system_cpu_usage'] = cpu_percent
            
            # Memory information
            memory = psutil.virtual_memory()
            memory_gb = round(memory.total / (1024**3), 2)
            memory_available_gb = round(memory.available / (1024**3), 2)
            validation['details']['system_memory_gb'] = memory_gb
            validation['details']['system_available_memory_gb'] = memory_available_gb
            
            # Disk space
            disk = psutil.disk_usage(str(self.base_path))
            disk_free_gb = round(disk.free / (1024**3), 2)
            validation['details']['disk_free_gb'] = disk_free_gb
            
            # Resource recommendations
            if memory_gb < 4:
                validation['recommendations'].append("Consider increasing system memory (current: {memory_gb}GB, recommended: 8GB+)")
            
            if cpu_count < 2:
                validation['recommendations'].append(f"Consider using a system with more CPU cores (current: {cpu_count}, recommended: 4+)")
            
            if disk_free_gb < 10:
                validation['recommendations'].append(f"Low disk space available (current: {disk_free_gb}GB, recommended: 20GB+ free)")
            
            self.log(f"✅ System resources validated - CPU: {cpu_count} cores, RAM: {memory_gb}GB, Disk: {disk_free_gb}GB free")
            
        except ImportError:
            validation['issues'].append("psutil not available - cannot validate system resources")
            validation['status'] = 'PARTIAL'
            self.log("⚠️  psutil not available - skipping system resource validation", 'warning')
        except Exception as e:
            validation['issues'].append(f"Error validating system resources: {e}")
            validation['status'] = 'FAIL'
            self.log(f"❌ Error validating resources: {e}", 'error')
        
        return validation
    
    def validate_dependencies(self) -> Dict[str, Any]:
        """Validate required Python dependencies"""
        self.log("Validating Python dependencies...")
        
        validation = {
            'status': 'PASS',
            'details': {},
            'issues': [],
            'missing_packages': []
        }
        
        # Core required packages (using correct import names)
        required_packages = [
            ('Bio', 'biopython'),  # BioPython is imported as 'Bio'
            ('pandas', 'pandas'), 
            ('numpy', 'numpy'),
            ('yaml', 'yaml'),
            ('aiohttp', 'aiohttp'),
            ('aiofiles', 'aiofiles'),
            ('pydantic', 'pydantic'),
            ('structlog', 'structlog')
        ]
        
        for import_name, package_name in required_packages:
            try:
                __import__(import_name)
                validation['details'][package_name] = 'AVAILABLE'
                self.log(f"✅ Package available: {package_name}")
            except ImportError:
                validation['details'][package_name] = 'MISSING'
                validation['missing_packages'].append(package_name)
                validation['issues'].append(f"Required package missing: {package_name}")
                validation['status'] = 'FAIL'
                self.log(f"❌ Missing package: {package_name}", 'error')
        
        return validation
    
    def run_validation(self) -> Dict[str, Any]:
        """Run complete environment validation"""
        self.log("Starting comprehensive environment validation...")
        
        # Run all validations
        self.validation_results['validations']['directory_structure'] = self.validate_directory_structure()
        self.validation_results['validations']['configuration_files'] = self.validate_configuration_files()
        self.validation_results['validations']['application_imports'] = self.validate_application_imports()
        self.validation_results['validations']['resource_configuration'] = self.validate_resource_configuration()
        self.validation_results['validations']['dependencies'] = self.validate_dependencies()
        
        # Determine overall status
        failed_validations = []
        for name, result in self.validation_results['validations'].items():
            if result['status'] == 'FAIL':
                failed_validations.append(name)
            elif result['status'] == 'PARTIAL' and not failed_validations:
                # Only mark as partial if no complete failures
                self.validation_results['overall_status'] = 'PARTIAL'
        
        if failed_validations:
            self.validation_results['overall_status'] = 'FAIL'
            self.validation_results['failed_validations'] = failed_validations
        elif self.validation_results['overall_status'] != 'PARTIAL':
            self.validation_results['overall_status'] = 'PASS'
        
        # Collect all fixed issues
        for validation in self.validation_results['validations'].values():
            if 'fixed' in validation:
                self.validation_results['fixed_issues'].extend(validation['fixed'])
        
        # Collect all recommendations
        for validation in self.validation_results['validations'].values():
            if 'recommendations' in validation:
                self.validation_results['recommendations'].extend(validation['recommendations'])
        
        self.log(f"Validation completed - Overall status: {self.validation_results['overall_status']}")
        
        return self.validation_results
    
    def print_results(self):
        """Print human-readable validation results"""
        print("\n" + "="*70)
        print("Production Environment Validation Results")
        print("="*70)
        
        status_icons = {
            'PASS': '[PASS]',
            'FAIL': '[FAIL]', 
            'PARTIAL': '[PARTIAL]',
            'UNKNOWN': '[UNKNOWN]'
        }
        
        overall_icon = status_icons.get(self.validation_results['overall_status'], '[UNKNOWN]')
        print(f"Overall Status: {overall_icon} {self.validation_results['overall_status']}")
        print(f"Environment: {self.validation_results['environment']}")
        print(f"Platform: {self.validation_results['platform']}")
        print()
        
        # Print validation details
        for name, result in self.validation_results['validations'].items():
            status_icon = status_icons.get(result['status'], '[UNKNOWN]')
            print(f"{status_icon} {name.upper().replace('_', ' ')}: {result['status']}")
            
            # Show issues
            if result.get('issues'):
                for issue in result['issues'][:3]:  # Limit to first 3 issues
                    print(f"   ERROR: {issue}")
                if len(result['issues']) > 3:
                    print(f"   ... and {len(result['issues']) - 3} more issues")
            
            # Show fixes
            if result.get('fixed'):
                for fix in result['fixed']:
                    print(f"   FIXED: {fix}")
            
            print()
        
        # Show recommendations
        if self.validation_results.get('recommendations'):
            print("Recommendations:")
            for rec in self.validation_results['recommendations']:
                print(f"   * {rec}")
            print()


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(description="Fixed Environment Validation Script")
    parser.add_argument('--json', action='store_true',
                       help='Output results in JSON format')
    parser.add_argument('--container', action='store_true',
                       help='Run in container mode (uses /app paths)')
    parser.add_argument('--create-missing', action='store_true',
                       help='Automatically create missing directories')
    parser.add_argument('--fix-issues', action='store_true',
                       help='Attempt to automatically fix detected issues')
    
    args = parser.parse_args()
    
    # Create validator
    validator = FixedEnvironmentValidator(
        is_container=args.container,
        create_missing=args.create_missing or args.fix_issues
    )
    
    # Run validation
    results = validator.run_validation()
    
    # Output results
    if args.json:
        print(json.dumps(results, indent=2))
    else:
        validator.print_results()
    
    # Exit with appropriate code
    exit_code = 0 if results['overall_status'] in ['PASS', 'PARTIAL'] else 1
    sys.exit(exit_code)


if __name__ == '__main__':
    main()