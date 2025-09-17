#!/usr/bin/env python3
"""
Production Environment Validation for ProductionWildTypeAligner
Validates production deployment readiness and configuration

This script validates:
1. Container environment setup
2. Configuration file integrity
3. Resource allocation correctness
4. Security configuration
5. Performance optimization settings

Usage:
    python validate_environment.py [--config-dir /path/to/config] [--fix-issues]

Author: GenomeAMRAnalyzer Pipeline
Version: 2.0 - Environment Validation
"""

import os
import sys
import yaml
import json
import logging
import argparse
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple


class EnvironmentValidator:
    """Comprehensive environment validation for production deployment"""
    
    def __init__(self, config_dir: str = '/app/config/production', fix_issues: bool = False):
        self.config_dir = Path(config_dir)
        self.fix_issues = fix_issues
        self.validation_results = {
            'timestamp': int(time.time()) if 'time' in sys.modules else 0,
            'overall_status': 'UNKNOWN',
            'validations': {},
            'recommendations': []
        }
        
        # Setup logging
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
        self.logger = logging.getLogger(__name__)
    
    def log(self, message: str, level: str = 'info'):
        """Log validation messages"""
        getattr(self.logger, level)(message)
    
    def validate_directory_structure(self) -> Dict[str, Any]:
        """Validate required directory structure"""
        self.log("Validating directory structure...")
        
        result = {
            'status': 'PASS',
            'details': {},
            'issues': [],
            'fixed': []
        }
        
        # Required directories
        required_dirs = [
            '/app/src',
            '/app/config',
            '/app/data',
            '/app/results',
            '/app/logs',
            '/app/cache'
        ]
        
        for dir_path in required_dirs:
            path = Path(dir_path)
            if path.exists():
                if path.is_dir():
                    # Check permissions
                    if os.access(dir_path, os.R_OK | os.W_OK):
                        result['details'][dir_path] = 'OK'
                    else:
                        result['details'][dir_path] = 'PERMISSION_ERROR'
                        result['issues'].append(f"No read/write permissions for {dir_path}")
                        result['status'] = 'FAIL'
                        
                        if self.fix_issues:
                            try:
                                os.chmod(dir_path, 0o755)
                                result['fixed'].append(f"Fixed permissions for {dir_path}")
                            except Exception as e:
                                result['issues'].append(f"Could not fix permissions for {dir_path}: {e}")
                else:
                    result['details'][dir_path] = 'NOT_DIRECTORY'
                    result['issues'].append(f"{dir_path} exists but is not a directory")
                    result['status'] = 'FAIL'
            else:
                result['details'][dir_path] = 'MISSING'
                result['issues'].append(f"Required directory missing: {dir_path}")
                result['status'] = 'FAIL'
                
                if self.fix_issues:
                    try:
                        path.mkdir(parents=True, exist_ok=True)
                        result['fixed'].append(f"Created missing directory: {dir_path}")
                        result['details'][dir_path] = 'CREATED'
                    except Exception as e:
                        result['issues'].append(f"Could not create directory {dir_path}: {e}")
        
        return result
    
    def validate_configuration_files(self) -> Dict[str, Any]:
        """Validate configuration file structure and content"""
        self.log("Validating configuration files...")
        
        result = {
            'status': 'PASS',
            'details': {},
            'issues': [],
            'fixed': []
        }
        
        # Required configuration files
        config_files = {
            'app.yml': 'Application configuration',
            '.env.template': 'Environment variables template'
        }
        
        for filename, description in config_files.items():
            file_path = self.config_dir / filename
            
            if file_path.exists():
                try:
                    if filename.endswith('.yml') or filename.endswith('.yaml'):
                        with open(file_path, 'r') as f:
                            config_data = yaml.safe_load(f)
                        
                        # Validate YAML structure
                        if isinstance(config_data, dict):
                            result['details'][filename] = 'VALID_YAML'
                            
                            # Validate required sections for app.yml
                            if filename == 'app.yml':
                                required_sections = ['application', 'logging', 'performance', 'security']
                                missing_sections = [section for section in required_sections 
                                                  if section not in config_data]
                                
                                if missing_sections:
                                    result['issues'].append(f"Missing sections in {filename}: {missing_sections}")
                                    result['status'] = 'FAIL'
                                else:
                                    result['details'][filename] = 'COMPLETE'
                        else:
                            result['details'][filename] = 'INVALID_YAML'
                            result['issues'].append(f"{filename} is not a valid YAML dictionary")
                            result['status'] = 'FAIL'
                    
                    elif filename.startswith('.env'):
                        # Validate environment file
                        with open(file_path, 'r') as f:
                            env_lines = f.readlines()
                        
                        env_vars = [line.strip() for line in env_lines 
                                   if line.strip() and not line.startswith('#')]
                        
                        result['details'][filename] = f"{len(env_vars)} variables defined"
                    
                except Exception as e:
                    result['details'][filename] = 'READ_ERROR'
                    result['issues'].append(f"Could not read {filename}: {e}")
                    result['status'] = 'FAIL'
            else:
                result['details'][filename] = 'MISSING'
                result['issues'].append(f"Required configuration file missing: {filename}")
                result['status'] = 'FAIL'
        
        return result
    
    def validate_resource_configuration(self) -> Dict[str, Any]:
        """Validate resource allocation and limits"""
        self.log("Validating resource configuration...")
        
        result = {
            'status': 'PASS',
            'details': {},
            'issues': [],
            'recommendations': []
        }
        
        # Check Docker resource limits
        try:
            import psutil
            
            # CPU information
            cpu_count = psutil.cpu_count()
            result['details']['system_cpu_cores'] = cpu_count
            
            # Memory information
            memory = psutil.virtual_memory()
            total_memory_gb = memory.total / (1024**3)
            result['details']['system_memory_gb'] = round(total_memory_gb, 2)
            
            # Load application configuration to check resource settings
            app_config_path = self.config_dir / 'app.yml'
            if app_config_path.exists():
                with open(app_config_path, 'r') as f:
                    app_config = yaml.safe_load(f)
                
                if 'performance' in app_config:
                    perf_config = app_config['performance']
                    
                    # Check CPU configuration
                    configured_cpu = perf_config.get('cpu_cores', 1)
                    if configured_cpu > cpu_count:
                        result['issues'].append(f"Configured CPU cores ({configured_cpu}) exceeds available ({cpu_count})")
                        result['status'] = 'FAIL'
                    elif configured_cpu > cpu_count * 0.8:
                        result['recommendations'].append(f"Consider reducing CPU allocation from {configured_cpu} to {int(cpu_count * 0.8)}")
                    
                    # Check memory configuration
                    configured_memory = perf_config.get('memory_limit_gb', 4)
                    if configured_memory > total_memory_gb:
                        result['issues'].append(f"Configured memory ({configured_memory}GB) exceeds available ({total_memory_gb:.1f}GB)")
                        result['status'] = 'FAIL'
                    elif configured_memory > total_memory_gb * 0.8:
                        result['recommendations'].append(f"Consider reducing memory allocation from {configured_memory}GB to {total_memory_gb * 0.8:.1f}GB")
                    
                    result['details']['configured_cpu_cores'] = configured_cpu
                    result['details']['configured_memory_gb'] = configured_memory
            
        except ImportError:
            result['details']['resource_check'] = 'PSUTIL_NOT_AVAILABLE'
            result['recommendations'].append("Install psutil for detailed resource validation")
        except Exception as e:
            result['issues'].append(f"Error checking system resources: {e}")
            result['status'] = 'FAIL'
        
        return result
    
    def validate_security_configuration(self) -> Dict[str, Any]:
        """Validate security configuration and settings"""
        self.log("Validating security configuration...")
        
        result = {
            'status': 'PASS',
            'details': {},
            'issues': [],
            'recommendations': []
        }
        
        # Check file permissions
        critical_files = [
            self.config_dir / 'app.yml',
            self.config_dir / '.env.template'
        ]
        
        for file_path in critical_files:
            if file_path.exists():
                stat = file_path.stat()
                permissions = oct(stat.st_mode)[-3:]
                
                result['details'][str(file_path)] = f"permissions: {permissions}"
                
                # Check if file is world-readable (security risk)
                if permissions.endswith('4') or permissions.endswith('6') or permissions.endswith('7'):
                    result['issues'].append(f"{file_path} is world-readable (permissions: {permissions})")
                    result['status'] = 'FAIL'
                    
                    if self.fix_issues:
                        try:
                            file_path.chmod(0o640)  # rw-r-----
                            result['fixed'].append(f"Fixed permissions for {file_path}")
                        except Exception as e:
                            result['issues'].append(f"Could not fix permissions for {file_path}: {e}")
        
        # Check for sensitive data in configuration
        try:
            app_config_path = self.config_dir / 'app.yml'
            if app_config_path.exists():
                with open(app_config_path, 'r') as f:
                    config_content = f.read()
                
                # Look for potential secrets (simple heuristic)
                sensitive_patterns = ['password', 'secret', 'key', 'token']
                found_patterns = []
                
                for pattern in sensitive_patterns:
                    if pattern in config_content.lower():
                        found_patterns.append(pattern)
                
                if found_patterns:
                    result['recommendations'].append(f"Consider moving sensitive data to environment variables: {found_patterns}")
                
                result['details']['security_scan'] = f"Checked for {len(sensitive_patterns)} patterns"
        
        except Exception as e:
            result['issues'].append(f"Error scanning configuration for security issues: {e}")
        
        return result
    
    def validate_application_configuration(self) -> Dict[str, Any]:
        """Validate application-specific configuration"""
        self.log("Validating application configuration...")
        
        result = {
            'status': 'PASS',
            'details': {},
            'issues': [],
            'recommendations': []
        }
        
        try:
            app_config_path = self.config_dir / 'app.yml'
            if app_config_path.exists():
                with open(app_config_path, 'r') as f:
                    app_config = yaml.safe_load(f)
                
                # Validate application section
                if 'application' in app_config:
                    app_section = app_config['application']
                    
                    # Check debug mode
                    debug_mode = app_section.get('debug', True)
                    if debug_mode:
                        result['issues'].append("Debug mode is enabled in production configuration")
                        result['status'] = 'FAIL'
                    
                    # Check environment setting
                    environment = app_section.get('environment', 'development')
                    if environment != 'production':
                        result['issues'].append(f"Environment is set to '{environment}', should be 'production'")
                        result['status'] = 'FAIL'
                    
                    result['details']['debug_mode'] = debug_mode
                    result['details']['environment'] = environment
                
                # Validate logging configuration
                if 'logging' in app_config:
                    logging_config = app_config['logging']
                    
                    log_level = logging_config.get('level', 'DEBUG')
                    if log_level == 'DEBUG':
                        result['recommendations'].append("Consider setting log level to INFO or WARNING in production")
                    
                    enable_console = logging_config.get('enable_console', True)
                    if enable_console:
                        result['recommendations'].append("Consider disabling console logging in production")
                    
                    result['details']['log_level'] = log_level
                    result['details']['console_logging'] = enable_console
                
                # Validate performance settings
                if 'performance' in app_config:
                    perf_config = app_config['performance']
                    
                    timeout = perf_config.get('timeout_seconds', 300)
                    if timeout < 300:
                        result['recommendations'].append(f"Timeout ({timeout}s) might be too low for large datasets")
                    
                    max_concurrent = perf_config.get('max_concurrent_jobs', 1)
                    if max_concurrent < 2:
                        result['recommendations'].append("Consider increasing max_concurrent_jobs for better performance")
                    
                    result['details']['timeout_seconds'] = timeout
                    result['details']['max_concurrent_jobs'] = max_concurrent
        
        except Exception as e:
            result['issues'].append(f"Error validating application configuration: {e}")
            result['status'] = 'FAIL'
        
        return result
    
    def run_all_validations(self) -> Dict[str, Any]:
        """Run all environment validations"""
        self.log("Starting comprehensive environment validation...")
        
        # Run all validation checks
        self.validation_results['validations']['directory_structure'] = self.validate_directory_structure()
        self.validation_results['validations']['configuration_files'] = self.validate_configuration_files()
        self.validation_results['validations']['resource_configuration'] = self.validate_resource_configuration()
        self.validation_results['validations']['security_configuration'] = self.validate_security_configuration()
        self.validation_results['validations']['application_configuration'] = self.validate_application_configuration()
        
        # Collect all recommendations
        for validation in self.validation_results['validations'].values():
            self.validation_results['recommendations'].extend(validation.get('recommendations', []))
        
        # Determine overall status
        failed_validations = [name for name, result in self.validation_results['validations'].items()
                            if result['status'] == 'FAIL']
        
        if failed_validations:
            self.validation_results['overall_status'] = 'FAIL'
            self.validation_results['failed_validations'] = failed_validations
        else:
            self.validation_results['overall_status'] = 'PASS'
        
        return self.validation_results
    
    def print_results(self):
        """Print validation results"""
        print("\n" + "="*60)
        print("🔧 Production Environment Validation Results")
        print("="*60)
        
        status_icon = "✅" if self.validation_results['overall_status'] == 'PASS' else "❌"
        print(f"Overall Status: {status_icon} {self.validation_results['overall_status']}")
        print()
        
        for validation_name, validation_result in self.validation_results['validations'].items():
            status_icon = "✅" if validation_result['status'] == 'PASS' else "❌"
            print(f"{status_icon} {validation_name.upper().replace('_', ' ')}: {validation_result['status']}")
            
            # Show issues
            if validation_result.get('issues'):
                for issue in validation_result['issues']:
                    print(f"   ❌ {issue}")
            
            # Show fixes applied
            if validation_result.get('fixed'):
                for fix in validation_result['fixed']:
                    print(f"   🔧 {fix}")
            
            # Show details
            if validation_result.get('details'):
                for key, value in validation_result['details'].items():
                    print(f"   📊 {key}: {value}")
            
            print()
        
        # Show recommendations
        if self.validation_results['recommendations']:
            print("💡 RECOMMENDATIONS:")
            for recommendation in set(self.validation_results['recommendations']):  # Remove duplicates
                print(f"   💡 {recommendation}")
            print()


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(description="Production Environment Validation")
    parser.add_argument('--config-dir', default='/app/config/production',
                       help='Path to configuration directory')
    parser.add_argument('--fix-issues', action='store_true',
                       help='Attempt to fix detected issues automatically')
    parser.add_argument('--json', action='store_true',
                       help='Output results in JSON format')
    
    args = parser.parse_args()
    
    # Initialize validator
    validator = EnvironmentValidator(
        config_dir=args.config_dir,
        fix_issues=args.fix_issues
    )
    
    # Run validations
    results = validator.run_all_validations()
    
    # Output results
    if args.json:
        print(json.dumps(results, indent=2))
    else:
        validator.print_results()
    
    # Exit with appropriate code
    exit_code = 0 if results['overall_status'] == 'PASS' else 1
    sys.exit(exit_code)


if __name__ == '__main__':
    import time  # Import here to avoid issues with timestamp
    main()