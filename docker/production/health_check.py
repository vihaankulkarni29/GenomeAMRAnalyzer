#!/usr/bin/env python3
"""
Production Health Check Script for ProductionWildTypeAligner
Comprehensive system health validation and monitoring

This script provides:
1. Application health validation
2. Dependency checking
3. Resource monitoring
4. Service connectivity testing
5. Performance benchmarking

Usage:
    python health_check.py [--full] [--json] [--quiet]

Author: GenomeAMRAnalyzer Pipeline
Version: 2.0 - Production Health Check
"""

import os
import sys
import json
import time
import logging
import argparse
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple
import psutil
import asyncio

# Add src directory to path
sys.path.insert(0, '/app/src')

class HealthChecker:
    """Comprehensive health checking for ProductionWildTypeAligner"""
    
    def __init__(self, verbose: bool = True, json_output: bool = False):
        self.verbose = verbose
        self.json_output = json_output
        self.results = {
            'timestamp': time.time(),
            'overall_status': 'UNKNOWN',
            'checks': {}
        }
        
        # Setup logging
        log_level = logging.DEBUG if verbose else logging.WARNING
        logging.basicConfig(level=log_level, format='%(asctime)s - %(levelname)s - %(message)s')
        self.logger = logging.getLogger(__name__)
    
    def log(self, message: str, level: str = 'info'):
        """Log message if verbose mode enabled"""
        if self.verbose and not self.json_output:
            getattr(self.logger, level)(message)
    
    def check_application_imports(self) -> Dict[str, Any]:
        """Check if core application modules can be imported"""
        self.log("Checking application imports...")
        
        check_result = {
            'status': 'PASS',
            'details': {},
            'errors': []
        }
        
        # Test core imports
        imports_to_test = [
            ('production_wildtype_aligner', 'ProductionWildTypeAligner'),
            ('sepi_configuration_manager', 'SEPIConfigurationManager'),
            ('sepi_integration', 'EnhancedSEPIReferenceManager')
        ]
        
        for module_name, class_name in imports_to_test:
            try:
                module = __import__(module_name)
                if hasattr(module, class_name):
                    check_result['details'][f'{module_name}.{class_name}'] = 'OK'
                    self.log(f"✅ {module_name}.{class_name} imported successfully")
                else:
                    check_result['details'][f'{module_name}.{class_name}'] = 'CLASS_NOT_FOUND'
                    check_result['errors'].append(f"Class {class_name} not found in {module_name}")
                    check_result['status'] = 'FAIL'
            except ImportError as e:
                check_result['details'][f'{module_name}.{class_name}'] = 'IMPORT_ERROR'
                check_result['errors'].append(f"Import error for {module_name}: {e}")
                check_result['status'] = 'FAIL'
        
        return check_result
    
    def check_dependencies(self) -> Dict[str, Any]:
        """Check external dependencies availability"""
        self.log("Checking dependencies...")
        
        check_result = {
            'status': 'PASS',
            'details': {},
            'errors': []
        }
        
        # Test Python packages
        required_packages = [
            'Bio', 'pandas', 'numpy', 'yaml', 'aiohttp', 'pydantic'
        ]
        
        for package in required_packages:
            try:
                __import__(package)
                check_result['details'][package] = 'OK'
                self.log(f"✅ {package} available")
            except ImportError:
                check_result['details'][package] = 'MISSING'
                check_result['errors'].append(f"Required package {package} not available")
                check_result['status'] = 'FAIL'
        
        # Check EMBOSS WATER
        try:
            result = subprocess.run(['water', '-help'], 
                                  capture_output=True, text=True, timeout=10)
            if result.returncode == 0:
                check_result['details']['emboss_water'] = 'OK'
                self.log("✅ EMBOSS WATER available")
            else:
                check_result['details']['emboss_water'] = 'NOT_FUNCTIONAL'
                self.log("⚠️  EMBOSS WATER not functional")
        except (subprocess.TimeoutExpired, FileNotFoundError):
            check_result['details']['emboss_water'] = 'NOT_FOUND'
            self.log("⚠️  EMBOSS WATER not found (BioPython fallback will be used)")
        
        return check_result
    
    def check_file_system(self) -> Dict[str, Any]:
        """Check file system and directory structure"""
        self.log("Checking file system...")
        
        check_result = {
            'status': 'PASS',
            'details': {},
            'errors': []
        }
        
        # Check required directories
        required_dirs = [
            '/app/data', '/app/results', '/app/logs', '/app/cache'
        ]
        
        for dir_path in required_dirs:
            path = Path(dir_path)
            if path.exists() and path.is_dir():
                # Check permissions
                if os.access(dir_path, os.R_OK | os.W_OK):
                    check_result['details'][dir_path] = 'OK'
                    self.log(f"✅ {dir_path} accessible")
                else:
                    check_result['details'][dir_path] = 'PERMISSION_ERROR'
                    check_result['errors'].append(f"No read/write access to {dir_path}")
                    check_result['status'] = 'FAIL'
            else:
                check_result['details'][dir_path] = 'NOT_FOUND'
                check_result['errors'].append(f"Required directory not found: {dir_path}")
                check_result['status'] = 'FAIL'
        
        # Check disk space
        try:
            disk_usage = psutil.disk_usage('/')
            free_gb = disk_usage.free / (1024**3)
            check_result['details']['disk_space_gb'] = round(free_gb, 2)
            
            if free_gb < 1.0:  # Less than 1GB free
                check_result['status'] = 'FAIL'
                check_result['errors'].append(f"Low disk space: {free_gb:.2f}GB free")
            elif free_gb < 5.0:  # Less than 5GB free
                self.log(f"⚠️  Low disk space warning: {free_gb:.2f}GB free")
            else:
                self.log(f"✅ Sufficient disk space: {free_gb:.2f}GB free")
                
        except Exception as e:
            check_result['details']['disk_space_gb'] = 'ERROR'
            check_result['errors'].append(f"Could not check disk space: {e}")
        
        return check_result
    
    def check_system_resources(self) -> Dict[str, Any]:
        """Check system resource availability"""
        self.log("Checking system resources...")
        
        check_result = {
            'status': 'PASS',
            'details': {},
            'errors': []
        }
        
        try:
            # CPU usage
            cpu_percent = psutil.cpu_percent(interval=1)
            check_result['details']['cpu_usage_percent'] = cpu_percent
            
            if cpu_percent > 90:
                check_result['status'] = 'FAIL'
                check_result['errors'].append(f"High CPU usage: {cpu_percent}%")
            elif cpu_percent > 70:
                self.log(f"⚠️  High CPU usage: {cpu_percent}%")
            else:
                self.log(f"✅ CPU usage normal: {cpu_percent}%")
            
            # Memory usage
            memory = psutil.virtual_memory()
            memory_percent = memory.percent
            memory_available_gb = memory.available / (1024**3)
            
            check_result['details']['memory_usage_percent'] = memory_percent
            check_result['details']['memory_available_gb'] = round(memory_available_gb, 2)
            
            if memory_percent > 90:
                check_result['status'] = 'FAIL'
                check_result['errors'].append(f"High memory usage: {memory_percent}%")
            elif memory_percent > 80:
                self.log(f"⚠️  High memory usage: {memory_percent}%")
            else:
                self.log(f"✅ Memory usage normal: {memory_percent}%")
            
            # Load average
            load_avg = os.getloadavg()
            check_result['details']['load_average'] = load_avg
            self.log(f"✅ Load average: {load_avg}")
            
        except Exception as e:
            check_result['status'] = 'FAIL'
            check_result['errors'].append(f"Could not check system resources: {e}")
        
        return check_result
    
    def check_application_functionality(self) -> Dict[str, Any]:
        """Test basic application functionality"""
        self.log("Checking application functionality...")
        
        check_result = {
            'status': 'PASS',
            'details': {},
            'errors': []
        }
        
        try:
            # Test configuration manager
            from sepi_configuration_manager import SEPIConfigurationManager
            
            config_manager = SEPIConfigurationManager(workspace_root='/app')
            gene_config = config_manager.get_gene_config('acrA')
            
            if gene_config:
                check_result['details']['configuration_manager'] = 'OK'
                self.log("✅ Configuration manager functional")
            else:
                check_result['details']['configuration_manager'] = 'CONFIG_ERROR'
                check_result['errors'].append("Configuration manager not returning gene configs")
                check_result['status'] = 'FAIL'
                
        except Exception as e:
            check_result['details']['configuration_manager'] = 'ERROR'
            check_result['errors'].append(f"Configuration manager error: {e}")
            check_result['status'] = 'FAIL'
        
        try:
            # Test aligner initialization (without actual processing)
            from production_wildtype_aligner import ProductionWildTypeAligner
            
            # Create temporary test directory
            import tempfile
            with tempfile.TemporaryDirectory() as temp_dir:
                aligner = ProductionWildTypeAligner(
                    output_dir=temp_dir,
                    sepi_path='/app/MetaDataHarvester/sepi2.0/sepi.py',
                    email='healthcheck@example.com',
                    max_concurrent=1
                )
                
                check_result['details']['aligner_initialization'] = 'OK'
                self.log("✅ Aligner initialization successful")
                
        except Exception as e:
            check_result['details']['aligner_initialization'] = 'ERROR'
            check_result['errors'].append(f"Aligner initialization error: {e}")
            check_result['status'] = 'FAIL'
        
        return check_result
    
    def check_network_connectivity(self) -> Dict[str, Any]:
        """Test network connectivity for external services"""
        self.log("Checking network connectivity...")
        
        check_result = {
            'status': 'PASS',
            'details': {},
            'errors': []
        }
        
        # Test URLs to check
        test_urls = [
            ('ncbi', 'https://www.ncbi.nlm.nih.gov'),
            ('eutils', 'https://eutils.ncbi.nlm.nih.gov'),
        ]
        
        for service_name, url in test_urls:
            try:
                import urllib.request
                response = urllib.request.urlopen(url, timeout=10)
                if response.status == 200:
                    check_result['details'][service_name] = 'OK'
                    self.log(f"✅ {service_name} connectivity OK")
                else:
                    check_result['details'][service_name] = f'HTTP_{response.status}'
                    self.log(f"⚠️  {service_name} returned HTTP {response.status}")
            except Exception as e:
                check_result['details'][service_name] = 'ERROR'
                self.log(f"⚠️  {service_name} connectivity error: {e}")
                # Don't fail on network issues as they might be temporary
        
        return check_result
    
    def run_all_checks(self, full_check: bool = False) -> Dict[str, Any]:
        """Run all health checks"""
        self.log("Starting comprehensive health check...")
        
        # Core checks (always run)
        self.results['checks']['imports'] = self.check_application_imports()
        self.results['checks']['dependencies'] = self.check_dependencies()
        self.results['checks']['filesystem'] = self.check_file_system()
        self.results['checks']['resources'] = self.check_system_resources()
        
        # Extended checks (run with --full flag)
        if full_check:
            self.results['checks']['functionality'] = self.check_application_functionality()
            self.results['checks']['connectivity'] = self.check_network_connectivity()
        
        # Determine overall status
        failed_checks = [name for name, result in self.results['checks'].items() 
                        if result['status'] == 'FAIL']
        
        if failed_checks:
            self.results['overall_status'] = 'FAIL'
            self.results['failed_checks'] = failed_checks
        else:
            self.results['overall_status'] = 'PASS'
        
        self.results['check_duration'] = time.time() - self.results['timestamp']
        
        return self.results
    
    def print_results(self):
        """Print health check results"""
        if self.json_output:
            print(json.dumps(self.results, indent=2))
        else:
            print("\n" + "="*60)
            print("🧬 ProductionWildTypeAligner Health Check Results")
            print("="*60)
            
            print(f"Overall Status: {'✅ PASS' if self.results['overall_status'] == 'PASS' else '❌ FAIL'}")
            print(f"Check Duration: {self.results['check_duration']:.2f} seconds")
            print()
            
            for check_name, check_result in self.results['checks'].items():
                status_icon = "✅" if check_result['status'] == 'PASS' else "❌"
                print(f"{status_icon} {check_name.upper()}: {check_result['status']}")
                
                if check_result['errors']:
                    for error in check_result['errors']:
                        print(f"   ❌ {error}")
                
                if self.verbose and check_result['details']:
                    for key, value in check_result['details'].items():
                        print(f"   📊 {key}: {value}")
                print()


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(description="ProductionWildTypeAligner Health Check")
    parser.add_argument('--full', action='store_true', help='Run full health check including functionality tests')
    parser.add_argument('--json', action='store_true', help='Output results in JSON format')
    parser.add_argument('--quiet', action='store_true', help='Suppress verbose output')
    
    args = parser.parse_args()
    
    # Initialize health checker
    health_checker = HealthChecker(
        verbose=not args.quiet,
        json_output=args.json
    )
    
    # Run health checks
    results = health_checker.run_all_checks(full_check=args.full)
    
    # Print results
    health_checker.print_results()
    
    # Exit with appropriate code
    exit_code = 0 if results['overall_status'] == 'PASS' else 1
    sys.exit(exit_code)


if __name__ == '__main__':
    main()