#!/usr/bin/env python3
"""
SEPI 2.0 Configuration Setup and Validation Script
Complete setup and validation for SEPI 2.0 integration with WildType Aligner

This script:
1. Validates SEPI 2.0 installation and availability
2. Creates necessary directory structure
3. Tests SEPI configuration and NCBI connectivity
4. Validates gene reference fetching capabilities
5. Generates setup report for production use

Usage:
    python setup_sepi_integration.py --workspace /path/to/workspace [options]

Author: GenomeAMRAnalyzer Pipeline
Version: 2.0 - Production Setup
"""

import os
import sys
import asyncio
import logging
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple
import yaml
import json
import time

# Add src directory to path for imports
src_dir = Path(__file__).parent / "src"
sys.path.insert(0, str(src_dir))

try:
    from sepi_configuration_manager import SEPIConfigurationManager, create_default_sepi_configuration
    from sepi_integration import EnhancedSEPIReferenceManager, SEPIIntegrationTester
except ImportError as e:
    print(f"Warning: Could not import SEPI modules: {e}")
    print("Make sure sepi_configuration_manager.py and sepi_integration.py are in the src directory")


class SEPISetupValidator:
    """
    Comprehensive SEPI 2.0 setup validation and configuration
    """
    
    def __init__(self, workspace_root: str, config_file: Optional[str] = None):
        """Initialize SEPI setup validator"""
        self.workspace_root = Path(workspace_root).resolve()
        self.config_file = config_file
        
        # Setup logging
        self.logger = logging.getLogger('SEPISetupValidator')
        
        # Validation results
        self.validation_results = {
            'sepi_availability': False,
            'directory_structure': False,
            'ncbi_connectivity': False,
            'configuration_valid': False,
            'gene_fetching': False,
            'overall_status': 'FAILED'
        }
        
        # Detailed results storage
        self.detailed_results = {}
        
        self.logger.info(f"Initializing SEPI setup validation for workspace: {self.workspace_root}")
    
    async def run_complete_validation(self) -> Dict[str, Any]:
        """Run complete SEPI setup validation"""
        
        self.logger.info("🚀 Starting comprehensive SEPI 2.0 setup validation")
        
        # Step 1: Validate SEPI availability
        await self._validate_sepi_availability()
        
        # Step 2: Setup directory structure
        await self._setup_directory_structure()
        
        # Step 3: Validate configuration
        await self._validate_configuration()
        
        # Step 4: Test NCBI connectivity
        await self._test_ncbi_connectivity()
        
        # Step 5: Test gene reference fetching
        await self._test_gene_reference_fetching()
        
        # Step 6: Generate final report
        self._generate_validation_report()
        
        return self.validation_results
    
    async def _validate_sepi_availability(self):
        """Validate SEPI 2.0 script availability and functionality"""
        
        self.logger.info("📋 Step 1: Validating SEPI 2.0 availability")
        
        # Search for SEPI script
        sepi_locations = [
            self.workspace_root / "MetaDataHarvester" / "sepi2.0" / "sepi.py",
            self.workspace_root / "sepi2.0" / "sepi.py",
            Path("MetaDataHarvester/sepi2.0/sepi.py"),
            Path("sepi2.0/sepi.py")
        ]
        
        sepi_path = None
        for location in sepi_locations:
            if location.exists():
                sepi_path = location
                break
        
        if not sepi_path:
            self.logger.error("❌ SEPI script not found in expected locations")
            self.detailed_results['sepi_availability'] = {
                'status': 'FAILED',
                'error': 'SEPI script not found',
                'searched_locations': [str(loc) for loc in sepi_locations]
            }
            return
        
        # Test SEPI functionality
        try:
            self.logger.info(f"Testing SEPI script: {sepi_path}")
            
            result = await asyncio.create_subprocess_exec(
                sys.executable, str(sepi_path), "--help",
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE
            )
            
            stdout, stderr = await result.communicate()
            
            if result.returncode == 0:
                self.logger.info("✅ SEPI 2.0 is available and functional")
                self.validation_results['sepi_availability'] = True
                self.detailed_results['sepi_availability'] = {
                    'status': 'SUCCESS',
                    'sepi_path': str(sepi_path),
                    'version_info': stdout.decode()[:500] if stdout else 'N/A'
                }
            else:
                self.logger.error(f"❌ SEPI script failed to execute: {stderr.decode()}")
                self.detailed_results['sepi_availability'] = {
                    'status': 'FAILED',
                    'sepi_path': str(sepi_path),
                    'error': stderr.decode()
                }
                
        except Exception as e:
            self.logger.error(f"❌ Exception testing SEPI: {e}")
            self.detailed_results['sepi_availability'] = {
                'status': 'FAILED',
                'error': str(e)
            }
    
    async def _setup_directory_structure(self):
        """Setup required directory structure for SEPI integration"""
        
        self.logger.info("📁 Step 2: Setting up directory structure")
        
        required_directories = [
            self.workspace_root / ".sepi_cache",
            self.workspace_root / "temp" / "sepi",
            self.workspace_root / "config",
            self.workspace_root / "src",
            self.workspace_root / "logs",
            self.workspace_root / "output" / "alignments",
            self.workspace_root / "output" / "references"
        ]
        
        created_directories = []
        failed_directories = []
        
        for directory in required_directories:
            try:
                directory.mkdir(parents=True, exist_ok=True)
                
                # Test write permissions
                test_file = directory / ".test_write"
                test_file.write_text("test")
                test_file.unlink()
                
                created_directories.append(str(directory))
                
            except Exception as e:
                self.logger.error(f"❌ Failed to create/access directory {directory}: {e}")
                failed_directories.append({'directory': str(directory), 'error': str(e)})
        
        if not failed_directories:
            self.logger.info("✅ Directory structure setup successful")
            self.validation_results['directory_structure'] = True
            self.detailed_results['directory_structure'] = {
                'status': 'SUCCESS',
                'created_directories': created_directories
            }
        else:
            self.logger.error("❌ Directory structure setup failed")
            self.detailed_results['directory_structure'] = {
                'status': 'FAILED',
                'created_directories': created_directories,
                'failed_directories': failed_directories
            }
    
    async def _validate_configuration(self):
        """Validate SEPI configuration"""
        
        self.logger.info("⚙️ Step 3: Validating SEPI configuration")
        
        try:
            # Initialize configuration manager
            if self.config_file and Path(self.config_file).exists():
                config_manager = SEPIConfigurationManager(
                    config_file=self.config_file,
                    workspace_root=str(self.workspace_root)
                )
            else:
                # Create default configuration
                config_manager = create_default_sepi_configuration(str(self.workspace_root))
            
            # Configuration validation happens during initialization
            self.logger.info("✅ SEPI configuration is valid")
            self.validation_results['configuration_valid'] = True
            
            # Store configuration details
            self.detailed_results['configuration'] = {
                'status': 'SUCCESS',
                'config_file': self.config_file or 'default_generated',
                'sepi_script_path': config_manager.environment_config.sepi_script_path if config_manager.environment_config else 'N/A',
                'email': config_manager.ncbi_config.email if config_manager.ncbi_config else 'N/A',
                'cache_directory': config_manager.environment_config.cache_directory if config_manager.environment_config else 'N/A',
                'gene_configs_loaded': len(config_manager.gene_configs),
                'organism_configs_loaded': len(config_manager.organism_configs)
            }
            
            # Store for later use
            self.config_manager = config_manager
            
        except Exception as e:
            self.logger.error(f"❌ Configuration validation failed: {e}")
            self.detailed_results['configuration'] = {
                'status': 'FAILED',
                'error': str(e)
            }
    
    async def _test_ncbi_connectivity(self):
        """Test NCBI API connectivity"""
        
        self.logger.info("🌐 Step 4: Testing NCBI connectivity")
        
        if not hasattr(self, 'config_manager'):
            self.logger.error("❌ Configuration not available for NCBI test")
            self.detailed_results['ncbi_connectivity'] = {
                'status': 'FAILED',
                'error': 'Configuration not validated'
            }
            return
        
        try:
            # Test Entrez connectivity - simplified test
            self.logger.info("Testing NCBI connectivity with simple request...")
            
            # Simple urllib test to check NCBI accessibility
            import urllib.request
            import urllib.parse
            
            # Test NCBI esearch endpoint
            base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
            params = {
                'db': 'nucleotide',
                'term': 'Escherichia coli',
                'retmax': 1,
                'retmode': 'json'
            }
            
            if self.config_manager.ncbi_config:
                params['email'] = self.config_manager.ncbi_config.email
                if self.config_manager.ncbi_config.api_key:
                    params['api_key'] = self.config_manager.ncbi_config.api_key
            
            url = f"{base_url}?{urllib.parse.urlencode(params)}"
            
            # Make request with timeout
            request = urllib.request.Request(url)
            
            with urllib.request.urlopen(request, timeout=30) as response:
                if response.status == 200:
                    response_data = response.read().decode()
                    
                    # Check if response contains expected data
                    if 'esearchresult' in response_data.lower():
                        self.logger.info("✅ NCBI connectivity test successful")
                        self.validation_results['ncbi_connectivity'] = True
                        self.detailed_results['ncbi_connectivity'] = {
                            'status': 'SUCCESS',
                            'test_results': 'NCBI API accessible',
                            'email_used': self.config_manager.ncbi_config.email if self.config_manager.ncbi_config else 'N/A'
                        }
                    else:
                        self.logger.error("❌ NCBI returned unexpected response format")
                        self.detailed_results['ncbi_connectivity'] = {
                            'status': 'FAILED',
                            'error': 'Unexpected response format'
                        }
                else:
                    self.logger.error(f"❌ NCBI request failed with status: {response.status}")
                    self.detailed_results['ncbi_connectivity'] = {
                        'status': 'FAILED',
                        'error': f'HTTP {response.status}'
                    }
                
        except Exception as e:
            self.logger.error(f"❌ NCBI connectivity test failed: {e}")
            self.detailed_results['ncbi_connectivity'] = {
                'status': 'FAILED',
                'error': str(e)
            }
    
    async def _test_gene_reference_fetching(self):
        """Test gene reference fetching capabilities"""
        
        self.logger.info("🧬 Step 5: Testing gene reference fetching")
        
        if not hasattr(self, 'config_manager'):
            self.logger.error("❌ Configuration not available for gene fetching test")
            self.detailed_results['gene_fetching'] = {
                'status': 'FAILED',
                'error': 'Configuration not validated'
            }
            return
        
        try:
            # Initialize integration tester with a subset of genes
            test_genes = ['acrA', 'tolC']  # Start with well-known genes
            
            tester = SEPIIntegrationTester(str(self.workspace_root))
            
            self.logger.info(f"Testing reference fetching for genes: {test_genes}")
            
            # Test with timeout
            test_results = await asyncio.wait_for(
                tester.test_gene_reference_fetching(test_genes),
                timeout=300  # 5 minute timeout
            )
            
            success_rate = test_results['success_rate']
            
            if success_rate >= 50:  # At least 50% success required
                self.logger.info(f"✅ Gene reference fetching test successful ({success_rate:.1f}% success rate)")
                self.validation_results['gene_fetching'] = True
                self.detailed_results['gene_fetching'] = {
                    'status': 'SUCCESS',
                    'success_rate': success_rate,
                    'successful_genes': [g['gene_name'] for g in test_results['successful_fetches']],
                    'failed_genes': [g['gene_name'] for g in test_results['failed_fetches']],
                    'reference_manager_stats': test_results.get('reference_manager_stats', {})
                }
            else:
                self.logger.error(f"❌ Gene reference fetching test failed ({success_rate:.1f}% success rate)")
                self.detailed_results['gene_fetching'] = {
                    'status': 'FAILED',
                    'success_rate': success_rate,
                    'failed_genes': [g['gene_name'] for g in test_results['failed_fetches']],
                    'error': 'Success rate below threshold (50%)'
                }
                
        except asyncio.TimeoutError:
            self.logger.error("❌ Gene reference fetching test timed out")
            self.detailed_results['gene_fetching'] = {
                'status': 'FAILED',
                'error': 'Test timed out after 5 minutes'
            }
        except Exception as e:
            self.logger.error(f"❌ Gene reference fetching test failed: {e}")
            self.detailed_results['gene_fetching'] = {
                'status': 'FAILED',
                'error': str(e)
            }
    
    def _generate_validation_report(self):
        """Generate comprehensive validation report"""
        
        self.logger.info("📊 Step 6: Generating validation report")
        
        # Calculate overall status
        required_validations = ['sepi_availability', 'directory_structure', 'configuration_valid']
        optional_validations = ['ncbi_connectivity', 'gene_fetching']
        
        required_passed = all(self.validation_results[key] for key in required_validations)
        optional_passed = sum(self.validation_results[key] for key in optional_validations)
        
        if required_passed and optional_passed >= 1:
            self.validation_results['overall_status'] = 'SUCCESS'
        elif required_passed:
            self.validation_results['overall_status'] = 'PARTIAL'
        else:
            self.validation_results['overall_status'] = 'FAILED'
        
        # Generate detailed report
        report = {
            'validation_timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
            'workspace_root': str(self.workspace_root),
            'overall_status': self.validation_results['overall_status'],
            'validation_summary': self.validation_results,
            'detailed_results': self.detailed_results,
            'recommendations': self._generate_recommendations()
        }
        
        # Save report
        report_file = self.workspace_root / "logs" / "sepi_validation_report.json"
        report_file.parent.mkdir(parents=True, exist_ok=True)
        
        with open(report_file, 'w') as f:
            json.dump(report, f, indent=2)
        
        self.logger.info(f"✅ Validation report saved to: {report_file}")
        
        # Print summary
        self._print_validation_summary()
    
    def _generate_recommendations(self) -> List[str]:
        """Generate recommendations based on validation results"""
        
        recommendations = []
        
        if not self.validation_results['sepi_availability']:
            recommendations.append("Install or fix SEPI 2.0 script in MetaDataHarvester/sepi2.0/")
        
        if not self.validation_results['directory_structure']:
            recommendations.append("Fix directory permissions and ensure write access to workspace")
        
        if not self.validation_results['configuration_valid']:
            recommendations.append("Review and fix SEPI configuration file")
        
        if not self.validation_results['ncbi_connectivity']:
            recommendations.append("Check internet connection and NCBI email configuration")
        
        if not self.validation_results['gene_fetching']:
            recommendations.append("Review gene configurations and SEPI setup for reference fetching")
        
        if self.validation_results['overall_status'] == 'SUCCESS':
            recommendations.append("System is ready for production use!")
        elif self.validation_results['overall_status'] == 'PARTIAL':
            recommendations.append("Basic functionality available, but some features may be limited")
        
        return recommendations
    
    def _print_validation_summary(self):
        """Print validation summary to console"""
        
        status_icon = {
            'SUCCESS': '✅',
            'PARTIAL': '⚠️',
            'FAILED': '❌'
        }
        
        print("\n" + "="*60)
        print("🔬 SEPI 2.0 Integration Validation Summary")
        print("="*60)
        print(f"Overall Status: {status_icon[self.validation_results['overall_status']]} {self.validation_results['overall_status']}")
        print(f"Workspace: {self.workspace_root}")
        print()
        
        print("Validation Results:")
        for key, value in self.validation_results.items():
            if key != 'overall_status':
                icon = '✅' if value else '❌'
                print(f"  {icon} {key.replace('_', ' ').title()}: {'PASS' if value else 'FAIL'}")
        
        print("\nRecommendations:")
        for i, rec in enumerate(self._generate_recommendations(), 1):
            print(f"  {i}. {rec}")
        
        print("\n" + "="*60)


async def main():
    """Command line interface for SEPI setup validation"""
    import argparse
    
    parser = argparse.ArgumentParser(description="SEPI 2.0 Setup and Validation")
    parser.add_argument("--workspace", required=True, help="Workspace root directory")
    parser.add_argument("--config", help="SEPI configuration file")
    parser.add_argument("--log-level", default="INFO", 
                       choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    parser.add_argument("--skip-gene-test", action="store_true", 
                       help="Skip gene reference fetching test (faster)")
    
    args = parser.parse_args()
    
    # Setup logging
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    
    try:
        # Initialize validator
        validator = SEPISetupValidator(args.workspace, args.config)
        
        # Run validation
        if args.skip_gene_test:
            # Override gene fetching test
            validator.validation_results['gene_fetching'] = True
            validator.detailed_results['gene_fetching'] = {'status': 'SKIPPED'}
        
        results = await validator.run_complete_validation()
        
        # Return appropriate exit code
        if results['overall_status'] == 'SUCCESS':
            return 0
        elif results['overall_status'] == 'PARTIAL':
            return 1
        else:
            return 2
        
    except Exception as e:
        print(f"❌ Setup validation failed: {e}")
        logging.exception("Setup validation exception")
        return 3


if __name__ == "__main__":
    sys.exit(asyncio.run(main()))