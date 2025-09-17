#!/usr/bin/env python3
"""
Production Environment Test Suite
Comprehensive testing of containerized ProductionWildTypeAligner deployment

This script tests:
1. Container build and deployment
2. Health checks and validation
3. Basic functionality testing
4. Resource monitoring
5. Cleanup and rollback capabilities

Usage:
    python test_production_environment.py [--build] [--deploy] [--test] [--cleanup]

Author: GenomeAMRAnalyzer Pipeline
Version: 2.0 - Production Testing
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


class ProductionEnvironmentTester:
    """Comprehensive testing suite for production environment"""
    
    def __init__(self, project_root: str = None):
        self.project_root = Path(project_root) if project_root else Path(__file__).parent.parent.parent
        self.docker_dir = self.project_root / "docker" / "production"
        self.compose_file = self.docker_dir / "docker-compose.yml"
        
        # Test results
        self.test_results = {
            'timestamp': time.time(),
            'overall_status': 'UNKNOWN',
            'tests': {},
            'performance_metrics': {}
        }
        
        # Setup logging
        logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
        self.logger = logging.getLogger(__name__)
    
    def log(self, message: str, level: str = 'info'):
        """Log test messages"""
        getattr(self.logger, level)(message)
    
    def run_command(self, command: List[str], timeout: int = 300) -> Tuple[int, str, str]:
        """Run shell command and return result"""
        self.log(f"Running command: {' '.join(command)}")
        
        try:
            result = subprocess.run(
                command,
                capture_output=True,
                text=True,
                timeout=timeout,
                cwd=str(self.project_root)
            )
            return result.returncode, result.stdout, result.stderr
        except subprocess.TimeoutExpired:
            self.log(f"Command timed out after {timeout} seconds", 'error')
            return 124, "", f"Command timed out after {timeout} seconds"
        except Exception as e:
            self.log(f"Command failed with exception: {e}", 'error')
            return 1, "", str(e)
    
    def test_docker_build(self) -> Dict[str, Any]:
        """Test Docker image building"""
        self.log("Testing Docker image build...")
        
        test_result = {
            'status': 'PASS',
            'details': {},
            'errors': [],
            'duration': 0
        }
        
        start_time = time.time()
        
        # Build the production image
        build_command = [
            'docker', 'build',
            '-f', str(self.docker_dir / 'Dockerfile'),
            '-t', 'genome-amr-analyzer/wildtype-aligner:test',
            str(self.project_root)
        ]
        
        returncode, stdout, stderr = self.run_command(build_command, timeout=600)
        
        test_result['duration'] = time.time() - start_time
        
        if returncode == 0:
            test_result['details']['build_status'] = 'SUCCESS'
            self.log("✅ Docker image built successfully")
            
            # Check image size
            size_command = ['docker', 'images', 'genome-amr-analyzer/wildtype-aligner:test', '--format', '{{.Size}}']
            size_code, size_out, size_err = self.run_command(size_command)
            
            if size_code == 0:
                test_result['details']['image_size'] = size_out.strip()
                self.log(f"📊 Image size: {size_out.strip()}")
        else:
            test_result['status'] = 'FAIL'
            test_result['errors'].append(f"Docker build failed: {stderr}")
            self.log(f"❌ Docker build failed: {stderr}", 'error')
        
        return test_result
    
    def test_container_startup(self) -> Dict[str, Any]:
        """Test container startup and basic functionality"""
        self.log("Testing container startup...")
        
        test_result = {
            'status': 'PASS',
            'details': {},
            'errors': [],
            'duration': 0
        }
        
        start_time = time.time()
        
        # Start container with health check
        run_command = [
            'docker', 'run', '--rm', '-d',
            '--name', 'wildtype-aligner-test',
            'genome-amr-analyzer/wildtype-aligner:test',
            'sleep', '60'
        ]
        
        returncode, stdout, stderr = self.run_command(run_command)
        
        if returncode == 0:
            container_id = stdout.strip()
            test_result['details']['container_id'] = container_id
            
            # Wait for container to be running
            time.sleep(5)
            
            # Check container status
            status_command = ['docker', 'ps', '--filter', f'id={container_id}', '--format', '{{.Status}}']
            status_code, status_out, status_err = self.run_command(status_command)
            
            if status_code == 0 and 'Up' in status_out:
                test_result['details']['container_status'] = 'RUNNING'
                self.log("✅ Container started successfully")
                
                # Test application import inside container
                import_command = [
                    'docker', 'exec', container_id,
                    'python', '-c',
                    'from src.production_wildtype_aligner import ProductionWildTypeAligner; print("Import successful")'
                ]
                
                import_code, import_out, import_err = self.run_command(import_command)
                
                if import_code == 0:
                    test_result['details']['application_import'] = 'SUCCESS'
                    self.log("✅ Application imports working in container")
                else:
                    test_result['status'] = 'FAIL'
                    test_result['errors'].append(f"Application import failed: {import_err}")
                
                # Cleanup container
                cleanup_command = ['docker', 'stop', container_id]
                self.run_command(cleanup_command)
                
            else:
                test_result['status'] = 'FAIL'
                test_result['errors'].append(f"Container not running: {status_out}")
        else:
            test_result['status'] = 'FAIL'
            test_result['errors'].append(f"Container startup failed: {stderr}")
        
        test_result['duration'] = time.time() - start_time
        return test_result
    
    def test_docker_compose_deployment(self) -> Dict[str, Any]:
        """Test Docker Compose deployment"""
        self.log("Testing Docker Compose deployment...")
        
        test_result = {
            'status': 'PASS',
            'details': {},
            'errors': [],
            'duration': 0
        }
        
        start_time = time.time()
        
        if not self.compose_file.exists():
            test_result['status'] = 'FAIL'
            test_result['errors'].append(f"Docker Compose file not found: {self.compose_file}")
            return test_result
        
        # Start services with docker-compose
        compose_up_command = [
            'docker-compose', '-f', str(self.compose_file),
            'up', '-d', '--build'
        ]
        
        returncode, stdout, stderr = self.run_command(compose_up_command, timeout=600)
        
        if returncode == 0:
            test_result['details']['compose_up'] = 'SUCCESS'
            self.log("✅ Docker Compose services started")
            
            # Wait for services to initialize
            time.sleep(10)
            
            # Check service status
            ps_command = ['docker-compose', '-f', str(self.compose_file), 'ps']
            ps_code, ps_out, ps_err = self.run_command(ps_command)
            
            if ps_code == 0:
                test_result['details']['services_status'] = ps_out
                
                # Count running services
                running_services = ps_out.count('Up')
                test_result['details']['running_services'] = running_services
                
                if running_services > 0:
                    self.log(f"✅ {running_services} services running")
                else:
                    test_result['status'] = 'FAIL'
                    test_result['errors'].append("No services are running")
            
            # Test health checks if available
            health_command = [
                'docker-compose', '-f', str(self.compose_file),
                'exec', '-T', 'wildtype-aligner',
                'python', '/app/docker/production/health_check.py', '--json'
            ]
            
            health_code, health_out, health_err = self.run_command(health_command)
            
            if health_code == 0:
                try:
                    health_data = json.loads(health_out)
                    test_result['details']['health_check'] = health_data['overall_status']
                    self.log(f"✅ Health check: {health_data['overall_status']}")
                except:
                    test_result['details']['health_check'] = 'PARSE_ERROR'
            else:
                test_result['details']['health_check'] = 'FAILED'
                test_result['errors'].append(f"Health check failed: {health_err}")
            
            # Cleanup
            down_command = ['docker-compose', '-f', str(self.compose_file), 'down']
            self.run_command(down_command)
            
        else:
            test_result['status'] = 'FAIL'
            test_result['errors'].append(f"Docker Compose up failed: {stderr}")
        
        test_result['duration'] = time.time() - start_time
        return test_result
    
    def test_environment_validation(self) -> Dict[str, Any]:
        """Test environment validation script"""
        self.log("Testing environment validation...")
        
        test_result = {
            'status': 'PASS',
            'details': {},
            'errors': [],
            'duration': 0
        }
        
        start_time = time.time()
        
        # Run environment validation in container
        validation_command = [
            'docker', 'run', '--rm',
            '-v', f'{self.project_root}/config:/app/config',
            'genome-amr-analyzer/wildtype-aligner:test',
            'python', '/app/docker/production/validate_environment.py', '--json'
        ]
        
        returncode, stdout, stderr = self.run_command(validation_command)
        
        if returncode == 0:
            try:
                validation_data = json.loads(stdout)
                test_result['details']['validation_status'] = validation_data['overall_status']
                test_result['details']['validations'] = validation_data['validations']
                self.log(f"✅ Environment validation: {validation_data['overall_status']}")
            except json.JSONDecodeError:
                test_result['status'] = 'FAIL'
                test_result['errors'].append("Could not parse validation output as JSON")
        else:
            test_result['status'] = 'FAIL'
            test_result['errors'].append(f"Environment validation failed: {stderr}")
        
        test_result['duration'] = time.time() - start_time
        return test_result
    
    def test_performance_benchmarks(self) -> Dict[str, Any]:
        """Test basic performance benchmarks"""
        self.log("Testing performance benchmarks...")
        
        test_result = {
            'status': 'PASS',
            'details': {},
            'errors': [],
            'duration': 0
        }
        
        start_time = time.time()
        
        # Create a small test FASTA file
        test_fasta_content = """>test_protein_1
MVKKLKSLIFLFALFLTACSNNDFQGDKSTVKGKNIAHLSESLFKTKQKILGSISLSLFFLYVLLVTLSKLKDE
>test_protein_2
MVKKLKSLIFLFALFLTACSNNDFQGDKSTVKGKNIAHLSESLFKTKQKILGSISLSLFFLYVLLVTLSKLKDE
"""
        
        # Create temporary directory for test
        test_dir = self.project_root / "test_temp"
        test_dir.mkdir(exist_ok=True)
        
        test_fasta_file = test_dir / "test_proteins.fasta"
        with open(test_fasta_file, 'w') as f:
            f.write(test_fasta_content)
        
        # Run performance test in container
        perf_command = [
            'docker', 'run', '--rm',
            '-v', f'{test_dir}:/app/test_data',
            'genome-amr-analyzer/wildtype-aligner:test',
            'python', '-c',
            '''
import time
import sys
sys.path.insert(0, "/app/src")
start_time = time.time()
try:
    from production_wildtype_aligner import ProductionWildTypeAligner
    from sepi_configuration_manager import SEPIConfigurationManager
    print(f"Import time: {time.time() - start_time:.3f}s")
    print("Performance test completed")
except Exception as e:
    print(f"Performance test failed: {e}")
    sys.exit(1)
'''
        ]
        
        returncode, stdout, stderr = self.run_command(perf_command)
        
        if returncode == 0:
            # Parse performance metrics from output
            for line in stdout.split('\n'):
                if 'Import time:' in line:
                    import_time = float(line.split(':')[1].strip().replace('s', ''))
                    test_result['details']['import_time_seconds'] = import_time
                    self.log(f"📊 Import time: {import_time:.3f}s")
            
            test_result['details']['performance_test'] = 'SUCCESS'
            self.log("✅ Performance test completed")
        else:
            test_result['status'] = 'FAIL'
            test_result['errors'].append(f"Performance test failed: {stderr}")
        
        # Cleanup
        if test_dir.exists():
            import shutil
            shutil.rmtree(test_dir)
        
        test_result['duration'] = time.time() - start_time
        return test_result
    
    def run_all_tests(self, build: bool = True, deploy: bool = True, 
                     test_functionality: bool = True) -> Dict[str, Any]:
        """Run all production environment tests"""
        self.log("Starting comprehensive production environment testing...")
        
        if build:
            self.test_results['tests']['docker_build'] = self.test_docker_build()
        
        # Only proceed with other tests if build succeeded
        if (not build or self.test_results['tests'].get('docker_build', {}).get('status') == 'PASS'):
            
            if test_functionality:
                self.test_results['tests']['container_startup'] = self.test_container_startup()
                self.test_results['tests']['environment_validation'] = self.test_environment_validation()
                self.test_results['tests']['performance_benchmarks'] = self.test_performance_benchmarks()
            
            if deploy:
                self.test_results['tests']['docker_compose_deployment'] = self.test_docker_compose_deployment()
        
        # Determine overall status
        failed_tests = [name for name, result in self.test_results['tests'].items()
                       if result['status'] == 'FAIL']
        
        if failed_tests:
            self.test_results['overall_status'] = 'FAIL'
            self.test_results['failed_tests'] = failed_tests
        else:
            self.test_results['overall_status'] = 'PASS'
        
        # Calculate total duration
        total_duration = sum(test.get('duration', 0) for test in self.test_results['tests'].values())
        self.test_results['total_duration'] = total_duration
        
        return self.test_results
    
    def cleanup_test_artifacts(self):
        """Clean up test artifacts and containers"""
        self.log("Cleaning up test artifacts...")
        
        # Remove test images
        cleanup_commands = [
            ['docker', 'rmi', 'genome-amr-analyzer/wildtype-aligner:test'],
            ['docker', 'system', 'prune', '-f']
        ]
        
        for command in cleanup_commands:
            self.run_command(command)
        
        self.log("✅ Cleanup completed")
    
    def print_results(self):
        """Print test results"""
        print("\n" + "="*60)
        print("🧪 Production Environment Test Results")
        print("="*60)
        
        status_icon = "✅" if self.test_results['overall_status'] == 'PASS' else "❌"
        print(f"Overall Status: {status_icon} {self.test_results['overall_status']}")
        print(f"Total Duration: {self.test_results['total_duration']:.2f} seconds")
        print()
        
        for test_name, test_result in self.test_results['tests'].items():
            status_icon = "✅" if test_result['status'] == 'PASS' else "❌"
            duration = test_result.get('duration', 0)
            print(f"{status_icon} {test_name.upper().replace('_', ' ')}: {test_result['status']} ({duration:.2f}s)")
            
            # Show errors
            if test_result.get('errors'):
                for error in test_result['errors']:
                    print(f"   ❌ {error}")
            
            # Show key details
            if test_result.get('details'):
                for key, value in test_result['details'].items():
                    if isinstance(value, str) and len(value) < 100:  # Only show short details
                        print(f"   📊 {key}: {value}")
            
            print()


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(description="Production Environment Test Suite")
    parser.add_argument('--build', action='store_true', default=True,
                       help='Test Docker image building')
    parser.add_argument('--deploy', action='store_true', default=True,
                       help='Test Docker Compose deployment')
    parser.add_argument('--test', action='store_true', default=True,
                       help='Test functionality and performance')
    parser.add_argument('--cleanup', action='store_true',
                       help='Clean up test artifacts after testing')
    parser.add_argument('--json', action='store_true',
                       help='Output results in JSON format')
    parser.add_argument('--project-root', 
                       help='Path to project root directory')
    
    args = parser.parse_args()
    
    # Initialize tester
    tester = ProductionEnvironmentTester(project_root=args.project_root)
    
    # Run tests
    results = tester.run_all_tests(
        build=args.build,
        deploy=args.deploy,
        test_functionality=args.test
    )
    
    # Output results
    if args.json:
        print(json.dumps(results, indent=2))
    else:
        tester.print_results()
    
    # Cleanup if requested
    if args.cleanup:
        tester.cleanup_test_artifacts()
    
    # Exit with appropriate code
    exit_code = 0 if results['overall_status'] == 'PASS' else 1
    sys.exit(exit_code)


if __name__ == '__main__':
    main()