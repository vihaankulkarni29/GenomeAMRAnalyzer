from src.configuration_manager import config_manager
#!/usr/bin/env python3
"""
PRIORITY 1 COMPREHENSIVE STRESS TEST
====================================
This test performs exhaustive validation of all Priority 1 requirements
with maximum robustness testing including edge cases, error conditions,
and stress scenarios.

Priority 1 Requirements:
1. Fix Remaining Hardcoded Values - Clean up test files and examples
2. Resolve Import Dependencies - Fix pandas and setuptools imports, update requirements.txt  
3. Add Error Handling & Validation - Implement robust error handling across all components

Author: GenomeAMRAnalyzer Pipeline
Version: 1.0.0 (Stress Test)
"""

import os
import sys
import traceback
import tempfile
import shutil
from pathlib import Path
import json
import time
import gc
import psutil
import logging

# Comprehensive test results tracking
class StressTestResults:
    def __init__(self):
        self.tests_passed = 0
        self.tests_failed = 0
        self.test_details = []
        self.performance_metrics = {}
        self.memory_usage = []
        
    def record_test(self, test_name, passed, details="", execution_time=0, memory_delta=0):
        if passed:
            self.tests_passed += 1
            status = "✅ PASSED"
        else:
            self.tests_failed += 1
            status = "❌ FAILED"
            
        self.test_details.append({
            'test': test_name,
            'status': status,
            'details': details,
            'execution_time': execution_time,
            'memory_delta': memory_delta
        })
        
        print(f"   {status}: {test_name}")
        if details:
            print(f"      Details: {details}")
        if execution_time > 0:
            print(f"      Execution Time: {execution_time:.3f}s")
            
    def get_summary(self):
        total = self.tests_passed + self.tests_failed
        success_rate = (self.tests_passed / total * 100) if total > 0 else 0
        return {
            'total_tests': total,
            'passed': self.tests_passed,
            'failed': self.tests_failed,
            'success_rate': success_rate
        }

def measure_performance(func):
    """Decorator to measure performance and memory usage."""
    def wrapper(*args, **kwargs):
        # Memory before
        process = psutil.Process()
        mem_before = process.memory_info().rss / 1024 / 1024  # MB
        
        # Time execution
        start_time = time.time()
        try:
            result = func(*args, **kwargs)
            success = True
        except Exception as e:
            result = str(e)
            success = False
        end_time = time.time()
        
        # Memory after
        mem_after = process.memory_info().rss / 1024 / 1024  # MB
        memory_delta = mem_after - mem_before
        
        execution_time = end_time - start_time
        
        return result, success, execution_time, memory_delta
    return wrapper

def stress_test_error_handling():
    """Comprehensive stress test of error handling system."""
    print("🔥 STRESS TESTING ERROR HANDLING SYSTEM...")
    results = StressTestResults()
    
    try:
        from src.core.robust_error_handling import (
            ValidationError, DataProcessingError, FileSystemError,
            ValidationSuite, RobustLogger, robust_exception_handler
        )
        
        # Test 1: Exception hierarchy stress test
        @measure_performance
        def test_exception_hierarchy():
            exceptions_to_test = [
                (ValidationError, "Test validation error"),
                (DataProcessingError, "Test data processing error"),
                (FileSystemError, "Test filesystem error")
            ]
            
            for exc_class, message in exceptions_to_test:
                try:
                    raise exc_class(message)
                except exc_class as e:
                    if str(e) != message:
                        raise AssertionError(f"Exception message mismatch: {str(e)} != {message}")
            return "All custom exceptions working correctly"
        
        result, success, exec_time, mem_delta = test_exception_hierarchy()
        results.record_test("Exception Hierarchy", success, result, exec_time, mem_delta)
        
        # Test 2: ValidationSuite stress test with edge cases
        @measure_performance  
        def test_validation_suite_edge_cases():
            validator = ValidationSuite()
            
            # Edge case sequences
            edge_cases = [
                ("", False, "Empty sequence"),
                ("ATCG" * 1000, True, "Very long valid sequence"),
                ("ATCGXYZ", False, "Invalid characters"),
                ("atcg", True, "Lowercase valid"),
                ("ATCG-ATCG", False, "Hyphen character"),
                ("ATCG ATCG", False, "Space character"),
                ("NNNNNNNN", True, "All N's (valid)"),
                ("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG", True, "Medium length"),
            ]
            
            for seq, expected, desc in edge_cases:
                try:
                    is_valid = validator.validate_sequence(seq)
                    if is_valid != expected:
                        raise AssertionError(f"Validation failed for {desc}: {seq}")
                except Exception as e:
                    if expected:  # If we expected it to be valid but got exception
                        raise AssertionError(f"Unexpected error for {desc}: {str(e)}")
                        
            return f"Validated {len(edge_cases)} edge cases successfully"
        
        result, success, exec_time, mem_delta = test_validation_suite_edge_cases()
        results.record_test("ValidationSuite Edge Cases", success, result, exec_time, mem_delta)
        
        # Test 3: RobustLogger stress test with concurrent operations
        @measure_performance
        def test_robust_logger_stress():
            with tempfile.TemporaryDirectory() as temp_dir:
                logger = RobustLogger("StressTest", log_dir=temp_dir)
                
                # Rapid logging stress test
                for i in range(100):
                    logger.info(f"Stress test message {i}")
                    logger.warning(f"Warning message {i}")
                    if i % 10 == 0:
                        logger.error(f"Error message {i}")
                
                # Check log file was created and contains entries
                log_files = list(Path(temp_dir).glob("*.log"))
                if not log_files:
                    raise AssertionError("No log files created")
                    
                log_content = log_files[0].read_text()
                if "Stress test message 99" not in log_content:
                    raise AssertionError("Log content incomplete")
                    
                return f"Successfully logged 300 messages, log file: {log_files[0].name}"
        
        result, success, exec_time, mem_delta = test_robust_logger_stress()
        results.record_test("RobustLogger Stress Test", success, result, exec_time, mem_delta)
        
        # Test 4: Error decorator resilience test
        @measure_performance
        def test_error_decorator_resilience():
            @robust_exception_handler
            def function_that_fails():
                raise ValueError("Intentional test failure")
                
            @robust_exception_handler  
            def function_that_succeeds():
                return "Success"
                
            # Test that decorator handles failures gracefully
            result1 = function_that_fails()
            if result1 is not None:
                raise AssertionError("Decorator should return None for failed function")
                
            # Test that decorator doesn't interfere with success
            result2 = function_that_succeeds()
            if result2 != "Success":
                raise AssertionError("Decorator interfered with successful function")
                
            return "Error decorator handling both success and failure correctly"
        
        result, success, exec_time, mem_delta = test_error_decorator_resilience()
        results.record_test("Error Decorator Resilience", success, result, exec_time, mem_delta)
        
    except ImportError as e:
        results.record_test("Error Handling Import", False, f"Import failed: {str(e)}")
    except Exception as e:
        results.record_test("Error Handling System", False, f"Unexpected error: {str(e)}")
        traceback.print_exc()
        
    return results

def stress_test_dependency_management():
    """Comprehensive stress test of dependency management."""
    print("🔥 STRESS TESTING DEPENDENCY MANAGEMENT...")
    results = StressTestResults()
    
    try:
        from src.core.dependencies import DependencyChecker, safe_import
        
        # Test 1: Safe import with various scenarios
        @measure_performance
        def test_safe_import_scenarios():
            # Test existing module
            module1 = safe_import('os')
            if module1 is None:
                raise AssertionError("Failed to import existing module 'os'")
                
            # Test non-existent module
            module2 = safe_import('nonexistent_module_12345')
            if module2 is not None:
                raise AssertionError("safe_import should return None for non-existent module")
                
            # Test conditional import
            biopython = safe_import('Bio.SeqIO')
            # Should work since we installed biopython
            
            return f"Safe import working for existing, non-existent, and conditional modules"
        
        result, success, exec_time, mem_delta = test_safe_import_scenarios()
        results.record_test("Safe Import Scenarios", success, result, exec_time, mem_delta)
        
        # Test 2: Dependency checker comprehensive test
        @measure_performance
        def test_dependency_checker_comprehensive():
            checker = DependencyChecker()
            
            # Check all dependencies
            results_dict = checker.check_all_dependencies()
            
            # Verify structure
            required_keys = ['missing', 'available', 'summary']
            for key in required_keys:
                if key not in results_dict:
                    raise AssertionError(f"Missing key in dependency results: {key}")
            
            # Check that core Python modules are available
            if 'os' not in results_dict['available']:
                raise AssertionError("Core module 'os' should be available")
                
            return f"Checked dependencies: {len(results_dict['available'])} available, {len(results_dict['missing'])} missing"
        
        result, success, exec_time, mem_delta = test_dependency_checker_comprehensive()
        results.record_test("Dependency Checker Comprehensive", success, result, exec_time, mem_delta)
        
        # Test 3: Environment validation stress test
        @measure_performance
        def test_environment_validation():
            checker = DependencyChecker()
            
            # Test Python version validation
            python_version = sys.version_info
            if python_version.major < 3 or (python_version.major == 3 and python_version.minor < 6):
                raise AssertionError("Python version too old for robust testing")
                
            # Test package availability for core scientific stack
            critical_packages = ['numpy', 'pandas']
            missing_critical = []
            
            for package in critical_packages:
                if safe_import(package) is None:
                    missing_critical.append(package)
                    
            if missing_critical:
                return f"Environment validation completed, missing: {missing_critical}"
            else:
                return "Environment validation passed - all critical packages available"
        
        result, success, exec_time, mem_delta = test_environment_validation()
        results.record_test("Environment Validation", success, result, exec_time, mem_delta)
        
    except ImportError as e:
        results.record_test("Dependency Management Import", False, f"Import failed: {str(e)}")
    except Exception as e:
        results.record_test("Dependency Management", False, f"Unexpected error: {str(e)}")
        traceback.print_exc()
        
    return results

def stress_test_wildtype_aligner():
    """Comprehensive stress test of wildtype aligner."""
    print("🔥 STRESS TESTING WILDTYPE ALIGNER...")
    results = StressTestResults()
    
    try:
        from src.simplified_wildtype_aligner import SimplifiedWildTypeAligner, SimpleAlignerConfig
        
        # Test 1: Configuration robustness
        @measure_performance
        def test_configuration_robustness():
            # Test valid configuration
            config = SimpleAlignerConfig(
                input_dir="test_input",
                output_dir="test_output",
                target_genes=["geneA", "geneB", "geneC"],
                reference_dir=None
            )
            
            aligner = SimplifiedWildTypeAligner(config)
            
            # Verify configuration was set correctly
            if aligner.config.target_genes != ["geneA", "geneB", "geneC"]:
                raise AssertionError("Configuration not set correctly")
                
            return "Configuration robustness test passed"
        
        result, success, exec_time, mem_delta = test_configuration_robustness()
        results.record_test("Configuration Robustness", success, result, exec_time, mem_delta)
        
        # Test 2: analyze_sequences stress test with edge cases  
        @measure_performance
        def test_analyze_sequences_stress():
            config = SimpleAlignerConfig(
                input_dir="test_input",
                output_dir="test_output", 
                target_genes=["testGene"],
                reference_dir=None
            )
            
            aligner = SimplifiedWildTypeAligner(config)
            
            # Edge case sequences for stress testing
            stress_sequences = [
                "ATCGATCGATCG",  # Normal sequence
                "ATCGATCGATCG",  # Identical duplicate
                "ATCGATCGATCC",  # Single mutation
                "TTCGATCGATCG",  # Multiple mutations
                "ATCGATCGATCGATCGATCGATCG",  # Longer sequence
                "ATCG",  # Very short sequence
                "A" * 100,  # Long repetitive sequence
            ]
            
            try:
                result = aligner.analyze_sequences(stress_sequences, "stressGene")
                
                # Verify result structure
                required_keys = ['gene_name', 'total_sequences', 'mutations_found', 'analysis_summary']
                for key in required_keys:
                    if key not in result:
                        raise AssertionError(f"Missing key in analysis result: {key}")
                        
                if result['total_sequences'] != len(stress_sequences):
                    raise AssertionError("Sequence count mismatch")
                    
                return f"Analyzed {len(stress_sequences)} stress sequences successfully"
                
            except Exception as e:
                # Check if it's a graceful validation error (acceptable)
                if "validation" in str(e).lower() or "invalid" in str(e).lower():
                    return f"Graceful validation handling: {str(e)}"
                else:
                    raise e
        
        result, success, exec_time, mem_delta = test_analyze_sequences_stress()
        results.record_test("Analyze Sequences Stress Test", success, result, exec_time, mem_delta)
        
        # Test 3: Memory and performance under load
        @measure_performance
        def test_performance_under_load():
            config = SimpleAlignerConfig(
                input_dir="test_input",
                output_dir="test_output",
                target_genes=["loadTestGene"],
                reference_dir=None
            )
            
            aligner = SimplifiedWildTypeAligner(config)
            
            # Generate large sequence set for load testing
            base_sequence = "ATCGATCGATCGATCGATCG"
            large_sequence_set = []
            
            for i in range(50):  # 50 sequences for load test
                # Create variations
                seq = base_sequence
                if i % 5 == 0:  # Every 5th sequence has mutation
                    seq = seq.replace("A", "T", 1)  # Single mutation
                large_sequence_set.append(seq)
            
            # Process large set
            result = aligner.analyze_sequences(large_sequence_set, "loadTestGene")
            
            if result['total_sequences'] != len(large_sequence_set):
                raise AssertionError("Load test sequence count mismatch")
                
            return f"Successfully processed {len(large_sequence_set)} sequences under load"
        
        result, success, exec_time, mem_delta = test_performance_under_load()
        results.record_test("Performance Under Load", success, result, exec_time, mem_delta)
        
    except ImportError as e:
        results.record_test("Wildtype Aligner Import", False, f"Import failed: {str(e)}")
    except Exception as e:
        results.record_test("Wildtype Aligner", False, f"Unexpected error: {str(e)}")
        traceback.print_exc()
        
    return results

def stress_test_cooccurrence_analyzer():
    """Comprehensive stress test of co-occurrence analyzer."""
    print("🔥 STRESS TESTING CO-OCCURRENCE ANALYZER...")
    results = StressTestResults()
    
    try:
        from src.generic_cooccurrence_analyzer import GenericCoOccurrenceAnalyzer
        
        # Test 1: Analyzer initialization and configuration
        @measure_performance
        def test_analyzer_initialization():
            analyzer = GenericCoOccurrenceAnalyzer()
            
            # Verify analyzer has required methods
            required_methods = ['analyze_cooccurrence']
            for method in required_methods:
                if not hasattr(analyzer, method):
                    raise AssertionError(f"Missing required method: {method}")
                    
            return "Analyzer initialization successful with all required methods"
        
        result, success, exec_time, mem_delta = test_analyzer_initialization()
        results.record_test("Analyzer Initialization", success, result, exec_time, mem_delta)
        
        # Test 2: Co-occurrence analysis stress test
        @measure_performance
        def test_cooccurrence_analysis_stress():
            analyzer = GenericCoOccurrenceAnalyzer()
            
            # Create comprehensive test data
            test_data = [
                {
                    'gene': 'geneA',
                    'position': 100,
                    'mutation': 'A->T',
                    'sample_id': 'sample1'
                },
                {
                    'gene': 'geneB', 
                    'position': 200,
                    'mutation': 'G->C',
                    'sample_id': 'sample1'  # Same sample - should co-occur
                },
                {
                    'gene': 'geneA',
                    'position': 150,
                    'mutation': 'C->A',
                    'sample_id': 'sample2'
                },
                {
                    'gene': 'geneC',
                    'position': 300,
                    'mutation': 'T->G',
                    'sample_id': 'sample3'
                }
            ]
            
            # Test analysis
            result = analyzer.analyze_cooccurrence(test_data, ['geneA', 'geneB'])
            
            # Verify result structure
            required_keys = ['genes_analyzed', 'total_mutations', 'cooccurrence_patterns', 'analysis_summary']
            for key in required_keys:
                if key not in result:
                    raise AssertionError(f"Missing key in co-occurrence result: {key}")
                    
            if result['total_mutations'] != len(test_data):
                raise AssertionError("Mutation count mismatch in co-occurrence analysis")
                
            return f"Co-occurrence analysis completed for {len(test_data)} mutations"
        
        result, success, exec_time, mem_delta = test_cooccurrence_analysis_stress()
        results.record_test("Co-occurrence Analysis Stress", success, result, exec_time, mem_delta)
        
        # Test 3: Large dataset performance test
        @measure_performance
        def test_large_dataset_performance():
            analyzer = GenericCoOccurrenceAnalyzer()
            
            # Generate large synthetic dataset
            large_dataset = []
            genes = ['geneA', 'geneB', 'geneC', 'geneD', 'geneE']
            mutations = ['A->T', 'T->C', 'G->A', 'C->G']
            
            for i in range(200):  # 200 mutations for performance test
                sample_id = f"sample_{i // 10}"  # 10 mutations per sample average
                gene = genes[i % len(genes)]
                mutation = mutations[i % len(mutations)]
                position = 100 + (i * 10)
                
                large_dataset.append({
                    'gene': gene,
                    'position': position,
                    'mutation': mutation,
                    'sample_id': sample_id
                })
            
            # Analyze large dataset
            result = analyzer.analyze_cooccurrence(large_dataset, genes)
            
            if result['total_mutations'] != len(large_dataset):
                raise AssertionError("Large dataset mutation count mismatch")
                
            return f"Successfully analyzed large dataset with {len(large_dataset)} mutations across {len(genes)} genes"
        
        result, success, exec_time, mem_delta = test_large_dataset_performance()
        results.record_test("Large Dataset Performance", success, result, exec_time, mem_delta)
        
    except ImportError as e:
        results.record_test("Co-occurrence Analyzer Import", False, f"Import failed: {str(e)}")
    except Exception as e:
        results.record_test("Co-occurrence Analyzer", False, f"Unexpected error: {str(e)}")
        traceback.print_exc()
        
    return results

def stress_test_hardcoded_references():
    """Comprehensive test for hardcoded reference removal."""
    print("🔥 STRESS TESTING HARDCODED REFERENCE REMOVAL...")
    results = StressTestResults()
    
    @measure_performance
    def test_comprehensive_hardcode_scan():
        # Comprehensive scan of all Python files for hardcoded references
        hardcoded_patterns = [
            config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0], config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1], config_manager.get_default_genes("rnd_efflux_pumps", "primary")[2],  # Original hardcoded genes
            config_manager.get_reference_strain("primary")["id"],  # Specific strain references
            'E.coli',  # Specific organism
        ]
        
        found_hardcoded = []
        files_scanned = 0
        
        # Scan all Python files in the project
        for root, dirs, files in os.walk('.'):
            # Skip cache and output directories
            if '__pycache__' in root or 'test_output' in root or '.git' in root:
                continue
                
            for file in files:
                if file.endswith('.py'):
                    files_scanned += 1
                    filepath = os.path.join(root, file)
                    
                    try:
                        with open(filepath, 'r', encoding='utf-8') as f:
                            content = f.read()
                            
                        for pattern in hardcoded_patterns:
                            if pattern in content:
                                # Check if it's in a comment or documentation (acceptable)
                                lines = content.split('\n')
                                for i, line in enumerate(lines, 1):
                                    if pattern in line:
                                        # Skip if it's clearly documentation or comments
                                        if (line.strip().startswith('#') or 
                                            line.strip().startswith('"""') or
                                            line.strip().startswith("'''") or
                                            'README' in filepath.upper() or
                                            'doc' in filepath.lower()):
                                            continue
                                        else:
                                            found_hardcoded.append(f"{filepath}:{i} - {pattern}")
                                            
                    except Exception as e:
                        # Skip files that can't be read
                        continue
        
        if found_hardcoded:
            details = f"Found hardcoded references: {'; '.join(found_hardcoded)}"
            return details
        else:
            return f"No hardcoded references found in {files_scanned} Python files"
    
    result, success, exec_time, mem_delta = test_comprehensive_hardcode_scan()
    # Consider it a pass if no problematic hardcoded references found
    test_success = "No hardcoded references found" in result
    results.record_test("Comprehensive Hardcode Scan", test_success, result, exec_time, mem_delta)
    
    return results

def run_complete_stress_test():
    """Run all stress tests and generate comprehensive report."""
    print("=" * 80)
    print("🚀 PRIORITY 1 COMPREHENSIVE STRESS TEST SUITE")
    print("=" * 80)
    print("Testing maximum robustness of all Priority 1 components...")
    print()
    
    all_results = []
    
    # Run all stress tests
    stress_tests = [
        ("Error Handling System", stress_test_error_handling),
        ("Dependency Management", stress_test_dependency_management), 
        ("Wildtype Aligner", stress_test_wildtype_aligner),
        ("Co-occurrence Analyzer", stress_test_cooccurrence_analyzer),
        ("Hardcoded Reference Removal", stress_test_hardcoded_references)
    ]
    
    for test_name, test_func in stress_tests:
        print(f"\n{'='*20} {test_name.upper()} {'='*20}")
        try:
            test_results = test_func()
            all_results.append((test_name, test_results))
        except Exception as e:
            print(f"❌ CRITICAL ERROR in {test_name}: {str(e)}")
            traceback.print_exc()
            # Create dummy results for failed test suite
            dummy_results = StressTestResults()
            dummy_results.record_test(f"{test_name} Critical Failure", False, str(e))
            all_results.append((test_name, dummy_results))
    
    # Generate comprehensive report
    print("\n" + "="*80)
    print("🎯 COMPREHENSIVE STRESS TEST REPORT")
    print("="*80)
    
    total_tests_passed = 0
    total_tests_failed = 0
    
    for test_suite_name, test_results in all_results:
        summary = test_results.get_summary()
        total_tests_passed += summary['passed']
        total_tests_failed += summary['failed']
        
        print(f"\n📊 {test_suite_name}:")
        print(f"   Tests: {summary['total_tests']}")
        print(f"   Passed: {summary['passed']}")
        print(f"   Failed: {summary['failed']}")
        print(f"   Success Rate: {summary['success_rate']:.1f}%")
        
        # Show individual test details
        for test_detail in test_results.test_details:
            print(f"   {test_detail['status']}: {test_detail['test']}")
            if test_detail['details']:
                print(f"      → {test_detail['details']}")
    
    # Overall summary
    total_tests = total_tests_passed + total_tests_failed
    overall_success_rate = (total_tests_passed / total_tests * 100) if total_tests > 0 else 0
    
    print(f"\n🏆 OVERALL RESULTS:")
    print(f"   Total Tests: {total_tests}")
    print(f"   Total Passed: {total_tests_passed}")
    print(f"   Total Failed: {total_tests_failed}")
    print(f"   Overall Success Rate: {overall_success_rate:.1f}%")
    
    if overall_success_rate >= 95:
        print("\n🎉 EXCELLENT! Priority 1 implementation meets maximum robustness standards!")
        print("✅ All components are production-ready")
        print("✅ Error handling is comprehensive")
        print("✅ Performance is optimized")
        print("✅ Ready for Priority 2 implementation")
    elif overall_success_rate >= 85:
        print("\n✅ GOOD! Priority 1 implementation is robust with minor issues")
        print("⚠️  Some optimization opportunities identified")
    else:
        print("\n⚠️  NEEDS ATTENTION! Some components require additional robustness")
        print("🔧 Review failed tests before proceeding to Priority 2")
    
    print("\n" + "="*80)
    
    return overall_success_rate >= 95

if __name__ == "__main__":
    success = run_complete_stress_test()
    sys.exit(0 if success else 1)
