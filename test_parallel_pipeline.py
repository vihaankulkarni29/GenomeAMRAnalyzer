#!/usr/bin/env python3
"""
Test script to validate parallel processing improvements in the GenomeAMRAnalyzer pipeline.
Tests both sequential and parallel modes to measure performance gains.
"""

import os
import sys
import time
import subprocess
import tempfile
import shutil
from pathlib import Path

# Add src directory to path to import our modules
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

def create_test_config(test_dir: Path, accessions: list, threads: int = None) -> Path:
    """Create a test configuration file with small batch of genomes"""
    
    # Create accession list file
    accession_file = test_dir / "test_accessions.txt"
    with open(accession_file, 'w') as f:
        for acc in accessions:
            f.write(f"{acc}\n")
    
    # Create minimal config
    config_content = f"""
database:
  path: "{test_dir / 'test.db'}"

logging:
  log_dir: "{test_dir / 'logs'}"

ncbi:
  email: "test@example.com"
  api_key: ""
  mock_mode: true  # Use mock mode to avoid real NCBI API calls during testing

input:
  accessions_file: "{accession_file}"

genome_output_dir: "{test_dir / 'genome_data'}"
card_output_dir: "{test_dir / 'card_results'}"
extractor_output_dir: "{test_dir / 'extracted_proteins'}"
alignment_output_dir: "{test_dir / 'alignments'}"
subscan_output_dir: "{test_dir / 'subscan_results'}"
cooccurrence_output_dir: "{test_dir / 'cooccurrence_results'}"

reference_dir: "reference_proteins"

subscan:
  min_confidence: "MEDIUM"

cooccurrence:
  min_count: 3
  significance_level: 0.01

run_id: "parallel_test_run"
"""
    
    config_file = test_dir / "test_config.yaml"
    with open(config_file, 'w') as f:
        f.write(config_content)
    
    return config_file

def run_pipeline_test(config_file: Path, threads: int = None, test_name: str = "test") -> dict:
    """Run the pipeline and measure performance"""
    
    cmd = [
        sys.executable, 
        "src/production_pipeline_orchestrator.py",
        "--config", str(config_file),
        "--stage", "card_scan"  # Just test the CARD scanning stage for speed
    ]
    
    if threads:
        cmd.extend(["--threads", str(threads)])
    
    print(f"\n{'='*60}")
    print(f"Running {test_name} with command: {' '.join(cmd)}")
    print(f"{'='*60}")
    
    start_time = time.time()
    
    try:
        # Run the pipeline
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300,  # 5 minute timeout
            cwd=os.path.dirname(__file__)
        )
        
        end_time = time.time()
        processing_time = end_time - start_time
        
        return {
            'success': result.returncode == 0,
            'processing_time': processing_time,
            'stdout': result.stdout,
            'stderr': result.stderr,
            'return_code': result.returncode,
            'threads': threads or 1
        }
        
    except subprocess.TimeoutExpired:
        end_time = time.time()
        return {
            'success': False,
            'processing_time': end_time - start_time,
            'stdout': "",
            'stderr': "Test timed out after 5 minutes",
            'return_code': -1,
            'threads': threads or 1
        }

def analyze_performance(results: list) -> None:
    """Analyze and report performance results"""
    
    print(f"\n{'='*80}")
    print("PARALLEL PIPELINE PERFORMANCE ANALYSIS")
    print(f"{'='*80}")
    
    successful_tests = [r for r in results if r['success']]
    failed_tests = [r for r in results if not r['success']]
    
    if not successful_tests:
        print("❌ No successful tests to analyze!")
        for i, test in enumerate(failed_tests):
            print(f"\nTest {i+1} (Failed - {test['threads']} threads):")
            print(f"  Return code: {test['return_code']}")
            print(f"  Error: {test['stderr'][:200]}...")
        return
    
    # Find sequential baseline (1 thread)
    sequential_test = next((r for r in successful_tests if r['threads'] == 1), None)
    
    print("Test Results Summary:")
    print(f"  Total tests: {len(results)}")
    print(f"  Successful: {len(successful_tests)}")
    print(f"  Failed: {len(failed_tests)}")
    
    print("\nPerformance Results:")
    for test in successful_tests:
        threads = test['threads']
        time_taken = test['processing_time']
        
        if sequential_test and threads > 1:
            speedup = sequential_test['processing_time'] / time_taken
            efficiency = speedup / threads * 100
            print(f"  {threads} threads: {time_taken:.2f}s (speedup: {speedup:.2f}x, efficiency: {efficiency:.1f}%)")
        else:
            print(f"  {threads} thread: {time_taken:.2f}s (baseline)")
    
    if sequential_test and len(successful_tests) > 1:
        parallel_tests = [r for r in successful_tests if r['threads'] > 1]
        if parallel_tests:
            best_test = min(parallel_tests, key=lambda x: x['processing_time'])
            max_speedup = sequential_test['processing_time'] / best_test['processing_time']
            print(f"\nBest Performance:")
            print(f"  Configuration: {best_test['threads']} threads")
            print(f"  Maximum speedup: {max_speedup:.2f}x")
            print(f"  Time reduction: {(1 - best_test['processing_time']/sequential_test['processing_time'])*100:.1f}%")
    
    # Show any failures
    if failed_tests:
        print(f"\n❌ Failed Tests ({len(failed_tests)}):")
        for test in failed_tests:
            print(f"  {test['threads']} threads: {test['stderr'][:100]}...")

def main():
    """Main test execution"""
    
    print("🧬 GenomeAMRAnalyzer Parallel Pipeline Performance Test")
    print("=" * 60)
    
    # Create temporary test directory
    test_dir = Path(tempfile.mkdtemp(prefix="genome_amr_test_"))
    print(f"Test directory: {test_dir}")
    
    try:
        # Use first few accessions for quick testing
        test_accessions = [
            "KWV17775.1",
            "KWV17821.1", 
            "KWV19057.1"
        ]
        
        print(f"Testing with {len(test_accessions)} genomes: {', '.join(test_accessions)}")
        
        # Create test configuration
        config_file = create_test_config(test_dir, test_accessions)
        
        # Test different thread configurations
        thread_configs = [1, 2, 4]  # Sequential, then parallel
        
        results = []
        
        for threads in thread_configs:
            test_name = f"Sequential Test" if threads == 1 else f"Parallel Test ({threads} threads)"
            
            result = run_pipeline_test(config_file, threads, test_name)
            results.append(result)
            
            if result['success']:
                print(f"✅ {test_name}: {result['processing_time']:.2f}s")
            else:
                print(f"❌ {test_name}: FAILED")
                print(f"   Error: {result['stderr'][:200]}...")
        
        # Analyze and report results
        analyze_performance(results)
        
    except Exception as e:
        print(f"❌ Test execution failed: {e}")
        import traceback
        traceback.print_exc()
    
    finally:
        # Cleanup test directory
        try:
            shutil.rmtree(test_dir)
            print(f"\n🧹 Cleaned up test directory: {test_dir}")
        except Exception as e:
            print(f"⚠️  Could not clean up test directory {test_dir}: {e}")

if __name__ == "__main__":
    main()