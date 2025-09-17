"""
Priority 2 Example Workflow: Basic AMR Analysis
----------------------------------------------
Demonstrates basic usage of the Priority 2 pipeline for AMR genome analysis.
"""

import os
import sys
import logging
from pathlib import Path

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent / "src"))

from priority2.config.pipeline_config import PipelineConfig, DatabaseConfig, ConfigManager
from priority2.core.pipeline_orchestrator import PipelineOrchestrator

def setup_logging():
    """Set up logging for the example."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )

def create_example_config(amr_db_path: str, output_dir: str) -> PipelineConfig:
    """Create example configuration for AMR analysis."""
    
    config_manager = ConfigManager()
    
    # Create a sensitive detection configuration
    config = config_manager.create_sensitive_config(amr_db_path, output_dir)
    
    # Customize for this example
    config.processing.threads = 4
    config.processing.min_identity = 80.0
    config.processing.min_coverage = 50.0
    config.output.generate_reports = True
    
    return config

def basic_amr_analysis_example():
    """
    Example: Basic AMR analysis workflow
    
    This example demonstrates:
    1. Configuration setup
    2. Pipeline initialization
    3. Sample processing
    4. Results interpretation
    """
    
    print("=" * 60)
    print("Priority 2 Example: Basic AMR Analysis")
    print("=" * 60)
    
    setup_logging()
    
    # Example paths (adjust for your system)
    amr_database = "/path/to/amr_database.fasta"  # Replace with actual path
    input_samples = [
        "/path/to/sample1.fasta",  # Replace with actual sample paths
        "/path/to/sample2.fasta",
        "/path/to/sample3.fasta"
    ]
    output_directory = "./example_output/basic_analysis"
    
    # Create output directory
    os.makedirs(output_directory, exist_ok=True)
    
    try:
        # Step 1: Create configuration
        print("Step 1: Creating pipeline configuration...")
        config = create_example_config(amr_database, output_directory)
        
        # Validate configuration
        config_manager = ConfigManager()
        if config_manager.validate_config(config):
            print("✓ Configuration validated successfully")
        else:
            print("⚠ Configuration validation warnings (check logs)")
        
        # Step 2: Initialize pipeline
        print("Step 2: Initializing pipeline...")
        orchestrator = PipelineOrchestrator(config)
        
        # Step 3: Define progress callback
        def progress_callback(message: str, progress: float):
            print(f"Progress: {progress:.1f}% - {message}")
        
        # Step 4: Run analysis
        print("Step 3: Running AMR analysis...")
        results = orchestrator.run_pipeline(
            input_samples=input_samples,
            progress_callback=progress_callback
        )
        
        # Step 5: Display results
        print("Step 4: Analysis Results")
        print("-" * 40)
        summary = orchestrator.get_results_summary()
        
        print(f"Total samples processed: {summary['total_samples']}")
        print(f"Successful analyses: {summary['successful']}")
        print(f"Failed analyses: {summary['failed']}")
        print(f"Success rate: {summary['success_rate']:.1f}%")
        print(f"Total processing time: {summary['total_time']}")
        
        if summary['errors'] > 0:
            print(f"Errors encountered: {summary['errors']}")
        
        print(f"\nDetailed results saved to: {output_directory}")
        
        return True
        
    except Exception as e:
        print(f"Error in basic analysis: {e}")
        return False

def high_throughput_example():
    """
    Example: High-throughput processing
    
    Demonstrates processing of large datasets with optimization.
    """
    
    print("=" * 60)
    print("Priority 2 Example: High-Throughput Processing")
    print("=" * 60)
    
    # Example for processing many samples efficiently
    amr_database = "/path/to/amr_database.fasta"
    output_directory = "./example_output/high_throughput"
    
    # Large sample list (100+ samples)
    input_samples = [f"/path/to/sample_{i}.fasta" for i in range(1, 101)]
    
    os.makedirs(output_directory, exist_ok=True)
    
    try:
        # Create high-throughput configuration
        config_manager = ConfigManager()
        config = config_manager.create_high_throughput_config(amr_database, output_directory)
        
        # Optimize for large datasets
        config.processing.max_parallel_samples = 8
        config.processing.chunk_size = 5000
        config.output.save_intermediate = False  # Save space
        config.output.compress_outputs = True
        
        print(f"Configuration: {config.processing.max_parallel_samples} parallel workers")
        print(f"Chunk size: {config.processing.chunk_size}")
        
        # Initialize and run
        orchestrator = PipelineOrchestrator(config)
        
        def high_throughput_progress(message: str, progress: float):
            if progress % 10 == 0:  # Report every 10%
                print(f"High-throughput progress: {progress:.0f}% - {message}")
        
        results = orchestrator.run_pipeline(
            input_samples=input_samples,
            progress_callback=high_throughput_progress
        )
        
        summary = orchestrator.get_results_summary()
        print(f"Processed {summary['total_samples']} samples")
        print(f"Throughput: {summary['total_samples'] / float(summary['total_time'].replace('s', '')):.2f} samples/second")
        
        return True
        
    except Exception as e:
        print(f"Error in high-throughput analysis: {e}")
        return False

def custom_analysis_example():
    """
    Example: Custom analysis with specific parameters
    
    Shows how to customize the pipeline for specific research questions.
    """
    
    print("=" * 60)
    print("Priority 2 Example: Custom Analysis Parameters")
    print("=" * 60)
    
    amr_database = "/path/to/amr_database.fasta"
    output_directory = "./example_output/custom_analysis"
    input_samples = ["/path/to/ecoli_samples.fasta"]
    
    os.makedirs(output_directory, exist_ok=True)
    
    try:
        # Create custom configuration
        from priority2.config.pipeline_config import ProcessingConfig, OutputConfig, AlignmentPreset
        
        # Custom processing parameters for E. coli RND efflux pump analysis
        custom_processing = ProcessingConfig(
            threads=6,
            preset=AlignmentPreset.MAP_ONT,
            min_identity=95.0,  # High stringency for exact matches
            min_coverage=80.0,  # Good coverage required
            quality_threshold=30,  # High quality alignments only
            enable_parallel=True,
            max_parallel_samples=2
        )
        
        # Custom output configuration
        custom_output = OutputConfig(
            base_output_dir=output_directory,
            generate_reports=True,
            report_format="html",
            save_intermediate=True
        )
        
        config = PipelineConfig(
            database=DatabaseConfig(amr_database_path=amr_database),
            processing=custom_processing,
            output=custom_output,
            config_name="custom_efflux_analysis"
        )
        
        print("Custom parameters:")
        print(f"- Identity threshold: {config.processing.min_identity}%")
        print(f"- Coverage threshold: {config.processing.min_coverage}%")
        print(f"- Quality threshold: {config.processing.quality_threshold}")
        
        # Run analysis
        orchestrator = PipelineOrchestrator(config)
        results = orchestrator.run_pipeline(input_samples)
        
        # Analyze results for efflux pumps specifically
        print("\nAnalyzing efflux pump patterns...")
        
        # This would typically involve post-processing the results
        # to identify co-occurrence patterns in efflux pump genes
        
        return True
        
    except Exception as e:
        print(f"Error in custom analysis: {e}")
        return False

def large_dataset_example():
    """
    Example: Processing very large datasets with streaming
    
    Demonstrates the large dataset processor for millions of sequences.
    """
    
    print("=" * 60)
    print("Priority 2 Example: Large Dataset Processing")
    print("=" * 60)
    
    from priority2.pipelines.large_dataset_processor import LargeDatasetProcessor
    
    amr_database = "/path/to/amr_database.fasta"
    large_dataset_file = "/path/to/million_sequences.fasta"  # Very large file
    output_directory = "./example_output/large_dataset"
    
    os.makedirs(output_directory, exist_ok=True)
    
    try:
        # Initialize large dataset processor
        processor = LargeDatasetProcessor(
            reference_path=amr_database,
            output_dir=output_directory,
            chunk_size=10000,  # Process 10K sequences per chunk
            max_workers=4,
            memory_limit_gb=16
        )
        
        def streaming_progress(stats):
            if stats.chunks_processed % 10 == 0:  # Every 10 chunks
                progress = (stats.processed_sequences / stats.total_sequences) * 100
                print(f"Streaming progress: {progress:.1f}% "
                      f"({stats.processed_sequences:,}/{stats.total_sequences:,} sequences)")
                print(f"  Throughput: {stats.current_throughput:.1f} seq/s")
                
                if stats.estimated_completion:
                    print(f"  ETA: {stats.estimated_completion/60:.1f} minutes")
        
        print(f"Processing large dataset: {large_dataset_file}")
        print(f"Chunk size: {processor.chunk_size:,} sequences")
        print(f"Max workers: {processor.max_workers}")
        
        # Process the large dataset
        results = processor.process_large_dataset(
            large_dataset_file,
            progress_callback=streaming_progress
        )
        
        # Display final statistics
        final_stats = processor.get_processing_stats()
        report = processor.create_processing_report()
        
        print("\nFinal Results:")
        print(f"Total sequences: {final_stats.total_sequences:,}")
        print(f"Processing rate: {report['dataset_info']['processing_rate']}")
        print(f"Peak memory usage: {report['performance_metrics']['peak_memory_mb']}")
        print(f"Total time: {report['performance_metrics']['total_time']}")
        
        return True
        
    except Exception as e:
        print(f"Error in large dataset processing: {e}")
        return False

def main():
    """Run all example workflows."""
    
    print("Priority 2 AMR Analysis Pipeline - Example Workflows")
    print("=" * 60)
    
    examples = [
        ("Basic AMR Analysis", basic_amr_analysis_example),
        ("High-Throughput Processing", high_throughput_example),
        ("Custom Analysis Parameters", custom_analysis_example),
        ("Large Dataset Processing", large_dataset_example)
    ]
    
    results = {}
    
    for name, example_func in examples:
        print(f"\nRunning: {name}")
        try:
            success = example_func()
            results[name] = "SUCCESS" if success else "FAILED"
        except Exception as e:
            print(f"Example failed with error: {e}")
            results[name] = "ERROR"
        
        print(f"Result: {results[name]}")
        print("-" * 60)
    
    # Summary
    print("\nExample Workflow Summary:")
    for name, result in results.items():
        status_icon = "✓" if result == "SUCCESS" else "✗"
        print(f"{status_icon} {name}: {result}")

if __name__ == "__main__":
    main()