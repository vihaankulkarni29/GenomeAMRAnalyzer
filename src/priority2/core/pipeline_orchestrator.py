"""
Pipeline Orchestrator
--------------------
High-level workflow orchestration for AMR genome analysis pipeline.
Manages sequence processing, analysis, and reporting with robust error handling
and performance optimization.
"""

import os
import json
import time
import logging
from pathlib import Path
from typing import List, Dict, Any, Optional, Callable
from dataclasses import dataclass, asdict, field
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading

from .enhanced_sequence_processor import EnhancedSequenceProcessor, ProcessingStats
from ..config.pipeline_config import PipelineConfig

@dataclass
class PipelineResults:
    """Results from pipeline execution."""
    total_samples: int = 0
    successful_samples: int = 0
    failed_samples: int = 0
    total_processing_time: float = 0.0
    alignment_stats: Dict[str, ProcessingStats] = field(default_factory=dict)
    output_files: Dict[str, List[str]] = field(default_factory=dict)
    errors: List[str] = field(default_factory=list)

    def __post_init__(self):
        # ensure types are proper even if deserialized/overridden
        if not isinstance(self.alignment_stats, dict):
            self.alignment_stats = {}
        if not isinstance(self.output_files, dict):
            self.output_files = {}
        if not isinstance(self.errors, list):
            self.errors = []

class PipelineOrchestrator:
    """
    Orchestrates the complete AMR analysis pipeline with robust error handling
    and performance optimization.
    """
    
    def __init__(self, config: PipelineConfig, logger: Optional[logging.Logger] = None):
        self.config = config
        self.logger = logger or self._setup_logger()
        self.results = PipelineResults()
        self._validate_config()
        
        # Initialize processor
        self.processor = EnhancedSequenceProcessor(
            reference_path=config.database.amr_database_path,
            threads=config.processing.threads,
            logger=self.logger,
            config=config.processing
        )
        
        # Thread safety
        self._lock = threading.Lock()

    def _setup_logger(self) -> logging.Logger:
        logger = logging.getLogger("PipelineOrchestrator")
        if not logger.handlers:
            handler = logging.StreamHandler()
            formatter = logging.Formatter('%(asctime)s %(levelname)s %(name)s: %(message)s')
            handler.setFormatter(formatter)
            logger.addHandler(handler)
        logger.setLevel(logging.INFO)
        return logger

    def _validate_config(self):
        """Validate pipeline configuration using central validator."""
        issues = self.config.validate()
        if issues:
            # non-fatal unless database missing
            missing_db = [i for i in issues if "AMR database" in i]
            if missing_db:
                raise ValueError("; ".join(issues))
            for i in issues:
                self.logger.warning(f"Configuration issue: {i}")
        # ensure output directory exists
        os.makedirs(self.config.output.base_output_dir, exist_ok=True)

    def run_pipeline(self, input_samples: List[str], 
                    progress_callback: Optional[Callable[[str, float], None]] = None) -> PipelineResults:
        """
        Run the complete AMR analysis pipeline.
        
        Args:
            input_samples: List of FASTA/FASTQ file paths
            progress_callback: Optional callback for progress updates
        
        Returns:
            PipelineResults with comprehensive execution statistics
        """
        start_time = time.time()
        self.logger.info(f"Starting pipeline with {len(input_samples)} samples")
        
        # Reset results
        self.results = PipelineResults()
        self.results.total_samples = len(input_samples)
        
        try:
            # Validate inputs
            valid_samples = self._validate_inputs(input_samples)
            
            if not valid_samples:
                raise ValueError("No valid input samples found")
            
            # Process samples
            if self.config.processing.enable_parallel and len(valid_samples) > 1:
                self._process_samples_parallel(valid_samples, progress_callback)
            else:
                self._process_samples_sequential(valid_samples, progress_callback)
            
            # Generate final report
            self._generate_pipeline_report()
            
        except Exception as e:
            self.logger.error(f"Pipeline execution failed: {e}")
            self.results.errors.append(str(e))
        finally:
            self.results.total_processing_time = time.time() - start_time
            self.logger.info(f"Pipeline completed in {self.results.total_processing_time:.2f}s")
        
        return self.results

    def _validate_inputs(self, input_samples: List[str]) -> List[str]:
        """Validate input sample files."""
        valid_samples = []
        for sample in input_samples:
            if os.path.exists(sample):
                # Basic format validation
                if self._is_valid_sequence_file(sample):
                    valid_samples.append(sample)
                else:
                    self.logger.warning(f"Invalid sequence file format: {sample}")
            else:
                self.logger.warning(f"Sample file not found: {sample}")
        
        self.logger.info(f"Validated {len(valid_samples)}/{len(input_samples)} input samples")
        return valid_samples

    def _is_valid_sequence_file(self, filepath: str) -> bool:
        """Check if file is a valid FASTA/FASTQ file."""
        try:
            with open(filepath, 'r') as f:
                first_line = f.readline().strip()
                return first_line.startswith('>') or first_line.startswith('@')
        except Exception:
            return False

    def _process_samples_parallel(self, samples: List[str], 
                                progress_callback: Optional[Callable[[str, float], None]]):
        """Process samples in parallel."""
        self.logger.info(f"Processing {len(samples)} samples in parallel (max workers: {self.config.processing.max_parallel_samples})")
        
        with ThreadPoolExecutor(max_workers=self.config.processing.max_parallel_samples) as executor:
            # Submit all jobs
            future_to_sample = {
                executor.submit(self._process_single_sample, sample, i): sample 
                for i, sample in enumerate(samples)
            }
            
            completed = 0
            for future in as_completed(future_to_sample):
                sample = future_to_sample[future]
                completed += 1
                
                try:
                    result = future.result()
                    if result:
                        with self._lock:
                            self.results.successful_samples += 1
                            sample_name = Path(sample).stem
                            self.results.alignment_stats[sample_name] = result
                    else:
                        with self._lock:
                            self.results.failed_samples += 1
                except Exception as e:
                    self.logger.error(f"Sample processing failed for {sample}: {e}")
                    with self._lock:
                        self.results.failed_samples += 1
                        self.results.errors.append(f"{sample}: {str(e)}")
                
                # Progress callback
                if progress_callback:
                    progress = (completed / len(samples)) * 100
                    progress_callback(f"Processed {completed}/{len(samples)} samples", progress)

    def _process_samples_sequential(self, samples: List[str], 
                                  progress_callback: Optional[Callable[[str, float], None]]):
        """Process samples sequentially."""
        self.logger.info(f"Processing {len(samples)} samples sequentially")
        
        for i, sample in enumerate(samples):
            try:
                result = self._process_single_sample(sample, i)
                if result:
                    self.results.successful_samples += 1
                    sample_name = Path(sample).stem
                    self.results.alignment_stats[sample_name] = result
                else:
                    self.results.failed_samples += 1
            except Exception as e:
                self.logger.error(f"Sample processing failed for {sample}: {e}")
                self.results.failed_samples += 1
                self.results.errors.append(f"{sample}: {str(e)}")
            
            # Progress callback
            if progress_callback:
                progress = ((i + 1) / len(samples)) * 100
                progress_callback(f"Processed {i + 1}/{len(samples)} samples", progress)

    def _process_single_sample(self, sample_path: str, sample_index: int) -> Optional[ProcessingStats]:
        """Process a single sample through the pipeline."""
        sample_name = Path(sample_path).stem
        sample_output_dir = os.path.join(self.config.output.base_output_dir, sample_name)
        os.makedirs(sample_output_dir, exist_ok=True)
        
        self.logger.info(f"Processing sample: {sample_name}")
        
        try:
            # Reset processor stats for this sample
            self.processor.reset_stats()
            
            # Run alignment
            alignment_files = self.processor.align_sequences(
                query_fastas=[sample_path],
                output_dir=sample_output_dir,
                validate_inputs=True
            )
            
            if alignment_files:
                # Collect quality metrics
                quality_metrics = self.processor.collect_quality_metrics(alignment_files)
                
                # Save sample results
                if self.config.output.save_intermediate:
                    self._save_sample_results(sample_name, sample_output_dir, quality_metrics)
                
                # Store output files
                with self._lock:
                    self.results.output_files[sample_name] = alignment_files
                
                return self.processor.get_processing_stats()
            else:
                self.logger.warning(f"No alignment files generated for {sample_name}")
                return None
                
        except Exception as e:
            self.logger.error(f"Failed to process sample {sample_name}: {e}")
            raise

    def _save_sample_results(self, sample_name: str, output_dir: str, quality_metrics: Dict[str, Any]):
        """Save intermediate results for a sample."""
        results_file = os.path.join(output_dir, f"{sample_name}_results.json")
        
        sample_results = {
            "sample_name": sample_name,
            "processing_stats": asdict(self.processor.get_processing_stats()),
            "quality_metrics": quality_metrics,
            "timestamp": time.time()
        }
        
        with open(results_file, 'w') as f:
            json.dump(sample_results, f, indent=2)

    def _generate_pipeline_report(self):
        """Generate comprehensive pipeline report."""
        report_path = os.path.join(self.config.output.base_output_dir, "pipeline_report.json")
        
        report_data = {
            "pipeline_config": self.config.to_dict(),
            "execution_results": asdict(self.results),
            "summary": {
                "success_rate": (self.results.successful_samples / self.results.total_samples * 100) 
                               if self.results.total_samples > 0 else 0,
                "average_processing_time": self.results.total_processing_time / self.results.total_samples
                                         if self.results.total_samples > 0 else 0
            },
            "timestamp": time.time()
        }
        
        with open(report_path, 'w') as f:
            json.dump(report_data, f, indent=2)
        
        self.logger.info(f"Pipeline report saved to: {report_path}")

    def get_results_summary(self) -> Dict[str, Any]:
        """Get a summary of pipeline results."""
        return {
            "total_samples": self.results.total_samples,
            "successful": self.results.successful_samples,
            "failed": self.results.failed_samples,
            "success_rate": (self.results.successful_samples / self.results.total_samples * 100) 
                           if self.results.total_samples > 0 else 0,
            "total_time": f"{self.results.total_processing_time:.2f}s",
            "errors": len(self.results.errors)
        }
