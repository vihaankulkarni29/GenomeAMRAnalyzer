"""
Large Dataset Processor
----------------------
High-performance processor for handling millions of sequences with streaming,
chunked processing, and memory optimization for massive genomic datasets.
"""

import os
import gc
import time
import logging
import threading
from pathlib import Path
from typing import List, Dict, Any, Optional, Iterator, Callable, Generator
from dataclasses import dataclass, field
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
import multiprocessing as mp
from queue import Queue, Empty
import tempfile
import shutil

# Expose AlignmentToolManager at module level for test patching
try:
    from ..integrations.external_alignment_tools import AlignmentToolManager as AlignmentToolManager  # type: ignore
except Exception:  # pragma: no cover - will be patched in tests
    AlignmentToolManager = None  # type: ignore

@dataclass
class ProcessingChunk:
    """Represents a chunk of sequences for processing."""
    chunk_id: int
    sequences: List[tuple]  # (name, sequence) pairs
    chunk_size: int
    temp_file_path: Optional[str] = None

@dataclass
class StreamingStats:
    """Statistics for streaming processing."""
    total_sequences: int = 0
    processed_sequences: int = 0
    chunks_processed: int = 0
    processing_time: float = 0.0
    peak_memory_mb: float = 0.0
    current_throughput: float = 0.0  # sequences per second
    estimated_completion: Optional[float] = None

class SequenceStreamReader:
    """Memory-efficient sequence file reader with streaming support."""
    
    def __init__(self, file_path: str, chunk_size: int = 1000):
        self.file_path = file_path
        self.chunk_size = chunk_size
        self.logger = logging.getLogger("SequenceStreamReader")
        
    def read_chunks(self) -> Generator[ProcessingChunk, None, None]:
        """Read sequences in chunks to minimize memory usage."""
        chunk_id = 0
        current_chunk = []
        
        try:
            with open(self.file_path, 'r') as f:
                current_name = None
                current_seq = []
                
                for line_num, line in enumerate(f):
                    line = line.strip()
                    
                    if line.startswith('>') or line.startswith('@'):
                        # Save previous sequence if exists
                        if current_name and current_seq:
                            sequence = ''.join(current_seq)
                            current_chunk.append((current_name, sequence))
                            
                            # Yield chunk if full
                            if len(current_chunk) >= self.chunk_size:
                                yield ProcessingChunk(
                                    chunk_id=chunk_id,
                                    sequences=current_chunk,
                                    chunk_size=len(current_chunk)
                                )
                                chunk_id += 1
                                current_chunk = []
                        
                        # Start new sequence
                        current_name = line[1:].split()[0]
                        current_seq = []
                        
                    elif line and not line.startswith('+'):  # Skip quality lines in FASTQ
                        current_seq.append(line)
                
                # Handle last sequence
                if current_name and current_seq:
                    sequence = ''.join(current_seq)
                    current_chunk.append((current_name, sequence))
                
                # Yield final chunk
                if current_chunk:
                    yield ProcessingChunk(
                        chunk_id=chunk_id,
                        sequences=current_chunk,
                        chunk_size=len(current_chunk)
                    )
                    
        except Exception as e:
            self.logger.error(f"Error reading sequence file {self.file_path}: {e}")
            raise

    def count_sequences(self) -> int:
        """Count total sequences in file efficiently."""
        count = 0
        try:
            with open(self.file_path, 'r') as f:
                for line in f:
                    if line.startswith('>') or line.startswith('@'):
                        count += 1
        except Exception as e:
            self.logger.error(f"Error counting sequences: {e}")
            return 0
        return count

class ChunkProcessor:
    """Processes individual chunks of sequences."""
    
    def __init__(self, reference_path: str, output_dir: str, processing_config: Dict[str, Any]):
        self.reference_path = reference_path
        self.output_dir = output_dir
        self.config = processing_config
        self.logger = logging.getLogger("ChunkProcessor")
        
    def process_chunk(self, chunk: ProcessingChunk) -> Dict[str, Any]:
        """Process a single chunk and return results."""
        start_time = time.time()
        
        # Create temporary files for this chunk
        temp_dir = os.path.join(self.output_dir, f"temp_chunk_{chunk.chunk_id}")
        os.makedirs(temp_dir, exist_ok=True)
        
        try:
            # Write chunk to temporary FASTA file
            chunk_fasta = os.path.join(temp_dir, f"chunk_{chunk.chunk_id}.fasta")
            self._write_chunk_to_file(chunk, chunk_fasta)
            
            # Process with alignment tool
            output_paf = os.path.join(temp_dir, f"chunk_{chunk.chunk_id}.paf")
            success = self._align_chunk(chunk_fasta, output_paf)
            
            if success:
                # Parse results
                results = self._parse_chunk_results(output_paf)
                
                # Clean up temporary files
                self._cleanup_temp_files(temp_dir)
                
                return {
                    'chunk_id': chunk.chunk_id,
                    'processed_sequences': chunk.chunk_size,
                    'processing_time': time.time() - start_time,
                    'results': results,
                    'success': True
                }
            else:
                return {
                    'chunk_id': chunk.chunk_id,
                    'processed_sequences': 0,
                    'processing_time': time.time() - start_time,
                    'results': [],
                    'success': False,
                    'error': 'Alignment failed'
                }
                
        except Exception as e:
            self.logger.error(f"Error processing chunk {chunk.chunk_id}: {e}")
            return {
                'chunk_id': chunk.chunk_id,
                'processed_sequences': 0,
                'processing_time': time.time() - start_time,
                'results': [],
                'success': False,
                'error': str(e)
            }
    
    def _write_chunk_to_file(self, chunk: ProcessingChunk, output_path: str):
        """Write chunk sequences to FASTA file."""
        with open(output_path, 'w') as f:
            for name, sequence in chunk.sequences:
                f.write(f">{name}\n{sequence}\n")
    
    def _align_chunk(self, input_fasta: str, output_paf: str) -> bool:
        """Align chunk using configured alignment tool."""
        try:
            # Use module-level AlignmentToolManager (allows test patching)
            if AlignmentToolManager is None:
                # Fallback import if not available at module import time
                from ..integrations.external_alignment_tools import AlignmentToolManager as _AlignmentToolManager  # type: ignore
                tool_manager = _AlignmentToolManager(self.reference_path)
            else:
                tool_manager = AlignmentToolManager(self.reference_path)
            aligner = tool_manager.get_best_aligner(
                dataset_size=1000,  # Chunk size
                sequence_type=self.config.get('sequence_type', 'long')
            )
            
            return aligner.align_file(input_fasta, output_paf)
            
        except Exception as e:
            self.logger.error(f"Chunk alignment failed: {e}")
            return False
    
    def _parse_chunk_results(self, paf_file: str) -> List[Dict[str, Any]]:
        """Parse alignment results from PAF file."""
        results = []
        try:
            with open(paf_file, 'r') as f:
                for line in f:
                    if line.strip():
                        fields = line.strip().split('\t')
                        if len(fields) >= 12:
                            result = {
                                'query_name': fields[0],
                                'target_name': fields[5],
                                'identity': (int(fields[9]) / int(fields[10]) * 100) if int(fields[10]) > 0 else 0,
                                'coverage': ((int(fields[3]) - int(fields[2])) / int(fields[1]) * 100) if int(fields[1]) > 0 else 0,
                                'mapq': int(fields[11])
                            }
                            results.append(result)
        except Exception as e:
            self.logger.error(f"Error parsing PAF file: {e}")
        
        return results
    
    def _cleanup_temp_files(self, temp_dir: str):
        """Clean up temporary files."""
        try:
            shutil.rmtree(temp_dir)
        except Exception as e:
            self.logger.warning(f"Failed to clean up temp directory {temp_dir}: {e}")

class ResultAggregator:
    """Aggregates results from processed chunks."""
    
    def __init__(self, output_dir: str):
        self.output_dir = output_dir
        self.logger = logging.getLogger("ResultAggregator")
        self.results_queue = Queue()
        self.aggregated_results = []
        self._lock = threading.Lock()
        
    def add_chunk_result(self, chunk_result: Dict[str, Any]):
        """Add result from a processed chunk."""
        with self._lock:
            self.aggregated_results.append(chunk_result)
            
    def finalize_results(self) -> Dict[str, Any]:
        """Finalize and save aggregated results."""
        total_sequences = sum(r.get('processed_sequences', 0) for r in self.aggregated_results)
        successful_chunks = len([r for r in self.aggregated_results if r.get('success', False)])
        total_processing_time = sum(r.get('processing_time', 0) for r in self.aggregated_results)
        
        # Combine all alignment results
        all_alignments = []
        for result in self.aggregated_results:
            if result.get('success', False):
                all_alignments.extend(result.get('results', []))
        
        # Save results
        summary = {
            'total_sequences_processed': total_sequences,
            'successful_chunks': successful_chunks,
            'total_chunks': len(self.aggregated_results),
            'total_processing_time': total_processing_time,
            'total_alignments': len(all_alignments),
            'average_throughput': total_sequences / total_processing_time if total_processing_time > 0 else 0
        }
        
        self._save_summary(summary)
        self._save_detailed_results(all_alignments)
        
        return summary
    
    def _save_summary(self, summary: Dict[str, Any]):
        """Save processing summary."""
        import json
        summary_path = os.path.join(self.output_dir, "processing_summary.json")
        with open(summary_path, 'w') as f:
            json.dump(summary, f, indent=2)
        self.logger.info(f"Processing summary saved to {summary_path}")
    
    def _save_detailed_results(self, alignments: List[Dict[str, Any]]):
        """Save detailed alignment results."""
        try:
            import pandas as pd  # type: ignore
        except ImportError:
            self.logger.warning("pandas not available; skipping detailed CSV results")
            return
        if alignments:
            df = pd.DataFrame(alignments)
            results_path = os.path.join(self.output_dir, "alignment_results.csv")
            df.to_csv(results_path, index=False)
            self.logger.info(f"Detailed results saved to {results_path}")

class LargeDatasetProcessor:
    """
    High-performance processor for massive genomic datasets with streaming,
    chunked processing, and memory optimization.
    """
    
    def __init__(self, reference_path: str, output_dir: str, 
                 chunk_size: int = 5000, max_workers: Optional[int] = None,
                 memory_limit_gb: int = 8):
        self.reference_path = reference_path
        self.output_dir = output_dir
        self.chunk_size = chunk_size
        self.max_workers = (max_workers or min((mp.cpu_count() or 1), 4))
        self.memory_limit_gb = memory_limit_gb
        self.logger = self._setup_logger()
        self.stats = StreamingStats()
        
        os.makedirs(self.output_dir, exist_ok=True)
        
    def _setup_logger(self) -> logging.Logger:
        logger = logging.getLogger("LargeDatasetProcessor")
        if not logger.handlers:
            handler = logging.StreamHandler()
            formatter = logging.Formatter('%(asctime)s %(levelname)s %(name)s: %(message)s')
            handler.setFormatter(formatter)
            logger.addHandler(handler)
        logger.setLevel(logging.INFO)
        return logger
    
    def process_large_dataset(self, input_file: str, 
                            progress_callback: Optional[Callable[[StreamingStats], None]] = None) -> Dict[str, Any]:
        """
        Process a large dataset with streaming and chunked processing.
        """
        start_time = time.time()
        self.logger.info(f"Starting large dataset processing: {input_file}")
        
        try:
            # Initialize components
            reader = SequenceStreamReader(input_file, self.chunk_size)
            aggregator = ResultAggregator(self.output_dir)
            
            # Count total sequences for progress tracking
            self.stats.total_sequences = reader.count_sequences()
            self.logger.info(f"Total sequences to process: {self.stats.total_sequences}")
            
            # Process chunks
            if self.max_workers == 1:
                self._process_chunks_sequential(reader, aggregator, progress_callback)
            else:
                self._process_chunks_parallel(reader, aggregator, progress_callback)
            
            # Finalize results
            final_results = aggregator.finalize_results()
            
            self.stats.processing_time = time.time() - start_time
            self.stats.current_throughput = self.stats.processed_sequences / self.stats.processing_time if self.stats.processing_time > 0 else 0
            
            self.logger.info(f"Large dataset processing completed in {self.stats.processing_time:.2f}s")
            self.logger.info(f"Throughput: {self.stats.current_throughput:.2f} sequences/second")
            
            return final_results
            
        except Exception as e:
            self.logger.error(f"Large dataset processing failed: {e}")
            raise
    
    def _process_chunks_sequential(self, reader: SequenceStreamReader, 
                                 aggregator: ResultAggregator,
                                 progress_callback: Optional[Callable[[StreamingStats], None]]):
        """Process chunks sequentially."""
        processor = ChunkProcessor(
            self.reference_path, 
            self.output_dir, 
            {'sequence_type': 'long'}
        )
        
        for chunk in reader.read_chunks():
            self.logger.info(f"Processing chunk {chunk.chunk_id} ({chunk.chunk_size} sequences)")
            
            result = processor.process_chunk(chunk)
            aggregator.add_chunk_result(result)
            
            # Update statistics
            self.stats.processed_sequences += result.get('processed_sequences', 0)
            self.stats.chunks_processed += 1
            
            # Progress callback
            if progress_callback:
                progress_callback(self.stats)
            
            # Memory management
            gc.collect()
    
    def _process_chunks_parallel(self, reader: SequenceStreamReader, 
                               aggregator: ResultAggregator,
                               progress_callback: Optional[Callable[[StreamingStats], None]]):
        """Process chunks in parallel."""
        
        with ProcessPoolExecutor(max_workers=self.max_workers) as executor:
            # Submit initial chunks
            futures = {}
            chunks_submitted = 0
            
            chunk_iterator = reader.read_chunks()
            
            # Submit initial batch
            for _ in range(self.max_workers * 2):  # Buffer some chunks
                try:
                    chunk = next(chunk_iterator)
                    future = executor.submit(self._process_chunk_wrapper, chunk)
                    futures[future] = chunk.chunk_id
                    chunks_submitted += 1
                except StopIteration:
                    break
            
            # Process results as they complete
            while futures:
                try:
                    for future in as_completed(futures, timeout=1.0):
                        chunk_id = futures.pop(future)
                        
                        try:
                            result = future.result()
                            aggregator.add_chunk_result(result)
                            
                            # Update statistics
                            self.stats.processed_sequences += result.get('processed_sequences', 0)
                            self.stats.chunks_processed += 1
                            
                            self.logger.info(f"Completed chunk {chunk_id} "
                                           f"({self.stats.processed_sequences}/{self.stats.total_sequences})")
                            
                            # Progress callback
                            if progress_callback:
                                progress_callback(self.stats)
                            
                            # Submit next chunk if available
                            try:
                                next_chunk = next(chunk_iterator)
                                next_future = executor.submit(self._process_chunk_wrapper, next_chunk)
                                futures[next_future] = next_chunk.chunk_id
                                chunks_submitted += 1
                            except StopIteration:
                                pass
                                
                        except Exception as e:
                            self.logger.error(f"Chunk {chunk_id} processing failed: {e}")
                            
                except TimeoutError:
                    # Check memory usage and potentially pause if needed
                    self._check_memory_usage()
                    continue
    
    def _process_chunk_wrapper(self, chunk: ProcessingChunk) -> Dict[str, Any]:
        """Wrapper for chunk processing in separate process."""
        processor = ChunkProcessor(
            self.reference_path, 
            self.output_dir, 
            {'sequence_type': 'long'}
        )
        return processor.process_chunk(chunk)
    
    def _check_memory_usage(self):
        """Check and manage memory usage."""
        try:
            import psutil  # type: ignore
            process = psutil.Process()
            memory_mb = process.memory_info().rss / 1024 / 1024
            
            self.stats.peak_memory_mb = max(self.stats.peak_memory_mb, memory_mb)
            
            if memory_mb > self.memory_limit_gb * 1024 * 0.8:  # 80% of limit
                self.logger.warning(f"High memory usage: {memory_mb:.1f}MB")
                gc.collect()
                
        except ImportError:
            pass  # psutil not available
    
    def get_processing_stats(self) -> StreamingStats:
        """Get current processing statistics."""
        return self.stats
    
    def estimate_completion_time(self) -> Optional[float]:
        """Estimate remaining processing time."""
        if self.stats.processed_sequences > 0 and self.stats.processing_time > 0:
            current_rate = self.stats.processed_sequences / self.stats.processing_time
            remaining_sequences = self.stats.total_sequences - self.stats.processed_sequences
            if remaining_sequences > 0:
                return remaining_sequences / current_rate
        return None
    
    def create_processing_report(self) -> Dict[str, Any]:
        """Create comprehensive processing report."""
        return {
            'dataset_info': {
                'total_sequences': self.stats.total_sequences,
                'processed_sequences': self.stats.processed_sequences,
                'processing_rate': f"{self.stats.current_throughput:.2f} seq/s"
            },
            'performance_metrics': {
                'total_time': f"{self.stats.processing_time:.2f}s",
                'chunks_processed': self.stats.chunks_processed,
                'peak_memory_mb': f"{self.stats.peak_memory_mb:.1f}MB",
                'average_chunk_time': f"{self.stats.processing_time / max(1, self.stats.chunks_processed):.2f}s"
            },
            'configuration': {
                'chunk_size': self.chunk_size,
                'max_workers': self.max_workers,
                'memory_limit_gb': self.memory_limit_gb
            }
        }