"""
Enhanced Sequence Processor
--------------------------
Core module for high-throughput, memory-optimized sequence processing with minimap2 integration.
Robust error handling, async support, and quality metrics collection.
"""

import os
import logging
import asyncio
import concurrent.futures
import subprocess
import json
from pathlib import Path
from typing import List, Optional, Dict, Any, Callable, Union, TYPE_CHECKING, cast
from dataclasses import dataclass

# Add import for configuration (only for typing)
if TYPE_CHECKING:
    from ..config.pipeline_config import ProcessingConfig, AlignmentPreset  # pragma: no cover
else:
    ProcessingConfig = Any  # type: ignore
    AlignmentPreset = Any  # type: ignore

try:
    import mappy as mp  # type: ignore
except ImportError:
    mp = None

class EnhancedSequenceProcessorError(Exception):
    """Custom exception for EnhancedSequenceProcessor errors."""
    pass

@dataclass
class ProcessingStats:
    """Statistics for processing operations."""
    total_files: int = 0
    successful_alignments: int = 0
    failed_alignments: int = 0
    total_sequences: int = 0
    processing_time: float = 0.0

class EnhancedSequenceProcessor:
    """
    High-throughput, memory-optimized sequence processor with minimap2 integration.
    Provides robust error handling and quality metrics collection.
    """
    def __init__(self, reference_path: str, threads: int = 4, logger: Optional[logging.Logger] = None, 
                 config: Optional[ProcessingConfig] = None, require_mappy: bool = True):
        self.reference_path = reference_path
        self.threads = threads
        self.logger = logger or self._setup_logger()
        self.stats = ProcessingStats()
        self._aligner = None  # Cache aligner instance
        self.config = config  # Store configuration
        self.require_mappy = require_mappy
        self._check_dependencies()

    def _setup_logger(self) -> logging.Logger:
        logger = logging.getLogger("EnhancedSequenceProcessor")
        if not logger.handlers:
            handler = logging.StreamHandler()
            formatter = logging.Formatter('%(asctime)s %(levelname)s %(name)s: %(message)s')
            handler.setFormatter(formatter)
            logger.addHandler(handler)
        logger.setLevel(logging.INFO)
        return logger

    def _check_dependencies(self):
        # Reference required; mappy optional (we can fallback to minimap2 CLI)
        if not os.path.exists(self.reference_path):
            raise EnhancedSequenceProcessorError(f"Reference file not found: {self.reference_path}")
        if mp is None:
            if self.require_mappy:
                raise EnhancedSequenceProcessorError("mappy is not available; install mappy or set require_mappy=False to enable minimap2 CLI fallback.")
            else:
                self.logger.warning("mappy not available; will fallback to minimap2 CLI if installed.")

    def _get_aligner(self):
        """Get cached aligner instance for better performance."""
        if self._aligner is None:
            try:
                # Use preset from config if available
                preset = "map-ont"  # default
                if self.config and hasattr(self.config, 'preset'):
                    preset = self.config.preset.value if hasattr(self.config.preset, 'value') else str(self.config.preset)
                # type safety for mappy usage
                assert mp is not None, "mappy is required but not available"
                self._aligner = mp.Aligner(self.reference_path, preset=preset, n_threads=self.threads)
                self.logger.info(f"Initialized minimap2 aligner with {self.threads} threads, preset: {preset}")
            except Exception as e:
                raise EnhancedSequenceProcessorError(f"Failed to initialize minimap2 aligner: {e}")
        return self._aligner

    def validate_inputs(self, query_fastas: List[str]) -> List[str]:
        """Validate input files and return list of valid files."""
        valid_files = []
        for fasta in query_fastas:
            if not os.path.exists(fasta):
                self.logger.warning(f"Input file not found: {fasta}")
                continue
            if not self._is_valid_fasta(fasta):
                self.logger.warning(f"Invalid FASTA format: {fasta}")
                continue
            valid_files.append(fasta)
        self.logger.info(f"Validated {len(valid_files)}/{len(query_fastas)} input files")
        return valid_files

    def _is_valid_fasta(self, filepath: str) -> bool:
        """Enhanced FASTA/FASTQ format validation."""
        try:
            with open(filepath, 'r') as f:
                first_line = f.readline().strip()
                # Support both FASTA and FASTQ
                return first_line.startswith('>') or first_line.startswith('@')
        except Exception:
            return False

    def align_sequences(
        self,
        query_fastas: List[str],
        output_dir: str,
        on_success: Optional[Callable[[str, str], None]] = None,
        on_error: Optional[Callable[[str, Exception], None]] = None,
        validate_inputs: bool = True
    ) -> List[str]:
        """
        Aligns a list of query FASTA files to the reference using minimap2 (mappy).
        Returns a list of output alignment file paths.
        Optional callbacks for success and error handling.
        """
        import time
        start_time = time.time()
        
        if validate_inputs:
            query_fastas = self.validate_inputs(query_fastas)
        
        os.makedirs(output_dir, exist_ok=True)
        alignments = []
        self.stats.total_files = len(query_fastas)
        
        for fasta in query_fastas:
            try:
                output_path = os.path.join(output_dir, Path(fasta).stem + ".paf")
                seq_count = self._run_alignment(fasta, output_path)
                alignments.append(output_path)
                self.stats.successful_alignments += 1
                self.stats.total_sequences += seq_count
                self.logger.info(f"Alignment successful for {fasta}, sequences: {seq_count}, output: {output_path}")
                if on_success:
                    on_success(fasta, output_path)
            except Exception as e:
                self.stats.failed_alignments += 1
                self.logger.error(f"Alignment failed for {fasta}: {e}")
                if on_error:
                    on_error(fasta, e)
        
        self.stats.processing_time = time.time() - start_time
        self.logger.info(f"Processing complete: {self.stats.successful_alignments}/{self.stats.total_files} files, "
                        f"{self.stats.total_sequences} sequences in {self.stats.processing_time:.2f}s")
        return alignments

    def _run_alignment(self, query_fasta: str, output_path: str) -> int:
        """
        Runs minimap2 alignment and writes output to file.
        Returns number of sequences processed.
        """
        # If mappy is available, use it; otherwise fallback to minimap2 CLI
        if mp is not None:
            aligner = self._get_aligner()
            seq_count = 0
            with open(query_fasta, 'r') as qf, open(output_path, 'w') as outf:
                for name, seq, qual in mp.fastx_read(qf):
                    seq_count += 1
                    for hit in aligner.map(seq):
                        target_len = getattr(hit, 'ctg_len', None)
                        if target_len is None:
                            target_len = len(hit.ctg) if hasattr(hit, 'ctg') else 0
                        outf.write(
                            f"{name}\t{len(seq)}\t{hit.q_st}\t{hit.q_en}\t{hit.strand}\t"
                            f"{hit.ctg}\t{target_len}\t{hit.r_st}\t{hit.r_en}\t"
                            f"{hit.mlen}\t{hit.blen}\t{hit.mapq}\n"
                        )
            return seq_count
        else:
            # Fallback: use minimap2 CLI via integrations
            try:
                from ..integrations.external_alignment_tools import Minimap2Aligner
                aligner = Minimap2Aligner(self.reference_path)
                success = aligner.align_file(query_fasta, output_path)
                if not success:
                    raise EnhancedSequenceProcessorError("minimap2 CLI alignment failed")
                # quick count sequences
                return sum(1 for line in open(query_fasta, 'r') if line.startswith('>') or line.startswith('@'))
            except Exception as e:
                raise EnhancedSequenceProcessorError(f"Alignment failed (no mappy, minimap2 fallback error): {e}")

    def process_batch_parallel(self, query_fastas: List[str], output_dir: str, max_workers: Optional[int] = None) -> List[str]:
        """
        Process multiple FASTA files in parallel using ThreadPoolExecutor.
        """
        max_workers = max_workers or min(len(query_fastas), self.threads)
        alignments = []
        
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_fasta = {
                executor.submit(self._process_single_file, fasta, output_dir): fasta 
                for fasta in self.validate_inputs(query_fastas)
            }
            
            for future in concurrent.futures.as_completed(future_to_fasta):
                fasta = future_to_fasta[future]
                try:
                    result = future.result()
                    if result:
                        alignments.append(result)
                        self.stats.successful_alignments += 1
                except Exception as e:
                    self.stats.failed_alignments += 1
                    self.logger.error(f"Parallel processing failed for {fasta}: {e}")
        
        return alignments

    def _process_single_file(self, fasta: str, output_dir: str) -> Optional[str]:
        """Process a single FASTA file for parallel execution."""
        try:
            output_path = os.path.join(output_dir, Path(fasta).stem + ".paf")
            self._run_alignment(fasta, output_path)
            return output_path
        except Exception as e:
            self.logger.error(f"Failed to process {fasta}: {e}")
            return None

    def collect_quality_metrics(self, alignment_files: List[str]) -> Dict[str, Any]:
        """
        Collects basic quality metrics from alignment files.
        Returns a dictionary mapping file paths to metrics.
        """
        metrics = {}
        for aln in alignment_files:
            try:
                with open(aln) as f:
                    lines = f.readlines()
                metrics[aln] = {"num_alignments": len(lines)}
            except Exception as e:
                self.logger.error(f"Failed to collect metrics for {aln}: {e}")
                metrics[aln] = {"num_alignments": 0, "error": str(e)}
        return metrics

    def apply_quality_filters(self, alignments: List[str]) -> List[str]:
        """Apply quality filters based on configuration."""
        if not self.config:
            return alignments
        
        filtered_alignments = []
        for aln_file in alignments:
            try:
                filtered_file = aln_file.replace('.paf', '_filtered.paf')
                with open(aln_file, 'r') as infile, open(filtered_file, 'w') as outfile:
                    for line in infile:
                        fields = line.strip().split('\t')
                        if len(fields) >= 12:
                            mapq = int(fields[11])
                            if mapq >= getattr(self.config, 'quality_threshold', 20):
                                outfile.write(line)
                filtered_alignments.append(filtered_file)
            except Exception as e:
                self.logger.error(f"Failed to filter {aln_file}: {e}")
                filtered_alignments.append(aln_file)  # Keep original if filtering fails
        
        return filtered_alignments

    def get_processing_stats(self) -> ProcessingStats:
        """Return current processing statistics."""
        return self.stats

    def reset_stats(self):
        """Reset processing statistics."""
        self.stats = ProcessingStats()

    def run_rgi(self, genome_fasta: str, output_dir: str, rgi_exec: str = "rgi", rgi_db: Optional[str] = None) -> Optional[str]:
        """
        Runs RGI on the given genome FASTA file to predict AMR genes and extract coordinates.
        Returns the path to the RGI output JSON file.
        """
        os.makedirs(output_dir, exist_ok=True)
        output_json = os.path.join(output_dir, f"{Path(genome_fasta).stem}_rgi.json")
        db_arg = f"--card_json {rgi_db}" if rgi_db else ""
        cmd = f'{rgi_exec} main --input_sequence "{genome_fasta}" --output_file "{output_json}" --input_type contig {db_arg} --local --clean --json'
        self.logger.info(f"Running RGI: {cmd}")
        try:
            result = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            self.logger.info(f"RGI completed for {genome_fasta}")
            return output_json
        except subprocess.CalledProcessError as e:
            self.logger.error(f"RGI failed for {genome_fasta}: {e.stderr.decode().strip()}")
            return None

    def extract_gene_coordinates(self, rgi_json: str, target_genes: Optional[List[str]] = None) -> List[Dict[str, Any]]:
        """
        Parses RGI output JSON and extracts coordinates for target genes.
        Returns a list of dicts with gene name, contig, start, end, and strand.
        """
        if not os.path.exists(rgi_json):
            self.logger.error(f"RGI output not found: {rgi_json}")
            return []
        try:
            with open(rgi_json, "r") as f:
                rgi_data = json.load(f)
            coords = []
            for entry in rgi_data if isinstance(rgi_data, list) else rgi_data.get("hits", []):
                gene = entry.get("best_hit_aro", "")
                if target_genes and gene not in target_genes:
                    continue
                coords.append({
                    "gene": gene,
                    "contig": entry.get("contig", ""),
                    "start": entry.get("start", 0),
                    "end": entry.get("stop", 0),
                    "strand": entry.get("strand", ""),
                    "product": entry.get("drug_class", ""),
                })
            self.logger.info(f"Extracted {len(coords)} gene coordinates from RGI output")
            return coords
        except Exception as e:
            self.logger.error(f"Failed to parse RGI output {rgi_json}: {e}")
            return []

    # Optional: Async support stub for future extension
    async def align_sequences_async(self, query_fastas: List[str], output_dir: str) -> List[str]:
        """
        Async version of align_sequences (stub for future implementation).
        """
        # ...implement async logic if needed...
        raise NotImplementedError("Async alignment not yet implemented.")
