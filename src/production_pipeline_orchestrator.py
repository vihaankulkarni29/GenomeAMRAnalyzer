#!/usr/bin/env python3
"""
Production Pipeline Orchestrator with Batch Integrity
Complete traceability and manifest system for large-scale genome processing

This module provides a senior bioinformatician-level orchestration system that:
1. Coordinates URL discovery → genome download → RGI analysis → protein extraction
2. Maintains strict accession-based naming across all pipeline stages
3. Generates comprehensive batch manifests with full provenance tracking
4. Provides integrity checking and validation across all processing steps
5. Handles 100-1000 genome batches with bulletproof traceability

Author: GenomeAMRAnalyzer Pipeline
Version: 2.0 - Senior Bioinformatician Grade
"""

import os
import sys
import time
import json
import logging
import hashlib
import asyncio
import multiprocessing
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Set, Any
from dataclasses import dataclass, asdict
import yaml

# Progress tracking
try:
    from tqdm import tqdm
    TQDM_AVAILABLE = True
except ImportError:
    TQDM_AVAILABLE = False
    # Fallback progress class
    class tqdm:
        def __init__(self, iterable=None, total=None, **kwargs):
            self.total = total
            self.n = 0
            
        def update(self, n=1):
            self.n += n
            
        def set_postfix(self, ordered_dict=None, **kwargs):
            pass
            
        def close(self):
            pass

# Import our production components
sys.path.append(str(Path(__file__).parent))
from core.url_to_genomes_workflow import URLToGenomesWorkflow
from core.genome_downloader import GenomeDownloader
from abricate_runner import run_abricate, ScanResults
from abricate_to_coords import convert_abricate_to_coords
from production_fasta_extractor import ProductionFastaExtractor


@dataclass
class PipelineManifest:
    """Complete pipeline manifest with multi-database traceability"""
    pipeline_version: str
    generated_timestamp: str
    pipeline_id: str
    source_url: str
    
    # Batch information
    total_genomes_discovered: int
    total_genomes_downloaded: int
    total_genomes_rgi_analyzed: int
    total_proteins_extracted: int
    
    # Multi-database tracking
    database_scan_results: Dict[str, Dict[str, Any]]  # database -> scan summary
    card_success_rate: float  # Critical database
    vfdb_success_rate: float  # Secondary database
    plasmidfinder_success_rate: float  # Secondary database
    
    # Processing statistics
    download_success_rate: float
    rgi_success_rate: float
    extraction_success_rate: float
    
    # File manifests
    genome_manifest_path: str
    coordinate_manifest_path: str
    extraction_manifest_path: str
    
    # Integrity checksums
    pipeline_checksum: str
    manifest_checksum: str
    
    # Accession tracking
    successful_accessions: List[str]
    failed_accessions: Dict[str, str]  # accession -> failure reason
    
    # Target genes
    target_genes: List[str]


@dataclass
class GenomeWorkItem:
    """Input work item for single-genome processing"""
    genome_id: str
    genome_file_path: str
    target_genes: List[str]
    output_dirs: Dict[str, str]  # Directory paths for outputs
    config: Dict[str, Any]       # Processing configuration
    pipeline_id: str


@dataclass 
class GenomeProcessingResult:
    """Result of processing a single genome"""
    genome_id: str
    success: bool
    error_message: Optional[str] = None
    
    # Processing stages completed
    abricate_card_success: bool = False
    abricate_vfdb_success: bool = False
    abricate_plasmidfinder_success: bool = False
    coordinates_conversion_success: bool = False
    protein_extraction_success: bool = False
    
    # Output file paths
    card_output_file: Optional[str] = None
    vfdb_output_file: Optional[str] = None
    plasmidfinder_output_file: Optional[str] = None
    coordinates_file: Optional[str] = None
    protein_file: Optional[str] = None
    
    # Statistics
    coordinates_found: int = 0
    proteins_extracted: int = 0
    processing_time: float = 0.0


def process_single_genome(work_item: GenomeWorkItem) -> Optional[GenomeProcessingResult]:
    """
    Process a single genome through the complete AMR analysis pipeline.
    
    This function encapsulates all per-genome processing steps:
    1. Abricate scanning against multiple databases (CARD, VFDB, PlasmidFinder)
    2. Coordinate conversion for CARD results
    3. Protein extraction using coordinates
    
    Args:
        work_item: GenomeWorkItem containing all necessary input data
        
    Returns:
        GenomeProcessingResult if successful, None if failed
    """
    start_time = time.time()
    genome_id = work_item.genome_id
    
    # Initialize result object
    result = GenomeProcessingResult(
        genome_id=genome_id,
        success=False
    )
    
    try:
        # Setup logging for this worker process
        logger = logging.getLogger(f"genome_worker_{genome_id}")
        logger.setLevel(logging.INFO)
        
        # Create a process-specific log handler if it doesn't exist
        if not logger.handlers:
            log_file = Path(work_item.output_dirs['logs']) / f"{genome_id}_worker.log"
            handler = logging.FileHandler(log_file)
            formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
            handler.setFormatter(formatter)
            logger.addHandler(handler)
        
        logger.info(f"Starting processing for genome {genome_id}")
        
        # Validate input files exist
        genome_file = Path(work_item.genome_file_path)
        if not genome_file.exists():
            raise FileNotFoundError(f"Genome file not found: {work_item.genome_file_path}")
        
        # Create temporary directory for this genome's intermediate files
        temp_genome_dir = Path(work_item.output_dirs['temp']) / genome_id
        temp_genome_dir.mkdir(parents=True, exist_ok=True)
        
        # Copy genome file to temp directory for isolated processing
        temp_genome_file = temp_genome_dir / genome_file.name
        import shutil
        shutil.copy2(genome_file, temp_genome_file)
        
        # Step 1: Abricate scanning against multiple databases
        databases_to_scan = [
            {'name': 'card', 'critical': True, 'description': 'CARD AMR Database'},
            {'name': 'vfdb', 'critical': False, 'description': 'VFDB Virulence Factors'},
            {'name': 'plasmidfinder', 'critical': False, 'description': 'PlasmidFinder Database'}
        ]
        
        abricate_results = {}
        
        for db_config in databases_to_scan:
            db_name = db_config['name']
            is_critical = db_config['critical']
            
            try:
                logger.info(f"Running abricate scan for {db_name} database")
                
                # Get database-specific output directory
                db_output_dir = work_item.output_dirs[f'abricate_{db_name}']
                
                # Run abricate for single genome against this database
                scan_results = run_abricate(str(temp_genome_dir), db_output_dir, db=db_name)
                
                # Store results
                abricate_results[db_name] = scan_results
                
                # Check if this genome was successfully processed
                if genome_id in [Path(f).stem.replace(f'_{db_name}', '') for f in scan_results.output_files]:
                    logger.info(f"Successfully scanned {genome_id} against {db_name}")
                    
                    # Set success flags and file paths
                    if db_name == 'card':
                        result.abricate_card_success = True
                        result.card_output_file = str(next(
                            f for f in scan_results.output_files 
                            if Path(f).stem.replace(f'_{db_name}', '') == genome_id
                        ))
                    elif db_name == 'vfdb':
                        result.abricate_vfdb_success = True
                        result.vfdb_output_file = str(next(
                            f for f in scan_results.output_files 
                            if Path(f).stem.replace(f'_{db_name}', '') == genome_id
                        ))
                    elif db_name == 'plasmidfinder':
                        result.abricate_plasmidfinder_success = True
                        result.plasmidfinder_output_file = str(next(
                            f for f in scan_results.output_files 
                            if Path(f).stem.replace(f'_{db_name}', '') == genome_id
                        ))
                        
                else:
                    logger.warning(f"No results found for {genome_id} in {db_name} database")
                    if is_critical:
                        raise RuntimeError(f"Critical database {db_name} scan failed for {genome_id}")
                        
            except Exception as e:
                logger.error(f"Abricate scan failed for {db_name}: {e}")
                if is_critical:
                    raise RuntimeError(f"Critical database {db_name} scan failed: {e}")
                else:
                    logger.warning(f"Non-critical database {db_name} scan failed, continuing: {e}")
        
        # Step 2: Convert CARD results to coordinates (only if CARD scan was successful)
        if result.abricate_card_success and result.card_output_file:
            try:
                logger.info(f"Converting CARD results to coordinates for {genome_id}")
                
                coords_output_file = Path(work_item.output_dirs['coordinates']) / f"{genome_id}_coordinates.csv"
                
                # Convert abricate results to coordinate format
                rows_converted = convert_abricate_to_coords(
                    Path(result.card_output_file), 
                    coords_output_file
                )
                
                if rows_converted > 0:
                    result.coordinates_conversion_success = True
                    result.coordinates_file = str(coords_output_file)
                    result.coordinates_found = rows_converted
                    logger.info(f"Successfully converted {rows_converted} coordinates for {genome_id}")
                else:
                    logger.warning(f"No coordinates found for {genome_id}")
                    result.coordinates_conversion_success = True  # Not an error, just no hits
                    result.coordinates_file = str(coords_output_file)
                    
            except Exception as e:
                logger.error(f"Coordinate conversion failed for {genome_id}: {e}")
                raise RuntimeError(f"Coordinate conversion failed: {e}")
        
        # Step 3: Protein extraction (only if coordinates were successfully generated)
        if result.coordinates_conversion_success and result.coordinates_file:
            try:
                logger.info(f"Extracting proteins for {genome_id}")
                
                # Initialize protein extractor for single genome
                extractor = ProductionFastaExtractor(work_item.output_dirs['proteins'])
                
                # Process single genome protein extraction
                # We need to create a mini-batch with just this genome
                single_genome_dir = temp_genome_dir
                single_coords_dir = Path(work_item.output_dirs['coordinates'])
                
                # Run extraction for this single genome
                extraction_success = extractor.process_batch_from_coordinates(
                    str(single_genome_dir),
                    str(single_coords_dir), 
                    work_item.target_genes
                )
                
                if extraction_success and genome_id in extractor.extraction_results:
                    result.protein_extraction_success = True
                    result.proteins_extracted = len([
                        r for r in extractor.extraction_results[genome_id] 
                        if r.extraction_success
                    ])
                    result.protein_file = str(
                        Path(work_item.output_dirs['proteins']) / "proteins" / f"{genome_id}_proteins.fasta"
                    )
                    logger.info(f"Successfully extracted {result.proteins_extracted} proteins for {genome_id}")
                else:
                    logger.warning(f"Protein extraction completed but no proteins found for {genome_id}")
                    
            except Exception as e:
                logger.error(f"Protein extraction failed for {genome_id}: {e}")
                # Protein extraction failure is not necessarily critical if coordinates exist
                logger.warning(f"Continuing despite protein extraction failure for {genome_id}")
        
        # Cleanup temporary files
        try:
            shutil.rmtree(temp_genome_dir)
        except Exception as e:
            logger.warning(f"Failed to cleanup temp directory for {genome_id}: {e}")
        
        # Calculate processing time
        result.processing_time = time.time() - start_time
        
        # Determine overall success
        # Success criteria: At least CARD scan and coordinate conversion must succeed
        result.success = result.abricate_card_success and result.coordinates_conversion_success
        
        if result.success:
            logger.info(f"Successfully completed processing for {genome_id} in {result.processing_time:.2f}s")
        else:
            logger.warning(f"Processing completed with limited success for {genome_id}")
            
        return result
        
    except Exception as e:
        # Handle any unexpected errors
        result.success = False
        result.error_message = str(e)
        result.processing_time = time.time() - start_time
        
        # Log error
        try:
            logger = logging.getLogger(f"genome_worker_{genome_id}")
            logger.error(f"Processing failed for {genome_id}: {e}")
        except:
            # Fallback logging if logger setup failed
            print(f"ERROR: Processing failed for {genome_id}: {e}")
        
        return result


class ProductionPipelineOrchestrator:
    """
    Senior bioinformatician-level pipeline orchestrator
    Handles complete URL-to-proteins workflow with bulletproof traceability
    """
    
    def __init__(self, output_base_dir: str, pipeline_config: Dict[str, Any]):
        """Initialize production pipeline orchestrator"""
        self.output_base_dir = Path(output_base_dir)
        self.output_base_dir.mkdir(parents=True, exist_ok=True)
        
        self.config = pipeline_config
        self.pipeline_id = f"pipeline_{int(time.time())}"
        
        # Create structured output directories with multi-database support
        self.directories = {
            'genomes': self.output_base_dir / "genomes",
            'coordinates': self.output_base_dir / "coordinates", 
            'abricate_card': self.output_base_dir / "abricate_card",
            'abricate_vfdb': self.output_base_dir / "abricate_vfdb", 
            'abricate_plasmidfinder': self.output_base_dir / "abricate_plasmidfinder",
            'proteins': self.output_base_dir / "proteins",
            'manifests': self.output_base_dir / "manifests",
            'logs': self.output_base_dir / "logs",
            'reports': self.output_base_dir / "reports",
            'temp': self.output_base_dir / "temp"  # Temporary directory for parallel processing
        }
        
        for dir_path in self.directories.values():
            dir_path.mkdir(parents=True, exist_ok=True)
        
        # Setup logging
        self._setup_logging()
        
        # Initialize tracking with multi-database support
        self.accession_registry: Dict[str, Dict[str, Any]] = {}
        self.processing_status: Dict[str, str] = {}  # accession -> status
        self.failure_registry: Dict[str, str] = {}  # accession -> failure reason
        self.database_results: Dict[str, ScanResults] = {}  # database -> scan results
        
        # Performance metrics
        self.performance_metrics = {
            'start_time': time.time(),
            'stage_times': {},
            'memory_usage': [],
            'processing_rates': {}
        }
        
        self.logger.info(f"Production pipeline orchestrator initialized: {self.pipeline_id}")
    
    def _setup_logging(self):
        """Setup comprehensive logging system"""
        log_file = self.directories['logs'] / f"{self.pipeline_id}.log"
        
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - [%(filename)s:%(lineno)d] - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        
        self.logger = logging.getLogger('ProductionPipeline')
        self.logger.setLevel(logging.INFO)
        
        # File handler
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(formatter)
        self.logger.addHandler(file_handler)
        
        # Console handler
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(formatter)
        self.logger.addHandler(console_handler)
    
    async def execute_complete_pipeline(self, source_url: str, target_genes: List[str]) -> PipelineManifest:
        """Execute complete URL-to-proteins pipeline with full traceability"""
        
        self.logger.info(f"Starting complete pipeline for: {source_url}")
        self.logger.info(f"Target genes: {', '.join(target_genes)}")
        
        pipeline_start = time.time()
        
        try:
            # Stage 1: URL Discovery and Genome Download
            stage_start = time.time()
            self.logger.info("=== STAGE 1: URL Discovery and Genome Download ===")
            
            genome_results = await self._execute_genome_discovery_and_download(source_url)
            
            self.performance_metrics['stage_times']['genome_download'] = time.time() - stage_start
            
            if not genome_results:
                raise RuntimeError("No genomes successfully downloaded")
            
            # Stage 2 & 3: Parallel Genome Processing (replaces sequential multi-database analysis and protein extraction)
            stage_start = time.time()
            self.logger.info("=== STAGE 2-3: Parallel Genome Processing ===")
            
            parallel_results = await self._execute_parallel_genome_processing(target_genes)
            
            self.performance_metrics['stage_times']['parallel_processing'] = time.time() - stage_start
            
            # Stage 4: Generate Master Manifest
            stage_start = time.time()
            self.logger.info("=== STAGE 4: Master Manifest Generation ===")
            
            master_manifest = self._generate_master_manifest(
                source_url, target_genes, genome_results, parallel_results
            )
            
            self.performance_metrics['stage_times']['manifest_generation'] = time.time() - stage_start
            
            # Final pipeline metrics
            total_time = time.time() - pipeline_start
            self.performance_metrics['total_pipeline_time'] = total_time
            
            self.logger.info(f"Pipeline completed successfully in {total_time:.2f} seconds")
            self.logger.info(f"Processed {len(master_manifest.successful_accessions)} accessions successfully")
            
            return master_manifest
            
        except Exception as e:
            self.logger.error(f"Pipeline failed: {e}")
            raise
    
    async def _execute_genome_discovery_and_download(self, source_url: str) -> Dict[str, Any]:
        """Execute genome discovery and download with accession tracking.

        This builds a temporary YAML config file for URLToGenomesWorkflow,
        then registers downloaded genomes into the accession registry.
        """
        
        try:
            # Build a temporary config file for the URLToGenomesWorkflow
            import yaml as _yaml
            tmp_cfg_path = self.directories['manifests'] / f"{self.pipeline_id}_url_to_genomes.yaml"
            tmp_cfg = {
                'directories': {
                    'genomes': str(self.directories['genomes']),
                    'logs': str(self.directories['logs']),
                },
                'ncbi_email': self.config.get('email', 'pipeline@example.com'),
                'ncbi_api_key': self.config.get('api_key'),
                'ncbi_search': {
                    'max_genomes': self.config.get('ncbi_search', {}).get('max_genomes', 100),
                    'retry_attempts': self.config.get('ncbi_search', {}).get('retry_attempts', 3),
                    'timeout_seconds': self.config.get('ncbi_search', {}).get('timeout_seconds', 60),
                }
            }
            with open(tmp_cfg_path, 'w') as _f:
                _yaml.safe_dump(tmp_cfg, _f)

            # Initialize URL workflow
            workflow = URLToGenomesWorkflow(str(tmp_cfg_path))
            
            # Execute workflow: returns (list_of_files, report)
            files, report = await workflow.run_complete_workflow(source_url)
            
            # Register accessions from downloaded files
            for file_path in files:
                fp = Path(file_path)
                genome_id = fp.stem
                self.accession_registry[genome_id] = {
                    'genome_file': str(fp),
                    'organism': '',
                    'genome_length': 0,
                    'download_timestamp': '',
                    'file_checksum': '',
                    'download_success': True
                }
                self.processing_status[genome_id] = 'genome_downloaded'
            
            self.logger.info(f"Downloaded {len(files)} genomes successfully")
            
            return {
                'files': files,
                'report': report,
            }
            
        except Exception as e:
            self.logger.error(f"Genome discovery/download failed: {e}")
            raise
    
    async def _execute_parallel_genome_processing(self, target_genes: List[str]) -> Dict[str, Any]:
        """Execute parallel genome processing using multiprocessing Pool.
        
        This method replaces the sequential stages 2 and 3 with parallel processing:
        - Multi-database analysis (CARD, VFDB, PlasmidFinder) 
        - Coordinate conversion
        - Protein extraction
        
        Each genome is processed completely in parallel, utilizing all available CPU cores.
        """
        
        try:
            self.logger.info("=== PARALLEL GENOME PROCESSING ===")
            self.logger.info(f"Processing genomes in parallel using {self.config.get('threads', 1)} threads")
            
            # Prepare work items for each downloaded genome
            work_items = []
            genomes_dir = Path(self.directories['genomes'])
            
            # Find all downloaded genome files
            genome_files = list(genomes_dir.glob("*.fasta")) + list(genomes_dir.glob("*.fa")) + list(genomes_dir.glob("*.fna"))
            
            if not genome_files:
                raise RuntimeError("No genome files found for processing")
            
            self.logger.info(f"Found {len(genome_files)} genome files for parallel processing")
            
            # Create work items for parallel processing
            for genome_file in genome_files:
                genome_id = genome_file.stem
                
                # Create work item with all necessary information
                work_item = GenomeWorkItem(
                    genome_id=genome_id,
                    genome_file_path=str(genome_file),
                    target_genes=target_genes,
                    output_dirs={
                        'abricate_card': str(self.directories['abricate_card']),
                        'abricate_vfdb': str(self.directories['abricate_vfdb']),
                        'abricate_plasmidfinder': str(self.directories['abricate_plasmidfinder']),
                        'coordinates': str(self.directories['coordinates']),
                        'proteins': str(self.directories['proteins']),
                        'logs': str(self.directories['logs']),
                        'temp': str(self.directories['temp'])
                    },
                    config=self.config.copy(),
                    pipeline_id=self.pipeline_id
                )
                work_items.append(work_item)
            
            # Execute parallel processing with progress tracking
            num_threads = self.config.get('threads', 1)
            processing_start = time.time()
            
            self.logger.info(f"Starting parallel processing of {len(work_items)} genomes with {num_threads} workers")
            
            # Use multiprocessing Pool with progress tracking
            successful_results = []
            failed_results = []
            
            with multiprocessing.Pool(processes=num_threads) as pool:
                # Use imap_unordered for real-time progress tracking
                if TQDM_AVAILABLE:
                    # Create progress bar with detailed information
                    progress_bar = tqdm(
                        total=len(work_items),
                        desc="Processing genomes",
                        unit="genome",
                        ncols=100,
                        bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}, {rate_fmt}]",
                        leave=True
                    )
                    
                    # Process with progress updates
                    for result in pool.imap_unordered(process_single_genome, work_items):
                        if result and result.success:
                            successful_results.append(result)
                            progress_bar.set_postfix({
                                'Success': len(successful_results), 
                                'Failed': len(failed_results),
                                'Current': f"{result.genome_id}"
                            })
                        else:
                            failed_results.append(result)
                            progress_bar.set_postfix({
                                'Success': len(successful_results), 
                                'Failed': len(failed_results),
                                'Current': f"{result.genome_id if result else 'Unknown'} (FAILED)"
                            })
                        progress_bar.update(1)
                    
                    progress_bar.close()
                    
                else:
                    # Fallback without progress bar
                    self.logger.info("Processing genomes (no progress bar - install tqdm for progress tracking)")
                    for i, result in enumerate(pool.imap_unordered(process_single_genome, work_items), 1):
                        if result and result.success:
                            successful_results.append(result)
                        else:
                            failed_results.append(result)
                        
                        # Log progress every 10% or every 5 genomes, whichever is smaller
                        log_interval = max(1, min(5, len(work_items) // 10))
                        if i % log_interval == 0 or i == len(work_items):
                            progress_pct = (i / len(work_items)) * 100
                            self.logger.info(f"Progress: {i}/{len(work_items)} ({progress_pct:.1f}%) - "
                                           f"Success: {len(successful_results)}, Failed: {len(failed_results)}")
            
            processing_time = time.time() - processing_start
            self.logger.info(f"Parallel processing completed in {processing_time:.2f} seconds")
            
            # Filter and aggregate results (already done above during processing)
            # successful_results and failed_results are already populated
            
            # Enhanced result analysis and reporting
            result_analysis = self._analyze_processing_results(successful_results, failed_results, processing_time)
            
            # Log comprehensive summary statistics
            self._log_processing_summary(result_analysis)
            
            # Update internal tracking with results
            self._update_internal_tracking(successful_results, failed_results)
            
            return result_analysis
            
        except Exception as e:
            self.logger.error(f"Parallel genome processing failed: {e}")
            raise
    
    def _analyze_processing_results(self, successful_results: List, failed_results: List, processing_time: float) -> Dict[str, Any]:
        """Analyze and aggregate parallel processing results with detailed statistics"""
        
        total_genomes = len(successful_results) + len(failed_results)
        
        if total_genomes == 0:
            return {
                'successful_results': [],
                'failed_results': [],
                'total_processing_time': processing_time,
                'success_rate': 0.0,
                'total_genomes': 0
            }
        
        # Calculate comprehensive statistics
        avg_processing_time = sum(r.processing_time for r in successful_results) / len(successful_results) if successful_results else 0
        total_coordinates = sum(r.coordinates_found for r in successful_results)
        total_proteins = sum(r.proteins_extracted for r in successful_results)
        
        # Database-specific success rates
        card_successes = sum(1 for r in successful_results if r.abricate_card_success)
        vfdb_successes = sum(1 for r in successful_results if r.abricate_vfdb_success)
        plasmidfinder_successes = sum(1 for r in successful_results if r.abricate_plasmidfinder_success)
        coordinates_successes = sum(1 for r in successful_results if r.coordinates_conversion_success)
        protein_successes = sum(1 for r in successful_results if r.protein_extraction_success)
        
        # Processing time statistics
        if successful_results:
            processing_times = [r.processing_time for r in successful_results]
            min_time = min(processing_times)
            max_time = max(processing_times)
            median_time = sorted(processing_times)[len(processing_times) // 2]
        else:
            min_time = max_time = median_time = 0
        
        # Failure analysis
        failure_reasons = {}
        for result in failed_results:
            if result and result.error_message:
                reason = result.error_message.split(':')[0] if ':' in result.error_message else result.error_message
                failure_reasons[reason] = failure_reasons.get(reason, 0) + 1
        
        # Parallelization efficiency
        theoretical_sequential_time = len(successful_results) * avg_processing_time
        speedup_factor = theoretical_sequential_time / processing_time if processing_time > 0 else 0
        
        return {
            'successful_results': successful_results,
            'failed_results': failed_results,
            'total_processing_time': processing_time,
            'avg_processing_time': avg_processing_time,
            'total_coordinates': total_coordinates,
            'total_proteins': total_proteins,
            'success_rate': len(successful_results) / total_genomes,
            'total_genomes': total_genomes,
            
            # Database-specific statistics
            'database_success_rates': {
                'card': card_successes / total_genomes,
                'vfdb': vfdb_successes / total_genomes,
                'plasmidfinder': plasmidfinder_successes / total_genomes,
                'coordinates': coordinates_successes / total_genomes,
                'proteins': protein_successes / total_genomes
            },
            
            # Processing time statistics
            'processing_time_stats': {
                'min': min_time,
                'max': max_time,
                'average': avg_processing_time,
                'median': median_time
            },
            
            # Failure analysis
            'failure_analysis': {
                'total_failures': len(failed_results),
                'failure_rate': len(failed_results) / total_genomes,
                'failure_reasons': failure_reasons
            },
            
            # Performance metrics
            'performance_metrics': {
                'speedup_factor': speedup_factor,
                'theoretical_sequential_time': theoretical_sequential_time,
                'efficiency': speedup_factor / self.config.get('threads', 1) if self.config.get('threads', 1) > 0 else 0
            }
        }
    
    def _log_processing_summary(self, analysis: Dict[str, Any]) -> None:
        """Log comprehensive processing summary with detailed statistics"""
        
        self.logger.info("=" * 60)
        self.logger.info("PARALLEL PROCESSING RESULTS SUMMARY")
        self.logger.info("=" * 60)
        
        # Overall statistics
        self.logger.info(f"Processing Overview:")
        self.logger.info(f"  Total genomes: {analysis['total_genomes']}")
        self.logger.info(f"  Successful: {len(analysis['successful_results'])}")
        self.logger.info(f"  Failed: {len(analysis['failed_results'])}")
        self.logger.info(f"  Success rate: {analysis['success_rate']*100:.1f}%")
        
        # Database-specific results
        db_rates = analysis['database_success_rates']
        self.logger.info(f"Database Success Rates:")
        self.logger.info(f"  CARD (critical): {db_rates['card']*100:.1f}%")
        self.logger.info(f"  VFDB (virulence): {db_rates['vfdb']*100:.1f}%")
        self.logger.info(f"  PlasmidFinder: {db_rates['plasmidfinder']*100:.1f}%")
        self.logger.info(f"  Coordinates: {db_rates['coordinates']*100:.1f}%")
        self.logger.info(f"  Protein extraction: {db_rates['proteins']*100:.1f}%")
        
        # Output statistics
        self.logger.info(f"Output Statistics:")
        self.logger.info(f"  Total coordinates found: {analysis['total_coordinates']}")
        self.logger.info(f"  Total proteins extracted: {analysis['total_proteins']}")
        self.logger.info(f"  Avg coordinates per genome: {analysis['total_coordinates']/len(analysis['successful_results']):.1f}" if analysis['successful_results'] else "  Avg coordinates per genome: 0")
        self.logger.info(f"  Avg proteins per genome: {analysis['total_proteins']/len(analysis['successful_results']):.1f}" if analysis['successful_results'] else "  Avg proteins per genome: 0")
        
        # Performance metrics
        perf = analysis['performance_metrics']
        time_stats = analysis['processing_time_stats']
        self.logger.info(f"Performance Metrics:")
        self.logger.info(f"  Total processing time: {analysis['total_processing_time']:.2f}s")
        self.logger.info(f"  Average time per genome: {time_stats['average']:.2f}s")
        self.logger.info(f"  Processing time range: {time_stats['min']:.2f}s - {time_stats['max']:.2f}s")
        self.logger.info(f"  Parallelization speedup: {perf['speedup_factor']:.1f}x")
        self.logger.info(f"  Parallel efficiency: {perf['efficiency']*100:.1f}%")
        
        # Failure analysis
        if analysis['failed_results']:
            failure = analysis['failure_analysis']
            self.logger.warning(f"Failure Analysis:")
            self.logger.warning(f"  Failure rate: {failure['failure_rate']*100:.1f}%")
            self.logger.warning(f"  Common failure reasons:")
            for reason, count in failure['failure_reasons'].items():
                self.logger.warning(f"    {reason}: {count} genomes")
        
        self.logger.info("=" * 60)
    
    def _update_internal_tracking(self, successful_results: List, failed_results: List) -> None:
        """Update internal tracking structures with parallel processing results"""
        
        for result in successful_results:
            genome_id = result.genome_id
            
            # Update processing status
            if result.protein_extraction_success:
                self.processing_status[genome_id] = 'proteins_extracted'
            elif result.coordinates_conversion_success:
                self.processing_status[genome_id] = 'rgi_analyzed'
            else:
                self.processing_status[genome_id] = 'partial_success'
            
            # Update accession registry
            if genome_id not in self.accession_registry:
                self.accession_registry[genome_id] = {}
            
            registry_entry = self.accession_registry[genome_id]
            registry_entry['rgi_coordinates_found'] = result.coordinates_found
            registry_entry['proteins_extracted'] = result.proteins_extracted
            registry_entry['processing_time'] = result.processing_time
            
            if result.coordinates_file:
                registry_entry['coordinate_file'] = result.coordinates_file
            if result.protein_file:
                registry_entry['protein_file'] = result.protein_file
        
        # Update failure registry
        for result in failed_results:
            if result:
                genome_id = result.genome_id
                self.processing_status[genome_id] = 'processing_failed'
                self.failure_registry[genome_id] = result.error_message or "Unknown processing failure"
    
    async def _execute_multi_database_analysis(self, target_genes: List[str]) -> Dict[str, Any]:
        """Execute comprehensive multi-database genomic profiling with failure-proof design.
        
        Scans genomes against CARD (critical), VFDB (virulence factors), and PlasmidFinder (plasmids).
        CARD failures are critical and stop the pipeline. Secondary database failures generate warnings only.
        """
        
        try:
            genomes_dir = str(self.directories['genomes'])
            coords_dir = self.directories['coordinates']
            
            # Define databases to scan with criticality levels
            databases_to_scan = [
                {'name': 'card', 'critical': True, 'description': 'CARD AMR Database'},
                {'name': 'vfdb', 'critical': False, 'description': 'VFDB Virulence Factors'},
                {'name': 'plasmidfinder', 'critical': False, 'description': 'PlasmidFinder Database'}
            ]
            
            self.logger.info(f"Starting multi-database scanning across {len(databases_to_scan)} databases")
            
            # Track overall processing results
            processed_accessions: List[str] = []
            failed_accessions: List[str] = []
            coord_manifest: Dict[str, str] = {}
            
            # Scan each database
            for db_config in databases_to_scan:
                db_name = db_config['name']
                is_critical = db_config['critical']
                description = db_config['description']
                
                self.logger.info(f"Scanning {description} ({'CRITICAL' if is_critical else 'SECONDARY'})")
                
                try:
                    # Get database-specific output directory
                    abricate_out = str(self.directories[f'abricate_{db_name}'])
                    
                    # Run abricate for this database
                    scan_results = run_abricate(genomes_dir, abricate_out, db=db_name)
                    
                    # Store scan results for reporting
                    self.database_results[db_name] = scan_results
                    
                    # Log scan summary
                    total_genomes = scan_results.total_genomes
                    successful = len(scan_results.successful_scans)
                    failed = len(scan_results.failed_genomes)
                    no_hits = len(scan_results.no_hits_genomes)
                    
                    self.logger.info(
                        f"{description}: {successful} hits, {no_hits} no hits, {failed} failed "
                        f"({total_genomes} total genomes)"
                    )
                    
                    # Handle critical database failures
                    if is_critical and scan_results.global_error:
                        error_msg = f"CRITICAL: {description} scan failed: {scan_results.global_error}"
                        self.logger.error(error_msg)
                        raise RuntimeError(error_msg)
                    
                    # For CARD database, convert to coordinates (only critical database needs coordinates)
                    if db_name == 'card':
                        self.logger.info("Converting CARD results to coordinate format for protein extraction")
                        
                        for output_file in scan_results.output_files:
                            genome_id = output_file.stem.replace(f'_{db_name}', '')
                            out_csv = coords_dir / f"{genome_id}_coordinates.csv"
                            
                            try:
                                rows = convert_abricate_to_coords(output_file, out_csv)
                                if rows > 0:
                                    if genome_id not in processed_accessions:
                                        processed_accessions.append(genome_id)
                                    coord_manifest[genome_id] = str(out_csv)
                                    self.processing_status[genome_id] = 'rgi_analyzed'
                                    self.accession_registry.setdefault(genome_id, {})['rgi_coordinates_found'] = rows
                                    self.accession_registry[genome_id]['coordinate_file'] = str(out_csv)
                                else:
                                    # No coordinates but still processed
                                    if genome_id not in processed_accessions:
                                        processed_accessions.append(genome_id)
                                    coord_manifest[genome_id] = str(out_csv)
                                    self.processing_status[genome_id] = 'rgi_analyzed'
                                    self.accession_registry.setdefault(genome_id, {})['rgi_coordinates_found'] = 0
                                    self.accession_registry[genome_id]['coordinate_file'] = str(out_csv)
                            except Exception as e:
                                self.logger.warning(f"Coordinate conversion failed for {genome_id}: {e}")
                                if genome_id not in failed_accessions:
                                    failed_accessions.append(genome_id)
                                self.processing_status[genome_id] = 'rgi_failed'
                                self.failure_registry[genome_id] = f'CARD coordinate conversion failed: {str(e)}'
                        
                        # Update accession tracking with CARD-specific failures
                        for genome_id, error in scan_results.failed_genomes.items():
                            if genome_id not in failed_accessions:
                                failed_accessions.append(genome_id)
                            self.processing_status[genome_id] = 'rgi_failed'
                            self.failure_registry[genome_id] = f'CARD scan failed: {error}'
                    
                except Exception as e:
                    if is_critical:
                        # Critical database failure - stop pipeline
                        error_msg = f"CRITICAL: {description} scan failed: {str(e)}"
                        self.logger.error(error_msg)
                        raise RuntimeError(error_msg)
                    else:
                        # Secondary database failure - log warning and continue
                        self.logger.warning(f"Secondary database {description} scan failed: {str(e)}")
                        # Create empty scan results for failed secondary database
                        self.database_results[db_name] = ScanResults(database=db_name)
                        self.database_results[db_name].global_error = str(e)
            
            successful_rgi = len(processed_accessions)
            self.logger.info(f"Multi-database analysis completed for {successful_rgi} accessions")

            # Write coordinate manifest for downstream reporting (only needed for CARD)
            try:
                import json as _json
                coord_manifest_file = coords_dir / "coordinate_manifest.json"
                with open(coord_manifest_file, 'w') as _f:
                    _json.dump(coord_manifest, _f, indent=2)
            except Exception as _e:
                self.logger.warning(f"Failed to write coordinate manifest: {_e}")

            return {
                'processed_accessions': processed_accessions,
                'failed_accessions': failed_accessions,
                'database_results': self.database_results,
                'coordinate_cache': coord_manifest,
                'stats': {
                    'total_databases_scanned': len(databases_to_scan),
                    'coordinates_generated': len(coord_manifest),
                }
            }
            
        except Exception as e:
            self.logger.error(f"Multi-database analysis failed: {e}")
            raise
    
    async def _execute_protein_extraction(self, target_genes: List[str]) -> Dict[str, Any]:
        """Execute protein extraction using RGI coordinates"""
        
        try:
            # Initialize production extractor
            extractor = ProductionFastaExtractor(str(self.directories['proteins']))
            
            # Process batch extraction
            success = extractor.process_batch_from_coordinates(
                str(self.directories['genomes']),
                str(self.directories['coordinates']),
                target_genes
            )
            
            if not success:
                raise RuntimeError("Protein extraction failed")
            
            # Update accession tracking with extraction results
            for accession in self.accession_registry:
                if accession in extractor.extraction_results:
                    self.processing_status[accession] = 'proteins_extracted'
                    
                    successful_extractions = len([r for r in extractor.extraction_results[accession] 
                                                if r.extraction_success])
                    self.accession_registry[accession]['proteins_extracted'] = successful_extractions
                    self.accession_registry[accession]['protein_file'] = str(
                        self.directories['proteins'] / "proteins" / f"{accession}_proteins.fasta"
                    )
                elif accession in extractor.failed_extractions:
                    self.processing_status[accession] = 'extraction_failed'
                    self.failure_registry[accession] = extractor.failed_extractions[accession]
            
            self.logger.info(f"Protein extraction completed for {extractor.stats['successful_accessions']} accessions")
            
            return {
                'extraction_results': extractor.extraction_results,
                'failed_extractions': extractor.failed_extractions,
                'stats': extractor.stats
            }
            
        except Exception as e:
            self.logger.error(f"Protein extraction failed: {e}")
            raise
    
    def _generate_master_manifest(self, source_url: str, target_genes: List[str],
                                genome_results: Dict, parallel_results: Dict) -> PipelineManifest:
        """Generate comprehensive master manifest for complete pipeline with parallel processing results"""
        
        try:
            # Calculate success rates
            total_discovered = len(self.accession_registry)
            downloaded_count = len([acc for acc, status in self.processing_status.items() 
                                  if 'downloaded' in status])
            rgi_count = len([acc for acc, status in self.processing_status.items() 
                           if status in ['rgi_analyzed', 'proteins_extracted']])
            extracted_count = len([acc for acc, status in self.processing_status.items() 
                                 if status == 'proteins_extracted'])
            
            download_rate = (downloaded_count / total_discovered) * 100 if total_discovered > 0 else 0
            rgi_rate = (rgi_count / downloaded_count) * 100 if downloaded_count > 0 else 0
            extraction_rate = (extracted_count / rgi_count) * 100 if rgi_count > 0 else 0
            
            # Identify successful accessions (completed all stages)
            successful_accessions = [acc for acc, status in self.processing_status.items() 
                                   if status == 'proteins_extracted']
            
            # Generate pipeline checksum
            pipeline_data = {
                'source_url': source_url,
                'target_genes': target_genes,
                'successful_accessions': sorted(successful_accessions),
                'pipeline_id': self.pipeline_id
            }
            pipeline_checksum = hashlib.sha256(
                json.dumps(pipeline_data, sort_keys=True).encode()
            ).hexdigest()
            
            # Calculate database success rates from parallel processing results
            successful_results = parallel_results.get('successful_results', [])
            total_processed = len(successful_results) + len(parallel_results.get('failed_results', []))
            
            if total_processed > 0:
                # Calculate database-specific success rates from parallel results
                card_successes = sum(1 for r in successful_results if r.abricate_card_success)
                vfdb_successes = sum(1 for r in successful_results if r.abricate_vfdb_success)
                plasmidfinder_successes = sum(1 for r in successful_results if r.abricate_plasmidfinder_success)
                
                card_success_rate = card_successes / total_processed
                vfdb_success_rate = vfdb_successes / total_processed
                plasmidfinder_success_rate = plasmidfinder_successes / total_processed
                
                database_scan_results = {
                    'card': {
                        'total_genomes': total_processed,
                        'successful_scans': card_successes,
                        'success_rate': card_success_rate,
                        'processing_method': 'parallel'
                    },
                    'vfdb': {
                        'total_genomes': total_processed,
                        'successful_scans': vfdb_successes,
                        'success_rate': vfdb_success_rate,
                        'processing_method': 'parallel'
                    },
                    'plasmidfinder': {
                        'total_genomes': total_processed,
                        'successful_scans': plasmidfinder_successes,
                        'success_rate': plasmidfinder_success_rate,
                        'processing_method': 'parallel'
                    }
                }
            else:
                # Fallback if no results available
                card_success_rate = 0.0
                vfdb_success_rate = 0.0
                plasmidfinder_success_rate = 0.0
                database_scan_results = {}
            
            # Create master manifest
            manifest = PipelineManifest(
                pipeline_version="2.0",
                generated_timestamp=time.strftime('%Y-%m-%d %H:%M:%S'),
                pipeline_id=self.pipeline_id,
                source_url=source_url,
                total_genomes_discovered=total_discovered,
                total_genomes_downloaded=downloaded_count,
                total_genomes_rgi_analyzed=rgi_count,
                total_proteins_extracted=parallel_results.get('total_proteins', 0),
                database_scan_results=database_scan_results,
                card_success_rate=card_success_rate,
                vfdb_success_rate=vfdb_success_rate,
                plasmidfinder_success_rate=plasmidfinder_success_rate,
                download_success_rate=download_rate,
                rgi_success_rate=rgi_rate,
                extraction_success_rate=extraction_rate,
                genome_manifest_path=str(self.directories['genomes'] / "genome_manifest.json"),
                coordinate_manifest_path=str(self.directories['coordinates'] / "coordinate_manifest.json"),
                extraction_manifest_path=str(self.directories['proteins'] / "manifests" / "extraction_manifest.json"),
                pipeline_checksum=pipeline_checksum,
                manifest_checksum="",  # Will be calculated after serialization
                successful_accessions=successful_accessions,
                failed_accessions=self.failure_registry,
                target_genes=target_genes
            )
            
            # Save master manifest
            manifest_file = self.directories['manifests'] / f"{self.pipeline_id}_master_manifest.json"
            manifest_dict = asdict(manifest)
            
            # Calculate manifest checksum
            manifest_checksum = hashlib.sha256(
                json.dumps(manifest_dict, sort_keys=True).encode()
            ).hexdigest()
            manifest.manifest_checksum = manifest_checksum
            manifest_dict['manifest_checksum'] = manifest_checksum
            
            with open(manifest_file, 'w') as f:
                json.dump(manifest_dict, f, indent=2)
            
            # Generate detailed accession report
            self._generate_accession_report()
            
            # Generate performance report
            self._generate_performance_report()
            
            self.logger.info(f"Master manifest generated: {manifest_file}")
            return manifest
            
        except Exception as e:
            self.logger.error(f"Failed to generate master manifest: {e}")
            raise
    
    def _generate_accession_report(self):
        """Generate detailed accession-by-accession report"""
        
        report_file = self.directories['reports'] / f"{self.pipeline_id}_accession_report.json"
        
        accession_report = {
            "report_version": "1.0",
            "generated_timestamp": time.strftime('%Y-%m-%d %H:%M:%S'),
            "pipeline_id": self.pipeline_id,
            "total_accessions": len(self.accession_registry),
            "accessions": {}
        }
        
        for accession, data in self.accession_registry.items():
            accession_report["accessions"][accession] = {
                **data,
                "final_status": self.processing_status.get(accession, 'unknown'),
                "failure_reason": self.failure_registry.get(accession, None)
            }
        
        with open(report_file, 'w') as f:
            json.dump(accession_report, f, indent=2)
        
        self.logger.info(f"Accession report generated: {report_file}")
    
    def _generate_performance_report(self):
        """Generate performance analysis report"""
        
        report_file = self.directories['reports'] / f"{self.pipeline_id}_performance_report.json"
        
        performance_report = {
            "report_version": "1.0",
            "generated_timestamp": time.strftime('%Y-%m-%d %H:%M:%S'),
            "pipeline_id": self.pipeline_id,
            "total_pipeline_time_seconds": self.performance_metrics['total_pipeline_time'],
            "stage_times": self.performance_metrics['stage_times'],
            "processing_rates": {
                "genomes_per_minute": (len(self.accession_registry) / 
                                     (self.performance_metrics['total_pipeline_time'] / 60)),
                "successful_accessions_per_minute": (len([acc for acc, status in self.processing_status.items() 
                                                         if status == 'proteins_extracted']) / 
                                                    (self.performance_metrics['total_pipeline_time'] / 60))
            }
        }
        
        with open(report_file, 'w') as f:
            json.dump(performance_report, f, indent=2)
        
        self.logger.info(f"Performance report generated: {report_file}")


async def main():
    """Command line interface for production pipeline orchestrator"""
    import argparse
    
    parser = argparse.ArgumentParser(description="Production Pipeline Orchestrator")
    parser.add_argument("--source-url", required=True, help="NCBI search URL for genome discovery")
    parser.add_argument("--target-genes", nargs="+", required=True, help="Target gene names")
    parser.add_argument("--output-dir", required=True, help="Base output directory")
    parser.add_argument("--config", help="YAML configuration file")
    parser.add_argument("--email", required=True, help="Email for NCBI API")
    parser.add_argument("--api-key", help="NCBI API key (optional but recommended)")
    parser.add_argument("--max-concurrent", type=int, default=10, help="Max concurrent downloads")
    parser.add_argument("--rgi-threads", type=int, default=4, help="RGI analysis threads")
    parser.add_argument("--threads", "-t", type=int, default=os.cpu_count() or 4, 
                       help=f"Number of parallel processing threads (default: {os.cpu_count() or 4})")
    parser.add_argument("--force", action="store_true", help="Force execution even with warnings")
    
    args = parser.parse_args()
    
    # Validate threads parameter
    max_cpu_cores = os.cpu_count() or 4  # Fallback to 4 if cpu_count returns None
    
    if args.threads < 1:
        print(f"ERROR: --threads must be at least 1, got {args.threads}")
        return 1
    elif args.threads > max_cpu_cores * 2:
        print(f"WARNING: --threads ({args.threads}) is higher than 2x CPU cores ({max_cpu_cores * 2})")
        print("This may cause performance degradation due to context switching overhead")
        
        if args.force:
            print("Proceeding due to --force flag")
        elif sys.stdin.isatty() and hasattr(sys.stdin, 'readable') and sys.stdin.readable():
            # Interactive terminal with readable stdin - prompt for confirmation
            try:
                response = input("Continue anyway? (y/N): ")
                if response.lower() != 'y':
                    print("Aborted by user")
                    return 1
            except (EOFError, OSError):
                # Input not available - treat as non-interactive
                print("ERROR: Running in non-interactive mode. Use --force to override this warning.")
                return 1
        else:
            # Non-interactive environment - abort without blocking
            print("ERROR: Running in non-interactive mode. Use --force to override this warning.")
            return 1
    
    print(f"Using {args.threads} parallel processing threads ({max_cpu_cores} CPU cores detected)")
    
    # Setup configuration
    config = {
        'email': args.email,
        'api_key': args.api_key,
        'max_concurrent_downloads': args.max_concurrent,
        'rgi_threads': args.rgi_threads,
        'threads': args.threads
    }
    
    if args.config and Path(args.config).exists():
        with open(args.config, 'r') as f:
            file_config = yaml.safe_load(f)
            config.update(file_config)
    
    # Initialize orchestrator
    orchestrator = ProductionPipelineOrchestrator(args.output_dir, config)
    
    try:
        # Execute complete pipeline
        manifest = await orchestrator.execute_complete_pipeline(args.source_url, args.target_genes)
        
        print(f"\nPipeline completed successfully!")
        print(f"Pipeline ID: {manifest.pipeline_id}")
        print(f"Successful accessions: {len(manifest.successful_accessions)}")
        print(f"Total proteins extracted: {manifest.total_proteins_extracted}")
        print(f"Master manifest: {orchestrator.directories['manifests'] / f'{manifest.pipeline_id}_master_manifest.json'}")
        
        return 0
        
    except Exception as e:
        print(f"Pipeline failed: {e}")
        return 1


if __name__ == "__main__":
    import asyncio
    sys.exit(asyncio.run(main()))