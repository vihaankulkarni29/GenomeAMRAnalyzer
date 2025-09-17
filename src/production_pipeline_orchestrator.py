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
from pathlib import Path
from typing import Dict, List, Optional, Set, Any
from dataclasses import dataclass, asdict
import yaml

# Import our production components
sys.path.append(str(Path(__file__).parent))
from core.url_to_genomes_workflow import URLToGenomesWorkflow
from core.genome_downloader import GenomeDownloadEngine
from card_runner import CARDRunner, RGIConfig
from production_fasta_extractor import ProductionFastaExtractor


@dataclass
class PipelineManifest:
    """Complete pipeline manifest with full traceability"""
    pipeline_version: str
    generated_timestamp: str
    pipeline_id: str
    source_url: str
    
    # Batch information
    total_genomes_discovered: int
    total_genomes_downloaded: int
    total_genomes_rgi_analyzed: int
    total_proteins_extracted: int
    
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
        
        # Create structured output directories
        self.directories = {
            'genomes': self.output_base_dir / "genomes",
            'coordinates': self.output_base_dir / "coordinates", 
            'proteins': self.output_base_dir / "proteins",
            'manifests': self.output_base_dir / "manifests",
            'logs': self.output_base_dir / "logs",
            'reports': self.output_base_dir / "reports"
        }
        
        for dir_path in self.directories.values():
            dir_path.mkdir(parents=True, exist_ok=True)
        
        # Setup logging
        self._setup_logging()
        
        # Initialize tracking
        self.accession_registry: Dict[str, Dict[str, Any]] = {}
        self.processing_status: Dict[str, str] = {}  # accession -> status
        self.failure_registry: Dict[str, str] = {}  # accession -> failure reason
        
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
            
            # Stage 2: CARD RGI Analysis
            stage_start = time.time()
            self.logger.info("=== STAGE 2: CARD RGI Analysis ===")
            
            rgi_results = await self._execute_rgi_analysis(target_genes)
            
            self.performance_metrics['stage_times']['rgi_analysis'] = time.time() - stage_start
            
            # Stage 3: Protein Extraction
            stage_start = time.time()
            self.logger.info("=== STAGE 3: Protein Extraction ===")
            
            extraction_results = await self._execute_protein_extraction(target_genes)
            
            self.performance_metrics['stage_times']['protein_extraction'] = time.time() - stage_start
            
            # Stage 4: Generate Master Manifest
            stage_start = time.time()
            self.logger.info("=== STAGE 4: Master Manifest Generation ===")
            
            master_manifest = self._generate_master_manifest(
                source_url, target_genes, genome_results, rgi_results, extraction_results
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
        """Execute genome discovery and download with accession tracking"""
        
        try:
            # Initialize URL workflow
            workflow = URLToGenomesWorkflow(
                output_dir=str(self.directories['genomes']),
                email=self.config.get('email', 'pipeline@example.com'),
                api_key=self.config.get('api_key'),
                max_concurrent=self.config.get('max_concurrent_downloads', 10)
            )
            
            # Execute workflow
            results = await workflow.run_complete_workflow(source_url)
            
            # Register accessions from download results
            for accession, result in results.items():
                self.accession_registry[accession] = {
                    'genome_file': result.file_path if result.success else None,
                    'organism': getattr(result, 'organism', ''),
                    'genome_length': getattr(result, 'genome_length', 0),
                    'download_timestamp': getattr(result, 'download_timestamp', ''),
                    'file_checksum': getattr(result, 'file_checksum', ''),
                    'download_success': result.success
                }
                
                if result.success:
                    self.processing_status[accession] = 'genome_downloaded'
                else:
                    self.processing_status[accession] = 'genome_download_failed'
                    self.failure_registry[accession] = result.error_message or 'Unknown download error'
            
            successful_downloads = len([r for r in results.values() if r.success])
            self.logger.info(f"Downloaded {successful_downloads}/{len(results)} genomes successfully")
            
            return results
            
        except Exception as e:
            self.logger.error(f"Genome discovery/download failed: {e}")
            raise
    
    async def _execute_rgi_analysis(self, target_genes: List[str]) -> Dict[str, Any]:
        """Execute CARD RGI analysis on downloaded genomes"""
        
        try:
            # Setup RGI configuration
            rgi_config = RGIConfig(
                input_dir=str(self.directories['genomes']),
                output_dir=str(self.directories['coordinates']),
                target_genes=target_genes,
                num_threads=self.config.get('rgi_threads', 4)
            )
            
            # Initialize CARD runner
            card_runner = CARDRunner(rgi_config)
            
            # Execute RGI analysis
            success = card_runner.run()
            
            if not success:
                raise RuntimeError("CARD RGI analysis failed")
            
            # Update accession tracking with RGI results
            for accession in self.accession_registry:
                if accession in card_runner.processed_accessions:
                    self.processing_status[accession] = 'rgi_analyzed'
                    
                    # Get coordinate count
                    if accession in card_runner.coordinate_cache:
                        coord_count = len(card_runner.coordinate_cache[accession])
                        self.accession_registry[accession]['rgi_coordinates_found'] = coord_count
                        self.accession_registry[accession]['coordinate_file'] = str(
                            card_runner._create_accession_output_path(accession, "coordinates")
                        )
                elif accession in card_runner.failed_accessions:
                    self.processing_status[accession] = 'rgi_failed'
                    self.failure_registry[accession] = 'RGI analysis failed'
            
            successful_rgi = len(card_runner.processed_accessions)
            self.logger.info(f"RGI analysis completed for {successful_rgi} accessions")
            
            return {
                'processed_accessions': card_runner.processed_accessions,
                'failed_accessions': card_runner.failed_accessions,
                'coordinate_cache': card_runner.coordinate_cache,
                'stats': card_runner.stats
            }
            
        except Exception as e:
            self.logger.error(f"RGI analysis failed: {e}")
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
                                genome_results: Dict, rgi_results: Dict, 
                                extraction_results: Dict) -> PipelineManifest:
        """Generate comprehensive master manifest for complete pipeline"""
        
        try:
            # Calculate success rates
            total_discovered = len(self.accession_registry)
            downloaded_count = len([acc for acc, status in self.processing_status.items() 
                                  if 'downloaded' in status])
            rgi_count = len([acc for acc, status in self.processing_status.items() 
                           if status == 'rgi_analyzed'])
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
            
            # Create master manifest
            manifest = PipelineManifest(
                pipeline_version="2.0",
                generated_timestamp=time.strftime('%Y-%m-%d %H:%M:%S'),
                pipeline_id=self.pipeline_id,
                source_url=source_url,
                total_genomes_discovered=total_discovered,
                total_genomes_downloaded=downloaded_count,
                total_genomes_rgi_analyzed=rgi_count,
                total_proteins_extracted=extraction_results['stats']['total_proteins_extracted'],
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
    
    args = parser.parse_args()
    
    # Setup configuration
    config = {
        'email': args.email,
        'api_key': args.api_key,
        'max_concurrent_downloads': args.max_concurrent,
        'rgi_threads': args.rgi_threads
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