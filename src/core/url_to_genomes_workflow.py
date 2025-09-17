"""
Complete URL-to-Genomes Workflow
===============================
Integrated workflow that takes an NCBI URL and produces downloaded FASTA files
ready for RGI processing.
"""

import asyncio
import logging
import time
from typing import List, Dict, Optional, Tuple
from pathlib import Path
import yaml

# Import our custom modules (adjust import paths as needed)
import sys
import os
sys.path.append(os.path.dirname(__file__))

from ncbi_genome_discovery import URLBasedGenomeDiscovery, GenomeMetadata
from genome_downloader import GenomeDownloader, DownloadResult


class URLToGenomesWorkflow:
    """
    Complete workflow: NCBI URL → Genome Discovery → Download → Ready for RGI
    """
    
    def __init__(self, config_file: str):
        """Initialize workflow with configuration"""
        self.config = self._load_config(config_file)
        self.logger = self._setup_logger()
        
        # Initialize components
        self.discovery_engine = URLBasedGenomeDiscovery(
            email=self.config['ncbi_email'],
            api_key=self.config.get('ncbi_api_key'),
            max_results=self.config['ncbi_search']['max_genomes'],
            logger=self.logger
        )
        
        self.downloader = GenomeDownloader(
            output_dir=self.config['directories']['genomes'],
            email=self.config['ncbi_email'],
            api_key=self.config.get('ncbi_api_key'),
            max_concurrent=3,
            retry_attempts=self.config['ncbi_search']['retry_attempts'],
            timeout_seconds=self.config['ncbi_search']['timeout_seconds'],
            logger=self.logger
        )
        
        # Results tracking
        self.discovered_genomes: List[GenomeMetadata] = []
        self.download_results: Dict[str, DownloadResult] = {}
        self.successful_files: List[str] = []
    
    def _load_config(self, config_file: str) -> dict:
        """Load configuration from YAML file"""
        try:
            with open(config_file, 'r') as f:
                return yaml.safe_load(f)
        except Exception as e:
            raise ValueError(f"Failed to load config file {config_file}: {e}")
    
    def _setup_logger(self) -> logging.Logger:
        """Setup comprehensive logging for the workflow"""
        logger = logging.getLogger("URLToGenomesWorkflow")
        if not logger.handlers:
            # Console handler
            console_handler = logging.StreamHandler()
            console_formatter = logging.Formatter('%(levelname)s %(name)s: %(message)s')
            console_handler.setFormatter(console_formatter)
            logger.addHandler(console_handler)
            
            # File handler
            log_dir = Path(self.config['directories']['logs'])
            log_dir.mkdir(parents=True, exist_ok=True)
            log_file = log_dir / "url_to_genomes_workflow.log"
            
            file_handler = logging.FileHandler(log_file)
            file_formatter = logging.Formatter('%(asctime)s %(levelname)s %(name)s: %(message)s')
            file_handler.setFormatter(file_formatter)
            logger.addHandler(file_handler)
        
        logger.setLevel(logging.INFO)
        return logger
    
    async def run_complete_workflow(self, ncbi_url: str) -> Tuple[List[str], str]:
        """
        Run the complete workflow from URL to downloaded genomes
        
        Args:
            ncbi_url: NCBI search URL
            
        Returns:
            Tuple of (list of downloaded file paths, workflow report)
        """
        self.logger.info("=== Starting Complete URL-to-Genomes Workflow ===")
        self.logger.info(f"Input URL: {ncbi_url}")
        
        workflow_start_time = time.time()
        
        try:
            # Phase 1: Discover genomes from URL
            self.logger.info("Phase 1: Discovering genomes from URL...")
            self.discovered_genomes, discovery_report = await self.discovery_engine.discover_from_url(ncbi_url)
            
            if not self.discovered_genomes:
                raise Exception("No genomes discovered from URL")
            
            self.logger.info(f"Discovered {len(self.discovered_genomes)} genomes")
            
            # Phase 2: Download discovered genomes
            self.logger.info("Phase 2: Downloading discovered genomes...")
            self.download_results = await self.downloader.download_genomes(self.discovered_genomes)
            
            # Phase 3: Collect successful downloads
            self.successful_files = self.downloader.get_downloaded_files()
            
            if not self.successful_files:
                raise Exception("No genomes were successfully downloaded")
            
            self.logger.info(f"Successfully downloaded {len(self.successful_files)} genomes")
            
            # Phase 4: Generate comprehensive report
            workflow_time = time.time() - workflow_start_time
            comprehensive_report = self._generate_comprehensive_report(
                discovery_report, 
                ncbi_url, 
                workflow_time
            )
            
            self.logger.info("=== Workflow Complete ===")
            self.logger.info(f"Ready for RGI processing: {len(self.successful_files)} FASTA files")
            
            return self.successful_files, comprehensive_report
            
        except Exception as e:
            self.logger.error(f"Workflow failed: {e}")
            raise
    
    def _generate_comprehensive_report(self, discovery_report: str, ncbi_url: str, workflow_time: float) -> str:
        """Generate comprehensive workflow report"""
        
        successful_downloads = [r for r in self.download_results.values() if r.success]
        failed_downloads = [r for r in self.download_results.values() if not r.success]
        
        report = f"""
URL-to-Genomes Workflow Report
=============================
Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}
Workflow Duration: {workflow_time:.1f} seconds

Input:
------
NCBI URL: {ncbi_url}
Max Genomes: {self.config['ncbi_search']['max_genomes']}
Target Genes: {', '.join(self.config['target_genes'])}

Results Summary:
---------------
Genomes Discovered: {len(self.discovered_genomes)}
Successful Downloads: {len(successful_downloads)}
Failed Downloads: {len(failed_downloads)}
Success Rate: {len(successful_downloads) / len(self.discovered_genomes) * 100:.1f}%
Total Data Downloaded: {sum(r.file_size_bytes for r in successful_downloads) / 1024 / 1024:.2f} MB

Ready for RGI Processing:
------------------------
"""
        
        for file_path in self.successful_files:
            file_size = Path(file_path).stat().st_size / 1024
            accession = Path(file_path).stem
            report += f"{accession}\t{file_path}\t{file_size:.1f} KB\n"
        
        if failed_downloads:
            report += "\nFailed Downloads:\n-----------------\n"
            for result in failed_downloads:
                report += f"{result.accession}\t{result.error_message}\n"
        
        report += f"\n{discovery_report}"
        
        # Save report to file
        report_file = Path(self.config['directories']['reports']) / "url_to_genomes_report.txt"
        report_file.parent.mkdir(parents=True, exist_ok=True)
        with open(report_file, 'w') as f:
            f.write(report)
        
        self.logger.info(f"Comprehensive report saved to: {report_file}")
        return report
    
    def get_genomes_for_rgi(self) -> List[Dict[str, str]]:
        """
        Get genome information formatted for RGI processing
        
        Returns:
            List of dictionaries with genome info for RGI
        """
        rgi_genomes = []
        
        for file_path in self.successful_files:
            accession = Path(file_path).stem
            
            # Find corresponding metadata
            genome_metadata = None
            for genome in self.discovered_genomes:
                if genome.accession == accession:
                    genome_metadata = genome
                    break
            
            rgi_genome = {
                'accession': accession,
                'file_path': file_path,
                'organism': genome_metadata.organism if genome_metadata else "Unknown",
                'length': genome_metadata.length if genome_metadata else 0,
                'ready_for_rgi': True
            }
            
            rgi_genomes.append(rgi_genome)
        
        return rgi_genomes
    
    def validate_for_rgi_processing(self) -> bool:
        """
        Validate that downloaded genomes are ready for RGI processing
        
        Returns:
            True if all files are valid for RGI
        """
        if not self.successful_files:
            self.logger.error("No successful downloads to validate")
            return False
        
        valid_count = 0
        
        for file_path in self.successful_files:
            try:
                path = Path(file_path)
                
                # Check file exists and has content
                if not path.exists() or path.stat().st_size == 0:
                    self.logger.error(f"File missing or empty: {file_path}")
                    continue
                
                # Check FASTA format
                with open(path, 'r') as f:
                    first_line = f.readline()
                    if not first_line.startswith('>'):
                        self.logger.error(f"Invalid FASTA format: {file_path}")
                        continue
                
                # Check minimum size (RGI needs substantial sequences)
                if path.stat().st_size < 100000:  # Less than 100KB
                    self.logger.warning(f"File may be too small for RGI: {file_path}")
                
                valid_count += 1
                self.logger.debug(f"Validated for RGI: {file_path}")
                
            except Exception as e:
                self.logger.error(f"Validation failed for {file_path}: {e}")
        
        success_rate = valid_count / len(self.successful_files)
        self.logger.info(f"RGI validation: {valid_count}/{len(self.successful_files)} files valid ({success_rate*100:.1f}%)")
        
        return success_rate >= 0.8  # At least 80% should be valid


async def main():
    """Example usage of the complete workflow"""
    
    # Test URL (E. coli with macrolide resistance)
    test_url = "https://www.ncbi.nlm.nih.gov/nuccore/?term=(Escherichia+coli+and+macrolide+resistance)+AND+%22Escherichia+coli%22%5Bporgn%3A__txid562%5D+and+complete+genome"
    
    try:
        # Initialize workflow
        workflow = URLToGenomesWorkflow("config/snakemake_config.yaml")
        
        # Run complete workflow
        genome_files, report = await workflow.run_complete_workflow(test_url)
        
        print(f"\n🎉 Workflow completed successfully!")
        print(f"Downloaded {len(genome_files)} genomes ready for RGI processing")
        
        # Validate for RGI
        if workflow.validate_for_rgi_processing():
            print("✅ All genomes validated and ready for RGI processing")
            
            # Show RGI-ready genomes
            rgi_genomes = workflow.get_genomes_for_rgi()
            print(f"\nGenomes ready for RGI:")
            for i, genome in enumerate(rgi_genomes[:5], 1):
                print(f"{i}. {genome['accession']} - {genome['organism'][:50]} ({genome['length']/1_000_000:.2f} MB)")
        else:
            print("⚠️  Some genomes may not be suitable for RGI processing")
        
    except Exception as e:
        print(f"❌ Workflow failed: {e}")


if __name__ == "__main__":
    asyncio.run(main())