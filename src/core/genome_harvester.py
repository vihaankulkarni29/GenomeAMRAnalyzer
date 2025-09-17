"""
Genome Harvester
---------------
Professional-grade genome downloader with robust error handling, validation, and logging.
Supports parallel downloads, retries, and comprehensive status tracking.
"""

import os
import re
import logging
import asyncio
import aiohttp
import aiofiles
from typing import List, Dict, Set, Optional, Tuple
from pathlib import Path
from dataclasses import dataclass
from concurrent.futures import ThreadPoolExecutor
import hashlib
import time


@dataclass
class AccessionStatus:
    """Track status of each accession throughout the process"""
    accession: str
    is_valid: bool = False
    is_duplicate: bool = False
    download_attempted: bool = False
    download_success: bool = False
    file_path: Optional[str] = None
    error_message: Optional[str] = None
    file_size_bytes: int = 0
    download_time_seconds: float = 0.0


class GenomeHarvester:
    """
    Professional genome downloader with:
    - Input validation and deduplication
    - Parallel downloads with rate limiting
    - Comprehensive error handling and logging
    - Quality control and file validation
    """
    
    def __init__(self, 
                 output_dir: str,
                 email: str,
                 api_key: Optional[str] = None,
                 max_concurrent: int = 5,
                 retry_attempts: int = 3,
                 logger: Optional[logging.Logger] = None):
        
        self.output_dir = Path(output_dir)
        self.email = email
        self.api_key = api_key
        self.max_concurrent = max_concurrent
        self.retry_attempts = retry_attempts
        self.logger = logger or self._setup_logger()
        
        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Status tracking
        self.status_log: Dict[str, AccessionStatus] = {}
        
        # Validation patterns
        self.valid_patterns = [
            r'^GCF_\d+\.\d+$',  # RefSeq assemblies
            r'^GCA_\d+\.\d+$',  # GenBank assemblies
            r'^NZ_\w+\.\d+$',   # RefSeq contigs/chromosomes
            r'^NC_\d+\.\d+$',   # RefSeq chromosomes
            r'^CP\d+\.\d+$',    # Complete genomes
        ]
        
    def _setup_logger(self) -> logging.Logger:
        """Setup comprehensive logging"""
        logger = logging.getLogger("GenomeHarvester")
        if not logger.handlers:
            # Console handler only for testing
            console_handler = logging.StreamHandler()
            console_formatter = logging.Formatter(
                '%(levelname)s %(name)s: %(message)s'
            )
            console_handler.setFormatter(console_formatter)
            logger.addHandler(console_handler)
            
        logger.setLevel(logging.INFO)
        return logger
    
    def validate_accession_list(self, accession_file: str) -> List[str]:
        """
        Validate, clean, and deduplicate accession list
        Returns: List of valid, unique accessions
        """
        self.logger.info(f"Validating accession list: {accession_file}")
        
        if not os.path.exists(accession_file):
            raise FileNotFoundError(f"Accession file not found: {accession_file}")
        
        valid_accessions = []
        seen_accessions: Set[str] = set()
        
        with open(accession_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                accession = line.strip()
                
                # Skip empty lines
                if not accession:
                    continue
                
                # Initialize status tracking
                status = AccessionStatus(accession=accession)
                
                # Validate format
                is_valid = any(re.match(pattern, accession) for pattern in self.valid_patterns)
                status.is_valid = is_valid
                
                if not is_valid:
                    status.error_message = f"Invalid accession format (line {line_num})"
                    self.logger.warning(f"Invalid accession format: {accession} (line {line_num})")
                    self.status_log[accession] = status
                    continue
                
                # Check for duplicates
                if accession in seen_accessions:
                    status.is_duplicate = True
                    status.error_message = f"Duplicate accession (line {line_num})"
                    self.logger.warning(f"Duplicate accession: {accession} (line {line_num})")
                    self.status_log[accession] = status
                    continue
                
                # Valid and unique
                seen_accessions.add(accession)
                valid_accessions.append(accession)
                self.status_log[accession] = status
        
        self.logger.info(f"Validation complete: {len(valid_accessions)} valid, "
                        f"{len(self.status_log) - len(valid_accessions)} invalid/duplicate")
        
        return valid_accessions
    
    def _construct_download_url(self, accession: str) -> str:
        """Construct NCBI FTP download URL for accession"""
        if accession.startswith(('GCF_', 'GCA_')):
            # Assembly accessions
            prefix = accession[:3]
            assembly_id = accession.replace('.', '_')
            return f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{prefix}/{accession[:7]}/{accession[:13]}/{assembly_id}/{assembly_id}_genomic.fna.gz"
        elif accession.startswith(('NZ_', 'NC_', 'CP')):
            # Individual sequences via NCBI E-utilities
            return f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession}&rettype=fasta&retmode=text"
        else:
            raise ValueError(f"Unsupported accession format: {accession}")
    
    async def _download_single_genome(self, 
                                    session: aiohttp.ClientSession, 
                                    accession: str) -> bool:
        """Download a single genome with retries and error handling"""
        
        status = self.status_log[accession]
        status.download_attempted = True
        
        output_file = self.output_dir / f"{accession}.fasta"
        
        # Skip if already exists and valid
        if output_file.exists() and output_file.stat().st_size > 0:
            status.download_success = True
            status.file_path = str(output_file)
            status.file_size_bytes = output_file.stat().st_size
            self.logger.info(f"Genome already exists: {accession}")
            return True
        
        url = self._construct_download_url(accession)
        start_time = time.time()
        
        for attempt in range(self.retry_attempts):
            try:
                self.logger.info(f"Downloading {accession} (attempt {attempt + 1}/{self.retry_attempts})")
                
                params = {'email': self.email}
                if self.api_key:
                    params['api_key'] = self.api_key
                
                async with session.get(url, params=params) as response:
                    if response.status == 200:
                        content = await response.read()
                        
                        # Write to file
                        async with aiofiles.open(output_file, 'wb') as f:
                            await f.write(content)
                        
                        # Validate downloaded file
                        if output_file.stat().st_size > 0:
                            status.download_success = True
                            status.file_path = str(output_file)
                            status.file_size_bytes = output_file.stat().st_size
                            status.download_time_seconds = time.time() - start_time
                            
                            self.logger.info(f"Successfully downloaded {accession} "
                                           f"({status.file_size_bytes} bytes)")
                            return True
                        else:
                            raise Exception("Downloaded file is empty")
                    
                    else:
                        raise Exception(f"HTTP {response.status}: {response.reason}")
            
            except Exception as e:
                self.logger.warning(f"Download attempt {attempt + 1} failed for {accession}: {e}")
                if attempt == self.retry_attempts - 1:
                    status.error_message = f"Download failed after {self.retry_attempts} attempts: {e}"
                    self.logger.error(f"Failed to download {accession}: {e}")
                    return False
                
                # Wait before retry
                await asyncio.sleep(2 ** attempt)
        
        return False
    
    async def download_genomes(self, accessions: List[str]) -> Dict[str, AccessionStatus]:
        """
        Download genomes with parallel processing and rate limiting
        Returns: Dictionary of accession -> status
        """
        self.logger.info(f"Starting download of {len(accessions)} genomes")
        
        # Rate limiting connector
        connector = aiohttp.TCPConnector(limit=self.max_concurrent)
        timeout = aiohttp.ClientTimeout(total=300)  # 5 minute timeout
        
        async with aiohttp.ClientSession(connector=connector, timeout=timeout) as session:
            semaphore = asyncio.Semaphore(self.max_concurrent)
            
            async def download_with_semaphore(accession: str) -> bool:
                async with semaphore:
                    return await self._download_single_genome(session, accession)
            
            # Execute downloads
            tasks = [download_with_semaphore(acc) for acc in accessions]
            results = await asyncio.gather(*tasks, return_exceptions=True)
        
        # Summary statistics
        successful = sum(1 for status in self.status_log.values() if status.download_success)
        failed = len(accessions) - successful
        total_size = sum(status.file_size_bytes for status in self.status_log.values() 
                        if status.download_success)
        
        self.logger.info(f"Download complete: {successful} successful, {failed} failed, "
                        f"{total_size / 1024 / 1024:.2f} MB total")
        
        return self.status_log
    
    def generate_status_report(self, output_file: Optional[str] = None) -> str:
        """Generate detailed status report"""
        
        if output_file is None:
            output_file = self.output_dir / "download_status_report.txt"
        
        successful = [s for s in self.status_log.values() if s.download_success]
        failed = [s for s in self.status_log.values() if not s.download_success and s.download_attempted]
        invalid = [s for s in self.status_log.values() if not s.is_valid]
        duplicates = [s for s in self.status_log.values() if s.is_duplicate]
        
        report = f"""
Genome Harvester Status Report
==============================
Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}

Summary:
--------
Total accessions processed: {len(self.status_log)}
Valid accessions: {len(self.status_log) - len(invalid)}
Invalid accessions: {len(invalid)}
Duplicate accessions: {len(duplicates)}
Successful downloads: {len(successful)}
Failed downloads: {len(failed)}
Total data downloaded: {sum(s.file_size_bytes for s in successful) / 1024 / 1024:.2f} MB

Successful Downloads:
--------------------
"""
        
        for status in successful:
            report += f"{status.accession}\t{status.file_size_bytes / 1024:.1f} KB\t{status.download_time_seconds:.1f}s\n"
        
        if failed:
            report += "\nFailed Downloads:\n-----------------\n"
            for status in failed:
                report += f"{status.accession}\t{status.error_message}\n"
        
        if invalid:
            report += "\nInvalid Accessions:\n-------------------\n"
            for status in invalid:
                report += f"{status.accession}\t{status.error_message}\n"
        
        with open(output_file, 'w') as f:
            f.write(report)
        
        self.logger.info(f"Status report saved to: {output_file}")
        return report


def main():
    """Example usage and testing"""
    
    # Configuration
    config = {
        'output_dir': 'genomes',
        'email': 'vihaankulkarni29@gmail.com',
        'api_key': None,
        'max_concurrent': 3,
        'retry_attempts': 3
    }
    
    # Initialize harvester
    harvester = GenomeHarvester(**config)
    
    # Process accession list
    try:
        valid_accessions = harvester.validate_accession_list('accession_list.txt')
        
        if valid_accessions:
            # Download genomes
            asyncio.run(harvester.download_genomes(valid_accessions))
            
            # Generate report
            harvester.generate_status_report()
        else:
            harvester.logger.error("No valid accessions found")
            
    except Exception as e:
        harvester.logger.error(f"Pipeline failed: {e}")
        raise


if __name__ == "__main__":
    main()