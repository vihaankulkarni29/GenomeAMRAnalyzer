"""
URL-Based Genome Harvester
--------------------------
Professional genome harvester that accepts NCBI search URLs and extracts accessions
from the search results. Much more flexible than static accession lists.
"""

import os
import re
import logging
import asyncio
import aiohttp
import aiofiles
import xml.etree.ElementTree as ET
from typing import List, Dict, Set, Optional, Tuple
from pathlib import Path
from dataclasses import dataclass
from urllib.parse import parse_qs, urlparse
import time
from bs4 import BeautifulSoup


@dataclass
class GenomeResult:
    """Track each genome result from NCBI query"""
    accession: str
    title: str = ""
    organism: str = ""
    length: int = 0
    download_attempted: bool = False
    download_success: bool = False
    file_path: Optional[str] = None
    error_message: Optional[str] = None
    file_size_bytes: int = 0
    download_time_seconds: float = 0.0


class URLGenomeHarvester:
    """
    Professional URL-based genome harvester with:
    - NCBI search URL parsing and accession extraction
    - Configurable result limits (e.g., first 10 genomes)
    - Parallel downloads with rate limiting
    - Comprehensive error handling and logging
    - Quality control and file validation
    """
    
    def __init__(self, 
                 output_dir: str,
                 email: str,
                 api_key: Optional[str] = None,
                 max_results: int = 10,
                 max_concurrent: int = 3,
                 retry_attempts: int = 3,
                 logger: Optional[logging.Logger] = None):
        
        self.output_dir = Path(output_dir)
        self.email = email
        self.api_key = api_key
        self.max_results = max_results
        self.max_concurrent = max_concurrent
        self.retry_attempts = retry_attempts
        self.logger = logger or self._setup_logger()
        
        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Results tracking
        self.genome_results: Dict[str, GenomeResult] = {}
        
        # NCBI E-utilities base URLs
        self.esearch_base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        self.efetch_base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        self.esummary_base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        
    def _setup_logger(self) -> logging.Logger:
        """Setup comprehensive logging"""
        logger = logging.getLogger("URLGenomeHarvester")
        if not logger.handlers:
            console_handler = logging.StreamHandler()
            console_formatter = logging.Formatter(
                '%(levelname)s %(name)s: %(message)s'
            )
            console_handler.setFormatter(console_formatter)
            logger.addHandler(console_handler)
            
        logger.setLevel(logging.INFO)
        return logger
    
    def parse_ncbi_url(self, ncbi_url: str) -> str:
        """
        Parse NCBI search URL and extract the search query
        """
        self.logger.info(f"Parsing NCBI URL: {ncbi_url}")
        
        try:
            parsed = urlparse(ncbi_url)
            
            # Handle different NCBI URL formats
            if 'nuccore' in parsed.path:
                # Extract the term parameter
                if '?term=' in ncbi_url:
                    # Direct term parameter
                    query_params = parse_qs(parsed.query)
                    if 'term' in query_params:
                        search_term = query_params['term'][0]
                        self.logger.info(f"Extracted search term: {search_term}")
                        return search_term
                
                # If no direct term, try to extract from URL fragments
                if 'term=' in ncbi_url:
                    term_start = ncbi_url.find('term=') + 5
                    term_end = ncbi_url.find('&', term_start)
                    if term_end == -1:
                        term_end = len(ncbi_url)
                    search_term = ncbi_url[term_start:term_end]
                    
                    # URL decode
                    import urllib.parse
                    search_term = urllib.parse.unquote_plus(search_term)
                    self.logger.info(f"Extracted and decoded search term: {search_term}")
                    return search_term
            
            raise ValueError("Could not extract search term from URL")
            
        except Exception as e:
            self.logger.error(f"Failed to parse NCBI URL: {e}")
            raise ValueError(f"Invalid NCBI URL format: {e}")
    
    async def search_ncbi(self, search_term: str) -> List[str]:
        """
        Search NCBI using E-utilities and return list of accession IDs
        """
        self.logger.info(f"Searching NCBI for: {search_term}")
        
        # Prepare search parameters
        params = {
            'db': 'nuccore',
            'term': search_term,
            'retmax': str(self.max_results),
            'retmode': 'xml',
            'tool': 'GenomeHarvester',
            'email': self.email
        }
        
        if self.api_key:
            params['api_key'] = self.api_key
        
        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(self.esearch_base, params=params) as response:
                    if response.status != 200:
                        raise Exception(f"NCBI search failed: HTTP {response.status}")
                    
                    xml_content = await response.text()
                    
                    # Parse XML response
                    root = ET.fromstring(xml_content)
                    
                    # Extract ID list
                    id_list = []
                    for id_elem in root.findall('.//Id'):
                        if id_elem.text:
                            id_list.append(id_elem.text)
                    
                    self.logger.info(f"Found {len(id_list)} genome IDs")
                    
                    if not id_list:
                        self.logger.warning("No genomes found for search term")
                        return []
                    
                    # Get accession numbers from IDs
                    accessions = await self._get_accessions_from_ids(id_list)
                    return accessions
                    
        except Exception as e:
            self.logger.error(f"NCBI search failed: {e}")
            raise
    
    async def _get_accessions_from_ids(self, id_list: List[str]) -> List[str]:
        """
        Convert NCBI IDs to accession numbers using esummary
        """
        self.logger.info(f"Converting {len(id_list)} IDs to accessions")
        
        params = {
            'db': 'nuccore',
            'id': ','.join(id_list),
            'retmode': 'xml',
            'tool': 'GenomeHarvester',
            'email': self.email
        }
        
        if self.api_key:
            params['api_key'] = self.api_key
        
        accessions = []
        
        try:
            async with aiohttp.ClientSession() as session:
                async with session.get(self.esummary_base, params=params) as response:
                    if response.status != 200:
                        raise Exception(f"NCBI esummary failed: HTTP {response.status}")
                    
                    xml_content = await response.text()
                    root = ET.fromstring(xml_content)
                    
                    # Parse summary information
                    for doc_sum in root.findall('.//DocSum'):
                        accession = None
                        title = ""
                        organism = ""
                        length = 0
                        
                        # Extract accession
                        for item in doc_sum.findall('.//Item'):
                            if item.get('Name') == 'AccessionVersion':
                                accession = item.text
                            elif item.get('Name') == 'Title':
                                title = item.text or ""
                            elif item.get('Name') == 'Organism':
                                organism = item.text or ""
                            elif item.get('Name') == 'Length':
                                try:
                                    length = int(item.text or 0)
                                except ValueError:
                                    length = 0
                        
                        if accession:
                            accessions.append(accession)
                            
                            # Store genome information
                            self.genome_results[accession] = GenomeResult(
                                accession=accession,
                                title=title,
                                organism=organism,
                                length=length
                            )
                    
                    self.logger.info(f"Successfully converted to {len(accessions)} accessions")
                    return accessions
                    
        except Exception as e:
            self.logger.error(f"Failed to get accessions from IDs: {e}")
            # Fallback: use IDs as accessions (less reliable)
            return id_list[:self.max_results]
    
    async def _download_single_genome(self, 
                                    session: aiohttp.ClientSession, 
                                    accession: str) -> bool:
        """Download a single genome with retries and error handling"""
        
        if accession not in self.genome_results:
            self.genome_results[accession] = GenomeResult(accession=accession)
        
        result = self.genome_results[accession]
        result.download_attempted = True
        
        output_file = self.output_dir / f"{accession}.fasta"
        
        # Skip if already exists and valid
        if output_file.exists() and output_file.stat().st_size > 0:
            result.download_success = True
            result.file_path = str(output_file)
            result.file_size_bytes = output_file.stat().st_size
            self.logger.info(f"Genome already exists: {accession}")
            return True
        
        start_time = time.time()
        
        for attempt in range(self.retry_attempts):
            try:
                self.logger.info(f"Downloading {accession} (attempt {attempt + 1}/{self.retry_attempts})")
                
                # Use efetch to download FASTA
                params = {
                    'db': 'nuccore',
                    'id': accession,
                    'rettype': 'fasta',
                    'retmode': 'text',
                    'tool': 'GenomeHarvester',
                    'email': self.email
                }
                
                if self.api_key:
                    params['api_key'] = self.api_key
                
                async with session.get(self.efetch_base, params=params) as response:
                    if response.status == 200:
                        content = await response.read()
                        
                        # Basic validation - should start with '>'
                        if content.startswith(b'>'):
                            # Write to file
                            async with aiofiles.open(output_file, 'wb') as f:
                                await f.write(content)
                            
                            # Validate downloaded file
                            if output_file.stat().st_size > 0:
                                result.download_success = True
                                result.file_path = str(output_file)
                                result.file_size_bytes = output_file.stat().st_size
                                result.download_time_seconds = time.time() - start_time
                                
                                self.logger.info(f"Successfully downloaded {accession} "
                                               f"({result.file_size_bytes} bytes)")
                                return True
                            else:
                                raise Exception("Downloaded file is empty")
                        else:
                            raise Exception("Downloaded content is not FASTA format")
                    
                    else:
                        raise Exception(f"HTTP {response.status}: {response.reason}")
            
            except Exception as e:
                self.logger.warning(f"Download attempt {attempt + 1} failed for {accession}: {e}")
                if attempt == self.retry_attempts - 1:
                    result.error_message = f"Download failed after {self.retry_attempts} attempts: {e}"
                    self.logger.error(f"Failed to download {accession}: {e}")
                    return False
                
                # Wait before retry
                await asyncio.sleep(2 ** attempt)
        
        return False
    
    async def download_genomes_from_url(self, ncbi_url: str) -> Dict[str, GenomeResult]:
        """
        Complete workflow: parse URL, search NCBI, download genomes
        Returns: Dictionary of accession -> GenomeResult
        """
        self.logger.info(f"Starting genome harvest from URL")
        
        try:
            # Parse search term from URL
            search_term = self.parse_ncbi_url(ncbi_url)
            
            # Search NCBI for accessions
            accessions = await self.search_ncbi(search_term)
            
            if not accessions:
                self.logger.warning("No accessions found from search")
                return self.genome_results
            
            self.logger.info(f"Found {len(accessions)} accessions, downloading first {min(len(accessions), self.max_results)}")
            
            # Limit to max_results
            accessions = accessions[:self.max_results]
            
            # Download genomes in parallel
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
            successful = sum(1 for result in self.genome_results.values() if result.download_success)
            failed = len(accessions) - successful
            total_size = sum(result.file_size_bytes for result in self.genome_results.values() 
                            if result.download_success)
            
            self.logger.info(f"Download complete: {successful} successful, {failed} failed, "
                            f"{total_size / 1024 / 1024:.2f} MB total")
            
            return self.genome_results
            
        except Exception as e:
            self.logger.error(f"Genome harvest failed: {e}")
            raise
    
    def generate_detailed_report(self, output_file: Optional[str] = None) -> str:
        """Generate detailed harvest report"""
        
        if output_file is None:
            output_file = self.output_dir / "harvest_report.txt"
        
        successful = [r for r in self.genome_results.values() if r.download_success]
        failed = [r for r in self.genome_results.values() if not r.download_success and r.download_attempted]
        
        report = f"""
URL-Based Genome Harvest Report
==============================
Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}

Summary:
--------
Total genomes found: {len(self.genome_results)}
Successful downloads: {len(successful)}
Failed downloads: {len(failed)}
Total data downloaded: {sum(r.file_size_bytes for r in successful) / 1024 / 1024:.2f} MB

Successful Downloads:
--------------------
"""
        
        for result in successful:
            report += f"{result.accession}\t{result.organism[:50]}\t{result.file_size_bytes / 1024:.1f} KB\t{result.download_time_seconds:.1f}s\n"
        
        if failed:
            report += "\nFailed Downloads:\n-----------------\n"
            for result in failed:
                report += f"{result.accession}\t{result.error_message}\n"
        
        with open(output_file, 'w') as f:
            f.write(report)
        
        self.logger.info(f"Detailed report saved to: {output_file}")
        return report


async def main():
    """Example usage and testing"""
    
    # Test URL (E. coli with macrolide resistance)
    test_url = "https://www.ncbi.nlm.nih.gov/nuccore/?term=(Escherichia+coli+and+macrolide+resistance)+AND+%22Escherichia+coli%22%5Bporgn%3A__txid562%5D+and+complete+genome"
    
    # Configuration
    harvester = URLGenomeHarvester(
        output_dir='url_genomes',
        email='vihaankulkarni29@gmail.com',
        max_results=10,  # First 10 genomes
        max_concurrent=3,
        retry_attempts=2
    )
    
    try:
        # Run the harvest
        results = await harvester.download_genomes_from_url(test_url)
        
        # Generate report
        harvester.generate_detailed_report()
        
        print(f"\n🎉 Harvest complete! Downloaded {len([r for r in results.values() if r.download_success])} genomes")
        
    except Exception as e:
        print(f"❌ Harvest failed: {e}")


if __name__ == "__main__":
    asyncio.run(main())