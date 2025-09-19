"""
Genome Download Engine
=====================
Professional genome downloader that takes discovered accessions and downloads
FASTA files with robust error handling and parallel processing.
"""

import os
import asyncio
import aiohttp
import aiofiles
import hashlib
import time
from typing import List, Dict, Optional, Tuple
from pathlib import Path
from dataclasses import dataclass
import logging
from ncbi_genome_discovery import GenomeMetadata


@dataclass
class DownloadResult:
    """Track download results for each genome with enhanced metadata"""
    accession: str
    success: bool = False
    file_path: Optional[str] = None
    file_size_bytes: int = 0
    download_time_seconds: float = 0.0
    error_message: Optional[str] = None
    retry_count: int = 0
    # Enhanced metadata for provenance
    organism: str = ""
    genome_length: int = 0
    download_timestamp: str = ""
    ncbi_url: str = ""
    file_checksum: Optional[str] = None


@dataclass
class GenomeDownloaderConfig:
    """Configuration settings for GenomeDownloader"""
    output_dir: str
    email: str
    api_key: Optional[str] = None
    max_concurrent: int = 3
    retry_attempts: int = 3
    timeout_seconds: int = 300


class GenomeDownloader:
    """
    Professional genome downloader with:
    - Parallel downloads with rate limiting
    - Comprehensive error handling and retries
    - FASTA validation and quality control
    - Progress tracking and detailed reporting
    """
    
    def __init__(self,
                 config: GenomeDownloaderConfig,
                 logger: Optional[logging.Logger] = None):
        
        self.output_dir = Path(config.output_dir)
        self.email = config.email
        self.api_key = config.api_key
        self.max_concurrent = config.max_concurrent
        self.retry_attempts = config.retry_attempts
        self.timeout_seconds = config.timeout_seconds
        self.logger = logger or self._setup_logger()
        
        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Download tracking
        self.download_results: Dict[str, DownloadResult] = {}
        
        # NCBI efetch endpoint
        self.efetch_base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    
    def _setup_logger(self) -> logging.Logger:
        """Setup comprehensive logging"""
        logger = logging.getLogger("GenomeDownloader")
        if not logger.handlers:
            # Console handler
            console_handler = logging.StreamHandler()
            console_formatter = logging.Formatter('%(levelname)s %(name)s: %(message)s')
            console_handler.setFormatter(console_formatter)
            logger.addHandler(console_handler)
            
            # File handler
            log_file = self.output_dir / "download.log"
            file_handler = logging.FileHandler(log_file)
            file_formatter = logging.Formatter('%(asctime)s %(levelname)s %(name)s: %(message)s')
            file_handler.setFormatter(file_formatter)
            logger.addHandler(file_handler)
        
        logger.setLevel(logging.INFO)
        return logger
    
    async def download_genomes(self, genomes: List[GenomeMetadata]) -> Dict[str, DownloadResult]:
        """
        Download genomes with parallel processing
        
        Args:
            genomes: List of GenomeMetadata objects to download
            
        Returns:
            Dictionary mapping accession to DownloadResult
        """
        self.logger.info(f"Starting download of {len(genomes)} genomes")
        
        # Initialize download results
        for genome in genomes:
            self.download_results[genome.accession] = DownloadResult(accession=genome.accession)
        
        # Setup connection limits and timeout
        connector = aiohttp.TCPConnector(limit=self.max_concurrent)
        timeout = aiohttp.ClientTimeout(total=self.timeout_seconds)
        
        async with aiohttp.ClientSession(connector=connector, timeout=timeout) as session:
            # Create semaphore for rate limiting
            semaphore = asyncio.Semaphore(self.max_concurrent)
            
            # Create download tasks
            tasks = []
            for genome in genomes:
                task = self._download_with_semaphore(session, semaphore, genome)
                tasks.append(task)
            
            # Execute downloads
            await asyncio.gather(*tasks, return_exceptions=True)
        
        # Generate summary and reports
        self._log_download_summary()
        
        # Generate detailed reports
        try:
            self.generate_download_report()
            self.generate_metadata_manifest()
            self.logger.info("Generated download report and metadata manifest")
        except Exception as e:
            self.logger.warning(f"Failed to generate reports: {e}")
        
        return self.download_results
    
    async def _download_with_semaphore(self, 
                                     session: aiohttp.ClientSession,
                                     semaphore: asyncio.Semaphore,
                                     genome: GenomeMetadata) -> bool:
        """Download single genome with semaphore control"""
        async with semaphore:
            return await self._download_single_genome(session, genome)
    
    async def _download_single_genome(self, 
                                    session: aiohttp.ClientSession,
                                    genome: GenomeMetadata) -> bool:
        """Download a single genome with retries and validation"""
        
        result = self.download_results[genome.accession]
        # Create output filename using accession-centric naming
        output_file = self.output_dir / self._create_accession_filename(genome.accession)
        
        # Skip if already exists and valid
        if await self._file_exists_and_valid(output_file):
            result.success = True
            result.file_path = str(output_file)
            result.file_size_bytes = output_file.stat().st_size
            self.logger.info(f"Genome already exists: {genome.accession}")
            return True
        
        start_time = time.time()
        
        # Retry loop
        for attempt in range(self.retry_attempts):
            try:
                result.retry_count = attempt
                self.logger.info(f"Downloading {genome.accession} (attempt {attempt + 1}/{self.retry_attempts})")
                
                # Prepare download parameters
                params = {
                    'db': 'nuccore',
                    'id': genome.accession,
                    'rettype': 'fasta',
                    'retmode': 'text',
                    'tool': 'GenomeHarvester',
                    'email': self.email
                }
                
                if self.api_key:
                    params['api_key'] = self.api_key
                
                # Download the genome
                async with session.get(self.efetch_base, params=params) as response:
                    if response.status == 200:
                        content = await response.read()
                        
                        # Validate content
                        if await self._validate_fasta_content(content, genome):
                            # Extract metadata from FASTA content
                            content_str = content.decode('utf-8', errors='ignore')
                            organism, genome_length = self._extract_metadata_from_fasta(content_str, genome.accession)
                            
                            # Update filename with organism if available
                            if organism:
                                output_file = self.output_dir / self._create_accession_filename(genome.accession, organism)
                            
                            # Write to file
                            async with aiofiles.open(output_file, 'wb') as f:
                                await f.write(content)
                            
                            # Calculate checksum
                            checksum = await self._calculate_file_checksum(output_file)
                            
                            # Final validation
                            if await self._validate_downloaded_file(output_file, genome):
                                result.success = True
                                result.file_path = str(output_file)
                                result.file_size_bytes = output_file.stat().st_size
                                result.download_time_seconds = time.time() - start_time
                                result.organism = organism
                                result.genome_length = genome_length
                                result.download_timestamp = time.strftime('%Y-%m-%d %H:%M:%S')
                                result.ncbi_url = f"{self.efetch_base}?{'&'.join(f'{k}={v}' for k, v in params.items())}"
                                result.file_checksum = checksum
                                
                                self.logger.info(f"Successfully downloaded {genome.accession} "
                                               f"({result.file_size_bytes / 1024:.1f} KB, {organism})")
                                return True
                            else:
                                raise Exception("Downloaded file failed validation")
                        else:
                            raise Exception("Downloaded content failed validation")
                    
                    else:
                        raise Exception(f"HTTP {response.status}: {response.reason}")
            
            except Exception as e:
                self.logger.warning(f"Download attempt {attempt + 1} failed for {genome.accession}: {e}")
                
                if attempt == self.retry_attempts - 1:
                    # Final failure
                    result.error_message = f"Download failed after {self.retry_attempts} attempts: {e}"
                    self.logger.error(f"Failed to download {genome.accession}: {e}")
                    return False
                
                # Wait before retry (exponential backoff)
                wait_time = 2 ** attempt
                await asyncio.sleep(wait_time)
        
        return False
    
    async def _file_exists_and_valid(self, file_path: Path) -> bool:
        """Check if file exists and is valid"""
        if not file_path.exists():
            return False
        
        if file_path.stat().st_size == 0:
            self.logger.warning(f"Removing empty file: {file_path}")
            file_path.unlink()
            return False
        
        # Quick FASTA validation
        try:
            async with aiofiles.open(file_path, 'rb') as f:
                first_bytes = await f.read(10)
                if not first_bytes.startswith(b'>'):
                    self.logger.warning(f"Removing invalid FASTA file: {file_path}")
                    file_path.unlink()
                    return False
        except Exception:
            return False
        
        return True
    
    async def _calculate_file_checksum(self, file_path: Path) -> str:
        """Calculate SHA256 checksum for downloaded file"""
        hash_sha256 = hashlib.sha256()
        try:
            async with aiofiles.open(file_path, 'rb') as f:
                chunk = await f.read(8192)
                while chunk:
                    hash_sha256.update(chunk)
                    chunk = await f.read(8192)
            return hash_sha256.hexdigest()
        except Exception as e:
            self.logger.warning(f"Failed to calculate checksum for {file_path}: {e}")
            return ""

    def _create_accession_filename(self, accession: str, organism: str = "") -> str:
        """Create standardized filename with accession as primary identifier"""
        # Clean organism name for safe filename
        clean_organism = "".join(c for c in organism if c.isalnum() or c in (' ', '-', '_')).strip()
        clean_organism = clean_organism.replace(' ', '_')
        
        if clean_organism and len(clean_organism) > 0:
            return f"{accession}_{clean_organism}.fasta"
        else:
            return f"{accession}.fasta"

    def _extract_metadata_from_fasta(self, content: str, accession: str) -> Tuple[str, int]:
        """Extract organism name and genome length from FASTA content"""
        organism = ""
        genome_length = 0
        
        lines = content.split('\n')
        for line in lines:
            if line.startswith('>'):
                # Extract organism from header
                if '[' in line and ']' in line:
                    start = line.find('[') + 1
                    end = line.find(']')
                    organism = line[start:end]
                    break
            elif not line.startswith('>') and line.strip():
                genome_length += len(line.strip())
        
        return organism, genome_length
    
    async def _validate_fasta_content(self, content: bytes, genome: GenomeMetadata) -> bool:
        """Validate downloaded FASTA content"""
        if len(content) < 100:
            return False
        
        # Must start with FASTA header
        if not content.startswith(b'>'):
            return False
        
        # Check for obvious error messages
        content_str = content.decode('utf-8', errors='ignore').lower()
        error_indicators = ['error', 'not found', 'invalid', 'failed']
        if any(indicator in content_str[:200] for indicator in error_indicators):
            return False
        
        # Should contain nucleotide sequences
        sequence_chars = set('atcgn')
        lines = content_str.split('\n')[1:]  # Skip header
        sequence_content = ''.join(lines).replace(' ', '').replace('\t', '')
        
        if len(sequence_content) < 1000:  # Very short sequence
            return False
        
        # Check if it looks like nucleotide sequence
        nucleotide_ratio = sum(1 for c in sequence_content.lower() if c in sequence_chars) / len(sequence_content)
        if nucleotide_ratio < 0.8:  # Less than 80% nucleotides
            return False
        
        return True
    
    async def _validate_downloaded_file(self, file_path: Path, genome: GenomeMetadata) -> bool:
        """Validate downloaded file on disk"""
        try:
            if not file_path.exists() or file_path.stat().st_size == 0:
                return False
            
            # Read first few lines to validate
            async with aiofiles.open(file_path, 'r') as f:
                first_line = await f.readline()
                if not first_line.startswith('>'):
                    return False
                
                # Read a bit of sequence
                second_line = await f.readline()
                if len(second_line.strip()) < 10:
                    return False
            
            return True
            
        except Exception as e:
            self.logger.warning(f"File validation failed for {file_path}: {e}")
            return False
    
    def _log_download_summary(self):
        """Log download summary statistics"""
        successful = [r for r in self.download_results.values() if r.success]
        failed = [r for r in self.download_results.values() if not r.success]
        
        total_size = sum(r.file_size_bytes for r in successful)
        avg_time = sum(r.download_time_seconds for r in successful) / len(successful) if successful else 0
        
        self.logger.info(f"Download Summary:")
        self.logger.info(f"  Successful: {len(successful)}")
        self.logger.info(f"  Failed: {len(failed)}")
        self.logger.info(f"  Total size: {total_size / 1024 / 1024:.2f} MB")
        self.logger.info(f"  Average time: {avg_time:.1f} seconds")
        
        if failed:
            self.logger.warning(f"Failed downloads:")
            for result in failed:
                self.logger.warning(f"  {result.accession}: {result.error_message}")
    
    def get_successful_downloads(self) -> List[str]:
        """Get list of successfully downloaded accessions"""
        return [r.accession for r in self.download_results.values() if r.success]
    
    def get_downloaded_files(self) -> List[str]:
        """Get list of successfully downloaded file paths"""
        return [r.file_path for r in self.download_results.values() if r.success and r.file_path]
    
    def generate_download_report(self, output_file: Optional[str] = None) -> str:
        """Generate detailed download report"""
        
        if output_file is None:
            output_file = str(self.output_dir / "download_report.txt")
        
        successful = [r for r in self.download_results.values() if r.success]
        failed = [r for r in self.download_results.values() if not r.success]
        
        report = f"""
Genome Download Report
=====================
Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}

Summary:
--------
Total genomes: {len(self.download_results)}
Successful downloads: {len(successful)}
Failed downloads: {len(failed)}
Success rate: {len(successful) / len(self.download_results) * 100:.1f}%
Total data downloaded: {sum(r.file_size_bytes for r in successful) / 1024 / 1024:.2f} MB

Successful Downloads:
--------------------
Accession\tOrganism\tFile Size\tGenome Length\tDownload Time\tRetries\tChecksum\tTimestamp
"""
        
        for result in successful:
            organism_display = result.organism[:30] + "..." if len(result.organism) > 30 else result.organism
            checksum_short = result.file_checksum[:8] if result.file_checksum else "N/A"
            report += f"{result.accession}\t{organism_display}\t{result.file_size_bytes / 1024:.1f} KB\t{result.genome_length:,} bp\t{result.download_time_seconds:.1f}s\t{result.retry_count}\t{checksum_short}\t{result.download_timestamp}\n"
        
        if failed:
            report += "\nFailed Downloads:\n-----------------\n"
            for result in failed:
                report += f"{result.accession}\t{result.error_message}\n"
        
        # Write report to file
        with open(output_file, 'w') as f:
            f.write(report)
        
        self.logger.info(f"Download report saved to: {output_file}")
        return report

    def generate_metadata_manifest(self, output_file: Optional[str] = None) -> str:
        """Generate JSON manifest with complete metadata for all downloads"""
        import json
        
        if output_file is None:
            output_file = str(self.output_dir / "genome_manifest.json")
        
        successful = [r for r in self.download_results.values() if r.success]
        
        manifest = {
            "manifest_version": "1.0",
            "generated_timestamp": time.strftime('%Y-%m-%d %H:%M:%S'),
            "total_genomes": len(successful),
            "genomes": []
        }
        
        for result in successful:
            genome_entry = {
                "accession": result.accession,
                "file_path": result.file_path,
                "organism": result.organism,
                "genome_length_bp": result.genome_length,
                "file_size_bytes": result.file_size_bytes,
                "download_timestamp": result.download_timestamp,
                "file_checksum_sha256": result.file_checksum,
                "ncbi_source_url": result.ncbi_url,
                "download_time_seconds": result.download_time_seconds,
                "retry_count": result.retry_count
            }
            manifest["genomes"].append(genome_entry)
        
        # Sort by accession for consistent ordering
        manifest["genomes"].sort(key=lambda x: x["accession"])
        
        # Write JSON manifest
        with open(output_file, 'w') as f:
            json.dump(manifest, f, indent=2)
        
        self.logger.info(f"Metadata manifest saved to: {output_file}")
        return output_file


# Integration test function
async def test_download_engine():
    """Test the download engine with a few real accessions"""
    from ncbi_genome_discovery import GenomeMetadata
    
    # Create test genomes (real E. coli accessions)
    test_genomes = [
        GenomeMetadata(accession="NC_000913.3", organism="Escherichia coli K-12", length=4641652),
        GenomeMetadata(accession="CP063194.1", organism="Escherichia coli", length=5176726),
        GenomeMetadata(accession="CP081855.1", organism="Escherichia coli", length=5098066)
    ]
    
    # Create configuration
    config = GenomeDownloaderConfig(
        output_dir="test_genomes",
        email="vihaankulkarni29@gmail.com",
        api_key="ef7622c2e716fa317fe04d24c42904211107",
        max_concurrent=2,
        retry_attempts=2
    )
    
    # Create downloader
    downloader = GenomeDownloader(config)
    
    try:
        # Download genomes
        results = await downloader.download_genomes(test_genomes)
        
        # Generate report
        report = downloader.generate_download_report()
        
        print(f"\n🎉 Download test complete!")
        print(f"Successful downloads: {len(downloader.get_successful_downloads())}")
        print(f"Downloaded files: {downloader.get_downloaded_files()}")
        
    except Exception as e:
        print(f"❌ Download test failed: {e}")


if __name__ == "__main__":
    asyncio.run(test_download_engine())