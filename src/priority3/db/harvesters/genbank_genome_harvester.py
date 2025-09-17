"""
GenBank Genome Harvester with E-utilities Integration
====================================================

Robust NCBI GenBank genome downloading with:
- E-utilities (esearch, esummary, efetch) integration
- Rate limiting and API key support
- Checkpointing and resume capability
- Comprehensive metadata extraction
- File integrity validation
- Mock mode for testing
"""

import os
import time
import logging
import requests
import hashlib
from datetime import datetime
from typing import List, Dict, Any, Optional, Union, Generator
from pathlib import Path
from dataclasses import dataclass
import xml.etree.ElementTree as ET
import json
from urllib.parse import urlencode
import gzip
import tempfile
import shutil

from ..repositories import GenomeRepository, GenomeRecord, calculate_file_hash

@dataclass
class SearchResult:
    """NCBI search result."""
    query: str
    total_count: int
    accessions: List[str]
    web_env: Optional[str] = None
    query_key: Optional[str] = None

@dataclass  
class GenomeMetadata:
    """Extended genome metadata from NCBI."""
    accession: str
    organism: str
    strain: Optional[str] = None
    biosample: Optional[str] = None
    bioproject: Optional[str] = None
    assembly_level: Optional[str] = None
    genome_size: Optional[int] = None
    gc_content: Optional[float] = None
    sequencing_tech: Optional[str] = None
    submitter: Optional[str] = None
    release_date: Optional[str] = None
    ftp_path: Optional[str] = None
    assembly_name: Optional[str] = None

class NCBIRateLimiter:
    """NCBI API rate limiting."""
    
    def __init__(self, requests_per_second: float = 3.0, api_key: Optional[str] = None):
        # With API key: 10 req/sec, without: 3 req/sec
        self.delay = 1.0 / (10.0 if api_key else requests_per_second)
        self.last_request = 0.0
        self.api_key = api_key
        
    def wait(self):
        """Wait if necessary to respect rate limits."""
        now = time.time()
        elapsed = now - self.last_request
        if elapsed < self.delay:
            time.sleep(self.delay - elapsed)
        self.last_request = time.time()

class GenBankGenomeHarvester:
    """
    Harvest genomes from NCBI GenBank with robust error handling.
    """
    
    def __init__(self, 
                 output_dir: str = "genomes",
                 db_path: str = "priority3.db",
                 api_key: Optional[str] = None,
                 email: Optional[str] = None,
                 mock_mode: bool = False):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        self.repository = GenomeRepository(db_path)
        self.rate_limiter = NCBIRateLimiter(api_key=api_key)
        self.api_key = api_key
        self.email = email or "user@example.com"
        self.mock_mode = mock_mode
        
        self.logger = logging.getLogger("GenBankHarvester")
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'GenomeAMRAnalyzer/1.0 (https://github.com/user/GenomeAMRAnalyzer)'
        })
        
        # NCBI base URLs
        self.eutils_base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
        self.ftp_base = "https://ftp.ncbi.nlm.nih.gov"
        
        # Checkpoint file for resuming
        self.checkpoint_file = self.output_dir / "harvest_checkpoint.json"
        
    def search_genomes(self, 
                      query: str, 
                      database: str = "assembly",
                      max_results: Optional[int] = None) -> SearchResult:
        """
        Search NCBI for genomes matching query.
        
        Args:
            query: NCBI search query (e.g., "Escherichia coli AND AMR")
            database: NCBI database (assembly, nucleotide)
            max_results: Maximum results to return
            
        Returns:
            SearchResult with accessions and metadata
        """
        if self.mock_mode:
            return self._mock_search_results(query, max_results or 100)
            
        try:
            # Step 1: esearch to get IDs
            search_params = {
                'db': database,
                'term': query,
                'usehistory': 'y',
                'tool': 'GenomeAMRAnalyzer',
                'email': self.email
            }
            
            if self.api_key:
                search_params['api_key'] = self.api_key
                
            if max_results is not None:
                search_params['retmax'] = str(max_results)
                
            self.rate_limiter.wait()
            response = self.session.get(
                f"{self.eutils_base}/esearch.fcgi",
                params=search_params
            )
            response.raise_for_status()
            
            # Parse XML response
            root = ET.fromstring(response.content)
            
            count_text = None
            count_node = root.find('Count')
            if count_node is not None and count_node.text is not None:
                count_text = count_node.text
            total_count = int(count_text) if count_text else 0

            web_env_node = root.find('WebEnv')
            web_env = web_env_node.text if web_env_node is not None else None

            qk_node = root.find('QueryKey')
            query_key = qk_node.text if qk_node is not None else None
            
            # Get ID list
            id_list = root.find('IdList')
            accessions = []
            if id_list is not None:
                for id_elem in id_list.findall('Id'):
                    if id_elem is not None and id_elem.text:
                        accessions.append(id_elem.text)
            
            self.logger.info(f"Found {total_count} genomes for query: {query}")
            self.logger.info(f"Retrieved {len(accessions)} accession IDs")
            
            return SearchResult(
                query=query,
                total_count=total_count,
                accessions=accessions,
                web_env=web_env,
                query_key=query_key
            )
            
        except Exception as e:
            self.logger.error(f"Failed to search genomes: {e}")
            raise
            
    def get_genome_metadata(self, accessions: List[str], database: str = "assembly") -> List[GenomeMetadata]:
        """Get detailed metadata for genome accessions."""
        if self.mock_mode:
            return self._mock_genome_metadata(accessions)
            
        try:
            metadata_list = []
            
            # Process in batches to avoid URL length limits
            batch_size = 100
            for i in range(0, len(accessions), batch_size):
                batch = accessions[i:i + batch_size]
                
                # esummary to get metadata
                summary_params = {
                    'db': database,
                    'id': ','.join(batch),
                    'tool': 'GenomeAMRAnalyzer',
                    'email': self.email
                }
                
                if self.api_key:
                    summary_params['api_key'] = self.api_key
                    
                self.rate_limiter.wait()
                response = self.session.get(
                    f"{self.eutils_base}/esummary.fcgi",
                    params=summary_params
                )
                response.raise_for_status()
                
                # Parse metadata from XML (handle eSummary XML structure)
                root = ET.fromstring(response.content)
                doc_summaries = root.findall('.//DocumentSummary')
                if not doc_summaries:
                    # Some endpoints/databases may use DocSum
                    doc_summaries = root.findall('.//DocSum')
                if not doc_summaries:
                    preview = response.text[:300].replace('\n', ' ')
                    self.logger.warning(
                        f"No DocumentSummary nodes found in eSummary response. Preview: {preview}"
                    )
                for doc_sum in doc_summaries:
                    metadata = self._parse_genome_metadata(doc_sum)
                    if metadata:
                        metadata_list.append(metadata)
                        
            self.logger.info(f"Retrieved metadata for {len(metadata_list)} genomes")
            return metadata_list
            
        except Exception as e:
            self.logger.error(f"Failed to get genome metadata: {e}")
            raise
            
    def _parse_genome_metadata(self, doc_sum: ET.Element) -> Optional[GenomeMetadata]:
        """Parse genome metadata from eSummary DocumentSummary/DocSum nodes."""
        try:
            # Accession resolution across possible layouts
            accession = (
                doc_sum.get('uid') or
                self._get_text_safe(doc_sum, './/AssemblyAccession') or
                self._get_text_safe(doc_sum, "./Item[@Name='AssemblyAccession']") or
                self._get_text_safe(doc_sum, './/Accession') or
                self._get_text_safe(doc_sum, "./Item[@Name='Accession']")
            )
            if not accession:
                return None
            # Extract basic info (support DocSum Item[@Name] style)
            get_item = lambda name: self._get_text_safe(doc_sum, f"./Item[@Name='{name}']")
            organism = (
                self._get_text_safe(doc_sum, './/Organism') or
                get_item('Organism')
            )
            strain = self._get_text_safe(doc_sum, './/Strain') or get_item('Strain')
            # BioSample and BioProject
            biosample = self._get_text_safe(doc_sum, './/BioSampleAccn') or get_item('BioSampleAccn')
            bioproject = self._get_text_safe(doc_sum, './/BioprojectAccn') or get_item('BioprojectAccn')
            # Assembly info
            assembly_level = self._get_text_safe(doc_sum, './/AssemblyLevel') or get_item('AssemblyLevel')
            # Size and composition
            genome_size = self._get_int_safe(doc_sum, './/TotalLength') or self._get_int_safe(doc_sum, "./Item[@Name='TotalLength']")
            gc_content = self._get_float_safe(doc_sum, './/GC') or self._get_float_safe(doc_sum, "./Item[@Name='GC']")
            # Additional metadata
            sequencing_tech = self._get_text_safe(doc_sum, './/SequencingTech') or get_item('SequencingTech')
            submitter = self._get_text_safe(doc_sum, './/SubmitterOrganization') or get_item('SubmitterOrganization')
            release_date = self._get_text_safe(doc_sum, './/AsmReleaseDate') or get_item('AsmReleaseDate')
            # FTP paths and assembly name (used to construct download URL)
            ftp_path = (
                self._get_text_safe(doc_sum, './/FtpPath_RefSeq') or
                self._get_text_safe(doc_sum, './/FtpPath_GenBank') or
                get_item('FtpPath_RefSeq') or
                get_item('FtpPath_GenBank')
            )
            assembly_name = self._get_text_safe(doc_sum, './/AssemblyName') or get_item('AssemblyName')
            return GenomeMetadata(
                accession=accession,
                organism=organism or "Unknown",
                strain=strain,
                biosample=biosample,
                bioproject=bioproject,
                assembly_level=assembly_level,
                genome_size=genome_size,
                gc_content=gc_content,
                sequencing_tech=sequencing_tech,
                submitter=submitter,
                release_date=release_date,
                ftp_path=ftp_path,
                assembly_name=assembly_name,
            )
            
        except Exception as e:
            self.logger.error(f"Failed to parse metadata: {e}")
            return None
            
    def download_genome(self, metadata: GenomeMetadata, retries: int = 3) -> bool:
        """
        Download genome FASTA file with retries and checkpointing.
        """
        if self.mock_mode:
            return self._mock_download_genome(metadata)
            
        output_file = self.output_dir / f"{metadata.accession}.fasta"
        
        # Check if already downloaded and valid
        if output_file.exists():
            existing_record = self.repository.get_genome(metadata.accession)
            if existing_record and existing_record.status == "downloaded":
                self.logger.info(f"Genome {metadata.accession} already downloaded")
                return True
                
        try:
            # Add to database as pending
            genome_record = GenomeRecord(
                accession=metadata.accession,
                organism=metadata.organism,
                biosample=metadata.biosample,
                bioproject=metadata.bioproject,
                assembly_level=metadata.assembly_level,
                genome_size=metadata.genome_size,
                status="downloading"
            )
            self.repository.add_genome(genome_record)
            
            # Get download URL (prefer from metadata FTP path when available)
            download_url = self._get_download_url_from_metadata(metadata) or self._get_genome_download_url(metadata.accession)
            
            for attempt in range(retries):
                try:
                    self.logger.info(f"Downloading {metadata.accession} (attempt {attempt + 1})")
                    
                    self.rate_limiter.wait()
                    response = self.session.get(download_url, stream=True, timeout=300)
                    response.raise_for_status()
                    
                    # Decide if content is gzipped based on URL or headers
                    is_gzip = download_url.endswith('.gz') or \
                        response.headers.get('Content-Encoding', '').lower() == 'gzip' or \
                        'gzip' in response.headers.get('Content-Type', '').lower()

                    # Download to temporary file first
                    suffix = '.fasta.gz' if is_gzip else '.fasta'
                    with tempfile.NamedTemporaryFile(delete=False, suffix=suffix) as tmp_file:
                        for chunk in response.iter_content(chunk_size=8192):
                            tmp_file.write(chunk)
                        tmp_path = tmp_file.name

                    # If gzipped, decompress to a second temp file for validation/output
                    if is_gzip:
                        decompressed_path = None
                        try:
                            with gzip.open(tmp_path, 'rb') as gz_in, tempfile.NamedTemporaryFile(delete=False, suffix='.fasta') as tmp_out:
                                shutil.copyfileobj(gz_in, tmp_out)
                                decompressed_path = tmp_out.name
                        finally:
                            # Remove the gz temp file regardless
                            try:
                                os.unlink(tmp_path)
                            except Exception:
                                pass
                        tmp_path = decompressed_path
                        
                    # Validate downloaded file
                    if self._validate_downloaded_genome(tmp_path):
                        # Move to final location
                        # Use shutil.move to be cross-device safe on Windows
                        shutil.move(tmp_path, str(output_file))
                        
                        # Calculate hash and update database
                        file_hash = calculate_file_hash(str(output_file))
                        self.repository.update_genome_status(
                            metadata.accession, 
                            "downloaded",
                            str(output_file),
                            file_hash
                        )
                        
                        self.logger.info(f"Successfully downloaded {metadata.accession}")
                        return True
                    else:
                        os.unlink(tmp_path)
                        self.logger.warning(f"Downloaded file validation failed: {metadata.accession}")
                        
                except Exception as e:
                    self.logger.warning(f"Download attempt {attempt + 1} failed: {e}")
                    if attempt < retries - 1:
                        time.sleep(2 ** attempt)  # Exponential backoff
                        
            # All retries failed
            self.repository.update_genome_status(metadata.accession, "error")
            self.logger.error(f"Failed to download {metadata.accession} after {retries} attempts")
            return False
            
        except Exception as e:
            self.repository.update_genome_status(metadata.accession, "error")
            self.logger.error(f"Download failed for {metadata.accession}: {e}")
            return False
            
    def harvest_genomes(self, 
                       query: str,
                       max_genomes: Optional[int] = None,
                       resume: bool = True) -> Dict[str, Any]:
        """
        Complete genome harvesting workflow.
        
        Args:
            query: NCBI search query
            max_genomes: Maximum genomes to download
            resume: Whether to resume from checkpoint
            
        Returns:
            Summary of harvesting results
        """
        start_time = datetime.now()
        
        # Load checkpoint if resuming
        checkpoint = self._load_checkpoint() if resume else {}
        
        try:
            # Search for genomes
            if 'search_result' not in checkpoint:
                self.logger.info(f"Searching for genomes: {query}")
                search_result = self.search_genomes(query, max_results=max_genomes)
                checkpoint['search_result'] = {
                    'query': search_result.query,
                    'total_count': search_result.total_count,
                    'accessions': search_result.accessions
                }
                self._save_checkpoint(checkpoint)
            else:
                search_result = SearchResult(**checkpoint['search_result'])
                self.logger.info("Resuming from checkpoint - search results loaded")
                
            # Get metadata
            if 'metadata' not in checkpoint:
                self.logger.info("Getting genome metadata...")
                metadata_list = self.get_genome_metadata(search_result.accessions)
                checkpoint['metadata'] = [
                    {
                        'accession': m.accession,
                        'organism': m.organism,
                        'biosample': m.biosample,
                        'bioproject': m.bioproject,
                        'assembly_level': m.assembly_level,
                        'genome_size': m.genome_size
                    } for m in metadata_list
                ]
                self._save_checkpoint(checkpoint)
            else:
                metadata_list = [GenomeMetadata(**m) for m in checkpoint['metadata']]
                self.logger.info("Resuming from checkpoint - metadata loaded")
                
            # Download genomes
            downloaded = checkpoint.get('downloaded', [])
            failed = checkpoint.get('failed', [])
            
            for metadata in metadata_list:
                if metadata.accession in downloaded:
                    continue
                    
                if metadata.accession in failed:
                    continue
                    
                success = self.download_genome(metadata)
                if success:
                    downloaded.append(metadata.accession)
                else:
                    failed.append(metadata.accession)
                    
                # Update checkpoint
                checkpoint['downloaded'] = downloaded
                checkpoint['failed'] = failed
                self._save_checkpoint(checkpoint)
                
            # Final summary
            end_time = datetime.now()
            duration = (end_time - start_time).total_seconds()
            
            summary = {
                'query': query,
                'total_found': search_result.total_count,
                'metadata_retrieved': len(metadata_list),
                'successfully_downloaded': len(downloaded),
                'failed_downloads': len(failed),
                'duration_seconds': duration,
                'start_time': start_time.isoformat(),
                'end_time': end_time.isoformat()
            }
            
            self.logger.info(f"Harvest complete: {len(downloaded)}/{len(metadata_list)} genomes downloaded")
            return summary
            
        except Exception as e:
            self.logger.error(f"Harvest failed: {e}")
            raise
            
    def _get_genome_download_url(self, accession: str) -> str:
        """Get download URL for genome (simplified)."""
        # In practice, this would parse NCBI FTP structure
        # For now, return a mock URL
        return f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{accession}/{accession}_genomic.fna.gz"

    def _get_download_url_from_metadata(self, metadata: GenomeMetadata) -> Optional[str]:
        """Construct a reliable HTTPS URL for genomic FASTA using metadata FTP path."""
        try:
            if not metadata.ftp_path:
                return None
            ftp_path = metadata.ftp_path
            # Normalize to HTTPS
            if ftp_path.startswith('ftp://'):
                ftp_path = 'https://' + ftp_path[len('ftp://'):]
            # Basename determines file prefix
            base = os.path.basename(ftp_path.rstrip('/'))
            # genomic FASTA file
            return f"{ftp_path}/{base}_genomic.fna.gz"
        except Exception:
            return None

    def harvest_by_accessions(self, accessions: List[str], resume: bool = True) -> Dict[str, Any]:
        """Harvest genomes for a given list of assembly accessions (GCF_*/GCA_*)."""
        start_time = datetime.now()
        checkpoint = self._load_checkpoint() if resume else {}
        try:
            # Get metadata for provided accessions
            if 'metadata' not in checkpoint or checkpoint.get('accession_mode') is not True:
                metadata_list = self.get_genome_metadata(accessions)
                checkpoint['metadata'] = [
                    {
                        'accession': m.accession,
                        'organism': m.organism,
                        'biosample': m.biosample,
                        'bioproject': m.bioproject,
                        'assembly_level': m.assembly_level,
                        'genome_size': m.genome_size,
                        'ftp_path': m.ftp_path,
                        'assembly_name': m.assembly_name,
                    } for m in metadata_list
                ]
                checkpoint['accession_mode'] = True
                self._save_checkpoint(checkpoint)
            else:
                metadata_list = [GenomeMetadata(**m) for m in checkpoint['metadata']]

            downloaded = checkpoint.get('downloaded', [])
            failed = checkpoint.get('failed', [])

            for metadata in metadata_list:
                if metadata.accession in downloaded or metadata.accession in failed:
                    continue
                if self.download_genome(metadata):
                    downloaded.append(metadata.accession)
                else:
                    failed.append(metadata.accession)
                checkpoint['downloaded'] = downloaded
                checkpoint['failed'] = failed
                self._save_checkpoint(checkpoint)

            end_time = datetime.now()
            duration = (end_time - start_time).total_seconds()
            return {
                'total_requested': len(accessions),
                'metadata_retrieved': len(metadata_list),
                'successfully_downloaded': len(downloaded),
                'failed_downloads': len(failed),
                'duration_seconds': duration,
                'start_time': start_time.isoformat(),
                'end_time': end_time.isoformat()
            }
        except Exception as e:
            self.logger.error(f"Harvest by accessions failed: {e}")
            raise
        
    def _validate_downloaded_genome(self, file_path: str) -> bool:
        """Validate downloaded genome file."""
        try:
            # Check file size
            if os.path.getsize(file_path) < 1000:  # Too small
                return False
                
            # Check if it's gzipped
            is_gzipped = file_path.endswith('.gz')
            
            # Basic FASTA validation
            opener = gzip.open if is_gzipped else open
            mode = 'rt' if is_gzipped else 'r'
            
            with opener(file_path, mode) as f:
                first_line = f.readline()
                if not str(first_line).startswith('>'):
                    return False
                    
                # Check for sequence content
                sequence_found = False
                for _ in range(10):  # Check first few lines
                    line = f.readline()
                    if not line:
                        break
                    if not str(line).startswith('>') and str(line).strip():
                        sequence_found = True
                        break
                        
                return sequence_found
                
        except Exception:
            return False
            
    def _load_checkpoint(self) -> Dict[str, Any]:
        """Load checkpoint from file."""
        try:
            if self.checkpoint_file.exists():
                with open(self.checkpoint_file, 'r') as f:
                    return json.load(f)
        except Exception as e:
            self.logger.warning(f"Failed to load checkpoint: {e}")
        return {}
        
    def _save_checkpoint(self, checkpoint: Dict[str, Any]):
        """Save checkpoint to file."""
        try:
            with open(self.checkpoint_file, 'w') as f:
                json.dump(checkpoint, f, indent=2)
        except Exception as e:
            self.logger.warning(f"Failed to save checkpoint: {e}")
            
    def _get_text_safe(self, element: ET.Element, xpath: str) -> Optional[str]:
        """Safely extract text from XML element."""
        try:
            found = element.find(xpath)
            return found.text if found is not None else None
        except Exception:
            return None
            
    def _get_int_safe(self, element: ET.Element, xpath: str) -> Optional[int]:
        """Safely extract integer from XML element."""
        try:
            text = self._get_text_safe(element, xpath)
            return int(text) if text else None
        except Exception:
            return None
            
    def _get_float_safe(self, element: ET.Element, xpath: str) -> Optional[float]:
        """Safely extract float from XML element."""
        try:
            text = self._get_text_safe(element, xpath)
            return float(text) if text else None
        except Exception:
            return None
            
    # Mock methods for testing
    def _mock_search_results(self, query: str, max_results: int) -> SearchResult:
        """Generate mock search results."""
        accessions = [f"GCF_{1000000 + i:06d}.1" for i in range(min(max_results, 50))]
        return SearchResult(
            query=query,
            total_count=len(accessions),
            accessions=accessions
        )
        
    def _mock_genome_metadata(self, accessions: List[str]) -> List[GenomeMetadata]:
        """Generate mock genome metadata."""
        return [
            GenomeMetadata(
                accession=acc,
                organism="Escherichia coli",
                strain=f"strain_{i}",
                biosample=f"SAMN{10000000 + i}",
                bioproject=f"PRJNA{100000 + i}",
                assembly_level="Complete Genome",
                genome_size=5000000 + i * 1000
            ) for i, acc in enumerate(accessions)
        ]
        
    def _mock_download_genome(self, metadata: GenomeMetadata) -> bool:
        """Mock genome download."""
        output_file = self.output_dir / f"{metadata.accession}.fasta"
        
        # Create a simple mock FASTA file
        with open(output_file, 'w') as f:
            f.write(f">{metadata.accession} {metadata.organism}\n")
            f.write("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n")
            
        # Update database
        genome_record = GenomeRecord(
            accession=metadata.accession,
            organism=metadata.organism,
            biosample=metadata.biosample,
            bioproject=metadata.bioproject,
            file_path=str(output_file),
            file_hash=calculate_file_hash(str(output_file)),
            download_date=datetime.now(),
            assembly_level=metadata.assembly_level,
            genome_size=metadata.genome_size,
            status="downloaded"
        )
        self.repository.add_genome(genome_record)
        
        return True
        
    def close(self):
        """Close harvester and cleanup."""
        self.repository.close()
        self.session.close()