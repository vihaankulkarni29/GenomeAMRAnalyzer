#!/usr/bin/env python3
"""
GenomeDownloader - Download genome sequences and annotations from NCBI

This script downloads complete bacterial genome sequences (FASTA) and annotations (GFF)
from NCBI using protein accession numbers. It finds the corresponding genome assemblies
and downloads the genome data for mutation analysis.

Similar to NCBIProteinHarvester but for genome data instead of protein sequences.

Usage:
    python GenomeDownloader.py --accession-list data/accession_list.txt --email your.email@example.com --output-dir genome_data/

Author: MetaDataHarvester Pipeline
Version: 1.0
"""

import os
import sys
import argparse
import logging
import time
import csv
from pathlib import Path
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass
import json
import requests
from urllib.parse import urlencode
import xml.etree.ElementTree as ET
from Bio import Entrez, SeqIO
import re

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


@dataclass
class GenomeDownloaderConfig:
    """Configuration for GenomeDownloader"""
    accession_file: Optional[str]
    output_dir: str
    email: str = "vihaankulkarni29@gmail.com"
    api_key: Optional[str] = None
    batch_size: int = 50  # Smaller batches for genome downloads
    retry_attempts: int = 3
    timeout: int = 60  # Longer timeout for large files
    delay_between_requests: float = 1.0  # Respect NCBI rate limits


class GenomeDownloader:
    """
    Downloads genome sequences and annotations from NCBI using protein accession numbers.
    Similar to NCBIProteinHarvester but downloads genome data instead of protein sequences.
    """

    def __init__(self, config: GenomeDownloaderConfig):
        """Initialize the genome downloader with configuration"""
        self.config = config
        self.setup_logging()
        self.logger.info(f"Attempting to connect to NCBI database with email: {self.config.email}")
        if self.config.api_key:
            self.logger.info(f"Using NCBI API key: {self.config.api_key[:8]}...")
        else:
            self.logger.info("No NCBI API key provided - using default rate limits")
        self.setup_ncbi()
        self.setup_directories()

        # Statistics
        self.stats = {
            'total_accessions': 0,
            'genomes_downloaded': 0,
            'metadata_extracted': 0,
            'fasta_files_created': 0,
            'gff_files_created': 0,
            'failed_downloads': 0
        }

    def setup_logging(self):
        """Setup logging configuration"""
        log_dir = Path(self.config.output_dir) / "logs"
        log_dir.mkdir(parents=True, exist_ok=True)

        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_dir / "genome_download.log"),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger('GenomeDownloader')

    def setup_ncbi(self):
        """Setup NCBI Entrez configuration"""
        Entrez.email = self.config.email
        if self.config.api_key:
            Entrez.api_key = self.config.api_key
            self.logger.info(f"Using NCBI API key: {self.config.api_key[:8]}...")
        else:
            self.logger.warning("No NCBI API key provided - using default rate limits")

    def setup_directories(self):
        """Create necessary output directories"""
        output_path = Path(self.config.output_dir)
        fasta_path = output_path / "fasta"
        gff_path = output_path / "gff"

        output_path.mkdir(parents=True, exist_ok=True)
        fasta_path.mkdir(parents=True, exist_ok=True)
        gff_path.mkdir(parents=True, exist_ok=True)

        self.logger.info(f"Output directory: {output_path.absolute()}")
        self.logger.info(f"FASTA directory: {fasta_path.absolute()}")
        self.logger.info(f"GFF directory: {gff_path.absolute()}")

    def parse_accession_file(self) -> List[str]:
        """
        Parse accession file to get protein IDs for genome download

        Returns:
            List of protein accession IDs
        """
        self.logger.info(f"Parsing accession file: {self.config.accession_file}")

        # Input validation
        if not os.path.exists(self.config.accession_file):
            raise FileNotFoundError(f"Accession file not found: {self.config.accession_file}")

        if os.path.getsize(self.config.accession_file) == 0:
            raise ValueError(f"Accession file is empty: {self.config.accession_file}")

        protein_ids = []

        try:
            with open(self.config.accession_file, 'r', encoding='utf-8') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if line and not line.startswith('#'):  # Skip empty lines and comments
                        # Clean the accession number
                        accession = line.split()[0]  # Take first word if multiple columns
                        if self._is_valid_protein_id(accession):
                            protein_ids.append(accession)
                        else:
                            self.logger.warning(f"Skipping invalid accession at line {line_num}: {accession}")

        except Exception as e:
            self.logger.error(f"Error reading accession file: {e}")
            raise

        # Remove duplicates
        original_count = len(protein_ids)
        protein_ids = list(set(protein_ids))
        if len(protein_ids) != original_count:
            self.logger.info(f"Removed {original_count - len(protein_ids)} duplicate accessions")

        # Apply batch size limit
        if len(protein_ids) > self.config.batch_size:
            limited_ids = protein_ids[:self.config.batch_size]
            self.logger.info(f"Limited to first {self.config.batch_size} proteins for processing")
            protein_ids = limited_ids

        self.logger.info(f"Successfully parsed {len(protein_ids)} accession numbers")
        if protein_ids:
            self.logger.info(f"Sample accessions: {protein_ids[:5]}")

        self.stats['total_accessions'] = len(protein_ids)
        return protein_ids

    def _is_valid_protein_id(self, protein_id):
        """Validate protein ID format"""
        if not protein_id or not isinstance(protein_id, str):
            return False

        protein_id = protein_id.strip()

        # Check for common protein ID patterns
        patterns = [
            r'^[A-Z]{1,2}_\d+$',  # NCBI accession format (NP_, XP_, etc.)
            r'^\d+$',             # GI number format
            r'^[A-Z]{2}\d{6,}$', # GenBank format
            r'^[A-Z]{3}\d{5,}$', # Another accession format
        ]

        for pattern in patterns:
            if re.match(pattern, protein_id):
                return True

        # Allow alphanumeric IDs that are reasonable length
        if 3 <= len(protein_id) <= 20 and protein_id.replace('_', '').replace('.', '').isalnum():
            return True

        return False

    def download_genomes_batch(self, protein_ids: List[str]) -> List[Dict]:
        """
        Download genome sequences and annotations for a batch of protein IDs

        Args:
            protein_ids: List of protein accession IDs

        Returns:
            List of genome metadata dictionaries
        """
        if not protein_ids:
            self.logger.warning("No protein IDs provided for genome download")
            return []

        genome_metadata = []
        self.logger.info(f"Starting genome download for {len(protein_ids)} proteins")

        for attempt in range(self.config.retry_attempts):
            try:
                self.logger.info(f"Genome download attempt {attempt + 1}/{self.config.retry_attempts}")

                # Process each protein ID to find and download its genome
                successful_downloads = 0

                for protein_id in protein_ids:
                    try:
                        self.logger.debug(f"Processing protein: {protein_id}")

                        # Find the genome assembly for this protein
                        genome_info = self._find_genome_for_protein(protein_id)
                        if not genome_info:
                            self.logger.warning(f"Could not find genome for protein {protein_id}")
                            continue

                        # Download genome FASTA and GFF
                        fasta_success = self._download_genome_fasta(genome_info)
                        gff_success = self._download_genome_gff(genome_info)

                        if fasta_success:
                            self.stats['genomes_downloaded'] += 1
                            self.stats['fasta_files_created'] += 1

                            # Create metadata record
                            metadata_record = {
                                'protein_accession': protein_id,
                                'genome_accession': genome_info.get('assembly_accession', ''),
                                'organism': genome_info.get('organism', 'Unknown'),
                                'assembly_level': genome_info.get('assembly_level', ''),
                                'genome_size': genome_info.get('genome_size', ''),
                                'fasta_file': str(self._get_fasta_path(protein_id)),
                                'gff_file': str(self._get_gff_path(protein_id)) if gff_success else '',
                                'biosample': genome_info.get('biosample', ''),
                                'bioproject': genome_info.get('bioproject', ''),
                                'date': genome_info.get('date', '')
                            }

                            genome_metadata.append(metadata_record)
                            successful_downloads += 1

                            if gff_success:
                                self.stats['gff_files_created'] += 1

                        else:
                            self.stats['failed_downloads'] += 1

                    except Exception as protein_error:
                        self.logger.warning(f"Failed to process protein {protein_id}: {protein_error}")
                        self.stats['failed_downloads'] += 1
                        continue

                self.logger.info(f"Successfully downloaded {successful_downloads}/{len(protein_ids)} genomes")
                break

            except Exception as e:
                error_msg = str(e).lower()
                self.logger.warning(f"Genome download attempt {attempt + 1} failed: {e}")

                if "timeout" in error_msg:
                    self.logger.warning("Network timeout detected")
                elif "rate limit" in error_msg:
                    self.logger.warning("Rate limit exceeded")
                elif "network" in error_msg:
                    self.logger.warning("Network connectivity issue")

                if attempt < self.config.retry_attempts - 1:
                    wait_time = (attempt + 1) * 5
                    self.logger.info(f"Waiting {wait_time} seconds before retry...")
                    time.sleep(wait_time)
                else:
                    self.logger.error(f"Failed to download genomes after {self.config.retry_attempts} attempts")

        # Add delay between requests
        if self.config.delay_between_requests > 0:
            self.logger.debug(f"Sleeping {self.config.delay_between_requests}s to respect rate limits")
            time.sleep(self.config.delay_between_requests)

        return genome_metadata

    def _find_genome_for_protein(self, accession: str) -> Optional[Dict]:
        """
        Find genome assembly information for a given accession.
        This handles both protein accessions and contig/scaffold accessions.

        Args:
            accession: NCBI accession (protein or nucleotide)

        Returns:
            Genome assembly information dictionary or None if not found
        """
        try:
            # Check if this is a contig/scaffold accession (common pattern)
            if accession.count('.') >= 2 or '0' in accession:
                # This looks like a contig/scaffold accession
                self.logger.debug(f"Treating {accession} as contig/scaffold accession")
                return self._handle_contig_accession(accession)
            else:
                # Try as protein accession
                self.logger.debug(f"Treating {accession} as protein accession")
                return self._handle_protein_accession(accession)

        except Exception as e:
            self.logger.warning(f"Failed to find genome for accession {accession}: {e}")

        return None

    def _handle_contig_accession(self, contig_accession: str) -> Optional[Dict]:
        """Handle contig/scaffold accessions by finding their parent assembly"""
        try:
            # Get contig information from nucleotide database
            handle = Entrez.efetch(db="nucleotide", id=contig_accession, rettype="gb", retmode="text")
            record_text = handle.read()
            handle.close()

            # Parse the record to find assembly information
            assembly_info = self._parse_contig_record_for_assembly(record_text, contig_accession)

            if assembly_info:
                return assembly_info
            else:
                # Try to find assembly by searching
                return self._search_assembly_for_contig(contig_accession)

        except Exception as e:
            self.logger.warning(f"Failed to handle contig accession {contig_accession}: {e}")
            return None

    def _handle_protein_accession(self, protein_accession: str) -> Optional[Dict]:
        """Handle protein accessions (original logic)"""
        try:
            # Get protein information to find the source genome
            handle = Entrez.efetch(db="protein", id=protein_accession, rettype="gb", retmode="text")
            record_text = handle.read()
            handle.close()

            # Parse the record to find assembly information
            assembly_info = self._parse_protein_record_for_assembly(record_text, protein_accession)

            if assembly_info:
                return assembly_info
            else:
                # Fallback: try to find assembly by searching
                return self._search_assembly_for_protein(protein_accession)

        except Exception as e:
            self.logger.warning(f"Failed to handle protein accession {protein_accession}: {e}")
            return None

    def _parse_protein_record_for_assembly(self, record_text: str, protein_accession: str) -> Optional[Dict]:
        """Parse protein record to extract assembly information"""
        try:
            # Look for assembly accession in the record
            lines = record_text.split('\n')

            assembly_accession = None
            organism = "Unknown"
            biosample = ""
            bioproject = ""

            for line in lines:
                line = line.strip()

                # Look for assembly accession
                if '/assembly="' in line or 'Assembly:' in line:
                    # Extract assembly accession
                    if '/assembly="' in line:
                        assembly_match = re.search(r'/assembly="([^"]+)"', line)
                        if assembly_match:
                            assembly_accession = assembly_match.group(1)
                    elif 'Assembly:' in line:
                        assembly_match = re.search(r'Assembly:\s*([^\s]+)', line)
                        if assembly_match:
                            assembly_accession = assembly_match.group(1)

                # Extract organism
                if '/organism="' in line:
                    organism_match = re.search(r'/organism="([^"]+)"', line)
                    if organism_match:
                        organism = organism_match.group(1)

                # Extract BioSample
                if '/db_xref="BioSample:' in line:
                    biosample_match = re.search(r'/db_xref="BioSample:([^"]+)"', line)
                    if biosample_match:
                        biosample = biosample_match.group(1)

                # Extract BioProject
                if '/db_xref="BioProject:' in line:
                    bioproject_match = re.search(r'/db_xref="BioProject:([^"]+)"', line)
                    if bioproject_match:
                        bioproject = bioproject_match.group(1)

            if assembly_accession:
                # Get additional assembly information
                return self._get_assembly_details(assembly_accession, organism, biosample, bioproject)
            else:
                # Try to find assembly from BioProject
                if bioproject:
                    return self._find_assembly_from_bioproject(bioproject, organism)

        except Exception as e:
            self.logger.debug(f"Failed to parse protein record: {e}")

        return None

    def _get_assembly_details(self, assembly_accession: str, organism: str, biosample: str, bioproject: str) -> Dict:
        """Get detailed assembly information"""
        try:
            # Query NCBI for assembly details
            handle = Entrez.esummary(db="assembly", id=assembly_accession)
            summary = Entrez.read(handle)
            handle.close()

            if summary['DocumentSummarySet']['DocumentSummary']:
                assembly = summary['DocumentSummarySet']['DocumentSummary'][0]

                return {
                    'assembly_accession': assembly_accession,
                    'organism': organism,
                    'assembly_level': assembly.get('AssemblyStatus', ''),
                    'genome_size': assembly.get('GenomeSize', ''),
                    'ftp_path': assembly.get('FtpPath_RefSeq', assembly.get('FtpPath_GenBank', '')),
                    'biosample': biosample,
                    'bioproject': bioproject,
                    'date': assembly.get('SubmissionDate', ''),
                    'assembly_name': assembly.get('AssemblyName', '')
                }

        except Exception as e:
            self.logger.debug(f"Failed to get assembly details for {assembly_accession}: {e}")

        # Return basic info if detailed query fails
        return {
            'assembly_accession': assembly_accession,
            'organism': organism,
            'assembly_level': 'Unknown',
            'genome_size': '',
            'ftp_path': '',
            'biosample': biosample,
            'bioproject': bioproject,
            'date': '',
            'assembly_name': ''
        }

    def _find_assembly_from_bioproject(self, bioproject: str, organism: str) -> Optional[Dict]:
        """Find assembly from BioProject ID"""
        try:
            handle = Entrez.esearch(db="assembly", term=f"bioproject:{bioproject}", retmax=1)
            record = Entrez.read(handle)
            handle.close()

            if record['IdList']:
                assembly_id = record['IdList'][0]
                handle = Entrez.esummary(db="assembly", id=assembly_id)
                summary = Entrez.read(handle)
                handle.close()

                assembly = summary['DocumentSummarySet']['DocumentSummary'][0]
                assembly_accession = assembly['AssemblyAccession']

                return {
                    'assembly_accession': assembly_accession,
                    'organism': organism,
                    'assembly_level': assembly.get('AssemblyStatus', ''),
                    'genome_size': assembly.get('GenomeSize', ''),
                    'ftp_path': assembly.get('FtpPath_RefSeq', assembly.get('FtpPath_GenBank', '')),
                    'biosample': '',
                    'bioproject': bioproject,
                    'date': assembly.get('SubmissionDate', ''),
                    'assembly_name': assembly.get('AssemblyName', '')
                }

        except Exception as e:
            self.logger.debug(f"Failed to find assembly from BioProject {bioproject}: {e}")

        return None

    def _search_assembly_for_protein(self, protein_accession: str) -> Optional[Dict]:
        """Search for assembly containing the protein"""
        try:
            # Search assemblies that contain this protein
            handle = Entrez.esearch(db="assembly", term=f"{protein_accession}[Protein]", retmax=1)
            record = Entrez.read(handle)
            handle.close()

            if record['IdList']:
                assembly_id = record['IdList'][0]
                handle = Entrez.esummary(db="assembly", id=assembly_id)
                summary = Entrez.read(handle)
                handle.close()

                assembly = summary['DocumentSummarySet']['DocumentSummary'][0]
                assembly_accession = assembly['AssemblyAccession']

                return {
                    'assembly_accession': assembly_accession,
                    'organism': assembly.get('SpeciesName', 'Unknown'),
                    'assembly_level': assembly.get('AssemblyStatus', ''),
                    'genome_size': assembly.get('GenomeSize', ''),
                    'ftp_path': assembly.get('FtpPath_RefSeq', assembly.get('FtpPath_GenBank', '')),
                    'biosample': assembly.get('BioSampleAccn', ''),
                    'bioproject': assembly.get('BioProjectAccn', ''),
                    'date': assembly.get('SubmissionDate', ''),
                    'assembly_name': assembly.get('AssemblyName', '')
                }

        except Exception as e:
            self.logger.debug(f"Failed to search assembly for protein {protein_accession}: {e}")

        return None

    def _parse_contig_record_for_assembly(self, record_text: str, contig_accession: str) -> Optional[Dict]:
        """Parse contig record to extract assembly information"""
        try:
            lines = record_text.split('\n')

            # Look for assembly information in the record
            for line in lines:
                line = line.strip()

                # Look for assembly accession in various formats
                if '/assembly="' in line or 'Assembly:' in line:
                    if '/assembly="' in line:
                        assembly_match = re.search(r'/assembly="([^"]+)"', line)
                        if assembly_match:
                            assembly_accession = assembly_match.group(1)
                    elif 'Assembly:' in line:
                        assembly_match = re.search(r'Assembly:\s*([^\s]+)', line)
                        if assembly_match:
                            assembly_accession = assembly_match.group(1)

                # Look for BioProject
                if '/db_xref="BioProject:' in line:
                    bioproject_match = re.search(r'/db_xref="BioProject:([^"]+)"', line)
                    if bioproject_match:
                        bioproject = bioproject_match.group(1)

                # Look for organism
                if '/organism="' in line:
                    organism_match = re.search(r'/organism="([^"]+)"', line)
                    if organism_match:
                        organism = organism_match.group(1)

            # If we found assembly info, get full details
            if 'assembly_accession' in locals() and assembly_accession:
                return self._get_assembly_details(assembly_accession,
                                                locals().get('organism', 'Unknown'),
                                                '',
                                                locals().get('bioproject', ''))

        except Exception as e:
            self.logger.debug(f"Failed to parse contig record: {e}")

        return None

    def _search_assembly_for_contig(self, contig_accession: str) -> Optional[Dict]:
        """Search for assembly containing the contig"""
        try:
            # Search assemblies that contain this contig
            handle = Entrez.esearch(db="assembly", term=f"{contig_accession}[Nucleotide]", retmax=1)
            record = Entrez.read(handle)
            handle.close()

            if record['IdList']:
                assembly_id = record['IdList'][0]
                handle = Entrez.esummary(db="assembly", id=assembly_id)
                summary = Entrez.read(handle)
                handle.close()

                assembly = summary['DocumentSummarySet']['DocumentSummary'][0]
                assembly_accession = assembly['AssemblyAccession']

                return {
                    'assembly_accession': assembly_accession,
                    'organism': assembly.get('SpeciesName', 'Unknown'),
                    'assembly_level': assembly.get('AssemblyStatus', ''),
                    'genome_size': assembly.get('GenomeSize', ''),
                    'ftp_path': assembly.get('FtpPath_RefSeq', assembly.get('FtpPath_GenBank', '')),
                    'biosample': assembly.get('BioSampleAccn', ''),
                    'bioproject': assembly.get('BioProjectAccn', ''),
                    'date': assembly.get('SubmissionDate', ''),
                    'assembly_name': assembly.get('AssemblyName', '')
                }

        except Exception as e:
            self.logger.debug(f"Failed to search assembly for contig {contig_accession}: {e}")

        return None

    def _download_genome_fasta(self, genome_info: Dict) -> bool:
        """Download genome FASTA file"""
        ftp_path = genome_info.get('ftp_path', '')
        if not ftp_path:
            self.logger.warning("No FTP path available for FASTA download")
            return False

    def run_download(self) -> bool:
        """
        Main method to run the genome downloading process

        Returns:
            True if successful, False otherwise
        """
        try:
            self.logger.info("Starting GenomeDownloader")
            self.logger.info(f"Accession file: {self.config.accession_file}")
            self.logger.info(f"Output directory: {self.config.output_dir}")
            self.logger.info(f"Batch size: {self.config.batch_size}")
            self.logger.info(f"Email: {self.config.email}")
            self.logger.info(f"API Key: {'Set' if self.config.api_key else 'Not set'}")

            # Parse accession file
            self.logger.info("Step 1: Parsing accession file")
            protein_ids = self.parse_accession_file()

            if not protein_ids:
                self.logger.error("No protein IDs found in accession file")
                return False

            # Process in batches
            all_genome_metadata = []
            batch_size = min(self.config.batch_size, 20)  # Smaller batches for genome downloads
            total_batches = (len(protein_ids) + batch_size - 1) // batch_size

            self.logger.info(f"Step 2: Processing {total_batches} batches of up to {batch_size} proteins each")

            for i in range(0, len(protein_ids), batch_size):
                batch_num = (i // batch_size) + 1
                batch_ids = protein_ids[i:i + batch_size]

                self.logger.info(f"Processing batch {batch_num}/{total_batches} ({len(batch_ids)} proteins)")

                # Download genomes for this batch
                batch_metadata = self.download_genomes_batch(batch_ids)
                if batch_metadata:
                    all_genome_metadata.extend(batch_metadata)
                    self.logger.info(f"Downloaded {len(batch_metadata)} genomes from batch {batch_num}")

            # Save results
            self.logger.info("Step 3: Saving results")
            success = self._save_results(all_genome_metadata)

            if not success:
                self.logger.error("Failed to save results")
                return False

            # Print final statistics
            self._print_statistics()

            success_rate = (self.stats['genomes_downloaded'] / self.stats['total_accessions'] * 100) if self.stats['total_accessions'] > 0 else 0
            self.logger.info("GENOME DOWNLOAD COMPLETED SUCCESSFULLY!")
            self.logger.info(f"Success rate: {success_rate:.1f}%")
            self.logger.info(f"Results saved to: {self.config.output_dir}")

            return True

        except KeyboardInterrupt:
            self.logger.warning("Genome download interrupted by user")
            return False
        except Exception as e:
            self.logger.error(f"Genome download failed with unexpected error: {e}")
            import traceback
            self.logger.error(f"Stack trace: {traceback.format_exc()}")
            return False

    def _save_results(self, genome_metadata: List[Dict]) -> bool:
        """Save genome metadata to CSV file"""
        try:
            if genome_metadata:
                # Save as CSV without pandas
                metadata_file = Path(self.config.output_dir) / "genome_metadata.csv"
                with open(metadata_file, 'w', newline='') as f:
                    if genome_metadata:
                        writer = csv.DictWriter(f, fieldnames=genome_metadata[0].keys())
                        writer.writeheader()
                        writer.writerows(genome_metadata)
                self.logger.info(f"Saved metadata for {len(genome_metadata)} genomes to {metadata_file}")

                # Save summary statistics
                summary_file = Path(self.config.output_dir) / "download_summary.json"
                with open(summary_file, 'w') as f:
                    json.dump(self.stats, f, indent=2)
                self.logger.info(f"Saved download summary to {summary_file}")

            return True

        except Exception as e:
            self.logger.error(f"Failed to save results: {e}")
            return False

    def _print_statistics(self):
        """Print download statistics"""
        self.logger.info("=" * 60)
        self.logger.info("GENOME DOWNLOAD STATISTICS")
        self.logger.info("=" * 60)
        self.logger.info(f"Total accessions: {self.stats['total_accessions']}")
        self.logger.info(f"Genomes downloaded: {self.stats['genomes_downloaded']}")
        self.logger.info(f"FASTA files created: {self.stats['fasta_files_created']}")
        self.logger.info(f"GFF files created: {self.stats['gff_files_created']}")
        self.logger.info(f"Failed downloads: {self.stats['failed_downloads']}")

        success_rate = (self.stats['genomes_downloaded'] / self.stats['total_accessions'] * 100) if self.stats['total_accessions'] > 0 else 0
        self.logger.info(f"Success rate: {success_rate:.1f}%")


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description="Download genome sequences and annotations from NCBI using protein accession numbers",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Download genomes from accession list
  python GenomeDownloader.py --accession-list data/accession_list.txt --email your.email@example.com

  # Download with API key for faster access
  python GenomeDownloader.py --accession-list data/accession_list.txt --email your.email@example.com --api-key YOUR_API_KEY

  # Download first 15 genomes only
  python GenomeDownloader.py --accession-list data/accession_list.txt --email your.email@example.com --batch-size 15

  # Specify custom output directory
  python GenomeDownloader.py --accession-list data/accession_list.txt --email your.email@example.com --output-dir genome_data
        """
    )

    parser.add_argument(
        '--accession-list',
        required=True,
        help='Path to file containing protein accession numbers (one per line)'
    )

    parser.add_argument(
        '--output-dir',
        required=True,
        help='Output directory for downloaded genomes (required for pipeline isolation)'
    )

    parser.add_argument(
        '--email',
        default="vihaankulkarni29@gmail.com",
        help="Email for NCBI Entrez (required by NCBI)"
    )

    parser.add_argument(
        '--api-key',
        help='NCBI API key for increased rate limits'
    )

    parser.add_argument(
        '--batch-size',
        type=int,
        default=50,
        help='Maximum number of genomes to download (default: 50)'
    )

    args = parser.parse_args()

    # Create configuration
    config = GenomeDownloaderConfig(
        accession_file=args.accession_list,
        output_dir=args.output_dir,
        email=args.email,
        api_key=args.api_key,
        batch_size=args.batch_size
    )

    # Run downloader
    downloader = GenomeDownloader(config)
    success = downloader.run_download()

    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()