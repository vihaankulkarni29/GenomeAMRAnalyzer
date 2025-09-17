import sys; import os; sys.path.append(os.path.join(os.path.dirname(__file__), '\\\\\\src')); from configuration_manager import config_manager
#!/usr/bin/env python3
"""
NCBIProteinHarvester.py - Download protein sequences and metadata from NCBI

This script downloads protein FASTA sequences and metadata for RND efflux pump proteins
from NCBI using protein accession numbers. It creates properly formatted FASTA files
and metadata CSV files ready for downstream analysis with WildTypeAligner.

Usage:
    python NCBIProteinHarvester.py --accession-list data/accession_list.txt --email your.email@example.com --output-dir protein_data

Author: MetaDataHarvester Pipeline
Version: 1.0
"""

import os
import sys
import argparse
import logging
import time
import re
from pathlib import Path
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass
import json
import pandas as pd
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


@dataclass
class ProteinHarvesterConfig:
    """Configuration for NCBIProteinHarvester"""
    accession_file: str
    output_dir: str
    email: str = "vihaankulkarni29@gmail.com"
    api_key: Optional[str] = None
    batch_size: int = 50
    retry_attempts: int = 3
    delay_between_requests: float = 0.5


class NCBIProteinHarvester:
    """
    Downloads protein sequences and metadata from NCBI using accession numbers.
    """

    def __init__(self, config: ProteinHarvesterConfig):
        """Initialize the protein harvester with configuration"""
        self.config = config
        self.setup_logging()
        self.setup_ncbi()
        self.setup_directories()

        # Statistics
        self.stats = {
            'total_accessions': 0,
            'proteins_downloaded': 0,
            'fasta_files_created': 0,
            'metadata_records': 0,
            'failed_downloads': 0
        }

        # RND efflux pump protein identifiers
        self.rnd_proteins = {
            config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0]: [config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0], 'AcrA', 'acrA_MG1655'],
            config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1]: [config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1], 'AcrB', 'acrB_MG1655'],
            config_manager.get_default_genes("rnd_efflux_pumps", "primary")[2]: [config_manager.get_default_genes("rnd_efflux_pumps", "primary")[2], 'TolC', 'tolC_MG1655'],
            'acrD': ['acrD', 'AcrD'],
            'acrE': ['acrE', 'AcrE'],
            'acrF': ['acrF', 'AcrF'],
            'mdtB': ['mdtB', 'MdtB'],
            'mdtC': ['mdtC', 'MdtC'],
            'oqxA': ['oqxA', 'OqxA'],
            'oqxB': ['oqxB', 'OqxB'],
            'acrR': ['acrR', 'AcrR'],
            'acrZ': ['acrZ', 'AcrZ'],
            'marA': ['marA', 'MarA'],
            'ramA': ['ramA', 'RamA'],
            'soxS': ['soxS', 'SoxS'],
            'rob': ['rob', 'Rob'],
            'eefA': ['eefA', 'EefA'],
            'eefB': ['eefB', 'EefB'],
            'eefC': ['eefC', 'EefC']
        }

    def setup_logging(self):
        """Setup logging configuration"""
        log_dir = Path(self.config.output_dir) / "logs"
        log_dir.mkdir(parents=True, exist_ok=True)

        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(name)s - %(message)s',
            handlers=[
                logging.FileHandler(log_dir / "protein_harvest.log"),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger('NCBIProteinHarvester')

    def setup_ncbi(self):
        """Setup NCBI Entrez configuration"""
        Entrez.email = self.config.email
        if self.config.api_key:
            Entrez.api_key = self.config.api_key
            self.logger.info(f"Using NCBI API key: {self.config.api_key[:8]}...")

    def setup_directories(self):
        """Create necessary output directories"""
        output_path = Path(self.config.output_dir)
        fasta_path = output_path / "fasta"
        metadata_path = output_path / "metadata"

        output_path.mkdir(parents=True, exist_ok=True)
        fasta_path.mkdir(parents=True, exist_ok=True)
        metadata_path.mkdir(parents=True, exist_ok=True)

        self.logger.info(f"Output directory: {output_path.absolute()}")
        self.logger.info(f"FASTA directory: {fasta_path.absolute()}")
        self.logger.info(f"Metadata directory: {metadata_path.absolute()}")

    def parse_accession_file(self) -> List[str]:
        """
        Parse accession file to get protein IDs

        Returns:
            List of protein accession IDs
        """
        self.logger.info(f"Parsing accession file: {self.config.accession_file}")

        if not os.path.exists(self.config.accession_file):
            raise FileNotFoundError(f"Accession file not found: {self.config.accession_file}")

        protein_ids = []

        try:
            with open(self.config.accession_file, 'r', encoding='utf-8') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if line and not line.startswith('#'):
                        accession = line.split()[0]
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

        self.logger.info(f"Successfully parsed {len(protein_ids)} accession numbers")
        if protein_ids:
            self.logger.info(f"Sample accessions: {protein_ids[:5]}")

        self.stats['total_accessions'] = len(protein_ids)
        return protein_ids

    def _is_valid_protein_id(self, protein_id: str) -> bool:
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
            r'^[A-Z]\d[A-Z]\d[A-Z]\d+$',  # UniProt format (P0AE06, Q8X8Z1, etc.)
        ]

        for pattern in patterns:
            if re.match(pattern, protein_id):
                return True

        # Allow alphanumeric IDs that are reasonable length
        if 3 <= len(protein_id) <= 20 and protein_id.replace('_', '').replace('.', '').isalnum():
            return True

        return False

    def download_proteins_batch(self, protein_ids: List[str]) -> List[Dict]:
        """
        Download protein sequences and metadata for a batch of accession IDs

        Args:
            protein_ids: List of protein accession IDs

        Returns:
            List of protein metadata dictionaries
        """
        if not protein_ids:
            self.logger.warning("No protein IDs provided for download")
            return []

        self.logger.info(f"Starting protein download for {len(protein_ids)} accessions")

        protein_metadata = []
        successful_downloads = 0

        for protein_id in protein_ids:
            try:
                self.logger.debug(f"Processing protein: {protein_id}")

                # Download protein sequence and metadata
                protein_info = self._download_single_protein(protein_id)

                if protein_info:
                    protein_metadata.append(protein_info)
                    successful_downloads += 1
                    self.stats['proteins_downloaded'] += 1

                    if protein_info.get('fasta_file'):
                        self.stats['fasta_files_created'] += 1

                else:
                    self.stats['failed_downloads'] += 1

            except Exception as protein_error:
                self.logger.warning(f"Failed to process protein {protein_id}: {protein_error}")
                self.stats['failed_downloads'] += 1
                continue

            # Respect NCBI rate limits
            if self.config.delay_between_requests > 0:
                time.sleep(self.config.delay_between_requests)

        self.logger.info(f"Successfully downloaded {successful_downloads}/{len(protein_ids)} proteins")
        return protein_metadata

    def _download_single_protein(self, accession: str) -> Optional[Dict]:
        """
        Download a single protein sequence and metadata

        Args:
            accession: Protein accession ID

        Returns:
            Protein metadata dictionary or None if failed
        """
        try:
            # Check if this is a UniProt accession (starts with a letter)
            if accession[0].isalpha():
                # Handle UniProt accession
                return self._download_uniprot_protein(accession)
            else:
                # Handle NCBI accession
                return self._download_ncbi_protein(accession)

        except Exception as e:
            self.logger.warning(f"Failed to download protein {accession}: {e}")
            return None

    def _download_ncbi_protein(self, accession: str) -> Optional[Dict]:
        """Download protein from NCBI"""
        try:
            # Fetch protein record from NCBI
            handle = Entrez.efetch(db="protein", id=accession, rettype="gb", retmode="text")
            record_text = handle.read()
            handle.close()

            # Parse the record
            protein_info = self._parse_protein_record(record_text, accession)

            if not protein_info:
                return None

            # Download FASTA sequence
            fasta_success = self._download_protein_fasta(accession, protein_info)

            if fasta_success:
                self.logger.info(f"Successfully downloaded NCBI protein: {accession}")
                return protein_info
            else:
                self.logger.warning(f"Failed to download FASTA for {accession}")
                return None

        except Exception as e:
            self.logger.warning(f"Failed to download NCBI protein {accession}: {e}")
            return None

    def _download_uniprot_protein(self, accession: str) -> Optional[Dict]:
        """Download protein from UniProt"""
        try:
            # Use UniProt REST API
            import requests

            # Get protein information from UniProt
            url = f"https://www.uniprot.org/uniprot/{accession}.txt"
            response = requests.get(url, timeout=30)

            if response.status_code != 200:
                self.logger.warning(f"UniProt returned status {response.status_code} for {accession}")
                return None

            record_text = response.text

            # Parse UniProt record
            protein_info = self._parse_uniprot_record(record_text, accession)

            if not protein_info:
                return None

            # Download FASTA sequence
            fasta_success = self._download_uniprot_fasta(accession, protein_info)

            if fasta_success:
                self.logger.info(f"Successfully downloaded UniProt protein: {accession}")
                return protein_info
            else:
                self.logger.warning(f"Failed to download FASTA for {accession}")
                return None

        except Exception as e:
            self.logger.warning(f"Failed to download UniProt protein {accession}: {e}")
            return None

    def _parse_protein_record(self, record_text: str, accession: str) -> Optional[Dict]:
        """
        Parse NCBI protein record to extract metadata

        Args:
            record_text: Raw NCBI record text
            accession: Protein accession

        Returns:
            Protein metadata dictionary
        """
        try:
            lines = record_text.split('\n')

            protein_info = {
                'Accession_Number': accession,
                'organism': 'Unknown',
                'protein_name': 'Unknown',
                'protein_family': 'Unknown',
                'gene_name': '',
                'locus_tag': '',
                'sequence_length': 0,
                'biosample': '',
                'bioproject': '',
                'taxonomy': '',
                'definition': '',
                'fasta_file': '',
                'metadata_file': ''
            }

            in_features = False
            current_feature = None

            for line in lines:
                line = line.strip()

                # Parse header information
                if line.startswith('DEFINITION'):
                    protein_info['definition'] = line.replace('DEFINITION', '').strip()

                elif line.startswith('ACCESSION'):
                    # Already have accession
                    pass

                elif line.startswith('VERSION'):
                    # Could extract version if needed
                    pass

                elif line.startswith('ORGANISM'):
                    organism = line.replace('ORGANISM', '').strip()
                    protein_info['organism'] = organism

                elif line.startswith('SOURCE'):
                    # Additional source information
                    pass

                elif '/db_xref="BioSample:' in line:
                    biosample_match = re.search(r'/db_xref="BioSample:([^"]+)"', line)
                    if biosample_match:
                        protein_info['biosample'] = biosample_match.group(1)

                elif '/db_xref="BioProject:' in line:
                    bioproject_match = re.search(r'/db_xref="BioProject:([^"]+)"', line)
                    if bioproject_match:
                        protein_info['bioproject'] = bioproject_match.group(1)

                elif '/gene=' in line:
                    gene_match = re.search(r'/gene="([^"]+)"', line)
                    if gene_match:
                        protein_info['gene_name'] = gene_match.group(1)

                elif '/locus_tag=' in line:
                    locus_match = re.search(r'/locus_tag="([^"]+)"', line)
                    if locus_match:
                        protein_info['locus_tag'] = locus_match.group(1)

                elif '/product=' in line:
                    product_match = re.search(r'/product="([^"]+)"', line)
                    if product_match:
                        protein_info['protein_name'] = product_match.group(1)

            # Determine protein family based on gene name and product
            protein_info['protein_family'] = self._determine_protein_family(protein_info)

            return protein_info

        except Exception as e:
            self.logger.warning(f"Failed to parse protein record for {accession}: {e}")
            return None

    def _determine_protein_family(self, protein_info: Dict) -> str:
        """Determine protein family based on gene name and product"""
        gene_name = protein_info.get('gene_name', '').lower()
        protein_name = protein_info.get('protein_name', '').lower()
        locus_tag = protein_info.get('locus_tag', '').lower()

        # Check against known RND proteins
        for family, identifiers in self.rnd_proteins.items():
            for identifier in identifiers:
                if (identifier.lower() in gene_name or
                    identifier.lower() in protein_name or
                    identifier.lower() in locus_tag):
                    return family

        # If no match found, try to infer from keywords
        if 'efflux' in protein_name or 'pump' in protein_name:
            if 'rnd' in protein_name or 'resistance' in protein_name:
                return 'rnd_efflux'
            else:
                return 'efflux_pump'
        elif 'regulator' in protein_name or 'transcriptional' in protein_name:
            return 'regulator'
        else:
            return 'unknown'

    def _download_protein_fasta(self, accession: str, protein_info: Dict) -> bool:
        """
        Download protein FASTA sequence

        Args:
            accession: Protein accession
            protein_info: Protein metadata dictionary

        Returns:
            True if successful, False otherwise
        """
        try:
            # Fetch FASTA sequence
            handle = Entrez.efetch(db="protein", id=accession, rettype="fasta", retmode="text")
            fasta_record = SeqIO.read(handle, "fasta")
            handle.close()

            # Create output filename
            protein_family = protein_info.get('protein_family', 'unknown')
            output_filename = f"{accession}_{protein_family}.fasta"
            output_path = Path(self.config.output_dir) / "fasta" / output_filename

            # Update protein info with sequence details
            protein_info['sequence_length'] = len(fasta_record.seq)
            protein_info['fasta_file'] = str(output_path)

            # Update record description with our metadata
            fasta_record.id = f"{accession}_{protein_family}"
            fasta_record.name = protein_family
            fasta_record.description = f"{protein_info.get('protein_name', 'Unknown protein')} from {protein_info.get('organism', 'Unknown organism')}"

            # Write FASTA file
            with open(output_path, 'w') as f:
                SeqIO.write(fasta_record, f, 'fasta')

            self.logger.debug(f"Downloaded FASTA: {accession} -> {output_path}")
            return True

        except Exception as e:
            self.logger.warning(f"Failed to download FASTA for {accession}: {e}")
            return False

    def _parse_uniprot_record(self, record_text: str, accession: str) -> Optional[Dict]:
        """
        Parse UniProt record to extract metadata

        Args:
            record_text: Raw UniProt record text
            accession: Protein accession

        Returns:
            Protein metadata dictionary
        """
        try:
            lines = record_text.split('\n')

            protein_info = {
                'Accession_Number': accession,
                'organism': 'Unknown',
                'protein_name': 'Unknown',
                'protein_family': 'Unknown',
                'gene_name': '',
                'locus_tag': '',
                'sequence_length': 0,
                'biosample': '',
                'bioproject': '',
                'taxonomy': '',
                'definition': '',
                'fasta_file': '',
                'metadata_file': ''
            }

            for line in lines:
                line = line.strip()

                # Parse UniProt format
                if line.startswith('ID'):
                    # ID line contains basic info
                    pass

                elif line.startswith('AC'):
                    # Accession numbers
                    pass

                elif line.startswith('DE'):
                    # Description/definition
                    if 'RecName:' in line:
                        # Extract protein name
                        name_match = re.search(r'RecName:\s*Full=(.+?);', line)
                        if name_match:
                            protein_info['protein_name'] = name_match.group(1).strip()

                elif line.startswith('GN'):
                    # Gene name
                    if 'Name=' in line:
                        gene_match = re.search(r'Name=(\w+)', line)
                        if gene_match:
                            protein_info['gene_name'] = gene_match.group(1)

                elif line.startswith('OS'):
                    # Organism
                    organism_match = re.search(r'OS\s+(.+)', line)
                    if organism_match:
                        protein_info['organism'] = organism_match.group(1).strip()

                elif line.startswith('OX'):
                    # Taxonomy ID
                    tax_match = re.search(r'NCBI_TaxID=(\d+)', line)
                    if tax_match:
                        protein_info['taxonomy'] = tax_match.group(1)

                elif line.startswith('DR'):
                    # Database cross-references
                    if 'BioSample;' in line:
                        bs_match = re.search(r'BioSample;\s*([^;]+)', line)
                        if bs_match:
                            protein_info['biosample'] = bs_match.group(1).strip()

            # Determine protein family based on gene name and product
            protein_info['protein_family'] = self._determine_protein_family(protein_info)

            return protein_info

        except Exception as e:
            self.logger.warning(f"Failed to parse UniProt record for {accession}: {e}")
            return None

    def _download_uniprot_fasta(self, accession: str, protein_info: Dict) -> bool:
        """
        Download protein FASTA sequence from UniProt

        Args:
            accession: Protein accession
            protein_info: Protein metadata dictionary

        Returns:
            True if successful, False otherwise
        """
        try:
            import requests

            # Get FASTA sequence from UniProt
            url = f"https://www.uniprot.org/uniprot/{accession}.fasta"
            response = requests.get(url, timeout=30)

            if response.status_code != 200:
                self.logger.warning(f"UniProt FASTA returned status {response.status_code} for {accession}")
                return False

            # Parse FASTA
            fasta_content = response.text
            lines = fasta_content.split('\n')

            if len(lines) < 2:
                self.logger.warning(f"Invalid FASTA format for {accession}")
                return False

            # Extract sequence
            sequence = ''.join(lines[1:]).replace('\n', '').replace(' ', '')

            # Create SeqRecord
            from Bio.SeqRecord import SeqRecord
            from Bio.Seq import Seq

            protein_family = protein_info.get('protein_family', 'unknown')
            record_name = protein_family
            record_description = f"{protein_info.get('protein_name', 'Unknown protein')} from {protein_info.get('organism', 'Unknown organism')}"

            seq_record = SeqRecord(
                Seq(sequence),
                id=accession,  # Use just accession for WildTypeAligner compatibility
                name=record_name,
                description=record_description
            )

            # Update protein info with sequence details
            protein_info['sequence_length'] = len(sequence)
            protein_info['fasta_file'] = str(self._get_fasta_path(accession))

            # Save FASTA file
            output_path = Path(self.config.output_dir) / "fasta" / f"{accession}_{protein_family}.fasta"
            with open(output_path, 'w') as f:
                from Bio import SeqIO
                SeqIO.write(seq_record, f, 'fasta')

            self.logger.debug(f"Downloaded UniProt FASTA: {accession} -> {output_path}")
            return True

        except Exception as e:
            self.logger.warning(f"Failed to download UniProt FASTA for {accession}: {e}")
            return False

    def _get_fasta_path(self, accession: str) -> Path:
        """Get the path for a FASTA file"""
        protein_family = 'unknown'  # This would be determined from protein_info
        return Path(self.config.output_dir) / "fasta" / f"{accession}_{protein_family}.fasta"

    def run_harvest(self) -> bool:
        """
        Main method to run the protein harvesting process

        Returns:
            True if successful, False otherwise
        """
        try:
            self.logger.info("Starting NCBI Protein Harvesting")
            self.logger.info(f"Accession file: {self.config.accession_file}")
            self.logger.info(f"Output directory: {self.config.output_dir}")
            self.logger.info(f"Batch size: {self.config.batch_size}")
            self.logger.info(f"Email: {self.config.email}")

            # Parse accession file
            self.logger.info("Step 1: Parsing accession file")
            protein_ids = self.parse_accession_file()

            if not protein_ids:
                self.logger.error("No protein IDs found in accession file")
                return False

            # Process in batches
            all_protein_metadata = []
            batch_size = min(self.config.batch_size, 20)  # Smaller batches for protein downloads
            total_batches = (len(protein_ids) + batch_size - 1) // batch_size

            self.logger.info(f"Step 2: Processing {total_batches} batches of up to {batch_size} proteins each")

            for i in range(0, len(protein_ids), batch_size):
                batch_num = (i // batch_size) + 1
                batch_ids = protein_ids[i:i + batch_size]

                self.logger.info(f"Processing batch {batch_num}/{total_batches} ({len(batch_ids)} proteins)")

                # Download proteins for this batch
                batch_metadata = self.download_proteins_batch(batch_ids)
                if batch_metadata:
                    all_protein_metadata.extend(batch_metadata)
                    self.logger.info(f"Downloaded {len(batch_metadata)} proteins from batch {batch_num}")

            # Save results
            self.logger.info("Step 3: Saving results")
            success = self._save_results(all_protein_metadata)

            if not success:
                self.logger.error("Failed to save results")
                return False

            # Print final statistics
            self._print_statistics()

            success_rate = (self.stats['proteins_downloaded'] / self.stats['total_accessions'] * 100) if self.stats['total_accessions'] > 0 else 0
            self.logger.info("PROTEIN HARVESTING COMPLETED SUCCESSFULLY!")
            self.logger.info(f"Success rate: {success_rate:.1f}%")
            self.logger.info(f"Results saved to: {self.config.output_dir}")

            return True

        except KeyboardInterrupt:
            self.logger.warning("Protein harvesting interrupted by user")
            return False
        except Exception as e:
            self.logger.error(f"Protein harvesting failed with unexpected error: {e}")
            import traceback
            self.logger.error(f"Stack trace: {traceback.format_exc()}")
            return False

    def _save_results(self, protein_metadata: List[Dict]) -> bool:
        """Save protein metadata to CSV file"""
        try:
            if protein_metadata:
                df = pd.DataFrame(protein_metadata)
                metadata_file = Path(self.config.output_dir) / "metadata" / "protein_metadata.csv"
                df.to_csv(metadata_file, index=False)
                self.logger.info(f"Saved metadata for {len(protein_metadata)} proteins to {metadata_file}")

                # Save summary statistics
                summary_file = Path(self.config.output_dir) / "metadata" / "harvest_summary.json"
                with open(summary_file, 'w') as f:
                    json.dump(self.stats, f, indent=2)
                self.logger.info(f"Saved harvest summary to {summary_file}")

            return True

        except Exception as e:
            self.logger.error(f"Failed to save results: {e}")
            return False

    def _print_statistics(self):
        """Print harvesting statistics"""
        self.logger.info("=" * 60)
        self.logger.info("PROTEIN HARVESTING STATISTICS")
        self.logger.info("=" * 60)
        self.logger.info(f"Total accessions: {self.stats['total_accessions']}")
        self.logger.info(f"Proteins downloaded: {self.stats['proteins_downloaded']}")
        self.logger.info(f"FASTA files created: {self.stats['fasta_files_created']}")
        self.logger.info(f"Failed downloads: {self.stats['failed_downloads']}")

        success_rate = (self.stats['proteins_downloaded'] / self.stats['total_accessions'] * 100) if self.stats['total_accessions'] > 0 else 0
        self.logger.info(f"Success rate: {success_rate:.1f}%")
        self.logger.info("=" * 60)


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description="Download protein sequences and metadata from NCBI using accession numbers",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Download proteins from accession list
  python NCBIProteinHarvester.py --accession-list data/accession_list.txt --email your.email@example.com

  # Download with API key for faster access
  python NCBIProteinHarvester.py --accession-list data/accession_list.txt --email your.email@example.com --api-key YOUR_API_KEY

  # Download first 25 proteins only
  python NCBIProteinHarvester.py --accession-list data/accession_list.txt --email your.email@example.com --batch-size 25

  # Specify custom output directory
  python NCBIProteinHarvester.py --accession-list data/accession_list.txt --email your.email@example.com --output-dir protein_data

The tool will:
  1. Read protein accession numbers from the input file
  2. Download protein sequences and metadata from NCBI
  3. Classify proteins into RND efflux pump families
  4. Save FASTA files and metadata CSV ready for WildTypeAligner
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
        help='Output directory for downloaded proteins and metadata'
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
        help='Maximum number of proteins to download (default: 50)'
    )

    args = parser.parse_args()

    # Create configuration
    config = ProteinHarvesterConfig(
        accession_file=args.accession_list,
        output_dir=args.output_dir,
        email=args.email,
        api_key=args.api_key,
        batch_size=args.batch_size
    )

    # Run harvester
    harvester = NCBIProteinHarvester(config)
    success = harvester.run_harvest()

    if success:
        print("\n" + "="*80)
        print("PROTEIN HARVESTING COMPLETED SUCCESSFULLY!")
        print("="*80)
        print(f"Output directory: {args.output_dir}")
        print(f"FASTA files: {args.output_dir}/fasta/")
        print(f"Metadata: {args.output_dir}/metadata/protein_metadata.csv")
        print("="*80)
        print("Ready for WildTypeAligner -> SubScan -> HTML Report Generator")
        print("="*80)

    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()