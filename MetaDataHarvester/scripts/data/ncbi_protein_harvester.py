import sys; import os; sys.path.append(os.path.join(os.path.dirname(__file__), '\\\\\\\\\src')); from configuration_manager import config_manager
#!/usr/bin/env python3
"""
NCBIProteinHarvester - Streamlined protein sequence and metadata extraction

This module processes NCBI protein accession lists to extract:
- Protein sequences in FASTA format
- Comprehensive metadata (organism, protein name, family, etc.)

Usage:
    python ncbi_protein_harvester.py --accession-file accessions.txt --output-dir results/

Author: MetaDataHarvester Pipeline
Version: 2.0 - Streamlined and Robust
"""

import os
import sys
import argparse
import logging
from pathlib import Path
from typing import List, Dict, Optional
import pandas as pd
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
import time
import re


class Config:
    """Configuration for NCBIProteinHarvester"""
    def __init__(self, accession_file: str, output_dir: str, email: str, api_key: Optional[str] = None):
        self.accession_file = accession_file
        self.output_dir = Path(output_dir)
        self.email = email
        self.api_key = api_key
        self.batch_size = 50  # Optimal batch size
        self.retry_attempts = 3
        self.timeout = 30


class NCBIProteinHarvester:
    """
    Streamlined protein sequence and metadata harvester
    """

    def __init__(self, config: Config):
        self.config = config
        self.setup_logging()
        self.setup_ncbi()
        self.setup_directories()

        # Statistics
        self.stats = {
            'total_ids': 0,
            'sequences_fetched': 0,
            'metadata_extracted': 0,
            'fasta_files_created': 0,
            'failed_sequences': 0
        }

    def setup_logging(self):
        """Setup logging configuration"""
        log_dir = self.config.output_dir / "logs"
        log_dir.mkdir(parents=True, exist_ok=True)

        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(log_dir / "harvest.log"),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger('NCBIProteinHarvester')

    def setup_ncbi(self):
        """Setup NCBI Entrez configuration"""
        Entrez.email = self.config.email
        if self.config.api_key:
            Entrez.api_key = self.config.api_key
            self.logger.info("Using NCBI API key for faster downloads")

    def setup_directories(self):
        """Create necessary output directories"""
        self.config.output_dir.mkdir(parents=True, exist_ok=True)
        (self.config.output_dir / "fasta").mkdir(parents=True, exist_ok=True)
        (self.config.output_dir / "logs").mkdir(parents=True, exist_ok=True)

    def parse_accession_file(self) -> List[str]:
        """Parse accession file to get protein IDs"""
        self.logger.info(f"Parsing accession file: {self.config.accession_file}")

        if not Path(self.config.accession_file).exists():
            raise FileNotFoundError(f"Accession file not found: {self.config.accession_file}")

        protein_ids = []
        with open(self.config.accession_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    accession = line.split()[0]  # Take first word
                    if self._is_valid_accession(accession):
                        protein_ids.append(accession)

        # Remove duplicates
        protein_ids = list(set(protein_ids))
        self.logger.info(f"Found {len(protein_ids)} unique valid accessions")
        self.stats['total_ids'] = len(protein_ids)
        return protein_ids

    def _is_valid_accession(self, accession: str) -> bool:
        """Validate protein accession format"""
        if not accession or len(accession) < 3:
            return False

        # Common NCBI accession patterns
        patterns = [
            r'^[A-Z]{1,2}_\d+(\.\d+)?$',  # Standard format (NP_123456, XP_123456, WP_123456.1)
            r'^\d+$',                      # GI numbers
            r'^[A-Z]{2,4}\d{5,}(\.\d+)?$', # GenBank format (KWV17775.1, APQ19878.1, etc.)
        ]

        return any(re.match(pattern, accession) for pattern in patterns)

    def fetch_sequences_batch(self, protein_ids: List[str]) -> List[SeqRecord]:
        """Fetch protein sequences in batches with robust error handling"""
        sequences = []

        for attempt in range(self.config.retry_attempts):
            try:
                self.logger.info(f"Fetching {len(protein_ids)} sequences (attempt {attempt + 1})")

                # Join IDs for batch request
                id_string = ",".join(protein_ids)

                # Fetch from NCBI
                handle = Entrez.efetch(
                    db="protein",
                    id=id_string,
                    rettype="fasta",
                    retmode="text",
                    timeout=self.config.timeout
                )

                # Parse sequences
                batch_sequences = list(SeqIO.parse(handle, "fasta"))
                handle.close()

                self.logger.info(f"Successfully fetched {len(batch_sequences)} sequences")
                sequences.extend(batch_sequences)
                break

            except Exception as e:
                self.logger.warning(f"Attempt {attempt + 1} failed: {e}")
                if attempt < self.config.retry_attempts - 1:
                    time.sleep(2 ** attempt)  # Exponential backoff
                else:
                    self.logger.error("All fetch attempts failed")
                    break

        return sequences

    def extract_metadata_batch(self, protein_ids: List[str]) -> pd.DataFrame:
        """Extract metadata for protein IDs"""
        metadata_records = []

        for attempt in range(self.config.retry_attempts):
            try:
                self.logger.info(f"Extracting metadata for {len(protein_ids)} proteins (attempt {attempt + 1})")

                # Join IDs for batch request
                id_string = ",".join(protein_ids)

                # Fetch metadata
                handle = Entrez.efetch(
                    db="protein",
                    id=id_string,
                    rettype="gb",
                    retmode="xml",
                    timeout=self.config.timeout
                )

                records = Entrez.read(handle)
                handle.close()

                # Process each record
                for record in records:
                    metadata = self._extract_record_metadata(record)
                    if metadata:
                        metadata_records.append(metadata)

                self.logger.info(f"Extracted metadata for {len(metadata_records)} proteins")
                break

            except Exception as e:
                self.logger.warning(f"Metadata extraction attempt {attempt + 1} failed: {e}")
                if attempt < self.config.retry_attempts - 1:
                    time.sleep(2 ** attempt)
                else:
                    self.logger.error("All metadata extraction attempts failed")
                    break

        return pd.DataFrame(metadata_records)

    def _extract_record_metadata(self, record: dict) -> Optional[Dict]:
        """Extract metadata from a single GenBank record"""
        try:
            accession = record.get('GBSeq_accession-version', record.get('GBSeq_primary-accession', 'Unknown'))
            organism = record.get('GBSeq_organism', 'Unknown')
            definition = record.get('GBSeq_definition', '')
            length = record.get('GBSeq_length', 0)

            # Extract protein name and family
            protein_name = self._extract_protein_name(definition)
            protein_family = self._determine_protein_family(definition)

            return {
                'Accession_Number': accession,
                'Organism': organism,
                'Protein_Name': protein_name,
                'Protein_Family': protein_family,
                'Protein_Length': int(length) if length else 0,
                'Full_Description': definition
            }

        except Exception as e:
            self.logger.debug(f"Failed to extract metadata: {e}")
            return None

    def _extract_protein_name(self, definition: str) -> str:
        """Extract protein name from definition"""
        if not definition:
            return 'Unknown'

        # Simple extraction patterns
        patterns = [
            r'(\w+(?:\s+\w+)*)\s+\[',  # Name before organism
            r'^([^,\[\]]+)',           # Everything before first comma/bracket
        ]

        for pattern in patterns:
            match = re.search(pattern, definition)
            if match:
                name = match.group(1).strip()
                if len(name) > 3:
                    return name

        return definition.split()[0] if definition.split() else 'Unknown'

    def _determine_protein_family(self, definition: str) -> str:
        """Determine protein family from definition"""
        text = definition.lower()

        # RND efflux pump families
        if 'acra' in text or 'periplasmic adaptor' in text:
            return 'acra'
        elif 'acrb' in text or 'transporter' in text:
            return 'acrb'
        elif 'tolc' in text or 'outer membrane' in text:
            return 'tolc'
        elif 'acr' in text and 'repressor' in text:
            return 'acrr'
        elif 'mar' in text and 'activator' in text:
            return 'mara'
        elif 'sox' in text and 'activator' in text:
            return 'soxs'
        elif 'oqx' in text:
            return 'oqxb'
        elif 'mdt' in text:
            return 'mdtb'
        elif 'eef' in text:
            return 'eefa'
        elif 'rnd' in text or 'efflux' in text:
            return 'RND_Efflux'
        else:
            return 'Other'

    def save_sequences_by_family(self, sequences: List[SeqRecord], metadata_df: pd.DataFrame):
        """Save sequences organized by protein family"""
        fasta_dir = self.config.output_dir / "fasta"
        family_sequences = {}

        for seq in sequences:
            # Extract accession from sequence ID
            seq_id = seq.id.split('|')[1] if '|' in seq.id else seq.id

            # Find corresponding metadata
            metadata_row = metadata_df[metadata_df['Accession_Number'] == seq_id]

            if not metadata_row.empty:
                family = metadata_row.iloc[0]['Protein_Family']
            else:
                family = 'Unknown'

            if family not in family_sequences:
                family_sequences[family] = []

            family_sequences[family].append(seq)

        # Save by family
        files_created = 0
        for family, family_seqs in family_sequences.items():
            if family_seqs:
                filename = f"{family.replace(' ', '_')}.fasta"
                filepath = fasta_dir / filename

                try:
                    SeqIO.write(family_seqs, filepath, "fasta")
                    files_created += 1
                    self.logger.info(f"Saved {len(family_seqs)} sequences to {filename}")
                except Exception as e:
                    self.logger.error(f"Failed to save {filename}: {e}")

        self.stats['fasta_files_created'] = files_created

    def run_harvest(self) -> bool:
        """Main harvesting process"""
        try:
            self.logger.info("Starting streamlined NCBI Protein Harvester")
            self.logger.info(f"Accession file: {self.config.accession_file}")
            self.logger.info(f"Output directory: {self.config.output_dir}")

            # Parse accession file
            self.logger.info("Step 1: Parsing accession file")
            protein_ids = self.parse_accession_file()

            if not protein_ids:
                self.logger.error("No valid protein IDs found")
                return False

            # Process in batches
            self.logger.info("Step 2: Fetching sequences and metadata")
            all_sequences = []
            all_metadata = []

            for i in range(0, len(protein_ids), self.config.batch_size):
                batch_ids = protein_ids[i:i + self.config.batch_size]
                batch_num = (i // self.config.batch_size) + 1
                total_batches = (len(protein_ids) + self.config.batch_size - 1) // self.config.batch_size

                self.logger.info(f"Processing batch {batch_num}/{total_batches} ({len(batch_ids)} proteins)")

                # Fetch sequences
                sequences = self.fetch_sequences_batch(batch_ids)
                if sequences:
                    all_sequences.extend(sequences)
                    self.stats['sequences_fetched'] += len(sequences)

                # Extract metadata
                metadata_df = self.extract_metadata_batch(batch_ids)
                if not metadata_df.empty:
                    all_metadata.append(metadata_df)
                    self.stats['metadata_extracted'] += len(metadata_df)

            # Combine metadata
            self.logger.info("Step 3: Combining results")
            if all_metadata:
                combined_metadata = pd.concat(all_metadata, ignore_index=True)
            else:
                self.logger.warning("No metadata extracted")
                combined_metadata = pd.DataFrame()

            # Save results
            self.logger.info("Step 4: Saving results")

            # Save metadata
            if not combined_metadata.empty:
                metadata_file = self.config.output_dir / "metadata.csv"
                combined_metadata.to_csv(metadata_file, index=False)
                self.logger.info(f"Saved metadata for {len(combined_metadata)} proteins")

            # Save sequences by family
            if all_sequences:
                self.save_sequences_by_family(all_sequences, combined_metadata)

            # Print statistics
            self.print_statistics()

            success_rate = (self.stats['sequences_fetched'] / self.stats['total_ids'] * 100) if self.stats['total_ids'] > 0 else 0
            self.logger.info("HARVEST COMPLETED SUCCESSFULLY!")
            self.logger.info(f"Success rate: {success_rate:.1f}%")
            return True

        except Exception as e:
            self.logger.error(f"Harvest failed: {e}")
            import traceback
            self.logger.error(f"Traceback: {traceback.format_exc()}")
            return False

    def print_statistics(self):
        """Print harvest statistics"""
        self.logger.info("=" * 50)
        self.logger.info("HARVEST STATISTICS")
        self.logger.info("=" * 50)
        self.logger.info(f"Total protein IDs: {self.stats['total_ids']}")
        self.logger.info(f"Sequences fetched: {self.stats['sequences_fetched']}")
        self.logger.info(f"Metadata records: {self.stats['metadata_extracted']}")
        self.logger.info(f"FASTA files created: {self.stats['fasta_files_created']}")
        self.logger.info(f"Failed sequences: {self.stats['failed_sequences']}")
        self.logger.info("=" * 50)


def main():
    """Main entry point"""
    parser = argparse.ArgumentParser(
        description="Streamlined NCBI protein sequence and metadata extraction"
    )

    parser.add_argument(
        "--accession-file",
        required=True,
        help="Path to text file with accession numbers (one per line)"
    )

    parser.add_argument(
        "--output-dir",
        required=True,
        help="Output directory for results"
    )

    parser.add_argument(
        "--email",
        required=True,
        help="Email for NCBI Entrez (required by NCBI)"
    )

    parser.add_argument(
        "--api-key",
        help="NCBI API key for increased rate limits"
    )

    args = parser.parse_args()

    # Create configuration
    config = Config(
        accession_file=args.accession_file,
        output_dir=args.output_dir,
        email=args.email,
        api_key=args.api_key
    )

    # Run harvester
    harvester = NCBIProteinHarvester(config)
    success = harvester.run_harvest()

    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
