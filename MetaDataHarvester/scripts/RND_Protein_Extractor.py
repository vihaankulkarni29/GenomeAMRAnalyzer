import sys; import os; sys.path.append(os.path.join(os.path.dirname(__file__), '\\\\\\src')); from configuration_manager import config_manager
#!/usr/bin/env python3
"""
RND_Protein_Extractor.py - Extract all RND efflux pump proteins from genomes

This tool extracts all RND (Resistance-Nodulation-Division) efflux pump proteins
from downloaded bacterial genomes. It searches for and extracts protein sequences
for the complete RND efflux pump system including:

- AcrAB-TolC system (AcrA, AcrB, TolC)
- Other RND pumps (AcrD, AcrE, AcrF, MdtB/C, OqxA/B, etc.)
- Regulatory proteins (MarA, RamA, SoxS, Rob, etc.)

Usage:
    python RND_Protein_Extractor.py --genome-dir data/genomes --output-dir extracted_rnd_proteins

Author: MetaDataHarvester Pipeline
Version: 1.0
"""

import os
import sys
import argparse
import logging
import re
from pathlib import Path
from typing import List, Dict, Optional, Set
import pandas as pd
from datetime import datetime

# BioPython imports
try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
except ImportError:
    print("ERROR: BioPython is required. Install with: pip install biopython")
    sys.exit(1)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


class RND_Protein_Extractor:
    """
    Extracts all RND efflux pump proteins from bacterial genomes.
    """

    def __init__(self, genome_dir: str, output_dir: str):
        """
        Initialize the RND protein extractor.

        Args:
            genome_dir: Directory containing downloaded genomes
            output_dir: Output directory for extracted proteins
        """
        self.genome_dir = Path(genome_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Define RND efflux pump proteins to extract
        self.rnd_proteins = {
            # AcrAB-TolC system
            config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0]: [config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0], 'AcrA', 'acrA_MG1655', 'EG11703'],
            config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1]: [config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1], 'AcrB', 'acrB_MG1655', 'ACRB'],
            config_manager.get_default_genes("rnd_efflux_pumps", "primary")[2]: [config_manager.get_default_genes("rnd_efflux_pumps", "primary")[2], 'TolC', 'tolC_MG1655'],

            # Other RND pumps
            'acrD': ['acrD', 'AcrD', 'CAD6000925'],
            'acrE': ['acrE', 'AcrE', 'CAD5995424'],
            'acrF': ['acrF', 'AcrF', 'CAD5995417'],
            'mdtB': ['mdtB', 'MdtB', 'CAD6003569'],
            'mdtC': ['mdtC', 'MdtC', 'CAD6003562'],
            'oqxA': ['oqxA', 'OqxA', 'WGZ82228'],
            'oqxB': ['oqxB', 'OqxB', 'XOG11541'],

            # Regulatory proteins
            'acrR': ['acrR', 'AcrR', 'CAM7858307', 'EG11741'],
            'acrZ': ['acrZ', 'AcrZ', 'VWQ01551'],
            'marA': ['marA', 'MarA', 'CAD6013869'],
            'marR': ['marR', 'MarR', 'CAD6013879'],
            'ramA': ['ramA', 'RamA', 'CAD6016219'],
            'rob': ['rob', 'Rob', 'CAD6017604'],
            'soxS': ['soxS', 'SoxS', 'AMW89049'],
            'envR': ['envR', 'EnvR', 'CAD6001772'],

            # EefABC system
            'eefA': ['eefA', 'EefA', 'AGT21931'],
            'eefB': ['eefB', 'EefB', 'ABV24995'],
            'eefC': ['eefC', 'EefC', 'OUH93039']
        }

        # Statistics
        self.stats = {
            'genomes_processed': 0,
            'total_proteins_extracted': 0,
            'proteins_by_type': {protein: 0 for protein in self.rnd_proteins.keys()},
            'failed_extractions': 0
        }

        logger.info("RND_Protein_Extractor initialized")
        logger.info(f"Genome directory: {self.genome_dir}")
        logger.info(f"Output directory: {self.output_dir}")
        logger.info(f"Target proteins: {len(self.rnd_proteins)} types")

    def extract_all_rnd_proteins(self) -> Dict:
        """
        Extract all RND proteins from all genomes in the genome directory.

        Returns:
            Dictionary with extraction results and statistics
        """
        logger.info("Starting RND protein extraction from all genomes")

        # Find all genome directories
        genome_dirs = self._find_genome_directories()
        if not genome_dirs:
            logger.error("No genome directories found")
            return {'success': False, 'error': 'No genome directories found'}

        logger.info(f"Found {len(genome_dirs)} genome directories to process")

        all_extracted_proteins = []

        for genome_dir in genome_dirs:
            logger.info(f"Processing genome: {genome_dir.name}")
            self.stats['genomes_processed'] += 1

            # Extract RND proteins from this genome
            extracted = self._extract_from_single_genome(genome_dir)
            if extracted:
                all_extracted_proteins.extend(extracted)

        # Save summary results
        self._save_extraction_summary(all_extracted_proteins)

        # Report final statistics
        self._report_final_statistics()

        logger.info(f"RND protein extraction complete. Extracted {len(all_extracted_proteins)} proteins.")

        return {
            'success': True,
            'total_proteins': len(all_extracted_proteins),
            'genomes_processed': self.stats['genomes_processed'],
            'output_directory': str(self.output_dir),
            'extraction_summary': str(self.output_dir / 'extraction_summary.csv')
        }

    def _find_genome_directories(self) -> List[Path]:
        """Find all genome directories in the genome directory."""
        genome_dirs = []

        if not self.genome_dir.exists():
            logger.error(f"Genome directory does not exist: {self.genome_dir}")
            return genome_dirs

        # Look for directories that contain genome files
        for item in self.genome_dir.iterdir():
            if item.is_dir():
                # Check if it contains FASTA and GFF files
                fasta_files = list(item.glob("*_genomic.fna"))
                gff_files = list(item.glob("*_genomic.gff"))

                if fasta_files and gff_files:
                    genome_dirs.append(item)

        return genome_dirs

    def _extract_from_single_genome(self, genome_dir: Path) -> List[Dict]:
        """
        Extract RND proteins from a single genome.

        Args:
            genome_dir: Path to genome directory

        Returns:
            List of extracted protein information
        """
        extracted_proteins = []

        try:
            # Find FASTA and GFF files
            fasta_files = list(genome_dir.glob("*_genomic.fna"))
            gff_files = list(genome_dir.glob("*_genomic.gff"))

            if not fasta_files or not gff_files:
                logger.warning(f"No FASTA or GFF files found in {genome_dir}")
                return extracted_proteins

            fasta_file = fasta_files[0]
            gff_file = gff_files[0]

            logger.info(f"Processing FASTA: {fasta_file.name}")
            logger.info(f"Processing GFF: {gff_file.name}")

            # Parse GFF file to find RND protein genes
            gff_features = self._parse_gff_file(gff_file)

            # Extract protein sequences
            for protein_type, search_terms in self.rnd_proteins.items():
                protein_info = self._extract_single_protein(
                    protein_type, search_terms, gff_features,
                    fasta_file, genome_dir
                )

                if protein_info:
                    extracted_proteins.append(protein_info)
                    self.stats['total_proteins_extracted'] += 1
                    self.stats['proteins_by_type'][protein_type] += 1

        except Exception as e:
            logger.error(f"Failed to extract proteins from {genome_dir}: {e}")
            self.stats['failed_extractions'] += 1

        return extracted_proteins

    def _parse_gff_file(self, gff_file: Path) -> List[Dict]:
        """Parse GFF file to extract gene features."""
        features = []

        try:
            with open(gff_file, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        continue

                    fields = line.strip().split('\t')
                    if len(fields) < 9:
                        continue

                    if fields[2] == 'CDS':  # Coding sequence
                        feature = {
                            'seqid': fields[0],
                            'source': fields[1],
                            'type': fields[2],
                            'start': int(fields[3]),
                            'end': int(fields[4]),
                            'score': fields[5],
                            'strand': fields[6],
                            'phase': fields[7],
                            'attributes': fields[8]
                        }

                        # Parse attributes
                        attr_dict = self._parse_gff_attributes(fields[8])
                        feature['attributes_dict'] = attr_dict

                        features.append(feature)

        except Exception as e:
            logger.warning(f"Failed to parse GFF file {gff_file}: {e}")

        return features

    def _parse_gff_attributes(self, attributes: str) -> Dict:
        """Parse GFF attribute string into dictionary."""
        attr_dict = {}

        # Split by semicolon and parse each attribute
        for attr in attributes.split(';'):
            if '=' in attr:
                key, value = attr.split('=', 1)
                attr_dict[key.strip()] = value.strip()
            elif ':' in attr:
                key, value = attr.split(':', 1)
                attr_dict[key.strip()] = value.strip()

        return attr_dict

    def _extract_single_protein(self, protein_type: str, search_terms: List[str],
                               gff_features: List[Dict], fasta_file: Path,
                               genome_dir: Path) -> Optional[Dict]:
        """
        Extract a single RND protein from genome.

        Args:
            protein_type: Type of protein (e.g., config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0])
            search_terms: List of terms to search for in gene names
            gff_features: Parsed GFF features
            fasta_file: Path to FASTA file
            genome_dir: Genome directory

        Returns:
            Protein information dictionary or None
        """
        # Find matching features in GFF
        matching_features = self._find_matching_features(search_terms, gff_features)

        if not matching_features:
            logger.debug(f"No {protein_type} features found")
            return None

        # Use the first match (could be improved to select best match)
        feature = matching_features[0]

        try:
            # Extract protein sequence
            protein_seq = self._extract_protein_sequence(feature, fasta_file)

            if not protein_seq:
                logger.warning(f"Could not extract sequence for {protein_type}")
                return None

            # Create protein record
            accession = genome_dir.name.split('_')[0]  # Extract accession from directory name
            protein_id = f"{accession}_{protein_type}"

            protein_record = SeqRecord(
                Seq(protein_seq),
                id=protein_id,
                name=protein_type,
                description=f"{protein_type} from {accession}"
            )

            # Save to FASTA file
            output_file = self.output_dir / f"{protein_id}.fasta"
            with open(output_file, 'w') as f:
                SeqIO.write(protein_record, f, 'fasta')

            logger.info(f"Extracted {protein_type} to {output_file}")

            return {
                'protein_type': protein_type,
                'genome_accession': accession,
                'protein_id': protein_id,
                'sequence_length': len(protein_seq),
                'gene_location': f"{feature['seqid']}:{feature['start']}-{feature['end']}",
                'strand': feature['strand'],
                'fasta_file': str(output_file),
                'search_terms_used': search_terms,
                'gff_attributes': feature['attributes']
            }

        except Exception as e:
            logger.warning(f"Failed to extract {protein_type}: {e}")
            return None

    def _find_matching_features(self, search_terms: List[str], gff_features: List[Dict]) -> List[Dict]:
        """Find GFF features that match the search terms."""
        matching_features = []

        for feature in gff_features:
            attr_dict = feature['attributes_dict']

            # Check various attribute fields for matches
            searchable_text = ''
            if 'gene' in attr_dict:
                searchable_text += attr_dict['gene'] + ' '
            if 'locus_tag' in attr_dict:
                searchable_text += attr_dict['locus_tag'] + ' '
            if 'product' in attr_dict:
                searchable_text += attr_dict['product'] + ' '
            if 'protein_id' in attr_dict:
                searchable_text += attr_dict['protein_id'] + ' '

            searchable_text = searchable_text.lower()

            # Check if any search term matches
            for term in search_terms:
                if term.lower() in searchable_text:
                    matching_features.append(feature)
                    break

        return matching_features

    def _extract_protein_sequence(self, feature: Dict, fasta_file: Path) -> Optional[str]:
        """Extract protein sequence from FASTA file using GFF coordinates."""
        try:
            # Parse FASTA file
            records = list(SeqIO.parse(fasta_file, 'fasta'))

            if not records:
                return None

            # Find the contig/chromosome that matches the feature
            target_seqid = feature['seqid']
            genome_record = None

            for record in records:
                if record.id == target_seqid or target_seqid in record.id:
                    genome_record = record
                    break

            if not genome_record:
                logger.warning(f"Could not find sequence {target_seqid} in FASTA file")
                return None

            # Extract DNA sequence
            start = feature['start'] - 1  # Convert to 0-based indexing
            end = feature['end']
            strand = feature['strand']

            if strand == '+':
                dna_seq = genome_record.seq[start:end]
            else:
                dna_seq = genome_record.seq[start:end].reverse_complement()

            # Translate to protein (assuming standard genetic code)
            protein_seq = dna_seq.translate(to_stop=True)

            return str(protein_seq)

        except Exception as e:
            logger.warning(f"Failed to extract protein sequence: {e}")
            return None

    def _save_extraction_summary(self, extracted_proteins: List[Dict]):
        """Save extraction summary to CSV file."""
        if not extracted_proteins:
            logger.warning("No proteins extracted - skipping summary save")
            return

        # Convert to DataFrame
        df = pd.DataFrame(extracted_proteins)

        # Save to CSV
        summary_file = self.output_dir / 'extraction_summary.csv'
        df.to_csv(summary_file, index=False)

        logger.info(f"Extraction summary saved to {summary_file}")

    def _report_final_statistics(self):
        """Report final extraction statistics."""
        logger.info("=" * 60)
        logger.info("RND PROTEIN EXTRACTION STATISTICS")
        logger.info("=" * 60)
        logger.info(f"Genomes processed: {self.stats['genomes_processed']}")
        logger.info(f"Total proteins extracted: {self.stats['total_proteins_extracted']}")
        logger.info(f"Failed extractions: {self.stats['failed_extractions']}")
        logger.info("")
        logger.info("Proteins extracted by type:")
        for protein_type, count in self.stats['proteins_by_type'].items():
            if count > 0:
                logger.info(f"  {protein_type}: {count}")
        logger.info("=" * 60)


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Extract all RND efflux pump proteins from bacterial genomes",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Extract RND proteins from downloaded genomes
  python RND_Protein_Extractor.py --genome-dir data/genomes --output-dir extracted_rnd_proteins

  # Use custom directories
  python RND_Protein_Extractor.py --genome-dir my_genomes --output-dir my_extracted_proteins

Target RND Proteins:
  - AcrAB-TolC system: AcrA, AcrB, TolC
  - Other RND pumps: AcrD, AcrE, AcrF, MdtB/C, OqxA/B
  - Regulatory proteins: MarA, RamA, SoxS, Rob, AcrR, AcrZ
  - EefABC system: EefA, EefB, EefC
        """
    )

    parser.add_argument(
        '--genome-dir',
        required=True,
        help='Directory containing downloaded genomes'
    )

    parser.add_argument(
        '--output-dir',
        required=True,
        help='Output directory for extracted proteins'
    )

    parser.add_argument(
        '--log-level',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
        default='INFO',
        help='Set logging level (default: INFO)'
    )

    args = parser.parse_args()

    # Set logging level
    logging.getLogger().setLevel(getattr(logging, args.log_level))

    # Create extractor
    extractor = RND_Protein_Extractor(
        genome_dir=args.genome_dir,
        output_dir=args.output_dir
    )

    # Extract proteins
    results = extractor.extract_all_rnd_proteins()

    if results['success']:
        logger.info("RND protein extraction completed successfully!")
        logger.info(f"Extracted proteins are available in: {results['output_directory']}")
        logger.info(f"Extraction summary: {results['extraction_summary']}")

        print("\n" + "="*80)
        print("RND PROTEIN EXTRACTION SUMMARY")
        print("="*80)
        print(f"Genomes processed: {results['genomes_processed']}")
        print(f"Total proteins extracted: {results['total_proteins']}")
        print(f"Output directory: {results['output_directory']}")
        print("="*80)

        sys.exit(0)
    else:
        logger.error(f"RND protein extraction failed: {results.get('error', 'Unknown error')}")
        sys.exit(1)


if __name__ == "__main__":
    main()