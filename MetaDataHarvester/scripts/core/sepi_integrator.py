import sys; import os; sys.path.append(os.path.join(os.path.dirname(__file__), '\\\\\\\\\src')); from configuration_manager import config_manager
#!/usr/bin/env python3
"""
SEPI 2.0 Integrator - Enhanced Protein Extraction Integration
Integrates SEPI 2.0 (https://github.com/vihaankulkarni29/sepi2.0) for efficient
RND efflux pump protein extraction and metadata collection.

Author: MetaDataHarvester Pipeline
Version: 2.0 - SEPI Integration
"""

import os
import sys
import logging
import subprocess
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import pandas as pd
import json
from datetime import datetime
import tempfile
import requests


class SEPIIntegrator:
    """
    Integrates SEPI 2.0 for automated RND efflux pump protein extraction
    """

    def __init__(self, sepi_path: Optional[str] = None, output_dir: str = "sepi_output"):
        self.sepi_path = Path(sepi_path) if sepi_path else self._find_sepi_installation()
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)

        # Setup logging
        self.logger = logging.getLogger(__name__)
        self._setup_logging()

        # Validate SEPI installation
        if not self._validate_sepi():
            raise RuntimeError("SEPI 2.0 not found or not properly installed")

    def _setup_logging(self):
        """Setup logging configuration"""
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )

    def _find_sepi_installation(self) -> Optional[Path]:
        """Find SEPI installation in common locations"""
        common_paths = [
            Path.home() / "sepi2.0",
            Path.home() / "Documents" / "sepi2.0",
            Path.home() / "Desktop" / "sepi2.0",
            Path("/opt/sepi2.0"),
            Path("/usr/local/sepi2.0")
        ]

        for path in common_paths:
            if path.exists() and (path / "sepi.py").exists():
                return path

        # Try to find in current environment
        try:
            result = subprocess.run(["which", "sepi"], capture_output=True, text=True)
            if result.returncode == 0:
                return Path(result.stdout.strip()).parent
        except:
            pass

        return None

    def _validate_sepi(self) -> bool:
        """Validate SEPI installation and functionality"""
        if not self.sepi_path or not self.sepi_path.exists():
            self.logger.error(f"SEPI path not found: {self.sepi_path}")
            return False

        sepi_script = self.sepi_path / "sepi.py"
        if not sepi_script.exists():
            self.logger.error(f"SEPI script not found: {sepi_script}")
            return False

        # Test SEPI functionality
        try:
            result = subprocess.run([
                sys.executable, str(sepi_script), "--help"
            ], capture_output=True, text=True, timeout=30)

            if result.returncode == 0 and "SEPI" in result.stdout:
                self.logger.info("SEPI validation successful")
                return True
            else:
                self.logger.error("SEPI validation failed")
                return False
        except Exception as e:
            self.logger.error(f"SEPI validation error: {e}")
            return False

    def extract_rnd_proteins(self, genome_files: List[str],
                           protein_family: str = "rnd_efflux",
                           output_format: str = "fasta") -> Dict[str, str]:
        """
        Extract RND efflux pump proteins using SEPI 2.0

        Args:
            genome_files: List of genome FASTA files
            protein_family: Target protein family
            output_format: Output format (fasta, gff, etc.)

        Returns:
            Dictionary mapping genome IDs to output file paths
        """
        self.logger.info(f"Extracting {protein_family} proteins from {len(genome_files)} genomes")

        results = {}

        for genome_file in genome_files:
            genome_path = Path(genome_file)
            if not genome_path.exists():
                self.logger.warning(f"Genome file not found: {genome_file}")
                continue

            genome_id = genome_path.stem
            output_prefix = self.output_dir / f"{genome_id}_{protein_family}"

            try:
                # Run SEPI extraction
                cmd = [
                    sys.executable,
                    str(self.sepi_path / "sepi.py"),
                    "--input", str(genome_path),
                    "--output", str(output_prefix),
                    "--family", protein_family,
                    "--format", output_format,
                    "--threads", "4"
                ]

                self.logger.info(f"Running SEPI for {genome_id}")
                result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)

                if result.returncode == 0:
                    # Find output files
                    output_files = list(self.output_dir.glob(f"{genome_id}_{protein_family}*"))
                    if output_files:
                        results[genome_id] = str(output_files[0])
                        self.logger.info(f"Successfully extracted proteins for {genome_id}")
                    else:
                        self.logger.warning(f"No output files found for {genome_id}")
                else:
                    self.logger.error(f"SEPI failed for {genome_id}: {result.stderr}")

            except subprocess.TimeoutExpired:
                self.logger.error(f"SEPI timeout for {genome_id}")
            except Exception as e:
                self.logger.error(f"SEPI error for {genome_id}: {e}")

        self.logger.info(f"SEPI extraction completed: {len(results)} successful")
        return results

    def batch_extract_from_accessions(self, accession_list: List[str],
                                    download_dir: str = "genome_downloads") -> Dict[str, str]:
        """
        Extract proteins from NCBI accessions using SEPI

        Args:
            accession_list: List of NCBI genome accessions
            download_dir: Directory to store downloaded genomes

        Returns:
            Dictionary mapping accessions to protein files
        """
        self.logger.info(f"Batch extracting from {len(accession_list)} accessions")

        download_path = Path(download_dir)
        download_path.mkdir(exist_ok=True)

        # Download genomes
        genome_files = self._download_genomes(accession_list, download_path)

        # Extract proteins
        protein_files = self.extract_rnd_proteins(genome_files)

        return protein_files

    def _download_genomes(self, accessions: List[str], download_dir: Path) -> List[str]:
        """Download genome sequences from NCBI"""
        from Bio import Entrez

        genome_files = []

        for accession in accessions:
            try:
                self.logger.info(f"Downloading {accession}")

                # Download genome
                handle = Entrez.efetch(db="nucleotide", id=accession,
                                     rettype="fasta", retmode="text")
                genome_data = handle.read()
                handle.close()

                # Save to file
                output_file = download_dir / f"{accession}.fasta"
                with open(output_file, 'w') as f:
                    f.write(genome_data)

                genome_files.append(str(output_file))
                self.logger.info(f"Downloaded {accession}")

            except Exception as e:
                self.logger.error(f"Failed to download {accession}: {e}")

        return genome_files

    def create_metadata_from_sepi_output(self, sepi_results: Dict[str, str]) -> pd.DataFrame:
        """
        Create comprehensive metadata from SEPI output

        Args:
            sepi_results: Dictionary from extract_rnd_proteins

        Returns:
            DataFrame with protein metadata
        """
        self.logger.info("Creating metadata from SEPI results")

        metadata_records = []

        for genome_id, protein_file in sepi_results.items():
            try:
                # Parse protein file
                proteins = self._parse_protein_file(protein_file)

                for protein in proteins:
                    metadata_records.append({
                        'accession': protein['accession'],
                        'genome_id': genome_id,
                        'protein_family': protein['family'],
                        'gene_name': protein['gene_name'],
                        'product': protein['product'],
                        'start_position': protein['start'],
                        'end_position': protein['end'],
                        'strand': protein['strand'],
                        'sequence_length': protein['length'],
                        'confidence_score': protein['confidence'],
                        'extraction_method': 'SEPI_2.0',
                        'timestamp': datetime.now().isoformat()
                    })

            except Exception as e:
                self.logger.error(f"Failed to parse {protein_file}: {e}")

        df = pd.DataFrame(metadata_records)
        self.logger.info(f"Created metadata for {len(df)} proteins")
        return df

    def _parse_protein_file(self, protein_file: str) -> List[Dict]:
        """Parse SEPI output protein file"""
        proteins = []

        try:
            with open(protein_file, 'r') as f:
                content = f.read()

            # Parse based on SEPI output format
            # This would need to be adapted based on SEPI's actual output format
            lines = content.split('\n')
            current_protein = {}

            for line in lines:
                if line.startswith('>'):
                    # New protein
                    if current_protein:
                        proteins.append(current_protein)

                    # Parse header
                    header = line[1:].strip()
                    parts = header.split('|')

                    current_protein = {
                        'accession': parts[0] if len(parts) > 0 else 'unknown',
                        'gene_name': parts[1] if len(parts) > 1 else 'unknown',
                        'family': parts[2] if len(parts) > 2 else 'unknown',
                        'start': int(parts[3]) if len(parts) > 3 else 0,
                        'end': int(parts[4]) if len(parts) > 4 else 0,
                        'strand': parts[5] if len(parts) > 5 else '+',
                        'confidence': float(parts[6]) if len(parts) > 6 else 1.0,
                        'product': parts[7] if len(parts) > 7 else 'hypothetical protein',
                        'sequence': '',
                        'length': 0
                    }
                else:
                    # Sequence line
                    if current_protein:
                        current_protein['sequence'] += line.strip()

            # Add last protein
            if current_protein:
                current_protein['length'] = len(current_protein['sequence'])
                proteins.append(current_protein)

        except Exception as e:
            self.logger.error(f"Error parsing protein file {protein_file}: {e}")

        return proteins

    def integrate_with_amr_pipeline(self, accession_list: List[str],
                                  protein_family: str = config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0]) -> Tuple[pd.DataFrame, Dict[str, str]]:
        """
        Complete integration with AMR pipeline

        Args:
            accession_list: List of genome accessions
            protein_family: Target protein family

        Returns:
            Tuple of (metadata DataFrame, protein file mapping)
        """
        self.logger.info("Starting SEPI-AMR pipeline integration")

        # Extract proteins using SEPI
        protein_files = self.batch_extract_from_accessions(accession_list)

        # Create metadata
        metadata_df = self.create_metadata_from_sepi_output(protein_files)

        # Filter for target protein family
        if protein_family != "all":
            metadata_df = metadata_df[metadata_df['protein_family'] == protein_family]

        # Save results
        metadata_file = self.output_dir / f"sepi_metadata_{protein_family}.csv"
        metadata_df.to_csv(metadata_file, index=False)

        self.logger.info(f"SEPI integration complete: {len(metadata_df)} proteins extracted")
        return metadata_df, protein_files

    def generate_sepi_report(self, metadata_df: pd.DataFrame,
                           protein_files: Dict[str, str]) -> None:
        """Generate SEPI integration report"""
        report_file = self.output_dir / "sepi_integration_report.txt"

        with open(report_file, 'w') as f:
            f.write("SEPI 2.0 INTEGRATION REPORT\n")
            f.write("=" * 40 + "\n\n")

            f.write("EXTRACTION SUMMARY:\n")
            f.write("-" * 20 + "\n")
            f.write(f"Total genomes processed: {len(protein_files)}\n")
            f.write(f"Total proteins extracted: {len(metadata_df)}\n")
            f.write(f"Protein families found: {metadata_df['protein_family'].nunique()}\n\n")

            if not metadata_df.empty:
                # Family distribution
                f.write("PROTEIN FAMILY DISTRIBUTION:\n")
                f.write("-" * 30 + "\n")
                family_counts = metadata_df['protein_family'].value_counts()
                for family, count in family_counts.items():
                    f.write(f"{family}: {count}\n")

                f.write("\nQUALITY METRICS:\n")
                f.write("-" * 20 + "\n")
                f.write(f"Average confidence: {metadata_df['confidence_score'].mean():.2f}\n")
                f.write(f"Average protein length: {metadata_df['sequence_length'].mean():.0f} aa\n")
                f.write(f"High confidence (>0.8): {(metadata_df['confidence_score'] > 0.8).sum()}\n")

        self.logger.info(f"SEPI report generated: {report_file}")


def main():
    """Main entry point"""
    import argparse

    parser = argparse.ArgumentParser(
        description="SEPI 2.0 integration for RND efflux pump extraction"
    )

    parser.add_argument(
        "--accession-list",
        required=True,
        help="File containing NCBI genome accessions (one per line)"
    )

    parser.add_argument(
        "--protein-family",
        default="rnd_efflux",
        help="Target protein family"
    )

    parser.add_argument(
        "--sepi-path",
        help="Path to SEPI 2.0 installation"
    )

    parser.add_argument(
        "--output-dir",
        default="sepi_integration_output",
        help="Output directory"
    )

    parser.add_argument(
        "--email",
        required=True,
        help="NCBI email for downloads"
    )

    args = parser.parse_args()

    # Setup NCBI
    from Bio import Entrez
    Entrez.email = args.email

    # Setup logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    # Initialize SEPI integrator
    integrator = SEPIIntegrator(args.sepi_path, args.output_dir)

    # Load accession list
    with open(args.accession_list, 'r') as f:
        accessions = [line.strip() for line in f if line.strip()]

    # Run integration
    metadata_df, protein_files = integrator.integrate_with_amr_pipeline(
        accessions, args.protein_family
    )

    # Generate report
    integrator.generate_sepi_report(metadata_df, protein_files)

    print(f"\nSEPI integration complete!")
    print(f"Processed {len(accessions)} genomes")
    print(f"Extracted {len(metadata_df)} proteins")
    print(f"Results saved to: {args.output_dir}")


if __name__ == "__main__":
    main()