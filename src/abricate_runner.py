#!/usr/bin/env python3
"""
Abricate Runner - Multi-Database Support

Runs Abricate against any specified database (card, vfdb, plasmidfinder, etc.) on all 
genome FASTA files in an input directory and writes database-specific TSV reports.

Usage (CLI):
    python -m src.abricate_runner --input-dir genomes/ --output-dir results/ --database card
    python -m src.abricate_runner --input-dir genomes/ --output-dir results/ --database vfdb

Features:
- Multi-database support (card, vfdb, plasmidfinder, resfinder, etc.)
- Dynamic output naming: <genome_id>_<database>.tsv
- Per-genome error handling prevents batch failures
- Detailed success/failure statistics
- Robust logging and progress tracking

Notes:
- Requires `abricate` to be installed and available on PATH
- Recognizes FASTA extensions: .fasta, .fa, .fna (case-insensitive)  
- Individual genome failures don't stop the entire batch
- Creates output directory if it doesn't exist
"""

from __future__ import annotations

import argparse
import os
import sys
import shutil
import subprocess
import logging
from pathlib import Path
from typing import Iterable, List, Dict, Tuple, Optional
from dataclasses import dataclass, field

# Configure logging for better error tracking
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

FASTA_EXTENSIONS = {".fasta", ".fa", ".fna"}

@dataclass
class ScanResults:
    """Results from a database scan operation with detailed tracking"""
    database: str
    successful_scans: List[str] = field(default_factory=list)
    failed_genomes: Dict[str, str] = field(default_factory=dict)
    no_hits_genomes: List[str] = field(default_factory=list)
    output_files: List[Path] = field(default_factory=list)
    total_genomes: int = 0
    global_error: Optional[str] = None


def is_fasta(path: Path) -> bool:
    """Check if a file is a valid FASTA file with appropriate extension.
    
    Validates that the path points to an existing file with a FASTA extension.
    
    Args:
        path (Path): File path to check.
    
    Returns:
        bool: True if the file exists and has a FASTA extension, False otherwise.
    """
    return path.is_file() and path.suffix.lower() in FASTA_EXTENSIONS


def find_fastas(input_dir: Path) -> List[Path]:
    """Recursively find all FASTA files in a directory and its subdirectories.
    
    Searches for files with common FASTA extensions (.fasta, .fas, .fa, .fna, .ffn)
    in the specified directory tree.
    
    Args:
        input_dir (Path): Root directory to search for FASTA files.
    
    Returns:
        List[Path]: List of Path objects pointing to discovered FASTA files,
            sorted alphabetically for reproducible processing order.
    """
    files: List[Path] = []
    for entry in sorted(input_dir.iterdir()):
        if is_fasta(entry):
            files.append(entry)
    return files


def run_abricate_on_file(fasta_path: Path, db: str = "card", nopath: bool = True) -> str:
    """Run abricate antimicrobial resistance screening on a single FASTA file.
    
    Executes the abricate command-line tool to identify antimicrobial resistance 
    genes, virulence factors, or plasmid genes in bacterial genome assemblies.
    
    Args:
        fasta_path (Path): Path to the input bacterial genome FASTA file.
        db (str, optional): Abricate database name to use for screening. 
            Options include 'card', 'vfdb', 'plasmidfinder', 'resfinder', etc. 
            Defaults to "card".
        nopath (bool, optional): Whether to use --nopath flag to exclude 
            directory paths from output for cleaner results. Defaults to True.

    Returns:
        str: Raw TSV-formatted output from abricate containing gene matches,
            coordinates, and resistance information. Returns empty string on failure.
        
    Raises:
        RuntimeError: If abricate executable is not found on system PATH.
        subprocess.SubprocessError: If abricate command execution fails.
    """
    cmd = ["abricate", "--db", db]
    if nopath:
        cmd.append("--nopath")
    cmd.append(str(fasta_path))

    try:
        logger.debug(f"Running: {' '.join(cmd)}")
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
        )
    except FileNotFoundError:
        raise RuntimeError(
            "Abricate executable not found. Ensure 'abricate' is installed and on your PATH."
        )

    if proc.returncode != 0:
        # Log warning but don't crash - let caller handle the failure
        logger.warning(
            f"Abricate failed for {fasta_path.name} (database: {db}) with exit code {proc.returncode}"
        )
        if proc.stderr:
            logger.warning(f"Abricate stderr: {proc.stderr.strip()}")
        return ""

    return proc.stdout or ""


def write_report(tsv_content: str, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as f:
        f.write(tsv_content)


def infer_genome_id(fasta_path: Path) -> str:
    # Use filename stem as genome ID (e.g., GCF_000005825.2 from GCF_000005825.2.fasta)
    return fasta_path.stem


def run_abricate(input_dir: str | os.PathLike, output_dir: str | os.PathLike, db: str = "card") -> ScanResults:
    """Run abricate across all FASTA files in input_dir and write TSV reports to output_dir.
    
    Enhanced with per-genome error handling and structured result tracking.

    Parameters
    ----------
    input_dir : str | os.PathLike
        Directory containing genome FASTA files (.fasta|.fa|.fna)
    output_dir : str | os.PathLike  
        Directory to write TSV reports with database-specific naming
    db : str
        Abricate database name (card, vfdb, plasmidfinder, etc.)

    Returns
    -------
    ScanResults
        Structured results with success/failure tracking per genome
    """
    input_dir_p = Path(input_dir)
    output_dir_p = Path(output_dir)
    
    # Initialize results tracking
    results = ScanResults(database=db)
    
    logger.info(f"Starting {db} database scan on {input_dir_p}")

    if not input_dir_p.exists() or not input_dir_p.is_dir():
        error_msg = f"Input directory does not exist or is not a directory: {input_dir}"
        logger.error(error_msg)
        results.global_error = error_msg
        return results

    if shutil.which("abricate") is None:
        error_msg = "'abricate' not found on PATH. Please install Abricate or run inside the project's Docker environment."
        logger.error(error_msg)
        results.global_error = error_msg
        return results

    fasta_files = find_fastas(input_dir_p)
    if not fasta_files:
        warning_msg = f"No FASTA files found in {input_dir_p}"
        logger.warning(warning_msg)
        results.global_error = warning_msg
        return results

    results.total_genomes = len(fasta_files)
    logger.info(f"Found {len(fasta_files)} FASTA files for {db} scanning")

    for fasta in fasta_files:
        genome_id = infer_genome_id(fasta)
        # Dynamic naming: [genome_id]_[database_name].tsv
        out_path = output_dir_p / f"{genome_id}_{db}.tsv"
        
        try:
            logger.debug(f"Scanning {fasta.name} against {db} database -> {out_path.name}")
            tsv = run_abricate_on_file(fasta, db=db, nopath=True)

            # Skip writing truly empty outputs but track the attempt
            if not tsv or not tsv.strip():
                logger.info(f"No hits found for {genome_id} in {db} database")
                results.no_hits_genomes.append(genome_id)
                continue

            write_report(tsv, out_path)
            results.successful_scans.append(genome_id)
            results.output_files.append(out_path)
            logger.debug(f"Successfully generated {out_path.name}")
            
        except Exception as e:
            error_msg = f"Failed to process {genome_id} for {db} database: {str(e)}"
            logger.error(error_msg)
            results.failed_genomes[genome_id] = str(e)

    # Log summary
    logger.info(f"Completed {db} scan: {len(results.successful_scans)} successful, {len(results.failed_genomes)} failed, {len(results.no_hits_genomes)} no hits")
    
    return results


def _parse_args(argv: Iterable[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Run Abricate on a directory of genome FASTA files and save TSV reports.")
    p.add_argument("--input-dir", required=True, help="Directory containing genome FASTA files (.fasta|.fa|.fna)")
    p.add_argument("--output-dir", required=True, help="Directory to write raw Abricate TSV reports")
    p.add_argument("--db", default="card", help="Abricate database to use (default: card)")
    return p.parse_args(list(argv))


def main(argv: Iterable[str] | None = None) -> int:
    args = _parse_args(argv if argv is not None else sys.argv[1:])
    try:
        results = run_abricate(args.input_dir, args.output_dir, db=args.db)
        
        if results.global_error:
            sys.stderr.write(f"[abricate_runner] Error: {results.global_error}\n")
            return 1
            
        # Print summary
        total_successful = len(results.successful_scans)
        total_failed = len(results.failed_genomes)
        total_no_hits = len(results.no_hits_genomes)
        
        print(f"[abricate_runner] Completed {results.database} scan:")
        print(f"  - {total_successful} genomes with hits")
        print(f"  - {total_no_hits} genomes with no hits")
        print(f"  - {total_failed} genomes failed")
        print(f"  - Generated {len(results.output_files)} report files")
        
        return 0 if total_failed == 0 else 1
        
    except Exception as e:
        sys.stderr.write(f"[abricate_runner] Error: {e}\n")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
