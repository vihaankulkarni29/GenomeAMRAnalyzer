#!/usr/bin/env python3
"""
URL → Genomes → Abricate (CARD) Orchestrator
======================================
End-to-end pipeline:
  1) Parse NCBI URL, discover genomes
  2) Download genomes (FASTA)
  3) Run RGI (CARD) on the downloaded FASTAs
  4) Emit coordinates and summary for downstream modules

Usage:
  python -m src.url_to_card_pipeline --url "<NCBI URL>" --config config/snakemake_config.yaml
"""

import argparse
import sys
import os
from pathlib import Path
import asyncio
import logging
from typing import List, Optional

# Ensure src imports
CURRENT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = CURRENT_DIR.parent
if str(CURRENT_DIR) not in sys.path:
    sys.path.insert(0, str(CURRENT_DIR))

from core.url_to_genomes_workflow import URLToGenomesWorkflow
from abricate_runner import run_abricate
from abricate_to_coords import convert_abricate_to_coords


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Run URL→Abricate (CARD) pipeline")
    p.add_argument("--url", required=True, help="NCBI nuccore search URL")
    p.add_argument("--config", default=str(PROJECT_ROOT / "config" / "snakemake_config.yaml"), help="Path to YAML config")
    # Kept for compatibility; no longer used, but accepted to avoid breaking scripts
    p.add_argument("--rgi-exec", default="rgi", help="(Deprecated) RGI executable path/name")
    p.add_argument("--threads", type=int, default=None, help="(Deprecated) RGI threads")
    return p


async def run_url_to_genomes(config_path: str, url: str):
    workflow = URLToGenomesWorkflow(config_path)
    files, report = await workflow.run_complete_workflow(url)
    return workflow, files, report


def run_rgi_on_files(genome_files: List[str], config_path: str, rgi_exec: str, threads_override: Optional[int] = None) -> bool:
    """Run Abricate on the discovered genomes and convert outputs to coordinate CSVs.

    Note: Function name kept for backward compatibility with callers.
    """
    import yaml
    from pathlib import Path
    import os

    with open(config_path, "r") as f:
        cfg = yaml.safe_load(f)

    # Reuse the existing directory entry in config for outputs
    base_out = Path(cfg['directories']['card_results']).resolve()
    coords_dir = base_out
    abricate_raw = base_out / "abricate_raw"
    abricate_raw.mkdir(parents=True, exist_ok=True)
    coords_dir.mkdir(parents=True, exist_ok=True)

    if not genome_files:
        return False

    input_dir = str(Path(genome_files[0]).parent)

    # Step 1: Run Abricate on input directory
    tsv_reports = run_abricate(input_dir, str(abricate_raw), db='card')

    # Step 2: Convert each TSV to legacy coordinate CSV
    generated = 0
    for tsv in tsv_reports:
        tsv_path = Path(tsv)
        genome_id = tsv_path.stem
        if genome_id.endswith('_abricate'):
            genome_id = genome_id[:-len('_abricate')]
        out_csv = coords_dir / f"{genome_id}_coordinates.csv"
        rows = convert_abricate_to_coords(tsv_path, out_csv)
        if rows >= 0:
            generated += 1

    return generated > 0


def main():
    parser = build_arg_parser()
    args = parser.parse_args()

    # Basic logging
    logging.basicConfig(level=logging.INFO, format='%(levelname)s %(message)s')

    # Run URL → genomes
    workflow, files, report = asyncio.run(run_url_to_genomes(args.config, args.url))
    logging.info(f"Discovered and downloaded {len(files)} genomes")
    
    # Validate for Abricate
    if not workflow.validate_for_rgi_processing():
        logging.warning("Some genomes may not be ideal for CARD scanning; continuing anyway")

    # Run Abricate → Coordinates
    ok = run_rgi_on_files(files, args.config, args.rgi_exec, args.threads)
    if not ok:
        sys.exit(1)

    logging.info("URL→Abricate pipeline completed. CARD-based coordinates ready.")


if __name__ == "__main__":
    main()
