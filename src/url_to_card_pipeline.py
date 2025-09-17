#!/usr/bin/env python3
"""
URL → Genomes → RGI (CARD) Orchestrator
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
from card_runner import CARDRunner, RGIConfig


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Run URL→RGI pipeline")
    p.add_argument("--url", required=True, help="NCBI nuccore search URL")
    p.add_argument("--config", default=str(PROJECT_ROOT / "config" / "snakemake_config.yaml"), help="Path to YAML config")
    p.add_argument("--rgi-exec", default="rgi", help="RGI executable path/name")
    p.add_argument("--threads", type=int, default=None, help="Override RGI threads")
    return p


async def run_url_to_genomes(config_path: str, url: str):
    workflow = URLToGenomesWorkflow(config_path)
    files, report = await workflow.run_complete_workflow(url)
    return workflow, files, report


def run_rgi_on_files(genome_files: List[str], config_path: str, rgi_exec: str, threads_override: Optional[int] = None) -> bool:
    import yaml
    with open(config_path, "r") as f:
        cfg = yaml.safe_load(f)

    output_dir = cfg['directories']['card_results']
    target_genes = cfg['target_genes']
    threads = threads_override or cfg.get('rgi', {}).get('threads', 1)

    rgi_cfg = RGIConfig(
        input_dir=str(Path(genome_files[0]).parent),
        output_dir=output_dir,
        target_genes=target_genes,
        rgi_executable=rgi_exec,
        num_threads=threads,
    )

    runner = CARDRunner(rgi_cfg)

    # Ensure only the discovered files are processed by creating a temp dir
    # pointing to the same directory is okay since CARDRunner scans *.fasta
    return runner.run_analysis()


def main():
    parser = build_arg_parser()
    args = parser.parse_args()

    # Basic logging
    logging.basicConfig(level=logging.INFO, format='%(levelname)s %(message)s')

    # Run URL → genomes
    workflow, files, report = asyncio.run(run_url_to_genomes(args.config, args.url))
    logging.info(f"Discovered and downloaded {len(files)} genomes")
    
    # Validate for RGI
    if not workflow.validate_for_rgi_processing():
        logging.warning("Some genomes may not be ideal for RGI; continuing anyway")

    # Run RGI
    ok = run_rgi_on_files(files, args.config, args.rgi_exec, args.threads)
    if not ok:
        sys.exit(1)

    logging.info("URL→RGI pipeline completed. CARD results ready.")


if __name__ == "__main__":
    main()
