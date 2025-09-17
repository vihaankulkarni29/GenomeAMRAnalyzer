"""
CARD RGI Integration Module
==========================

- Runs CARD RGI on downloaded genome FASTA files
- Parses RGI output for gene coordinates and resistance annotations
- Stores results and provenance for downstream extraction
- Robust error handling and status tracking

Engineering Principles:
- Fail-fast on missing input or RGI errors
- Modular, testable, and extensible
- Full provenance and logging
"""

import os
import subprocess
import logging
import time
from pathlib import Path
from typing import List, Dict, Any, Optional, Callable
import json

class CARDRGIRunner:
    """
    Runs CARD RGI on genome FASTA files and parses output.
    """
    def __init__(self,
                 rgi_path: str = "rgi",
                 card_db: Optional[str] = None,
                 output_dir: str = "card_results",
                 timeout_seconds: int = 600,
                 retries: int = 1,
                 backoff_factor: float = 2.0,
                 executor: Optional[Callable[..., subprocess.CompletedProcess]] = None):
        self.rgi_path = rgi_path
        self.card_db = card_db
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True, parents=True)
        self.logger = logging.getLogger("CARDRGIRunner")
        self.timeout_seconds = timeout_seconds
        self.retries = max(1, retries)
        self.backoff_factor = backoff_factor
        # Executor is primarily for testing; default to subprocess.run
        self._exec: Callable[..., subprocess.CompletedProcess] = executor or subprocess.run

    def run_rgi(self, fasta_path: str, sample_id: str, output_format: str = "txt") -> Optional[str]:
        """Run RGI on a single FASTA file. Returns path to output .txt/.csv or None on error."""
        # output_format: "txt" (tabular), "csv", or "json"
        if output_format not in ("txt", "csv", "json"):
            raise ValueError("output_format must be 'txt', 'csv', or 'json'")
        # Validate input
        if not fasta_path or not os.path.exists(fasta_path):
            self.logger.error(f"Input FASTA not found: {fasta_path}")
            return None
        if os.path.getsize(fasta_path) == 0:
            self.logger.error(f"Input FASTA is empty: {fasta_path}")
            return None

        ext = output_format if output_format != "txt" else "txt"
        output_file = self.output_dir / f"{sample_id}_rgi.{ext}"
        cmd = [
            self.rgi_path, "main",
            "--input_sequence", fasta_path,
            "--output_file", str(output_file),
            "--input_type", "contig",
            "--local",
            "--output_format", output_format
        ]
        if self.card_db:
            cmd += ["--card_json", self.card_db]

        attempt = 0
        delay = 0.0
        while attempt < self.retries:
            if delay > 0:
                time.sleep(delay)
            try:
                result = self._exec(cmd, capture_output=True, text=True, check=True, timeout=self.timeout_seconds)
                # Validate output file exists and looks sane
                if not os.path.exists(output_file):
                    raise RuntimeError("RGI did not produce output file")
                if os.path.getsize(output_file) == 0:
                    raise RuntimeError("RGI output file is empty")
                # Basic validation by format
                if output_format in ("txt", "csv"):
                    if not self._validate_tabular(str(output_file)):
                        raise RuntimeError("RGI tabular output failed validation")
                elif output_format == "json":
                    if not self._validate_json(str(output_file)):
                        raise RuntimeError("RGI JSON output failed validation")
                self.logger.info(f"RGI run complete for {sample_id}: {output_file}")
                return str(output_file)
            except subprocess.TimeoutExpired as e:
                self.logger.error(f"RGI timeout for {sample_id} after {self.timeout_seconds}s")
            except subprocess.CalledProcessError as e:
                # Include a snippet of stderr for debugging
                stderr_snippet = (e.stderr or "").strip().splitlines()[-1:] if e.stderr else []
                self.logger.error(f"RGI failed for {sample_id} (code {e.returncode}): {' '.join(stderr_snippet)}")
            except FileNotFoundError:
                self.logger.error(f"RGI binary not found: {self.rgi_path}")
                break
            except Exception as e:
                self.logger.error(f"RGI error for {sample_id}: {e}")
            attempt += 1
            delay = self.backoff_factor ** attempt
        return None

    def parse_rgi_output(self, rgi_json_path: str) -> List[Dict[str, Any]]:
        """Parse RGI JSON output for gene coordinates and annotations."""
        try:
            with open(rgi_json_path, 'r') as f:
                data = json.load(f)
            # Each entry is a hit with coordinates and annotation
            return data if isinstance(data, list) else []
        except Exception as e:
            self.logger.error(f"Failed to parse RGI output {rgi_json_path}: {e}")
            return []

    def run_batch(self, fasta_files: List[str], sample_ids: List[str], output_format: str = "txt") -> Dict[str, str]:
        """Run RGI on a batch of FASTA files. Returns mapping of sample_id to output file path."""
        results = {}
        for fasta, sid in zip(fasta_files, sample_ids):
            out = self.run_rgi(fasta, sid, output_format=output_format)
            if out:
                results[sid] = out
        return results

    def extract_gene_hits(self, rgi_tabular_path: str, gene_list: List[str]) -> List[Dict[str, str]]:
        """Extract only user-specified genes from RGI .txt/.csv output for FastaAAExtractor compatibility."""
        import csv
        hits = []
        try:
            with open(rgi_tabular_path, 'r', newline='') as f:
                reader = csv.DictReader(f, delimiter='\t' if rgi_tabular_path.endswith('.txt') else ',')
                for row in reader:
                    if row.get('Best_Hit_ARO') and row['Best_Hit_ARO'] in gene_list:
                        hits.append(row)
            return hits
        except Exception as e:
            self.logger.error(f"Failed to extract gene hits from {rgi_tabular_path}: {e}")
            return []

    # ---------------- Internal validators ----------------
    def _validate_tabular(self, file_path: str) -> bool:
        try:
            with open(file_path, 'r', newline='') as f:
                header = f.readline().strip()
                if not header:
                    return False
                # Expect at least these columns (RGI v5 typical columns include Best_Hit_ARO)
                required = ["Best_Hit_ARO"]
                return all(col in header for col in required)
        except Exception:
            return False

    def _validate_json(self, file_path: str) -> bool:
        try:
            with open(file_path, 'r') as f:
                data = json.load(f)
            # Accept list with at least 0 elements; stricter checks can be added
            return isinstance(data, list)
        except Exception:
            return False
