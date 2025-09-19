#!/usr/bin/env python3
"""
Abricate → Coordinates Adapter

Converts a raw Abricate TSV output into the coordinate CSV format expected by the
FastaAAExtractor integrator (drop-in replacement for legacy RGI coordinates).

Expected Abricate columns (tab-separated, comment lines starting with '#'):
  FILE, SEQUENCE, START, END, STRAND, GENE, COVERAGE, COVERAGE_MAP, GAPS,
  PERCENT_COVERAGE, PERCENT_IDENTITY, DATABASE, ACCESSION, PRODUCT, RESISTANCE

Output CSV columns (superset for compatibility):
  genome_id, accession, contig_id, start, end, strand, gene_name,
  cut_off, pass_bitscore, best_hit_aro, model_type, drug_class,
  resistance_mechanism, amr_gene_family, analysis_timestamp,
  rgi_version, card_version

Notes:
- genome_id is derived from the input TSV filename stem, with trailing
  "_abricate" removed when present.
- Unknown legacy RGI fields are populated with placeholder values to match
  downstream expectations without changing FastaAAExtractor.
- Empty inputs result in an empty CSV (with header only). Downstream
  orchestrator should avoid passing empty coordinate files to the extractor.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable, List
import sys
import io
import pandas as pd

ABR_COLS = [
    "FILE", "SEQUENCE", "START", "END", "STRAND", "GENE",
    "COVERAGE", "COVERAGE_MAP", "GAPS", "PERCENT_COVERAGE",
    "PERCENT_IDENTITY", "DATABASE", "ACCESSION", "PRODUCT", "RESISTANCE"
]

OUT_COLS = [
    "genome_id", "accession", "contig_id", "start", "end", "strand", "gene_name",
    "cut_off", "pass_bitscore", "best_hit_aro", "model_type", "drug_class",
    "resistance_mechanism", "amr_gene_family", "analysis_timestamp",
    "rgi_version", "card_version"
]


def _derive_genome_id_from_tsv(tsv_path: Path) -> str:
    stem = tsv_path.stem
    # Strip a trailing "_abricate" if present
    if stem.endswith("_abricate"):
        stem = stem[: -len("_abricate")]
    return stem


def _read_abricate_tsv(tsv_path: Path) -> pd.DataFrame:
    # Read, skipping comment lines starting with '#'
    with tsv_path.open("r", encoding="utf-8", errors="ignore") as f:
        lines = [ln for ln in f if not ln.lstrip().startswith("#")]
    if not lines:
        # return empty DataFrame with expected columns to keep downstream consistent
        return pd.DataFrame(columns=ABR_COLS)

    # Pandas can read from a string buffer
    buf = io.StringIO("".join(lines))
    try:
        df = pd.read_csv(buf, sep="\t", header=0, dtype=str)
    except Exception as e:
        raise ValueError(f"Failed to read Abricate TSV '{tsv_path}': {e}")

    # Normalize columns: upper-case, strip spaces, map known aliases
    col_map = {c: c.strip().upper() for c in df.columns}
    df.columns = [col_map[c] for c in df.columns]

    # Ensure required columns exist (fill if missing)
    for req in ABR_COLS:
        if req not in df.columns:
            df[req] = pd.Series([None] * len(df), dtype="object")

    return df[ABR_COLS]


def convert_abricate_to_coords(in_tsv_path: str | Path, out_csv_path: str | Path) -> int:
    """Convert a raw Abricate TSV file to a coordinate CSV.

    Returns the number of rows written to the output CSV (0 if empty).
    """
    in_path = Path(in_tsv_path)
    out_path = Path(out_csv_path)

    if not in_path.exists():
        raise FileNotFoundError(f"Input TSV not found: {in_path}")

    genome_id = _derive_genome_id_from_tsv(in_path)
    abr = _read_abricate_tsv(in_path)

    if abr.empty:
        # Write empty CSV with header to keep behavior explicit
        empty_out = pd.DataFrame(columns=OUT_COLS)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        empty_out.to_csv(out_path, index=False)
        print(f"[abricate_to_coords] No entries in {in_path.name}; wrote empty coordinates CSV.")
        return 0

    # Build output DataFrame with exact schema
    out = pd.DataFrame({
        "genome_id": genome_id,
        "accession": genome_id,  # legacy name; align with expectations
        "contig_id": abr["SEQUENCE"].fillna("").astype(str),
        "start": pd.to_numeric(abr["START"], errors="coerce").astype("Int64"),
        "end": pd.to_numeric(abr["END"], errors="coerce").astype("Int64"),
        "strand": abr["STRAND"].fillna("").astype(str),
        "gene_name": abr["GENE"].fillna("").astype(str),
        # Placeholders to match legacy RGI schema
        "cut_off": "Perfect",
        "pass_bitscore": 500.0,
        "best_hit_aro": "unknown",
        "model_type": "unknown",
        "drug_class": "unknown",
        "resistance_mechanism": "unknown",
        "amr_gene_family": "unknown",
        "analysis_timestamp": "unknown",
        "rgi_version": "unknown",
        "card_version": "unknown",
    })

    # Ensure correct column order
    out = out[OUT_COLS]

    # Drop rows missing critical coordinates
    before = len(out)
    out = out.dropna(subset=["start", "end", "gene_name", "contig_id"]).copy()
    after = len(out)
    if after < before:
        print(f"[abricate_to_coords] Dropped {before - after} row(s) with incomplete coordinates from {in_path.name}.")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_path, index=False)
    print(f"[abricate_to_coords] Wrote {after} coordinate row(s) -> {out_path}")
    return int(after)


def _parse_args(argv: Iterable[str]) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Convert Abricate TSV to coordinate CSV (FastaAAExtractor-compatible)")
    p.add_argument("--in-tsv", required=True, help="Path to input *_abricate.tsv file")
    p.add_argument("--out-csv", required=True, help="Path to output coordinates .csv file")
    return p.parse_args(list(argv))


def main(argv: Iterable[str] | None = None) -> int:
    args = _parse_args(sys.argv[1:] if argv is None else argv)
    try:
        convert_abricate_to_coords(args.in_tsv, args.out_csv)
        return 0
    except SystemExit:
        raise
    except Exception as e:
        sys.stderr.write(f"[abricate_to_coords] Error: {e}\n")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
