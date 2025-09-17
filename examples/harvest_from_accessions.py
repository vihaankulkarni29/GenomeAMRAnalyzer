"""
Example: Harvest genomes from a list of assembly accessions (GCF_/GCA_)
This script reads accessions from a text file (one per line, comments with # allowed)
 and runs GenBankGenomeHarvester to download genomic FASTAs ready for CARD RGI.
"""
from pathlib import Path
import sys

# Ensure project root is on sys.path to import src.priority3
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
import argparse

from src.priority3.db.harvesters.genbank_genome_harvester import GenBankGenomeHarvester


def read_accessions(file_path: Path):
    accs = []
    for line in file_path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        accs.append(line)
    return accs


def main():
    p = argparse.ArgumentParser(description="Harvest genomes from assembly accessions for CARD RGI")
    p.add_argument("accessions_file", type=str, help="Path to text file with GCF_/GCA_ accessions")
    p.add_argument("--out", dest="out", default="genomes", help="Output directory for FASTA files")
    p.add_argument("--db", dest="db", default="priority3.db", help="SQLite DB path")
    p.add_argument("--api-key", dest="api_key", default=None, help="NCBI API key for higher rate limits")
    p.add_argument("--email", dest="email", default=None, help="Contact email for NCBI eutils")
    p.add_argument("--resume", dest="resume", action="store_true", help="Resume from checkpoint if present")
    p.add_argument("--mock", dest="mock", action="store_true", help="Run in mock mode (no network)")
    args = p.parse_args()

    accessions = read_accessions(Path(args.accessions_file))
    if not accessions:
        print("No accessions found in file")
        return 1

    harv = GenBankGenomeHarvester(output_dir=args.out, db_path=args.db, api_key=args.api_key, email=args.email, mock_mode=args.mock)
    try:
        summary = harv.harvest_by_accessions(accessions, resume=args.resume)
        print("Harvest summary:\n", summary)
        print(f"FASTAs saved in: {args.out}")
    finally:
        harv.close()


if __name__ == "__main__":
    raise SystemExit(main())
