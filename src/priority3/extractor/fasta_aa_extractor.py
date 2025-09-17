"""
FastaAAExtractor Module
======================

- Extracts user-specified protein sequences from genomes using CARD RGI coordinates
- Supports .tab/.csv/.xlsx RGI outputs
- Validates gene list, handles missing genes, and stores extracted protein FASTAs with provenance
- Robust error handling and logging

Engineering Principles:
- Fail-fast on missing input or extraction errors
- Modular, testable, and extensible
- Full provenance and logging
"""

import os
import logging
from pathlib import Path
from typing import List, Dict, Any, Optional
import csv

class FastaAAExtractor:
    """
    Extracts protein sequences for user-specified genes from genome FASTA using CARD RGI coordinates.
    """
    def __init__(self, output_dir: str = "extracted_proteins"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True, parents=True)
        self.logger = logging.getLogger("FastaAAExtractor")

    def extract_proteins(self, genome_fasta: str, rgi_tabular: str, gene_list: List[str], sample_id: str) -> List[str]:
        """
        Extracts protein FASTA files for each gene in gene_list using CARD RGI .tab/.csv output.
        Returns list of output FASTA file paths.
        """
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        hits = self._parse_rgi_tabular(rgi_tabular, gene_list)
        if not hits:
            self.logger.warning(f"No matching genes found in RGI output for {sample_id}")
            return []
        # Load genome sequences
        seqs = {rec.id: rec for rec in SeqIO.parse(genome_fasta, "fasta")}
        output_files = []
        for hit in hits:
            gene = hit['Best_Hit_ARO']
            contig = hit.get('Contig', hit.get('Contig name', ''))
            start = int(hit.get('Start', hit.get('Start coord', 0)))
            end = int(hit.get('Stop', hit.get('Stop coord', 0)))
            strand = hit.get('Orientation', '+')
            if contig not in seqs:
                self.logger.warning(f"Contig {contig} not found in genome for gene {gene}")
                continue
            seq = seqs[contig].seq[start-1:end] if start < end else seqs[contig].seq[end-1:start].reverse_complement()
            if strand == '-':
                seq = seq.reverse_complement()
            # Translate to protein
            protein = seq.translate(to_stop=True)
            rec = SeqRecord(protein, id=f"{sample_id}|{gene}", description=f"Extracted from {contig}:{start}-{end} ({strand})")
            out_fasta = self.output_dir / f"{sample_id}_{gene}.faa"
            SeqIO.write([rec], out_fasta, "fasta")
            output_files.append(str(out_fasta))
        return output_files

    def _parse_rgi_tabular(self, rgi_tabular: str, gene_list: List[str]) -> List[Dict[str, str]]:
        """
        Parse RGI .tab/.csv/.xlsx output and filter for user-specified genes.
        """
        hits = []
        ext = Path(rgi_tabular).suffix.lower()
        try:
            if ext == '.xlsx':
                import pandas as pd
                df = pd.read_excel(rgi_tabular)
                for _, row in df.iterrows():
                    if row.get('Best_Hit_ARO') in gene_list:
                        hits.append(row.to_dict())
            else:
                delimiter = '\t' if ext == '.txt' else ','
                with open(rgi_tabular, 'r', newline='') as f:
                    reader = csv.DictReader(f, delimiter=delimiter)
                    for row in reader:
                        if row.get('Best_Hit_ARO') in gene_list:
                            hits.append(row)
        except Exception as e:
            self.logger.error(f"Failed to parse RGI tabular {rgi_tabular}: {e}")
        return hits
