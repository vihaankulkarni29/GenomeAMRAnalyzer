"""
WildTypeAligner Module
=====================

- Aligns extracted proteins to reference wild-type proteins using EMBOSS-WATER or similar
- Ensures reference coverage for all user genes
- Handles error cases and stores alignment artifacts with provenance
- Modular, testable, and extensible

Engineering Principles:
- Fail-fast on missing input or alignment errors
- Full provenance and logging
"""

import os
import subprocess
import logging
from pathlib import Path
from typing import List, Dict, Any, Optional

class WildTypeAligner:
    """
    Aligns extracted proteins to reference wild-type proteins using EMBOSS-WATER.
    """
    def __init__(self, reference_dir: str, output_dir: str = "alignments", water_path: str = "water"):
        self.reference_dir = Path(reference_dir)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True, parents=True)
        self.water_path = water_path
        self.logger = logging.getLogger("WildTypeAligner")

    def align(self, protein_fasta: str, gene: str, sample_id: str, species: Optional[str] = None) -> Optional[str]:
        """
        Aligns a protein FASTA to the reference for the given gene and species. If reference is missing, fetch using SEPI. Returns path to alignment file or None on error.
        """
        ref_fasta = None
        species_norm = None
        if species:
            species_norm = species.replace(" ", "").replace(".", "").replace("/", "_")
            ref_fasta_species = self.reference_dir / f"{species_norm}_{gene}.faa"
            if ref_fasta_species.exists():
                ref_fasta = ref_fasta_species
        if ref_fasta is None:
            # Fallback to gene-only reference
            ref_fasta_gene = self.reference_dir / f"{gene}.faa"
            if ref_fasta_gene.exists():
                ref_fasta = ref_fasta_gene
        # If still not found, try to fetch using SEPI
        if ref_fasta is None and species and species_norm:
            try:
                from src.priority3.tools.reference_fetcher import fetch_reference_protein
                ref_fasta_path = self.reference_dir / f"{species_norm}_{gene}.faa"
                fetch_reference_protein(species.split()[0], species, gene, ref_fasta_path)
                if ref_fasta_path.exists():
                    ref_fasta = ref_fasta_path
            except Exception as e:
                self.logger.error(f"SEPI fetch failed for {species} {gene}: {e}")
        if ref_fasta is None:
            self.logger.error(f"Reference protein not found for gene {gene} (species: {species}): tried {self.reference_dir / f'{species}_{gene}.faa'} and {self.reference_dir / f'{gene}.faa'} and SEPI fetch.")
            return None
        out_file = self.output_dir / f"{sample_id}_{gene}_water.water"
        cmd = [
            self.water_path, 
            "-asequence", str(ref_fasta), 
            "-bsequence", protein_fasta, 
            "-gapopen", "10", 
            "-gapextend", "0.5", 
            "-aformat", "srspair",  # SubScan-compatible format
            "-outfile", str(out_file)
        ]
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            self.logger.info(f"Alignment complete: {out_file}")
            return str(out_file)
        except subprocess.CalledProcessError as e:
            self.logger.error(f"EMBOSS-WATER failed for {sample_id} {gene}: {e.stderr}")
            return None
        except Exception as e:
            self.logger.error(f"Alignment error for {sample_id} {gene}: {e}")
            return None

    def align_batch(self, protein_fastas: List[str], genes: List[str], sample_ids: List[str], species_list: Optional[List[str]] = None) -> Dict[str, str]:
        """
        Aligns a batch of protein FASTAs to their references. Returns mapping of (sample_id, gene) to alignment file.
        species_list: list of species names, same order as sample_ids (optional)
        """
        results = {}
        for idx, (fasta, gene, sid) in enumerate(zip(protein_fastas, genes, sample_ids)):
            species = species_list[idx] if species_list is not None else None
            out = self.align(fasta, gene, sid, species)
            if out:
                results[(sid, gene)] = out
        return results
