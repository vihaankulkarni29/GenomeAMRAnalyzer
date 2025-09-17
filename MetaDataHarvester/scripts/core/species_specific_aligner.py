import sys; import os; sys.path.append(os.path.join(os.path.dirname(__file__), '\\\\\\\\\src')); from configuration_manager import config_manager
#!/usr/bin/env python3
"""
SpeciesSpecificAligner - Eliminates false positives by using species-appropriate references
Addresses critical flaw: comparing all bacteria to E. coli references

Author: MetaDataHarvester Pipeline
Version: 2.0 - Species-Specific Alignment
"""

import os
import sys
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
from dataclasses import dataclass
from Bio import SeqIO, AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import PairwiseAligner
import pandas as pd
import numpy as np
from collections import defaultdict


@dataclass
class AlignmentResult:
    """Results of species-specific alignment"""
    query_accession: str
    reference_species: str
    reference_protein: str
    alignment_score: float
    sequence_identity: float
    coverage: float
    phylogenetic_distance: float
    quality_score: float
    warnings: List[str]


@dataclass
class ReferenceDatabase:
    """Species-specific reference database"""
    genus: str
    references: Dict[str, SeqRecord]  # protein_name -> sequence
    metadata: Dict[str, Dict]  # protein_name -> metadata
    phylogenetic_relatives: List[str]  # related genera for fallback


class SpeciesSpecificAligner:
    """
    Aligns proteins to species-appropriate references instead of universal E. coli
    """

    def __init__(self, reference_db_path: str = "references", output_dir: str = "species_alignments"):
        self.reference_db_path = Path(reference_db_path)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)

        # Setup logging
        self.logger = logging.getLogger(__name__)

        # Load reference database
        self.reference_db = self._load_reference_database()

        # Setup aligner with optimized parameters
        self.aligner = self._setup_aligner()

    def _load_reference_database(self) -> Dict[str, ReferenceDatabase]:
        """Load species-specific reference database"""
        reference_db = {}

        if not self.reference_db_path.exists():
            self.logger.warning(f"Reference database path {self.reference_db_path} does not exist")
            return reference_db

        # Load each genus directory
        for genus_dir in self.reference_db_path.iterdir():
            if genus_dir.is_dir():
                genus = genus_dir.name
                references = {}
                metadata = {}

                # Load reference sequences
                for fasta_file in genus_dir.glob("*.fasta"):
                    protein_name = fasta_file.stem.replace("_reference", "")

                    try:
                        records = list(SeqIO.parse(fasta_file, 'fasta'))
                        if records:
                            references[protein_name] = records[0]

                            # Load metadata if available
                            metadata_file = fasta_file.with_suffix('.json')
                            if metadata_file.exists():
                                import json
                                with open(metadata_file, 'r') as f:
                                    metadata[protein_name] = json.load(f)

                    except Exception as e:
                        self.logger.warning(f"Could not load reference {fasta_file}: {e}")

                if references:
                    # Determine phylogenetic relatives
                    phylogenetic_relatives = self._get_phylogenetic_relatives(genus)

                    reference_db[genus] = ReferenceDatabase(
                        genus=genus,
                        references=references,
                        metadata=metadata,
                        phylogenetic_relatives=phylogenetic_relatives
                    )

        self.logger.info(f"Loaded reference database for {len(reference_db)} genera")
        return reference_db

    def _get_phylogenetic_relatives(self, genus: str) -> List[str]:
        """Get phylogenetically related genera for fallback alignment"""
        phylogenetic_map = {
            'Escherichia': ['Shigella', 'Salmonella', 'Klebsiella', 'Enterobacter'],
            'Pseudomonas': ['Acinetobacter', 'Burkholderia', 'Stenotrophomonas'],
            'Klebsiella': ['Escherichia', 'Shigella', 'Salmonella', 'Enterobacter'],
            'Acinetobacter': ['Pseudomonas', 'Burkholderia', 'Stenotrophomonas'],
            'Salmonella': ['Escherichia', 'Shigella', 'Klebsiella', 'Enterobacter'],
            'Shigella': ['Escherichia', 'Salmonella', 'Klebsiella', 'Enterobacter']
        }

        return phylogenetic_map.get(genus, ['Escherichia'])  # Default fallback

    def _setup_aligner(self) -> PairwiseAligner:
        """Setup aligner with optimized parameters for protein alignment"""
        aligner = PairwiseAligner()
        aligner.mode = 'global'
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.open_gap_score = -2
        aligner.extend_gap_score = -0.5

        # Use BLOSUM62-like scoring for proteins
        aligner.substitution_matrix = self._get_blosum62_matrix()

        return aligner

    def _get_blosum62_matrix(self) -> Dict[Tuple[str, str], float]:
        """Get simplified BLOSUM62-like substitution matrix"""
        # Simplified matrix for common amino acids
        matrix = {}

        # High similarity groups
        similar_groups = [
            ['I', 'L', 'V', 'M'],  # Hydrophobic
            ['F', 'Y', 'W'],      # Aromatic
            ['D', 'E'],           # Acidic
            ['K', 'R'],           # Basic
            ['S', 'T'],           # Hydroxyl
            ['A', 'G']            # Small
        ]

        # Self matches
        for aa in 'ACDEFGHIKLMNPQRSTVWY':
            matrix[(aa, aa)] = 4

        # Similar amino acids
        for group in similar_groups:
            for aa1 in group:
                for aa2 in group:
                    if aa1 != aa2:
                        matrix[(aa1, aa2)] = 1
                        matrix[(aa2, aa1)] = 1

        # Dissimilar amino acids (negative scores)
        dissimilar_pairs = [
            ('A', 'C'), ('A', 'D'), ('A', 'E'), ('A', 'F'), ('A', 'H'),
            ('A', 'I'), ('A', 'K'), ('A', 'L'), ('A', 'M'), ('A', 'N'),
            ('A', 'P'), ('A', 'Q'), ('A', 'R'), ('A', 'S'), ('A', 'T'),
            ('A', 'V'), ('A', 'W'), ('A', 'Y'), ('C', 'D'), ('C', 'E'),
            ('C', 'F'), ('C', 'G'), ('C', 'H'), ('C', 'I'), ('C', 'K'),
            ('C', 'L'), ('C', 'M'), ('C', 'N'), ('C', 'P'), ('C', 'Q'),
            ('C', 'R'), ('C', 'S'), ('C', 'T'), ('C', 'V'), ('C', 'W'),
            ('C', 'Y'), ('D', 'F'), ('D', 'G'), ('D', 'H'), ('D', 'I'),
            ('D', 'K'), ('D', 'L'), ('D', 'M'), ('D', 'P'), ('D', 'R'),
            ('D', 'S'), ('D', 'T'), ('D', 'V'), ('D', 'W'), ('D', 'Y'),
            ('E', 'F'), ('E', 'G'), ('E', 'H'), ('E', 'I'), ('E', 'K'),
            ('E', 'L'), ('E', 'M'), ('E', 'P'), ('E', 'R'), ('E', 'S'),
            ('E', 'T'), ('E', 'V'), ('E', 'W'), ('E', 'Y'), ('F', 'G'),
            ('F', 'H'), ('F', 'I'), ('F', 'K'), ('F', 'L'), ('F', 'M'),
            ('F', 'N'), ('F', 'P'), ('F', 'Q'), ('F', 'R'), ('F', 'S'),
            ('F', 'T'), ('F', 'V'), ('F', 'Y'), ('G', 'H'), ('G', 'I'),
            ('G', 'K'), ('G', 'L'), ('G', 'M'), ('G', 'N'), ('G', 'P'),
            ('G', 'Q'), ('G', 'R'), ('G', 'S'), ('G', 'T'), ('G', 'V'),
            ('G', 'W'), ('G', 'Y'), ('H', 'I'), ('H', 'K'), ('H', 'L'),
            ('H', 'M'), ('H', 'N'), ('H', 'P'), ('H', 'Q'), ('H', 'R'),
            ('H', 'S'), ('H', 'T'), ('H', 'V'), ('H', 'W'), ('H', 'Y'),
            ('I', 'K'), ('I', 'N'), ('I', 'P'), ('I', 'Q'), ('I', 'R'),
            ('I', 'S'), ('I', 'T'), ('I', 'W'), ('I', 'Y'), ('K', 'L'),
            ('K', 'M'), ('K', 'N'), ('K', 'P'), ('K', 'Q'), ('K', 'S'),
            ('K', 'T'), ('K', 'V'), ('K', 'W'), ('K', 'Y'), ('L', 'N'),
            ('L', 'P'), ('L', 'Q'), ('L', 'R'), ('L', 'S'), ('L', 'T'),
            ('L', 'W'), ('L', 'Y'), ('M', 'N'), ('M', 'P'), ('M', 'Q'),
            ('M', 'R'), ('M', 'S'), ('M', 'T'), ('M', 'W'), ('M', 'Y'),
            ('N', 'P'), ('N', 'Q'), ('N', 'R'), ('N', 'V'), ('N', 'W'),
            ('N', 'Y'), ('P', 'Q'), ('P', 'R'), ('P', 'S'), ('P', 'T'),
            ('P', 'V'), ('P', 'W'), ('P', 'Y'), ('Q', 'S'), ('Q', 'T'),
            ('Q', 'V'), ('Q', 'W'), ('Q', 'Y'), ('R', 'S'), ('R', 'T'),
            ('R', 'V'), ('R', 'W'), ('R', 'Y'), ('S', 'V'), ('S', 'W'),
            ('S', 'Y'), ('T', 'V'), ('T', 'W'), ('T', 'Y'), ('V', 'W'),
            ('V', 'Y'), ('W', 'Y')
        ]

        for aa1, aa2 in dissimilar_pairs:
            matrix[(aa1, aa2)] = -1
            matrix[(aa2, aa1)] = -1

        return matrix

    def align_protein_to_species_reference(self, protein_sequence: SeqRecord,
                                          metadata: Dict) -> AlignmentResult:
        """
        Align protein to species-appropriate reference

        Args:
            protein_sequence: Query protein sequence
            metadata: Metadata with genus, species, and protein family info

        Returns:
            Alignment result with quality metrics
        """
        genus = metadata.get('genus', 'unknown')
        protein_family = metadata.get('protein_family', 'unknown')

        # Find appropriate reference
        reference_info = self._find_reference(genus, protein_family)

        if not reference_info:
            return AlignmentResult(
                query_accession=protein_sequence.id,
                reference_species="unknown",
                reference_protein="unknown",
                alignment_score=0.0,
                sequence_identity=0.0,
                coverage=0.0,
                phylogenetic_distance=1.0,
                quality_score=0.0,
                warnings=["No suitable reference found"]
            )

        # Perform alignment
        reference_seq = reference_info['sequence']
        alignments = self.aligner.align(protein_sequence.seq, reference_seq.seq)

        if not alignments:
            return AlignmentResult(
                query_accession=protein_sequence.id,
                reference_species=reference_info['genus'],
                reference_protein=reference_info['protein_name'],
                alignment_score=0.0,
                sequence_identity=0.0,
                coverage=0.0,
                phylogenetic_distance=reference_info['phylogenetic_distance'],
                quality_score=0.0,
                warnings=["Alignment failed"]
            )

        # Get best alignment
        best_alignment = alignments[0]

        # Calculate alignment metrics
        alignment_score = best_alignment.score
        sequence_identity = self._calculate_identity(best_alignment)
        coverage = self._calculate_coverage(best_alignment, protein_sequence.seq, reference_seq.seq)
        phylogenetic_distance = reference_info['phylogenetic_distance']

        # Assess alignment quality
        quality_score, warnings = self._assess_alignment_quality(
            best_alignment, sequence_identity, coverage, phylogenetic_distance
        )

        return AlignmentResult(
            query_accession=protein_sequence.id,
            reference_species=reference_info['genus'],
            reference_protein=reference_info['protein_name'],
            alignment_score=alignment_score,
            sequence_identity=sequence_identity,
            coverage=coverage,
            phylogenetic_distance=phylogenetic_distance,
            quality_score=quality_score,
            warnings=warnings
        )

    def _find_reference(self, genus: str, protein_family: str) -> Optional[Dict]:
        """Find best reference for given genus and protein family"""
        # Priority 1: Exact genus match
        if genus in self.reference_db:
            genus_db = self.reference_db[genus]
            # Look for any reference from this genus that matches the protein family
            for ref_protein, sequence in genus_db.references.items():
                if self._protein_matches_family(ref_protein, protein_family):
                    return {
                        'genus': genus,
                        'protein_name': ref_protein,
                        'sequence': sequence,
                        'phylogenetic_distance': 0.0,
                        'selection_method': 'exact_match'
                    }

        # Priority 2: Same protein family from phylogenetically close genus
        for close_genus in self._get_phylogenetic_relatives(genus):
            if close_genus in self.reference_db:
                close_genus_db = self.reference_db[close_genus]
                for ref_protein, sequence in close_genus_db.references.items():
                    if self._protein_matches_family(ref_protein, protein_family):
                        distance = self._calculate_phylogenetic_distance(genus, close_genus)
                        return {
                            'genus': close_genus,
                            'protein_name': ref_protein,
                            'sequence': sequence,
                            'phylogenetic_distance': distance,
                            'selection_method': 'phylogenetic_match'
                        }

        # Priority 3: Any reference from same genus (different protein)
        if genus in self.reference_db:
            genus_db = self.reference_db[genus]
            if genus_db.references:
                # Use first available reference
                ref_protein = list(genus_db.references.keys())[0]
                return {
                    'genus': genus,
                    'protein_name': ref_protein,
                    'sequence': genus_db.references[ref_protein],
                    'phylogenetic_distance': 0.5,  # Same genus, different protein
                    'selection_method': 'genus_fallback'
                }

        # Priority 4: No suitable reference found
        return None

    def _calculate_phylogenetic_distance(self, genus1: str, genus2: str) -> float:
        """Calculate phylogenetic distance between genera"""
        # Simplified distance matrix
        distance_matrix = {
            ('Escherichia', 'Shigella'): 0.1,
            ('Escherichia', 'Salmonella'): 0.2,
            ('Escherichia', 'Klebsiella'): 0.3,
            ('Escherichia', 'Enterobacter'): 0.3,
            ('Pseudomonas', 'Acinetobacter'): 0.4,
            ('Pseudomonas', 'Burkholderia'): 0.5,
        }

        key = tuple(sorted([genus1, genus2]))
        return distance_matrix.get(key, 0.8)  # Default distance

    def _protein_matches_family(self, protein_name: str, target_family: str) -> bool:
        """Check if a protein belongs to the target RND family"""
        protein_lower = protein_name.lower()
        family_lower = target_family.lower()

        # RND family mapping
        rnd_families = {
            'acra': ['acra', 'acr_a', 'acra_MG1655'],
            'acrb': ['acrb', 'acr_b', 'acrb_MG1655'],
            'tolc': ['tolc', config_manager.get_default_genes("rnd_efflux_pumps", "primary")[2], 'tolc_MG1655'],
            'acrd': ['acrd', 'acr_d'],
            'acre': ['acre', 'acr_e'],
            'acrf': ['acrf', 'acr_f'],
            'mdtb': ['mdtb', 'mdtB'],
            'mdtc': ['mdtc', 'mdtC'],
            'oqxa': ['oqxa', 'oqxA'],
            'oqxb': ['oqxb', 'oqxB'],
            'acrr': ['acrr', 'acrR'],
            'acrz': ['acrz', 'acrZ'],
            'mara': ['mara', 'marA'],
            'ramA': ['rama', 'ramA'],
            'soxs': ['soxs', 'soxS'],
            'rob': ['rob', 'rob'],
            'eefa': ['eefa', 'eefA'],
            'eefb': ['eefb', 'eefB'],
            'eefc': ['eefc', 'eefC']
        }

        if family_lower in rnd_families:
            return any(identifier in protein_lower for identifier in rnd_families[family_lower])

        return False

    def _calculate_identity(self, alignment) -> float:
        """Calculate sequence identity from alignment"""
        aligned_seq1 = str(alignment[0])
        aligned_seq2 = str(alignment[1])

        matches = 0
        total = 0

        for aa1, aa2 in zip(aligned_seq1, aligned_seq2):
            if aa1 != '-' and aa2 != '-':
                total += 1
                if aa1 == aa2:
                    matches += 1

        return matches / total if total > 0 else 0.0

    def _calculate_coverage(self, alignment, query_seq, reference_seq) -> float:
        """Calculate alignment coverage"""
        aligned_query = str(alignment[0])
        aligned_ref = str(alignment[1])

        # Count non-gap positions in query
        query_aligned_length = sum(1 for aa in aligned_query if aa != '-')
        query_total_length = len(query_seq)

        return query_aligned_length / query_total_length if query_total_length > 0 else 0.0

    def _assess_alignment_quality(self, alignment, identity: float, coverage: float,
                                phylogenetic_distance: float) -> Tuple[float, List[str]]:
        """Assess overall alignment quality"""
        warnings = []
        quality_score = 0.0

        # Identity score (0-40 points)
        if identity >= 0.8:
            quality_score += 40
        elif identity >= 0.6:
            quality_score += 30
        elif identity >= 0.4:
            quality_score += 20
        elif identity >= 0.2:
            quality_score += 10
        else:
            warnings.append(".1f")

        # Coverage score (0-40 points)
        if coverage >= 0.9:
            quality_score += 40
        elif coverage >= 0.7:
            quality_score += 30
        elif coverage >= 0.5:
            quality_score += 20
        else:
            warnings.append(".1f")

        # Phylogenetic distance penalty (0-20 points)
        if phylogenetic_distance <= 0.2:
            quality_score += 20
        elif phylogenetic_distance <= 0.5:
            quality_score += 15
        elif phylogenetic_distance <= 0.8:
            quality_score += 10
        else:
            warnings.append(".1f")

        # Normalize to 0-1 scale
        quality_score = quality_score / 100.0

        return quality_score, warnings

    def align_protein_dataset(self, fasta_file: str, metadata_file: str) -> List[AlignmentResult]:
        """
        Align entire protein dataset using species-specific references

        Args:
            fasta_file: FASTA file with protein sequences
            metadata_file: CSV file with enhanced metadata

        Returns:
            List of alignment results
        """
        # Load sequences
        sequences = list(SeqIO.parse(fasta_file, 'fasta'))

        # Load metadata
        metadata_df = pd.read_csv(metadata_file)
        metadata_dict = metadata_df.set_index('accession').to_dict('index')

        results = []

        for seq in sequences:
            accession = seq.id.split('.')[0]  # Remove version number
            seq_metadata = metadata_dict.get(accession, {})

            result = self.align_protein_to_species_reference(seq, seq_metadata)
            results.append(result)

            self.logger.info(
                f"Aligned {accession}: {result.reference_species} {result.reference_protein} "
                f"(quality: {result.quality_score:.2f})"
            )

        return results

    def save_alignment_results(self, results: List[AlignmentResult],
                              output_file: str = "alignment_results.csv") -> None:
        """Save alignment results to CSV"""
        output_path = self.output_dir / output_file

        data = []
        for result in results:
            row = {
                'query_accession': result.query_accession,
                'reference_species': result.reference_species,
                'reference_protein': result.reference_protein,
                'alignment_score': result.alignment_score,
                'sequence_identity': result.sequence_identity,
                'coverage': result.coverage,
                'phylogenetic_distance': result.phylogenetic_distance,
                'quality_score': result.quality_score,
                'warnings': '; '.join(result.warnings)
            }
            data.append(row)

        df = pd.DataFrame(data)
        df.to_csv(output_path, index=False)

        self.logger.info(f"Saved alignment results for {len(results)} proteins to {output_path}")

    def generate_alignment_report(self, results: List[AlignmentResult]) -> None:
        """Generate comprehensive alignment quality report"""
        report_file = self.output_dir / "alignment_quality_report.txt"

        with open(report_file, 'w') as f:
            f.write("SPECIES-SPECIFIC ALIGNMENT QUALITY REPORT\n")
            f.write("=" * 50 + "\n\n")

            f.write(f"Total alignments: {len(results)}\n\n")

            # Quality distribution
            quality_scores = [r.quality_score for r in results]
            high_quality = sum(1 for q in quality_scores if q >= 0.8)
            medium_quality = sum(1 for q in quality_scores if 0.6 <= q < 0.8)
            low_quality = sum(1 for q in quality_scores if q < 0.6)

            f.write("ALIGNMENT QUALITY DISTRIBUTION:\n")
            f.write("-" * 35 + "\n")
            f.write(f"High quality (≥0.8): {high_quality}\n")
            f.write(f"Medium quality (0.6-0.8): {medium_quality}\n")
            f.write(f"Low quality (<0.6): {low_quality}\n")
            f.write(f"Average quality score: {sum(quality_scores) / len(quality_scores):.1f}\n")
            f.write("\n")

            # Reference species distribution
            species_counts = {}
            for result in results:
                species = result.reference_species
                species_counts[species] = species_counts.get(species, 0) + 1

            f.write("REFERENCE SPECIES DISTRIBUTION:\n")
            f.write("-" * 35 + "\n")
            for species, count in sorted(species_counts.items()):
                f.write(f"{species}: {count}\n")
            f.write("\n")

            # Warnings summary
            all_warnings = []
            for result in results:
                all_warnings.extend(result.warnings)

            warning_counts = {}
            for warning in all_warnings:
                warning_counts[warning] = warning_counts.get(warning, 0) + 1

            if warning_counts:
                f.write("ALIGNMENT WARNINGS:\n")
                f.write("-" * 20 + "\n")
                for warning, count in sorted(warning_counts.items()):
                    f.write(f"{warning}: {count}\n")

        self.logger.info(f"Generated alignment report: {report_file}")


def main():
    """Main entry point"""
    import argparse

    parser = argparse.ArgumentParser(
        description="Species-specific protein alignment system"
    )

    parser.add_argument(
        "--fasta-file",
        required=True,
        help="FASTA file containing protein sequences"
    )

    parser.add_argument(
        "--metadata-file",
        required=True,
        help="CSV file with enhanced metadata"
    )

    parser.add_argument(
        "--reference-db",
        default="references",
        help="Path to reference database directory"
    )

    parser.add_argument(
        "--output-dir",
        default="species_alignments",
        help="Output directory for results"
    )

    args = parser.parse_args()

    # Setup logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    # Initialize aligner
    aligner = SpeciesSpecificAligner(args.reference_db, args.output_dir)

    # Perform alignments
    results = aligner.align_protein_dataset(args.fasta_file, args.metadata_file)

    # Save results
    aligner.save_alignment_results(results)
    aligner.generate_alignment_report(results)

    print(f"\nSpecies-specific alignment complete!")
    print(f"Processed {len(results)} proteins")
    print(f"Results saved to: {args.output_dir}")


if __name__ == "__main__":
    main()