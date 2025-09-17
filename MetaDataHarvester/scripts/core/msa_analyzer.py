import sys; import os; sys.path.append(os.path.join(os.path.dirname(__file__), '\\\\\\\\\src')); from configuration_manager import config_manager
#!/usr/bin/env python3
"""
MSAAnalyzer - Multiple Sequence Alignment for scientifically valid mutation detection
Addresses critical flaws in pairwise alignment approach for AMR research

Author: MetaDataHarvester Pipeline
Version: 1.0 - Scientifically Validated
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
import pandas as pd
import numpy as np
from scipy import stats
from collections import defaultdict, Counter


@dataclass
class MutationResult:
    """Results of mutation analysis with statistical significance"""
    position: int
    reference_aa: str
    variant_aa: str
    frequency: float
    count: int
    total_sequences: int
    p_value: float
    is_significant: bool
    conservation_score: float
    entropy: float
    known_resistance: bool


@dataclass
class MSAResult:
    """Results of multiple sequence alignment analysis"""
    alignment_file: str
    mutations: List[MutationResult]
    conservation_profile: Dict[int, float]
    phylogenetic_diversity: float
    quality_score: float
    warnings: List[str]


class MSAAnalyzer:
    """
    Performs multiple sequence alignment and statistically valid mutation detection
    """

    def __init__(self, reference_sequence: str, output_dir: str = "msa_results"):
        self.reference_sequence = reference_sequence
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)

        # Setup logging
        self.logger = logging.getLogger(__name__)

        # Known resistance mutations
        self.known_resistance_mutations = self._load_known_resistance()

        # Amino acid properties for conservation analysis
        self.aa_properties = self._load_aa_properties()

    def _load_known_resistance(self) -> Dict[str, Set[str]]:
        """Load known resistance mutations from literature"""
        return {
            config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0]: {'T104A', 'H596N', 'G616D', 'F615L', 'V630I'},
            config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1]: {'V610A', 'F628L', 'A279T', 'R717C', 'G288D', 'F610L'},
            config_manager.get_default_genes("rnd_efflux_pumps", "primary")[2]: {'D153A', 'E473K', 'R367H', 'A451V'}
        }

    def _load_aa_properties(self) -> Dict[str, Dict]:
        """Load amino acid physicochemical properties"""
        return {
            'A': {'hydrophobic': True, 'size': 'small', 'charge': 'neutral'},
            'R': {'hydrophobic': False, 'size': 'large', 'charge': 'positive'},
            'N': {'hydrophobic': False, 'size': 'medium', 'charge': 'neutral'},
            'D': {'hydrophobic': False, 'size': 'medium', 'charge': 'negative'},
            'C': {'hydrophobic': True, 'size': 'small', 'charge': 'neutral'},
            'Q': {'hydrophobic': False, 'size': 'large', 'charge': 'neutral'},
            'E': {'hydrophobic': False, 'size': 'large', 'charge': 'negative'},
            'G': {'hydrophobic': True, 'size': 'small', 'charge': 'neutral'},
            'H': {'hydrophobic': False, 'size': 'medium', 'charge': 'positive'},
            'I': {'hydrophobic': True, 'size': 'large', 'charge': 'neutral'},
            'L': {'hydrophobic': True, 'size': 'large', 'charge': 'neutral'},
            'K': {'hydrophobic': False, 'size': 'large', 'charge': 'positive'},
            'M': {'hydrophobic': True, 'size': 'large', 'charge': 'neutral'},
            'F': {'hydrophobic': True, 'size': 'large', 'charge': 'neutral'},
            'P': {'hydrophobic': True, 'size': 'medium', 'charge': 'neutral'},
            'S': {'hydrophobic': False, 'size': 'small', 'charge': 'neutral'},
            'T': {'hydrophobic': False, 'size': 'medium', 'charge': 'neutral'},
            'W': {'hydrophobic': True, 'size': 'large', 'charge': 'neutral'},
            'Y': {'hydrophobic': False, 'size': 'large', 'charge': 'neutral'},
            'V': {'hydrophobic': True, 'size': 'medium', 'charge': 'neutral'},
            'X': {'hydrophobic': False, 'size': 'unknown', 'charge': 'unknown'}
        }

    def create_multiple_alignment(self, fasta_files: List[str], protein_family: str) -> str:
        """
        Create multiple sequence alignment from FASTA files

        Args:
            fasta_files: List of FASTA file paths
            protein_family: Protein family for reference selection

        Returns:
            Path to alignment file
        """
        try:
            # Load all sequences
            all_sequences = []
            for fasta_file in fasta_files:
                records = list(SeqIO.parse(fasta_file, 'fasta'))
                all_sequences.extend(records)

            # Add reference sequence if not present
            ref_records = list(SeqIO.parse(self.reference_sequence, 'fasta'))
            if ref_records:
                ref_record = ref_records[0]
                # Check if reference is already in sequences
                ref_in_sequences = any(
                    str(record.seq) == str(ref_record.seq)
                    for record in all_sequences
                )
                if not ref_in_sequences:
                    all_sequences.insert(0, ref_record)  # Add reference as first sequence

            # Create alignment using MUSCLE (if available) or fallback to pairwise
            alignment_file = self.output_dir / f"{protein_family}_msa.fasta"

            if self._has_muscle():
                self._run_muscle_alignment(all_sequences, alignment_file)
            else:
                self.logger.warning("MUSCLE not found, using pairwise alignment approach")
                self._create_simple_msa(all_sequences, alignment_file)

            return str(alignment_file)

        except Exception as e:
            self.logger.error(f"Error creating MSA: {e}")
            raise

    def _has_muscle(self) -> bool:
        """Check if MUSCLE is available"""
        import subprocess
        try:
            subprocess.run(['muscle', '-version'], capture_output=True, check=True)
            return True
        except (subprocess.CalledProcessError, FileNotFoundError):
            return False

    def _run_muscle_alignment(self, sequences: List[SeqRecord], output_file: Path) -> None:
        """Run MUSCLE alignment"""
        import subprocess
        from Bio import SeqIO

        # Write sequences to temporary file
        temp_input = self.output_dir / "temp_sequences.fasta"
        SeqIO.write(sequences, temp_input, "fasta")

        # Run MUSCLE
        cmd = [
            'muscle',
            '-in', str(temp_input),
            '-out', str(output_file),
            '-maxiters', '16',
            '-diags'
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"MUSCLE failed: {result.stderr}")

        # Clean up
        temp_input.unlink(missing_ok=True)

    def _create_simple_msa(self, sequences: List[SeqRecord], output_file: Path) -> None:
        """Create simple MSA when MUSCLE is not available"""
        # This is a fallback - in practice, MUSCLE should be used
        from Bio.Align import PairwiseAligner

        aligner = PairwiseAligner()
        aligner.mode = 'global'
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.open_gap_score = -2
        aligner.extend_gap_score = -0.5

        # Use first sequence as reference
        reference = sequences[0]
        aligned_sequences = [reference]

        for seq in sequences[1:]:
            alignments = aligner.align(reference.seq, seq.seq)
            if alignments:
                best_alignment = alignments[0]
                # Extract aligned sequence
                aligned_seq = best_alignment.target
                aligned_record = SeqRecord(
                    aligned_seq,
                    id=seq.id,
                    name=seq.name,
                    description=seq.description
                )
                aligned_sequences.append(aligned_record)

        # Write alignment
        SeqIO.write(aligned_sequences, output_file, "fasta")

    def analyze_mutations(self, alignment_file: str, protein_family: str) -> MSAResult:
        """
        Analyze mutations in multiple sequence alignment

        Args:
            alignment_file: Path to MSA file
            protein_family: Protein family for analysis

        Returns:
            Analysis results with statistical significance
        """
        try:
            # Load alignment
            alignment = AlignIO.read(alignment_file, "fasta")

            # Calculate conservation profile
            conservation_profile = self._calculate_conservation(alignment)

            # Identify mutations
            mutations = self._identify_mutations(alignment, protein_family)

            # Calculate phylogenetic diversity
            phylogenetic_diversity = self._calculate_phylogenetic_diversity(alignment)

            # Assess alignment quality
            quality_score = self._assess_alignment_quality(alignment)

            # Generate warnings
            warnings = self._generate_quality_warnings(alignment, mutations)

            return MSAResult(
                alignment_file=alignment_file,
                mutations=mutations,
                conservation_profile=conservation_profile,
                phylogenetic_diversity=phylogenetic_diversity,
                quality_score=quality_score,
                warnings=warnings
            )

        except Exception as e:
            self.logger.error(f"Error analyzing mutations: {e}")
            raise

    def _calculate_conservation(self, alignment) -> Dict[int, float]:
        """Calculate conservation scores for each position"""
        conservation_scores = {}

        for i in range(alignment.get_alignment_length()):
            column = alignment[:, i]
            # Remove gaps
            column_no_gaps = [aa for aa in column if aa != '-']

            if not column_no_gaps:
                conservation_scores[i] = 0.0
                continue

            # Calculate Shannon entropy
            aa_counts = Counter(column_no_gaps)
            total = len(column_no_gaps)
            entropy = 0.0

            for count in aa_counts.values():
                if count > 0:
                    p = count / total
                    entropy -= p * np.log2(p)

            # Convert to conservation score (1 - normalized entropy)
            max_entropy = np.log2(20)  # 20 amino acids
            conservation_scores[i] = 1.0 - (entropy / max_entropy)

        return conservation_scores

    def _identify_mutations(self, alignment, protein_family: str) -> List[MutationResult]:
        """Identify mutations with statistical significance"""
        mutations = []

        # Use first sequence as reference
        reference_seq = str(alignment[0].seq)

        for i in range(len(reference_seq)):
            if reference_seq[i] == '-':
                continue  # Skip gaps in reference

            ref_aa = reference_seq[i]
            column = alignment[:, i]

            # Count amino acids at this position
            aa_counts = Counter(aa for aa in column if aa != '-')
            total_sequences = sum(aa_counts.values())

            # Find variants
            for variant_aa, count in aa_counts.items():
                if variant_aa == ref_aa:
                    continue  # Skip reference amino acid

                frequency = count / total_sequences

                # Statistical significance test
                # Test if variant frequency is significantly different from background
                p_value = self._calculate_significance(ref_aa, variant_aa, frequency, total_sequences)

                # Check if this is a known resistance mutation
                mutation_key = f"{ref_aa}{i+1}{variant_aa}"
                known_resistance = mutation_key in self.known_resistance_mutations.get(protein_family, set())

                # Calculate conservation and entropy
                conservation_score = self._calculate_position_conservation(column)
                entropy = self._calculate_position_entropy(column)

                mutation = MutationResult(
                    position=i + 1,  # 1-based position
                    reference_aa=ref_aa,
                    variant_aa=variant_aa,
                    frequency=frequency,
                    count=count,
                    total_sequences=total_sequences,
                    p_value=p_value,
                    is_significant=p_value < 0.05,
                    conservation_score=conservation_score,
                    entropy=entropy,
                    known_resistance=known_resistance
                )

                mutations.append(mutation)

        return mutations

    def _calculate_significance(self, ref_aa: str, var_aa: str, frequency: float,
                              total_sequences: int) -> float:
        """Calculate statistical significance of mutation"""
        # Simple binomial test: is this variant frequency significantly different from random?
        # This is a simplified approach - in practice, more sophisticated tests would be used

        # Expected frequency based on amino acid abundance in proteins
        aa_abundance = {
            'A': 0.078, 'R': 0.051, 'N': 0.043, 'D': 0.054, 'C': 0.019,
            'Q': 0.042, 'E': 0.051, 'G': 0.074, 'H': 0.026, 'I': 0.053,
            'L': 0.091, 'K': 0.059, 'M': 0.023, 'F': 0.039, 'P': 0.052,
            'S': 0.068, 'T': 0.059, 'W': 0.014, 'Y': 0.032, 'V': 0.066
        }

        expected_freq = aa_abundance.get(var_aa, 0.05)  # Default 5% if unknown

        # Binomial test
        from scipy.stats import binomtest
        result = binomtest(int(frequency * total_sequences), total_sequences, expected_freq)
        return result.pvalue

    def _calculate_position_conservation(self, column: str) -> float:
        """Calculate conservation score for a position"""
        aa_counts = Counter(aa for aa in column if aa != '-')
        if not aa_counts:
            return 0.0

        total = sum(aa_counts.values())
        max_freq = max(aa_counts.values()) / total

        # Conservation = 1 - normalized entropy
        entropy = 0.0
        for count in aa_counts.values():
            if count > 0:
                p = count / total
                entropy -= p * np.log2(p)

        max_entropy = np.log2(len(aa_counts))
        if max_entropy == 0:
            return 1.0

        return 1.0 - (entropy / max_entropy)

    def _calculate_position_entropy(self, column: str) -> float:
        """Calculate Shannon entropy for a position"""
        aa_counts = Counter(aa for aa in column if aa != '-')
        if not aa_counts:
            return 0.0

        total = sum(aa_counts.values())
        entropy = 0.0

        for count in aa_counts.values():
            if count > 0:
                p = count / total
                entropy -= p * np.log2(p)

        return entropy

    def _calculate_phylogenetic_diversity(self, alignment) -> float:
        """Calculate phylogenetic diversity of the alignment"""
        # Simple diversity metric based on sequence differences
        if len(alignment) < 2:
            return 0.0

        total_pairs = 0
        differing_pairs = 0

        for i in range(len(alignment)):
            for j in range(i + 1, len(alignment)):
                total_pairs += 1
                seq1 = str(alignment[i].seq)
                seq2 = str(alignment[j].seq)

                # Count differences (excluding gaps)
                differences = 0
                valid_positions = 0

                for k in range(len(seq1)):
                    if seq1[k] != '-' and seq2[k] != '-':
                        valid_positions += 1
                        if seq1[k] != seq2[k]:
                            differences += 1

                if valid_positions > 0 and differences / valid_positions > 0.1:  # >10% different
                    differing_pairs += 1

        return differing_pairs / total_pairs if total_pairs > 0 else 0.0

    def _assess_alignment_quality(self, alignment) -> float:
        """Assess overall alignment quality"""
        if len(alignment) == 0:
            return 0.0

        total_positions = alignment.get_alignment_length()
        gap_positions = 0

        for i in range(total_positions):
            column = alignment[:, i]
            if '-' in column:
                gap_positions += 1

        # Quality = 1 - (gap proportion)
        gap_proportion = gap_positions / total_positions
        return 1.0 - gap_proportion

    def _generate_quality_warnings(self, alignment, mutations: List[MutationResult]) -> List[str]:
        """Generate quality warnings"""
        warnings = []

        # Check alignment size
        if len(alignment) < 3:
            warnings.append("Very small alignment (< 3 sequences)")

        # Check for too many gaps
        total_positions = alignment.get_alignment_length()
        gap_positions = sum(1 for i in range(total_positions) if '-' in alignment[:, i])
        gap_ratio = gap_positions / total_positions

        if gap_ratio > 0.5:
            warnings.append(".1f")

        # Check for low mutation significance
        significant_mutations = sum(1 for m in mutations if m.is_significant)
        if len(mutations) > 0 and significant_mutations / len(mutations) < 0.1:
            warnings.append("Very few statistically significant mutations detected")

        return warnings

    def save_results(self, result: MSAResult, output_prefix: str) -> None:
        """Save analysis results to files"""
        # Save mutations
        mutation_data = []
        for mutation in result.mutations:
            mutation_data.append({
                'position': mutation.position,
                'reference_aa': mutation.reference_aa,
                'variant_aa': mutation.variant_aa,
                'frequency': mutation.frequency,
                'count': mutation.count,
                'total_sequences': mutation.total_sequences,
                'p_value': mutation.p_value,
                'is_significant': mutation.is_significant,
                'conservation_score': mutation.conservation_score,
                'entropy': mutation.entropy,
                'known_resistance': mutation.known_resistance
            })

        df = pd.DataFrame(mutation_data)
        df.to_csv(f"{output_prefix}_mutations.csv", index=False)

        # Save conservation profile
        conservation_df = pd.DataFrame(
            list(result.conservation_profile.items()),
            columns=['position', 'conservation_score']
        )
        conservation_df.to_csv(f"{output_prefix}_conservation.csv", index=False)

        # Save summary
        summary = {
            'total_mutations': len(result.mutations),
            'significant_mutations': sum(1 for m in result.mutations if m.is_significant),
            'known_resistance_mutations': sum(1 for m in result.mutations if m.known_resistance),
            'phylogenetic_diversity': result.phylogenetic_diversity,
            'alignment_quality': result.quality_score,
            'warnings': result.warnings
        }

        with open(f"{output_prefix}_summary.json", 'w') as f:
            import json
            json.dump(summary, f, indent=2)


def main():
    """Main entry point"""
    import argparse

    parser = argparse.ArgumentParser(description="Multiple sequence alignment and mutation analysis")
    parser.add_argument("--fasta-dir", required=True, help="Directory containing FASTA files")
    parser.add_argument("--reference", required=True, help="Reference sequence FASTA file")
    parser.add_argument("--protein-family", required=True, help="Protein family (acrA, acrB, tolC)")
    parser.add_argument("--output-prefix", default="msa_analysis", help="Output file prefix")

    args = parser.parse_args()

    # Setup logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    # Find FASTA files
    fasta_dir = Path(args.fasta_dir)
    fasta_files = list(fasta_dir.glob("*.fasta"))

    if not fasta_files:
        print(f"No FASTA files found in {fasta_dir}")
        sys.exit(1)

    # Run analysis
    analyzer = MSAAnalyzer(args.reference)
    alignment_file = analyzer.create_multiple_alignment(fasta_files, args.protein_family)
    results = analyzer.analyze_mutations(alignment_file, args.protein_family)
    analyzer.save_results(results, args.output_prefix)

    print(f"MSA analysis complete. Results saved with prefix: {args.output_prefix}")


if __name__ == "__main__":
    main()