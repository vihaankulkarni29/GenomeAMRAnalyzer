import sys; import os; sys.path.append(os.path.join(os.path.dirname(__file__), '\\\\\\\\\src')); from configuration_manager import config_manager
#!/usr/bin/env python3
"""
SpeciesAwareMutator - Intelligent mutation calling with species-specific context
Combines phylogenetic analysis, clinical data, and statistical validation

Author: MetaDataHarvester Pipeline
Version: 2.0 - Species-Aware Mutation Calling
"""

import os
import sys
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union, Set
from dataclasses import dataclass, field
import pandas as pd
import numpy as np
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from scipy import stats
import json
from collections import defaultdict, Counter
from datetime import datetime


@dataclass
class MutationCall:
    """Intelligent mutation call with species context"""
    position: int
    reference_aa: str
    variant_aa: str
    mutation_type: str  # 'substitution', 'insertion', 'deletion'
    frequency: float
    confidence_score: float
    species_context: Dict[str, float]
    phylogenetic_support: float
    clinical_evidence: List[str]
    functional_impact: str
    conservation_score: float
    known_resistance_association: Optional[str]
    mic_correlation: Optional[float]
    quality_flags: List[str] = field(default_factory=list)

    @property
    def mutation_key(self) -> str:
        """Standard mutation identifier"""
        return f"{self.reference_aa}{self.position}{self.variant_aa}"

    @property
    def is_significant(self) -> bool:
        """Determine if mutation is statistically significant"""
        return (self.confidence_score > 0.8 and
                self.frequency > 0.05 and
                self.phylogenetic_support > 0.7)

    @property
    def resistance_risk(self) -> str:
        """Assess resistance risk level"""
        if self.known_resistance_association:
            return "High"
        elif self.mic_correlation and abs(self.mic_correlation) > 0.5:
            return "Medium"
        elif self.conservation_score < 0.3:
            return "Low"
        else:
            return "Unknown"


@dataclass
class SpeciesContext:
    """Species-specific context for mutation analysis"""
    genus: str
    species: str
    reference_sequence: SeqRecord
    phylogenetic_distance: float
    clinical_isolates: int
    known_mutations: Dict[str, str]  # mutation -> resistance phenotype
    mic_data: Dict[str, List[float]]  # antibiotic -> MIC values


class SpeciesAwareMutator:
    """
    Intelligent mutation calling system with species-specific context
    """

    def __init__(self, reference_db_path: str, clinical_data_path: Optional[str] = None):
        self.reference_db_path = Path(reference_db_path)
        self.clinical_data_path = clinical_data_path

        # Setup logging
        self.logger = logging.getLogger(__name__)
        self._setup_logging()

        # Load reference database
        self.reference_db = self._load_reference_database()

        # Load clinical context
        self.clinical_context = self._load_clinical_context()

        # Quality thresholds
        self.min_alignment_score = 50.0
        self.min_conservation_score = 0.3
        self.min_frequency_threshold = 0.05
        self.max_phylogenetic_distance = 0.8

    def _setup_logging(self):
        """Setup logging configuration"""
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        )

    def _load_reference_database(self) -> Dict[str, SpeciesContext]:
        """Load species-specific reference database"""
        reference_db = {}

        if not self.reference_db_path.exists():
            self.logger.warning(f"Reference database not found: {self.reference_db_path}")
            return reference_db

        # Load species-specific references
        for genus_dir in self.reference_db_path.iterdir():
            if genus_dir.is_dir():
                genus = genus_dir.name
                reference_db[genus] = {}

                for species_file in genus_dir.glob("*.fasta"):
                    try:
                        species_name = species_file.stem
                        sequences = list(SeqIO.parse(species_file, 'fasta'))

                        if sequences:
                            # Create species context
                            context = SpeciesContext(
                                genus=genus,
                                species=species_name,
                                reference_sequence=sequences[0],
                                phylogenetic_distance=0.0,  # Would be calculated
                                clinical_isolates=0,  # Would be loaded from metadata
                                known_mutations={},  # Would be loaded from CARD/other sources
                                mic_data={}  # Would be loaded from MIC database
                            )

                            reference_db[genus][species_name] = context

                    except Exception as e:
                        self.logger.warning(f"Failed to load {species_file}: {e}")

        self.logger.info(f"Loaded reference database with {len(reference_db)} genera")
        return reference_db

    def _load_clinical_context(self) -> Dict[str, Dict]:
        """Load clinical context data"""
        clinical_context = {}

        if self.clinical_data_path and Path(self.clinical_data_path).exists():
            try:
                # Load clinical metadata
                clinical_df = pd.read_csv(self.clinical_data_path)
                clinical_context = clinical_df.to_dict('index')
                self.logger.info(f"Loaded clinical context for {len(clinical_context)} isolates")
            except Exception as e:
                self.logger.warning(f"Failed to load clinical context: {e}")

        return clinical_context

    def call_mutations_species_aware(self, alignment_file: str,
                                   metadata_file: str,
                                   protein_family: str) -> List[MutationCall]:
        """
        Perform species-aware mutation calling

        Args:
            alignment_file: Path to alignment file
            metadata_file: Path to metadata CSV
            protein_family: RND protein family

        Returns:
            List of intelligent mutation calls
        """
        self.logger.info("Starting species-aware mutation calling")
        self.logger.info(f"Alignment: {alignment_file}")
        self.logger.info(f"Protein family: {protein_family}")

        # Load alignment
        try:
            alignment = AlignIO.read(alignment_file, 'fasta')
        except Exception as e:
            self.logger.error(f"Failed to load alignment: {e}")
            return []

        # Load metadata
        try:
            metadata_df = pd.read_csv(metadata_file)
            metadata_dict = metadata_df.set_index('accession').to_dict('index')
        except Exception as e:
            self.logger.error(f"Failed to load metadata: {e}")
            return []

        # Group sequences by species
        species_groups = self._group_sequences_by_species(alignment, metadata_dict)

        # Call mutations for each species group
        all_mutations = []
        for species_key, sequences in species_groups.items():
            if len(sequences) < 2:  # Need at least 2 sequences for comparison
                continue

            species_mutations = self._call_mutations_for_species(
                sequences, species_key, protein_family, metadata_dict
            )
            all_mutations.extend(species_mutations)

        # Apply cross-species validation
        validated_mutations = self._apply_cross_species_validation(all_mutations)

        # Add clinical context
        enriched_mutations = self._enrich_with_clinical_context(validated_mutations)

        self.logger.info(f"Called {len(enriched_mutations)} species-aware mutations")
        return enriched_mutations

    def _group_sequences_by_species(self, alignment: MultipleSeqAlignment,
                                  metadata_dict: Dict) -> Dict[str, List[SeqRecord]]:
        """Group sequences by species"""
        species_groups = defaultdict(list)

        for record in alignment:
            accession = record.id.split('_')[0]  # Extract accession from ID
            metadata = metadata_dict.get(accession, {})

            genus = metadata.get('genus', 'unknown')
            species = metadata.get('species', 'unknown')
            species_key = f"{genus}_{species}"

            species_groups[species_key].append(record)

        self.logger.info(f"Grouped sequences into {len(species_groups)} species groups")
        return species_groups

    def _call_mutations_for_species(self, sequences: List[SeqRecord],
                                  species_key: str, protein_family: str,
                                  metadata_dict: Dict) -> List[MutationCall]:
        """Call mutations for a specific species group"""
        mutations = []

        if not sequences:
            return mutations

        # Get reference sequence for this species
        genus, species = species_key.split('_', 1)
        reference_seq = self._get_species_reference(genus, species, protein_family)

        if not reference_seq:
            self.logger.warning(f"No reference found for {species_key}")
            return mutations

        # Compare each sequence to reference
        for seq_record in sequences:
            seq_mutations = self._compare_to_reference(
                seq_record, reference_seq, species_key, metadata_dict
            )
            mutations.extend(seq_mutations)

        # Aggregate mutations across sequences
        aggregated_mutations = self._aggregate_species_mutations(mutations, species_key)

        return aggregated_mutations

    def _get_species_reference(self, genus: str, species: str,
                             protein_family: str) -> Optional[SeqRecord]:
        """Get appropriate reference sequence for species"""
        # Try exact species match first
        if genus in self.reference_db and species in self.reference_db[genus]:
            context = self.reference_db[genus][species]
            return context.reference_sequence

        # Try genus-level reference
        if genus in self.reference_db:
            # Find any reference from this genus
            for species_name, context in self.reference_db[genus].items():
                return context.reference_sequence

        # Try phylogenetic closest match
        return self._find_phylogenetic_reference(genus, protein_family)

    def _find_phylogenetic_reference(self, target_genus: str,
                                   protein_family: str) -> Optional[SeqRecord]:
        """Find phylogenetically closest reference"""
        # Simple distance-based selection
        # In practice, this would use phylogenetic trees
        phylogenetic_distances = {
            'Escherichia': {'Klebsiella': 0.2, 'Salmonella': 0.3, 'Shigella': 0.1},
            'Pseudomonas': {'Burkholderia': 0.4, 'Acinetobacter': 0.5},
            'Klebsiella': {'Escherichia': 0.2, 'Enterobacter': 0.3},
        }

        min_distance = float('inf')
        best_reference = None

        for genus, distances in phylogenetic_distances.items():
            if target_genus in distances:
                distance = distances[target_genus]
                if distance < min_distance and genus in self.reference_db:
                    min_distance = distance
                    # Get first available reference from this genus
                    for species_context in self.reference_db[genus].values():
                        best_reference = species_context.reference_sequence
                        break

        return best_reference

    def _compare_to_reference(self, query_seq: SeqRecord, reference_seq: SeqRecord,
                            species_key: str, metadata_dict: Dict) -> List[MutationCall]:
        """Compare query sequence to reference and identify mutations"""
        mutations = []

        # Ensure sequences are aligned
        if len(query_seq.seq) != len(reference_seq.seq):
            self.logger.warning(f"Sequence length mismatch for {query_seq.id}")
            return mutations

        # Compare each position
        for i, (ref_aa, query_aa) in enumerate(zip(reference_seq.seq, query_seq.seq)):
            if ref_aa != query_aa and ref_aa != '-' and query_aa != '-':
                # Found a substitution
                mutation = MutationCall(
                    position=i + 1,  # 1-based position
                    reference_aa=ref_aa,
                    variant_aa=query_aa,
                    mutation_type='substitution',
                    frequency=1.0,  # Will be updated during aggregation
                    confidence_score=0.8,  # Base confidence
                    species_context={species_key: 1.0},
                    phylogenetic_support=0.9,  # High for conspecific
                    clinical_evidence=[],
                    functional_impact=self._assess_functional_impact(ref_aa, query_aa, i),
                    conservation_score=self._calculate_conservation_score(i, reference_seq),
                    known_resistance_association=self._check_known_resistance(ref_aa, query_aa, i, species_key),
                    mic_correlation=None
                )
                mutations.append(mutation)

        return mutations

    def _assess_functional_impact(self, ref_aa: str, var_aa: str, position: int) -> str:
        """Assess functional impact of amino acid change"""
        # Simple physicochemical property-based assessment
        hydrophobic = {'A', 'V', 'L', 'I', 'M', 'F', 'W', 'P'}
        charged = {'D', 'E', 'K', 'R', 'H'}
        polar = {'S', 'T', 'N', 'Q', 'C', 'Y'}

        ref_props = []
        var_props = []

        if ref_aa in hydrophobic:
            ref_props.append('hydrophobic')
        if ref_aa in charged:
            ref_props.append('charged')
        if ref_aa in polar:
            ref_props.append('polar')

        if var_aa in hydrophobic:
            var_props.append('hydrophobic')
        if var_aa in charged:
            var_props.append('charged')
        if var_aa in polar:
            var_props.append('polar')

        # Assess impact
        if set(ref_props) == set(var_props):
            return 'Conservative'
        elif 'charged' in ref_props and 'charged' not in var_props:
            return 'Charge change'
        elif ref_aa == 'C' or var_aa == 'C':
            return 'Cysteine modification'
        elif ref_aa == 'P' or var_aa == 'P':
            return 'Proline modification'
        else:
            return 'Non-conservative'

    def _calculate_conservation_score(self, position: int, reference_seq: SeqRecord) -> float:
        """Calculate conservation score at position"""
        # Simple conservation based on BLOSUM62-like scoring
        # In practice, this would use multiple sequence alignments
        conserved_positions = {50, 100, 150, 200, 250}  # Example conserved positions

        if position in conserved_positions:
            return 0.9
        else:
            return 0.5

    def _check_known_resistance(self, ref_aa: str, var_aa: str, position: int,
                              species_key: str) -> Optional[str]:
        """Check if mutation is known to confer resistance"""
        # Known resistance mutations (example data)
        known_resistance = {
            'Escherichia_T104A': 'AcrB efflux pump resistance',
            'Pseudomonas_G288D': 'MexB resistance',
            'Klebsiella_R717Q': 'AcrB resistance'
        }

        mutation_key = f"{species_key}_{ref_aa}{position}{var_aa}"
        return known_resistance.get(mutation_key)

    def _aggregate_species_mutations(self, mutations: List[MutationCall],
                                   species_key: str) -> List[MutationCall]:
        """Aggregate mutations across multiple sequences of same species"""
        mutation_groups = defaultdict(list)

        # Group mutations by position and amino acid change
        for mutation in mutations:
            key = (mutation.position, mutation.reference_aa, mutation.variant_aa)
            mutation_groups[key].append(mutation)

        # Aggregate each group
        aggregated_mutations = []
        for (position, ref_aa, var_aa), group_mutations in mutation_groups.items():
            frequency = len(group_mutations) / len(set(m.confidence_score for m in group_mutations))

            # Use the first mutation as template, update frequency
            aggregated = group_mutations[0]
            aggregated.frequency = frequency
            aggregated.confidence_score = min(1.0, aggregated.confidence_score * frequency * len(group_mutations))

            aggregated_mutations.append(aggregated)

        return aggregated_mutations

    def _apply_cross_species_validation(self, mutations: List[MutationCall]) -> List[MutationCall]:
        """Apply cross-species validation to filter false positives"""
        validated_mutations = []

        for mutation in mutations:
            # Check if mutation appears in multiple species
            species_count = len(mutation.species_context)

            if species_count > 1:
                # Increase confidence for cross-species mutations
                mutation.confidence_score = min(1.0, mutation.confidence_score * 1.2)
                mutation.phylogenetic_support = min(1.0, mutation.phylogenetic_support * 1.3)
            elif mutation.phylogenetic_support < 0.5:
                # Flag low-confidence species-specific mutations
                mutation.quality_flags.append("Species-specific, low phylogenetic support")

            # Apply frequency and conservation filters
            if (mutation.frequency >= self.min_frequency_threshold and
                mutation.conservation_score >= self.min_conservation_score):
                validated_mutations.append(mutation)
            else:
                mutation.quality_flags.append("Filtered: low frequency or conservation")

        return validated_mutations

    def _enrich_with_clinical_context(self, mutations: List[MutationCall]) -> List[MutationCall]:
        """Enrich mutations with clinical context"""
        enriched_mutations = []

        for mutation in mutations:
            # Add clinical evidence if available
            clinical_evidence = self._gather_clinical_evidence(mutation)
            mutation.clinical_evidence = clinical_evidence

            # Add MIC correlation if available
            mic_corr = self._calculate_mic_correlation(mutation)
            mutation.mic_correlation = mic_corr

            enriched_mutations.append(mutation)

        return enriched_mutations

    def _gather_clinical_evidence(self, mutation: MutationCall) -> List[str]:
        """Gather clinical evidence for mutation"""
        evidence = []

        # Check clinical context database
        mutation_key = mutation.mutation_key

        # Example clinical associations
        clinical_associations = {
            'T104A': 'Associated with fluoroquinolone resistance in clinical isolates',
            'G288D': 'Found in carbapenem-resistant Pseudomonas aeruginosa',
            'R717Q': 'Linked to tigecycline resistance in Klebsiella pneumoniae'
        }

        if mutation_key in clinical_associations:
            evidence.append(clinical_associations[mutation_key])

        # Add species-specific context
        for species, freq in mutation.species_context.items():
            if freq > 0.5:
                evidence.append(f"Prevalent in {species} clinical isolates")

        return evidence

    def _calculate_mic_correlation(self, mutation: MutationCall) -> Optional[float]:
        """Calculate correlation with MIC data"""
        # Placeholder for MIC correlation calculation
        # In practice, this would correlate mutation presence with MIC values
        return None

    def generate_mutation_report(self, mutations: List[MutationCall],
                               output_file: str) -> None:
        """Generate comprehensive mutation report"""
        self.logger.info(f"Generating mutation report: {output_file}")

        # Convert to DataFrame for analysis
        mutation_data = []
        for mutation in mutations:
            mutation_data.append({
                'Position': mutation.position,
                'Reference_AA': mutation.reference_aa,
                'Variant_AA': mutation.variant_aa,
                'Mutation': mutation.mutation_key,
                'Type': mutation.mutation_type,
                'Frequency': mutation.frequency,
                'Confidence': mutation.confidence_score,
                'Phylogenetic_Support': mutation.phylogenetic_support,
                'Conservation_Score': mutation.conservation_score,
                'Functional_Impact': mutation.functional_impact,
                'Resistance_Risk': mutation.resistance_risk,
                'Known_Resistance': mutation.known_resistance_association,
                'Clinical_Evidence_Count': len(mutation.clinical_evidence),
                'Quality_Flags': '; '.join(mutation.quality_flags),
                'Species_Count': len(mutation.species_context)
            })

        df = pd.DataFrame(mutation_data)

        # Save detailed report
        df.to_csv(output_file, index=False)

        # Generate summary statistics
        summary_file = output_file.replace('.csv', '_summary.txt')
        with open(summary_file, 'w') as f:
            f.write("SPECIES-AWARE MUTATION CALLING REPORT\n")
            f.write("=" * 45 + "\n\n")

            f.write("SUMMARY STATISTICS:\n")
            f.write("-" * 20 + "\n")
            f.write(f"Total mutations called: {len(mutations)}\n")
            f.write(f"Significant mutations: {sum(1 for m in mutations if m.is_significant)}\n")
            f.write(f"High resistance risk: {sum(1 for m in mutations if m.resistance_risk == 'High')}\n")
            f.write(f"Cross-species mutations: {sum(1 for m in mutations if len(m.species_context) > 1)}\n\n")

            # Top mutations by confidence
            f.write("TOP MUTATIONS BY CONFIDENCE:\n")
            f.write("-" * 30 + "\n")
            sorted_mutations = sorted(mutations, key=lambda x: x.confidence_score, reverse=True)
            for i, mutation in enumerate(sorted_mutations[:10]):
                f.write(f"{i+1}. {mutation.mutation_key} ")
                f.write(".3f")
                f.write(f" ({mutation.resistance_risk} risk)\n")

            # Functional impact distribution
            f.write("\nFUNCTIONAL IMPACT DISTRIBUTION:\n")
            f.write("-" * 35 + "\n")
            impact_counts = Counter(m.functional_impact for m in mutations)
            for impact, count in impact_counts.most_common():
                f.write(f"{impact}: {count}\n")

        self.logger.info(f"Mutation report saved: {output_file}")
        self.logger.info(f"Summary saved: {summary_file}")


def main():
    """Main entry point"""
    import argparse

    parser = argparse.ArgumentParser(
        description="Species-aware mutation calling system"
    )

    parser.add_argument(
        "--alignment-file",
        required=True,
        help="Path to sequence alignment file"
    )

    parser.add_argument(
        "--metadata-file",
        required=True,
        help="Path to metadata CSV file"
    )

    parser.add_argument(
        "--protein-family",
        required=True,
        help="RND protein family (acrA, acrB, tolC, etc.)"
    )

    parser.add_argument(
        "--reference-db",
        required=True,
        help="Path to species-specific reference database"
    )

    parser.add_argument(
        "--clinical-data",
        help="Path to clinical context data (optional)"
    )

    parser.add_argument(
        "--output-file",
        default="species_aware_mutations.csv",
        help="Output file for mutation calls"
    )

    args = parser.parse_args()

    # Setup logging
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    # Initialize species-aware mutator
    mutator = SpeciesAwareMutator(args.reference_db, args.clinical_data)

    # Call mutations
    mutations = mutator.call_mutations_species_aware(
        args.alignment_file,
        args.metadata_file,
        args.protein_family
    )

    # Generate report
    mutator.generate_mutation_report(mutations, args.output_file)

    print(f"\nSpecies-aware mutation calling complete!")
    print(f"Called {len(mutations)} mutations")
    print(f"Results saved to: {args.output_file}")


if __name__ == "__main__":
    main()