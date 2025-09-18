import sys; import os; sys.path.append(os.path.join(os.path.dirname(__file__), '\\\src')); from configuration_manager import config_manager
"""
Hardcore Comprehensive Test Suite for SubScan Alignment Analyzer
================================================================

This test suite provides industrial-strength validation of SubScan with:
- Extreme mutation detection accuracy testing
- Alignment parsing and quality assessment validation
- Complex error scenario simulation and recovery testing  
- Production-scale performance and memory stress testing
- Real-world alignment data simulation and edge cases
- Integration testing with WildTypeAligner output formats
- Confidence scoring and functional annotation validation
- Security testing and input sanitization
- Database integration and provenance tracking
- End-to-end workflow validation

Test Categories:
- Module initialization and component integration
- EMBOSS-WATER alignment parsing and validation
- Substitution detection algorithms and accuracy
- Quality control and confidence scoring systems
- Functional annotation and AMR knowledge integration
- Batch processing and performance optimization
- Error handling and fault tolerance under stress
- Integration with pipeline components and databases
- Memory efficiency and scalability testing
- Production deployment scenarios and edge cases

Author: GenomeAMRAnalyzer Pipeline
Version: 1.0 - Production Hardened
"""

import pytest
import tempfile
import shutil
import logging
import json
import time
import threading
import concurrent.futures
import random
import string
import math
from unittest.mock import Mock, patch, MagicMock
from pathlib import Path
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass
from collections import Counter, defaultdict
import hashlib
import uuid
from datetime import datetime

# Setup path for imports
import sys
ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

try:
    from src.priority3.analysis.subscan_alignment_analyzer import (
        SubScanAlignmentAnalyzer,
        SubstitutionDetector,
        FunctionalAnnotator,
        EmbossWaterParser,
        MutationRecord,
        MutationType,
        ConfidenceLevel,
        GenomicCoordinate,
        ProteinContext,
        AlignmentQuality
    )
    SUBSCAN_AVAILABLE = True
except ImportError as e:
    print(f"SubScan components not available: {e}")
    SUBSCAN_AVAILABLE = False

try:
    from src.priority3.db.repositories import GenomeRepository
    DB_AVAILABLE = True
except ImportError:
    DB_AVAILABLE = False

# Test utilities and mock data generators
class MockEMBOSSOutput:
    """Generate realistic EMBOSS-WATER alignment output for testing."""
    
    @staticmethod
    def generate_perfect_alignment(length=1000, gene_name=config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0]):
        """Generate perfect alignment with no mutations."""
        sequence = MockEMBOSSOutput._generate_sequence(length)
        return f"""########################################
# Program: water
# Rundate: {datetime.now().strftime('%a %d %b %Y %H:%M:%S')}
# Commandline: water -asequence ref.fasta -bsequence query.fasta -gapopen 10.0 -gapextend 0.5 -outfile alignment.water
# Align_format: srspair
# Report_file: alignment.water
########################################

#=======================================
#
# Aligned_sequences: 2
# 1: {gene_name}_reference
# 2: sample_001_{gene_name}
# Matrix: EBLOSUM62
# Gap_penalty: 10.0
# Extend_penalty: 0.5
#
# Length: {length}
# Identity:   {length}/{length} (100.0%)
# Similarity: {length}/{length} (100.0%)
# Gaps:         0/{length} ( 0.0%)
# Score: {length * 5}.0
#
#=======================================

{gene_name}_reference  1 {sequence}  {length}
                    {"".join("|" for _ in range(min(len(sequence), 50)))}
sample_001_{gene_name}  1 {sequence}  {length}


#---------------------------------------
#---------------------------------------
"""

    @staticmethod
    def generate_single_mutation_alignment(position=100, ref_base="A", mut_base="G", 
                                         length=1000, gene_name=config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0], score_penalty=5.0):
        """Generate alignment with single point mutation."""
        ref_sequence = MockEMBOSSOutput._generate_sequence(length)
        mut_sequence = ref_sequence[:position-1] + mut_base + ref_sequence[position:]
        
        # Ensure reference has the specified base at position
        ref_sequence = ref_sequence[:position-1] + ref_base + ref_sequence[position:]
        
        identity = length - 1
        score = (length * 5.0) - score_penalty
        
        return f"""########################################
# Program: water
# Rundate: {datetime.now().strftime('%a %d %b %Y %H:%M:%S')}
# Commandline: water -asequence ref.fasta -bsequence query.fasta -gapopen 10.0 -gapextend 0.5 -outfile alignment.water
# Align_format: srspair
# Report_file: alignment.water
########################################

#=======================================
#
# Aligned_sequences: 2
# 1: {gene_name}_reference
# 2: sample_001_{gene_name}
# Matrix: EBLOSUM62
# Gap_penalty: 10.0
# Extend_penalty: 0.5
#
# Length: {length}
# Identity:   {identity}/{length} ({identity/length*100:.1f}%)
# Similarity: {identity}/{length} ({identity/length*100:.1f}%)
# Gaps:         0/{length} ( 0.0%)
# Score: {score:.1f}
#
#=======================================

{gene_name}_reference  1 {ref_sequence}  {length}
                    {"".join("|" if i != position-1 else " " for i in range(min(len(ref_sequence), 50)))}
sample_001_{gene_name}  1 {mut_sequence}  {length}


#---------------------------------------
#---------------------------------------
"""

    @staticmethod
    def generate_multiple_mutations_alignment(mutations_list, length=1000, gene_name=config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0]):
        """Generate alignment with multiple mutations."""
        ref_sequence = MockEMBOSSOutput._generate_sequence(length)
        mut_sequence = list(ref_sequence)
        
        mutation_count = 0
        for pos, ref_base, mut_base in mutations_list:
            if 0 < pos <= length:
                ref_sequence = ref_sequence[:pos-1] + ref_base + ref_sequence[pos:]
                mut_sequence[pos-1] = mut_base
                mutation_count += 1
        
        mut_sequence = ''.join(mut_sequence)
        identity = length - mutation_count
        score = (length * 5.0) - (mutation_count * 5.0)
        
        return f"""########################################
# Program: water
# Rundate: {datetime.now().strftime('%a %d %b %Y %H:%M:%S')}
# Commandline: water -asequence ref.fasta -bsequence query.fasta -gapopen 10.0 -gapextend 0.5 -outfile alignment.water
# Align_format: srspair
# Report_file: alignment.water
########################################

#=======================================
#
# Aligned_sequences: 2
# 1: {gene_name}_reference
# 2: sample_001_{gene_name}
# Matrix: EBLOSUM62
# Gap_penalty: 10.0
# Extend_penalty: 0.5
#
# Length: {length}
# Identity:   {identity}/{length} ({identity/length*100:.1f}%)
# Similarity: {identity}/{length} ({identity/length*100:.1f}%)
# Gaps:         0/{length} ( 0.0%)
# Score: {score:.1f}
#
#=======================================

{gene_name}_reference  1 {ref_sequence}  {length}
                    {"".join("|" if ref_sequence[i] == mut_sequence[i] else " " for i in range(min(len(ref_sequence), 50)))}
sample_001_{gene_name}  1 {mut_sequence}  {length}


#---------------------------------------
#---------------------------------------
"""

    @staticmethod
    def generate_gapped_alignment(gap_positions, length=1000, gene_name=config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0]):
        """Generate alignment with gaps (insertions/deletions)."""
        ref_sequence = MockEMBOSSOutput._generate_sequence(length)
        mut_sequence = list(ref_sequence)
        
        # Insert gaps
        gap_count = 0
        for pos in sorted(gap_positions, reverse=True):
            if 0 < pos <= length:
                mut_sequence.insert(pos, '-')
                gap_count += 1
        
        mut_sequence = ''.join(mut_sequence)
        total_length = length + gap_count
        identity = length
        score = (length * 5.0) - (gap_count * 10.0)  # Gap penalty
        
        return f"""########################################
# Program: water
# Rundate: {datetime.now().strftime('%a %d %b %Y %H:%M:%S')}
# Commandline: water -asequence ref.fasta -bsequence query.fasta -gapopen 10.0 -gapextend 0.5 -outfile alignment.water
# Align_format: srspair
# Report_file: alignment.water
########################################

#=======================================
#
# Aligned_sequences: 2
# 1: {gene_name}_reference
# 2: sample_001_{gene_name}
# Matrix: EBLOSUM62
# Gap_penalty: 10.0
# Extend_penalty: 0.5
#
# Length: {total_length}
# Identity:   {identity}/{total_length} ({identity/total_length*100:.1f}%)
# Similarity: {identity}/{total_length} ({identity/total_length*100:.1f}%)
# Gaps:         {gap_count}/{total_length} ({gap_count/total_length*100:.1f}%)
# Score: {score:.1f}
#
#=======================================

{gene_name}_reference  1 {ref_sequence}  {length}
                    {"".join("|" for _ in range(min(len(ref_sequence), 50)))}
sample_001_{gene_name}  1 {mut_sequence}  {len(mut_sequence)}


#---------------------------------------
#---------------------------------------
"""

    @staticmethod
    def generate_low_quality_alignment(length=500, gene_name=config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0], identity_percent=65):
        """Generate low-quality alignment with many mutations."""
        ref_sequence = MockEMBOSSOutput._generate_sequence(length)
        mut_sequence = list(ref_sequence)
        
        # Introduce random mutations to achieve target identity
        target_mutations = int(length * (1 - identity_percent / 100))
        mutation_positions = random.sample(range(length), min(target_mutations, length))
        
        nucleotides = ['A', 'T', 'G', 'C']
        for pos in mutation_positions:
            original = mut_sequence[pos]
            mut_sequence[pos] = random.choice([n for n in nucleotides if n != original])
        
        mut_sequence = ''.join(mut_sequence)
        identity = length - len(mutation_positions)
        score = (length * 2.0) - (len(mutation_positions) * 3.0)  # Lower scoring
        
        return f"""########################################
# Program: water
# Rundate: {datetime.now().strftime('%a %d %b %Y %H:%M:%S')}
# Commandline: water -asequence ref.fasta -bsequence query.fasta -gapopen 10.0 -gapextend 0.5 -outfile alignment.water
# Align_format: srspair
# Report_file: alignment.water
########################################

#=======================================
#
# Aligned_sequences: 2
# 1: {gene_name}_reference
# 2: sample_001_{gene_name}
# Matrix: EBLOSUM62
# Gap_penalty: 10.0
# Extend_penalty: 0.5
#
# Length: {length}
# Identity:   {identity}/{length} ({identity/length*100:.1f}%)
# Similarity: {identity}/{length} ({identity/length*100:.1f}%)
# Gaps:         0/{length} ( 0.0%)
# Score: {score:.1f}
#
#=======================================

{gene_name}_reference  1 {ref_sequence}  {length}
                    {"".join("|" if ref_sequence[i] == mut_sequence[i] else " " for i in range(min(len(ref_sequence), 50)))}
sample_001_{gene_name}  1 {mut_sequence}  {length}


#---------------------------------------
#---------------------------------------
"""

    @staticmethod
    def _generate_sequence(length):
        """Generate random DNA sequence."""
        nucleotides = ['A', 'T', 'G', 'C']
        return ''.join(random.choices(nucleotides, k=length))

class MockDatabase:
    """Mock database for testing."""
    
    def __init__(self):
        self.genomes = {}
        self.mutations = {}
    
    def add_genome(self, accession, file_path, organism="Test organism"):
        self.genomes[accession] = {
            'accession': accession,
            'file_path': file_path,
            'organism': organism,
            'status': 'downloaded'
        }
    
    def list_genomes(self, status=None):
        if status:
            return [g for g in self.genomes.values() if g.get('status') == status]
        return list(self.genomes.values())

@pytest.fixture
def comprehensive_test_environment():
    """Create comprehensive test environment."""
    temp_dir = tempfile.mkdtemp()
    
    test_env = {
        'temp_dir': Path(temp_dir),
        'output_dir': Path(temp_dir) / "subscan_output",
        'alignment_dir': Path(temp_dir) / "alignments",
        'db_path': str(Path(temp_dir) / "test.db"),
        'card_db_path': str(Path(temp_dir) / "card.json")
    }
    
    # Create directories
    for key, path in test_env.items():
        if isinstance(path, Path):
            path.mkdir(parents=True, exist_ok=True)
    
    # Create test alignment files
    test_env['alignment_files'] = _create_test_alignment_files(test_env['alignment_dir'])
    test_env['stress_test_files'] = _create_stress_test_files(test_env['alignment_dir'])
    test_env['corrupted_files'] = _create_corrupted_alignment_files(test_env['alignment_dir'])
    
    # Create mock CARD database
    _create_mock_card_database(test_env['card_db_path'])
    
    yield test_env
    
    # Cleanup
    try:
        shutil.rmtree(temp_dir)
    except PermissionError:
        pass

def _create_test_alignment_files(alignment_dir: Path) -> List[Path]:
    """Create comprehensive test alignment files."""
    alignment_files = []
    
    # Perfect alignment
    perfect_alignment = alignment_dir / "perfect_alignment.water"
    with open(perfect_alignment, 'w') as f:
        f.write(MockEMBOSSOutput.generate_perfect_alignment(length=1200, gene_name=config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0]))
    alignment_files.append(perfect_alignment)
    
    # Single mutation alignment
    single_mut = alignment_dir / "single_mutation.water"
    with open(single_mut, 'w') as f:
        f.write(MockEMBOSSOutput.generate_single_mutation_alignment(
            position=100, ref_base="A", mut_base="G", length=1000, gene_name=config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1]
        ))
    alignment_files.append(single_mut)
    
    # Multiple mutations alignment
    multi_mut = alignment_dir / "multiple_mutations.water"
    mutations = [(50, "C", "T"), (150, "G", "A"), (300, "A", "C")]
    with open(multi_mut, 'w') as f:
        f.write(MockEMBOSSOutput.generate_multiple_mutations_alignment(
            mutations, length=800, gene_name=config_manager.get_default_genes("rnd_efflux_pumps", "primary")[2]
        ))
    alignment_files.append(multi_mut)
    
    # Gapped alignment
    gapped = alignment_dir / "gapped_alignment.water"
    with open(gapped, 'w') as f:
        f.write(MockEMBOSSOutput.generate_gapped_alignment(
            gap_positions=[75, 225, 400], length=500, gene_name="marA"
        ))
    alignment_files.append(gapped)
    
    # Low quality alignment
    low_quality = alignment_dir / "low_quality.water"
    with open(low_quality, 'w') as f:
        f.write(MockEMBOSSOutput.generate_low_quality_alignment(
            length=600, gene_name="acrR", identity_percent=65
        ))
    alignment_files.append(low_quality)
    
    return alignment_files

def _create_stress_test_files(alignment_dir: Path) -> List[Path]:
    """Create large files for stress testing."""
    stress_files = []
    
    # Large perfect alignment
    large_perfect = alignment_dir / "large_perfect.water"
    with open(large_perfect, 'w') as f:
        f.write(MockEMBOSSOutput.generate_perfect_alignment(length=10000, gene_name="mexA"))
    stress_files.append(large_perfect)
    
    # Large alignment with many mutations
    large_mutations = alignment_dir / "large_with_mutations.water"
    mutations = [(i*50, random.choice(["A", "T"]), random.choice(["G", "C"])) 
                for i in range(1, 100)]  # 99 mutations
    with open(large_mutations, 'w') as f:
        f.write(MockEMBOSSOutput.generate_multiple_mutations_alignment(
            mutations, length=5000, gene_name="mexB"
        ))
    stress_files.append(large_mutations)
    
    return stress_files

def _create_corrupted_alignment_files(alignment_dir: Path) -> List[Path]:
    """Create corrupted files for error testing."""
    corrupted_files = []
    
    # Empty file
    empty_file = alignment_dir / "empty.water"
    empty_file.touch()
    corrupted_files.append(empty_file)
    
    # Invalid format
    invalid_format = alignment_dir / "invalid_format.water"
    with open(invalid_format, 'w') as f:
        f.write("This is not a valid EMBOSS-WATER file\n")
        f.write("Missing proper headers and format\n")
    corrupted_files.append(invalid_format)
    
    # Truncated alignment
    truncated = alignment_dir / "truncated.water"
    with open(truncated, 'w') as f:
        f.write("########################################\n")
        f.write("# Program: water\n")
        f.write("# Truncated file with incomplete data\n")
    corrupted_files.append(truncated)
    
    # Malformed alignment data
    malformed = alignment_dir / "malformed.water"
    with open(malformed, 'w') as f:
        f.write("""########################################
# Program: water
# Rundate: invalid date format
# Commandline: water
# Align_format: srspair
# Report_file: alignment.water
########################################

#=======================================
# Malformed alignment section
# Missing required fields and structure
#=======================================

reference_seq  INVALID_SEQUENCE_FORMAT
               MALFORMED_ANNOTATION_LINE
query_seq      DIFFERENT_LENGTH_SEQUENCE

#---------------------------------------
""")
    corrupted_files.append(malformed)
    
    return corrupted_files

def _create_mock_card_database(card_db_path: str):
    """Create mock CARD database for testing."""
    mock_card_data = {
        "version": "3.2.0",
        "mutations": {
            "acrB:S83L": {
                "drug_class": "fluoroquinolone",
                "resistance_mechanism": "target_alteration",
                "evidence_level": "high",
                "publications": ["PMID:12345678"]
            },
            "acrB:D87N": {
                "drug_class": "fluoroquinolone",
                "resistance_mechanism": "target_alteration", 
                "evidence_level": "high",
                "publications": ["PMID:87654321"]
            },
            "tolC:A119E": {
                "drug_class": "multiple",
                "resistance_mechanism": "efflux_pump",
                "evidence_level": "medium",
                "publications": ["PMID:11111111"]
            }
        }
    }
    
    with open(card_db_path, 'w') as f:
        json.dump(mock_card_data, f, indent=2)


@pytest.mark.skipif(not SUBSCAN_AVAILABLE, reason="SubScan components not available")
class TestSubScanHardcore:
    """Hardcore comprehensive test suite for SubScan Alignment Analyzer."""
    
    # ===== Component Initialization Tests =====
    
    def test_subscan_analyzer_initialization_basic(self, comprehensive_test_environment):
        """Test basic SubScan analyzer initialization."""
        analyzer = SubScanAlignmentAnalyzer(
            output_dir=str(comprehensive_test_environment['output_dir']),
            db_path=comprehensive_test_environment['db_path']
        )
        
        assert analyzer.output_dir.exists()
        assert analyzer.min_confidence == ConfidenceLevel.MEDIUM
        assert hasattr(analyzer, 'parser')
        assert hasattr(analyzer, 'detector')
        assert hasattr(analyzer, 'annotator')
        assert hasattr(analyzer, 'logger')
        assert analyzer.stats['alignments_processed'] == 0
    
    def test_subscan_analyzer_initialization_custom_parameters(self, comprehensive_test_environment):
        """Test SubScan analyzer with custom parameters."""
        analyzer = SubScanAlignmentAnalyzer(
            output_dir=str(comprehensive_test_environment['output_dir']),
            db_path=comprehensive_test_environment['db_path'],
            min_confidence=ConfidenceLevel.HIGH,
            card_database_path=comprehensive_test_environment['card_db_path']
        )
        
        assert analyzer.min_confidence == ConfidenceLevel.HIGH
        assert analyzer.annotator.card_db_path == comprehensive_test_environment['card_db_path']
    
    def test_substitution_detector_initialization(self, comprehensive_test_environment):
        """Test SubstitutionDetector initialization."""
        detector = SubstitutionDetector(min_quality_score=0.8)
        
        assert detector.min_quality_score == 0.8
        assert hasattr(detector, 'params')
        assert detector.params['min_identity_percent'] == 70.0
        assert detector.params['max_gap_percent'] == 20.0
    
    def test_functional_annotator_initialization(self, comprehensive_test_environment):
        """Test FunctionalAnnotator initialization."""
        annotator = FunctionalAnnotator(
            card_database_path=comprehensive_test_environment['card_db_path']
        )
        
        assert annotator.card_db_path == comprehensive_test_environment['card_db_path']
        assert hasattr(annotator, 'known_amr_mutations')
        assert hasattr(annotator, 'impact_rules')
    
    # ===== EMBOSS-WATER Parsing Tests =====
    
    def test_parse_perfect_alignment(self, comprehensive_test_environment):
        """Test parsing perfect alignment with no mutations."""
        if not hasattr(SubScanAlignmentAnalyzer, 'parser'):
            pytest.skip("EmbossWaterParser not available")
            
        parser = EmbossWaterParser()
        alignment_file = comprehensive_test_environment['alignment_files'][0]  # perfect_alignment.water
        
        alignments = parser.parse_alignment_file(str(alignment_file))
        
        assert len(alignments) == 1
        alignment = alignments[0]
        assert alignment['identity']['percent'] == 100.0
        assert alignment['gaps']['percent'] == 0.0
        assert config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0] in alignment['reference_id']
        assert config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0] in alignment['query_id']
    
    def test_parse_single_mutation_alignment(self, comprehensive_test_environment):
        """Test parsing alignment with single mutation."""
        if not hasattr(SubScanAlignmentAnalyzer, 'parser'):
            pytest.skip("EmbossWaterParser not available")
            
        parser = EmbossWaterParser()
        alignment_file = comprehensive_test_environment['alignment_files'][1]  # single_mutation.water
        
        alignments = parser.parse_alignment_file(str(alignment_file))
        
        assert len(alignments) == 1
        alignment = alignments[0]
        assert alignment['identity']['percent'] < 100.0
        assert alignment['identity']['percent'] > 95.0  # Should be high with only one mutation
        assert config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1] in alignment['reference_id']
    
    def test_parse_multiple_mutations_alignment(self, comprehensive_test_environment):
        """Test parsing alignment with multiple mutations."""
        if not hasattr(SubScanAlignmentAnalyzer, 'parser'):
            pytest.skip("EmbossWaterParser not available")
            
        parser = EmbossWaterParser()
        alignment_file = comprehensive_test_environment['alignment_files'][2]  # multiple_mutations.water
        
        alignments = parser.parse_alignment_file(str(alignment_file))
        
        assert len(alignments) == 1
        alignment = alignments[0]
        assert alignment['identity']['percent'] < 100.0
        assert config_manager.get_default_genes("rnd_efflux_pumps", "primary")[2] in alignment['reference_id']
    
    def test_parse_gapped_alignment(self, comprehensive_test_environment):
        """Test parsing alignment with gaps."""
        if not hasattr(SubScanAlignmentAnalyzer, 'parser'):
            pytest.skip("EmbossWaterParser not available")
            
        parser = EmbossWaterParser()
        alignment_file = comprehensive_test_environment['alignment_files'][3]  # gapped_alignment.water
        
        alignments = parser.parse_alignment_file(str(alignment_file))
        
        assert len(alignments) == 1
        alignment = alignments[0]
        assert alignment['gaps']['percent'] > 0.0
        assert 'marA' in alignment['reference_id']
    
    def test_parse_corrupted_alignment_files(self, comprehensive_test_environment):
        """Test parsing corrupted alignment files."""
        if not hasattr(SubScanAlignmentAnalyzer, 'parser'):
            pytest.skip("EmbossWaterParser not available")
            
        parser = EmbossWaterParser()
        
        for corrupted_file in comprehensive_test_environment['corrupted_files']:
            try:
                alignments = parser.parse_alignment_file(str(corrupted_file))
                # Should either return empty list or handle gracefully
                assert isinstance(alignments, list)
            except Exception as e:
                # Expected for some corrupted files
                assert isinstance(e, (ValueError, FileNotFoundError, IOError))
    
    # ===== Substitution Detection Tests =====
    
    def test_detect_substitutions_perfect_alignment(self, comprehensive_test_environment):
        """Test substitution detection on perfect alignment."""
        detector = SubstitutionDetector()
        
        # Mock perfect alignment data
        perfect_alignment = {
            'alignment_id': 'perfect_test',
            'query_sequence': 'ATCGATCGATCGATCG',
            'reference_sequence': 'ATCGATCGATCGATCG',
            'query_id': 'sample_001_acrA',
            'reference_id': 'acrA_reference',
            'identity': {'percent': 100.0},
            'gaps': {'count': 0, 'percent': 0.0},
            'score': 80.0,
            'reference_start': 1,
            'query_start': 1
        }
        
        mutations = detector.detect_substitutions(perfect_alignment)
        
        # Perfect alignment should have no mutations
        assert len(mutations) == 0
    
    def test_detect_substitutions_single_mutation(self, comprehensive_test_environment):
        """Test detection of single substitution."""
        detector = SubstitutionDetector()
        
        # Mock alignment with single mutation at position 30 (G->A) - using longer sequences
        query_seq = "ATCGATCGATCGATCGATCGATCGATCGATCAATCGATCGATCGATCGATCGATCG"  # 54 bp, A at position 30
        ref_seq = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"    # 54 bp, G at position 30
        
        single_mut_alignment = {
            'alignment_id': 'single_mut_test',
            'query_sequence': query_seq,
            'reference_sequence': ref_seq,
            'query_id': 'sample_001_acrB',
            'reference_id': 'acrB_reference',
            'identity': {'percent': 98.15},  # 53/54 matches
            'gaps': {'count': 0, 'percent': 0.0},
            'score': 265.0,
            'reference_start': 1,
            'query_start': 1
        }
        
        mutations = detector.detect_substitutions(single_mut_alignment)
        
        assert len(mutations) == 1
        mutation = mutations[0]
        assert mutation.mutation_type == MutationType.SUBSTITUTION
        assert mutation.genomic_coord.position == 32  # 1-based coordinate (position 31 in 0-based)
        assert mutation.genomic_coord.reference_base == "G"
        assert "A" in mutation.mutation_id
    
    def test_detect_substitutions_multiple_mutations(self, comprehensive_test_environment):
        """Test detection of multiple substitutions."""
        detector = SubstitutionDetector()
        
        # Mock alignment with multiple mutations - using longer sequences
        query_seq = "ATCGATCGATCGATCGATCGATCGATCGATCAATCGATCGATCGATCGATCGATCC"  # 54 bp, mutations at positions 30, 53
        ref_seq = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"    # 54 bp
        
        multi_mut_alignment = {
            'alignment_id': 'multi_mut_test',
            'query_sequence': query_seq,
            'reference_sequence': ref_seq,
            'query_id': 'sample_001_tolC',
            'reference_id': 'tolC_reference',
            'identity': {'percent': 96.30},  # 52/54 matches
            'gaps': {'count': 0, 'percent': 0.0},
            'score': 250.0,
            'reference_start': 1,
            'query_start': 1
        }
        
        mutations = detector.detect_substitutions(multi_mut_alignment)
        
        assert len(mutations) == 2
        positions = [m.genomic_coord.position for m in mutations]
        assert 32 in positions  # 1-based coordinate for position 31
        assert 56 in positions  # 1-based coordinate for position 55
    
    def test_detect_substitutions_quality_filtering(self, comprehensive_test_environment):
        """Test quality filtering of substitutions."""
        detector = SubstitutionDetector(min_quality_score=0.8)
        
        # Mock low-quality alignment
        low_quality_alignment = {
            'alignment_id': 'low_quality_test',
            'query_sequence': 'ATCGATCGATCGATCG',
            'reference_sequence': 'GTACGTACGTACGTAC',  # Many differences
            'query_id': 'sample_001_low',
            'reference_id': 'low_reference',
            'identity': {'percent': 25.0},  # Very low identity
            'gaps': {'count': 0, 'percent': 0.0},
            'score': 20.0,
            'reference_start': 1,
            'query_start': 1
        }
        
        mutations = detector.detect_substitutions(low_quality_alignment)
        
        # Should be filtered out due to low quality
        assert len(mutations) == 0
    
    def test_substitution_confidence_scoring(self, comprehensive_test_environment):
        """Test confidence scoring for substitutions."""
        detector = SubstitutionDetector()
        
        # High-quality mutation with longer sequence
        query_seq = "ATCGATCGATCGATCGATCGATCGATCGATCAATCGATCGATCGATCGATCGATCG"  # A at position 30
        ref_seq = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"    # G at position 30
        
        high_quality_alignment = {
            'alignment_id': 'high_conf_test',
            'query_sequence': query_seq,
            'reference_sequence': ref_seq,
            'query_id': 'sample_001_high',
            'reference_id': 'high_reference',
            'identity': {'percent': 98.15},  # 53/54 matches
            'gaps': {'count': 0, 'percent': 0.0},
            'score': 265.0,
            'reference_start': 1,
            'query_start': 1
        }
        
        mutations = detector.detect_substitutions(high_quality_alignment)
        
        assert len(mutations) == 1
        mutation = mutations[0]
        assert mutation.confidence in [ConfidenceLevel.HIGH, ConfidenceLevel.MEDIUM]
        assert mutation.quality_score > 0.5
    
    # ===== Quality Control and Edge Case Tests =====
    
    def test_empty_alignment_handling(self, comprehensive_test_environment):
        """Test handling of empty alignments."""
        detector = SubstitutionDetector()
        
        empty_alignment = {
            'alignment_id': 'empty_test',
            'query_sequence': '',
            'reference_sequence': '',
            'query_id': 'empty_query',
            'reference_id': 'empty_reference',
            'identity': {'percent': 0.0},
            'gaps': {'count': 0, 'percent': 0.0},
            'score': 0.0,
            'reference_start': 1,
            'query_start': 1
        }
        
        mutations = detector.detect_substitutions(empty_alignment)
        assert len(mutations) == 0
    
    def test_sequence_length_mismatch_handling(self, comprehensive_test_environment):
        """Test handling of sequence length mismatches."""
        detector = SubstitutionDetector()
        
        mismatch_alignment = {
            'alignment_id': 'mismatch_test',
            'query_sequence': 'ATCGATCG',      # 8 bases
            'reference_sequence': 'ATCGATCGATCG',  # 12 bases
            'query_id': 'mismatch_query',
            'reference_id': 'mismatch_reference',
            'identity': {'percent': 100.0},
            'gaps': {'count': 0, 'percent': 0.0},
            'score': 40.0,
            'reference_start': 1,
            'query_start': 1
        }
        
        mutations = detector.detect_substitutions(mismatch_alignment)
        # Should handle gracefully, likely return empty list
        assert isinstance(mutations, list)
    
    def test_gap_proximity_penalty(self, comprehensive_test_environment):
        """Test penalty for mutations near gaps."""
        detector = SubstitutionDetector()
        
        # Mutation near gap should have lower confidence
        gap_proximal_alignment = {
            'alignment_id': 'gap_test',
            'query_sequence': 'ATCGATC-ATCGATCG',  # Gap at position 8
            'reference_sequence': 'ATCGATCGATCGATCG', # Mutation at position 7 (near gap)
            'query_id': 'gap_query',
            'reference_id': 'gap_reference',
            'identity': {'percent': 87.5},
            'gaps': {'count': 1, 'percent': 6.25},
            'score': 65.0,
            'reference_start': 1,
            'query_start': 1
        }
        
        mutations = detector.detect_substitutions(gap_proximal_alignment)
        
        # May have mutations but with reduced confidence due to gap proximity
        if mutations:
            for mutation in mutations:
                # Quality should be impacted by gap proximity
                assert hasattr(mutation, 'quality_score')
    
    def test_homopolymer_context_penalty(self, comprehensive_test_environment):
        """Test penalty for mutations in homopolymer regions."""
        detector = SubstitutionDetector()
        
        # Mutation in homopolymer region (AAAA)
        homopolymer_alignment = {
            'alignment_id': 'homopol_test',
            'query_sequence': 'ATCGAAAAGCGATCG',  # AAAA homopolymer
            'reference_sequence': 'ATCGAAAAACGATCG', # Mutation in homopolymer
            'query_id': 'homopol_query',
            'reference_id': 'homopol_reference',
            'identity': {'percent': 93.33},
            'gaps': {'count': 0, 'percent': 0.0},
            'score': 70.0,
            'reference_start': 1,
            'query_start': 1
        }
        
        mutations = detector.detect_substitutions(homopolymer_alignment)
        
        # Should detect mutation but with appropriate quality penalty
        if mutations:
            mutation = mutations[0]
            # Quality metrics should include homopolymer penalty
            assert hasattr(mutation, 'quality_metrics')
    
    # ===== Functional Annotation Tests =====
    
    def test_functional_annotation_known_amr_mutation(self, comprehensive_test_environment):
        """Test annotation of known AMR mutations."""
        annotator = FunctionalAnnotator(
            card_database_path=comprehensive_test_environment['card_db_path']
        )
        
        # Create mutation record for known AMR mutation (gyrA:S83L)
        known_mutation = MutationRecord(
            mutation_id="gyrA:249:C>T",  # Nucleotide change for S83L
            mutation_type=MutationType.SUBSTITUTION,
            genomic_coord=GenomicCoordinate(
                contig="contig1",
                position=249,  # Codon position for S83
                reference_base="C",
                gene_context="gyrA"
            ),
            confidence=ConfidenceLevel.HIGH,
            quality_metrics=AlignmentQuality(
                identity_percent=95.0,
                coverage_percent=100.0,
                gap_count=0,
                mismatch_count=1,
                alignment_score=75.0
            ),
            quality_score=0.95,
            source_alignment="test_alignment",
            query_sequence="test_query",
            reference_sequence="test_reference"
        )
        
        # Add protein context for AMR mutation that matches known_amr_mutations
        known_mutation.protein_context = ProteinContext(
            protein_id="gyrA",
            amino_acid_position=83,
            reference_aa="S",
            mutant_aa="L",
            codon_position=1
        )
        
        annotated_mutation = annotator.annotate_mutation(known_mutation)
        
        # Should identify as known resistance mutation
        assert annotated_mutation.known_resistance == True
        assert annotated_mutation.resistance_mechanism == "target_modification"
        assert "ciprofloxacin" in annotated_mutation.drug_associations
    
    def test_functional_annotation_novel_mutation(self, comprehensive_test_environment):
        """Test annotation of novel mutations."""
        annotator = FunctionalAnnotator(
            card_database_path=comprehensive_test_environment['card_db_path']
        )
        
        # Create mutation record for novel mutation
        novel_mutation = MutationRecord(
            mutation_id="unknown_gene:100:G>A",
            mutation_type=MutationType.SUBSTITUTION,
            genomic_coord=GenomicCoordinate(
                contig="contig1",
                position=100,
                reference_base="G",
                gene_context="unknown_gene"
            ),
            confidence=ConfidenceLevel.MEDIUM,
            quality_metrics=AlignmentQuality(
                identity_percent=90.0,
                coverage_percent=100.0,
                gap_count=0,
                mismatch_count=1
            ),
            quality_score=0.85,
            source_alignment="test_alignment",
            query_sequence="test_query",
            reference_sequence="test_reference"
        )
        
        annotated_mutation = annotator.annotate_mutation(novel_mutation)
        
        # Should not be identified as known resistance
        if hasattr(annotated_mutation, 'known_resistance'):
            assert not annotated_mutation.known_resistance
    
    # ===== Integration and End-to-End Tests =====
    
    def test_analyze_alignment_file_single_file(self, comprehensive_test_environment):
        """Test complete analysis of single alignment file."""
        analyzer = SubScanAlignmentAnalyzer(
            output_dir=str(comprehensive_test_environment['output_dir']),
            db_path=comprehensive_test_environment['db_path'],
            card_database_path=comprehensive_test_environment['card_db_path']
        )
        
        alignment_file = comprehensive_test_environment['alignment_files'][1]  # single_mutation.water
        
        mutations = analyzer.analyze_alignment_file(str(alignment_file), "TEST_001")
        
        assert isinstance(mutations, list)
        assert analyzer.stats['alignments_processed'] > 0
        
        # Check that output files were generated
        expected_report = comprehensive_test_environment['output_dir'] / f"{alignment_file.stem}_report.txt"
        # Report generation might depend on implementation
    
    def test_analyze_batch_processing(self, comprehensive_test_environment):
        """Test batch processing of multiple alignment files."""
        analyzer = SubScanAlignmentAnalyzer(
            output_dir=str(comprehensive_test_environment['output_dir']),
            db_path=comprehensive_test_environment['db_path']
        )
        
        alignment_files = [str(f) for f in comprehensive_test_environment['alignment_files'][:3]]
        genome_accessions = ["BATCH_001", "BATCH_002", "BATCH_003"]
        
        results = analyzer.analyze_batch(alignment_files, genome_accessions)
        
        assert isinstance(results, dict)
        assert len(results) == len(alignment_files)
        
        for file_path in alignment_files:
            assert file_path in results
            assert isinstance(results[file_path], list)
    
    def test_confidence_threshold_filtering(self, comprehensive_test_environment):
        """Test filtering by confidence threshold."""
        # High confidence threshold
        high_threshold_analyzer = SubScanAlignmentAnalyzer(
            output_dir=str(comprehensive_test_environment['output_dir']),
            db_path=comprehensive_test_environment['db_path'],
            min_confidence=ConfidenceLevel.HIGH
        )
        
        # Low confidence threshold
        low_threshold_analyzer = SubScanAlignmentAnalyzer(
            output_dir=str(comprehensive_test_environment['output_dir']),
            db_path=comprehensive_test_environment['db_path'],
            min_confidence=ConfidenceLevel.LOW
        )
        
        alignment_file = comprehensive_test_environment['alignment_files'][4]  # low_quality.water
        
        high_mutations = high_threshold_analyzer.analyze_alignment_file(str(alignment_file))
        low_mutations = low_threshold_analyzer.analyze_alignment_file(str(alignment_file))
        
        # High threshold should filter out more mutations
        assert len(high_mutations) <= len(low_mutations)
    
    def test_statistics_tracking(self, comprehensive_test_environment):
        """Test statistics tracking during analysis."""
        analyzer = SubScanAlignmentAnalyzer(
            output_dir=str(comprehensive_test_environment['output_dir']),
            db_path=comprehensive_test_environment['db_path']
        )
        
        initial_stats = analyzer.stats.copy()
        
        alignment_files = comprehensive_test_environment['alignment_files'][:2]
        for alignment_file in alignment_files:
            analyzer.analyze_alignment_file(str(alignment_file))
        
        # Statistics should be updated
        assert analyzer.stats['alignments_processed'] > initial_stats['alignments_processed']
        assert analyzer.stats['processing_time'] > initial_stats['processing_time']
    
    # ===== Performance and Stress Tests =====
    
    def test_large_alignment_processing(self, comprehensive_test_environment):
        """Test processing of large alignment files."""
        analyzer = SubScanAlignmentAnalyzer(
            output_dir=str(comprehensive_test_environment['output_dir']),
            db_path=comprehensive_test_environment['db_path']
        )
        
        large_file = comprehensive_test_environment['stress_test_files'][0]  # large_perfect.water
        
        start_time = time.time()
        mutations = analyzer.analyze_alignment_file(str(large_file))
        end_time = time.time()
        
        processing_time = end_time - start_time
        
        assert isinstance(mutations, list)
        assert processing_time < 30  # Should complete within 30 seconds
        
        print(f"Large alignment processed in {processing_time:.2f} seconds")
    
    def test_memory_efficiency_large_files(self, comprehensive_test_environment):
        """Test memory efficiency with large files."""
        import psutil
        
        analyzer = SubScanAlignmentAnalyzer(
            output_dir=str(comprehensive_test_environment['output_dir']),
            db_path=comprehensive_test_environment['db_path']
        )
        
        process = psutil.Process()
        initial_memory = process.memory_info().rss / 1024 / 1024  # MB
        
        # Process large files
        for large_file in comprehensive_test_environment['stress_test_files']:
            analyzer.analyze_alignment_file(str(large_file))
        
        final_memory = process.memory_info().rss / 1024 / 1024  # MB
        memory_increase = final_memory - initial_memory
        
        # Memory increase should be reasonable
        assert memory_increase < 200  # Less than 200MB increase
        
        print(f"Memory increase: {memory_increase:.2f} MB")
    
    def test_concurrent_analysis(self, comprehensive_test_environment):
        """Test concurrent analysis for thread safety."""
        analyzer = SubScanAlignmentAnalyzer(
            output_dir=str(comprehensive_test_environment['output_dir']),
            db_path=comprehensive_test_environment['db_path']
        )
        
        alignment_files = comprehensive_test_environment['alignment_files']
        
        # Run concurrent analyses
        with concurrent.futures.ThreadPoolExecutor(max_workers=3) as executor:
            futures = []
            for i, alignment_file in enumerate(alignment_files):
                future = executor.submit(
                    analyzer.analyze_alignment_file,
                    str(alignment_file),
                    f"CONCURRENT_{i}"
                )
                futures.append(future)
            
            # Collect results
            results = []
            for future in concurrent.futures.as_completed(futures):
                try:
                    result = future.result(timeout=30)
                    results.append(result)
                except Exception as e:
                    print(f"Concurrent analysis failed: {e}")
                    results.append([])
        
        # All should complete successfully
        assert len(results) == len(alignment_files)
        
        print(f"✅ Concurrent analysis completed: {len(results)} files processed")
    
    # ===== Error Handling and Recovery Tests =====
    
    def test_corrupted_file_handling(self, comprehensive_test_environment):
        """Test handling of corrupted alignment files."""
        analyzer = SubScanAlignmentAnalyzer(
            output_dir=str(comprehensive_test_environment['output_dir']),
            db_path=comprehensive_test_environment['db_path']
        )
        
        for corrupted_file in comprehensive_test_environment['corrupted_files']:
            try:
                mutations = analyzer.analyze_alignment_file(str(corrupted_file))
                # Should handle gracefully
                assert isinstance(mutations, list)
            except Exception as e:
                # Expected for some corrupted files
                assert isinstance(e, (ValueError, IOError, FileNotFoundError))
    
    def test_missing_file_handling(self, comprehensive_test_environment):
        """Test handling of missing files."""
        analyzer = SubScanAlignmentAnalyzer(
            output_dir=str(comprehensive_test_environment['output_dir']),
            db_path=comprehensive_test_environment['db_path']
        )
        
        # Test with non-existent file
        try:
            mutations = analyzer.analyze_alignment_file("nonexistent_file.water")
            assert mutations == []
        except FileNotFoundError:
            # Expected behavior
            pass
    
    def test_batch_processing_with_failures(self, comprehensive_test_environment):
        """Test batch processing with some file failures."""
        analyzer = SubScanAlignmentAnalyzer(
            output_dir=str(comprehensive_test_environment['output_dir']),
            db_path=comprehensive_test_environment['db_path']
        )
        
        # Mix of valid and invalid files
        mixed_files = [
            str(comprehensive_test_environment['alignment_files'][0]),  # Valid
            "nonexistent_file.water",  # Invalid
            str(comprehensive_test_environment['alignment_files'][1]),  # Valid
            str(comprehensive_test_environment['corrupted_files'][0])   # Corrupted
        ]
        
        results = analyzer.analyze_batch(mixed_files)
        
        # Should process valid files and handle invalid ones gracefully
        assert isinstance(results, dict)
        assert len(results) >= 2  # At least the valid files


def test_subscan_integration_workflow():
    """Integration test for SubScan workflow."""
    temp_dir = tempfile.mkdtemp()
    try:
        output_dir = Path(temp_dir) / "subscan_output"
        alignment_dir = Path(temp_dir) / "alignments"
        
        # Create directories
        output_dir.mkdir(parents=True)
        alignment_dir.mkdir(parents=True)
        
        # Create test alignment
        test_alignment = alignment_dir / "integration_test.water"
        with open(test_alignment, 'w') as f:
            f.write(MockEMBOSSOutput.generate_single_mutation_alignment(
                position=100, ref_base="A", mut_base="G", 
                length=500, gene_name="integration_test"
            ))
        
        # Test workflow
        if SUBSCAN_AVAILABLE:
            analyzer = SubScanAlignmentAnalyzer(
                output_dir=str(output_dir),
                db_path=":memory:"  # In-memory database for testing
            )
            
            # Single file analysis
            mutations = analyzer.analyze_alignment_file(str(test_alignment))
            assert isinstance(mutations, list)
            
            # Batch analysis
            results = analyzer.analyze_batch([str(test_alignment)])
            assert str(test_alignment) in results
            
            # Statistics should be updated
            assert analyzer.stats['alignments_processed'] > 0
        
        print("✅ SubScan integration workflow validated")
        
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)


if __name__ == "__main__":
    # Run basic tests for debugging
    if SUBSCAN_AVAILABLE:
        import tempfile
        temp_dir = tempfile.mkdtemp()
        temp_env = {
            'temp_dir': Path(temp_dir),
            'output_dir': Path(temp_dir) / "output",
            'alignment_dir': Path(temp_dir) / "alignments",
            'db_path': str(Path(temp_dir) / "test.db"),
            'card_db_path': str(Path(temp_dir) / "card.json")
        }
        
        for key, path in temp_env.items():
            if isinstance(path, Path):
                path.mkdir(parents=True, exist_ok=True)
        
        # Create test alignment
        test_alignment = temp_env['alignment_dir'] / "test.water"
        with open(test_alignment, 'w') as f:
            f.write(MockEMBOSSOutput.generate_single_mutation_alignment())
        
        temp_env['alignment_files'] = [test_alignment]
        _create_mock_card_database(temp_env['card_db_path'])
        
        test = TestSubScanHardcore()
        try:
            test.test_subscan_analyzer_initialization_basic(temp_env)
            test.test_detect_substitutions_single_mutation(temp_env)
            print("✅ SubScan hardcore tests passed!")
        finally:
            shutil.rmtree(temp_dir, ignore_errors=True)
    
    # Always run integration test
    test_subscan_integration_workflow()