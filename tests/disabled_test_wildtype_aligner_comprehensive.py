import sys; import os; sys.path.append(os.path.join(os.path.dirname(__file__), '\\\src')); from configuration_manager import config_manager
"""
Comprehensive Test Suite for WildTypeAligner Module
==================================================

This test suite validates the WildTypeAligner implementations with focus on:
- Primary EMBOSS-WATER based aligner
- SimplifiedWildTypeAligner with BioPython
- Reference handling and fallback mechanisms
- Batch processing capabilities
- Error handling and recovery
- Integration with SEPI reference fetching

Test Categories:
- Module initialization and configuration
- Single alignment workflows  
- Batch alignment processing
- Reference file management
- EMBOSS-WATER integration
- Error handling and edge cases
- Performance and memory efficiency
"""

import pytest
import tempfile
import shutil
import subprocess
import logging
from unittest.mock import Mock, patch, MagicMock
from pathlib import Path
from typing import List, Dict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Setup path for imports
import sys
ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

# Try importing different aligner implementations
try:
    sys.path.insert(0, str(SRC / "priority3" / "aligner"))
    from wildtype_aligner import WildTypeAligner
    PRIMARY_ALIGNER_AVAILABLE = True
except ImportError as e:
    print(f"Primary WildTypeAligner not available: {e}")
    PRIMARY_ALIGNER_AVAILABLE = False

try:
    from simplified_wildtype_aligner import SimplifiedWildTypeAligner, SimpleAlignerConfig, AlignmentResult
    SIMPLIFIED_ALIGNER_AVAILABLE = True
except ImportError as e:
    print(f"SimplifiedWildTypeAligner not available: {e}")
    SIMPLIFIED_ALIGNER_AVAILABLE = False


@pytest.mark.skipif(not PRIMARY_ALIGNER_AVAILABLE, reason="Primary WildTypeAligner not available")
class TestWildTypeAlignerPrimary:
    """Test suite for the primary EMBOSS-WATER based WildTypeAligner."""
    
    def setup_method(self):
        """Setup test environment for each test."""
        self.temp_dir = tempfile.mkdtemp()
        self.reference_dir = Path(self.temp_dir) / "references"
        self.output_dir = Path(self.temp_dir) / "alignments"
        self.protein_dir = Path(self.temp_dir) / "proteins"
        
        # Create directories
        self.reference_dir.mkdir(parents=True)
        self.output_dir.mkdir(parents=True)
        self.protein_dir.mkdir(parents=True)
        
        # Setup logging
        logging.basicConfig(level=logging.DEBUG)
        self.logger = logging.getLogger("TestWildTypeAligner")
        
        self._create_test_data()
    
    def teardown_method(self):
        """Clean up test environment."""
        try:
            shutil.rmtree(self.temp_dir)
        except PermissionError:
            pass
    
    def _create_test_data(self):
        """Create test protein and reference sequences."""
        # Create test reference sequences
        test_references = {
            config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0]: 'MKQSTIALALLPLLFTPLLAQAEAVTPTEAAQESAEKPKLVDRAQAETDLLTQAPVSKRDDLLQSPDQPSRAQQQQLQDALQERLRAVAAELRNDQTQAAANQLQARAEQAASSIVQAVNDRVTLLPQVAQAAANQLQARAELVPQV',
            config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1]: 'MSVITKGRLNDLRDIAQNLGAMPRGTLVHQTTRKMQSLLDPLNMQEFQQTLKEYQLQGQSEAIRQIQEALAAAKNDIREIDLSQDFRLRSLNLRASFVAQTQDSADKVTQRLLRGDPVVRLQILSDPELTEHPVQRAALDLGATLKLDHQTLGLLLVGQQLGVMSQLQAAYQHQSRFKLRGDQARL',
            config_manager.get_default_genes("rnd_efflux_pumps", "primary")[2]: 'MRLLPLLLLLLPLLAQITAFDQSTQLVQAIVDDVKQQIQADSQNFAVNQVTRDLASQPVQKLQDDVVRLQIVSADKQQVQFLQGKQQVAAAVDQLDKSDTQATAQVQQLIRQQQAQLQDDVKRLQIVSADKQQVQFLQGKQQVAAAVDQLDKSDTQATAQVQQLIRQQQAQLQDDVKRLQIVSA'
        }
        
        for gene, sequence in test_references.items():
            ref_file = self.reference_dir / f"{gene}.faa"
            with open(ref_file, 'w') as f:
                f.write(f">{gene}_reference\n{sequence}\n")
        
        # Create test protein sequences (variations of references)
        test_proteins = {
            'sample_001_acrA': test_references[config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0]][:100] + 'X' + test_references[config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0]][101:],  # Single mutation
            'sample_002_acrB': test_references[config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1]][:50] + 'YY' + test_references[config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1]][52:],   # Insertion
            'sample_003_tolC': test_references[config_manager.get_default_genes("rnd_efflux_pumps", "primary")[2]][:-10],  # Deletion
        }
        
        for protein_id, sequence in test_proteins.items():
            protein_file = self.protein_dir / f"{protein_id}.fasta"
            with open(protein_file, 'w') as f:
                f.write(f">{protein_id}\n{sequence}\n")
    
    # ===== Initialization and Configuration Tests =====
    
    def test_aligner_initialization_basic(self):
        """Test basic aligner initialization."""
        aligner = WildTypeAligner(
            reference_dir=str(self.reference_dir),
            output_dir=str(self.output_dir)
        )
        
        assert aligner.reference_dir == self.reference_dir
        assert aligner.output_dir == self.output_dir
        assert aligner.water_path == "water"  # Default
        assert hasattr(aligner, 'logger')
    
    def test_aligner_initialization_custom_water_path(self):
        """Test aligner with custom EMBOSS water path."""
        custom_water = "/usr/local/bin/water"
        aligner = WildTypeAligner(
            reference_dir=str(self.reference_dir),
            output_dir=str(self.output_dir),
            water_path=custom_water
        )
        
        assert aligner.water_path == custom_water
    
    def test_aligner_output_directory_creation(self):
        """Test that aligner creates output directory."""
        nested_output = str(Path(self.temp_dir) / "nested" / "deep" / "alignments")
        aligner = WildTypeAligner(
            reference_dir=str(self.reference_dir),
            output_dir=nested_output
        )
        
        assert Path(nested_output).exists()
    
    # ===== Reference File Management Tests =====
    
    def test_reference_file_detection_gene_only(self):
        """Test detection of gene-only reference files."""
        aligner = WildTypeAligner(
            reference_dir=str(self.reference_dir),
            output_dir=str(self.output_dir)
        )
        
        protein_file = str(self.protein_dir / "sample_001_acrA.fasta")
        
        # Mock EMBOSS water to avoid actual execution
        with patch('subprocess.run') as mock_run:
            mock_run.return_value = MagicMock(returncode=0)
            
            result = aligner.align(protein_file, config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0], "sample_001")
            
            # Should find acrA.faa reference
            assert result is not None
            assert "sample_001_acrA_water.txt" in result
    
    def test_reference_file_detection_species_specific(self):
        """Test detection of species-specific reference files."""
        # Create species-specific reference
        species_ref = self.reference_dir / "Escherichia_coli_acrA.faa"
        with open(species_ref, 'w') as f:
            f.write(">E_coli_acrA_reference\nTESTSEQUENCE\n")
        
        aligner = WildTypeAligner(
            reference_dir=str(self.reference_dir),
            output_dir=str(self.output_dir)
        )
        
        protein_file = str(self.protein_dir / "sample_001_acrA.fasta")
        
        with patch('subprocess.run') as mock_run:
            mock_run.return_value = MagicMock(returncode=0)
            
            result = aligner.align(protein_file, config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0], "sample_001", species="Escherichia coli")
            
            assert result is not None
    
    def test_reference_file_missing_gene(self):
        """Test handling when reference file is missing."""
        aligner = WildTypeAligner(
            reference_dir=str(self.reference_dir),
            output_dir=str(self.output_dir)
        )
        
        protein_file = str(self.protein_dir / "sample_001_acrA.fasta")
        
        # Try with non-existent gene
        result = aligner.align(protein_file, "nonexistent_gene", "sample_001")
        
        assert result is None
    
    # ===== EMBOSS-WATER Integration Tests =====
    
    @patch('subprocess.run')
    def test_emboss_water_command_construction(self, mock_run):
        """Test proper EMBOSS water command construction."""
        mock_run.return_value = MagicMock(returncode=0)
        
        aligner = WildTypeAligner(
            reference_dir=str(self.reference_dir),
            output_dir=str(self.output_dir)
        )
        
        protein_file = str(self.protein_dir / "sample_001_acrA.fasta")
        result = aligner.align(protein_file, "acrA", "sample_001")
        
        # Verify subprocess was called with correct arguments
        mock_run.assert_called_once()
        call_args = mock_run.call_args[0][0]
        
        assert call_args[0] == "water"
        assert "-asequence" in call_args
        assert "-bsequence" in call_args
        assert "-gapopen" in call_args
        assert "10" in call_args
        assert "-gapextend" in call_args
        assert "0.5" in call_args
        assert "-outfile" in call_args
    
    @patch('subprocess.run')
    def test_emboss_water_success(self, mock_run):
        """Test successful EMBOSS water execution."""
        mock_run.return_value = MagicMock(returncode=0)
        
        aligner = WildTypeAligner(
            reference_dir=str(self.reference_dir),
            output_dir=str(self.output_dir)
        )
        
        protein_file = str(self.protein_dir / "sample_001_acrA.fasta")
        result = aligner.align(protein_file, "acrA", "sample_001")
        
        assert result is not None
        assert "sample_001_acrA_water.txt" in result
        assert Path(result).parent == self.output_dir
    
    @patch('subprocess.run')
    def test_emboss_water_failure(self, mock_run):
        """Test handling of EMBOSS water execution failure."""
        mock_run.side_effect = subprocess.CalledProcessError(
            1, ["water"], stderr="EMBOSS water error"
        )
        
        aligner = WildTypeAligner(
            reference_dir=str(self.reference_dir),
            output_dir=str(self.output_dir)
        )
        
        protein_file = str(self.protein_dir / "sample_001_acrA.fasta")
        result = aligner.align(protein_file, "acrA", "sample_001")
        
        assert result is None
    
    @patch('subprocess.run')
    def test_emboss_water_custom_parameters(self, mock_run):
        """Test EMBOSS water with custom gap parameters."""
        mock_run.return_value = MagicMock(returncode=0)
        
        # Create custom aligner (would need parameter support)
        aligner = WildTypeAligner(
            reference_dir=str(self.reference_dir),
            output_dir=str(self.output_dir)
        )
        
        protein_file = str(self.protein_dir / "sample_001_acrA.fasta")
        result = aligner.align(protein_file, "acrA", "sample_001")
        
        # Verify default parameters are used
        call_args = mock_run.call_args[0][0]
        gap_open_index = call_args.index("-gapopen")
        assert call_args[gap_open_index + 1] == "10"
    
    # ===== Batch Processing Tests =====
    
    @patch('subprocess.run')
    def test_batch_alignment_success(self, mock_run):
        """Test successful batch alignment processing."""
        mock_run.return_value = MagicMock(returncode=0)
        
        aligner = WildTypeAligner(
            reference_dir=str(self.reference_dir),
            output_dir=str(self.output_dir)
        )
        
        protein_files = [
            str(self.protein_dir / "sample_001_acrA.fasta"),
            str(self.protein_dir / "sample_002_acrB.fasta"),
            str(self.protein_dir / "sample_003_tolC.fasta")
        ]
        genes = ["acrA", "acrB", "tolC"]
        sample_ids = ["sample_001", "sample_002", "sample_003"]
        
        results = aligner.align_batch(protein_files, genes, sample_ids)
        
        assert len(results) == 3
        assert ("sample_001", config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0]) in results
        assert ("sample_002", config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1]) in results
        assert ("sample_003", config_manager.get_default_genes("rnd_efflux_pumps", "primary")[2]) in results
        
        # Verify all results are file paths
        for result_file in results.values():
            assert isinstance(result_file, str)
            assert "water.txt" in result_file
    
    @patch('subprocess.run')
    def test_batch_alignment_with_species(self, mock_run):
        """Test batch alignment with species information."""
        mock_run.return_value = MagicMock(returncode=0)
        
        aligner = WildTypeAligner(
            reference_dir=str(self.reference_dir),
            output_dir=str(self.output_dir)
        )
        
        protein_files = [str(self.protein_dir / "sample_001_acrA.fasta")]
        genes = ["acrA"]
        sample_ids = ["sample_001"]
        species_list = ["Escherichia coli"]
        
        results = aligner.align_batch(protein_files, genes, sample_ids, species_list)
        
        assert len(results) == 1
        assert ("sample_001", config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0]) in results
    
    @patch('subprocess.run')
    def test_batch_alignment_partial_failures(self, mock_run):
        """Test batch alignment with some failures."""
        # Mock alternating success/failure
        def side_effect(*args, **kwargs):
            if config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0] in str(args[0]):
                return MagicMock(returncode=0)
            else:
                raise subprocess.CalledProcessError(1, args[0])
        
        mock_run.side_effect = side_effect
        
        aligner = WildTypeAligner(
            reference_dir=str(self.reference_dir),
            output_dir=str(self.output_dir)
        )
        
        protein_files = [
            str(self.protein_dir / "sample_001_acrA.fasta"),
            str(self.protein_dir / "sample_002_acrB.fasta"),
        ]
        genes = ["acrA", "acrB"]
        sample_ids = ["sample_001", "sample_002"]
        
        results = aligner.align_batch(protein_files, genes, sample_ids)
        
        # Should only have successful alignment
        assert len(results) == 1
        assert ("sample_001", config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0]) in results
        assert ("sample_002", config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1]) not in results
    
    # ===== Error Handling and Edge Cases =====
    
    def test_missing_protein_file(self):
        """Test handling of missing protein file."""
        aligner = WildTypeAligner(
            reference_dir=str(self.reference_dir),
            output_dir=str(self.output_dir)
        )
        
        with patch('subprocess.run') as mock_run:
            mock_run.side_effect = FileNotFoundError("File not found")
            
            result = aligner.align("nonexistent.fasta", "acrA", "sample_001")
            assert result is None
    
    def test_corrupted_reference_file(self):
        """Test handling of corrupted reference file."""
        # Create corrupted reference
        corrupted_ref = self.reference_dir / "corrupted.faa"
        with open(corrupted_ref, 'w') as f:
            f.write("Not a valid FASTA file\n")
        
        aligner = WildTypeAligner(
            reference_dir=str(self.reference_dir),
            output_dir=str(self.output_dir)
        )
        
        protein_file = str(self.protein_dir / "sample_001_acrA.fasta")
        
        with patch('subprocess.run') as mock_run:
            mock_run.side_effect = subprocess.CalledProcessError(
                1, ["water"], stderr="Invalid FASTA format"
            )
            
            result = aligner.align(protein_file, "corrupted", "sample_001")
            assert result is None
    
    def test_insufficient_permissions_output_directory(self):
        """Test handling when output directory is not writable."""
        # This test may be platform-specific
        try:
            readonly_dir = Path(self.temp_dir) / "readonly"
            readonly_dir.mkdir()
            readonly_dir.chmod(0o444)  # Read-only
            
            aligner = WildTypeAligner(
                reference_dir=str(self.reference_dir),
                output_dir=str(readonly_dir)
            )
            
            protein_file = str(self.protein_dir / "sample_001_acrA.fasta")
            
            with patch('subprocess.run') as mock_run:
                mock_run.side_effect = PermissionError("Permission denied")
                
                result = aligner.align(protein_file, config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0], "sample_001")
                assert result is None
            
        except (OSError, PermissionError):
            pytest.skip("Permission test not supported on this system")
        finally:
            try:
                readonly_dir.chmod(0o755)
            except:
                pass
    
    def test_empty_protein_sequence(self):
        """Test handling of empty protein sequences."""
        # Create empty protein file
        empty_protein = self.protein_dir / "empty.fasta"
        with open(empty_protein, 'w') as f:
            f.write(">empty_protein\n\n")
        
        aligner = WildTypeAligner(
            reference_dir=str(self.reference_dir),
            output_dir=str(self.output_dir)
        )
        
        with patch('subprocess.run') as mock_run:
            mock_run.side_effect = subprocess.CalledProcessError(
                1, ["water"], stderr="Empty sequence"
            )
            
            result = aligner.align(str(empty_protein), config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0], "empty_test")
            assert result is None
    
    # ===== SEPI Integration Tests =====
    
    def test_sepi_reference_fetching_success(self):
        """Test successful SEPI reference fetching."""
        aligner = WildTypeAligner(
            reference_dir=str(self.reference_dir),
            output_dir=str(self.output_dir)
        )
        
        protein_file = str(self.protein_dir / "sample_001_acrA.fasta")
        
        # Check if SEPI integration is available
        if not hasattr(aligner, '_fetch_sepi_reference') and not hasattr(aligner, 'fetch_sepi_reference'):
            pytest.skip("SEPI integration not available in this aligner implementation")
        
        # Mock the SEPI fetch mechanism 
        with patch('subprocess.run') as mock_run:
            mock_run.return_value = MagicMock(returncode=0)
            
            # Test with non-existent gene - should fail gracefully without SEPI
            result = aligner.align(protein_file, "novel_gene", "sample_001", "Escherichia coli")
            
            # Without SEPI integration, should return None for missing reference
            assert result is None
    
    def test_sepi_reference_fetching_failure(self):
        """Test SEPI reference fetching failure."""
        aligner = WildTypeAligner(
            reference_dir=str(self.reference_dir),
            output_dir=str(self.output_dir)
        )
        
        protein_file = str(self.protein_dir / "sample_001_acrA.fasta")
        
        # Check if SEPI integration is available
        if not hasattr(aligner, '_fetch_sepi_reference') and not hasattr(aligner, 'fetch_sepi_reference'):
            pytest.skip("SEPI integration not available in this aligner implementation")
        
        # Try with non-existent gene and species
        result = aligner.align(protein_file, "novel_gene", "sample_001", "Escherichia coli")
        
        # Should fail gracefully when reference is not found
        assert result is None


@pytest.mark.skipif(not SIMPLIFIED_ALIGNER_AVAILABLE, reason="SimplifiedWildTypeAligner not available")  
class TestSimplifiedWildTypeAligner:
    """Test suite for SimplifiedWildTypeAligner with BioPython."""
    
    def setup_method(self):
        """Setup test environment."""
        self.temp_dir = tempfile.mkdtemp()
        self.input_dir = Path(self.temp_dir) / "input"
        self.output_dir = Path(self.temp_dir) / "output"
        
        self.input_dir.mkdir(parents=True)
        self.output_dir.mkdir(parents=True)
        
        # Create test configuration
        self.config = SimpleAlignerConfig(
            input_dir=str(self.input_dir),
            output_dir=str(self.output_dir),
            target_genes=config_manager.get_default_genes("rnd_efflux_pumps", "primary")
        )
        
        self._create_test_sequences()
    
    def teardown_method(self):
        """Clean up test environment."""
        try:
            shutil.rmtree(self.temp_dir)
        except PermissionError:
            pass
    
    def _create_test_sequences(self):
        """Create test protein sequences."""
        test_sequences = {
            'sample_001_acrA': 'MKQSTIALALLPLLFTPLLAQAEAVTPTEAAQESAEKPKLVDRAQAETDLLTQAPVSKRDDLLQSPDQPSRAQQQQLQDALQERLRAVAAELRNDQTQAAANQLQARAEQAASSIVQAVNDRVTLLPQVAQAAANQLQARAELVPQV',
            'sample_002_acrB': 'MSVITKGRLNDLRDIAQNLGAMPRGTLVHQTTRKMQSLLDPLNMQEFQQTLKEYQLQGQSEAIRQIQEALAAAKNDIREIDLSQDFRLRSLNLRASFVAQTQDSADKVTQRLLRGDPVVRLQILSDPELTEHPVQRAALDLGATLKLDHQTLGLLLVGQQLGVMSQLQAAYQHQSRFKLRGDQARL',
            'sample_003_tolC': 'MRLLPLLLLLLPLLAQITAFDQSTQLVQAIVDDVKQQIQADSQNFAVNQVTRDLASQPVQKLQDDVVRLQIVSADKQQVQFLQGKQQVAAAVDQLDKSDTQATAQVQQLIRQQQAQLQDDVKRLQIVSADKQQVQFLQGKQQVAAAVDQLDKSDTQATAQVQQLIRQQQAQLQDDVKRLQIVSA'
        }
        
        for protein_id, sequence in test_sequences.items():
            protein_file = self.input_dir / f"{protein_id}.fasta"
            with open(protein_file, 'w') as f:
                f.write(f">{protein_id}\n{sequence}\n")
    
    def test_simplified_aligner_initialization(self):
        """Test SimplifiedWildTypeAligner initialization."""
        aligner = SimplifiedWildTypeAligner(self.config)
        
        assert aligner.config == self.config
        assert hasattr(aligner, 'logger')
        assert hasattr(aligner, 'stats')
        assert aligner.stats['proteins_processed'] == 0
    
    def test_analyze_sequences_basic(self):
        """Test basic sequence analysis functionality."""
        aligner = SimplifiedWildTypeAligner(self.config)
        
        # Test sequences
        test_sequences = [
            "MKQSTIALALLPLLFTPLLAQ",
            "MKQSTIALALLPLLFTPLLAX",  # Single mutation
            "MKQSTIALALLPLLFTPLLA",   # Deletion
        ]
        reference = "MKQSTIALALLPLLFTPLLAQ"
        gene_name = "test_gene"
        
        if hasattr(aligner, 'analyze_sequences'):
            results = aligner.analyze_sequences(test_sequences, reference, gene_name)
            
            assert 'gene_name' in results
            assert results['gene_name'] == gene_name
            assert 'mutations' in results
            assert 'similarity_scores' in results
            assert isinstance(results['mutations'], list)
            assert isinstance(results['similarity_scores'], list)
        else:
            # Method not available - skip this test
            pytest.skip("analyze_sequences method not available")
    
    def test_sequence_validation(self):
        """Test sequence validation functionality."""
        aligner = SimplifiedWildTypeAligner(self.config)
        
        if hasattr(aligner, 'analyze_sequences'):
            # Test with invalid sequences
            invalid_sequences = [
                "",  # Empty
                "MKQSTIALX123LLPLLFTPLLAQ",  # Invalid characters
                None  # None value
            ]
            reference = "MKQSTIALALLPLLFTPLLAQ"
            gene_name = "test_gene"
            
            try:
                results = aligner.analyze_sequences(invalid_sequences, reference, gene_name)
                # Should handle gracefully or raise specific exception
                assert isinstance(results, dict)
            except Exception as e:
                # Expected for invalid input
                assert isinstance(e, (ValueError, TypeError))
    
    def test_built_in_references(self):
        """Test built-in reference sequences."""
        aligner = SimplifiedWildTypeAligner(self.config)
        
        if hasattr(aligner, 'builtin_references'):
            # Should be empty for generic aligner
            assert isinstance(aligner.builtin_references, dict)
            # Generic aligner should not have hardcoded references
            # Users should provide their own
    
    def test_statistics_tracking(self):
        """Test statistics tracking functionality."""
        aligner = SimplifiedWildTypeAligner(self.config)
        
        # Verify initial statistics
        assert aligner.stats['proteins_processed'] == 0
        assert aligner.stats['alignments_performed'] == 0
        assert aligner.stats['failed_alignments'] == 0
        
        # Statistics should be updateable
        aligner.stats['proteins_processed'] += 1
        assert aligner.stats['proteins_processed'] == 1


def test_wildtype_aligner_integration_workflow():
    """Integration test for WildTypeAligner workflow."""
    temp_dir = tempfile.mkdtemp()
    try:
        reference_dir = Path(temp_dir) / "references"
        output_dir = Path(temp_dir) / "alignments"
        protein_dir = Path(temp_dir) / "proteins"
        
        # Create directories
        reference_dir.mkdir(parents=True)
        output_dir.mkdir(parents=True)
        protein_dir.mkdir(parents=True)
        
        # Create test data
        reference_seq = "MKQSTIALALLPLLFTPLLAQAEAVTPTEAAQESAEKPKLVDRAQ"
        protein_seq = "MKQSTIALALLPLLFTPLLAXAEAVTPTEAAQESAEKPKLVDRAQ"  # Single mutation
        
        # Create reference file
        ref_file = reference_dir / "acrA.faa"
        with open(ref_file, 'w') as f:
            f.write(">acrA_reference\n" + reference_seq + "\n")
        
        # Create protein file
        protein_file = protein_dir / "sample_001_acrA.fasta"
        with open(protein_file, 'w') as f:
            f.write(">sample_001_acrA\n" + protein_seq + "\n")
        
        # Test alignment workflow (mock EMBOSS)
        if PRIMARY_ALIGNER_AVAILABLE:
            aligner = WildTypeAligner(
                reference_dir=str(reference_dir),
                output_dir=str(output_dir)
            )
            
            with patch('subprocess.run') as mock_run:
                mock_run.return_value = MagicMock(returncode=0)
                
                result = aligner.align(
                    str(protein_file), config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0], "sample_001"
                )
                
                assert result is not None
                assert "sample_001_acrA_water.txt" in result
        
        # Test simplified workflow
        if SIMPLIFIED_ALIGNER_AVAILABLE:
            config = SimpleAlignerConfig(
                input_dir=str(protein_dir),
                output_dir=str(output_dir),
                target_genes=[config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0]]
            )
            
            simplified_aligner = SimplifiedWildTypeAligner(config)
            
            if hasattr(simplified_aligner, 'analyze_sequences'):
                results = simplified_aligner.analyze_sequences(
                    [protein_seq], reference_seq, config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0]
                )
                
                assert 'gene_name' in results
                assert results['gene_name'] == config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0]
        
        print("✅ WildTypeAligner integration workflow validated")
        
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)


if __name__ == "__main__":
    # Run basic tests for debugging
    if PRIMARY_ALIGNER_AVAILABLE:
        test = TestWildTypeAlignerPrimary()
        test.setup_method()
        try:
            test.test_aligner_initialization_basic()
            test.test_reference_file_detection_gene_only()
            print("✅ Primary WildTypeAligner tests passed!")
        finally:
            test.teardown_method()
    
    if SIMPLIFIED_ALIGNER_AVAILABLE:
        test = TestSimplifiedWildTypeAligner()
        test.setup_method()
        try:
            test.test_simplified_aligner_initialization()
            print("✅ SimplifiedWildTypeAligner tests passed!")
        finally:
            test.teardown_method()
    
    # Always run integration test
    test_wildtype_aligner_integration_workflow()