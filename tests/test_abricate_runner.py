#!/usr/bin/env python3
"""
Unit tests for src/abricate_runner.py

This test suite validates the abricate_runner module's core functionality including:
- Command construction and subprocess calls
- File discovery and validation  
- Output filename generation
- Error handling and statistics
- Multi-database support

Uses pytest-mock to patch subprocess calls and file system operations for 
isolated unit testing without requiring actual abricate installation.
"""

import pytest
import subprocess
import tempfile
from pathlib import Path
from unittest.mock import Mock, call
import sys
import os

# Add src directory to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from abricate_runner import (
    run_abricate_on_file,
    find_fastas,
    is_fasta,
    infer_genome_id,
    run_abricate,
    main,
    ScanResults,
    FASTA_EXTENSIONS
)


class TestAbricateRunnerCore:
    """Test core functionality of abricate_runner module"""
    
    def test_is_fasta_valid_extensions(self, tmp_path):
        """Test that is_fasta correctly identifies FASTA files by extension"""
        # Create test files with various extensions
        valid_files = []
        for ext in ['.fasta', '.fa', '.fna', '.FASTA', '.FA']:  # Include uppercase
            test_file = tmp_path / f"genome{ext}"
            test_file.touch()
            valid_files.append(test_file)
            
        invalid_files = []
        for ext in ['.txt', '.tsv', '.csv', '.fastq', '.fq']:
            test_file = tmp_path / f"genome{ext}"
            test_file.touch()
            invalid_files.append(test_file)
        
        # Test valid FASTA files
        for fasta_file in valid_files:
            assert is_fasta(fasta_file), f"Should recognize {fasta_file.suffix} as FASTA"
            
        # Test invalid files
        for invalid_file in invalid_files:
            assert not is_fasta(invalid_file), f"Should not recognize {invalid_file.suffix} as FASTA"
            
        # Test non-existent file
        non_existent = tmp_path / "missing.fasta"
        assert not is_fasta(non_existent), "Should return False for non-existent files"

    def test_find_fastas_discovery(self, tmp_path):
        """Test that find_fastas correctly discovers FASTA files in directory"""
        # Create mix of FASTA and non-FASTA files
        fasta_files = []
        for i, ext in enumerate(['.fasta', '.fa', '.fna']):
            fasta_file = tmp_path / f"genome_{i}{ext}"
            fasta_file.touch()
            fasta_files.append(fasta_file)
            
        # Create non-FASTA files that should be ignored
        (tmp_path / "readme.txt").touch()
        (tmp_path / "results.tsv").touch()
        
        # Test discovery
        discovered = find_fastas(tmp_path)
        
        assert len(discovered) == 3, "Should find exactly 3 FASTA files"
        for fasta_file in fasta_files:
            assert fasta_file in discovered, f"Should discover {fasta_file.name}"
            
        # Test sorting (should be alphabetical)
        assert discovered == sorted(discovered), "Results should be sorted alphabetically"

    def test_infer_genome_id_extraction(self):
        """Test genome ID extraction from various filename patterns"""
        test_cases = [
            ("genome_123.fasta", "genome_123"),
            ("E_coli_K12.fa", "E_coli_K12"),
            ("GCF_000005825.2_ASM582v2_genomic.fna", "GCF_000005825.2_ASM582v2_genomic"),
            ("sample.FASTA", "sample"),
            ("/path/to/complex_name.v2.fasta", "complex_name.v2")
        ]
        
        for filename, expected_id in test_cases:
            path = Path(filename)
            result = infer_genome_id(path)
            assert result == expected_id, f"Expected '{expected_id}' for '{filename}', got '{result}'"


class TestAbricateSubprocessIntegration:
    """Test subprocess integration with proper mocking"""
    
    def test_run_abricate_command_construction(self, mocker):
        """Test that abricate commands are constructed correctly"""
        # Mock subprocess.run
        mock_run = mocker.patch('abricate_runner.subprocess.run')
        mock_run.return_value = Mock(
            returncode=0,
            stdout="SEQUENCE\tSTART\tEND\tSTRAND\tGENE\tCOVERAGE\tIDENTITY\tDATABASE\n",
            stderr=""
        )
        
        test_fasta = Path("/test/genome.fasta")
        
        # Test default CARD database
        result = run_abricate_on_file(test_fasta)
        
        expected_cmd = ["abricate", "--db", "card", "--nopath", str(test_fasta)]
        mock_run.assert_called_once_with(
            expected_cmd,
            capture_output=True,
            text=True
        )
        assert "SEQUENCE\tSTART\tEND" in result
        
    def test_run_abricate_different_databases(self, mocker):
        """Test command construction for different databases"""
        mock_run = mocker.patch('abricate_runner.subprocess.run')
        mock_run.return_value = Mock(returncode=0, stdout="mock_output", stderr="")
        
        test_fasta = Path("/test/genome.fasta")
        
        # Test different databases
        databases = ["vfdb", "plasmidfinder", "resfinder", "ncbi"]
        for db in databases:
            mock_run.reset_mock()
            run_abricate_on_file(test_fasta, db=db)
            
            expected_cmd = ["abricate", "--db", db, "--nopath", str(test_fasta)]
            mock_run.assert_called_once_with(
                expected_cmd,
                capture_output=True,
                text=True
            )
    
    def test_run_abricate_nopath_flag(self, mocker):
        """Test nopath flag behavior"""
        mock_run = mocker.patch('abricate_runner.subprocess.run')
        mock_run.return_value = Mock(returncode=0, stdout="output", stderr="")
        
        test_fasta = Path("/test/genome.fasta")
        
        # Test with nopath=True (default)
        run_abricate_on_file(test_fasta, nopath=True)
        args, kwargs = mock_run.call_args
        assert "--nopath" in args[0]
        
        # Test with nopath=False
        mock_run.reset_mock()
        run_abricate_on_file(test_fasta, nopath=False)
        args, kwargs = mock_run.call_args
        assert "--nopath" not in args[0]

    def test_run_abricate_error_handling(self, mocker):
        """Test error handling for subprocess failures"""
        # Test FileNotFoundError (abricate not installed)
        mocker.patch('abricate_runner.subprocess.run', side_effect=FileNotFoundError())
        
        with pytest.raises(RuntimeError, match="Abricate executable not found"):
            run_abricate_on_file(Path("/test/genome.fasta"))
            
        # Test non-zero return code
        mock_run = mocker.patch('abricate_runner.subprocess.run')
        mock_run.return_value = Mock(
            returncode=1,
            stdout="",
            stderr="Database error"
        )
        
        result = run_abricate_on_file(Path("/test/genome.fasta"))
        assert result == "", "Should return empty string on failure"


class TestAbricateBatchProcessing:
    """Test batch processing functionality"""
    
    def test_run_abricate_batch_success(self, mocker, tmp_path):
        """Test successful batch processing of multiple genomes"""
        # Create test FASTA files
        fasta_files = []
        for i in range(3):
            fasta_file = tmp_path / f"genome_{i}.fasta"
            fasta_file.write_text(f">genome_{i}\nATCGATCG\n")
            fasta_files.append(fasta_file)
            
        output_dir = tmp_path / "output"
        
        # Mock successful abricate runs
        mock_run = mocker.patch('abricate_runner.subprocess.run')
        mock_run.return_value = Mock(
            returncode=0,
            stdout="SEQUENCE\tSTART\tEND\tSTRAND\tGENE\tCOVERAGE\tIDENTITY\tDATABASE\n" +
                   "genome_0\t100\t200\t+\tacrA\t95.5\t98.2\tcard\n",
            stderr=""
        )
        
        # Mock file writing
        mocker.patch('abricate_runner.write_report')
        
        # Run batch processing
        results = run_abricate(str(tmp_path), str(output_dir), db="card")
        
        # Verify results
        assert isinstance(results, ScanResults)
        assert results.database == "card"
        assert results.total_genomes == 3
        assert len(results.successful_scans) == 3
        assert len(results.failed_genomes) == 0
        assert len(results.no_hits_genomes) == 0
        assert results.global_error is None
        
        # Verify subprocess calls
        assert mock_run.call_count == 3
        
    def test_run_abricate_mixed_results(self, mocker, tmp_path):
        """Test batch processing with mixed success/failure results"""
        # Create test files
        (tmp_path / "success.fasta").write_text(">success\nATCG\n")
        (tmp_path / "no_hits.fasta").write_text(">no_hits\nATCG\n") 
        (tmp_path / "failure.fasta").write_text(">failure\nATCG\n")
        
        output_dir = tmp_path / "output"
        
        # Mock different abricate outcomes
        def mock_run_side_effect(cmd, **kwargs):
            if "success.fasta" in str(cmd):
                return Mock(returncode=0, stdout="gene_hit_data\n", stderr="")
            elif "no_hits.fasta" in str(cmd):
                return Mock(returncode=0, stdout="", stderr="")
            else:  # failure.fasta
                return Mock(returncode=1, stdout="", stderr="error")
                
        mocker.patch('abricate_runner.subprocess.run', side_effect=mock_run_side_effect)
        mocker.patch('abricate_runner.write_report')
        
        results = run_abricate(str(tmp_path), str(output_dir), db="vfdb")
        
        assert results.database == "vfdb"
        assert len(results.successful_scans) == 1
        assert len(results.no_hits_genomes) == 1  
        assert len(results.failed_genomes) == 1
        assert "failure" in results.failed_genomes

    def test_output_filename_generation(self, mocker, tmp_path):
        """Test that output filenames are generated correctly"""
        # Create test file
        (tmp_path / "E_coli_K12.fasta").write_text(">test\nATCG\n")
        output_dir = tmp_path / "output"
        
        # Mock successful run
        mocker.patch('abricate_runner.subprocess.run', return_value=Mock(
            returncode=0, stdout="gene_data\n", stderr=""
        ))
        
        # Mock write_report to capture the output path
        written_files = []
        def mock_write(content, path):
            written_files.append(path)
            
        mocker.patch('abricate_runner.write_report', side_effect=mock_write)
        
        # Test different databases
        for db in ["card", "vfdb", "plasmidfinder"]:
            written_files.clear()
            run_abricate(str(tmp_path), str(output_dir), db=db)
            
            assert len(written_files) == 1
            output_file = written_files[0]
            expected_name = f"E_coli_K12_{db}.tsv"
            assert output_file.name == expected_name, f"Expected {expected_name}, got {output_file.name}"


class TestAbricateMainFunction:
    """Test the main CLI function"""
    
    def test_main_argument_parsing(self, mocker):
        """Test main function argument parsing and execution"""
        # Mock the run_abricate function
        mock_results = ScanResults(database="card")
        mock_results.successful_scans = ["genome1", "genome2"]
        mock_results.failed_genomes = {}
        mock_results.no_hits_genomes = []
        mock_results.output_files = [Path("out1.tsv"), Path("out2.tsv")]
        
        mock_run = mocker.patch('abricate_runner.run_abricate', return_value=mock_results)
        
        # Test successful execution
        argv = ["--input-dir", "/input", "--output-dir", "/output", "--db", "vfdb"]
        result = main(argv)
        
        assert result == 0  # Success
        mock_run.assert_called_once_with("/input", "/output", db="vfdb")
        
    def test_main_with_failures(self, mocker, capsys):
        """Test main function with processing failures"""
        # Mock results with failures
        mock_results = ScanResults(database="card")
        mock_results.successful_scans = ["genome1"]
        mock_results.failed_genomes = {"genome2": "processing error"}
        mock_results.no_hits_genomes = []
        mock_results.output_files = [Path("out1.tsv")]
        
        mocker.patch('abricate_runner.run_abricate', return_value=mock_results)
        
        argv = ["--input-dir", "/input", "--output-dir", "/output"]
        result = main(argv)
        
        assert result == 1  # Failure due to failed genomes
        
        # Check output contains summary
        captured = capsys.readouterr()
        assert "1 genomes with hits" in captured.out
        assert "1 genomes failed" in captured.out
        
    def test_main_global_error(self, mocker, capsys):
        """Test main function with global errors"""
        # Mock results with global error
        mock_results = ScanResults(database="card")
        mock_results.global_error = "Input directory not found"
        
        mocker.patch('abricate_runner.run_abricate', return_value=mock_results)
        
        argv = ["--input-dir", "/nonexistent", "--output-dir", "/output"]
        result = main(argv)
        
        assert result == 1
        captured = capsys.readouterr()
        assert "Input directory not found" in captured.err

    def test_main_exception_handling(self, mocker, capsys):
        """Test main function exception handling"""
        # Mock function to raise exception
        mocker.patch('abricate_runner.run_abricate', side_effect=Exception("Unexpected error"))
        
        argv = ["--input-dir", "/input", "--output-dir", "/output"]
        result = main(argv)
        
        assert result == 1
        captured = capsys.readouterr()
        assert "Unexpected error" in captured.err


class TestScanResults:
    """Test ScanResults dataclass functionality"""
    
    def test_scan_results_initialization(self):
        """Test ScanResults initialization and default values"""
        results = ScanResults(database="test_db")
        
        assert results.database == "test_db"
        assert results.successful_scans == []
        assert results.failed_genomes == {}
        assert results.no_hits_genomes == []
        assert results.output_files == []
        assert results.total_genomes == 0
        assert results.global_error is None
        
    def test_scan_results_data_tracking(self):
        """Test ScanResults data tracking capabilities"""
        results = ScanResults(database="card")
        
        # Add test data
        results.successful_scans.append("genome1")
        results.failed_genomes["genome2"] = "error message"
        results.no_hits_genomes.append("genome3")
        results.output_files.append(Path("output.tsv"))
        results.total_genomes = 3
        
        # Verify data is tracked correctly
        assert len(results.successful_scans) == 1
        assert len(results.failed_genomes) == 1
        assert len(results.no_hits_genomes) == 1
        assert len(results.output_files) == 1
        assert results.total_genomes == 3
        assert "genome2" in results.failed_genomes


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
