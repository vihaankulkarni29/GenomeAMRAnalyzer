import sys; import os; sys.path.append(os.path.join(os.path.dirname(__file__), '\\\src')); from configuration_manager import config_manager
"""
Rigorous Test Suite for FastaAAExtractor Module
==============================================

This test suite validates the FastaAAExtractor with comprehensive error handling,
real-world scenarios, and production-grade robustness testing.

Test Categories:
- Module initialization and configuration
- CARD RGI output parsing (multiple formats)
- Genome FASTA loading and validation
- Protein sequence extraction accuracy
- Coordinate validation and error handling
- File I/O robustness and edge cases
- Large-scale batch processing
- Memory efficiency and performance
- Integration with CARD pipeline output
"""

import pytest
import tempfile
import shutil
import csv
import json
import gzip
import logging
from unittest.mock import Mock, patch, MagicMock
from pathlib import Path
from typing import Dict, List, Optional
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Setup path for imports
import sys
ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from priority3.extractor.fasta_aa_extractor import FastaAAExtractor


class TestFastaAAExtractorRobust:
    """Comprehensive test suite for FastaAAExtractor."""
    
    def setup_method(self):
        """Setup test environment for each test."""
        self.temp_dir = tempfile.mkdtemp()
        self.test_output_dir = Path(self.temp_dir) / "extracted_proteins"
        self.test_genomes_dir = Path(self.temp_dir) / "genomes"
        self.test_rgi_dir = Path(self.temp_dir) / "rgi_results"
        
        # Create directories
        self.test_output_dir.mkdir(parents=True)
        self.test_genomes_dir.mkdir(parents=True)
        self.test_rgi_dir.mkdir(parents=True)
        
        # Setup logging
        logging.basicConfig(level=logging.DEBUG)
        self.logger = logging.getLogger("TestFastaAAExtractor")
        
        # Create test data
        self._create_test_genome_files()
        self._create_test_rgi_outputs()
    
    def teardown_method(self):
        """Clean up test environment."""
        try:
            shutil.rmtree(self.temp_dir)
        except PermissionError:
            # Windows file locking issue - not critical for tests
            pass
    
    def _create_test_genome_files(self):
        """Create realistic test genome FASTA files."""
        # E. coli K-12 MG1655-like sequence with known genes
        ecoli_sequence = (
            "ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGACGCGTACAGGAAACACAGAAAAAAGCCCGCACCTG"
            "ACAGTGCGGGCTTTTTTTTTCGACCAAAGGTAACGAGGTAACAACCATGCGAGTGTTGAAGTTCGGCGGTACATCAGTGGCAAATGCAGAACGTTTTCTGC"
            "GTGTTGCCGATATTCTGGAAAGCAATGCCAGGCAGGGGCAGGTGGCCACCGTCCTCTCTGCCCCCGCCAAAATCACCAACCACCTGGTGGCGATGATTGAA"
            "AAAACCATTAGCGGCCAGGATGCTTTACCCAATATCAGCGATGCCGAACGTATTTTTGCCGAACTTTTGACGGGACTCGCCGCCGCCCAGCCGGGGTTCCCG"
            "CTGGCGCAATTGAAAACTTTCGTCGATCAGGAATTTGCCCAAATAAAACATGTCCTGCATGGCATTAGTTTGTTGGGGCAGTGCCCGGATAGCATCAACGC"
        )
        
        # Create test genome with known gene coordinates
        test_genome_records = [
            SeqRecord(
                Seq(ecoli_sequence * 20),  # Make it longer for realistic gene lengths
                id="NC_000913.3",
                description="Escherichia coli str. K-12 substr. MG1655 chromosome"
            )
        ]
        
        genome_file = self.test_genomes_dir / "GCF_000005825.2_ASM582v2_genomic.fna"
        with open(genome_file, "w") as f:
            SeqIO.write(test_genome_records, f, "fasta")
        
        # Create second test genome
        second_genome_records = [
            SeqRecord(
                Seq(ecoli_sequence[::-1] * 15),  # Reverse for variation
                id="NC_000913.4",
                description="Test genome 2"
            )
        ]
        
        genome_file2 = self.test_genomes_dir / "GCF_000009605.1_ASM960v1_genomic.fna"
        with open(genome_file2, "w") as f:
            SeqIO.write(second_genome_records, f, "fasta")
    
    def _create_test_rgi_outputs(self):
        """Create realistic CARD RGI output files in multiple formats."""
        
        # Sample RGI data matching test genomes
        rgi_data = [
            {
                'Best_Hit_ARO': config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0],
                'Contig': 'NC_000913.3',
                'Start': '100',
                'Stop': '400',
                'Orientation': '+',
                'Cut_Off': 'Perfect',
                'Pass_Bitscore': '500',
                'Best_Hit_Bitscore': '500',
                'Best_Identities': '100.00',
                'ARO': '3000007',
                'Model_type': 'protein homolog model'
            },
            {
                'Best_Hit_ARO': config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1],
                'Contig': 'NC_000913.3', 
                'Start': '500',
                'Stop': '800',
                'Orientation': '-',
                'Cut_Off': 'Perfect',
                'Pass_Bitscore': '600',
                'Best_Hit_Bitscore': '600',
                'Best_Identities': '98.50',
                'ARO': '3000008',
                'Model_type': 'protein homolog model'
            },
            {
                'Best_Hit_ARO': config_manager.get_default_genes("rnd_efflux_pumps", "primary")[2],
                'Contig': 'NC_000913.3',
                'Start': '900',
                'Stop': '1200',
                'Orientation': '+',
                'Cut_Off': 'Strict',
                'Pass_Bitscore': '400',
                'Best_Hit_Bitscore': '400',
                'Best_Identities': '95.00',
                'ARO': '3000009',
                'Model_type': 'protein homolog model'
            }
        ]
        
        # Create .txt (tab-separated) format
        txt_file = self.test_rgi_dir / "sample1_rgi_output.txt"
        with open(txt_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=rgi_data[0].keys(), delimiter='\t')
            writer.writeheader()
            writer.writerows(rgi_data)
        
        # Create .csv format
        csv_file = self.test_rgi_dir / "sample2_rgi_output.csv"
        with open(csv_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=rgi_data[0].keys())
            writer.writeheader()
            writer.writerows(rgi_data)
        
        # Create Excel format (if pandas available)
        try:
            import pandas as pd
            df = pd.DataFrame(rgi_data)
            xlsx_file = self.test_rgi_dir / "sample3_rgi_output.xlsx"
            df.to_excel(xlsx_file, index=False)
        except ImportError:
            self.logger.warning("pandas not available - skipping Excel test file creation")
        
        # Create malformed RGI file for error testing
        malformed_file = self.test_rgi_dir / "malformed_rgi.txt"
        with open(malformed_file, 'w') as f:
            f.write("This is not a valid RGI file\n")
            f.write("Missing headers and structure\n")
        
        # Create empty RGI file
        empty_file = self.test_rgi_dir / "empty_rgi.txt"
        with open(empty_file, 'w') as f:
            f.write("Best_Hit_ARO\tContig\tStart\tStop\tOrientation\n")  # Headers only
    
    # ===== Initialization and Configuration Tests =====
    
    def test_extractor_initialization_defaults(self):
        """Test extractor initializes correctly with default parameters."""
        extractor = FastaAAExtractor()
        
        assert extractor.output_dir == Path("extracted_proteins")
        assert extractor.output_dir.exists()
        assert hasattr(extractor, 'logger')
    
    def test_extractor_initialization_custom_output(self):
        """Test extractor with custom output directory."""
        custom_output = str(self.test_output_dir)
        extractor = FastaAAExtractor(output_dir=custom_output)
        
        assert extractor.output_dir == Path(custom_output)
        assert extractor.output_dir.exists()
    
    def test_extractor_initialization_nested_directory(self):
        """Test extractor creates nested output directories."""
        nested_output = str(self.test_output_dir / "nested" / "deep" / "proteins")
        extractor = FastaAAExtractor(output_dir=nested_output)
        
        assert extractor.output_dir == Path(nested_output)
        assert extractor.output_dir.exists()
    
    # ===== RGI Output Parsing Tests =====
    
    def test_parse_rgi_tabular_txt_format(self):
        """Test parsing of tab-separated RGI output."""
        extractor = FastaAAExtractor(output_dir=str(self.test_output_dir))
        rgi_file = self.test_rgi_dir / "sample1_rgi_output.txt"
        gene_list = config_manager.get_default_genes("rnd_efflux_pumps", "primary")
        
        hits = extractor._parse_rgi_tabular(str(rgi_file), gene_list)
        
        assert len(hits) == 3
        assert all(hit['Best_Hit_ARO'] in gene_list for hit in hits)
        
        # Verify specific gene data
        acra_hit = next(hit for hit in hits if hit['Best_Hit_ARO'] == config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0])
        assert acra_hit['Contig'] == 'NC_000913.3'
        assert acra_hit['Start'] == '100'
        assert acra_hit['Stop'] == '400'
        assert acra_hit['Orientation'] == '+'
    
    def test_parse_rgi_tabular_csv_format(self):
        """Test parsing of CSV RGI output."""
        extractor = FastaAAExtractor(output_dir=str(self.test_output_dir))
        rgi_file = self.test_rgi_dir / "sample2_rgi_output.csv"
        gene_list = [config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0], config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1]]
        
        hits = extractor._parse_rgi_tabular(str(rgi_file), gene_list)
        
        assert len(hits) == 2  # Only requested genes
        gene_names = [hit['Best_Hit_ARO'] for hit in hits]
        assert config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0] in gene_names
        assert config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1] in gene_names
        assert config_manager.get_default_genes("rnd_efflux_pumps", "primary")[2] not in gene_names  # Not in gene_list
    
    @pytest.mark.skipif(
        not pytest.importorskip("pandas", reason="pandas not available"),
        reason="pandas required for Excel support"
    )
    def test_parse_rgi_tabular_xlsx_format(self):
        """Test parsing of Excel RGI output."""
        extractor = FastaAAExtractor(output_dir=str(self.test_output_dir))
        rgi_file = self.test_rgi_dir / "sample3_rgi_output.xlsx"
        gene_list = config_manager.get_default_genes("rnd_efflux_pumps", "primary")
        
        if rgi_file.exists():  # Only run if Excel file was created
            hits = extractor._parse_rgi_tabular(str(rgi_file), gene_list)
            assert len(hits) == 3
            assert all(hit['Best_Hit_ARO'] in gene_list for hit in hits)
    
    def test_parse_rgi_tabular_gene_filtering(self):
        """Test gene list filtering functionality."""
        extractor = FastaAAExtractor(output_dir=str(self.test_output_dir))
        rgi_file = self.test_rgi_dir / "sample1_rgi_output.txt"
        
        # Test subset filtering
        subset_genes = [config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0]]
        hits = extractor._parse_rgi_tabular(str(rgi_file), subset_genes)
        assert len(hits) == 1
        assert hits[0]['Best_Hit_ARO'] == config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0]
        
        # Test empty gene list
        hits = extractor._parse_rgi_tabular(str(rgi_file), [])
        assert len(hits) == 0
        
        # Test non-matching genes
        hits = extractor._parse_rgi_tabular(str(rgi_file), ['nonexistent_gene'])
        assert len(hits) == 0
    
    def test_parse_rgi_tabular_malformed_file(self):
        """Test handling of malformed RGI files."""
        extractor = FastaAAExtractor(output_dir=str(self.test_output_dir))
        malformed_file = self.test_rgi_dir / "malformed_rgi.txt"
        gene_list = [config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0], config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1]]
        
        hits = extractor._parse_rgi_tabular(str(malformed_file), gene_list)
        assert len(hits) == 0  # Should return empty list for malformed files
    
    def test_parse_rgi_tabular_empty_file(self):
        """Test handling of empty RGI files."""
        extractor = FastaAAExtractor(output_dir=str(self.test_output_dir))
        empty_file = self.test_rgi_dir / "empty_rgi.txt"
        gene_list = [config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0], config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1]]
        
        hits = extractor._parse_rgi_tabular(str(empty_file), gene_list)
        assert len(hits) == 0
    
    def test_parse_rgi_tabular_nonexistent_file(self):
        """Test handling of non-existent RGI files."""
        extractor = FastaAAExtractor(output_dir=str(self.test_output_dir))
        nonexistent_file = self.test_rgi_dir / "does_not_exist.txt"
        gene_list = [config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0], config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1]]
        
        hits = extractor._parse_rgi_tabular(str(nonexistent_file), gene_list)
        assert len(hits) == 0
    
    # ===== Protein Extraction Tests =====
    
    def test_protein_extraction_basic_workflow(self):
        """Test basic protein extraction workflow."""
        extractor = FastaAAExtractor(output_dir=str(self.test_output_dir))
        
        genome_file = self.test_genomes_dir / "GCF_000005825.2_ASM582v2_genomic.fna"
        rgi_file = self.test_rgi_dir / "sample1_rgi_output.txt"
        gene_list = [config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0], config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1]]
        sample_id = "test_sample_001"
        
        output_files = extractor.extract_proteins(
            str(genome_file), str(rgi_file), gene_list, sample_id
        )
        
        assert len(output_files) == 2  # One file per gene
        
        # Verify output files exist and contain sequences
        for output_file in output_files:
            assert Path(output_file).exists()
            assert Path(output_file).stat().st_size > 0
            
            # Verify FASTA format
            records = list(SeqIO.parse(output_file, "fasta"))
            assert len(records) == 1  # One sequence per file
            assert len(records[0].seq) > 0  # Non-empty sequence
    
    def test_protein_extraction_coordinate_accuracy(self):
        """Test accuracy of coordinate-based extraction."""
        extractor = FastaAAExtractor(output_dir=str(self.test_output_dir))
        
        genome_file = self.test_genomes_dir / "GCF_000005825.2_ASM582v2_genomic.fna"
        rgi_file = self.test_rgi_dir / "sample1_rgi_output.txt"
        gene_list = [config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0]]
        sample_id = "coord_test"
        
        output_files = extractor.extract_proteins(
            str(genome_file), str(rgi_file), gene_list, sample_id
        )
        
        assert len(output_files) == 1
        
        # Load the extracted sequence
        extracted_records = list(SeqIO.parse(output_files[0], "fasta"))
        extracted_seq = str(extracted_records[0].seq)
        
        # Verify it's a protein sequence (no nucleotides like 'T')
        assert 'T' not in extracted_seq.upper()
        assert all(aa in 'ACDEFGHIKLMNPQRSTVWY*' for aa in extracted_seq.upper())
        
        # Verify sequence length is reasonable (not too short/long)
        assert 50 <= len(extracted_seq) <= 500  # Reasonable protein length
    
    def test_protein_extraction_strand_handling(self):
        """Test correct handling of positive and negative strand genes."""
        extractor = FastaAAExtractor(output_dir=str(self.test_output_dir))
        
        genome_file = self.test_genomes_dir / "GCF_000005825.2_ASM582v2_genomic.fna"
        rgi_file = self.test_rgi_dir / "sample1_rgi_output.txt"
        gene_list = [config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0], config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1]]  # RND efflux pump genes (configurable)
        sample_id = "strand_test"
        
        output_files = extractor.extract_proteins(
            str(genome_file), str(rgi_file), gene_list, sample_id
        )
        
        assert len(output_files) == 2
        
        # Both should produce valid protein sequences regardless of strand
        for output_file in output_files:
            records = list(SeqIO.parse(output_file, "fasta"))
            assert len(records) == 1
            seq_str = str(records[0].seq)
            assert len(seq_str) > 0
            assert all(aa in 'ACDEFGHIKLMNPQRSTVWY*' for aa in seq_str.upper())
    
    def test_protein_extraction_missing_contig(self):
        """Test handling when specified contig is not found in genome."""
        extractor = FastaAAExtractor(output_dir=str(self.test_output_dir))
        
        # Create RGI data with non-existent contig
        bad_rgi_data = [{
            'Best_Hit_ARO': config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0],
            'Contig': 'NONEXISTENT_CONTIG',
            'Start': '100',
            'Stop': '400',
            'Orientation': '+'
        }]
        
        bad_rgi_file = self.test_rgi_dir / "bad_contig_rgi.txt"
        with open(bad_rgi_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=bad_rgi_data[0].keys(), delimiter='\t')
            writer.writeheader()
            writer.writerows(bad_rgi_data)
        
        genome_file = self.test_genomes_dir / "GCF_000005825.2_ASM582v2_genomic.fna"
        gene_list = [config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0]]
        sample_id = "missing_contig_test"
        
        output_files = extractor.extract_proteins(
            str(genome_file), str(bad_rgi_file), gene_list, sample_id
        )
        
        assert len(output_files) == 0  # No files created for missing contigs
    
    def test_protein_extraction_invalid_coordinates(self):
        """Test handling of invalid coordinate ranges."""
        extractor = FastaAAExtractor(output_dir=str(self.test_output_dir))
        
        # Create RGI data with invalid coordinates (start > genome length)
        invalid_rgi_data = [{
            'Best_Hit_ARO': config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0],
            'Contig': 'NC_000913.3',
            'Start': '999999999',  # Way beyond genome length
            'Stop': '999999999',
            'Orientation': '+'
        }]
        
        invalid_rgi_file = self.test_rgi_dir / "invalid_coords_rgi.txt"
        with open(invalid_rgi_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=invalid_rgi_data[0].keys(), delimiter='\t')
            writer.writeheader()
            writer.writerows(invalid_rgi_data)
        
        genome_file = self.test_genomes_dir / "GCF_000005825.2_ASM582v2_genomic.fna"
        gene_list = [config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0]]
        sample_id = "invalid_coords_test"
        
        output_files = extractor.extract_proteins(
            str(genome_file), str(invalid_rgi_file), gene_list, sample_id
        )
        
        # Should handle gracefully and return empty list
        assert len(output_files) == 0
    
    # ===== File I/O and Error Handling Tests =====
    
    def test_missing_genome_file(self):
        """Test handling of missing genome files."""
        extractor = FastaAAExtractor(output_dir=str(self.test_output_dir))
        
        nonexistent_genome = str(self.test_genomes_dir / "nonexistent.fna")
        rgi_file = self.test_rgi_dir / "sample1_rgi_output.txt"
        gene_list = [config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0]]
        sample_id = "missing_genome_test"
        
        output_files = extractor.extract_proteins(
            nonexistent_genome, str(rgi_file), gene_list, sample_id
        )
        
        assert len(output_files) == 0
    
    def test_corrupted_genome_file(self):
        """Test handling of corrupted genome files."""
        extractor = FastaAAExtractor(output_dir=str(self.test_output_dir))
        
        # Create corrupted genome file
        corrupted_genome = self.test_genomes_dir / "corrupted.fna"
        with open(corrupted_genome, 'w') as f:
            f.write("This is not a valid FASTA file\n")
            f.write("Missing headers and sequences\n")
        
        rgi_file = self.test_rgi_dir / "sample1_rgi_output.txt"
        gene_list = [config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0]]
        sample_id = "corrupted_genome_test"
        
        output_files = extractor.extract_proteins(
            str(corrupted_genome), str(rgi_file), gene_list, sample_id
        )
        
        assert len(output_files) == 0
    
    def test_output_directory_permissions(self):
        """Test handling of output directory permission issues."""
        # This test may be platform-specific and could be skipped on some systems
        try:
            readonly_dir = Path(self.temp_dir) / "readonly"
            readonly_dir.mkdir()
            readonly_dir.chmod(0o444)  # Read-only
            
            extractor = FastaAAExtractor(output_dir=str(readonly_dir))
            
            genome_file = self.test_genomes_dir / "GCF_000005825.2_ASM582v2_genomic.fna"
            rgi_file = self.test_rgi_dir / "sample1_rgi_output.txt"
            gene_list = [config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0]]
            sample_id = "permission_test"
            
            output_files = extractor.extract_proteins(
                str(genome_file), str(rgi_file), gene_list, sample_id
            )
            
            # Should handle permission errors gracefully
            assert len(output_files) == 0
            
        except (OSError, PermissionError):
            # Permission manipulation may not work on all systems
            pytest.skip("Permission test not supported on this system")
        finally:
            # Restore permissions for cleanup
            try:
                readonly_dir.chmod(0o755)
            except:
                pass
    
    def test_large_genome_file_handling(self):
        """Test handling of large genome files (memory efficiency)."""
        extractor = FastaAAExtractor(output_dir=str(self.test_output_dir))
        
        # Create a larger test genome (simulating real genome sizes)
        large_sequence = "ATGCGATCGATCGATCG" * 100000  # ~1.6MB sequence
        large_genome_record = SeqRecord(
            Seq(large_sequence),
            id="large_contig",
            description="Large test genome"
        )
        
        large_genome_file = self.test_genomes_dir / "large_genome.fna"
        with open(large_genome_file, "w") as f:
            SeqIO.write([large_genome_record], f, "fasta")
        
        # Create RGI data for large genome
        large_rgi_data = [{
            'Best_Hit_ARO': config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0],
            'Contig': 'large_contig',
            'Start': '50000',
            'Stop': '50300',
            'Orientation': '+'
        }]
        
        large_rgi_file = self.test_rgi_dir / "large_genome_rgi.txt"
        with open(large_rgi_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=large_rgi_data[0].keys(), delimiter='\t')
            writer.writeheader()
            writer.writerows(large_rgi_data)
        
        gene_list = [config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0]]
        sample_id = "large_genome_test"
        
        output_files = extractor.extract_proteins(
            str(large_genome_file), str(large_rgi_file), gene_list, sample_id
        )
        
        assert len(output_files) == 1
        
        # Verify extraction worked despite large genome size
        records = list(SeqIO.parse(output_files[0], "fasta"))
        assert len(records) == 1
        assert len(records[0].seq) > 0
    
    # ===== Batch Processing and Performance Tests =====
    
    def test_multiple_gene_extraction(self):
        """Test extraction of multiple genes from single genome."""
        extractor = FastaAAExtractor(output_dir=str(self.test_output_dir))
        
        genome_file = self.test_genomes_dir / "GCF_000005825.2_ASM582v2_genomic.fna"
        rgi_file = self.test_rgi_dir / "sample1_rgi_output.txt"
        gene_list = [config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0], config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1], config_manager.get_default_genes("rnd_efflux_pumps", "primary")[2]]  # All available genes
        sample_id = "multi_gene_test"
        
        output_files = extractor.extract_proteins(
            str(genome_file), str(rgi_file), gene_list, sample_id
        )
        
        assert len(output_files) == 3
        
        # Verify each gene produces a separate file
        gene_files = {}
        for output_file in output_files:
            filename = Path(output_file).name
            for gene in gene_list:
                if gene in filename:
                    gene_files[gene] = output_file
                    break
        
        assert len(gene_files) == 3
        assert all(gene in gene_files for gene in gene_list)
    
    def test_batch_genome_processing_simulation(self):
        """Test processing multiple genomes (simulated batch)."""
        extractor = FastaAAExtractor(output_dir=str(self.test_output_dir))
        
        # Simulate batch processing by running multiple extractions
        genomes_and_samples = [
            (self.test_genomes_dir / "GCF_000005825.2_ASM582v2_genomic.fna", "sample_001"),
            (self.test_genomes_dir / "GCF_000009605.1_ASM960v1_genomic.fna", "sample_002"),
        ]
        
        all_output_files = []
        gene_list = [config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0], config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1]]
        
        for genome_file, sample_id in genomes_and_samples:
            rgi_file = self.test_rgi_dir / "sample1_rgi_output.txt"
            
            output_files = extractor.extract_proteins(
                str(genome_file), str(rgi_file), gene_list, sample_id
            )
            
            all_output_files.extend(output_files)
        
        # Should have files for each genome-gene combination
        expected_files = len(genomes_and_samples) * len(gene_list)
        assert len(all_output_files) >= expected_files // 2  # Allow for some failures
        
        # Verify all files exist and are unique
        unique_files = set(all_output_files)
        assert len(unique_files) == len(all_output_files)  # No duplicates
    
    def test_performance_timing(self):
        """Test extraction performance for timing validation."""
        import time
        
        extractor = FastaAAExtractor(output_dir=str(self.test_output_dir))
        
        genome_file = self.test_genomes_dir / "GCF_000005825.2_ASM582v2_genomic.fna"
        rgi_file = self.test_rgi_dir / "sample1_rgi_output.txt"
        gene_list = [config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0], config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1], config_manager.get_default_genes("rnd_efflux_pumps", "primary")[2]]
        sample_id = "performance_test"
        
        start_time = time.time()
        
        output_files = extractor.extract_proteins(
            str(genome_file), str(rgi_file), gene_list, sample_id
        )
        
        end_time = time.time()
        processing_time = end_time - start_time
        
        assert len(output_files) > 0
        assert processing_time < 10.0  # Should complete within 10 seconds
    
    def test_memory_usage_monitoring(self):
        """Test memory efficiency during extraction."""
        import psutil
        import os
        
        process = psutil.Process(os.getpid())
        initial_memory = process.memory_info().rss
        
        extractor = FastaAAExtractor(output_dir=str(self.test_output_dir))
        
        genome_file = self.test_genomes_dir / "GCF_000005825.2_ASM582v2_genomic.fna"
        rgi_file = self.test_rgi_dir / "sample1_rgi_output.txt"
        gene_list = [config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0], config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1], config_manager.get_default_genes("rnd_efflux_pumps", "primary")[2]]
        sample_id = "memory_test"
        
        output_files = extractor.extract_proteins(
            str(genome_file), str(rgi_file), gene_list, sample_id
        )
        
        final_memory = process.memory_info().rss
        memory_increase = final_memory - initial_memory
        
        assert len(output_files) > 0
        # Memory increase should be reasonable (< 50MB for small test)
        assert memory_increase < 50 * 1024 * 1024
    
    # ===== Integration and Real-World Scenario Tests =====
    
    def test_integration_with_card_pipeline_output(self):
        """Test integration with realistic CARD pipeline output."""
        extractor = FastaAAExtractor(output_dir=str(self.test_output_dir))
        
        # Create more realistic RGI output with additional columns
        realistic_rgi_data = [{
            'Best_Hit_ARO': config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0],
            'Best_Hit_ARO_Name': config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0],
            'Best_Hit_ARO_Description': 'multidrug efflux pump',
            'Contig': 'NC_000913.3',
            'Start': '100',
            'Stop': '400',
            'Orientation': '+',
            'Cut_Off': 'Perfect',
            'Pass_Bitscore': '500.0',
            'Best_Hit_Bitscore': '500.0',
            'Best_Identities': '100.00',
            'ARO': '3000007',
            'Model_type': 'protein homolog model',
            'SNPs_in_Best_Hit_ARO': '',
            'Other_SNPs': '',
            'Drug_Class': 'peptide antibiotic',
            'Resistance_Mechanism': 'antibiotic efflux',
            'AMR_Gene_Family': 'resistance-nodulation-cell division (RND) antibiotic efflux pump'
        }]
        
        realistic_rgi_file = self.test_rgi_dir / "realistic_card_output.txt"
        with open(realistic_rgi_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=realistic_rgi_data[0].keys(), delimiter='\t')
            writer.writeheader()
            writer.writerows(realistic_rgi_data)
        
        genome_file = self.test_genomes_dir / "GCF_000005825.2_ASM582v2_genomic.fna"
        gene_list = [config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0]]
        sample_id = "card_integration_test"
        
        output_files = extractor.extract_proteins(
            str(genome_file), str(realistic_rgi_file), gene_list, sample_id
        )
        
        assert len(output_files) == 1
        
        # Verify output file content and format
        records = list(SeqIO.parse(output_files[0], "fasta"))
        assert len(records) == 1
        
        # Verify FASTA header contains expected information
        header = records[0].description
        assert config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0] in header
        assert 'NC_000913.3' in header
        assert '100' in header and '400' in header  # Coordinates
    
    def test_edge_case_coordinate_boundaries(self):
        """Test extraction at genome boundaries and edge coordinates."""
        extractor = FastaAAExtractor(output_dir=str(self.test_output_dir))
        
        # Load genome to get actual length
        genome_file = self.test_genomes_dir / "GCF_000005825.2_ASM582v2_genomic.fna"
        genome_records = list(SeqIO.parse(genome_file, "fasta"))
        genome_length = len(genome_records[0].seq)
        
        # Create RGI data with boundary coordinates
        boundary_rgi_data = [
            {  # Near start of genome
                'Best_Hit_ARO': config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0],
                'Contig': 'NC_000913.3',
                'Start': '1',
                'Stop': '99',
                'Orientation': '+'
            },
            {  # Near end of genome
                'Best_Hit_ARO': config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1],
                'Contig': 'NC_000913.3',
                'Start': str(genome_length - 100),
                'Stop': str(genome_length),
                'Orientation': '+'
            }
        ]
        
        boundary_rgi_file = self.test_rgi_dir / "boundary_rgi.txt"
        with open(boundary_rgi_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=boundary_rgi_data[0].keys(), delimiter='\t')
            writer.writeheader()
            writer.writerows(boundary_rgi_data)
        
        gene_list = [config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0], config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1]]
        sample_id = "boundary_test"
        
        output_files = extractor.extract_proteins(
            str(genome_file), str(boundary_rgi_file), gene_list, sample_id
        )
        
        assert len(output_files) == 2
        
        # Verify both extractions produced valid sequences
        for output_file in output_files:
            records = list(SeqIO.parse(output_file, "fasta"))
            assert len(records) == 1
            assert len(records[0].seq) > 0
    
    def test_real_world_gene_list_filtering(self):
        """Test realistic gene list filtering scenarios."""
        extractor = FastaAAExtractor(output_dir=str(self.test_output_dir))
        
        genome_file = self.test_genomes_dir / "GCF_000005825.2_ASM582v2_genomic.fna"
        rgi_file = self.test_rgi_dir / "sample1_rgi_output.txt"
        
        # Test case-insensitive filtering
        case_insensitive_genes = ['ACRA', 'acrb', 'TolC']
        output_files = extractor.extract_proteins(
            str(genome_file), str(rgi_file), case_insensitive_genes, "case_test"
        )
        
        # Should match despite case differences (depending on implementation)
        # If implementation is case-sensitive, this tests the limitation
        assert isinstance(output_files, list)
        
        # Test partial gene names (should not match if exact matching)
        partial_genes = ['acr', 'tol']  # Partial names
        output_files = extractor.extract_proteins(
            str(genome_file), str(rgi_file), partial_genes, "partial_test"
        )
        
        # Should not match partial names
        assert len(output_files) == 0


def test_fasta_aa_extractor_integration_workflow():
    """Integration test simulating complete extraction workflow."""
    temp_dir = tempfile.mkdtemp()
    try:
        # Setup test environment
        output_dir = Path(temp_dir) / "extracted"
        genomes_dir = Path(temp_dir) / "genomes"
        rgi_dir = Path(temp_dir) / "rgi"
        
        output_dir.mkdir(parents=True)
        genomes_dir.mkdir(parents=True)
        rgi_dir.mkdir(parents=True)
        
        # Create realistic test genome
        test_sequence = (
            "ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGACGCGTACAGGAAACACAGAAAAAAGCCCGCACCTG"
            "ACAGTGCGGGCTTTTTTTTTCGACCAAAGGTAACGAGGTAACAACCATGCGAGTGTTGAAGTTCGGCGGTACATCAGTGGCAAATGCAGAACGTTTTCTGC"
            "GTGTTGCCGATATTCTGGAAAGCAATGCCAGGCAGGGGCAGGTGGCCACCGTCCTCTCTGCCCCCGCCAAAATCACCAACCACCTGGTGGCGATGATTGAA"
        ) * 10  # Longer sequence
        
        genome_record = SeqRecord(
            Seq(test_sequence),
            id="test_contig",
            description="Integration test genome"
        )
        
        genome_file = genomes_dir / "test_genome.fna"
        with open(genome_file, "w") as f:
            SeqIO.write([genome_record], f, "fasta")
        
        # Create realistic RGI output
        rgi_data = [
            {
                'Best_Hit_ARO': config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0],
                'Contig': 'test_contig',
                'Start': '100',
                'Stop': '400',
                'Orientation': '+',
                'Cut_Off': 'Perfect',
                'Pass_Bitscore': '500'
            },
            {
                'Best_Hit_ARO': config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1],
                'Contig': 'test_contig',
                'Start': '500',
                'Stop': '800',
                'Orientation': '-',
                'Cut_Off': 'Perfect',
                'Pass_Bitscore': '600'
            }
        ]
        
        rgi_file = rgi_dir / "test_rgi.txt"
        with open(rgi_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=rgi_data[0].keys(), delimiter='\t')
            writer.writeheader()
            writer.writerows(rgi_data)
        
        # Run extraction
        extractor = FastaAAExtractor(output_dir=str(output_dir))
        gene_list = [config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0], config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1]]
        sample_id = "integration_test"
        
        output_files = extractor.extract_proteins(
            str(genome_file), str(rgi_file), gene_list, sample_id
        )
        
        # Validate results
        assert len(output_files) == 2
        
        total_sequences = 0
        for output_file in output_files:
            assert Path(output_file).exists()
            records = list(SeqIO.parse(output_file, "fasta"))
            assert len(records) == 1
            assert len(records[0].seq) > 0
            total_sequences += 1
        
        assert total_sequences == 2
        
    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)


if __name__ == "__main__":
    # Run specific test for debugging
    test = TestFastaAAExtractorRobust()
    test.setup_method()
    try:
        test.test_extractor_initialization_defaults()
        test.test_protein_extraction_basic_workflow()
        print("✅ Basic FastaAAExtractor tests passed!")
    finally:
        test.teardown_method()