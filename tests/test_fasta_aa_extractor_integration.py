import sys; import os; sys.path.append(os.path.join(os.path.dirname(__file__), '\\\src')); from configuration_manager import config_manager
#!/usr/bin/env python3
"""
Comprehensive Test Suite for FastaAAExtractor Integration
Tests all functionality of the FastaAAExtractor pipeline component

This test suite validates:
1. CARD coordinate loading from multiple formats
2. Protein extraction workflows (internal and external)
3. WildTypeAligner preparation
4. Error handling and edge cases
5. Output file generation and quality

Author: GenomeAMRAnalyzer Pipeline
Version: 1.0 - Production Ready
"""

import os
import sys
import tempfile
import shutil
import json
import csv
from pathlib import Path
from unittest.mock import patch, MagicMock
import traceback

# Import the module we're testing
try:
    from fasta_aa_extractor_integration import (
        FastaAAExtractorIntegrator,
        GeneCoordinate,
        ExtractionResult
    )
    INTEGRATION_MODULE_AVAILABLE = True
except ImportError:
    INTEGRATION_MODULE_AVAILABLE = False
    print("Warning: FastaAAExtractor integration module not available")


class TestFastaAAExtractorIntegration:
    """Comprehensive test suite for FastaAAExtractor integration"""
    
    def __init__(self):
        """Initialize test environment"""
        self.test_dir = None
        self.integrator = None
        self.test_results = []
        
    def setup_test_environment(self):
        """Create temporary test environment"""
        self.test_dir = Path(tempfile.mkdtemp(prefix="fasta_extractor_test_"))
        
        # Create directory structure
        self.genomes_dir = self.test_dir / "genomes"
        self.genomes_dir.mkdir()
        
        self.references_dir = self.test_dir / "references"
        self.references_dir.mkdir()
        
        self.output_dir = self.test_dir / "output"
        self.output_dir.mkdir()
        
        # Create mock genome files
        self._create_mock_genomes()
        
        # Create mock reference files
        self._create_mock_references()
        
        # Create mock coordinate files
        self._create_mock_coordinates()
        
        print(f"Test environment created at: {self.test_dir}")
    
    def cleanup_test_environment(self):
        """Clean up test environment"""
        if self.test_dir and self.test_dir.exists():
            shutil.rmtree(self.test_dir)
            print("Test environment cleaned up")
    
    def _create_mock_genomes(self):
        """Create mock genome FASTA files"""
        # Simple nucleotide sequences for testing
        mock_genomes = {
            "GCF_000005825.2": {
                "sequence": "ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGACGCGTACAGGAAACACAGAAAAAAGCCCGCACCTGACAGTGCGGGCTTTTTTTTTCGACCAAAGGTAACGAGGTAACAACCATGCGAGTGTTGAAGTTCGGCGGTACATCAGTGGCAAATGCAGAACGTTTTCTGCGTGTTGCCGATATTCTGGAAAGCAATGCCAGGCAGGGGCAGGTGGCCACCGTCCTCTCTGCCCCCGCCAAAATCACCAACCACCTGGTGGCGATGATTGAAAAAACCATTAGCGGCCAGGATGCTTTACCCAATATCAGCGATGCCGAACGTATTTTTGCCGAACTTTTGACGGGACTCGCCGCCGCCCAGCCGGGGTTCCCGCTGGCGCAATTGAAAACTTTCGTCGATCAGGAATTTGCCCAAATAAAACATGTCCTGCATGGCATTAGTTTGTTGGGGCAGTGCCCGGATAGCATCAACGCTGCGCTGATTTGCCGTGGCGAGAAAATGTCGATCGCCATTATGGCCGGCGTATTAGAAGCGCGCGGTCACAACGTTACTGTTATCGATCCGGTCGTTCTGAATGGCGGTTTCCGGGGCTGGCCGCATCATGGCGGCTTGGGTTGGTTTGGGTTTGGGTTTGGGTTTGGGTTTGGGTTTGGGTTTGGGTTTGGGTTTGGGTTTGGGTTTGGGTTTGGGTTTGGGTTTG",
                "genes": {
                    "mdtF": {"start": 100, "end": 450, "strand": "+"},
                    config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0]: {"start": 500, "end": 850, "strand": "+"},
                    config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1]: {"start": 900, "end": 1400, "strand": "+"}
                }
            },
            "GCF_000006945.2": {
                "sequence": "ATGAAACGCATTAGCACCACCATTACCACCACCATCACCATTACCACAGGTAACGGTGCGGGCTGACGCGTACAGGAAACACAGAAAAAAGCCCGCACCTGACAGTGCGGGCTTTTTTTTTCGACCAAAGGTAACGAGGTAACAACCATGCGAGTGTTGAAGTTCGGCGGTACATCAGTGGCAAATGCAGAACGTTTTCTGCGTGTTGCCGATATTCTGGAAAGCAATGCCAGGCAGGGGCAGGTGGCCACCGTCCTCTCTGCCCCCGCCAAAATCACCAACCACCTGGTGGCGATGATTGAAAAAACCATTAGCGGCCAGGATGCTTTACCCAATATCAGCGATGCCGAACGTATTTTTGCCGAACTTTTGACGGGACTCGCCGCCGCCCAGCCGGGGTTCCCGCTGGCGCAATTGAAAACTTTCGTCGATCAGGAATTTGCCCAAATAAAACATGTCCTGCATGGCATTAGTTTGTTGGGGCAGTGCCCGGATAGCATCAACGCTGCGCTGATTTGCCGTGGCGAGAAAATGTCGATCGCCATTATGGCCGGCGTATTAGAAGCGCGCGGTCACAACGTTACTGTTATCGATCCGGTCGTTCTGAATGGCGGTTTCCGGGGCTGGCCGCATCATGGCGGCTTGGGTTGGTTTGGGTTTGGGTTTGGGTTTGGGTTTGGGTTTGGGTTTGGGTTTGGGTTTGGGTTTGGGTTTGGGTTTGGGTTTGGGTTTG",
                "genes": {
                    "mdtF": {"start": 120, "end": 470, "strand": "+"},
                    config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0]: {"start": 520, "end": 870, "strand": "+"},
                    config_manager.get_default_genes("rnd_efflux_pumps", "primary")[2]: {"start": 920, "end": 1320, "strand": "+"}
                }
            }
        }
        
        for genome_id, data in mock_genomes.items():
            genome_file = self.genomes_dir / f"{genome_id}_genomic.fna"
            with open(genome_file, 'w') as f:
                f.write(f">{genome_id}_chromosome\n")
                # Write sequence in lines of 80 characters
                seq = data["sequence"]
                for i in range(0, len(seq), 80):
                    f.write(f"{seq[i:i+80]}\n")
    
    def _create_mock_references(self):
        """Create mock reference sequence files"""
        reference_sequences = {
            "mdtF": "MAPLKTLLSLLQFRSLLLLGLLLLPSVREEEERERERERMLRTNSLTLRKRTKRTKRTN",
            config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0]: "MARAPLKTLLSLLQFRSLLLLGLLLLPSVREEEERERERMLRTNSLTLRKRTKRTKRTN",
            config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1]: "MPLKTLLSLLQFRSLLLLGLLLLPSVREEEERERERMLRTNSLTLRKRTKRTKRTNERE",
            config_manager.get_default_genes("rnd_efflux_pumps", "primary")[2]: "MKTLLSLLQFRSLLLLGLLLLPSVREEEERERERMLRTNSLTLRKRTKRTKRTNERERE"
        }
        
        for gene, sequence in reference_sequences.items():
            ref_file = self.references_dir / f"{gene}_reference.faa"
            with open(ref_file, 'w') as f:
                f.write(f">{gene}_MG1655_reference\n")
                f.write(f"{sequence}\n")
    
    def _create_mock_coordinates(self):
        """Create mock coordinate files in different formats"""
        # CSV format
        csv_file = self.test_dir / "coordinates.csv"
        with open(csv_file, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['genome_id', 'gene_name', 'start', 'end', 'strand', 'contig'])
            writer.writerow(['GCF_000005825.2', 'mdtF', '100', '450', '+', 'GCF_000005825.2_chromosome'])
            writer.writerow(['GCF_000005825.2', config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0], '500', '850', '+', 'GCF_000005825.2_chromosome'])
            writer.writerow(['GCF_000005825.2', config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1], '900', '1400', '+', 'GCF_000005825.2_chromosome'])
            writer.writerow(['GCF_000006945.2', 'mdtF', '120', '470', '+', 'GCF_000006945.2_chromosome'])
            writer.writerow(['GCF_000006945.2', config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0], '520', '870', '+', 'GCF_000006945.2_chromosome'])
            writer.writerow(['GCF_000006945.2', config_manager.get_default_genes("rnd_efflux_pumps", "primary")[2], '920', '1320', '+', 'GCF_000006945.2_chromosome'])
        
        # TSV format
        tsv_file = self.test_dir / "coordinates.tsv"
        with open(tsv_file, 'w') as f:
            f.write("genome_id\tgene_name\tstart\tend\tstrand\tcontig\n")
            f.write("GCF_000005825.2\tmdtF\t100\t450\t+\tGCF_000005825.2_chromosome\n")
            f.write("GCF_000005825.2\tacrA\t500\t850\t+\tGCF_000005825.2_chromosome\n")
            f.write("GCF_000005825.2\tacrB\t900\t1400\t+\tGCF_000005825.2_chromosome\n")
            f.write("GCF_000006945.2\tmdtF\t120\t470\t+\tGCF_000006945.2_chromosome\n")
            f.write("GCF_000006945.2\tacrA\t520\t870\t+\tGCF_000006945.2_chromosome\n")
            f.write("GCF_000006945.2\ttolC\t920\t1320\t+\tGCF_000006945.2_chromosome\n")
        
        # JSON format
        json_file = self.test_dir / "coordinates.json"
        json_data = {
            "GCF_000005825.2": {
                "mdtF": {"start": 100, "end": 450, "strand": "+", "contig": "GCF_000005825.2_chromosome"},
                config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0]: {"start": 500, "end": 850, "strand": "+", "contig": "GCF_000005825.2_chromosome"},
                config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1]: {"start": 900, "end": 1400, "strand": "+", "contig": "GCF_000005825.2_chromosome"}
            },
            "GCF_000006945.2": {
                "mdtF": {"start": 120, "end": 470, "strand": "+", "contig": "GCF_000006945.2_chromosome"},
                config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0]: {"start": 520, "end": 870, "strand": "+", "contig": "GCF_000006945.2_chromosome"},
                config_manager.get_default_genes("rnd_efflux_pumps", "primary")[2]: {"start": 920, "end": 1320, "strand": "+", "contig": "GCF_000006945.2_chromosome"}
            }
        }
        
        with open(json_file, 'w') as f:
            json.dump(json_data, f, indent=2)
    
    def run_test(self, test_name, test_function):
        """Run individual test and capture results"""
        print(f"\n{'='*20} Running {test_name} {'='*20}")
        
        try:
            result = test_function()
            status = "PASS" if result else "FAIL"
            self.test_results.append({
                'test': test_name,
                'status': status,
                'result': result,
                'error': None
            })
            print(f"✓ {test_name}: {status}")
            return result
            
        except Exception as e:
            self.test_results.append({
                'test': test_name,
                'status': 'ERROR',
                'result': False,
                'error': str(e)
            })
            print(f"✗ {test_name}: ERROR - {e}")
            print(f"Traceback: {traceback.format_exc()}")
            return False
    
    def test_gene_coordinate_validation(self):
        """Test GeneCoordinate class validation"""
        print("Testing GeneCoordinate validation...")
        
        # Valid coordinate
        try:
            coord = GeneCoordinate(
                gene_name=config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0],
                genome_id="test_genome",
                start=100,
                end=400,
                strand="+"
            )
            assert coord.length == 301
            assert coord.is_forward_strand == True
            print("  ✓ Valid coordinate creation")
        except Exception as e:
            print(f"  ✗ Valid coordinate creation failed: {e}")
            return False
        
        # Invalid coordinates
        try:
            GeneCoordinate(
                gene_name="invalid",
                genome_id="test",
                start=400,
                end=100,  # Invalid: end < start
                strand="+"
            )
            print("  ✗ Should have failed with invalid coordinates")
            return False
        except ValueError:
            print("  ✓ Correctly rejected invalid coordinates")
        
        # Invalid strand
        try:
            GeneCoordinate(
                gene_name="invalid",
                genome_id="test",
                start=100,
                end=400,
                strand="invalid"  # Invalid strand
            )
            print("  ✗ Should have failed with invalid strand")
            return False
        except ValueError:
            print("  ✓ Correctly rejected invalid strand")
        
        return True
    
    def test_integrator_initialization(self):
        """Test FastaAAExtractorIntegrator initialization"""
        print("Testing integrator initialization...")
        
        try:
            integrator = FastaAAExtractorIntegrator(
                output_dir=str(self.output_dir / "init_test"),
                log_level="INFO"
            )
            
            # Check if output directory was created
            if not (self.output_dir / "init_test").exists():
                print("  ✗ Output directory not created")
                return False
            
            # Check if logging is set up
            if not hasattr(integrator, 'logger'):
                print("  ✗ Logger not initialized")
                return False
            
            print("  ✓ Integrator initialized successfully")
            return True
            
        except Exception as e:
            print(f"  ✗ Initialization failed: {e}")
            return False
    
    def test_coordinate_loading_csv(self):
        """Test loading coordinates from CSV file"""
        print("Testing CSV coordinate loading...")
        
        try:
            integrator = FastaAAExtractorIntegrator(
                output_dir=str(self.output_dir / "csv_test")
            )
            
            # Load CSV coordinates
            csv_file = self.test_dir / "coordinates.csv"
            integrator.load_card_coordinates(csv_file, file_format="csv")
            
            # Check if coordinates were loaded
            if len(integrator.gene_coordinates) != 2:
                print(f"  ✗ Expected 2 genomes, got {len(integrator.gene_coordinates)}")
                return False
            
            # Check specific genome
            gcf_825_coords = integrator.gene_coordinates.get("GCF_000005825.2", [])
            if len(gcf_825_coords) != 3:
                print(f"  ✗ Expected 3 genes for GCF_000005825.2, got {len(gcf_825_coords)}")
                return False
            
            gene_names = [coord.gene_name for coord in gcf_825_coords]
            expected_genes = ['mdtF', config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0], config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1]]
            if set(gene_names) != set(expected_genes):
                print(f"  ✗ Expected genes {expected_genes}, got {gene_names}")
                return False
            
            print("  ✓ CSV coordinates loaded correctly")
            return True
            
        except Exception as e:
            print(f"  ✗ CSV loading failed: {e}")
            return False
    
    def test_coordinate_loading_json(self):
        """Test loading coordinates from JSON file"""
        print("Testing JSON coordinate loading...")
        
        try:
            integrator = FastaAAExtractorIntegrator(
                output_dir=str(self.output_dir / "json_test")
            )
            
            # Load JSON coordinates
            json_file = self.test_dir / "coordinates.json"
            integrator.load_card_coordinates(json_file, file_format="json")
            
            # Check if coordinates were loaded
            if len(integrator.gene_coordinates) != 2:
                print(f"  ✗ Expected 2 genomes, got {len(integrator.gene_coordinates)}")
                return False
            
            # Check specific coordinate
            gcf_945_coords = integrator.gene_coordinates.get("GCF_000006945.2", [])
            tolc_coord = None
            for coord in gcf_945_coords:
                if coord.gene_name == config_manager.get_default_genes("rnd_efflux_pumps", "primary")[2]:
                    tolc_coord = coord
                    break
            
            if not tolc_coord:
                print("  ✗ tolC coordinate not found")
                return False
            
            if tolc_coord.start != 920 or tolc_coord.end != 1320:
                print(f"  ✗ tolC coordinates incorrect: {tolc_coord.start}-{tolc_coord.end}")
                return False
            
            print("  ✓ JSON coordinates loaded correctly")
            return True
            
        except Exception as e:
            print(f"  ✗ JSON loading failed: {e}")
            return False
    
    def test_genome_file_discovery(self):
        """Test finding genome files"""
        print("Testing genome file discovery...")
        
        try:
            integrator = FastaAAExtractorIntegrator(
                output_dir=str(self.output_dir / "genome_test")
            )
            
            # Load coordinates first
            csv_file = self.test_dir / "coordinates.csv"
            integrator.load_card_coordinates(csv_file)
            
            # Test genome file discovery
            genome_files = integrator._find_genome_files(self.genomes_dir)
            
            if len(genome_files) != 2:
                print(f"  ✗ Expected 2 genome files, found {len(genome_files)}")
                return False
            
            expected_genomes = {"GCF_000005825.2", "GCF_000006945.2"}
            found_genomes = set(genome_files.keys())
            
            if found_genomes != expected_genomes:
                print(f"  ✗ Expected genomes {expected_genomes}, found {found_genomes}")
                return False
            
            # Check file paths
            for genome_id, file_path in genome_files.items():
                if not file_path.exists():
                    print(f"  ✗ Genome file not found: {file_path}")
                    return False
            
            print("  ✓ Genome files discovered correctly")
            return True
            
        except Exception as e:
            print(f"  ✗ Genome discovery failed: {e}")
            return False
    
    def test_gene_filtering(self):
        """Test filtering coordinates by gene list"""
        print("Testing gene filtering...")
        
        try:
            integrator = FastaAAExtractorIntegrator(
                output_dir=str(self.output_dir / "filter_test")
            )
            
            # Load all coordinates
            csv_file = self.test_dir / "coordinates.csv"
            integrator.load_card_coordinates(csv_file)
            
            # Count initial coordinates
            initial_count = sum(len(coords) for coords in integrator.gene_coordinates.values())
            
            # Filter to specific genes
            integrator._filter_coordinates_by_genes(['mdtF', config_manager.get_default_genes("rnd_efflux_pumps", "primary")[2]])
            
            # Check filtering
            filtered_count = sum(len(coords) for coords in integrator.gene_coordinates.values())
            
            if filtered_count >= initial_count:
                print(f"  ✗ Filtering didn't reduce coordinates: {initial_count} -> {filtered_count}")
                return False
            
            # Check that only desired genes remain
            remaining_genes = set()
            for coords in integrator.gene_coordinates.values():
                remaining_genes.update(coord.gene_name for coord in coords)
            
            expected_genes = {'mdtF', config_manager.get_default_genes("rnd_efflux_pumps", "primary")[2]}
            if remaining_genes != expected_genes:
                print(f"  ✗ Expected genes {expected_genes}, got {remaining_genes}")
                return False
            
            print("  ✓ Gene filtering worked correctly")
            return True
            
        except Exception as e:
            print(f"  ✗ Gene filtering failed: {e}")
            return False
    
    def test_internal_extraction_simulation(self):
        """Test internal protein extraction (simulated without BioPython)"""
        print("Testing internal extraction simulation...")
        
        try:
            integrator = FastaAAExtractorIntegrator(
                output_dir=str(self.output_dir / "internal_test")
            )
            
            # Load coordinates
            csv_file = self.test_dir / "coordinates.csv"
            integrator.load_card_coordinates(csv_file)
            
            # Create mock extraction results (simulating what would happen)
            mock_results = [
                ExtractionResult(
                    gene_name="mdtF",
                    genome_id="GCF_000005825.2",
                    sequence="MKRLAPPITTTITITTQVTVRADAYRKHRKPPPDSAGFFFFFFFFDQRVTVTACEPTFGKLSSANAMQNFLRCCDILESNARGRGAVATCLCCTKITLGLAIDEKTISGRMLYPISDPERISFELDRPAPAGFPLAQEKRHYLIMFPKKCSHGILFL",
                    coordinates=GeneCoordinate("mdtF", "GCF_000005825.2", 100, 450, "+"),
                    extraction_method="simulated_internal"
                ),
                ExtractionResult(
                    gene_name=config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0],
                    genome_id="GCF_000005825.2",
                    sequence="MKRLAPPITTTITITTQVTVRADAYRKHRKPPPDSAGFFFFFFFFDQRVTVTACEPTFGKLSSANAMQNFLRCCDILESNARGRGAVATCLCCTKITLGLAIDEKTISGRMLYPISDPERISFELDRPAPAGFPLAQEKRHYLIMFPKKCSHGILFL",
                    coordinates=GeneCoordinate(config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0], "GCF_000005825.2", 500, 850, "+"),
                    extraction_method="simulated_internal"
                )
            ]
            
            integrator.extraction_results = mock_results
            integrator.stats['sequences_generated'] = len(mock_results)
            
            # Test output file generation
            output_files = integrator._generate_output_files()
            
            # Check output files
            required_outputs = ['individual_fasta', 'combined_fasta', 'metadata', 'summary']
            for output_type in required_outputs:
                if output_type not in output_files:
                    print(f"  ✗ Missing output: {output_type}")
                    return False
            
            # Check combined FASTA file
            combined_fasta = Path(output_files['combined_fasta'])
            if not combined_fasta.exists():
                print("  ✗ Combined FASTA file not created")
                return False
            
            with open(combined_fasta, 'r') as f:
                content = f.read()
                if ">GCF_000005825.2_mdtF" not in content:
                    print("  ✗ Expected sequence header not found")
                    return False
            
            print("  ✓ Internal extraction simulation successful")
            return True
            
        except Exception as e:
            print(f"  ✗ Internal extraction simulation failed: {e}")
            return False
    
    def test_reference_sequence_discovery(self):
        """Test finding reference sequences"""
        print("Testing reference sequence discovery...")
        
        try:
            integrator = FastaAAExtractorIntegrator(
                output_dir=str(self.output_dir / "ref_test")
            )
            
            # Test finding reference for each gene
            test_genes = ['mdtF', config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0], config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1], config_manager.get_default_genes("rnd_efflux_pumps", "primary")[2]]
            
            for gene in test_genes:
                ref_file = integrator._find_reference_sequence(self.references_dir, gene)
                
                if not ref_file:
                    print(f"  ✗ Reference not found for {gene}")
                    return False
                
                if not ref_file.exists():
                    print(f"  ✗ Reference file doesn't exist: {ref_file}")
                    return False
                
                expected_name = f"{gene}_reference.faa"
                if ref_file.name != expected_name:
                    print(f"  ✗ Unexpected reference file: {ref_file.name}")
                    return False
            
            # Test non-existent gene
            nonexistent_ref = integrator._find_reference_sequence(self.references_dir, "nonexistent")
            if nonexistent_ref is not None:
                print("  ✗ Should not have found reference for nonexistent gene")
                return False
            
            print("  ✓ Reference sequence discovery working correctly")
            return True
            
        except Exception as e:
            print(f"  ✗ Reference discovery failed: {e}")
            return False
    
    def test_wild_type_aligner_preparation(self):
        """Test preparation for WildTypeAligner"""
        print("Testing WildTypeAligner preparation...")
        
        try:
            integrator = FastaAAExtractorIntegrator(
                output_dir=str(self.output_dir / "aligner_test")
            )
            
            # Add mock extraction results
            mock_results = [
                ExtractionResult(
                    gene_name="mdtF",
                    genome_id="GCF_000005825.2",
                    sequence="MKRLAPPITTTITITTQVTVRADAYRKHRKPPPDSAGFFFFFFFFDQRVTVTACEPTFGKLSSANAMQNFLRCCDILESNARGRGAVATCLCCTKITLGLAIDEKTISGRMLYPISDPERISFELDRPAPAGFPLAQEKRHYLIMFPKKCSHGILFL",
                    coordinates=GeneCoordinate("mdtF", "GCF_000005825.2", 100, 450, "+")
                ),
                ExtractionResult(
                    gene_name="mdtF",
                    genome_id="GCF_000006945.2",
                    sequence="MKRLAPPITTTITITTQVTVRADAYRKHRKPPPDSAGFFFFFFFFDQRVTVTACEPTFGKLSSANAMQNFLRCCDILESNARGRGAVATCLCCTKITLGLAIDEKTISGRMLYPISDPERISFELDRPAPAGFPLAQEKRHYLIMFPKKCSHGILFL",
                    coordinates=GeneCoordinate("mdtF", "GCF_000006945.2", 120, 470, "+")
                ),
                ExtractionResult(
                    gene_name=config_manager.get_default_genes("rnd_efflux_pumps", "primary")[2],
                    genome_id="GCF_000006945.2",
                    sequence="MKRLAPPITTTITITTQVTVRADAYRKHRKPPPDSAGFFFFFFFFDQRVTVTACEPTFGKLSSANAMQNFLRCCDILESNARGRGAVATCLCCTKITLGLAIDEKTISGRMLYPISDPERISFELDRPAPAGFPLAQEKRHYLIMFPKKCSHGILFL",
                    coordinates=GeneCoordinate(config_manager.get_default_genes("rnd_efflux_pumps", "primary")[2], "GCF_000006945.2", 920, 1320, "+")
                )
            ]
            
            integrator.extraction_results = mock_results
            
            # Prepare for WildTypeAligner
            prepared_files = integrator.prepare_for_wild_type_aligner(
                reference_dir=self.references_dir,
                output_structure="by_gene"
            )
            
            # Check that files were prepared for each gene
            expected_genes = {'mdtF', config_manager.get_default_genes("rnd_efflux_pumps", "primary")[2]}
            prepared_genes = set(prepared_files.keys())
            
            if prepared_genes != expected_genes:
                print(f"  ✗ Expected genes {expected_genes}, prepared {prepared_genes}")
                return False
            
            # Check file structure for one gene
            mdtf_files = prepared_files.get('mdtF', {})
            required_keys = ['query_sequences', 'reference_sequence', 'output_directory']
            
            for key in required_keys:
                if key not in mdtf_files:
                    print(f"  ✗ Missing key {key} in mdtF preparation")
                    return False
                
                file_path = Path(mdtf_files[key])
                if not file_path.exists():
                    print(f"  ✗ File not created: {file_path}")
                    return False
            
            # Check query sequences file content
            query_file = Path(mdtf_files['query_sequences'])
            with open(query_file, 'r') as f:
                content = f.read()
                if content.count('>') != 2:  # Should have 2 mdtF sequences
                    print(f"  ✗ Expected 2 mdtF sequences, found {content.count('>')}")
                    return False
            
            print("  ✓ WildTypeAligner preparation successful")
            return True
            
        except Exception as e:
            print(f"  ✗ WildTypeAligner preparation failed: {e}")
            return False
    
    def test_error_handling(self):
        """Test error handling and edge cases"""
        print("Testing error handling...")
        
        try:
            # Test with non-existent coordinates file
            integrator = FastaAAExtractorIntegrator(
                output_dir=str(self.output_dir / "error_test")
            )
            
            try:
                integrator.load_card_coordinates("nonexistent_file.csv")
                print("  ✗ Should have failed with non-existent file")
                return False
            except FileNotFoundError:
                print("  ✓ Correctly handled non-existent coordinates file")
            
            # Test with empty coordinates
            empty_file = self.test_dir / "empty.csv"
            with open(empty_file, 'w') as f:
                f.write("")  # Empty file
            
            try:
                integrator.load_card_coordinates(empty_file)
                print("  ✗ Should have failed with empty file")
                return False
            except ValueError:
                print("  ✓ Correctly handled empty coordinates file")
            
            # Test with malformed coordinates
            malformed_file = self.test_dir / "malformed.csv"
            with open(malformed_file, 'w') as f:
                f.write("genome_id,gene_name,start,end,strand\n")
                f.write("test,gene1,invalid,400,+\n")  # Invalid start position
            
            integrator.load_card_coordinates(malformed_file)
            # Should handle malformed entries gracefully (warning, not error)
            print("  ✓ Gracefully handled malformed coordinates")
            
            return True
            
        except Exception as e:
            print(f"  ✗ Error handling test failed: {e}")
            return False
    
    def run_all_tests(self):
        """Run all tests and report results"""
        if not INTEGRATION_MODULE_AVAILABLE:
            print("ERROR: FastaAAExtractor integration module not available")
            return False
        
        print("Starting comprehensive test suite for FastaAAExtractor Integration")
        print(f"Test environment: {self.test_dir}")
        
        # Setup test environment
        self.setup_test_environment()
        
        # Define all tests
        tests = [
            ("Gene Coordinate Validation", self.test_gene_coordinate_validation),
            ("Integrator Initialization", self.test_integrator_initialization),
            ("CSV Coordinate Loading", self.test_coordinate_loading_csv),
            ("JSON Coordinate Loading", self.test_coordinate_loading_json),
            ("Genome File Discovery", self.test_genome_file_discovery),
            ("Gene Filtering", self.test_gene_filtering),
            ("Internal Extraction Simulation", self.test_internal_extraction_simulation),
            ("Reference Sequence Discovery", self.test_reference_sequence_discovery),
            ("WildTypeAligner Preparation", self.test_wild_type_aligner_preparation),
            ("Error Handling", self.test_error_handling)
        ]
        
        # Run all tests
        for test_name, test_function in tests:
            self.run_test(test_name, test_function)
        
        # Print summary
        self.print_test_summary()
        
        # Cleanup
        self.cleanup_test_environment()
        
        # Return overall success
        failed_tests = [r for r in self.test_results if r['status'] != 'PASS']
        return len(failed_tests) == 0
    
    def print_test_summary(self):
        """Print comprehensive test summary"""
        print("\n" + "="*80)
        print("FASTAAEXTRACTOR INTEGRATION TEST SUMMARY")
        print("="*80)
        
        total_tests = len(self.test_results)
        passed_tests = len([r for r in self.test_results if r['status'] == 'PASS'])
        failed_tests = len([r for r in self.test_results if r['status'] == 'FAIL'])
        error_tests = len([r for r in self.test_results if r['status'] == 'ERROR'])
        
        print(f"Total Tests: {total_tests}")
        print(f"Passed: {passed_tests}")
        print(f"Failed: {failed_tests}")
        print(f"Errors: {error_tests}")
        print(f"Success Rate: {(passed_tests/total_tests)*100:.1f}%")
        
        if failed_tests > 0 or error_tests > 0:
            print(f"\nFailed/Error Tests:")
            for result in self.test_results:
                if result['status'] != 'PASS':
                    print(f"  {result['test']}: {result['status']}")
                    if result['error']:
                        print(f"    Error: {result['error']}")
        
        print("\nDetailed Results:")
        for result in self.test_results:
            status_symbol = "✓" if result['status'] == 'PASS' else "✗"
            print(f"  {status_symbol} {result['test']}: {result['status']}")
        
        print("\n" + "="*80)
        
        if passed_tests == total_tests:
            print("🎉 ALL TESTS PASSED! FastaAAExtractor Integration is ready for production.")
        else:
            print("⚠️  Some tests failed. Please review and fix issues before production use.")
        
        print("="*80)


def main():
    """Main function to run tests"""
    tester = TestFastaAAExtractorIntegration()
    success = tester.run_all_tests()
    
    # Exit with appropriate code
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()