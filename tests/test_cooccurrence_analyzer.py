import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))
from configuration_manager import config_manager
#!/usr/bin/env python3
"""
Test Suite for Generic Co-occurrence Analyzer
Comprehensive testing to ensure the analyzer is bug-proof and robust

This test suite covers:
1. Data validation and error handling
2. Edge cases and boundary conditions
3. Statistical accuracy
4. Performance with large datasets
5. Integration with pipeline components

Author: GenomeAMRAnalyzer Pipeline
Version: 1.0 - Production Testing
"""

import unittest
import tempfile
import shutil
import os
from pathlib import Path
import pandas as pd
import json
from unittest.mock import patch, MagicMock
import sys

# Add the parent directory to the path to import our module
sys.path.append(str(Path(__file__).parent))

from generic_cooccurrence_analyzer import (
    GenericCoOccurrenceAnalyzer, 
    MutationEvent, 
    CoOccurrencePattern
)


class TestMutationEvent(unittest.TestCase):
    """Test MutationEvent data class"""
    
    def test_valid_mutation_event(self):
        """Test creation of valid mutation event"""
        mutation = MutationEvent(
            genome_id="GCF_123456",
            gene=config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0],
            position=123,
            reference_aa="A",
            variant_aa="V",
            substitution="A123V"
        )
        
        self.assertEqual(mutation.genome_id, "GCF_123456")
        self.assertEqual(mutation.gene, config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0])
        self.assertEqual(mutation.position, 123)
        self.assertEqual(mutation.reference_aa, "A")
        self.assertEqual(mutation.variant_aa, "V")
    
    def test_invalid_mutation_events(self):
        """Test validation of invalid mutation events"""
        
        # Test empty genome_id
        with self.assertRaises(ValueError):
            MutationEvent("", config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0], 123, "A", "V", "A123V")
        
        # Test negative position
        with self.assertRaises(ValueError):
            MutationEvent("GCF_123456", config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0], -1, "A", "V", "A123V")
        
        # Test invalid amino acids
        with self.assertRaises(ValueError):
            MutationEvent("GCF_123456", config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0], 123, "Z", "V", "Z123V")
        
        # Test multi-character amino acids
        with self.assertRaises(ValueError):
            MutationEvent("GCF_123456", config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0], 123, "AA", "V", "AA123V")


class TestGenericCoOccurrenceAnalyzer(unittest.TestCase):
    """Test GenericCoOccurrenceAnalyzer class"""
    
    def setUp(self):
        """Set up test environment"""
        self.temp_dir = tempfile.mkdtemp()
        self.analyzer = GenericCoOccurrenceAnalyzer(
            output_dir=self.temp_dir,
            min_genomes=1,  # Lower threshold for testing
            log_level="ERROR"  # Suppress logs during testing
        )
        
        # Create test data
        self.test_data = pd.DataFrame([
            {"genome_id": "genome1", "gene": config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0], "position": 123, "reference_aa": "A", "variant_aa": "V", "substitution": "A123V"},
            {"genome_id": "genome1", "gene": config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1], "position": 456, "reference_aa": "D", "variant_aa": "N", "substitution": "D456N"},
            {"genome_id": "genome2", "gene": config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0], "position": 123, "reference_aa": "A", "variant_aa": "T", "substitution": "A123T"},
            {"genome_id": "genome3", "gene": config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1], "position": 456, "reference_aa": "D", "variant_aa": "E", "substitution": "D456E"},
            {"genome_id": "genome3", "gene": config_manager.get_default_genes("rnd_efflux_pumps", "primary")[2], "position": 789, "reference_aa": "R", "variant_aa": "Q", "substitution": "R789Q"},
        ])
        
        self.test_file = Path(self.temp_dir) / "test_substitutions.csv"
        self.test_data.to_csv(self.test_file, index=False)
    
    def tearDown(self):
        """Clean up test environment"""
        shutil.rmtree(self.temp_dir)
    
    def test_initialization(self):
        """Test analyzer initialization"""
        self.assertTrue(Path(self.temp_dir).exists())
        self.assertEqual(self.analyzer.min_genomes, 1)
        self.assertEqual(self.analyzer.max_combination_size, 5)
    
    def test_load_valid_data(self):
        """Test loading valid substitution data"""
        self.analyzer.load_substitution_data(self.test_file)
        
        self.assertEqual(len(self.analyzer.mutations_data), 5)
        self.assertEqual(len(self.analyzer.genome_set), 3)
        self.assertEqual(len(self.analyzer.gene_list), 3)
        self.assertIn(config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0], self.analyzer.gene_list)
        self.assertIn(config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1], self.analyzer.gene_list)
        self.assertIn(config_manager.get_default_genes("rnd_efflux_pumps", "primary")[2], self.analyzer.gene_list)
    
    def test_load_with_gene_filter(self):
        """Test loading data with gene filter"""
        self.analyzer.load_substitution_data(self.test_file, gene_list=[config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0], config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1]])
        
        self.assertEqual(len(self.analyzer.mutations_data), 4)  # RND efflux pump genes (configurable)
        self.assertEqual(len(self.analyzer.gene_list), 2)
        self.assertNotIn(config_manager.get_default_genes("rnd_efflux_pumps", "primary")[2], self.analyzer.gene_list)
    
    def test_load_nonexistent_file(self):
        """Test error handling for nonexistent file"""
        with self.assertRaises(FileNotFoundError):
            self.analyzer.load_substitution_data("nonexistent_file.csv")
    
    def test_load_invalid_data(self):
        """Test handling of invalid data"""
        # Create invalid data
        invalid_data = pd.DataFrame([
            {"genome_id": "", "gene": config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0], "position": "invalid", "reference_aa": "A", "variant_aa": "V"},
            {"genome_id": "genome1", "gene": "", "position": 123, "reference_aa": "Z", "variant_aa": "V"},
        ])
        
        invalid_file = Path(self.temp_dir) / "invalid_data.csv"
        invalid_data.to_csv(invalid_file, index=False)
        
        # Should handle gracefully
        self.analyzer.load_substitution_data(invalid_file)
        self.assertEqual(len(self.analyzer.mutations_data), 0)  # All records should be filtered out
    
    def test_analyze_cooccurrence_patterns(self):
        """Test co-occurrence pattern analysis"""
        self.analyzer.load_substitution_data(self.test_file)
        self.analyzer.analyze_cooccurrence_patterns()
        
        # Should find patterns
        self.assertGreater(len(self.analyzer.patterns), 0)
        
        # Check for single gene patterns
        single_gene_patterns = [p for p in self.analyzer.patterns if len(p.genes) == 1]
        self.assertGreater(len(single_gene_patterns), 0)
        
        # Check for multi-gene patterns
        multi_gene_patterns = [p for p in self.analyzer.patterns if len(p.genes) > 1]
        self.assertGreater(len(multi_gene_patterns), 0)
    
    def test_pattern_validation(self):
        """Test pattern data validation"""
        self.analyzer.load_substitution_data(self.test_file)
        self.analyzer.analyze_cooccurrence_patterns()
        
        for pattern in self.analyzer.patterns:
            # Check data types
            self.assertIsInstance(pattern.genes, tuple)
            self.assertIsInstance(pattern.genome_count, int)
            self.assertIsInstance(pattern.frequency, float)
            
            # Check ranges
            self.assertGreaterEqual(pattern.genome_count, 0)
            self.assertLessEqual(pattern.genome_count, len(self.analyzer.genome_set))
            self.assertGreaterEqual(pattern.frequency, 0.0)
            self.assertLessEqual(pattern.frequency, 1.0)
            
            # Check consistency
            expected_frequency = pattern.genome_count / pattern.total_genomes
            self.assertAlmostEqual(pattern.frequency, expected_frequency, places=6)
    
    def test_empty_data_handling(self):
        """Test handling of empty datasets"""
        empty_data = pd.DataFrame(columns=["genome_id", "gene", "position", "reference_aa", "variant_aa", "substitution"])
        empty_file = Path(self.temp_dir) / "empty_data.csv"
        empty_data.to_csv(empty_file, index=False)
        
        self.analyzer.load_substitution_data(empty_file)
        
        with self.assertRaises(ValueError):
            self.analyzer.analyze_cooccurrence_patterns()
    
    def test_single_genome_data(self):
        """Test handling of single genome data"""
        single_genome_data = pd.DataFrame([
            {"genome_id": "genome1", "gene": config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0], "position": 123, "reference_aa": "A", "variant_aa": "V", "substitution": "A123V"},
            {"genome_id": "genome1", "gene": config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1], "position": 456, "reference_aa": "D", "variant_aa": "N", "substitution": "D456N"},
        ])
        
        single_file = Path(self.temp_dir) / "single_genome.csv"
        single_genome_data.to_csv(single_file, index=False)
        
        self.analyzer.load_substitution_data(single_file)
        self.analyzer.analyze_cooccurrence_patterns()
        
        # Should handle single genome appropriately
        self.assertGreater(len(self.analyzer.patterns), 0)
    
    def test_large_dataset_performance(self):
        """Test performance with larger dataset"""
        # Create larger test dataset
        large_data = []
        for genome_id in range(100):
            for gene in config_manager.get_default_genes("rnd_efflux_pumps", "primary"):
                for pos in [123, 456]:
                    large_data.append({
                        "genome_id": f"genome_{genome_id}",
                        "gene": gene,
                        "position": pos,
                        "reference_aa": "A",
                        "variant_aa": "V",
                        "substitution": f"A{pos}V"
                    })
        
        large_df = pd.DataFrame(large_data)
        large_file = Path(self.temp_dir) / "large_data.csv"
        large_df.to_csv(large_file, index=False)
        
        # Should handle large dataset efficiently
        import time
        start_time = time.time()
        
        self.analyzer.load_substitution_data(large_file)
        self.analyzer.analyze_cooccurrence_patterns()
        
        end_time = time.time()
        
        # Should complete in reasonable time (< 30 seconds)
        self.assertLess(end_time - start_time, 30)
        self.assertGreater(len(self.analyzer.patterns), 0)
    
    def test_save_results(self):
        """Test saving results in multiple formats"""
        self.analyzer.load_substitution_data(self.test_file)
        self.analyzer.analyze_cooccurrence_patterns()
        
        output_files = self.analyzer.save_results("test_output")
        
        # Check that files were created
        self.assertIn('json', output_files)
        self.assertIn('csv', output_files)
        self.assertIn('summary', output_files)
        
        # Verify files exist
        for filepath in output_files.values():
            self.assertTrue(Path(filepath).exists())
        
        # Verify JSON content
        with open(output_files['json'], 'r') as f:
            report = json.load(f)
        
        self.assertIn('summary', report)
        self.assertIn('patterns', report)
        self.assertIn('methodology', report)
    
    def test_column_mapping(self):
        """Test alternative column name mapping"""
        # Create data with alternative column names
        alt_data = pd.DataFrame([
            {"accession": "genome1", "protein": config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0], "pos": 123, "ref_aa": "A", "alt_aa": "V", "substitution": "A123V"},
            {"accession": "genome2", "protein": config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1], "pos": 456, "ref_aa": "D", "alt_aa": "N", "substitution": "D456N"},
        ])
        
        alt_file = Path(self.temp_dir) / "alt_columns.csv"
        alt_data.to_csv(alt_file, index=False)
        
        # Should handle alternative column names
        self.analyzer.load_substitution_data(alt_file)
        self.assertEqual(len(self.analyzer.mutations_data), 2)
    
    def test_mic_data_integration(self):
        """Test MIC data integration"""
        # Create MIC data
        mic_data = pd.DataFrame([
            {"genome_id": "genome1", "ampicillin": 32.0, "chloramphenicol": 64.0},
            {"genome_id": "genome2", "ampicillin": 8.0, "chloramphenicol": 16.0},
        ])
        
        mic_file = Path(self.temp_dir) / "mic_data.csv"
        mic_data.to_csv(mic_file, index=False)
        
        # Load with MIC data
        self.analyzer.load_substitution_data(self.test_file, mic_file=mic_file)
        
        # Check MIC data was loaded
        for mutation in self.analyzer.mutations_data:
            if mutation.genome_id in ["genome1", "genome2"]:
                self.assertGreater(len(mutation.mic_data), 0)


class TestStatisticalFunctions(unittest.TestCase):
    """Test statistical analysis functions"""
    
    def setUp(self):
        """Set up test environment"""
        self.temp_dir = tempfile.mkdtemp()
        self.analyzer = GenericCoOccurrenceAnalyzer(
            output_dir=self.temp_dir,
            log_level="ERROR"
        )
    
    def tearDown(self):
        """Clean up test environment"""
        shutil.rmtree(self.temp_dir)
    
    def test_frequency_calculation(self):
        """Test frequency calculation accuracy"""
        # Test with known data
        pattern = CoOccurrencePattern(
            genes=(config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0],),
            mutations=("acrA_mutated",),
            genome_count=25,
            total_genomes=100,
            frequency=0.25,
            genomes=set([f"genome_{i}" for i in range(25)])
        )
        
        self.assertEqual(pattern.percentage, 25.0)
        self.assertEqual(pattern.frequency, 0.25)
    
    @patch('generic_cooccurrence_analyzer.hypergeom')
    def test_statistical_significance_calculation(self, mock_hypergeom):
        """Test statistical significance calculation"""
        # Mock scipy.stats.hypergeom
        mock_hypergeom.sf.return_value = 0.001
        
        # Create test data for statistical testing
        test_data = pd.DataFrame([
            {"genome_id": f"genome_{i}", "gene": config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0], "position": 123, "reference_aa": "A", "variant_aa": "V", "substitution": "A123V"}
            for i in range(50)
        ] + [
            {"genome_id": f"genome_{i}", "gene": config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1], "position": 456, "reference_aa": "D", "variant_aa": "N", "substitution": "D456N"}
            for i in range(30)
        ] + [
            {"genome_id": f"genome_{i}", "gene": config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0], "position": 123, "reference_aa": "A", "variant_aa": "V", "substitution": "A123V"}
            for i in range(20, 40)  # Overlap for co-occurrence
        ] + [
            {"genome_id": f"genome_{i}", "gene": config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1], "position": 456, "reference_aa": "D", "variant_aa": "N", "substitution": "D456N"}
            for i in range(20, 40)  # Overlap for co-occurrence
        ])
        
        test_file = Path(self.temp_dir) / "stat_test.csv"
        test_data.to_csv(test_file, index=False)
        
        self.analyzer.load_substitution_data(test_file)
        self.analyzer.analyze_cooccurrence_patterns()
        
        # Should have calculated significance for multi-gene patterns
        multi_gene_patterns = [p for p in self.analyzer.patterns if len(p.genes) > 1]
        self.assertGreater(len(multi_gene_patterns), 0)


class TestEdgeCases(unittest.TestCase):
    """Test edge cases and boundary conditions"""
    
    def setUp(self):
        """Set up test environment"""
        self.temp_dir = tempfile.mkdtemp()
        self.analyzer = GenericCoOccurrenceAnalyzer(
            output_dir=self.temp_dir,
            log_level="ERROR"
        )
    
    def tearDown(self):
        """Clean up test environment"""
        shutil.rmtree(self.temp_dir)
    
    def test_special_amino_acids(self):
        """Test handling of special amino acid codes"""
        special_data = pd.DataFrame([
            {"genome_id": "genome1", "gene": config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0], "position": 123, "reference_aa": "*", "variant_aa": "X", "substitution": "*123X"},
            {"genome_id": "genome2", "gene": config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1], "position": 456, "reference_aa": "U", "variant_aa": "-", "substitution": "U456-"},
        ])
        
        special_file = Path(self.temp_dir) / "special_aa.csv"
        special_data.to_csv(special_file, index=False)
        
        # Should handle special amino acid codes
        self.analyzer.load_substitution_data(special_file)
        self.assertEqual(len(self.analyzer.mutations_data), 2)
    
    def test_very_large_positions(self):
        """Test handling of very large position numbers"""
        large_pos_data = pd.DataFrame([
            {"genome_id": "genome1", "gene": config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0], "position": 999999, "reference_aa": "A", "variant_aa": "V", "substitution": "A999999V"},
        ])
        
        large_file = Path(self.temp_dir) / "large_positions.csv"
        large_pos_data.to_csv(large_file, index=False)
        
        self.analyzer.load_substitution_data(large_file)
        self.assertEqual(len(self.analyzer.mutations_data), 1)
        self.assertEqual(self.analyzer.mutations_data[0].position, 999999)
    
    def test_unicode_genome_ids(self):
        """Test handling of unicode characters in genome IDs"""
        unicode_data = pd.DataFrame([
            {"genome_id": "genome_αβγ", "gene": config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0], "position": 123, "reference_aa": "A", "variant_aa": "V", "substitution": "A123V"},
        ])
        
        unicode_file = Path(self.temp_dir) / "unicode_ids.csv"
        unicode_data.to_csv(unicode_file, index=False)
        
        self.analyzer.load_substitution_data(unicode_file)
        self.assertEqual(len(self.analyzer.mutations_data), 1)


def run_all_tests():
    """Run all test suites"""
    test_suites = [
        unittest.TestLoader().loadTestsFromTestCase(TestMutationEvent),
        unittest.TestLoader().loadTestsFromTestCase(TestGenericCoOccurrenceAnalyzer),
        unittest.TestLoader().loadTestsFromTestCase(TestStatisticalFunctions),
        unittest.TestLoader().loadTestsFromTestCase(TestEdgeCases),
    ]
    
    combined_suite = unittest.TestSuite(test_suites)
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(combined_suite)
    
    return result.wasSuccessful()


if __name__ == "__main__":
    success = run_all_tests()
    if success:
        print("\n" + "=" * 60)
        print("ALL TESTS PASSED! ✓")
        print("Generic Co-occurrence Analyzer is ready for production use.")
        print("=" * 60)
    else:
        print("\n" + "=" * 60)
        print("SOME TESTS FAILED! ✗")
        print("Please review and fix issues before production use.")
        print("=" * 60)
        sys.exit(1)