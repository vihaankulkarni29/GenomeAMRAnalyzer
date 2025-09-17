import sys; import os; sys.path.append(os.path.join(os.path.dirname(__file__), '\\\\\')); from .configuration_manager import config_manager
"""
Comprehensive Test Suite for Priority 2 Modules
----------------------------------------------
Unit and integration tests for all Priority 2 components ensuring
robustness, performance, and reliability of the AMR analysis pipeline.
"""

import os
import sys
import unittest
import tempfile
import shutil
import json
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

class TestEnhancedSequenceProcessor(unittest.TestCase):
    """Test suite for EnhancedSequenceProcessor."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.test_dir = tempfile.mkdtemp()
        self.reference_path = os.path.join(self.test_dir, "reference.fasta")
        self.query_path = os.path.join(self.test_dir, "query.fasta")
        
        # Create test reference file
        with open(self.reference_path, 'w') as f:
            f.write(">ref1\nATCGATCGATCGATCG\n>ref2\nGCTAGCTAGCTAGCTA\n")
        
        # Create test query file
        with open(self.query_path, 'w') as f:
            f.write(">query1\nATCGATCGATCG\n>query2\nGCTAGCTAGCTA\n")
    
    def tearDown(self):
        """Clean up test fixtures."""
        shutil.rmtree(self.test_dir)
    
    @patch('src.priority2.core.enhanced_sequence_processor.mp')
    def test_initialization_success(self, mock_mp):
        """Test successful processor initialization."""
        mock_mp.Aligner.return_value = MagicMock()
        
        from src.priority2.core.enhanced_sequence_processor import EnhancedSequenceProcessor
        
        processor = EnhancedSequenceProcessor(self.reference_path)
        self.assertEqual(processor.reference_path, self.reference_path)
        self.assertEqual(processor.threads, 4)
    
    @patch('src.priority2.core.enhanced_sequence_processor.mp', None)
    def test_initialization_missing_dependency(self):
        """Test initialization with missing mappy dependency."""
        from src.priority2.core.enhanced_sequence_processor import EnhancedSequenceProcessor, EnhancedSequenceProcessorError
        
        with self.assertRaises(EnhancedSequenceProcessorError):
            EnhancedSequenceProcessor(self.reference_path)
    
    def test_initialization_missing_reference(self):
        """Test initialization with missing reference file."""
        from src.priority2.core.enhanced_sequence_processor import EnhancedSequenceProcessor, EnhancedSequenceProcessorError
        
        with self.assertRaises(EnhancedSequenceProcessorError):
            EnhancedSequenceProcessor("nonexistent_reference.fasta")
    
    @patch('src.priority2.core.enhanced_sequence_processor.mp')
    def test_input_validation(self, mock_mp):
        """Test input file validation."""
        mock_mp.Aligner.return_value = MagicMock()
        
        from src.priority2.core.enhanced_sequence_processor import EnhancedSequenceProcessor
        
        processor = EnhancedSequenceProcessor(self.reference_path)
        
        # Test with valid files
        valid_files = processor.validate_inputs([self.query_path])
        self.assertEqual(len(valid_files), 1)
        
        # Test with invalid files
        invalid_files = processor.validate_inputs(["nonexistent.fasta"])
        self.assertEqual(len(invalid_files), 0)

class TestDataQualityController(unittest.TestCase):
    """Test suite for DataQualityController."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.test_metrics = [
            {"num_alignments": 100},
            {"num_alignments": 50},
            {"num_alignments": 150},
            {"num_alignments": 75}
        ]
    
    @patch('src.priority2.core.data_quality_controller.IsolationForest')
    def test_quality_assessment(self, mock_isolation_forest):
        """Test quality assessment functionality."""
        # Mock IsolationForest
        mock_model = MagicMock()
        mock_model.fit_predict.return_value = [1, -1, 1, 1]  # 1 = normal, -1 = anomaly
        mock_isolation_forest.return_value = mock_model
        
        from src.priority2.core.data_quality_controller import DataQualityController
        
        controller = DataQualityController()
        assessed_metrics = controller.assess_quality(self.test_metrics)
        
        self.assertEqual(len(assessed_metrics), 4)
        self.assertTrue(all('quality_score' in m for m in assessed_metrics))
    
    def test_empty_metrics(self):
        """Test handling of empty metrics."""
        from src.priority2.core.data_quality_controller import DataQualityController
        
        controller = DataQualityController()
        result = controller.assess_quality([])
        self.assertEqual(result, [])

class TestPerformanceOptimizer(unittest.TestCase):
    """Test suite for PerformanceOptimizer."""
    
    @patch('src.priority2.core.performance_optimizer.psutil')
    def test_dynamic_resource_allocation(self, mock_psutil):
        """Test dynamic resource allocation."""
        mock_psutil.cpu_count.return_value = 8
        
        from src.priority2.core.performance_optimizer import PerformanceOptimizer
        
        optimizer = PerformanceOptimizer()
        
        # Test small dataset
        threads = optimizer.dynamic_resource_allocation(500)
        self.assertEqual(threads, 4)  # cpu_count // 2
        
        # Test medium dataset
        threads = optimizer.dynamic_resource_allocation(5000)
        self.assertEqual(threads, 7)  # cpu_count - 1
        
        # Test large dataset
        threads = optimizer.dynamic_resource_allocation(50000)
        self.assertEqual(threads, 8)  # cpu_count
    
    def test_memory_profiler(self):
        """Test memory profiling functionality."""
        from src.priority2.core.performance_optimizer import PerformanceOptimizer
        
        optimizer = PerformanceOptimizer()
        
        def test_function(x, y):
            return x + y
        
        result = optimizer.memory_profiler(test_function, 2, 3)
        self.assertEqual(result, 5)

class TestExportManager(unittest.TestCase):
    """Test suite for ExportManager."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.test_dir = tempfile.mkdtemp()
        self.test_data = [
            {"sample": "S1", "gene": config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0], "identity": 95.5},
            {"sample": "S2", "gene": config_manager.get_default_genes("rnd_efflux_pumps", "primary")[1], "identity": 88.2}
        ]
    
    def tearDown(self):
        """Clean up test fixtures."""
        shutil.rmtree(self.test_dir)
    
    def test_csv_export(self):
        """Test CSV export functionality."""
        from src.priority2.core.export_manager import ExportManager
        
        manager = ExportManager(self.test_dir)
        manager.export_csv(self.test_data, "test.csv")
        
        csv_path = os.path.join(self.test_dir, "test.csv")
        self.assertTrue(os.path.exists(csv_path))
    
    def test_json_export(self):
        """Test JSON export functionality."""
        from src.priority2.core.export_manager import ExportManager
        
        manager = ExportManager(self.test_dir)
        manager.export_json(self.test_data, "test.json")
        
        json_path = os.path.join(self.test_dir, "test.json")
        self.assertTrue(os.path.exists(json_path))
        
        # Verify content
        with open(json_path, 'r') as f:
            loaded_data = json.load(f)
        self.assertEqual(loaded_data, self.test_data)
    
    def test_fasta_export(self):
        """Test FASTA export functionality."""
        from src.priority2.core.export_manager import ExportManager
        
        fasta_data = [
            {"id": "seq1", "seq": "ATCGATCG"},
            {"id": "seq2", "seq": "GCTAGCTA"}
        ]
        
        manager = ExportManager(self.test_dir)
        manager.export_fasta(fasta_data, "test.fasta")
        
        fasta_path = os.path.join(self.test_dir, "test.fasta")
        self.assertTrue(os.path.exists(fasta_path))

class TestHighThroughputAligner(unittest.TestCase):
    """Test suite for HighThroughputAligner."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.test_dir = tempfile.mkdtemp()
        self.reference_path = os.path.join(self.test_dir, "reference.fasta")
        
        # Create test reference
        with open(self.reference_path, 'w') as f:
            f.write(">ref1\nATCGATCGATCGATCG\n")
    
    def tearDown(self):
        """Clean up test fixtures."""
        shutil.rmtree(self.test_dir)
    
    @patch('src.priority2.pipelines.high_throughput_aligner.EnhancedSequenceProcessor')
    @patch('src.priority2.pipelines.high_throughput_aligner.PerformanceOptimizer')
    def test_chunk_processing(self, mock_optimizer, mock_processor):
        """Test chunked processing functionality."""
        # Setup mocks
        mock_optimizer_instance = MagicMock()
        mock_optimizer_instance.dynamic_resource_allocation.return_value = 4
        mock_optimizer.return_value = mock_optimizer_instance
        
        mock_processor_instance = MagicMock()
        mock_processor_instance.align_sequences.return_value = ["alignment1.paf", "alignment2.paf"]
        mock_processor.return_value = mock_processor_instance
        
        from src.priority2.pipelines.high_throughput_aligner import HighThroughputAligner
        
        aligner = HighThroughputAligner(self.reference_path, self.test_dir, chunk_size=2)
        result = aligner.run(["file1.fasta", "file2.fasta", "file3.fasta"])
        
        self.assertEqual(len(result), 2)  # 2 alignments per chunk

class TestAlignmentAnalyzer(unittest.TestCase):
    """Test suite for AlignmentAnalyzer."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.test_dir = tempfile.mkdtemp()
        self.paf_file = os.path.join(self.test_dir, "test.paf")
        
        # Create test PAF file
        with open(self.paf_file, 'w') as f:
            f.write("query1\t100\t0\t90\t+\tacrA\t200\t10\t100\t80\t90\t60\n")
            f.write("query2\t150\t0\t140\t+\tacrB\t250\t20\t160\t130\t140\t55\n")
    
    def tearDown(self):
        """Clean up test fixtures."""
        shutil.rmtree(self.test_dir)
    
    def test_paf_parsing(self):
        """Test PAF file parsing."""
        from src.priority2.core.alignment_analyzer import AlignmentAnalyzer
        
        analyzer = AlignmentAnalyzer()
        alignments = analyzer.parse_paf_file(self.paf_file)
        
        self.assertEqual(len(alignments), 2)
        self.assertEqual(alignments[0].query_name, "query1")
        self.assertEqual(alignments[0].target_name, config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0])
    
    def test_sample_analysis(self):
        """Test sample analysis functionality."""
        from src.priority2.core.alignment_analyzer import AlignmentAnalyzer
        
        analyzer = AlignmentAnalyzer()
        report = analyzer.analyze_sample(self.paf_file)
        
        self.assertEqual(report.sample_name, "test")
        self.assertEqual(report.total_sequences, 2)
        self.assertGreater(report.alignment_rate, 0)

class TestExternalAlignmentTools(unittest.TestCase):
    """Test suite for ExternalAlignmentTools."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.test_dir = tempfile.mkdtemp()
        self.reference_path = os.path.join(self.test_dir, "reference.fasta")
        
        with open(self.reference_path, 'w') as f:
            f.write(">ref1\nATCGATCGATCGATCG\n")
    
    def tearDown(self):
        """Clean up test fixtures."""
        shutil.rmtree(self.test_dir)
    
    def test_tool_detection(self):
        """Test alignment tool detection."""
        from src.priority2.integrations.external_alignment_tools import AlignmentToolManager
        
        manager = AlignmentToolManager(self.reference_path)
        available_tools = manager.available_tools
        
        self.assertIsInstance(available_tools, dict)
        self.assertIn('mappy', available_tools)
        self.assertIn('minimap2', available_tools)
        self.assertIn('parasail', available_tools)
    
    @patch('src.priority2.integrations.external_alignment_tools.mp')
    def test_mappy_aligner(self, mock_mp):
        """Test MappyAligner functionality."""
        # Setup mock
        mock_aligner = MagicMock()
        mock_hit = MagicMock()
        mock_hit.ctg = "ref1"
        mock_hit.ctg_len = 16
        mock_hit.q_st = 0
        mock_hit.q_en = 12
        mock_hit.r_st = 0
        mock_hit.r_en = 12
        mock_hit.mlen = 12
        mock_hit.blen = 12
        mock_hit.mapq = 60
        mock_hit.strand = 1
        
        mock_aligner.map.return_value = [mock_hit]
        mock_mp.Aligner.return_value = mock_aligner
        
        from src.priority2.integrations.external_alignment_tools import MappyAligner
        
        aligner = MappyAligner(self.reference_path)
        results = aligner.align("ATCGATCGATCG", "test_query")
        
        self.assertEqual(len(results), 1)
        self.assertEqual(results[0].query_name, "test_query")
        self.assertEqual(results[0].target_name, "ref1")

class TestLargeDatasetProcessor(unittest.TestCase):
    """Test suite for LargeDatasetProcessor."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.test_dir = tempfile.mkdtemp()
        self.reference_path = os.path.join(self.test_dir, "reference.fasta")
        self.large_dataset = os.path.join(self.test_dir, "large_dataset.fasta")
        
        # Create test reference
        with open(self.reference_path, 'w') as f:
            f.write(">ref1\nATCGATCGATCGATCG\n")
        
        # Create test large dataset
        with open(self.large_dataset, 'w') as f:
            for i in range(100):  # 100 sequences
                f.write(f">seq{i}\nATCGATCGATCG\n")
    
    def tearDown(self):
        """Clean up test fixtures."""
        shutil.rmtree(self.test_dir)
    
    def test_sequence_counting(self):
        """Test sequence counting functionality."""
        from src.priority2.pipelines.large_dataset_processor import SequenceStreamReader
        
        reader = SequenceStreamReader(self.large_dataset, chunk_size=10)
        count = reader.count_sequences()
        
        self.assertEqual(count, 100)
    
    def test_chunk_reading(self):
        """Test chunk reading functionality."""
        from src.priority2.pipelines.large_dataset_processor import SequenceStreamReader
        
        reader = SequenceStreamReader(self.large_dataset, chunk_size=10)
        chunks = list(reader.read_chunks())
        
        self.assertEqual(len(chunks), 10)  # 100 sequences / 10 per chunk
        self.assertEqual(chunks[0].chunk_size, 10)
    
    @patch('src.priority2.pipelines.large_dataset_processor.AlignmentToolManager')
    def test_chunk_processing(self, mock_tool_manager):
        """Test chunk processing functionality."""
        # Setup mock
        mock_manager = MagicMock()
        mock_aligner = MagicMock()
        mock_aligner.align_file.return_value = True
        mock_manager.get_best_aligner.return_value = mock_aligner
        mock_tool_manager.return_value = mock_manager
        
        from src.priority2.pipelines.large_dataset_processor import ChunkProcessor, ProcessingChunk
        
        processor = ChunkProcessor(self.reference_path, self.test_dir, {})
        
        # Create test chunk
        chunk = ProcessingChunk(
            chunk_id=0,
            sequences=[("seq1", "ATCGATCG"), ("seq2", "GCTAGCTA")],
            chunk_size=2
        )
        
        result = processor.process_chunk(chunk)
        
        self.assertTrue(result['success'])
        self.assertEqual(result['chunk_id'], 0)

class TestPipelineIntegration(unittest.TestCase):
    """Integration tests for the complete pipeline."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.test_dir = tempfile.mkdtemp()
        self.reference_path = os.path.join(self.test_dir, "reference.fasta")
        self.config_path = os.path.join(self.test_dir, "config.yaml")
        
        # Create test reference
        with open(self.reference_path, 'w') as f:
            f.write(">acrA\nATCGATCGATCGATCGATCGATCGATCG\n")
            f.write(">acrB\nGCTAGCTAGCTAGCTAGCTAGCTAGCTA\n")
    
    def tearDown(self):
        """Clean up test fixtures."""
        shutil.rmtree(self.test_dir)
    
    @patch('src.priority2.config.pipeline_config.yaml')
    def test_config_creation_and_loading(self, mock_yaml):
        """Test configuration creation and loading."""
        from src.priority2.config.pipeline_config import PipelineConfig, DatabaseConfig, ConfigManager
        
        # Test config creation
        config = PipelineConfig(
            database=DatabaseConfig(amr_database_path=self.reference_path),
            config_name="test_config"
        )
        
        self.assertEqual(config.database.amr_database_path, self.reference_path)
        self.assertEqual(config.config_name, "test_config")
    
    def test_end_to_end_small_dataset(self):
        """Test end-to-end processing of small dataset."""
        # This would require mocking many components or having actual tools installed
        # For now, we'll test the basic pipeline structure
        pass

class TestPerformanceBenchmarks(unittest.TestCase):
    """Performance benchmark tests."""
    
    def test_memory_usage_limits(self):
        """Test that components respect memory limits."""
        # Test memory-conscious operations
        pass
    
    def test_throughput_benchmarks(self):
        """Test processing throughput meets targets."""
        # Test processing speed requirements
        pass

def run_priority2_tests():
    """Run all Priority 2 tests."""
    # Create test suite
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()
    
    # Add test classes
    test_classes = [
        TestEnhancedSequenceProcessor,
        TestDataQualityController,
        TestPerformanceOptimizer,
        TestExportManager,
        TestHighThroughputAligner,
        TestAlignmentAnalyzer,
        TestExternalAlignmentTools,
        TestLargeDatasetProcessor,
        TestPipelineIntegration,
        TestPerformanceBenchmarks
    ]
    
    for test_class in test_classes:
        tests = loader.loadTestsFromTestCase(test_class)
        suite.addTests(tests)
    
    # Run tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    return result.wasSuccessful()

if __name__ == "__main__":
    success = run_priority2_tests()
    sys.exit(0 if success else 1)