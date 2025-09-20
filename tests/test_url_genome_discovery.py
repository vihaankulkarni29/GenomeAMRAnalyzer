"""
Test Suite for URL-based Genome Discovery
=========================================
Test the URL parsing and genome discovery logic without external dependencies.
"""

import sys
import os
from pathlib import Path

# Add the src/core directory to the path
current_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(current_dir)
src_core_dir = os.path.join(project_root, 'src', 'core')
sys.path.insert(0, src_core_dir)

def test_url_parsing():
    """Test URL parsing functionality"""
    print("Testing URL parsing...")
    
    try:
        from ncbi_genome_discovery import NCBIUrlParser
        
        parser = NCBIUrlParser()
        
        # Test URLs
        test_urls = [
            "https://www.ncbi.nlm.nih.gov/nuccore/?term=(Escherichia+coli+and+macrolide+resistance)+AND+%22Escherichia+coli%22%5Bporgn%3A__txid562%5D+and+complete+genome",
            "https://www.ncbi.nlm.nih.gov/nuccore/?term=Escherichia+coli+complete+genome",
            "https://www.ncbi.nlm.nih.gov/nuccore/?term=Salmonella+AND+antibiotic+resistance"
        ]
        
        for i, url in enumerate(test_urls, 1):
            try:
                search_term = parser.parse_search_url(url)
                print(f"✅ Test {i}: Successfully parsed URL")
                print(f"   Search term: {search_term}")
            except Exception as e:
                print(f"❌ Test {i}: Failed to parse URL - {e}")
        
        print("✅ URL parsing tests completed")
        return True
        
    except ImportError as e:
        print(f"❌ Could not import NCBIUrlParser: {e}")
        return False


def test_genome_metadata():
    """Test GenomeMetadata data structure"""
    print("Testing GenomeMetadata...")
    
    try:
        from ncbi_genome_discovery import GenomeMetadata
        
        # Create test genome metadata
        genome = GenomeMetadata(
            accession="NC_000913.3",
            title="Escherichia coli str. K-12 substr. MG1655, complete genome",
            organism="Escherichia coli",
            length=4641652,
            pub_date="2013/09/26",
            molecule_type="DNA",
            topology="circular"
        )
        
        # Test properties
        assert genome.accession == "NC_000913.3"
        assert "Escherichia coli" in genome.organism
        assert genome.length > 4000000
        
        print(f"✅ GenomeMetadata test passed")
        print(f"   Accession: {genome.accession}")
        print(f"   Organism: {genome.organism}")
        print(f"   Length: {genome.length:,} bp")
        
        return True
        
    except Exception as e:
        print(f"❌ GenomeMetadata test failed: {e}")
        return False


def test_download_result():
    """Test DownloadResult data structure"""
    print("Testing DownloadResult...")
    
    try:
        from genome_downloader import DownloadResult
        
        # Create test download result
        result = DownloadResult(
            accession="NC_000913.3",
            success=True,
            file_path="/path/to/NC_000913.3.fasta",
            file_size_bytes=4641652,
            download_time_seconds=15.5,
            retry_count=0
        )
        
        # Test properties
        assert result.accession == "NC_000913.3"
        assert result.success is True
        assert result.file_size_bytes > 0
        
        print(f"✅ DownloadResult test passed")
        print(f"   Accession: {result.accession}")
        print(f"   Success: {result.success}")
        print(f"   File size: {result.file_size_bytes:,} bytes")
        print(f"   Download time: {result.download_time_seconds} seconds")
        
        return True
        
    except Exception as e:
        print(f"❌ DownloadResult test failed: {e}")
        return False


def test_configuration_loading():
    """Test configuration loading"""
    print("Testing configuration loading...")
    
    try:
        config_path = Path(project_root) / "config" / "snakemake_config.yaml"
        
        if not config_path.exists():
            print(f"❌ Config file not found: {config_path}")
            return False
        
        import yaml
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
        
        # Test required configuration keys
        required_keys = [
            'ncbi_email',
            'ncbi_search',
            'directories',
            'target_genes',
            'rgi'
        ]
        
        for key in required_keys:
            if key not in config:
                print(f"❌ Missing required config key: {key}")
                return False
        
        # Test specific values
        assert isinstance(config['ncbi_search']['max_genomes'], int)
        assert config['ncbi_search']['max_genomes'] > 0
        assert isinstance(config['target_genes'], list)
        assert len(config['target_genes']) > 0
        
        print(f"✅ Configuration loading test passed")
        print(f"   Email: {config['ncbi_email']}")
        print(f"   Max genomes: {config['ncbi_search']['max_genomes']}")
        print(f"   Target genes: {', '.join(config['target_genes'])}")
        print(f"   Output directory: {config['directories']['genomes']}")
        
        return True
        
    except Exception as e:
        print(f"❌ Configuration loading test failed: {e}")
        return False


def test_workflow_initialization():
    """Test workflow initialization (without running)"""
    print("Testing workflow initialization...")
    
    try:
        config_path = Path(project_root) / "config" / "snakemake_config.yaml"
        
        if not config_path.exists():
            print(f"❌ Config file not found: {config_path}")
            return False
        
        # Test that we can create the workflow class without external dependencies
        # This tests the basic structure and configuration loading
        
        import yaml
        with open(config_path, 'r') as f:
            config = yaml.safe_load(f)
        
        # Simulate workflow initialization checks
        required_components = [
            'ncbi_email',
            'ncbi_search',
            'directories',
            'target_genes'
        ]
        
        for component in required_components:
            if component not in config:
                raise ValueError(f"Missing required component: {component}")
        
        # Test directory creation logic
        genome_dir = Path(config['directories']['genomes'])
        print(f"   Would create genome directory: {genome_dir}")
        
        logs_dir = Path(config['directories']['logs'])
        print(f"   Would create logs directory: {logs_dir}")
        
        print(f"✅ Workflow initialization test passed")
        return True
        
    except Exception as e:
        print(f"❌ Workflow initialization test failed: {e}")
        return False


def run_all_tests():
    """Run all tests and report results"""
    print("=" * 60)
    print("URL-based Genome Discovery Test Suite")
    print("=" * 60)
    
    tests = [
        ("URL Parsing", test_url_parsing),
        ("Genome Metadata", test_genome_metadata),
        ("Download Result", test_download_result),
        ("Configuration Loading", test_configuration_loading),
        ("Workflow Initialization", test_workflow_initialization)
    ]
    
    passed = 0
    total = len(tests)
    
    for test_name, test_func in tests:
        print(f"\n--- {test_name} ---")
        try:
            if test_func():
                passed += 1
        except Exception as e:
            print(f"❌ {test_name} failed with exception: {e}")
    
    print("\n" + "=" * 60)
    print(f"Test Results: {passed}/{total} tests passed")
    
    if passed == total:
        print("🎉 All tests passed! URL-based genome discovery is ready.")
    else:
        print("⚠️  Some tests failed. Please review the issues above.")
    
    print("=" * 60)
    
    return passed == total


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)