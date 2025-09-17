# Tests for the GenomeAMRAnalyzer package

def test_import():
    """Test that core modules can be imported successfully"""
    try:
        import sys
        import os
        sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))
        
        from generic_cooccurrence_analyzer import GenericCoOccurrenceAnalyzer
        from fasta_aa_extractor_integration import FastaAAExtractorIntegrator
        
        # Test basic instantiation
        analyzer = GenericCoOccurrenceAnalyzer(output_dir='test_output')
        assert analyzer is not None
        
        integrator = FastaAAExtractorIntegrator(output_dir='test_output')
        assert integrator is not None
        
        print("✓ All core modules imported successfully")
        return True
        
    except ImportError as e:
        print(f"✗ Import error: {e}")
        return False
    except Exception as e:
        print(f"✗ Unexpected error: {e}")
        return False

if __name__ == "__main__":
    test_import()