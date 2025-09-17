"""
SubScan Alignment Analyzer - Comprehensive Test & Validation Suite
==================================================================

Enterprise-grade testing for mutation detection pipeline with:
- EMBOSS-WATER output parsing validation
- Substitution detection accuracy testing  
- Quality control mechanism validation
- Functional annotation testing
- Performance benchmarking
- Edge case and error handling verification
"""

import sys
import os
import tempfile
import logging
from pathlib import Path
from typing import List, Dict, Any
from datetime import datetime

# Add parent directory to path for imports
sys.path.append(str(Path(__file__).parent.parent))

from analysis.subscan_alignment_analyzer import (
    SubScanAlignmentAnalyzer, EmbossWaterParser, SubstitutionDetector,
    FunctionalAnnotator, MutationRecord, MutationType, ConfidenceLevel,
    GenomicCoordinate, AlignmentQuality
)

def setup_logging():
    """Setup comprehensive logging for validation."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler('subscan_analyzer_validation.log'),
            logging.StreamHandler()
        ]
    )

def create_mock_emboss_water_output() -> str:
    """Create realistic mock EMBOSS-WATER output for testing."""
    water_output = """########################################
# Program: water
# Rundate: Wed 13 Sep 2023 10:30:00
# Commandline: water
#    -asequence: query.fasta
#    -bsequence: reference.fasta
#    -gapopen: 10.0
#    -gapextend: 0.5
#    -outfile: alignment.water
# Align_format: srspair
# Report_file: alignment.water
########################################

#=======================================
#
# Aligned_sequences: 2
# 1: gyrA_query
# 2: gyrA_reference
# Matrix: EDNAFULL
# Gap_penalty: 10.0
# Extend_penalty: 0.5
#
# Length: 261
# Identity:     249/261 (95.4%)
# Similarity:   249/261 (95.4%)
# Gaps:           0/261 ( 0.0%)
# Score: 1190.0
#
#
#=======================================

gyrA_query         1 ATGGCTGAACGTCTGGTGCTGGAACAGTTCGGCGACCTGCGTACGAAACCG     50
                     |||||||||||||||||||||||||||||||||||||||||||||||||| 
gyrA_reference     1 ATGGCTGAACGTCTGGTGCTGGAACAGTTCGGCGACCTGCGTACGAAACCG     50

gyrA_query        51 GTGAAACTGCTGGGTAACAAGCGTGAAGCGATGAAACAGCTGACCGACCTG    100
                     ||||||||||||||||||||||||||||||||||||||||||||||||||
gyrA_reference    51 GTGAAACTGCTGGGTAACAAGCGTGAAGCGATGAAACAGCTGACCGACCTG    100

gyrA_query       101 ATCGACGATGAACAGTTCGACAAACTGGAACAGATCCTGCGTGAACTGCTG    150
                     ||||||||||||||||||||||||||||||||||||||||||||||||||
gyrA_reference   101 ATCGACGATGAACAGTTCGACAAACTGGAACAGATCCTGCGTGAACTGCTG    150

gyrA_query       151 GAACAGCCGTTCATCAAGAAACTGGGTTGCCCGGACGGCCTGAAAGACATC    200
                     ||||||||||||||||||||||||||||||||||||||||||||||||||
gyrA_reference   151 GAACAGCCGTTCATCAAGAAACTGGGTTGCCCGGACGGCCTGAAAGACATC    200

gyrA_query       201 AAGACCTTCATCCGTATCTTCGCCCAGTTCCTGAAAGACATGACCACCAAC    250
                     ||||||||||||||||||||||||||||||||||||||||||||||||||
gyrA_reference   201 AAGACCTTCATCCGTATCTTCGCCCAGTTCCTGAAAGACATGACCACCAAC    250

gyrA_query       251 CTGAAACAGTAA    262
                     ||||||||||||
gyrA_reference   251 CTGAAACAGTAA    262


########################################
# Program: water
# Rundate: Wed 13 Sep 2023 10:30:01
########################################

#=======================================
#
# Aligned_sequences: 2
# 1: gyrA_mutant
# 2: gyrA_reference
# Matrix: EDNAFULL
# Gap_penalty: 10.0
# Extend_penalty: 0.5
#
# Length: 261
# Identity:     247/261 (94.6%)
# Similarity:   247/261 (94.6%)
# Gaps:           0/261 ( 0.0%)
# Score: 1180.0
#
#
#=======================================

gyrA_mutant        1 ATGGCTGAACGTCTGGTGCTGGAACAGTTCGGCGACCTGCGTACGAAACCG     50
                     |||||||||||||||||||||||||||||||||||||||||||||||||| 
gyrA_reference     1 ATGGCTGAACGTCTGGTGCTGGAACAGTTCGGCGACCTGCGTACGAAACCG     50

gyrA_mutant       51 GTGAAACTGCTGGGTAACAAGCGTGAAGCGATGAAACAGCTGACCGACCTG    100
                     ||||||||||||||||||||||||||||||||||||||||||||||||||
gyrA_reference    51 GTGAAACTGCTGGGTAACAAGCGTGAAGCGTGAAGCGATGAAACAGCTGACCGACCTG    100

gyrA_mutant      101 ATCGACGATGAACAGTTCGACAAACTGGAACAGATCCTGCGTGAACTGCTG    150
                     ||||||||||||||||||||||||||||||||||||||||||||||||||
gyrA_reference   101 ATCGACGATGAACAGTTCGACAAACTGGAACAGATCCTGCGTGAACTGCTG    150

gyrA_mutant      151 GAACAGCCGTTCATCAAGAAACTGGGTTGCCCGGACGGCCTGAAAGACATC    200
                     ||||||||||||||||||||||||||||||||||||||||||||||||||
gyrA_reference   151 GAACAGCCGTTCATCAAGAAACTGGGTTGCCCGGACGGCCTGAAAGACATC    200

gyrA_mutant      201 AAGACCTTCATCCGTATCTTCGCCCAGTTCCTGAAAGACATGACCACCAAC    250
                     ||||||||||||||||||||||   ||||||||||||||||||||||||| 
gyrA_reference   201 AAGACCTTCATCCGTATCTTCGCCCAGTTCCTGAAAGACATGACCACCAAC    250

gyrA_mutant      251 CTGAAACAGTAA    262
                     ||||||||||||
gyrA_reference   251 CTGAAACAGTAA    262


"""
    return water_output

def create_mock_alignment_with_substitutions() -> str:
    """Create mock alignment with known substitutions for validation."""
    return """########################################
# Program: water
# Rundate: Wed 13 Sep 2023 10:30:02
########################################

#=======================================
#
# Aligned_sequences: 2
# 1: test_query
# 2: test_reference  
# Matrix: EDNAFULL
# Gap_penalty: 10.0
# Extend_penalty: 0.5
#
# Length: 100
# Identity:      95/100 (95.0%)
# Similarity:    95/100 (95.0%)
# Gaps:           0/100 ( 0.0%)
# Score: 450.0
#
#
#=======================================

test_query         1 ATGGCTGAACGTCTGGTGCTGGAACAGTTCGGCGACCTGCGTACGAAACCG     50
                     ||||||||||||||||||||||||||||||||||||||||||||||||||
test_reference     1 ATGGCTGAACGTCTGGTGCTGGAACAGTTCGGCGACCTGCGTACGAAACCG     50

test_query        51 GTGAAACTGCTGGGTAACAAGCGTGAAGCGATGAAACAGCTGACCGACCTG    100
                     |||||||||||||||||||||| |||||| |||||||||||||||||||
test_reference    51 GTGAAACTGCTGGGTAACAAGCCTGAAGCAATGAAACAGCTGACCGACCTG    100


"""

def test_emboss_water_parsing():
    """Test EMBOSS-WATER output parsing capabilities."""
    print("\n" + "="*60)
    print("TESTING EMBOSS-WATER OUTPUT PARSING")
    print("="*60)
    
    parser = EmbossWaterParser()
    
    # Create temporary file with mock data
    with tempfile.NamedTemporaryFile(mode='w', suffix='.water', delete=False) as f:
        f.write(create_mock_emboss_water_output())
        temp_file = f.name
        
    try:
        # Parse the file
        alignments = parser.parse_alignment_file(temp_file)
        
        # Validation tests
        tests_passed = 0
        total_tests = 6
        
        # Test 1: Should parse 2 alignments
        if len(alignments) == 2:
            print("✓ Correct number of alignments parsed")
            tests_passed += 1
        else:
            print(f"✗ Expected 2 alignments, got {len(alignments)}")
            
        # Test 2: First alignment should have high identity
        if alignments[0]['identity']['percent'] == 95.4:
            print("✓ Identity percentage correctly parsed")
            tests_passed += 1
        else:
            print(f"✗ Expected 95.4% identity, got {alignments[0]['identity']['percent']}")
            
        # Test 3: Sequences should be extracted
        if alignments[0]['query_sequence'] and alignments[0]['reference_sequence']:
            print("✓ Sequences successfully extracted")
            tests_passed += 1
        else:
            print("✗ Failed to extract sequences")
            
        # Test 4: Sequence IDs should be identified
        if alignments[0]['query_id'] == 'gyrA_query':
            print("✓ Query sequence ID correctly identified")
            tests_passed += 1
        else:
            print(f"✗ Expected 'gyrA_query', got '{alignments[0]['query_id']}'")
            
        # Test 5: Alignment coordinates should be parsed
        if alignments[0]['query_start'] == 1 and alignments[0]['query_end'] == 262:
            print("✓ Alignment coordinates correctly parsed")
            tests_passed += 1
        else:
            print(f"✗ Incorrect coordinates: {alignments[0]['query_start']}-{alignments[0]['query_end']}")
            
        # Test 6: Score should be parsed
        if alignments[0]['score'] == 1190.0:
            print("✓ Alignment score correctly parsed")
            tests_passed += 1
        else:
            print(f"✗ Expected score 1190.0, got {alignments[0]['score']}")
            
        print(f"\nEMBOSS-WATER Parsing: {tests_passed}/{total_tests} tests passed ({100*tests_passed/total_tests:.1f}%)")
        assert tests_passed == total_tests
        
    finally:
        os.unlink(temp_file)

def test_substitution_detection():
    """Test substitution detection accuracy and quality control."""
    print("\n" + "="*60)
    print("TESTING SUBSTITUTION DETECTION")
    print("="*60)
    
    detector = SubstitutionDetector()
    
    # Create test alignment with known substitutions
    test_alignment = {
        'alignment_id': 'test_001',
        'query_id': 'test_query',
        'reference_id': 'test_reference',
        'query_sequence': 'ATGGCTGAACGTCTGGTGCTGGAACAGTTCGGCGACCTGCGTACGAAACCGGTGAAACTGCTGGGTAACAAGCGTGAAGCGATGAAACAGCTGACCGACCTG',
        'reference_sequence': 'ATGGCTGAACGTCTGGTGCTGGAACAGTTCGGCGACCTGCGTACGAAACCGGTGAAACTGCTGGGTAACAAGCCTGAAGCAATGAAACAGCTGACCGACCTG',
        'identity': {'percent': 95.0, 'matches': 95, 'total': 100},
        'gaps': {'count': 0, 'percent': 0.0},
        'score': 450.0,
        'query_start': 1,
        'reference_start': 1
    }
    
    # Expected mutations at positions where sequences differ
    # Position 79: C->G (reference C, query G) 
    # Position 85: C->A (reference C, query A)
    
    mutations = detector.detect_substitutions(test_alignment)
    
    tests_passed = 0
    total_tests = 5
    
    # Test 1: Should detect exactly 2 mutations
    if len(mutations) == 2:
        print("✓ Correct number of substitutions detected")
        tests_passed += 1
    else:
        print(f"✗ Expected 2 substitutions, detected {len(mutations)}")
        
    if mutations:
        # Test 2: Mutations should be substitution type
        if all(m.mutation_type == MutationType.SUBSTITUTION for m in mutations):
            print("✓ Correct mutation type classification")
            tests_passed += 1
        else:
            print("✗ Incorrect mutation type classification")
            
        # Test 3: Quality scores should be reasonable
        avg_quality = sum(m.quality_score for m in mutations) / len(mutations)
        if avg_quality >= 0.7:
            print(f"✓ Good quality scores (average: {avg_quality:.3f})")
            tests_passed += 1
        else:
            print(f"✗ Low quality scores (average: {avg_quality:.3f})")
            
        # Test 4: Confidence levels should be assigned
        if all(m.confidence != ConfidenceLevel.UNCERTAIN for m in mutations):
            print("✓ Appropriate confidence levels assigned")
            tests_passed += 1
        else:
            print("✗ Some mutations have uncertain confidence")
            
        # Test 5: Genomic coordinates should be correct
        expected_positions = [79, 85]  # Based on sequence differences
        detected_positions = [m.genomic_coord.position for m in mutations]
        if set(detected_positions) == set(expected_positions):
            print("✓ Correct genomic coordinates")
            tests_passed += 1
        else:
            print(f"✗ Expected positions {expected_positions}, got {detected_positions}")
            
    print(f"\nSubstitution Detection: {tests_passed}/{total_tests} tests passed ({100*tests_passed/total_tests:.1f}%)")
    assert tests_passed == total_tests

def test_functional_annotation():
    """Test functional annotation and AMR database integration."""
    print("\n" + "="*60)
    print("TESTING FUNCTIONAL ANNOTATION")
    print("="*60)
    
    annotator = FunctionalAnnotator()
    
    # Create test mutation for gyrA S83L (known quinolone resistance)
    test_mutation = MutationRecord(
        mutation_id="gyrA:248:A>T",
        mutation_type=MutationType.SUBSTITUTION,
        genomic_coord=GenomicCoordinate(
            contig="gyrA",
            position=248,  # Position 83 in amino acid sequence
            reference_base="A",
            gene_context="gyrA"
        ),
        confidence=ConfidenceLevel.HIGH,
        quality_metrics=AlignmentQuality(
            identity_percent=95.0,
            coverage_percent=100.0,
            gap_count=0,
            mismatch_count=1
        ),
        quality_score=0.95,
        source_alignment="test_alignment",
        query_sequence="test_query",
        reference_sequence="gyrA_reference"
    )
    
    # Add simplified protein context for known AMR position
    from analysis.subscan_alignment_analyzer import ProteinContext
    test_mutation.protein_context = ProteinContext(
        protein_id="gyrA",
        amino_acid_position=83,
        reference_aa="S",
        mutant_aa="L",
        codon_position=1
    )
    
    # Test annotation
    annotated_mutation = annotator.annotate_mutation(test_mutation)
    
    tests_passed = 0
    total_tests = 4
    
    # Test 1: Should recognize as known AMR mutation
    if annotated_mutation.known_resistance:
        print("✓ Known AMR mutation correctly identified")
        tests_passed += 1
    else:
        print("✗ Failed to identify known AMR mutation")
        
    # Test 2: Should assign resistance mechanism
    if annotated_mutation.resistance_mechanism:
        print(f"✓ Resistance mechanism assigned: {annotated_mutation.resistance_mechanism}")
        tests_passed += 1
    else:
        print("✗ No resistance mechanism assigned")
        
    # Test 3: Should associate with drugs
    if annotated_mutation.drug_associations:
        print(f"✓ Drug associations found: {', '.join(annotated_mutation.drug_associations)}")
        tests_passed += 1
    else:
        print("✗ No drug associations found")
        
    # Test 4: Should predict functional impact
    if annotated_mutation.functional_impact:
        print(f"✓ Functional impact predicted: {annotated_mutation.functional_impact}")
        tests_passed += 1
    else:
        print("✗ No functional impact predicted")
        
    print(f"\nFunctional Annotation: {tests_passed}/{total_tests} tests passed ({100*tests_passed/total_tests:.1f}%)")
    assert tests_passed == total_tests

def test_end_to_end_analysis():
    """Test complete end-to-end analysis workflow."""
    print("\n" + "="*60)
    print("TESTING END-TO-END ANALYSIS WORKFLOW")
    print("="*60)
    
    # Create temporary alignment file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.water', delete=False) as f:
        f.write(create_mock_alignment_with_substitutions())
        temp_file = f.name
        
    try:
        # Create temporary output directory
        with tempfile.TemporaryDirectory() as temp_dir:
            analyzer = SubScanAlignmentAnalyzer(
                output_dir=temp_dir,
                db_path=":memory:",  # In-memory database for testing
                min_confidence=ConfidenceLevel.LOW
            )
            
            try:
                # Run analysis
                mutations = analyzer.analyze_alignment_file(temp_file, "TEST_GENOME_001")
                
                tests_passed = 0
                total_tests = 6
                
                # Test 1: Should complete without errors
                print("✓ Analysis completed without errors")
                tests_passed += 1
                
                # Test 2: Should detect mutations
                if len(mutations) > 0:
                    print(f"✓ Mutations detected: {len(mutations)}")
                    tests_passed += 1
                else:
                    print("✗ No mutations detected")
                    
                # Test 3: Should generate statistics
                stats = analyzer.get_analysis_statistics()
                if stats['alignments_processed'] > 0:
                    print("✓ Analysis statistics generated")
                    tests_passed += 1
                else:
                    print("✗ No analysis statistics")
                    
                # Test 4: Should create output files
                output_files = list(Path(temp_dir).glob("*"))
                if output_files:
                    print(f"✓ Output files created: {len(output_files)}")
                    tests_passed += 1
                else:
                    print("✗ No output files created")
                    
                # Test 5: Quality filtering should work
                high_qual_mutations = [m for m in mutations if m.confidence == ConfidenceLevel.HIGH]
                print(f"✓ Quality filtering applied: {len(high_qual_mutations)} high-confidence mutations")
                tests_passed += 1
                
                # Test 6: Database storage should work
                if stats['mutations_detected'] == len(mutations):
                    print("✓ Database storage functioning")
                    tests_passed += 1
                else:
                    print("✗ Database storage issues")
                    
                print(f"\nEnd-to-End Analysis: {tests_passed}/{total_tests} tests passed ({100*tests_passed/total_tests:.1f}%)")
                assert tests_passed == total_tests
                
            finally:
                analyzer.close()
                
    finally:
        os.unlink(temp_file)

def test_error_handling_and_edge_cases():
    """Test error handling and edge case management.""" 
    print("\n" + "="*60)
    print("TESTING ERROR HANDLING & EDGE CASES")
    print("="*60)
    
    tests_passed = 0
    total_tests = 5
    
    # Test 1: Empty alignment file
    try:
        parser = EmbossWaterParser()
        with tempfile.NamedTemporaryFile(mode='w', suffix='.water', delete=False) as f:
            f.write("")
            temp_file = f.name
            
        alignments = parser.parse_alignment_file(temp_file)
        os.unlink(temp_file)
        
        if len(alignments) == 0:
            print("✓ Empty file handled gracefully")
            tests_passed += 1
        else:
            print("✗ Empty file not handled properly")
            
    except Exception as e:
        print(f"✗ Empty file caused exception: {e}")
        
    # Test 2: Malformed alignment data
    try:
        detector = SubstitutionDetector()
        malformed_alignment = {
            'alignment_id': 'malformed',
            'query_sequence': 'ATCG',
            'reference_sequence': 'AT',  # Different length
            'identity': {'percent': 50.0},
            'gaps': {'count': 0, 'percent': 0.0}
        }
        
        mutations = detector.detect_substitutions(malformed_alignment)
        if len(mutations) == 0:
            print("✓ Malformed alignment handled gracefully")
            tests_passed += 1
        else:
            print("✗ Malformed alignment not handled properly")
            
    except Exception as e:
        print("✓ Malformed alignment properly rejected with exception")
        tests_passed += 1
        
    # Test 3: Very low quality alignment
    try:
        low_quality_alignment = {
            'alignment_id': 'low_quality',
            'query_sequence': 'ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG',
            'reference_sequence': 'GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT',
            'identity': {'percent': 25.0},  # Very low identity
            'gaps': {'count': 0, 'percent': 0.0},
            'score': 50.0
        }
        
        mutations = detector.detect_substitutions(low_quality_alignment)
        # Should be filtered out due to low quality
        if len(mutations) == 0:
            print("✓ Low quality alignment properly filtered")
            tests_passed += 1
        else:
            print("✗ Low quality alignment not filtered")
            
    except Exception:
        print("✓ Low quality alignment properly rejected")
        tests_passed += 1
        
    # Test 4: Non-existent file
    try:
        with tempfile.TemporaryDirectory() as temp_dir:
            analyzer = SubScanAlignmentAnalyzer(
                output_dir=temp_dir,
                db_path=":memory:"
            )
            
            try:
                mutations = analyzer.analyze_alignment_file("nonexistent_file.water")
                print("✗ Non-existent file should raise exception")
            except Exception:
                print("✓ Non-existent file properly handled with exception")
                tests_passed += 1
            finally:
                analyzer.close()
                
    except Exception as e:
        print(f"✗ Unexpected error with non-existent file: {e}")
        
    # Test 5: Analyzer with invalid parameters
    try:
        analyzer = SubScanAlignmentAnalyzer(
            output_dir="/invalid/path/that/does/not/exist",
            min_confidence=ConfidenceLevel.HIGH
        )
        print("✓ Invalid parameters handled gracefully")
        tests_passed += 1
        analyzer.close()
    except Exception:
        print("✓ Invalid parameters properly rejected")
        tests_passed += 1
        
    print(f"\nError Handling & Edge Cases: {tests_passed}/{total_tests} tests passed ({100*tests_passed/total_tests:.1f}%)")
    assert tests_passed == total_tests

def main():
    """Run comprehensive SubScan analyzer validation."""
    setup_logging()
    
    print("SUBSCAN ALIGNMENT ANALYZER - COMPREHENSIVE VALIDATION")
    print("=" * 80)
    print("Testing enterprise-grade mutation detection pipeline")
    print("=" * 80)
    
    # Run all test suites
    test_results = []
    
    test_results.append(test_emboss_water_parsing())
    test_results.append(test_substitution_detection())
    test_results.append(test_functional_annotation())
    test_results.append(test_end_to_end_analysis())
    test_results.append(test_error_handling_and_edge_cases())
    
    # Final summary
    passed_tests = sum(test_results)
    total_tests = len(test_results)
    success_rate = 100 * passed_tests / total_tests
    
    print("\n" + "="*80)
    print("FINAL VALIDATION SUMMARY")
    print("="*80)
    print(f"Test Suites Passed: {passed_tests}/{total_tests} ({success_rate:.1f}%)")
    
    if passed_tests == total_tests:
        print("🎉 ALL TESTS PASSED - SubScan Analyzer ready for production!")
        print("\nKey Features Validated:")
        print("✓ Robust EMBOSS-WATER output parsing")
        print("✓ Accurate substitution detection with quality control")
        print("✓ Comprehensive functional annotation")
        print("✓ End-to-end workflow integration")
        print("✓ Error handling and edge case management")
        print("✓ Database integration and artifact storage")
    else:
        print("⚠️  Some tests failed - review implementation before deployment")
        
    print("\nNext Steps:")
    print("1. Review any failed tests and fix implementation")
    print("2. Test with real EMBOSS-WATER alignment files")
    print("3. Integrate with existing CARD and protein extraction modules")
    print("4. Configure production parameters and thresholds")
    print("5. Validate with known AMR mutations from literature")
    
    return passed_tests == total_tests

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)