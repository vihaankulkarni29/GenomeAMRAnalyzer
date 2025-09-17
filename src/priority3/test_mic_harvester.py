"""
MIC Metadata Harvester Test & Validation
========================================

Comprehensive testing and validation for the MIC metadata harvester.
Demonstrates enterprise-grade reliability and data quality assurance.
"""

import sys
import os
import logging
from pathlib import Path
from typing import List, Dict, Any

# Add parent directory to path for imports
sys.path.append(str(Path(__file__).parent.parent))

from metadata.mic_metadata_harvester import (
    NCBIMICHarvester, AntibioticStandardizer, MICNormalizer,
    MICRecord, BioSampleMetadata, MICUnit, ResistanceProfile
)

def setup_logging():
    """Setup comprehensive logging for validation."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler('mic_harvester_validation.log'),
            logging.StreamHandler()
        ]
    )

def test_antibiotic_standardization():
    """Test antibiotic name standardization with edge cases."""
    print("\n" + "="*60)
    print("TESTING ANTIBIOTIC STANDARDIZATION")
    print("="*60)
    
    standardizer = AntibioticStandardizer()
    
    test_cases = [
        # Standard names
        ("ampicillin", "ampicillin", 1.0),
        ("ciprofloxacin", "ciprofloxacin", 1.0),
        
        # Common abbreviations
        ("amp", "ampicillin", 0.9),
        ("cipro", "ciprofloxacin", 0.9),
        ("gent", "gentamicin", 0.9),
        
        # WHONET codes
        ("AMP", "ampicillin", 1.0),
        ("CIP", "ciprofloxacin", 1.0),
        
        # Variations and misspellings
        ("ampic", "ampicillin", 0.9),
        ("ciproflox", "ciprofloxacin", 0.8),
        
        # Complex names
        ("trimethoprim-sulfamethoxazole", "trimethoprim-sulfamethoxazole", 0.9),
        
        # Unknown antibiotics
        ("unknowndrugXYZ", "unknowndrugXYZ", 0.3),
    ]
    
    passed = 0
    total = len(test_cases)
    
    for input_name, expected_output, min_confidence in test_cases:
        standardized, confidence = standardizer.standardize(input_name)
        
        success = (standardized == expected_output and confidence >= min_confidence)
        status = "✓ PASS" if success else "✗ FAIL"
        
        print(f"{status} | {input_name:25} → {standardized:25} | confidence: {confidence:.2f}")
        
        if success:
            passed += 1
        else:
            print(f"      Expected: {expected_output}, Got: {standardized}")
            
    print(f"\nAntibiotic Standardization: {passed}/{total} tests passed ({100*passed/total:.1f}%)")
    assert passed == total

def test_mic_normalization():
    """Test MIC value normalization with various formats."""
    print("\n" + "="*60)
    print("TESTING MIC VALUE NORMALIZATION")
    print("="*60)
    
    normalizer = MICNormalizer()
    
    test_cases = [
        # Standard cases
        ("16", "mg/L", "ampicillin", 16.0, MICUnit.MG_L, 1.0),
        ("8", "μg/mL", "ciprofloxacin", 8.0, MICUnit.MG_L, 1.0),
        
        # With operators
        ("<=4", "mg/L", "gentamicin", 4.0, MICUnit.MG_L, 0.8),
        (">32", "mg/L", "vancomycin", 32.0, MICUnit.MG_L, 0.8),
        ("≤2", "mg/L", "tetracycline", 2.0, MICUnit.MG_L, 0.8),
        
        # Ranges
        ("2-4", "mg/L", "ampicillin", 2.83, MICUnit.MG_L, 0.7),  # sqrt(2*4)
        ("0.5/1.0", "mg/L", "ciprofloxacin", 0.71, MICUnit.MG_L, 0.7),  # sqrt(0.5*1.0)
        
        # Decimal values
        ("0.25", "mg/L", "gentamicin", 0.25, MICUnit.MG_L, 1.0),
        ("1.5", "μg/mL", "ciprofloxacin", 1.5, MICUnit.MG_L, 1.0),
        
        # Edge cases
        ("", "mg/L", "ampicillin", None, MICUnit.UNKNOWN, 0.0),
        ("invalid", "mg/L", "ampicillin", None, MICUnit.UNKNOWN, 0.0),
    ]
    
    passed = 0
    total = len(test_cases)
    
    for value_str, unit_str, antibiotic, expected_value, expected_unit, min_quality in test_cases:
        result_value, result_unit, quality = normalizer.normalize_mic_value(value_str, unit_str, antibiotic)
        
        # Check if result matches expectation (with tolerance for floats)
        value_match = (
            (expected_value is None and result_value is None) or
            (expected_value is not None and result_value is not None and abs(result_value - expected_value) < 0.1)
        )
        
        unit_match = result_unit == expected_unit
        quality_match = quality >= min_quality
        
        success = value_match and unit_match and quality_match
        status = "✓ PASS" if success else "✗ FAIL"
        
        print(f"{status} | {value_str:10} {unit_str:8} → {result_value:8} {result_unit.value:8} | Q: {quality:.2f}")
        
        if success:
            passed += 1
        else:
            print(f"      Expected: {expected_value} {expected_unit.value}, Got: {result_value} {result_unit.value}")
            
    print(f"\nMIC Normalization: {passed}/{total} tests passed ({100*passed/total:.1f}%)")
    assert passed == total

def test_mock_mic_harvesting():
    """Test MIC harvesting with mock data."""
    print("\n" + "="*60)
    print("TESTING MOCK MIC HARVESTING")
    print("="*60)
    
    # Test accessions
    test_accessions = [
        "GCF_000001405.39",  # Human reference
        "GCF_000005825.2",   # E. coli K-12
        "GCF_000009605.1",   # S. aureus
        "GCF_000195955.2",   # Mycobacterium tuberculosis
        "GCF_000006825.1",   # Pseudomonas aeruginosa
    ]
    
    # Initialize harvester in mock mode
    harvester = NCBIMICHarvester(
        db_path="test_priority3.db",
        mock_mode=True
    )
    
    try:
        # Harvest MIC data
        print(f"Harvesting MIC data for {len(test_accessions)} genomes...")
        results = harvester.harvest_mic_data(test_accessions)
        
        print(f"\nResults Summary:")
        print(f"- Genomes processed: {len(test_accessions)}")
        print(f"- Genomes with MIC data: {len(results)}")
        
        # Analyze results
        total_mic_records = 0
        antibiotics_found = set()
        quality_scores = []
        
        for accession, metadata in results.items():
            print(f"\n{accession}:")
            print(f"  - Organism: {metadata.organism}")
            print(f"  - Strain: {metadata.strain}")
            print(f"  - BioSample: {metadata.biosample_id}")
            print(f"  - MIC Records: {len(metadata.mic_records)}")
            
            total_mic_records += len(metadata.mic_records)
            
            for mic_record in metadata.mic_records:
                antibiotics_found.add(mic_record.antibiotic_standardized)
                quality_scores.append(mic_record.quality_score)
                
                print(f"    • {mic_record.antibiotic_standardized}: {mic_record.mic_value} {mic_record.mic_unit.value} ({mic_record.resistance_profile.value})")
                
        # Summary statistics
        avg_quality = sum(quality_scores) / len(quality_scores) if quality_scores else 0
        
        print(f"\n" + "-"*50)
        print(f"MOCK HARVESTING SUMMARY:")
        print(f"- Total MIC records: {total_mic_records}")
        print(f"- Unique antibiotics: {len(antibiotics_found)}")
        print(f"- Average quality score: {avg_quality:.3f}")
        print(f"- Antibiotics found: {', '.join(sorted(antibiotics_found))}")
        
        # Validation checks
        checks_passed = 0
        total_checks = 4
        
        # Check 1: All accessions should have results
        if len(results) == len(test_accessions):
            print("✓ All accessions have results")
            checks_passed += 1
        else:
            print("✗ Missing results for some accessions")
            
        # Check 2: All results should have MIC records
        all_have_mics = all(len(metadata.mic_records) > 0 for metadata in results.values())
        if all_have_mics:
            print("✓ All results have MIC records")
            checks_passed += 1
        else:
            print("✗ Some results missing MIC records")
            
        # Check 3: Quality scores should be reasonable
        if avg_quality >= 0.8:
            print("✓ Good average quality score")
            checks_passed += 1
        else:
            print("✗ Low average quality score")
            
        # Check 4: Should find multiple antibiotics
        if len(antibiotics_found) >= 3:
            print("✓ Multiple antibiotics found")
            checks_passed += 1
        else:
            print("✗ Too few antibiotics found")
            
        print(f"\nMock Harvesting Validation: {checks_passed}/{total_checks} checks passed")
        assert checks_passed == total_checks
        
    finally:
        harvester.close()

def test_data_quality_assessment():
    """Test data quality assessment features."""
    print("\n" + "="*60)
    print("TESTING DATA QUALITY ASSESSMENT")
    print("="*60)
    
    # Create test metadata with various quality issues
    test_cases = [
        {
            'name': 'High Quality',
            'attributes': {
                'strain': 'ATCC 25922',
                'isolation_source': 'blood',
                'collection_date': '2023-01-15',
                'geographic_location': 'USA: California'
            },
            'mic_count': 6,
            'avg_quality': 0.95,
            'expected_flags': []
        },
        {
            'name': 'Missing Strain',
            'attributes': {
                'isolation_source': 'urine',
                'collection_date': '2023-02-01'
            },
            'mic_count': 4,
            'avg_quality': 0.85,
            'expected_flags': ['missing_strain_info']
        },
        {
            'name': 'Limited MIC Data',
            'attributes': {
                'strain': 'clinical_isolate_1',
                'isolation_source': 'wound'
            },
            'mic_count': 2,
            'avg_quality': 0.9,
            'expected_flags': ['missing_collection_date', 'limited_mic_data']
        },
        {
            'name': 'Poor Quality MICs',
            'attributes': {
                'strain': 'unknown',
                'collection_date': '2023-03-01'
            },
            'mic_count': 5,
            'avg_quality': 0.4,
            'expected_flags': ['missing_isolation_source', 'poor_mic_quality']
        }
    ]
    
    harvester = NCBIMICHarvester(mock_mode=True)
    
    try:
        passed = 0
        total = len(test_cases)
        
        for test_case in test_cases:
            # Create mock MIC records with specified quality
            mic_records = []
            for i in range(test_case['mic_count']):
                mic_record = MICRecord(
                    accession="TEST_001",
                    antibiotic="ampicillin",
                    antibiotic_standardized="ampicillin", 
                    mic_value=16.0,
                    mic_unit=MICUnit.MG_L,
                    mic_unit_original="mg/L",
                    resistance_profile=ResistanceProfile.RESISTANT,
                    quality_score=test_case['avg_quality']
                )
                mic_records.append(mic_record)
                
            # Assess quality
            flags = harvester._assess_metadata_quality(test_case['attributes'], mic_records)
            
            # Check if expected flags are present
            expected_set = set(test_case['expected_flags'])
            actual_set = set(flags)
            
            success = expected_set.issubset(actual_set)
            status = "✓ PASS" if success else "✗ FAIL"
            
            print(f"{status} | {test_case['name']:20} | Flags: {', '.join(flags) if flags else 'none'}")
            
            if success:
                passed += 1
            else:
                print(f"      Expected flags: {test_case['expected_flags']}")
                print(f"      Actual flags: {flags}")
                
        print(f"\nData Quality Assessment: {passed}/{total} tests passed ({100*passed/total:.1f}%)")
        assert passed == total
        
    finally:
        harvester.close()

def main():
    """Run comprehensive MIC harvester validation."""
    setup_logging()
    
    print("MIC METADATA HARVESTER - COMPREHENSIVE VALIDATION")
    print("=" * 80)
    print("Testing enterprise-grade AMR phenotype data collection")
    print("=" * 80)
    
    # Run all test suites
    test_results = []
    
    test_results.append(test_antibiotic_standardization())
    test_results.append(test_mic_normalization()) 
    test_results.append(test_mock_mic_harvesting())
    test_results.append(test_data_quality_assessment())
    
    # Final summary
    passed_tests = sum(test_results)
    total_tests = len(test_results)
    success_rate = 100 * passed_tests / total_tests
    
    print("\n" + "="*80)
    print("FINAL VALIDATION SUMMARY")
    print("="*80)
    print(f"Test Suites Passed: {passed_tests}/{total_tests} ({success_rate:.1f}%)")
    
    if passed_tests == total_tests:
        print("🎉 ALL TESTS PASSED - MIC Harvester ready for production!")
        print("\nKey Features Validated:")
        print("✓ Robust antibiotic name standardization")
        print("✓ Accurate MIC value normalization")
        print("✓ Comprehensive data quality assessment")
        print("✓ Mock data generation for testing")
        print("✓ Error handling and edge cases")
    else:
        print("⚠️  Some tests failed - review implementation before deployment")
        
    print("\nNext Steps:")
    print("1. Review any failed tests and fix implementation")
    print("2. Test with real NCBI data (requires API access)")
    print("3. Integrate with GenBank harvester for complete workflow")
    print("4. Configure production database and logging")
    
    return passed_tests == total_tests

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)