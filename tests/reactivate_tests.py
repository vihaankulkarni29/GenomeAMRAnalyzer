#!/usr/bin/env python3
"""
Test Reactivation Script for GenomeAMRAnalyzer
==============================================

This script reactivates disabled test files that are relevant to the new Abricate workflow.
It copies content from disabled_test_*.py files to their corresponding test_*.py files,
with updates for Abricate compatibility where needed.

Priority Tests for Reactivation:
1. FASTA AA Extractor tests - Core pipeline functionality
2. Workflow/pipeline integration tests  
3. Database repository tests
4. Genome harvester tests
5. URL discovery tests

Skip:
- Pure RGI tests (replaced by Abricate)
- CARD-specific tests (replaced by Abricate databases)
"""

import os
import shutil
from pathlib import Path


def get_file_size(filepath):
    """Get file size, return 0 if file doesn't exist."""
    try:
        return os.path.getsize(filepath)
    except OSError:
        return 0


def should_reactivate_test(test_name):
    """Determine if a test should be reactivated based on relevance to Abricate workflow."""
    # Core pipeline components - High priority
    high_priority = [
        "fasta_aa_extractor_simplified",
        "fasta_aa_extractor_robust", 
        "workflow_hardcore",
        "db_repository",
        "genome_harvester",
        "genome_harvester_simple",
        "genbank_genome_harvester"
    ]
    
    # Integration and infrastructure - Medium priority
    medium_priority = [
        "production_wildtype_aligner_integration",
        "url_genome_discovery",
        "url_discovery_hardcore",
        "downloader_hardcore",
        "performance_hardcore",
        "html_report_generator_hardcore"
    ]
    
    # Analysis components - Medium priority
    analysis_priority = [
        "cooccurrence_analyzer",
        "cooccurrence_hardcore",
        "subscan_hardcore",
        "mic_metadata_harvester_corrected",
        "mic_metadata_harvester_robust"
    ]
    
    # Skip RGI/CARD specific tests
    skip_rgi = [
        "card_rgi_runner",
        "wildtype_aligner_comprehensive"  # Likely RGI-specific
    ]
    
    base_name = test_name.replace("disabled_test_", "").replace("test_", "").replace(".py", "")
    
    if base_name in skip_rgi:
        return False, "RGI/CARD specific - replaced by Abricate"
    elif base_name in high_priority:
        return True, "High priority - core pipeline"
    elif base_name in medium_priority:
        return True, "Medium priority - integration/infrastructure"
    elif base_name in analysis_priority:
        return True, "Medium priority - analysis components"
    else:
        return False, "Low priority or unknown"


def reactivate_test_file(disabled_path, active_path):
    """Copy content from disabled test to active test file."""
    try:
        # Read content from disabled file
        with open(disabled_path, 'r', encoding='utf-8') as f:
            content = f.read()
        
        # Write to active file
        with open(active_path, 'w', encoding='utf-8') as f:
            f.write(content)
        
        return True, f"Successfully copied {disabled_path.name} -> {active_path.name}"
    except Exception as e:
        return False, f"Error copying {disabled_path.name}: {e}"


def analyze_test_files():
    """Analyze all test files and provide reactivation recommendations."""
    tests_dir = Path(__file__).parent
    disabled_files = list(tests_dir.glob("disabled_test_*.py"))
    
    results = {
        'to_reactivate': [],
        'to_skip': [],
        'already_active': [],
        'errors': []
    }
    
    for disabled_file in disabled_files:
        # Get corresponding active filename
        active_name = disabled_file.name.replace("disabled_", "")
        active_file = tests_dir / active_name
        
        # Check sizes
        disabled_size = get_file_size(disabled_file)
        active_size = get_file_size(active_file)
        
        # Determine if should reactivate
        should_reactivate, reason = should_reactivate_test(disabled_file.name)
        
        status = {
            'disabled_file': disabled_file,
            'active_file': active_file,
            'disabled_size': disabled_size,
            'active_size': active_size,
            'reason': reason
        }
        
        if not should_reactivate:
            results['to_skip'].append(status)
        elif active_size > 0 and active_size >= disabled_size * 0.8:  # Active file has substantial content
            results['already_active'].append(status)
        elif disabled_size > 0:
            results['to_reactivate'].append(status)
        else:
            status['reason'] = "Disabled file is empty"
            results['errors'].append(status)
    
    return results


def generate_reactivation_commands():
    """Generate shell commands for reactivating tests."""
    analysis = analyze_test_files()
    commands = []
    
    print("=" * 60)
    print("TEST REACTIVATION ANALYSIS")
    print("=" * 60)
    
    print(f"\nHIGH PRIORITY - TO REACTIVATE ({len(analysis['to_reactivate'])} files):")
    for item in analysis['to_reactivate']:
        disabled_file = item['disabled_file']
        active_file = item['active_file']
        print(f"  ✓ {disabled_file.name} -> {active_file.name}")
        print(f"    Size: {item['disabled_size']} -> {item['active_size']} bytes")
        print(f"    Reason: {item['reason']}")
        
        # Generate copy command
        cmd = f'Copy-Item "{disabled_file}" "{active_file}" -Force'
        commands.append(cmd)
    
    print(f"\nALREADY ACTIVE ({len(analysis['already_active'])} files):")
    for item in analysis['already_active']:
        print(f"  ⚡ {item['active_file'].name} (Size: {item['active_size']} bytes)")
        print(f"     Reason: {item['reason']}")
    
    print(f"\nSKIPPING ({len(analysis['to_skip'])} files):")
    for item in analysis['to_skip']:
        print(f"  ⏭️  {item['disabled_file'].name}")
        print(f"     Reason: {item['reason']}")
    
    if analysis['errors']:
        print(f"\nERRORS ({len(analysis['errors'])} files):")
        for item in analysis['errors']:
            print(f"  ❌ {item['disabled_file'].name}")
            print(f"     Reason: {item['reason']}")
    
    print("\n" + "=" * 60)
    print("POWERSHELL COMMANDS TO REACTIVATE TESTS:")
    print("=" * 60)
    
    if commands:
        print("# Navigate to tests directory")
        print('cd "c:\\Users\\Vihaan\\Documents\\GenomeAMRAnalyzer\\tests"')
        print()
        for cmd in commands:
            print(cmd)
        print()
        print("# Verify reactivation")
        print("Get-ChildItem test_*.py | Where-Object { $_.Length -gt 0 } | Select-Object Name, Length | Sort-Object Name")
    else:
        print("# No files need reactivation - all relevant tests are already active!")
    
    return commands


if __name__ == "__main__":
    generate_reactivation_commands()