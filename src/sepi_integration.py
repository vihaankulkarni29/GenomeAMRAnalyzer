#!/usr/bin/env python3
"""
SEPI 2.0 Integration for Production WildType Aligner
Seamless integration between SEPI 2.0 and ProductionWildTypeAligner

This module provides:
1. Enhanced SEPIReferenceManager with configuration-driven operation
2. Automated SEPI configuration generation for target genes
3. Robust caching and fallback mechanisms
4. Production-grade error handling and logging

Author: GenomeAMRAnalyzer Pipeline  
Version: 2.0 - Production SEPI Integration
"""

import os
import sys
import time
import asyncio
import logging
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any, Union
import yaml
import json

# Import our configuration manager
from sepi_configuration_manager import SEPIConfigurationManager, GeneTargetConfig, OrganismPreferences


class EnhancedSEPIReferenceManager:
    """
    Enhanced SEPI reference manager with configuration-driven operation
    Integrates with SEPIConfigurationManager for robust reference fetching
    """
    
    def __init__(self, config_manager: SEPIConfigurationManager):
        """Initialize enhanced SEPI reference manager"""
        self.config_manager = config_manager
        self.logger = logging.getLogger('EnhancedSEPIReferenceManager')
        
        # Validate configuration
        if not all([config_manager.environment_config, config_manager.ncbi_config]):
            raise ValueError("Configuration manager not properly initialized")
        
        # Configuration shortcuts
        self.env_config = config_manager.environment_config
        self.ncbi_config = config_manager.ncbi_config
        self.quality_config = config_manager.quality_config
        
        # Cache and state
        self.reference_cache: Dict[str, Dict[str, Any]] = {}
        self.failed_fetches: Dict[str, str] = {}
        self.sepi_availability = self._check_sepi_availability()
        
        # Statistics
        self.stats = {
            'cache_hits': 0,
            'sepi_fetches': 0,
            'sepi_failures': 0,
            'local_fallbacks': 0
        }
        
        self.logger.info("Enhanced SEPI Reference Manager initialized")
    
    def _check_sepi_availability(self) -> bool:
        """Check if SEPI 2.0 is available and functional"""
        try:
            if not self.env_config.sepi_script_path_obj.exists():
                self.logger.warning(f"SEPI script not found: {self.env_config.sepi_script_path}")
                return False
            
            # Test SEPI with help command
            result = subprocess.run(
                [sys.executable, self.env_config.sepi_script_path, "--help"],
                capture_output=True, text=True, timeout=30
            )
            
            if result.returncode == 0:
                self.logger.info("SEPI 2.0 availability confirmed")
                return True
            else:
                self.logger.warning(f"SEPI test failed: {result.stderr}")
                return False
                
        except Exception as e:
            self.logger.warning(f"SEPI availability check failed: {e}")
            return False
    
    def _generate_cache_key(self, gene_name: str, organism: str) -> str:
        """Generate cache key for reference sequences"""
        import hashlib
        key_string = f"{gene_name}_{organism}".lower().replace(" ", "_")
        return hashlib.sha256(key_string.encode()).hexdigest()[:16]
    
    async def fetch_reference_with_sepi(self, gene_name: str, organism: str) -> Optional[Dict[str, Any]]:
        """
        Fetch reference sequence using SEPI 2.0 with configuration-driven approach
        """
        
        if not self.sepi_availability:
            self.logger.warning("SEPI not available, cannot fetch references")
            return None
        
        cache_key = self._generate_cache_key(gene_name, organism)
        
        # Check if we've already failed this fetch
        if cache_key in self.failed_fetches:
            self.logger.info(f"Skipping previously failed fetch: {gene_name} from {organism}")
            return None
        
        try:
            self.logger.info(f"Attempting SEPI fetch: {gene_name} from {organism}")
            
            # Generate SEPI configuration file
            sepi_config_file = self.config_manager.generate_sepi_config_file(gene_name, organism)
            
            # Create temporary output directory
            current_timestamp = int(time.time())
            temp_output = self.env_config.temp_directory_obj / f"sepi_fetch_{gene_name}_{current_timestamp}"
            temp_output.mkdir(exist_ok=True)
            
            # Construct SEPI command using configuration file
            sepi_cmd = [
                sys.executable, self.env_config.sepi_script_path,
                "--config", sepi_config_file,
                "--output", str(temp_output)
            ]
            
            # Execute SEPI command with timeout
            self.logger.info(f"Executing SEPI command: {' '.join(sepi_cmd)}")
            
            process = await asyncio.create_subprocess_exec(
                *sepi_cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
                cwd=temp_output
            )
            
            try:
                stdout, stderr = await asyncio.wait_for(
                    process.communicate(), 
                    timeout=self.env_config.timeout_seconds
                )
            except asyncio.TimeoutError:
                process.kill()
                await process.wait()
                self.logger.error(f"SEPI fetch timed out for {gene_name} from {organism}")
                self.failed_fetches[cache_key] = "Timeout"
                return None
            
            if process.returncode == 0:
                # Parse SEPI output
                fasta_files = list(temp_output.glob("**/*.fasta"))
                json_files = list(temp_output.glob("**/*.json"))
                
                if fasta_files:
                    reference_data = await self._parse_sepi_output(
                        fasta_files[0], json_files[0] if json_files else None,
                        gene_name, organism, sepi_config_file
                    )
                    
                    if reference_data:
                        # Cache the reference
                        self.reference_cache[cache_key] = reference_data
                        self.stats['sepi_fetches'] += 1
                        
                        # Move reference to permanent cache
                        cache_filename = f"{gene_name}_{organism.replace(' ', '_')}_reference.fasta"
                        cache_path = self.env_config.cache_directory_obj / cache_filename
                        
                        import shutil
                        shutil.copy2(fasta_files[0], cache_path)
                        reference_data['cached_file_path'] = str(cache_path)
                        
                        self.logger.info(f"Successfully fetched and cached {gene_name} from {organism}")
                        
                        # Cleanup temporary files
                        import shutil
                        shutil.rmtree(temp_output, ignore_errors=True)
                        
                        return reference_data
                else:
                    self.logger.warning(f"No FASTA files found in SEPI output for {gene_name}")
            else:
                error_msg = stderr.decode() if stderr else "Unknown SEPI error"
                self.logger.error(f"SEPI fetch failed for {gene_name} from {organism}: {error_msg}")
                self.failed_fetches[cache_key] = error_msg
                
            # Cleanup on failure
            import shutil
            shutil.rmtree(temp_output, ignore_errors=True)
            
            # Clean up config file
            try:
                os.unlink(sepi_config_file)
            except:
                pass
            
        except Exception as e:
            self.logger.error(f"SEPI fetch exception for {gene_name} from {organism}: {e}")
            self.failed_fetches[cache_key] = str(e)
        
        self.stats['sepi_failures'] += 1
        return None
    
    async def _parse_sepi_output(self, fasta_file: Path, json_file: Optional[Path],
                               gene_name: str, organism: str, config_file: str) -> Optional[Dict[str, Any]]:
        """Parse SEPI output files and extract reference data"""
        
        try:
            # Parse FASTA file
            try:
                from Bio import SeqIO
                BIOPYTHON_AVAILABLE = True
            except ImportError:
                BIOPYTHON_AVAILABLE = False
                self.logger.warning("BioPython not available for FASTA parsing")
                return None
            
            if BIOPYTHON_AVAILABLE:
                records = list(SeqIO.parse(fasta_file, "fasta"))
                if not records:
                    return None
                
                # Use the first record (best match)
                record = records[0]
                sequence = str(record.seq)
                
                # Parse metadata from JSON if available
                metadata = {}
                if json_file and json_file.exists():
                    try:
                        with open(json_file, 'r') as f:
                            metadata = json.load(f)
                    except Exception as e:
                        self.logger.warning(f"Failed to parse JSON metadata: {e}")
                
                # Calculate sequence checksum
                import hashlib
                checksum = hashlib.sha256(sequence.encode()).hexdigest()[:16]
                
                reference_data = {
                    'gene_name': gene_name,
                    'organism': organism,
                    'accession': record.id,
                    'sequence': sequence,
                    'sequence_length': len(sequence),
                    'source': 'SEPI_2.0',
                    'fetch_timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
                    'sequence_checksum': checksum,
                    'sepi_metadata': {
                        'config_file': config_file,
                        'fasta_file': str(fasta_file),
                        'json_file': str(json_file) if json_file else None,
                        'additional_metadata': metadata,
                        'organism_queried': organism
                    },
                    'quality_score': 1.0,  # High quality from SEPI
                    'reference_file_path': str(fasta_file)
                }
                
                # Validate against quality configuration
                if self.quality_config:
                    if not self._validate_sequence_quality(reference_data):
                        self.logger.warning(f"Reference sequence failed quality validation: {gene_name}")
                        return None
                
                return reference_data
            
        except Exception as e:
            self.logger.error(f"Failed to parse SEPI output: {e}")
            return None
        
        return None
    
    def _validate_sequence_quality(self, reference_data: Dict[str, Any]) -> bool:
        """Validate reference sequence against quality configuration"""
        
        if not self.quality_config:
            return True
        
        sequence = reference_data['sequence']
        sequence_length = len(sequence)
        
        # Length validation
        if sequence_length < self.quality_config.min_sequence_length:
            self.logger.warning(f"Sequence too short: {sequence_length} < {self.quality_config.min_sequence_length}")
            return False
        
        if sequence_length > self.quality_config.max_sequence_length:
            self.logger.warning(f"Sequence too long: {sequence_length} > {self.quality_config.max_sequence_length}")
            return False
        
        # Protein sequence validation
        if self.quality_config.validate_protein_sequence:
            valid_aa = set('ACDEFGHIKLMNPQRSTVWY*')
            invalid_chars = set(sequence.upper()) - valid_aa
            if invalid_chars:
                self.logger.warning(f"Invalid amino acid characters found: {invalid_chars}")
                return False
        
        # Check for partial sequences
        if self.quality_config.exclude_partial_sequences:
            if sequence.upper().startswith('M') == False:
                self.logger.warning("Sequence does not start with methionine (potentially partial)")
                return False
        
        return True
    
    async def get_reference_sequence(self, gene_name: str, query_organism: Optional[str] = None) -> Optional[Dict[str, Any]]:
        """
        Get reference sequence with intelligent organism selection and caching
        """
        
        # Get gene configuration for organism preferences
        gene_config = self.config_manager.get_gene_config(gene_name)
        
        # Determine target organisms with configuration-driven preferences
        target_organisms = []
        
        # Priority 1: Gene-specific organism preferences
        if gene_config.organism_preferences:
            target_organisms.extend(gene_config.organism_preferences)
        
        # Priority 2: Query organism if provided
        if query_organism and query_organism not in target_organisms:
            target_organisms.insert(0, query_organism)
        
        # Priority 3: General organism preferences from configuration
        general_organisms = ['Escherichia coli K-12 MG1655', 'Pseudomonas aeruginosa PAO1', 
                           'Klebsiella pneumoniae', 'Acinetobacter baumannii']
        for org in general_organisms:
            if org not in target_organisms:
                target_organisms.append(org)
        
        # Try each organism in order
        for organism in target_organisms:
            cache_key = self._generate_cache_key(gene_name, organism)
            
            # Check cache first
            if cache_key in self.reference_cache:
                self.logger.info(f"Using cached reference for {gene_name} from {organism}")
                self.stats['cache_hits'] += 1
                return self.reference_cache[cache_key]
            
            # Try SEPI fetch
            reference = await self.fetch_reference_with_sepi(gene_name, organism)
            if reference:
                return reference
        
        # Fallback: Check local reference files
        local_reference = self._check_local_references(gene_name, target_organisms)
        if local_reference:
            self.stats['local_fallbacks'] += 1
            return local_reference
        
        self.logger.error(f"Failed to obtain reference for {gene_name} from any organism")
        return None
    
    def _check_local_references(self, gene_name: str, organisms: List[str]) -> Optional[Dict[str, Any]]:
        """Check for local reference files as fallback"""
        
        # Check common reference directories
        reference_dirs = [
            self.env_config.cache_directory_obj,
            Path("references"),
            Path("../references"),
            Path("../MetaDataHarvester/references")
        ]
        
        for ref_dir in reference_dirs:
            if not ref_dir.exists():
                continue
                
            # Look for gene-specific reference files
            patterns = [
                f"{gene_name}.fasta",
                f"{gene_name}.faa", 
                f"{gene_name}_reference.fasta",
                f"*{gene_name}*.fasta"
            ]
            
            for pattern in patterns:
                matches = list(ref_dir.glob(pattern))
                if matches:
                    ref_file = matches[0]
                    
                    try:
                        from Bio import SeqIO
                        records = list(SeqIO.parse(ref_file, "fasta"))
                        if records:
                            record = records[0]
                            sequence = str(record.seq)
                            
                            import hashlib
                            checksum = hashlib.sha256(sequence.encode()).hexdigest()[:16]
                            
                            reference_data = {
                                'gene_name': gene_name,
                                'organism': "local_reference",
                                'accession': record.id,
                                'sequence': sequence,
                                'sequence_length': len(sequence),
                                'source': "LOCAL_FILE",
                                'fetch_timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
                                'sequence_checksum': checksum,
                                'sepi_metadata': {"local_file_path": str(ref_file)},
                                'reference_file_path': str(ref_file),
                                'quality_score': 0.8  # Lower quality, unknown provenance
                            }
                            
                            self.logger.info(f"Using local reference file for {gene_name}: {ref_file}")
                            return reference_data
                            
                    except Exception as e:
                        self.logger.warning(f"Failed to parse local reference {ref_file}: {e}")
                        continue
        
        return None
    
    def get_statistics(self) -> Dict[str, Any]:
        """Get reference manager statistics"""
        return {
            'sepi_availability': self.sepi_availability,
            'cached_references': len(self.reference_cache),
            'failed_fetches': len(self.failed_fetches),
            **self.stats
        }
    
    def clear_cache(self):
        """Clear reference cache"""
        self.reference_cache.clear()
        self.failed_fetches.clear()
        self.logger.info("Reference cache cleared")


class SEPIIntegrationTester:
    """
    Test suite for SEPI 2.0 integration
    """
    
    def __init__(self, workspace_root: str):
        """Initialize SEPI integration tester"""
        self.workspace_root = Path(workspace_root)
        self.logger = logging.getLogger('SEPIIntegrationTester')
        
        # Initialize configuration manager
        self.config_manager = SEPIConfigurationManager(workspace_root=str(self.workspace_root))
        
        # Initialize enhanced reference manager
        self.reference_manager = EnhancedSEPIReferenceManager(self.config_manager)
    
    async def test_gene_reference_fetching(self, gene_names: List[str]) -> Dict[str, Any]:
        """Test reference fetching for multiple genes"""
        
        results = {
            'successful_fetches': [],
            'failed_fetches': [],
            'total_tested': len(gene_names),
            'success_rate': 0.0
        }
        
        self.logger.info(f"Testing reference fetching for {len(gene_names)} genes")
        
        for gene_name in gene_names:
            try:
                self.logger.info(f"Testing fetch for gene: {gene_name}")
                
                reference = await self.reference_manager.get_reference_sequence(gene_name)
                
                if reference:
                    results['successful_fetches'].append({
                        'gene_name': gene_name,
                        'organism': reference['organism'],
                        'sequence_length': reference['sequence_length'],
                        'source': reference['source'],
                        'quality_score': reference['quality_score']
                    })
                    self.logger.info(f"✅ Successfully fetched {gene_name}")
                else:
                    results['failed_fetches'].append({
                        'gene_name': gene_name,
                        'error': 'No reference found'
                    })
                    self.logger.warning(f"❌ Failed to fetch {gene_name}")
                    
            except Exception as e:
                results['failed_fetches'].append({
                    'gene_name': gene_name,
                    'error': str(e)
                })
                self.logger.error(f"❌ Exception fetching {gene_name}: {e}")
        
        # Calculate success rate
        success_count = len(results['successful_fetches'])
        results['success_rate'] = (success_count / len(gene_names)) * 100 if gene_names else 0
        
        # Add statistics
        results['reference_manager_stats'] = self.reference_manager.get_statistics()
        
        self.logger.info(f"Testing complete: {success_count}/{len(gene_names)} successful ({results['success_rate']:.1f}%)")
        
        return results
    
    def save_test_results(self, results: Dict[str, Any], output_file: str):
        """Save test results to file"""
        
        output_path = Path(output_file)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_path, 'w') as f:
            json.dump(results, f, indent=2)
        
        self.logger.info(f"Test results saved to: {output_path}")


async def main():
    """Command line interface for SEPI integration testing"""
    import argparse
    
    parser = argparse.ArgumentParser(description="SEPI 2.0 Integration Tester")
    parser.add_argument("--workspace", required=True, help="Workspace root directory")
    parser.add_argument("--genes", nargs="+", default=['acrA', 'acrB', 'tolC', 'mexA', 'mexB'],
                       help="Gene names to test")
    parser.add_argument("--output", help="Output file for test results")
    parser.add_argument("--log-level", default="INFO", 
                       choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    
    args = parser.parse_args()
    
    # Setup logging
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    
    try:
        # Initialize tester
        tester = SEPIIntegrationTester(args.workspace)
        
        # Test gene reference fetching
        results = await tester.test_gene_reference_fetching(args.genes)
        
        # Save results
        if args.output:
            tester.save_test_results(results, args.output)
        
        # Print summary
        print("\n🧬 SEPI 2.0 Integration Test Results:")
        print(f"   Tested genes: {results['total_tested']}")
        print(f"   Successful: {len(results['successful_fetches'])}")
        print(f"   Failed: {len(results['failed_fetches'])}")
        print(f"   Success rate: {results['success_rate']:.1f}%")
        
        if results['successful_fetches']:
            print("\n✅ Successful fetches:")
            for fetch in results['successful_fetches']:
                print(f"   {fetch['gene_name']} from {fetch['organism']} ({fetch['source']})")
        
        if results['failed_fetches']:
            print("\n❌ Failed fetches:")
            for fetch in results['failed_fetches']:
                print(f"   {fetch['gene_name']}: {fetch['error']}")
        
        return 0 if results['success_rate'] > 50 else 1
        
    except Exception as e:
        print(f"❌ Integration test failed: {e}")
        return 1


if __name__ == "__main__":
    import asyncio
    sys.exit(asyncio.run(main()))