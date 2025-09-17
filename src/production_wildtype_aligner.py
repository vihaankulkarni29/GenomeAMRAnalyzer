#!/usr/bin/env python3
"""
Production WildType Aligner - Senior Bioinformatician Grade
Advanced protein alignment system with intelligent reference management

This module provides a production-ready alignment system that:
1. Integrates seamlessly with ProductionFastaExtractor outputs
2. Uses SEPI 2.0 for dynamic wild-type reference fetching
3. Provides multi-algorithm alignment with quality assessment
4. Maintains complete provenance tracking from extraction to alignment
5. Generates comprehensive manifests for SubScan integration

Author: GenomeAMRAnalyzer Pipeline
Version: 2.0 - Senior Bioinformatician Grade
"""

import os
import sys
import time
import hashlib
import logging
import asyncio
import subprocess
import tempfile
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Set, Any
from dataclasses import dataclass, asdict
import json
import csv
from collections import defaultdict
import re

# BioPython imports for sequence handling
try:
    from Bio import SeqIO, Align
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Align import PairwiseAligner
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False
    print("Warning: BioPython not available. Some features may be limited.")


@dataclass
class ReferenceSequence:
    """Reference sequence with complete provenance"""
    gene_name: str
    organism: str
    accession: str
    sequence: str
    sequence_length: int
    source: str  # SEPI_2.0, MANUAL, CACHED
    fetch_timestamp: str
    sequence_checksum: str
    sepi_metadata: Dict[str, Any]
    reference_file_path: str = ""
    quality_score: float = 0.0


@dataclass
class AlignmentQualityMetrics:
    """Comprehensive alignment quality assessment"""
    identity_percentage: float
    similarity_percentage: float
    coverage_query: float
    coverage_reference: float
    alignment_length: int
    query_length: int
    reference_length: int
    gap_count: int
    gap_percentage: float
    score: float
    normalized_score: float
    confidence_level: str  # HIGH, MEDIUM, LOW, FAILED
    algorithm_used: str


@dataclass
class AlignmentResult:
    """Complete alignment result with provenance"""
    accession: str
    gene_name: str
    query_sequence: str
    reference_sequence: str
    query_checksum: str
    reference_checksum: str
    reference_organism: str
    alignment_file_path: str
    quality_metrics: AlignmentQualityMetrics
    alignment_success: bool = True
    error_message: str = ""
    alignment_timestamp: str = ""
    alignment_algorithm: str = ""
    reference_source: str = ""
    # Provenance chain
    extraction_source_file: str = ""
    coordinate_source_file: str = ""


class SEPIReferenceManager:
    """
    Senior bioinformatician-grade reference management with SEPI 2.0 integration
    """
    
    def __init__(self, sepi_path: str, cache_dir: str, email: str):
        """Initialize SEPI reference manager"""
        self.sepi_path = Path(sepi_path)
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.email = email
        
        # Reference cache and metadata
        self.reference_cache: Dict[str, ReferenceSequence] = {}
        self.organism_cache: Dict[str, str] = {}
        
        # Setup logging
        self.logger = logging.getLogger('SEPIReferenceManager')
        
        # SEPI availability check
        self.sepi_available = self._check_sepi_availability()
        
    def _check_sepi_availability(self) -> bool:
        """Check if SEPI 2.0 is available and functional"""
        try:
            if not self.sepi_path.exists():
                self.logger.warning(f"SEPI path not found: {self.sepi_path}")
                return False
            
            # Test SEPI with help command
            result = subprocess.run(
                [sys.executable, str(self.sepi_path), "--help"],
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
        key_string = f"{gene_name}_{organism}".lower().replace(" ", "_")
        return hashlib.sha256(key_string.encode()).hexdigest()[:16]
    
    def _extract_organism_from_protein_header(self, protein_header: str) -> Optional[str]:
        """Extract organism information from protein FASTA header"""
        # Parse header format: >accession|gene_name|coordinates|strand|contig|checksum
        parts = protein_header.split('|')
        if len(parts) >= 1:
            accession = parts[0].lstrip('>')
            
            # Use organism mapping cache or try to infer from accession
            if accession in self.organism_cache:
                return self.organism_cache[accession]
            
            # Try to infer organism from accession patterns
            organism_patterns = {
                r'^GCF_.*': 'inferred_from_ncbi',
                r'^GCA_.*': 'inferred_from_ncbi'
            }
            
            for pattern, organism_type in organism_patterns.items():
                if re.match(pattern, accession):
                    return organism_type
        
        return None
    
    def _determine_reference_organisms(self, gene_name: str, query_organism: Optional[str] = None) -> List[str]:
        """Determine optimal reference organisms for gene fetching"""
        
        # Gene-specific organism preferences (senior bioinformatician knowledge)
        gene_organism_preferences = {
            'mexA': ['Pseudomonas aeruginosa PAO1', 'Pseudomonas aeruginosa', 'Pseudomonas'],
            'mexB': ['Pseudomonas aeruginosa PAO1', 'Pseudomonas aeruginosa', 'Pseudomonas'],
            'mexC': ['Pseudomonas aeruginosa', 'Pseudomonas'],
            'mexD': ['Pseudomonas aeruginosa', 'Pseudomonas'],
            'mexE': ['Pseudomonas aeruginosa', 'Pseudomonas'],
            'mexF': ['Pseudomonas aeruginosa', 'Pseudomonas'],
            'oprM': ['Pseudomonas aeruginosa PAO1', 'Pseudomonas aeruginosa', 'Pseudomonas'],
            'oprN': ['Pseudomonas aeruginosa', 'Pseudomonas'],
            'acrA': ['Escherichia coli K-12 MG1655', 'Escherichia coli', 'Escherichia'],
            'acrB': ['Escherichia coli K-12 MG1655', 'Escherichia coli', 'Escherichia'],
            'tolC': ['Escherichia coli K-12 MG1655', 'Escherichia coli', 'Escherichia']
        }
        
        gene_lower = gene_name.lower()
        
        # Priority 1: Gene-specific preferences
        if gene_lower in gene_organism_preferences:
            organisms = gene_organism_preferences[gene_lower].copy()
        else:
            # Priority 2: Query organism if provided
            organisms = []
            if query_organism:
                organisms.append(query_organism)
        
        # Priority 3: Add query organism to preferences if not already there
        if query_organism and query_organism not in organisms:
            organisms.insert(0, query_organism)
        
        # Priority 4: Generic fallbacks
        fallback_organisms = [
            'Escherichia coli',
            'Pseudomonas aeruginosa', 
            'Klebsiella pneumoniae',
            'Acinetobacter baumannii'
        ]
        
        for fallback in fallback_organisms:
            if fallback not in organisms:
                organisms.append(fallback)
        
        return organisms
    
    async def fetch_reference_with_sepi(self, gene_name: str, organisms: List[str]) -> Optional[ReferenceSequence]:
        """Fetch reference sequence using SEPI 2.0 with fallback strategy"""
        
        if not self.sepi_available:
            self.logger.warning("SEPI not available, cannot fetch references")
            return None
        
        current_timestamp = time.strftime('%Y-%m-%d %H:%M:%S')
        
        for organism in organisms:
            try:
                self.logger.info(f"Attempting SEPI fetch: {gene_name} from {organism}")
                
                # Create temporary output directory
                temp_output = self.cache_dir / f"sepi_temp_{gene_name}_{int(time.time())}"
                temp_output.mkdir(exist_ok=True)
                
                # Construct SEPI command
                sepi_cmd = [
                    sys.executable, str(self.sepi_path),
                    "--organism", organism,
                    "--proteins", gene_name,
                    "--assembly_level", "complete_genome",
                    "--output", str(temp_output / f"{gene_name}_{organism.replace(' ', '_')}"),
                    "--email", self.email
                ]
                
                # Execute SEPI command
                result = await asyncio.create_subprocess_exec(
                    *sepi_cmd,
                    stdout=asyncio.subprocess.PIPE,
                    stderr=asyncio.subprocess.PIPE,
                    cwd=temp_output
                )
                
                stdout, stderr = await result.communicate()
                
                if result.returncode == 0:
                    # Parse SEPI output
                    fasta_files = list(temp_output.glob("**/*.fasta"))
                    
                    if fasta_files:
                        # Read the first FASTA file found
                        fasta_file = fasta_files[0]
                        
                        if BIOPYTHON_AVAILABLE:
                            records = list(SeqIO.parse(fasta_file, "fasta"))
                            if records:
                                record = records[0]
                                sequence = str(record.seq)
                                accession = record.id
                                
                                # Calculate checksum
                                checksum = hashlib.sha256(sequence.encode()).hexdigest()[:16]
                                
                                # Move reference to permanent cache
                                cache_filename = f"{gene_name}_{organism.replace(' ', '_')}_reference.fasta"
                                cache_path = self.cache_dir / cache_filename
                                shutil.copy2(fasta_file, cache_path)
                                
                                reference = ReferenceSequence(
                                    gene_name=gene_name,
                                    organism=organism,
                                    accession=accession,
                                    sequence=sequence,
                                    sequence_length=len(sequence),
                                    source="SEPI_2.0",
                                    fetch_timestamp=current_timestamp,
                                    sequence_checksum=checksum,
                                    sepi_metadata={
                                        "sepi_output_dir": str(temp_output),
                                        "sepi_command": " ".join(sepi_cmd),
                                        "sepi_stdout": stdout.decode() if stdout else "",
                                        "organism_attempted": organism
                                    },
                                    reference_file_path=str(cache_path),
                                    quality_score=1.0  # High quality from SEPI
                                )
                                
                                # Cache the reference
                                cache_key = self._generate_cache_key(gene_name, organism)
                                self.reference_cache[cache_key] = reference
                                
                                # Cleanup temporary directory
                                shutil.rmtree(temp_output, ignore_errors=True)
                                
                                self.logger.info(f"Successfully fetched {gene_name} from {organism} via SEPI")
                                return reference
                
                # Cleanup on failure
                shutil.rmtree(temp_output, ignore_errors=True)
                
            except Exception as e:
                self.logger.warning(f"SEPI fetch failed for {gene_name} from {organism}: {e}")
                continue
        
        self.logger.warning(f"All SEPI attempts failed for {gene_name}")
        return None
    
    async def get_reference_sequence(self, gene_name: str, query_organism: Optional[str] = None) -> Optional[ReferenceSequence]:
        """Get reference sequence with intelligent fallback strategy"""
        
        # Check cache first
        organisms = self._determine_reference_organisms(gene_name, query_organism)
        
        for organism in organisms:
            cache_key = self._generate_cache_key(gene_name, organism)
            if cache_key in self.reference_cache:
                self.logger.info(f"Using cached reference for {gene_name} from {organism}")
                return self.reference_cache[cache_key]
        
        # Try SEPI fetch
        reference = await self.fetch_reference_with_sepi(gene_name, organisms)
        if reference:
            return reference
        
        # Fallback: Check local reference files (if any)
        local_refs = self._check_local_references(gene_name, organisms)
        if local_refs:
            return local_refs
        
        self.logger.error(f"Failed to obtain reference for {gene_name}")
        return None
    
    def _check_local_references(self, gene_name: str, organisms: List[str]) -> Optional[ReferenceSequence]:
        """Check for local reference files as fallback"""
        
        # Check common reference directories
        reference_dirs = [
            self.cache_dir,
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
                        if BIOPYTHON_AVAILABLE:
                            records = list(SeqIO.parse(ref_file, "fasta"))
                            if records:
                                record = records[0]
                                sequence = str(record.seq)
                                
                                reference = ReferenceSequence(
                                    gene_name=gene_name,
                                    organism="local_reference",
                                    accession=record.id,
                                    sequence=sequence,
                                    sequence_length=len(sequence),
                                    source="LOCAL_FILE",
                                    fetch_timestamp=time.strftime('%Y-%m-%d %H:%M:%S'),
                                    sequence_checksum=hashlib.sha256(sequence.encode()).hexdigest()[:16],
                                    sepi_metadata={"local_file_path": str(ref_file)},
                                    reference_file_path=str(ref_file),
                                    quality_score=0.8  # Lower quality, unknown provenance
                                )
                                
                                self.logger.info(f"Using local reference file for {gene_name}: {ref_file}")
                                return reference
                                
                    except Exception as e:
                        self.logger.warning(f"Failed to parse local reference {ref_file}: {e}")
                        continue
        
        return None


class AlignmentEngine:
    """
    Production-grade alignment engine with multiple algorithms
    """
    
    def __init__(self, output_dir: str, emboss_water_path: str = "water"):
        """Initialize alignment engine"""
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        self.emboss_water_path = emboss_water_path
        self.logger = logging.getLogger('AlignmentEngine')
        
        # Check algorithm availability
        self.emboss_available = self._check_emboss_availability()
        self.biopython_available = BIOPYTHON_AVAILABLE
        
        # Algorithm priorities
        self.algorithm_priority = ['EMBOSS_WATER', 'BIOPYTHON_PAIRWISE']
        
    def _check_emboss_availability(self) -> bool:
        """Check if EMBOSS WATER is available"""
        try:
            result = subprocess.run(
                [self.emboss_water_path, "-help"],
                capture_output=True, text=True, timeout=10
            )
            
            if result.returncode == 0 or "water" in result.stderr.lower():
                self.logger.info("EMBOSS WATER availability confirmed")
                return True
            else:
                self.logger.warning("EMBOSS WATER not available")
                return False
                
        except Exception as e:
            self.logger.warning(f"EMBOSS WATER check failed: {e}")
            return False
    
    def _select_algorithm(self, query_length: int, reference_length: int) -> str:
        """Select optimal alignment algorithm based on sequence properties"""
        
        # Algorithm selection logic (senior bioinformatician approach)
        max_length = max(query_length, reference_length)
        length_ratio = max(query_length, reference_length) / min(query_length, reference_length)
        
        # For very long sequences or highly divergent lengths, prefer EMBOSS
        if max_length > 2000 or length_ratio > 3.0:
            if self.emboss_available:
                return 'EMBOSS_WATER'
        
        # For standard protein alignments, use availability priority
        for algorithm in self.algorithm_priority:
            if algorithm == 'EMBOSS_WATER' and self.emboss_available:
                return 'EMBOSS_WATER'
            elif algorithm == 'BIOPYTHON_PAIRWISE' and self.biopython_available:
                return 'BIOPYTHON_PAIRWISE'
        
        raise RuntimeError("No alignment algorithms available")
    
    async def align_with_emboss_water(self, query_seq: str, reference_seq: str, 
                                    output_file: Path) -> Tuple[bool, Dict[str, Any]]:
        """Perform alignment using EMBOSS WATER"""
        
        try:
            # Create temporary files for sequences
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as query_tmp:
                query_tmp.write(f">query\n{query_seq}\n")
                query_tmp_path = query_tmp.name
            
            with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as ref_tmp:
                ref_tmp.write(f">reference\n{reference_seq}\n")
                ref_tmp_path = ref_tmp.name
            
            # EMBOSS WATER command
            water_cmd = [
                self.emboss_water_path,
                "-asequence", query_tmp_path,
                "-bsequence", ref_tmp_path,
                "-outfile", str(output_file),
                "-gapopen", "10.0",
                "-gapextend", "0.5",
                "-auto"
            ]
            
            # Execute WATER
            result = await asyncio.create_subprocess_exec(
                *water_cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE
            )
            
            stdout, stderr = await result.communicate()
            
            # Cleanup temporary files
            os.unlink(query_tmp_path)
            os.unlink(ref_tmp_path)
            
            if result.returncode == 0 and output_file.exists():
                # Parse WATER output for statistics
                alignment_stats = self._parse_water_output(output_file)
                return True, alignment_stats
            else:
                error_msg = stderr.decode() if stderr else "Unknown WATER error"
                self.logger.error(f"WATER alignment failed: {error_msg}")
                return False, {"error": error_msg}
                
        except Exception as e:
            self.logger.error(f"EMBOSS WATER execution failed: {e}")
            return False, {"error": str(e)}
    
    def _parse_water_output(self, water_file: Path) -> Dict[str, Any]:
        """Parse EMBOSS WATER output file for alignment statistics"""
        
        stats = {
            "identity_percentage": 0.0,
            "similarity_percentage": 0.0,
            "gaps": 0,
            "score": 0.0,
            "alignment_length": 0
        }
        
        try:
            with open(water_file, 'r') as f:
                content = f.read()
            
            # Parse identity
            identity_match = re.search(r'# Identity:\s+(\d+)/(\d+)\s+\(([0-9.]+)%\)', content)
            if identity_match:
                stats["identity_percentage"] = float(identity_match.group(3))
                stats["alignment_length"] = int(identity_match.group(2))
            
            # Parse similarity
            similarity_match = re.search(r'# Similarity:\s+(\d+)/(\d+)\s+\(([0-9.]+)%\)', content)
            if similarity_match:
                stats["similarity_percentage"] = float(similarity_match.group(3))
            
            # Parse gaps
            gaps_match = re.search(r'# Gaps:\s+(\d+)/(\d+)\s+\(([0-9.]+)%\)', content)
            if gaps_match:
                stats["gaps"] = int(gaps_match.group(1))
            
            # Parse score
            score_match = re.search(r'# Score:\s+([0-9.]+)', content)
            if score_match:
                stats["score"] = float(score_match.group(1))
                
        except Exception as e:
            self.logger.warning(f"Failed to parse WATER output: {e}")
        
        return stats
    
    def align_with_biopython(self, query_seq: str, reference_seq: str) -> Tuple[bool, Dict[str, Any]]:
        """Perform alignment using BioPython as fallback"""
        
        if not BIOPYTHON_AVAILABLE:
            return False, {"error": "BioPython not available"}
        
        try:
            # Configure pairwise aligner
            aligner = PairwiseAligner()
            aligner.match_score = 2
            aligner.mismatch_score = -1
            aligner.open_gap_score = -2
            aligner.extend_gap_score = -0.5
            
            # Perform alignment
            alignments = aligner.align(query_seq, reference_seq)
            
            if alignments:
                best_alignment = alignments[0]
                
                # Calculate statistics
                alignment_length = len(best_alignment.query)
                matches = sum(1 for i, j in zip(best_alignment.query, best_alignment.target) 
                            if i == j and i != '-' and j != '-')
                gaps = alignment_length - len(query_seq.replace('-', ''))
                
                identity_percentage = (matches / alignment_length) * 100 if alignment_length > 0 else 0
                
                stats = {
                    "identity_percentage": identity_percentage,
                    "similarity_percentage": identity_percentage,  # Simplified
                    "gaps": gaps,
                    "score": best_alignment.score,
                    "alignment_length": alignment_length
                }
                
                return True, stats
            else:
                return False, {"error": "No alignments generated"}
                
        except Exception as e:
            self.logger.error(f"BioPython alignment failed: {e}")
            return False, {"error": str(e)}
    
    async def perform_alignment(self, query_seq: str, reference_seq: str, 
                              output_file: Path, algorithm: Optional[str] = None) -> Tuple[bool, Dict[str, Any]]:
        """Perform alignment with specified or auto-selected algorithm"""
        
        if not algorithm:
            algorithm = self._select_algorithm(len(query_seq), len(reference_seq))
        
        self.logger.info(f"Using alignment algorithm: {algorithm}")
        
        if algorithm == 'EMBOSS_WATER':
            return await self.align_with_emboss_water(query_seq, reference_seq, output_file)
        elif algorithm == 'BIOPYTHON_PAIRWISE':
            # For BioPython, create a simple alignment file
            success, stats = self.align_with_biopython(query_seq, reference_seq)
            if success:
                with open(output_file, 'w') as f:
                    f.write(f"# BioPython Pairwise Alignment\n")
                    f.write(f"# Identity: {stats['identity_percentage']:.1f}%\n")
                    f.write(f"# Score: {stats['score']}\n")
                    f.write(f"Query: {query_seq}\n")
                    f.write(f"Reference: {reference_seq}\n")
            return success, stats
        else:
            return False, {"error": f"Unsupported algorithm: {algorithm}"}


class QualityAssessment:
    """
    Comprehensive alignment quality assessment
    """
    
    def __init__(self):
        """Initialize quality assessment"""
        self.logger = logging.getLogger('QualityAssessment')
        
        # Quality thresholds (senior bioinformatician standards)
        self.quality_thresholds = {
            'HIGH': {'identity': 95.0, 'coverage': 90.0, 'score_ratio': 0.8},
            'MEDIUM': {'identity': 80.0, 'coverage': 70.0, 'score_ratio': 0.5},
            'LOW': {'identity': 60.0, 'coverage': 50.0, 'score_ratio': 0.3}
        }
    
    def calculate_coverage(self, alignment_length: int, query_length: int, reference_length: int) -> Tuple[float, float]:
        """Calculate query and reference coverage"""
        query_coverage = (alignment_length / query_length) * 100 if query_length > 0 else 0
        reference_coverage = (alignment_length / reference_length) * 100 if reference_length > 0 else 0
        return query_coverage, reference_coverage
    
    def normalize_score(self, score: float, query_length: int, reference_length: int) -> float:
        """Normalize alignment score by sequence length"""
        max_length = max(query_length, reference_length)
        return score / max_length if max_length > 0 else 0
    
    def assess_alignment_quality(self, alignment_stats: Dict[str, Any], 
                               query_length: int, reference_length: int,
                               algorithm: str) -> AlignmentQualityMetrics:
        """Perform comprehensive quality assessment"""
        
        identity = alignment_stats.get('identity_percentage', 0.0)
        similarity = alignment_stats.get('similarity_percentage', identity)
        alignment_length = alignment_stats.get('alignment_length', 0)
        gap_count = alignment_stats.get('gaps', 0)
        score = alignment_stats.get('score', 0.0)
        
        # Calculate coverage
        query_coverage, reference_coverage = self.calculate_coverage(
            alignment_length, query_length, reference_length
        )
        
        # Calculate gap percentage
        gap_percentage = (gap_count / alignment_length) * 100 if alignment_length > 0 else 0
        
        # Normalize score
        normalized_score = self.normalize_score(score, query_length, reference_length)
        
        # Determine confidence level
        confidence_level = self._determine_confidence_level(
            identity, query_coverage, reference_coverage, normalized_score
        )
        
        return AlignmentQualityMetrics(
            identity_percentage=identity,
            similarity_percentage=similarity,
            coverage_query=query_coverage,
            coverage_reference=reference_coverage,
            alignment_length=alignment_length,
            query_length=query_length,
            reference_length=reference_length,
            gap_count=gap_count,
            gap_percentage=gap_percentage,
            score=score,
            normalized_score=normalized_score,
            confidence_level=confidence_level,
            algorithm_used=algorithm
        )
    
    def _determine_confidence_level(self, identity: float, query_coverage: float, 
                                  reference_coverage: float, normalized_score: float) -> str:
        """Determine overall confidence level based on multiple metrics"""
        
        min_coverage = min(query_coverage, reference_coverage)
        
        # Check against thresholds
        for level in ['HIGH', 'MEDIUM', 'LOW']:
            thresholds = self.quality_thresholds[level]
            
            if (identity >= thresholds['identity'] and 
                min_coverage >= thresholds['coverage'] and 
                normalized_score >= thresholds['score_ratio']):
                return level
        
        return 'FAILED'


class ProductionWildTypeAligner:
    """
    Senior bioinformatician-grade wild-type alignment system
    Complete integration with ProductionFastaExtractor and SEPI 2.0
    """
    
    def __init__(self, output_dir: str, sepi_path: str, email: str, 
                 emboss_water_path: str = "water", max_concurrent: int = 10):
        """Initialize production wild-type aligner"""
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Create output subdirectories
        self.directories = {
            'alignments': self.output_dir / "alignments",
            'references': self.output_dir / "references", 
            'quality_reports': self.output_dir / "quality_reports",
            'logs': self.output_dir / "logs",
            'manifests': self.output_dir / "manifests"
        }
        
        for dir_path in self.directories.values():
            dir_path.mkdir(parents=True, exist_ok=True)
        
        self.email = email
        self.max_concurrent = max_concurrent
        
        # Initialize components
        self.reference_manager = SEPIReferenceManager(
            sepi_path, str(self.directories['references']), email
        )
        self.alignment_engine = AlignmentEngine(
            str(self.directories['alignments']), emboss_water_path
        )
        self.quality_assessment = QualityAssessment()
        
        # Setup logging
        self._setup_logging()
        
        # Data storage
        self.alignment_results: Dict[str, List[AlignmentResult]] = {}
        self.failed_alignments: Dict[str, str] = {}
        self.processed_accessions: Set[str] = set()
        
        # Statistics
        self.stats = {
            'total_accessions': 0,
            'successful_accessions': 0,
            'failed_accessions': 0,
            'total_alignments': 0,
            'successful_alignments': 0,
            'failed_alignments': 0,
            'high_quality_alignments': 0,
            'medium_quality_alignments': 0,
            'low_quality_alignments': 0
        }
        
        self.logger.info("Production WildType Aligner initialized")
    
    def _setup_logging(self):
        """Setup comprehensive logging"""
        log_file = self.directories['logs'] / "wildtype_alignment.log"
        
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - [%(filename)s:%(lineno)d] - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        
        self.logger = logging.getLogger('ProductionWildTypeAligner')
        self.logger.setLevel(logging.INFO)
        
        # File handler
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(formatter)
        self.logger.addHandler(file_handler)
    
    def _parse_protein_header(self, header: str) -> Dict[str, str]:
        """Parse protein FASTA header for metadata"""
        # Format: >accession|gene_name|start-end|strand|contig|checksum
        parts = header.lstrip('>').split('|')
        
        if len(parts) >= 6:
            return {
                'accession': parts[0],
                'gene_name': parts[1],
                'coordinates': parts[2],
                'strand': parts[3],
                'contig': parts[4],
                'checksum': parts[5]
            }
        else:
            return {'accession': header.lstrip('>'), 'gene_name': 'unknown'}
    
    async def process_accession_proteins(self, accession: str, protein_file: Path, 
                                       target_genes: List[str]) -> List[AlignmentResult]:
        """Process all proteins for a single accession"""
        
        alignment_results = []
        current_timestamp = time.strftime('%Y-%m-%d %H:%M:%S')
        
        try:
            if not BIOPYTHON_AVAILABLE:
                raise RuntimeError("BioPython required for protein processing")
            
            # Parse protein file
            protein_records = list(SeqIO.parse(protein_file, "fasta"))
            
            self.logger.info(f"Processing {len(protein_records)} proteins for {accession}")
            
            for record in protein_records:
                header_info = self._parse_protein_header(record.description)
                gene_name = header_info.get('gene_name', 'unknown')
                
                # Filter by target genes
                if gene_name.lower() not in [g.lower() for g in target_genes]:
                    continue
                
                try:
                    # Get reference sequence
                    reference = await self.reference_manager.get_reference_sequence(gene_name)
                    
                    if not reference:
                        result = AlignmentResult(
                            accession=accession,
                            gene_name=gene_name,
                            query_sequence=str(record.seq),
                            reference_sequence="",
                            query_checksum=header_info.get('checksum', ''),
                            reference_checksum="",
                            reference_organism="",
                            alignment_file_path="",
                            quality_metrics=AlignmentQualityMetrics(
                                identity_percentage=0, similarity_percentage=0,
                                coverage_query=0, coverage_reference=0,
                                alignment_length=0, query_length=len(record.seq),
                                reference_length=0, gap_count=0, gap_percentage=0,
                                score=0, normalized_score=0, confidence_level="FAILED",
                                algorithm_used="NONE"
                            ),
                            alignment_success=False,
                            error_message="Reference sequence not found",
                            alignment_timestamp=current_timestamp
                        )
                        alignment_results.append(result)
                        continue
                    
                    # Perform alignment
                    alignment_file = self.directories['alignments'] / f"{accession}_{gene_name}_alignment.water"
                    
                    success, alignment_stats = await self.alignment_engine.perform_alignment(
                        str(record.seq), reference.sequence, alignment_file
                    )
                    
                    if success:
                        # Quality assessment
                        quality_metrics = self.quality_assessment.assess_alignment_quality(
                            alignment_stats, len(record.seq), reference.sequence_length,
                            alignment_stats.get('algorithm', 'UNKNOWN')
                        )
                        
                        result = AlignmentResult(
                            accession=accession,
                            gene_name=gene_name,
                            query_sequence=str(record.seq),
                            reference_sequence=reference.sequence,
                            query_checksum=header_info.get('checksum', ''),
                            reference_checksum=reference.sequence_checksum,
                            reference_organism=reference.organism,
                            alignment_file_path=str(alignment_file),
                            quality_metrics=quality_metrics,
                            alignment_success=True,
                            alignment_timestamp=current_timestamp,
                            alignment_algorithm=quality_metrics.algorithm_used,
                            reference_source=reference.source
                        )
                        
                        # Update statistics
                        self.stats['successful_alignments'] += 1
                        if quality_metrics.confidence_level == 'HIGH':
                            self.stats['high_quality_alignments'] += 1
                        elif quality_metrics.confidence_level == 'MEDIUM':
                            self.stats['medium_quality_alignments'] += 1
                        elif quality_metrics.confidence_level == 'LOW':
                            self.stats['low_quality_alignments'] += 1
                        
                    else:
                        error_msg = alignment_stats.get('error', 'Alignment failed')
                        result = AlignmentResult(
                            accession=accession,
                            gene_name=gene_name,
                            query_sequence=str(record.seq),
                            reference_sequence=reference.sequence,
                            query_checksum=header_info.get('checksum', ''),
                            reference_checksum=reference.sequence_checksum,
                            reference_organism=reference.organism,
                            alignment_file_path="",
                            quality_metrics=AlignmentQualityMetrics(
                                identity_percentage=0, similarity_percentage=0,
                                coverage_query=0, coverage_reference=0,
                                alignment_length=0, query_length=len(record.seq),
                                reference_length=reference.sequence_length,
                                gap_count=0, gap_percentage=0,
                                score=0, normalized_score=0, confidence_level="FAILED",
                                algorithm_used="NONE"
                            ),
                            alignment_success=False,
                            error_message=error_msg,
                            alignment_timestamp=current_timestamp,
                            reference_source=reference.source
                        )
                        
                        self.stats['failed_alignments'] += 1
                    
                    self.stats['total_alignments'] += 1
                    alignment_results.append(result)
                    
                except Exception as e:
                    self.logger.error(f"Error processing {gene_name} for {accession}: {e}")
                    self.stats['failed_alignments'] += 1
                    self.stats['total_alignments'] += 1
            
            return alignment_results
            
        except Exception as e:
            self.logger.error(f"Failed to process proteins for {accession}: {e}")
            return []
    
    async def process_batch_alignments(self, protein_dir: str, 
                                     target_genes: List[str]) -> bool:
        """Process batch of protein files for alignment"""
        
        protein_path = Path(protein_dir)
        
        if not protein_path.exists():
            self.logger.error(f"Protein directory not found: {protein_path}")
            return False
        
        # Find protein files
        protein_files = list(protein_path.glob("*_proteins.fasta"))
        
        if not protein_files:
            self.logger.error(f"No protein files found in {protein_path}")
            return False
        
        self.logger.info(f"Processing {len(protein_files)} protein files")
        self.stats['total_accessions'] = len(protein_files)
        
        # Process files concurrently
        semaphore = asyncio.Semaphore(self.max_concurrent)
        
        async def process_single_file(protein_file: Path):
            async with semaphore:
                accession = protein_file.stem.replace('_proteins', '')
                
                if accession in self.processed_accessions:
                    self.logger.warning(f"Accession {accession} already processed")
                    return
                
                try:
                    results = await self.process_accession_proteins(
                        accession, protein_file, target_genes
                    )
                    
                    if results:
                        self.alignment_results[accession] = results
                        successful_alignments = len([r for r in results if r.alignment_success])
                        
                        if successful_alignments > 0:
                            self.stats['successful_accessions'] += 1
                            self.logger.info(f"Processed {accession}: {successful_alignments} alignments")
                        else:
                            self.stats['failed_accessions'] += 1
                            self.failed_alignments[accession] = "All alignments failed"
                    else:
                        self.stats['failed_accessions'] += 1
                        self.failed_alignments[accession] = "Protein processing failed"
                    
                    self.processed_accessions.add(accession)
                    
                except Exception as e:
                    self.logger.error(f"Error processing {accession}: {e}")
                    self.failed_alignments[accession] = str(e)
                    self.stats['failed_accessions'] += 1
        
        # Execute all processing tasks
        tasks = [process_single_file(pf) for pf in protein_files]
        await asyncio.gather(*tasks, return_exceptions=True)
        
        # Generate final manifest
        self.generate_alignment_manifest()
        
        # Log final statistics
        self.logger.info(f"Alignment processing complete:")
        self.logger.info(f"  Successful accessions: {self.stats['successful_accessions']}")
        self.logger.info(f"  Failed accessions: {self.stats['failed_accessions']}")
        self.logger.info(f"  Total alignments: {self.stats['total_alignments']}")
        self.logger.info(f"  High quality: {self.stats['high_quality_alignments']}")
        self.logger.info(f"  Medium quality: {self.stats['medium_quality_alignments']}")
        self.logger.info(f"  Low quality: {self.stats['low_quality_alignments']}")
        
        return self.stats['successful_accessions'] > 0
    
    def generate_alignment_manifest(self, output_file: Optional[str] = None) -> str:
        """Generate comprehensive alignment manifest"""
        
        if output_file is None:
            output_file = str(self.directories['manifests'] / "alignment_manifest.json")
        
        manifest = {
            "manifest_version": "1.0",
            "generated_timestamp": time.strftime('%Y-%m-%d %H:%M:%S'),
            "aligner_version": "2.0",
            "total_accessions_processed": len(self.alignment_results),
            "successful_alignments": self.stats['successful_alignments'],
            "failed_alignments": self.stats['failed_alignments'],
            "quality_distribution": {
                "high_quality": self.stats['high_quality_alignments'],
                "medium_quality": self.stats['medium_quality_alignments'],
                "low_quality": self.stats['low_quality_alignments']
            },
            "accessions": {}
        }
        
        # Add detailed data for each accession
        for accession, results in self.alignment_results.items():
            successful_results = [r for r in results if r.alignment_success]
            
            manifest["accessions"][accession] = {
                "alignment_files": [r.alignment_file_path for r in successful_results],
                "proteins_aligned": len(successful_results),
                "alignment_failures": len(results) - len(successful_results),
                "alignments": [
                    {
                        "gene_name": r.gene_name,
                        "query_sequence_checksum": r.query_checksum,
                        "reference_sequence_checksum": r.reference_checksum,
                        "reference_organism": r.reference_organism,
                        "alignment_file": r.alignment_file_path,
                        "quality_metrics": asdict(r.quality_metrics),
                        "alignment_timestamp": r.alignment_timestamp,
                        "alignment_algorithm": r.alignment_algorithm,
                        "reference_source": r.reference_source
                    }
                    for r in successful_results
                ]
            }
        
        # Add failed accessions
        if self.failed_alignments:
            manifest["failed_accessions"] = self.failed_alignments
        
        # Write JSON manifest
        Path(output_file).parent.mkdir(parents=True, exist_ok=True)
        with open(output_file, 'w') as f:
            json.dump(manifest, f, indent=2)
        
        self.logger.info(f"Alignment manifest saved to: {output_file}")
        return output_file


async def main():
    """Command line interface for production wild-type aligner"""
    import argparse
    
    parser = argparse.ArgumentParser(description="Production WildType Aligner")
    parser.add_argument("--protein-dir", required=True, 
                       help="Directory containing protein FASTA files from ProductionFastaExtractor")
    parser.add_argument("--target-genes", nargs="+", required=True,
                       help="Target gene names for alignment")
    parser.add_argument("--output-dir", required=True, 
                       help="Output directory for alignments")
    parser.add_argument("--sepi-path", required=True,
                       help="Path to SEPI 2.0 script")
    parser.add_argument("--email", required=True,
                       help="Email for NCBI API (required by SEPI)")
    parser.add_argument("--emboss-water-path", default="water",
                       help="Path to EMBOSS WATER executable")
    parser.add_argument("--max-concurrent", type=int, default=10,
                       help="Maximum concurrent alignments")
    parser.add_argument("--log-level", default="INFO", 
                       choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    
    args = parser.parse_args()
    
    # Setup logging
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    
    # Initialize aligner
    aligner = ProductionWildTypeAligner(
        args.output_dir, args.sepi_path, args.email,
        args.emboss_water_path, args.max_concurrent
    )
    
    try:
        # Process batch alignments
        success = await aligner.process_batch_alignments(
            args.protein_dir, args.target_genes
        )
        
        if success:
            print("Alignment processing completed successfully")
            print(f"Results saved to: {args.output_dir}")
            print(f"Manifest: {aligner.directories['manifests'] / 'alignment_manifest.json'}")
            return 0
        else:
            print("Alignment processing failed")
            return 1
            
    except Exception as e:
        print(f"Alignment processing failed: {e}")
        return 1


if __name__ == "__main__":
    import asyncio
    sys.exit(asyncio.run(main()))