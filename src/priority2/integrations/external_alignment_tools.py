"""
External Alignment Tools Integration
-----------------------------------
Wrappers and integrations for external alignment tools (minimap2, mappy, parasail).
Provides unified interface and automatic tool selection based on dataset characteristics.
"""

import os
import subprocess
import logging
from pathlib import Path
from typing import List, Dict, Any, Optional, Union
from dataclasses import dataclass
from abc import ABC, abstractmethod

# Optional imports for different tools
try:
    import mappy as mp
except ImportError:
    mp = None

try:
    import parasail
except ImportError:
    parasail = None

@dataclass
class AlignmentResult:
    """Standardized alignment result format."""
    query_name: str
    target_name: str
    query_length: int
    target_length: int
    query_start: int
    query_end: int
    target_start: int
    target_end: int
    matches: int
    alignment_length: int
    mapq: int
    strand: str
    identity: float
    coverage: float

class BaseAligner(ABC):
    """Base class for alignment tool wrappers."""
    
    def __init__(self, reference_path: str, **kwargs):
        self.reference_path = reference_path
        self.logger = logging.getLogger(self.__class__.__name__)
        
    @abstractmethod
    def align(self, query_sequence: str, query_name: str) -> List[AlignmentResult]:
        """Align a single sequence and return results."""
        pass
    
    @abstractmethod
    def align_file(self, query_file: str, output_file: str) -> bool:
        """Align sequences from file and save results."""
        pass
    
    @abstractmethod
    def is_available(self) -> bool:
        """Check if the tool is available."""
        pass

class MappyAligner(BaseAligner):
    """Wrapper for mappy (minimap2 Python bindings)."""
    
    def __init__(self, reference_path: str, preset: str = "map-ont", threads: int = 4, **kwargs):
        super().__init__(reference_path, **kwargs)
        self.preset = preset
        self.threads = threads
        self._aligner = None
        
        if not self.is_available():
            raise RuntimeError("mappy is not available")
            
        self._initialize_aligner()
    
    def _initialize_aligner(self):
        """Initialize the mappy aligner."""
        try:
            self._aligner = mp.Aligner(
                self.reference_path, 
                preset=self.preset, 
                n_threads=self.threads
            )
            self.logger.info(f"Initialized mappy aligner with preset: {self.preset}")
        except Exception as e:
            raise RuntimeError(f"Failed to initialize mappy aligner: {e}")
    
    def align(self, query_sequence: str, query_name: str) -> List[AlignmentResult]:
        """Align a single sequence using mappy."""
        results = []
        
        try:
            for hit in self._aligner.map(query_sequence):
                result = AlignmentResult(
                    query_name=query_name,
                    target_name=hit.ctg,
                    query_length=len(query_sequence),
                    target_length=hit.ctg_len,
                    query_start=hit.q_st,
                    query_end=hit.q_en,
                    target_start=hit.r_st,
                    target_end=hit.r_en,
                    matches=hit.mlen,
                    alignment_length=hit.blen,
                    mapq=hit.mapq,
                    strand="+" if hit.strand == 1 else "-",
                    identity=(hit.mlen / hit.blen * 100) if hit.blen > 0 else 0,
                    coverage=((hit.q_en - hit.q_st) / len(query_sequence) * 100) if len(query_sequence) > 0 else 0
                )
                results.append(result)
        except Exception as e:
            self.logger.error(f"Alignment failed for {query_name}: {e}")
        
        return results
    
    def align_file(self, query_file: str, output_file: str) -> bool:
        """Align sequences from FASTA/FASTQ file."""
        try:
            with open(query_file, 'r') as qf, open(output_file, 'w') as of:
                for name, seq, qual in mp.fastx_read(qf):
                    alignments = self.align(seq, name)
                    for aln in alignments:
                        # Write in PAF format
                        of.write(f"{aln.query_name}\t{aln.query_length}\t{aln.query_start}\t{aln.query_end}\t"
                                f"{aln.strand}\t{aln.target_name}\t{aln.target_length}\t{aln.target_start}\t"
                                f"{aln.target_end}\t{aln.matches}\t{aln.alignment_length}\t{aln.mapq}\n")
            return True
        except Exception as e:
            self.logger.error(f"File alignment failed: {e}")
            return False
    
    def is_available(self) -> bool:
        """Check if mappy is available."""
        return mp is not None

class Minimap2Aligner(BaseAligner):
    """Wrapper for standalone minimap2 executable."""
    
    def __init__(self, reference_path: str, preset: str = "map-ont", threads: int = 4, 
                 minimap2_path: Optional[str] = None, **kwargs):
        super().__init__(reference_path, **kwargs)
        self.preset = preset
        self.threads = threads
        self.minimap2_path = minimap2_path or "minimap2"
        
        if not self.is_available():
            raise RuntimeError("minimap2 executable not found")
    
    def align(self, query_sequence: str, query_name: str) -> List[AlignmentResult]:
        """Align single sequence (not efficient for standalone minimap2)."""
        # For single sequences, write to temp file and use align_file
        import tempfile
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as temp_query:
            temp_query.write(f">{query_name}\n{query_sequence}\n")
            temp_query_path = temp_query.name
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.paf', delete=False) as temp_output:
            temp_output_path = temp_output.name
        
        try:
            success = self.align_file(temp_query_path, temp_output_path)
            if success:
                return self._parse_paf_file(temp_output_path)
            return []
        finally:
            os.unlink(temp_query_path)
            os.unlink(temp_output_path)
    
    def align_file(self, query_file: str, output_file: str) -> bool:
        """Align sequences using standalone minimap2."""
        try:
            cmd = [
                self.minimap2_path,
                "-x", self.preset,
                "-t", str(self.threads),
                "--paf-no-hit",
                self.reference_path,
                query_file
            ]
            
            with open(output_file, 'w') as of:
                result = subprocess.run(cmd, stdout=of, stderr=subprocess.PIPE, 
                                      text=True, check=True)
            
            self.logger.info(f"minimap2 alignment completed: {query_file} -> {output_file}")
            return True
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"minimap2 failed: {e.stderr}")
            return False
        except Exception as e:
            self.logger.error(f"minimap2 execution failed: {e}")
            return False
    
    def _parse_paf_file(self, paf_file: str) -> List[AlignmentResult]:
        """Parse PAF file and return alignment results."""
        results = []
        try:
            with open(paf_file, 'r') as f:
                for line in f:
                    if line.strip():
                        result = self._parse_paf_line(line.strip())
                        if result:
                            results.append(result)
        except Exception as e:
            self.logger.error(f"Failed to parse PAF file: {e}")
        return results
    
    def _parse_paf_line(self, line: str) -> Optional[AlignmentResult]:
        """Parse a single PAF line."""
        try:
            fields = line.split('\t')
            if len(fields) >= 12:
                return AlignmentResult(
                    query_name=fields[0],
                    target_name=fields[5],
                    query_length=int(fields[1]),
                    target_length=int(fields[6]),
                    query_start=int(fields[2]),
                    query_end=int(fields[3]),
                    target_start=int(fields[7]),
                    target_end=int(fields[8]),
                    matches=int(fields[9]),
                    alignment_length=int(fields[10]),
                    mapq=int(fields[11]),
                    strand=fields[4],
                    identity=(int(fields[9]) / int(fields[10]) * 100) if int(fields[10]) > 0 else 0,
                    coverage=((int(fields[3]) - int(fields[2])) / int(fields[1]) * 100) if int(fields[1]) > 0 else 0
                )
        except (ValueError, IndexError):
            pass
        return None
    
    def is_available(self) -> bool:
        """Check if minimap2 executable is available."""
        try:
            result = subprocess.run([self.minimap2_path, "--version"], 
                                  capture_output=True, text=True, timeout=5)
            return result.returncode == 0
        except Exception:
            return False

class ParasailAligner(BaseAligner):
    """Wrapper for parasail (SIMD-accelerated alignment)."""
    
    def __init__(self, reference_path: str, match_score: int = 2, mismatch_score: int = -1, 
                 gap_open: int = -1, gap_extend: int = -1, **kwargs):
        super().__init__(reference_path, **kwargs)
        self.match_score = match_score
        self.mismatch_score = mismatch_score
        self.gap_open = gap_open
        self.gap_extend = gap_extend
        
        if not self.is_available():
            raise RuntimeError("parasail is not available")
        
        self._load_reference()
    
    def _load_reference(self):
        """Load reference sequences."""
        self.reference_seqs = {}
        try:
            with open(self.reference_path, 'r') as f:
                current_name = None
                current_seq = []
                
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        if current_name:
                            self.reference_seqs[current_name] = ''.join(current_seq)
                        current_name = line[1:].split()[0]
                        current_seq = []
                    else:
                        current_seq.append(line)
                
                if current_name:
                    self.reference_seqs[current_name] = ''.join(current_seq)
                    
        except Exception as e:
            raise RuntimeError(f"Failed to load reference: {e}")
    
    def align(self, query_sequence: str, query_name: str) -> List[AlignmentResult]:
        """Align using parasail."""
        results = []
        
        for ref_name, ref_seq in self.reference_seqs.items():
            try:
                alignment = parasail.sw_trace_striped_16(
                    query_sequence, ref_seq,
                    self.gap_open, self.gap_extend,
                    parasail.blosum62
                )
                
                if alignment.score > 0:
                    # Calculate metrics
                    identity = (alignment.matches / alignment.length * 100) if alignment.length > 0 else 0
                    coverage = (alignment.end_query / len(query_sequence) * 100) if len(query_sequence) > 0 else 0
                    
                    result = AlignmentResult(
                        query_name=query_name,
                        target_name=ref_name,
                        query_length=len(query_sequence),
                        target_length=len(ref_seq),
                        query_start=alignment.beg_query,
                        query_end=alignment.end_query,
                        target_start=alignment.beg_ref,
                        target_end=alignment.end_ref,
                        matches=alignment.matches,
                        alignment_length=alignment.length,
                        mapq=min(60, int(alignment.score / 10)),  # Approximate MAPQ
                        strand="+",
                        identity=identity,
                        coverage=coverage
                    )
                    results.append(result)
                    
            except Exception as e:
                self.logger.error(f"Parasail alignment failed for {ref_name}: {e}")
        
        return results
    
    def align_file(self, query_file: str, output_file: str) -> bool:
        """Align sequences from file using parasail."""
        try:
            with open(output_file, 'w') as of:
                # Read FASTA file manually
                with open(query_file, 'r') as qf:
                    current_name = None
                    current_seq = []
                    
                    for line in qf:
                        line = line.strip()
                        if line.startswith('>'):
                            if current_name:
                                seq = ''.join(current_seq)
                                alignments = self.align(seq, current_name)
                                for aln in alignments:
                                    of.write(f"{aln.query_name}\t{aln.query_length}\t{aln.query_start}\t{aln.query_end}\t"
                                            f"{aln.strand}\t{aln.target_name}\t{aln.target_length}\t{aln.target_start}\t"
                                            f"{aln.target_end}\t{aln.matches}\t{aln.alignment_length}\t{aln.mapq}\n")
                            current_name = line[1:].split()[0]
                            current_seq = []
                        else:
                            current_seq.append(line)
                    
                    if current_name:
                        seq = ''.join(current_seq)
                        alignments = self.align(seq, current_name)
                        for aln in alignments:
                            of.write(f"{aln.query_name}\t{aln.query_length}\t{aln.query_start}\t{aln.query_end}\t"
                                    f"{aln.strand}\t{aln.target_name}\t{aln.target_length}\t{aln.target_start}\t"
                                    f"{aln.target_end}\t{aln.matches}\t{aln.alignment_length}\t{aln.mapq}\n")
            return True
        except Exception as e:
            self.logger.error(f"Parasail file alignment failed: {e}")
            return False
    
    def is_available(self) -> bool:
        """Check if parasail is available."""
        return parasail is not None

class AlignmentToolManager:
    """Manager for selecting and using alignment tools."""
    
    def __init__(self, reference_path: str):
        self.reference_path = reference_path
        self.logger = logging.getLogger("AlignmentToolManager")
        self.available_tools = self._detect_available_tools()
    
    def _detect_available_tools(self) -> Dict[str, bool]:
        """Detect which alignment tools are available."""
        tools = {}
        
        # Test mappy
        try:
            aligner = MappyAligner(self.reference_path)
            tools['mappy'] = True
            self.logger.info("mappy is available")
        except Exception:
            tools['mappy'] = False
        
        # Test minimap2
        try:
            aligner = Minimap2Aligner(self.reference_path)
            tools['minimap2'] = True
            self.logger.info("minimap2 is available")
        except Exception:
            tools['minimap2'] = False
        
        # Test parasail
        try:
            aligner = ParasailAligner(self.reference_path)
            tools['parasail'] = True
            self.logger.info("parasail is available")
        except Exception:
            tools['parasail'] = False
        
        return tools
    
    def get_best_aligner(self, dataset_size: int, sequence_type: str = "long") -> BaseAligner:
        """Select the best aligner based on dataset characteristics."""
        
        # Decision logic based on dataset size and type
        if dataset_size > 10000:
            # Large datasets - prefer minimap2 for speed
            if self.available_tools.get('minimap2', False):
                return Minimap2Aligner(self.reference_path)
            elif self.available_tools.get('mappy', False):
                return MappyAligner(self.reference_path)
        
        elif dataset_size < 100:
            # Small datasets - parasail for accuracy
            if self.available_tools.get('parasail', False):
                return ParasailAligner(self.reference_path)
        
        # Default to mappy if available
        if self.available_tools.get('mappy', False):
            return MappyAligner(self.reference_path)
        elif self.available_tools.get('minimap2', False):
            return Minimap2Aligner(self.reference_path)
        elif self.available_tools.get('parasail', False):
            return ParasailAligner(self.reference_path)
        
        raise RuntimeError("No alignment tools are available")
    
    def benchmark_tools(self, test_query: str, test_name: str = "test") -> Dict[str, Dict[str, Any]]:
        """Benchmark available alignment tools."""
        import time
        
        results = {}
        
        for tool_name, available in self.available_tools.items():
            if not available:
                continue
                
            try:
                # Create aligner
                if tool_name == 'mappy':
                    aligner = MappyAligner(self.reference_path)
                elif tool_name == 'minimap2':
                    aligner = Minimap2Aligner(self.reference_path)
                elif tool_name == 'parasail':
                    aligner = ParasailAligner(self.reference_path)
                else:
                    continue
                
                # Benchmark
                start_time = time.time()
                alignments = aligner.align(test_query, test_name)
                end_time = time.time()
                
                results[tool_name] = {
                    'execution_time': end_time - start_time,
                    'num_alignments': len(alignments),
                    'average_identity': sum(aln.identity for aln in alignments) / len(alignments) if alignments else 0,
                    'success': True
                }
                
            except Exception as e:
                results[tool_name] = {
                    'execution_time': None,
                    'num_alignments': 0,
                    'error': str(e),
                    'success': False
                }
        
        return results
    
    def get_tool_recommendations(self, dataset_size: int) -> List[str]:
        """Get tool recommendations based on dataset characteristics."""
        recommendations = []
        
        if dataset_size > 50000:
            recommendations.extend(['minimap2', 'mappy'])
        elif dataset_size > 1000:
            recommendations.extend(['mappy', 'minimap2'])
        else:
            recommendations.extend(['parasail', 'mappy'])
        
        # Filter by availability
        return [tool for tool in recommendations if self.available_tools.get(tool, False)]