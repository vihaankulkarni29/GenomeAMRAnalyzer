"""
Internal EMBOSS-style Alignment Generator
Replaces EMBOSS water dependency with internal implementation that produces
EMBOSS-compatible output format for SubScan compatibility.
"""

import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from Bio import SeqIO, Align
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np

logger = logging.getLogger(__name__)

class EmbossStyleAligner:
    """
    Internal aligner that produces EMBOSS water-compatible output.
    Uses BioPython's pairwise alignment with EMBOSS-style formatting.
    """
    
    def __init__(self, gap_open: float = 10.0, gap_extend: float = 0.5, matrix: str = "BLOSUM62"):
        self.gap_open = gap_open
        self.gap_extend = gap_extend
        self.matrix = matrix
        
        # Setup BioPython aligner
        self.aligner = Align.PairwiseAligner()
        self.aligner.substitution_matrix = Align.substitution_matrices.load(matrix)
        self.aligner.open_gap_score = -gap_open
        self.aligner.extend_gap_score = -gap_extend
        self.aligner.mode = 'local'  # Local alignment like EMBOSS water
    
    def align_sequences(self, seq1: str, seq2: str, seq1_id: str = "seq1", seq2_id: str = "seq2") -> str:
        """
        Align two sequences and return EMBOSS-style formatted output.
        
        Args:
            seq1: First sequence (usually reference)
            seq2: Second sequence (usually query)
            seq1_id: Identifier for first sequence
            seq2_id: Identifier for second sequence
        
        Returns:
            EMBOSS-style alignment text
        """
        try:
            # Perform alignment
            alignments = self.aligner.align(seq1, seq2)
            if not alignments:
                return self._format_no_alignment(seq1_id, seq2_id)
            
            # Get best alignment
            best_alignment = alignments[0]
            
            # Format as EMBOSS-style output
            return self._format_emboss_style(best_alignment, seq1_id, seq2_id, seq1, seq2)
            
        except Exception as e:
            logger.error(f"Alignment failed for {seq1_id} vs {seq2_id}: {e}")
            return self._format_error_alignment(seq1_id, seq2_id, str(e))
    
    def _format_emboss_style(self, alignment, seq1_id: str, seq2_id: str, seq1: str, seq2: str) -> str:
        """Format alignment in EMBOSS water style."""
        
        # Get alignment details
        score = alignment.score
        aligned_seq1, aligned_seq2 = alignment.aligned
        
        # Calculate statistics
        length = len(str(alignment).split('\n')[0])
        identity = self._calculate_identity(str(alignment))
        similarity = self._calculate_similarity(str(alignment))
        gaps = self._calculate_gaps(str(alignment))
        
        # Create EMBOSS-style header
        header = f"""########################################
# Program: water (GenomeAMRAnalyzer internal)
# Rundate: {self._get_timestamp()}
# Commandline: water -asequence {seq1_id} -bsequence {seq2_id}
# Align_format: pair
# Report_file: stdout
########################################

#=======================================
#
# Aligned_sequences: 2
# 1: {seq1_id}
# 2: {seq2_id}
# Matrix: {self.matrix}
# Gap_penalty: {self.gap_open}
# Extend_penalty: {self.gap_extend}
#
# Length: {length}
# Identity: {identity['matches']}/{identity['total']} ({identity['percent']:.1f}%)
# Similarity: {similarity['matches']}/{similarity['total']} ({similarity['percent']:.1f}%)
# Gaps: {gaps['gaps']}/{gaps['total']} ({gaps['percent']:.1f}%)
# Score: {score:.1f}
#
#
#=======================================
"""
        
        # Format alignment body
        alignment_body = self._format_alignment_body(alignment, seq1_id, seq2_id)
        
        # EMBOSS-style footer
        footer = f"""
#---------------------------------------
#---------------------------------------
"""
        
        return header + alignment_body + footer
    
    def _format_alignment_body(self, alignment, seq1_id: str, seq2_id: str) -> str:
        """Format the alignment sequences in EMBOSS style."""
        
        # Get aligned sequences as strings
        alignment_str = str(alignment)
        lines = alignment_str.split('\n')
        
        if len(lines) < 3:
            return "# No alignment found\n"
        
        aligned1 = lines[0]
        match_line = lines[1] 
        aligned2 = lines[2]
        
        # Break into chunks of 50 characters (EMBOSS standard)
        chunk_size = 50
        body = ""
        
        start1 = 1
        start2 = 1
        
        for i in range(0, len(aligned1), chunk_size):
            chunk1 = aligned1[i:i+chunk_size]
            chunk_match = match_line[i:i+chunk_size]
            chunk2 = aligned2[i:i+chunk_size]
            
            # Calculate end positions (accounting for gaps)
            end1 = start1 + len(chunk1.replace('-', '')) - 1
            end2 = start2 + len(chunk2.replace('-', '')) - 1
            
            # Format chunk
            body += f"\n{seq1_id:<15} {start1:>6} {chunk1} {end1:<6}"
            body += f"\n{'':<15} {'':<6} {chunk_match} {'':<6}"
            body += f"\n{seq2_id:<15} {start2:>6} {chunk2} {end2:<6}\n"
            
            # Update start positions for next chunk
            start1 = end1 + 1 if end1 >= start1 else start1
            start2 = end2 + 1 if end2 >= start2 else start2
        
        return body
    
    def _calculate_identity(self, alignment_str: str) -> Dict[str, float]:
        """Calculate identity statistics from alignment string."""
        lines = alignment_str.split('\n')
        if len(lines) < 3:
            return {"matches": 0, "total": 0, "percent": 0.0}
        
        match_line = lines[1]
        matches = match_line.count('|') + match_line.count(':')
        total = len(match_line.replace(' ', ''))
        percent = (matches / total * 100) if total > 0 else 0.0
        
        return {"matches": matches, "total": total, "percent": percent}
    
    def _calculate_similarity(self, alignment_str: str) -> Dict[str, float]:
        """Calculate similarity statistics (identity + similar residues)."""
        lines = alignment_str.split('\n')
        if len(lines) < 3:
            return {"matches": 0, "total": 0, "percent": 0.0}
        
        match_line = lines[1]
        # Count exact matches (|) and similarities (:, +)
        matches = match_line.count('|') + match_line.count(':') + match_line.count('+')
        total = len(match_line.replace(' ', ''))
        percent = (matches / total * 100) if total > 0 else 0.0
        
        return {"matches": matches, "total": total, "percent": percent}
    
    def _calculate_gaps(self, alignment_str: str) -> Dict[str, float]:
        """Calculate gap statistics."""
        lines = alignment_str.split('\n')
        if len(lines) < 3:
            return {"gaps": 0, "total": 0, "percent": 0.0}
        
        seq1 = lines[0]
        seq2 = lines[2]
        gaps = seq1.count('-') + seq2.count('-')
        total = len(seq1) + len(seq2)
        percent = (gaps / total * 100) if total > 0 else 0.0
        
        return {"gaps": gaps, "total": total, "percent": percent}
    
    def _get_timestamp(self) -> str:
        """Get current timestamp in EMBOSS format."""
        from datetime import datetime
        return datetime.now().strftime("%a %b %d %H:%M:%S %Y")
    
    def _format_no_alignment(self, seq1_id: str, seq2_id: str) -> str:
        """Format output when no alignment is found."""
        return f"""########################################
# Program: water (GenomeAMRAnalyzer internal)
# No significant alignment found between {seq1_id} and {seq2_id}
########################################
"""
    
    def _format_error_alignment(self, seq1_id: str, seq2_id: str, error: str) -> str:
        """Format output when alignment fails."""
        return f"""########################################
# Program: water (GenomeAMRAnalyzer internal)
# Error aligning {seq1_id} and {seq2_id}: {error}
########################################
"""
    
    def align_fasta_files(self, query_fasta: Path, reference_fasta: Path, output_file: Path) -> bool:
        """
        Align sequences from FASTA files and write EMBOSS-style output.
        
        Args:
            query_fasta: Path to query FASTA file
            reference_fasta: Path to reference FASTA file  
            output_file: Path to write alignment output
        
        Returns:
            True if successful, False otherwise
        """
        try:
            # Read sequences
            query_records = list(SeqIO.parse(query_fasta, "fasta"))
            ref_records = list(SeqIO.parse(reference_fasta, "fasta"))
            
            if not query_records or not ref_records:
                logger.error(f"No sequences found in input files")
                return False
            
            # Use first sequence from each file
            query_seq = str(query_records[0].seq)
            ref_seq = str(ref_records[0].seq)
            query_id = query_records[0].id
            ref_id = ref_records[0].id
            
            # Perform alignment
            alignment_output = self.align_sequences(ref_seq, query_seq, ref_id, query_id)
            
            # Write output
            output_file.parent.mkdir(parents=True, exist_ok=True)
            with output_file.open("w") as f:
                f.write(alignment_output)
            
            logger.info(f"Alignment written to {output_file}")
            return True
            
        except Exception as e:
            logger.error(f"FASTA alignment failed: {e}")
            return False


def create_mock_emboss_alignment(query_id: str, ref_id: str, output_file: Path) -> bool:
    """Create a mock EMBOSS-style alignment for testing."""
    
    mock_output = f"""########################################
# Program: water (GenomeAMRAnalyzer internal - MOCK)
# Rundate: {EmbossStyleAligner()._get_timestamp()}
# Commandline: water -asequence {ref_id} -bsequence {query_id}
# Align_format: pair
# Report_file: {output_file}
########################################

#=======================================
#
# Aligned_sequences: 2
# 1: {ref_id}
# 2: {query_id}
# Matrix: BLOSUM62
# Gap_penalty: 10.0
# Extend_penalty: 0.5
#
# Length: 300
# Identity: 285/300 (95.0%)
# Similarity: 295/300 (98.3%)
# Gaps: 5/300 (1.7%)
# Score: 1420.5
#
#
#=======================================

{ref_id:<15}      1 MVKIEFFKRTGQPANFDPVITGTQSIFGNHSPLTAMGQSDFKQKFRVLFI     50
                     ||||||||||||||||||||||||||||||||||||||||||||||||||
{query_id:<15}      1 MVKIEFFKRTGQPANFDPVITGTQSIFGNHSPLTAMGQSDFKQKFRVLFI     50

{ref_id:<15}     51 TSRRSVRFNRQHLIRFVFNHDENGTTFTPRLEQAGKAIGVAYEFSLLFAD    100
                     ||||||||||||||||||||||||||||||||||||||||||||||||||
{query_id:<15}     51 TSRRSVRFNRQHLIRFVFNHDENGTTFTPRLEQAGKAIGVAYEFSLLFAD    100

#---------------------------------------
#---------------------------------------
"""
    
    try:
        output_file.parent.mkdir(parents=True, exist_ok=True)
        with output_file.open("w") as f:
            f.write(mock_output)
        return True
    except Exception as e:
        logger.error(f"Failed to create mock alignment: {e}")
        return False


if __name__ == "__main__":
    # Test the aligner
    aligner = EmbossStyleAligner()
    
    # Test sequences
    seq1 = "MVKIEFFKRTGQPANFDPVITGTQSIFGNHSPLTAMGQSDFKQKFRVLFI"
    seq2 = "MVKIEFFKRTGQPANFDPVITGTQSIFGNHSPLTAMGQSDFKQKFRVLFI"
    
    result = aligner.align_sequences(seq1, seq2, "reference", "query")
    print(result)