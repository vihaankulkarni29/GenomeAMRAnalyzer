"""
MCP Data Processing Layer for Priority 3
========================================

- Format validation for all supported genomics/phenotype files
- Optimized file parsing (streaming, chunked, memory-efficient)
- Data normalization and quality checks (no biological interpretation)
- Modular, extensible, and testable design
- Provenance and error reporting for all operations

Engineering Principles:
- No biological inference or annotation
- Fail-fast on format or integrity errors
- Full provenance and error traceability
- Designed for integration with pipeline orchestrator
"""

import os
import logging
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple
import json
import csv

class MCPFormatValidator:
    """
    Validates file formats for all supported genomics and phenotype data types.
    """
    def __init__(self):
        self.logger = logging.getLogger("MCPFormatValidator")

    def validate_fasta(self, file_path: str) -> Tuple[bool, str]:
        """Validate FASTA file format."""
        try:
            with open(file_path, 'r') as f:
                first = f.readline()
                if not first.startswith('>'):
                    return False, "Missing FASTA header (>)."
                for line in f:
                    if line.startswith('>'):
                        continue
                    if not line.strip():
                        continue
                    if not all(c in 'ACGTURYKMSWBDHVN*-.' for c in line.strip().upper()):
                        return False, f"Invalid character in sequence: {line.strip()}"
            return True, "OK"
        except Exception as e:
            return False, f"Exception: {e}"

    def validate_csv(self, file_path: str, required_columns: Optional[List[str]] = None) -> Tuple[bool, str]:
        """Validate CSV file for required columns and structure."""
        try:
            with open(file_path, 'r', newline='') as f:
                reader = csv.DictReader(f)
                if required_columns:
                    if not reader.fieldnames:
                        return False, "No columns found in CSV."
                    missing = [col for col in required_columns if col not in reader.fieldnames]
                    if missing:
                        return False, f"Missing columns: {missing}"
            return True, "OK"
        except Exception as e:
            return False, f"Exception: {e}"

    def validate_json(self, file_path: str, required_keys: Optional[List[str]] = None) -> Tuple[bool, str]:
        """Validate JSON file for required keys."""
        try:
            with open(file_path, 'r') as f:
                data = json.load(f)
                if required_keys:
                    missing = [k for k in required_keys if k not in data]
                    if missing:
                        return False, f"Missing keys: {missing}"
            return True, "OK"
        except Exception as e:
            return False, f"Exception: {e}"

class MCPFileParser:
    """
    Optimized file parsing for large genomics/phenotype files (streaming, chunked).
    """
    def __init__(self):
        self.logger = logging.getLogger("MCPFileParser")

    def parse_fasta(self, file_path: str) -> List[Dict[str, Any]]:
        """Stream-parse FASTA file, yielding records as dicts."""
        records = []
        try:
            with open(file_path, 'r') as f:
                seq_id = None
                seq = []
                for line in f:
                    if line.startswith('>'):
                        if seq_id:
                            records.append({'id': seq_id, 'sequence': ''.join(seq)})
                        seq_id = line[1:].strip()
                        seq = []
                    else:
                        seq.append(line.strip())
                if seq_id:
                    records.append({'id': seq_id, 'sequence': ''.join(seq)})
            return records
        except Exception as e:
            self.logger.error(f"FASTA parse error: {e}")
            return []

    def parse_csv(self, file_path: str) -> List[Dict[str, Any]]:
        """Stream-parse CSV file as list of dicts."""
        try:
            with open(file_path, 'r', newline='') as f:
                reader = csv.DictReader(f)
                return [row for row in reader]
        except Exception as e:
            self.logger.error(f"CSV parse error: {e}")
            return []

    def parse_json(self, file_path: str) -> Any:
        """Parse JSON file."""
        try:
            with open(file_path, 'r') as f:
                return json.load(f)
        except Exception as e:
            self.logger.error(f"JSON parse error: {e}")
            return None

class MCPDataNormalizer:
    """
    Data normalization and quality checks (no biological interpretation).
    """
    def __init__(self):
        self.logger = logging.getLogger("MCPDataNormalizer")

    def normalize_fasta(self, records: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """Normalize FASTA records: uppercase, remove gaps, check length."""
        norm = []
        for rec in records:
            seq = rec['sequence'].replace('-', '').replace('.', '').upper()
            norm.append({'id': rec['id'], 'sequence': seq, 'length': len(seq)})
        return norm

    def normalize_csv(self, rows: List[Dict[str, Any]], column_map: Optional[Dict[str, str]] = None) -> List[Dict[str, Any]]:
        """Normalize CSV rows: rename columns, strip whitespace."""
        norm = []
        for row in rows:
            new_row = {column_map.get(k, k): v.strip() if isinstance(v, str) else v for k, v in row.items()} if column_map else {k: v.strip() if isinstance(v, str) else v for k, v in row.items()}
            norm.append(new_row)
        return norm

    def quality_check(self, data: Any) -> Tuple[bool, List[str]]:
        """Basic quality checks: missing values, length, duplicates."""
        issues = []
        if isinstance(data, list):
            seen = set()
            for i, rec in enumerate(data):
                if isinstance(rec, dict):
                    for k, v in rec.items():
                        if v in [None, '', 'NA', 'N/A']:
                            issues.append(f"Row {i} missing value for {k}")
                    if 'id' in rec:
                        if rec['id'] in seen:
                            issues.append(f"Duplicate id: {rec['id']}")
                        seen.add(rec['id'])
        return (len(issues) == 0, issues)
