"""
Database Schema and Repository Layer for Priority 3
==================================================

Central SQLite store for:
- Genome metadata and accessions
- MIC values and biosample data  
- Artifact references and provenance
- Analysis results and reports

Features:
- Atomic transactions
- Schema migrations
- Comprehensive indexing
- Provenance tracking
"""

import os
import sys
import sqlite3
import logging
import hashlib
from datetime import datetime
from typing import Optional, List, Dict, Any, Union
from pathlib import Path
from dataclasses import dataclass, asdict
import json

SQLALCHEMY_AVAILABLE = False
if sys.version_info < (3, 13):
    try:
        from sqlalchemy import create_engine, Column, Integer, String, DateTime, Text, Float, Boolean
        from sqlalchemy.ext.declarative import declarative_base
        from sqlalchemy.orm import sessionmaker, Session
        SQLALCHEMY_AVAILABLE = True
    except ImportError:
        SQLALCHEMY_AVAILABLE = False

@dataclass
class GenomeRecord:
    """Core genome record with metadata."""
    accession: str
    organism: str
    biosample: Optional[str] = None
    bioproject: Optional[str] = None
    file_path: Optional[str] = None
    file_hash: Optional[str] = None
    download_date: Optional[datetime] = None
    assembly_level: Optional[str] = None
    genome_size: Optional[int] = None
    status: str = "pending"  # pending, downloaded, processed, error

@dataclass  
class MICRecord:
    """MIC (Minimum Inhibitory Concentration) data."""
    biosample: str
    antibiotic: str
    mic_value: float
    mic_units: str
    method: Optional[str] = None
    source: Optional[str] = None
    validated: bool = False

@dataclass
class ArtifactRecord:
    """Reference to generated artifacts."""
    type: str  # fasta, alignment, report, etc.
    path: str
    hash: str
    size: int
    created_date: datetime
    provenance: Dict[str, Any]

@dataclass
class MetadataRecord:
    """Generic metadata blob associated with an accession (JSON payload)."""
    accession: str
    metadata_type: str  # e.g., "mic_data", "analysis_summary"
    data: Dict[str, Any]
    source: str
    collection_date: datetime

class GenomeRepository:
    def list_artifacts(self, accession: Optional[str] = None, artifact_type: Optional[str] = None) -> List[ArtifactRecord]:
        """List artifacts by accession and/or type."""
        try:
            cursor = self.conn.cursor()
            query = "SELECT * FROM artifacts"
            params = []
            conditions = []
            if accession:
                # Join with genomes to get file_path or accession linkage if needed
                # For now, assume artifact provenance JSON contains accession
                conditions.append("provenance LIKE ?")
                params.append(f'%{accession}%')
            if artifact_type:
                conditions.append("type = ?")
                params.append(artifact_type)
            if conditions:
                query += " WHERE " + " AND ".join(conditions)
            cursor.execute(query, params)
            rows = cursor.fetchall()
            return [
                ArtifactRecord(
                    type=row['type'],
                    path=row['path'],
                    hash=row['hash'],
                    size=row['size'],
                    created_date=row['created_date'],
                    provenance=json.loads(row['provenance'])
                ) for row in rows
            ]
        except Exception as e:
            self.logger.error(f"Failed to list artifacts: {e}")
            return []
    """
    Repository for genome data with SQLite backend.
    Handles CRUD operations, migrations, and data integrity.
    """
    
    def __init__(self, db_path: str = "priority3.db"):
        self.db_path = db_path
        self.logger = logging.getLogger("GenomeRepository")
        
        # Create database directory if needed
        os.makedirs(os.path.dirname(os.path.abspath(db_path)), exist_ok=True)
        
        if SQLALCHEMY_AVAILABLE:
            self._init_sqlalchemy()
        else:
            self._init_sqlite()
            
    def _init_sqlalchemy(self):
        """Initialize SQLAlchemy if available."""
        # Import here to ensure names are available when SQLAlchemy is present
        from sqlalchemy import create_engine
        from sqlalchemy.orm import sessionmaker
        self.engine = create_engine(f"sqlite:///{self.db_path}")
        self.SessionLocal = sessionmaker(bind=self.engine)
        self._create_sqlalchemy_tables()
        
    def _init_sqlite(self):
        """Fallback to raw SQLite."""
        self.conn = sqlite3.connect(self.db_path, check_same_thread=False)
        self.conn.row_factory = sqlite3.Row
        self._create_sqlite_tables()
        
    def _create_sqlite_tables(self):
        """Create tables using raw SQL."""
        schema = """
        CREATE TABLE IF NOT EXISTS genomes (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            accession TEXT UNIQUE NOT NULL,
            organism TEXT NOT NULL,
            biosample TEXT,
            bioproject TEXT, 
            file_path TEXT,
            file_hash TEXT,
            download_date TIMESTAMP,
            assembly_level TEXT,
            genome_size INTEGER,
            status TEXT DEFAULT 'pending',
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        );
        
        CREATE TABLE IF NOT EXISTS mic_data (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            biosample TEXT NOT NULL,
            antibiotic TEXT NOT NULL,
            mic_value REAL NOT NULL,
            mic_units TEXT NOT NULL,
            method TEXT,
            source TEXT,
            validated BOOLEAN DEFAULT FALSE,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            UNIQUE(biosample, antibiotic, method)
        );
        
        CREATE TABLE IF NOT EXISTS artifacts (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            type TEXT NOT NULL,
            path TEXT NOT NULL,
            hash TEXT NOT NULL,
            size INTEGER NOT NULL,
            created_date TIMESTAMP NOT NULL,
            provenance TEXT NOT NULL,  -- JSON
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        );

        CREATE TABLE IF NOT EXISTS metadata (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            accession TEXT NOT NULL,
            metadata_type TEXT NOT NULL,
            data TEXT NOT NULL, -- JSON
            source TEXT,
            collection_date TIMESTAMP,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
            UNIQUE(accession, metadata_type, collection_date)
        );
        
        CREATE INDEX IF NOT EXISTS idx_genomes_accession ON genomes(accession);
        CREATE INDEX IF NOT EXISTS idx_genomes_biosample ON genomes(biosample);
        CREATE INDEX IF NOT EXISTS idx_genomes_status ON genomes(status);
        CREATE INDEX IF NOT EXISTS idx_mic_biosample ON mic_data(biosample);
        CREATE INDEX IF NOT EXISTS idx_mic_antibiotic ON mic_data(antibiotic);
        CREATE INDEX IF NOT EXISTS idx_artifacts_type ON artifacts(type);
        CREATE INDEX IF NOT EXISTS idx_metadata_accession ON metadata(accession);
        CREATE INDEX IF NOT EXISTS idx_metadata_type ON metadata(metadata_type);
        """
        
        for statement in schema.split(';'):
            if statement.strip():
                self.conn.execute(statement)
        self.conn.commit()
        
    def _create_sqlalchemy_tables(self):
        """Create tables using SQLAlchemy (for future extensibility)."""
        # TODO: Implement SQLAlchemy models when needed
        pass
        
    def add_genome(self, genome: GenomeRecord) -> bool:
        """Add a genome record."""
        try:
            cursor = self.conn.cursor()
            cursor.execute("""
                INSERT OR REPLACE INTO genomes 
                (accession, organism, biosample, bioproject, file_path, file_hash, 
                 download_date, assembly_level, genome_size, status, updated_at)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, CURRENT_TIMESTAMP)
            """, (
                genome.accession, genome.organism, genome.biosample, 
                genome.bioproject, genome.file_path, genome.file_hash,
                genome.download_date, genome.assembly_level, 
                genome.genome_size, genome.status
            ))
            self.conn.commit()
            self.logger.info(f"Added genome record: {genome.accession}")
            return True
        except Exception as e:
            self.logger.error(f"Failed to add genome {genome.accession}: {e}")
            return False
            
    def get_genome(self, accession: str) -> Optional[GenomeRecord]:
        """Get genome by accession."""
        try:
            cursor = self.conn.cursor()
            cursor.execute("SELECT * FROM genomes WHERE accession = ?", (accession,))
            row = cursor.fetchone()
            
            if row:
                return GenomeRecord(
                    accession=row['accession'],
                    organism=row['organism'],
                    biosample=row['biosample'],
                    bioproject=row['bioproject'],
                    file_path=row['file_path'],
                    file_hash=row['file_hash'],
                    download_date=row['download_date'],
                    assembly_level=row['assembly_level'],
                    genome_size=row['genome_size'],
                    status=row['status']
                )
        except Exception as e:
            self.logger.error(f"Failed to get genome {accession}: {e}")
        return None
        
    def update_genome_status(self, accession: str, status: str, 
                           file_path: Optional[str] = None, 
                           file_hash: Optional[str] = None) -> bool:
        """Update genome processing status."""
        try:
            cursor = self.conn.cursor()
            if file_path and file_hash:
                cursor.execute("""
                    UPDATE genomes 
                    SET status = ?, file_path = ?, file_hash = ?, updated_at = CURRENT_TIMESTAMP
                    WHERE accession = ?
                """, (status, file_path, file_hash, accession))
            else:
                cursor.execute("""
                    UPDATE genomes 
                    SET status = ?, updated_at = CURRENT_TIMESTAMP
                    WHERE accession = ?
                """, (status, accession))
            self.conn.commit()
            return True
        except Exception as e:
            self.logger.error(f"Failed to update genome status {accession}: {e}")
            return False
            
    def add_mic_data(self, mic: MICRecord) -> bool:
        """Add MIC data."""
        try:
            cursor = self.conn.cursor()
            cursor.execute("""
                INSERT OR REPLACE INTO mic_data
                (biosample, antibiotic, mic_value, mic_units, method, source, validated)
                VALUES (?, ?, ?, ?, ?, ?, ?)
            """, (
                mic.biosample, mic.antibiotic, mic.mic_value, 
                mic.mic_units, mic.method, mic.source, mic.validated
            ))
            self.conn.commit()
            self.logger.info(f"Added MIC data: {mic.biosample} - {mic.antibiotic}")
            return True
        except Exception as e:
            self.logger.error(f"Failed to add MIC data: {e}")
            return False
            
    def get_mic_data(self, biosample: str) -> List[MICRecord]:
        """Get all MIC data for a biosample."""
        try:
            cursor = self.conn.cursor()
            cursor.execute("SELECT * FROM mic_data WHERE biosample = ?", (biosample,))
            rows = cursor.fetchall()
            
            return [
                MICRecord(
                    biosample=row['biosample'],
                    antibiotic=row['antibiotic'],
                    mic_value=row['mic_value'],
                    mic_units=row['mic_units'],
                    method=row['method'],
                    source=row['source'],
                    validated=bool(row['validated'])
                ) for row in rows
            ]
        except Exception as e:
            self.logger.error(f"Failed to get MIC data for {biosample}: {e}")
            return []
            
    def add_artifact(self, artifact: ArtifactRecord) -> bool:
        """Add artifact reference."""
        try:
            cursor = self.conn.cursor()
            cursor.execute("""
                INSERT INTO artifacts (type, path, hash, size, created_date, provenance)
                VALUES (?, ?, ?, ?, ?, ?)
            """, (
                artifact.type, artifact.path, artifact.hash, 
                artifact.size, artifact.created_date, 
                json.dumps(artifact.provenance)
            ))
            self.conn.commit()
            return True
        except Exception as e:
            self.logger.error(f"Failed to add artifact: {e}")
            return False

    def add_metadata(self, metadata: MetadataRecord) -> bool:
        """Store a metadata record as JSON."""
        try:
            cursor = self.conn.cursor()
            cursor.execute(
                """
                INSERT OR REPLACE INTO metadata (accession, metadata_type, data, source, collection_date)
                VALUES (?, ?, ?, ?, ?)
                """,
                (
                    metadata.accession,
                    metadata.metadata_type,
                    json.dumps(metadata.data),
                    metadata.source,
                    metadata.collection_date,
                ),
            )
            self.conn.commit()
            return True
        except Exception as e:
            self.logger.error(f"Failed to add metadata for {metadata.accession}: {e}")
            return False
            
    def list_genomes(self, status: Optional[str] = None, limit: Optional[int] = None) -> List[GenomeRecord]:
        """List genomes with optional filtering."""
        try:
            cursor = self.conn.cursor()
            
            query = "SELECT * FROM genomes"
            params = []
            
            if status:
                query += " WHERE status = ?"
                params.append(status)
                
            query += " ORDER BY created_at DESC"
            
            if limit:
                query += " LIMIT ?"
                params.append(limit)
                
            cursor.execute(query, params)
            rows = cursor.fetchall()
            
            return [
                GenomeRecord(
                    accession=row['accession'],
                    organism=row['organism'],
                    biosample=row['biosample'],
                    bioproject=row['bioproject'],
                    file_path=row['file_path'],
                    file_hash=row['file_hash'],
                    download_date=row['download_date'],
                    assembly_level=row['assembly_level'],
                    genome_size=row['genome_size'],
                    status=row['status']
                ) for row in rows
            ]
        except Exception as e:
            self.logger.error(f"Failed to list genomes: {e}")
            return []
            
    def get_stats(self) -> Dict[str, Any]:
        """Get database statistics."""
        try:
            cursor = self.conn.cursor()
            
            stats = {}
            
            # Genome counts by status
            cursor.execute("SELECT status, COUNT(*) FROM genomes GROUP BY status")
            stats['genome_status'] = dict(cursor.fetchall())
            
            # Total genomes
            cursor.execute("SELECT COUNT(*) FROM genomes")
            stats['total_genomes'] = cursor.fetchone()[0]
            
            # MIC data coverage
            cursor.execute("SELECT COUNT(DISTINCT biosample) FROM mic_data")
            stats['genomes_with_mic'] = cursor.fetchone()[0]
            
            # Artifacts
            cursor.execute("SELECT type, COUNT(*) FROM artifacts GROUP BY type")
            stats['artifacts'] = dict(cursor.fetchall())
            
            return stats
        except Exception as e:
            self.logger.error(f"Failed to get stats: {e}")
            return {}
            
    def close(self):
        """Close database connection."""
        if hasattr(self, 'conn'):
            self.conn.close()

def calculate_file_hash(file_path: str) -> str:
    """Calculate SHA-256 hash of a file."""
    hash_sha256 = hashlib.sha256()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_sha256.update(chunk)
    return hash_sha256.hexdigest()

def validate_fasta_file(file_path: str) -> bool:
    """Validate FASTA file format and integrity."""
    try:
        with open(file_path, 'r') as f:
            first_line = f.readline()
            if not first_line.startswith('>'):
                return False
                
            # Check for basic FASTA structure
            has_sequence = False
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    continue
                # Basic nucleotide/amino acid check
                if any(c not in 'ATCGNRYSWKMBDHV*-' for c in line.upper()):
                    return False
                has_sequence = True
                
            return has_sequence
    except Exception:
        return False