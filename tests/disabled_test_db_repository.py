import os
import tempfile
import pytest
from src.priority3.db.repositories import GenomeRepository, GenomeRecord, MICRecord, ArtifactRecord
from datetime import datetime

def test_genome_crud():
    db_fd, db_path = tempfile.mkstemp()
    repo = GenomeRepository(db_path)
    genome = GenomeRecord(
        accession="TEST123",
        organism="Escherichia coli",
        biosample="BS1",
        bioproject="BP1",
        file_path="/tmp/test.fna",
        file_hash="abc123",
        download_date=datetime.now(),
        assembly_level="Complete",
        genome_size=5000000,
        status="downloaded"
    )
    # Add genome
    assert repo.add_genome(genome)
    # Retrieve genome
    g2 = repo.get_genome("TEST123")
    assert g2 is not None
    assert g2.organism == "Escherichia coli"
    # Update genome
    assert repo.update_genome_status("TEST123", "processed")
    g3 = repo.get_genome("TEST123")
    assert g3.status == "processed"
    # Delete genome (simulate by direct SQL for now)
    repo.conn.execute("DELETE FROM genomes WHERE accession = ?", ("TEST123",))
    repo.conn.commit()
    assert repo.get_genome("TEST123") is None
    os.close(db_fd)
    os.remove(db_path)

def test_mic_crud():
    db_fd, db_path = tempfile.mkstemp()
    repo = GenomeRepository(db_path)
    mic = MICRecord(
        biosample="BS1",
        antibiotic="ampicillin",
        mic_value=2.0,
        mic_units="mg/L",
        method="broth",
        source="test",
        validated=True
    )
    # Add MIC (simulate by direct SQL for now)
    repo.conn.execute("INSERT INTO mic_data (biosample, antibiotic, mic_value, mic_units, method, source, validated) VALUES (?, ?, ?, ?, ?, ?, ?)",
        (mic.biosample, mic.antibiotic, mic.mic_value, mic.mic_units, mic.method, mic.source, mic.validated))
    repo.conn.commit()
    # Retrieve MIC
    cur = repo.conn.execute("SELECT * FROM mic_data WHERE biosample = ?", ("BS1",))
    row = cur.fetchone()
    assert row is not None
    assert row[1] == "BS1"
    assert row[2] == "ampicillin"
    os.close(db_fd)
    os.remove(db_path)

def test_artifact_crud():
    db_fd, db_path = tempfile.mkstemp()
    repo = GenomeRepository(db_path)
    artifact = ArtifactRecord(
        type="alignment",
        path="/tmp/align.txt",
        hash="def456",
        size=1234,
        created_date=datetime.now(),
        provenance={"accession": "TEST123"}
    )
    # Add artifact (simulate by direct SQL for now)
    repo.conn.execute("INSERT INTO artifacts (type, path, hash, size, created_date, provenance) VALUES (?, ?, ?, ?, ?, ?)",
        (artifact.type, artifact.path, artifact.hash, artifact.size, artifact.created_date, str(artifact.provenance)))
    repo.conn.commit()
    # Retrieve artifact
    cur = repo.conn.execute("SELECT * FROM artifacts WHERE path = ?", ("/tmp/align.txt",))
    row = cur.fetchone()
    assert row is not None
    assert row[1] == "alignment"
    os.close(db_fd)
    os.remove(db_path)
