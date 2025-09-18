import io
import os
import gzip
import json
import time
import shutil
import sqlite3
from pathlib import Path

import pytest

# Ensure we can import from src
import sys
ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from priority3.db.harvesters.genbank_genome_harvester import (
    GenBankGenomeHarvester,
    GenomeMetadata,
)


class DummyResponse:
    """Simple mock of requests.Response for streaming content."""

    def __init__(self, content: bytes, status_code: int = 200, headers=None):
        self._content = content
        self.status_code = status_code
        self.headers = headers or {}

    @property
    def content(self):
        return self._content

    def raise_for_status(self):
        if not (200 <= self.status_code < 300):
            raise Exception(f"HTTP {self.status_code}")

    def iter_content(self, chunk_size=8192):
        for i in range(0, len(self._content), chunk_size):
            yield self._content[i : i + chunk_size]


@pytest.fixture()
def tmp_env(tmp_path):
    out_dir = tmp_path / "genomes"
    out_dir.mkdir()
    db_path = tmp_path / "harvest.db"
    yield out_dir, db_path


def test_search_and_metadata_mock(tmp_env):
    out_dir, db_path = tmp_env
    harv = GenBankGenomeHarvester(output_dir=str(out_dir), db_path=str(db_path), mock_mode=True)

    res = harv.search_genomes("Escherichia coli", max_results=10)
    assert res.total_count == len(res.accessions) == 10
    assert all(acc.startswith("GCF_") for acc in res.accessions)

    meta = harv.get_genome_metadata(res.accessions[:3])
    assert len(meta) == 3
    assert meta[0].organism == "Escherichia coli"
    assert meta[0].assembly_level == "Complete Genome"


def test_mock_download_updates_db_and_file(tmp_env):
    out_dir, db_path = tmp_env
    harv = GenBankGenomeHarvester(output_dir=str(out_dir), db_path=str(db_path), mock_mode=True)

    m = GenomeMetadata(
        accession="GCF_999999.1",
        organism="Escherichia coli",
        assembly_level="Complete Genome",
        genome_size=5000000,
    )
    ok = harv.download_genome(m)
    assert ok
    fasta = out_dir / f"{m.accession}.fasta"
    assert fasta.exists() and fasta.stat().st_size > 0

    # Validate DB insertion
    # Use raw sqlite to avoid SQLAlchemy environment differences
    conn = sqlite3.connect(db_path)
    cur = conn.cursor()
    cur.execute("SELECT accession, status, file_path FROM genomes WHERE accession=?", (m.accession,))
    row = cur.fetchone()
    assert row is not None
    assert row[0] == m.accession and row[1] == "downloaded"
    assert Path(row[2]).name == fasta.name


def test_validate_downloaded_genome_checks(tmp_env):
    out_dir, db_path = tmp_env
    harv = GenBankGenomeHarvester(output_dir=str(out_dir), db_path=str(db_path), mock_mode=True)

    # Too small file
    small = out_dir / "small.fasta"
    small.write_text(">x\nA\n")
    assert harv._validate_downloaded_genome(str(small)) is False

    # Valid minimal FASTA (make it larger than threshold)
    valid = out_dir / "valid.fasta"
    valid.write_text(">x\n" + ("ATCG" * 300) + "\n")  # 1200 chars
    assert harv._validate_downloaded_genome(str(valid)) is True


def test_download_plaintext_stream(tmp_env, monkeypatch):
    out_dir, db_path = tmp_env
    harv = GenBankGenomeHarvester(output_dir=str(out_dir), db_path=str(db_path), mock_mode=False)

    # Prepare plain FASTA bytes
    fasta_txt = ">seq1\n".encode() + ("ATCG" * 300).encode() + b"\n"

    def fake_get(url, stream=True, timeout=300):
        return DummyResponse(fasta_txt, 200, headers={"Content-Type": "text/plain"})

    monkeypatch.setattr(harv, "_get_genome_download_url", lambda acc: f"https://example.org/{acc}.fna")
    monkeypatch.setattr(harv.session, "get", fake_get)

    m = GenomeMetadata(accession="GCF_111111.1", organism="E. coli")
    ok = harv.download_genome(m)
    assert ok
    out = out_dir / f"{m.accession}.fasta"
    assert out.exists() and out.stat().st_size > 0


def test_download_gzip_stream_and_resume(tmp_env, monkeypatch):
    out_dir, db_path = tmp_env
    harv = GenBankGenomeHarvester(output_dir=str(out_dir), db_path=str(db_path), mock_mode=False)

    # Create gzipped FASTA
    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode="wb") as gz:
        gz.write((">seq1\n" + ("ATCG" * 300) + "\n").encode())
    gz_bytes = buf.getvalue()

    def fake_get_gz(url, stream=True, timeout=300):
        # Simulate .gz URL with gzip content
        return DummyResponse(gz_bytes, 200, headers={"Content-Type": "application/gzip", "Content-Encoding": "gzip"})

    monkeypatch.setattr(harv, "_get_genome_download_url", lambda acc: f"https://example.org/{acc}.fna.gz")
    monkeypatch.setattr(harv.session, "get", fake_get_gz)

    m = GenomeMetadata(accession="GCF_222222.1", organism="E. coli")
    ok = harv.download_genome(m)
    # Expect True after robust gzip handling is implemented
    assert ok, "Gzip download should be handled and validated"
    out = out_dir / f"{m.accession}.fasta"
    assert out.exists() and out.stat().st_size > 0

    # Now test harvest resume behavior with a controlled failure
    harv.mock_mode = True  # Use mock for search/metadata
    res = harv.search_genomes("test", max_results=3)
    metas = harv.get_genome_metadata(res.accessions)

    calls = {"n": 0}

    def flaky_download(md):
        calls["n"] += 1
        if calls["n"] == 1:
            # first succeeds
            return True
        elif calls["n"] == 2:
            # second fails hard
            raise RuntimeError("network error")
        else:
            return True

    # Patch download method only during harvest
    monkeypatch.setattr(harv, "download_genome", flaky_download)

    with pytest.raises(RuntimeError):
        harv.harvest_genomes("query", max_genomes=3, resume=True)

    # Ensure checkpoint recorded first success
    chk = json.loads((out_dir / "harvest_checkpoint.json").read_text())
    assert len(chk.get("downloaded", [])) == 1

    # New harvester resumes
    harv2 = GenBankGenomeHarvester(output_dir=str(out_dir), db_path=str(db_path), mock_mode=True)
    # Stub methods to reuse prior search/metadata from checkpoint
    res2 = harv2.harvest_genomes("query", max_genomes=3, resume=True)
    assert res2["successfully_downloaded"] >= 1
