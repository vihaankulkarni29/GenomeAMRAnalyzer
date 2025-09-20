import asyncio
from pathlib import Path
import re
import pytest
from aioresponses import aioresponses

from genome_downloader import GenomeDownloader
from ncbi_genome_discovery import GenomeMetadata


@pytest.mark.asyncio
async def test_downloader_success(tmp_output_dir):
    downloader = GenomeDownloader(
        output_dir=str(tmp_output_dir),
        email="test@example.com",
        max_concurrent=2,
        retry_attempts=1,
        timeout_seconds=30,
    )

    genomes = [
        GenomeMetadata(accession="NC_000913.3", organism="E. coli", length=1_000_000, title="complete genome")
    ]

    # Well-formed FASTA content
    fasta = b">NC_000913.3\n" + b"ATGC" * 300 + b"\n"

    efetch_pattern = re.compile(r"https://eutils\.ncbi\.nlm\.nih\.gov/entrez/eutils/efetch\.fcgi.*")

    with aioresponses() as mocked:
        mocked.get(efetch_pattern, status=200, body=fasta, repeat=True)
        results = await downloader.download_genomes(genomes)

    result = results["NC_000913.3"]
    assert result.success is True
    assert result.file_path is not None
    assert Path(result.file_path).exists()


@pytest.mark.asyncio
async def test_downloader_retry_then_fail(tmp_output_dir):
    downloader = GenomeDownloader(
        output_dir=str(tmp_output_dir),
        email="test@example.com",
        max_concurrent=1,
        retry_attempts=2,
        timeout_seconds=10,
    )

    genomes = [GenomeMetadata(accession="BAD123", organism="X", length=1_500_000, title="complete genome")]

    with aioresponses() as mocked:
        mocked.get(downloader.efetch_base, status=500)
        results = await downloader.download_genomes(genomes)

    result = results["BAD123"]
    assert result.success is False
    assert result.error_message is not None
