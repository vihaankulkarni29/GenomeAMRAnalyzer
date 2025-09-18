import asyncio
import time
import re
import pytest
from aioresponses import aioresponses

from genome_downloader import GenomeDownloader
from ncbi_genome_discovery import GenomeMetadata


@pytest.mark.asyncio
async def test_parallel_download_performance(tmp_output_dir):
    downloader = GenomeDownloader(
        output_dir=str(tmp_output_dir),
        email="test@example.com",
        max_concurrent=4,
        retry_attempts=1,
        timeout_seconds=30,
    )

    # Create 8 fake genomes
    genomes = [
        GenomeMetadata(accession=f"ACC{i}", organism="E. coli", length=1_000_000, title="complete genome")
        for i in range(8)
    ]

    # 30KB-like fasta
    fasta = b">ACC\n" + (b"ATGC" * 7000) + b"\n"

    with aioresponses() as mocked:
        mocked.get(re.compile(r"https://eutils\.ncbi\.nlm\.nih\.gov/entrez/eutils/efetch\.fcgi.*"), status=200, body=fasta, repeat=True)
        start = time.time()
        await downloader.download_genomes(genomes)
        elapsed = time.time() - start

    # With concurrency 4 and small payloads, 8 downloads should finish well under 10s
    assert elapsed < 10
    assert len(downloader.get_successful_downloads()) == 8
