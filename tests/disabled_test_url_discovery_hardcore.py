import asyncio
from pathlib import Path
import pytest
from aioresponses import aioresponses

from ncbi_genome_discovery import NCBIUrlParser, NCBIGenomeDiscovery, URLBasedGenomeDiscovery


@pytest.mark.parametrize("url,expected", [
    (
        "https://www.ncbi.nlm.nih.gov/nuccore/?term=Escherichia+coli+complete+genome",
        "Escherichia coli complete genome",
    ),
    (
        "https://www.ncbi.nlm.nih.gov/nuccore/?term=(Escherichia+coli+and+macrolide+resistance)+AND+%22Escherichia+coli%22%5Bporgn%3A__txid562%5D+and+complete+genome",
        '(Escherichia coli and macrolide resistance) AND "Escherichia coli"[porgn:__txid562] and complete genome',
    ),
])
def test_parser_extracts_term(url, expected):
    parser = NCBIUrlParser()
    term = parser.parse_search_url(url)
    assert term == expected


@pytest.mark.asyncio
def test_discovery_handles_empty_results(monkeypatch):
    # Mock esearch returning no IDs
    discovery = NCBIGenomeDiscovery(email="test@example.com", max_results=5)

    async def mock_search(_):
        return []

    monkeypatch.setattr(discovery, "_search_genomes", mock_search)
    genomes = asyncio.get_event_loop().run_until_complete(discovery.discover_genomes("anything"))
    assert genomes == []


@pytest.mark.asyncio
async def test_discovery_parses_metadata(monkeypatch):
    discovery = NCBIGenomeDiscovery(email="test@example.com", max_results=5)

    async def mock_search(_):
        return ["1", "2"]

    async def mock_get_batch(ids):
        from ncbi_genome_discovery import GenomeMetadata
        return [
            GenomeMetadata(accession="NC_000913.3", organism="Escherichia coli", length=4641652, title="complete genome"),
            GenomeMetadata(accession="CP000000.1", organism="Escherichia coli", length=1800000, title="chromosome")
        ]

    monkeypatch.setattr(discovery, "_search_genomes", mock_search)
    monkeypatch.setattr(discovery, "_get_batch_metadata", mock_get_batch)

    genomes = await discovery.discover_genomes("Escherichia coli complete genome")
    # Should filter to max_results and quality criteria
    assert all(g.length >= 1_000_000 for g in genomes)
    assert any("complete" in g.title or "chromosome" in g.title for g in genomes)


@pytest.mark.asyncio
async def test_url_based_discovery_pipeline(monkeypatch):
    # End-to-end of URLBasedGenomeDiscovery with mocks underneath
    url_engine = URLBasedGenomeDiscovery(email="test@example.com", max_results=3)

    async def mock_discover(term):
        from ncbi_genome_discovery import GenomeMetadata
        return [
            GenomeMetadata(accession="NC_000913.3", organism="Escherichia coli", length=4641652, title="complete genome"),
        ]

    def mock_report():
        return "report"

    monkeypatch.setattr(url_engine.genome_discovery, "discover_genomes", mock_discover)
    monkeypatch.setattr(url_engine.genome_discovery, "generate_discovery_report", mock_report)

    genomes, report = await url_engine.discover_from_url(
        "https://www.ncbi.nlm.nih.gov/nuccore/?term=Escherichia+coli+complete+genome"
    )
    assert len(genomes) == 1
    assert report == "report"
