import asyncio
from pathlib import Path
import re
import pytest
from aioresponses import aioresponses

from url_to_genomes_workflow import URLToGenomesWorkflow
from ncbi_genome_discovery import GenomeMetadata


@pytest.mark.asyncio
async def test_workflow_end_to_end(tmp_output_dir, monkeypatch, project_root):
    # Use a copy of config with tmp directories
    config_path = project_root / "config" / "snakemake_config.yaml"
    assert config_path.exists()

    # Patch directories in the loaded config at runtime via monkeypatch
    import yaml
    with open(config_path, 'r') as f:
        cfg = yaml.safe_load(f)

    cfg['directories']['genomes'] = str(tmp_output_dir)
    cfg['directories']['logs'] = str(tmp_output_dir / 'logs')
    cfg['directories']['reports'] = str(tmp_output_dir / 'reports')

    # Create a temp config file
    tmp_cfg_path = tmp_output_dir / 'config.yaml'
    with open(tmp_cfg_path, 'w') as f:
        yaml.safe_dump(cfg, f)

    workflow = URLToGenomesWorkflow(str(tmp_cfg_path))

    # Mock discovery to return one genome
    async def mock_discover_from_url(url):
        return [GenomeMetadata(accession="NC_000913.3", organism="E. coli", length=1_000_000, title="complete genome")], "report"

    workflow.discovery_engine.discover_from_url = mock_discover_from_url

    # Mock efetch to return valid FASTA
    fasta = b">NC_000913.3\n" + b"ATGC" * 300 + b"\n"

    with aioresponses() as mocked:
        mocked.get(re.compile(r"https://eutils\.ncbi\.nlm\.nih\.gov/entrez/eutils/efetch\.fcgi.*"), status=200, body=fasta, repeat=True)
        files, report = await workflow.run_complete_workflow("https://www.ncbi.nlm.nih.gov/nuccore/?term=Escherichia+coli")

    # Validate outputs
    assert len(files) == 1
    assert Path(files[0]).exists()
    # Report contains a Ready for RGI section header indirectly via log; in the saved report verify key lines
    assert "URL-to-Genomes Workflow Report" in report
    assert "Ready for RGI Processing:" in report
    
    assert workflow.validate_for_rgi_processing() is True
