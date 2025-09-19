import pytest
import pandas as pd
from pathlib import Path
import src.abricate_to_coords as abricate_to_coords

@pytest.fixture
def sample_tsv(tmp_path):
    asset = tmp_path / "sample_abricate_vfdb.tsv"
    asset.write_text("contig_id\tstart\tend\tstrand\tgene_name\tcut_off\tbest_hit_aro\ncontigA\t100\t200\t+\tvfgA\t99\tARO:1234\n")
    return asset

def test_abricate_to_coords_transformation(sample_tsv, tmp_path):
    out_csv = tmp_path / "coords.csv"
    abricate_to_coords.main([str(sample_tsv), str(out_csv)])
    df = pd.read_csv(out_csv)
    expected_cols = ["contig_id", "start", "end", "strand", "gene_name", "cut_off", "best_hit_aro"]
    for col in expected_cols:
        assert col in df.columns
    # Check key data matches
    assert df.loc[0, "contig_id"] == "contigA"
    assert df.loc[0, "start"] == 100
    assert df.loc[0, "end"] == 200
