import pytest
from pytest_mock import MockerFixture
from pathlib import Path
import src.abricate_runner as abricate_runner

@pytest.fixture
def mock_input_dir(tmp_path):
    # Create mock input directory with fake FASTA files
    fasta1 = tmp_path / "genome1.fasta"
    fasta1.write_text(">contig1\nATGC")
    return tmp_path

def test_run_abricate_main_invokes_subprocess(mocker: MockerFixture, mock_input_dir):
    # Patch subprocess.run to avoid running actual abricate
    mock_run = mocker.patch("subprocess.run")
    mock_run.return_value.stdout = "MOCK_TSV_OUTPUT"
    # Call main with mock input and db
    argv = [str(mock_input_dir), "--db", "vfdb"]
    abricate_runner.main(argv)
    # Check subprocess.run called with correct args
    expected_cmd = ["abricate", "--db", "vfdb", "--nopath", str(mock_input_dir / "genome1.fasta")]
    mock_run.assert_called()
    # Check output filename logic
    output_files = abricate_runner.find_fastas(mock_input_dir)
    for fasta in output_files:
        out_name = abricate_runner.infer_genome_id(fasta) + "_vfdb.tsv"
        assert out_name.startswith("genome1_vfdb.tsv")
