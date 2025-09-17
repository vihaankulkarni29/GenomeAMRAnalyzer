import io
import os
from pathlib import Path
import json
import tempfile
import pytest

from src.priority3.card.rgi_runner import CARDRGIRunner
import sys; import os; sys.path.append(os.path.join(os.path.dirname(__file__), '\\\src')); from configuration_manager import config_manager


class DummyCompleted:
    def __init__(self, returncode=0, stdout="", stderr=""):
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr


def write_text(path: Path, text: str):
    path.write_text(text, encoding="utf-8")


@pytest.fixture()
def tmp_env(tmp_path):
    out = tmp_path / "card_out"
    out.mkdir()
    fasta = tmp_path / "genome.fasta"
    write_text(fasta, ">x\n" + ("ATCG" * 300) + "\n")
    yield out, fasta


def test_missing_fasta(tmp_env):
    out, fasta = tmp_env
    runner = CARDRGIRunner(output_dir=str(out))
    assert runner.run_rgi(str(fasta) + ".missing", "S1") is None


def test_empty_fasta(tmp_env):
    out, fasta = tmp_env
    empty = fasta.with_name("empty.fasta")
    empty.write_text("")
    runner = CARDRGIRunner(output_dir=str(out))
    assert runner.run_rgi(str(empty), "S1") is None


def test_binary_not_found(tmp_env):
    out, fasta = tmp_env
    runner = CARDRGIRunner(rgi_path="__no_such_binary__", output_dir=str(out), retries=1)
    assert runner.run_rgi(str(fasta), "S1") is None


def test_timeout(tmp_env):
    out, fasta = tmp_env

    def slow_exec(*args, **kwargs):
        raise pytest.raises(subprocess.TimeoutExpired, match="timeout")

    # Instead of raising via pytest.raises, simulate TimeoutExpired explicitly
    import subprocess

    def raise_timeout(*args, **kwargs):
        raise subprocess.TimeoutExpired(cmd=args[0], timeout=0.01)

    runner = CARDRGIRunner(output_dir=str(out), timeout_seconds=1, retries=1, executor=raise_timeout)
    assert runner.run_rgi(str(fasta), "S1") is None


def test_happy_path_tabular(tmp_env):
    out, fasta = tmp_env
    calls = {"n": 0}

    def exec_ok(cmd, capture_output=True, text=True, check=True, timeout=None):
        calls["n"] += 1
        # emulate producing an output file
        for idx, token in enumerate(cmd):
            if token == "--output_file":
                out_path = Path(cmd[idx + 1])
                out_path.write_text("Best_Hit_ARO\tother\nacrA\tx\n")
        return DummyCompleted(0, stdout="ok")

    runner = CARDRGIRunner(output_dir=str(out), executor=exec_ok)
    outpath = runner.run_rgi(str(fasta), "S1", output_format="txt")
    assert outpath is not None
    assert Path(outpath).exists()
    assert Path(outpath).read_text().startswith("Best_Hit_ARO")


def test_retry_then_success(tmp_env):
    out, fasta = tmp_env
    state = {"attempt": 0}

    import subprocess

    def flaky_exec(cmd, capture_output=True, text=True, check=True, timeout=None):
        state["attempt"] += 1
        if state["attempt"] == 1:
            raise subprocess.CalledProcessError(1, cmd, stderr="boom")
        # success on second
        for idx, token in enumerate(cmd):
            if token == "--output_file":
                Path(cmd[idx + 1]).write_text("Best_Hit_ARO\tother\nacrA\tx\n")
        return DummyCompleted(0, stdout="ok")

    runner = CARDRGIRunner(output_dir=str(out), retries=2, executor=flaky_exec)
    outpath = runner.run_rgi(str(fasta), "S1", output_format="txt")
    assert outpath is not None


def test_json_validation_and_parse(tmp_env):
    out, fasta = tmp_env
    calls = {"n": 0}

    def exec_json(cmd, capture_output=True, text=True, check=True, timeout=None):
        for idx, token in enumerate(cmd):
            if token == "--output_file":
                Path(cmd[idx + 1]).write_text(json.dumps([{"gene": config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0]}]))
        return DummyCompleted(0, stdout="ok")

    runner = CARDRGIRunner(output_dir=str(out), executor=exec_json)
    outpath = runner.run_rgi(str(fasta), "S1", output_format="json")
    assert outpath is not None
    parsed = runner.parse_rgi_output(outpath)
    assert isinstance(parsed, list)
    assert parsed and parsed[0].get("gene") == config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0]


def test_extract_gene_hits(tmp_env):
    out, fasta = tmp_env
    tab = out / "S1_rgi.txt"
    tab.write_text("Best_Hit_ARO\tcol\nacrA\t1\nacrB\t2\n")
    runner = CARDRGIRunner(output_dir=str(out))
    hits = runner.extract_gene_hits(str(tab), [config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0]]) 
    assert len(hits) == 1 and hits[0]["Best_Hit_ARO"] == config_manager.get_default_genes("rnd_efflux_pumps", "primary")[0]


def test_batch_mixed(tmp_env):
    out, fasta = tmp_env
    ok_fasta = fasta
    bad_fasta = out / "missing.fasta"

    import subprocess

    def exec_ok(cmd, **kwargs):
        for idx, token in enumerate(cmd):
            if token == "--output_file":
                Path(cmd[idx + 1]).write_text("Best_Hit_ARO\tcol\nacrA\t1\n")
        return DummyCompleted(0)

    runner = CARDRGIRunner(output_dir=str(out), executor=exec_ok)
    res = runner.run_batch([str(ok_fasta), str(bad_fasta)], ["S1", "S2"], output_format="txt")
    assert "S1" in res and "S2" not in res


def test_large_batch_stress(tmp_env):
    out, fasta = tmp_env

    # Create 120 fake FASTA files (copy of valid one); simulate every 10th failing
    fastas = []
    sids = []
    for i in range(120):
        fp = out / f"G_{i}.fasta"
        fp.write_text(fasta.read_text())
        fastas.append(str(fp))
        sids.append(f"S{i}")

    import subprocess

    import subprocess

    def exec_mixed(cmd, **kwargs):
        # detect sid from output file name
        out_path = None
        sid = None
        for idx, token in enumerate(cmd):
            if token == "--output_file":
                out_path = Path(cmd[idx + 1])
                sid = out_path.name.split("_rgi")[0]
                break
        assert out_path is not None and sid is not None
        # fail every 10th sample
        idxnum = int(sid[1:])
        if idxnum % 10 == 0:
            raise subprocess.CalledProcessError(1, cmd, stderr="simulated error")
        out_path.write_text("Best_Hit_ARO\tcol\nacrA\t1\n")
        return subprocess.CompletedProcess(cmd, 0, stdout="ok", stderr="")

    runner = CARDRGIRunner(output_dir=str(out), executor=exec_mixed, retries=1)
    results = runner.run_batch(fastas, sids, output_format="txt")

    # Expect ~90% success (since every 10th failed)
    assert len(results) == 108  # 120 - 12 failures
    # spot check an expected success and failure
    assert "S1" in results
    assert "S10" not in results
