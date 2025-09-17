import os
import sys
from pathlib import Path
import pytest

# Ensure src/core is importable
PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_CORE = PROJECT_ROOT / 'src' / 'core'
if str(SRC_CORE) not in sys.path:
    sys.path.insert(0, str(SRC_CORE))

@pytest.fixture(scope='session')
def project_root() -> Path:
    return PROJECT_ROOT

@pytest.fixture(scope='session')
def tmp_output_dir(tmp_path_factory):
    d = tmp_path_factory.mktemp('url_genomes')
    return d
