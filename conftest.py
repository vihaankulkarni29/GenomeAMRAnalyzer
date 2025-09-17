"""
Pytest configuration: ensure repository modules are importable in tests.

This adds the project root, 'src', and 'src/priority3' to sys.path so tests
can import packages like 'metadata', 'analysis', and top-level modules without
per-test sys.path hacks.
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SRC = ROOT / "src"
PRIORITY3 = SRC / "priority3"

for p in (ROOT, SRC, PRIORITY3):
    ps = str(p)
    if ps not in sys.path:
        sys.path.insert(0, ps)
