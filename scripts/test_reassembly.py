#!/usr/bin/env python3
"""Run the affected atlas, preservation and publication tests, not old campaigns."""
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
PATTERNS = ("test_nsc_atlas*.py", "test_nsc_cosmological_atlas*.py", "test_reassembl*.py",
            "test_replay_history.py", "test_lab_consolidation.py", "test_current_check_routing.py",
            "test_paper*.py", "test_main_paper*.py", "test_publication*.py", "test_claim*.py", "test_v070_release.py",
            "test_provenance*.py", "test_atlas_root_reassembly.py", "test_root_runtime.py",
            "test_reassembly_retirement.py")

if __name__ == "__main__":
    tests = sorted({str(p.relative_to(ROOT)) for pattern in PATTERNS for p in (ROOT/"tests").glob(pattern)
                    if p.name != "test_publication.py"})  # frozen release contract runs in explicit history
    raise SystemExit(subprocess.call([sys.executable, "-m", "pytest", "-q", *tests], cwd=ROOT))
