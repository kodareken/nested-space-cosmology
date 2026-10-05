"""Make the recursive_horizons package importable for every lab test.

Several modules under lab/src are only reachable when that directory is on
the path, and tests that import them otherwise fail at collection for a
reason that has nothing to do with the test.
"""
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
for _path in (str(LAB / "src"), str(LAB / "scripts")):
    if _path not in sys.path:
        sys.path.insert(0, _path)
