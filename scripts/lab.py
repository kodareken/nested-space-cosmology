#!/usr/bin/env python3
"""Run Python commands in the active lab without mixing publication owners."""
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def environment(root=ROOT):
    env = os.environ.copy()
    lab = root / "lab"
    env["PYTHONPATH"] = os.pathsep.join([str(lab / "src"), str(lab / "scripts")])
    # Only exact pinned source Merkle paths are retained; never add this partial
    # historical object store to the live repository's permanent alternates.
    env["GIT_ALTERNATE_OBJECT_DIRECTORIES"] = str(lab / ".source-history/objects")
    return env


def main():
    if len(sys.argv) == 1:
        print("Usage: python scripts/lab.py <script.py | -m module> [arguments]", file=sys.stderr)
        return 2
    return subprocess.call([sys.executable, *sys.argv[1:]], cwd=ROOT / "lab", env=environment())


if __name__ == "__main__":
    raise SystemExit(main())
