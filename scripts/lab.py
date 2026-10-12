#!/usr/bin/env python3
"""Compatibility launcher for root Python commands; historical replay is explicit."""
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def environment(root=ROOT):
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join([str(root / "src"), str(root / "scripts"),
        *[p for p in env.get("PYTHONPATH", "").split(os.pathsep) if "site-packages" in p]])
    # Only exact pinned source Merkle paths are retained; never add this partial
    # historical object store to the live repository's permanent alternates.
    env["GIT_ALTERNATE_OBJECT_DIRECTORIES"] = str(root / "archive/source-history/objects")
    return env


def main():
    if len(sys.argv) == 1:
        print("Usage: python scripts/lab.py <script.py | -m module> [arguments]", file=sys.stderr)
        return 2
    return subprocess.call([sys.executable, *sys.argv[1:]], cwd=ROOT, env=environment())


if __name__ == "__main__":
    raise SystemExit(main())
