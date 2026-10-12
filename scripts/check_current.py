#!/usr/bin/env python3
"""Check the canonical root; historical evidence has an explicit replay command."""
import argparse
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--local-evidence", action="store_true")
    parser.add_argument("--historical-archive", action="store_true")
    args = parser.parse_args(argv)
    if args.local_evidence or args.historical_archive:
        parser.error("the old archive route is retired; use python scripts/replay_history.py all")
    for name in ("check_lab_snapshot.py", "check_publication.py"):
        result = subprocess.call([sys.executable, str(ROOT / "scripts" / name)], cwd=ROOT)
        if result:
            return result
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
