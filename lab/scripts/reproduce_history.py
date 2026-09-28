#!/usr/bin/env python3
"""Run selected non-campaign reproductions in the original Git checkpoint."""

from __future__ import annotations

import argparse
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile

ROOT = Path(__file__).resolve().parents[1]
ANCHOR = "5f38712ca01ddd71e715fd265088925a73369aba"
COMMANDS = {
    "scale-closure": [sys.executable, "-B", "scripts/check_nsc_scale_closure.py", "--check"],
    "core-identities": [sys.executable, "-B", "scripts/reproduce_core.py", "--output", "build/historical-core.json"],
    "repository-check": [sys.executable, "-B", "scripts/check_repo.py"],
}


def materialize(destination: Path) -> None:
    subprocess.run(["git", "clone", "--quiet", "--no-checkout", "--shared", str(ROOT), str(destination)], check=True)
    subprocess.run(["git", "update-ref", "--no-deref", "HEAD", ANCHOR], cwd=destination, check=True)
    subprocess.run(["git", "read-tree", ANCHOR], cwd=destination, check=True)
    data = subprocess.check_output(["git", "--no-replace-objects", "archive", ANCHOR], cwd=ROOT)
    with tarfile.open(fileobj=io.BytesIO(data)) as archive:
        # Historical PLAN.md is a machine-specific symlink. Keep the original
        # target as text, never follow it or import the host's private plan.
        for member in archive:
            target = destination / member.name
            if not target.resolve().is_relative_to(destination.resolve()):
                raise RuntimeError("unsafe archive member")
            if member.isdir():
                target.mkdir(parents=True, exist_ok=True)
            elif member.issym():
                target.parent.mkdir(parents=True, exist_ok=True)
                target.with_name(target.name + ".symlink-target.txt").write_text(member.linkname)
            elif member.isfile():
                target.parent.mkdir(parents=True, exist_ok=True)
                source = archive.extractfile(member)
                assert source is not None
                target.write_bytes(source.read())
                target.chmod(member.mode)
            else:
                raise RuntimeError("unsupported archive member")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=COMMANDS)
    parser.add_argument("--describe", action="store_true")
    args = parser.parse_args()
    if args.describe:
        print(json.dumps({"checkpoint": ANCHOR, "command": COMMANDS[args.case],
                          "raw_stores_available": False,
                          "historical_PLAN_symlink_materialized": False}))
        return 0
    with tempfile.TemporaryDirectory(prefix="nsc-history-") as directory:
        work = Path(directory)
        materialize(work)
        (work / "build").mkdir(exist_ok=True)
        environment = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", PYTHONPATH=str(work / "src"))
        result = subprocess.run(COMMANDS[args.case], cwd=work, env=environment, check=False)
        print(json.dumps({"checkpoint": ANCHOR, "case": args.case, "exit_code": result.returncode,
                          "historical_failure_not_repaired": result.returncode != 0}))
        return result.returncode


if __name__ == "__main__":
    raise SystemExit(main())
