"""Produce a clean backup archive of the checkout.

The 2026-10-05 backup of nester-space-cosmology was 9.1 GB and 94,630 files,
of which 14,045 were a virtual environment and the rest included build
caches, pytest caches and macOS metadata. This writes the same tree without
those, so the archive is small enough to look at and quick to verify.

The Git directory is kept. On this machine the branch is fifty commits ahead
of origin, so the archive is one of the few copies of that history.

Usage:
    python backup_lab_tree.py [--root PATH] [--out PATH] [--verify]

Default root is the parent of this repository, default output is
<root>.clean.7z beside it.
"""

from __future__ import annotations

import argparse
import pathlib
import subprocess
import sys

#: Directories never worth archiving.
EXCLUDE_DIRS = {
    ".venv",
    "venv",
    "__pycache__",
    ".pytest_cache",
    ".mypy_cache",
    ".ruff_cache",
    ".build",
    "node_modules",
    ".hypa",
    ".hypa-pi",
    ".codex",
    ".cua-driver",
    "_build",
    "dist",
    "build",
}

#: File suffixes never worth archiving.
EXCLUDE_SUFFIXES = (".pyc", ".pyo", ".pyd")

#: macOS metadata that this tree carried from its original machine.
EXCLUDE_NAMES = {".DS_Store"}

SEVEN_ZIP_CANDIDATES = (
    r"C:/Program Files/7-Zip/7z.exe",
    "7z",
)


def find_seven_zip() -> str:
    for candidate in SEVEN_ZIP_CANDIDATES:
        if pathlib.Path(candidate).exists():
            return candidate
    return "7z"


def excluded(relative: pathlib.PurePath) -> bool:
    for part in relative.parts:
        if part in EXCLUDE_DIRS:
            return True
        if part.startswith("._"):
            return True
        if part in EXCLUDE_NAMES:
            return True
    return relative.name.endswith(EXCLUDE_SUFFIXES)


def survey(root: pathlib.Path) -> tuple[int, int, dict[str, int]]:
    kept = 0
    dropped = 0
    dropped_by: dict[str, int] = {}
    for path in root.rglob("*"):
        rel = path.relative_to(root)
        if excluded(rel):
            if path.is_file():
                dropped += 1
                reason = next(
                    (p for p in rel.parts if p in EXCLUDE_DIRS),
                    "<suffix or metadata>",
                )
                dropped_by[reason] = dropped_by.get(reason, 0) + 1
            continue
        if path.is_file():
            kept += 1
    return kept, dropped, dropped_by


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    default_root = pathlib.Path(__file__).resolve().parents[2].parent
    parser.add_argument("--root", default=str(default_root))
    parser.add_argument("--out", default=None)
    parser.add_argument("--verify", action="store_true",
                        help="run 7z t on the archive after writing it")
    args = parser.parse_args(argv)

    root = pathlib.Path(args.root).resolve()
    if not root.is_dir():
        raise SystemExit(f"root is not a directory: {root}")
    out = pathlib.Path(args.out) if args.out else root.with_suffix(".clean.7z")
    if out.exists():
        raise SystemExit(f"refusing to overwrite {out}")

    kept, dropped, dropped_by = survey(root)
    print(f"root          {root}")
    print(f"files kept    {kept}")
    print(f"files skipped {dropped}")
    for reason, count in sorted(dropped_by.items(), key=lambda kv: -kv[1]):
        print(f"    {count:7d}  {reason}")

    seven = find_seven_zip()
    command = [
        seven, "a", "-t7z", "-mx=5",
        "-xr!.venv", "-xr!__pycache__", "-xr!.pytest_cache",
        "-xr!.mypy_cache", "-xr!.ruff_cache", "-xr!.build",
        "-xr!node_modules", "-xr!._*", "-xr!.DS_Store",
        "-xr!*.pyc", "-xr!*.pyo",
        str(out), str(root),
    ]
    print(f"\nrunning {seven} ...", flush=True)
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        print(result.stdout[-2000:])
        print(result.stderr[-2000:])
        raise SystemExit(f"7z failed with code {result.returncode}")

    size = out.stat().st_size
    print(f"wrote {out}")
    print(f"size  {size / 1e9:.2f} GB")

    if args.verify:
        print("verifying ...", flush=True)
        check = subprocess.run([seven, "t", str(out)],
                               capture_output=True, text=True)
        tail = [line for line in check.stdout.splitlines() if line.strip()][-4:]
        for line in tail:
            print("   ", line)
        return check.returncode
    return 0


if __name__ == "__main__":
    sys.exit(main())