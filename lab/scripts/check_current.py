#!/usr/bin/env python3
"""Verify archive preservation and portable current research entrypoints."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
from urllib.parse import unquote

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "archive/pre-closure-2026-09-08/manifest.json"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--local-evidence", action="store_true")
    args = parser.parse_args()
    manifest = json.loads(MANIFEST.read_text())
    for row in manifest["entries"]:
        if row["disposition"] == "retained":
            data = subprocess.check_output(
                ["git", "cat-file", "blob", row["git_blob"]], cwd=ROOT
            )
            if hashlib.sha256(data).hexdigest() != row["sha256"] or len(data) != row["bytes"]:
                raise RuntimeError(f"historical Git blob changed: {row['original_path']}")
            continue
        path = ROOT / row["storage_path"]
        if path.is_symlink() or not path.resolve().is_relative_to(ROOT):
            raise RuntimeError(f"invalid archive path: {path}")
        if hashlib.sha256(path.read_bytes()).hexdigest() != row["sha256"]:
            raise RuntimeError(f"historical bytes changed: {row['original_path']}")
    if args.local_evidence:
        for row in manifest["protected_local_evidence"]:
            path = ROOT / row["path"]
            with path.open("rb") as stream:
                digest = hashlib.file_digest(stream, "sha256").hexdigest()
            if digest != row["sha256"] or path.stat().st_size != row["bytes"]:
                raise RuntimeError(f"local evidence changed: {path}")
    if (ROOT / "PLAN.md").exists():
        raise RuntimeError("the execution plan belongs to the active Codex task, not the repository")
    for name in ("README.md", "docs/active-code-map.md", "docs/claim-ledger.md", "results/README.md"):
        if not (ROOT / name).is_file():
            raise RuntimeError(f"missing current entrypoint: {name}")
    for document in [ROOT / "README.md", *sorted((ROOT / "docs").glob("*.md")), ROOT / "results/README.md"]:
        for target in re.findall(r"\[[^\]]*\]\(([^)]+)\)", document.read_text()):
            if "://" in target or target.startswith("#"):
                continue
            relative = unquote(target.split("#", 1)[0].strip("<>"))
            if not (document.parent / relative).exists():
                # Douglas removed local plans. Keep the authenticated historical
                # derivation byte-identical; only its retired plan link is exempt.
                if (document.relative_to(ROOT).as_posix() == "docs/nsc-compatible-prepared-history.md"
                        and relative == "nsc-compatible-history-plan.md"):
                    receipt = json.loads((ROOT / "results/development/nsc-compatible-prepared-history.json").read_text())
                    expected = receipt["source_sha256"]["docs/nsc-compatible-prepared-history.md"]
                    if hashlib.sha256(document.read_bytes()).hexdigest() == expected:
                        continue
                raise RuntimeError(f"broken current link in {document.name}: {target}")
    print(f"Archive: {len(manifest['entries'])} original file hashes preserved; current entrypoints portable.")
    if args.local_evidence:
        print(f"Protected local evidence: {len(manifest['protected_local_evidence'])} files unchanged.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
