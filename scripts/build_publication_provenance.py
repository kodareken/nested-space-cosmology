#!/usr/bin/env python3
"""Write the deterministic provenance envelope for the curated repository."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "PUBLICATION-PROVENANCE.json"
FRONTIER = ROOT / "results" / "nsc-2-zeta1-recursion-map.json"
IGNORED_PARTS = {
    ".git",
    ".pytest_cache",
    ".ruff_cache",
    ".venv",
    "__pycache__",
    ".build",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def included_files() -> list[Path]:
    paths = []
    for path in ROOT.rglob("*"):
        if not path.is_file() or path == OUTPUT:
            continue
        relative = path.relative_to(ROOT)
        if any(
            part in IGNORED_PARTS or part.endswith(".egg-info")
            for part in relative.parts
        ):
            continue
        if path.name == ".DS_Store" or path.suffix == ".pyc":
            continue
        paths.append(path)
    return sorted(paths, key=lambda item: str(item.relative_to(ROOT)))


def main() -> int:
    record = {
        "schema": "NSC-PUBLICATION-PROVENANCE-v1",
        "public_version": "0.1.0",
        "source_project": "BlackHoles-Infinity",
        "source_commit": "ff2cf2722b966589b98a61accdbb6cee819a58c7",
        "follow_up_source_commit": "5f38712ca01ddd71e715fd265088925a73369aba",
        "nsc3_source_commit": "445d5b069adea8b5641384b0ee02f07fe9cfcc0a",
        "nsc3_introduced_commit": "3490f19164eb9303915db41b8c90eca0f40a836e",
        "source_branch": "doug/commit-branch",
        "export_timestamp_utc": "2026-09-07T00:00:00Z",
        "source_worktree_clean": True,
        "frontier_artifact_id": "NSC-2-ZETA1-RECURSION-MAP",
        "frontier_sha256": sha256(FRONTIER),
        "follow_up_artifact_id": "NSC-2-ZETA1-UNIT-CLOSURE-CHECK",
        "follow_up_sha256": sha256(ROOT / "results" / "nsc-2-zeta1-unit-closure-check.json"),
        "historical_archive_location": "not_public",
        "excluded_category_summaries": [
            "historical FGC configurations, runners, certificates, and tests",
            "raw campaign stores and local scientific data",
            "superseded manuscripts, generated previews, and caches",
            "private plans, host paths, refs, reflogs, and development history",
            "NSC experiments outside the current frontier dependency closure",
        ],
        "included_paths": [
            {
                "path": str(path.relative_to(ROOT)),
                "sha256": sha256(path),
                "bytes": path.stat().st_size,
            }
            for path in included_files()
        ],
    }
    OUTPUT.write_text(
        json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"wrote {OUTPUT.relative_to(ROOT)} with {len(record['included_paths'])} paths")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
