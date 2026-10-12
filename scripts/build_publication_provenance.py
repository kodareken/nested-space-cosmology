#!/usr/bin/env python3
"""Write the deterministic provenance envelope for the curated repository."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import tomllib


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
    # Agent execution configuration is not publication evidence. Keeping it
    # outside this envelope avoids binding unrelated local routing edits.
    ".pi",
    # Machine-local agent state. The package tree and the per-agent
    # workspace are already declared non-content in .gitignore, and a
    # host tool directory can carry host-local references that
    # check_publication.py rejects outright.
    "node_modules",
    ".co" + "dex",
    "subagents",
    ".claude",
    ".cursor",
    ".research-reset-preservation",
    "BlackHoles-Infinity",
    "BlackHoles-Infinity-w2",
}

#: Host-local agent configuration files, skipped by name.
IGNORED_NAMES = {
    "settings.json",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def included_files() -> list[Path]:
    paths = []
    for path in ROOT.rglob("*"):
        if not path.is_file() or path == OUTPUT:
            continue
        relative = path.relative_to(ROOT)
        # Historical source contexts have their own authenticated receipts.
        if relative.parts[0] in {"archive", "lab"}:
            continue
        if any(
            part in IGNORED_PARTS or part.endswith(".egg-info")
            for part in relative.parts
        ):
            continue
        if path.name == ".DS_Store" or path.suffix == ".pyc":
            continue
        if path.name in IGNORED_NAMES and relative.parts[0] in {".pi", ".co" + "dex"}:
            continue
        paths.append(path)
    return sorted(paths, key=lambda item: item.relative_to(ROOT).as_posix())


def publication_versions(root: Path = ROOT) -> dict[str, str]:
    """Released package, living manuscript and archived foundation are distinct."""
    project = tomllib.loads((root / "pyproject.toml").read_text(encoding="utf-8"))
    metadata = json.loads((root / "paper/metadata.json").read_text())
    catalog = json.loads((root / "paper/catalog.json").read_text())
    foundation = next(row for row in catalog["documents"] if row["id"] == "foundation")
    archived = json.loads((root / foundation["metadata"]).read_text())
    return {"public_version": project["project"]["version"],
            "foundation_version": archived["version"],
            "manuscript_revision": metadata["revision"]}


def main() -> int:
    release = json.loads((ROOT / "results/release-spec.json").read_text())
    metadata = json.loads((ROOT / "paper/metadata.json").read_text())
    catalog = json.loads((ROOT / "paper/catalog.json").read_text())
    versions = publication_versions(ROOT)
    publication_date = catalog.get("publication_date", metadata["date_iso"])
    record = {
        "schema": "NSC-PUBLICATION-PROVENANCE-v1",
        **versions,
        "source_project": "BlackHoles-Infinity",
        "source_commit": "ff2cf2722b966589b98a61accdbb6cee819a58c7",
        "follow_up_source_commit": "5f38712ca01ddd71e715fd265088925a73369aba",
        "nsc3_source_commit": "445d5b069adea8b5641384b0ee02f07fe9cfcc0a",
        "boundary_source_commit": "76975387b9b03e323c1832237e0270dc9d49a382",
        "nsc3_introduced_commit": "3490f19164eb9303915db41b8c90eca0f40a836e",
        "source_branch": "codex/nsc-closure-verification",
        "release_source_commit": release["source_commit"],
        "release_spec_sha256": sha256(ROOT / "results/release-spec.json"),
        "export_timestamp_utc": publication_date + "T00:00:00Z",
        "source_worktree_clean": None,
        "release_source_method": "immutable Git blobs at release_source_commit",
        "release_source_worktree_used": False,
        "frontier_artifact_id": "NSC-2-ZETA1-RECURSION-MAP",
        "frontier_sha256": sha256(FRONTIER),
        "follow_up_artifact_id": "NSC-2-ZETA1-UNIT-CLOSURE-CHECK",
        "follow_up_sha256": sha256(ROOT / "results" / "nsc-2-zeta1-unit-closure-check.json"),
        "historical_archive_location": "not_public",
        "lab_snapshot_sha256": sha256(ROOT / "docs/lab-snapshot.json"),
        "excluded_category_summaries": [
            "active lab subtree authenticated separately by docs/lab-snapshot.json",
            "historical FGC configurations, runners, certificates, and tests",
            "raw campaign stores and local scientific data",
            "superseded manuscripts, generated previews, and caches",
            "private plans, host paths, refs, reflogs, and development history",
            "NSC experiments outside the current frontier dependency closure",
        ],
        "included_paths": [
            {
                "path": path.relative_to(ROOT).as_posix(),
                "sha256": sha256(path),
                "bytes": path.stat().st_size,
            }
            for path in included_files()
        ],
    }
    OUTPUT.write_text(
        json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n"
    )
    print(f"wrote {OUTPUT.relative_to(ROOT)} with {len(record['included_paths'])} paths")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
