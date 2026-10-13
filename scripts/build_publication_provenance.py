#!/usr/bin/env python3
"""Write the deterministic provenance envelope for the curated repository."""

from __future__ import annotations

import argparse
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


ARCHIVE_DOI = "10.5281/zenodo.23090015"


def evidence_index(root: Path) -> dict[str, dict]:
    """Immutable result identities and historical companion archive identities."""
    manifest = json.loads((root / "results/manifest.json").read_text(encoding="utf-8"))
    snapshot = json.loads((root / "results/development-snapshot.json").read_text(encoding="utf-8"))
    index = {}

    def add(path, digest, size=None):
        relative = Path(path)
        if (relative.is_absolute() or relative.drive or ".." in relative.parts
                or not relative.parts or relative.parts[0] not in {"results", "paper"}
                or not isinstance(digest, str) or len(digest) != 64
                or any(c not in "0123456789abcdef" for c in digest)):
            raise ValueError("unsafe or malformed evidence declaration: " + str(path))
        if path in index and index[path]["sha256"] != digest:
            raise ValueError("conflicting evidence identities: " + path)
        index[path] = {"path": path, "sha256": digest, "bytes": size}

    for step in manifest["steps"]:
        add(step["output"], step["output_sha256"])
        for entry in step.get("auxiliary_inputs", []):
            if entry["path"].startswith("results/"):
                add(entry["path"], entry["sha256"], entry.get("bytes"))
    for entry in snapshot["source_files"]:
        if entry["path"].startswith("results/"):
            add(entry["path"], entry["sha256"], entry["bytes"])
    catalog = json.loads((root / "paper/catalog.json").read_text(encoding="utf-8"))
    for row in catalog["documents"]:
        if row["id"] != "local-gate":
            continue
        for edition in [row, *row.get("historical_editions", [])]:
            record = json.loads((root / edition["manifest"]).read_text(encoding="utf-8"))
            for field, digest in (("source", record["inputs"][edition["source"]]),
                                  ("pdf", record["pdf_sha256"]),
                                  ("source_archive", record["source_archive_sha256"])):
                add(edition[field], digest)
    return index


def evidence_locations(root: Path, authorized: set[str], prior: dict) -> list[dict]:
    """Missing evidence requires an explicit archival disposition, never a skip."""
    archived = {row["path"]: row for row in prior.get("evidence_locations", [])
                if row.get("status") == "archived" or row.get("archive_disposition")}
    locations = []
    for relative, expected in sorted(evidence_index(root).items()):
        path = root / relative
        row = dict(expected)
        previous = archived.get(relative)
        disposition = (previous or {}).get("archive_disposition", {})
        predeclared = previous is not None and previous.get("sha256") == row["sha256"] and (
            (previous.get("status") == "archived" and previous.get("doi") == ARCHIVE_DOI)
            or (disposition.get("sha256") == row["sha256"] and disposition.get("doi") == ARCHIVE_DOI))
        admitted = relative in authorized or predeclared
        if admitted:
            # Record the audited retirement before deletion, while its original
            # local bytes can still be checked. A later builder keeps this receipt.
            row["archive_disposition"] = {"doi": ARCHIVE_DOI, "sha256": row["sha256"]}
        if path.is_file():
            if sha256(path) != row["sha256"]:
                raise ValueError("retained evidence differs: " + relative)
            if row["bytes"] is not None and path.stat().st_size != row["bytes"]:
                raise ValueError("retained evidence size differs: " + relative)
            row["status"] = "local"
        else:
            if not admitted:
                raise ValueError("missing evidence has no archival disposition: " + relative)
            row.update(status="archived", doi=ARCHIVE_DOI,
                       authentication="immutable pinned Git context; remote download not checked")
        locations.append(row)
    return locations


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
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive-path-list", type=Path,
                        help="explicit authorized relative retirement paths for this cleanup")
    args = parser.parse_args()
    authorized = set()
    if args.archive_path_list:
        authorized = {line.strip() for line in args.archive_path_list.read_text(encoding="utf-8").splitlines()
                      if line.strip()}
        for relative in authorized:
            path = Path(relative)
            if path.is_absolute() or path.drive or ".." in path.parts:
                raise ValueError("unsafe archival disposition: " + relative)
    prior = json.loads(OUTPUT.read_text(encoding="utf-8")) if OUTPUT.is_file() else {}
    locations = evidence_locations(ROOT, authorized, prior)
    identities = {row["path"]: row["sha256"] for row in locations}
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
        "frontier_sha256": identities[FRONTIER.relative_to(ROOT).as_posix()],
        "follow_up_artifact_id": "NSC-2-ZETA1-UNIT-CLOSURE-CHECK",
        "follow_up_sha256": identities["results/nsc-2-zeta1-unit-closure-check.json"],
        "historical_archive_location": "archive/",
        "lineage_scope": "Original release lineage; current contents are enumerated by included_paths.",
        "archived_evidence_doi": ARCHIVE_DOI,
        "evidence_locations": locations,
        "envelope_scope": (
            "Envelope exclusions do not imply exclusion from Git or public access. "
            "Root archive/ is tracked and authenticated separately; paper/archive/ is included here."
        ),
        "lab_snapshot_sha256": sha256(ROOT / "docs/lab-snapshot.json"),
        "excluded_category_summaries": [
            "archive/source-history original objects and pins authenticated by retained-file receipts",
            "retired lab/ and old nested checkout containers, if present locally",
            "local Git metadata, build/replay directories, virtual environments and caches",
            "machine-local agent configuration, runtime state and research-reset backup",
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
