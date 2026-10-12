#!/usr/bin/env python3
"""Plan and preserve the relocated research tree without overwriting owners.

The receipt authenticates original bytes and explicitly records duplicate and
variant destinations. Applying it moves only files inside the named workspace;
existing different files are never replaced. It does not delete directories,
rewrite scientific records, stage Git changes, or retire a worktree.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import zipfile


ROOT = Path(__file__).resolve().parents[1]
BASELINE = "04a0cb9378977772243ebf7c9be3c74182b030f2"
ARCHIVE = Path("archive/reassembly-20261012")
MANIFEST = Path("results/reassembly/manifest.json")
ACTIVE = {"src", "scripts", "tests", "docs", "results"}
HELPERS = (
    ".atlas-v2-finalize.py", ".atlas-v3-review.py", ".atlas-v5-inspect.py",
    ".causal-observer-summary.py", ".connection-slice-summary.py",
    ".moving-wall-summary.py",
)


def confined(root: Path, relative: str | Path) -> Path:
    """Resolve before any move and reject escapes or filesystem indirection."""
    root = root.resolve()
    value = root / relative
    resolved = value.resolve()
    if not resolved.is_relative_to(root) or resolved == root:
        raise ValueError(f"path leaves the workspace: {relative}")
    if value.is_symlink() or value.is_junction():
        raise ValueError(f"refuse filesystem indirection: {relative}")
    return resolved


def digest(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def runtime(path: Path) -> bool:
    return any(p in {"__pycache__", ".pytest_cache", ".pytest-tmp"}
               or p.endswith(".egg-info") or p.startswith("pytest-of-")
               for p in path.parts)


def plan(root: Path = ROOT) -> dict:
    rows = []

    def add(source: Path, destination: Path, origin: str, logical: str,
            *, copy: bool = False, disposition: str = "retained") -> None:
        source_path = confined(root, source)
        destination_path = confined(root, destination)
        sha = digest(source_path)
        action = "copy" if copy else "move"
        if destination_path.exists():
            if not destination_path.is_file() or digest(destination_path) != sha:
                raise ValueError(f"unresolved destination collision: {destination}")
            action = "existing_identical"
        rows.append({"origin": origin, "original_path": logical,
                     "source": source.as_posix(), "destination": destination.as_posix(),
                     "sha256": sha, "bytes": source_path.stat().st_size,
                     "action": action, "disposition": disposition})

    lab = confined(root, "BlackHoles-Infinity/lab")
    for source in sorted(lab.rglob("*")):
        if not source.is_file():
            continue
        rel = source.relative_to(lab)
        if rel.parts[0] == ".source-history":
            target = Path("archive/source-history") / Path(*rel.parts[1:])
            disposition = "immutable_source_proof"
        elif rel.parts[0] in ACTIVE and not runtime(rel):
            target, disposition = rel, "active_owner"
            existing = confined(root, target)
            if existing.exists() and digest(existing) != digest(source):
                target = ARCHIVE / "lab-variants" / rel
                disposition = "preserved_variant_root_owner_kept"
        else:
            target = ARCHIVE / "lab-context" / rel
            disposition = "historical_context" if not runtime(rel) else "runtime_preserved"
        add(source.relative_to(root), target, "lab", "lab/" + rel.as_posix(),
            disposition=disposition)

    receipt_path = root / ".research-reset-preservation/20261010T112227Z/BlackHoles-Infinity-w2/receipt.json"
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    for entry in sorted(receipt["files"], key=lambda item: item["path"]):
        rel = Path(entry["path"])
        source = Path("BlackHoles-Infinity-w2") / rel
        live = confined(root, source)
        if digest(live) != entry["sha256"] or live.stat().st_size != entry["bytes"]:
            raise ValueError(f"w2 has changed since its preservation receipt: {rel}")
        add(source, ARCHIVE / "w2/files" / rel, "w2", rel.as_posix(), copy=True,
            disposition="runtime_preserved" if runtime(rel) else "exploratory_snapshot")
        active = Path(*rel.parts[1:]) if rel.parts[0] == "lab" else rel
        if runtime(rel) or active.parts[0] not in ACTIVE | {"paper"}:
            continue
        if any(row["destination"] == active.as_posix() for row in rows) or (root / active).exists():
            continue
        add(source, active, "w2", rel.as_posix(), copy=True,
            disposition="exploratory_candidate_not_admitted_result")
    add(receipt_path.relative_to(root), ARCHIVE / "w2/receipt.json", "w2", "receipt.json", copy=True)

    for directory in sorted(root.glob(".atlas-*-preview*")):
        for source in sorted(directory.rglob("*")):
            if source.is_file():
                rel = source.relative_to(root)
                add(rel, ARCHIVE / "previews-and-helpers" / rel, "preview", rel.as_posix(),
                    disposition="superseded_preview")
    for name in HELPERS:
        if (root / name).exists():
            add(Path(name), ARCHIVE / "previews-and-helpers" / name, "helper", name,
                disposition="retired_execution_helper")
    destinations = [row["destination"] for row in rows]
    if len(destinations) != len(set(destinations)):
        raise ValueError("two original entries target the same destination")
    return {
        "schema": "NSC-REASSEMBLY-v1", "baseline_commit": BASELINE,
        "canonical_root": ".", "original_layout": "lab/* under relocated BlackHoles-Infinity/lab",
        "historical_replay": "Original commit paths materialize under ignored .build/replay; records are not rebound.",
        "w2_base_commit": receipt["head"],
        "rows": rows,
        "counts": dict(sorted(Counter(row["action"] for row in rows).items())),
        "dispositions": dict(sorted(Counter(row["disposition"] for row in rows).items())),
    }


def apply(record: dict, root: Path = ROOT) -> dict:
    """Apply an authenticated plan, resumably, without replacing any file."""
    completed = Counter()
    for row in record["rows"]:
        source = confined(root, row["source"])
        destination = confined(root, row["destination"])
        if destination.exists():
            if not destination.is_file() or digest(destination) != row["sha256"]:
                raise ValueError(f"destination drift: {row['destination']}")
            completed["already_preserved"] += 1
            continue
        if row["action"] == "existing_identical":
            raise ValueError(f"former identical owner disappeared: {row['destination']}")
        if not source.is_file() or digest(source) != row["sha256"]:
            raise ValueError(f"source drift: {row['source']}")
        destination.parent.mkdir(parents=True, exist_ok=True)
        if row["action"] == "copy":
            shutil.copy2(source, destination)
        elif row["action"] == "move":
            # Both final absolute paths were resolved and confined above.
            shutil.move(str(source), str(destination))
        else:
            raise ValueError(f"unknown action: {row['action']}")
        if digest(destination) != row["sha256"]:
            raise ValueError(f"preservation check failed: {row['destination']}")
        completed[row["action"]] += 1
    return dict(completed)


def original_bytes_available(record: dict, row: dict, root: Path) -> bool:
    if row["origin"] == "w2":
        snapshot = confined(root, ARCHIVE / "w2/files" / row["original_path"])
        return snapshot.is_file() and digest(snapshot) == row["sha256"]
    if row["origin"] == "lab":
        result = subprocess.run(["git", "show", record["baseline_commit"] + ":" + row["original_path"]],
                                cwd=root, capture_output=True, check=False)
        return result.returncode == 0 and hashlib.sha256(result.stdout).hexdigest() == row["sha256"]
    return False


def record_successors(record: dict, paths: list[str], root: Path = ROOT) -> None:
    """Declare deliberate current edits without replacing their original identity."""
    requested = {Path(path).as_posix() for path in paths}
    found = set()
    for row in record["rows"]:
        if row["destination"] not in requested:
            continue
        current = confined(root, row["destination"])
        if not current.is_file():
            raise ValueError(f"successor is missing: {row['destination']}")
        sha = digest(current)
        if sha != row["sha256"]:
            if not original_bytes_available(record, row, root):
                raise ValueError(f"original source unavailable for successor: {row['original_path']}")
            row["current_sha256"] = sha
            row["current_bytes"] = current.stat().st_size
            row["successor_reason"] = "Explicit canonical-root migration or current documentation edit; original identity retained."
        found.add(row["destination"])
    for destination in sorted(requested - found):
        variant = next((row for row in record["rows"]
                        if row["origin"] == "lab" and row["original_path"] == "lab/" + destination), None)
        if variant is None:
            raise ValueError("not a migrated destination or merged owner: " + destination)
        original = subprocess.run(["git", "show", record["baseline_commit"] + ":" + destination],
                                  cwd=root, capture_output=True, check=False)
        if original.returncode:
            raise ValueError("original root owner unavailable: " + destination)
        current = confined(root, destination)
        merged = {"destination": destination, "current_sha256": digest(current),
                  "current_bytes": current.stat().st_size,
                  "original_root_path": destination,
                  "original_root_sha256": hashlib.sha256(original.stdout).hexdigest(),
                  "preserved_lab_variant": variant["destination"],
                  "preserved_lab_sha256": variant["sha256"]}
        record.setdefault("merges", {})[destination] = merged


def pack_local_storage(record: dict, root: Path = ROOT) -> dict:
    """Preserve ignored/runtime bytes in one deterministic, portable ZIP."""
    destinations = [row["destination"] for row in record["rows"]]
    ignored_run = subprocess.run(["git", "check-ignore", "--no-index", "--stdin", "-z"],
        input=("\0".join(destinations) + "\0").encode(), cwd=root, capture_output=True)
    if ignored_run.returncode not in (0, 1):
        raise RuntimeError(ignored_run.stderr.decode())
    ignored = {path.decode() for path in ignored_run.stdout.split(b"\0") if path}
    rows = sorted((row for row in record["rows"]
                   if row["destination"] in ignored or row["disposition"] == "runtime_preserved"),
                  key=lambda row: row["destination"])
    relative = ARCHIVE / "local-runtime-and-ignored.zip"
    destination = confined(root, relative)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(destination, "x", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for row in rows:
            source = confined(root, row["destination"])
            raw = source.read_bytes()
            if hashlib.sha256(raw).hexdigest() != row.get("current_sha256", row["sha256"]):
                raise ValueError("packed source drift: " + row["destination"])
            info = zipfile.ZipInfo(row["destination"], date_time=(1980, 1, 1, 0, 0, 0))
            info.create_system = 3
            info.external_attr = 0o100644 << 16
            info.compress_type = zipfile.ZIP_DEFLATED
            archive.writestr(info, raw, compress_type=zipfile.ZIP_DEFLATED, compresslevel=9)
            row["storage_archive"] = relative.as_posix()
            row["archive_member"] = row["destination"]
    receipt = {"path": relative.as_posix(), "sha256": digest(destination),
               "bytes": destination.stat().st_size, "members": len(rows)}
    record.setdefault("storage_archives", []).append(receipt)
    return receipt


def check_storage_archives(record: dict, root: Path) -> set[str]:
    """Authenticate every packed member even when its local copy exists."""
    verified = set()
    seen = set()
    for receipt in record.get("storage_archives", []):
        if receipt["path"] in seen:
            raise ValueError("duplicate storage archive")
        seen.add(receipt["path"])
        path = confined(root, receipt["path"])
        if (not path.is_file() or path.stat().st_size != receipt["bytes"]
                or digest(path) != receipt["sha256"]):
            raise ValueError("storage archive mismatch: " + receipt["path"])
        rows = [row for row in record["rows"] if row.get("storage_archive") == receipt["path"]]
        expected = {row["archive_member"] for row in rows}
        with zipfile.ZipFile(path) as archive:
            names = archive.namelist()
            if len(names) != len(set(names)) or set(names) != expected or len(names) != receipt["members"]:
                raise ValueError("storage archive member inventory mismatch")
            for row in rows:
                confined(root, row["archive_member"])
                raw = archive.read(row["archive_member"])
                if hashlib.sha256(raw).hexdigest() != row.get("current_sha256", row["sha256"]):
                    raise ValueError("storage archive member mismatch: " + row["archive_member"])
                verified.add(row["destination"])
    if any(row.get("storage_archive") not in seen for row in record["rows"] if "storage_archive" in row):
        raise ValueError("undeclared storage archive")
    return verified


def check(record: dict, root: Path = ROOT) -> dict:
    """Authenticate every current destination and every declared original separately."""
    verified = 0
    historical = 0
    packed = check_storage_archives(record, root)
    packed_only = 0
    errors = []
    for row in record["rows"]:
        target = confined(root, row["destination"])
        expected = row.get("current_sha256", row["sha256"])
        if not target.exists() and row["destination"] in packed:
            packed_only += 1
        elif not target.is_file() or digest(target) != expected:
            errors.append(row["destination"])
            continue
        if "current_sha256" in row:
            if not original_bytes_available(record, row, root):
                errors.append("original:" + row["original_path"])
                continue
            historical += 1
        verified += 1
    for merge in record.get("merges", {}).values():
        current = confined(root, merge["destination"])
        variant = confined(root, merge["preserved_lab_variant"])
        original = subprocess.run(["git", "show", record["baseline_commit"] + ":" + merge["original_root_path"]],
                                  cwd=root, capture_output=True, check=False)
        if (not current.is_file() or digest(current) != merge["current_sha256"]
                or not variant.is_file() or digest(variant) != merge["preserved_lab_sha256"]
                or original.returncode
                or hashlib.sha256(original.stdout).hexdigest() != merge["original_root_sha256"]):
            errors.append("merge:" + merge["destination"])
    if errors:
        raise ValueError("original bytes unavailable: " + ", ".join(errors))
    return {"preserved_current": verified, "declared_originals_verified": historical,
            "merged_owners_verified": len(record.get("merges", {})),
            "packed_members_verified": len(packed), "packed_only_destinations": packed_only,
            "total": verified}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=MANIFEST)
    parser.add_argument("--apply", action="store_true")
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--record-successors", nargs="+", metavar="PATH")
    parser.add_argument("--pack-local-storage", action="store_true")
    args = parser.parse_args()
    destination = confined(ROOT, args.manifest)
    if args.pack_local_storage:
        record = json.loads(destination.read_text(encoding="utf-8"))
        receipt = pack_local_storage(record)
        destination.write_text(json.dumps(record, sort_keys=True, indent=2) + "\n", encoding="utf-8", newline="\n")
        print(json.dumps(receipt, sort_keys=True))
    elif args.record_successors:
        record = json.loads(destination.read_text(encoding="utf-8"))
        record_successors(record, args.record_successors)
        destination.write_text(json.dumps(record, sort_keys=True, indent=2) + "\n", encoding="utf-8", newline="\n")
        print(json.dumps({"explicit_successors": args.record_successors}))
    elif args.apply or args.check:
        record = json.loads(destination.read_text(encoding="utf-8"))
        print(json.dumps(check(record) if args.check else apply(record), sort_keys=True))
    else:
        record = plan()
        destination.parent.mkdir(parents=True, exist_ok=True)
        with destination.open("x", encoding="utf-8", newline="\n") as stream:
            json.dump(record, stream, sort_keys=True, indent=2)
            stream.write("\n")
        print(json.dumps({"manifest": str(args.manifest), "files": len(record["rows"]),
                          "counts": record["counts"], "dispositions": record["dispositions"]}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
