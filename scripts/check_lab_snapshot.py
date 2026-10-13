#!/usr/bin/env python3
"""Authenticate retained current files and unchanged historical source objects.

This verifies file integrity, not a closed scientific gate. Historical import
provenance is metadata; deliberate later edits require a current successor receipt.
"""
import hashlib
import json
from pathlib import Path
import re
import subprocess
import zlib

from lab import ROOT, environment

SECRET_PATTERNS = (
    re.compile(rb"gh[pousr]_[A-Za-z0-9_]{20,}"),
    re.compile(rb"sk-(?:proj-)?[A-Za-z0-9_-]{20,}"),
    re.compile(rb"AKIA[0-9A-Z]{16}"),
    re.compile(rb"BEGIN (?:RSA|OPENSSH|EC|PGP) PRIVATE KEY"),
)


def safe_path(root, relative):
    path = Path(relative)
    if path.is_absolute() or ".." in path.parts or path.as_posix() != relative:
        raise ValueError("invalid snapshot path")
    full = root / path
    if full.is_symlink() or not full.resolve().is_relative_to(root.resolve()):
        raise ValueError("escaping snapshot path")
    return full


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def check_objects(root=ROOT):
    cache = root / "archive/source-history"
    manifest = json.loads((cache / "manifest.json").read_text())
    if manifest["schema"] != "NSC-PINNED-SOURCE-OBJECTS-v1":
        raise ValueError("unknown source object schema")
    seen = set()
    for row in manifest["objects"]:
        oid = row["oid"]
        if not re.fullmatch(r"[0-9a-f]{40}", oid) or oid in seen:
            raise ValueError("invalid source object identity")
        seen.add(oid)
        framed = zlib.decompress((cache / "objects" / oid[:2] / oid[2:]).read_bytes())
        header, raw = framed.split(b"\0", 1)
        if (hashlib.sha1(framed).hexdigest() != oid
                or hashlib.sha256(framed).hexdigest() != row["sha256"]
                or header != f"{row['type']} {row['bytes']}".encode()
                or len(raw) != row["bytes"]):
            raise ValueError("corrupt pinned source object: " + oid)
        if any(pattern.search(raw) for pattern in SECRET_PATTERNS):
            raise ValueError("possible credential in pinned source object: " + oid)
    actual = {p.parent.name + p.name for p in (cache / "objects").glob("*/*") if p.is_file()}
    if actual != seen:
        raise ValueError("source object inventory mismatch")
    for pin in manifest["pins"]:
        if not re.fullmatch(r"[0-9a-f]{40}", pin["commit"]):
            raise ValueError("invalid pinned commit")
        safe_path(root, pin["path"])
        raw = subprocess.check_output(
            ["git", "-C", str(root), "show", pin["commit"] + ":" + pin["path"]],
            env=environment(root), stderr=subprocess.PIPE,
        )
        if hashlib.sha256(raw).hexdigest() != pin["sha256"]:
            raise ValueError("pinned source replay mismatch")
    return len(seen)


def check_retained_files(manifest, root=ROOT):
    """Authenticate only explicitly retained current-root files."""
    if manifest.get("schema") != "NSC-CONSOLIDATED-LAB-v2":
        raise ValueError("unknown lab snapshot schema")
    if manifest.get("canonical_root") != ".":
        raise ValueError("snapshot must use the current repository root")
    seen = set()
    total_bytes = 0
    for row in manifest["files"]:
        name = row["path"]
        if name in seen:
            raise ValueError("duplicated retained snapshot path: " + name)
        seen.add(name)
        if name.startswith("archive/") and not name.startswith("archive/source-history/"):
            raise ValueError("retired archive is not a current snapshot input: " + name)
        full = safe_path(root, name)
        expected_bytes = row["bytes"]
        if type(expected_bytes) is not int or expected_bytes < 0:
            raise ValueError("invalid retained byte count: " + name)
        expected_hash = row["sha256"]
        if not isinstance(expected_hash, str) or not re.fullmatch(r"[0-9a-f]{64}", expected_hash):
            raise ValueError("invalid retained SHA256 identity: " + name)
        if not full.is_file():
            raise FileNotFoundError("retained snapshot file missing: " + name)
        if full.stat().st_size != expected_bytes or digest(full) != expected_hash:
            raise ValueError("retained snapshot bytes changed: " + name)
        total_bytes += expected_bytes
    return {"files": len(seen), "bytes": total_bytes}


def check_snapshot(root=ROOT):
    """Check current retained bytes and the unchanged historical Merkle pins."""
    manifest = json.loads((root / "docs/lab-snapshot.json").read_text(encoding="utf-8"))
    checked = check_retained_files(manifest, root)
    count = check_objects(root)
    print(f"root preservation passed: {checked['files']} retained files; "
          f"{checked['bytes']} bytes; {count} pinned source objects; scientific gate OPEN")
    return {**checked, "pinned_objects": count}

if __name__ == "__main__":
    check_snapshot()
