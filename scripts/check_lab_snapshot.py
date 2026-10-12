#!/usr/bin/env python3
"""Authenticate the imported lab and retained historical source objects.

This verifies file integrity, not a closed scientific gate. The manifest is
an import receipt; deliberate later edits require an explicit successor receipt.
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


def check_snapshot(root=ROOT):
    """Verify the relocation and the unchanged original import identities."""
    from reassemble_repository import check as check_reassembly

    original = json.loads((root / "docs/lab-snapshot.json").read_text(encoding="utf-8"))
    migration = json.loads((root / "results/reassembly/manifest.json").read_text(encoding="utf-8"))
    if original["schema"] != "NSC-CONSOLIDATED-LAB-v1":
        raise ValueError("unknown lab snapshot schema")
    if migration["schema"] != "NSC-REASSEMBLY-v1":
        raise ValueError("unknown reassembly schema")
    originals = {row["original_path"]: row for row in migration["rows"] if row["origin"] == "lab"}
    seen = set()
    for row in original["files"]:
        name = row["path"]
        if name in seen or not name.startswith("lab/"):
            raise ValueError("invalid or duplicated lab path")
        seen.add(name)
        migrated = originals.get(name)
        if migrated is None or migrated["sha256"] != row["sha256"] or migrated["bytes"] != row["bytes"]:
            raise ValueError("original snapshot identity changed: " + name)
    for omitted in original["intentionally_omitted"]:
        if "lab/" + omitted in originals:
            raise ValueError("historical payload unexpectedly imported: " + omitted)
    checked = check_reassembly(migration, root)
    count = check_objects(root)
    print(f"root preservation passed: {checked['total']} destinations; "
          f"{len(seen)} original import identities; {count} pinned source objects; scientific gate OPEN")

if __name__ == "__main__":
    check_snapshot()
