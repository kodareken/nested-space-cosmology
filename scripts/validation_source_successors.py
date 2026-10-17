"""Authenticate explicitly reviewed successors of saved implementation sources.

This module never executes historical code. A registered successor authenticates
the original bytes from the pinned Git commit and the exact current bytes. Data,
scientific documents, unregistered changes and changed inventories remain strict.
Creation-time source inventories must continue hashing the current files.
Implementation successors preserve provenance; they do not assert numerical
equivalence or transfer an old result to a changed calculation.
"""
from __future__ import annotations

import hashlib
from functools import lru_cache
import json
from pathlib import Path, PurePosixPath
import re
import subprocess


ORIGIN_COMMIT = "fc7170922ec1e2d9ccec59357f3c68e6a577fb3d"
REGISTRY = "docs/storage/validation-source-successors.json"
SCHEMA = "nsc-validation-source-successors-v1"
RESOLVER_SOURCES = ("scripts/validation_source_successors.py",
                    "tests/test_validation_source_successors.py",
                    "src/recursive_horizons/nsc_atlas_validation.py",
                    "tests/test_nsc_atlas_validation.py",
                    "tests/control_budget.py", "tests/test_control_budget.py")
REVIEWED_PACKAGE_PATHS = frozenset({
    "src/recursive_horizons/nsc_atlas_characteristic_trace.py",
    "src/recursive_horizons/nsc_atlas_refined_cache.py",
    "src/recursive_horizons/nsc_atlas_forcing_cache.py",
    "src/recursive_horizons/nsc_atlas_cap_field.py",
    "src/recursive_horizons/nsc_atlas_background_profile.py",
    "src/recursive_horizons/nsc_atlas_prepared_response.py",
    "src/recursive_horizons/nsc_atlas_characteristic_background.py",
    "src/recursive_horizons/nsc_atlas_cap_mass_derivatives.py",
    "src/recursive_horizons/nsc_atlas_cap_response.py",
    "src/recursive_horizons/nsc_atlas_connection_slice.py",
})
CLASSIFICATIONS = frozenset({"validation-only", "test-only", "implementation-successor"})
_SHA256 = re.compile(r"[0-9a-f]{64}\Z")


def _digest(raw):
    return hashlib.sha256(raw).hexdigest()


def _path(root, name):
    if (not isinstance(name, str) or not name or "\\" in name
            or ":" in name or name != PurePosixPath(name).as_posix()
            or PurePosixPath(name).is_absolute()
            or any(part in {"", ".", "..", ".git"} for part in name.split("/"))):
        raise ValueError("source binding requires a safe repository-relative path")
    path = root / name
    if any(parent.is_symlink() for parent in (path, *path.parents) if parent != root):
        raise ValueError("source binding cannot traverse a symlink")
    if not path.resolve().is_relative_to(root):
        raise ValueError("source binding exits repository")
    return path


def _source_path(name):
    parts = PurePosixPath(name).parts
    return (len(parts) == 2 and parts[0] in {"scripts", "tests"} and name.endswith(".py")) or name in REVIEWED_PACKAGE_PATHS


def _pairs(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError("duplicate validation-successor registry key")
        result[key] = value
    return result


def _registry(root):
    document = json.loads(_path(root, REGISTRY).read_bytes(), object_pairs_hook=_pairs)
    if (not isinstance(document, dict) or document.get("schema") != SCHEMA
            or document.get("origin_commit") != ORIGIN_COMMIT
            or set(document.get("resolver_bindings", {})) != set(RESOLVER_SOURCES)
            or not isinstance(document.get("successors"), dict)):
        raise ValueError("invalid validation-successor registry")
    for name, expected in document["resolver_bindings"].items():
        if not isinstance(expected, str) or not _SHA256.fullmatch(expected) or _digest(_path(root, name).read_bytes()) != expected:
            raise ValueError("validation-successor resolver drift: " + name)
    for name, row in document["successors"].items():
        _path(root, name)
        if (not _source_path(name) or not isinstance(row, dict)
                or set(row) not in ({"original_sha256", "successor_sha256", "classification", "reason"},
                                    {"original_sha256", "successor_sha256", "classification", "reason", "original_commit"})
                or ("original_commit" in row and (not isinstance(row["original_commit"],str)
                    or re.fullmatch(r"[0-9a-f]{40}",row["original_commit"]) is None))
                or row["classification"] not in CLASSIFICATIONS
                or not isinstance(row["reason"], str) or not row["reason"].strip()
                or any(not isinstance(row[k], str) or not _SHA256.fullmatch(row[k]) for k in ("original_sha256", "successor_sha256"))
                or (row["classification"] == "test-only" and not name.startswith("tests/"))
                or (row["classification"] == "implementation-successor" and name not in REVIEWED_PACKAGE_PATHS)):
            raise ValueError("invalid validation-successor entry: " + name)
    return document["successors"]


@lru_cache(maxsize=128)
def _original_bytes(root, commit, name, expected):
    """Cache only successfully authenticated immutable Git bytes.

    Current files and the registry are deliberately checked outside this cache
    on every call. Exceptions are not cached by lru_cache.
    """
    process = subprocess.run(["git", "-C", root, "show", commit + ":" + name],
                             capture_output=True, check=False)
    if process.returncode or _digest(process.stdout) != expected:
        raise ValueError("saved source is not authenticated by the pinned origin: " + name)
    return process.stdout


def verify_saved_source(name, expected, root):
    """Return authenticated ORIGINAL bytes; current bytes when already exact.

    Returning the original bytes preserves saved byte-count receipts. Current
    validation code is still what runs; Git is used only for byte authentication.
    A changed data file is rejected even if someone adds it to the registry.
    """
    root = Path(root).resolve()
    if not isinstance(expected, str) or not _SHA256.fullmatch(expected):
        raise ValueError("invalid saved source digest")
    raw = _path(root, name).read_bytes()
    actual = _digest(raw)
    if actual == expected:
        return raw
    if not _source_path(name):
        raise ValueError("saved data/source binding changed: " + name)
    try:
        row = _registry(root).get(name)
    except FileNotFoundError as error:
        raise ValueError("unregistered saved source binding change: " + name) from error
    if (row is None or row["original_sha256"] != expected
            or row["successor_sha256"] != actual):
        raise ValueError("unregistered saved source binding change: " + name)
    return _original_bytes(str(root), row.get("original_commit",ORIGIN_COMMIT), name, expected)


def verify_saved_bindings(expected, root, *, expected_paths=None):
    """Verify a saved inventory against exact files or named source successors."""
    if not isinstance(expected, dict) or not expected:
        raise ValueError("saved source bindings must be a nonempty mapping")
    if expected_paths is not None and set(expected) != set(expected_paths):
        raise ValueError("saved source inventory changed")
    for name, digest in expected.items():
        verify_saved_source(name, digest, root)
    return True
