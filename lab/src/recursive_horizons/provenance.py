"""Resolve an original scientific input without changing its recorded identity."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path, PurePosixPath

ARCHIVE = "archive/pre-closure-2026-09-08"


def resolve_source(root: Path, original_path: str, expected_sha256: str) -> Path:
    """Find the exact recorded bytes, either active or in the frozen archive.

    A differing active file is an error, not permission to hide drift by
    falling back to history. The caller keeps the ORIGINAL path in its result.
    """
    relative = PurePosixPath(original_path)
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError("source path must stay within the repository")
    root = root.resolve()
    candidate = root / relative
    if not candidate.exists():
        manifest = json.loads((root / ARCHIVE / "manifest.json").read_text())
        entries = {e["original_path"]: e for e in manifest["entries"]}
        item = entries.get(original_path)
        if item is None or item["sha256"] != expected_sha256:
            raise RuntimeError(f"no authenticated archived source: {original_path}")
        candidate = root / item["storage_path"]
    if not candidate.resolve().is_relative_to(root) or candidate.is_symlink():
        raise RuntimeError(f"source escapes repository: {original_path}")
    if hashlib.sha256(candidate.read_bytes()).hexdigest() != expected_sha256:
        raise RuntimeError(f"source drift: {original_path}")
    return candidate
