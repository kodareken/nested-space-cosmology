"""Resolve an original scientific input without changing its recorded identity."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import subprocess

ARCHIVE = "archive/pre-closure-2026-09-08"


def resolve_pinned_source_bytes(
        root: Path, original_path: str, expected_sha256: str, *, commit: str) -> bytes:
    """Authenticate a source at an explicitly named immutable Git commit.

    This is a historical replay operation, separate from ``resolve_source``'s
    current-file drift checks. It never searches history or reads a result
    artifact in place of the pinned source. Retained content carriers are
    available through process-local alternates without changing the cache.
    """
    relative = PurePosixPath(original_path)
    if (relative.is_absolute() or ".." in relative.parts
            or relative.as_posix() != original_path):
        raise ValueError("source path must stay within the repository")
    source_parts = relative.parts[1:] if relative.parts[:1] == ("lab",) else relative.parts
    if not source_parts or source_parts[0] not in {"src", "scripts", "tests", "docs"}:
        raise ValueError("pinned replay requires an implementation source")
    if not isinstance(commit, str) or not re.fullmatch(r"[0-9a-f]{40}", commit):
        raise ValueError("full immutable commit required")
    if not isinstance(expected_sha256, str) or not re.fullmatch(r"[0-9a-f]{64}", expected_sha256):
        raise ValueError("source SHA-256 required")
    root = root.resolve()
    env = os.environ.copy()
    cache = root / "lab/.source-history/objects"
    if cache.is_dir():
        env["GIT_ALTERNATE_OBJECT_DIRECTORIES"] = str(cache)
    try:
        raw = subprocess.check_output(
            ["git", "-C", str(root), "show", f"{commit}:{original_path}"],
            env=env, stderr=subprocess.PIPE, timeout=30)
    except (OSError, subprocess.SubprocessError) as exc:
        raise RuntimeError(f"pinned source unavailable: {original_path}") from exc
    if hashlib.sha256(raw).hexdigest() != expected_sha256:
        raise RuntimeError(f"pinned source drift: {original_path}")
    return raw


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
