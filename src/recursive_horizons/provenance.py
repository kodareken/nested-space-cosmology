"""Authenticate current root sources or explicitly pinned original Git bytes."""

from __future__ import annotations

from pathlib import Path

from .nsc_frozen_sources import authenticate_source


def resolve_pinned_source_bytes(
        root: Path, original_path: str, expected_sha256: str, *, commit: str) -> bytes:
    """Read one explicit original source context without executing or copying it."""
    try:
        return authenticate_source(root, original_path, expected_sha256, context=commit)
    except RuntimeError as error:
        # Preserve the public provenance exception wording for its callers.
        raise RuntimeError(str(error).replace("frozen source", "pinned source", 1)) from error


def resolve_source(root: Path, original_path: str, expected_sha256: str) -> Path:
    """Require exact current root bytes; missing or changed inputs never use history.

    Retired migration receipts and archived source copies are not current-source
    substitutes. Original source authentication uses resolve_pinned_source_bytes
    with a full explicit commit and returns bytes rather than a current path.
    """
    root = Path(root).resolve()
    try:
        authenticate_source(root, original_path, expected_sha256, context=None)
    except FileNotFoundError as error:
        raise RuntimeError("current source unavailable: " + original_path) from error
    return root / original_path
