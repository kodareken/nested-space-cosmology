"""Authenticate current files or explicitly selected original source bytes.

Historical reads use local Git objects only. They never import old code,
materialize a checkout, or substitute historical numerical outputs. Callers
record current successor bindings separately from a saved source identity.
"""

from __future__ import annotations

from contextlib import contextmanager
from contextvars import ContextVar
import hashlib
import os
from pathlib import Path, PurePath, PurePosixPath, PureWindowsPath
import re
import subprocess
from typing import Iterator


_SELECTED_COMMIT: ContextVar[str | None] = ContextVar(
    "nsc_explicit_frozen_source_commit", default=None)
_SOURCE_SUFFIXES = {
    ".py", ".pyi", ".cpp", ".c", ".h", ".hpp", ".md", ".rst",
    ".tex", ".bib", ".toml", ".txt", ".yaml", ".yml",
    ".cfg", ".ini", ".cff",
}
_ROOT_SOURCES = {
    "pyproject.toml", "requirements-history.txt", "README.md", "THEORY.md",
    "AGENTS.md", "CITATION.cff",
}


def _relative_path(relative_path: str | PurePath) -> str:
    name = relative_path.as_posix() if isinstance(relative_path, PurePath) else relative_path
    if not isinstance(name, str):
        raise ValueError("source path must be repository-relative")
    relative = PurePosixPath(name)
    if (not name or "\\" in name or ":" in name or relative.is_absolute()
            or PureWindowsPath(name).drive or ".." in relative.parts
            or relative.as_posix() != name or name == "."):
        raise ValueError("source path must be safe and repository-relative")
    return name


def _source_path(name: str) -> bool:
    parts = PurePosixPath(name).parts
    logical = parts[1:] if parts[:1] == ("lab",) else parts
    # Receipts can live outside the top-level results directory. JSON is not
    # an implementation source in any retained contract; publication metadata
    # has its own explicit authenticator. Evidence subtrees stay local even
    # when they contain a text file with an otherwise source-like extension.
    if any(part in {"results", "artifacts", "evidence"} or part.endswith("-evidence")
           for part in logical):
        return False
    return bool(logical) and (
        "/".join(logical) in _ROOT_SOURCES
        or (logical[0] in {"src", "scripts", "tests", "docs", "paper"}
            and PurePosixPath(name).suffix in _SOURCE_SUFFIXES))


def _commit(context: str) -> str:
    if not isinstance(context, str) or not re.fullmatch(r"[0-9a-f]{40}", context):
        raise ValueError("full immutable source commit required")
    return context


def is_implementation_source(relative_path: str | PurePath) -> bool:
    """Classify a safe relative source path; receipts and evidence are excluded.

    This predicate performs no lookup and selects no historical context.
    Callers with legacy Git adapters must apply it before historical reads.
    """
    return _source_path(_relative_path(relative_path))


def selected_source_context() -> str | None:
    """Return the explicitly selected scope for source-authentication cache keys.

    An authenticated cache must partition by this value so an original-source
    replay cannot satisfy a later current-root lookup. Direct authentication
    with context=None continues to ignore this selection and use local files.
    """
    return _SELECTED_COMMIT.get()


def _current_bytes(root: Path, name: str) -> bytes:
    candidate = root / name
    if candidate.is_symlink() or not candidate.resolve().is_relative_to(root):
        raise ValueError("current source escapes repository: " + name)
    if not candidate.is_file():
        raise FileNotFoundError("current source unavailable: " + name)
    return candidate.read_bytes()


def _source_bytes(root: str | Path, relative_path: str | PurePath,
                  context: str | None) -> bytes:
    root = Path(root).resolve()
    name = _relative_path(relative_path)
    if context is None:
        return _current_bytes(root, name)
    context = _commit(context)
    if not _source_path(name):
        raise ValueError("frozen authentication requires an implementation source")
    environment = os.environ.copy()
    # The explicit repository owns this lookup, independently of a caller's
    # temporary Git work-tree configuration.
    for key in ("GIT_DIR", "GIT_WORK_TREE", "GIT_INDEX_FILE", "GIT_COMMON_DIR"):
        environment.pop(key, None)
    objects = root / "archive/source-history/objects"
    if objects.is_dir():
        inherited = environment.get("GIT_ALTERNATE_OBJECT_DIRECTORIES", "")
        environment["GIT_ALTERNATE_OBJECT_DIRECTORIES"] = (
            str(objects) + (os.pathsep + inherited if inherited else ""))
    try:
        return subprocess.check_output(
            ["git", "--no-pager", "-C", str(root), "show", "--no-textconv",
             f"{context}:{name}"], env=environment, stderr=subprocess.PIPE,
            timeout=30,
        )
    except (OSError, subprocess.SubprocessError) as error:
        raise RuntimeError(f"frozen source unavailable at {context}: {name}") from error


def source_sha256(root: str | Path, relative_path: str | PurePath, *,
                  context: str | None = None) -> str:
    """Hash current bytes by default; only a supplied commit selects history."""
    return hashlib.sha256(_source_bytes(root, relative_path, context)).hexdigest()


def authenticate_source(root: str | Path, relative_path: str | PurePath,
                        expected_sha256: str, *, context: str | None = None) -> bytes:
    """Return exact authenticated bytes, preserving the original path identity."""
    if (not isinstance(expected_sha256, str)
            or not re.fullmatch(r"[0-9a-f]{64}", expected_sha256)):
        raise ValueError("source SHA-256 required")
    raw = _source_bytes(root, relative_path, context)
    if hashlib.sha256(raw).hexdigest() != expected_sha256:
        origin = "current" if context is None else "frozen"
        raise RuntimeError(f"{origin} source drift: {relative_path}")
    return raw


@contextmanager
def frozen_source_context(full_commit: str) -> Iterator[str]:
    """Explicitly select original sources for legacy-signature hash helpers."""
    commit = _commit(full_commit)
    token = _SELECTED_COMMIT.set(commit)
    try:
        yield commit
    finally:
        _SELECTED_COMMIT.reset(token)


def selected_source_sha256(root: str | Path, relative_path: str | PurePath) -> str:
    """Use the selected source pin; result and binary output paths stay local.

    This opt-in adapter is distinct from source_sha256(context=None), which
    always reads current files even inside a frozen_source_context.
    """
    name = _relative_path(relative_path)
    context = _SELECTED_COMMIT.get() if _source_path(name) else None
    return source_sha256(root, name, context=context)
