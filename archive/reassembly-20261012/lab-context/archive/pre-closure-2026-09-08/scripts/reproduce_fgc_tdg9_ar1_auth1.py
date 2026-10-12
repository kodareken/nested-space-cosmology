#!/usr/bin/env python3
"""Write once or compactly verify FGC-1-TDG9-AR1-AUTH1.

``--write`` performs AR1's sole live, read-only PREF2 authentication,
environment check, exact prospective Git-delta check, and one-time
diagnostic-namespace absence observation. ``--verify-compact`` reads only
tracked configuration/result bytes, so it remains valid after the separately
authorized diagnostic publishes its raw result. The eventual authority commit
is necessarily supplied at runtime rather than self-bound here.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import secrets
import stat
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.fgc.evolution import tdg9_ar1_authority as authority  # noqa: E402


REJECTED_RESULT_SHA256S = (
    "c98d9b53eb08f388fdd620bc04b73311da02a497ad1594886bd44c7acc73a6f7",
    "1c889bd63f78313a820fcb8d4d298fcaa7d56e06ced4ce0c99c12b859da43503",
)
_REWRITE_STAGE_PREFIX = ".fgc-1-tdg9-ar1-auth1.rewrite-"
_MAX_RESULT_BYTES = 8 * 1024 * 1024


def _config() -> bytes:
    return authority._read_regular(  # type: ignore[attr-defined]
        ROOT, authority.CONFIG_PATH, 2 * 1024 * 1024
    )


def _result() -> bytes:
    return authority._read_regular(  # type: ignore[attr-defined]
        ROOT, authority.RESULT_PATH, 8 * 1024 * 1024
    )


def _publish_new(raw: bytes, root: Path = ROOT) -> None:
    path = root / authority.RESULT_PATH
    directory_fd = os.open(
        path.parent,
        os.O_RDONLY
        | getattr(os, "O_DIRECTORY", 0)
        | getattr(os, "O_NOFOLLOW", 0),
    )
    leaf_fd = -1
    created = False
    try:
        leaf_fd = os.open(
            path.name,
            os.O_WRONLY
            | os.O_CREAT
            | os.O_EXCL
            | getattr(os, "O_NOFOLLOW", 0),
            0o644,
            dir_fd=directory_fd,
        )
        created = True
        view = memoryview(raw)
        while view:
            written = os.write(leaf_fd, view)
            if written <= 0:
                raise OSError("short TDG9 AR1 authority result write")
            view = view[written:]
        os.fsync(leaf_fd)
        os.close(leaf_fd)
        leaf_fd = -1
        os.fsync(directory_fd)
    except Exception:
        if leaf_fd != -1:
            os.close(leaf_fd)
        if created:
            try:
                os.unlink(path.name, dir_fd=directory_fd)
            except OSError:
                pass
        raise
    finally:
        os.close(directory_fd)


def _read_rejected_identity(directory_fd: int, leaf: str) -> tuple[int, int, int, int]:
    try:
        before = os.stat(leaf, dir_fd=directory_fd, follow_symlinks=False)
    except OSError as exc:
        raise ValueError("rejected TDG9 AR1 result is absent") from exc
    if stat.S_ISLNK(before.st_mode) or not stat.S_ISREG(before.st_mode):
        raise ValueError("rejected TDG9 AR1 result is not a regular leaf")
    if before.st_size > _MAX_RESULT_BYTES:
        raise ValueError("rejected TDG9 AR1 result is too large")
    descriptor = os.open(
        leaf,
        os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
        dir_fd=directory_fd,
    )
    try:
        active = os.fstat(descriptor)
        identity = (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns)
        if identity != (active.st_dev, active.st_ino, active.st_size, active.st_mtime_ns):
            raise ValueError("rejected TDG9 AR1 result raced before read")
        digest = sha256()
        size = 0
        while size <= _MAX_RESULT_BYTES:
            block = os.read(descriptor, min(1 << 20, _MAX_RESULT_BYTES + 1 - size))
            if not block:
                break
            digest.update(block)
            size += len(block)
        after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    if identity != (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns):
        raise ValueError("rejected TDG9 AR1 result changed during read")
    if size != before.st_size or digest.hexdigest() not in REJECTED_RESULT_SHA256S:
        raise ValueError("rejected TDG9 AR1 result has unknown bytes")
    return identity


def _replace_rejected(raw: bytes, root: Path = ROOT) -> None:
    """Atomically replace only the exact rejected uncommitted compact result."""

    path = root / authority.RESULT_PATH
    directory_fd = os.open(
        path.parent,
        os.O_RDONLY
        | getattr(os, "O_DIRECTORY", 0)
        | getattr(os, "O_NOFOLLOW", 0),
    )
    stage_name = f"{_REWRITE_STAGE_PREFIX}{os.getpid()}-{secrets.token_hex(8)}"
    stage_created = False
    stage_fd = -1
    try:
        original_identity = _read_rejected_identity(directory_fd, path.name)
        stage_fd = os.open(
            stage_name,
            os.O_WRONLY
            | os.O_CREAT
            | os.O_EXCL
            | getattr(os, "O_NOFOLLOW", 0),
            0o644,
            dir_fd=directory_fd,
        )
        stage_created = True
        view = memoryview(raw)
        while view:
            written = os.write(stage_fd, view)
            if written <= 0:
                raise OSError("short TDG9 AR1 rewrite stage write")
            view = view[written:]
        os.fsync(stage_fd)
        os.close(stage_fd)
        stage_fd = -1
        # This is intentionally the final operation before the atomic replace.
        if _read_rejected_identity(directory_fd, path.name) != original_identity:
            raise ValueError("rejected TDG9 AR1 result identity changed")
        os.replace(
            stage_name,
            path.name,
            src_dir_fd=directory_fd,
            dst_dir_fd=directory_fd,
        )
        stage_created = False
        os.fsync(directory_fd)
    finally:
        if stage_fd != -1:
            os.close(stage_fd)
        if stage_created:
            try:
                os.unlink(stage_name, dir_fd=directory_fd)
            except OSError:
                pass
        os.close(directory_fd)


def write() -> bytes:
    config = _config()
    raw = authority.canonical_pretty(authority.build_prelaunch(config, ROOT))
    _publish_new(raw)
    observed = _result()
    authority.validate_compact(config, observed)
    if observed != raw:
        raise ValueError("published TDG9 AR1 authority result bytes differ")
    return raw


def rewrite_rejected() -> bytes:
    config = _config()
    raw = authority.canonical_pretty(authority.build_prelaunch(config, ROOT))
    _replace_rejected(raw)
    observed = _result()
    authority.validate_compact(config, observed)
    if observed != raw:
        raise ValueError("rewritten TDG9 AR1 authority result bytes differ")
    return raw


def verify_compact() -> bytes:
    config, result = _config(), _result()
    authority.validate_compact(config, result)
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--rewrite-rejected", action="store_true")
    mode.add_argument("--verify-compact", action="store_true")
    args = parser.parse_args(argv)
    if args.write:
        raw = write()
        selected = "write"
        live_observation = True
    elif args.rewrite_rejected:
        raw = rewrite_rejected()
        selected = "rewrite_rejected"
        live_observation = True
    else:
        raw = verify_compact()
        selected = "verify_compact"
        live_observation = False
    print(
        json.dumps(
            {
                "artifact_id": authority.ARTIFACT_ID,
                "mode": selected,
                "result_sha256": sha256(raw).hexdigest(),
                "live_observation": live_observation,
                "environment_bound": True,
                "selection_bound": True,
                "committed_image_bound": True,
                "authority_commit_self_bound": False,
                "campaign_store_mutated": False,
                "diagnostic_executed": False,
                "continuation_authorized": False,
                "candidate_branch_opened": False,
                "physical_result_earned": False,
            },
            sort_keys=True,
            separators=(",", ":"),
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
