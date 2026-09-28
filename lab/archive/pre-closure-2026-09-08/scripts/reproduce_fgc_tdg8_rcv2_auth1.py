#!/usr/bin/env python3
"""Write or verify the compact FGC-1-TDG8-RCV2-AUTH1 authority.

``--write`` performs RCV2's sole read-only observation of the exact live
entry. ``--verify-compact`` reads tracked bytes only and never opens a run
store.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.fgc.evolution import tdg8_bounded_retry_authority as authority  # noqa: E402


def _config() -> bytes:
    return authority.rcv1._read_leaf(
        ROOT, authority.CONFIG_PATH, "RCV2 config", 1024 * 1024,
    )


def _result() -> bytes:
    return authority.rcv1._read_leaf(
        ROOT, authority.RESULT_PATH, "RCV2 result", 4 * 1024 * 1024,
    )


def _publish_new(raw: bytes) -> None:
    path = ROOT / authority.RESULT_PATH
    descriptor = os.open(
        path.parent,
        os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0),
    )
    try:
        flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0)
        file_descriptor = os.open(path.name, flags, 0o644, dir_fd=descriptor)
        try:
            view = memoryview(raw)
            while view:
                written = os.write(file_descriptor, view)
                if written <= 0:
                    raise OSError("short RCV2 result write")
                view = view[written:]
            os.fsync(file_descriptor)
        except Exception:
            os.close(file_descriptor)
            try:
                os.unlink(path.name, dir_fd=descriptor)
            except OSError:
                pass
            raise
        else:
            os.close(file_descriptor)
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def write() -> bytes:
    config = _config()
    raw = authority.canonical(authority.build_prelaunch(config, ROOT))
    _publish_new(raw)
    observed = _result()
    authority.validate_compact(config, observed)
    if observed != raw:
        raise ValueError("published RCV2 result bytes differ")
    return raw


def verify_compact() -> bytes:
    config, result = _config(), _result()
    authority.validate_compact(config, result)
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--verify-compact", action="store_true")
    args = parser.parse_args(argv)
    if args.write:
        result = write()
        selected = "write"
    else:
        result = verify_compact()
        selected = "verify_compact"
    print(json.dumps({
        "artifact_id": authority.ARTIFACT_ID,
        "mode": selected,
        "result_sha256": sha256(result).hexdigest(),
        "store_mutated": False,
        "candidate_branch_opened": False,
        "calibration_result_earned": False,
        "physical_result_earned": False,
    }, sort_keys=True, separators=(",", ":")))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
