#!/usr/bin/env python3
"""Write once or compactly verify FGC-1-TDG8-RCV3-FRZ1.

``--write`` is the only mode that observes the ignored source campaign and the
one-time destination-absence premise.  ``--verify-compact`` is store-blind, so
it remains valid after a separately authorized bootstrap creates the frozen
projection.
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

from recursive_horizons.fgc.evolution import tdg8_rcv3_freeze as freeze  # noqa: E402


def _config() -> bytes:
    return freeze.read_nofollow(ROOT, freeze.CONFIG_PATH, "RCV3 FRZ1 config")


def _result() -> bytes:
    return freeze.read_nofollow(ROOT, freeze.RESULT_PATH, "RCV3 FRZ1 result")


def _publish_new(raw: bytes) -> None:
    path = ROOT / freeze.RESULT_PATH
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
                raise OSError("short RCV3 FRZ1 result write")
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


def write() -> bytes:
    config = _config()
    evidence = freeze.derive_live_evidence(ROOT)
    raw = freeze.canonical_result(freeze.build_result(config, evidence))
    _publish_new(raw)
    observed = _result()
    freeze.validate_compact(config, observed)
    if observed != raw:
        raise ValueError("published RCV3 FRZ1 result bytes differ")
    return raw


def verify_compact() -> bytes:
    config, result = _config(), _result()
    freeze.validate_compact(config, result)
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write", action="store_true")
    mode.add_argument("--verify-compact", action="store_true")
    args = parser.parse_args(argv)
    if args.write:
        raw = write()
        selected = "write"
        live_observation = True
    else:
        raw = verify_compact()
        selected = "verify_compact"
        live_observation = False
    print(json.dumps({
        "artifact_id": freeze.ARTIFACT_ID,
        "mode": selected,
        "result_sha256": sha256(raw).hexdigest(),
        "live_observation": live_observation,
        "source_store_mutated": False,
        "destination_store_created": False,
        "bootstrap_authorized": True,
        "bounded_execution_authorized": False,
        "candidate_branch_opened": False,
        "physical_result_earned": False,
    }, sort_keys=True, separators=(",", ":")))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
