#!/usr/bin/env python3
"""Refresh or verify the external-image inventory for PRO19 event one.

This command only reads source/configuration files.  It does not open a run
namespace, construct a solver, or authorize a trajectory.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src")]

from recursive_horizons.fgc.evolution.proto19_launch_manifest import (  # noqa: E402
    Proto19LaunchManifestError, refresh_manifest, render_manifest, verify_manifest,
)


MANIFEST = ROOT / "configs/fgc/fgc-1-pro19-launch-authority.toml"
RELATIVE_MANIFEST = MANIFEST.relative_to(ROOT).as_posix()


def main() -> int:
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--refresh", action="store_true", help="write no-follow SHA-256 values for every existing binding")
    mode.add_argument("--verify", action="store_true", help="reject placeholders and verify every binding")
    args = parser.parse_args()
    raw = MANIFEST.read_bytes()
    try:
        if args.refresh:
            refreshed = render_manifest(refresh_manifest(ROOT, raw, manifest_path=RELATIVE_MANIFEST))
            temporary = MANIFEST.with_name(f".{MANIFEST.name}.tmp")
            if temporary.exists():
                raise Proto19LaunchManifestError("stale launch-manifest temporary file exists")
            temporary.write_bytes(refreshed)
            temporary.replace(MANIFEST)
            print("PRO19 launch manifest refreshed; it remains uncommitted external-image evidence")
        else:
            verify_manifest(ROOT, raw, manifest_path=RELATIVE_MANIFEST)
            print("PRO19 launch manifest verified")
    except Proto19LaunchManifestError as error:
        raise SystemExit(str(error)) from error
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
