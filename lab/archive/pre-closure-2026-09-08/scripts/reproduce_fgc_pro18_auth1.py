#!/usr/bin/env python3
"""Reproduce the content-only FGC-1-PRO18-AUTH1 authority artifact.

This command deliberately has no launch mode.  It reads only tracked compact
inputs plus the prospective source files, validates the pinned evolution
environment, and writes/verifies the one canonical AUTH1 JSON result.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import platform
import sys
import tomllib

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src")]

from recursive_horizons.fgc.evolution.proto17_pure_construction import digest  # noqa: E402
from recursive_horizons.fgc.evolution.proto18_authority import (  # noqa: E402
    Proto18AuthorityError,
    build_authority,
    verify_authority_payload,
)

CONFIG = ROOT / "configs/fgc/fgc-1-pro18-auth1.toml"
RESULT = ROOT / "results/fgc-1-pro18-auth1.json"


def canonical(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, indent=2, ensure_ascii=True,
                       allow_nan=False) + "\n").encode("utf-8")


def read_canonical(path: Path) -> dict[str, object]:
    def no_duplicates(pairs: list[tuple[str, object]]) -> dict[str, object]:
        result: dict[str, object] = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(key)
            result[key] = value
        return result
    raw = path.read_bytes()
    try:
        value = json.loads(raw.decode("utf-8"), object_pairs_hook=no_duplicates)
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as error:
        raise Proto18AuthorityError("AUTH1_RESULT_INVALID", "stored result is malformed") from error
    if not isinstance(value, dict) or canonical(value) != raw:
        raise Proto18AuthorityError("AUTH1_RESULT_INVALID", "stored result is noncanonical")
    return value


def observed_environment() -> dict[str, object]:
    return {
        "contract_id": "FGC-1-PRO18-AUTH1-evolution-environment-v1",
        "implementation": "CPython" if sys.implementation.name == "cpython" else sys.implementation.name,
        "python_version": platform.python_version(),
        "python_build": " | ".join(platform.python_build()),
        "python_compiler": platform.python_compiler(),
        "numpy_version": np.__version__,
        "numpy_build_configuration": json.dumps(
            np.show_config(mode="dicts"), sort_keys=True, separators=(",", ":"), ensure_ascii=True,
        ),
        "blas_lapack": "Accelerate" if "accelerate" in json.dumps(np.show_config(mode="dicts")).lower() else "other",
        "operating_system": platform.platform(),
        "system": platform.system(),
        "machine": platform.machine(),
        "processor": platform.processor(),
        "platform": sys.platform,
        "byteorder": sys.byteorder,
    }


def require_pinned_environment(config: dict[str, object]) -> None:
    declared = config.get("numerical_environment")
    if not isinstance(declared, dict) or not isinstance(declared.get("contract"), dict):
        raise Proto18AuthorityError("AUTH1_ENVIRONMENT_DRIFT", "environment contract is absent")
    actual = observed_environment()
    if actual != declared["contract"]:
        raise Proto18AuthorityError("AUTH1_ENVIRONMENT_DRIFT", "current interpreter/environment differs from AUTH1")
    if digest(actual) != declared.get("canonical_sha256"):
        raise Proto18AuthorityError("AUTH1_ENVIRONMENT_DRIFT", "current environment digest differs from AUTH1")


def atomic_write(path: Path, payload: bytes) -> None:
    temporary = path.with_name(f".{path.name}.tmp")
    if temporary.exists():
        raise Proto18AuthorityError("AUTH1_RESULT_INVALID", "stale temporary result exists")
    try:
        with temporary.open("xb") as handle:
            handle.write(payload)
            handle.flush()
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--write", action="store_true", help="write the fixed AUTH1 result path")
    group.add_argument("--verify", action="store_true", help="verify the fixed AUTH1 result path")
    args = parser.parse_args()
    config = tomllib.loads(CONFIG.read_text("utf-8"))
    require_pinned_environment(config)
    expected = build_authority(ROOT, config)
    if args.write:
        atomic_write(RESULT, canonical(expected))
    else:
        verify_authority_payload(ROOT, config, read_canonical(RESULT))
    print(json.dumps({
        "artifact_id": "FGC-1-PRO18-AUTH1",
        "result_sha256": sha256(RESULT.read_bytes()).hexdigest() if RESULT.exists() else None,
        "verified": bool(args.verify),
        "written": bool(args.write),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
