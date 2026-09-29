#!/usr/bin/env python3
"""Record or check the finite channel-turnover witness. The physical gate stays open."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

OUTPUT = ROOT / "results/development/nsc-finite-turnover-v1.json"
SCHEMA = "NSC-FINITE-CHANNEL-TURNOVER-v1"
PASS = "PASS_FINITE_CHANNEL_TURNOVER"
# These three sources only. Not this record, and not a Git HEAD.
SOURCES = (
    "src/recursive_horizons/nsc_finite_turnover.py",
    "scripts/derive_nsc_finite_turnover.py",
    "src/recursive_horizons/nsc_nested_qualities.py",
)
METRICS = (
    "G_cut_exact",
    "total_inventory_exact",
    "sign_control_value_exact",
    "probe_response_change_fro",
    "omitted_memory_error_fro",
    "omitted_drive_error_fro",
    "dropped_cross_covariance_error_fro",
    "half_identity_response_change_fro",
    "far_probe_response_change_fro",
)


def source_bindings() -> dict[str, str]:
    bindings = {}
    for relative in SOURCES:
        path = ROOT / relative
        if not path.is_file():
            raise FileNotFoundError(path)
        bindings[relative] = hashlib.sha256(path.read_bytes()).hexdigest()
    return bindings


def require_api(result: dict) -> None:
    if not isinstance(result, dict):
        raise TypeError("compute_result() must return a dict")
    required = (
        "schema",
        "verdict",
        "physical_local_gate",
        "inputs",
        "checks",
        "metrics",
        "assumptions",
    )
    missing = [key for key in required if key not in result]
    if missing:
        raise ValueError(f"compute_result() missing {missing}")
    if result["schema"] != SCHEMA:
        raise ValueError(f"schema must be {SCHEMA}")
    if result["physical_local_gate"] != "OPEN":
        raise ValueError("physical_local_gate must stay OPEN")
    checks = result["checks"]
    if not isinstance(checks, dict) or not checks:
        raise ValueError("checks must be a non-empty object of booleans")
    if any(type(value) is not bool for value in checks.values()):
        raise ValueError("checks must be booleans only")
    if result["verdict"] == PASS and not all(checks.values()):
        raise ValueError("PASS_FINITE_CHANNEL_TURNOVER requires every check to pass")
    metrics = result["metrics"]
    if not isinstance(metrics, dict):
        raise ValueError("metrics must be an object")
    missing_metrics = [key for key in METRICS if key not in metrics]
    if missing_metrics:
        raise ValueError(f"metrics missing {missing_metrics}")


def payload() -> dict:
    from recursive_horizons.nsc_finite_turnover import compute_result

    result = compute_result()
    require_api(result)
    assembled = dict(result)
    assembled["source_bindings"] = source_bindings()
    return assembled


def canonical(value: dict) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n").encode("utf-8")


def write_new(path: Path, body: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    temporary = Path(name)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(body)
            handle.flush()
            os.fsync(handle.fileno())
        try:
            os.link(temporary, path)
        except FileExistsError as error:
            raise FileExistsError(f"{path} exists; use --check") from error
    finally:
        temporary.unlink(missing_ok=True)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--record", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = parser.parse_args(argv)
    if args.check and not OUTPUT.is_file():
        raise FileNotFoundError(OUTPUT)
    body = canonical(payload())
    if args.record:
        write_new(OUTPUT, body)
    elif OUTPUT.read_bytes() != body:
        raise ValueError("finite channel-turnover record differs from recomputation")
    summary = json.loads(body)
    print(json.dumps({
        "schema": summary["schema"],
        "verdict": summary["verdict"],
        "physical_local_gate": summary["physical_local_gate"],
        "checks": len(summary["checks"]),
        "output": OUTPUT.relative_to(ROOT).as_posix(),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
