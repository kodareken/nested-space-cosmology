#!/usr/bin/env python3
"""Record or check the imbalance mean-circulation witness. The physical gate stays open."""
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

OUTPUT = ROOT / "results/development/nsc-imbalance-turnover-v1.json"
SCHEMA = "NSC-IMBALANCE-TURNOVER-v1"
PASS = "PASS_IMBALANCE_MEAN_CIRCULATION"
# Actual imbalance code, this recorder, the original finite window, and the helpers it calls.
# Not this record, and not a Git HEAD.
SOURCES = (
    "src/recursive_horizons/nsc_imbalance_turnover.py",
    "scripts/derive_nsc_imbalance_turnover.py",
    "src/recursive_horizons/nsc_nested_qualities.py",
    "src/recursive_horizons/nsc_finite_turnover.py",
)
METRICS = (
    "Gbar_region0_cut_exact",
    "Gbar_region0_cut_lower_exact",
    "Gbar_region0_cut_positive",
    "net_region0_cut_exact",
    "Gbar_region1_to_region2_exact",
    "Gbar_region1_to_region2_lower_exact",
    "persistent_mean_mode_circulation",
    "nonzero_mean_channels",
    "mean_adjacent_currents_exact",
    "initial_adjacent_current_derivatives_exact",
    "projection_coefficients_exact",
    "gram_determinant_exact",
    "total_content_exact",
    "total_energy_exact",
    "initial_region_populations_exact",
    "mean_region_populations_exact",
    "initial_region_population_rates_exact",
    "mean_region_population_rates_exact",
    "uniform_adjacent_currents_max_abs_exact",
    "delta_scale_factor_exact",
    "spectral_dephasing_frobenius_error",
    "spectral_dephasing_digits",
    "projection_entry_bitlength",
    "numeric_eigenvalue_gap_control",
    "filtered_mean_response_frobenius",
    "full_versus_reduced_frobenius",
    "self_energy_versus_covariance_frobenius",
    "perturbed_filtered_response_frobenius",
    "perturbed_full_versus_reduced_frobenius",
    "perturbed_versus_baseline_response_frobenius",
    "frozen_cbar_shortcut_versus_recomputed_frobenius",
    "frozen_cbar_shortcut_versus_baseline_frobenius",
    "mode_filling_to_mean_current_map",
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
        "verdict_scope",
        "question",
        "answer",
        "physical_local_gate",
        "inputs",
        "checks",
        "metrics",
        "finite_time_bound",
        "finite_time_bound_error",
        "filling_control_error",
        "preparation_conditions",
        "assumptions",
        "remaining_connection",
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
        raise ValueError("PASS_IMBALANCE_MEAN_CIRCULATION requires every check to pass")
    metrics = result["metrics"]
    if not isinstance(metrics, dict):
        raise ValueError("metrics must be an object")
    missing_metrics = [key for key in METRICS if key not in metrics]
    if missing_metrics:
        raise ValueError(f"metrics missing {missing_metrics}")


def payload() -> dict:
    from recursive_horizons.nsc_imbalance_turnover import compute_result

    result = compute_result()
    require_api(result)
    assembled = dict(result)
    assembled["source_bindings"] = source_bindings()
    return assembled


def canonical(value: dict) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n").encode("utf-8")


def reject_mismatch(recorded: bytes, fresh: bytes) -> None:
    previous = json.loads(recorded)
    current = json.loads(fresh)
    recorded_sources = previous.get("source_bindings")
    fresh_sources = current.get("source_bindings")
    if not isinstance(recorded_sources, dict) or not isinstance(fresh_sources, dict):
        raise ValueError("source hash mismatch")
    for relative in SOURCES:
        if relative not in recorded_sources or relative not in fresh_sources:
            raise FileNotFoundError(f"missing bound source: {relative}")
        if recorded_sources[relative] != fresh_sources[relative]:
            raise ValueError(f"source hash mismatch: {relative}")
    if recorded != fresh:
        raise ValueError("imbalance-turnover record differs from recomputation")


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
    else:
        reject_mismatch(OUTPUT.read_bytes(), body)
    summary = json.loads(body)
    try:
        output = OUTPUT.resolve().relative_to(ROOT.resolve()).as_posix()
    except ValueError:
        output = OUTPUT.as_posix()
    print(json.dumps({
        "schema": summary["schema"],
        "verdict": summary["verdict"],
        "physical_local_gate": summary["physical_local_gate"],
        "checks": len(summary["checks"]),
        "output": output,
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
