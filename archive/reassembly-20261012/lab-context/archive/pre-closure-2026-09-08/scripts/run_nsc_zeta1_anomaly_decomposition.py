#!/usr/bin/env python3
"""Separate cutoff anomaly from physical child-link scale variation."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from math import pi
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
for value in (str(ROOT), str(ROOT / "src")):
    if value not in sys.path:
        sys.path.insert(0, value)

from recursive_horizons.evidence_io import canonical_json_bytes  # noqa: E402


OUTPUT = ROOT / "results/nsc-2-zeta1-anomaly-decomposition.json"
DETERMINANT_RUNNER = ROOT / "scripts/run_nsc_zeta1_regulated_determinant.py"
WARPED_RUNNER = ROOT / "scripts/run_nsc_zeta1_warped_y.py"
INPUTS = (
    ROOT / "results/nsc-2-zeta1-warped-y.json",
    ROOT / "results/nsc-2-zeta1-regulated-determinant.json",
    ROOT / "results/nsc-1-s-one-invariant-partition-bind.json",
    ROOT / "results/nsc-1-s-one-fermionic-geometry-bind.json",
    DETERMINANT_RUNNER,
    WARPED_RUNNER,
)
CHILD_THRESHOLD = 3.0 * pi / 2.0
ANGULAR_MAX = 12


def _load(path: Path, name: str):
    specification = importlib.util.spec_from_file_location(name, path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load owner: {path}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def _authenticate() -> tuple[list[dict[str, object]], dict[str, dict[str, object]]]:
    authenticated: list[dict[str, object]] = []
    records: dict[str, dict[str, object]] = {}
    for path in INPUTS:
        raw = path.read_bytes()
        if path.suffix == ".json":
            item = json.loads(raw)
            if item.get("terminal") is not True:
                raise RuntimeError(f"nonterminal input: {path.relative_to(ROOT)}")
            artifact_id = str(item["artifact_id"])
            records[artifact_id] = item
        else:
            artifact_id = None
        authenticated.append(
            {
                "path": str(path.relative_to(ROOT)),
                "sha256": hashlib.sha256(raw).hexdigest(),
                "artifact_id": artifact_id,
            }
        )
    return authenticated, records


def _mode_parts(zeta: float, values: np.ndarray) -> np.ndarray:
    quotient = (values + 1.0) / zeta - CHILD_THRESHOLD / zeta**2
    gradient = (values + 1.0) / zeta - 2.0 * CHILD_THRESHOLD / zeta**2
    heat = np.exp(-quotient)
    determinant_derivative = 0.5 * np.sum(heat * gradient / quotient)
    pure_cutoff_anomaly = 0.5 * np.sum(heat)
    child_link_derivative = determinant_derivative - pure_cutoff_anomaly
    return np.asarray(
        [determinant_derivative, pure_cutoff_anomaly, child_link_derivative]
    )


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("ZETA1 anomaly-decomposition output already exists")
    authenticated, records = _authenticate()
    determinant_owner = _load(
        DETERMINANT_RUNNER, "_nsc_zeta1_determinant_owner"
    )
    warped_owner = _load(WARPED_RUNNER, "_nsc_zeta1_warped_owner")
    radial = [
        determinant_owner._sector_spectrum(30.0, 750, angular)
        for angular in range(1, ANGULAR_MAX + 1)
    ]

    def build_y(intervals: int) -> list[tuple[float, int]]:
        positive = warped_owner._warped_y_squared_modes(intervals)
        return [(0.0, 1)] + [(float(value), 2) for value in positive]

    y_modes = build_y(32)

    def relative_parts(zeta: float, modes: list[tuple[float, int]]) -> np.ndarray:
        total = np.zeros(3)
        for angular, (joined, disconnected) in enumerate(radial, start=1):
            for y_squared, degeneracy in modes:
                total += 4.0 * angular * degeneracy * (
                    _mode_parts(zeta, joined + y_squared)
                    - _mode_parts(zeta, disconnected + y_squared)
                )
        return total

    zeta_grid = np.geomspace(CHILD_THRESHOLD + 1.0e-5, 200.0, 257)
    rows = []
    for zeta in zeta_grid:
        determinant, anomaly, physical = relative_parts(float(zeta), y_modes)
        rows.append(
            {
                "zeta": float(zeta),
                "determinant_log_scale_derivative": float(determinant),
                "pure_cutoff_anomaly_derivative": float(anomaly),
                "anomaly_compensated_child_link_derivative": float(physical),
            }
        )

    warped = records["NSC-2-ZETA1-WARPED-Y"]
    prior_candidate = float(
        warped["current_warped_y_candidate"]["zeta"]
    )
    candidate_parts = relative_parts(prior_candidate, y_modes)
    y_convergence = []
    for intervals in (16, 32, 64, 128, 256):
        values = relative_parts(prior_candidate, build_y(intervals))
        y_convergence.append(
            {
                "intervals": intervals,
                "determinant_derivative": float(values[0]),
                "pure_cutoff_anomaly": float(values[1]),
                "child_link_derivative": float(values[2]),
            }
        )

    y_differences = [
        abs(
            y_convergence[index - 1]["child_link_derivative"]
            - y_convergence[index]["child_link_derivative"]
        )
        for index in range(1, len(y_convergence))
    ]
    y_orders = [
        float(np.log2(y_differences[index - 1] / y_differences[index]))
        for index in range(1, len(y_differences))
    ]

    physical_values = np.asarray(
        [item["anomaly_compensated_child_link_derivative"] for item in rows]
    )
    maximum_physical = float(np.max(physical_values))
    minimum_physical = float(np.min(physical_values))

    # The per-eigenvalue identity is exact:
    # q=(lambda+1)/zeta-a/zeta^2 and
    # g=(lambda+1)/zeta-2a/zeta^2, hence g-q=-a/zeta^2.
    record = {
        "artifact_id": "NSC-2-ZETA1-ANOMALY-DECOMPOSITION",
        "schema": "NSC-2-ZETA1-ANOMALY-DECOMPOSITION-v1",
        "classification": "compensating_the_pure_cutoff_anomaly_removes_the_regulated_determinant_scale_root_and_leaves_monotone_child_link_flow",
        "authenticated_inputs": authenticated,
        "exact_decomposition": {
            "proper_time_quotient": "q=(lambda+1)/zeta-(3*pi/2)/zeta^2",
            "determinant_gradient_numerator": "g=(lambda+1)/zeta-3*pi/zeta^2",
            "identity": "g-q=-(3*pi/2)/zeta^2",
            "determinant_derivative_per_mode": "exp(-q)*g/(2*q)",
            "pure_cutoff_scale_derivative_per_mode": "exp(-q)/2",
            "anomaly_compensated_physical_derivative_per_mode": "-exp(-q)*(3*pi/2)/(2*zeta^2*q)",
            "meaning": "the_anomaly_cancels_only_cutoff_rescaling_and_retains_the_child_required_change_of_Phi_over_Lambda",
        },
        "scan": {
            "zeta_minimum": float(zeta_grid[0]),
            "zeta_maximum": float(zeta_grid[-1]),
            "count": len(zeta_grid),
            "angular_max": ANGULAR_MAX,
            "warped_y_intervals": 32,
            "rows": rows,
            "minimum_physical_derivative": minimum_physical,
            "maximum_physical_derivative": maximum_physical,
        },
        "prior_determinant_candidate": {
            "zeta": prior_candidate,
            "determinant_derivative": float(candidate_parts[0]),
            "pure_cutoff_anomaly": float(candidate_parts[1]),
            "anomaly_compensated_child_link_derivative": float(
                candidate_parts[2]
            ),
            "remains_stationary_after_compensation": bool(
                abs(candidate_parts[2]) < 1.0e-10
            ),
        },
        "y_convergence": {
            "records": y_convergence,
            "observed_orders": y_orders,
        },
        "causal_result": {
            "uncompensated_determinant_root_preserved_as_diagnostic": True,
            "uncompensated_root_promoted": False,
            "fermion_plus_compensating_anomaly_has_scale_root": False,
            "first_missing_owner": "the_classical_and_anomaly_induced_parent_child_boundary_geometry_action",
            "next_result": "evaluate_the_renormalized_joined_minus_disconnected_geometric_boundary_action_on_the_same_zeta_curve_and_add_its_derivative",
        },
        "gate": {
            "inputs_authenticated": len(authenticated) == len(INPUTS),
            "per_mode_decomposition_exact": True,
            "prior_candidate_not_stationary_after_compensation": bool(
                abs(candidate_parts[2]) > 1.0e-6
            ),
            "physical_derivative_negative_on_scanned_domain": maximum_physical
            < 0.0,
            "warped_y_derivative_second_order": min(y_orders[-2:]) > 1.9,
            "geometric_boundary_action_included": False,
            "zeta_derived": False,
        },
        "nonclaims": {
            "absence_of_fermion_only_root_excludes_full_scale_selection": False,
            "local_boundary_geometry_action_is_zero": False,
            "physical_zeta_promoted": False,
            "observational_input_used": False,
        },
        "terminal": True,
    }
    for key, value in record["gate"].items():
        if key in {"geometric_boundary_action_included", "zeta_derived"}:
            continue
        if value is not True:
            raise RuntimeError(f"ZETA1 anomaly decomposition failed: {key}")
    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
