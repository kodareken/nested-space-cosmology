#!/usr/bin/env python3
"""Insert the exact conformal compact-y Dirac spectrum into ZETA1."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from math import log, pi, sqrt
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import svdvals
from scipy.optimize import brentq


ROOT = Path(__file__).resolve().parents[1]
for value in (str(ROOT), str(ROOT / "src")):
    if value not in sys.path:
        sys.path.insert(0, value)

from recursive_horizons.evidence_io import canonical_json_bytes  # noqa: E402


OUTPUT = ROOT / "results/nsc-2-zeta1-warped-y.json"
DETERMINANT_RUNNER = ROOT / "scripts/run_nsc_zeta1_regulated_determinant.py"
INPUTS = (
    ROOT / "results/nsc-2-zeta1-orbifold-parity.json",
    ROOT / "results/nsc-2-zeta1-regulated-determinant.json",
    ROOT / "results/nsc-2-zeta1-foliation.json",
    DETERMINANT_RUNNER,
)
CHILD_THRESHOLD = 3.0 * pi / 2.0
ANGULAR_MAX = 12


def _load_owner():
    specification = importlib.util.spec_from_file_location(
        "_nsc_zeta1_determinant_owner", DETERMINANT_RUNNER
    )
    if specification is None or specification.loader is None:
        raise RuntimeError("cannot load determinant owner")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def _authenticate() -> list[dict[str, object]]:
    authenticated = []
    for path in INPUTS:
        raw = path.read_bytes()
        if path.suffix == ".json":
            item = json.loads(raw)
            if item.get("terminal") is not True:
                raise RuntimeError(f"nonterminal input: {path.relative_to(ROOT)}")
            artifact_id = item["artifact_id"]
        else:
            artifact_id = None
        authenticated.append(
            {
                "path": str(path.relative_to(ROOT)),
                "sha256": hashlib.sha256(raw).hexdigest(),
                "artifact_id": artifact_id,
            }
        )
    return authenticated


def _warped_y_squared_modes(intervals: int) -> np.ndarray:
    spacing = 2.0 / intervals
    nodes = np.linspace(-1.0, 1.0, intervals + 1)
    edges = 0.5 * (nodes[:-1] + nodes[1:])
    sigma_nodes = -18.0 * nodes**2 / 1015.0
    sigma_edges = -18.0 * edges**2 / 1015.0

    incidence = np.zeros((intervals, intervals + 1))
    index = np.arange(intervals)
    incidence[index, index] = -1.0 / spacing
    incidence[index, index + 1] = 1.0 / spacing
    conformal_dirac = (
        np.diag(np.exp(-sigma_edges / 2.0))
        @ incidence
        @ np.diag(np.exp(-sigma_nodes / 2.0))
    )

    node_weights = np.full(intervals + 1, spacing)
    node_weights[[0, -1]] = spacing / 2.0
    edge_weights = np.full(intervals, spacing)
    orthonormal_dirac = (
        np.diag(np.sqrt(edge_weights))
        @ conformal_dirac
        @ np.diag(1.0 / np.sqrt(node_weights))
    )
    singular = np.sort(svdvals(orthonormal_dirac))
    return singular**2


def _mode_derivative(zeta: float, values: np.ndarray) -> float:
    quotient = (values + 1.0) / zeta - CHILD_THRESHOLD / zeta**2
    gradient = (values + 1.0) / zeta - 2.0 * CHILD_THRESHOLD / zeta**2
    return float(np.sum(0.5 * np.exp(-quotient) * gradient / quotient))


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("ZETA1 warped-y output already exists")
    authenticated = _authenticate()
    owner = _load_owner()
    radial_sectors = [
        owner._sector_spectrum(30.0, 750, angular)
        for angular in range(1, ANGULAR_MAX + 1)
    ]

    convergence = []
    for intervals in (16, 32, 64, 128, 256):
        positive_y_squared = _warped_y_squared_modes(intervals)
        y_modes = [(0.0, 1)] + [
            (float(value), 2) for value in positive_y_squared
        ]

        def derivative(zeta: float) -> float:
            total = 0.0
            for angular, (joined, disconnected) in enumerate(
                radial_sectors, start=1
            ):
                for y_squared, degeneracy in y_modes:
                    total += 4.0 * angular * degeneracy * (
                        _mode_derivative(zeta, joined + y_squared)
                        - _mode_derivative(
                            zeta, disconnected + y_squared
                        )
                    )
            return total

        lower = CHILD_THRESHOLD + 1.0e-5
        upper = 8.0
        lower_derivative = derivative(lower)
        upper_derivative = derivative(upper)
        if not lower_derivative < 0.0 < upper_derivative:
            raise RuntimeError("warped-y determinant root is not bracketed")
        root = brentq(derivative, lower, upper, xtol=1.0e-13, rtol=1.0e-13)
        convergence.append(
            {
                "intervals": intervals,
                "spacing": 2.0 / intervals,
                "first_positive_y_eigenvalue": float(
                    positive_y_squared[0]
                ),
                "zeta": root,
                "mu_squared": 1.0 - CHILD_THRESHOLD / root,
                "mu": sqrt(1.0 - CHILD_THRESHOLD / root),
                "root_residual": derivative(root),
                "lower_derivative": lower_derivative,
                "upper_derivative": upper_derivative,
            }
        )

    roots = [item["zeta"] for item in convergence]
    orders = [
        log(
            abs(
                (roots[index - 2] - roots[index - 1])
                / (roots[index - 1] - roots[index])
            ),
            2.0,
        )
        for index in range(2, len(roots))
    ]
    extrapolated = (4.0 * roots[-1] - roots[-2]) / 3.0
    error = abs(extrapolated - roots[-1])
    mu_squared = 1.0 - CHILD_THRESHOLD / extrapolated
    record = {
        "artifact_id": "NSC-2-ZETA1-WARPED-Y",
        "schema": "NSC-2-ZETA1-WARPED-Y-v1",
        "classification": "the_exact_conformal_compact_y_Dirac_spectrum_preserves_the_orbifold_scale_minimum",
        "authenticated_inputs": authenticated,
        "operator": {
            "sigma": "-18*y^2/1015",
            "flat_measure_Dirac": "Q_y=exp(-sigma_edge/2)*d_y*exp(-sigma_node/2)",
            "discretization": "weighted_node_edge_incidence_with_trapezoid_node_mass",
            "orbifold_spectrum": "one_kernel_zero_mode_plus_two_copies_of_each_positive_singular_value",
            "factorized_radial_approximation": True,
        },
        "convergence": convergence,
        "observed_orders": orders,
        "current_warped_y_candidate": {
            "zeta": extrapolated,
            "error_estimate": error,
            "zeta_minus_3pi_over_2": extrapolated - CHILD_THRESHOLD,
            "mu_squared": mu_squared,
            "mu": sqrt(mu_squared),
            "root_remains_on_physical_side": extrapolated > CHILD_THRESHOLD,
        },
        "causal_result": {
            "orbifold_zero_mode_retained": True,
            "exact_y_warp_changes_but_does_not_remove_root": True,
            "next_owner": "nonseparable_warp_times_radial_operator_plus_horizon_lapse_and_shift",
            "next_result": "construct_the_sparse_two_dimensional_Dirac_operator_and_repeat_the_relative_determinant_root",
        },
        "gate": {
            "inputs_authenticated": len(authenticated) == len(INPUTS),
            "all_roots_bracketed": all(
                item["lower_derivative"] < 0.0 < item["upper_derivative"]
                for item in convergence
            ),
            "root_residuals_below_1e_9": all(
                abs(item["root_residual"]) < 1.0e-9
                for item in convergence
            ),
            "second_order_y_convergence": min(orders[-2:]) > 1.9,
            "warped_root_above_child_threshold": extrapolated
            > CHILD_THRESHOLD,
            "full_nonseparable_transition_zeta_derived": False,
        },
        "nonclaims": {
            "factorized_radial_y_sum_is_the_full_Dirac_operator": False,
            "lapse_and_shift_included": False,
            "local_matrix_anomaly_included": False,
            "physical_zeta_promoted": False,
        },
        "terminal": True,
    }
    for key, value in record["gate"].items():
        if key == "full_nonseparable_transition_zeta_derived":
            continue
        if value is not True:
            raise RuntimeError(f"ZETA1 warped-y gate failed: {key}")
    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
