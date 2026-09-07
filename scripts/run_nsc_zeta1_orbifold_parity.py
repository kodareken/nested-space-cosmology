#!/usr/bin/env python3
"""Select the compact-y parity and compute its unwarped ZETA1 root."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from math import pi, sqrt
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import brentq


ROOT = Path(__file__).resolve().parents[1]
for value in (str(ROOT), str(ROOT / "src")):
    if value not in sys.path:
        sys.path.insert(0, value)

from recursive_horizons.evidence_io import canonical_json_bytes  # noqa: E402


OUTPUT = ROOT / "results/nsc-2-zeta1-orbifold-parity.json"
DETERMINANT_RUNNER = ROOT / "scripts/run_nsc_zeta1_regulated_determinant.py"
INPUTS = (
    ROOT / "results/nsc-2-zeta1-y-boundary-sensitivity.json",
    ROOT / "results/nsc-2-zeta1-regulated-determinant.json",
    ROOT / "results/nsc-1-s-one-warped-resolution-bind.json",
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


def _mode_derivative(zeta: float, values: np.ndarray) -> float:
    quotient = (values + 1.0) / zeta - CHILD_THRESHOLD / zeta**2
    gradient = (values + 1.0) / zeta - 2.0 * CHILD_THRESHOLD / zeta**2
    return float(np.sum(0.5 * np.exp(-quotient) * gradient / quotient))


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("ZETA1 orbifold-parity output already exists")
    authenticated = _authenticate()
    owner = _load_owner()
    sectors = [
        owner._sector_spectrum(30.0, 750, angular)
        for angular in range(1, ANGULAR_MAX + 1)
    ]

    convergence = []
    for y_max in (4, 8, 16, 32):
        # One even chirality has q_0=0 and Neumann modes.  The opposite odd
        # chirality has Dirichlet modes.  Nonzero q_n=n*pi/2 are therefore
        # paired, while the zero mode is single.
        y_modes = [(0.0, 1)] + [
            (n * pi / 2.0, 2) for n in range(1, y_max + 1)
        ]

        def derivative(zeta: float) -> float:
            total = 0.0
            for angular, (joined, disconnected) in enumerate(
                sectors, start=1
            ):
                for y_value, y_degeneracy in y_modes:
                    total += 4.0 * angular * y_degeneracy * (
                        _mode_derivative(zeta, joined + y_value**2)
                        - _mode_derivative(
                            zeta, disconnected + y_value**2
                        )
                    )
            return total

        lower = CHILD_THRESHOLD + 1.0e-5
        upper = 8.0
        lower_derivative = derivative(lower)
        upper_derivative = derivative(upper)
        if not lower_derivative < 0.0 < upper_derivative:
            raise RuntimeError("orbifold determinant root is not bracketed")
        root = brentq(derivative, lower, upper, xtol=1.0e-13, rtol=1.0e-13)
        convergence.append(
            {
                "y_max": y_max,
                "zeta": root,
                "mu_squared": 1.0 - CHILD_THRESHOLD / root,
                "mu": sqrt(1.0 - CHILD_THRESHOLD / root),
                "lower_derivative": lower_derivative,
                "upper_derivative": upper_derivative,
                "root_residual": derivative(root),
            }
        )

    selected = convergence[-1]
    y_tail = abs(convergence[-1]["zeta"] - convergence[-2]["zeta"])
    record = {
        "artifact_id": "NSC-2-ZETA1-ORBIFOLD-PARITY",
        "schema": "NSC-2-ZETA1-ORBIFOLD-PARITY-v1",
        "classification": "joint_fermion_metric_orbifold_parity_selects_one_even_zero_mode_and_an_unwarped_child_allowed_scale_root",
        "authenticated_inputs": authenticated,
        "parity_derivation": {
            "geometry": "S1_over_Z2_reflection_in_the_compact_outside_direction",
            "metric_tangent_components": "even_with_Neumann_zero_mode",
            "metric_mixed_normal_components": "odd_with_Dirichlet_boundary_values",
            "fermion_condition": "Psi(-y)=plus_or_minus_gamma_y*Psi(y)",
            "even_chirality": "Neumann_tower_including_one_zero_mode",
            "odd_chirality": "Dirichlet_tower_without_zero_mode",
            "combined_spectrum": "q_0=0_once; q_n=n*pi/2_twice_for_n>=1",
            "selection_reason": "self_adjoint_orbifold_parity_and_the_existing_even_local_gravity_zero_mode_not_preservation_of_the_scale_root",
        },
        "convergence": convergence,
        "current_unwarped_orbifold_candidate": {
            "zeta": selected["zeta"],
            "y_tail": y_tail,
            "zeta_minus_3pi_over_2": selected["zeta"] - CHILD_THRESHOLD,
            "mu_squared": selected["mu_squared"],
            "mu": selected["mu"],
        },
        "next_result": "replace_the_factorized_q_y_spectrum_with_the_exact_warped_conformal_Dirac_operator_and_retest_the_regulated_scale_minimum",
        "gate": {
            "inputs_authenticated": len(authenticated) == len(INPUTS),
            "one_chiral_zero_mode": True,
            "nonzero_modes_paired": True,
            "root_bracketed_at_all_y_cutoffs": all(
                item["lower_derivative"] < 0.0 < item["upper_derivative"]
                for item in convergence
            ),
            "root_residuals_below_1e_10": all(
                abs(item["root_residual"]) < 1.0e-10
                for item in convergence
            ),
            "y_tower_converged_below_1e_10": y_tail < 1.0e-10,
            "candidate_above_child_threshold": selected["zeta"]
            > CHILD_THRESHOLD,
            "exact_warped_zeta_derived": False,
        },
        "nonclaims": {
            "unwarped_factorization_is_the_exact_y_operator": False,
            "lapse_and_shift_included": False,
            "local_matrix_anomaly_included": False,
            "physical_zeta_promoted": False,
        },
        "terminal": True,
    }
    for key, value in record["gate"].items():
        if key == "exact_warped_zeta_derived":
            continue
        if value is not True:
            raise RuntimeError(f"ZETA1 orbifold parity failed: {key}")
    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
