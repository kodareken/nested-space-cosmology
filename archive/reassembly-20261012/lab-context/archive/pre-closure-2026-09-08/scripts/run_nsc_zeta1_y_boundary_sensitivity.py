#!/usr/bin/env python3
"""Measure ZETA1 determinant sensitivity to compact-y boundary spectra."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from math import pi
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import brentq


ROOT = Path(__file__).resolve().parents[1]
for value in (str(ROOT), str(ROOT / "src")):
    if value not in sys.path:
        sys.path.insert(0, value)

from recursive_horizons.evidence_io import canonical_json_bytes  # noqa: E402


OUTPUT = ROOT / "results/nsc-2-zeta1-y-boundary-sensitivity.json"
DETERMINANT_RUNNER = ROOT / "scripts/run_nsc_zeta1_regulated_determinant.py"
INPUTS = (
    ROOT / "results/nsc-2-zeta1-regulated-determinant.json",
    ROOT / "results/nsc-2-zeta1-self-adjoint-domain.json",
    ROOT / "results/nsc-1-s-one-warped-resolution-bind.json",
    DETERMINANT_RUNNER,
)
CHILD_THRESHOLD = 3.0 * pi / 2.0
ANGULAR_MAX = 12
Y_MAX = 16


def _load_determinant_owner():
    specification = importlib.util.spec_from_file_location(
        "_nsc_zeta1_determinant_owner", DETERMINANT_RUNNER
    )
    if specification is None or specification.loader is None:
        raise RuntimeError("cannot load ZETA1 determinant owner")
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


def _y_modes(kind: str, maximum: int) -> list[tuple[float, int]]:
    if kind == "periodic":
        return [(0.0, 1)] + [(n * pi, 2) for n in range(1, maximum + 1)]
    if kind == "neumann":
        return [(0.0, 1)] + [
            (n * pi / 2.0, 1) for n in range(1, maximum + 1)
        ]
    if kind == "dirichlet":
        return [(n * pi / 2.0, 1) for n in range(1, maximum + 1)]
    if kind == "MIT_bag":
        return [
            ((n + 0.5) * pi / 2.0, 1) for n in range(maximum)
        ]
    raise ValueError(f"unknown y boundary kind: {kind}")


def _mode_derivative(zeta: float, values: np.ndarray) -> float:
    quotient = (values + 1.0) / zeta - CHILD_THRESHOLD / zeta**2
    gradient = (values + 1.0) / zeta - 2.0 * CHILD_THRESHOLD / zeta**2
    return float(np.sum(0.5 * np.exp(-quotient) * gradient / quotient))


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("ZETA1 y-boundary sensitivity output already exists")
    authenticated = _authenticate()
    owner = _load_determinant_owner()
    sectors = [
        owner._sector_spectrum(30.0, 750, angular)
        for angular in range(1, ANGULAR_MAX + 1)
    ]

    boundary_results = []
    for kind in ("periodic", "neumann", "dirichlet", "MIT_bag"):
        convergence = []
        for y_max in (4, 8, Y_MAX):
            y_modes = _y_modes(kind, y_max)

            def derivative(zeta: float) -> float:
                total = 0.0
                for angular, (joined, disconnected) in enumerate(
                    sectors, start=1
                ):
                    for y_value, y_degeneracy in y_modes:
                        joined_value = _mode_derivative(
                            zeta, joined + y_value**2
                        )
                        disconnected_value = _mode_derivative(
                            zeta, disconnected + y_value**2
                        )
                        total += (
                            4.0
                            * angular
                            * y_degeneracy
                            * (joined_value - disconnected_value)
                        )
                return total

            zeta_grid = np.geomspace(
                CHILD_THRESHOLD + 1.0e-5, 20.0, 401
            )
            derivatives = np.asarray(
                [derivative(float(zeta)) for zeta in zeta_grid]
            )
            changes = np.flatnonzero(derivatives[:-1] * derivatives[1:] < 0.0)
            roots = [
                brentq(
                    derivative,
                    float(zeta_grid[index]),
                    float(zeta_grid[index + 1]),
                )
                for index in changes
            ]
            minimum_index = int(np.argmin(derivatives))
            convergence.append(
                {
                    "y_max": y_max,
                    "zero_mode_present": y_modes[0][0] == 0.0,
                    "root_count": len(roots),
                    "roots": roots,
                    "minimum_derivative": float(derivatives[minimum_index]),
                    "minimum_derivative_zeta": float(zeta_grid[minimum_index]),
                    "maximum_derivative": float(np.max(derivatives)),
                }
            )
        boundary_results.append(
            {
                "boundary_kind": kind,
                "convergence": convergence,
                "final": convergence[-1],
            }
        )

    by_kind = {item["boundary_kind"]: item for item in boundary_results}
    periodic = by_kind["periodic"]["final"]
    neumann = by_kind["neumann"]["final"]
    dirichlet = by_kind["dirichlet"]["final"]
    mit = by_kind["MIT_bag"]["final"]
    record = {
        "artifact_id": "NSC-2-ZETA1-Y-BOUNDARY-SENSITIVITY",
        "schema": "NSC-2-ZETA1-Y-BOUNDARY-SENSITIVITY-v1",
        "classification": "the_regulated_scale_minimum_survives_only_compact_y_spectra_with_an_even_zero_mode_in_the_factorized_control",
        "authenticated_inputs": authenticated,
        "method": {
            "radial_angular_owner": "NSC-2-ZETA1-REGULATED-DETERMINANT",
            "angular_max": ANGULAR_MAX,
            "factorized_control": "lambda_total=lambda_radial_partner+q_y^2",
            "warp_included": False,
            "boundary_spectra": {
                "periodic": "q_0=0_then_q_n=n*pi_with_double_nonzero_degeneracy",
                "neumann": "q_n=n*pi/2_including_n=0",
                "dirichlet": "q_n=n*pi/2_for_n>=1",
                "MIT_bag": "q_n=(n+1/2)*pi/2",
            },
            "y_mode_cutoffs": [4, 8, Y_MAX],
        },
        "boundary_results": boundary_results,
        "causal_result": {
            "periodic_zero_mode_root": periodic["roots"],
            "neumann_zero_mode_root": neumann["roots"],
            "dirichlet_root": dirichlet["roots"],
            "MIT_bag_root": mit["roots"],
            "root_existence_tracks_zero_mode": True,
            "existing_local_gravity_requirement": "the_warped_resolution_result_requires_an_even_normalizable_zero_mode_for_local_4D_gravity",
            "next_owner": "derive_the_joint_fermion_metric_orbifold_parities_and_APS_Calderon_boundary_data_then_run_the_exact_warped_y_operator",
        },
        "gate": {
            "inputs_authenticated": len(authenticated) == len(INPUTS),
            "periodic_root_present": len(periodic["roots"]) == 1,
            "neumann_root_present": len(neumann["roots"]) == 1,
            "dirichlet_root_absent": not dirichlet["roots"],
            "MIT_bag_root_absent": not mit["roots"],
            "all_y_sums_converged_by_16": all(
                item["convergence"][-1]["roots"]
                == item["convergence"][-2]["roots"]
                for item in boundary_results
            ),
            "physical_y_boundary_selected": False,
            "full_warped_zeta_derived": False,
        },
        "nonclaims": {
            "zero_mode_selected_to_force_a_root": False,
            "factorized_unwarped_result_is_the_full_y_operator": False,
            "fermion_and_metric_parities_already_identical": False,
            "physical_zeta_promoted": False,
        },
        "terminal": True,
    }
    for key, value in record["gate"].items():
        if key in {"physical_y_boundary_selected", "full_warped_zeta_derived"}:
            continue
        if value is not True:
            raise RuntimeError(f"ZETA1 y-boundary sensitivity failed: {key}")
    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
