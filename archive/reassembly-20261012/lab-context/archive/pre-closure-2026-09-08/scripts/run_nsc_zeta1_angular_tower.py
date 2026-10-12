#!/usr/bin/env python3
"""Test the lowest ZETA1 scale minimum against the spinor angular tower."""

from __future__ import annotations

import hashlib
import json
from math import pi
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import brentq


ROOT = Path(__file__).resolve().parents[1]
for value in (str(ROOT), str(ROOT / "src")):
    if value not in sys.path:
        sys.path.insert(0, value)

from recursive_horizons.evidence_io import canonical_json_bytes  # noqa: E402


OUTPUT = ROOT / "results/nsc-2-zeta1-angular-tower.json"
INPUTS = (
    ROOT / "results/nsc-2-zeta1-lowest-mode.json",
    ROOT / "results/nsc-2-zeta1-self-adjoint-domain.json",
    ROOT / "results/nsc-1-s-one-child-scale-correction.json",
)
BOX_RADIUS = 30.0
HALF_INTERVALS = 750
ANGULAR_MAX = 12
CHILD_THRESHOLD = 3.0 * pi / 2.0


def _authenticate() -> list[dict[str, object]]:
    authenticated = []
    for path in INPUTS:
        raw = path.read_bytes()
        item = json.loads(raw)
        if item.get("terminal") is not True:
            raise RuntimeError(f"nonterminal input: {path.relative_to(ROOT)}")
        authenticated.append(
            {
                "path": str(path.relative_to(ROOT)),
                "sha256": hashlib.sha256(raw).hexdigest(),
                "artifact_id": item["artifact_id"],
            }
        )
    return authenticated


def _sector_spectrum(angular: int) -> tuple[np.ndarray, np.ndarray]:
    spacing = BOX_RADIUS / HALF_INTERVALS
    joined_parts = []
    disconnected_parts = []

    def spectrum(
        left: float, right: float, count: int, partner_sign: int
    ) -> np.ndarray:
        rho = np.linspace(left + spacing, right - spacing, count)
        superpotential = angular / np.sqrt(1.0 + rho**2)
        derivative = -angular * rho / (1.0 + rho**2) ** 1.5
        potential = superpotential**2 + partner_sign * derivative
        diagonal = 2.0 / spacing**2 + potential
        off_diagonal = np.full(count - 1, -1.0 / spacing**2)
        return eigh_tridiagonal(diagonal, off_diagonal, eigvals_only=True)

    for partner_sign in (1, -1):
        joined_parts.append(
            spectrum(
                -BOX_RADIUS,
                BOX_RADIUS,
                2 * HALF_INTERVALS - 1,
                partner_sign,
            )
        )
        disconnected_parts.extend(
            (
                spectrum(
                    -BOX_RADIUS,
                    0.0,
                    HALF_INTERVALS - 1,
                    partner_sign,
                ),
                spectrum(
                    0.0,
                    BOX_RADIUS,
                    HALF_INTERVALS - 1,
                    partner_sign,
                ),
            )
        )
    return np.concatenate(joined_parts), np.concatenate(disconnected_parts)


def _relative_derivative(
    zeta: float, joined: np.ndarray, disconnected: np.ndarray
) -> float:
    mu_squared = 1.0 - CHILD_THRESHOLD / zeta

    def contribution(values: np.ndarray) -> float:
        heat = np.exp(-(values + mu_squared) / zeta)
        logarithmic_gradient = (
            (values + 1.0) / zeta
            - 2.0 * CHILD_THRESHOLD / zeta**2
        )
        return float(np.sum(heat * logarithmic_gradient))

    return contribution(joined) - contribution(disconnected)


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("ZETA1 angular-tower output already exists")
    authenticated = _authenticate()
    sectors = [_sector_spectrum(angular) for angular in range(1, ANGULAR_MAX + 1)]
    zeta_grid = np.geomspace(CHILD_THRESHOLD + 1.0e-6, 200.0, 1025)

    cumulative = []
    for angular_cutoff in range(1, ANGULAR_MAX + 1):
        derivatives = []
        for zeta in zeta_grid:
            value = 0.0
            for angular, (joined, disconnected) in enumerate(
                sectors[:angular_cutoff], start=1
            ):
                # Dirac eigenvalues on S2 are +/-angular, each with
                # multiplicity 2*angular.  The two signs exchange the radial
                # partner pair, giving total multiplicity 4*angular.
                value += 4.0 * angular * _relative_derivative(
                    float(zeta), joined, disconnected
                )
            derivatives.append(value)
        values = np.asarray(derivatives)
        sign_changes = np.flatnonzero(values[:-1] * values[1:] < 0.0)
        roots = []
        for index in sign_changes:
            def cumulative_derivative(zeta: float) -> float:
                return sum(
                    4.0
                    * angular
                    * _relative_derivative(zeta, joined, disconnected)
                    for angular, (joined, disconnected) in enumerate(
                        sectors[:angular_cutoff], start=1
                    )
                )

            roots.append(
                brentq(
                    cumulative_derivative,
                    float(zeta_grid[index]),
                    float(zeta_grid[index + 1]),
                )
            )
        minimum_index = int(np.argmin(values))
        cumulative.append(
            {
                "angular_cutoff": angular_cutoff,
                "included_spinor_degeneracy": sum(
                    4 * angular for angular in range(1, angular_cutoff + 1)
                ),
                "minimum_derivative": float(values[minimum_index]),
                "minimum_derivative_zeta": float(zeta_grid[minimum_index]),
                "maximum_derivative": float(np.max(values)),
                "sign_change_count": int(len(sign_changes)),
                "roots": roots,
            }
        )

    selected_zeta = (5.0, 10.0, 20.0, 50.0, 100.0, 200.0)
    sector_contributions = []
    for angular, (joined, disconnected) in enumerate(sectors, start=1):
        sector_contributions.append(
            {
                "angular": angular,
                "degeneracy": 4 * angular,
                "derivatives": {
                    str(zeta): 4.0
                    * angular
                    * _relative_derivative(zeta, joined, disconnected)
                    for zeta in selected_zeta
                },
            }
        )

    first = cumulative[0]
    second = cumulative[1]
    full = cumulative[-1]
    record = {
        "artifact_id": "NSC-2-ZETA1-ANGULAR-TOWER",
        "schema": "NSC-2-ZETA1-ANGULAR-TOWER-v1",
        "classification": "the_first_omitted_spinor_angular_sector_removes_the_lowest_mode_heat_only_scale_minimum",
        "authenticated_inputs": authenticated,
        "method": {
            "box_radius": BOX_RADIUS,
            "half_intervals": HALF_INTERVALS,
            "spacing": BOX_RADIUS / HALF_INTERVALS,
            "spinor_sphere_spectrum": "eigenvalues_plus_or_minus_n_with_multiplicity_2n_per_sign",
            "angular_range": [1, ANGULAR_MAX],
            "zeta_scan": {
                "minimum": float(zeta_grid[0]),
                "maximum": float(zeta_grid[-1]),
                "count": len(zeta_grid),
                "spacing": "logarithmic",
            },
            "same_child_constraint": "mu_squared=1-3*pi/(2*zeta)",
        },
        "cumulative_results": cumulative,
        "sector_contributions": sector_contributions,
        "causal_result": {
            "lowest_sector_root_reproduced": bool(first["roots"]),
            "root_survives_first_omitted_angular_sector": bool(second["roots"]),
            "heat_only_angular_tower_has_stationary_root": bool(full["roots"]),
            "first_owner": "the_scale_derivative_requires_the_anomaly_consistent_fermionic_determinant_and_warped_y_spectrum_not_the_radial_heat_trace_alone",
            "next_result": "add_the_y_warp_and_compute_the_regulated_joined_disconnected_fermionic_determinant_before_retesting_scale_stationarity",
        },
        "gate": {
            "inputs_authenticated": len(authenticated) == len(INPUTS),
            "lowest_sector_root_present": len(first["roots"]) == 1,
            "first_omitted_sector_removes_root": not second["roots"],
            "angular_12_derivative_positive_on_scanned_domain": full[
                "minimum_derivative"
            ]
            > 0.0,
            "full_anomaly_consistent_zeta_derived": False,
        },
        "nonclaims": {
            "finite_angular_scan_is_an_all_mode_theorem": False,
            "absence_of_heat_only_root_excludes_full_stationarity": False,
            "y_warp_lapse_shift_included": False,
            "fermionic_determinant_and_anomaly_included": False,
        },
        "terminal": True,
    }
    for key, value in record["gate"].items():
        if key == "full_anomaly_consistent_zeta_derived":
            continue
        if value is not True:
            raise RuntimeError(f"ZETA1 angular-tower gate failed: {key}")
    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
