#!/usr/bin/env python3
"""Test whether the induced spectral Einstein block restores critical interface rank."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
for value in (str(ROOT), str(ROOT / "src")):
    if value not in sys.path:
        sys.path.insert(0, value)

from recursive_horizons.evidence_io import canonical_json_bytes  # noqa: E402


OUTPUT = ROOT / "results/nsc-1-s-one-induced-interface-rank.json"
KINETIC = ROOT / "results/nsc-1-s-one-gauss-bonnet-kinetic-rank.json"
SPECTRAL_BOUNDARY = ROOT / "results/nsc-1-s-one-spectral-boundary-constraint.json"
KINETIC_RUNNER = ROOT / "scripts/run_s_one_gauss_bonnet_kinetic_rank.py"


def _load_kinetic_owner():
    specification = importlib.util.spec_from_file_location(
        "_nsc_gauss_bonnet_kinetic_rank", KINETIC_RUNNER
    )
    if specification is None or specification.loader is None:
        raise RuntimeError("cannot load kinetic-rank implementation owner")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def _einstein_principal(
    metric: np.ndarray, covector: np.ndarray, perturbation: np.ndarray
) -> np.ndarray:
    dimension = len(covector)
    riemann = np.zeros((dimension, dimension, dimension, dimension))
    for first in range(dimension):
        for second in range(dimension):
            for third in range(dimension):
                for fourth in range(dimension):
                    riemann[first, second, third, fourth] = 0.5 * (
                        covector[third]
                        * covector[second]
                        * perturbation[first, fourth]
                        + covector[fourth]
                        * covector[first]
                        * perturbation[second, third]
                        - covector[fourth]
                        * covector[second]
                        * perturbation[first, third]
                        - covector[third]
                        * covector[first]
                        * perturbation[second, fourth]
                    )
    inverse = np.linalg.inv(metric)
    ricci = np.einsum("ac,abcd->bd", inverse, riemann)
    ricci_scalar = np.einsum("ab,ab", inverse, ricci)
    return ricci - 0.5 * metric * ricci_scalar


def _interface_record(owner, background, rho: float) -> dict[str, object]:
    diagonal, riemann = background(rho, 0.0, np.pi / 2)
    metric = np.diag(np.asarray(diagonal, dtype=float).reshape(5))
    riemann = np.asarray(riemann, dtype=float)
    positive_directions = [
        index for index in range(5) if metric[index, index] > 0
    ]
    if len(positive_directions) != 1:
        raise RuntimeError("interface point does not have one time direction")
    time_index = positive_directions[0]
    covector = np.zeros(5)
    covector[time_index] = 1.0
    spatial = [index for index in range(5) if index != time_index]
    pairs = [
        (left, right)
        for offset, left in enumerate(spatial)
        for right in spatial[offset:]
    ]
    bulk = np.zeros((10, 10))
    induced = np.zeros((10, 10))
    alpha = 1015.0 / 144.0
    induced_metric = metric[:4, :4]
    induced_covector = covector[:4]
    for column, (left, right) in enumerate(pairs):
        perturbation = np.zeros((5, 5))
        perturbation[left, right] = 1.0
        perturbation[right, left] = 1.0
        bulk_output = owner._principal_action(
            metric, riemann, covector, perturbation, alpha
        )
        if left < 4 and right < 4:
            induced_perturbation = perturbation[:4, :4]
            induced_output = _einstein_principal(
                induced_metric, induced_covector, induced_perturbation
            )
        else:
            induced_output = np.zeros((4, 4))
        for row, (first, second) in enumerate(pairs):
            bulk[row, column] = bulk_output[first, second]
            if first < 4 and second < 4:
                induced[row, column] = induced_output[first, second]

    combined = bulk + induced

    def rank_data(matrix: np.ndarray) -> tuple[int, list[float], float]:
        singular_values = np.linalg.svd(matrix, compute_uv=False)
        tolerance = float(np.max(singular_values) * 1.0e-10)
        rank = int(np.count_nonzero(singular_values > tolerance))
        return rank, [float(value) for value in singular_values], tolerance

    bulk_rank, bulk_singular, bulk_tolerance = rank_data(bulk)
    induced_rank, induced_singular, induced_tolerance = rank_data(induced)
    combined_rank, combined_singular, combined_tolerance = rank_data(combined)

    gauge_residual = 0.0
    for direction in range(5):
        gauge_vector = np.zeros(5)
        gauge_vector[direction] = 1.0
        gauge_perturbation = np.outer(covector, gauge_vector) + np.outer(
            gauge_vector, covector
        )
        bulk_gauge = owner._principal_action(
            metric, riemann, covector, gauge_perturbation, alpha
        )
        if direction < 4:
            induced_gauge = _einstein_principal(
                induced_metric,
                induced_covector,
                gauge_perturbation[:4, :4],
            )
            embedded = np.zeros((5, 5))
            embedded[:4, :4] = induced_gauge
        else:
            embedded = np.zeros((5, 5))
        gauge_residual = max(
            gauge_residual,
            float(np.max(np.abs(bulk_gauge + embedded))),
        )

    return {
        "rho": rho,
        "time_coordinate_index": time_index,
        "spatial_component_order": [f"{left}{right}" for left, right in pairs],
        "critical_bulk": {
            "rank": bulk_rank,
            "singular_values": bulk_singular,
            "rank_tolerance": bulk_tolerance,
        },
        "induced_room_Einstein": {
            "rank": induced_rank,
            "singular_values": induced_singular,
            "rank_tolerance": induced_tolerance,
        },
        "combined_interface": {
            "rank": combined_rank,
            "singular_values": combined_singular,
            "rank_tolerance": combined_tolerance,
            "condition_number": combined_singular[0] / combined_singular[-1],
        },
        "gauge_principal_residual": gauge_residual,
    }


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("induced-interface rank output already exists")

    authenticated = []
    for path in (KINETIC, SPECTRAL_BOUNDARY, KINETIC_RUNNER):
        raw = path.read_bytes()
        if path.suffix == ".json":
            item = json.loads(raw)
            if item.get("terminal") is not True:
                raise RuntimeError(f"nonterminal input: {path.relative_to(ROOT)}")
            artifact_id = item.get("artifact_id")
        else:
            artifact_id = None
        authenticated.append(
            {
                "path": str(path.relative_to(ROOT)),
                "sha256": hashlib.sha256(raw).hexdigest(),
                "artifact_id": artifact_id,
            }
        )

    owner = _load_kinetic_owner()
    background = owner._background_function()
    samples = [
        _interface_record(owner, background, rho)
        for rho in (-2.0, -1.0, 0.0, 0.5, 1.0, 2.0)
    ]

    record = {
        "artifact_id": "NSC-1-S-ONE-INDUCED-INTERFACE-RANK",
        "schema": "NSC-1-S-ONE-INDUCED-INTERFACE-RANK-v1",
        "classification": "the_induced_local_spectral_Einstein_block_restores_full_time_principal_rank_on_the_critical_bulk_sheet",
        "authenticated_inputs": authenticated,
        "interface_action": {
            "bulk": "integral_M5 sqrt(g5)*[R5+(1015/144)*L_star^2*GB5]",
            "local_room": "positive_overall_normalization_times_integral_y=0 sqrt(g4)*R4",
            "support": "the_induced_block_acts_only_on_metric_components_tangent_to_the_room",
            "normalization_used": 1.0,
            "normalization_role": "the_single_allowed_overall_local_room_unit_not_a_fitted_dimensionless_ratio",
        },
        "samples": samples,
        "result": {
            "critical_bulk_ranks": sorted(
                {item["critical_bulk"]["rank"] for item in samples}
            ),
            "induced_room_ranks": sorted(
                {item["induced_room_Einstein"]["rank"] for item in samples}
            ),
            "combined_interface_ranks": sorted(
                {item["combined_interface"]["rank"] for item in samples}
            ),
            "maximum_gauge_residual": max(
                item["gauge_principal_residual"] for item in samples
            ),
            "maximum_condition_number": max(
                item["combined_interface"]["condition_number"]
                for item in samples
            ),
        },
        "physical_interpretation": {
            "bulk_criticality_retained": True,
            "local_propagation_restored_by_room_boundary": True,
            "ordinary_gravity_is_a_boundary_localized_mode": True,
            "dark_black_geometric_invariant_unchanged": True,
            "new_phantom_or_dark_substance_added": False,
        },
        "next_result": "evaluate_the_full_covariant_characteristic_polynomial_of_the_combined_bulk_boundary_system_and_test_common_hyperbolic_cones",
        "nonclaims": {
            "combined_kinetic_signs_positive": False,
            "complete_characteristic_roots_real": False,
            "full_hyperbolicity_proved": False,
            "relative_bulk_boundary_normalization_predicted": False,
        },
        "gate": {
            "inputs_authenticated": len(authenticated) == 3,
            "critical_bulk_rank_five": all(
                item["critical_bulk"]["rank"] == 5 for item in samples
            ),
            "induced_tangential_rank_six": all(
                item["induced_room_Einstein"]["rank"] == 6
                for item in samples
            ),
            "combined_interface_rank_ten": all(
                item["combined_interface"]["rank"] == 10 for item in samples
            ),
            "gauge_principal_residual_zero": all(
                item["gauge_principal_residual"] == 0.0 for item in samples
            ),
            "full_characteristic_and_sign_test_passed": False,
        },
        "terminal": True,
    }
    for key, value in record["gate"].items():
        if key != "full_characteristic_and_sign_test_passed" and value is not True:
            raise RuntimeError(f"induced-interface rank gate failed: {key}")
    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
