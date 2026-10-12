#!/usr/bin/env python3
"""Scan the induced-room normalization for real interface characteristics."""

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


OUTPUT = ROOT / "results/nsc-1-s-one-interface-characteristic-scan.json"
RANK_RESULT = ROOT / "results/nsc-1-s-one-induced-interface-rank.json"
CRITICAL_RESULT = ROOT / "results/nsc-1-s-one-gauss-bonnet-criticality.json"
KINETIC_RUNNER = ROOT / "scripts/run_s_one_gauss_bonnet_kinetic_rank.py"
INTERFACE_RUNNER = ROOT / "scripts/run_s_one_induced_interface_rank.py"


def _load(path: Path, name: str):
    specification = importlib.util.spec_from_file_location(name, path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"cannot load implementation owner: {path}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def _symbol_matrix(
    kinetic_owner,
    interface_owner,
    background,
    rho: float,
    covector: np.ndarray,
    alpha: float,
) -> tuple[list[tuple[int, int]], np.ndarray, np.ndarray, int]:
    diagonal, riemann = background(rho, 0.0, np.pi / 2)
    metric = np.diag(np.asarray(diagonal, dtype=float).reshape(5))
    riemann = np.asarray(riemann, dtype=float)
    time_index = [index for index in range(5) if metric[index, index] > 0][0]
    spatial = [index for index in range(5) if index != time_index]
    pairs = [
        (left, right)
        for offset, left in enumerate(spatial)
        for right in spatial[offset:]
    ]
    bulk = np.zeros((10, 10))
    induced = np.zeros((10, 10))
    for column, (left, right) in enumerate(pairs):
        perturbation = np.zeros((5, 5))
        perturbation[left, right] = 1.0
        perturbation[right, left] = 1.0
        bulk_output = kinetic_owner._principal_action(
            metric, riemann, covector, perturbation, alpha
        )
        if left < 4 and right < 4:
            induced_output = interface_owner._einstein_principal(
                metric[:4, :4], covector[:4], perturbation[:4, :4]
            )
        else:
            induced_output = np.zeros((4, 4))
        for row, (first, second) in enumerate(pairs):
            bulk[row, column] = bulk_output[first, second]
            if first < 4 and second < 4:
                induced[row, column] = induced_output[first, second]
    return pairs, bulk, induced, time_index


def _pencil(
    kinetic_owner,
    interface_owner,
    background,
    rho: float,
    spatial_offset: int,
    alpha: float,
):
    diagonal, _ = background(rho, 0.0, np.pi / 2)
    metric = np.diag(np.asarray(diagonal, dtype=float).reshape(5))
    time_index = [index for index in range(5) if metric[index, index] > 0][0]
    spatial = [index for index in range(5) if index != time_index]
    time_covector = np.zeros(5)
    time_covector[time_index] = 1.0
    space_covector = np.zeros(5)
    space_covector[spatial[spatial_offset]] = 1.0
    _, bulk_a, induced_a, _ = _symbol_matrix(
        kinetic_owner,
        interface_owner,
        background,
        rho,
        time_covector,
        alpha,
    )
    _, bulk_c, induced_c, _ = _symbol_matrix(
        kinetic_owner,
        interface_owner,
        background,
        rho,
        space_covector,
        alpha,
    )
    _, bulk_sum, induced_sum, _ = _symbol_matrix(
        kinetic_owner,
        interface_owner,
        background,
        rho,
        time_covector + space_covector,
        alpha,
    )
    return {
        "rho": rho,
        "time_index": time_index,
        "space_index": spatial[spatial_offset],
        "bulk_a": bulk_a,
        "bulk_b": bulk_sum - bulk_a - bulk_c,
        "bulk_c": bulk_c,
        "induced_a": induced_a,
        "induced_b": induced_sum - induced_a - induced_c,
        "induced_c": induced_c,
    }


def _roots(pencil, normalization: float, *, include_bulk_gauss_bonnet: bool):
    if include_bulk_gauss_bonnet:
        matrix_a = pencil["bulk_a"] + normalization * pencil["induced_a"]
        matrix_b = pencil["bulk_b"] + normalization * pencil["induced_b"]
        matrix_c = pencil["bulk_c"] + normalization * pencil["induced_c"]
    else:
        matrix_a = pencil["bulk_a"]
        matrix_b = pencil["bulk_b"]
        matrix_c = pencil["bulk_c"]
    inverse_a_b = np.linalg.solve(matrix_a, matrix_b)
    inverse_a_c = np.linalg.solve(matrix_a, matrix_c)
    companion = np.block(
        [
            [np.zeros((10, 10)), np.eye(10)],
            [-inverse_a_c, -inverse_a_b],
        ]
    )
    return np.linalg.eigvals(companion)


def _physical_roots(roots: np.ndarray) -> np.ndarray:
    """Select the ten 5D graviton roots, excluding the gauge-chart zero sector."""

    order = np.argsort(np.abs(roots))[::-1]
    return roots[order[:10]]


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("interface characteristic-scan output already exists")

    authenticated = []
    for path in (
        RANK_RESULT,
        CRITICAL_RESULT,
        KINETIC_RUNNER,
        INTERFACE_RUNNER,
    ):
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

    kinetic_owner = _load(KINETIC_RUNNER, "_nsc_kinetic_owner")
    interface_owner = _load(INTERFACE_RUNNER, "_nsc_interface_owner")
    background = kinetic_owner._background_function()
    rho_samples = (-2.0, -1.0, 0.0, 0.5, 1.0, 2.0)
    candidate_pencils = [
        _pencil(
            kinetic_owner,
            interface_owner,
            background,
            rho,
            direction,
            1015.0 / 144.0,
        )
        for rho in rho_samples
        for direction in range(4)
    ]
    gr_pencils = [
        _pencil(
            kinetic_owner,
            interface_owner,
            background,
            rho,
            direction,
            0.0,
        )
        for rho in rho_samples
        for direction in range(4)
    ]

    def score(pencil, normalization: float, candidate: bool) -> float:
        roots = _physical_roots(
            _roots(
                pencil,
                normalization,
                include_bulk_gauss_bonnet=candidate,
            )
        )
        return float(np.max(np.abs(roots.imag)))

    gr_scores = [score(pencil, 0.0, False) for pencil in gr_pencils]
    normalizations = np.logspace(-4, 6, 401)
    scan = []
    best = None
    best_tangent = None
    for normalization in normalizations:
        scores = [
            score(pencil, float(normalization), True)
            for pencil in candidate_pencils
        ]
        tangent_scores = [
            value
            for value, pencil in zip(scores, candidate_pencils, strict=True)
            if pencil["space_index"] != 4
        ]
        worst_index = int(np.argmax(scores))
        item = {
            "normalization": float(normalization),
            "maximum_imaginary_root": max(scores),
            "maximum_tangent_imaginary_root": max(tangent_scores),
            "worst_rho": candidate_pencils[worst_index]["rho"],
            "worst_space_index": candidate_pencils[worst_index]["space_index"],
        }
        scan.append(item)
        if best is None or item["maximum_imaginary_root"] < best[
            "maximum_imaginary_root"
        ]:
            best = item
        if best_tangent is None or item[
            "maximum_tangent_imaginary_root"
        ] < best_tangent["maximum_tangent_imaginary_root"]:
            best_tangent = item
    if best is None or best_tangent is None:
        raise RuntimeError("empty characteristic scan")

    selected = []
    for target in (1.0e-4, 1.0e-2, 1.0, 10.0, 1.0e2, 1.0e4, 1.0e6):
        selected.append(min(scan, key=lambda item: abs(item["normalization"] - target)))

    record = {
        "artifact_id": "NSC-1-S-ONE-INTERFACE-CHARACTERISTIC-SCAN",
        "schema": "NSC-1-S-ONE-INTERFACE-CHARACTERISTIC-SCAN-v1",
        "classification": "no_positive_induced_Einstein_normalization_restores_real_sampled_characteristics_without_the_off_diagonal_sheet_mode",
        "authenticated_inputs": authenticated,
        "method": {
            "principal_symbol": "direct_covariant_EGB_curvature_variation_plus_tangential_induced_Einstein_block",
            "gauge": "nonzero_time_frequency_gauge_with_ten_spatial_symmetric_metric_variables",
            "quadratic_pencil": "P(omega,k)=A*omega^2+B*omega*k+C*k^2",
            "root_solver": "20_by_20_first_companion_matrix",
            "physical_root_selection": "ten_largest_magnitude_roots_remove_the_ten_GR_gauge_chart_constraint_roots_at_zero",
            "rho_samples": list(rho_samples),
            "spatial_directions_per_point": 4,
            "normalization_grid": {
                "minimum": float(normalizations[0]),
                "maximum": float(normalizations[-1]),
                "count": len(normalizations),
                "spacing": "logarithmic",
            },
        },
        "GR_control": {
            "maximum_imaginary_root": max(gr_scores),
            "all_roots_real_tolerance": 1.0e-7,
            "passes": max(gr_scores) < 1.0e-7,
            "root_structure": "five_real_plus_minus_metric_null_pairs_after_removing_ten_constraint_zero_roots",
        },
        "candidate_scan": {
            "best_all_directions": best,
            "best_tangent_directions": best_tangent,
            "selected_grid_points": selected,
            "any_all_direction_pass": any(
                item["maximum_imaginary_root"] < 1.0e-7 for item in scan
            ),
            "any_tangent_direction_pass": any(
                item["maximum_tangent_imaginary_root"] < 1.0e-7
                for item in scan
            ),
        },
        "causal_result": {
            "time_principal_rank_restored": True,
            "real_characteristics_restored": False,
            "normalization_rescaling_is_not_the_fix": True,
            "first_missing_operator": "off_diagonal_parent_child_Phi_kinetic_and_gradient_principal_block",
            "reason": "the_one_operator_formula_contains_Phi_but_the_tested_interface_pencil_included_only_diagonal_bulk_and_room_metric_blocks",
        },
        "scope": {
            "development_scan_not_final_hyperbolicity_proof": True,
            "mixed_dimensional_normal_characteristics_require_interface_boundary_analysis": True,
            "tangent_complex_roots_already_exclude_the_tested_diagonal_only_pencil": True,
            "Nested_Space_theory_rejected": False,
        },
        "next_result": "derive_the_Phi_principal_block_from_the_two_sheet_spectral_action_and_solve_the_coupled_metric_Phi_characteristic_pencil",
        "gate": {
            "inputs_authenticated": len(authenticated) == 4,
            "GR_control_real": max(gr_scores) < 1.0e-7,
            "normalization_grid_complete": len(scan) == 401,
            "no_all_direction_normalization_passes": not any(
                item["maximum_imaginary_root"] < 1.0e-7 for item in scan
            ),
            "no_tangent_normalization_passes": not any(
                item["maximum_tangent_imaginary_root"] < 1.0e-7
                for item in scan
            ),
            "off_diagonal_Phi_characteristics_tested": False,
        },
        "terminal": True,
    }
    for key, value in record["gate"].items():
        if key != "off_diagonal_Phi_characteristics_tested" and value is not True:
            raise RuntimeError(f"interface characteristic scan failed: {key}")
    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
