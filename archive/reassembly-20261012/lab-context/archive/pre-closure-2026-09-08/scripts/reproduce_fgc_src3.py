#!/usr/bin/env python3
"""Reproduce SRC3's reference-covariant GR-0 source certificate."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import struct
import subprocess
import tempfile
import tomllib
from typing import Any, Mapping, Sequence

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
import sys

sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.gr0_direct_source import (  # noqa: E402
    gr0_ref1_residual_batch,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    array_content_sha256,
)
from recursive_horizons.fgc.evolution.src2_exact_oracle import (  # noqa: E402
    exact_affine_system,
    exact_gaussian_solve,
    exact_ref1_residual,
)
from recursive_horizons.fgc.evolution.src3_reference_balanced_source import (  # noqa: E402
    diagnose_gr0_reference_balanced_accelerations,
    gr0_reference_balanced_ref1_residual_batch,
)


ARTIFACT_ID = "FGC-1-SRC3"
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-src3.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-src3.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-src3.md"
SOURCE_MODULE = (
    REPOSITORY
    / "src/recursive_horizons/fgc/evolution/src3_reference_balanced_source.py"
)
IMPLEMENTATION = (
    Path(__file__).resolve(),
    SOURCE_MODULE,
    REPOSITORY / "src/recursive_horizons/fgc/evolution/src2_exact_oracle.py",
    REPOSITORY / "src/recursive_horizons/fgc/evolution/gr0_direct_source.py",
    REPOSITORY / "src/recursive_horizons/fgc/spherical_reduction.py",
    REPOSITORY / "src/recursive_horizons/fgc/modified_harmonic_reference.py",
    REPOSITORY / "src/recursive_horizons/fgc/reference_connection.py",
)
FIELD_NAMES = ("u", "p", "q", "p_r", "q_r")
FIELD_ORDER = ["alpha", "shift", "lambda", "R", "phi", "chi"]
EQUATION_ORDER = [
    "metric_tt_mhg",
    "metric_tr_mhg",
    "metric_rr_mhg",
    "metric_theta_theta_mhg",
    "scalar_phi",
    "scalar_chi",
]


EXPECTED_SCOPE = {
    "role": "prospectively_frozen_reference_covariant_binary64_GR0_REF1_source_evaluator_gate",
    "branch": "GR-0",
    "captured_amplitude": "5/2",
    "captured_member": "RK4-8193",
    "captured_grid_points_away_from_center": 8192,
    "existing_CAP1_capture_used_as_development_evidence": True,
    "fresh_GR0_trajectory_read": False,
    "SGBL_trajectory_read": False,
    "FGCQR_trajectory_read": False,
    "collapse_or_trapped_outcome_classified": False,
    "mechanism_question_answered": False,
}
EXPECTED_EVALUATOR = {
    "continuum_equations": "unredefined_ACT1_VAR1_GR0_REF1",
    "arithmetic": "IEEE_754_binary64",
    "identity": "C^a_bc=one_half_g^ad_times_bar_nabla_b_h_dc_plus_bar_nabla_c_h_db_minus_bar_nabla_d_h_bc",
    "perturbation": "h_ab=g_ab-minus-flat_spherical_reference_ab",
    "connection_difference_derivative": "differentiate_the_reference_covariant_identity_before_rounding",
    "ricci_assembly": "flat_reference_connection_plus_C_and_partial_C",
    "reference_connection": "flat_spherical_annulus",
    "tilde_normal_factor": 4,
    "hat_normal_factor": 9,
    "raw_complete_residual_maximum": "1e-12",
    "kinetic_condition_number_maximum": "1e10",
    "maximum_complete_refinement_iterations": 16,
    "strict_raw_gate": True,
    "continuum_equation_change_forbidden": True,
    "reference_geometry_change_forbidden": True,
    "reference_derivative_map_change_forbidden": True,
    "source_threshold_change_forbidden": True,
    "grid_method_CFL_retry_and_physical_input_change_forbidden": True,
    "higher_precision_runtime_fallback_used": False,
}
EXPECTED_PROOF = {
    "exact_reference_must_be_bitwise_zero_at_all_frozen_radii": True,
    "compact_point_stable_root_residual_maximum": "1e-24",
    "compact_point_stable_root_exact_dyadic_residual_maximum": "1e-24",
    "compact_point_stable_root_maximum_ULP_distance_from_rounded_exact_root": 8,
    "exact_zero_unit_seed_normalized_error_maximum": "3.5527136788005009e-15",
    "independent_nontrivial_control_count": 2,
    "independent_control_stable_root_residual_maximum": "1e-12",
    "independent_controls_must_use_exact_zero_unit_oracle": True,
    "captured_full_grid_stable_residual_must_pass_unchanged_strict_gate": True,
    "captured_full_grid_condition_must_pass_unchanged_limit": True,
    "captured_full_grid_raw_bundle_is_optional_for_clean_clone_reproduction": True,
    "captured_full_grid_raw_bundle_must_match_all_frozen_hashes_when_present": True,
    "captured_full_grid_partial_presence_must_fail_closed": True,
    "legacy_evaluator_must_reproduce_the_PREF11_false_rejection": True,
    "one_bit_input_mutation_must_change_the_stable_result": True,
    "invalid_shapes_factors_radii_and_condition_limits_must_fail_closed": True,
    "result_must_bind_config_controls_implementation_and_authority_hashes": True,
    "amplitude_three_direct_coarse_phi_derivative_veto_remains_binding": True,
    "no_successor_protocol_or_trajectory_is_authorized": True,
    "no_SGBL_or_FGCQR_outcome_may_be_read": True,
    "canonical_hash_bound_result_required": True,
}
EXPECTED_CAPTURE = {
    "point_count": 8192,
    "stable_complete_residual_infinity": "3.296668493746324e-14",
    "stable_kinetic_condition_infinity_maximum": "2859227.196975944",
    "maximum_residual_point_index": 872,
    "maximum_residual_row_index": 3,
    "maximum_residual_radius_binary64_hex": "0x1.b480000000000p+3",
    "stable_acceleration_array_sha256": "98360dd52ae44e8cf5e43861d0fe9e51cef2fa53e7ecf1324390574f3b4ee744",
    "stable_residual_array_sha256": "1def670bbdab9fdc518a7944bf319f03a4925d3111f8e341e9a53a3c4072c1fc",
    "stable_raw_gate_passed": True,
    "legacy_complete_residual_infinity": "1.0659145473163184e-12",
    "legacy_raw_gate_passed": False,
    "fresh_trajectory_result": False,
}
EXPECTED_SUCCESSOR = {
    "next_gate": "FGC-1-RSP1",
    "purpose": "prospectively_freeze_and_adjudicate_the_independent_amplitude_three_resolution_spectrum_veto",
    "SRC3_evaluator_may_be_consumed_only_after_its_canonical_certificate_passes": True,
    "PROTO12_may_be_considered_only_after_SRC3_and_the_independent_resolution_study_pass": True,
    "fresh_GR0_campaign_must_use_a_separately_frozen_protocol_and_namespace": True,
    "fresh_GR0_outcomes_must_not_be_tuned_after_inspection": True,
}
EXPECTED_CLAIMS = {
    "SRC3_reference_covariant_binary64_evaluator_derived": True,
    "SRC3_exact_reference_bitwise_preserved": True,
    "SRC3_compact_and_independent_exact_oracle_controls_passed": True,
    "SRC3_captured_full_grid_source_gate_passed": True,
    "SRC3_unredefined_equations_and_raw_threshold_preserved": True,
    "SRC3_runtime_source_instrument_authorized_for_a_separately_frozen_GR0_protocol": True,
    "amplitude_three_spectral_veto_cleared": False,
    "PROTO12_frozen": False,
    "fresh_GR0_calibration_completed": False,
    "classical_spherical_diagnostic_authorized": False,
    "SGBL_execution_authorized": False,
    "FGCQR_holdout_execution_authorized": False,
    "FGCQR_mechanism_rejected": False,
    "retained_EFT_evolution_authorized": False,
    "physical_transition_claim_authorized": False,
    "general_gradient_route_rejected": False,
    "singularity_resolution_derived": False,
    "child_domain_or_topology_derived": False,
    "dark_sector_mechanism_derived": False,
    "varying_locally_measured_c_derived": False,
}


def _rel(path: Path) -> str:
    return path.resolve().relative_to(REPOSITORY).as_posix()


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _canonical(value: Mapping[str, Any]) -> str:
    return json.dumps(
        value,
        sort_keys=True,
        indent=2,
        ensure_ascii=True,
        allow_nan=False,
    ) + "\n"


def _atomic_write(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        Path(temporary).replace(path)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise


def _reject_duplicate_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    answer: dict[str, Any] = {}
    for key, value in pairs:
        if key in answer:
            raise ValueError(f"duplicate JSON key: {key}")
        answer[key] = value
    return answer


def _load_canonical_json(path: Path, artifact_id: str | None = None) -> dict[str, Any]:
    raw = path.read_bytes()
    value = json.loads(raw, object_pairs_hook=_reject_duplicate_pairs)
    if not isinstance(value, dict) or raw != _canonical(value).encode("utf-8"):
        raise ValueError(f"{_rel(path)} must contain canonical JSON")
    if artifact_id is not None and value.get("artifact_id") != artifact_id:
        raise ValueError(f"{_rel(path)} artifact identifier differs")
    return value


def _strict_keys(name: str, value: Mapping[str, Any], expected: set[str]) -> None:
    if set(value) != expected:
        raise ValueError(f"{name} keys differ")


def _fraction_float(value: object) -> float:
    if not isinstance(value, str):
        raise TypeError("frozen numerical bounds must be rational strings")
    return float(Fraction(value))


def _hex_float(value: object) -> float:
    if not isinstance(value, str):
        raise TypeError("control values must be binary64 hex strings")
    answer = float.fromhex(value)
    if not np.isfinite(answer) or answer.hex() != value:
        raise ValueError("control value is not canonical finite binary64 hex")
    return answer


def _hex_row(value: object, *, name: str) -> np.ndarray:
    if not isinstance(value, list) or len(value) != 6:
        raise ValueError(f"{name} must contain six binary64 hex values")
    return np.asarray([_hex_float(item) for item in value], dtype=np.float64)


def _ordered_float_bits(value: float) -> int:
    signed = struct.unpack(">q", struct.pack(">d", value))[0]
    return 0x8000000000000000 - signed if signed < 0 else signed


def _ulp_distance(left: float, right: float) -> int:
    if left == right:
        return 0
    return abs(_ordered_float_bits(left) - _ordered_float_bits(right))


def _git(*arguments: str) -> subprocess.CompletedProcess[bytes]:
    return subprocess.run(
        ["git", *arguments],
        cwd=REPOSITORY,
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    with path.open("rb") as handle:
        config = tomllib.load(handle)
    _strict_keys(
        "SRC3 config",
        config,
        {
            "schema_version",
            "artifact_id",
            "project_version",
            "metric_signature",
            "riemann_convention",
            "scope",
            "lineage",
            "evaluator",
            "proof_contract",
            "captured_full_grid_evidence",
            "successor_boundary",
            "claims",
        },
    )
    if (
        config["schema_version"] != 1
        or config["artifact_id"] != ARTIFACT_ID
        or config["project_version"] != "0.11.0"
        or config["metric_signature"] != "-+++"
        or config["riemann_convention"] != "plus_partial_mu_gamma_nu"
        or config["scope"] != EXPECTED_SCOPE
        or config["evaluator"] != EXPECTED_EVALUATOR
        or config["proof_contract"] != EXPECTED_PROOF
        or config["captured_full_grid_evidence"] != EXPECTED_CAPTURE
        or config["successor_boundary"] != EXPECTED_SUCCESSOR
        or config["claims"] != EXPECTED_CLAIMS
    ):
        raise ValueError("SRC3 scope, evaluator, proof, successor, or claim contract differs")
    lineage = config["lineage"]
    _strict_keys(
        "SRC3 lineage",
        lineage,
        {
            "src3_parent_git_commit",
            "pref11_result",
            "pref11_result_sha256",
            "pref11_config",
            "pref11_config_sha256",
            "pref11_point_fixture",
            "pref11_point_fixture_sha256",
            "src3_control_fixture",
            "src3_control_fixture_sha256",
            "cap1_manifest",
            "cap1_manifest_sha256",
            "cap1_replay_result",
            "cap1_replay_result_sha256",
            "cap1_raw_fixture",
            "cap1_raw_fixture_sha256",
            "cap1_arithmetic_result",
            "cap1_arithmetic_result_sha256",
        },
    )
    if lineage["src3_parent_git_commit"] != "e519de824bfc4a8ad5265b6d8a209f4f791ff218":
        raise ValueError("SRC3 predecessor commit differs")
    return config


def _validate_tracked_lineage(config: Mapping[str, Any]) -> None:
    lineage = config["lineage"]
    for path_key, hash_key, artifact_id in (
        ("pref11_result", "pref11_result_sha256", "FGC-1-SRC2-PREF11"),
        ("pref11_config", "pref11_config_sha256", None),
        ("pref11_point_fixture", "pref11_point_fixture_sha256", "FGC-1-SRC2-PREF11-POINT0"),
        ("src3_control_fixture", "src3_control_fixture_sha256", "FGC-1-SRC3-CONTROLS"),
    ):
        path = REPOSITORY / lineage[path_key]
        if not path.is_file() or _sha(path) != lineage[hash_key]:
            raise ValueError(f"SRC3 tracked lineage hash differs: {lineage[path_key]}")
        if artifact_id is not None:
            _load_canonical_json(path, artifact_id)
    commit = lineage["src3_parent_git_commit"]
    if _git("cat-file", "-e", f"{commit}^{{commit}}").returncode != 0:
        raise ValueError("SRC3 predecessor commit is unavailable")
    if _git("merge-base", "--is-ancestor", commit, "HEAD").returncode != 0:
        raise ValueError("SRC3 predecessor commit is not an ancestor of HEAD")


def _load_point(config: Mapping[str, Any]) -> tuple[dict[str, Any], tuple[np.ndarray, ...], float]:
    path = REPOSITORY / config["lineage"]["pref11_point_fixture"]
    fixture = _load_canonical_json(path, "FGC-1-SRC2-PREF11-POINT0")
    if (
        fixture.get("field_order") != FIELD_ORDER
        or fixture.get("equation_order") != EQUATION_ORDER
        or fixture.get("point_index") != 0
        or fixture.get("coordinate_radius") != "0x1.0000000000000p-6"
    ):
        raise ValueError("SRC3 PREF11 compact point identity differs")
    lower = fixture.get("lower_jet")
    if not isinstance(lower, dict) or set(lower) != set(FIELD_NAMES):
        raise ValueError("SRC3 PREF11 lower-jet inventory differs")
    rows = tuple(_hex_row(lower[name], name=f"compact {name}") for name in FIELD_NAMES)
    return fixture, rows, _hex_float(fixture["coordinate_radius"])


def _load_controls(config: Mapping[str, Any]) -> dict[str, Any]:
    path = REPOSITORY / config["lineage"]["src3_control_fixture"]
    fixture = _load_canonical_json(path, "FGC-1-SRC3-CONTROLS")
    _strict_keys(
        "SRC3 controls",
        fixture,
        {
            "schema_version",
            "artifact_id",
            "exact_reference_radii_binary64_hex",
            "field_order",
            "independent_nontrivial_controls",
        },
    )
    controls = fixture["independent_nontrivial_controls"]
    radii = fixture["exact_reference_radii_binary64_hex"]
    if (
        fixture["schema_version"] != 1
        or fixture["field_order"] != FIELD_ORDER
        or not isinstance(radii, list)
        or len(radii) != 7
        or not isinstance(controls, list)
        or len(controls) != EXPECTED_PROOF["independent_nontrivial_control_count"]
    ):
        raise ValueError("SRC3 control fixture scope differs")
    if len({_hex_float(value) for value in radii}) != len(radii):
        raise ValueError("SRC3 exact-reference radii are not unique")
    identifiers: set[str] = set()
    for control in controls:
        if not isinstance(control, dict):
            raise ValueError("SRC3 nontrivial control must be an object")
        _strict_keys("SRC3 nontrivial control", control, {"control_id", "lower_jet", "radius_binary64_hex"})
        if not isinstance(control["control_id"], str) or control["control_id"] in identifiers:
            raise ValueError("SRC3 control identifier is invalid or repeated")
        identifiers.add(control["control_id"])
        lower = control["lower_jet"]
        if not isinstance(lower, dict) or set(lower) != set(FIELD_NAMES):
            raise ValueError("SRC3 nontrivial lower-jet inventory differs")
        for name in FIELD_NAMES:
            _hex_row(lower[name], name=f"{control['control_id']} {name}")
        _hex_float(control["radius_binary64_hex"])
    return fixture


def _exact_reference_record(fixture: Mapping[str, Any]) -> dict[str, Any]:
    radii = np.asarray(
        [_hex_float(value) for value in fixture["exact_reference_radii_binary64_hex"]],
        dtype=np.float64,
    )
    points = radii.size
    u = np.zeros((points, 6), dtype=np.float64)
    p = np.zeros_like(u)
    q = np.zeros_like(u)
    p_r = np.zeros_like(u)
    q_r = np.zeros_like(u)
    u[:, 0] = 1.0
    u[:, 2] = 1.0
    u[:, 3] = radii
    q[:, 3] = 1.0
    evaluated = gr0_reference_balanced_ref1_residual_batch(
        u,
        p,
        q,
        np.zeros((1, points, 6), dtype=np.float64),
        p_r,
        q_r,
        radii,
    )
    zero_surfaces = {
        "full_residual": evaluated.full_residual,
        "unredefined_metric_residual": evaluated.unredefined_metric_residual,
        "scalar_residual": evaluated.scalar_residual,
        "gauge_constraint": evaluated.gauge_constraint,
        "hamiltonian_constraint": evaluated.hamiltonian_constraint,
        "momentum_constraint": evaluated.momentum_constraint,
        "ricci_scalar": evaluated.ricci_scalar,
        "ricci_squared": evaluated.ricci_squared,
    }
    if not all(np.array_equal(value, np.zeros_like(value)) for value in zero_surfaces.values()):
        raise ValueError("SRC3 exact reference is not bitwise stationary")
    solved = diagnose_gr0_reference_balanced_accelerations(
        u, p, q, p_r, q_r, radii
    )
    if (
        not solved.raw_gate_passed
        or solved.residual_infinity != 0.0
        or not np.array_equal(solved.accelerations, np.zeros_like(solved.accelerations))
    ):
        raise ValueError("SRC3 exact-reference affine solve is not bitwise zero")
    return {
        "radii_binary64_hex": list(fixture["exact_reference_radii_binary64_hex"]),
        "all_complete_residual_constraint_and_curvature_surfaces_bitwise_zero": True,
        "affine_acceleration_bitwise_zero": True,
        "verified_residual_infinity": 0.0,
        "kinetic_condition_infinity_maximum": solved.kinetic_condition_infinity_maximum,
    }


def _point_record(rows: tuple[np.ndarray, ...], radius: float) -> dict[str, Any]:
    constant, matrix = exact_affine_system(*rows, radius)
    exact_root = exact_gaussian_solve(matrix, tuple(-value for value in constant))
    rounded_exact = np.asarray([float(value) for value in exact_root], dtype=np.float64)
    solved = diagnose_gr0_reference_balanced_accelerations(
        *(value[None, :] for value in rows), np.asarray((radius,), dtype=np.float64)
    )
    stable_exact = exact_ref1_residual(*rows, radius, solved.accelerations[0])
    exact_infinity = float(max((abs(value) for value in stable_exact), default=Fraction(0)))
    ulps = [
        _ulp_distance(left, right)
        for left, right in zip(solved.accelerations[0], rounded_exact, strict=True)
    ]
    legacy = gr0_ref1_residual_batch(
        *(value[None, :] for value in rows[:3]),
        solved.accelerations[None, :, :],
        *(value[None, :] for value in rows[3:]),
        np.asarray((radius,), dtype=np.float64),
    ).full_residual[0, 0]
    legacy_infinity = float(np.max(np.abs(legacy), initial=0.0))
    proof = EXPECTED_PROOF
    if (
        not solved.raw_gate_passed
        or solved.residual_infinity >= _fraction_float(proof["compact_point_stable_root_residual_maximum"])
        or exact_infinity >= _fraction_float(proof["compact_point_stable_root_exact_dyadic_residual_maximum"])
        or max(ulps) > proof["compact_point_stable_root_maximum_ULP_distance_from_rounded_exact_root"]
        or legacy_infinity <= _fraction_float(EXPECTED_EVALUATOR["raw_complete_residual_maximum"])
    ):
        raise ValueError("SRC3 compact-point proof burden failed")
    return {
        "radius_binary64_hex": radius.hex(),
        "stable_root_binary64_hex": [value.hex() for value in solved.accelerations[0]],
        "rounded_exact_root_binary64_hex": [value.hex() for value in rounded_exact],
        "stable_root_component_ULP_distance_from_rounded_exact_root": ulps,
        "stable_binary64_residual_infinity": solved.residual_infinity,
        "stable_root_exact_dyadic_residual_infinity": exact_infinity,
        "legacy_binary64_residual_at_stable_root_infinity": legacy_infinity,
        "unchanged_strict_raw_gate_passed": solved.raw_gate_passed,
        "legacy_false_rejection_preserved": True,
    }


def _nontrivial_control_record(control: Mapping[str, Any]) -> dict[str, Any]:
    lower = tuple(
        _hex_row(control["lower_jet"][name], name=f"{control['control_id']} {name}")
        for name in FIELD_NAMES
    )
    radius = _hex_float(control["radius_binary64_hex"])
    seeds = np.zeros((7, 1, 6), dtype=np.float64)
    for field in range(6):
        seeds[field + 1, 0, field] = 1.0
    stable_seed_residuals = gr0_reference_balanced_ref1_residual_batch(
        *(value[None, :] for value in lower[:3]),
        seeds,
        *(value[None, :] for value in lower[3:]),
        np.asarray((radius,), dtype=np.float64),
    ).full_residual[:, 0]
    exact_seed_residuals: list[tuple[Fraction, ...]] = [
        exact_ref1_residual(*lower, radius, acceleration)
        for acceleration in seeds[:, 0]
    ]
    exact_seed_float = np.asarray(
        [[float(value) for value in row] for row in exact_seed_residuals],
        dtype=np.float64,
    )
    absolute_error = float(
        np.max(np.abs(stable_seed_residuals - exact_seed_float), initial=0.0)
    )
    exact_scale = max(1.0, float(np.max(np.abs(exact_seed_float), initial=0.0)))
    normalized_error = absolute_error / exact_scale
    constant = exact_seed_residuals[0]
    matrix = tuple(
        tuple(
            exact_seed_residuals[column + 1][row] - constant[row]
            for column in range(6)
        )
        for row in range(6)
    )
    exact_root = exact_gaussian_solve(matrix, tuple(-value for value in constant))
    solved = diagnose_gr0_reference_balanced_accelerations(
        *(value[None, :] for value in lower), np.asarray((radius,), dtype=np.float64)
    )
    stable_exact = exact_ref1_residual(*lower, radius, solved.accelerations[0])
    stable_exact_infinity = float(
        max((abs(value) for value in stable_exact), default=Fraction(0))
    )
    if (
        normalized_error > _fraction_float(EXPECTED_PROOF["exact_zero_unit_seed_normalized_error_maximum"])
        or not solved.raw_gate_passed
        or solved.residual_infinity >= _fraction_float(EXPECTED_PROOF["independent_control_stable_root_residual_maximum"])
    ):
        raise ValueError(f"SRC3 independent control failed: {control['control_id']}")
    return {
        "control_id": control["control_id"],
        "radius_binary64_hex": radius.hex(),
        "exact_zero_unit_seed_scale": exact_scale,
        "stable_vs_exact_seed_absolute_error_infinity": absolute_error,
        "stable_vs_exact_seed_normalized_error_infinity": normalized_error,
        "rounded_exact_root_binary64_hex": [float(value).hex() for value in exact_root],
        "stable_root_binary64_hex": [value.hex() for value in solved.accelerations[0]],
        "stable_binary64_residual_infinity": solved.residual_infinity,
        "stable_root_exact_dyadic_residual_infinity": stable_exact_infinity,
        "unchanged_strict_raw_gate_passed": solved.raw_gate_passed,
    }


def _mutation_and_failure_record(rows: tuple[np.ndarray, ...], radius: float) -> dict[str, Any]:
    baseline = diagnose_gr0_reference_balanced_accelerations(
        *(value[None, :] for value in rows), np.asarray((radius,))
    )
    mutated = list(rows)
    mutated_u = mutated[0].copy()
    mutated_u[0] = np.nextafter(mutated_u[0], np.inf)
    mutated[0] = mutated_u
    changed = diagnose_gr0_reference_balanced_accelerations(
        *(value[None, :] for value in mutated), np.asarray((radius,))
    )
    one_bit_visible = not np.array_equal(changed.accelerations, baseline.accelerations)

    failures: dict[str, bool] = {}
    try:
        diagnose_gr0_reference_balanced_accelerations(
            *(value[None, :] for value in rows), np.asarray((0.0,))
        )
    except ValueError:
        failures["nonpositive_radius"] = True
    try:
        diagnose_gr0_reference_balanced_accelerations(
            *(value[None, :] for value in rows),
            np.asarray((radius,)),
            condition_number_maximum=1.0,
        )
    except ValueError:
        failures["condition_limit"] = True
    try:
        gr0_reference_balanced_ref1_residual_batch(
            *(value[None, :] for value in rows[:3]),
            np.zeros((1, 1, 6), dtype=np.float64),
            *(value[None, :] for value in rows[3:]),
            np.asarray((radius,)),
            tilde_normal_factor=4.0,
            hat_normal_factor=4.0,
        )
    except ValueError:
        failures["invalid_auxiliary_factors"] = True
    if not one_bit_visible or failures != {
        "nonpositive_radius": True,
        "condition_limit": True,
        "invalid_auxiliary_factors": True,
    }:
        raise ValueError("SRC3 mutation or typed-failure controls failed")
    return {
        "one_bit_alpha_mutation_visible": one_bit_visible,
        "baseline_acceleration_array_sha256": array_content_sha256(baseline.accelerations),
        "mutated_acceleration_array_sha256": array_content_sha256(changed.accelerations),
        "typed_failure_controls": failures,
    }


def _validate_raw_bundle_if_present(
    config: Mapping[str, Any],
    point_rows: tuple[np.ndarray, ...],
    point_radius: float,
) -> None:
    lineage = config["lineage"]
    bindings = (
        ("cap1_manifest", "cap1_manifest_sha256"),
        ("cap1_replay_result", "cap1_replay_result_sha256"),
        ("cap1_raw_fixture", "cap1_raw_fixture_sha256"),
        ("cap1_arithmetic_result", "cap1_arithmetic_result_sha256"),
    )
    paths = [REPOSITORY / lineage[path_key] for path_key, _hash_key in bindings]
    present = [path.is_file() for path in paths]
    if any(present) and not all(present):
        raise ValueError("CAP1 raw bundle is only partially present")
    if not all(present):
        return
    for path, (_path_key, hash_key) in zip(paths, bindings, strict=True):
        if _sha(path) != lineage[hash_key]:
            raise ValueError(f"CAP1 raw bundle hash differs: {_rel(path)}")
    manifest = _load_canonical_json(paths[0])
    replay = _load_canonical_json(paths[1])
    arithmetic = _load_canonical_json(paths[3])
    if (
        manifest.get("FGCQR_outcome_read") is not False
        or replay.get("complete_rejection_trace_reproduced") is not True
        or replay.get("source_only_rejection_count") != 164
        or replay.get("FGCQR_outcome_read") is not False
        or arithmetic.get("best_candidate_strict_raw_gate_passed") is not False
        or arithmetic.get("PROTO12_frozen") is not False
        or arithmetic.get("FGCQR_holdout_execution_authorized") is not False
    ):
        raise ValueError("CAP1 raw bundle scope or outcome differs")
    with np.load(paths[2], allow_pickle=False) as raw:
        required = {
            "metadata_utf8",
            "u",
            "p",
            "q",
            "p_r",
            "q_r",
            "radii",
            "forward_constant",
            "forward_jacobian",
            "baseline_acceleration",
            "baseline_complete_residual",
            "best_acceleration",
            "best_complete_residual",
        }
        if set(raw.files) != required:
            raise ValueError("CAP1 raw fixture inventory differs")
        if raw["radii"].shape != (EXPECTED_CAPTURE["point_count"],):
            raise ValueError("CAP1 raw grid size differs")
        for name, expected in zip(FIELD_NAMES, point_rows, strict=True):
            if not np.array_equal(raw[name][0], expected):
                raise ValueError(f"SRC3 compact point differs from raw {name}")
        if raw["radii"][0] != point_radius:
            raise ValueError("SRC3 compact radius differs from raw fixture")
        legacy_infinity = float(
            np.max(np.abs(raw["baseline_complete_residual"]), initial=0.0)
        )
        solved = diagnose_gr0_reference_balanced_accelerations(
            *(raw[name] for name in (*FIELD_NAMES, "radii")),
            raw_tolerance=_fraction_float(EXPECTED_EVALUATOR["raw_complete_residual_maximum"]),
            condition_number_maximum=_fraction_float(EXPECTED_EVALUATOR["kinetic_condition_number_maximum"]),
            maximum_refinement_iterations=EXPECTED_EVALUATOR["maximum_complete_refinement_iterations"],
        )
    last = solved.iterations[-1]
    if (
        not solved.raw_gate_passed
        or solved.residual_infinity != _fraction_float(EXPECTED_CAPTURE["stable_complete_residual_infinity"])
        or solved.kinetic_condition_infinity_maximum
        != _fraction_float(EXPECTED_CAPTURE["stable_kinetic_condition_infinity_maximum"])
        or last.maximum_residual_point_index != EXPECTED_CAPTURE["maximum_residual_point_index"]
        or last.maximum_residual_row_index != EXPECTED_CAPTURE["maximum_residual_row_index"]
        or last.maximum_residual_radius.hex() != EXPECTED_CAPTURE["maximum_residual_radius_binary64_hex"]
        or array_content_sha256(solved.accelerations)
        != EXPECTED_CAPTURE["stable_acceleration_array_sha256"]
        or array_content_sha256(solved.residuals)
        != EXPECTED_CAPTURE["stable_residual_array_sha256"]
        or legacy_infinity != _fraction_float(EXPECTED_CAPTURE["legacy_complete_residual_infinity"])
        or legacy_infinity < _fraction_float(EXPECTED_EVALUATOR["raw_complete_residual_maximum"])
    ):
        raise ValueError("SRC3 captured full-grid evidence differs")


def reproduce(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    config = load_config(config_path)
    _validate_tracked_lineage(config)
    point_fixture, point_rows, point_radius = _load_point(config)
    controls = _load_controls(config)
    reference = _exact_reference_record(controls)
    point = _point_record(point_rows, point_radius)
    independent = [
        _nontrivial_control_record(control)
        for control in controls["independent_nontrivial_controls"]
    ]
    attacks = _mutation_and_failure_record(point_rows, point_radius)
    _validate_raw_bundle_if_present(config, point_rows, point_radius)
    lineage = config["lineage"]
    control_path = REPOSITORY / lineage["src3_control_fixture"]
    point_path = REPOSITORY / lineage["pref11_point_fixture"]
    predecessor_path = REPOSITORY / lineage["pref11_result"]
    capture = dict(config["captured_full_grid_evidence"])
    capture.update(
        {
            "stable_complete_residual_infinity": _fraction_float(
                capture["stable_complete_residual_infinity"]
            ),
            "stable_kinetic_condition_infinity_maximum": _fraction_float(
                capture["stable_kinetic_condition_infinity_maximum"]
            ),
            "legacy_complete_residual_infinity": _fraction_float(
                capture["legacy_complete_residual_infinity"]
            ),
            "raw_bundle_optional_for_clean_clone_reproduction": True,
            "partial_or_hash_mismatched_raw_bundle_fails_closed": True,
        }
    )
    return {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": config["project_version"],
        "classification": "prospectively_frozen_reference_covariant_binary64_GR0_source_evaluator",
        "generated_by": _rel(Path(__file__)),
        "derivation_document": _rel(OWNER_DOCUMENT),
        "derivation_document_sha256": _sha(OWNER_DOCUMENT),
        "source_config_sha256": {
            _rel(config_path): _sha(config_path),
            _rel(control_path): _sha(control_path),
            _rel(point_path): _sha(point_path),
        },
        "predecessor_sha256": {_rel(predecessor_path): _sha(predecessor_path)},
        "implementation_sha256": {_rel(path): _sha(path) for path in IMPLEMENTATION},
        "scope_bindings": dict(config["scope"]),
        "gate_status": dict(EXPECTED_CLAIMS),
        "artifact_payload": {
            "reference_covariant_identity": {
                "metric_perturbation": config["evaluator"]["perturbation"],
                "connection_difference": config["evaluator"]["identity"],
                "derivative_order": config["evaluator"]["connection_difference_derivative"],
                "ricci_assembly": config["evaluator"]["ricci_assembly"],
                "continuum_equations": config["evaluator"]["continuum_equations"],
                "continuum_equations_changed": False,
                "reference_geometry_changed": False,
                "raw_residual_threshold_changed": False,
            },
            "exact_reference_control": reference,
            "pref11_compact_point_control": point,
            "independent_nontrivial_exact_controls": independent,
            "mutation_and_typed_failure_controls": attacks,
            "captured_full_grid_development_evidence": capture,
            "successor_boundary": dict(config["successor_boundary"]),
            "epistemic_boundary": {
                "established": [
                    "the unchanged GR-0 REF1 equations admit a reference-covariant binary64 evaluation order that is bitwise exact on spherical Minkowski",
                    "the SRC3 evaluator converges to the exact dyadic oracle on the PREF11 point and two independent nontrivial controls",
                    "the SRC3 evaluator clears the unchanged strict raw source gate across the captured 8192-point CAP1 call",
                    "the specific CAP1 source wall was numerical evaluation-order cancellation rather than a demonstrated equation obstruction",
                    "SRC3 is fit only for consumption by a separately frozen prospective GR-0 protocol",
                ],
                "not_established": [
                    "resolution convergence or removal of amplitude 3's independent coarse derivative-tail veto",
                    "a completed or eligible fresh GR-0 calibration",
                    "SGB-L or FGC-QR dynamics, activation, collapse, or defocusing",
                    "retained-EFT validity or a physical transition",
                    "a singularity resolution, child domain, dark sector, or variable local light speed",
                ],
            },
        },
        "nonclaims": {
            key: value for key, value in EXPECTED_CLAIMS.items() if value is False
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    record = reproduce(args.config.resolve())
    payload = _canonical(record).encode("utf-8")
    output = args.output.resolve()
    if args.check:
        if not output.is_file() or output.read_bytes() != payload:
            raise SystemExit(f"stored result differs: {output}")
        print(f"verified {output}")
        return
    _atomic_write(output, payload)
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
