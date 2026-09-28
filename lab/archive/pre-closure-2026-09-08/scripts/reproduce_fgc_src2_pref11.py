#!/usr/bin/env python3
"""Reproduce PREF11's exact diagnosis of the captured SRC2 source wall."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import tomllib
from typing import Any, Mapping

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
import sys

sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    array_content_sha256,
)
from recursive_horizons.fgc.evolution.src2_affine_arithmetic import (  # noqa: E402
    assemble_forward_affine_system,
)
from recursive_horizons.fgc.evolution.src2_exact_oracle import (  # noqa: E402
    diagnose_exact_oracle,
    exact_affine_system,
)


ARTIFACT_ID = "FGC-1-SRC2-PREF11"
FIXTURE_ARTIFACT_ID = "FGC-1-SRC2-PREF11-POINT0"
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-src2-pref11.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-src2-pref11.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-src2-pref11.md"
ORACLE_MODULE = (
    REPOSITORY
    / "src/recursive_horizons/fgc/evolution/src2_exact_oracle.py"
)
IMPLEMENTATION = (
    Path(__file__).resolve(),
    ORACLE_MODULE,
    REPOSITORY / "src/recursive_horizons/fgc/spherical_reduction.py",
    REPOSITORY / "src/recursive_horizons/fgc/modified_harmonic_reference.py",
    REPOSITORY / "src/recursive_horizons/fgc/reference_connection.py",
    REPOSITORY / "src/recursive_horizons/fgc/evolution/gr0_direct_source.py",
    REPOSITORY / "src/recursive_horizons/fgc/evolution/vectorized_source.py",
)
Q = Fraction


EXPECTED_CLAIMS = {
    "PROTO11_source_wall_replayed": True,
    "CAP1_frozen_arithmetic_routes_cleared_the_wall": False,
    "captured_GR0_point0_exact_affine_system_nonsingular": True,
    "captured_GR0_point0_exact_root_exists": True,
    "captured_GR0_point0_binary64_residual_evaluator_incomplete": True,
    "captured_GR0_point0_wall_is_a_floating_cancellation_instrument_failure": True,
    "captured_GR0_point0_physical_or_structural_breakdown_demonstrated": False,
    "cancellation_resistant_source_evaluator_required": True,
    "amplitude_three_spectral_veto_cleared": False,
    "PROTO12_frozen": False,
    "GR0_calibration_completed": False,
    "classical_spherical_diagnostic_authorized": False,
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


def _bytes_sha(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


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


def _fraction_float(value: str) -> float:
    return float(Q(value))


def _hex_float(value: object) -> float:
    if not isinstance(value, str):
        raise TypeError("compact point values must be binary64 hex strings")
    answer = float.fromhex(value)
    if not np.isfinite(answer) or answer.hex() != value:
        raise ValueError("compact point value is not canonical finite binary64 hex")
    return answer


def _hex_row(value: object, *, name: str) -> np.ndarray:
    if not isinstance(value, list) or len(value) != 6:
        raise ValueError(f"{name} must contain six binary64 hex values")
    return np.asarray([_hex_float(item) for item in value], dtype=np.float64)


def _git(*args: str, check: bool = True) -> subprocess.CompletedProcess[bytes]:
    return subprocess.run(
        ["git", *args],
        cwd=REPOSITORY,
        check=check,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    with path.open("rb") as handle:
        config = tomllib.load(handle)
    _strict_keys(
        "PREF11 config",
        config,
        {
            "schema_version",
            "artifact_id",
            "project_version",
            "metric_signature",
            "riemann_convention",
            "scope",
            "lineage",
            "cap1_outcome",
            "exact_oracle",
            "proof_contract",
            "successor_boundary",
            "claims",
        },
    )
    scope = config["scope"]
    lineage = config["lineage"]
    cap1 = config["cap1_outcome"]
    oracle = config["exact_oracle"]
    proof = config["proof_contract"]
    successor = config["successor_boundary"]
    _strict_keys(
        "PREF11 scope",
        scope,
        {
            "role",
            "branch",
            "amplitude",
            "member",
            "captured_point_index",
            "captured_coordinate_radius",
            "GR0_trajectory_read",
            "SGBL_trajectory_read",
            "FGCQR_trajectory_read",
            "collapse_or_trapped_outcome_classified",
            "mechanism_question_answered",
        },
    )
    _strict_keys(
        "PREF11 lineage",
        lineage,
        {
            "cal8_result",
            "cal8_result_sha256",
            "cap1_config",
            "cap1_config_sha256",
            "cap1_git_commit",
            "cap1_manifest",
            "cap1_manifest_sha256",
            "cap1_replay_result",
            "cap1_replay_result_sha256",
            "cap1_raw_fixture",
            "cap1_raw_fixture_sha256",
            "cap1_arithmetic_result",
            "cap1_arithmetic_result_sha256",
            "compact_point_fixture",
            "compact_point_fixture_sha256",
        },
    )
    _strict_keys(
        "PREF11 CAP1 outcome",
        cap1,
        {
            "complete_rejection_trace_reproduced",
            "source_only_rejection_count",
            "target_source_call_matches",
            "baseline_complete_residual_infinity",
            "raw_complete_residual_maximum",
            "best_frozen_route",
            "best_frozen_route_passed",
            "power_two_equilibration_passed",
            "symmetric_affine_extraction_passed",
            "captured_coefficient_exact_binary_refinement_passed",
            "maximum_captured_kinetic_condition_infinity",
        },
    )
    _strict_keys(
        "PREF11 exact oracle",
        oracle,
        {
            "binary64_literals_are_reconstructed_as_exact_dyadic_rationals",
            "continuum_equations",
            "affine_extraction",
            "linear_solve",
            "root_rounding",
            "reference_connection",
            "tilde_normal_factor",
            "hat_normal_factor",
            "exact_root_residual_required_to_be_zero",
            "exact_rounded_root_residual_maximum",
            "floating_evaluator_false_rejection_minimum",
            "raw_threshold_change_forbidden",
            "continuum_equation_change_forbidden",
            "reference_derivative_map_change_forbidden",
        },
    )
    _strict_keys(
        "PREF11 proof contract",
        proof,
        {
            "compact_fixture_must_be_canonical_JSON",
            "compact_fixture_must_match_the_raw_fixture_point_when_raw_bundle_is_present",
            "raw_bundle_is_optional_for_clean_clone_reproduction",
            "raw_bundle_must_match_all_frozen_hashes_when_present",
            "raw_bundle_partial_presence_must_fail_closed",
            "CAP1_manifest_must_bind_the_immutable_git_commit_and_config",
            "CAP1_replay_must_bind_the_raw_fixture_and_reproduce_the_complete_trace",
            "CAP1_arithmetic_result_must_report_no_frozen_route_pass",
            "exact_affine_matrix_must_be_nonsingular",
            "exact_root_must_zero_all_six_exact_rows",
            "binary64_rounded_exact_root_must_pass_the_exact_rational_residual_bound",
            "direct_and_generic_binary64_evaluators_must_both_false_reject_the_rounded_exact_root",
            "captured_baseline_and_forward_coefficients_must_match_the_compact_fixture_bits",
            "result_must_bind_config_fixture_implementation_and_authority_hashes",
            "amplitude_three_direct_coarse_phi_derivative_veto_remains_binding",
            "no_successor_protocol_or_trajectory_is_authorized",
            "no_SGBL_or_FGCQR_outcome_may_be_read",
            "canonical_hash_bound_result_required",
        },
    )
    _strict_keys(
        "PREF11 successor boundary",
        successor,
        {
            "next_gate",
            "purpose",
            "must_match_exact_oracle_on_compact_point_and_independent_nontrivial_controls",
            "must_preserve_unredefined_equations_reference_map_thresholds_and_protocol_inputs",
            "must_be_frozen_before_any_new_GR0_campaign",
            "PROTO12_may_be_considered_only_after_SRC3_passes",
            "amplitude_three_resolution_spectrum_veto_requires_an_independent_prospective_resolution_study",
        },
    )
    if (
        config["schema_version"] != 1
        or config["artifact_id"] != ARTIFACT_ID
        or config["metric_signature"] != "-+++"
        or config["riemann_convention"] != "plus_partial_mu_gamma_nu"
        or scope.get("role")
        != "hash_bound_exact_rational_diagnosis_of_the_captured_PROTO11_GR0_point_zero_source_wall"
        or scope.get("branch") != "GR-0"
        or scope.get("amplitude") != "5/2"
        or scope.get("member") != "RK4-8193"
        or scope.get("captured_point_index") != 0
        or scope.get("captured_coordinate_radius") != "1/64"
        or scope.get("GR0_trajectory_read") is not True
        or scope.get("SGBL_trajectory_read") is not False
        or scope.get("FGCQR_trajectory_read") is not False
        or scope.get("collapse_or_trapped_outcome_classified") is not False
        or scope.get("mechanism_question_answered") is not False
        or lineage.get("cap1_git_commit")
        != "241e9b99a44985202ad35bdfaf9e6e1cf2d2c3d6"
        or cap1.get("complete_rejection_trace_reproduced") is not True
        or cap1.get("source_only_rejection_count") != 164
        or cap1.get("target_source_call_matches") != 2
        or cap1.get("best_frozen_route") != "PROTO11_binary64_baseline"
        or cap1.get("best_frozen_route_passed") is not False
        or cap1.get("power_two_equilibration_passed") is not False
        or cap1.get("symmetric_affine_extraction_passed") is not False
        or cap1.get("captured_coefficient_exact_binary_refinement_passed") is not False
        or cap1.get("baseline_complete_residual_infinity")
        != "1.0659145473163184e-12"
        or cap1.get("raw_complete_residual_maximum") != "1e-12"
        or cap1.get("maximum_captured_kinetic_condition_infinity")
        != "2859227.196975944"
        or oracle.get("binary64_literals_are_reconstructed_as_exact_dyadic_rationals")
        is not True
        or oracle.get("continuum_equations") != "unredefined_ACT1_VAR1_GR0_REF1"
        or oracle.get("affine_extraction")
        != "exact_zero_plus_six_exact_unit_accelerations"
        or oracle.get("linear_solve")
        != "exact_fraction_gaussian_elimination"
        or oracle.get("root_rounding")
        != "correctly_rounded_binary64_per_component"
        or oracle.get("reference_connection") != "flat_spherical_annulus"
        or oracle.get("tilde_normal_factor") != 4
        or oracle.get("hat_normal_factor") != 9
        or oracle.get("exact_root_residual_required_to_be_zero") is not True
        or oracle.get("exact_rounded_root_residual_maximum") != "1e-26"
        or oracle.get("floating_evaluator_false_rejection_minimum") != "1e-12"
        or oracle.get("raw_threshold_change_forbidden") is not True
        or oracle.get("continuum_equation_change_forbidden") is not True
        or oracle.get("reference_derivative_map_change_forbidden") is not True
        or any(value is not True for value in proof.values())
        or successor.get("next_gate") != "FGC-1-SRC3"
        or successor.get("purpose")
        != "derive_and_predeclare_a_cancellation_resistant_runtime_REF1_source_evaluator"
        or any(
            successor.get(key) is not True
            for key in (
                "must_match_exact_oracle_on_compact_point_and_independent_nontrivial_controls",
                "must_preserve_unredefined_equations_reference_map_thresholds_and_protocol_inputs",
                "must_be_frozen_before_any_new_GR0_campaign",
                "amplitude_three_resolution_spectrum_veto_requires_an_independent_prospective_resolution_study",
            )
        )
        or successor.get("PROTO12_may_be_considered_only_after_SRC3_passes")
        is not True
        or config["claims"] != EXPECTED_CLAIMS
    ):
        raise ValueError("PREF11 scope, proof, successor, or claim contract differs")
    for path_key, hash_key in (
        ("cal8_result", "cal8_result_sha256"),
        ("cap1_config", "cap1_config_sha256"),
        ("compact_point_fixture", "compact_point_fixture_sha256"),
    ):
        target = REPOSITORY / lineage[path_key]
        if not target.is_file() or _sha(target) != lineage[hash_key]:
            raise ValueError(f"PREF11 tracked lineage differs: {lineage[path_key]}")
    commit = lineage["cap1_git_commit"]
    _git("cat-file", "-e", f"{commit}^{{commit}}")
    if _git("merge-base", "--is-ancestor", commit, "HEAD", check=False).returncode:
        raise ValueError("CAP1 commit is not an ancestor of HEAD")
    cap1_blob = _git("show", f"{commit}:{lineage['cap1_config']}").stdout
    if _bytes_sha(cap1_blob) != lineage["cap1_config_sha256"]:
        raise ValueError("CAP1 config differs at its immutable commit")
    return config


def _validate_compact_fixture(
    config: Mapping[str, Any],
) -> tuple[dict[str, Any], dict[str, np.ndarray | float]]:
    lineage = config["lineage"]
    path = REPOSITORY / lineage["compact_point_fixture"]
    fixture = _load_canonical_json(path, FIXTURE_ARTIFACT_ID)
    _strict_keys(
        "PREF11 compact fixture",
        fixture,
        {
            "schema_version",
            "artifact_id",
            "action_parameters",
            "baseline_acceleration",
            "baseline_complete_residual",
            "coordinate_radius",
            "equation_order",
            "field_order",
            "forward_constant",
            "forward_jacobian",
            "hat_normal_factor",
            "lower_jet",
            "point_index",
            "raw_fixture_array_sha256",
            "raw_fixture_sha256",
            "source_call",
            "tilde_normal_factor",
        },
    )
    expected_fields = ["alpha", "shift", "lambda", "R", "phi", "chi"]
    expected_equations = [
        "metric_tt_mhg",
        "metric_tr_mhg",
        "metric_rr_mhg",
        "metric_theta_theta_mhg",
        "scalar_phi",
        "scalar_chi",
    ]
    if (
        fixture["schema_version"] != 1
        or fixture["action_parameters"]
        != {"g4": "1/2", "model_id": "GR-0", "mu": "3", "planck_mass": "2"}
        or fixture["field_order"] != expected_fields
        or fixture["equation_order"] != expected_equations
        or fixture["point_index"] != 0
        or fixture["raw_fixture_sha256"] != lineage["cap1_raw_fixture_sha256"]
        or fixture["tilde_normal_factor"] != 4
        or fixture["hat_normal_factor"] != 9
        or fixture["source_call"]
        != {
            "amplitude": "5/2",
            "branch": "GR-0",
            "member": "RK4-8193",
            "stage_time": "0x1.21526722367a9p-4",
        }
    ):
        raise ValueError("PREF11 compact fixture identity or scope differs")
    lower = fixture["lower_jet"]
    _strict_keys("PREF11 compact lower jet", lower, {"u", "p", "q", "p_r", "q_r"})
    arrays: dict[str, np.ndarray | float] = {
        name: _hex_row(lower[name], name=name)
        for name in ("u", "p", "q", "p_r", "q_r")
    }
    arrays["radius"] = _hex_float(fixture["coordinate_radius"])
    arrays["baseline_acceleration"] = _hex_row(
        fixture["baseline_acceleration"], name="baseline_acceleration"
    )
    arrays["baseline_complete_residual"] = _hex_row(
        fixture["baseline_complete_residual"], name="baseline_complete_residual"
    )
    arrays["forward_constant"] = _hex_row(
        fixture["forward_constant"], name="forward_constant"
    )
    jacobian = fixture["forward_jacobian"]
    if not isinstance(jacobian, list) or len(jacobian) != 6:
        raise ValueError("forward_jacobian must contain six rows")
    arrays["forward_jacobian"] = np.stack(
        [_hex_row(row, name="forward_jacobian row") for row in jacobian]
    )
    if arrays["radius"] != 1.0 / 64.0:
        raise ValueError("PREF11 compact coordinate radius differs")
    return fixture, arrays


def _validate_raw_bundle_if_present(
    config: Mapping[str, Any],
    fixture: Mapping[str, Any],
    compact: Mapping[str, np.ndarray | float],
) -> None:
    lineage = config["lineage"]
    bindings = (
        ("cap1_manifest", "cap1_manifest_sha256"),
        ("cap1_replay_result", "cap1_replay_result_sha256"),
        ("cap1_raw_fixture", "cap1_raw_fixture_sha256"),
        ("cap1_arithmetic_result", "cap1_arithmetic_result_sha256"),
    )
    paths = [REPOSITORY / lineage[path_key] for path_key, _hash_key in bindings]
    present = [path.exists() for path in paths]
    if any(present) and not all(present):
        raise ValueError("CAP1 raw bundle is only partially present")
    if not all(present):
        return
    for path, (_path_key, hash_key) in zip(paths, bindings, strict=True):
        if not path.is_file() or _sha(path) != lineage[hash_key]:
            raise ValueError(f"CAP1 raw bundle hash differs: {_rel(path)}")
    manifest = _load_canonical_json(paths[0])
    replay = _load_canonical_json(paths[1])
    arithmetic = _load_canonical_json(paths[3])
    if (
        manifest.get("git_commit") != lineage["cap1_git_commit"]
        or manifest.get("config_sha256") != lineage["cap1_config_sha256"]
        or manifest.get("FGCQR_outcome_read") is not False
        or replay.get("fixture_sha256") != lineage["cap1_raw_fixture_sha256"]
        or replay.get("complete_rejection_trace_reproduced") is not True
        or replay.get("source_only_rejection_count") != 164
        or replay.get("target_source_call_matches") != 2
        or replay.get("FGCQR_outcome_read") is not False
        or replay.get("mechanism_question_answered") is not False
        or arithmetic.get("fixture_sha256") != lineage["cap1_raw_fixture_sha256"]
        or arithmetic.get("best_candidate_label") != "PROTO11_binary64_baseline"
        or arithmetic.get("best_candidate_strict_raw_gate_passed") is not False
        or arithmetic.get("PROTO12_frozen") is not False
        or arithmetic.get("FGCQR_holdout_execution_authorized") is not False
        or arithmetic.get("mechanism_question_answered") is not False
    ):
        raise ValueError("CAP1 raw replay or arithmetic outcome differs")
    with np.load(paths[2], allow_pickle=False) as raw:
        expected_keys = {
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
        if set(raw.files) != expected_keys:
            raise ValueError("CAP1 raw fixture array inventory differs")
        metadata_bytes = bytes(raw["metadata_utf8"])
        metadata = json.loads(metadata_bytes, object_pairs_hook=_reject_duplicate_pairs)
        if metadata_bytes != _canonical(metadata).encode("utf-8"):
            raise ValueError("CAP1 raw fixture metadata is not canonical JSON")
        metadata_hashes = metadata.get("array_sha256", {})
        if (
            not isinstance(metadata_hashes, dict)
            or any(
                metadata_hashes.get(name) != expected
                for name, expected in fixture["raw_fixture_array_sha256"].items()
            )
            or metadata.get("captured_source_residual")
            != _fraction_float(config["cap1_outcome"]["baseline_complete_residual_infinity"])
            or metadata.get("FGCQR_outcome_read") is not False
            or metadata.get("mechanism_question_answered") is not False
        ):
            raise ValueError("CAP1 raw fixture metadata differs")
        for name, expected_hash in fixture["raw_fixture_array_sha256"].items():
            if array_content_sha256(raw[name]) != expected_hash:
                raise ValueError(f"CAP1 raw fixture array hash differs: {name}")
        for name in ("u", "p", "q", "p_r", "q_r"):
            if not np.array_equal(raw[name][0], compact[name]):
                raise ValueError(f"compact point differs from CAP1 raw {name}")
        if (
            raw["radii"][0] != compact["radius"]
            or not np.array_equal(raw["forward_constant"][0], compact["forward_constant"])
            or not np.array_equal(raw["forward_jacobian"][0], compact["forward_jacobian"])
            or not np.array_equal(
                raw["baseline_acceleration"][0], compact["baseline_acceleration"]
            )
            or not np.array_equal(
                raw["baseline_complete_residual"][0],
                compact["baseline_complete_residual"],
            )
        ):
            raise ValueError("compact point differs from the CAP1 raw source records")


def _coefficient_diagnosis(
    compact: Mapping[str, np.ndarray | float],
) -> dict[str, Any]:
    lower = tuple(
        compact[name] for name in ("u", "p", "q", "p_r", "q_r")
    )
    radius = compact["radius"]
    assert all(isinstance(item, np.ndarray) for item in lower)
    assert isinstance(radius, float)
    exact_constant, exact_matrix = exact_affine_system(*lower, radius)
    exact_constant_float = np.asarray(
        [float(value) for value in exact_constant], dtype=np.float64
    )
    exact_matrix_float = np.asarray(
        [[float(value) for value in row] for row in exact_matrix],
        dtype=np.float64,
    )
    captured_constant = compact["forward_constant"]
    captured_matrix = compact["forward_jacobian"]
    assert isinstance(captured_constant, np.ndarray)
    assert isinstance(captured_matrix, np.ndarray)
    return {
        "exact_constant_binary64_hex": [
            value.hex() for value in exact_constant_float
        ],
        "exact_matrix_binary64_hex": [
            [value.hex() for value in row] for row in exact_matrix_float
        ],
        "captured_forward_constant_minus_exact_constant_infinity": float(
            np.max(np.abs(captured_constant - exact_constant_float), initial=0.0)
        ),
        "captured_forward_matrix_minus_exact_matrix_infinity": float(
            np.max(np.abs(captured_matrix - exact_matrix_float), initial=0.0)
        ),
        "exact_matrix_binary64_infinity_condition": float(
            np.linalg.cond(exact_matrix_float, np.inf)
        ),
        "captured_matrix_binary64_infinity_condition": float(
            np.linalg.cond(captured_matrix, np.inf)
        ),
        "constant_term_cancellation_dominates_captured_coefficient_error": bool(
            np.max(np.abs(captured_constant - exact_constant_float), initial=0.0)
            > np.max(np.abs(captured_matrix - exact_matrix_float), initial=0.0)
        ),
    }


def _validate_outcome(
    config: Mapping[str, Any],
    compact: Mapping[str, np.ndarray | float],
    oracle: Mapping[str, Any],
    coefficients: Mapping[str, Any],
) -> None:
    cap1 = config["cap1_outcome"]
    contract = config["exact_oracle"]
    baseline = compact["baseline_complete_residual"]
    assert isinstance(baseline, np.ndarray)
    direct_baseline = np.asarray(
        oracle["floating_residuals"]["baseline"]["direct"], dtype=np.float64
    )
    if not np.array_equal(baseline, direct_baseline):
        raise ValueError("compact CAP1 baseline residual differs from the direct evaluator")
    raw_limit = _fraction_float(cap1["raw_complete_residual_maximum"])
    floating_minimum = _fraction_float(
        contract["floating_evaluator_false_rejection_minimum"]
    )
    if (
        float(np.max(np.abs(baseline), initial=0.0))
        != _fraction_float(cap1["baseline_complete_residual_infinity"])
        or oracle.get("exact_affine_nonsingular") is not True
        or oracle.get("exact_root_residual_is_zero") is not True
        or oracle.get("rounded_exact_root_residual_infinity")
        >= _fraction_float(contract["exact_rounded_root_residual_maximum"])
        or oracle["floating_residual_infinities"]["rounded_exact_root_direct"]
        <= floating_minimum
        or oracle["floating_residual_infinities"]["rounded_exact_root_generic"]
        <= floating_minimum
        or oracle["floating_residual_infinities"]["rounded_exact_root_direct"]
        <= raw_limit
        or oracle["floating_residual_infinities"]["rounded_exact_root_generic"]
        <= raw_limit
        or coefficients["constant_term_cancellation_dominates_captured_coefficient_error"]
        is not True
        or coefficients["captured_matrix_binary64_infinity_condition"]
        >= _fraction_float(cap1["maximum_captured_kinetic_condition_infinity"])
    ):
        raise ValueError("PREF11 exact-oracle outcome differs from the frozen burden")
    lower = tuple(
        compact[name] for name in ("u", "p", "q", "p_r", "q_r")
    )
    radius = compact["radius"]
    assert all(isinstance(item, np.ndarray) for item in lower)
    assert isinstance(radius, float)
    forward = assemble_forward_affine_system(
        *(np.asarray(item, dtype=np.float64)[None, :] for item in lower),
        np.asarray((radius,), dtype=np.float64),
    )
    captured_constant = compact["forward_constant"]
    captured_matrix = compact["forward_jacobian"]
    assert isinstance(captured_constant, np.ndarray)
    assert isinstance(captured_matrix, np.ndarray)
    if (
        not np.array_equal(forward.constant[0], captured_constant)
        or not np.array_equal(forward.jacobian[0], captured_matrix)
    ):
        raise ValueError("compact captured affine coefficients do not reproduce")


def reproduce(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    config = load_config(config_path)
    fixture, compact = _validate_compact_fixture(config)
    _validate_raw_bundle_if_present(config, fixture, compact)
    oracle = diagnose_exact_oracle(
        *(compact[name] for name in ("u", "p", "q", "p_r", "q_r")),
        compact["radius"],
        compact["baseline_acceleration"],
    )
    coefficients = _coefficient_diagnosis(compact)
    _validate_outcome(config, compact, oracle, coefficients)
    lineage = config["lineage"]
    compact_path = REPOSITORY / lineage["compact_point_fixture"]
    cal8_path = REPOSITORY / lineage["cal8_result"]
    cap1_config_path = REPOSITORY / lineage["cap1_config"]
    return {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": config["project_version"],
        "classification": "captured_GR0_point0_binary64_evaluator_cancellation_diagnosis",
        "generated_by": _rel(Path(__file__)),
        "derivation_document": _rel(OWNER_DOCUMENT),
        "derivation_document_sha256": _sha(OWNER_DOCUMENT),
        "source_config_sha256": {
            _rel(config_path): _sha(config_path),
            _rel(compact_path): _sha(compact_path),
        },
        "predecessor_sha256": {
            _rel(cal8_path): _sha(cal8_path),
            _rel(cap1_config_path): _sha(cap1_config_path),
        },
        "implementation_sha256": {
            _rel(path): _sha(path) for path in IMPLEMENTATION
        },
        "scope_bindings": dict(config["scope"]),
        "gate_status": dict(EXPECTED_CLAIMS),
        "artifact_payload": {
            "cap1_raw_bundle_binding": {
                "git_commit": lineage["cap1_git_commit"],
                "manifest_sha256": lineage["cap1_manifest_sha256"],
                "replay_result_sha256": lineage["cap1_replay_result_sha256"],
                "fixture_sha256": lineage["cap1_raw_fixture_sha256"],
                "arithmetic_result_sha256": lineage["cap1_arithmetic_result_sha256"],
                "raw_bundle_optional_for_clean_clone_reproduction": True,
                "partial_or_hash_mismatched_raw_bundle_fails_closed": True,
            },
            "cap1_frozen_outcome": dict(config["cap1_outcome"]),
            "compact_point_fixture": {
                "artifact_id": fixture["artifact_id"],
                "sha256": _sha(compact_path),
                "point_index": fixture["point_index"],
                "coordinate_radius_binary64_hex": fixture["coordinate_radius"],
                "source_call": dict(fixture["source_call"]),
                "literal_binary64_bits_are_the_only_point_input": True,
            },
            "exact_oracle": oracle,
            "coefficient_diagnosis": coefficients,
            "successor_boundary": dict(config["successor_boundary"]),
            "epistemic_boundary": {
                "established": [
                    "CAP1 reproduced the complete implicated GR-0 source-wall trace and none of its frozen arithmetic routes cleared the unchanged direct residual gate",
                    "the captured dyadic REF1 system at point zero is exactly affine and nonsingular",
                    "its exact root zeros all six rows and its binary64 rounding has an exact residual below 1e-26",
                    "both current binary64 residual implementations falsely reject that rounded exact root above 1e-12",
                    "the specific captured point-zero wall is an evaluator-cancellation failure rather than a demonstrated equation or physical obstruction",
                ],
                "not_established": [
                    "a cancellation-resistant full-grid runtime source evaluator",
                    "an eligible GR-0 calibration amplitude",
                    "resolution convergence for the independent amplitude-3 spectral veto",
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
