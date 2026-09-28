#!/usr/bin/env python3
"""Bind the post-freeze TDG5 stage-complete refinement theorem."""

from __future__ import annotations

import argparse
from copy import deepcopy
from fractions import Fraction
from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import tomllib
from typing import Any, Callable, Mapping


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.tdg5_stage_complete_refinement import (  # noqa: E402
    stage_complete_temporal_refinement_preflight,
)
from recursive_horizons.fgc.evolution.tdg5_stage_complete_refinement_theorem import (  # noqa: E402
    exact_derivative_root_candidates,
    independent_hermite_coefficients,
    numerical_engine_shape_certificate,
    restrict_cubic_to_half,
    stage_complete_refinement_theorem_certificate,
)


ARTIFACT_ID = "FGC-1-TDG5-PREF20"
PROJECT_VERSION = "0.11.0"
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-tdg5-pref20.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-tdg5-pref20.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-tdg5-pref20.md"
FREEZE_RESULT = REPOSITORY / "results/fgc-1-tdg5-frz1.json"
FREEZE_CONFIG = REPOSITORY / "configs/fgc/fgc-1-tdg5-frz1.toml"
FREEZE_REPRODUCER = REPOSITORY / "scripts/reproduce_fgc_tdg5_frz1.py"
FREEZE_MODULE = (
    REPOSITORY
    / "src/recursive_horizons/fgc/evolution/tdg5_stage_complete_refinement.py"
)
FREEZE_DOCUMENT = REPOSITORY / "docs/fgc-tdg5-frz1.md"
NUMERICAL_ENGINE = (
    REPOSITORY / "src/recursive_horizons/fgc/evolution/numerical_engine.py"
)
THEOREM_MODULE = (
    REPOSITORY
    / "src/recursive_horizons/fgc/evolution/tdg5_stage_complete_refinement_theorem.py"
)
IMPLEMENTATION = (Path(__file__).resolve(), THEOREM_MODULE)


EXPECTED_CLAIMS = {
    "TDG4_sampling_identifiability_theorem_completed": True,
    "finite_samples_alone_proved_insufficient_for_continuum_top_band_bound": True,
    "current_temporal_spectral_rule_retired_as_continuum_admission": True,
    "current_temporal_spectral_rule_retained_as_sampled_diagnostic": True,
    "replacement_temporal_gate_design_authorized": True,
    "TDG5_replacement_design_frozen": True,
    "TDG5_exact_no_history_controls_pass": True,
    "TDG5_stage_complete_refinement_theorem_completed": True,
    "runtime_refinement_pair_implementation_authorized": True,
    "runtime_refinement_pair_implemented": False,
    "runtime_thresholds_frozen": False,
    "replacement_temporal_admission_defined": False,
    "PROTO14_frozen": False,
    "fresh_GR0_dynamic_calibration_completed": False,
    "GR0_case_eligible": False,
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
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _bytes_sha(payload: bytes) -> str:
    return sha256(payload).hexdigest()


def _canonical(value: Mapping[str, Any]) -> str:
    return json.dumps(
        value, indent=2, sort_keys=True, ensure_ascii=True, allow_nan=False
    ) + "\n"


def _canonical_bytes(value: Mapping[str, Any]) -> bytes:
    return _canonical(value).encode("utf-8")


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


def _git(*args: str, check: bool = True) -> subprocess.CompletedProcess[bytes]:
    return subprocess.run(
        ["git", *args],
        cwd=REPOSITORY,
        check=check,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )


def _strict_keys(name: str, value: Mapping[str, Any], expected: set[str]) -> None:
    if set(value) != expected:
        raise ValueError(f"{name} keys differ")


def _expected_scope() -> dict[str, Any]:
    return {
        "target": "FGC-2-SF1-TDG5",
        "role": "post_freeze_stage_complete_refinement_theorem_binder",
        "calibration_branch": "GR-0",
        "freeze_compact_result_read": True,
        "actual_terminal_history_arrays_consumed": False,
        "campaign_checkpoint_loaded": False,
        "campaign_checkpoint_resumed": False,
        "state_advanced": False,
        "SGBL_trajectory_read": False,
        "FGCQR_trajectory_read": False,
        "runtime_refinement_pair_implemented": False,
        "runtime_thresholds_frozen": False,
        "replacement_temporal_admission_defined": False,
        "physical_or_candidate_question_answered": False,
    }


def _expected_lineage() -> dict[str, Any]:
    return {
        "freeze_commit": "1716cc4acbc95a5b11e6e162e405c43330d389a6",
        "freeze_config_sha256": "90f2d548a0cc9062615e69b16c175ad9ecf2db6d49d97ecf35697b9173a50878",
        "freeze_result_sha256": "2af573edda929424542129a8391c8bcf34c081206ec73d8c4dd040fbe47a4c04",
        "freeze_reproducer_sha256": "c976d553d9bbf8c730f0fe0c7f5eb2ef0f5fcbb28b59b10e48b7e5ebb1a2d303",
        "freeze_module_sha256": "56cf8607adf6aa769eeade329feadf1dd0fc963002de7cc970961b1aa9ba6d55",
        "freeze_document_sha256": "f549406b381e1be70d04ea741c7c86d9c0e09140fdc104f876697bae4389c795",
        "numerical_engine_sha256": "8ea99604a85b4e1acbf869ed2462ea4aa8301c06e068330a1d0ce66a51a0bebf",
        "freeze_commit_must_be_ancestor_of_HEAD": True,
        "all_bound_blobs_must_match_freeze_commit_and_worktree": True,
        "freeze_result_must_be_canonical_and_prospective": True,
        "raw_campaign_bundle_must_not_be_loaded": True,
    }


def _expected_theorem() -> dict[str, Any]:
    return {
        "declared_interval_class": (
            "P3_on_each_accepted_normalized_time_interval"
        ),
        "data_order": [
            "left_value",
            "left_scaled_slope",
            "right_value",
            "right_scaled_slope",
        ],
        "sampling_determinant": "1",
        "basis_determinant": "1",
        "sampling_infinity_norm": "6",
        "basis_infinity_norm": "9",
        "condition_number_infinity_upper_bound": "54",
        "left_half_restriction_determinant": "1/64",
        "right_half_restriction_determinant": "1/64",
        "maximum_candidate_count_per_half_interval": 4,
        "complete_candidate_set": (
            "endpoints_plus_every_real_interior_derivative_root"
        ),
        "endpoint_only_bump_absolute_maximum": "1",
        "RK4_formal_order": 4,
        "RK4_Richardson_denominator": 15,
        "SSPRK3_formal_order": 3,
        "SSPRK3_Richardson_denominator": 7,
        "immutable_engine_exposes_complete_stage_and_endpoint_records": True,
        "fresh_candidate_endpoint_RHS_record_is_available": True,
        "runtime_pair_can_be_built_without_mutating_existing_step_proposal_shape": True,
        "declared_cubic_is_exact_PDE_history": False,
        "trajectory_asymptotic_regime_proved": False,
        "rigorous_global_PDE_error_bound_proved": False,
        "binary64_extremum_evaluator_implemented": False,
        "atomic_coarse_fine_transaction_implemented": False,
    }


def _expected_authorization() -> dict[str, Any]:
    return {
        "authorized_object": (
            "implementation_and_synthetic_validation_of_the_frozen_same_grid_refinement_pair"
        ),
        "complete_dimensionless_u_p_q_owned_state_required": True,
        "same_bitwise_initial_state_required": True,
        "same_source_projector_and_spatial_operator_required": True,
        "coarse_one_full_step_is_evidence_only": True,
        "fine_two_half_step_path_is_the_only_eventual_committable_path": True,
        "fine_two_half_steps_must_commit_atomically_or_neither_commits": True,
        "unchanged_stage_source_health_CFL_boundary_scale_and_transaction_guards_required": True,
        "continuous_difference_must_include_all_endpoints_and_real_interior_derivative_roots": True,
        "binary64_root_and_value_evaluation_must_add_an_outward_arithmetic_debit": True,
        "thresholds_require_a_separate_prospective_freeze": True,
        "production_trajectory_authorized": False,
        "historical_PROTO13_stop_preserved": True,
        "historical_CAL10_result_reclassified": False,
    }


def _expected_limitations() -> dict[str, Any]:
    return {
        "Hermite_extension_is_a_declared_numerical_object_not_the_exact_PDE_history": True,
        "finite_dimensional_injectivity_does_not_prove_PDE_accuracy": True,
        "Richardson_factor_assumes_a_shared_leading_error_coefficient": True,
        "Richardson_factor_assumes_the_tested_step_is_in_the_asymptotic_regime": True,
        "Richardson_debit_is_not_a_rigorous_global_PDE_bound": True,
        "binary64_continuous_extremum_implementation_remains_open": True,
        "spatial_ladders_and_cross_method_agreement_remain_required": True,
        "complete_DEF1_error_ledger_remains_required": True,
        "observed_trapped_or_Raychaudhuri_margin_may_not_set_the_threshold": True,
    }


def _expected_successor() -> dict[str, Any]:
    return {
        "TDG4_sampling_identifiability_theorem_completed": True,
        "current_temporal_spectral_rule_retired_as_continuum_admission": True,
        "current_temporal_spectral_rule_retained_as_sampled_diagnostic": True,
        "replacement_temporal_gate_design_authorized": True,
        "TDG5_replacement_design_frozen": True,
        "TDG5_exact_no_history_controls_pass": True,
        "TDG5_stage_complete_refinement_theorem_completed": True,
        "runtime_refinement_pair_implementation_authorized": True,
        "runtime_refinement_pair_implemented": False,
        "runtime_thresholds_frozen": False,
        "replacement_temporal_admission_defined": False,
        "PROTO14_frozen": False,
        "fresh_GR0_dynamic_calibration_completed": False,
        "GR0_case_eligible": False,
        "SGBL_execution_authorized": False,
        "FGCQR_holdout_execution_authorized": False,
        "DEF1_execution_authorized": False,
    }


def validate_config_data(config: Mapping[str, Any]) -> None:
    _strict_keys(
        "PREF20 config",
        config,
        {
            "schema_version",
            "artifact_id",
            "project_version",
            "metric_signature",
            "riemann_convention",
            "freeze_result",
            "freeze_config",
            "freeze_reproducer",
            "freeze_module",
            "freeze_document",
            "numerical_engine",
            "scope",
            "immutable_lineage",
            "theorem_result",
            "runtime_authorization",
            "limitation_audit",
            "proof_contract",
            "successor_boundary",
            "claims",
        },
    )
    expected_paths = {
        "freeze_result": _rel(FREEZE_RESULT),
        "freeze_config": _rel(FREEZE_CONFIG),
        "freeze_reproducer": _rel(FREEZE_REPRODUCER),
        "freeze_module": _rel(FREEZE_MODULE),
        "freeze_document": _rel(FREEZE_DOCUMENT),
        "numerical_engine": _rel(NUMERICAL_ENGINE),
    }
    if (
        config.get("schema_version") != 1
        or config.get("artifact_id") != ARTIFACT_ID
        or config.get("project_version") != PROJECT_VERSION
        or config.get("metric_signature") != "-+++"
        or config.get("riemann_convention") != "plus_partial_mu_gamma_nu"
        or any(config.get(key) != value for key, value in expected_paths.items())
        or config.get("scope") != _expected_scope()
        or config.get("immutable_lineage") != _expected_lineage()
        or config.get("theorem_result") != _expected_theorem()
        or config.get("runtime_authorization") != _expected_authorization()
        or config.get("limitation_audit") != _expected_limitations()
        or config.get("successor_boundary") != _expected_successor()
        or config.get("claims") != EXPECTED_CLAIMS
    ):
        raise ValueError("PREF20 identity, lineage, theorem, or boundary differs")
    if set(config.get("proof_contract", {}).values()) != {True}:
        raise ValueError("PREF20 proof contract must be all-of")


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    with path.open("rb") as handle:
        config = tomllib.load(handle)
    validate_config_data(config)
    return config


def _verify_lineage(config: Mapping[str, Any]) -> dict[str, Any]:
    lineage = config["immutable_lineage"]
    commit = lineage["freeze_commit"]
    _git("cat-file", "-e", f"{commit}^{{commit}}")
    if _git("merge-base", "--is-ancestor", commit, "HEAD", check=False).returncode:
        raise ValueError("PREF20 freeze commit is not an ancestor of HEAD")
    tracked = {
        _rel(FREEZE_CONFIG): lineage["freeze_config_sha256"],
        _rel(FREEZE_RESULT): lineage["freeze_result_sha256"],
        _rel(FREEZE_REPRODUCER): lineage["freeze_reproducer_sha256"],
        _rel(FREEZE_MODULE): lineage["freeze_module_sha256"],
        _rel(FREEZE_DOCUMENT): lineage["freeze_document_sha256"],
        _rel(NUMERICAL_ENGINE): lineage["numerical_engine_sha256"],
    }
    for relative, expected_hash in tracked.items():
        committed = _git("show", f"{commit}:{relative}").stdout
        current = (REPOSITORY / relative).read_bytes()
        if (
            committed != current
            or _bytes_sha(committed) != expected_hash
            or _bytes_sha(current) != expected_hash
        ):
            raise ValueError(f"PREF20 bound blob differs: {relative}")

    freeze_raw = FREEZE_RESULT.read_bytes()
    freeze = json.loads(freeze_raw)
    status = freeze.get("gate_status", {})
    boundary = freeze.get("artifact_payload", {}).get("claim_boundary", {})
    if (
        freeze_raw != _canonical_bytes(freeze)
        or freeze.get("artifact_id") != "FGC-1-TDG5-FRZ1"
        or freeze.get("classification")
        != "prospective_TDG5_stage_complete_temporal_refinement_design_freeze"
        or status.get("TDG5_replacement_design_frozen") is not True
        or status.get("TDG5_exact_no_history_controls_pass") is not True
        or status.get("TDG5_stage_complete_refinement_theorem_completed")
        is not False
        or status.get("runtime_refinement_pair_implementation_authorized")
        is not False
        or status.get("runtime_refinement_pair_implemented") is not False
        or status.get("runtime_thresholds_frozen") is not False
        or status.get("replacement_temporal_admission_defined") is not False
        or boundary.get("actual_terminal_history_arrays_consumed") is not False
        or boundary.get("campaign_checkpoint_loaded") is not False
        or boundary.get("state_advanced") is not False
        or boundary.get("historical_PROTO13_stop_preserved") is not True
        or boundary.get("historical_CAL10_result_reclassified") is not False
    ):
        raise ValueError("PREF20 prospective freeze record differs")

    freeze_config = tomllib.loads(FREEZE_CONFIG.read_text(encoding="utf-8"))
    frozen_design = freeze_config.get("frozen_replacement_design", {})
    same_grid = frozen_design.get("same_grid_step_refinement", {})
    if (
        freeze_config.get("artifact_id") != "FGC-1-TDG5-FRZ1"
        or frozen_design.get("monitored_state") != "dimensionless_u_p_q"
        or frozen_design.get("same_spatial_operator_required") is not True
        or frozen_design.get("same_initial_state_required") is not True
        or same_grid.get("RK4_formal_order") != 4
        or same_grid.get("RK4_Richardson_denominator") != 15
        or same_grid.get("SSPRK3_formal_order") != 3
        or same_grid.get("SSPRK3_Richardson_denominator") != 7
        or same_grid.get("fine_path_is_the_only_committable_path") is not True
        or same_grid.get("coarse_path_is_evidence_only") is not True
    ):
        raise ValueError("PREF20 frozen design contract differs")
    return {
        "freeze_commit": commit,
        "freeze_commit_is_ancestor_of_HEAD": True,
        "bound_blob_sha256": tracked,
        "freeze_result_verified_canonical_and_prospective": True,
        "frozen_complete_state": frozen_design["monitored_state"],
        "frozen_same_initial_state_required": True,
        "frozen_same_spatial_operator_required": True,
        "frozen_method_orders": {"RK4": 4, "SSPRK3": 3},
        "frozen_Richardson_denominators": {"RK4": 15, "SSPRK3": 7},
        "historical_PROTO13_stop_preserved": True,
        "historical_CAL10_result_reclassified": False,
        "actual_terminal_history_arrays_consumed": False,
        "campaign_checkpoint_loaded": False,
        "campaign_checkpoint_resumed": False,
        "state_advanced": False,
    }


def _freeze_crosscheck(theorem: Mapping[str, Any]) -> dict[str, Any]:
    preflight = stage_complete_temporal_refinement_preflight()
    if (
        preflight.get("all_exact_controls_pass") is not True
        or preflight.get("actual_terminal_histories_consumed") is not False
        or preflight.get("campaign_checkpoint_loaded") is not False
        or preflight.get("state_advanced") is not False
        or preflight.get("runtime_refinement_pair_implemented") is not False
        or preflight.get("runtime_thresholds_frozen") is not False
        or preflight.get("replacement_temporal_admission_defined") is not False
    ):
        raise ValueError("PREF20 freeze preflight boundary differs")
    frozen_map = preflight["Hermite_data_map"]
    theorem_map = theorem["independent_Hermite_map"]
    alias = preflight["endpoint_only_alias_control"]
    extrema = theorem["continuous_extremum_theorem"]
    frozen_methods = preflight["method_controls"]
    theorem_methods = theorem["method_owned_Richardson_debits"]
    if (
        frozen_map.get("determinant")
        != theorem_map.get("sampling_determinant_by_fraction_Gaussian_elimination")
        or frozen_map.get("sampling_infinity_norm")
        != theorem_map.get("sampling_infinity_norm")
        or frozen_map.get("inverse_infinity_norm")
        != theorem_map.get("basis_infinity_norm")
        or frozen_map.get("condition_number_infinity_upper_bound")
        != theorem_map.get("condition_number_infinity_upper_bound")
        or alias.get("adversarial_midpoint_value")
        != extrema.get("endpoint_only_bump_control", {}).get("absolute_maximum")
    ):
        raise ValueError("PREF20 Hermite or alias crosscheck differs")
    for method in ("RK4", "SSPRK3"):
        frozen = frozen_methods[method]
        independent = theorem_methods[method]
        if (
            frozen["order"]["formal_order"]
            != independent["formal_order"]
            or frozen["step_doubling"]["Richardson_denominator"]
            != independent["Richardson_denominator"]
        ):
            raise ValueError(f"PREF20 {method} crosscheck differs")
    return {
        "freeze_exact_preflight_reexecuted": True,
        "independent_basis_expansion_matches_frozen_Hermite_map": True,
        "independent_elimination_matches_frozen_determinant_and_norms": True,
        "endpoint_alias_maximum_matches_frozen_midpoint_control": True,
        "both_method_orders_and_Richardson_denominators_match": True,
        "freeze_preflight_is_a_crosscheck_not_the_theorem_authority": True,
        "actual_terminal_histories_consumed": False,
    }


def _expect_config_rejection(
    config: Mapping[str, Any], mutate: Callable[[dict[str, Any]], None]
) -> bool:
    attacked = deepcopy(config)
    mutate(attacked)
    try:
        validate_config_data(attacked)
    except (TypeError, ValueError):
        return True
    return False


def _mutation_controls(config: Mapping[str, Any]) -> dict[str, bool]:
    short_data_rejected = False
    invalid_half_rejected = False
    irrational_root_requires_future_evaluator = False
    try:
        independent_hermite_coefficients((0, 1, 2))
    except (TypeError, ValueError):
        short_data_rejected = True
    try:
        restrict_cubic_to_half((0, 1, 2, 3), half=2)
    except (TypeError, ValueError):
        invalid_half_rejected = True
    try:
        exact_derivative_root_candidates((0, Fraction(-1, 2), 0, 1))
    except (TypeError, ValueError):
        irrational_root_requires_future_evaluator = True

    engine_source = NUMERICAL_ENGINE.read_text(encoding="utf-8")
    attacked_source = engine_source.replace(
        'evaluate("candidate_endpoint", final_time, candidate)',
        'evaluate("candidate_shadow", final_time, candidate)',
    )
    attacked_engine = numerical_engine_shape_certificate(attacked_source)
    controls = {
        "short_Hermite_data_mutation_rejected": short_data_rejected,
        "invalid_half_interval_mutation_rejected": invalid_half_rejected,
        "irrational_root_requires_future_outward_evaluator": (
            irrational_root_requires_future_evaluator
        ),
        "candidate_endpoint_stage_mutation_rejected": (
            attacked_engine[
                "runtime_pair_can_be_built_without_mutating_the_existing_step_proposal_shape"
            ]
            is False
        ),
        "freeze_commit_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["immutable_lineage"].__setitem__(
                "freeze_commit", "0" * 40
            ),
        ),
        "freeze_result_hash_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["immutable_lineage"].__setitem__(
                "freeze_result_sha256", "0" * 64
            ),
        ),
        "Hermite_determinant_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["theorem_result"].__setitem__(
                "sampling_determinant", "0"
            ),
        ),
        "restriction_determinant_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["theorem_result"].__setitem__(
                "left_half_restriction_determinant", "1/32"
            ),
        ),
        "endpoint_RHS_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["theorem_result"].__setitem__(
                "fresh_candidate_endpoint_RHS_record_is_available", False
            ),
        ),
        "implemented_runtime_promotion_rejected": _expect_config_rejection(
            config,
            lambda value: value["claims"].__setitem__(
                "runtime_refinement_pair_implemented", True
            ),
        ),
        "threshold_promotion_rejected": _expect_config_rejection(
            config,
            lambda value: value["claims"].__setitem__(
                "runtime_thresholds_frozen", True
            ),
        ),
        "replacement_admission_promotion_rejected": _expect_config_rejection(
            config,
            lambda value: value["claims"].__setitem__(
                "replacement_temporal_admission_defined", True
            ),
        ),
        "GR0_eligibility_promotion_rejected": _expect_config_rejection(
            config,
            lambda value: value["claims"].__setitem__(
                "GR0_case_eligible", True
            ),
        ),
        "candidate_promotion_rejected": _expect_config_rejection(
            config,
            lambda value: value["claims"].__setitem__(
                "FGCQR_holdout_execution_authorized", True
            ),
        ),
        "historical_reclassification_mutation_rejected": (
            _expect_config_rejection(
                config,
                lambda value: value["runtime_authorization"].__setitem__(
                    "historical_CAL10_result_reclassified", True
                ),
            )
        ),
    }
    if set(controls.values()) != {True}:
        raise ValueError(f"PREF20 mutation control failed: {controls}")
    return controls


def record(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    config = load_config(config_path)
    lineage = _verify_lineage(config)
    engine_source = NUMERICAL_ENGINE.read_text(encoding="utf-8")
    theorem = stage_complete_refinement_theorem_certificate(engine_source)
    if (
        theorem.get("all_independent_theorem_controls_pass") is not True
        or theorem.get("TDG5_stage_complete_refinement_theorem_completed")
        is not True
        or theorem.get("runtime_refinement_pair_implementation_authorized")
        is not True
        or theorem.get("runtime_refinement_pair_implemented") is not False
        or theorem.get("runtime_thresholds_frozen") is not False
        or theorem.get("replacement_temporal_admission_defined") is not False
        or theorem.get("PROTO14_frozen") is not False
        or theorem.get("GR0_case_eligible") is not False
        or theorem.get("actual_terminal_history_arrays_consumed") is not False
        or theorem.get("campaign_checkpoint_loaded") is not False
        or theorem.get("state_advanced") is not False
        or theorem.get("SGBL_trajectory_read") is not False
        or theorem.get("FGCQR_trajectory_read") is not False
        or theorem.get("physical_or_candidate_question_answered") is not False
    ):
        raise ValueError("PREF20 theorem certificate differs")
    crosscheck = _freeze_crosscheck(theorem)
    mutations = _mutation_controls(config)
    result = {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "metric_signature": "-+++",
        "riemann_convention": "plus_partial_mu_gamma_nu",
        "classification": (
            "completed_TDG5_stage_complete_refinement_theorem_"
            "runtime_pair_implementation_authorized"
        ),
        "generated_by": _rel(Path(__file__)),
        "gate_status": dict(EXPECTED_CLAIMS),
        "nonclaims": {
            key: False
            for key, value in sorted(EXPECTED_CLAIMS.items())
            if value is False
        },
        "source_config_sha256": {_rel(config_path): _sha(config_path)},
        "implementation_sha256": {
            _rel(path): _sha(path) for path in IMPLEMENTATION
        },
        "artifact_payload": {
            "immutable_lineage": lineage,
            "theorem_certificate": theorem,
            "freeze_crosscheck": crosscheck,
            "runtime_authorization": deepcopy(config["runtime_authorization"]),
            "limitation_audit": deepcopy(config["limitation_audit"]),
            "mutation_controls": mutations,
            "proof_contract": dict(config["proof_contract"]),
            "successor_boundary": dict(config["successor_boundary"]),
            "claim_boundary": {
                "actual_terminal_history_arrays_consumed": False,
                "campaign_checkpoint_loaded": False,
                "state_advanced": False,
                "historical_PROTO13_stop_preserved": True,
                "historical_CAL10_result_reclassified": False,
                "declared_cubic_is_exact_PDE_history": False,
                "trajectory_asymptotic_regime_proved": False,
                "rigorous_global_PDE_error_bound_proved": False,
                "stage_complete_refinement_theorem_completed": True,
                "runtime_refinement_pair_implementation_authorized": True,
                "runtime_refinement_pair_implemented": False,
                "runtime_thresholds_frozen": False,
                "replacement_temporal_admission_defined": False,
                "new_production_trajectory_authorized": False,
                "GR0_case_eligible": False,
                "candidate_branch_opened": False,
                "physical_question_answered": False,
            },
        },
        "derivation_document": _rel(OWNER_DOCUMENT),
        "derivation_document_sha256": _sha(OWNER_DOCUMENT),
    }
    return json.loads(
        json.dumps(result, sort_keys=True, ensure_ascii=True, allow_nan=False)
    )


def _validate_stored_record(
    value: Mapping[str, Any], config: Mapping[str, Any]
) -> None:
    payload = value.get("artifact_payload", {})
    theorem = payload.get("theorem_certificate", {})
    boundary = payload.get("claim_boundary", {})
    if (
        value.get("schema_version") != 1
        or value.get("artifact_id") != ARTIFACT_ID
        or value.get("project_version") != PROJECT_VERSION
        or value.get("classification")
        != "completed_TDG5_stage_complete_refinement_theorem_runtime_pair_implementation_authorized"
        or value.get("generated_by") != "scripts/reproduce_fgc_tdg5_pref20.py"
        or value.get("gate_status") != EXPECTED_CLAIMS
        or value.get("nonclaims")
        != {
            key: False
            for key, item in sorted(EXPECTED_CLAIMS.items())
            if item is False
        }
        or value.get("source_config_sha256")
        != {_rel(DEFAULT_CONFIG): _sha(DEFAULT_CONFIG)}
        or payload.get("runtime_authorization")
        != config["runtime_authorization"]
        or payload.get("limitation_audit") != config["limitation_audit"]
        or payload.get("successor_boundary") != config["successor_boundary"]
        or theorem.get("TDG5_stage_complete_refinement_theorem_completed")
        is not True
        or theorem.get("runtime_refinement_pair_implementation_authorized")
        is not True
        or theorem.get("runtime_refinement_pair_implemented") is not False
        or theorem.get("runtime_thresholds_frozen") is not False
        or boundary.get("historical_PROTO13_stop_preserved") is not True
        or boundary.get("historical_CAL10_result_reclassified") is not False
        or boundary.get("declared_cubic_is_exact_PDE_history") is not False
        or boundary.get("trajectory_asymptotic_regime_proved") is not False
        or boundary.get("rigorous_global_PDE_error_bound_proved") is not False
        or boundary.get("runtime_refinement_pair_implementation_authorized")
        is not True
        or boundary.get("runtime_refinement_pair_implemented") is not False
        or boundary.get("runtime_thresholds_frozen") is not False
        or boundary.get("replacement_temporal_admission_defined") is not False
        or boundary.get("new_production_trajectory_authorized") is not False
        or boundary.get("candidate_branch_opened") is not False
    ):
        raise ValueError("PREF20 stored result boundary differs")
    implementation = value.get("implementation_sha256", {})
    if (
        not isinstance(implementation, Mapping)
        or set(implementation) != {_rel(path) for path in IMPLEMENTATION}
        or any(
            _sha(REPOSITORY / relative) != expected
            for relative, expected in implementation.items()
        )
    ):
        raise ValueError("PREF20 implementation hash ledger differs")
    document = REPOSITORY / str(value.get("derivation_document", ""))
    if (
        value.get("derivation_document") != "docs/fgc-tdg5-pref20.md"
        or not document.is_file()
        or _sha(document) != value.get("derivation_document_sha256")
    ):
        raise ValueError("PREF20 derivation document differs")
    if (
        not isinstance(payload.get("mutation_controls"), Mapping)
        or set(payload["mutation_controls"].values()) != {True}
        or not isinstance(payload.get("proof_contract"), Mapping)
        or set(payload["proof_contract"].values()) != {True}
    ):
        raise ValueError("PREF20 proof or mutation ledger differs")


def verify_canonical(
    config_path: Path = DEFAULT_CONFIG,
    output_path: Path = DEFAULT_OUTPUT,
) -> dict[str, Any]:
    config = load_config(config_path)
    raw = output_path.read_bytes()
    stored = json.loads(raw)
    if not isinstance(stored, dict) or raw != _canonical_bytes(stored):
        raise ValueError("PREF20 stored result is not canonical JSON")
    _validate_stored_record(stored, config)
    observed = record(config_path)
    if stored != observed:
        raise ValueError("PREF20 stored result differs from reproduction")
    return stored


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    if args.write and args.verify:
        raise SystemExit("choose either --write or --verify")
    if args.write:
        value = record(args.config)
        _atomic_write(args.output, _canonical_bytes(value))
        print(f"wrote {_rel(args.output)}")
    elif args.verify:
        value = verify_canonical(args.config, args.output)
        print(
            f"PASS {ARTIFACT_ID}: theorem="
            f"{value['gate_status']['TDG5_stage_complete_refinement_theorem_completed']} "
            "runtime_implementation_authorized="
            f"{value['gate_status']['runtime_refinement_pair_implementation_authorized']} "
            "runtime_implemented="
            f"{value['gate_status']['runtime_refinement_pair_implemented']}"
        )
    else:
        print(_canonical(record(args.config)), end="")


if __name__ == "__main__":
    main()
