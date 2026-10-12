#!/usr/bin/env python3
"""Reproduce the prospective FGC-1-TDG5-FRZ1 replacement-design freeze."""

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
    TDG5_HERMITE_DATA_ORDER,
    TDG5_METHOD_ORDERS,
    TDG5_OWNED_ROW_POLICY,
    TDG5_REQUIRED_RUNTIME_LAYERS,
    TDG5_SELECTED_SUFFICIENT_ROUTE,
    butcher_order_certificate,
    evaluate_cubic,
    hermite_coefficients,
    stage_complete_temporal_refinement_preflight,
    step_doubling_certificate,
)


ARTIFACT_ID = "FGC-1-TDG5-FRZ1"
PROJECT_VERSION = "0.11.0"
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-tdg5-frz1.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-tdg5-frz1.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-tdg5-frz1.md"
SOURCE_RESULT = REPOSITORY / "results/fgc-1-tdg4-pref19.json"
SOURCE_CONFIG = REPOSITORY / "configs/fgc/fgc-1-tdg4-pref19.toml"
SOURCE_REPRODUCER = REPOSITORY / "scripts/reproduce_fgc_tdg4_pref19.py"
SOURCE_THEOREM_MODULE = (
    REPOSITORY
    / "src/recursive_horizons/fgc/evolution/tdg4_sampling_identifiability_theorem.py"
)
SOURCE_DOCUMENT = REPOSITORY / "docs/fgc-tdg4-pref19.md"
NUMERICAL_ENGINE = (
    REPOSITORY / "src/recursive_horizons/fgc/evolution/numerical_engine.py"
)
PROTO13_RUNNER = REPOSITORY / "scripts/run_fgc_gr0_calibration_v13.py"
PROTOCOL_V13 = REPOSITORY / "configs/fgc/fgc-2-sf1-protocol-v13.toml"
CAL10_CONFIG = REPOSITORY / "configs/fgc/fgc-1-cal10-pref15.toml"
CAL10_RESULT = REPOSITORY / "results/fgc-1-cal10-pref15.json"
DESIGN_MODULE = (
    REPOSITORY
    / "src/recursive_horizons/fgc/evolution/tdg5_stage_complete_refinement.py"
)
IMPLEMENTATION = (Path(__file__).resolve(), DESIGN_MODULE)


EXPECTED_CLAIMS = {
    "TDG4_sampling_identifiability_theorem_completed": True,
    "finite_samples_alone_proved_insufficient_for_continuum_top_band_bound": True,
    "current_temporal_spectral_rule_retired_as_continuum_admission": True,
    "current_temporal_spectral_rule_retained_as_sampled_diagnostic": True,
    "replacement_temporal_gate_design_authorized": True,
    "TDG5_replacement_design_frozen": True,
    "TDG5_exact_no_history_controls_pass": True,
    "TDG5_stage_complete_refinement_theorem_completed": False,
    "runtime_refinement_pair_implementation_authorized": False,
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
        value,
        indent=2,
        sort_keys=True,
        ensure_ascii=True,
        allow_nan=False,
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
        "role": (
            "prospective_stage_complete_continuous_extension_and_same_grid_"
            "temporal_refinement_design"
        ),
        "calibration_branch": "GR-0",
        "source_compact_result_read": True,
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
        "checkpoint_commit": "b866d24fd253ca56db509a2e011323394ae4147b",
        "source_config_sha256": (
            "22727f144dd1666699b8c4a2d0d1b58d7a482b4611020b71fe7a021cfa39cf98"
        ),
        "source_result_sha256": (
            "0db51678fa0fb32280eb3a352fa8bfacb613de19e0c14bdd99d9e2ed08e9e11b"
        ),
        "source_reproducer_sha256": (
            "ae660d7ad166a833b4f63d58abdb0ad827ae7ae0657b8ea6f922bd7a2a38c31f"
        ),
        "source_theorem_module_sha256": (
            "e703a3db4a7a09dee68a9f74cdddc18f3183e1f890aff67f196367a26e229b42"
        ),
        "source_document_sha256": (
            "cafdea686218647a06286f0a8f324ea40b3fd97ab7434509bdeb4af19299ce0a"
        ),
        "numerical_engine_sha256": (
            "8ea99604a85b4e1acbf869ed2462ea4aa8301c06e068330a1d0ce66a51a0bebf"
        ),
        "PROTO13_runner_sha256": (
            "bdc31c8ab4205ff86a134e9744aaaf8d085344aa1ba6a19cebb3fdbad94286d7"
        ),
        "protocol_v13_sha256": (
            "fc9398bad6ea4448f659188641e6da6bd7c2e295a443fe73bf73eb0308480692"
        ),
        "CAL10_result_sha256": (
            "1ebbdc7222baee8134a6d466589770ee445aee594f10630fc09aa6264cb4ac18"
        ),
        "CAL10_config_sha256": (
            "3661ba8cb6b933ab47daf8897c577b4631544605175dbe15d0c084542cf6c9eb"
        ),
        "checkpoint_commit_must_be_ancestor_of_HEAD": True,
        "all_bound_blobs_must_match_commit_and_worktree": True,
        "source_result_must_be_canonical_and_complete": True,
        "raw_campaign_bundle_must_not_be_loaded": True,
    }


def _expected_frozen_design() -> dict[str, Any]:
    return {
        "selected_sufficient_route": TDG5_SELECTED_SUFFICIENT_ROUTE,
        "replaced_object": (
            "PROTO4_inherited_64_point_temporal_spectrum_as_continuum_admission"
        ),
        "retained_object": (
            "PROTO4_inherited_temporal_spectrum_as_sampled_diagnostic_only"
        ),
        "monitored_state": "dimensionless_u_p_q",
        "owned_rows": "noncentre_rows_excluding_four_projector_owned_outer_rows",
        "time_coordinate": "fixed_modified_harmonic_coordinate_time",
        "same_spatial_operator_required": True,
        "same_initial_state_required": True,
        "normal_flow_tracer_spectrum_required_for_admission": False,
        "normal_flow_tracer_spectrum_remains_public_diagnostic": True,
        "actual_PDE_continuum_solution_identified": False,
        "Hermite_continuous_extension": {
            "declared_interval_class": (
                "P3_on_each_accepted_normalized_time_interval"
            ),
            "data_order": list(TDG5_HERMITE_DATA_ORDER),
            "endpoint_physical_slopes_are_multiplied_by_interval_width": True,
            "sampling_matrix_determinant": "1",
            "sampling_infinity_norm": "6",
            "inverse_infinity_norm": "9",
            "condition_number_infinity_upper_bound": "54",
            "data_map_is_injective": True,
            "continuous_coarse_fine_difference_is_piecewise_cubic": True,
            "maximum_requires_endpoints_and_every_real_derivative_root_in_each_half_interval": True,
            "binary64_extremum_evaluation_requires_an_outward_roundoff_debit": True,
        },
        "same_grid_step_refinement": {
            "coarse_path": "one_full_step_of_width_Delta_t",
            "fine_path": "two_consecutive_half_steps_of_width_Delta_t_over_2",
            "RK4_formal_order": 4,
            "RK4_Richardson_denominator": 15,
            "SSPRK3_formal_order": 3,
            "SSPRK3_Richardson_denominator": 7,
            "fine_path_is_the_only_committable_path": True,
            "coarse_path_is_evidence_only": True,
            "each_method_owns_its_own_refinement_debit": True,
            "each_method_owned_spatial_ladder_remains_required": True,
            "cross_method_agreement_remains_required": True,
            "formal_order_and_scalar_contraction_controls_do_not_prove_trajectory_asymptotics": True,
        },
        "transaction_ownership": {
            "coarse_and_fine_paths_must_start_from_one_bitwise_identical_accepted_state": True,
            "coarse_and_fine_paths_must_use_one_source_projector_and_spatial_operator": True,
            "both_paths_must_pass_unchanged_source_health_CFL_boundary_and_scale_premises": True,
            "shadow_path_may_not_mutate_runtime_monitor_causal_tracer_or_retry_state": True,
            "fine_two_half_steps_commit_atomically_or_neither_commits": True,
            "rejected_proposals_preserve_complete_typed_evidence": True,
            "projector_owned_rows_remain_under_existing_exact_projector_and_constraint_gates": True,
        },
        "error_semantics": {
            "continuous_extension_removes_hidden_between_record_freedom_only_inside_the_declared_numerical_class": True,
            "same_grid_difference_is_a_temporal_refinement_measure_not_a_spatial_measure": True,
            "Richardson_factor_is_a_numerical_error_debit_not_a_rigorous_PDE_error_bound": True,
            "debit_thresholds_must_be_frozen_from_independent_controls_before_any_successor_trajectory": True,
            "debit_may_not_be_normalized_by_the_observed_trapped_or_Raychaudhuri_margin": True,
            "complete_DEF1_error_ledger_must_keep_temporal_spatial_source_constraint_affine_null_interpolation_boundary_and_roundoff_terms_separate": True,
            "successful_TDG5_design_does_not_reclassify_CAL10_or_make_GR0_eligible": True,
        },
    }


def _expected_alternatives() -> dict[str, Any]:
    return {
        "full_semidiscrete_Gronwall_route": (
            "deferred_because_a_global_infinity_norm_bound_scales_with_inverse_"
            "grid_spacing_and_is_expected_to_be_non_discriminating_without_a_"
            "new_energy_estimate"
        ),
        "independent_continuum_derivative_bound_route": (
            "deferred_because_the_current_GR0_health_contract_does_not_bound_"
            "every_required_state_derivative_between_steps"
        ),
        "additional_finite_sampling_route": "forbidden_by_TDG4",
        "endpoint_only_dense_output_route": (
            "rejected_by_exact_interior_alias_control"
        ),
        "selected_route_is_purpose_built_on_existing_RK_stage_and_transaction_records": True,
        "no_new_runtime_dependency": True,
    }


def _expected_successor() -> dict[str, Any]:
    return {
        "TDG4_sampling_identifiability_theorem_completed": True,
        "current_temporal_spectral_rule_retired_as_continuum_admission": True,
        "current_temporal_spectral_rule_retained_as_sampled_diagnostic": True,
        "replacement_temporal_gate_design_authorized": True,
        "TDG5_replacement_design_frozen": True,
        "TDG5_exact_no_history_controls_pass": True,
        "TDG5_stage_complete_refinement_theorem_completed": False,
        "runtime_refinement_pair_implementation_authorized": False,
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


PROOF_KEYS = {
    "source_commit_and_every_bound_blob_must_match",
    "TDG4_theorem_and_narrow_retirement_must_remain_complete",
    "no_raw_history_checkpoint_or_candidate_state_may_be_loaded",
    "Hermite_sampling_and_inverse_matrices_must_multiply_to_identity_exactly",
    "Hermite_inverse_condition_upper_bound_must_equal_54",
    "endpoint_only_alias_must_be_detected_by_endpoint_slopes",
    "RK4_order_conditions_through_order_four_must_pass_exactly",
    "SSPRK3_order_conditions_through_order_three_must_pass_exactly",
    "same_grid_step_doubling_scalar_controls_must_contract_exactly",
    "lineage_method_map_transaction_nonclaim_and_claim_mutations_must_fail_closed",
    "canonical_hash_bound_result_required",
}


def validate_config_data(config: Mapping[str, Any]) -> None:
    _strict_keys(
        "TDG5 freeze config",
        config,
        {
            "schema_version",
            "artifact_id",
            "project_version",
            "metric_signature",
            "riemann_convention",
            "source_result",
            "source_config",
            "source_reproducer",
            "source_theorem_module",
            "source_document",
            "numerical_engine",
            "PROTO13_runner",
            "protocol_v13",
            "CAL10_config",
            "CAL10_result",
            "scope",
            "immutable_lineage",
            "frozen_replacement_design",
            "alternatives",
            "proof_contract",
            "successor_boundary",
            "claims",
        },
    )
    expected_paths = {
        "source_result": _rel(SOURCE_RESULT),
        "source_config": _rel(SOURCE_CONFIG),
        "source_reproducer": _rel(SOURCE_REPRODUCER),
        "source_theorem_module": _rel(SOURCE_THEOREM_MODULE),
        "source_document": _rel(SOURCE_DOCUMENT),
        "numerical_engine": _rel(NUMERICAL_ENGINE),
        "PROTO13_runner": _rel(PROTO13_RUNNER),
        "protocol_v13": _rel(PROTOCOL_V13),
        "CAL10_config": _rel(CAL10_CONFIG),
        "CAL10_result": _rel(CAL10_RESULT),
    }
    if (
        config.get("schema_version") != 1
        or config.get("artifact_id") != ARTIFACT_ID
        or config.get("project_version") != PROJECT_VERSION
        or config.get("metric_signature") != "-+++"
        or config.get("riemann_convention") != "plus_partial_mu_gamma_nu"
        or any(config.get(key) != value for key, value in expected_paths.items())
    ):
        raise ValueError("TDG5 freeze identity differs")
    if config.get("scope") != _expected_scope():
        raise ValueError("TDG5 freeze scope differs")
    if config.get("immutable_lineage") != _expected_lineage():
        raise ValueError("TDG5 freeze lineage differs")
    if config.get("frozen_replacement_design") != _expected_frozen_design():
        raise ValueError("TDG5 frozen replacement design differs")
    if config.get("alternatives") != _expected_alternatives():
        raise ValueError("TDG5 alternative ledger differs")
    proof = config.get("proof_contract", {})
    if set(proof) != PROOF_KEYS or set(proof.values()) != {True}:
        raise ValueError("TDG5 proof contract differs")
    if config.get("successor_boundary") != _expected_successor():
        raise ValueError("TDG5 successor boundary differs")
    if config.get("claims") != EXPECTED_CLAIMS:
        raise ValueError("TDG5 claims differ")


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    with path.open("rb") as handle:
        config = tomllib.load(handle)
    validate_config_data(config)
    return config


def _load_toml(path: Path) -> dict[str, Any]:
    with path.open("rb") as handle:
        return tomllib.load(handle)


def _verify_source_lineage(config: Mapping[str, Any]) -> dict[str, Any]:
    lineage = config["immutable_lineage"]
    commit = lineage["checkpoint_commit"]
    _git("cat-file", "-e", f"{commit}^{{commit}}")
    if _git("merge-base", "--is-ancestor", commit, "HEAD", check=False).returncode:
        raise ValueError("TDG5 source commit is not an ancestor of HEAD")
    tracked = {
        _rel(SOURCE_CONFIG): lineage["source_config_sha256"],
        _rel(SOURCE_RESULT): lineage["source_result_sha256"],
        _rel(SOURCE_REPRODUCER): lineage["source_reproducer_sha256"],
        _rel(SOURCE_THEOREM_MODULE): lineage["source_theorem_module_sha256"],
        _rel(SOURCE_DOCUMENT): lineage["source_document_sha256"],
        _rel(NUMERICAL_ENGINE): lineage["numerical_engine_sha256"],
        _rel(PROTO13_RUNNER): lineage["PROTO13_runner_sha256"],
        _rel(PROTOCOL_V13): lineage["protocol_v13_sha256"],
        _rel(CAL10_CONFIG): lineage["CAL10_config_sha256"],
        _rel(CAL10_RESULT): lineage["CAL10_result_sha256"],
    }
    for relative, expected_hash in tracked.items():
        committed = _git("show", f"{commit}:{relative}").stdout
        current = (REPOSITORY / relative).read_bytes()
        if (
            committed != current
            or _bytes_sha(committed) != expected_hash
            or _bytes_sha(current) != expected_hash
        ):
            raise ValueError(f"TDG5 source tracked blob differs: {relative}")

    source_raw = SOURCE_RESULT.read_bytes()
    source = json.loads(source_raw)
    status = source.get("gate_status", {})
    boundary = source.get("artifact_payload", {}).get("claim_boundary", {})
    if (
        source_raw != _canonical_bytes(source)
        or source.get("artifact_id") != "FGC-1-TDG4-PREF19"
        or source.get("classification")
        != "completed_TDG4_sampling_identifiability_theorem_current_sampled_rule_retired_as_continuum_admission"
        or status.get("TDG4_sampling_identifiability_theorem_completed")
        is not True
        or status.get(
            "finite_samples_alone_proved_insufficient_for_continuum_top_band_bound"
        )
        is not True
        or status.get(
            "current_temporal_spectral_rule_retired_as_continuum_admission"
        )
        is not True
        or status.get(
            "current_temporal_spectral_rule_retained_as_sampled_diagnostic"
        )
        is not True
        or status.get("replacement_temporal_gate_design_authorized") is not True
        or status.get("replacement_temporal_admission_defined") is not False
        or status.get("PROTO14_frozen") is not False
        or status.get("GR0_case_eligible") is not False
        or status.get("SGBL_execution_authorized") is not False
        or status.get("FGCQR_holdout_execution_authorized") is not False
        or boundary.get("actual_terminal_history_arrays_consumed") is not False
        or boundary.get("campaign_checkpoint_loaded") is not False
        or boundary.get("state_advanced") is not False
        or boundary.get("historical_PROTO13_stop_preserved") is not True
        or boundary.get("historical_CAL10_result_reclassified") is not False
        or boundary.get("replacement_temporal_admission_defined") is not False
        or boundary.get("new_trajectory_authorized") is not False
    ):
        raise ValueError("TDG5 source PREF19 compact result differs")

    protocol = _load_toml(PROTOCOL_V13)
    if (
        protocol.get("artifact_id") != "FGC-2-SF1-PROTO13"
        or protocol.get("amendment", {}).get(
            "absolute_spectral_budgets_and_direct_ratio_ceiling_unchanged"
        )
        is not True
    ):
        raise ValueError("TDG5 PROTO13 historical contract differs")

    cal10_raw = CAL10_RESULT.read_bytes()
    cal10 = json.loads(cal10_raw)
    if (
        cal10_raw != _canonical_bytes(cal10)
        or cal10.get("artifact_id") != "FGC-1-CAL10-PREF15"
        or cal10.get("gate_status", {}).get("PROTO13_temporal_admission_failed")
        is not True
        or cal10.get("gate_status", {}).get(
            "PROTO13_temporal_failure_cause_derived"
        )
        is not False
        or cal10.get("artifact_payload", {})
        .get("claim_boundary", {})
        .get("temporal_stop_is_not_a_candidate_action_result")
        is not True
    ):
        raise ValueError("TDG5 CAL10 historical boundary differs")

    return {
        "checkpoint_commit": commit,
        "checkpoint_commit_is_ancestor_of_HEAD": True,
        "tracked_blob_sha256": tracked,
        "source_compact_result_verified": True,
        "source_compact_result_read": True,
        "historical_PROTO13_stop_preserved": True,
        "historical_CAL10_result_reclassified": False,
        "actual_terminal_history_arrays_consumed": False,
        "campaign_checkpoint_loaded": False,
        "campaign_checkpoint_resumed": False,
        "state_advanced": False,
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
    invalid_method_rejected = False
    inexact_data_rejected = False
    invalid_interval_rejected = False
    try:
        butcher_order_certificate("Euler")
    except (TypeError, ValueError):
        invalid_method_rejected = True
    try:
        hermite_coefficients((0, 1, 2, 3.0))
    except (TypeError, ValueError):
        inexact_data_rejected = True
    try:
        evaluate_cubic((0, 1, 2, 3), Fraction(3, 2))
    except (TypeError, ValueError):
        invalid_interval_rejected = True

    attacks: dict[str, Callable[[dict[str, Any]], None]] = {
        "source_hash_mutation_rejected": lambda value: value[
            "immutable_lineage"
        ].__setitem__("source_result_sha256", "0" * 64),
        "selected_route_mutation_rejected": lambda value: value[
            "frozen_replacement_design"
        ].__setitem__("selected_sufficient_route", "additional_finite_sampling"),
        "Hermite_determinant_mutation_rejected": lambda value: value[
            "frozen_replacement_design"
        ]["Hermite_continuous_extension"].__setitem__(
            "sampling_matrix_determinant", "0"
        ),
        "Hermite_condition_mutation_rejected": lambda value: value[
            "frozen_replacement_design"
        ]["Hermite_continuous_extension"].__setitem__(
            "condition_number_infinity_upper_bound", "53"
        ),
        "owned_state_mutation_rejected": lambda value: value[
            "frozen_replacement_design"
        ].__setitem__("monitored_state", "64_normal_flow_tracers"),
        "method_order_mutation_rejected": lambda value: value[
            "frozen_replacement_design"
        ]["same_grid_step_refinement"].__setitem__("SSPRK3_formal_order", 2),
        "Richardson_denominator_mutation_rejected": lambda value: value[
            "frozen_replacement_design"
        ]["same_grid_step_refinement"].__setitem__("RK4_Richardson_denominator", 16),
        "transaction_ownership_mutation_rejected": lambda value: value[
            "frozen_replacement_design"
        ]["transaction_ownership"].__setitem__(
            "fine_two_half_steps_commit_atomically_or_neither_commits", False
        ),
        "margin_normalization_mutation_rejected": lambda value: value[
            "frozen_replacement_design"
        ]["error_semantics"].__setitem__(
            "debit_may_not_be_normalized_by_the_observed_trapped_or_Raychaudhuri_margin",
            False,
        ),
        "raw_history_scope_promotion_rejected": lambda value: value[
            "scope"
        ].__setitem__("actual_terminal_history_arrays_consumed", True),
        "theorem_completion_promotion_rejected": lambda value: value[
            "claims"
        ].__setitem__("TDG5_stage_complete_refinement_theorem_completed", True),
        "runtime_implementation_authorization_promotion_rejected": lambda value: value[
            "claims"
        ].__setitem__("runtime_refinement_pair_implementation_authorized", True),
        "replacement_admission_promotion_rejected": lambda value: value[
            "claims"
        ].__setitem__("replacement_temporal_admission_defined", True),
        "PROTO14_promotion_rejected": lambda value: value["claims"].__setitem__(
            "PROTO14_frozen", True
        ),
        "GR0_eligibility_promotion_rejected": lambda value: value[
            "claims"
        ].__setitem__("GR0_case_eligible", True),
        "candidate_promotion_rejected": lambda value: value["claims"].__setitem__(
            "FGCQR_holdout_execution_authorized", True
        ),
        "physical_claim_promotion_rejected": lambda value: value[
            "claims"
        ].__setitem__("physical_transition_claim_authorized", True),
    }
    controls = {
        "invalid_method_mutation_rejected": invalid_method_rejected,
        "inexact_Hermite_data_mutation_rejected": inexact_data_rejected,
        "invalid_Hermite_interval_mutation_rejected": invalid_interval_rejected,
        **{
            name: _expect_config_rejection(config, attack)
            for name, attack in attacks.items()
        },
    }
    if set(controls.values()) != {True}:
        raise ValueError(f"TDG5 mutation control failed: {controls}")
    return controls


def _validate_preflight(preflight: Mapping[str, Any]) -> None:
    hermite = preflight.get("Hermite_data_map", {})
    alias = preflight.get("endpoint_only_alias_control", {})
    methods = preflight.get("method_controls", {})
    if (
        preflight.get("selected_sufficient_route")
        != TDG5_SELECTED_SUFFICIENT_ROUTE
        or preflight.get("owned_row_policy") != TDG5_OWNED_ROW_POLICY
        or tuple(preflight.get("required_runtime_layers", ()))
        != TDG5_REQUIRED_RUNTIME_LAYERS
        or preflight.get("all_exact_controls_pass") is not True
        or hermite.get("data_order") != list(TDG5_HERMITE_DATA_ORDER)
        or hermite.get("determinant") != "1"
        or hermite.get("left_inverse_exact") is not True
        or hermite.get("right_inverse_exact") is not True
        or hermite.get("condition_number_infinity_upper_bound") != "54"
        or alias.get("endpoint_values_equal") is not True
        or alias.get("adversarial_midpoint_value") != "1"
        or alias.get("stage_complete_Hermite_data_detects_it") is not True
        or set(methods) != {method for method, _ in TDG5_METHOD_ORDERS}
        or preflight.get("actual_terminal_histories_consumed") is not False
        or preflight.get("campaign_checkpoint_loaded") is not False
        or preflight.get("campaign_checkpoint_resumed") is not False
        or preflight.get("state_advanced") is not False
        or preflight.get("runtime_refinement_pair_implemented") is not False
        or preflight.get("runtime_thresholds_frozen") is not False
        or preflight.get("replacement_temporal_admission_defined") is not False
        or preflight.get("PROTO14_frozen") is not False
        or preflight.get("GR0_case_eligible") is not False
    ):
        raise ValueError("TDG5 exact preflight differs")
    for method, formal_order in TDG5_METHOD_ORDERS:
        order = methods[method].get("order", {})
        doubling = methods[method].get("step_doubling", {})
        if (
            order.get("formal_order") != formal_order
            or order.get("all_required_order_conditions_pass") is not True
            or doubling.get("formal_order") != formal_order
            or doubling.get("Richardson_denominator") != 2**formal_order - 1
            or doubling.get(
                "contraction_is_stricter_than_two_to_minus_formal_order"
            )
            is not True
            or doubling.get("Richardson_factor_is_not_a_rigorous_PDE_error_bound")
            is not True
        ):
            raise ValueError(f"TDG5 {method} exact control differs")


def record(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    config = load_config(config_path)
    lineage = _verify_source_lineage(config)
    preflight = stage_complete_temporal_refinement_preflight()
    _validate_preflight(preflight)
    mutations = _mutation_controls(config)
    result = {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "metric_signature": "-+++",
        "riemann_convention": "plus_partial_mu_gamma_nu",
        "classification": (
            "prospective_TDG5_stage_complete_temporal_refinement_design_freeze"
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
            "frozen_replacement_design": deepcopy(
                config["frozen_replacement_design"]
            ),
            "alternatives": deepcopy(config["alternatives"]),
            "exact_no_history_preflight": preflight,
            "mutation_controls": mutations,
            "proof_contract": dict(config["proof_contract"]),
            "successor_boundary": dict(config["successor_boundary"]),
            "claim_boundary": {
                "actual_terminal_history_arrays_consumed": False,
                "campaign_checkpoint_loaded": False,
                "state_advanced": False,
                "historical_PROTO13_stop_preserved": True,
                "historical_CAL10_result_reclassified": False,
                "replacement_design_frozen": True,
                "stage_complete_refinement_theorem_completed": False,
                "runtime_refinement_pair_implementation_authorized": False,
                "runtime_refinement_pair_implemented": False,
                "runtime_thresholds_frozen": False,
                "replacement_temporal_admission_defined": False,
                "new_trajectory_authorized": False,
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
    expected_nonclaims = {
        key: False
        for key, item in sorted(EXPECTED_CLAIMS.items())
        if item is False
    }
    expected_boundary = {
        "actual_terminal_history_arrays_consumed": False,
        "campaign_checkpoint_loaded": False,
        "state_advanced": False,
        "historical_PROTO13_stop_preserved": True,
        "historical_CAL10_result_reclassified": False,
        "replacement_design_frozen": True,
        "stage_complete_refinement_theorem_completed": False,
        "runtime_refinement_pair_implementation_authorized": False,
        "runtime_refinement_pair_implemented": False,
        "runtime_thresholds_frozen": False,
        "replacement_temporal_admission_defined": False,
        "new_trajectory_authorized": False,
        "GR0_case_eligible": False,
        "candidate_branch_opened": False,
        "physical_question_answered": False,
    }
    if (
        value.get("schema_version") != 1
        or value.get("artifact_id") != ARTIFACT_ID
        or value.get("project_version") != PROJECT_VERSION
        or value.get("classification")
        != "prospective_TDG5_stage_complete_temporal_refinement_design_freeze"
        or value.get("generated_by") != "scripts/reproduce_fgc_tdg5_frz1.py"
        or value.get("gate_status") != EXPECTED_CLAIMS
        or value.get("nonclaims") != expected_nonclaims
        or value.get("source_config_sha256")
        != {_rel(DEFAULT_CONFIG): _sha(DEFAULT_CONFIG)}
        or payload.get("frozen_replacement_design")
        != config["frozen_replacement_design"]
        or payload.get("alternatives") != config["alternatives"]
        or payload.get("successor_boundary") != config["successor_boundary"]
        or payload.get("claim_boundary") != expected_boundary
    ):
        raise ValueError("TDG5 stored result boundary differs")
    implementation = value.get("implementation_sha256", {})
    if (
        not isinstance(implementation, Mapping)
        or set(implementation) != {_rel(path) for path in IMPLEMENTATION}
        or any(
            _sha(REPOSITORY / relative) != expected
            for relative, expected in implementation.items()
        )
    ):
        raise ValueError("TDG5 implementation hash ledger differs")
    document = REPOSITORY / str(value.get("derivation_document", ""))
    if (
        value.get("derivation_document") != "docs/fgc-tdg5-frz1.md"
        or not document.is_file()
        or _sha(document) != value.get("derivation_document_sha256")
    ):
        raise ValueError("TDG5 derivation document differs")
    preflight = payload.get("exact_no_history_preflight", {})
    _validate_preflight(preflight)
    mutations = payload.get("mutation_controls", {})
    if not isinstance(mutations, Mapping) or set(mutations.values()) != {True}:
        raise ValueError("TDG5 mutation controls differ")


def verify_canonical(
    config_path: Path = DEFAULT_CONFIG,
    output_path: Path = DEFAULT_OUTPUT,
) -> dict[str, Any]:
    config = load_config(config_path)
    raw = output_path.read_bytes()
    stored = json.loads(raw)
    if not isinstance(stored, dict) or raw != _canonical_bytes(stored):
        raise ValueError("TDG5 stored result is not canonical JSON")
    _validate_stored_record(stored, config)
    observed = record(config_path)
    if stored != observed:
        raise ValueError("TDG5 stored result differs from reproduction")
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
            f"PASS {ARTIFACT_ID}: exact_preflight="
            f"{value['artifact_payload']['exact_no_history_preflight']['all_exact_controls_pass']}"
        )
    else:
        print(_canonical(record(args.config)), end="")


if __name__ == "__main__":
    main()
