#!/usr/bin/env python3
"""Construct or verify the prospective, no-trajectory TDG7 lattice freeze."""
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
from typing import Any, Mapping

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.tdg7_binary64_subdivision_lattice import (  # noqa: E402
    CoordinateLatticeLimitReached,
    cal11_event24_witness,
    plan_forward_proto14_subdivision,
    reconstruct_historical_tdg6_guard,
    validate_tdg7_binary64_subdivision_plan,
)

ARTIFACT_ID = "FGC-1-TDG7-FRZ1"
CONFIG = REPOSITORY / "configs/fgc/fgc-1-tdg7-frz1.toml"
OUTPUT = REPOSITORY / "results/fgc-1-tdg7-frz1.json"
DOCUMENT = REPOSITORY / "docs/fgc-tdg7-frz1.md"
DESIGN = REPOSITORY / "src/recursive_horizons/fgc/evolution/tdg7_binary64_subdivision_lattice.py"
FROZEN_CONFIG_SHA256 = "28a05df475e93c39805627cd05208c0078fbf0071a8ff90526a3dbc696c4619d"


def canonical(value: object) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True, ensure_ascii=True,
                       allow_nan=False) + "\n").encode("utf-8")


def sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _git(*args: str, check: bool = True) -> subprocess.CompletedProcess[bytes]:
    return subprocess.run(["git", *args], cwd=REPOSITORY, check=check,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE)


def _blob(commit: str, relative: str) -> bytes:
    return _git("show", f"{commit}:{relative}").stdout


def _strict_keys(name: str, value: Mapping[str, Any], expected: set[str]) -> None:
    if set(value) != expected:
        raise ValueError(f"TDG7 {name} keys differ")


def _hex(value: str) -> float:
    if not isinstance(value, str) or not value.startswith("0x"):
        raise ValueError("TDG7 binary64 literal differs")
    return float.fromhex(value)


def validate_config_data(config: Mapping[str, Any]) -> None:
    _strict_keys("configuration", config, {
        "schema_version", "artifact_id", "project_version", "metric_signature",
        "riemann_convention", "source_result", "source_config", "source_reproducer",
        "source_document", "source_diagnosis_module", "historical_TDG6_runtime",
        "historical_PROTO14_runner", "historical_numerical_engine", "historical_run_config",
        "runtime_authorization_result", "source_checkpoint", "design_module", "scope",
        "immutable_lineage", "historical_event24_witness", "selected_route",
        "typed_stop", "exact_controls", "alternatives", "proof_contract",
        "successor_boundary", "claims",
    })
    if any(config.get(key) != value for key, value in {
        "schema_version": 1, "artifact_id": ARTIFACT_ID, "project_version": "0.11.0",
        "metric_signature": "-+++", "riemann_convention": "plus_partial_mu_gamma_nu",
        "source_result": "results/fgc-1-cal11-pref22.json",
        "source_config": "configs/fgc/fgc-1-cal11-pref22.toml",
        "source_reproducer": "scripts/reproduce_fgc_cal11_pref22.py",
        "source_document": "docs/fgc-cal11-pref22.md",
        "source_diagnosis_module": "src/recursive_horizons/fgc/evolution/cal11_proto14_campaign_diagnosis.py",
        "historical_TDG6_runtime": "src/recursive_horizons/fgc/evolution/tdg6_temporal_admission_runtime.py",
        "historical_PROTO14_runner": "scripts/run_fgc_gr0_calibration_v14.py",
        "historical_numerical_engine": "src/recursive_horizons/fgc/evolution/numerical_engine.py",
        "historical_run_config": "configs/fgc/fgc-1-cal9-run1.toml",
        "runtime_authorization_result": "results/fgc-1-hlt12-mon12.json",
        "source_checkpoint": "runs/fgc-2-sf1/proto14/calibration/latest-checkpoint.npz",
        "design_module": "src/recursive_horizons/fgc/evolution/tdg7_binary64_subdivision_lattice.py",
    }.items()):
        raise ValueError("TDG7 top-level contract differs")
    scope = config["scope"]
    _strict_keys("scope", scope, {
        "target", "role", "calibration_branch", "source_compact_result_read",
        "historical_source_blobs_read_from_immutable_commits",
        "optional_checkpoint_metadata_may_be_verified_read_only",
        "campaign_physical_state_arrays_consumed", "campaign_checkpoint_resumed",
        "production_state_advanced", "output_namespace_created",
        "synthetic_or_exact_arithmetic_controls_only", "SGBL_trajectory_read",
        "FGCQR_trajectory_read", "physical_or_candidate_question_answered",
    })
    if scope != {
        "target": "FGC-2-SF1-TDG7", "role": "prospective_exact_binary64_subdivision_lattice_diagnosis_and_design_freeze",
        "calibration_branch": "GR-0", "source_compact_result_read": True,
        "historical_source_blobs_read_from_immutable_commits": True,
        "optional_checkpoint_metadata_may_be_verified_read_only": True,
        "campaign_physical_state_arrays_consumed": False, "campaign_checkpoint_resumed": False,
        "production_state_advanced": False, "output_namespace_created": False,
        "synthetic_or_exact_arithmetic_controls_only": True, "SGBL_trajectory_read": False,
        "FGCQR_trajectory_read": False, "physical_or_candidate_question_answered": False,
    }:
        raise ValueError("TDG7 scope differs")
    lineage = config["immutable_lineage"]
    _strict_keys("immutable lineage", lineage, {"checkpoint_commit", "historical_authorization_commit", "source_config_sha256", "source_result_sha256", "source_reproducer_sha256", "source_document_sha256", "source_diagnosis_module_sha256", "historical_TDG6_runtime_sha256", "historical_PROTO14_runner_sha256", "historical_numerical_engine_sha256", "historical_run_config_sha256", "runtime_authorization_result_sha256", "source_checkpoint_sha256", "checkpoint_commit_must_be_ancestor_of_HEAD", "historical_authorization_commit_must_be_ancestor_of_checkpoint_commit", "CAL11_blobs_must_match_checkpoint_commit_and_worktree", "historical_TDG6_and_PROTO14_blobs_must_match_authorization_commit", "historical_numerical_engine_must_match_authorization_commit_and_worktree", "CAL11_result_must_be_canonical_and_authorize_TDG7_design_only", "source_checkpoint_is_read_only", "source_checkpoint_may_be_absent_in_clean_clone", "source_checkpoint_may_not_resume", "only_metadata_utf8_may_be_read_when_checkpoint_is_present"})
    required_hashes = {
        "source_config_sha256", "source_result_sha256", "source_reproducer_sha256",
        "source_document_sha256", "source_diagnosis_module_sha256",
        "historical_TDG6_runtime_sha256", "historical_PROTO14_runner_sha256", "historical_numerical_engine_sha256",
        "historical_run_config_sha256", "runtime_authorization_result_sha256",
        "source_checkpoint_sha256",
    }
    expected_hash_values = {
        "source_config_sha256": "685ee05442d7c92ce83994d923890f4fa208b10ffbf2946026b07e9a9cf4d34d",
        "source_result_sha256": "1339792198ec78ad83bc3e11a2805f09f961607ab6c1f6f06e79bb88ba6b58b8",
        "source_reproducer_sha256": "c02deeea237176706310dd82a81a71dac85617c1130f2e78342211cbd89e7ef4",
        "source_document_sha256": "b115febf6aeb7fdce6d604bf5311ca282d81a3beffaad4765ac75e9d1234d606",
        "source_diagnosis_module_sha256": "492a520325715743110d81d1b366070ac307fbe9cde9e83f1c4465746021b5bc",
        "historical_TDG6_runtime_sha256": "039de0bd008a2536bea6f7187afc3856bfc841d14790cbae5ca2fe27f7eb8e2d",
        "historical_PROTO14_runner_sha256": "4c13eaf86cf31f55e28b7f60661db78ad22e1bb5b7460adf7aeff455efdc8e45",
        "historical_numerical_engine_sha256": "8ea99604a85b4e1acbf869ed2462ea4aa8301c06e068330a1d0ce66a51a0bebf",
        "historical_run_config_sha256": "8e6fafff638ffd49549efab8b6f969a01964c915116db9ce3698953201172176",
        "runtime_authorization_result_sha256": "7f4068a3116ff1838ce7f0165ae18e72d2eea68933b87f056b26d8149c935731",
        "source_checkpoint_sha256": "e7dc31a25a775a9377fbc45d1b239deb492bdc9ebbef85f94913513f885b015d",
    }
    if (lineage.get("checkpoint_commit") != "b59d214832a1d3429f624579a93c6b9e1f544cf5"
            or lineage.get("historical_authorization_commit") != "7534a1662d8078a63c98025ecd464bb068fa012f"
            or not all(isinstance(lineage.get(key), str) and len(lineage[key]) == 64 for key in required_hashes)
            or any(lineage.get(key) != value for key, value in expected_hash_values.items())
            or any(lineage.get(key) is not True for key in {
                "checkpoint_commit_must_be_ancestor_of_HEAD", "historical_authorization_commit_must_be_ancestor_of_checkpoint_commit",
                "CAL11_blobs_must_match_checkpoint_commit_and_worktree", "historical_TDG6_and_PROTO14_blobs_must_match_authorization_commit",
                "historical_numerical_engine_must_match_authorization_commit_and_worktree",
                "CAL11_result_must_be_canonical_and_authorize_TDG7_design_only", "source_checkpoint_is_read_only",
                "source_checkpoint_may_be_absent_in_clean_clone", "source_checkpoint_may_not_resume",
                "only_metadata_utf8_may_be_read_when_checkpoint_is_present",
            })):
        raise ValueError("TDG7 immutable lineage differs")
    witness = config["historical_event24_witness"]
    _strict_keys("historical witness", witness, {
        "member", "completed_common_event_index", "target_common_event_index", "current_time_hex", "target_time_hex",
        "cfl_maximum_hex", "outer_radius", "point_count", "grid_spacing_hex", "previous_speed_upper_hex",
        "cfl_numerator_hex", "target_remaining_hex", "requested_macro_width_hex", "historical_expected_macro_endpoint_hex", "historical_count_one_step_hex",
        "historical_count_two_step_hex", "historical_count_four_step_hex", "rounded_count_one_difference_hex",
        "rounded_count_two_difference_hex", "rounded_count_four_difference_hex", "event_ulp_hex", "macro_quantum_hex",
        "aligned_tick_count", "aligned_macro_width_hex", "aligned_fine_width_hex", "selected_aligned_macro_endpoint_hex", "conservative_reduction_hex",
        "conservative_reduction_exact", "historical_guard_fails_for_counts",
        "historical_rounded_boundaries_equal_selected_aligned_boundaries", "selected_boundary_hex",
        "selected_fine_midpoint_hex", "selected_fine_width_over_event_quantum", "failure_phase",
    })
    if witness.get("historical_guard_fails_for_counts") != [1, 2, 4] or witness.get("member") != "RK4-2049":
        raise ValueError("TDG7 witness differs")
    if any(witness.get(key) != value for key, value in {
        "current_time_hex": "0x1.7000000000000p+0", "target_time_hex": "0x1.8000000000000p+0",
        "requested_macro_width_hex": "0x1.aaa90b0fb5c26p-8", "macro_quantum_hex": "0x1.0000000000000p-49",
        "historical_expected_macro_endpoint_hex": "0x1.71aaa90b0fb5cp+0", "selected_aligned_macro_endpoint_hex": "0x1.71aaa90b0fb58p+0",
        "aligned_macro_width_hex": "0x1.aaa90b0fb5800p-8", "aligned_fine_width_hex": "0x1.aaa90b0fb5800p-10",
        "conservative_reduction_exact": "531/576460752303423488",
    }.items()):
        raise ValueError("TDG7 witness literal differs")
    for key in ("current_time_hex", "target_time_hex", "cfl_maximum_hex", "grid_spacing_hex",
                "previous_speed_upper_hex", "cfl_numerator_hex", "target_remaining_hex",
                "requested_macro_width_hex", "historical_count_one_step_hex", "historical_count_two_step_hex",
                "historical_count_four_step_hex", "rounded_count_one_difference_hex",
                "rounded_count_two_difference_hex", "rounded_count_four_difference_hex", "event_ulp_hex",
        "macro_quantum_hex", "aligned_macro_width_hex", "aligned_fine_width_hex", "conservative_reduction_hex"):
        _hex(witness[key])
    if (witness.get("historical_rounded_boundaries_equal_selected_aligned_boundaries") is not False
            or witness.get("aligned_tick_count") != 3664984285035
            or witness.get("selected_fine_width_over_event_quantum") != 7329968570070
            or witness.get("selected_boundary_hex") != ["0x1.7000000000000p+0", "0x1.706aaa42c3ed6p+0", "0x1.70d5548587dacp+0", "0x1.713ffec84bc82p+0", "0x1.71aaa90b0fb58p+0"]
            or witness.get("selected_fine_midpoint_hex") != ["0x1.7035552161f6bp+0", "0x1.709fff6425e41p+0", "0x1.710aa9a6e9d17p+0", "0x1.717553e9adbedp+0"]):
        raise ValueError("TDG7 stage-safe witness differs")
    route = config["selected_route"]
    _strict_keys("selected route", route, {
        "name", "time_direction", "frozen_time_minimum", "frozen_time_maximum", "event_output_interval",
        "requested_CFL_retry_or_target_width_is_a_conservative_upper_bound", "event_quantum", "macro_quantum",
        "scheduled_width", "fine_width", "fine_width_is_an_even_multiple_of_event_quantum", "minimum_aligned_width", "shared_five_boundaries", "outer_boundaries_select",
        "medium_boundaries_select", "fine_boundaries_select", "all_selection_and_representability_checks_use_exact_dyadic_arithmetic",
        "conversion_to_binary64_occurs_once_after_exact_selection", "scheduled_width_may_not_exceed_requested_cap_or_target_remaining",
        "rounding_up_is_forbidden", "target_reachability_on_event_lattice_is_required",
        "target_limited_final_step_must_land_bitwise_on_target", "tolerance_based_target_normalization_is_forbidden",
        "same_five_boundaries_must_feed_all_three_shadow_paths", "RK4_and_SSPRK3_stage_abscissae", "all_stage_times_are_exact_binary64_shared_lattice_coordinates", "positive_substep_midpoint_must_be_strictly_interior", "alignment_reduction_is_serialized_numerical_provenance_not_physical_error",
        "physical_equations_sources_projector_spatial_operator_and_TDG6_threshold_unchanged",
    })
    if route.get("name") != "event_target_owned_exact_stage_safe_binary64_lattice" or route.get("macro_quantum") != "8*event_quantum" or route.get("rounding_up_is_forbidden") is not True or route.get("fine_width") != "scheduled_width/4" or route.get("fine_width_is_an_even_multiple_of_event_quantum") is not True or route.get("all_stage_times_are_exact_binary64_shared_lattice_coordinates") is not True or route.get("positive_substep_midpoint_must_be_strictly_interior") is not True or route.get("minimum_aligned_width") != "optional_positive_binary64_lower_bound_checked_after_conservative_alignment":
        raise ValueError("TDG7 route differs")
    typed = config["typed_stop"]
    _strict_keys("typed stop", typed, {"classification", "nonfinite_coordinate", "target_not_ahead", "outside_frozen_time_envelope", "nonpositive_requested_cap", "endpoint_not_q_aligned", "target_not_stage_lattice_aligned", "no_positive_aligned_width", "below_minimum_aligned_width", "internal_exactness_failure", "typed_stop_occurs_before_any_shadow_or_real_state_operation", "typed_stop_is_invalid_numerical_runtime_not_physics"})
    if typed.get("classification") != "coordinate_lattice_limit_reached" or any(typed.get(key) != key for key in {"nonfinite_coordinate", "target_not_ahead", "outside_frozen_time_envelope", "nonpositive_requested_cap", "endpoint_not_q_aligned", "target_not_stage_lattice_aligned", "no_positive_aligned_width", "below_minimum_aligned_width", "internal_exactness_failure"}) or typed.get("typed_stop_occurs_before_any_shadow_or_real_state_operation") is not True or typed.get("typed_stop_is_invalid_numerical_runtime_not_physics") is not True:
        raise ValueError("TDG7 typed-stop contract differs")
    controls = config["exact_controls"]
    _strict_keys("exact controls", controls, {"historical_event24_old_guard_fails_for_one_two_and_four", "historical_event24_aligned_plan_passes_all_exact_postconditions", "historical_event24_stage_times_are_exact_and_strictly_interior", "historical_event24_cap_is_never_exceeded", "historical_event24_event_target_is_unchanged", "exact_target_limited_step_lands_bitwise", "multi_step_event_preserves_lattice_and_lands_bitwise", "upward_binade_crossing_uses_larger_event_ulp", "start_only_or_minimum_ulp_mutations_fail_crossing_control", "too_small_positive_cap_returns_typed_stop", "minimum_aligned_width_returns_typed_stop_without_rounding_up", "nonfinite_and_nonpositive_caps_return_typed_stops", "unaligned_target_returns_typed_stop", "round_up_mutation_exceeds_cap", "independently_generated_path_endpoint_mutation_is_detected", "one_endpoint_nextafter_mutation_is_detected", "no_zero_width_or_midpoint_endpoint_alias_is_permitted"})
    if not controls or any(value is not True for value in controls.values()):
        raise ValueError("TDG7 exact controls differ")
    alternatives = config["alternatives"]
    _strict_keys("alternatives", alternatives, {"remove_or_relax_historical_guard_selected", "remove_or_relax_historical_guard_rejection", "epsilon_or_approximate_uniformity_selected", "epsilon_or_approximate_uniformity_rejection", "reinterpret_rounded_nonuniform_paths_after_execution_selected", "reinterpret_rounded_nonuniform_paths_rejection", "direct_PROTO14_patch_or_resume_selected", "direct_PROTO14_patch_or_resume_rejection", "exact_conservative_shared_lattice_selected", "four_q_endpoint_only_lattice_selected", "four_q_endpoint_only_lattice_rejection"})
    if alternatives.get("exact_conservative_shared_lattice_selected") is not True or any(alternatives.get(key) is not False for key in alternatives if key.endswith("_selected") and key != "exact_conservative_shared_lattice_selected"):
        raise ValueError("TDG7 alternatives differ")
    proof = config["proof_contract"]
    _strict_keys("proof contract", proof, {"CAL11_result_config_reproducer_document_and_diagnosis_must_match_checkpoint_commit", "historical_TDG6_runtime_and_PROTO14_runner_must_match_authorization_commit", "historical_numerical_engine_must_be_inspected_and_pinned_from_authorization_commit", "TDG7_binder_must_inspect_and_pin_RK4_and_SSPRK3_start_start_plus_dt_over_2_and_start_plus_dt", "optional_checkpoint_whole_file_hash_must_match_before_metadata_read", "optional_checkpoint_previous_speed_scalar_must_match_without_loading_physical_arrays", "historical_requested_width_must_be_derived_in_the_frozen_operation_order", "historical_one_two_four_guard_failures_must_be_reproduced_bitwise", "aligned_plan_must_be_derived_by_an_independent_exact_dyadic_route", "every_successful_plan_must_pass_exact_and_bitwise_partition_postconditions", "typed_stop_and_no_mutation_contract_must_be_exercised", "binade_target_retry_cap_quantum_and_claim_mutations_must_fail_closed", "canonical_hash_bound_result_required", "historical_PROTO14_files_may_not_be_modified", "runtime_repair_must_be_a_separate_post_commit_artifact", "independent_TDG7_binder_must_precede_runtime_implementation", "new_protocol_namespace_and_runtime_authorization_must_follow_runtime_qualification"})
    if any(value is not True for value in proof.values()):
        raise ValueError("TDG7 proof contract differs")
    successor = config["successor_boundary"]
    claims = config["claims"]
    expected_true = {"CAL11_terminal_invalid_runtime_result_preserved", "PROTO14_runtime_fault_cause_fully_derived", "TDG7_subdivision_lattice_design_frozen"}
    expected_claim_true = expected_true | {
        "historical_requested_width_and_rounded_boundary_difference_are_bitwise_distinct",
        "historical_failure_is_binary64_scheduler_contract_not_TDG6_temporal_admission",
        "historical_failure_is_not_constraint_spatial_source_health_or_physical_model",
    }
    expected_claim_keys = expected_claim_true | {"TDG7_independent_binder_completed", "TDG7_runtime_repair_implementation_authorized", "TDG7_runtime_repair_implemented", "PROTO14_runtime_mutation_authorized", "PROTO14_terminal_checkpoint_may_resume", "fresh_replacement_protocol_frozen", "fresh_GR0_calibration_authorized", "fresh_GR0_dynamic_calibration_completed", "GR0_case_eligible", "classical_spherical_diagnostic_authorized", "SGBL_execution_authorized", "FGCQR_holdout_execution_authorized", "FGCQR_mechanism_rejected", "DEF1_execution_authorized", "retained_EFT_evolution_authorized", "physical_transition_claim_authorized", "general_gradient_route_rejected", "singularity_resolution_derived", "child_domain_or_topology_derived", "dark_sector_mechanism_derived", "varying_locally_measured_c_derived"}
    expected_false = set(claims) - expected_true
    if (set(successor) != {"CAL11_terminal_invalid_runtime_result_preserved", "PROTO14_runtime_fault_cause_fully_derived", "TDG7_subdivision_lattice_design_frozen", "TDG7_independent_binder_authorized", "TDG7_independent_binder_completed", "TDG7_runtime_repair_implementation_authorized", "TDG7_runtime_repair_implemented", "PROTO14_runtime_mutation_authorized", "PROTO14_terminal_checkpoint_may_resume", "fresh_replacement_protocol_frozen", "fresh_GR0_calibration_authorized"}
            or successor.get("TDG7_independent_binder_authorized") is not True
            or any(successor.get(key) is not False for key in set(successor) - expected_true - {"TDG7_independent_binder_authorized"})
            or set(claims) != expected_claim_keys
            or any(claims.get(key) is not True for key in expected_claim_true)
            or any(claims.get(key) is not False for key in set(claims) - expected_claim_true)):
        raise ValueError("TDG7 claim boundary differs")


def load_config(path: Path = CONFIG) -> dict[str, Any]:
    if path.resolve() == CONFIG.resolve() and sha(path) != FROZEN_CONFIG_SHA256:
        raise ValueError("TDG7 frozen TOML bytes differ")
    with path.open("rb") as handle:
        config = tomllib.load(handle)
    validate_config_data(config)
    return config


def _bind_blobs(config: Mapping[str, Any]) -> dict[str, dict[str, str]]:
    lineage = config["immutable_lineage"]
    checkpoint = lineage["checkpoint_commit"]
    historic = lineage["historical_authorization_commit"]
    if _git("merge-base", "--is-ancestor", checkpoint, "HEAD", check=False).returncode:
        raise ValueError("TDG7 checkpoint commit is not an ancestor of HEAD")
    if _git("merge-base", "--is-ancestor", historic, checkpoint, check=False).returncode:
        raise ValueError("TDG7 historical authorization is not an ancestor")
    current = {
        "source_config": config["source_config"], "source_result": config["source_result"],
        "source_reproducer": config["source_reproducer"], "source_document": config["source_document"],
        "source_diagnosis_module": config["source_diagnosis_module"],
    }
    historical = {
        "historical_TDG6_runtime": config["historical_TDG6_runtime"],
        "historical_PROTO14_runner": config["historical_PROTO14_runner"],
        "historical_numerical_engine": config["historical_numerical_engine"],
        "historical_run_config": config["historical_run_config"],
        "runtime_authorization_result": config["runtime_authorization_result"],
    }
    bound: dict[str, dict[str, str]] = {"checkpoint": {}, "historical": {}}
    for kind, paths, commit in (("checkpoint", current, checkpoint), ("historical", historical, historic)):
        for name, relative in paths.items():
            digest = sha256(_blob(commit, relative)).hexdigest()
            if digest != lineage[f"{name}_sha256"]:
                raise ValueError(f"TDG7 immutable blob differs: {relative}")
            if kind == "checkpoint" and sha(REPOSITORY / relative) != digest:
                raise ValueError(f"TDG7 CAL11 worktree blob differs: {relative}")
            if name == "historical_numerical_engine" and sha(REPOSITORY / relative) != digest:
                raise ValueError("TDG7 historical numerical-engine worktree blob differs")
            bound[kind][relative] = digest
    old_runtime = _blob(historic, config["historical_TDG6_runtime"])
    old_runner = _blob(historic, config["historical_PROTO14_runner"])
    old_engine = _blob(historic, config["historical_numerical_engine"])
    for marker in (b"def _subdivision_boundaries", b"step = width / count", b"boundaries = tuple(start + index * step", b"TDG6 subdivision is not bitwise uniform"):
        if marker not in old_runtime:
            raise ValueError("TDG7 historical TDG6 arithmetic marker differs")
    # The runner owns the historical call boundary; the frozen witness binds
    # the inner CFL operation order independently, without importing it.
    for marker in (b"active.advance_to(", b"retry_factor=retry_factor", b"maximum_CFL_retries=32"):
        if marker not in old_runner:
            raise ValueError("TDG7 historical PROTO14 CFL operation order differs")
    # This pins the real stage-time expressions rather than inferring their
    # abscissae through HLT12 or a later runtime wrapper.
    for marker in (
        b'evaluate("rk4_k1", start, state)',
        b'evaluate("rk4_k2", start + dt / 2, y2)',
        b'evaluate("rk4_k4", start + dt, y4)',
        b'evaluate("ssprk3_s0", start, state)',
        b'evaluate("ssprk3_s1", start + dt, y1)',
        b'evaluate("ssprk3_s2", start + dt / 2, y2)',
        b'evaluate("candidate_endpoint", final_time, candidate)',
    ):
        if marker not in old_engine:
            raise ValueError("TDG7 historical numerical-engine stage marker differs")
    source = _blob(checkpoint, config["source_result"])
    parsed = json.loads(source)
    if canonical(parsed) != source or parsed.get("artifact_id") != "FGC-1-CAL11-PREF22":
        raise ValueError("TDG7 CAL11 result is not canonical")
    boundary = parsed.get("artifact_payload", {}).get("successor_boundary", {})
    if boundary.get("TDG7_subdivision_uniformity_diagnosis_design_authorized") is not True or boundary.get("TDG7_runtime_repair_implemented") is not False:
        raise ValueError("TDG7 is not authorized by CAL11")
    return bound


def _old_failures(witness: Mapping[str, Any]) -> dict[str, Any]:
    current = _hex(witness["current_time_hex"])
    # Frozen PROTO14 operation order: CFL cap first, then target minimum.
    cfl = _hex(witness["cfl_maximum_hex"])
    spacing = _hex(witness["grid_spacing_hex"])
    speed = _hex(witness["previous_speed_upper_hex"])
    target_remaining = _hex(witness["target_remaining_hex"])
    requested = min(target_remaining, cfl * spacing / speed)
    if requested.hex() != witness["requested_macro_width_hex"]:
        raise ValueError("TDG7 frozen requested width differs")
    reconstruction = reconstruct_historical_tdg6_guard(current, requested)
    if reconstruction.failed_counts != (1, 2, 4):
        raise ValueError("TDG7 historical guard failures differ")
    failures: dict[str, Any] = {}
    for count, step, boundaries, differences in zip(
        (1, 2, 4), reconstruction.nominal_steps, reconstruction.rounded_boundaries,
        reconstruction.rounded_differences, strict=True,
    ):
        failures[str(count)] = {
            "nominal_step_hex": step.hex(),
            "rounded_boundaries_hex": [item.hex() for item in boundaries],
            "rounded_adjacent_differences_hex": [item.hex() for item in differences],
            "guard_passes": False,
        }
    return {
        "requested_macro_width_hex": requested.hex(),
        "expected_final_hex": reconstruction.expected_final.hex(),
        "old_guard_failures": failures,
    }


def _four_q_midpoint_alias_attack(witness: Mapping[str, Any]) -> dict[str, Any]:
    """Prove why the historical endpoint-only 4Q plan is insufficient.

    Its four-quarter boundary values can be representable while an odd-Q fine
    width gives a half-Q exact midpoint.  The conversion back to binary64 is
    necessarily not that exact stage coordinate, so 4Q cannot certify the
    actual c=1/2 stage abscissa.
    """

    current = Fraction.from_float(_hex(witness["current_time_hex"]))
    q = Fraction.from_float(_hex(witness["event_ulp_hex"]))
    endpoint_only_macro = Fraction.from_float(
        _hex(witness["rounded_count_one_difference_hex"])
    )
    endpoint_only_fine = endpoint_only_macro / 4
    if endpoint_only_fine / q != 7_329_968_570_071:
        raise ValueError("TDG7 4Q midpoint attack width differs")
    exact_midpoint = current + endpoint_only_fine / 2
    rounded_midpoint = float(exact_midpoint)
    if Fraction.from_float(rounded_midpoint) == exact_midpoint:
        raise ValueError("TDG7 4Q midpoint attack unexpectedly exact")
    return {
        "endpoint_only_macro_quantum_hex": float(4 * q).hex(),
        "endpoint_only_fine_width_over_Q": 7_329_968_570_071,
        "exact_first_fine_midpoint": str(exact_midpoint),
        "rounded_first_fine_midpoint_hex": rounded_midpoint.hex(),
        "midpoint_exact_binary64": False,
        "four_q_endpoint_only_route_rejected": True,
    }


def _optional_checkpoint(config: Mapping[str, Any]) -> dict[str, Any]:
    path = REPOSITORY / config["source_checkpoint"]
    if not path.is_file():
        return {"whole_file_hash_required_before_metadata_read": True, "metadata_utf8_only": True, "physical_state_arrays_consumed": False, "clean_clone_absence_permitted": True}
    if sha(path) != config["immutable_lineage"]["source_checkpoint_sha256"]:
        raise ValueError("TDG7 checkpoint hash differs")
    # The ZIP member is read directly: no state/tracer/debit array is opened.
    import zipfile
    with zipfile.ZipFile(path) as archive:
        member = "metadata_utf8.npy"
        if member not in archive.namelist():
            raise ValueError("TDG7 checkpoint lacks metadata_utf8")
        raw = archive.read(member)
    # NPY v1/v2 header then uint8 payload; metadata is JSON and no ndarray state is loaded.
    if raw[:6] != b"\x93NUMPY":
        raise ValueError("TDG7 metadata encoding differs")
    header_length = int.from_bytes(raw[8:10] if raw[6] == 1 else raw[8:12], "little")
    offset = 10 + header_length if raw[6] == 1 else 12 + header_length
    metadata = json.loads(raw[offset:].decode("utf-8"))
    member = metadata.get("members", {}).get("RK4-2049", {})
    causal = member.get("causal_state", {})
    if float(causal.get("previous_speed_upper")).hex() != config["historical_event24_witness"]["previous_speed_upper_hex"]:
        raise ValueError("TDG7 checkpoint metadata speed differs")
    return {"whole_file_hash_required_before_metadata_read": True, "metadata_utf8_only": True,
            "physical_state_arrays_consumed": False, "clean_clone_absence_permitted": True}


def build(config: Mapping[str, Any]) -> dict[str, Any]:
    bindings = _bind_blobs(config)
    witness = config["historical_event24_witness"]
    historical = _old_failures(witness)
    plan = cal11_event24_witness()
    validate_tdg7_binary64_subdivision_plan(plan)
    if (plan.macro_width.hex() != witness["aligned_macro_width_hex"]
            or plan.fine_width.hex() != witness["aligned_fine_width_hex"]
            or plan.quantum.hex() != witness["event_ulp_hex"]
            or [item.hex() for item in plan.boundaries] != witness["selected_boundary_hex"]
            or [item[1].hex() for item in plan.four_quarter_stage_times] != witness["selected_fine_midpoint_hex"]
            or Fraction.from_float(plan.fine_width) / Fraction.from_float(plan.quantum)
            != witness["selected_fine_width_over_event_quantum"]
            or plan.conservative_reduction.hex() != witness["conservative_reduction_hex"]):
        raise ValueError("TDG7 selected lattice plan differs")
    if (historical["expected_final_hex"] != witness["historical_expected_macro_endpoint_hex"]
            or plan.endpoint.hex() != witness["selected_aligned_macro_endpoint_hex"]
            or plan.endpoint.hex() == historical["expected_final_hex"]):
        raise ValueError("TDG7 historical and selected macro endpoints differ")
    try:
        plan_forward_proto14_subdivision(23.0 / 16.0, 23.0 / 16.0 + 2.0**-52, 1.0)
    except CoordinateLatticeLimitReached as stopped:
        typed = stopped.reason
    else:
        raise ValueError("TDG7 typed stop control unexpectedly passed")
    claims = dict(config["claims"])
    # Execute representative mutation controls rather than serializing flags alone.
    if plan.macro_width > plan.requested_cap:
        raise ValueError("TDG7 cap control differs")
    try:
        plan_forward_proto14_subdivision(23.0 / 16.0, 24.0 / 16.0, 0.0)
    except CoordinateLatticeLimitReached as stopped:
        if stopped.reason != "nonpositive_requested_cap":
            raise ValueError("TDG7 cap mutation control differs")
    else:
        raise ValueError("TDG7 cap mutation unexpectedly passed")
    four_q_attack = _four_q_midpoint_alias_attack(witness)
    if historical["old_guard_failures"]["4"]["rounded_boundaries_hex"] == [
        item.hex() for item in plan.boundaries
    ]:
        raise ValueError("TDG7 historical 4Q path was mistaken for stage-safe plan")
    payload = {
        "immutable_lineage": {"checkpoint_commit": config["immutable_lineage"]["checkpoint_commit"], "historical_authorization_commit": config["immutable_lineage"]["historical_authorization_commit"], "checkpoint_commit_is_ancestor_of_HEAD": True, "historical_authorization_commit_is_ancestor_of_checkpoint_commit": True, "bound_source_sha256": bindings},
        "scope": dict(config["scope"]),
        "historical_event24_witness": historical,
        "selected_lattice_plan": {
            "event_quantum_hex": plan.quantum.hex(), "macro_quantum_hex": "0x1.0000000000000p-49",
            "macro_width_exact": "3664984285035/562949953421312", "macro_width_hex": plan.macro_width.hex(),
            "macro_eight_Q_tick_count": 3664984285035,
            "fine_width_exact": "3664984285035/2251799813685248", "fine_width_hex": plan.fine_width.hex(),
            "fine_width_over_Q": 7329968570070, "selection_budget_hex": plan.selection_budget.hex(),
            "conservative_reduction_exact": "531/576460752303423488", "conservative_reduction_hex": plan.conservative_reduction.hex(),
            "boundaries_hex": [item.hex() for item in plan.boundaries],
            "fine_stage_midpoints_hex": [item[1].hex() for item in plan.four_quarter_stage_times],
            "stage_times_exact_and_strictly_interior": True, "target_limited": plan.target_limited,
            "cap_not_exceeded": plan.macro_width <= plan.requested_cap,
        },
        "historical_vs_selected_first_macro_endpoint": {
            "historical_expected_endpoint_hex": historical["expected_final_hex"],
            "selected_aligned_endpoint_hex": plan.endpoint.hex(),
            "equal": False,
            "event_target_unchanged": True,
        },
        "four_q_endpoint_only_midpoint_attack": four_q_attack,
        "optional_checkpoint_metadata": _optional_checkpoint(config),
        "exact_design_preflight": {**dict(config["exact_controls"]), "typed_stop_reason": typed},
        "alternatives": dict(config["alternatives"]), "proof_contract": dict(config["proof_contract"]),
        "claim_boundary": claims, "successor_boundary": dict(config["successor_boundary"]),
        "mutation_controls": {"lineage_hash_mutation_rejected": True, "route_mutation_rejected": True, "quantum_or_threshold_mutation_rejected": True, "cap_round_up_mutation_rejected": True, "candidate_or_claim_promotion_rejected": True},
    }
    return {
        "artifact_id": ARTIFACT_ID, "artifact_payload": payload,
        "classification": "prospective_exact_binary64_subdivision_lattice_diagnosis_design_freeze",
        "derivation_document": "docs/fgc-tdg7-frz1.md", "derivation_document_sha256": sha(DOCUMENT),
        "gate_status": "PASS_EXACT_BINARY64_SUBDIVISION_LATTICE_DESIGN_FROZEN",
        "generated_by": "scripts/reproduce_fgc_tdg7_frz1.py",
        "implementation_sha256": {"configs/fgc/fgc-1-tdg7-frz1.toml": sha(CONFIG), "docs/fgc-tdg7-frz1.md": sha(DOCUMENT), "scripts/reproduce_fgc_tdg7_frz1.py": sha(Path(__file__)), "src/recursive_horizons/fgc/evolution/tdg7_binary64_subdivision_lattice.py": sha(DESIGN)},
        "metric_signature": config["metric_signature"], "riemann_convention": config["riemann_convention"],
        "nonclaims": {key: False for key, value in claims.items() if value is False},
        "project_version": config["project_version"], "schema_version": 1,
        "source_config_sha256": {"configs/fgc/fgc-1-tdg7-frz1.toml": sha(CONFIG)},
    }


def _atomic_write(path: Path, data: bytes) -> None:
    descriptor, temporary = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(data); handle.flush(); os.fsync(handle.fileno())
        Path(temporary).replace(path)
    except BaseException:
        Path(temporary).unlink(missing_ok=True); raise


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--verify", "--check", action="store_true")
    args = parser.parse_args()
    record = build(load_config())
    payload = canonical(record)
    if args.verify:
        if not args.output.is_file() or args.output.read_bytes() != payload:
            raise SystemExit("TDG7 tracked result differs from canonical reproduction")
    else:
        _atomic_write(args.output, payload)
    print(payload.decode("utf-8"), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
