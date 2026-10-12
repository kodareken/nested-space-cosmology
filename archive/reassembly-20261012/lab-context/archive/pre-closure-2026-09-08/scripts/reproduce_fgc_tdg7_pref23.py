#!/usr/bin/env python3
"""Construct or verify the independent TDG7 stage-safe lattice binder."""
from __future__ import annotations

import argparse
from copy import deepcopy
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

from recursive_horizons.fgc.evolution import tdg7_stage_safe_lattice_theorem as theorem  # noqa: E402

ARTIFACT_ID = "FGC-1-TDG7-PREF23"
PROJECT_VERSION = "0.11.0"
CONFIG = REPOSITORY / "configs/fgc/fgc-1-tdg7-pref23.toml"
OUTPUT = REPOSITORY / "results/fgc-1-tdg7-pref23.json"
DOCUMENT = REPOSITORY / "docs/fgc-tdg7-pref23.md"
BINDER_MODULE = REPOSITORY / "src/recursive_horizons/fgc/evolution/tdg7_stage_safe_lattice_theorem.py"

FREEZE_COMMIT = "ec403c3dbcec27babb861cfc8587d53a0136cfa0"
HISTORICAL_COMMIT = "7534a1662d8078a63c98025ecd464bb068fa012f"
CAL11_COMMIT = "b59d214832a1d3429f624579a93c6b9e1f544cf5"

FREEZE_PATHS = {
    "freeze_config": "configs/fgc/fgc-1-tdg7-frz1.toml",
    "freeze_result": "results/fgc-1-tdg7-frz1.json",
    "freeze_reproducer": "scripts/reproduce_fgc_tdg7_frz1.py",
    "freeze_module": "src/recursive_horizons/fgc/evolution/tdg7_binary64_subdivision_lattice.py",
    "freeze_document": "docs/fgc-tdg7-frz1.md",
}
HISTORICAL_PATHS = {
    "historical_TDG6_runtime": "src/recursive_horizons/fgc/evolution/tdg6_temporal_admission_runtime.py",
    "historical_numerical_engine": "src/recursive_horizons/fgc/evolution/numerical_engine.py",
    "historical_PROTO14_runner": "scripts/run_fgc_gr0_calibration_v14.py",
    "historical_run_config": "configs/fgc/fgc-1-cal9-run1.toml",
    "runtime_authorization_result": "results/fgc-1-hlt12-mon12.json",
}
CAL11_PATHS = {
    "cal11_result": "results/fgc-1-cal11-pref22.json",
    "cal11_config": "configs/fgc/fgc-1-cal11-pref22.toml",
}

EXPECTED_HASHES = {
    "freeze_config_sha256": "28a05df475e93c39805627cd05208c0078fbf0071a8ff90526a3dbc696c4619d",
    "freeze_result_sha256": "9d6e4426653fe159338a12513e9e5e0e5d6b57e59dae210b008e8a6ea6928ff7",
    "freeze_reproducer_sha256": "72bfc07a3f949429e406eede93d4056e1ac5b1f37568c4be3ea4328560ab7b4f",
    "freeze_module_sha256": "2bc86880d46a3b4c6cd743559d0bc4b51b09058382341dfffe6358b4565e44a9",
    "freeze_document_sha256": "3a4a8ae70d418b2d9d32a8209a0235051f3d1df8b776bfb4a9d4616a48b2b75d",
    "historical_TDG6_runtime_sha256": "039de0bd008a2536bea6f7187afc3856bfc841d14790cbae5ca2fe27f7eb8e2d",
    "historical_numerical_engine_sha256": "8ea99604a85b4e1acbf869ed2462ea4aa8301c06e068330a1d0ce66a51a0bebf",
    "historical_PROTO14_runner_sha256": "4c13eaf86cf31f55e28b7f60661db78ad22e1bb5b7460adf7aeff455efdc8e45",
    "historical_run_config_sha256": "8e6fafff638ffd49549efab8b6f969a01964c915116db9ce3698953201172176",
    "runtime_authorization_result_sha256": "7f4068a3116ff1838ce7f0165ae18e72d2eea68933b87f056b26d8149c935731",
    "cal11_result_sha256": "1339792198ec78ad83bc3e11a2805f09f961607ab6c1f6f06e79bb88ba6b58b8",
    "cal11_config_sha256": "685ee05442d7c92ce83994d923890f4fa208b10ffbf2946026b07e9a9cf4d34d",
}

EXPECTED_SCOPE = {
    "target": "FGC-2-SF1-TDG7",
    "role": "post_freeze_independent_stage_safe_exact_binary64_shared_lattice_theorem_and_runtime_repair_feasibility_binder",
    "calibration_branch": "GR-0",
    "freeze_compact_result_read": True,
    "historical_source_blobs_read_from_immutable_commits": True,
    "CAL11_compact_result_read": True,
    "actual_campaign_history_arrays_consumed": False,
    "campaign_checkpoint_loaded": False,
    "campaign_checkpoint_resumed": False,
    "production_state_advanced": False,
    "output_namespace_created": False,
    "synthetic_exact_and_source_shape_controls_only": True,
    "SGBL_trajectory_read": False,
    "FGCQR_trajectory_read": False,
    "physical_or_candidate_question_answered": False,
}

CLAIM_TRUE = {
    "CAL11_terminal_invalid_runtime_result_preserved",
    "PROTO14_runtime_fault_cause_fully_derived",
    "historical_requested_width_and_rounded_boundary_difference_are_bitwise_distinct",
    "historical_failure_is_binary64_scheduler_contract_not_TDG6_temporal_admission",
    "historical_failure_is_not_constraint_spatial_source_health_or_physical_model",
    "TDG7_subdivision_lattice_design_frozen",
    "TDG7_independent_binder_completed",
    "TDG7_runtime_repair_implementation_authorized",
}
CLAIM_FALSE = {
    "TDG7_runtime_repair_implemented", "PROTO14_runtime_mutation_authorized",
    "PROTO14_terminal_checkpoint_may_resume", "fresh_replacement_protocol_frozen",
    "fresh_GR0_calibration_authorized", "fresh_GR0_dynamic_calibration_completed",
    "GR0_case_eligible", "classical_spherical_diagnostic_authorized",
    "SGBL_execution_authorized", "FGCQR_holdout_execution_authorized",
    "FGCQR_mechanism_rejected", "DEF1_execution_authorized",
    "retained_EFT_evolution_authorized", "physical_transition_claim_authorized",
    "general_gradient_route_rejected", "singularity_resolution_derived",
    "child_domain_or_topology_derived", "dark_sector_mechanism_derived",
    "varying_locally_measured_c_derived",
}


def canonical(value: object) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True, ensure_ascii=True,
                       allow_nan=False) + "\n").encode("utf-8")


def _sha_bytes(payload: bytes) -> str:
    return sha256(payload).hexdigest()


def sha(path: Path) -> str:
    return _sha_bytes(path.read_bytes())


def _git(*args: str, check: bool = True) -> subprocess.CompletedProcess[bytes]:
    return subprocess.run(["git", *args], cwd=REPOSITORY, check=check,
                          stdout=subprocess.PIPE, stderr=subprocess.PIPE)


def _blob(commit: str, relative: str) -> bytes:
    return _git("show", f"{commit}:{relative}").stdout


def _strict_keys(name: str, mapping: Mapping[str, Any], expected: set[str]) -> None:
    if set(mapping) != expected:
        raise ValueError(f"PREF23 {name} keys differ")


def _all_true(name: str, mapping: Mapping[str, Any]) -> None:
    if not mapping or set(mapping.values()) != {True}:
        raise ValueError(f"PREF23 {name} must be all true")


def _expected_lineage() -> dict[str, Any]:
    return {
        "freeze_commit": FREEZE_COMMIT,
        "historical_authorization_commit": HISTORICAL_COMMIT,
        "cal11_checkpoint_commit": CAL11_COMMIT,
        **EXPECTED_HASHES,
        "freeze_commit_must_be_ancestor_of_HEAD": True,
        "historical_authorization_commit_must_be_ancestor_of_freeze_commit": True,
        "cal11_checkpoint_commit_must_be_ancestor_of_freeze_commit": True,
        "all_freeze_blobs_must_match_freeze_commit_and_worktree": True,
        "historical_source_blobs_must_match_historical_commit": True,
        "historical_numerical_engine_must_be_read_from_historical_commit_not_live_worktree": True,
        "CAL11_blobs_must_match_CAL11_commit_and_worktree": True,
        "freeze_result_must_be_canonical_and_authorize_independent_binder_only": True,
        "raw_campaign_bundle_must_not_be_loaded": True,
    }


def _expected_theorem_contract() -> dict[str, Any]:
    return {
        "theorem_must_not_import_TDG7_design_helper": True,
        "event_quantum": "max_binary64_ulp_of_current_time_and_event_target",
        "endpoint_only_macro_quantum": "4*Q",
        "stage_safe_macro_quantum": "8*Q",
        "selected_width": "floor(min(requested_cap,target-current)/(8*Q))*(8*Q)",
        "fine_width": "selected_width/4",
        "fine_width_must_be_even_multiple_of_Q": True,
        "stage_abscissae": ["0", "1/2", "1"],
        "all_stage_times_must_be_exact_binary64_shared_lattice_coordinates": True,
        "positive_fine_midpoints_must_be_strictly_interior": True,
        "four_Q_endpoint_only_midpoint_alias_must_be_independently_reproduced": True,
        "rounding_up_is_forbidden": True,
        "target_stage_lattice_alignment_required": True,
        "typed_coordinate_lattice_stop_precedes_shadow_or_real_state_operation": True,
    }


def _expected_witness() -> dict[str, Any]:
    return {
        "current_time_hex": "0x1.7000000000000p+0",
        "target_time_hex": "0x1.8000000000000p+0",
        "requested_macro_width_hex": "0x1.aaa90b0fb5c26p-8",
        "historical_expected_macro_endpoint_hex": "0x1.71aaa90b0fb5cp+0",
        "event_quantum_hex": "0x1.0000000000000p-52",
        "stage_safe_macro_quantum_hex": "0x1.0000000000000p-49",
        "stage_safe_eight_Q_tick_count": 3664984285035,
        "selected_macro_width_exact": "3664984285035/562949953421312",
        "selected_macro_width_hex": "0x1.aaa90b0fb5800p-8",
        "selected_fine_width_exact": "3664984285035/2251799813685248",
        "selected_fine_width_hex": "0x1.aaa90b0fb5800p-10",
        "selected_fine_width_over_Q": 7329968570070,
        "selected_aligned_macro_endpoint_hex": "0x1.71aaa90b0fb58p+0",
        "conservative_reduction_exact": "531/576460752303423488",
        "conservative_reduction_hex": "0x1.0980000000000p-50",
        "selected_boundaries_hex": ["0x1.7000000000000p+0", "0x1.706aaa42c3ed6p+0", "0x1.70d5548587dacp+0", "0x1.713ffec84bc82p+0", "0x1.71aaa90b0fb58p+0"],
        "selected_fine_midpoints_hex": ["0x1.7035552161f6bp+0", "0x1.709fff6425e41p+0", "0x1.710aa9a6e9d17p+0", "0x1.717553e9adbedp+0"],
        "historical_guard_fails_for_counts": [1, 2, 4],
        "historical_and_selected_first_macro_endpoints_are_equal": False,
        "event_target_is_unchanged": True,
    }


def _expected_authorization() -> dict[str, Any]:
    return {
        "authorized_object": "implementation_and_synthetic_qualification_of_the_frozen_TDG7_stage_safe_exact_binary64_shared_lattice_runtime_repair",
        "independent_binder_completed": True,
        "TDG7_runtime_repair_implementation_authorized": True,
        "one_shared_stage_safe_plan_must_own_all_one_two_four_paths": True,
        "actual_RK4_and_SSPRK3_stage_times_must_use_shared_plan_coordinates": True,
        "typed_lattice_stop_must_precede_shadow_source_projector_tracer_monitor_or_state_work": True,
        "physical_equations_source_projector_spatial_operator_and_TDG6_threshold_unchanged": True,
        "historical_PROTO14_runtime_mutation_authorized": False,
        "historical_terminal_checkpoint_may_resume": False,
        "new_protocol_or_output_namespace_authorized": False,
        "fresh_GR0_calibration_authorized": False,
        "candidate_execution_authorized": False,
        "production_trajectory_authorized": False,
    }


def validate_config_data(config: Mapping[str, Any]) -> None:
    _strict_keys("config", config, {
        "schema_version", "artifact_id", "project_version", "metric_signature", "riemann_convention",
        *FREEZE_PATHS, *HISTORICAL_PATHS, *CAL11_PATHS, "binder_module", "scope",
        "immutable_lineage", "independent_stage_safe_theorem", "historical_event24_witness",
        "runtime_repair_authorization", "limitation_audit", "proof_contract",
        "successor_boundary", "claims",
    })
    expected_paths = {**FREEZE_PATHS, **HISTORICAL_PATHS, **CAL11_PATHS,
                      "binder_module": "src/recursive_horizons/fgc/evolution/tdg7_stage_safe_lattice_theorem.py"}
    if (config.get("schema_version") != 1 or config.get("artifact_id") != ARTIFACT_ID
            or config.get("project_version") != PROJECT_VERSION
            or config.get("metric_signature") != "-+++"
            or config.get("riemann_convention") != "plus_partial_mu_gamma_nu"
            or any(config.get(key) != value for key, value in expected_paths.items())
            or config.get("scope") != EXPECTED_SCOPE
            or config.get("immutable_lineage") != _expected_lineage()
            or config.get("independent_stage_safe_theorem") != _expected_theorem_contract()
            or config.get("historical_event24_witness") != _expected_witness()
            or config.get("runtime_repair_authorization") != _expected_authorization()):
        raise ValueError("PREF23 identity, lineage, theorem, witness, or authorization differs")
    _strict_keys("limitation audit", config["limitation_audit"], {
        "independent_theorem_is_not_runtime_repair_implementation", "historical_CAL11_abort_remains_invalid_runtime_not_physics",
        "stage_safe_coordinate_plan_is_not_TDG6_temporal_admission_result", "source_shape_feasibility_is_not_synthetic_qualification",
        "new_runtime_must_be_qualified_on_synthetic_controls", "historical_checkpoint_may_not_resume",
        "new_protocol_namespace_and_runtime_authorization_remain_required_before_calibration", "no_GR0_SGBL_or_FGCQR_outcome_is_classified",
    })
    _all_true("limitation audit", config["limitation_audit"])
    _strict_keys("proof contract", config["proof_contract"], {
        "freeze_commit_and_all_freeze_blobs_must_match", "CAL11_and_historical_source_blobs_must_match_their_immutable_commits",
        "theorem_must_rederive_eight_Q_without_importing_freeze_design_helper", "historical_one_two_four_guard_failures_must_be_rederived_from_source_expression",
        "historical_numerical_engine_stage_abscissae_must_be_inspected_directly", "four_Q_midpoint_alias_and_eight_Q_stage_exactness_must_be_mutation_attacked",
        "all_witness_fractions_boundaries_and_midpoints_must_match_exactly", "runtime_repair_authorization_must_not_promote_calibration_or_candidate_execution",
        "no_history_checkpoint_state_tracer_or_debit_array_may_be_loaded", "lineage_theorem_authorization_nonclaim_and_claim_mutations_must_fail_closed",
        "canonical_hash_bound_result_required",
    })
    _all_true("proof contract", config["proof_contract"])
    expected_successor = {
        "CAL11_terminal_invalid_runtime_result_preserved": True, "PROTO14_runtime_fault_cause_fully_derived": True,
        "TDG7_subdivision_lattice_design_frozen": True, "TDG7_independent_binder_completed": True,
        "TDG7_runtime_repair_implementation_authorized": True, "TDG7_runtime_repair_implemented": False,
        "PROTO14_runtime_mutation_authorized": False, "PROTO14_terminal_checkpoint_may_resume": False,
        "fresh_replacement_protocol_frozen": False, "fresh_GR0_calibration_authorized": False,
    }
    if config.get("successor_boundary") != expected_successor:
        raise ValueError("PREF23 successor boundary differs")
    claims = config.get("claims", {})
    if set(claims) != CLAIM_TRUE | CLAIM_FALSE or any(claims.get(key) is not True for key in CLAIM_TRUE) or any(claims.get(key) is not False for key in CLAIM_FALSE):
        raise ValueError("PREF23 claim boundary differs")


def load_config(path: Path = CONFIG) -> dict[str, Any]:
    with path.open("rb") as handle:
        config = tomllib.load(handle)
    validate_config_data(config)
    return config


def _canonical_object(payload: bytes, name: str) -> dict[str, Any]:
    parsed = json.loads(payload)
    if not isinstance(parsed, dict) or canonical(parsed) != payload:
        raise ValueError(f"PREF23 {name} is not canonical JSON")
    return parsed


def _verify_lineage(config: Mapping[str, Any]) -> tuple[dict[str, Any], bytes, bytes]:
    lineage = config["immutable_lineage"]
    for commit in (FREEZE_COMMIT, HISTORICAL_COMMIT, CAL11_COMMIT):
        _git("cat-file", "-e", f"{commit}^{{commit}}")
    if _git("merge-base", "--is-ancestor", FREEZE_COMMIT, "HEAD", check=False).returncode:
        raise ValueError("PREF23 freeze commit is not an ancestor of HEAD")
    for older in (HISTORICAL_COMMIT, CAL11_COMMIT):
        if _git("merge-base", "--is-ancestor", older, FREEZE_COMMIT, check=False).returncode:
            raise ValueError("PREF23 predecessor commit is not an ancestor of freeze")
    ledgers: dict[str, dict[str, str]] = {"freeze": {}, "historical": {}, "cal11": {}}
    for label, commit, paths in (("freeze", FREEZE_COMMIT, FREEZE_PATHS), ("historical", HISTORICAL_COMMIT, HISTORICAL_PATHS), ("cal11", CAL11_COMMIT, CAL11_PATHS)):
        for field, relative in paths.items():
            expected = lineage[f"{field}_sha256"]
            payload = _blob(commit, relative)
            if _sha_bytes(payload) != expected:
                raise ValueError(f"PREF23 immutable source differs: {relative}")
            # Freeze and CAL11 are sealed inputs whose current bytes must still agree.
            if label in {"freeze", "cal11"} and sha(REPOSITORY / relative) != expected:
                raise ValueError(f"PREF23 sealed worktree source differs: {relative}")
            ledgers[label][relative] = expected
    freeze = _canonical_object(_blob(FREEZE_COMMIT, FREEZE_PATHS["freeze_result"]), "freeze result")
    if (freeze.get("artifact_id") != "FGC-1-TDG7-FRZ1"
            or freeze.get("gate_status") != "PASS_EXACT_BINARY64_SUBDIVISION_LATTICE_DESIGN_FROZEN"
            or freeze.get("artifact_payload", {}).get("successor_boundary", {}).get("TDG7_independent_binder_authorized") is not True
            or freeze.get("artifact_payload", {}).get("successor_boundary", {}).get("TDG7_runtime_repair_implementation_authorized") is not False):
        raise ValueError("PREF23 freeze boundary differs")
    cal11 = _canonical_object(_blob(CAL11_COMMIT, CAL11_PATHS["cal11_result"]), "CAL11 result")
    if (cal11.get("artifact_id") != "FGC-1-CAL11-PREF22"
            or cal11.get("gate_status", {}).get("PROTO14_invalid_runtime_or_nonconverged_terminal_observed") is not True):
        raise ValueError("PREF23 CAL11 boundary differs")
    engine = _blob(HISTORICAL_COMMIT, HISTORICAL_PATHS["historical_numerical_engine"])
    tdg6_runtime = _blob(HISTORICAL_COMMIT, HISTORICAL_PATHS["historical_TDG6_runtime"])
    return ({
        "freeze_commit": FREEZE_COMMIT, "freeze_commit_is_ancestor_of_HEAD": True,
        "historical_authorization_commit": HISTORICAL_COMMIT,
        "historical_authorization_commit_is_ancestor_of_freeze_commit": True,
        "cal11_checkpoint_commit": CAL11_COMMIT,
        "cal11_checkpoint_commit_is_ancestor_of_freeze_commit": True,
        "bound_source_sha256": ledgers,
    }, engine, tdg6_runtime)


def _historical_guard_payload() -> dict[str, Any]:
    witness = theorem.cal11_stage_safe_witness_certificate()
    guard = theorem.reconstruct_historical_tdg6_guard(
        float.fromhex("0x1.7000000000000p+0"), float.fromhex(witness["requested_cap_hex"]),
    )
    controls: dict[str, Any] = {}
    for count, step, boundaries, differences in zip((1, 2, 4), guard.nominal_steps, guard.boundaries, guard.adjacent_differences, strict=True):
        controls[str(count)] = {
            "nominal_step_hex": step.hex(), "rounded_boundaries_hex": [value.hex() for value in boundaries],
            "rounded_adjacent_differences_hex": [value.hex() for value in differences], "guard_passes": False,
        }
    return {"requested_cap_hex": witness["requested_cap_hex"], "historical_expected_endpoint_hex": guard.expected_final.hex(), "historical_failed_counts": list(guard.failed_counts), "old_guard_failures": controls}


def _independent_binder(
    config: Mapping[str, Any], engine: bytes, tdg6_runtime: bytes
) -> dict[str, Any]:
    source = BINDER_MODULE.read_text(encoding="utf-8")
    forbidden = "tdg7_binary64_subdivision_lattice"
    if forbidden in source:
        raise ValueError("PREF23 independent theorem imports TDG7 design helper")
    witness = theorem.cal11_stage_safe_witness_certificate()
    properties = theorem.deterministic_stage_safe_property_certificate()
    historical = _historical_guard_payload()
    engine_contract = theorem.parse_immutable_engine_stage_contract(engine.decode("utf-8"))
    historical_source_contract = theorem.parse_historical_tdg6_subdivision_contract(
        tdg6_runtime.decode("utf-8")
    )
    expected_witness = config["historical_event24_witness"]
    expected = {
        "requested_cap_hex": expected_witness["requested_macro_width_hex"],
        "historical_failed_counts": expected_witness["historical_guard_fails_for_counts"],
        "historical_expected_endpoint_hex": expected_witness["historical_expected_macro_endpoint_hex"],
        "event_quantum_hex": expected_witness["event_quantum_hex"],
        "macro_quantum_hex": expected_witness["stage_safe_macro_quantum_hex"],
        "macro_width_hex": expected_witness["selected_macro_width_hex"],
        "fine_width_hex": expected_witness["selected_fine_width_hex"],
        "fine_width_over_Q": expected_witness["selected_fine_width_over_Q"],
        "conservative_reduction": expected_witness["conservative_reduction_exact"],
        "selected_boundaries_hex": expected_witness["selected_boundaries_hex"],
        "selected_fine_midpoints_hex": expected_witness["selected_fine_midpoints_hex"],
        "endpoint_only_four_q_counterexample_derived_from_historical_one_step": True,
        "endpoint_only_four_q_midpoint_exact": False,
        "endpoint_only_four_q_midpoint_binary64_hex": "0x1.7035552161f6cp+0",
        "stage_triplet_count": 7,
    }
    if witness != expected or historical["historical_failed_counts"] != [1, 2, 4]:
        raise ValueError("PREF23 independent CAL11 witness differs")
    return {
        "theorem_module_imports_TDG7_design_helper": False,
        "independent_CAL11_event24_witness": witness,
        "independent_historical_TDG6_guard_reconstruction": historical,
        "independent_historical_TDG6_source_contract": {
            "step_is_width_over_count": historical_source_contract.step_is_width_over_count,
            "boundaries_are_start_plus_index_times_step": historical_source_contract.boundaries_are_start_plus_index_times_step,
            "expected_final_is_start_plus_width": historical_source_contract.expected_final_is_start_plus_width,
            "macro_endpoint_guard_present": historical_source_contract.macro_endpoint_guard_present,
            "adjacent_uniformity_guard_present": historical_source_contract.adjacent_uniformity_guard_present,
        },
        "independent_engine_stage_abscissa_contract": {
            "stage_time_kinds": [list(item) for item in engine_contract.stage_time_kinds],
            "candidate_endpoint_uses_final_time": engine_contract.candidate_endpoint_uses_final_time,
        },
        "deterministic_stage_safe_property_controls": properties,
        "all_independent_binder_controls_pass": True,
        "TDG7_independent_binder_completed": True,
        "TDG7_runtime_repair_implementation_authorized": True,
        "TDG7_runtime_repair_implemented": False,
    }


def _config_mutation_rejected(config: Mapping[str, Any], mutate: Any) -> bool:
    candidate = deepcopy(config)
    mutate(candidate)
    try:
        validate_config_data(candidate)
    except ValueError:
        return True
    return False


def _executed_mutation_controls(
    config: Mapping[str, Any], engine: bytes, tdg6_runtime: bytes
) -> dict[str, bool]:
    """Execute every attack advertised by the canonical certificate."""

    lineage = _config_mutation_rejected(
        config,
        lambda value: value["immutable_lineage"].__setitem__(
            "freeze_result_sha256", "0" * 64
        ),
    )
    four_q = _config_mutation_rejected(
        config,
        lambda value: value["independent_stage_safe_theorem"].__setitem__(
            "stage_safe_macro_quantum", "4*Q"
        ),
    )
    witness = _config_mutation_rejected(
        config,
        lambda value: value["historical_event24_witness"].__setitem__(
            "selected_fine_midpoints_hex", ["0x1.0p+0"] * 4
        ),
    )
    claim = _config_mutation_rejected(
        config,
        lambda value: value["claims"].__setitem__(
            "FGCQR_holdout_execution_authorized", True
        ),
    )
    broken_engine = engine.decode("utf-8").replace(
        'evaluate("rk4_k2", start + dt / 2, y2)',
        'evaluate("rk4_k2", start + dt, y2)',
    )
    try:
        theorem.parse_immutable_engine_stage_contract(broken_engine)
    except theorem.IndependentTDG7TheoremError:
        engine_stage = True
    else:
        engine_stage = False
    broken_runtime = tdg6_runtime.decode("utf-8").replace(
        "step = width / count", "step = width / (count + 1)"
    )
    try:
        theorem.parse_historical_tdg6_subdivision_contract(broken_runtime)
    except theorem.IndependentTDG7TheoremError:
        historical_expression = True
    else:
        historical_expression = False
    broken_guard_action = tdg6_runtime.decode("utf-8").replace(
        "TDG6 subdivision is not bitwise uniform", "TDG6 nonuniformity ignored"
    )
    try:
        theorem.parse_historical_tdg6_subdivision_contract(broken_guard_action)
    except theorem.IndependentTDG7TheoremError:
        historical_guard_action = True
    else:
        historical_guard_action = False
    try:
        theorem.derive_stage_safe_lattice(23.0 / 16.0, 24.0 / 16.0, 0.0)
    except theorem.IndependentTDG7TheoremError as error:
        cap = error.reason == "nonpositive_requested_cap"
    else:
        cap = False
    return {
        "lineage_hash_mutation_rejected": lineage,
        "four_Q_route_mutation_rejected": four_q,
        "witness_width_or_midpoint_mutation_rejected": witness,
        "historical_engine_stage_mutation_rejected": engine_stage,
        "historical_source_boundary_expression_mutation_rejected": (
            historical_expression and historical_guard_action
        ),
        "cap_typed_stop_control_passed": cap,
        "candidate_or_claim_promotion_rejected": claim,
    }


def build(config: Mapping[str, Any]) -> dict[str, Any]:
    lineage, engine, tdg6_runtime = _verify_lineage(config)
    binder = _independent_binder(config, engine, tdg6_runtime)
    claims = dict(config["claims"])
    authorization = dict(config["runtime_repair_authorization"])
    mutation_controls = _executed_mutation_controls(config, engine, tdg6_runtime)
    if set(mutation_controls.values()) != {True}:
        raise ValueError("PREF23 executable mutation control differs")
    payload = {
        "immutable_lineage": lineage,
        "freeze_crosscheck": {
            "freeze_artifact_id": "FGC-1-TDG7-FRZ1", "freeze_is_canonical": True,
            "freeze_authorized_independent_binder_only": True,
        },
        "independent_binder": binder,
        "runtime_repair_authorization": authorization,
        "limitation_audit": dict(config["limitation_audit"]),
        "proof_contract": dict(config["proof_contract"]),
        "mutation_controls": mutation_controls,
        "claim_boundary": claims,
        "successor_boundary": dict(config["successor_boundary"]),
    }
    return {
        "artifact_id": ARTIFACT_ID,
        "artifact_payload": payload,
        "classification": "independent_stage_safe_binary64_subdivision_binder_and_runtime_repair_implementation_authorization",
        "derivation_document": "docs/fgc-tdg7-pref23.md",
        "derivation_document_sha256": sha(DOCUMENT),
        "gate_status": "PASS_TDG7_RUNTIME_REPAIR_IMPLEMENTATION_AUTHORIZED",
        "generated_by": "scripts/reproduce_fgc_tdg7_pref23.py",
        "implementation_sha256": {
            "configs/fgc/fgc-1-tdg7-pref23.toml": sha(CONFIG),
            "docs/fgc-tdg7-pref23.md": sha(DOCUMENT),
            "scripts/reproduce_fgc_tdg7_pref23.py": sha(Path(__file__)),
            "src/recursive_horizons/fgc/evolution/tdg7_stage_safe_lattice_theorem.py": sha(BINDER_MODULE),
        },
        "metric_signature": config["metric_signature"],
        "nonclaims": {key: False for key in sorted(CLAIM_FALSE)},
        "project_version": PROJECT_VERSION,
        "riemann_convention": config["riemann_convention"],
        "schema_version": 1,
        "source_config_sha256": sha(CONFIG),
    }


def _atomic_write(path: Path, payload: bytes) -> None:
    descriptor, temporary = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        Path(temporary).replace(path)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--verify", "--check", action="store_true")
    args = parser.parse_args()
    record = build(load_config())
    output = args.output.resolve()
    payload = canonical(record)
    if args.verify:
        if not output.is_file() or output.read_bytes() != payload:
            raise SystemExit("PREF23 tracked result differs from canonical reproduction")
    else:
        _atomic_write(output, payload)
    print(payload.decode("utf-8"), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
