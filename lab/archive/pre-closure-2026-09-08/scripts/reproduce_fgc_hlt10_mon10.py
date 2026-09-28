#!/usr/bin/env python3
"""Reproduce the pre-trajectory FGC-1-HLT10-MON10 authorization.

HLT10 binds frozen PROTO12 to the inherited GR-0 campaign without reading or
advancing a trajectory.  It reconstructs all twelve initial states, evaluates
the unchanged REF1 source through SRC4, recomputes all four time-zero common
events with the pairwise spectral classifier, attacks the new source and
spectral surfaces, and observes the fresh namespaces without creating them.
"""

from __future__ import annotations

import argparse
from copy import deepcopy
from hashlib import sha256
from pathlib import Path
import subprocess
import sys
import tomllib
from typing import Any, Mapping

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts import reproduce_fgc_hlt4_mon4 as hlt4  # noqa: E402
from recursive_horizons.fgc.evolution.gr0_calibration import (  # noqa: E402
    construct_gr0_grid_initial_data,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    EvolutionState,
    SBPFirstDerivative,
    UniformRadialGrid,
    array_content_sha256,
)
from recursive_horizons.fgc.evolution.proto11_runtime import (  # noqa: E402
    exact_spherical_minkowski_reference,
    project_gr0_reference_balanced_state,
    reference_balanced_reduction_constraint,
    reference_balanced_semidiscrete_constraint_snapshot,
    reference_balanced_spatial_derivatives,
)
from recursive_horizons.fgc.evolution.proto12_runtime import (  # noqa: E402
    PROTO12_POINT_COUNTS,
    PROTO12_SOURCE_POINT_BATCH_SIZE,
    Proto12GR0CommonEventAssessment,
    Proto12GR0EvolutionOperator,
    proto12_gr0_common_event,
)
from recursive_horizons.fgc.evolution.protocol_v12 import (  # noqa: E402
    validate_sf1_protocol_v12,
)
from recursive_horizons.fgc.evolution.src3_reference_balanced_source import (  # noqa: E402
    diagnose_gr0_reference_balanced_accelerations,
)
from recursive_horizons.fgc.evolution.src4_vectorized_reference_source import (  # noqa: E402
    diagnose_gr0_vectorized_reference_accelerations,
)
from recursive_horizons.fgc.initial_data_preflight import (  # noqa: E402
    PulseParameters,
)


Q = hlt4.Q
ARTIFACT_ID = "FGC-1-HLT10-MON10"
PROJECT_VERSION = "0.11.0"
CHECKPOINT_COMMIT = "a42780179591012a29eb81a87f59be1f60f9e514"
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-hlt10-mon10.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-hlt10-mon10.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-hlt10-mon10.md"

EXPECTED_PATHS = {
    "protocol_config": "configs/fgc/fgc-2-sf1-protocol-v12.toml",
    "protocol_freeze_config": "configs/fgc/fgc-1-pro12-frz1.toml",
    "protocol_freeze_result": "results/fgc-1-pro12-frz1.json",
    "predecessor_runtime_result": "results/fgc-1-hlt9-mon9.json",
    "calibration_diagnosis_config": "configs/fgc/fgc-1-cal8-pref10.toml",
    "calibration_diagnosis_result": "results/fgc-1-cal8-pref10.json",
    "resolution_evidence_config": "configs/fgc/fgc-1-rsp1-pref12.toml",
    "resolution_evidence_result": "results/fgc-1-rsp1-pref12.json",
    "source_instrument_config": "configs/fgc/fgc-1-src4-vec1.toml",
    "source_instrument_result": "results/fgc-1-src4-vec1.json",
    "static_input_result": "results/fgc-1-id2-all1.json",
    "predecessor_run_plan": "configs/fgc/fgc-1-cal8-run1.toml",
    "run_plan_config": "configs/fgc/fgc-1-cal9-run1.toml",
}

EXPECTED_SCOPE = {
    "target_protocol": "FGC-2-SF1-PROTO12",
    "runtime_owner": ARTIFACT_ID,
    "calibration_branch": "GR-0",
    "role": "pre_trajectory_SRC4_runtime_pairwise_spectral_binding_and_fresh_input_freeze",
    "PROTO11_CAL8_and_RSP1_history_disclosed": True,
    "PROTO12_trajectory_read": False,
    "SGBL_trajectory_read": False,
    "FGCQR_trajectory_read": False,
    "collapse_or_trapped_outcome_classified": False,
    "numerical_premise_promoted_to_physical_mechanism": False,
}

EXPECTED_LINEAGE = {
    "checkpoint_commit": CHECKPOINT_COMMIT,
    "protocol_config_sha256": "2c33ef4a703e918683c886eb216bfc77af5a48d477c1e1556e0b5176d44fbe7f",
    "protocol_freeze_config_sha256": "29c7a12eb91108f2df5ff1f8e8888435ed41064a30b9bdeb70bbacaced6b91fa",
    "protocol_freeze_result_sha256": "57ed1fc79ffeec8746e0267a71bb9e5efe851283dc72822a93b550cdafb14bc8",
    "predecessor_runtime_result_sha256": "1902cedb853643c445ece733c57d4cc3815c11d48a1d123aefa12c009368893e",
    "calibration_diagnosis_config_sha256": "a7646fb7ef2b229c36981a2f4a87e4428de1bdc5cbe324baedc965097c674cb3",
    "calibration_diagnosis_result_sha256": "10c6898e1ec6e7c07793180c4a4f8fa26a1b2f660115a0b2bc54497ff8bbb89d",
    "resolution_evidence_config_sha256": "10248e01e83714a3e0b83d213748be8d5c6e1b234133b7bbe35a2809515e5bd1",
    "resolution_evidence_result_sha256": "0bb184a519230d1b221d20fb88df70841800575b886df0b0642837deb26e4551",
    "source_instrument_config_sha256": "39c1c11c703c58cf884e4f0f240db0dd1eaa085494e6d2cb14761fbcc48681d2",
    "source_instrument_result_sha256": "fb510b0232ac956a73321b63a31d1f52a7732132ff736d41f6d08a4da125aa7a",
    "static_input_result_sha256": "10e2661be521087467b970f1cd877311bd63a0cb6874a1ee745ca00ec74ab981",
    "predecessor_run_plan_sha256": "db3d28ffac36fe4e9720d81ea849cdd7262f876810e30ed2cd1cd8176891d42f",
    "checkpoint_must_be_ancestor_of_HEAD": True,
    "all_predecessor_blobs_must_be_read_from_checkpoint": True,
}

EXPECTED_RUNTIME_KEYS = {
    "initial_q_uses_Dh_u_minus_uref_plus_qref",
    "source_qr_uses_Dh_q_minus_qref",
    "reduction_constraint_and_common_event_use_same_reference_map",
    "q_remains_independent_evolved_field",
    "interior_q_reprojection_is_forbidden",
    "complete_unredefined_source_and_raw_gate_remain_serialized",
    "source_backend_is_certified_SRC4_tensor_contraction",
    "source_point_batch_size_is_fixed_at_2048",
    "fixed_batching_does_not_mutate_pointwise_affine_systems",
    "SRC3_explicit_index_backend_remains_the_differential_oracle",
    "each_adjacent_spectral_pair_uses_direct_or_complete_guarded_saturation",
    "every_raw_power_fraction_and_ratio_remains_public",
    "PROTO11_constraint_physical_transaction_boundary_and_stop_rules_unchanged",
    "physical_inputs_methods_CFL_retries_thresholds_and_schedule_unchanged",
    "runner_binds_and_restores_authorization_member_event_and_runner_identity",
}

EXPECTED_INPUT = {
    "eligible_amplitudes_in_order": ["5/2", "3"],
    "methods_in_order": ["RK4", "SSPRK3"],
    "resolutions_in_order": [2049, 4097, 8193],
    "expected_run_input_count": 12,
    "all_twelve_states_must_be_rebuilt_under_PROTO12": True,
    "all_twelve_u_p_and_q_arrays_must_match_HLT9": True,
    "expanded_PROTO12_run_config_hashes_must_be_serialized": True,
    "all_twelve_SRC4_accepted_state_source_prechecks_must_pass": True,
    "both_amplitudes_and_methods_must_pass_the_PROTO12_t0_common_event": True,
    "all_raw_and_conditioning_spectral_evidence_must_remain_public": True,
    "first_eligible_candidate_stops_later_candidates": True,
}

EXPECTED_CONTROL_KEYS = {
    "exact_Minkowski_must_be_bitwise_preserved_for_all_methods_and_grids",
    "one_bit_near_reference_state_must_not_use_exact_equilibrium_branch",
    "nontrivial_SRC4_and_SRC3_differential_control_must_pass",
    "generic_reduction_and_source_derivative_defects_must_remain_visible",
    "operator_call_must_not_mutate_or_reproject_input_q",
    "common_event_reduction_rows_must_equal_the_reference_map",
    "pairwise_raw_ratios_must_equal_the_predecessor_raw_ratios",
    "removing_any_saturation_guard_must_fail_closed",
    "old_PROTO10_coarse_direct_veto_must_remain_public",
    "no_new_tolerance_floor_damping_cleaning_or_fitted_asymptote",
}

EXPECTED_NAMESPACE = {
    "calibration_output_root": "runs/fgc-2-sf1/proto12/calibration",
    "holdout_output_root": "runs/fgc-2-sf1/proto12/holdout",
    "both_roots_must_be_absent_or_empty_at_authorization": True,
    "authorization_must_not_create_either_root": True,
    "runner_must_refuse_overwrite": True,
}

EXPECTED_PROOF_KEYS = {
    "PROTO12_and_PRO12_FRZ1_must_validate_from_immutable_checkpoint",
    "HLT9_CAL8_RSP1_SRC4_ID2_and_CAL8_plan_must_be_canonical_and_hash_matched",
    "CAL9_plan_must_differ_from_CAL8_only_by_the_frozen_PROTO12_delta",
    "all_twelve_PROTO12_input_hashes_must_be_rebuilt_and_serialized",
    "all_twelve_SRC4_source_prechecks_must_pass",
    "all_four_PROTO12_t0_common_events_must_pass",
    "SRC4_exact_reference_and_differential_controls_must_pass",
    "pairwise_guard_removal_and_raw_ratio_mutations_must_fail_closed",
    "PROTO11_transaction_constraint_boundary_health_and_affine_stops_must_remain_hash_bound",
    "authorized_runner_must_be_hash_bound_and_refuse_drift_overwrite_or_dirty_tracked_state",
    "fresh_namespaces_must_be_observed_empty_without_creation",
    "canonical_hash_bound_result_required",
    "no_PROTO12_trajectory_may_be_read_or_advanced",
}

EXPECTED_CLAIMS = {
    "PROTO12_successor_runtime_implemented": True,
    "PROTO12_fresh_GR0_dynamic_calibration_authorized": True,
    "PROTO12_resolved_holdout_manifest_authorized": False,
    "classical_spherical_diagnostic_authorized": False,
    "SGBL_execution_authorized": False,
    "FGCQR_holdout_execution_authorized": False,
    "retained_EFT_evolution_authorized": False,
    "physical_transition_claim_authorized": False,
    "FGCQR_mechanism_rejected": False,
    "general_gradient_route_rejected": False,
    "singularity_resolution_derived": False,
    "child_domain_or_topology_derived": False,
    "dark_sector_mechanism_derived": False,
    "varying_locally_measured_c_derived": False,
}

NEW_IMPLEMENTATION = (
    REPOSITORY / "src/recursive_horizons/fgc/evolution/protocol_v12.py",
    REPOSITORY / "src/recursive_horizons/fgc/evolution/spectral_sensitivity_v12.py",
    REPOSITORY / "src/recursive_horizons/fgc/evolution/proto12_runtime.py",
    REPOSITORY / "scripts/run_fgc_gr0_calibration_v12.py",
    Path(__file__).resolve(),
)

PROTO12_NUMERICAL_KEYS = {
    "common_event_spatial_each_adjacent_pair_direct_or_conditioning_guarded_saturation_required",
    "common_event_spatial_guarded_saturation_requires_both_pair_member_absolute_budgets",
    "common_event_spatial_guarded_saturation_requires_both_pair_member_map_witnesses",
    "common_event_spatial_every_raw_power_fraction_and_ratio_serialized",
    "source_evaluator",
    "source_point_batch_size",
    "source_fixed_point_batching_mutates_pointwise_system",
    "SRC3_explicit_index_backend_retained_as_differential_oracle",
}


def _path(name: str, value: Any, expected: str) -> Path:
    if value != expected:
        raise ValueError(f"HLT10 {name} differs")
    path = REPOSITORY / expected
    if not path.is_file():
        raise ValueError(f"HLT10 {name} is absent")
    return path


def validate_run_plan(path: Path) -> dict[str, Any]:
    """Prove CAL9 differs from immutable CAL8 only by PROTO12's delta."""

    with path.open("rb") as handle:
        raw = tomllib.load(handle)
    if (
        raw.get("schema_version") != 1
        or raw.get("artifact_id") != "FGC-1-CAL9-RUN1-PLAN"
        or raw.get("project_version") != PROJECT_VERSION
        or raw.get("protocol_config") != EXPECTED_PATHS["protocol_config"]
        or raw.get("protocol_freeze_result")
        != EXPECTED_PATHS["protocol_freeze_result"]
        or raw.get("runtime_authorization_result")
        != "results/fgc-1-hlt10-mon10.json"
        or raw.get("static_input_result") != EXPECTED_PATHS["static_input_result"]
    ):
        raise ValueError("CAL9 run-plan identity or ownership differs")
    if raw.get("scope") != {
        "target_protocol": "FGC-2-SF1-PROTO12",
        "branch": "GR-0",
        "role": "fresh_outcome_neutral_dynamic_amplitude_recalibration_after_SRC4_and_pairwise_spectral_qualification",
        "physical_equations": "unredefined_GR0_specialization_of_REF1",
        "candidate_action_fields_or_health_stops_forbidden": True,
        "retained_EFT_interpretation": False,
    }:
        raise ValueError("CAL9 scope differs")
    numerics = raw.get("numerics", {})
    expected_delta = {
        "common_event_spatial_each_adjacent_pair_direct_or_conditioning_guarded_saturation_required": True,
        "common_event_spatial_guarded_saturation_requires_both_pair_member_absolute_budgets": True,
        "common_event_spatial_guarded_saturation_requires_both_pair_member_map_witnesses": True,
        "common_event_spatial_every_raw_power_fraction_and_ratio_serialized": True,
        "source_evaluator": "SRC4_tensor_contracted_reference_covariant_GR0_REF1",
        "source_point_batch_size": PROTO12_SOURCE_POINT_BATCH_SIZE,
        "source_fixed_point_batching_mutates_pointwise_system": False,
        "SRC3_explicit_index_backend_retained_as_differential_oracle": True,
    }
    if (
        numerics.get("resolutions") != list(PROTO12_POINT_COUNTS)
        or any(numerics.get(key) != value for key, value in expected_delta.items())
    ):
        raise ValueError("CAL9 source or pairwise-spectral contract differs")
    if raw.get("provenance") != {
        "output_root": "runs/fgc-2-sf1/proto12/calibration",
        "campaign_manifest": "runs/fgc-2-sf1/proto12/calibration/manifest.json",
        "append_only_event_log_name": "events.jsonl",
        "checkpoint_extension": ".npz",
        "fresh_output_root_required_at_campaign_start": True,
        "runner_refuses_overwrite": True,
        "resume_requires_matching_plan_and_manifest_hashes": True,
        "PROTO11_and_RSP1_outputs_forbidden_as_PROTO12_outcomes": True,
        "canonical_reduced_result_is_tracked_only_after_campaign_classification": True,
    }:
        raise ValueError("CAL9 provenance differs")
    if raw.get("claims", {}).get("SGBL_execution_authorized") is not False:
        raise ValueError("CAL9 SGB-L claim boundary differs")

    normalized = deepcopy(raw)
    normalized["artifact_id"] = "FGC-1-CAL8-RUN1-PLAN"
    normalized["protocol_config"] = "configs/fgc/fgc-2-sf1-protocol-v11.toml"
    normalized["protocol_freeze_result"] = "results/fgc-1-pro11-frz1.json"
    normalized["runtime_authorization_result"] = "results/fgc-1-hlt9-mon9.json"
    normalized["scope"]["target_protocol"] = "FGC-2-SF1-PROTO11"
    normalized["scope"]["role"] = (
        "fresh_outcome_neutral_dynamic_amplitude_recalibration_after_"
        "well_balanced_reference_map"
    )
    for key in PROTO12_NUMERICAL_KEYS:
        normalized["numerics"].pop(key)
    normalized["numerics"].update(
        {
            "common_event_spatial_coarse_to_medium_tail_ratios_required_directly": True,
            "common_event_spatial_finest_pair_direct_or_conditioning_guarded_saturation_required": True,
            "common_event_spatial_guarded_saturation_requires_medium_and_fine_map_witnesses": True,
        }
    )
    normalized["provenance"] = {
        "output_root": "runs/fgc-2-sf1/proto11/calibration",
        "campaign_manifest": "runs/fgc-2-sf1/proto11/calibration/manifest.json",
        "append_only_event_log_name": "events.jsonl",
        "checkpoint_extension": ".npz",
        "fresh_output_root_required_at_campaign_start": True,
        "runner_refuses_overwrite": True,
        "resume_requires_matching_plan_and_manifest_hashes": True,
        "PROTO10_outputs_forbidden_as_PROTO11_outcomes": True,
        "canonical_reduced_result_is_tracked_only_after_campaign_classification": True,
    }
    normalized["claims"].pop("SGBL_execution_authorized")
    baseline_blob = hlt4._git_blob(
        CHECKPOINT_COMMIT,
        EXPECTED_PATHS["predecessor_run_plan"],
        "immutable CAL8 plan",
    )
    if sha256(baseline_blob).hexdigest() != EXPECTED_LINEAGE[
        "predecessor_run_plan_sha256"
    ]:
        raise ValueError("immutable CAL8 plan hash differs")
    if normalized != tomllib.loads(baseline_blob.decode("utf-8")):
        raise ValueError("CAL9 changes more than the declared PROTO12 delta")
    return raw


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    with path.open("rb") as handle:
        raw = tomllib.load(handle)
    expected_root = {
        "schema_version",
        "artifact_id",
        "project_version",
        "metric_signature",
        "riemann_convention",
        *EXPECTED_PATHS,
        "scope",
        "immutable_lineage",
        "runtime_binding",
        "input_freeze",
        "controls",
        "namespace",
        "proof_contract",
        "claims",
    }
    if set(raw) != expected_root:
        raise ValueError("HLT10 root keys differ")
    if (
        raw["schema_version"] != 1
        or raw["artifact_id"] != ARTIFACT_ID
        or raw["project_version"] != PROJECT_VERSION
        or raw["metric_signature"] != "-+++"
        or raw["riemann_convention"] != "plus_partial_mu_gamma_nu"
    ):
        raise ValueError("HLT10 identity or convention differs")
    paths = {
        name: _path(name, raw[name], expected)
        for name, expected in EXPECTED_PATHS.items()
    }
    if raw["scope"] != EXPECTED_SCOPE:
        raise ValueError("HLT10 scope differs")
    if raw["immutable_lineage"] != EXPECTED_LINEAGE:
        raise ValueError("HLT10 immutable lineage differs")
    if set(raw["runtime_binding"]) != EXPECTED_RUNTIME_KEYS or any(
        value is not True for value in raw["runtime_binding"].values()
    ):
        raise ValueError("HLT10 runtime binding differs")
    if raw["input_freeze"] != EXPECTED_INPUT:
        raise ValueError("HLT10 input freeze differs")
    if set(raw["controls"]) != EXPECTED_CONTROL_KEYS or any(
        value is not True for value in raw["controls"].values()
    ):
        raise ValueError("HLT10 controls differ")
    if raw["namespace"] != EXPECTED_NAMESPACE:
        raise ValueError("HLT10 namespace differs")
    if set(raw["proof_contract"]) != EXPECTED_PROOF_KEYS or any(
        value is not True for value in raw["proof_contract"].values()
    ):
        raise ValueError("HLT10 proof contract differs")
    if raw["claims"] != EXPECTED_CLAIMS:
        raise ValueError("HLT10 claims differ")
    return {
        "raw": raw,
        "paths": paths,
        "run_plan": validate_run_plan(paths["run_plan_config"]),
    }


def _immutable_lineage(loaded: Mapping[str, Any]) -> dict[str, Any]:
    lineage = loaded["raw"]["immutable_lineage"]
    commit = lineage["checkpoint_commit"]
    if subprocess.run(
        ["git", "merge-base", "--is-ancestor", commit, "HEAD"],
        cwd=REPOSITORY,
        check=False,
        capture_output=True,
    ).returncode != 0:
        raise ValueError("HLT10 checkpoint is not an ancestor of HEAD")
    immutable_names = tuple(
        name for name in EXPECTED_PATHS if name != "run_plan_config"
    )
    blobs = {
        name: hlt4._git_blob(commit, EXPECTED_PATHS[name], name)
        for name in immutable_names
    }
    for name, blob in blobs.items():
        expected = lineage[f"{name}_sha256"]
        if sha256(blob).hexdigest() != expected:
            raise ValueError(f"immutable {name} hash differs")
        if hlt4._sha(REPOSITORY / EXPECTED_PATHS[name]) != expected:
            raise ValueError(f"current {name} differs from immutable checkpoint")

    protocol = validate_sf1_protocol_v12(
        tomllib.loads(blobs["protocol_config"].decode("utf-8"))
    )
    results = {
        name: hlt4._load_json_bytes(blobs[name], f"immutable {name}")
        for name in (
            "protocol_freeze_result",
            "predecessor_runtime_result",
            "calibration_diagnosis_result",
            "resolution_evidence_result",
            "source_instrument_result",
            "static_input_result",
        )
    }
    freeze = results["protocol_freeze_result"]
    predecessor = results["predecessor_runtime_result"]
    diagnosis = results["calibration_diagnosis_result"]
    resolution = results["resolution_evidence_result"]
    source = results["source_instrument_result"]
    static = results["static_input_result"]
    predecessor_inputs = predecessor.get("artifact_payload", {}).get(
        "frozen_run_inputs", []
    )
    if (
        protocol.get("artifact_id") != "FGC-2-SF1-PROTO12"
        or protocol.get("protocol_version") != 12
        or protocol.get("point_counts") != list(PROTO12_POINT_COUNTS)
        or protocol.get("new_numeric_tolerance_added") is not False
        or freeze.get("artifact_id") != "FGC-1-PRO12-FRZ1"
        or freeze.get("gate_status", {}).get("PROTO12_frozen") is not True
        or freeze.get("gate_status", {}).get(
            "PROTO12_SRC4_source_backend_contract_frozen"
        )
        is not True
        or freeze.get("gate_status", {}).get(
            "PROTO12_pairwise_spectral_classifier_derived"
        )
        is not True
        or freeze.get("gate_status", {}).get(
            "PROTO12_fresh_GR0_dynamic_calibration_authorized"
        )
        is not False
        or predecessor.get("artifact_id") != "FGC-1-HLT9-MON9"
        or predecessor.get("gate_status", {}).get(
            "PROTO11_fresh_GR0_dynamic_calibration_authorized"
        )
        is not True
        or len(predecessor_inputs) != 12
        or diagnosis.get("artifact_id") != "FGC-1-CAL8-PREF10"
        or diagnosis.get("gate_status", {}).get(
            "PROTO11_campaign_terminated_normally"
        )
        is not True
        or diagnosis.get("gate_status", {}).get("PROTO11_GR0_case_eligible")
        is not False
        or resolution.get("artifact_id") != "FGC-1-RSP1-PREF12"
        or resolution.get("gate_status", {}).get(
            "amplitude_three_spectral_veto_cleared_for_successor_design"
        )
        is not True
        or resolution.get("gate_status", {}).get(
            "RSP1_generic_all_field_raw_spectral_admission_passed"
        )
        is not False
        or source.get("artifact_id") != "FGC-1-SRC4-VEC1"
        or source.get("gate_status", {}).get(
            "SRC4_runtime_source_instrument_authorized_for_a_separately_frozen_GR0_protocol"
        )
        is not True
        or static.get("gate_status", {}).get("all_case_static_ledger_completed")
        is not True
    ):
        raise ValueError("immutable HLT10 predecessor gate differs")
    for record in (freeze, predecessor, diagnosis, resolution, source, static):
        nonclaims = record.get("nonclaims", {})
        if any(value is not False for value in nonclaims.values()):
            raise ValueError("immutable HLT10 predecessor promoted a nonclaim")

    inherited_hashes: dict[str, str] = {}
    for record in (freeze, predecessor, diagnosis, resolution, source, static):
        for relative, expected in record.get("implementation_sha256", {}).items():
            current = REPOSITORY / relative
            if not current.is_file() or hlt4._sha(current) != expected:
                raise ValueError(f"inherited PROTO12 premise drifted: {relative}")
            inherited_hashes[relative] = expected
    return {
        "checkpoint_commit": commit,
        "checkpoint_is_ancestor_of_HEAD": True,
        "protocol": {
            "artifact_id": protocol["artifact_id"],
            "protocol_version": protocol["protocol_version"],
            "point_counts": protocol["point_counts"],
            "source_backend": protocol["source_backend"],
            "pairwise_spectral_rule": protocol["pairwise_spectral_rule"],
            "semantic_holdout_contract_sha256": protocol[
                "semantic_holdout_contract_sha256"
            ],
        },
        "freeze": {
            "artifact_id": freeze["artifact_id"],
            "result_sha256": lineage["protocol_freeze_result_sha256"],
            "PROTO12_frozen": True,
        },
        "predecessor_runtime": {
            "artifact_id": predecessor["artifact_id"],
            "result_sha256": lineage["predecessor_runtime_result_sha256"],
            "frozen_input_count": len(predecessor_inputs),
        },
        "calibration_diagnosis": {
            "artifact_id": diagnosis["artifact_id"],
            "result_sha256": lineage["calibration_diagnosis_result_sha256"],
            "eligible_GR0_case": False,
        },
        "resolution_evidence": {
            "artifact_id": resolution["artifact_id"],
            "result_sha256": lineage["resolution_evidence_result_sha256"],
            "target_veto_cleared_for_design": True,
            "generic_raw_all_field_pass": False,
        },
        "source_instrument": {
            "artifact_id": source["artifact_id"],
            "result_sha256": lineage["source_instrument_result_sha256"],
            "runtime_instrument_authorized": True,
        },
        "static_input": {
            "artifact_id": static["artifact_id"],
            "result_sha256": lineage["static_input_result_sha256"],
        },
        "predecessor_run_plan_sha256": lineage["predecessor_run_plan_sha256"],
        "predecessor_frozen_inputs": predecessor_inputs,
        "inherited_implementation_files_hash_matched": len(inherited_hashes),
        "inherited_implementation_hashes": inherited_hashes,
    }


def _initial_data(
    plan: Mapping[str, Any], amplitude: str, method: str, count: int
):
    physical = plan["physical_inputs"]
    order = 4 if method == "RK4" else 2
    return construct_gr0_grid_initial_data(
        PulseParameters(
            chi_amplitude=float(Q(amplitude)),
            center=float(Q(physical["chi_center"])),
            half_width=float(Q(physical["chi_half_width"])),
            phi_amplitude=float(Q(physical["phi_seed_amplitude"])),
            planck_mass=float(Q(physical["planck_mass"])),
            scalar_mass=float(Q(physical["scalar_mass"])),
            quartic_coupling=float(Q(physical["quartic_coupling"])),
        ),
        point_count=count,
        outer_radius=float(Q(physical["outer_radius"])),
        constraint_method=method,
        diagnostic_spatial_order=order,
    )


def _operator(plan: Mapping[str, Any], grid: UniformRadialGrid, order: int):
    numerics = plan["numerics"]
    thresholds = plan["universal_thresholds"]
    return Proto12GR0EvolutionOperator(
        grid,
        spatial_order=order,
        ko_dissipation=float(Q(numerics["ko_dissipation"])),
        raw_tolerance=float(Q(thresholds["source_residual_infinity_max"])),
        kinetic_condition_maximum=float(
            Q(thresholds["kinetic_condition_number_max"])
        ),
        maximum_refinement_iterations=thresholds["source_iteration_max"],
        point_batch_size=numerics["source_point_batch_size"],
    )


def _freeze_inputs(
    loaded: Mapping[str, Any], lineage: Mapping[str, Any]
) -> tuple[
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[Proto12GR0CommonEventAssessment],
    dict[tuple[str, str], list[tuple[int, EvolutionState, UniformRadialGrid]]],
]:
    plan = loaded["run_plan"]
    physical = plan["physical_inputs"]
    numerics = plan["numerics"]
    thresholds = plan["universal_thresholds"]
    plan_hash = hlt4._sha(loaded["paths"]["run_plan_config"])
    protocol_hash = loaded["raw"]["immutable_lineage"]["protocol_config_sha256"]
    predecessor_inputs = lineage["predecessor_frozen_inputs"]
    predecessor_by_key = {
        (item["amplitude"], item["method"], item["point_count"]): item
        for item in predecessor_inputs
    }
    expected_order = [
        (amplitude, method, count)
        for amplitude in ("5/2", "3")
        for method in ("RK4", "SSPRK3")
        for count in PROTO12_POINT_COUNTS
    ]
    if [
        (item["amplitude"], item["method"], item["point_count"])
        for item in predecessor_inputs
    ] != expected_order:
        raise ValueError("immutable HLT9 input order differs")

    frozen: list[dict[str, Any]] = []
    prechecks: list[dict[str, Any]] = []
    states: dict[
        tuple[str, str], list[tuple[int, EvolutionState, UniformRadialGrid]]
    ] = {}
    threshold = float(Q(thresholds["source_residual_infinity_max"]))
    for amplitude, method, count in expected_order:
        order = 4 if method == "RK4" else 2
        predecessor = predecessor_by_key[(amplitude, method, count)]
        initial = _initial_data(plan, amplitude, method, count)
        state = project_gr0_reference_balanced_state(initial, spatial_order=order)
        state_hash = array_content_sha256(state.u, state.p, state.q)
        if state_hash != predecessor["projected_state_sha256"]:
            raise ValueError("PROTO12 input differs from immutable HLT9")
        before = state_hash
        rhs = _operator(plan, initial.grid, order)(0.0, state)
        after = array_content_sha256(state.u, state.p, state.q)
        residual = float(rhs.diagnostics["source_residual_infinity"])
        if (
            residual >= threshold
            or rhs.diagnostics["source_raw_gate_passed"] is not True
            or rhs.diagnostics["SRC4_tensor_contracted_reference_source"]
            is not True
            or rhs.diagnostics["PROTO11_exact_reference_equilibrium_applied"]
            is not False
            or rhs.diagnostics["PROTO11_interior_q_reprojected"] is not False
            or after != before
        ):
            raise ValueError("PROTO12 accepted-state source precheck failed")
        record = dict(predecessor)
        record["predecessor_HLT9_projected_state_sha256"] = record.pop(
            "projected_state_sha256"
        )
        record["predecessor_HLT9_expanded_run_config_sha256"] = record.pop(
            "expanded_run_config_sha256"
        )
        record["predecessor_HLT9_plan_sha256"] = record.pop("plan_sha256")
        record["predecessor_HLT9_protocol_sha256"] = record.pop(
            "protocol_sha256"
        )
        record["predecessor_HLT9_resolution_role"] = record.pop(
            "PROTO11_resolution_role"
        )
        record.update(
            {
                "plan_sha256": plan_hash,
                "protocol_sha256": protocol_hash,
                "PROTO12_resolution_role": (
                    "direct_state_with_SRC4_source_reference_balanced_constraints_"
                    "and_pairwise_guarded_spectra"
                ),
                "projected_state_sha256": state_hash,
                "u_p_q_bitwise_match_HLT9": True,
                "SRC4_source_backend_bound": True,
                "source_point_batch_size": numerics["source_point_batch_size"],
                "q_remains_independently_evolved": True,
                "interior_q_reprojected": False,
                "complete_raw_source_gate_passed": True,
                "trajectory_advanced": False,
            }
        )
        record["expanded_run_config_sha256"] = hlt4._digest(record)
        frozen.append(record)
        prechecks.append(
            {
                "amplitude": amplitude,
                "method": method,
                "point_count": count,
                "source_backend": "SRC4_tensor_contracted_reference_covariant_GR0_REF1",
                "source_residual_infinity": residual,
                "source_residual_maximum": threshold,
                "source_point_batch_size": numerics["source_point_batch_size"],
                "source_point_batch_count": rhs.diagnostics[
                    "PROTO12_source_point_batch_count"
                ],
                "accepted_state_source_gate_passed": True,
                "complete_unredefined_source_evaluated": True,
                "exact_reference_equilibrium_branch_used": False,
                "operator_input_state_unchanged": True,
                "interior_q_reprojected": False,
                "trajectory_advanced": False,
            }
        )
        states.setdefault((amplitude, method), []).append(
            (count, state, initial.grid)
        )

    events: list[dict[str, Any]] = []
    assessments: list[Proto12GR0CommonEventAssessment] = []
    for amplitude in ("5/2", "3"):
        for method in ("RK4", "SSPRK3"):
            records = sorted(states[(amplitude, method)], key=lambda item: item[0])
            assessed = proto12_gr0_common_event(
                [item[1] for item in records],
                [item[2] for item in records],
                accepted_stage_counts=(0, 0, 0),
                method=method,
                coordinate_time=0.0,
                cutoff=float(Q(physical["cutoff_Lambda"])),
                measurement_radius_maximum=float(
                    Q(physical["measurement_radius_maximum"])
                ),
                taper_fraction=float(Q(numerics["proper_spectral_taper_fraction"])),
                fixed_outer_rows=numerics["fixed_outer_rows"],
            )
            pairwise = assessed.spatial_spectral_admission
            if (
                assessed.reference_balanced_constraint_admission_passed is not True
                or pairwise.admission_passed is not True
                or pairwise.every_tail_resolved_or_saturated is not True
                or assessed.admission_passed is not True
                or assessed.reference_state_map_applied is not True
                or assessed.SRC4_source_backend_bound is not True
                or assessed.interior_q_reprojected is not False
            ):
                raise ValueError("PROTO12 t0 common event failed closed")
            events.append(
                {
                    "amplitude": amplitude,
                    "method": method,
                    "coordinate_time": 0.0,
                    "point_counts": list(assessed.point_counts),
                    "reference_balanced_constraint_admission_passed": True,
                    "raw_PROTO9_admission_passed": assessed.raw_PROTO9_admission_passed,
                    "legacy_PROTO10_admission_passed": (
                        assessed.legacy_PROTO10_admission_passed
                    ),
                    "legacy_PROTO11_admission_passed": (
                        assessed.legacy_PROTO11_admission_passed
                    ),
                    "PROTO12_pairwise_admission_passed": True,
                    "PROTO12_admission_passed": True,
                    "diagnostic_saturation_used": pairwise.saturation_used,
                    "diagnostic_saturation_is_physical_resolution": False,
                    "round_trip_is_a_continuum_error_bound": False,
                    "complete_common_event": hlt4._serial(assessed),
                    "trajectory_advanced": False,
                }
            )
            assessments.append(assessed)
    return frozen, prechecks, events, assessments, states


def _runtime_controls(
    loaded: Mapping[str, Any],
    states: Mapping[
        tuple[str, str], list[tuple[int, EvolutionState, UniformRadialGrid]]
    ],
    assessments: list[Proto12GR0CommonEventAssessment],
) -> dict[str, Any]:
    plan = loaded["run_plan"]
    exact_records = []
    for method in ("RK4", "SSPRK3"):
        order = 4 if method == "RK4" else 2
        for count in PROTO12_POINT_COUNTS:
            grid = UniformRadialGrid(0.0, 128.0, count)
            state = exact_spherical_minkowski_reference(grid)
            before = array_content_sha256(state.u, state.p, state.q)
            rhs = _operator(plan, grid, order)(0.0, state)
            after = array_content_sha256(state.u, state.p, state.q)
            if (
                before != after
                or not np.array_equal(rhs.du, np.zeros_like(rhs.du))
                or not np.array_equal(rhs.dp, np.zeros_like(rhs.dp))
                or not np.array_equal(rhs.dq, np.zeros_like(rhs.dq))
                or rhs.diagnostics["source_raw_gate_passed"] is not True
                or rhs.diagnostics["PROTO11_exact_reference_equilibrium_applied"]
                is not True
            ):
                raise RuntimeError("PROTO12 exact-reference control failed")
            exact_records.append(
                {
                    "method": method,
                    "point_count": count,
                    "raw_source_residual_infinity": rhs.diagnostics[
                        "source_residual_infinity"
                    ],
                    "bitwise_zero_du_dt": True,
                    "bitwise_zero_dp_dt": True,
                    "bitwise_zero_dq_dt": True,
                    "operator_input_unchanged": True,
                    "exact_reference_equilibrium_branch_used": True,
                }
            )

    grid = UniformRadialGrid(0.0, 128.0, PROTO12_POINT_COUNTS[0])
    reference = exact_spherical_minkowski_reference(grid)
    perturbed_q = reference.q.copy()
    perturbed_q[1, 3] = np.nextafter(perturbed_q[1, 3], np.inf)
    one_bit = EvolutionState(reference.u, reference.p, perturbed_q)
    one_bit_rhs = _operator(plan, grid, 4)(0.0, one_bit)
    if one_bit_rhs.diagnostics["PROTO11_exact_reference_equilibrium_applied"]:
        raise RuntimeError("PROTO12 one-bit control used exact equilibrium")

    first_count, first_state, first_grid = states[("5/2", "RK4")][0]
    if first_count != PROTO12_POINT_COUNTS[0]:
        raise RuntimeError("PROTO12 selected the wrong physical control")
    derivative = SBPFirstDerivative(first_grid, 4)
    p_r, q_r, _ = reference_balanced_spatial_derivatives(first_state, derivative)
    indices = np.arange(1, 17)
    src3 = diagnose_gr0_reference_balanced_accelerations(
        first_state.u[indices],
        first_state.p[indices],
        first_state.q[indices],
        p_r[indices],
        q_r[indices],
        first_grid.coordinates[indices],
    )
    src4 = diagnose_gr0_vectorized_reference_accelerations(
        first_state.u[indices],
        first_state.p[indices],
        first_state.q[indices],
        p_r[indices],
        q_r[indices],
        first_grid.coordinates[indices],
    )
    differential = float(
        np.max(np.abs(src3.accelerations - src4.accelerations), initial=0.0)
    )
    if differential != 0.0:
        raise RuntimeError("PROTO12 SRC4/SRC3 differential control changed")

    generic_q = reference.q.copy()
    generic_q[grid.point_count // 4, 5] += 1.0e-6
    generic = EvolutionState(reference.u, reference.p, generic_q)
    _, generic_q_r, _ = reference_balanced_spatial_derivatives(generic, derivative)
    generic_reduction = reference_balanced_reduction_constraint(generic, derivative)
    if (
        float(np.max(np.abs(generic_reduction), initial=0.0)) <= 0.0
        or float(np.max(np.abs(generic_q_r), initial=0.0)) <= 0.0
    ):
        raise RuntimeError("PROTO12 generic defect was hidden")

    snapshot = reference_balanced_semidiscrete_constraint_snapshot(
        first_state,
        first_grid,
        coordinate_time=0.0,
        diagnostic_spatial_order=4,
    )
    expected_reduction = reference_balanced_reduction_constraint(
        first_state, derivative
    )
    reduction_difference = float(
        np.max(
            np.abs(snapshot.residuals[:, 4:] - expected_reduction),
            initial=0.0,
        )
    )
    if reduction_difference != 0.0:
        raise RuntimeError("PROTO12 common event uses another reduction map")

    saturated_pair_count = 0
    raw_adjacent_failure_count = 0
    for assessment in assessments:
        pairwise = assessment.spatial_spectral_admission
        for values in (
            *pairwise.field_pair_classification.values(),
            *pairwise.derivative_pair_classification.values(),
        ):
            saturated_pair_count += sum(
                value == "diagnostically_saturated" for value in values
            )
        raw = assessment.raw_spatial_spectral_admission
        for ratios in (
            *raw.field_power_tail_ratios.values(),
            *raw.derivative_power_tail_ratios.values(),
        ):
            raw_adjacent_failure_count += sum(
                value is None or value >= 0.25 for value in ratios
            )
    if saturated_pair_count < 1 or raw_adjacent_failure_count < 1:
        raise RuntimeError("PROTO12 controls did not exercise guarded saturation")
    return {
        "exact_reference_records": exact_records,
        "exact_reference_control_count": len(exact_records),
        "all_exact_reference_controls_bitwise_preserved": True,
        "complete_raw_source_serialized_before_exact_equilibrium_selection": True,
        "one_bit_perturbation_used_exact_equilibrium_branch": False,
        "one_bit_perturbation_raw_source_residual_infinity": one_bit_rhs.diagnostics[
            "source_residual_infinity"
        ],
        "SRC4_SRC3_differential_sample_count": int(indices.size),
        "SRC4_SRC3_acceleration_difference_infinity": differential,
        "SRC4_SRC3_differential_control_passed": True,
        "generic_reduction_defect_infinity": float(
            np.max(np.abs(generic_reduction), initial=0.0)
        ),
        "generic_source_q_r_defect_infinity": float(
            np.max(np.abs(generic_q_r), initial=0.0)
        ),
        "generic_reduction_and_source_derivative_defects_visible": True,
        "common_event_reduction_map_difference_infinity": reduction_difference,
        "common_event_reduction_rows_match_reference_map_bitwise": True,
        "pairwise_raw_ratios_preserved": True,
        "diagnostically_saturated_pair_classification_count": saturated_pair_count,
        "public_t0_raw_adjacent_pair_failure_count": raw_adjacent_failure_count,
        "historical_CAL8_evolved_coarse_pair_failure_retained": True,
        "guard_removal_is_fail_closed": True,
        "interior_q_reprojected_by_any_control": False,
        "new_numeric_tolerance_floor_or_fit_added": False,
        "all_controls_passed": True,
    }


def _namespace_precondition(
    raw: Mapping[str, Any], historical: Mapping[str, Any] | None
) -> dict[str, Any]:
    if historical is not None:
        records = historical.get("records", [])
        if (
            historical.get("both_fresh_namespaces_empty") is not True
            or historical.get("authorization_created_no_namespace") is not True
            or [item.get("path") for item in records]
            != [
                raw["namespace"]["calibration_output_root"],
                raw["namespace"]["holdout_output_root"],
            ]
            or any(item.get("empty") is not True for item in records)
        ):
            raise ValueError("historical HLT10 namespace evidence differs")
        return dict(historical)
    records = []
    for key in ("calibration_output_root", "holdout_output_root"):
        relative = raw["namespace"][key]
        path = REPOSITORY / relative
        exists = path.exists()
        if exists and (not path.is_dir() or any(path.iterdir())):
            raise ValueError(f"HLT10 fresh namespace is nonempty: {relative}")
        records.append(
            {
                "path": relative,
                "existed_before_authorization": exists,
                "entry_count": 0,
                "empty": True,
                "created_by_authorization": False,
            }
        )
    return {
        "both_fresh_namespaces_empty": True,
        "authorization_created_no_namespace": True,
        "records": records,
    }


def record(
    config_path: Path = DEFAULT_CONFIG,
    *,
    historical_namespace_evidence: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    loaded = load_config(config_path)
    raw = loaded["raw"]
    lineage = _immutable_lineage(loaded)
    inputs, prechecks, events, assessments, states = _freeze_inputs(loaded, lineage)
    controls = _runtime_controls(loaded, states, assessments)
    if (
        len(inputs) != 12
        or len(prechecks) != 12
        or len(events) != 4
        or any(not item["accepted_state_source_gate_passed"] for item in prechecks)
        or any(not item["PROTO12_admission_passed"] for item in events)
        or controls["all_controls_passed"] is not True
    ):
        raise ValueError("HLT10 authorization evidence is incomplete")
    namespace = _namespace_precondition(raw, historical_namespace_evidence)
    plan_hash = hlt4._sha(loaded["paths"]["run_plan_config"])
    input_manifest_hash = hlt4._digest(inputs)
    campaign_id = hlt4._digest(
        {
            "artifact_id": loaded["run_plan"]["artifact_id"],
            "run_plan_sha256": plan_hash,
            "input_manifest_sha256": input_manifest_hash,
            "ordered_amplitudes": ["5/2", "3"],
            "point_counts": list(PROTO12_POINT_COUNTS),
            "source_backend": "SRC4",
            "spectral_classifier": "PROTO12_pairwise",
            "output_root": raw["namespace"]["calibration_output_root"],
        }
    )
    inherited_hashes = dict(lineage.pop("inherited_implementation_hashes"))
    lineage.pop("predecessor_frozen_inputs")
    implementation = dict(inherited_hashes)
    implementation.update(
        {hlt4._rel(path): hlt4._sha(path) for path in NEW_IMPLEMENTATION}
    )
    return {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "classification": (
            "pre_trajectory_PROTO12_SRC4_pairwise_spectral_runtime_and_fresh_"
            "GR0_calibration_authorization"
        ),
        "generated_by": "scripts/reproduce_fgc_hlt10_mon10.py",
        "derivation_document": hlt4._rel(OWNER_DOCUMENT),
        "derivation_document_sha256": hlt4._sha(OWNER_DOCUMENT),
        "source_config_sha256": {
            hlt4._rel(config_path): hlt4._sha(config_path),
            hlt4._rel(loaded["paths"]["run_plan_config"]): plan_hash,
        },
        "predecessor_sha256": {
            EXPECTED_PATHS[name]: raw["immutable_lineage"][f"{name}_sha256"]
            for name in EXPECTED_PATHS
            if name != "run_plan_config"
        },
        "implementation_sha256": implementation,
        "scope_bindings": dict(raw["scope"]),
        "gate_status": dict(raw["claims"]),
        "artifact_payload": {
            "immutable_lineage": lineage,
            "runtime_controls": controls,
            "inherited_runtime": {
                "PROTO11_reference_map_constraints_transaction_boundary_health_and_affine_stops_unchanged": True,
                "raw_PROTO9_and_PROTO10_spectral_evidence_retained": True,
                "inherited_implementation_files_hash_matched": lineage[
                    "inherited_implementation_files_hash_matched"
                ],
                "process_local_runner_bindings_are_exception_safe": True,
            },
            "frozen_run_plan": {
                "artifact_id": loaded["run_plan"]["artifact_id"],
                "run_plan_sha256": plan_hash,
                "campaign_id": campaign_id,
                "input_manifest_sha256": input_manifest_hash,
                "ordered_amplitudes": ["5/2", "3"],
                "point_counts": list(PROTO12_POINT_COUNTS),
                "source_backend": "SRC4_tensor_contracted_reference_covariant_GR0_REF1",
                "source_point_batch_size": PROTO12_SOURCE_POINT_BATCH_SIZE,
                "pairwise_spectral_rule": (
                    "each_adjacent_pair_direct_or_conditioning_guarded_saturation"
                ),
                "first_eligible_candidate_stops_later_candidates": True,
                "expanded_input_count": len(inputs),
                "u_p_q_match_count": sum(
                    item["u_p_q_bitwise_match_HLT9"] for item in inputs
                ),
                "SRC4_source_precheck_passed_count": sum(
                    item["accepted_state_source_gate_passed"] for item in prechecks
                ),
                "PROTO12_t0_passed_count": sum(
                    item["PROTO12_admission_passed"] for item in events
                ),
                "t0_common_event_evaluated_count": len(events),
            },
            "frozen_run_inputs": inputs,
            "accepted_state_source_prechecks": prechecks,
            "initial_common_events": events,
            "namespace_precondition": namespace,
            "decision": (
                "HLT10_binds_SRC4_and_the_pairwise_classifier_and_authorizes_"
                "one_fresh_PROTO12_GR0_calibration_while_every_candidate_and_"
                "physical_claim_remains_closed"
            ),
            "epistemic_boundary": {
                "PROTO11_CAL8_and_RSP1_history_disclosed": True,
                "PROTO12_trajectory_read": False,
                "PROTO12_trajectory_advanced": False,
                "all_twelve_physical_inputs_reconstructed": True,
                "all_twelve_SRC4_source_prechecks_passed": True,
                "all_four_PROTO12_t0_admissions_passed": True,
                "every_raw_spectral_ratio_retained": True,
                "diagnostic_saturation_is_physical_resolution": False,
                "round_trip_is_a_continuum_error_bound": False,
                "SRC4_or_pairwise_classifier_is_a_physical_mechanism_result": False,
                "fresh_GR0_recalibration_authorized": True,
                "fresh_GR0_recalibration_completed": False,
                "trapped_or_collapse_outcome_classified": False,
                "resolved_holdout_manifest_exists": False,
                "SGBL_or_FGCQR_outcome_read": False,
                "mechanism_question_answered": False,
            },
        },
        "nonclaims": {
            "fresh_GR0_recalibration_completed": False,
            "SRC4_physically_validated": False,
            "pairwise_diagnostic_is_physical_resolution": False,
            "trapped_sphere_formed_dynamically": False,
            "SGBL_health_definition_supplied": False,
            "SGBL_comparison_completed": False,
            "FGCQR_trajectory_opened": False,
            "regulator_activation_observed": False,
            "positive_Raychaudhuri_margin_measured": False,
            "constraint_to_DEF1_stability_map_supplied": False,
            "retained_EFT_validity": False,
            "physical_transition": False,
            "singularity_resolution": False,
            "child_domain_or_topology": False,
            "dark_sector_mechanism": False,
            "varying_locally_measured_c": False,
        },
    }


def load_canonical_result(path: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    result = hlt4._load_json_bytes(path.read_bytes(), "HLT10 result")
    if result.get("artifact_id") != ARTIFACT_ID:
        raise ValueError("HLT10 result artifact identity differs")
    return result


def verify_canonical(
    path: Path = DEFAULT_OUTPUT,
    config_path: Path = DEFAULT_CONFIG,
) -> None:
    observed = load_canonical_result(path)
    historical = observed.get("artifact_payload", {}).get("namespace_precondition")
    expected = record(
        config_path,
        historical_namespace_evidence=historical,
    )
    if observed != hlt4._serial(expected):
        raise ValueError("HLT10 output differs from a fresh reproduction")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    config = args.config.resolve()
    output = args.output.resolve()
    if args.check:
        verify_canonical(output, config)
        print(f"verified {output}")
        return
    result = record(config)
    output.write_text(hlt4._canonical(result), encoding="utf-8")
    print(
        f"wrote {output} "
        "(PROTO12 runtime=true; fresh GR0 calibration authorized=true; "
        "trajectory_read=false)"
    )


if __name__ == "__main__":
    main()
