#!/usr/bin/env python3
"""Reproduce the pre-trajectory FGC-1-HLT9-MON9 authorization record.

HLT9 binds the outcome-neutral PROTO11 reference-state differentiation map to
the existing GR-0 campaign machinery.  It rebuilds every frozen input, runs
the unchanged raw source gate, recomputes all four time-zero common events,
and attacks the exact-reference exception and the independently evolved
``q`` contract.  It never advances a trajectory or creates a run namespace.
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
from recursive_horizons.fgc.evolution.calibration_runtime import (  # noqa: E402
    project_gr0_semidiscrete_state,
)
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
    Proto11GR0CommonEventAssessment,
    Proto11GR0EvolutionOperator,
    exact_spherical_minkowski_reference,
    project_gr0_reference_balanced_state,
    proto11_gr0_common_event,
    reference_balanced_reduction_constraint,
    reference_balanced_semidiscrete_constraint_snapshot,
    reference_balanced_spatial_derivatives,
)
from recursive_horizons.fgc.evolution.proto9_runtime import (  # noqa: E402
    PROTO9_POINT_COUNTS,
)
from recursive_horizons.fgc.evolution.protocol_v11 import (  # noqa: E402
    validate_sf1_protocol_v11,
)
from recursive_horizons.fgc.evolution.vectorized_source import (  # noqa: E402
    ADM_CENTER_PARITIES,
)
from recursive_horizons.fgc.initial_data_preflight import (  # noqa: E402
    PulseParameters,
)


Q = hlt4.Q
ARTIFACT_ID = "FGC-1-HLT9-MON9"
PROJECT_VERSION = "0.11.0"
CHECKPOINT_COMMIT = "b53c0f1a38823b84590e3946fb31146be99ed7b5"
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-hlt9-mon9.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-hlt9-mon9.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-hlt9-mon9.md"

EXPECTED_PATHS = {
    "protocol_config": "configs/fgc/fgc-2-sf1-protocol-v11.toml",
    "protocol_freeze_config": "configs/fgc/fgc-1-pro11-frz1.toml",
    "protocol_freeze_result": "results/fgc-1-pro11-frz1.json",
    "predecessor_runtime_result": "results/fgc-1-hlt8-mon8.json",
    "diagnosis_result": "results/fgc-1-cal7-pref9.json",
    "static_input_result": "results/fgc-1-id2-all1.json",
    "predecessor_run_plan": "configs/fgc/fgc-1-cal7-run1.toml",
    "run_plan_config": "configs/fgc/fgc-1-cal8-run1.toml",
}

EXPECTED_SCOPE = {
    "target_protocol": "FGC-2-SF1-PROTO11",
    "runtime_owner": ARTIFACT_ID,
    "calibration_branch": "GR-0",
    "role": "pre_trajectory_reference_balanced_runtime_binding_and_fresh_input_freeze",
    "PROTO10_campaign_and_CAL7_diagnosis_history_disclosed": True,
    "PROTO11_trajectory_read": False,
    "SGBL_trajectory_read": False,
    "FGCQR_trajectory_read": False,
    "collapse_or_trapped_outcome_classified": False,
    "reference_map_promoted_to_physical_mechanism": False,
}

EXPECTED_LINEAGE = {
    "checkpoint_commit": CHECKPOINT_COMMIT,
    "protocol_config_sha256": "31d21bf6257d6afe033ea8a0d48b59360187a1d9baf2bc95c2fcb086e39ce440",
    "protocol_freeze_config_sha256": "4e3b84e51ffae75a97ba81fe7ed22d9af0c41c526f74f3f51f593a4f66ec9bbc",
    "protocol_freeze_result_sha256": "a6d3958859111fe967db493333e2c4df9146e96913a91de6813e0a4e5f6e260f",
    "predecessor_runtime_result_sha256": "108150a8a63f952ff4e873a5723ed3ecee23f5c7cda718b91621bcade840f1b5",
    "diagnosis_result_sha256": "aab4ea62d569e224bf0bff0f7c3368b2d8dd07007735c1d1993ea3f9ef5b5d2a",
    "static_input_result_sha256": "10e2661be521087467b970f1cd877311bd63a0cb6874a1ee745ca00ec74ab981",
    "predecessor_run_plan_sha256": "1f60e7134873edcbbcd80e149eeab01308cd900008c4753118674d7dfa3ff71d",
    "checkpoint_must_be_ancestor_of_HEAD": True,
    "all_predecessor_blobs_must_be_read_from_checkpoint": True,
}

EXPECTED_RUNTIME_KEYS = {
    "initial_q_uses_Dh_u_minus_uref_plus_qref",
    "source_qr_uses_Dh_q_minus_qref",
    "reduction_constraint_uses_same_reference_map",
    "common_event_constraints_use_same_reference_map",
    "p_r_and_q_time_evolution_remain_Dh_p",
    "q_remains_independent_evolved_field",
    "interior_q_reprojection_is_forbidden",
    "exact_reference_equilibrium_uses_bitwise_identity_not_tolerance",
    "complete_unredefined_source_and_raw_gate_remain_serialized",
    "inherited_PROTO10_spatial_spectral_compositor_byte_unchanged",
    "inherited_PROTO7_evolution_transaction_and_source_retry_byte_unchanged",
    "physical_inputs_equations_methods_CFL_retries_thresholds_and_stops_unchanged",
    "runner_binds_and_restores_authorization_member_event_and_runner_identity",
}

EXPECTED_INPUT = {
    "eligible_amplitudes_in_order": ["5/2", "3"],
    "methods_in_order": ["RK4", "SSPRK3"],
    "resolutions_in_order": [2049, 4097, 8193],
    "expected_run_input_count": 12,
    "all_twelve_states_must_be_rebuilt_under_PROTO11": True,
    "all_twelve_u_and_p_arrays_must_match_HLT8": True,
    "all_twelve_q_arrays_must_match_the_frozen_reference_map": True,
    "expanded_PROTO11_run_config_hashes_must_be_serialized": True,
    "all_twelve_accepted_state_source_prechecks_must_pass": True,
    "both_amplitudes_and_methods_must_pass_the_PROTO11_t0_common_event": True,
    "all_PROTO10_spectral_rules_and_fresh_raw_evidence_must_remain_public": True,
    "first_eligible_candidate_stops_later_candidates": True,
}

EXPECTED_CONTROL_KEYS = {
    "exact_Minkowski_must_be_bitwise_preserved_for_all_methods_and_grids",
    "one_bit_near_reference_state_must_not_use_exact_equilibrium_branch",
    "nontrivial_REF1_source_controls_must_pass_the_unchanged_raw_gate",
    "generic_reduction_defect_must_remain_visible",
    "source_derivative_defect_must_remain_visible",
    "operator_call_must_not_mutate_or_reproject_input_q",
    "common_event_reduction_rows_must_equal_the_reference_map",
    "old_direct_q_derivative_map_must_not_be_operative",
}

EXPECTED_NAMESPACE = {
    "calibration_output_root": "runs/fgc-2-sf1/proto11/calibration",
    "holdout_output_root": "runs/fgc-2-sf1/proto11/holdout",
    "both_roots_must_be_absent_or_empty_at_authorization": True,
    "authorization_must_not_create_either_root": True,
    "runner_must_refuse_overwrite": True,
}

EXPECTED_PROOF_KEYS = {
    "PROTO11_and_PRO11_FRZ1_must_validate_from_immutable_checkpoint",
    "HLT8_CAL7_ID2_and_CAL7_plan_must_be_canonical_and_hash_matched_from_immutable_checkpoint",
    "CAL8_plan_must_differ_from_CAL7_only_by_the_frozen_reference_map_and_namespace",
    "all_twelve_PROTO11_input_hashes_must_be_rebuilt_and_serialized",
    "all_twelve_t0_accepted_source_prechecks_must_pass",
    "all_four_PROTO11_t0_common_events_must_pass",
    "exact_reference_and_generic_defect_controls_must_pass",
    "interior_q_nonreprojection_must_be_attacked",
    "PROTO10_spectral_and_direct_coarse_veto_contract_must_remain_hash_bound",
    "inherited_transaction_source_boundary_health_and_affine_stops_must_remain_hash_bound",
    "authorized_runner_must_be_hash_bound_and_refuse_drift_overwrite_or_dirty_tracked_state",
    "fresh_namespaces_must_be_observed_empty_without_creation",
    "canonical_hash_bound_result_required",
    "no_PROTO11_trajectory_may_be_read_or_advanced",
}

EXPECTED_CLAIMS = {
    "PROTO11_successor_runtime_implemented": True,
    "PROTO11_fresh_GR0_dynamic_calibration_authorized": True,
    "PROTO11_resolved_holdout_manifest_authorized": False,
    "classical_spherical_diagnostic_authorized": False,
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
    REPOSITORY / "src/recursive_horizons/fgc/evolution/protocol_v11.py",
    REPOSITORY / "src/recursive_horizons/fgc/evolution/proto11_runtime.py",
    REPOSITORY / "scripts/run_fgc_gr0_calibration_v11.py",
    Path(__file__).resolve(),
)

REFERENCE_NUMERICS = {
    "initial_q_reference_state_rule": "D_h(u-u_ref)+q_ref",
    "source_q_r_reference_state_rule": "D_h(q-q_ref)",
    "reduction_constraint_reference_state_rule": "q-(D_h(u-u_ref)+q_ref)",
    "common_event_constraints_use_same_reference_state_rules": True,
    "q_remains_independently_evolved": True,
    "interior_q_reprojection_after_initialization_or_Runge_Kutta_stage": False,
    "exact_reference_equilibrium_requires_bitwise_state_identity": True,
    "exact_reference_equilibrium_adds_no_numeric_tolerance": True,
    "complete_unredefined_raw_source_residual_remains_serialized": True,
}


def _path(name: str, value: Any, expected: str) -> Path:
    if value != expected:
        raise ValueError(f"HLT9 {name} differs")
    path = REPOSITORY / expected
    if not path.is_file():
        raise ValueError(f"HLT9 {name} is absent")
    return path


def validate_run_plan(path: Path) -> dict[str, Any]:
    """Prove CAL8 differs from immutable CAL7 only by PROTO11's map."""

    with path.open("rb") as handle:
        raw = tomllib.load(handle)
    if (
        raw.get("schema_version") != 1
        or raw.get("artifact_id") != "FGC-1-CAL8-RUN1-PLAN"
        or raw.get("project_version") != PROJECT_VERSION
        or raw.get("protocol_config") != EXPECTED_PATHS["protocol_config"]
        or raw.get("protocol_freeze_result")
        != EXPECTED_PATHS["protocol_freeze_result"]
        or raw.get("runtime_authorization_result")
        != "results/fgc-1-hlt9-mon9.json"
        or raw.get("static_input_result") != EXPECTED_PATHS["static_input_result"]
    ):
        raise ValueError("CAL8 run-plan identity or ownership differs")
    if raw.get("scope") != {
        "target_protocol": "FGC-2-SF1-PROTO11",
        "branch": "GR-0",
        "role": "fresh_outcome_neutral_dynamic_amplitude_recalibration_after_well_balanced_reference_map",
        "physical_equations": "unredefined_GR0_specialization_of_REF1",
        "candidate_action_fields_or_health_stops_forbidden": True,
        "retained_EFT_interpretation": False,
    }:
        raise ValueError("CAL8 scope differs")
    numerics = raw.get("numerics", {})
    if (
        numerics.get("resolutions") != list(PROTO9_POINT_COUNTS)
        or any(numerics.get(key) != value for key, value in REFERENCE_NUMERICS.items())
    ):
        raise ValueError("CAL8 reference-map contract differs")
    if raw.get("provenance") != {
        "output_root": "runs/fgc-2-sf1/proto11/calibration",
        "campaign_manifest": "runs/fgc-2-sf1/proto11/calibration/manifest.json",
        "append_only_event_log_name": "events.jsonl",
        "checkpoint_extension": ".npz",
        "fresh_output_root_required_at_campaign_start": True,
        "runner_refuses_overwrite": True,
        "resume_requires_matching_plan_and_manifest_hashes": True,
        "PROTO10_outputs_forbidden_as_PROTO11_outcomes": True,
        "canonical_reduced_result_is_tracked_only_after_campaign_classification": True,
    }:
        raise ValueError("CAL8 provenance differs")

    normalized = deepcopy(raw)
    normalized["artifact_id"] = "FGC-1-CAL7-RUN1-PLAN"
    normalized["protocol_config"] = "configs/fgc/fgc-2-sf1-protocol-v10.toml"
    normalized["protocol_freeze_result"] = "results/fgc-1-pro10-frz1.json"
    normalized["runtime_authorization_result"] = "results/fgc-1-hlt8-mon8.json"
    normalized["scope"]["target_protocol"] = "FGC-2-SF1-PROTO10"
    normalized["scope"]["role"] = (
        "fresh_outcome_neutral_dynamic_amplitude_recalibration_after_guarded_finest_pair_spectral_interpretation"
    )
    for key in REFERENCE_NUMERICS:
        normalized["numerics"].pop(key)
    normalized["provenance"] = {
        "output_root": "runs/fgc-2-sf1/proto10/calibration",
        "campaign_manifest": "runs/fgc-2-sf1/proto10/calibration/manifest.json",
        "append_only_event_log_name": "events.jsonl",
        "checkpoint_extension": ".npz",
        "fresh_output_root_required_at_campaign_start": True,
        "runner_refuses_overwrite": True,
        "resume_requires_matching_plan_and_manifest_hashes": True,
        "PROTO9_outputs_forbidden_as_PROTO10_outcomes": True,
        "canonical_reduced_result_is_tracked_only_after_campaign_classification": True,
    }
    baseline_blob = hlt4._git_blob(
        CHECKPOINT_COMMIT,
        EXPECTED_PATHS["predecessor_run_plan"],
        "immutable CAL7 plan",
    )
    if sha256(baseline_blob).hexdigest() != EXPECTED_LINEAGE[
        "predecessor_run_plan_sha256"
    ]:
        raise ValueError("immutable CAL7 plan hash differs")
    if normalized != tomllib.loads(baseline_blob.decode("utf-8")):
        raise ValueError("CAL8 changes more than the declared PROTO11 delta")
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
        raise ValueError("HLT9 root keys differ")
    if (
        raw["schema_version"] != 1
        or raw["artifact_id"] != ARTIFACT_ID
        or raw["project_version"] != PROJECT_VERSION
        or raw["metric_signature"] != "-+++"
        or raw["riemann_convention"] != "plus_partial_mu_gamma_nu"
    ):
        raise ValueError("HLT9 identity or convention differs")
    paths = {
        name: _path(name, raw[name], expected)
        for name, expected in EXPECTED_PATHS.items()
    }
    if raw["scope"] != EXPECTED_SCOPE:
        raise ValueError("HLT9 scope differs")
    if raw["immutable_lineage"] != EXPECTED_LINEAGE:
        raise ValueError("HLT9 immutable lineage differs")
    if set(raw["runtime_binding"]) != EXPECTED_RUNTIME_KEYS or any(
        value is not True for value in raw["runtime_binding"].values()
    ):
        raise ValueError("HLT9 runtime binding differs")
    if raw["input_freeze"] != EXPECTED_INPUT:
        raise ValueError("HLT9 input freeze differs")
    if set(raw["controls"]) != EXPECTED_CONTROL_KEYS or any(
        value is not True for value in raw["controls"].values()
    ):
        raise ValueError("HLT9 controls differ")
    if raw["namespace"] != EXPECTED_NAMESPACE:
        raise ValueError("HLT9 namespace differs")
    if set(raw["proof_contract"]) != EXPECTED_PROOF_KEYS or any(
        value is not True for value in raw["proof_contract"].values()
    ):
        raise ValueError("HLT9 proof contract differs")
    if raw["claims"] != EXPECTED_CLAIMS:
        raise ValueError("HLT9 claims differ")
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
        raise ValueError("HLT9 checkpoint is not an ancestor of HEAD")
    immutable_names = (
        "protocol_config",
        "protocol_freeze_config",
        "protocol_freeze_result",
        "predecessor_runtime_result",
        "diagnosis_result",
        "static_input_result",
        "predecessor_run_plan",
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

    protocol = validate_sf1_protocol_v11(
        tomllib.loads(blobs["protocol_config"].decode("utf-8"))
    )
    freeze = hlt4._load_json_bytes(
        blobs["protocol_freeze_result"], "immutable PRO11-FRZ1"
    )
    predecessor = hlt4._load_json_bytes(
        blobs["predecessor_runtime_result"], "immutable HLT8"
    )
    diagnosis = hlt4._load_json_bytes(
        blobs["diagnosis_result"], "immutable CAL7/PREF9"
    )
    static = hlt4._load_json_bytes(blobs["static_input_result"], "immutable ID2")
    predecessor_inputs = predecessor.get("artifact_payload", {}).get(
        "frozen_run_inputs", []
    )
    predecessor_events = predecessor.get("artifact_payload", {}).get(
        "initial_common_events", []
    )
    diagnosis_payload = diagnosis.get("artifact_payload", {})
    diagnosis_gates = diagnosis.get("gate_status", {})
    if (
        protocol.get("artifact_id") != "FGC-2-SF1-PROTO11"
        or protocol.get("protocol_version") != 11
        or protocol.get("point_counts") != list(PROTO9_POINT_COUNTS)
        or protocol.get("interior_q_reprojection_permitted") is not False
        or protocol.get("new_numeric_tolerance_added") is not False
        or freeze.get("artifact_id") != "FGC-1-PRO11-FRZ1"
        or any(value is not False for value in freeze.get("gate_status", {}).values())
        or predecessor.get("artifact_id") != "FGC-1-HLT8-MON8"
        or predecessor.get("gate_status", {}).get(
            "PROTO10_fresh_GR0_dynamic_calibration_authorized"
        )
        is not True
        or predecessor.get("gate_status", {}).get(
            "PROTO10_resolved_holdout_manifest_authorized"
        )
        is not False
        or diagnosis.get("artifact_id") != "FGC-1-CAL7-PREF9"
        or diagnosis_gates.get("PROTO10_campaign_terminated_normally") is not True
        or diagnosis_gates.get("PROTO10_GR0_case_eligible") is not False
        or diagnosis_gates.get(
            "PROTO10_center_adjacent_binary64_source_floor_localized"
        )
        is not True
        or "SGB-L or FGC-QR dynamics"
        not in diagnosis_payload.get("epistemic_boundary", {}).get(
            "not_observed", []
        )
        or static.get("gate_status", {}).get("all_case_static_ledger_completed")
        is not True
        or len(predecessor_inputs) != 12
        or len(predecessor_events) != 4
        or any(
            item.get("guarded_PROTO10_admission_passed") is not True
            for item in predecessor_events
        )
        or any(value is not False for value in freeze.get("nonclaims", {}).values())
        or any(value is not False for value in predecessor.get("nonclaims", {}).values())
        or any(value is not False for value in diagnosis.get("nonclaims", {}).values())
    ):
        raise ValueError("immutable HLT9 predecessor gate differs")

    hash_sets = (
        freeze.get("implementation_sha256", {}),
        predecessor.get("implementation_sha256", {}),
        diagnosis.get("implementation_sha256", {}),
        static.get("implementation_sha256", {}),
    )
    checked: set[str] = set()
    inherited_hashes: dict[str, str] = {}
    for hashes in hash_sets:
        for relative, expected in hashes.items():
            current = REPOSITORY / relative
            if not current.is_file() or hlt4._sha(current) != expected:
                raise ValueError(f"inherited PROTO11 premise drifted: {relative}")
            checked.add(relative)
            inherited_hashes[relative] = expected
    return {
        "checkpoint_commit": commit,
        "checkpoint_is_ancestor_of_HEAD": True,
        "protocol": {
            "artifact_id": protocol["artifact_id"],
            "protocol_version": protocol["protocol_version"],
            "point_counts": protocol["point_counts"],
            "initial_q_rule": protocol["initial_q_rule"],
            "source_q_r_rule": protocol["source_q_r_rule"],
            "reduction_constraint_rule": protocol["reduction_constraint_rule"],
            "semantic_holdout_contract_sha256": protocol[
                "semantic_holdout_contract_sha256"
            ],
        },
        "freeze": {
            "artifact_id": freeze["artifact_id"],
            "result_sha256": lineage["protocol_freeze_result_sha256"],
            "claims_all_false": True,
        },
        "predecessor_runtime": {
            "artifact_id": predecessor["artifact_id"],
            "result_sha256": lineage["predecessor_runtime_result_sha256"],
            "frozen_input_count": len(predecessor_inputs),
            "t0_event_count": len(predecessor_events),
        },
        "diagnosis": {
            "artifact_id": diagnosis["artifact_id"],
            "result_sha256": lineage["diagnosis_result_sha256"],
            "campaign_terminated_normally": True,
            "eligible_GR0_case": False,
            "localized_reference_map_requirement": True,
        },
        "static_input": {
            "artifact_id": static["artifact_id"],
            "result_sha256": lineage["static_input_result_sha256"],
        },
        "predecessor_run_plan_sha256": lineage["predecessor_run_plan_sha256"],
        "inherited_implementation_files_hash_matched": len(checked),
        "predecessor_frozen_inputs": predecessor_inputs,
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
    return Proto11GR0EvolutionOperator(
        grid,
        spatial_order=order,
        ko_dissipation=float(Q(numerics["ko_dissipation"])),
        raw_tolerance=float(Q(thresholds["source_residual_infinity_max"])),
        kinetic_condition_maximum=float(
            Q(thresholds["kinetic_condition_number_max"])
        ),
    )


def _freeze_inputs(
    loaded: Mapping[str, Any], lineage: Mapping[str, Any]
) -> tuple[
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[Proto11GR0CommonEventAssessment],
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
        for count in PROTO9_POINT_COUNTS
    ]
    if [
        (item["amplitude"], item["method"], item["point_count"])
        for item in predecessor_inputs
    ] != expected_order:
        raise ValueError("immutable HLT8 input order differs")

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
        old_state = project_gr0_semidiscrete_state(initial, spatial_order=order)
        if array_content_sha256(old_state.u, old_state.p, old_state.q) != predecessor[
            "projected_state_sha256"
        ]:
            raise ValueError("reconstructed input differs from immutable HLT8")
        state = project_gr0_reference_balanced_state(initial, spatial_order=order)
        reference = exact_spherical_minkowski_reference(initial.grid)
        expected_q = SBPFirstDerivative(initial.grid, order).differentiate(
            state.u - reference.u,
            center_parities=ADM_CENTER_PARITIES,
        ) + reference.q
        if (
            not np.array_equal(state.u, old_state.u)
            or not np.array_equal(state.p, old_state.p)
            or not np.array_equal(state.q, expected_q)
        ):
            raise ValueError("PROTO11 projected input violates the reference map")
        state_hash_before = array_content_sha256(state.u, state.p, state.q)
        rhs = _operator(plan, initial.grid, order)(0.0, state)
        state_hash_after = array_content_sha256(state.u, state.p, state.q)
        residual = float(rhs.diagnostics["source_residual_infinity"])
        if (
            residual >= threshold
            or rhs.diagnostics["source_raw_gate_passed"] is not True
            or rhs.diagnostics["PROTO11_exact_reference_equilibrium_applied"]
            is not False
            or rhs.diagnostics["PROTO11_interior_q_reprojected"] is not False
            or state_hash_after != state_hash_before
        ):
            raise ValueError("PROTO11 accepted-state source precheck failed")

        record = dict(predecessor)
        record["predecessor_HLT8_projected_state_sha256"] = record.pop(
            "projected_state_sha256"
        )
        record["predecessor_HLT8_expanded_run_config_sha256"] = record.pop(
            "expanded_run_config_sha256"
        )
        record["predecessor_HLT8_plan_sha256"] = record.pop("plan_sha256")
        record["predecessor_HLT8_protocol_sha256"] = record.pop(
            "protocol_sha256"
        )
        record["predecessor_HLT8_resolution_role"] = record.pop(
            "PROTO10_resolution_role"
        )
        record.update(
            {
                "plan_sha256": plan_hash,
                "protocol_sha256": protocol_hash,
                "PROTO11_resolution_role": (
                    "direct_state_with_reference_balanced_constraints_and_"
                    "guarded_PROTO10_spectra"
                ),
                "projected_state_sha256": state_hash_before,
                "u_bitwise_matches_HLT8": True,
                "p_bitwise_matches_HLT8": True,
                "q_matches_frozen_reference_map": True,
                "q_differs_from_HLT8_native_map": bool(
                    not np.array_equal(state.q, old_state.q)
                ),
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
                "source_residual_infinity": residual,
                "source_residual_maximum": threshold,
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
    assessments: list[Proto11GR0CommonEventAssessment] = []
    for amplitude in ("5/2", "3"):
        for method in ("RK4", "SSPRK3"):
            records = sorted(states[(amplitude, method)], key=lambda item: item[0])
            assessed = proto11_gr0_common_event(
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
            if (
                assessed.reference_balanced_constraint_admission_passed is not True
                or assessed.legacy_PROTO10_admission_passed is not True
                or assessed.spatial_spectral_admission.admission_passed is not True
                or assessed.spatial_spectral_admission.saturation_used is not True
                or assessed.admission_passed is not True
                or assessed.reference_state_map_applied is not True
                or assessed.interior_q_reprojected is not False
            ):
                raise ValueError("PROTO11 t0 common event failed closed")
            events.append(
                {
                    "amplitude": amplitude,
                    "method": method,
                    "coordinate_time": 0.0,
                    "point_counts": list(assessed.point_counts),
                    "reference_balanced_constraint_admission_passed": True,
                    "legacy_PROTO10_admission_passed": True,
                    "raw_PROTO9_admission_passed": (
                        assessed.raw_PROTO9_admission_passed
                    ),
                    "guarded_PROTO10_spatial_admission_passed": True,
                    "PROTO11_admission_passed": True,
                    "diagnostic_saturation_used": (
                        assessed.spatial_spectral_admission.saturation_used
                    ),
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
) -> dict[str, Any]:
    plan = loaded["run_plan"]
    exact_records = []
    any_old_direct_difference = False
    for method in ("RK4", "SSPRK3"):
        order = 4 if method == "RK4" else 2
        for count in PROTO9_POINT_COUNTS:
            grid = UniformRadialGrid(0.0, 128.0, count)
            state = exact_spherical_minkowski_reference(grid)
            before = array_content_sha256(state.u, state.p, state.q)
            rhs = _operator(plan, grid, order)(0.0, state)
            after = array_content_sha256(state.u, state.p, state.q)
            old_direct_q = SBPFirstDerivative(grid, order).differentiate(
                state.u,
                center_parities=ADM_CENTER_PARITIES,
            )
            old_difference = float(np.max(np.abs(old_direct_q - state.q), initial=0.0))
            any_old_direct_difference = any_old_direct_difference or old_difference > 0.0
            if (
                before != after
                or not np.array_equal(rhs.du, np.zeros_like(rhs.du))
                or not np.array_equal(rhs.dp, np.zeros_like(rhs.dp))
                or not np.array_equal(rhs.dq, np.zeros_like(rhs.dq))
                or rhs.diagnostics["source_raw_gate_passed"] is not True
                or rhs.diagnostics["PROTO11_exact_reference_equilibrium_applied"]
                is not True
            ):
                raise RuntimeError("PROTO11 exact-reference control failed")
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
                    "old_direct_q_map_difference_infinity": old_difference,
                }
            )
    if not any_old_direct_difference:
        raise RuntimeError("PROTO11 control did not distinguish the old direct map")

    grid = UniformRadialGrid(0.0, 128.0, PROTO9_POINT_COUNTS[0])
    reference = exact_spherical_minkowski_reference(grid)
    perturbed_q = reference.q.copy()
    perturbed_q[1, 3] = np.nextafter(perturbed_q[1, 3], np.inf)
    one_bit = EvolutionState(reference.u, reference.p, perturbed_q)
    one_bit_rhs = _operator(plan, grid, 4)(0.0, one_bit)
    if one_bit_rhs.diagnostics["PROTO11_exact_reference_equilibrium_applied"]:
        raise RuntimeError("PROTO11 one-bit control used the exact-reference branch")

    generic_q = reference.q.copy()
    generic_q[grid.point_count // 4, 5] += 1.0e-6
    generic = EvolutionState(reference.u, reference.p, generic_q)
    derivative = SBPFirstDerivative(grid, 4)
    _, generic_q_r, _ = reference_balanced_spatial_derivatives(generic, derivative)
    generic_reduction = reference_balanced_reduction_constraint(generic, derivative)
    if (
        float(np.max(np.abs(generic_reduction), initial=0.0)) <= 0.0
        or float(np.max(np.abs(generic_q_r), initial=0.0)) <= 0.0
    ):
        raise RuntimeError("PROTO11 generic defect was hidden")

    first_count, first_state, first_grid = states[("5/2", "RK4")][0]
    if first_count != PROTO9_POINT_COUNTS[0]:
        raise RuntimeError("PROTO11 control selected the wrong physical state")
    snapshot = reference_balanced_semidiscrete_constraint_snapshot(
        first_state,
        first_grid,
        coordinate_time=0.0,
        diagnostic_spatial_order=4,
    )
    expected_reduction = reference_balanced_reduction_constraint(
        first_state,
        SBPFirstDerivative(first_grid, 4),
    )
    reduction_difference = float(
        np.max(
            np.abs(snapshot.residuals[:, 4:] - expected_reduction),
            initial=0.0,
        )
    )
    if reduction_difference != 0.0:
        raise RuntimeError("PROTO11 common-event reduction rows use another map")
    return {
        "exact_reference_records": exact_records,
        "exact_reference_control_count": len(exact_records),
        "all_exact_reference_controls_bitwise_preserved": True,
        "complete_raw_source_serialized_before_exact_equilibrium_selection": True,
        "one_bit_perturbation_field": "q_R_at_first_positive_node",
        "one_bit_perturbation_used_exact_equilibrium_branch": False,
        "one_bit_perturbation_raw_source_residual_infinity": one_bit_rhs.diagnostics[
            "source_residual_infinity"
        ],
        "generic_reduction_defect_infinity": float(
            np.max(np.abs(generic_reduction), initial=0.0)
        ),
        "generic_source_q_r_defect_infinity": float(
            np.max(np.abs(generic_q_r), initial=0.0)
        ),
        "generic_reduction_defect_visible": True,
        "generic_source_derivative_defect_visible": True,
        "common_event_reduction_map_difference_infinity": reduction_difference,
        "common_event_reduction_rows_match_reference_map_bitwise": True,
        "old_direct_map_distinguished_on_at_least_one_frozen_operator": True,
        "interior_q_reprojected_by_any_control": False,
        "new_numeric_tolerance_added": False,
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
            raise ValueError("historical HLT9 namespace evidence differs")
        return dict(historical)
    records = []
    for key in ("calibration_output_root", "holdout_output_root"):
        relative = raw["namespace"][key]
        path = REPOSITORY / relative
        exists = path.exists()
        if exists and (not path.is_dir() or any(path.iterdir())):
            raise ValueError(f"HLT9 fresh namespace is nonempty: {relative}")
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
    inputs, prechecks, events, _assessments, states = _freeze_inputs(
        loaded, lineage
    )
    controls = _runtime_controls(loaded, states)
    if (
        len(inputs) != 12
        or len(prechecks) != 12
        or len(events) != 4
        or any(not item["accepted_state_source_gate_passed"] for item in prechecks)
        or any(not item["PROTO11_admission_passed"] for item in events)
        or controls["all_controls_passed"] is not True
    ):
        raise ValueError("HLT9 authorization evidence is incomplete")
    namespace = _namespace_precondition(raw, historical_namespace_evidence)
    plan_hash = hlt4._sha(loaded["paths"]["run_plan_config"])
    input_manifest_hash = hlt4._digest(inputs)
    campaign_id = hlt4._digest(
        {
            "artifact_id": loaded["run_plan"]["artifact_id"],
            "run_plan_sha256": plan_hash,
            "input_manifest_sha256": input_manifest_hash,
            "ordered_amplitudes": ["5/2", "3"],
            "point_counts": list(PROTO9_POINT_COUNTS),
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
            "pre_trajectory_PROTO11_reference_balanced_runtime_and_fresh_GR0_"
            "calibration_authorization"
        ),
        "generated_by": "scripts/reproduce_fgc_hlt9_mon9.py",
        "derivation_document": hlt4._rel(OWNER_DOCUMENT),
        "derivation_document_sha256": hlt4._sha(OWNER_DOCUMENT),
        "source_config_sha256": {
            hlt4._rel(config_path): hlt4._sha(config_path),
            hlt4._rel(loaded["paths"]["run_plan_config"]): plan_hash,
        },
        "predecessor_sha256": {
            EXPECTED_PATHS[name]: raw["immutable_lineage"][f"{name}_sha256"]
            for name in (
                "protocol_config",
                "protocol_freeze_config",
                "protocol_freeze_result",
                "predecessor_runtime_result",
                "diagnosis_result",
                "static_input_result",
                "predecessor_run_plan",
            )
        },
        "implementation_sha256": implementation,
        "scope_bindings": dict(raw["scope"]),
        "gate_status": dict(raw["claims"]),
        "artifact_payload": {
            "immutable_lineage": lineage,
            "reference_map_controls": controls,
            "inherited_runtime": {
                "PROTO10_spectral_compositor_and_direct_coarse_veto_unchanged": True,
                "PROTO7_transaction_source_retry_boundary_health_and_affine_stops_unchanged": True,
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
                "point_counts": list(PROTO9_POINT_COUNTS),
                "first_eligible_candidate_stops_later_candidates": True,
                "expanded_input_count": len(inputs),
                "u_p_match_count": sum(
                    item["u_bitwise_matches_HLT8"]
                    and item["p_bitwise_matches_HLT8"]
                    for item in inputs
                ),
                "reference_mapped_q_count": sum(
                    item["q_matches_frozen_reference_map"] for item in inputs
                ),
                "source_precheck_passed_count": sum(
                    item["accepted_state_source_gate_passed"] for item in prechecks
                ),
                "PROTO11_t0_passed_count": sum(
                    item["PROTO11_admission_passed"] for item in events
                ),
                "t0_common_event_evaluated_count": len(events),
            },
            "frozen_run_inputs": inputs,
            "accepted_state_source_prechecks": prechecks,
            "initial_common_events": events,
            "namespace_precondition": namespace,
            "decision": (
                "HLT9_binds_the_PROTO11_reference_map_and_authorizes_one_fresh_"
                "GR0_calibration_but_keeps_every_candidate_and_physical_claim_closed"
            ),
            "epistemic_boundary": {
                "PROTO10_campaign_and_CAL7_diagnosis_history_disclosed": True,
                "PROTO11_trajectory_read": False,
                "PROTO11_trajectory_advanced": False,
                "all_twelve_physical_inputs_reconstructed": True,
                "all_twelve_source_prechecks_passed": True,
                "all_four_PROTO11_t0_admissions_passed": True,
                "complete_PROTO10_spectral_evidence_retained": True,
                "diagnostic_saturation_is_physical_resolution": False,
                "round_trip_is_a_continuum_error_bound": False,
                "reference_map_is_a_physical_mechanism_result": False,
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
            "reference_map_physically_validated": False,
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
    result = hlt4._load_json_bytes(path.read_bytes(), "HLT9 result")
    if result.get("artifact_id") != ARTIFACT_ID:
        raise ValueError("HLT9 result artifact identity differs")
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
        raise ValueError("HLT9 output differs from a fresh reproduction")


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
        "(PROTO11 runtime=true; fresh GR0 calibration authorized=true; "
        "trajectory_read=false)"
    )


if __name__ == "__main__":
    main()
