#!/usr/bin/env python3
"""Regenerate the pre-trajectory FGC-1-HLT7-MON7 authorization record."""

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
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-hlt7-mon7.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-hlt7-mon7.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-hlt7-mon7.md"

from scripts import reproduce_fgc_hlt4_mon4 as hlt4  # noqa: E402
from scripts import reproduce_fgc_hlt6_mon6 as hlt6  # noqa: E402
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
from recursive_horizons.fgc.evolution.proto7_runtime import (  # noqa: E402
    Proto7GR0EvolutionOperator,
)
from recursive_horizons.fgc.evolution.proto9_runtime import (  # noqa: E402
    PROTO9_POINT_COUNTS,
    proto9_gr0_common_event,
)
from recursive_horizons.fgc.evolution.protocol_v9 import (  # noqa: E402
    validate_sf1_protocol_v9,
)
from recursive_horizons.fgc.evolution.vectorized_source import (  # noqa: E402
    ADM_CENTER_PARITIES,
)
from recursive_horizons.fgc.initial_data_preflight import (  # noqa: E402
    PulseParameters,
)


Q = hlt4.Q
ARTIFACT_ID = "FGC-1-HLT7-MON7"
PROJECT_VERSION = "0.11.0"
CHECKPOINT_COMMIT = "833910d1b31ed820745d61e0ce7976965e37b487"
EXPECTED_PATHS = {
    "protocol_config": "configs/fgc/fgc-2-sf1-protocol-v9.toml",
    "protocol_freeze_config": "configs/fgc/fgc-1-pro9-frz1.toml",
    "protocol_freeze_result": "results/fgc-1-pro9-frz1.json",
    "predecessor_runtime_result": "results/fgc-1-hlt6-mon6.json",
    "static_input_result": "results/fgc-1-id2-all1.json",
    "predecessor_run_plan": "configs/fgc/fgc-1-cal5-run1.toml",
    "run_plan_config": "configs/fgc/fgc-1-cal6-run1.toml",
}
EXPECTED_SCOPE = {
    "target_protocol": "FGC-2-SF1-PROTO9",
    "runtime_owner": ARTIFACT_ID,
    "calibration_branch": "GR-0",
    "role": "pre_trajectory_resolution_ladder_runtime_binding_and_fresh_input_freeze",
    "PROTO8_calibration_history_disclosed": True,
    "PROTO9_trajectory_read": False,
    "SGBL_trajectory_read": False,
    "FGCQR_trajectory_read": False,
    "collapse_or_trapped_outcome_classified": False,
    "conditional_8193_projection_used_as_evidence": False,
}
EXPECTED_LINEAGE = {
    "checkpoint_commit": CHECKPOINT_COMMIT,
    "protocol_config_sha256": "81ce8426b006518d58accae50e43fcbfe942d013dfd010a716ff265b57afb48e",
    "protocol_freeze_config_sha256": "92f00e53585dc74e7438c353d4fdbcfc0b3b2e634e134cc15e26d7700dea7f1e",
    "protocol_freeze_result_sha256": "cccda72803bd4b6bbdbc091a90e9746e1cc4b7952ef9f4dd9b684aee962b892d",
    "predecessor_runtime_result_sha256": "cda19a43fc08d6100511be046fb12a178adb29848429aa656a07160cf6787cbc",
    "static_input_result_sha256": "10e2661be521087467b970f1cd877311bd63a0cb6874a1ee745ca00ec74ab981",
    "predecessor_run_plan_sha256": "7967c1a785e9b6afcac3d516a426db4ed3652d58ab87fddb8b8874d824129e75",
    "checkpoint_must_be_ancestor_of_HEAD": True,
    "all_predecessor_blobs_must_be_read_from_checkpoint": True,
}
EXPECTED_RUNTIME = {
    "inherited_PROTO8_common_event_compositor_byte_unchanged": True,
    "inherited_PROTO7_evolution_transaction_and_source_retry_byte_unchanged": True,
    "exact_point_counts_required": [2049, 4097, 8193],
    "new_8193_state_required_for_each_amplitude_and_method": True,
    "conditional_projection_is_not_runtime_evidence": True,
    "common_event_constraints_spectra_and_thresholds_unchanged": True,
    "physical_inputs_equations_methods_CFL_and_stops_unchanged": True,
    "runner_binds_and_restores_immutable_predecessor_engine": True,
}
EXPECTED_INPUT = {
    "eligible_amplitudes_in_order": ["5/2", "3"],
    "methods_in_order": ["RK4", "SSPRK3"],
    "resolutions_in_order": [2049, 4097, 8193],
    "expected_run_input_count": 12,
    "overlapping_2049_and_4097_state_hashes_must_match_HLT6": True,
    "four_new_8193_states_must_be_constructed_directly": True,
    "physical_u_p_and_initial_profiles_unchanged": True,
    "q_must_be_initialized_by_native_SBP_operator": True,
    "expanded_PROTO9_run_config_hashes_must_be_serialized": True,
    "all_twelve_accepted_state_source_prechecks_must_pass": True,
    "both_amplitudes_and_methods_must_be_evaluated_at_the_PROTO9_t0_common_event": True,
    "first_eligible_candidate_stops_later_candidates": True,
}
EXPECTED_NAMESPACE = {
    "calibration_output_root": "runs/fgc-2-sf1/proto9/calibration",
    "holdout_output_root": "runs/fgc-2-sf1/proto9/holdout",
    "both_roots_must_be_absent_or_empty_at_authorization": True,
    "authorization_must_not_create_either_root": True,
    "runner_must_refuse_overwrite": True,
}
EXPECTED_PROOF_KEYS = {
    "PROTO9_and_PRO9_FRZ1_must_validate_from_immutable_checkpoint",
    "HLT6_ID2_and_CAL5_plan_must_be_canonical_and_hash_matched_from_immutable_checkpoint",
    "CAL6_plan_must_differ_from_CAL5_only_by_the_frozen_ladder_protocol_and_namespace",
    "old_PROTO8_ladder_must_fail_the_PROTO9_runtime_adapter",
    "all_overlapping_projected_inputs_must_match_HLT6",
    "all_four_new_8193_inputs_must_be_constructed_without_projection",
    "all_twelve_t0_accepted_source_prechecks_must_pass",
    "all_four_PROTO9_t0_common_events_must_be_serialized_without_reinterpretation",
    "inherited_constraint_and_spectral_thresholds_must_remain_unchanged",
    "inherited_transaction_source_boundary_health_and_affine_stops_must_remain_hash_bound",
    "authorized_runner_must_be_hash_bound_and_refuse_drift_overwrite_or_dirty_tracked_state",
    "fresh_namespaces_must_be_observed_empty_without_creation",
    "canonical_hash_bound_result_required",
    "no_PROTO9_trajectory_may_be_read_or_advanced",
}
EXPECTED_CLAIMS = {
    "PROTO9_successor_runtime_implemented": True,
    "PROTO9_fresh_GR0_dynamic_calibration_authorized": False,
    "PROTO9_resolved_holdout_manifest_authorized": False,
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
IMPLEMENTATION = (
    REPOSITORY / "src/recursive_horizons/fgc/evolution/calibration_runtime.py",
    REPOSITORY / "src/recursive_horizons/fgc/evolution/cal2_source_diagnosis.py",
    REPOSITORY / "src/recursive_horizons/fgc/evolution/proto5_runtime.py",
    REPOSITORY / "src/recursive_horizons/fgc/evolution/proto6_runtime.py",
    REPOSITORY / "src/recursive_horizons/fgc/evolution/proto7_runtime.py",
    REPOSITORY / "src/recursive_horizons/fgc/evolution/cal4_common_event_diagnosis.py",
    REPOSITORY / "src/recursive_horizons/fgc/evolution/proto8_runtime.py",
    REPOSITORY / "src/recursive_horizons/fgc/evolution/proto9_runtime.py",
    REPOSITORY / "src/recursive_horizons/fgc/evolution/protocol_v9.py",
    REPOSITORY / "scripts/run_fgc_gr0_calibration.py",
    REPOSITORY / "scripts/run_fgc_gr0_calibration_v6.py",
    REPOSITORY / "scripts/run_fgc_gr0_calibration_v7.py",
    REPOSITORY / "scripts/run_fgc_gr0_calibration_v8.py",
    REPOSITORY / "scripts/run_fgc_gr0_calibration_v9.py",
    Path(__file__).resolve(),
)


def _path(name: str, value: Any, expected: str) -> Path:
    if value != expected:
        raise ValueError(f"HLT7 {name} differs")
    path = REPOSITORY / expected
    if not path.is_file():
        raise ValueError(f"HLT7 {name} is absent")
    return path


def validate_run_plan(path: Path) -> dict[str, Any]:
    """Prove CAL6 differs from immutable CAL5 only by PROTO9's freeze."""

    with path.open("rb") as handle:
        raw = tomllib.load(handle)
    if (
        raw.get("schema_version") != 1
        or raw.get("artifact_id") != "FGC-1-CAL6-RUN1-PLAN"
        or raw.get("project_version") != PROJECT_VERSION
        or raw.get("protocol_config") != EXPECTED_PATHS["protocol_config"]
        or raw.get("protocol_freeze_result")
        != EXPECTED_PATHS["protocol_freeze_result"]
        or raw.get("runtime_authorization_result")
        != "results/fgc-1-hlt7-mon7.json"
        or raw.get("static_input_result") != EXPECTED_PATHS["static_input_result"]
    ):
        raise ValueError("CAL6 run-plan identity or ownership differs")
    if raw.get("scope") != {
        "target_protocol": "FGC-2-SF1-PROTO9",
        "branch": "GR-0",
        "role": "fresh_outcome_neutral_dynamic_amplitude_recalibration_after_resolution_ladder_shift",
        "physical_equations": "unredefined_GR0_specialization_of_REF1",
        "candidate_action_fields_or_health_stops_forbidden": True,
        "retained_EFT_interpretation": False,
    }:
        raise ValueError("CAL6 scope differs")
    if raw.get("numerics", {}).get("resolutions") != list(PROTO9_POINT_COUNTS):
        raise ValueError("CAL6 resolution ladder differs")
    if raw.get("provenance") != {
        "output_root": "runs/fgc-2-sf1/proto9/calibration",
        "campaign_manifest": "runs/fgc-2-sf1/proto9/calibration/manifest.json",
        "append_only_event_log_name": "events.jsonl",
        "checkpoint_extension": ".npz",
        "fresh_output_root_required_at_campaign_start": True,
        "runner_refuses_overwrite": True,
        "resume_requires_matching_plan_and_manifest_hashes": True,
        "PROTO8_outputs_forbidden_as_PROTO9_outcomes": True,
        "canonical_reduced_result_is_tracked_only_after_campaign_classification": True,
    }:
        raise ValueError("CAL6 provenance differs")

    normalized = deepcopy(raw)
    normalized["artifact_id"] = "FGC-1-CAL5-RUN1-PLAN"
    normalized["protocol_config"] = "configs/fgc/fgc-2-sf1-protocol-v8.toml"
    normalized["protocol_freeze_result"] = "results/fgc-1-pro8-frz1.json"
    normalized["runtime_authorization_result"] = "results/fgc-1-hlt6-mon6.json"
    normalized["scope"]["target_protocol"] = "FGC-2-SF1-PROTO8"
    normalized["scope"]["role"] = (
        "fresh_outcome_neutral_dynamic_amplitude_recalibration_after_common_event_ownership_repair"
    )
    normalized["numerics"]["resolutions"] = [1025, 2049, 4097]
    normalized["provenance"] = {
        "output_root": "runs/fgc-2-sf1/proto8/calibration",
        "campaign_manifest": "runs/fgc-2-sf1/proto8/calibration/manifest.json",
        "append_only_event_log_name": "events.jsonl",
        "checkpoint_extension": ".npz",
        "fresh_output_root_required_at_campaign_start": True,
        "runner_refuses_overwrite": True,
        "resume_requires_matching_plan_and_manifest_hashes": True,
        "PROTO7_outputs_forbidden_as_PROTO8_outcomes": True,
        "canonical_reduced_result_is_tracked_only_after_campaign_classification": True,
    }
    baseline_blob = hlt4._git_blob(
        CHECKPOINT_COMMIT,
        EXPECTED_PATHS["predecessor_run_plan"],
        "immutable CAL5 plan",
    )
    if normalized != tomllib.loads(baseline_blob.decode("utf-8")):
        raise ValueError("CAL6 changes more than the declared PROTO9 delta")
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
        "namespace",
        "proof_contract",
        "claims",
    }
    if set(raw) != expected_root:
        raise ValueError("HLT7 root keys differ")
    if (
        raw["schema_version"] != 1
        or raw["artifact_id"] != ARTIFACT_ID
        or raw["project_version"] != PROJECT_VERSION
        or raw["metric_signature"] != "-+++"
        or raw["riemann_convention"] != "plus_partial_mu_gamma_nu"
    ):
        raise ValueError("HLT7 identity or convention differs")
    paths = {
        name: _path(name, raw[name], expected)
        for name, expected in EXPECTED_PATHS.items()
    }
    if raw["scope"] != EXPECTED_SCOPE:
        raise ValueError("HLT7 scope differs")
    if raw["immutable_lineage"] != EXPECTED_LINEAGE:
        raise ValueError("HLT7 immutable lineage differs")
    if raw["runtime_binding"] != EXPECTED_RUNTIME:
        raise ValueError("HLT7 runtime binding differs")
    if raw["input_freeze"] != EXPECTED_INPUT:
        raise ValueError("HLT7 input freeze differs")
    if raw["namespace"] != EXPECTED_NAMESPACE:
        raise ValueError("HLT7 namespace differs")
    if set(raw["proof_contract"]) != EXPECTED_PROOF_KEYS or any(
        value is not True for value in raw["proof_contract"].values()
    ):
        raise ValueError("HLT7 proof contract differs")
    if raw["claims"] != EXPECTED_CLAIMS:
        raise ValueError("HLT7 claims differ")
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
        raise ValueError("HLT7 checkpoint is not an ancestor of HEAD")
    immutable_names = (
        "protocol_config",
        "protocol_freeze_config",
        "protocol_freeze_result",
        "predecessor_runtime_result",
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
        current = REPOSITORY / EXPECTED_PATHS[name]
        if sha256(current.read_bytes()).hexdigest() != expected:
            raise ValueError(f"current {name} differs from immutable checkpoint")

    protocol = validate_sf1_protocol_v9(
        tomllib.loads(blobs["protocol_config"].decode("utf-8"))
    )
    freeze = hlt4._load_json_bytes(
        blobs["protocol_freeze_result"], "immutable PRO9-FRZ1"
    )
    predecessor = hlt4._load_json_bytes(
        blobs["predecessor_runtime_result"], "immutable HLT6"
    )
    static = hlt4._load_json_bytes(blobs["static_input_result"], "immutable ID2")
    inputs = predecessor.get("artifact_payload", {}).get("frozen_run_inputs", [])
    if (
        protocol.get("artifact_id") != "FGC-2-SF1-PROTO9"
        or protocol.get("outcome_neutral_contract_validated") is not True
        or protocol.get("point_counts") != list(PROTO9_POINT_COUNTS)
        or protocol.get("new_8193_data_required") is not True
        or freeze.get("artifact_id") != "FGC-1-PRO9-FRZ1"
        or freeze.get("artifact_payload", {})
        .get("protocol_validation", {})
        .get("frozen")
        is not True
        or freeze.get("gate_status", {}).get("PROTO9_successor_runtime_implemented")
        is not False
        or predecessor.get("artifact_id") != "FGC-1-HLT6-MON6"
        or predecessor.get("gate_status", {}).get(
            "PROTO8_successor_runtime_compositor_implemented"
        )
        is not True
        or predecessor.get("gate_status", {}).get(
            "PROTO8_fresh_GR0_dynamic_calibration_authorized"
        )
        is not True
        or predecessor.get("gate_status", {}).get("FGCQR_holdout_execution_authorized")
        is not False
        or static.get("gate_status", {}).get("all_case_static_ledger_completed")
        is not True
        or len(inputs) != 12
        or any(value is not False for value in freeze.get("nonclaims", {}).values())
        or any(value is not False for value in predecessor.get("nonclaims", {}).values())
    ):
        raise ValueError("immutable HLT7 predecessor gate differs")
    for relative, expected in predecessor.get("implementation_sha256", {}).items():
        current = REPOSITORY / relative
        if not current.is_file() or hlt4._sha(current) != expected:
            raise ValueError(f"inherited PROTO8 runtime drifted: {relative}")
    return {
        "checkpoint_commit": commit,
        "checkpoint_is_ancestor_of_HEAD": True,
        "protocol": {
            "artifact_id": protocol["artifact_id"],
            "protocol_version": protocol["protocol_version"],
            "point_counts": protocol["point_counts"],
            "semantic_holdout_contract_sha256": protocol[
                "semantic_holdout_contract_sha256"
            ],
        },
        "freeze": {
            "artifact_id": freeze["artifact_id"],
            "result_sha256": lineage["protocol_freeze_result_sha256"],
            "claims_all_false": all(
                value is False for value in freeze["gate_status"].values()
            ),
        },
        "predecessor_runtime": {
            "artifact_id": predecessor["artifact_id"],
            "result_sha256": lineage["predecessor_runtime_result_sha256"],
            "implementation_files_hash_matched": len(
                predecessor["implementation_sha256"]
            ),
        },
        "static_input": {
            "artifact_id": static["artifact_id"],
            "result_sha256": lineage["static_input_result_sha256"],
        },
        "predecessor_run_plan_sha256": lineage["predecessor_run_plan_sha256"],
        "predecessor_frozen_inputs": inputs,
    }


def _flat_state(grid: UniformRadialGrid, order: int) -> EvolutionState:
    u = np.zeros((grid.point_count, 6), dtype=np.float64)
    u[:, 0] = 1.0
    u[:, 2] = 1.0
    u[:, 3] = grid.coordinates
    p = np.zeros_like(u)
    q = SBPFirstDerivative(grid, order).differentiate(
        u, center_parities=ADM_CENTER_PARITIES
    )
    return EvolutionState(u, p, q)


def _adapter_controls() -> dict[str, Any]:
    inherited = hlt6._adapter_controls()
    old_grids = tuple(
        UniformRadialGrid(0.0, 128.0, count) for count in (1025, 2049, 4097)
    )
    old_states = tuple(_flat_state(grid, 4) for grid in old_grids)
    old_rejected = False
    try:
        proto9_gr0_common_event(
            old_states,
            old_grids,
            accepted_stage_counts=(0, 0, 0),
            method="RK4",
            coordinate_time=0.0,
        )
    except ValueError:
        old_rejected = True
    reordered_grids = tuple(reversed(old_grids))
    reordered_states = tuple(reversed(old_states))
    reordered_rejected = False
    try:
        proto9_gr0_common_event(
            reordered_states,
            reordered_grids,
            accepted_stage_counts=(0, 0, 0),
            method="RK4",
            coordinate_time=0.0,
        )
    except ValueError:
        reordered_rejected = True
    if not old_rejected or not reordered_rejected:
        raise RuntimeError("PROTO9 resolution adapter controls failed")
    return {
        "inherited_PROTO8_adapter_controls": inherited,
        "old_PROTO8_ladder_rejected": True,
        "reordered_ladder_rejected": True,
        "conditional_projection_path_exists": False,
        "exact_new_ladder_tested_on_physical_t0_inputs": True,
        "all_controls_passed": True,
    }


def _freeze_inputs(
    loaded: Mapping[str, Any],
    lineage: Mapping[str, Any],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
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
    expected_old_order = [
        (amplitude, method, count)
        for amplitude in ("5/2", "3")
        for method in ("RK4", "SSPRK3")
        for count in (1025, 2049, 4097)
    ]
    if [
        (item["amplitude"], item["method"], item["point_count"])
        for item in predecessor_inputs
    ] != expected_old_order:
        raise ValueError("immutable HLT6 input order differs")

    frozen: list[dict[str, Any]] = []
    source_prechecks: list[dict[str, Any]] = []
    states: dict[tuple[str, str], list[tuple[int, EvolutionState, Any]]] = {}
    for amplitude in ("5/2", "3"):
        for method in ("RK4", "SSPRK3"):
            order = 4 if method == "RK4" else 2
            seed = predecessor_by_key[(amplitude, method, 4097)]
            old_hashes = {
                predecessor_by_key[(amplitude, method, count)][
                    "projected_state_sha256"
                ]
                for count in (1025, 2049, 4097)
            }
            for count in PROTO9_POINT_COUNTS:
                initial = construct_gr0_grid_initial_data(
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
                state = project_gr0_semidiscrete_state(initial, spatial_order=order)
                state_hash = array_content_sha256(state.u, state.p, state.q)
                u_preserved = np.array_equal(state.u, initial.state.u)
                p_preserved = np.array_equal(state.p, initial.state.p)
                native_q = SBPFirstDerivative(initial.grid, order).differentiate(
                    state.u,
                    center_parities=ADM_CENTER_PARITIES,
                )
                q_native = np.array_equal(state.q, native_q)
                overlap = count in (2049, 4097)
                if overlap and state_hash != predecessor_by_key[
                    (amplitude, method, count)
                ]["projected_state_sha256"]:
                    raise ValueError("PROTO9 overlapping projected state differs from HLT6")
                if count == 8193 and state_hash in old_hashes:
                    raise ValueError("PROTO9 finest state was not newly constructed")
                if not (u_preserved and p_preserved and q_native):
                    raise ValueError("PROTO9 semidiscrete input projection differs")

                operator = Proto7GR0EvolutionOperator(
                    initial.grid,
                    spatial_order=order,
                    ko_dissipation=float(Q(numerics["ko_dissipation"])),
                    raw_tolerance=float(
                        Q(thresholds["source_residual_infinity_max"])
                    ),
                    kinetic_condition_maximum=float(
                        Q(thresholds["kinetic_condition_number_max"])
                    ),
                )
                rhs = operator(0.0, state)
                residual = float(rhs.diagnostics["source_residual_infinity"])
                threshold = float(Q(thresholds["source_residual_infinity_max"]))
                if residual >= threshold:
                    raise ValueError("PROTO9 t0 accepted-state source precheck failed")
                source_prechecks.append(
                    {
                        "amplitude": amplitude,
                        "method": method,
                        "point_count": count,
                        "source_residual_infinity": residual,
                        "source_residual_maximum": threshold,
                        "accepted_state_source_gate_passed": True,
                        "trajectory_advanced": False,
                    }
                )

                record = dict(seed)
                record.update(
                    {
                        "point_count": count,
                        "grid_spacing": initial.grid.spacing,
                        "projected_state_sha256": state_hash,
                        "plan_sha256": plan_hash,
                        "protocol_sha256": protocol_hash,
                        "u_bitwise_preserved": True,
                        "p_bitwise_preserved": True,
                        "q_initialized_by_native_SBP_operator": True,
                        "PROTO9_resolution_role": (
                            "new_observed_finest_grid"
                            if count == 8193
                            else "reconstructed_overlap_grid"
                        ),
                        "overlap_projected_state_matches_HLT6": overlap,
                        "new_8193_state_constructed": count == 8193,
                        "conditional_projection_used": False,
                        "trajectory_advanced": False,
                    }
                )
                record.pop("expanded_run_config_sha256", None)
                record["expanded_run_config_sha256"] = hlt4._digest(record)
                frozen.append(record)
                states.setdefault((amplitude, method), []).append(
                    (count, state, initial.grid)
                )

    common_events: list[dict[str, Any]] = []
    for amplitude in ("5/2", "3"):
        for method in ("RK4", "SSPRK3"):
            records = sorted(states[(amplitude, method)], key=lambda item: item[0])
            assessed = proto9_gr0_common_event(
                [item[1] for item in records],
                [item[2] for item in records],
                accepted_stage_counts=[0, 0, 0],
                method=method,
                coordinate_time=0.0,
                cutoff=float(Q(physical["cutoff_Lambda"])),
                measurement_radius_maximum=float(
                    Q(physical["measurement_radius_maximum"])
                ),
                taper_fraction=float(Q(numerics["proper_spectral_taper_fraction"])),
                fixed_outer_rows=numerics["fixed_outer_rows"],
            )
            common_events.append(
                {
                    "amplitude": amplitude,
                    "method": method,
                    "coordinate_time": 0.0,
                    "point_counts": list(assessed.point_counts),
                    "accepted_stage_counts": list(assessed.accepted_stage_counts),
                    "excluded_outer_rows": list(
                        assessed.constraint_admission.excluded_outer_rows
                    ),
                    "owned_raw_global_norms": list(
                        assessed.constraint_admission.owned_raw_global_norms
                    ),
                    "constraint_admission_passed": (
                        assessed.constraint_admission.admission_passed
                    ),
                    "constraint_coarsest_guard_passed": (
                        assessed.constraint_admission.coarsest_guard_passed
                    ),
                    "constraint_finest_guard_passed": (
                        assessed.constraint_admission.finest_guard_passed
                    ),
                    "constraint_monotone_refinement_passed": (
                        assessed.constraint_admission.monotone_refinement_passed
                    ),
                    "constraint_finest_pair_order_passed": (
                        assessed.constraint_admission.finest_pair_order_passed
                    ),
                    "constraint_minimum_finite_finest_pair_order": (
                        assessed.constraint_admission.minimum_finite_finest_pair_order
                    ),
                    "coarsest_spatial_individual_budget_passed": (
                        assessed.spatial_spectral_admission.coarsest_individual_budget_passed
                    ),
                    "finest_pair_spatial_individual_budgets_passed": (
                        assessed.spatial_spectral_admission.finest_pair_individual_budgets_passed
                    ),
                    "all_adjacent_spatial_tail_ratios_passed": (
                        assessed.spatial_spectral_admission.every_nested_tail_ratio_passed
                    ),
                    "field_power_tail_ratios": {
                        name: list(values)
                        for name, values in assessed.spatial_spectral_admission.field_power_tail_ratios.items()
                    },
                    "derivative_power_tail_ratios": {
                        name: list(values)
                        for name, values in assessed.spatial_spectral_admission.derivative_power_tail_ratios.items()
                    },
                    "maximum_spatial_round_trip_interpolation_infinity": (
                        assessed.maximum_spatial_round_trip_interpolation_infinity
                    ),
                    "new_8193_state_present": True,
                    "conditional_projection_used": False,
                    "admission_passed": assessed.admission_passed,
                    "trajectory_advanced": False,
                }
            )
    return frozen, source_prechecks, common_events


def _namespace_precondition(
    raw: Mapping[str, Any],
    historical: Mapping[str, Any] | None,
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
            raise ValueError("historical HLT7 namespace evidence differs")
        return dict(historical)
    records = []
    for key in ("calibration_output_root", "holdout_output_root"):
        relative = raw["namespace"][key]
        path = REPOSITORY / relative
        exists = path.exists()
        if exists and (not path.is_dir() or any(path.iterdir())):
            raise ValueError(f"HLT7 fresh namespace is nonempty: {relative}")
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
    controls = _adapter_controls()
    inputs, source_prechecks, common_events = _freeze_inputs(loaded, lineage)
    if (
        len(common_events) != 4
        or any(not item["constraint_admission_passed"] for item in common_events)
        or any(
            not item["finest_pair_spatial_individual_budgets_passed"]
            for item in common_events
        )
        or any(
            item["all_adjacent_spatial_tail_ratios_passed"]
            for item in common_events
        )
        or any(item["admission_passed"] for item in common_events)
    ):
        raise ValueError("HLT7 t0 obstruction no longer has its diagnosed shape")
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
    public_lineage = dict(lineage)
    public_lineage.pop("predecessor_frozen_inputs")
    return {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "classification": "pre_trajectory_PROTO9_resolution_ladder_runtime_with_t0_nested_tail_obstruction",
        "generated_by": "scripts/reproduce_fgc_hlt7_mon7.py",
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
                "static_input_result",
                "predecessor_run_plan",
            )
        },
        "implementation_sha256": {
            hlt4._rel(path): hlt4._sha(path) for path in IMPLEMENTATION
        },
        "scope_bindings": dict(raw["scope"]),
        "gate_status": dict(raw["claims"]),
        "artifact_payload": {
            "immutable_lineage": public_lineage,
            "resolution_adapter_controls": controls,
            "inherited_runtime": {
                "PROTO8_common_event_and_PROTO7_evolution_unchanged": True,
                "inherited_implementation_files_hash_matched": public_lineage[
                    "predecessor_runtime"
                ]["implementation_files_hash_matched"],
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
                "new_8193_input_count": sum(
                    item["point_count"] == 8193 for item in inputs
                ),
                "t0_common_event_passed_count": sum(
                    item["admission_passed"] for item in common_events
                ),
                "t0_common_event_evaluated_count": len(common_events),
            },
            "frozen_run_inputs": inputs,
            "accepted_state_source_prechecks": source_prechecks,
            "initial_common_events": common_events,
            "namespace_precondition": namespace,
            "decision": "HLT7_implements_the_PROTO9_resolution_binding_but_keeps_calibration_closed_by_the_t0_nested_tail_veto",
            "epistemic_boundary": {
                "PROTO8_calibration_history_disclosed": True,
                "PROTO9_trajectory_read": False,
                "PROTO9_trajectory_advanced": False,
                "old_ladder_used_as_PROTO9_evidence": False,
                "conditional_8193_projection_used_as_evidence": False,
                "four_new_8193_t0_states_constructed": True,
                "t0_physical_inputs_reconstructed": True,
                "trapped_or_collapse_outcome_classified": False,
                "SGBL_or_FGCQR_outcome_read": False,
                "fresh_GR0_recalibration_authorized": False,
                "t0_constraint_admission_passed_for_all_four_cases": True,
                "t0_absolute_spatial_budgets_passed_for_all_four_cases": True,
                "t0_nested_tail_admission_passed_for_any_case": False,
                "resolved_holdout_manifest_exists": False,
                "mechanism_question_answered": False,
            },
        },
        "nonclaims": {
            "fresh_GR0_recalibration_completed": False,
            "trapped_sphere_formed_dynamically": False,
            "SGBL_health_definition_supplied": False,
            "SGBL_comparison_completed": False,
            "FGCQR_trajectory_opened": False,
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
    return hlt4._load_json_bytes(path.read_bytes(), "HLT7 result")


def verify_canonical(path: Path, config_path: Path = DEFAULT_CONFIG) -> None:
    observed = load_canonical_result(path)
    loaded = load_config(config_path)
    calibration_root = REPOSITORY / loaded["raw"]["namespace"][
        "calibration_output_root"
    ]
    if calibration_root.exists() and any(calibration_root.iterdir()):
        raise ValueError(
            "HLT7 left PROTO9 calibration unauthorized but its output root is nonempty"
        )
    expected = record(config_path)
    if observed != hlt4._serial(expected):
        raise ValueError("HLT7 output differs from a fresh reproduction")


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
        "(PROTO9 runtime=true; fresh calibration authorized=false; trajectory_read=false)"
    )


if __name__ == "__main__":
    main()
