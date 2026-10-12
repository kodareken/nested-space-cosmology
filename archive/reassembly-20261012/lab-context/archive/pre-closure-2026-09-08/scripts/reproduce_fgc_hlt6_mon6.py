#!/usr/bin/env python3
"""Regenerate the pre-trajectory FGC-1-HLT6-MON6 authorization record."""

from __future__ import annotations

import argparse
from copy import deepcopy
from dataclasses import replace
from hashlib import sha256
from pathlib import Path
import subprocess
import sys
import tomllib
from typing import Any, Mapping

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-hlt6-mon6.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-hlt6-mon6.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-hlt6-mon6.md"

from scripts import reproduce_fgc_hlt4_mon4 as hlt4  # noqa: E402
from scripts import reproduce_fgc_hlt5_mon5 as hlt5  # noqa: E402
from recursive_horizons.fgc.evolution.cal4_common_event_diagnosis import (  # noqa: E402
    finest_pair_nested_spectral_admission,
    owned_semidiscrete_constraint_snapshot,
)
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
from recursive_horizons.fgc.evolution.proto4_admission import (  # noqa: E402
    PROTO4_CONSTRAINT_ORDER,
    SpectralPowerBudget,
)
from recursive_horizons.fgc.evolution.proto7_runtime import (  # noqa: E402
    Proto7GR0EvolutionOperator,
)
from recursive_horizons.fgc.evolution.proto8_runtime import (  # noqa: E402
    proto8_gr0_common_event,
)
from recursive_horizons.fgc.evolution.protocol_v8 import (  # noqa: E402
    validate_sf1_protocol_v8,
)
from recursive_horizons.fgc.evolution.vectorized_source import (  # noqa: E402
    ADM_CENTER_PARITIES,
)
from recursive_horizons.fgc.initial_data_preflight import (  # noqa: E402
    PulseParameters,
)


Q = hlt4.Q
ARTIFACT_ID = "FGC-1-HLT6-MON6"
PROJECT_VERSION = "0.11.0"
CHECKPOINT_COMMIT = "b4c3cdcd216e305ef55e2f2700fbee4111afc95c"
EXPECTED_PATHS = {
    "protocol_config": "configs/fgc/fgc-2-sf1-protocol-v8.toml",
    "protocol_freeze_config": "configs/fgc/fgc-1-pro8-frz1.toml",
    "protocol_freeze_result": "results/fgc-1-pro8-frz1.json",
    "predecessor_runtime_result": "results/fgc-1-hlt5-mon5.json",
    "static_input_result": "results/fgc-1-id2-all1.json",
    "predecessor_run_plan": "configs/fgc/fgc-1-cal4-run1.toml",
    "run_plan_config": "configs/fgc/fgc-1-cal5-run1.toml",
}
EXPECTED_SCOPE = {
    "target_protocol": "FGC-2-SF1-PROTO8",
    "runtime_owner": ARTIFACT_ID,
    "calibration_branch": "GR-0",
    "role": "pre_trajectory_common_event_diagnostic_ownership_runtime_compositor_and_fresh_input_freeze",
    "PROTO7_calibration_history_disclosed": True,
    "PROTO8_trajectory_read": False,
    "SGBL_trajectory_read": False,
    "FGCQR_trajectory_read": False,
    "collapse_or_trapped_outcome_classified": False,
    "candidate_action_health_values_invented_for_GR0": False,
}
EXPECTED_LINEAGE = {
    "checkpoint_commit": CHECKPOINT_COMMIT,
    "protocol_config_sha256": "5de88456ed041f8a4293150f49b4caea057560dd796c8ea2b9a8110ecf7770ab",
    "protocol_freeze_config_sha256": "ec62cc24dc6127520bccb296488798602fa83bbe1af3b8044b4193584c122ae9",
    "protocol_freeze_result_sha256": "e6fc8a94518a3458fb139773f20f0f3009940006c9571c30e75ba4a8f45c5c64",
    "predecessor_runtime_result_sha256": "d32104f266feb81b1f9ca00a7d3e080dccb7a362d1775e252e307fe0e8c0daf9",
    "static_input_result_sha256": "10e2661be521087467b970f1cd877311bd63a0cb6874a1ee745ca00ec74ab981",
    "predecessor_run_plan_sha256": "c5f289a9be96134364f028cd808c24d418c428614ad5ce6a3878258d5fcab8ac",
    "checkpoint_must_be_ancestor_of_HEAD": True,
    "all_predecessor_blobs_must_be_read_from_checkpoint": True,
}
EXPECTED_COMMON_EVENT = {
    "applies_only_after_all_six_members_commit_to_one_bitwise_common_time": True,
    "projector_fixed_outer_rows": 4,
    "RK4_derivative_stencil_reach_rows": 3,
    "SSPRK3_derivative_stencil_reach_rows": 1,
    "full_domain_and_owned_domain_raw_constraint_residuals_are_public": True,
    "owned_raw_constraint_norms_decide_unchanged_magnitude_guards": True,
    "minimum_finest_pair_constraint_order": "3/2",
    "accumulated_roundoff_formula": "4096*binary64_epsilon*(1+accepted_stage_count)",
    "accumulated_roundoff_changes_order_classification_only": True,
    "coarsest_spatial_spectrum_is_public_convergence_witness": True,
    "medium_and_finest_spatial_spectra_each_pass_every_unchanged_absolute_budget": True,
    "all_adjacent_spatial_tail_ratios_pass_unchanged_ceiling": True,
    "causal_past_temporal_spectral_contract_inherited_unchanged": True,
    "trapped_sign_and_error_contract_inherited_unchanged": True,
    "source_solver_transaction_retry_boundary_and_health_stops_inherited_unchanged": True,
}
EXPECTED_INPUT = {
    "eligible_amplitudes_in_order": ["5/2", "3"],
    "methods_in_order": ["RK4", "SSPRK3"],
    "resolutions_in_order": [1025, 2049, 4097],
    "expected_run_input_count": 12,
    "projected_state_hashes_must_match_HLT5": True,
    "physical_u_p_q_and_initial_profiles_unchanged": True,
    "expanded_PROTO8_run_config_hashes_must_be_serialized": True,
    "all_twelve_accepted_state_source_prechecks_must_pass": True,
    "both_amplitudes_and_methods_must_pass_the_PROTO8_t0_common_event": True,
    "first_eligible_candidate_stops_later_candidates": True,
}
EXPECTED_NAMESPACE = {
    "calibration_output_root": "runs/fgc-2-sf1/proto8/calibration",
    "holdout_output_root": "runs/fgc-2-sf1/proto8/holdout",
    "both_roots_must_be_absent_or_empty_at_authorization": True,
    "authorization_must_not_create_either_root": True,
    "runner_must_refuse_overwrite": True,
}
EXPECTED_PROOF_KEYS = {
    "PROTO8_and_PRO8_FRZ1_must_validate_from_immutable_checkpoint",
    "HLT5_ID2_and_CAL4_plan_must_be_canonical_and_hash_matched_from_immutable_checkpoint",
    "PROTO8_run_plan_must_differ_from_PROTO7_only_by_frozen_common_event_ownership_and_namespace",
    "outer_projector_defect_must_remain_public_but_not_owned",
    "interior_constraint_defect_above_enclosure_must_fail",
    "accumulated_roundoff_must_not_change_raw_magnitude_guards",
    "coarse_only_spectral_miss_requires_resolved_finest_pair_and_contracting_tails",
    "medium_or_finest_absolute_spectral_miss_must_fail",
    "noncontracting_adjacent_tail_must_fail",
    "inherited_PROTO7_transaction_and_source_implementations_must_remain_hash_bound",
    "all_twelve_projected_inputs_must_reconstruct_deterministically",
    "all_twelve_t0_accepted_source_prechecks_must_pass",
    "all_four_PROTO8_t0_common_events_must_pass",
    "authorized_runner_must_be_hash_bound_and_refuse_drift_overwrite_or_dirty_tracked_state",
    "fresh_namespaces_must_be_observed_empty_without_creation",
    "canonical_hash_bound_result_required",
    "no_PROTO8_trajectory_may_be_read_or_advanced",
}
EXPECTED_CLAIMS = {
    "PROTO8_successor_runtime_compositor_implemented": True,
    "PROTO8_fresh_GR0_dynamic_calibration_authorized": True,
    "PROTO8_resolved_holdout_manifest_authorized": False,
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
    REPOSITORY / "scripts/run_fgc_gr0_calibration.py",
    REPOSITORY / "scripts/run_fgc_gr0_calibration_v6.py",
    REPOSITORY / "scripts/run_fgc_gr0_calibration_v7.py",
    REPOSITORY / "scripts/run_fgc_gr0_calibration_v8.py",
    Path(__file__).resolve(),
)


def _path(name: str, value: Any, expected: str) -> Path:
    if value != expected:
        raise ValueError(f"HLT6 {name} differs")
    path = REPOSITORY / expected
    if not path.is_file():
        raise ValueError(f"HLT6 {name} is absent")
    return path


def validate_run_plan(path: Path) -> dict[str, Any]:
    with path.open("rb") as handle:
        raw = tomllib.load(handle)
    if (
        raw.get("schema_version") != 1
        or raw.get("artifact_id") != "FGC-1-CAL5-RUN1-PLAN"
        or raw.get("project_version") != PROJECT_VERSION
        or raw.get("protocol_config") != EXPECTED_PATHS["protocol_config"]
        or raw.get("protocol_freeze_result") != EXPECTED_PATHS["protocol_freeze_result"]
        or raw.get("runtime_authorization_result") != "results/fgc-1-hlt6-mon6.json"
        or raw.get("static_input_result") != EXPECTED_PATHS["static_input_result"]
    ):
        raise ValueError("CAL5 run-plan identity or ownership differs")
    if raw.get("scope") != {
        "target_protocol": "FGC-2-SF1-PROTO8",
        "branch": "GR-0",
        "role": "fresh_outcome_neutral_dynamic_amplitude_recalibration_after_common_event_ownership_repair",
        "physical_equations": "unredefined_GR0_specialization_of_REF1",
        "candidate_action_fields_or_health_stops_forbidden": True,
        "retained_EFT_interpretation": False,
    }:
        raise ValueError("CAL5 scope differs")
    numerics = raw.get("numerics", {})
    expected_delta = {
        "common_event_constraint_domain": "evolution_owned_rows_excluding_fixed_projector_and_stencil_reach",
        "common_event_full_and_owned_constraint_residuals_serialized": True,
        "common_event_constraint_roundoff_formula": "4096*binary64_epsilon*(1+accepted_stage_count)",
        "common_event_constraint_roundoff_changes_order_only": True,
        "common_event_spatial_coarse_role": "public_convergence_witness",
        "common_event_spatial_medium_and_finest_absolute_budgets_required": True,
        "common_event_spatial_all_adjacent_tail_ratios_required": True,
    }
    if any(numerics.get(key) != value for key, value in expected_delta.items()):
        raise ValueError("CAL5 PROTO8 numerical ownership differs")
    if raw.get("provenance") != {
        "output_root": "runs/fgc-2-sf1/proto8/calibration",
        "campaign_manifest": "runs/fgc-2-sf1/proto8/calibration/manifest.json",
        "append_only_event_log_name": "events.jsonl",
        "checkpoint_extension": ".npz",
        "fresh_output_root_required_at_campaign_start": True,
        "runner_refuses_overwrite": True,
        "resume_requires_matching_plan_and_manifest_hashes": True,
        "PROTO7_outputs_forbidden_as_PROTO8_outcomes": True,
        "canonical_reduced_result_is_tracked_only_after_campaign_classification": True,
    }:
        raise ValueError("CAL5 provenance differs")

    normalized = deepcopy(raw)
    normalized["artifact_id"] = "FGC-1-CAL4-RUN1-PLAN"
    normalized["protocol_config"] = "configs/fgc/fgc-2-sf1-protocol-v7.toml"
    normalized["protocol_freeze_result"] = "results/fgc-1-pro7-frz1.json"
    normalized["runtime_authorization_result"] = "results/fgc-1-hlt5-mon5.json"
    normalized["scope"]["target_protocol"] = "FGC-2-SF1-PROTO7"
    normalized["scope"]["role"] = (
        "fresh_outcome_neutral_dynamic_amplitude_recalibration_after_general_transaction_repair"
    )
    for key in expected_delta:
        normalized["numerics"].pop(key)
    normalized["provenance"] = {
        "output_root": "runs/fgc-2-sf1/proto7/calibration",
        "campaign_manifest": "runs/fgc-2-sf1/proto7/calibration/manifest.json",
        "append_only_event_log_name": "events.jsonl",
        "checkpoint_extension": ".npz",
        "fresh_output_root_required_at_campaign_start": True,
        "runner_refuses_overwrite": True,
        "resume_requires_matching_plan_and_manifest_hashes": True,
        "PROTO6_outputs_forbidden_as_PROTO7_outcomes": True,
        "canonical_reduced_result_is_tracked_only_after_campaign_classification": True,
    }
    baseline_blob = hlt4._git_blob(
        CHECKPOINT_COMMIT,
        EXPECTED_PATHS["predecessor_run_plan"],
        "immutable CAL4 plan",
    )
    baseline = tomllib.loads(baseline_blob.decode("utf-8"))
    if normalized != baseline:
        raise ValueError("CAL5 changes more than the declared PROTO8 delta")
    return raw


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    with path.open("rb") as handle:
        raw = tomllib.load(handle)
    expected_root = {
        "schema_version", "artifact_id", "project_version", "metric_signature",
        "riemann_convention", *EXPECTED_PATHS, "scope", "immutable_lineage",
        "common_event_runtime", "input_freeze", "namespace", "proof_contract",
        "claims",
    }
    if set(raw) != expected_root:
        raise ValueError("HLT6 root keys differ")
    if (
        raw["schema_version"] != 1
        or raw["artifact_id"] != ARTIFACT_ID
        or raw["project_version"] != PROJECT_VERSION
        or raw["metric_signature"] != "-+++"
        or raw["riemann_convention"] != "plus_partial_mu_gamma_nu"
    ):
        raise ValueError("HLT6 identity or convention differs")
    paths = {
        name: _path(name, raw[name], expected)
        for name, expected in EXPECTED_PATHS.items()
    }
    if raw["scope"] != EXPECTED_SCOPE:
        raise ValueError("HLT6 scope differs")
    if raw["immutable_lineage"] != EXPECTED_LINEAGE:
        raise ValueError("HLT6 immutable lineage differs")
    if raw["common_event_runtime"] != EXPECTED_COMMON_EVENT:
        raise ValueError("HLT6 common-event runtime differs")
    if raw["input_freeze"] != EXPECTED_INPUT:
        raise ValueError("HLT6 input freeze differs")
    if raw["namespace"] != EXPECTED_NAMESPACE:
        raise ValueError("HLT6 namespace differs")
    if set(raw["proof_contract"]) != EXPECTED_PROOF_KEYS or any(
        value is not True for value in raw["proof_contract"].values()
    ):
        raise ValueError("HLT6 proof contract differs")
    if raw["claims"] != EXPECTED_CLAIMS:
        raise ValueError("HLT6 claims differ")
    return {"raw": raw, "paths": paths, "run_plan": validate_run_plan(paths["run_plan_config"])}


def _immutable_lineage(loaded: Mapping[str, Any]) -> dict[str, Any]:
    lineage = loaded["raw"]["immutable_lineage"]
    commit = lineage["checkpoint_commit"]
    if subprocess.run(
        ["git", "merge-base", "--is-ancestor", commit, "HEAD"],
        cwd=REPOSITORY,
        check=False,
        capture_output=True,
    ).returncode != 0:
        raise ValueError("HLT6 checkpoint is not an ancestor of HEAD")
    immutable_names = (
        "protocol_config", "protocol_freeze_config", "protocol_freeze_result",
        "predecessor_runtime_result", "static_input_result", "predecessor_run_plan",
    )
    blobs = {
        name: hlt4._git_blob(commit, EXPECTED_PATHS[name], name)
        for name in immutable_names
    }
    for name, blob in blobs.items():
        if sha256(blob).hexdigest() != lineage[f"{name}_sha256"]:
            raise ValueError(f"immutable {name} hash differs")
        current = REPOSITORY / EXPECTED_PATHS[name]
        if sha256(current.read_bytes()).hexdigest() != lineage[f"{name}_sha256"]:
            raise ValueError(f"current {name} differs from immutable checkpoint")

    protocol = validate_sf1_protocol_v8(
        tomllib.loads(blobs["protocol_config"].decode("utf-8"))
    )
    freeze = hlt4._load_json_bytes(blobs["protocol_freeze_result"], "immutable PRO8-FRZ1")
    predecessor = hlt4._load_json_bytes(
        blobs["predecessor_runtime_result"], "immutable HLT5"
    )
    static = hlt4._load_json_bytes(blobs["static_input_result"], "immutable ID2")
    inputs = predecessor.get("artifact_payload", {}).get("frozen_run_inputs", [])
    if (
        protocol.get("artifact_id") != "FGC-2-SF1-PROTO8"
        or protocol.get("outcome_neutral_contract_validated") is not True
        or freeze.get("artifact_id") != "FGC-1-PRO8-FRZ1"
        or freeze.get("gate_status", {}).get("PROTO8_outcome_neutral_protocol_frozen")
        is not True
        or freeze.get("gate_status", {}).get("PROTO8_fresh_GR0_dynamic_calibration_authorized")
        is not False
        or predecessor.get("artifact_id") != "FGC-1-HLT5-MON5"
        or predecessor.get("gate_status", {}).get("PROTO7_successor_runtime_compositor_implemented")
        is not True
        or predecessor.get("gate_status", {}).get("PROTO7_fresh_GR0_dynamic_calibration_authorized")
        is not True
        or predecessor.get("gate_status", {}).get("FGCQR_holdout_execution_authorized")
        is not False
        or static.get("gate_status", {}).get("all_case_static_ledger_completed") is not True
        or len(inputs) != 12
        or any(value is not False for value in freeze.get("nonclaims", {}).values())
        or any(value is not False for value in predecessor.get("nonclaims", {}).values())
    ):
        raise ValueError("immutable HLT6 predecessor gate differs")
    for relative, expected in predecessor.get("implementation_sha256", {}).items():
        current = REPOSITORY / relative
        if not current.is_file() or hlt4._sha(current) != expected:
            raise ValueError(f"inherited PROTO7 runtime drifted: {relative}")
    return {
        "checkpoint_commit": commit,
        "checkpoint_is_ancestor_of_HEAD": True,
        "protocol": {
            "artifact_id": protocol["artifact_id"],
            "protocol_version": protocol["protocol_version"],
            "semantic_holdout_contract_sha256": protocol["semantic_holdout_contract_sha256"],
        },
        "freeze": {
            "artifact_id": freeze["artifact_id"],
            "result_sha256": lineage["protocol_freeze_result_sha256"],
            "claims_all_false": all(
                value is False for value in freeze["artifact_payload"]["protocol_certificate"]["claims"].values()
            ),
        },
        "predecessor_runtime": {
            "artifact_id": predecessor["artifact_id"],
            "result_sha256": lineage["predecessor_runtime_result_sha256"],
            "implementation_files_hash_matched": len(predecessor["implementation_sha256"]),
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


def _budget(
    sample_count: int,
    *,
    field_fraction: float,
    derivative_fraction: float,
    admitted: bool,
) -> SpectralPowerBudget:
    nyquist = sample_count // 2
    return SpectralPowerBudget(
        sample_count=sample_count,
        top_band_first_bin=max(0, nyquist - 2),
        nyquist_bin=nyquist,
        total_field_power=1.0,
        top_band_field_power_fraction=field_fraction,
        total_derivative_weighted_power=1.0,
        top_band_derivative_weighted_power_fraction=derivative_fraction,
        rms_angular_scale=1.0,
        rms_scale_over_cutoff=0.1,
        last_occupied_bin=nyquist,
        last_bin_alias_diagnostic=True,
        individual_admission_passed=admitted,
    )


def _adapter_controls() -> dict[str, Any]:
    grid = UniformRadialGrid(0.0, 8.0, 65)
    state = _flat_state(grid, 4)
    component = PROTO4_CONSTRAINT_ORDER.index("reduction_alpha")
    outer_q = state.q.copy()
    outer_q[-5, 0] += 1.0e-4
    outer = owned_semidiscrete_constraint_snapshot(
        EvolutionState(state.u, state.p, outer_q),
        grid,
        coordinate_time=0.25,
        diagnostic_spatial_order=4,
        fixed_outer_rows=4,
        accepted_stage_count=20,
    )
    interior_q = state.q.copy()
    interior_q[20, 0] += 1.0e-6
    interior = owned_semidiscrete_constraint_snapshot(
        EvolutionState(state.u, state.p, interior_q),
        grid,
        coordinate_time=0.25,
        diagnostic_spatial_order=4,
        fixed_outer_rows=4,
        accepted_stage_count=20,
    )
    outer_pass = (
        outer.full_domain_norms.component_infinity[component] > 1.0e-6
        and outer.owned_domain_norms.component_infinity[component]
        < outer.accumulated_roundoff_enclosure
        and outer.full_component_maximum_radii[component] > outer.last_owned_radius
    )
    interior_veto = (
        interior.owned_domain_norms.component_infinity[component]
        > 1000.0 * interior.accumulated_roundoff_enclosure
        and interior.owned_component_maximum_radii[component]
        <= interior.last_owned_radius
    )

    counts = (33, 65, 129)
    resolved = (
        {"field": _budget(33, field_fraction=1.0e-4, derivative_fraction=1.0e-2, admitted=False)},
        {"field": _budget(65, field_fraction=1.0e-6, derivative_fraction=1.0e-4, admitted=True)},
        {"field": _budget(129, field_fraction=1.0e-8, derivative_fraction=1.0e-6, admitted=True)},
    )
    coarse_witness = finest_pair_nested_spectral_admission(counts, resolved)
    medium_miss = list(resolved)
    medium_miss[1] = {"field": replace(medium_miss[1]["field"], individual_admission_passed=False)}
    medium_veto = finest_pair_nested_spectral_admission(counts, medium_miss)
    tail_miss = list(resolved)
    tail_miss[2] = {
        "field": _budget(129, field_fraction=5.0e-7, derivative_fraction=5.0e-5, admitted=True)
    }
    tail_veto = finest_pair_nested_spectral_admission(counts, tail_miss)
    all_passed = (
        outer_pass
        and interior_veto
        and not coarse_witness.coarsest_individual_budget_passed
        and coarse_witness.finest_pair_individual_budgets_passed
        and coarse_witness.every_nested_tail_ratio_passed
        and coarse_witness.admission_passed
        and not medium_veto.admission_passed
        and not tail_veto.admission_passed
    )
    if not all_passed:
        raise RuntimeError("PROTO8 common-event adapter controls failed")
    return {
        "outer_projector_defect": {
            "full_raw_component_infinity": outer.full_domain_norms.component_infinity[component],
            "owned_raw_component_infinity": outer.owned_domain_norms.component_infinity[component],
            "accumulated_roundoff_enclosure": outer.accumulated_roundoff_enclosure,
            "full_maximum_radius": outer.full_component_maximum_radii[component],
            "last_owned_radius": outer.last_owned_radius,
            "public_but_not_owned_control_passed": True,
        },
        "interior_defect": {
            "owned_raw_component_infinity": interior.owned_domain_norms.component_infinity[component],
            "accumulated_roundoff_enclosure": interior.accumulated_roundoff_enclosure,
            "owned_maximum_radius": interior.owned_component_maximum_radii[component],
            "interior_veto_control_passed": True,
        },
        "spatial_spectrum": {
            "coarse_witness_finest_pair_control_passed": True,
            "medium_absolute_miss_vetoed": True,
            "noncontracting_tail_vetoed": True,
        },
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
    expected_order = [
        (amplitude, method, count)
        for amplitude in ("5/2", "3")
        for method in ("RK4", "SSPRK3")
        for count in (1025, 2049, 4097)
    ]
    if [
        (item["amplitude"], item["method"], item["point_count"])
        for item in predecessor_inputs
    ] != expected_order:
        raise ValueError("immutable HLT5 input order differs")

    frozen: list[dict[str, Any]] = []
    source_prechecks: list[dict[str, Any]] = []
    states: dict[tuple[str, str], list[tuple[int, EvolutionState, Any]]] = {}
    for predecessor in predecessor_inputs:
        amplitude = predecessor["amplitude"]
        method = predecessor["method"]
        count = predecessor["point_count"]
        order = 4 if method == "RK4" else 2
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
        if state_hash != predecessor["projected_state_sha256"]:
            raise ValueError("PROTO8 projected state differs from immutable HLT5")
        operator = Proto7GR0EvolutionOperator(
            initial.grid,
            spatial_order=order,
            ko_dissipation=float(Q(numerics["ko_dissipation"])),
            raw_tolerance=float(Q(thresholds["source_residual_infinity_max"])),
            kinetic_condition_maximum=float(Q(thresholds["kinetic_condition_number_max"])),
        )
        rhs = operator(0.0, state)
        residual = float(rhs.diagnostics["source_residual_infinity"])
        threshold = float(Q(thresholds["source_residual_infinity_max"]))
        if residual >= threshold:
            raise ValueError("PROTO8 t0 accepted-state source precheck failed")
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
        record = dict(predecessor)
        record["plan_sha256"] = plan_hash
        record["protocol_sha256"] = protocol_hash
        record["common_event_constraint_domain"] = (
            "evolution_owned_rows_excluding_projector_and_stencil_reach"
        )
        record["common_event_constraint_roundoff"] = (
            "4096*binary64_epsilon*(1+accepted_stage_count)_order_only"
        )
        record["common_event_spatial_spectral_roles"] = (
            "coarse_witness_medium_fine_absolute_all_adjacent_contraction"
        )
        record["trajectory_advanced"] = False
        record.pop("expanded_run_config_sha256", None)
        record["expanded_run_config_sha256"] = hlt4._digest(record)
        frozen.append(record)
        states.setdefault((amplitude, method), []).append((count, state, initial.grid))

    common_events: list[dict[str, Any]] = []
    for amplitude in ("5/2", "3"):
        for method in ("RK4", "SSPRK3"):
            records = sorted(states[(amplitude, method)], key=lambda item: item[0])
            assessed = proto8_gr0_common_event(
                [item[1] for item in records],
                [item[2] for item in records],
                accepted_stage_counts=[0, 0, 0],
                method=method,
                coordinate_time=0.0,
                cutoff=float(Q(physical["cutoff_Lambda"])),
                measurement_radius_maximum=float(Q(physical["measurement_radius_maximum"])),
                taper_fraction=float(Q(numerics["proper_spectral_taper_fraction"])),
                fixed_outer_rows=numerics["fixed_outer_rows"],
            )
            if not assessed.admission_passed:
                raise ValueError("PROTO8 initial common event failed")
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
                    "coarsest_spatial_individual_budget_passed": (
                        assessed.spatial_spectral_admission.coarsest_individual_budget_passed
                    ),
                    "finest_pair_spatial_individual_budgets_passed": (
                        assessed.spatial_spectral_admission.finest_pair_individual_budgets_passed
                    ),
                    "all_adjacent_spatial_tail_ratios_passed": (
                        assessed.spatial_spectral_admission.every_nested_tail_ratio_passed
                    ),
                    "admission_passed": True,
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
            != [raw["namespace"]["calibration_output_root"], raw["namespace"]["holdout_output_root"]]
            or any(item.get("empty") is not True for item in records)
        ):
            raise ValueError("historical HLT6 namespace evidence differs")
        return dict(historical)
    records = []
    for key in ("calibration_output_root", "holdout_output_root"):
        relative = raw["namespace"][key]
        path = REPOSITORY / relative
        exists = path.exists()
        if exists and (not path.is_dir() or any(path.iterdir())):
            raise ValueError(f"HLT6 fresh namespace is nonempty: {relative}")
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
    namespace = _namespace_precondition(raw, historical_namespace_evidence)
    plan_hash = hlt4._sha(loaded["paths"]["run_plan_config"])
    input_manifest_hash = hlt4._digest(inputs)
    campaign_id = hlt4._digest(
        {
            "artifact_id": loaded["run_plan"]["artifact_id"],
            "run_plan_sha256": plan_hash,
            "input_manifest_sha256": input_manifest_hash,
            "ordered_amplitudes": ["5/2", "3"],
            "output_root": raw["namespace"]["calibration_output_root"],
        }
    )
    public_lineage = dict(lineage)
    public_lineage.pop("predecessor_frozen_inputs")
    return {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "classification": "pre_trajectory_PROTO8_GR0_common_event_ownership_runtime_and_fresh_calibration_authorization",
        "generated_by": "scripts/reproduce_fgc_hlt6_mon6.py",
        "derivation_document": hlt4._rel(OWNER_DOCUMENT),
        "derivation_document_sha256": hlt4._sha(OWNER_DOCUMENT),
        "source_config_sha256": {
            hlt4._rel(config_path): hlt4._sha(config_path),
            hlt4._rel(loaded["paths"]["run_plan_config"]): plan_hash,
        },
        "predecessor_sha256": {
            EXPECTED_PATHS[name]: raw["immutable_lineage"][f"{name}_sha256"]
            for name in (
                "protocol_config", "protocol_freeze_config", "protocol_freeze_result",
                "predecessor_runtime_result", "static_input_result", "predecessor_run_plan",
            )
        },
        "implementation_sha256": {
            hlt4._rel(path): hlt4._sha(path) for path in IMPLEMENTATION
        },
        "scope_bindings": dict(raw["scope"]),
        "gate_status": dict(raw["claims"]),
        "artifact_payload": {
            "immutable_lineage": public_lineage,
            "common_event_adapter_controls": controls,
            "inherited_PROTO7_transaction": {
                "source_solver_transaction_retry_and_non_source_stops_unchanged": True,
                "inherited_implementation_files_hash_matched": public_lineage[
                    "predecessor_runtime"
                ]["implementation_files_hash_matched"],
            },
            "frozen_run_plan": {
                "artifact_id": loaded["run_plan"]["artifact_id"],
                "run_plan_sha256": plan_hash,
                "campaign_id": campaign_id,
                "input_manifest_sha256": input_manifest_hash,
                "ordered_amplitudes": ["5/2", "3"],
                "first_eligible_candidate_stops_later_candidates": True,
                "expanded_input_count": len(inputs),
            },
            "frozen_run_inputs": inputs,
            "accepted_state_source_prechecks": source_prechecks,
            "initial_common_events": common_events,
            "namespace_precondition": namespace,
            "decision": "HLT6_implements_only_PROTO8_common_event_ownership_and_authorizes_one_fresh_GR0_recalibration",
            "epistemic_boundary": {
                "PROTO7_calibration_history_disclosed": True,
                "PROTO8_trajectory_read": False,
                "PROTO8_trajectory_advanced": False,
                "synthetic_failures_used_only_for_adapter_controls": True,
                "t0_physical_inputs_reconstructed": True,
                "trapped_or_collapse_outcome_classified": False,
                "SGBL_or_FGCQR_outcome_read": False,
                "fresh_GR0_recalibration_is_now_authorized_but_not_completed": True,
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
    return hlt4._load_json_bytes(path.read_bytes(), "HLT6 result")


def _verify_postlaunch(
    loaded: Mapping[str, Any],
    observed: Mapping[str, Any],
    result_path: Path,
) -> None:
    raw = loaded["raw"]
    historical = observed["artifact_payload"]["namespace_precondition"]
    _namespace_precondition(raw, historical)
    calibration_root = REPOSITORY / raw["namespace"]["calibration_output_root"]
    holdout_root = REPOSITORY / raw["namespace"]["holdout_output_root"]
    if not calibration_root.is_dir():
        raise ValueError("postlaunch HLT6 calibration directory is absent")
    if holdout_root.exists() and (
        not holdout_root.is_dir() or any(holdout_root.iterdir())
    ):
        raise ValueError("HLT6 does not authorize a nonempty holdout namespace")
    manifest = hlt4._load_json_bytes(
        (calibration_root / "manifest.json").read_bytes(), "CAL5 launch manifest"
    )
    frozen = observed["artifact_payload"]["frozen_run_plan"]
    expected = {
        "schema_version": 1,
        "runner_id": "FGC-1-CAL5-RUN1-RUNNER",
        "campaign_id": frozen["campaign_id"],
        "plan_path": EXPECTED_PATHS["run_plan_config"],
        "plan_sha256": frozen["run_plan_sha256"],
        "authorization_path": hlt4._rel(result_path),
        "authorization_sha256": hlt4._sha(result_path),
        "input_manifest_sha256": frozen["input_manifest_sha256"],
        "ordered_amplitudes": frozen["ordered_amplitudes"],
        "implementation_sha256": observed["implementation_sha256"],
        "outcome_fields_present_at_creation": False,
    }
    for key, value in expected.items():
        if manifest.get(key) != value:
            raise ValueError(f"postlaunch CAL5 manifest field differs: {key}")
    commit = manifest.get("git_commit")
    if not isinstance(commit, str) or subprocess.run(
        ["git", "merge-base", "--is-ancestor", commit, "HEAD"],
        cwd=REPOSITORY,
        check=False,
        capture_output=True,
    ).returncode != 0:
        raise ValueError("postlaunch CAL5 Git checkpoint is invalid")


def verify_canonical(path: Path, config_path: Path = DEFAULT_CONFIG) -> None:
    observed = load_canonical_result(path)
    loaded = load_config(config_path)
    calibration_root = REPOSITORY / loaded["raw"]["namespace"]["calibration_output_root"]
    if calibration_root.exists() and any(calibration_root.iterdir()):
        _verify_postlaunch(loaded, observed, path)
        expected = record(
            config_path,
            historical_namespace_evidence=observed["artifact_payload"]["namespace_precondition"],
        )
    else:
        expected = record(config_path)
    if observed != hlt4._serial(expected):
        raise ValueError("HLT6 output differs from a fresh reproduction")


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
        "(PROTO8 runtime=true; fresh calibration authorized=true; trajectory_read=false)"
    )


if __name__ == "__main__":
    main()
