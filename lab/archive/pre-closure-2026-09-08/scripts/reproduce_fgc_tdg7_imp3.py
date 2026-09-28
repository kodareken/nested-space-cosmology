#!/usr/bin/env python3
"""Reproduce the synthetic FGC-1-TDG7-IMP3 runtime-repair certificate.

This artifact binds the sealed PREF23 authorization and the current TDG7
adapter/test bytes, then executes only synthetic states.  It never opens a
CAL11 or PROTO14 history array, checkpoint, campaign member, or runs namespace.
"""

from __future__ import annotations

import argparse
from collections import Counter
from copy import deepcopy
from dataclasses import asdict, replace
from hashlib import sha256
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import tomllib
from typing import Any, Callable, Mapping

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts.run_fgc_gr0_calibration import NormalFlowTracers  # noqa: E402
from recursive_horizons.fgc.evolution.boundary_domain import (  # noqa: E402
    BoundaryGeometry,
    CausalBudgetState,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    COMPARATOR_METHOD,
    PRIMARY_METHOD,
    EvolutionRHS,
    EvolutionState,
    array_content_sha256,
)
from recursive_horizons.fgc.evolution.proto5_runtime import (  # noqa: E402
    GR0RuntimeStageTransaction,
    GR0UniversalThresholds,
)
from recursive_horizons.fgc.evolution.tdg6_temporal_admission_design import (  # noqa: E402
    TDG6_COMPLETE_STATE_CHANNELS,
    TDG6_MINIMUM_MACRO_STEP,
)
from recursive_horizons.fgc.evolution import (  # noqa: E402
    tdg6_temporal_admission_runtime as tdg6,
    tdg7_stage_safe_runtime as runtime,
)


ARTIFACT_ID = "FGC-1-TDG7-IMP3"
PROJECT_VERSION = "0.11.0"
CONFIG = REPOSITORY / "configs/fgc/fgc-1-tdg7-imp3.toml"
OUTPUT = REPOSITORY / "results/fgc-1-tdg7-imp3.json"
DOCUMENT = REPOSITORY / "docs/fgc-tdg7-imp3.md"
PREDECESSOR_COMMIT = "e9289ce479a2cd720d6c58f2ca3d653f573d8eec"

PREF23_PATHS = {
    "PREF23_config": "configs/fgc/fgc-1-tdg7-pref23.toml",
    "PREF23_document": "docs/fgc-tdg7-pref23.md",
    "PREF23_result": "results/fgc-1-tdg7-pref23.json",
    "PREF23_reproducer": "scripts/reproduce_fgc_tdg7_pref23.py",
    "PREF23_reproduction_test": "tests/test_fgc_tdg7_pref23_reproduction.py",
    "PREF23_binder_module": (
        "src/recursive_horizons/fgc/evolution/tdg7_stage_safe_lattice_theorem.py"
    ),
}
PREF23_HASHES = {
    "PREF23_config_sha256": "410e6e3451f53dd6d262a2e1c409773955f5d27621b3d95c5e977d7e244643e7",
    "PREF23_document_sha256": "ef0673f179139ff10c6f3b906d3fd362ac2fd53fd8d1ae265fbb47a70b91bfc7",
    "PREF23_result_sha256": "09b40197f043026c59fb1efff2b6d24197ad7aeca03b784da32b602165e97843",
    "PREF23_reproducer_sha256": "4055de98048d7f40cee67e9768fd3b2d453406e26d85f982a205f95a882ddf47",
    "PREF23_reproduction_test_sha256": "32bf4f2d3a932937d08cd926654e57b78accade4ad1c8f79dc782d0a6b269534",
    "PREF23_binder_module_sha256": "f9d3ec171d9250b8998d53459de9ee2fe01ec1c4315c25782d97f4b85056454f",
}
CORE_PATHS = {
    "core_runtime": "src/recursive_horizons/fgc/evolution/tdg7_stage_safe_runtime.py",
    "core_runtime_test": "tests/test_fgc_tdg7_stage_safe_runtime.py",
}
CORE_HASHES = {
    "core_runtime_sha256": "751ceb0e209d65ee9e215188326cd20add7a189730b2fadbe1974bcb0d0150a1",
    "core_runtime_test_sha256": "2665b88f3549739fc97f1d838d82abacb3da0b52cbfa7497f02b14f06724e929",
}

START = 23.0 / 16.0
TARGET = 24.0 / 16.0
CAP = 0.125

EXPECTED_SCOPE = {
    "target": "FGC-2-SF1-TDG7",
    "role": "authorized_stage_safe_runtime_repair_implementation_and_synthetic_qualification",
    "calibration_branch": "GR-0",
    "sealed_predecessor_read_only": True,
    "PREF23_compact_artifacts_read_only": True,
    "CAL11_history_arrays_loaded": False,
    "PROTO14_history_arrays_loaded": False,
    "campaign_checkpoint_loaded": False,
    "campaign_checkpoint_resumed": False,
    "production_state_advanced": False,
    "runs_namespace_created": False,
    "synthetic_controls_only": True,
    "SGBL_trajectory_read": False,
    "FGCQR_trajectory_read": False,
    "physical_or_candidate_question_answered": False,
}
EXPECTED_RUNTIME_CONTRACT = {
    "one_immutable_shared_TDG7_plan_owns_one_two_four_paths": True,
    "actual_RK4_and_SSPRK3_stage_names_and_times_rechecked_against_plan": True,
    "RK4_total_shadow_stage_records": 35,
    "RK4_fine_shadow_stage_records": 20,
    "SSPRK3_total_shadow_stage_records": 28,
    "SSPRK3_fine_shadow_stage_records": 16,
    "typed_coordinate_lattice_stops_precede_shadow_source_projector_tracer_monitor_or_state_work": True,
    "all_typed_lattice_stops_are_exercised_with_zero_hooks": True,
    "TDG6_minimum_width_is_frozen_and_not_caller_overridable": True,
    "fresh_entrypoint_rejects_active_retry_ledger": True,
    "only_complete_fine_path_commits_atomically": True,
    "outer_and_medium_paths_are_evidence_only": True,
}
EXPECTED_TDG6_INHERITANCE = {
    "complete_channel_count": 18,
    "exact_zero_classification_preserved": True,
    "one_channel_resolved_order_failure_vetoes_all_of_admission": True,
    "large_convergent_difference_passes_without_magnitude_tolerance": True,
    "large_convergent_debit_remains_outward_and_greater_than_nine": True,
    "source_only_late_shadow_stop_remains_PROTO7_owned": True,
    "non_source_late_shadow_stop_remains_nonretry_terminal_evidence": True,
    "tracer_commit_exception_rolls_back_real_monitor_causal_and_tracer_boundaries": True,
    "valid_TDG6_retry_is_replanned_through_TDG7": True,
    "retry_successor_binds_predecessor_plan_retry_and_exact_updated_ledger": True,
    "retry_successor_rejects_mismatched_method_over_limit_and_forged_evidence_and_prefix": True,
    "two_minimum_to_minimum_successor_is_valid": True,
    "TDG6_owns_normal_minimum_macro_step_exhaustion": True,
    "TDG7_lattice_exhaustion_is_defensive_for_forged_required_wrapper_only": True,
}
EXPECTED_SYNTHETIC = {
    "both_methods_execute_all_seven_shadow_proposals": True,
    "all_stage_times_are_recorded_as_exact_binary64_hex": True,
    "all_typed_lattice_stops_precede_work": True,
    "fine_only_atomic_commit_required": True,
    "TDG6_interval_and_retry_controls_required": True,
    "late_shadow_ownership_and_tracer_rollback_required": True,
    "lineage_hash_claim_and_result_mutations_must_fail_closed": True,
}
EXPECTED_LIMITATIONS = {
    "synthetic_controls_are_not_PROTO14_trajectory_execution": True,
    "historical_CAL11_invalid_runtime_is_not_physics": True,
    "current_runtime_repair_does_not_mutate_or_resume_historical_PROTO14": True,
    "stateless_adapter_cannot_detect_discarded_external_retry": True,
    "future_protocol_cursor_required_before_calibration": True,
    "fresh_vs_retry_protocol_discriminator_implemented": False,
    "prospective_replacement_protocol_design_only_no_protocol_file_runtime_or_run": True,
    "no_new_protocol_or_output_namespace_is_created": True,
    "no_calibration_or_candidate_eligibility_is_derived": True,
    "no_SGBL_or_FGCQR_or_DEF1_execution_is_authorized": True,
    "no_retained_EFT_or_physical_transition_or_cosmological_claim_is_derived": True,
}
EXPECTED_SUCCESSOR = {
    "TDG7_independent_binder_completed": True,
    "TDG7_runtime_repair_implementation_authorized": True,
    "TDG7_runtime_repair_implemented": True,
    "TDG7_runtime_repair_synthetic_qualification_passed": True,
    "replacement_protocol_freeze_design_authorized": True,
    "fresh_replacement_protocol_frozen": False,
    "new_protocol_authorized": False,
    "new_output_namespace_authorized": False,
    "PROTO14_runtime_mutation_authorized": False,
    "PROTO14_terminal_checkpoint_may_resume": False,
    "new_protocol_or_output_namespace_authorized": False,
    "fresh_GR0_calibration_authorized": False,
    "GR0_case_eligible": False,
    "SGBL_execution_authorized": False,
    "FGCQR_holdout_execution_authorized": False,
    "DEF1_execution_authorized": False,
    "retained_EFT_evolution_authorized": False,
    "physical_transition_claim_authorized": False,
}
CLAIM_TRUE = {
    "CAL11_terminal_invalid_runtime_result_preserved",
    "TDG7_independent_binder_completed",
    "TDG7_runtime_repair_implementation_authorized",
    "TDG7_runtime_repair_implemented",
    "TDG7_runtime_repair_synthetic_qualification_passed",
    "one_shared_stage_safe_plan_owns_all_paths",
    "typed_lattice_stops_precede_work",
    "TDG6_threshold_and_all_of_admission_unchanged",
    "stateless_adapter_cannot_detect_discarded_external_retry",
    "future_protocol_cursor_required_before_calibration",
    "replacement_protocol_freeze_design_authorized",
}
CLAIM_FALSE = {
    "fresh_vs_retry_protocol_discriminator_implemented",
    "fresh_replacement_protocol_frozen",
    "PROTO14_runtime_mutation_authorized",
    "PROTO14_terminal_checkpoint_may_resume",
    "new_protocol_authorized",
    "new_output_namespace_authorized",
    "fresh_GR0_calibration_authorized",
    "fresh_GR0_dynamic_calibration_completed",
    "GR0_case_eligible",
    "SGBL_execution_authorized",
    "FGCQR_holdout_execution_authorized",
    "DEF1_execution_authorized",
    "retained_EFT_evolution_authorized",
    "physical_transition_claim_authorized",
    "physical_or_cosmological_claim_derived",
    "singularity_resolution_derived",
    "child_domain_or_topology_derived",
    "dark_sector_mechanism_derived",
    "varying_locally_measured_c_derived",
}


def canonical(value: object) -> bytes:
    return (
        json.dumps(value, indent=2, sort_keys=True, ensure_ascii=True, allow_nan=False)
        + "\n"
    ).encode("utf-8")


def sha(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _sha_bytes(value: bytes) -> str:
    return sha256(value).hexdigest()


def _rel(path: Path) -> str:
    return path.resolve().relative_to(REPOSITORY).as_posix()


def _git(*args: str, check: bool = True) -> subprocess.CompletedProcess[bytes]:
    return subprocess.run(
        ["git", *args],
        cwd=REPOSITORY,
        check=check,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )


def _strict_keys(name: str, mapping: Mapping[str, Any], expected: set[str]) -> None:
    if set(mapping) != expected:
        raise ValueError(f"IMP3 {name} keys differ")


def _all_true(name: str, mapping: Mapping[str, Any]) -> None:
    if not mapping or set(mapping.values()) != {True}:
        raise ValueError(f"IMP3 {name} must be all true")


def _expected_lineage() -> dict[str, Any]:
    return {
        "predecessor_commit": PREDECESSOR_COMMIT,
        **PREF23_HASHES,
        **CORE_HASHES,
        "predecessor_commit_must_be_ancestor_of_HEAD": True,
        "every_PREF23_blob_must_match_predecessor_and_worktree": True,
        "PREF23_compact_result_must_remain_canonical": True,
        "live_core_runtime_and_test_must_match_bound_hashes": True,
        "CAL11_invalid_runtime_status_preserved_without_loading_history": True,
    }


def validate_config_data(config: Mapping[str, Any]) -> None:
    _strict_keys(
        "config",
        config,
        {
            "schema_version",
            "artifact_id",
            "project_version",
            "metric_signature",
            "riemann_convention",
            "predecessor_commit",
            *PREF23_PATHS,
            *CORE_PATHS,
            "scope",
            "immutable_lineage",
            "stage_safe_runtime_contract",
            "tdg6_inheritance",
            "synthetic_qualification",
            "limitations",
            "proof_contract",
            "successor_boundary",
            "claims",
        },
    )
    expected_paths = {**PREF23_PATHS, **CORE_PATHS}
    if (
        config.get("schema_version") != 1
        or config.get("artifact_id") != ARTIFACT_ID
        or config.get("project_version") != PROJECT_VERSION
        or config.get("metric_signature") != "-+++"
        or config.get("riemann_convention") != "plus_partial_mu_gamma_nu"
        or config.get("predecessor_commit") != PREDECESSOR_COMMIT
        or any(config.get(key) != value for key, value in expected_paths.items())
        or config.get("scope") != EXPECTED_SCOPE
        or config.get("immutable_lineage") != _expected_lineage()
        or config.get("stage_safe_runtime_contract") != EXPECTED_RUNTIME_CONTRACT
        or config.get("tdg6_inheritance") != EXPECTED_TDG6_INHERITANCE
        or config.get("synthetic_qualification") != EXPECTED_SYNTHETIC
        or config.get("limitations") != EXPECTED_LIMITATIONS
        or config.get("successor_boundary") != EXPECTED_SUCCESSOR
    ):
        raise ValueError("IMP3 identity, lineage, contract, or boundary differs")
    _strict_keys(
        "proof contract",
        config["proof_contract"],
        {
            "predecessor_and_all_PREF23_hashes_must_match_git_and_worktree",
            "current_core_runtime_and_test_hashes_must_match",
            "canonical_PREF23_boundary_must_authorize_repair_only",
            "all_actual_method_stages_must_match_one_shared_plan",
            "typed_stops_must_precede_every_runtime_hook",
            "TDG6_all_of_thresholds_and_retry_evidence_must_remain_unchanged",
            "synthetic_controls_must_not_load_CAL11_or_PROTO14_arrays_or_create_runs_namespace",
            "fine_only_commit_late_shadow_ownership_and_tracer_rollback_must_pass",
            "valid_retry_replan_TDG6_minimum_exhaustion_and_defensive_TDG7_stop_must_pass",
            "retry_successor_lineage_and_frozen_minimum_width_must_fail_closed",
            "disjoint_entrypoints_and_protocol_discriminator_limitation_must_be_explicit",
            "replacement_protocol_design_authorization_must_not_promote_freeze_runtime_or_calibration",
            "mutation_nonpromotion_lineage_and_hash_attacks_must_fail_closed",
            "canonical_hash_bound_result_required",
        },
    )
    _all_true("proof contract", config["proof_contract"])
    claims = config.get("claims", {})
    if (
        set(claims) != CLAIM_TRUE | CLAIM_FALSE
        or any(claims.get(key) is not True for key in CLAIM_TRUE)
        or any(claims.get(key) is not False for key in CLAIM_FALSE)
    ):
        raise ValueError("IMP3 claims differ")


def load_config(path: Path = CONFIG) -> dict[str, Any]:
    with path.open("rb") as handle:
        config = tomllib.load(handle)
    validate_config_data(config)
    return config


def _canonical_object(raw: bytes, name: str) -> dict[str, Any]:
    value = json.loads(raw)
    if not isinstance(value, dict) or canonical(value) != raw:
        raise ValueError(f"IMP3 {name} is not canonical JSON")
    return value


def _verify_lineage(config: Mapping[str, Any]) -> dict[str, Any]:
    lineage = config["immutable_lineage"]
    _git("cat-file", "-e", f"{PREDECESSOR_COMMIT}^{{commit}}")
    if _git("merge-base", "--is-ancestor", PREDECESSOR_COMMIT, "HEAD", check=False).returncode:
        raise ValueError("IMP3 sealed predecessor is not an ancestor of HEAD")
    bound: dict[str, str] = {}
    for field, relative in PREF23_PATHS.items():
        expected = lineage[f"{field}_sha256"]
        committed = _git("show", f"{PREDECESSOR_COMMIT}:{relative}").stdout
        worktree = (REPOSITORY / relative).read_bytes()
        if committed != worktree or _sha_bytes(committed) != expected:
            raise ValueError(f"IMP3 PREF23 predecessor differs: {relative}")
        bound[relative] = expected
    core: dict[str, str] = {}
    for field, relative in CORE_PATHS.items():
        expected = lineage[f"{field}_sha256"]
        observed = sha(REPOSITORY / relative)
        if observed != expected:
            raise ValueError(f"IMP3 bound core hash differs: {relative}")
        core[relative] = observed

    pref23 = _canonical_object(
        (REPOSITORY / PREF23_PATHS["PREF23_result"]).read_bytes(), "PREF23 result"
    )
    payload = pref23.get("artifact_payload", {})
    boundary = payload.get("claim_boundary", {})
    authorization = payload.get("runtime_repair_authorization", {})
    if (
        pref23.get("artifact_id") != "FGC-1-TDG7-PREF23"
        or pref23.get("gate_status")
        != "PASS_TDG7_RUNTIME_REPAIR_IMPLEMENTATION_AUTHORIZED"
        or boundary.get("CAL11_terminal_invalid_runtime_result_preserved") is not True
        or boundary.get("TDG7_independent_binder_completed") is not True
        or boundary.get("TDG7_runtime_repair_implementation_authorized") is not True
        or boundary.get("TDG7_runtime_repair_implemented") is not False
        or authorization.get("one_shared_stage_safe_plan_must_own_all_one_two_four_paths")
        is not True
        or authorization.get("historical_PROTO14_runtime_mutation_authorized")
        is not False
        or authorization.get("historical_terminal_checkpoint_may_resume") is not False
        or authorization.get("new_protocol_or_output_namespace_authorized") is not False
        or authorization.get("fresh_GR0_calibration_authorized") is not False
        or authorization.get("candidate_execution_authorized") is not False
        or authorization.get("production_trajectory_authorized") is not False
    ):
        raise ValueError("IMP3 PREF23 authorization boundary differs")
    return {
        "sealed_predecessor_commit": PREDECESSOR_COMMIT,
        "sealed_predecessor_is_ancestor_of_HEAD": True,
        "PREF23_git_and_worktree_sha256": bound,
        "current_core_sha256": core,
        "PREF23_result_canonical": True,
        "PREF23_independent_binder_completed": True,
        "PREF23_authorized_runtime_repair_only": True,
        "CAL11_invalid_runtime_status_preserved": True,
        "CAL11_history_arrays_loaded": False,
        "PROTO14_history_arrays_loaded": False,
        "campaign_checkpoint_loaded": False,
        "campaign_checkpoint_resumed": False,
        "runs_namespace_created": False,
    }


def _state() -> EvolutionState:
    radii = np.linspace(0.0, 1.0, 9)
    u = np.zeros((9, 6))
    u[:, 0] = 1.0
    u[:, 1] = 0.025
    u[:, 2] = 1.0
    u[:, 3] = radii
    q = np.zeros_like(u)
    q[:, 3] = 1.0
    return EvolutionState(u, np.zeros_like(u), q)


def _rhs(*, failed_time: float | None = None, source_only: bool = False):
    def synthetic_rhs(time: float, state: EvolutionState) -> EvolutionRHS:
        failed = failed_time is not None and (
            np.float64(time).tobytes() == np.float64(failed_time).tobytes()
        )
        return EvolutionRHS(
            state.u,
            state.p,
            state.q,
            {
                "source_residual_infinity": 2.0 if failed and source_only else 0.0,
                "source_refinement_iterations": 0,
                "source_residual_decreased_monotonically": True,
                "kinetic_condition_infinity": 1.0,
                "coordinate_speed_upper": 1.0,
                "minimum_lapse": 0.0 if failed and not source_only else 1.0,
                "minimum_radial_metric": 1.0,
                "minimum_areal_radius_away_from_center": 0.125,
            },
        )

    return synthetic_rhs


def _transaction() -> GR0RuntimeStageTransaction:
    return GR0RuntimeStageTransaction(
        thresholds=GR0UniversalThresholds(),
        causal_state=CausalBudgetState(accepted_time=START, previous_speed_upper=1.0),
        boundary_geometry=BoundaryGeometry(128.0, 24.0, 16.0, 0.375),
        grid_spacing=1.0,
        cfl_maximum=1.0,
    )


def _tracers(state: EvolutionState) -> NormalFlowTracers:
    return NormalFlowTracers.create(
        minimum=0.25,
        maximum=0.75,
        spacing=0.25,
        state=state,
        coordinates=np.linspace(0.0, 1.0, 9),
        cutoff=1.0,
        outer_radius=128.0,
    )


def _prepare(
    method: str = PRIMARY_METHOD,
    *,
    state: EvolutionState | None = None,
    transaction: GR0RuntimeStageTransaction | None = None,
    tracers: NormalFlowTracers | None = None,
    ledger: tdg6.TDG6TemporalLedger | None = None,
    rhs: Callable[[float, EvolutionState], EvolutionRHS] | None = None,
    requested_cap: float = CAP,
):
    current = _state() if state is None else state
    tx = _transaction() if transaction is None else transaction
    flow = _tracers(current) if tracers is None else tracers
    temporal = (
        tdg6.TDG6TemporalLedger.zero(initial_time=START) if ledger is None else ledger
    )
    prepared = runtime.prepare_tdg7_initial_stage_safe_runtime(
        method=method,
        time=START,
        event_target=TARGET,
        requested_cap=requested_cap,
        state=current,
        rhs=_rhs() if rhs is None else rhs,
        projector=None,
        transaction=tx,
        tracers=flow,
        coordinates=np.linspace(0.0, 1.0, 9),
        temporal_ledger=temporal,
        previous_step_index=0,
        previous_transaction_serial=0,
    )
    return prepared, current, tx, flow, temporal


def _interval(value: float, *, subintervals: int):
    coefficients = np.tile(
        np.asarray((value, 0.0, 0.0, 0.0), dtype=np.float64), (subintervals, 1)
    )
    return tdg6._magnitude_interval(
        coefficients,
        np.zeros_like(coefficients),
        subinterval_count=subintervals,
        owned_row_count=1,
    )


def _synthetic_admission(
    *, method: str = PRIMARY_METHOD, failed_channel: str | None = None
) -> tdg6.TDG6ContinuousAdmissionEvidence:
    outer = {
        channel: _interval(16.0, subintervals=2) for channel in TDG6_COMPLETE_STATE_CHANNELS
    }
    fine = {
        channel: _interval(1.0, subintervals=4) for channel in TDG6_COMPLETE_STATE_CHANNELS
    }
    if failed_channel is not None:
        outer[failed_channel] = _interval(1.0e-12, subintervals=2)
        fine[failed_channel] = _interval(9.0e-13, subintervals=4)
    return tdg6.classify_tdg6_runtime_intervals(
        method=method, outer_intervals=outer, finest_intervals=fine
    )


def _plan_hex(plan: Any) -> dict[str, Any]:
    return {
        "current_hex": plan.current.hex(),
        "target_hex": plan.event_target.hex(),
        "requested_cap_hex": plan.requested_cap.hex(),
        "quantum_hex": plan.quantum.hex(),
        "macro_width_hex": plan.macro_width.hex(),
        "fine_width_hex": plan.fine_width.hex(),
        "minimum_width_hex": plan.minimum_width.hex(),
        "boundaries_hex": [value.hex() for value in plan.boundaries],
        "one_full_stage_times_hex": [[value.hex() for value in item] for item in plan.one_full_stage_times],
        "two_half_stage_times_hex": [[value.hex() for value in item] for item in plan.two_half_stage_times],
        "four_quarter_stage_times_hex": [[value.hex() for value in item] for item in plan.four_quarter_stage_times],
    }


def _path_stage_records(prepared: Any) -> dict[str, list[list[dict[str, str]]]]:
    result: dict[str, list[list[dict[str, str]]]] = {}
    for path in (
        prepared.prepared_tdg6.outer,
        prepared.prepared_tdg6.medium,
        prepared.prepared_tdg6.fine,
    ):
        result[path.level] = [
            [
                {"stage_name": record.stage_name, "time_hex": record.time.hex()}
                for record in attempt.proposal.stages
            ]
            for attempt in path.attempts
        ]
    return result


def _method_control(method: str) -> dict[str, Any]:
    prepared, state, transaction, tracers, ledger = _prepare(method)
    inner = prepared.prepared_tdg6
    paths = (inner.outer, inner.medium, inner.fine)
    monitor_before = transaction.state
    causal_before = transaction.causal_state
    tracer_before = runtime._tracer_snapshot(tracers)
    state_before = array_content_sha256(state.u, state.p, state.q)
    if (
        transaction.state != monitor_before
        or transaction.causal_state != causal_before
        or runtime._tracer_snapshot(tracers) != tracer_before
        or array_content_sha256(state.u, state.p, state.q) != state_before
        or ledger != inner.original_temporal_ledger
    ):
        raise ValueError("IMP3 synthetic preparation mutated its accepted boundary")
    runtime.validate_tdg7_stage_safe_prepared(prepared)
    committed = runtime.commit_tdg7_stage_safe_runtime(
        prepared,
        transaction=transaction,
        tracers=tracers,
        temporal_ledger=ledger,
        current_time=START,
        current_state=state,
        current_step_index=0,
        current_transaction_serial=0,
    )
    counts = Counter(item.classification for item in inner.continuous_admission.channel_admissions)
    return {
        "method": method,
        "plan": _plan_hex(prepared.plan),
        "actual_shadow_stage_records": _path_stage_records(prepared),
        "shadow_proposal_count": sum(len(path.attempts) for path in paths),
        "outer_stage_records": inner.outer.stage_record_count,
        "medium_stage_records": inner.medium.stage_record_count,
        "fine_stage_records": inner.fine.stage_record_count,
        "total_shadow_stage_records": sum(path.stage_record_count for path in paths),
        "classification_counts": dict(sorted(counts.items())),
        "complete_admission_passed": inner.continuous_admission.complete_admission_passed,
        "failed_channels": list(inner.continuous_admission.failed_channels),
        "real_boundary_unchanged_during_prepare": True,
        "one_shared_plan_owns_all_paths": True,
        "outer_and_medium_committed": False,
        "fine_quarter_steps_committed_atomically": committed.fine_quarter_steps_committed_atomically,
        "committed_time_hex": committed.time.hex(),
        "committed_state_sha256": array_content_sha256(
            committed.state.u, committed.state.p, committed.state.q
        ),
        "fine_shadow_state_sha256": array_content_sha256(
            inner.fine.final_accepted.state.u,
            inner.fine.final_accepted.state.p,
            inner.fine.final_accepted.state.q,
        ),
        "committed_tracer_snapshot": asdict(runtime._tracer_snapshot(tracers)),
        "fine_shadow_tracer_snapshot": asdict(inner.fine.final_tracer_snapshot),
    }


def _lattice_stop_controls() -> dict[str, dict[str, Any]]:
    q = math.ulp(START)
    cases = (
        ("nonfinite_time", math.nan, TARGET, CAP, "nonfinite_coordinate"),
        ("target_not_ahead", START, START, CAP, "target_not_ahead"),
        ("outside_envelope", 1.0, TARGET, CAP, "outside_frozen_time_envelope"),
        ("endpoint_not_q_aligned", np.nextafter(2.0, -math.inf), 2.0, CAP, "endpoint_not_q_aligned"),
        ("target_not_stage_aligned", START, START + q, CAP, "target_not_stage_lattice_aligned"),
        ("no_positive_width", START, TARGET, q, "no_positive_aligned_width"),
        ("below_frozen_minimum", START, TARGET, 8.0 * q, "below_minimum_aligned_width"),
        ("zero_cap", START, TARGET, 0.0, "nonpositive_requested_cap"),
        ("negative_cap", START, TARGET, -1.0, "nonpositive_requested_cap"),
        ("nonfinite_cap_nan", START, TARGET, math.nan, "nonfinite_coordinate"),
        ("nonfinite_cap_inf", START, TARGET, math.inf, "nonfinite_coordinate"),
    )
    outcomes: dict[str, dict[str, Any]] = {}
    for case, time, target, cap, expected_reason in cases:
        hooks: list[str] = []

        class Poison:
            def __getattr__(self, name: str) -> object:
                hooks.append(name)
                raise AssertionError(f"runtime hook touched: {name}")

        try:
            runtime.prepare_tdg7_initial_stage_safe_runtime(
                method=PRIMARY_METHOD,
                time=time,
                event_target=target,
                requested_cap=cap,
                state=Poison(),
                rhs=Poison(),
                projector=Poison(),
                transaction=Poison(),
                tracers=Poison(),
                coordinates=Poison(),
                temporal_ledger=Poison(),
                previous_step_index=0,
                previous_transaction_serial=0,
            )
        except runtime.TDG7CoordinateLatticeStop as error:
            if error.reason != expected_reason or hooks:
                raise ValueError("IMP3 typed lattice stop touched a runtime hook") from error
        else:
            raise ValueError(f"IMP3 typed lattice stop did not raise: {expected_reason}")
        outcomes[case] = {
            "case": case,
            "typed_stop": "TDG7CoordinateLatticeStop",
            "reason": expected_reason,
            "hook_count": 0,
        }
    return outcomes


def _frozen_minimum_width_control() -> dict[str, Any]:
    hooks: list[str] = []

    class Poison:
        def __getattr__(self, name: str) -> object:
            hooks.append(name)
            raise AssertionError(f"runtime hook touched: {name}")

    try:
        runtime.prepare_tdg7_initial_stage_safe_runtime(
            method=PRIMARY_METHOD,
            time=START,
            event_target=TARGET,
            requested_cap=CAP,
            state=Poison(),
            rhs=Poison(),
            projector=Poison(),
            transaction=Poison(),
            tracers=Poison(),
            coordinates=Poison(),
            temporal_ledger=Poison(),
            previous_step_index=0,
            previous_transaction_serial=0,
            minimum_width=0.0,
        )
    except TypeError:
        rejected = True
    else:
        rejected = False
    if not rejected or hooks:
        raise ValueError("IMP3 frozen TDG6 minimum-width lock differs")
    return {
        "caller_minimum_width_override_rejected": True,
        "runtime_hook_count": 0,
        "frozen_minimum_width_hex": float(TDG6_MINIMUM_MACRO_STEP).hex(),
    }


def _interval_controls() -> dict[str, Any]:
    zero = tdg6.classify_tdg6_runtime_intervals(
        method=PRIMARY_METHOD,
        outer_intervals={
            channel: _interval(0.0, subintervals=2) for channel in TDG6_COMPLETE_STATE_CHANNELS
        },
        finest_intervals={
            channel: _interval(0.0, subintervals=4) for channel in TDG6_COMPLETE_STATE_CHANNELS
        },
    )
    failed_channel = TDG6_COMPLETE_STATE_CHANNELS[7]
    veto = _synthetic_admission(failed_channel=failed_channel)
    large = tdg6.classify_tdg6_runtime_intervals(
        method=PRIMARY_METHOD,
        outer_intervals={
            channel: _interval(80.0, subintervals=2) for channel in TDG6_COMPLETE_STATE_CHANNELS
        },
        finest_intervals={
            channel: _interval(10.0, subintervals=4) for channel in TDG6_COMPLETE_STATE_CHANNELS
        },
    )
    if not (
        zero.complete_admission_passed
        and {item.classification for item in zero.channel_admissions} == {"exact_zero"}
        and not veto.complete_admission_passed
        and veto.failed_channels == (failed_channel,)
        and veto.channel_admissions[7].classification == "resolved_order_failure"
        and large.complete_admission_passed
        and min(large.finest_pair_debit_vector) > 9.0
    ):
        raise ValueError("IMP3 inherited TDG6 interval controls failed")
    return {
        "channel_count": len(TDG6_COMPLETE_STATE_CHANNELS),
        "exact_zero_classification": "exact_zero",
        "exact_zero_complete_admission_passed": True,
        "one_channel_veto_channel": failed_channel,
        "one_channel_veto_classification": "resolved_order_failure",
        "one_channel_veto_complete_admission_passed": False,
        "large_convergent_complete_admission_passed": True,
        "large_convergent_minimum_retained_debit": min(large.finest_pair_debit_vector),
        "absolute_or_magnitude_tolerance_used": False,
        "complete_admission_is_all_of": True,
    }


class _FailOnceTracers(NormalFlowTracers):
    fail_next_commit: bool = False

    def commit_advance(self, positions: np.ndarray, proper_times: np.ndarray) -> None:
        if self.fail_next_commit:
            self.fail_next_commit = False
            raise RuntimeError("injected TDG7 tracer commit failure")
        super().commit_advance(positions, proper_times)


def _failure_controls() -> dict[str, Any]:
    outcomes: dict[str, Any] = {}
    for source_only, label in ((False, "non_source"), (True, "source_only")):
        state = _state()
        transaction = _transaction()
        tracers = _tracers(state)
        ledger = tdg6.TDG6TemporalLedger.zero(initial_time=START)
        before = (
            array_content_sha256(state.u, state.p, state.q),
            transaction.state,
            transaction.causal_state,
            runtime._tracer_snapshot(tracers),
            ledger,
        )
        try:
            _prepare(
                state=state,
                transaction=transaction,
                tracers=tracers,
                ledger=ledger,
                rhs=_rhs(failed_time=1.4765625, source_only=source_only),
            )
        except tdg6.TDG6RefinementPathStop as error:
            if (
                error.path != "fine_2"
                or error.source_retry_owned_by_PROTO7 is not source_only
                or (source_only and error.source_retry is None)
                or (not source_only and error.source_retry is not None)
                or error.temporal_retry
            ):
                raise ValueError("IMP3 late-shadow ownership changed") from error
        else:
            raise ValueError("IMP3 late-shadow synthetic stop did not raise")
        after = (
            array_content_sha256(state.u, state.p, state.q),
            transaction.state,
            transaction.causal_state,
            runtime._tracer_snapshot(tracers),
            ledger,
        )
        if after != before:
            raise ValueError("IMP3 late shadow stop mutated a real boundary")
        outcomes[label] = {
            "path": "fine_2",
            "source_retry_owned_by_PROTO7": source_only,
            "real_boundaries_unchanged": True,
        }

    prepared, state, transaction, tracers, ledger = _prepare()
    tracer_before = runtime._tracer_snapshot(tracers)
    tracers.positions[0] = np.nextafter(tracers.positions[0], np.inf)
    try:
        runtime.commit_tdg7_stage_safe_runtime(
            prepared,
            transaction=transaction,
            tracers=tracers,
            temporal_ledger=ledger,
            current_time=START,
            current_state=state,
            current_step_index=0,
            current_transaction_serial=0,
        )
    except RuntimeError:
        drift_rejected = True
    else:
        drift_rejected = False
    if not drift_rejected or runtime._tracer_snapshot(tracers) == tracer_before:
        # A rejected drift must remain an observed drift; no hidden restoration is claimed.
        raise ValueError("IMP3 tracer precommit drift control differs")

    state = _state()
    base = _tracers(state)
    tracers = _FailOnceTracers(
        base.labels,
        base.positions,
        base.proper_times,
        base.event_proper_times,
        base.event_fields,
        base.cutoff,
        base.outer_radius,
    )
    transaction = _transaction()
    ledger = tdg6.TDG6TemporalLedger.zero(initial_time=START)
    prepared, _, _, _, _ = _prepare(
        state=state, transaction=transaction, tracers=tracers, ledger=ledger
    )
    monitor_before = transaction.state
    causal_before = transaction.causal_state
    tracer_before = runtime._tracer_snapshot(tracers)
    tracers.fail_next_commit = True
    try:
        runtime.commit_tdg7_stage_safe_runtime(
            prepared,
            transaction=transaction,
            tracers=tracers,
            temporal_ledger=ledger,
            current_time=START,
            current_state=state,
            current_step_index=0,
            current_transaction_serial=0,
        )
    except RuntimeError as error:
        rollback = (
            str(error) == "injected TDG7 tracer commit failure"
            and transaction.state == monitor_before
            and transaction.causal_state == causal_before
            and runtime._tracer_snapshot(tracers) == tracer_before
        )
    else:
        rollback = False
    if not rollback:
        raise ValueError("IMP3 tracer rollback control failed")
    outcomes["precommit_tracer_drift_rejected"] = drift_rejected
    outcomes["late_tracer_commit_rollback_passed"] = rollback
    return outcomes


def _retry_controls() -> dict[str, Any]:
    prepared, state, transaction, tracers, ledger = _prepare()
    failed_channel = TDG6_COMPLETE_STATE_CHANNELS[0]
    failed_inner = replace(
        prepared.prepared_tdg6,
        continuous_admission=_synthetic_admission(failed_channel=failed_channel),
    )
    failed = runtime.TDG7PreparedStageSafeRuntime(prepared.plan, failed_inner)
    sink: list[dict[str, object]] = []
    try:
        runtime.require_tdg7_stage_safe_admission(
            failed,
            transaction=transaction,
            tracers=tracers,
            temporal_ledger=ledger,
            current_time=START,
            current_state=state,
            current_step_index=0,
            current_transaction_serial=0,
            durable_rejection_sink=lambda item: sink.append(dict(item)),
        )
    except tdg6.TDG6TemporalRetryRequired as error:
        retry = error
    else:
        raise ValueError("IMP3 failed TDG6 admission did not issue a retry")
    successor = runtime.replan_tdg7_after_tdg6_retry(failed, retry)
    try:
        runtime.prepare_tdg7_initial_stage_safe_runtime(
            method=PRIMARY_METHOD,
            time=successor.plan.current,
            event_target=successor.plan.event_target,
            requested_cap=successor.plan.requested_cap,
            state=state,
            rhs=_rhs(),
            projector=None,
            transaction=transaction,
            tracers=tracers,
            coordinates=np.linspace(0.0, 1.0, 9),
            temporal_ledger=successor.temporal_ledger,
            previous_step_index=0,
            previous_transaction_serial=0,
        )
    except ValueError as error:
        active_retry_ledger_rejected_by_fresh_entrypoint = (
            "active retry ledger" in str(error)
        )
    else:
        active_retry_ledger_rejected_by_fresh_entrypoint = False
    resumed = runtime.prepare_tdg7_retry_successor(
        successor=successor,
        method=PRIMARY_METHOD,
        state=state,
        rhs=_rhs(),
        projector=None,
        transaction=transaction,
        tracers=tracers,
        coordinates=np.linspace(0.0, 1.0, 9),
        previous_step_index=0,
        previous_transaction_serial=0,
    )
    second_wrapped = runtime.TDG7PreparedStageSafeRuntime(
        resumed.plan,
        replace(
            resumed.prepared_tdg6,
            continuous_admission=failed.prepared_tdg6.continuous_admission,
        ),
    )
    try:
        runtime.require_tdg7_stage_safe_admission(
            second_wrapped,
            transaction=transaction,
            tracers=tracers,
            temporal_ledger=successor.temporal_ledger,
            current_time=START,
            current_state=state,
            current_step_index=0,
            current_transaction_serial=0,
            durable_rejection_sink=lambda _item: None,
        )
    except tdg6.TDG6TemporalRetryRequired as error:
        second_retry = error
    else:
        raise ValueError("IMP3 second TDG6 retry did not issue")
    prefix_retained = (
        second_retry.updated_ledger.serialized_temporal_rejections[:-1]
        == successor.temporal_ledger.serialized_temporal_rejections
    )
    last_accepted_retry_retained = (
        second_retry.updated_ledger.last_accepted_macro_step_temporal_retry_count
        == successor.temporal_ledger.last_accepted_macro_step_temporal_retry_count
    )
    latest_rejection = second_retry.updated_ledger.serialized_temporal_rejections[-1]
    forged_prefix = replace(
        second_retry.updated_ledger,
        serialized_temporal_rejections=(latest_rejection, latest_rejection),
    )
    try:
        runtime.replan_tdg7_after_tdg6_retry(
            second_wrapped,
            tdg6.TDG6TemporalRetryRequired(second_retry.evidence, forged_prefix),
        )
    except ValueError as error:
        forged_prefix_rejected = "durable TDG6 evidence" in str(error)
    else:
        forged_prefix_rejected = False

    minimum_state = _state()
    minimum_transaction = _transaction()
    minimum_tracers = _tracers(minimum_state)
    minimum_ledger = tdg6.TDG6TemporalLedger.zero(initial_time=START)
    twice_minimum = 2.0 * float(TDG6_MINIMUM_MACRO_STEP)
    minimum_prepared = runtime.prepare_tdg7_initial_stage_safe_runtime(
        method=PRIMARY_METHOD,
        time=START,
        event_target=TARGET,
        requested_cap=twice_minimum,
        state=minimum_state,
        rhs=_rhs(),
        projector=None,
        transaction=minimum_transaction,
        tracers=minimum_tracers,
        coordinates=np.linspace(0.0, 1.0, 9),
        temporal_ledger=minimum_ledger,
        previous_step_index=0,
        previous_transaction_serial=0,
    )
    minimum_failed = runtime.TDG7PreparedStageSafeRuntime(
        minimum_prepared.plan,
        replace(
            minimum_prepared.prepared_tdg6,
            continuous_admission=failed.prepared_tdg6.continuous_admission,
        ),
    )
    try:
        runtime.require_tdg7_stage_safe_admission(
            minimum_failed,
            transaction=minimum_transaction,
            tracers=minimum_tracers,
            temporal_ledger=minimum_ledger,
            current_time=START,
            current_state=minimum_state,
            current_step_index=0,
            current_transaction_serial=0,
            durable_rejection_sink=lambda _item: None,
        )
    except tdg6.TDG6TemporalRetryRequired as error:
        minimum_retry = error
    else:
        raise ValueError("IMP3 two-minimum TDG6 retry did not issue")
    minimum_successor = runtime.replan_tdg7_after_tdg6_retry(
        minimum_failed, minimum_retry
    )
    minimum_resumed = runtime.prepare_tdg7_retry_successor(
        successor=minimum_successor,
        method=PRIMARY_METHOD,
        state=minimum_state,
        rhs=_rhs(),
        projector=None,
        transaction=minimum_transaction,
        tracers=minimum_tracers,
        coordinates=np.linspace(0.0, 1.0, 9),
        previous_step_index=0,
        previous_transaction_serial=0,
    )
    minimum_wrapped = runtime.TDG7PreparedStageSafeRuntime(
        minimum_resumed.plan,
        replace(
            minimum_resumed.prepared_tdg6,
            continuous_admission=failed.prepared_tdg6.continuous_admission,
        ),
    )
    try:
        runtime.require_tdg7_stage_safe_admission(
            minimum_wrapped,
            transaction=minimum_transaction,
            tracers=minimum_tracers,
            temporal_ledger=minimum_successor.temporal_ledger,
            current_time=START,
            current_state=minimum_state,
            current_step_index=0,
            current_transaction_serial=0,
            durable_rejection_sink=lambda _item: None,
        )
    except tdg6.TDG6TemporalRetryExhausted as error:
        minimum_stop = error
    else:
        raise ValueError("IMP3 minimum TDG6 retry did not exhaust")
    impossible_required = tdg6.TDG6TemporalRetryRequired(
        minimum_stop.evidence, minimum_stop.updated_ledger
    )
    try:
        runtime.replan_tdg7_after_tdg6_retry(minimum_wrapped, impossible_required)
    except runtime.TDG7RetryLatticeExhausted as error:
        exhausted = error
    else:
        raise ValueError("IMP3 forged minimum retry did not reach TDG7 defense")
    over_limit = replace(
        retry.evidence,
        retry_count_for_current_macro_step=33,
        cumulative_temporal_retry_count=33,
    )
    try:
        runtime.replan_tdg7_after_tdg6_retry(
            failed, tdg6.TDG6TemporalRetryRequired(over_limit, retry.updated_ledger)
        )
    except ValueError as error:
        over_limit_rejected = "exceeds TDG6 retry limit" in str(error)
    else:
        over_limit_rejected = False
    mismatched_method = replace(retry.evidence, method=COMPARATOR_METHOD)
    try:
        runtime.replan_tdg7_after_tdg6_retry(
            failed,
            tdg6.TDG6TemporalRetryRequired(mismatched_method, retry.updated_ledger),
        )
    except ValueError as error:
        mismatched_method_rejected = "differs" in str(error)
    else:
        mismatched_method_rejected = False
    forged_mapping = json.loads(retry.updated_ledger.serialized_temporal_rejections[-1])
    forged_mapping["failed_channels"] = []
    forged_ledger = replace(
        retry.updated_ledger,
        serialized_temporal_rejections=(
            json.dumps(forged_mapping, sort_keys=True, separators=(",", ":")),
        ),
    )
    try:
        runtime.replan_tdg7_after_tdg6_retry(
            failed, tdg6.TDG6TemporalRetryRequired(retry.evidence, forged_ledger)
        )
    except ValueError as error:
        forged_evidence_rejected = "durable TDG6 evidence" in str(error)
    else:
        forged_evidence_rejected = False
    if not (
        len(sink) == 1
        and len(sink[0].get("channel_admissions", [])) == 18
        and successor.predecessor == failed
        and successor.temporal_ledger == retry.updated_ledger
        and successor.plan.macro_width < failed.plan.macro_width
        and active_retry_ledger_rejected_by_fresh_entrypoint
        and resumed.plan == successor.plan
        and resumed.prepared_tdg6.original_temporal_ledger == successor.temporal_ledger
        and prefix_retained
        and last_accepted_retry_retained
        and forged_prefix_rejected
        and minimum_prepared.plan.macro_width == twice_minimum
        and minimum_successor.plan.macro_width == float(TDG6_MINIMUM_MACRO_STEP)
        and minimum_stop.reason == "minimum_macro_step"
        and minimum_stop.physical_classification is False
        and len(minimum_stop.updated_ledger.serialized_temporal_rejections) == 2
        and exhausted.evidence == minimum_stop.evidence
        and exhausted.updated_ledger == minimum_stop.updated_ledger
        and exhausted.requested_half_cap == float(TDG6_MINIMUM_MACRO_STEP) / 2.0
        and exhausted.lattice_reason == "below_minimum_aligned_width"
        and exhausted.physical_classification is False
        and over_limit_rejected
        and mismatched_method_rejected
        and forged_evidence_rejected
    ):
        raise ValueError("IMP3 retry/replan evidence differs")
    return {
        "first_retry_channel_decision_count": len(sink[0]["channel_admissions"]),
        "retry_is_physical_classification": False,
        "original_macro_width_hex": failed.plan.macro_width.hex(),
        "replanned_macro_width_hex": successor.plan.macro_width.hex(),
        "replan_strictly_reduced": True,
        "retry_successor_binds_predecessor_plan_and_exact_updated_ledger": True,
        "canonical_durable_TDG6_evidence_retained": True,
        "active_retry_ledger_rejected_by_fresh_entrypoint": (
            active_retry_ledger_rejected_by_fresh_entrypoint
        ),
        "retry_successor_resume_uses_exact_updated_ledger": True,
        "full_ledger_prefix_retained": prefix_retained,
        "last_accepted_retry_count_retained": last_accepted_retry_retained,
        "forged_full_ledger_prefix_rejected": forged_prefix_rejected,
        "two_minimum_to_minimum_successor": True,
        "TDG6_minimum_macro_step_exhaustion_reason": minimum_stop.reason,
        "TDG6_minimum_macro_step_exhaustion_is_physical": False,
        "TDG7_lattice_exhaustion_is_defensive_forged_wrapper_only": True,
        "over_limit_retry_lineage_rejected": over_limit_rejected,
        "mismatched_method_retry_lineage_rejected": mismatched_method_rejected,
        "forged_canonical_retry_evidence_rejected": forged_evidence_rejected,
        "defensive_lattice_stop": "TDG7RetryLatticeExhausted",
        "defensive_lattice_reason": exhausted.lattice_reason,
        "defensive_lattice_stop_retains_TDG6_evidence": True,
        "defensive_lattice_stop_is_physical_classification": False,
    }


def _synthetic_certificate() -> dict[str, Any]:
    methods = {
        "RK4": _method_control(PRIMARY_METHOD),
        "SSPRK3": _method_control(COMPARATOR_METHOD),
    }
    stops = _lattice_stop_controls()
    frozen_minimum = _frozen_minimum_width_control()
    intervals = _interval_controls()
    failures = _failure_controls()
    retries = _retry_controls()
    protocol_discriminator = {
        "fresh_vs_retry_protocol_discriminator_implemented": (
            runtime.FRESH_VS_RETRY_PROTOCOL_DISCRIMINATOR_IMPLEMENTED
        ),
        "stateless_adapter_cannot_detect_discarded_external_retry": (
            runtime.STATELESS_ADAPTER_CANNOT_DETECT_DISCARDED_EXTERNAL_RETRY
        ),
        "future_protocol_cursor_required_before_calibration": True,
    }
    if not (
        methods["RK4"]["total_shadow_stage_records"] == 35
        and methods["RK4"]["fine_stage_records"] == 20
        and methods["SSPRK3"]["total_shadow_stage_records"] == 28
        and methods["SSPRK3"]["fine_stage_records"] == 16
        and len(stops) == 11
        and all(item["hook_count"] == 0 for item in stops.values())
        and frozen_minimum["caller_minimum_width_override_rejected"] is True
        and protocol_discriminator
        == {
            "fresh_vs_retry_protocol_discriminator_implemented": False,
            "stateless_adapter_cannot_detect_discarded_external_retry": True,
            "future_protocol_cursor_required_before_calibration": True,
        }
        and all(
            item["shadow_proposal_count"] == 7
            and item["complete_admission_passed"] is True
            and item["one_shared_plan_owns_all_paths"] is True
            and item["outer_and_medium_committed"] is False
            and item["fine_quarter_steps_committed_atomically"] is True
            and item["committed_state_sha256"] == item["fine_shadow_state_sha256"]
            and item["committed_tracer_snapshot"] == item["fine_shadow_tracer_snapshot"]
            for item in methods.values()
        )
    ):
        raise ValueError("IMP3 synthetic method controls failed")
    return {
        "runtime_environment": {
            "python_implementation": sys.implementation.name,
            "python_major_minor": f"{sys.version_info.major}.{sys.version_info.minor}",
            "numpy": np.__version__,
            "binary64_epsilon": float(np.finfo(np.float64).eps),
        },
        "synthetic_states_only": True,
        "CAL11_history_arrays_loaded": False,
        "PROTO14_history_arrays_loaded": False,
        "runs_namespace_created": False,
        "method_controls": methods,
        "typed_lattice_stop_controls": stops,
        "frozen_minimum_width_control": frozen_minimum,
        "protocol_discriminator_limitation": protocol_discriminator,
        "TDG6_all_of_interval_controls": intervals,
        "late_shadow_and_tracer_controls": failures,
        "retry_replan_controls": retries,
        "all_synthetic_controls_pass": True,
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
    complete = _synthetic_admission()
    outer = {item.channel: item.outer_difference for item in complete.channel_admissions}
    fine = {item.channel: item.finest_difference for item in complete.channel_admissions}
    try:
        tdg6.classify_tdg6_runtime_intervals(
            method=PRIMARY_METHOD,
            outer_intervals=dict(reversed(tuple(outer.items()))),
            finest_intervals=fine,
        )
    except ValueError:
        reordered_channels = True
    else:
        reordered_channels = False
    controls = {
        "reordered_TDG6_channel_map_rejected": reordered_channels,
        "sealed_predecessor_commit_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["immutable_lineage"].__setitem__(
                "predecessor_commit", "0" * 40
            ),
        ),
        "PREF23_result_hash_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["immutable_lineage"].__setitem__(
                "PREF23_result_sha256", "0" * 64
            ),
        ),
        "core_runtime_hash_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["immutable_lineage"].__setitem__(
                "core_runtime_sha256", "0" * 64
            ),
        ),
        "shared_plan_contract_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["stage_safe_runtime_contract"].__setitem__(
                "one_immutable_shared_TDG7_plan_owns_one_two_four_paths", False
            ),
        ),
        "TDG6_all_of_contract_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["tdg6_inheritance"].__setitem__(
                "one_channel_resolved_order_failure_vetoes_all_of_admission", False
            ),
        ),
        "PROTO14_mutation_promotion_rejected": _expect_config_rejection(
            config,
            lambda value: value["claims"].__setitem__(
                "PROTO14_runtime_mutation_authorized", True
            ),
        ),
        "replacement_protocol_freeze_promotion_rejected": _expect_config_rejection(
            config,
            lambda value: value["claims"].__setitem__(
                "fresh_replacement_protocol_frozen", True
            ),
        ),
        "new_protocol_promotion_rejected": _expect_config_rejection(
            config,
            lambda value: value["claims"].__setitem__(
                "new_protocol_authorized", True
            ),
        ),
        "namespace_promotion_rejected": _expect_config_rejection(
            config,
            lambda value: value["claims"].__setitem__(
                "new_output_namespace_authorized", True
            ),
        ),
        "calibration_promotion_rejected": _expect_config_rejection(
            config,
            lambda value: value["claims"].__setitem__(
                "fresh_GR0_calibration_authorized", True
            ),
        ),
        "FGCQR_promotion_rejected": _expect_config_rejection(
            config,
            lambda value: value["claims"].__setitem__(
                "FGCQR_holdout_execution_authorized", True
            ),
        ),
        "physical_claim_promotion_rejected": _expect_config_rejection(
            config,
            lambda value: value["claims"].__setitem__(
                "physical_transition_claim_authorized", True
            ),
        ),
    }
    if set(controls.values()) != {True}:
        raise ValueError(f"IMP3 mutation controls failed: {controls}")
    return controls


def build(config: Mapping[str, Any]) -> dict[str, Any]:
    validate_config_data(config)
    lineage = _verify_lineage(config)
    synthetic = _synthetic_certificate()
    mutations = _mutation_controls(config)
    return {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "metric_signature": "-+++",
        "riemann_convention": "plus_partial_mu_gamma_nu",
        "classification": (
            "authorized_TDG7_stage_safe_runtime_repair_implemented_and_"
            "synthetically_qualified_without_history_execution"
        ),
        "generated_by": _rel(Path(__file__)),
        "gate_status": dict(config["claims"]),
        "nonclaims": {key: False for key in sorted(CLAIM_FALSE)},
        "source_config_sha256": {_rel(CONFIG): sha(CONFIG)},
        "implementation_sha256": {
            _rel(CONFIG): sha(CONFIG),
            _rel(DOCUMENT): sha(DOCUMENT),
            _rel(Path(__file__)): sha(Path(__file__)),
            **{relative: sha(REPOSITORY / relative) for relative in CORE_PATHS.values()},
        },
        "derivation_document": _rel(DOCUMENT),
        "derivation_document_sha256": sha(DOCUMENT),
        "artifact_payload": {
            "immutable_lineage": lineage,
            "stage_safe_runtime_contract": dict(config["stage_safe_runtime_contract"]),
            "tdg6_inheritance": dict(config["tdg6_inheritance"]),
            "synthetic_qualification": synthetic,
            "limitations": dict(config["limitations"]),
            "proof_contract": dict(config["proof_contract"]),
            "mutation_controls": mutations,
            "successor_boundary": dict(config["successor_boundary"]),
            "claim_boundary": dict(config["claims"]),
        },
    }


def _validate_stored_record(value: Mapping[str, Any], config: Mapping[str, Any]) -> None:
    payload = value.get("artifact_payload", {})
    synthetic = payload.get("synthetic_qualification", {})
    if (
        value.get("schema_version") != 1
        or value.get("artifact_id") != ARTIFACT_ID
        or value.get("project_version") != PROJECT_VERSION
        or value.get("classification")
        != "authorized_TDG7_stage_safe_runtime_repair_implemented_and_synthetically_qualified_without_history_execution"
        or value.get("generated_by") != "scripts/reproduce_fgc_tdg7_imp3.py"
        or value.get("gate_status") != config["claims"]
        or value.get("nonclaims") != {key: False for key in sorted(CLAIM_FALSE)}
        or value.get("source_config_sha256") != {_rel(CONFIG): sha(CONFIG)}
        or payload.get("stage_safe_runtime_contract")
        != config["stage_safe_runtime_contract"]
        or payload.get("tdg6_inheritance") != config["tdg6_inheritance"]
        or payload.get("limitations") != config["limitations"]
        or payload.get("proof_contract") != config["proof_contract"]
        or payload.get("successor_boundary") != config["successor_boundary"]
        or payload.get("claim_boundary") != config["claims"]
        or synthetic.get("all_synthetic_controls_pass") is not True
        or synthetic.get("synthetic_states_only") is not True
        or synthetic.get("CAL11_history_arrays_loaded") is not False
        or synthetic.get("PROTO14_history_arrays_loaded") is not False
        or synthetic.get("runs_namespace_created") is not False
    ):
        raise ValueError("IMP3 stored result boundary differs")
    implementation = value.get("implementation_sha256", {})
    expected_implementation_paths = {
        _rel(CONFIG),
        _rel(DOCUMENT),
        _rel(Path(__file__)),
        *CORE_PATHS.values(),
    }
    if (
        not isinstance(implementation, Mapping)
        or set(implementation) != expected_implementation_paths
        or any(sha(REPOSITORY / relative) != expected for relative, expected in implementation.items())
    ):
        raise ValueError("IMP3 implementation hash ledger differs")
    if (
        value.get("derivation_document") != "docs/fgc-tdg7-imp3.md"
        or value.get("derivation_document_sha256") != sha(DOCUMENT)
        or not isinstance(payload.get("mutation_controls"), Mapping)
        or set(payload["mutation_controls"].values()) != {True}
    ):
        raise ValueError("IMP3 stored documentation or mutations differ")


def verify_canonical(
    config_path: Path = CONFIG, output_path: Path = OUTPUT
) -> dict[str, Any]:
    config = load_config(config_path)
    raw = output_path.read_bytes()
    stored = json.loads(raw)
    if not isinstance(stored, dict) or canonical(stored) != raw:
        raise ValueError("IMP3 stored result is not canonical JSON")
    _validate_stored_record(stored, config)
    observed = build(config)
    if stored != observed:
        raise ValueError("IMP3 stored result differs from canonical reproduction")
    return stored


def _atomic_write(path: Path, payload: bytes) -> None:
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


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--verify", "--check", action="store_true")
    args = parser.parse_args()
    if args.write and args.verify:
        raise SystemExit("choose either --write or --verify")
    if args.write:
        _atomic_write(args.output, canonical(build(load_config())))
        print(f"wrote {_rel(args.output)}")
    elif args.verify:
        value = verify_canonical(output_path=args.output)
        print(
            f"PASS {ARTIFACT_ID}: implemented="
            f"{value['gate_status']['TDG7_runtime_repair_implemented']} synthetic="
            f"{value['gate_status']['TDG7_runtime_repair_synthetic_qualification_passed']} "
            f"PROTO14_mutation={value['gate_status']['PROTO14_runtime_mutation_authorized']}"
        )
    else:
        print(canonical(build(load_config())).decode("utf-8"), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
