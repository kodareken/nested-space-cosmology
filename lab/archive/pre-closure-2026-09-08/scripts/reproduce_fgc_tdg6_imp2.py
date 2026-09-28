#!/usr/bin/env python3
"""Bind the TDG6 production compositor and synthetic qualification."""

from __future__ import annotations

import argparse
from collections import Counter
from copy import deepcopy
from dataclasses import asdict, replace
from hashlib import sha256
import json
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
    TDG6_MAXIMUM_TEMPORAL_RETRIES_PER_STEP,
    TDG6_MINIMUM_MACRO_STEP,
)
from recursive_horizons.fgc.evolution import (  # noqa: E402
    tdg6_temporal_admission_runtime as runtime,
)


ARTIFACT_ID = "FGC-1-TDG6-IMP2"
PROJECT_VERSION = "0.11.0"
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-tdg6-imp2.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-tdg6-imp2.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-tdg6-imp2.md"
PREF21_CONFIG = REPOSITORY / "configs/fgc/fgc-1-tdg6-pref21.toml"
PREF21_RESULT = REPOSITORY / "results/fgc-1-tdg6-pref21.json"
PREF21_REPRODUCER = REPOSITORY / "scripts/reproduce_fgc_tdg6_pref21.py"
PREF21_MODULE = (
    REPOSITORY
    / "src/recursive_horizons/fgc/evolution/tdg6_temporal_admission_theorem.py"
)
PREF21_DOCUMENT = REPOSITORY / "docs/fgc-tdg6-pref21.md"
TDG6_DESIGN = (
    REPOSITORY
    / "src/recursive_horizons/fgc/evolution/tdg6_temporal_admission_design.py"
)
TDG5_RUNTIME = (
    REPOSITORY
    / "src/recursive_horizons/fgc/evolution/tdg5_stage_complete_refinement_runtime.py"
)
PROTO7_RUNTIME = REPOSITORY / "src/recursive_horizons/fgc/evolution/proto7_runtime.py"
TRACER_OWNER = REPOSITORY / "scripts/run_fgc_gr0_calibration.py"
RUNTIME_MODULE = (
    REPOSITORY
    / "src/recursive_horizons/fgc/evolution/tdg6_temporal_admission_runtime.py"
)
RUNTIME_TEST = REPOSITORY / "tests/test_fgc_tdg6_temporal_admission_runtime.py"
IMPLEMENTATION = (Path(__file__).resolve(), RUNTIME_MODULE, RUNTIME_TEST)


EXPECTED_CLAIMS = {
    "TDG4_sampling_identifiability_theorem_completed": True,
    "current_temporal_spectral_rule_retired_as_continuum_admission": True,
    "current_temporal_spectral_rule_retained_as_sampled_diagnostic": True,
    "TDG5_runtime_refinement_pair_implemented": True,
    "TDG6_threshold_and_admission_design_frozen": True,
    "runtime_thresholds_frozen": True,
    "replacement_temporal_admission_defined": True,
    "TDG6_independent_binder_completed": True,
    "production_compositor_implementation_authorized": True,
    "production_compositor_implemented": True,
    "production_compositor_synthetic_qualification_passed": True,
    "PROTO14_freeze_design_authorized": True,
    "production_trajectory_authorized": False,
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
        "target": "FGC-2-SF1-TDG6",
        "role": (
            "production_shaped_three_level_temporal_compositor_implementation_"
            "and_synthetic_qualification"
        ),
        "calibration_branch": "GR-0",
        "PREF21_compact_result_read": True,
        "actual_PROTO13_history_arrays_consumed": False,
        "campaign_checkpoint_loaded": False,
        "campaign_checkpoint_resumed": False,
        "production_state_advanced": False,
        "synthetic_runtime_controls_only": True,
        "SGBL_trajectory_read": False,
        "FGCQR_trajectory_read": False,
        "physical_or_candidate_question_answered": False,
    }


def _expected_lineage() -> dict[str, Any]:
    return {
        "PREF21_commit": "b42f9548715ece3cd6ed5751c4f0b9aa8ead3ce3",
        "PREF21_config_sha256": (
            "2138c4064126b7da9910e86d5bdf466f7a90391030e1d25b97779258b006d313"
        ),
        "PREF21_result_sha256": (
            "88adb24e3b2fe32d29052defe9eed6481c7a20462386ef70a8fc5cd633aa1c1b"
        ),
        "PREF21_reproducer_sha256": (
            "f7a37603dd78988b776019a1740e16e8e034498f37687b5731bd328c84196988"
        ),
        "PREF21_module_sha256": (
            "629c7b832abd5967ee522074b4f31496be97d92b8a13a43447c3fc4d171a84d5"
        ),
        "PREF21_document_sha256": (
            "d2a5d97b11ddbc7939dc81210d8fd44fb3f7d40ed87412fea7147b7cada40ef5"
        ),
        "TDG6_design_module_sha256": (
            "87d9a90e272a89479690c0d432f2d5cd979d88841615aaaeb9b49a2c0cca651e"
        ),
        "TDG5_runtime_sha256": (
            "75fc11943f459ca7549f75a812fe961e55ae04bbdb13450c1fdd509b5b53a080"
        ),
        "PROTO7_runtime_sha256": (
            "ffa75767c1fdbae42a09c8ef50276183861fc9255ad0061e5a1ee3673f142734"
        ),
        "normal_flow_tracer_owner_sha256": (
            "4e3734495edc0a289e5c61dc68e193cbe19629582234dd04c0c20a46d52774b3"
        ),
        "PREF21_commit_must_be_ancestor_of_HEAD": True,
        "all_bound_predecessor_blobs_must_match_commit_and_worktree": True,
        "PREF21_result_must_be_canonical_and_authorize_only_implementation": True,
        "raw_campaign_bundle_must_not_be_loaded": True,
    }


def _expected_compositor() -> dict[str, Any]:
    return {
        "level_names": ["outer", "medium", "fine"],
        "level_step_counts": [1, 2, 4],
        "all_levels_start_from_one_bitwise_accepted_state": True,
        "all_levels_use_same_method_source_projector_spatial_grid_and_macro_interval": True,
        "outer_and_medium_are_evidence_only": True,
        "only_complete_fine_four_quarter_path_is_committable": True,
        "all_seven_proposals_use_unchanged_PROTO7_transaction_gates": True,
        "RK4_total_shadow_stage_records": 35,
        "RK4_committed_fine_stage_records": 20,
        "SSPRK3_total_shadow_stage_records": 28,
        "SSPRK3_committed_fine_stage_records": 16,
        "production_trajectory_executed": False,
    }


def _expected_interval() -> dict[str, Any]:
    return {
        "channel_order": list(TDG6_COMPLETE_STATE_CHANNELS),
        "channel_count": 18,
        "D01_subinterval_count": 2,
        "D12_subinterval_count": 4,
        "endpoint_and_all_numerically_real_stationary_candidates_retained": True,
        "Bernstein_convex_hull_certification_retained": True,
        "coefficient_and_evaluation_debits_retained": True,
        "lower_bound_subtracts_complete_coefficient_and_evaluation_debits_outward": True,
        "upper_bound_is_complete_certified_continuous_upper": True,
        "classification_uses_exact_Fraction_images_of_binary64_bounds": True,
        "exact_squared_admission": "8*U12^2<=L01^2",
        "complete_admission_is_all_of": True,
        "per_channel_debit_is_full_D12_upper": True,
        "debit_divided_by_Richardson_denominator": False,
        "absolute_or_outcome_scaled_tolerance_used": False,
    }


def _expected_transaction() -> dict[str, Any]:
    return {
        "outer_medium_and_fine_have_independent_monitor_ledgers": True,
        "outer_medium_and_fine_have_independent_causal_ledgers": True,
        "outer_medium_and_fine_have_independent_normal_flow_tracers": True,
        "tracer_preview_runs_before_every_shadow_acceptance": True,
        "shadow_tracer_commit_mutates_only_its_level_clone": True,
        "prepare_revalidates_real_state_monitor_causal_tracer_and_debit_boundaries": True,
        "commit_revalidates_real_time_state_step_serial_monitor_causal_tracer_and_debit_boundaries": True,
        "late_commit_exception_rolls_back_monitor_causal_and_complete_tracer_state": True,
        "fine_state_monitor_causal_tracer_and_debit_are_one_return_boundary": True,
        "source_only_shadow_retry_remains_owned_by_PROTO7": True,
        "non_source_shadow_failure_remains_inherited_terminal_evidence": True,
    }


def _expected_retry() -> dict[str, Any]:
    return {
        "temporal_retry_factor": "1/2",
        "maximum_temporal_retries_per_macro_step": 32,
        "minimum_macro_step": "1/1073741824",
        "temporal_retry_counter_is_separate_from_source_and_CFL_counters": True,
        "rejection_evidence_contains_all_18_channel_decisions": True,
        "durable_sink_must_return_before_retry_ledger_exists": True,
        "failed_durable_sink_leaves_retry_ledger_unchanged": True,
        "accepted_step_resets_only_current_retry_counter": True,
        "cumulative_retry_count_and_complete_rejection_history_are_retained": True,
        "debit_accumulation_is_componentwise_nonnegative_and_outward": True,
        "checkpoint_round_trips_debit_vector_retry_counts_and_complete_rejection_history": True,
        "checkpoint_debit_vector_is_content_hashed": True,
        "checkpoint_claim_promotion_fails_closed": True,
    }


def _expected_synthetic() -> dict[str, Any]:
    return {
        "both_methods_must_prepare_all_seven_paths": True,
        "both_methods_must_commit_only_the_fine_path": True,
        "exact_zero_and_resolved_pass_classes_required": True,
        "one_tiny_nonconvergent_channel_must_veto_all_of_admission": True,
        "large_convergent_difference_must_pass_and_retain_large_debit": True,
        "late_fine_non_source_failure_must_leave_every_real_ledger_unchanged": True,
        "source_only_fine_failure_must_retain_PROTO7_ownership": True,
        "tracer_or_debit_boundary_drift_must_fail_before_commit": True,
        "tracer_commit_failure_must_roll_back_every_mutable_real_ledger": True,
        "durable_temporal_retry_and_exhaustion_controls_required": True,
        "checkpoint_round_trip_and_mutation_controls_required": True,
        "dense_continuous_polynomial_audit_required": True,
        "implementation_and_claim_mutations_must_fail_closed": True,
    }


def _expected_limitations() -> dict[str, Any]:
    return {
        "piecewise_cubic_extension_is_not_exact_PDE_history": True,
        "three_level_pass_does_not_prove_trajectory_asymptotics": True,
        "local_accumulated_debit_is_not_global_PDE_error_bound": True,
        "future_constraint_to_Raychaudhuri_stability_map_required": True,
        "spatial_ladders_and_cross_method_agreement_remain_required": True,
        "synthetic_runtime_qualification_is_not_production_trajectory": True,
        "historical_CAL10_and_PROTO13_results_are_not_reclassified": True,
        "no_GR0_SGBL_or_FGCQR_outcome_is_classified": True,
    }


def _expected_successor() -> dict[str, Any]:
    return {
        key: EXPECTED_CLAIMS[key]
        for key in (
            "TDG6_independent_binder_completed",
            "production_compositor_implementation_authorized",
            "production_compositor_implemented",
            "production_compositor_synthetic_qualification_passed",
            "PROTO14_freeze_design_authorized",
            "production_trajectory_authorized",
            "PROTO14_frozen",
            "fresh_GR0_dynamic_calibration_completed",
            "GR0_case_eligible",
            "SGBL_execution_authorized",
            "FGCQR_holdout_execution_authorized",
        )
    } | {"DEF1_execution_authorized": False}


def validate_config_data(config: Mapping[str, Any]) -> None:
    _strict_keys(
        "IMP2 config",
        config,
        {
            "schema_version",
            "artifact_id",
            "project_version",
            "metric_signature",
            "riemann_convention",
            "PREF21_config",
            "PREF21_result",
            "PREF21_reproducer",
            "PREF21_module",
            "PREF21_document",
            "TDG6_design_module",
            "TDG5_runtime",
            "PROTO7_runtime",
            "normal_flow_tracer_owner",
            "implementation_module",
            "scope",
            "immutable_lineage",
            "three_level_compositor",
            "continuous_interval",
            "transaction_and_tracer",
            "temporal_retry_and_checkpoint",
            "synthetic_qualification",
            "limitations",
            "proof_contract",
            "successor_boundary",
            "claims",
        },
    )
    paths = {
        "PREF21_config": _rel(PREF21_CONFIG),
        "PREF21_result": _rel(PREF21_RESULT),
        "PREF21_reproducer": _rel(PREF21_REPRODUCER),
        "PREF21_module": _rel(PREF21_MODULE),
        "PREF21_document": _rel(PREF21_DOCUMENT),
        "TDG6_design_module": _rel(TDG6_DESIGN),
        "TDG5_runtime": _rel(TDG5_RUNTIME),
        "PROTO7_runtime": _rel(PROTO7_RUNTIME),
        "normal_flow_tracer_owner": _rel(TRACER_OWNER),
        "implementation_module": _rel(RUNTIME_MODULE),
    }
    if (
        config.get("schema_version") != 1
        or config.get("artifact_id") != ARTIFACT_ID
        or config.get("project_version") != PROJECT_VERSION
        or config.get("metric_signature") != "-+++"
        or config.get("riemann_convention") != "plus_partial_mu_gamma_nu"
        or any(config.get(key) != value for key, value in paths.items())
        or config.get("scope") != _expected_scope()
        or config.get("immutable_lineage") != _expected_lineage()
        or config.get("three_level_compositor") != _expected_compositor()
        or config.get("continuous_interval") != _expected_interval()
        or config.get("transaction_and_tracer") != _expected_transaction()
        or config.get("temporal_retry_and_checkpoint") != _expected_retry()
        or config.get("synthetic_qualification") != _expected_synthetic()
        or config.get("limitations") != _expected_limitations()
        or config.get("successor_boundary") != _expected_successor()
        or config.get("claims") != EXPECTED_CLAIMS
    ):
        raise ValueError("IMP2 identity, runtime contract, or claim boundary differs")
    proof = config.get("proof_contract", {})
    if not isinstance(proof, Mapping) or set(proof.values()) != {True}:
        raise ValueError("IMP2 proof contract must be all-of")


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    with path.open("rb") as handle:
        config = tomllib.load(handle)
    validate_config_data(config)
    return config


def _verify_lineage(config: Mapping[str, Any]) -> dict[str, Any]:
    lineage = config["immutable_lineage"]
    commit = lineage["PREF21_commit"]
    _git("cat-file", "-e", f"{commit}^{{commit}}")
    if _git("merge-base", "--is-ancestor", commit, "HEAD", check=False).returncode:
        raise ValueError("IMP2 PREF21 commit is not an ancestor of HEAD")
    tracked = {
        _rel(PREF21_CONFIG): lineage["PREF21_config_sha256"],
        _rel(PREF21_RESULT): lineage["PREF21_result_sha256"],
        _rel(PREF21_REPRODUCER): lineage["PREF21_reproducer_sha256"],
        _rel(PREF21_MODULE): lineage["PREF21_module_sha256"],
        _rel(PREF21_DOCUMENT): lineage["PREF21_document_sha256"],
        _rel(TDG6_DESIGN): lineage["TDG6_design_module_sha256"],
        _rel(TDG5_RUNTIME): lineage["TDG5_runtime_sha256"],
        _rel(PROTO7_RUNTIME): lineage["PROTO7_runtime_sha256"],
        _rel(TRACER_OWNER): lineage["normal_flow_tracer_owner_sha256"],
    }
    for relative, expected in tracked.items():
        committed = _git("show", f"{commit}:{relative}").stdout
        current = (REPOSITORY / relative).read_bytes()
        if (
            committed != current
            or _bytes_sha(committed) != expected
            or _bytes_sha(current) != expected
        ):
            raise ValueError(f"IMP2 bound predecessor differs: {relative}")
    raw = PREF21_RESULT.read_bytes()
    source = json.loads(raw)
    boundary = source.get("artifact_payload", {}).get("claim_boundary", {})
    binder = source.get("artifact_payload", {}).get("independent_binder", {})
    authorization = source.get("artifact_payload", {}).get(
        "implementation_authorization", {}
    )
    if (
        raw != _canonical_bytes(source)
        or source.get("artifact_id") != "FGC-1-TDG6-PREF21"
        or source.get("classification")
        != "independent_temporal_admission_binder_and_runtime_implementation_authorization"
        or source.get("gate_status")
        != "PASS_PRODUCTION_COMPOSITOR_IMPLEMENTATION_AUTHORIZED"
        or boundary.get("TDG6_independent_binder_completed") is not True
        or boundary.get("production_compositor_implementation_authorized") is not True
        or boundary.get("production_compositor_implemented") is not False
        or authorization.get("production_trajectory_authorized") is not False
        or boundary.get("PROTO14_frozen") is not False
        or boundary.get("GR0_case_eligible") is not False
        or binder.get("actual_PROTO13_history_arrays_consumed") is not False
        or binder.get("production_state_advanced") is not False
        or binder.get("physical_or_candidate_question_answered") is not False
        or authorization.get("historical_CAL10_result_reclassified") is not False
        or authorization.get("historical_PROTO13_stop_preserved") is not True
    ):
        raise ValueError("IMP2 source authorization differs")
    return {
        "PREF21_commit": commit,
        "PREF21_commit_is_ancestor_of_HEAD": True,
        "bound_predecessor_sha256": tracked,
        "PREF21_result_verified_canonical_and_implementation_only": True,
        "actual_PROTO13_history_arrays_consumed": False,
        "campaign_checkpoint_loaded": False,
        "campaign_checkpoint_resumed": False,
        "production_state_advanced": False,
        "historical_PROTO13_stop_preserved": True,
        "historical_CAL10_result_reclassified": False,
    }


def _state() -> EvolutionState:
    radii = np.linspace(0.0, 1.0, 9)
    u = np.zeros((9, 6))
    u[:, 0] = 1.0
    u[:, 1] = 0.025
    u[:, 2] = 1.0
    u[:, 3] = radii
    p = np.zeros_like(u)
    q = np.zeros_like(u)
    q[:, 3] = 1.0
    return EvolutionState(u, p, q)


def _rhs(*, failed_time: float | None = None, source_only: bool = False):
    def rhs(time: float, state: EvolutionState) -> EvolutionRHS:
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

    return rhs


def _transaction() -> GR0RuntimeStageTransaction:
    return GR0RuntimeStageTransaction(
        thresholds=GR0UniversalThresholds(),
        causal_state=CausalBudgetState(previous_speed_upper=1.0),
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
    ledger: runtime.TDG6TemporalLedger | None = None,
    rhs=None,
    step_size: float = 0.125,
) -> tuple[
    runtime.TDG6PreparedGR0Compositor,
    EvolutionState,
    GR0RuntimeStageTransaction,
    NormalFlowTracers,
    runtime.TDG6TemporalLedger,
]:
    current = _state() if state is None else state
    tx = _transaction() if transaction is None else transaction
    flow = _tracers(current) if tracers is None else tracers
    temporal = (
        runtime.TDG6TemporalLedger.zero(initial_time=0.0)
        if ledger is None
        else ledger
    )
    prepared = runtime.prepare_tdg6_gr0_compositor(
        method=method,
        time=0.0,
        step_size=step_size,
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
        np.asarray((value, 0.0, 0.0, 0.0), dtype=np.float64),
        (subintervals, 1),
    )
    return runtime._magnitude_interval(
        coefficients,
        np.zeros_like(coefficients),
        subinterval_count=subintervals,
        owned_row_count=1,
    )


def _synthetic_admission(
    *, method: str = PRIMARY_METHOD, failed_channel: str | None = None
) -> runtime.TDG6ContinuousAdmissionEvidence:
    outer = {
        channel: _interval(16.0, subintervals=2)
        for channel in TDG6_COMPLETE_STATE_CHANNELS
    }
    fine = {
        channel: _interval(1.0, subintervals=4)
        for channel in TDG6_COMPLETE_STATE_CHANNELS
    }
    if failed_channel is not None:
        outer[failed_channel] = _interval(1.0e-12, subintervals=2)
        fine[failed_channel] = _interval(9.0e-13, subintervals=4)
    return runtime.classify_tdg6_runtime_intervals(
        method=method, outer_intervals=outer, finest_intervals=fine
    )


def _method_control(method: str) -> dict[str, Any]:
    prepared, state, transaction, tracers, ledger = _prepare(method)
    monitor_before = transaction.state
    causal_before = transaction.causal_state
    tracer_before = runtime._tracer_snapshot(tracers)
    state_before = array_content_sha256(state.u, state.p, state.q)
    if (
        transaction.state != monitor_before
        or transaction.causal_state != causal_before
        or runtime._tracer_snapshot(tracers) != tracer_before
        or array_content_sha256(state.u, state.p, state.q) != state_before
        or ledger != prepared.original_temporal_ledger
    ):
        raise ValueError("TDG6 prepare mutated the accepted boundary")
    committed = runtime.commit_tdg6_gr0_compositor(
        prepared,
        transaction=transaction,
        tracers=tracers,
        temporal_ledger=ledger,
        current_time=0.0,
        current_state=state,
        current_step_index=0,
        current_transaction_serial=0,
    )
    counts = Counter(
        item.classification for item in prepared.continuous_admission.channel_admissions
    )
    return {
        "method": method,
        "shadow_proposal_count": sum(
            len(path.attempts)
            for path in (prepared.outer, prepared.medium, prepared.fine)
        ),
        "outer_stage_records": prepared.outer.stage_record_count,
        "medium_stage_records": prepared.medium.stage_record_count,
        "fine_stage_records": prepared.fine.stage_record_count,
        "total_shadow_stage_records": sum(
            path.stage_record_count
            for path in (prepared.outer, prepared.medium, prepared.fine)
        ),
        "classification_counts": dict(sorted(counts.items())),
        "complete_admission_passed": (
            prepared.continuous_admission.complete_admission_passed
        ),
        "failed_channels": list(prepared.continuous_admission.failed_channels),
        "maximum_fine_debit": max(
            prepared.continuous_admission.finest_pair_debit_vector
        ),
        "real_boundary_unchanged_during_prepare": True,
        "outer_and_medium_committed": False,
        "fine_quarter_steps_committed_atomically": (
            committed.fine_quarter_steps_committed_atomically
        ),
        "committed_time": committed.time,
        "committed_step_index": committed.step_index,
        "committed_transaction_serial": committed.transaction_serial,
        "committed_monitor_stage_count": transaction.state.accepted_stage_count,
        "accepted_macro_step_count": (
            committed.temporal_ledger.accepted_macro_step_count
        ),
        "committed_state_sha256": array_content_sha256(
            committed.state.u, committed.state.p, committed.state.q
        ),
        "fine_shadow_state_sha256": array_content_sha256(
            prepared.fine.final_accepted.state.u,
            prepared.fine.final_accepted.state.p,
            prepared.fine.final_accepted.state.q,
        ),
        "committed_tracer_snapshot": asdict(runtime._tracer_snapshot(tracers)),
        "fine_shadow_tracer_snapshot": asdict(
            prepared.fine.final_tracer_snapshot
        ),
    }


def _interval_controls() -> dict[str, Any]:
    zero = runtime.classify_tdg6_runtime_intervals(
        method=PRIMARY_METHOD,
        outer_intervals={
            channel: _interval(0.0, subintervals=2)
            for channel in TDG6_COMPLETE_STATE_CHANNELS
        },
        finest_intervals={
            channel: _interval(0.0, subintervals=4)
            for channel in TDG6_COMPLETE_STATE_CHANNELS
        },
    )
    failed_channel = TDG6_COMPLETE_STATE_CHANNELS[7]
    tiny = _synthetic_admission(failed_channel=failed_channel)
    large = runtime.classify_tdg6_runtime_intervals(
        method=PRIMARY_METHOD,
        outer_intervals={
            channel: _interval(80.0, subintervals=2)
            for channel in TDG6_COMPLETE_STATE_CHANNELS
        },
        finest_intervals={
            channel: _interval(10.0, subintervals=4)
            for channel in TDG6_COMPLETE_STATE_CHANNELS
        },
    )
    coefficients = np.asarray(
        ((0.0, 4.0, -4.0, 0.0), (0.0, 9.0 / 16.0, -1.5, 1.0))
    )
    dense_interval = runtime._magnitude_interval(
        coefficients,
        np.zeros_like(coefficients),
        subinterval_count=2,
        owned_row_count=1,
    )
    theta = np.linspace(0.0, 1.0, 131_073)
    dense_maximum = float(
        np.max(
            np.abs(
                coefficients[:, 0, None]
                + theta
                * (
                    coefficients[:, 1, None]
                    + theta
                    * (
                        coefficients[:, 2, None]
                        + theta * coefficients[:, 3, None]
                    )
                )
            )
        )
    )
    if not (
        zero.complete_admission_passed
        and {item.classification for item in zero.channel_admissions}
        == {"exact_zero"}
        and not tiny.complete_admission_passed
        and tiny.failed_channels == (failed_channel,)
        and tiny.channel_admissions[7].classification == "resolved_order_failure"
        and large.complete_admission_passed
        and min(large.finest_pair_debit_vector) > 9.0
        and dense_interval.lower_bound <= dense_maximum <= dense_interval.upper_bound
    ):
        raise ValueError("TDG6 interval controls failed")
    return {
        "exact_zero_classification": "exact_zero",
        "exact_zero_complete_admission_passed": True,
        "tiny_nonconvergent_failed_channel": failed_channel,
        "tiny_nonconvergent_classification": "resolved_order_failure",
        "tiny_nonconvergent_complete_admission_passed": False,
        "large_convergent_complete_admission_passed": True,
        "large_convergent_minimum_retained_debit": min(
            large.finest_pair_debit_vector
        ),
        "absolute_state_tolerance_used": False,
        "physical_signal_normalization_used": False,
        "dense_audit_sample_count": theta.size,
        "dense_audit_maximum": dense_maximum,
        "dense_audit_lower_bound": dense_interval.lower_bound,
        "dense_audit_upper_bound": dense_interval.upper_bound,
        "dense_audit_enclosed": True,
    }


class _FailOnceTracers(NormalFlowTracers):
    fail_next_commit: bool = False

    def commit_advance(self, positions, proper_times) -> None:
        if self.fail_next_commit:
            self.fail_next_commit = False
            raise RuntimeError("injected tracer commit failure")
        super().commit_advance(positions, proper_times)


def _failure_controls() -> dict[str, Any]:
    state = _state()
    transaction = _transaction()
    tracers = _tracers(state)
    ledger = runtime.TDG6TemporalLedger.zero(initial_time=0.0)
    monitor_before = transaction.state
    causal_before = transaction.causal_state
    tracer_before = runtime._tracer_snapshot(tracers)
    state_before = array_content_sha256(state.u, state.p, state.q)
    non_source_path = ""
    try:
        _prepare(
            state=state,
            transaction=transaction,
            tracers=tracers,
            ledger=ledger,
            rhs=_rhs(failed_time=0.078125),
        )
    except runtime.TDG6RefinementPathStop as error:
        if error.source_retry is not None or error.temporal_retry:
            raise ValueError("non-source failure ownership changed") from error
        non_source_path = error.path
    if (
        non_source_path != "fine_2"
        or transaction.state != monitor_before
        or transaction.causal_state != causal_before
        or runtime._tracer_snapshot(tracers) != tracer_before
        or array_content_sha256(state.u, state.p, state.q) != state_before
        or ledger != runtime.TDG6TemporalLedger.zero(initial_time=0.0)
    ):
        raise ValueError("late fine failure mutated a real boundary")

    source_path = ""
    source_owned = False
    try:
        _prepare(rhs=_rhs(failed_time=0.078125, source_only=True))
    except runtime.TDG6RefinementPathStop as error:
        source_path = error.path
        source_owned = error.source_retry_owned_by_PROTO7
        if error.source_retry is None or error.temporal_retry:
            raise ValueError("source failure ownership changed") from error

    prepared, state, transaction, tracers, ledger = _prepare()
    tracers.positions[0] = np.nextafter(tracers.positions[0], np.inf)
    drift_rejected = False
    try:
        runtime.commit_tdg6_gr0_compositor(
            prepared,
            transaction=transaction,
            tracers=tracers,
            temporal_ledger=ledger,
            current_time=0.0,
            current_state=state,
            current_step_index=0,
            current_transaction_serial=0,
        )
    except RuntimeError:
        drift_rejected = True

    state = _state()
    base = _tracers(state)
    tracers = _FailOnceTracers(**base.__dict__)
    transaction = _transaction()
    ledger = runtime.TDG6TemporalLedger.zero(initial_time=0.0)
    prepared, _, _, _, _ = _prepare(
        state=state,
        transaction=transaction,
        tracers=tracers,
        ledger=ledger,
    )
    monitor_before = transaction.state
    causal_before = transaction.causal_state
    tracer_before = runtime._tracer_snapshot(tracers)
    tracers.fail_next_commit = True
    rollback_passed = False
    try:
        runtime.commit_tdg6_gr0_compositor(
            prepared,
            transaction=transaction,
            tracers=tracers,
            temporal_ledger=ledger,
            current_time=0.0,
            current_state=state,
            current_step_index=0,
            current_transaction_serial=0,
        )
    except RuntimeError:
        rollback_passed = (
            transaction.state == monitor_before
            and transaction.causal_state == causal_before
            and runtime._tracer_snapshot(tracers) == tracer_before
        )
    if not (source_path == "fine_2" and source_owned and drift_rejected and rollback_passed):
        raise ValueError("TDG6 failure/atomicity controls failed")
    return {
        "late_non_source_stop_path": non_source_path,
        "late_non_source_stop_preserved_all_real_boundaries": True,
        "source_only_stop_path": source_path,
        "source_only_retry_owned_by_PROTO7": source_owned,
        "tracer_boundary_drift_rejected_before_commit": drift_rejected,
        "late_tracer_commit_exception_rolled_back_all_mutable_ledgers": rollback_passed,
    }


def _retry_checkpoint_controls() -> dict[str, Any]:
    prepared, state, transaction, tracers, ledger = _prepare()
    failed_channel = TDG6_COMPLETE_STATE_CHANNELS[7]
    failed = replace(
        prepared,
        continuous_admission=_synthetic_admission(failed_channel=failed_channel),
    )
    sink: list[dict[str, object]] = []
    try:
        runtime.require_tdg6_temporal_admission(
            failed,
            transaction=transaction,
            tracers=tracers,
            temporal_ledger=ledger,
            current_time=0.0,
            current_state=state,
            current_step_index=0,
            current_transaction_serial=0,
            durable_rejection_sink=lambda item: sink.append(dict(item)),
        )
    except runtime.TDG6TemporalRetryRequired as error:
        retried = error.updated_ledger
        retry_size = error.retry_step_size
        if error.physical_classification:
            raise ValueError("TDG6 temporal retry became physical") from error
    else:
        raise ValueError("TDG6 failed admission did not request retry")
    if len(sink) != 1 or len(sink[0].get("channel_admissions", [])) != 18:
        raise ValueError("TDG6 durable retry omitted complete channel evidence")

    metadata, debit = runtime.tdg6_checkpoint_extension(retried)
    restored = runtime.restore_tdg6_checkpoint_extension(metadata, debit)
    if restored != retried:
        raise ValueError("TDG6 checkpoint round trip differs")
    debit_mutation_rejected = False
    changed = debit.copy()
    changed[0] = np.nextafter(changed[0], np.inf)
    try:
        runtime.restore_tdg6_checkpoint_extension(metadata, changed)
    except ValueError:
        debit_mutation_rejected = True
    claim_mutation_rejected = False
    promoted = dict(metadata)
    promoted["physical_classification"] = True
    try:
        runtime.restore_tdg6_checkpoint_extension(promoted, debit)
    except ValueError:
        claim_mutation_rejected = True

    sink_failure_preserved_ledger = False
    try:
        runtime.require_tdg6_temporal_admission(
            failed,
            transaction=transaction,
            tracers=tracers,
            temporal_ledger=ledger,
            current_time=0.0,
            current_state=state,
            current_step_index=0,
            current_transaction_serial=0,
            durable_rejection_sink=lambda _item: (_ for _ in ()).throw(
                OSError("injected durable sink failure")
            ),
        )
    except OSError:
        sink_failure_preserved_ledger = (
            ledger.cumulative_temporal_retry_count == 0
            and ledger.serialized_temporal_rejections == ()
        )

    accepted, _, _, _, _ = _prepare(
        state=state,
        transaction=transaction,
        tracers=tracers,
        ledger=retried,
        step_size=retry_size,
    )
    committed = runtime.commit_tdg6_gr0_compositor(
        accepted,
        transaction=transaction,
        tracers=tracers,
        temporal_ledger=retried,
        current_time=0.0,
        current_state=state,
        current_step_index=0,
        current_transaction_serial=0,
    )
    reset_only_current = (
        committed.temporal_ledger.cumulative_temporal_retry_count == 1
        and committed.temporal_ledger.current_macro_step_temporal_retry_count == 0
        and committed.temporal_ledger.last_accepted_macro_step_temporal_retry_count
        == 1
        and len(committed.temporal_ledger.serialized_temporal_rejections) == 1
    )

    prepared, state, transaction, tracers, ledger = _prepare()
    failed = replace(
        prepared,
        continuous_admission=_synthetic_admission(failed_channel=failed_channel),
    )
    exhaustion_sink: list[dict[str, object]] = []
    exhausted_reason = ""
    for _ in range(TDG6_MAXIMUM_TEMPORAL_RETRIES_PER_STEP + 1):
        failed = replace(failed, original_temporal_ledger=ledger)
        try:
            runtime.require_tdg6_temporal_admission(
                failed,
                transaction=transaction,
                tracers=tracers,
                temporal_ledger=ledger,
                current_time=0.0,
                current_state=state,
                current_step_index=0,
                current_transaction_serial=0,
                durable_rejection_sink=lambda item: exhaustion_sink.append(dict(item)),
            )
        except runtime.TDG6TemporalRetryRequired as error:
            ledger = error.updated_ledger
        except runtime.TDG6TemporalRetryExhausted as error:
            ledger = error.updated_ledger
            exhausted_reason = error.reason
            break
    if not (
        debit_mutation_rejected
        and claim_mutation_rejected
        and sink_failure_preserved_ledger
        and reset_only_current
        and exhausted_reason == "maximum_temporal_retries"
        and len(exhaustion_sink) == TDG6_MAXIMUM_TEMPORAL_RETRIES_PER_STEP + 1
        and ledger.cumulative_temporal_retry_count
        == TDG6_MAXIMUM_TEMPORAL_RETRIES_PER_STEP + 1
    ):
        raise ValueError("TDG6 retry/checkpoint controls failed")
    return {
        "first_retry_step_size": retry_size,
        "durable_record_count_before_first_retry": 1,
        "first_retry_channel_decision_count": 18,
        "first_retry_physical_classification": False,
        "failed_sink_preserved_input_ledger": sink_failure_preserved_ledger,
        "accepted_step_reset_only_current_retry_counter": reset_only_current,
        "checkpoint_round_trip_passed": True,
        "checkpoint_debit_sha256": metadata["accumulated_debit_sha256"],
        "checkpoint_debit_mutation_rejected": debit_mutation_rejected,
        "checkpoint_claim_promotion_rejected": claim_mutation_rejected,
        "retry_exhaustion_reason": exhausted_reason,
        "retry_exhaustion_record_count": len(exhaustion_sink),
        "retry_exhaustion_is_physical_classification": False,
        "minimum_macro_step": float(TDG6_MINIMUM_MACRO_STEP),
    }


def _synthetic_certificate() -> dict[str, Any]:
    methods = {
        "RK4": _method_control(PRIMARY_METHOD),
        "SSPRK3": _method_control(COMPARATOR_METHOD),
    }
    interval = _interval_controls()
    failures = _failure_controls()
    retry = _retry_checkpoint_controls()
    if (
        methods["RK4"]["total_shadow_stage_records"] != 35
        or methods["RK4"]["fine_stage_records"] != 20
        or methods["SSPRK3"]["total_shadow_stage_records"] != 28
        or methods["SSPRK3"]["fine_stage_records"] != 16
        or any(
            item["shadow_proposal_count"] != 7
            or item["complete_admission_passed"] is not True
            or item["outer_and_medium_committed"] is not False
            or item["fine_quarter_steps_committed_atomically"] is not True
            or item["committed_state_sha256"] != item["fine_shadow_state_sha256"]
            or item["committed_tracer_snapshot"]
            != item["fine_shadow_tracer_snapshot"]
            for item in methods.values()
        )
    ):
        raise ValueError("TDG6 method controls failed")
    return {
        "runtime_environment": {
            "python_implementation": sys.implementation.name,
            "python_major_minor": f"{sys.version_info.major}.{sys.version_info.minor}",
            "numpy": np.__version__,
            "binary64_epsilon": float(np.finfo(np.float64).eps),
        },
        "method_controls": methods,
        "continuous_interval_controls": interval,
        "failure_and_atomicity_controls": failures,
        "retry_and_checkpoint_controls": retry,
        "complete_owned_state": "dimensionless_u_p_q",
        "channel_order": list(TDG6_COMPLETE_STATE_CHANNELS),
        "synthetic_states_only": True,
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
    malformed_channel_order = False
    complete = _synthetic_admission()
    outer = {item.channel: item.outer_difference for item in complete.channel_admissions}
    fine = {item.channel: item.finest_difference for item in complete.channel_admissions}
    try:
        runtime.classify_tdg6_runtime_intervals(
            method=PRIMARY_METHOD,
            outer_intervals=dict(reversed(tuple(outer.items()))),
            finest_intervals=fine,
        )
    except ValueError:
        malformed_channel_order = True
    controls = {
        "runtime_reordered_channel_map_rejected": malformed_channel_order,
        "PREF21_commit_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["immutable_lineage"].__setitem__(
                "PREF21_commit", "0" * 40
            ),
        ),
        "PREF21_result_hash_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["immutable_lineage"].__setitem__(
                "PREF21_result_sha256", "0" * 64
            ),
        ),
        "fine_only_commit_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["three_level_compositor"].__setitem__(
                "outer_and_medium_are_evidence_only", False
            ),
        ),
        "threshold_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["continuous_interval"].__setitem__(
                "exact_squared_admission", "7*U12^2<=L01^2"
            ),
        ),
        "retry_limit_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["temporal_retry_and_checkpoint"].__setitem__(
                "maximum_temporal_retries_per_macro_step", 33
            ),
        ),
        "trajectory_promotion_rejected": _expect_config_rejection(
            config,
            lambda value: value["claims"].__setitem__(
                "production_trajectory_authorized", True
            ),
        ),
        "PROTO14_promotion_rejected": _expect_config_rejection(
            config,
            lambda value: value["claims"].__setitem__("PROTO14_frozen", True),
        ),
        "GR0_eligibility_promotion_rejected": _expect_config_rejection(
            config,
            lambda value: value["claims"].__setitem__("GR0_case_eligible", True),
        ),
        "candidate_promotion_rejected": _expect_config_rejection(
            config,
            lambda value: value["claims"].__setitem__(
                "FGCQR_holdout_execution_authorized", True
            ),
        ),
        "PDE_history_promotion_rejected": _expect_config_rejection(
            config,
            lambda value: value["limitations"].__setitem__(
                "piecewise_cubic_extension_is_not_exact_PDE_history", False
            ),
        ),
        "historical_reclassification_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["limitations"].__setitem__(
                "historical_CAL10_and_PROTO13_results_are_not_reclassified", False
            ),
        ),
    }
    if set(controls.values()) != {True}:
        raise ValueError(f"IMP2 mutation control failed: {controls}")
    return controls


def record(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    config = load_config(config_path)
    lineage = _verify_lineage(config)
    synthetic = _synthetic_certificate()
    mutations = _mutation_controls(config)
    result = {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "metric_signature": "-+++",
        "riemann_convention": "plus_partial_mu_gamma_nu",
        "classification": (
            "production_TDG6_three_level_compositor_implemented_and_"
            "synthetically_qualified_PROTO14_design_only_authorized"
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
            "three_level_compositor": deepcopy(config["three_level_compositor"]),
            "continuous_interval": deepcopy(config["continuous_interval"]),
            "transaction_and_tracer": deepcopy(config["transaction_and_tracer"]),
            "temporal_retry_and_checkpoint": deepcopy(
                config["temporal_retry_and_checkpoint"]
            ),
            "synthetic_qualification": synthetic,
            "limitations": deepcopy(config["limitations"]),
            "mutation_controls": mutations,
            "proof_contract": dict(config["proof_contract"]),
            "successor_boundary": dict(config["successor_boundary"]),
            "claim_boundary": {
                "actual_PROTO13_history_arrays_consumed": False,
                "campaign_checkpoint_loaded": False,
                "campaign_checkpoint_resumed": False,
                "production_state_advanced": False,
                "synthetic_states_advanced": True,
                "historical_PROTO13_stop_preserved": True,
                "historical_CAL10_result_reclassified": False,
                "piecewise_cubic_extension_is_exact_PDE_history": False,
                "trajectory_asymptotic_regime_proved": False,
                "rigorous_global_PDE_error_bound_proved": False,
                "production_compositor_implemented": True,
                "production_compositor_synthetic_qualification_passed": True,
                "PROTO14_freeze_design_authorized": True,
                "production_trajectory_authorized": False,
                "PROTO14_frozen": False,
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
    boundary = payload.get("claim_boundary", {})
    synthetic = payload.get("synthetic_qualification", {})
    if (
        value.get("schema_version") != 1
        or value.get("artifact_id") != ARTIFACT_ID
        or value.get("project_version") != PROJECT_VERSION
        or value.get("classification")
        != "production_TDG6_three_level_compositor_implemented_and_synthetically_qualified_PROTO14_design_only_authorized"
        or value.get("generated_by") != "scripts/reproduce_fgc_tdg6_imp2.py"
        or value.get("gate_status") != EXPECTED_CLAIMS
        or value.get("nonclaims")
        != {
            key: False
            for key, item in sorted(EXPECTED_CLAIMS.items())
            if item is False
        }
        or value.get("source_config_sha256")
        != {_rel(DEFAULT_CONFIG): _sha(DEFAULT_CONFIG)}
        or payload.get("three_level_compositor")
        != config["three_level_compositor"]
        or payload.get("continuous_interval") != config["continuous_interval"]
        or payload.get("transaction_and_tracer")
        != config["transaction_and_tracer"]
        or payload.get("temporal_retry_and_checkpoint")
        != config["temporal_retry_and_checkpoint"]
        or payload.get("limitations") != config["limitations"]
        or payload.get("successor_boundary") != config["successor_boundary"]
        or synthetic.get("all_synthetic_controls_pass") is not True
        or synthetic.get("synthetic_states_only") is not True
        or boundary.get("production_compositor_implemented") is not True
        or boundary.get("production_compositor_synthetic_qualification_passed")
        is not True
        or boundary.get("PROTO14_freeze_design_authorized") is not True
        or boundary.get("production_trajectory_authorized") is not False
        or boundary.get("PROTO14_frozen") is not False
        or boundary.get("GR0_case_eligible") is not False
        or boundary.get("candidate_branch_opened") is not False
        or boundary.get("physical_question_answered") is not False
        or boundary.get("actual_PROTO13_history_arrays_consumed") is not False
        or boundary.get("production_state_advanced") is not False
        or boundary.get("historical_CAL10_result_reclassified") is not False
    ):
        raise ValueError("IMP2 stored result boundary differs")
    implementation = value.get("implementation_sha256", {})
    if (
        not isinstance(implementation, Mapping)
        or set(implementation) != {_rel(path) for path in IMPLEMENTATION}
        or any(
            _sha(REPOSITORY / relative) != expected
            for relative, expected in implementation.items()
        )
    ):
        raise ValueError("IMP2 implementation hash ledger differs")
    document = REPOSITORY / str(value.get("derivation_document", ""))
    if (
        value.get("derivation_document") != "docs/fgc-tdg6-imp2.md"
        or not document.is_file()
        or _sha(document) != value.get("derivation_document_sha256")
    ):
        raise ValueError("IMP2 derivation document differs")
    if (
        not isinstance(payload.get("mutation_controls"), Mapping)
        or set(payload["mutation_controls"].values()) != {True}
        or not isinstance(payload.get("proof_contract"), Mapping)
        or set(payload["proof_contract"].values()) != {True}
    ):
        raise ValueError("IMP2 proof or mutation ledger differs")


def verify_canonical(
    config_path: Path = DEFAULT_CONFIG,
    output_path: Path = DEFAULT_OUTPUT,
) -> dict[str, Any]:
    config = load_config(config_path)
    raw = output_path.read_bytes()
    stored = json.loads(raw)
    if not isinstance(stored, dict) or raw != _canonical_bytes(stored):
        raise ValueError("IMP2 stored result is not canonical JSON")
    _validate_stored_record(stored, config)
    observed = record(config_path)
    if stored != observed:
        raise ValueError("IMP2 stored result differs from reproduction")
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
        status = value["gate_status"]
        print(
            f"PASS {ARTIFACT_ID}: compositor="
            f"{status['production_compositor_implemented']} synthetic="
            f"{status['production_compositor_synthetic_qualification_passed']} "
            f"PROTO14_design={status['PROTO14_freeze_design_authorized']} "
            f"trajectory={status['production_trajectory_authorized']}"
        )
    else:
        print(_canonical(record(args.config)), end="")


if __name__ == "__main__":
    main()
