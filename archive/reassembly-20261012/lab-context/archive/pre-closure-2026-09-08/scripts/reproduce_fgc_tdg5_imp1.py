#!/usr/bin/env python3
"""Bind the TDG5 atomic runtime-pair implementation and synthetic controls."""

from __future__ import annotations

import argparse
from copy import deepcopy
from dataclasses import asdict
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
from recursive_horizons.fgc.evolution.proto7_runtime import (  # noqa: E402
    Proto7TerminalStop,
)
from recursive_horizons.fgc.evolution.tdg5_stage_complete_refinement_runtime import (  # noqa: E402
    TDG5RefinementPathStop,
    binary64_cubic_absolute_envelope,
    commit_tdg5_gr0_refinement_pair,
    prepare_tdg5_gr0_refinement_pair,
)


ARTIFACT_ID = "FGC-1-TDG5-IMP1"
PROJECT_VERSION = "0.11.0"
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-tdg5-imp1.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-tdg5-imp1.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-tdg5-imp1.md"
SOURCE_RESULT = REPOSITORY / "results/fgc-1-tdg5-pref20.json"
SOURCE_CONFIG = REPOSITORY / "configs/fgc/fgc-1-tdg5-pref20.toml"
SOURCE_REPRODUCER = REPOSITORY / "scripts/reproduce_fgc_tdg5_pref20.py"
SOURCE_THEOREM = (
    REPOSITORY
    / "src/recursive_horizons/fgc/evolution/tdg5_stage_complete_refinement_theorem.py"
)
SOURCE_DOCUMENT = REPOSITORY / "docs/fgc-tdg5-pref20.md"
NUMERICAL_ENGINE = (
    REPOSITORY / "src/recursive_horizons/fgc/evolution/numerical_engine.py"
)
PROTO7_RUNTIME = (
    REPOSITORY / "src/recursive_horizons/fgc/evolution/proto7_runtime.py"
)
GR0_TRANSACTION = (
    REPOSITORY / "src/recursive_horizons/fgc/evolution/proto5_runtime.py"
)
RUNTIME_MODULE = (
    REPOSITORY
    / "src/recursive_horizons/fgc/evolution/tdg5_stage_complete_refinement_runtime.py"
)
IMPLEMENTATION = (Path(__file__).resolve(), RUNTIME_MODULE)


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
    "runtime_refinement_pair_implemented": True,
    "runtime_refinement_pair_synthetic_validation_passed": True,
    "runtime_threshold_design_authorized": True,
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
        "role": "post_theorem_runtime_pair_implementation_and_synthetic_validation",
        "calibration_branch": "GR-0",
        "source_compact_result_read": True,
        "actual_terminal_history_arrays_consumed": False,
        "campaign_checkpoint_loaded": False,
        "campaign_checkpoint_resumed": False,
        "production_state_advanced": False,
        "synthetic_states_advanced": True,
        "SGBL_trajectory_read": False,
        "FGCQR_trajectory_read": False,
        "runtime_thresholds_frozen": False,
        "replacement_temporal_admission_defined": False,
        "physical_or_candidate_question_answered": False,
    }


def _expected_lineage() -> dict[str, Any]:
    return {
        "checkpoint_commit": "76da6873f4625c2db930f861fae5f4361d70d2c2",
        "source_config_sha256": "cb22ee7ca18ddd07b9b0f1ad411428c030d4a4d4a203caebe3c853251389f082",
        "source_result_sha256": "abb81aa4a6df47a90a9011eb5c1c9f1065a495077935d6164d5bd7c2f93df55a",
        "source_reproducer_sha256": "c2f1e3b61ac51ddcf21840bc0cabbb6b88d76bcb7b8f9d83c6d0155b39dc6220",
        "source_theorem_module_sha256": "07bd73402b750f51961cf580072144b332c2a384f6240b9f28f61ef0f613bff7",
        "source_document_sha256": "7f9fbf8022ae3680bcf6fdee73ee7d150d33cfdccfc3c8afae9eeef61f63e45a",
        "numerical_engine_sha256": "8ea99604a85b4e1acbf869ed2462ea4aa8301c06e068330a1d0ce66a51a0bebf",
        "PROTO7_runtime_sha256": "ffa75767c1fdbae42a09c8ef50276183861fc9255ad0061e5a1ee3673f142734",
        "GR0_transaction_runtime_sha256": "607c9d5c353e37f8c68c032ec70df773c326335ba723a7db7317191bc7260545",
        "checkpoint_commit_must_be_ancestor_of_HEAD": True,
        "all_bound_blobs_must_match_checkpoint_commit_and_worktree": True,
        "source_result_must_be_canonical_and_authorize_implementation_only": True,
        "raw_campaign_bundle_must_not_be_loaded": True,
    }


def _expected_runtime_contract() -> dict[str, Any]:
    return {
        "complete_state": "dimensionless_u_p_q",
        "owned_rows": "noncentre_rows_excluding_four_projector_owned_outer_rows",
        "same_bitwise_initial_state_required": True,
        "same_method_source_projector_and_spatial_operator_required": True,
        "coarse_path": "one_full_step_of_width_Delta_t_on_isolated_shadow_transaction",
        "fine_path": "two_half_steps_of_width_Delta_t_over_2_on_one_isolated_shadow_transaction",
        "coarse_path_is_evidence_only": True,
        "fine_path_is_the_only_committable_path": True,
        "fine_half_steps_commit_atomically_or_neither_commits": True,
        "commit_token_must_revalidate_state_time_step_serial_monitor_and_causal_ledgers": True,
        "unchanged_PROTO7_source_health_CFL_boundary_and_transaction_guards_required": True,
        "fresh_endpoint_RHS_required_for_all_three_proposals": True,
        "shared_initial_and_fine_midpoint_RHS_fields_must_match_bitwise": True,
        "RK4_formal_order": 4,
        "RK4_Richardson_denominator": 15,
        "SSPRK3_formal_order": 3,
        "SSPRK3_Richardson_denominator": 7,
    }


def _expected_envelope() -> dict[str, Any]:
    return {
        "declared_interval_class": "P3_on_each_accepted_normalized_time_interval",
        "both_endpoints_and_every_numerically_real_interior_derivative_root_evaluated": True,
        "scaled_quadratic_formula_avoids_derivative_overflow": True,
        "near_degenerate_discriminants_are_counted_explicitly": True,
        "coefficient_construction_roundoff_is_separate": True,
        "root_and_value_binary64_roundoff_is_separate": True,
        "independent_Bernstein_convex_hull_envelope_is_computed": True,
        "Bernstein_certification_slack_is_separate": True,
        "certified_upper_bound_dominates_candidate_and_Bernstein_routes": True,
        "exact_zero_must_remain_zero": True,
        "conditional_fine_path_Richardson_debit_uses_method_denominator": True,
    }


def _expected_synthetic() -> dict[str, Any]:
    return {
        "methods": ["RK4", "SSPRK3"],
        "grid_point_count": 9,
        "field_count": 6,
        "owned_row_count": 4,
        "step_size": "1/8",
        "endpoint_only_bump_raw_maximum": "1",
        "endpoint_only_bump_real_interior_root_count": 1,
        "genuine_cubic_real_interior_root_count": 2,
        "dense_audit_sample_count": 131073,
        "zero_control_has_exact_zero_envelope": True,
        "fine_right_injected_non_source_stop_time": "3/32",
        "fine_right_injected_stop_must_leave_real_state_and_ledgers_unchanged": True,
        "RK4_committed_stage_count": 10,
        "RK4_member_transaction_serial": 10,
        "SSPRK3_committed_stage_count": 8,
        "SSPRK3_member_transaction_serial": 8,
        "coarse_endpoint_must_not_equal_committed_fine_endpoint_on_nontrivial_control": True,
        "commit_drift_mutation_must_fail_closed": True,
        "one_bit_coefficient_mutation_must_change_evidence": True,
    }


def _expected_limitations() -> dict[str, Any]:
    return {
        "declared_cubic_is_exact_PDE_history": False,
        "finite_dimensional_envelope_proves_PDE_accuracy": False,
        "trajectory_asymptotic_regime_proved": False,
        "rigorous_global_PDE_error_bound_proved": False,
        "Bernstein_slack_is_arithmetic_roundoff": False,
        "conditional_Richardson_debit_is_an_admission_threshold": False,
        "runtime_thresholds_frozen": False,
        "production_compositor_integrated": False,
        "production_trajectory_authorized": False,
        "historical_PROTO13_stop_preserved": True,
        "historical_CAL10_result_reclassified": False,
        "spatial_ladders_cross_method_constraints_and_DEF1_ledger_remain_required": True,
    }


def _expected_successor() -> dict[str, Any]:
    return {
        key: EXPECTED_CLAIMS[key]
        for key in (
            "TDG4_sampling_identifiability_theorem_completed",
            "current_temporal_spectral_rule_retired_as_continuum_admission",
            "current_temporal_spectral_rule_retained_as_sampled_diagnostic",
            "replacement_temporal_gate_design_authorized",
            "TDG5_replacement_design_frozen",
            "TDG5_exact_no_history_controls_pass",
            "TDG5_stage_complete_refinement_theorem_completed",
            "runtime_refinement_pair_implementation_authorized",
            "runtime_refinement_pair_implemented",
            "runtime_refinement_pair_synthetic_validation_passed",
            "runtime_threshold_design_authorized",
            "runtime_thresholds_frozen",
            "replacement_temporal_admission_defined",
            "PROTO14_frozen",
            "fresh_GR0_dynamic_calibration_completed",
            "GR0_case_eligible",
            "SGBL_execution_authorized",
            "FGCQR_holdout_execution_authorized",
        )
    } | {"DEF1_execution_authorized": False}


def validate_config_data(config: Mapping[str, Any]) -> None:
    _strict_keys(
        "IMP1 config",
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
            "PROTO7_runtime",
            "GR0_transaction_runtime",
            "runtime_module",
            "scope",
            "immutable_lineage",
            "runtime_contract",
            "continuous_envelope",
            "synthetic_validation",
            "limitation_audit",
            "proof_contract",
            "successor_boundary",
            "claims",
        },
    )
    expected_paths = {
        "source_result": _rel(SOURCE_RESULT),
        "source_config": _rel(SOURCE_CONFIG),
        "source_reproducer": _rel(SOURCE_REPRODUCER),
        "source_theorem_module": _rel(SOURCE_THEOREM),
        "source_document": _rel(SOURCE_DOCUMENT),
        "numerical_engine": _rel(NUMERICAL_ENGINE),
        "PROTO7_runtime": _rel(PROTO7_RUNTIME),
        "GR0_transaction_runtime": _rel(GR0_TRANSACTION),
        "runtime_module": _rel(RUNTIME_MODULE),
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
        or config.get("runtime_contract") != _expected_runtime_contract()
        or config.get("continuous_envelope") != _expected_envelope()
        or config.get("synthetic_validation") != _expected_synthetic()
        or config.get("limitation_audit") != _expected_limitations()
        or config.get("successor_boundary") != _expected_successor()
        or config.get("claims") != EXPECTED_CLAIMS
    ):
        raise ValueError("IMP1 identity, runtime contract, or claim boundary differs")
    if set(config.get("proof_contract", {}).values()) != {True}:
        raise ValueError("IMP1 proof contract must be all-of")


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    with path.open("rb") as handle:
        config = tomllib.load(handle)
    validate_config_data(config)
    return config


def _verify_lineage(config: Mapping[str, Any]) -> dict[str, Any]:
    lineage = config["immutable_lineage"]
    commit = lineage["checkpoint_commit"]
    _git("cat-file", "-e", f"{commit}^{{commit}}")
    if _git("merge-base", "--is-ancestor", commit, "HEAD", check=False).returncode:
        raise ValueError("IMP1 checkpoint commit is not an ancestor of HEAD")
    tracked = {
        _rel(SOURCE_CONFIG): lineage["source_config_sha256"],
        _rel(SOURCE_RESULT): lineage["source_result_sha256"],
        _rel(SOURCE_REPRODUCER): lineage["source_reproducer_sha256"],
        _rel(SOURCE_THEOREM): lineage["source_theorem_module_sha256"],
        _rel(SOURCE_DOCUMENT): lineage["source_document_sha256"],
        _rel(NUMERICAL_ENGINE): lineage["numerical_engine_sha256"],
        _rel(PROTO7_RUNTIME): lineage["PROTO7_runtime_sha256"],
        _rel(GR0_TRANSACTION): lineage["GR0_transaction_runtime_sha256"],
    }
    for relative, expected_hash in tracked.items():
        committed = _git("show", f"{commit}:{relative}").stdout
        current = (REPOSITORY / relative).read_bytes()
        if (
            committed != current
            or _bytes_sha(committed) != expected_hash
            or _bytes_sha(current) != expected_hash
        ):
            raise ValueError(f"IMP1 bound blob differs: {relative}")

    raw = SOURCE_RESULT.read_bytes()
    source = json.loads(raw)
    status = source.get("gate_status", {})
    boundary = source.get("artifact_payload", {}).get("claim_boundary", {})
    if (
        raw != _canonical_bytes(source)
        or source.get("artifact_id") != "FGC-1-TDG5-PREF20"
        or source.get("classification")
        != "completed_TDG5_stage_complete_refinement_theorem_runtime_pair_implementation_authorized"
        or status.get("TDG5_stage_complete_refinement_theorem_completed") is not True
        or status.get("runtime_refinement_pair_implementation_authorized") is not True
        or status.get("runtime_refinement_pair_implemented") is not False
        or status.get("runtime_thresholds_frozen") is not False
        or status.get("replacement_temporal_admission_defined") is not False
        or boundary.get("actual_terminal_history_arrays_consumed") is not False
        or boundary.get("campaign_checkpoint_loaded") is not False
        or boundary.get("state_advanced") is not False
        or boundary.get("historical_PROTO13_stop_preserved") is not True
        or boundary.get("historical_CAL10_result_reclassified") is not False
    ):
        raise ValueError("IMP1 source authorization differs")
    return {
        "checkpoint_commit": commit,
        "checkpoint_commit_is_ancestor_of_HEAD": True,
        "bound_blob_sha256": tracked,
        "source_result_verified_canonical_and_implementation_only": True,
        "actual_terminal_history_arrays_consumed": False,
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
    u[:, 2] = 1.0
    u[:, 3] = radii
    p = np.zeros_like(u)
    q = np.zeros_like(u)
    q[:, 3] = 1.0
    return EvolutionState(u, p, q)


def _rhs(*, failed_time: float | None = None):
    def rhs(time: float, state: EvolutionState) -> EvolutionRHS:
        failed = failed_time is not None and (
            np.float64(time).tobytes() == np.float64(failed_time).tobytes()
        )
        return EvolutionRHS(
            state.u,
            state.p,
            state.q,
            {
                "source_residual_infinity": 0.0,
                "source_refinement_iterations": 0,
                "source_residual_decreased_monotonically": True,
                "kinetic_condition_infinity": 1.0,
                "coordinate_speed_upper": 1.0,
                "minimum_lapse": 0.0 if failed else 1.0,
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


def _method_control(method: str) -> dict[str, Any]:
    initial = _state()
    transaction = _transaction()
    monitor_before = transaction.state
    causal_before = transaction.causal_state
    prepared = prepare_tdg5_gr0_refinement_pair(
        method=method,
        time=0.0,
        step_size=0.125,
        state=initial,
        rhs=_rhs(),
        projector=None,
        transaction=transaction,
        previous_step_index=0,
        previous_transaction_serial=0,
    )
    if transaction.state != monitor_before or transaction.causal_state != causal_before:
        raise ValueError("IMP1 prepare mutated the real transaction")
    committed = commit_tdg5_gr0_refinement_pair(
        prepared,
        transaction=transaction,
        current_time=0.0,
        current_state=initial,
        current_step_index=0,
        current_transaction_serial=0,
    )
    coarse = prepared.coarse.accepted
    fine = prepared.fine_right.accepted
    if coarse is None or fine is None:
        raise ValueError("IMP1 synthetic pair omitted an accepted shadow path")
    coarse_hash = array_content_sha256(
        coarse.state.u, coarse.state.p, coarse.state.q
    )
    fine_hash = array_content_sha256(fine.state.u, fine.state.p, fine.state.q)
    committed_hash = array_content_sha256(
        committed.state.u, committed.state.p, committed.state.q
    )
    evidence = prepared.continuous_difference
    if (
        coarse_hash == fine_hash
        or committed_hash != fine_hash
        or transaction.state.accepted_stage_count
        != len(prepared.fine_left.proposal.stages)
        + len(prepared.fine_right.proposal.stages)
        or committed.coarse_path_committed is not False
        or committed.fine_half_steps_committed_atomically is not True
    ):
        raise ValueError("IMP1 synthetic commit ownership differs")
    return {
        "method": evidence.method_label,
        "formal_order": evidence.formal_order,
        "Richardson_denominator": evidence.richardson_denominator,
        "coarse_stage_count": len(prepared.coarse.proposal.stages),
        "fine_left_stage_count": len(prepared.fine_left.proposal.stages),
        "fine_right_stage_count": len(prepared.fine_right.proposal.stages),
        "committed_stage_count": transaction.state.accepted_stage_count,
        "committed_member_transaction_serial": committed.transaction_serial,
        "transaction_last_stage_serial": transaction.state.last_transaction_serial,
        "committed_step_index": committed.step_index,
        "coarse_endpoint_sha256": coarse_hash,
        "fine_endpoint_sha256": fine_hash,
        "committed_endpoint_sha256": committed_hash,
        "coarse_endpoint_differs_from_fine_endpoint": coarse_hash != fine_hash,
        "prepared_real_monitor_unchanged": True,
        "prepared_real_causal_ledger_unchanged": True,
        "coarse_path_committed": False,
        "fine_half_steps_committed_atomically": True,
        "continuous_difference": asdict(evidence),
    }


def _continuous_controls() -> dict[str, Any]:
    bump = binary64_cubic_absolute_envelope(
        np.asarray((0.0, 4.0, -4.0, 0.0))
    )
    cubic = binary64_cubic_absolute_envelope(
        np.asarray((0.0, 9.0 / 16.0, -1.5, 1.0))
    )
    zero = binary64_cubic_absolute_envelope(np.zeros((3, 4)))
    controls = np.asarray(
        (
            (0.0, 4.0, -4.0, 0.0),
            (0.0, 9.0 / 16.0, -1.5, 1.0),
            (0.25, -0.75, 0.125, 0.5),
            (-1.0, 2.0, -3.0, 1.25),
        ),
        dtype=np.float64,
    )
    combined = binary64_cubic_absolute_envelope(controls)
    theta = np.linspace(0.0, 1.0, 131_073)
    dense = float(
        np.max(
            np.abs(
                controls[:, 0, None]
                + theta
                * (
                    controls[:, 1, None]
                    + theta
                    * (controls[:, 2, None] + theta * controls[:, 3, None])
                )
            )
        )
    )
    mutation = controls[0].copy()
    mutation[1] = np.nextafter(mutation[1], np.inf)
    mutated = binary64_cubic_absolute_envelope(mutation)
    if (
        bump.raw_candidate_maximum != 1.0
        or bump.real_interior_root_count != 1
        or cubic.real_interior_root_count != 2
        or any(asdict(zero).get(name) != 0.0 for name in (
            "raw_candidate_maximum",
            "coefficient_construction_debit",
            "outward_arithmetic_debit",
            "bernstein_certification_slack",
            "certified_continuous_upper_bound",
        ))
        or dense > combined.certified_continuous_upper_bound
        or asdict(mutated) == asdict(bump)
    ):
        raise ValueError("IMP1 continuous synthetic control failed")
    return {
        "endpoint_only_bump": asdict(bump),
        "genuine_cubic_two_root_control": asdict(cubic),
        "exact_zero_control": asdict(zero),
        "dense_audit": {
            "polynomial_count": controls.shape[0],
            "sample_count_per_polynomial": theta.size,
            "dense_observed_maximum": dense,
            "certified_continuous_upper_bound": (
                combined.certified_continuous_upper_bound
            ),
            "certified_bound_dominates_dense_audit": (
                combined.certified_continuous_upper_bound >= dense
            ),
        },
        "one_bit_mutation": {
            "baseline_coefficient": float(controls[0, 1]),
            "mutated_coefficient": float(mutation[1]),
            "evidence_changed": asdict(mutated) != asdict(bump),
            "mutated_evidence": asdict(mutated),
        },
    }


def _atomic_failure_controls() -> dict[str, Any]:
    initial = _state()
    transaction = _transaction()
    monitor_before = transaction.state
    causal_before = transaction.causal_state
    initial_hash = array_content_sha256(initial.u, initial.p, initial.q)
    late_failure = False
    failure_reason = None
    try:
        prepare_tdg5_gr0_refinement_pair(
            method=PRIMARY_METHOD,
            time=0.0,
            step_size=0.125,
            state=initial,
            rhs=_rhs(failed_time=0.09375),
            projector=None,
            transaction=transaction,
            previous_step_index=0,
            previous_transaction_serial=0,
        )
    except TDG5RefinementPathStop as error:
        late_failure = (
            error.path == "fine_right"
            and isinstance(error.cause, Proto7TerminalStop)
            and error.cause.reason == "nonpositive_lapse"
            and error.real_transaction_advanced is False
            and error.accepted_state_advanced is False
        )
        failure_reason = None if error.cause is None else error.cause.reason
    unchanged = (
        transaction.state == monitor_before
        and transaction.causal_state == causal_before
        and array_content_sha256(initial.u, initial.p, initial.q) == initial_hash
    )

    drift_transaction = _transaction()
    prepared = prepare_tdg5_gr0_refinement_pair(
        method=PRIMARY_METHOD,
        time=0.0,
        step_size=0.125,
        state=initial,
        rhs=_rhs(),
        projector=None,
        transaction=drift_transaction,
        previous_step_index=0,
        previous_transaction_serial=0,
    )
    drift_transaction.causal_state = CausalBudgetState(previous_speed_upper=1.25)
    drifted = drift_transaction.causal_state
    drift_rejected = False
    try:
        commit_tdg5_gr0_refinement_pair(
            prepared,
            transaction=drift_transaction,
            current_time=0.0,
            current_state=initial,
            current_step_index=0,
            current_transaction_serial=0,
        )
    except RuntimeError:
        drift_rejected = (
            drift_transaction.causal_state == drifted
            and drift_transaction.state.accepted_stage_count == 0
        )
    if not late_failure or not unchanged or not drift_rejected:
        raise ValueError("IMP1 atomic failure control failed")
    return {
        "fine_right_injected_stop_time": 0.09375,
        "fine_right_injected_stop_path": "fine_right",
        "fine_right_injected_stop_reason": failure_reason,
        "real_state_monitor_and_causal_ledgers_unchanged": unchanged,
        "post_prepare_causal_ledger_drift_rejected": drift_rejected,
        "no_coarse_or_partial_fine_commit_on_failure": True,
    }


def _synthetic_certificate() -> dict[str, Any]:
    continuous = _continuous_controls()
    methods = {
        "RK4": _method_control(PRIMARY_METHOD),
        "SSPRK3": _method_control(COMPARATOR_METHOD),
    }
    atomic = _atomic_failure_controls()
    if (
        methods["RK4"]["committed_stage_count"] != 10
        or methods["RK4"]["committed_member_transaction_serial"] != 10
        or methods["SSPRK3"]["committed_stage_count"] != 8
        or methods["SSPRK3"]["committed_member_transaction_serial"] != 8
        or any(
            item["coarse_endpoint_differs_from_fine_endpoint"] is not True
            or item["fine_half_steps_committed_atomically"] is not True
            for item in methods.values()
        )
    ):
        raise ValueError("IMP1 method-owned synthetic controls differ")
    return {
        "runtime_environment": {
            "python_implementation": sys.implementation.name,
            "python_major_minor": f"{sys.version_info.major}.{sys.version_info.minor}",
            "numpy": np.__version__,
            "binary64_epsilon": float(np.finfo(np.float64).eps),
        },
        "continuous_envelope_controls": continuous,
        "method_owned_pair_controls": methods,
        "atomic_failure_controls": atomic,
        "complete_owned_state": "dimensionless_u_p_q",
        "owned_row_policy": (
            "noncentre_rows_excluding_four_projector_owned_outer_rows"
        ),
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
    invalid_coefficients = False
    negative_radius = False
    try:
        binary64_cubic_absolute_envelope(np.ones(3))
    except (TypeError, ValueError):
        invalid_coefficients = True
    try:
        binary64_cubic_absolute_envelope(
            np.ones(4), np.asarray((0.0, 0.0, -1.0, 0.0))
        )
    except (TypeError, ValueError):
        negative_radius = True
    controls = {
        "invalid_cubic_shape_rejected": invalid_coefficients,
        "negative_coefficient_radius_rejected": negative_radius,
        "checkpoint_commit_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["immutable_lineage"].__setitem__(
                "checkpoint_commit", "0" * 40
            ),
        ),
        "source_result_hash_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["immutable_lineage"].__setitem__(
                "source_result_sha256", "0" * 64
            ),
        ),
        "coarse_commit_promotion_rejected": _expect_config_rejection(
            config,
            lambda value: value["runtime_contract"].__setitem__(
                "coarse_path_is_evidence_only", False
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
        "PROTO14_promotion_rejected": _expect_config_rejection(
            config,
            lambda value: value["claims"].__setitem__("PROTO14_frozen", True),
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
        "PDE_history_promotion_rejected": _expect_config_rejection(
            config,
            lambda value: value["limitation_audit"].__setitem__(
                "declared_cubic_is_exact_PDE_history", True
            ),
        ),
        "historical_reclassification_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["limitation_audit"].__setitem__(
                "historical_CAL10_result_reclassified", True
            ),
        ),
    }
    if set(controls.values()) != {True}:
        raise ValueError(f"IMP1 mutation control failed: {controls}")
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
            "implemented_and_synthetically_validated_TDG5_atomic_runtime_pair_"
            "threshold_design_only_authorized"
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
            "runtime_contract": deepcopy(config["runtime_contract"]),
            "continuous_envelope_contract": deepcopy(
                config["continuous_envelope"]
            ),
            "synthetic_validation": synthetic,
            "limitation_audit": deepcopy(config["limitation_audit"]),
            "mutation_controls": mutations,
            "proof_contract": dict(config["proof_contract"]),
            "successor_boundary": dict(config["successor_boundary"]),
            "claim_boundary": {
                "actual_terminal_history_arrays_consumed": False,
                "campaign_checkpoint_loaded": False,
                "campaign_checkpoint_resumed": False,
                "production_state_advanced": False,
                "synthetic_states_advanced": True,
                "historical_PROTO13_stop_preserved": True,
                "historical_CAL10_result_reclassified": False,
                "declared_cubic_is_exact_PDE_history": False,
                "trajectory_asymptotic_regime_proved": False,
                "rigorous_global_PDE_error_bound_proved": False,
                "runtime_refinement_pair_implemented": True,
                "runtime_refinement_pair_synthetic_validation_passed": True,
                "runtime_threshold_design_authorized": True,
                "runtime_thresholds_frozen": False,
                "replacement_temporal_admission_defined": False,
                "production_compositor_integrated": False,
                "new_production_trajectory_authorized": False,
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
    synthetic = payload.get("synthetic_validation", {})
    if (
        value.get("schema_version") != 1
        or value.get("artifact_id") != ARTIFACT_ID
        or value.get("project_version") != PROJECT_VERSION
        or value.get("classification")
        != "implemented_and_synthetically_validated_TDG5_atomic_runtime_pair_threshold_design_only_authorized"
        or value.get("generated_by") != "scripts/reproduce_fgc_tdg5_imp1.py"
        or value.get("gate_status") != EXPECTED_CLAIMS
        or value.get("nonclaims")
        != {
            key: False
            for key, item in sorted(EXPECTED_CLAIMS.items())
            if item is False
        }
        or value.get("source_config_sha256")
        != {_rel(DEFAULT_CONFIG): _sha(DEFAULT_CONFIG)}
        or payload.get("runtime_contract") != config["runtime_contract"]
        or payload.get("continuous_envelope_contract")
        != config["continuous_envelope"]
        or payload.get("limitation_audit") != config["limitation_audit"]
        or payload.get("successor_boundary") != config["successor_boundary"]
        or synthetic.get("all_synthetic_controls_pass") is not True
        or synthetic.get("synthetic_states_only") is not True
        or boundary.get("runtime_refinement_pair_implemented") is not True
        or boundary.get("runtime_refinement_pair_synthetic_validation_passed")
        is not True
        or boundary.get("runtime_threshold_design_authorized") is not True
        or boundary.get("runtime_thresholds_frozen") is not False
        or boundary.get("replacement_temporal_admission_defined") is not False
        or boundary.get("production_compositor_integrated") is not False
        or boundary.get("new_production_trajectory_authorized") is not False
        or boundary.get("PROTO14_frozen") is not False
        or boundary.get("GR0_case_eligible") is not False
        or boundary.get("candidate_branch_opened") is not False
        or boundary.get("physical_question_answered") is not False
        or boundary.get("actual_terminal_history_arrays_consumed") is not False
        or boundary.get("production_state_advanced") is not False
        or boundary.get("historical_CAL10_result_reclassified") is not False
    ):
        raise ValueError("IMP1 stored result boundary differs")
    implementation = value.get("implementation_sha256", {})
    if (
        not isinstance(implementation, Mapping)
        or set(implementation) != {_rel(path) for path in IMPLEMENTATION}
        or any(
            _sha(REPOSITORY / relative) != expected
            for relative, expected in implementation.items()
        )
    ):
        raise ValueError("IMP1 implementation hash ledger differs")
    document = REPOSITORY / str(value.get("derivation_document", ""))
    if (
        value.get("derivation_document") != "docs/fgc-tdg5-imp1.md"
        or not document.is_file()
        or _sha(document) != value.get("derivation_document_sha256")
    ):
        raise ValueError("IMP1 derivation document differs")
    if (
        not isinstance(payload.get("mutation_controls"), Mapping)
        or set(payload["mutation_controls"].values()) != {True}
        or not isinstance(payload.get("proof_contract"), Mapping)
        or set(payload["proof_contract"].values()) != {True}
    ):
        raise ValueError("IMP1 proof or mutation ledger differs")


def verify_canonical(
    config_path: Path = DEFAULT_CONFIG,
    output_path: Path = DEFAULT_OUTPUT,
) -> dict[str, Any]:
    config = load_config(config_path)
    raw = output_path.read_bytes()
    stored = json.loads(raw)
    if not isinstance(stored, dict) or raw != _canonical_bytes(stored):
        raise ValueError("IMP1 stored result is not canonical JSON")
    _validate_stored_record(stored, config)
    observed = record(config_path)
    if stored != observed:
        raise ValueError("IMP1 stored result differs from reproduction")
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
            f"PASS {ARTIFACT_ID}: runtime_pair="
            f"{status['runtime_refinement_pair_implemented']} synthetic="
            f"{status['runtime_refinement_pair_synthetic_validation_passed']} "
            f"thresholds={status['runtime_thresholds_frozen']} "
            f"PROTO14={status['PROTO14_frozen']}"
        )
    else:
        print(_canonical(record(args.config)), end="")


if __name__ == "__main__":
    main()
