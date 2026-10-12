#!/usr/bin/env python3
"""Reproduce the FGC-1-TDG6-PREF21 independent binder."""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import os
import subprocess
import sys
import tempfile
import tomllib
from pathlib import Path
from typing import Any, Mapping


REPOSITORY = Path(__file__).resolve().parents[1]
SOURCE = REPOSITORY / "src"
if str(SOURCE) not in sys.path:
    sys.path.insert(0, str(SOURCE))

from recursive_horizons.fgc.evolution.tdg6_temporal_admission_theorem import (  # noqa: E402
    TDG6_BINDER_CHANNELS,
    direct_restriction_matrix,
    runtime_extension_certificate,
    tdg6_independent_binder_certificate,
)


ARTIFACT_ID = "FGC-1-TDG6-PREF21"
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-tdg6-pref21.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-tdg6-pref21.json"
DERIVATION_DOCUMENT = REPOSITORY / "docs/fgc-tdg6-pref21.md"
BINDER_MODULE = (
    REPOSITORY
    / "src/recursive_horizons/fgc/evolution/tdg6_temporal_admission_theorem.py"
)


TOP_LEVEL_PATHS = {
    "freeze_result": "results/fgc-1-tdg6-frz1.json",
    "freeze_config": "configs/fgc/fgc-1-tdg6-frz1.toml",
    "freeze_reproducer": "scripts/reproduce_fgc_tdg6_frz1.py",
    "freeze_module": "src/recursive_horizons/fgc/evolution/tdg6_temporal_admission_design.py",
    "freeze_document": "docs/fgc-tdg6-frz1.md",
    "TDG5_runtime": "src/recursive_horizons/fgc/evolution/tdg5_stage_complete_refinement_runtime.py",
    "PROTO7_runtime": "src/recursive_horizons/fgc/evolution/proto7_runtime.py",
    "PROTO7_runner": "scripts/run_fgc_gr0_calibration_v7.py",
    "PROTO13_runner": "scripts/run_fgc_gr0_calibration_v13.py",
    "binder_module": "src/recursive_horizons/fgc/evolution/tdg6_temporal_admission_theorem.py",
}

EXPECTED_SCOPE = {
    "target": "FGC-2-SF1-TDG6",
    "role": "post_freeze_independent_temporal_admission_and_runtime_feasibility_binder",
    "calibration_branch": "GR-0",
    "freeze_compact_result_read": True,
    "actual_PROTO13_history_arrays_consumed": False,
    "campaign_checkpoint_loaded": False,
    "campaign_checkpoint_resumed": False,
    "production_state_advanced": False,
    "synthetic_exact_and_source_shape_controls_only": True,
    "SGBL_trajectory_read": False,
    "FGCQR_trajectory_read": False,
    "physical_or_candidate_question_answered": False,
}

EXPECTED_LINEAGE = {
    "freeze_commit": "dff13cfc55e93124aa1fec3e350eae12dd0c6213",
    "freeze_config_sha256": "407143881d2b15ca1837da8f2b7830d5bcaba74e6a7283b7519eec6262303349",
    "freeze_result_sha256": "3723532961c9d10bce9df2b96e0068437652d657a50d1d60036360c7acb417eb",
    "freeze_reproducer_sha256": "274a8e0157021a5ce496b385f101337ff8f6457c88b7d83ab0e83e00c2535a63",
    "freeze_module_sha256": "87d9a90e272a89479690c0d432f2d5cd979d88841615aaaeb9b49a2c0cca651e",
    "freeze_document_sha256": "c97e224956a150a750ae83391c2e2bb94cef180c92e342210ff7fc8c3a05c422",
    "TDG5_runtime_sha256": "75fc11943f459ca7549f75a812fe961e55ae04bbdb13450c1fdd509b5b53a080",
    "PROTO7_runtime_sha256": "ffa75767c1fdbae42a09c8ef50276183861fc9255ad0061e5a1ee3673f142734",
    "PROTO7_runner_sha256": "04324e7066dc5a3c4ea3f6200612cc6c1a51fae5adba0c28ecb10ebec2f5d2e7",
    "PROTO13_runner_sha256": "bdc31c8ab4205ff86a134e9744aaaf8d085344aa1ba6a19cebb3fdbad94286d7",
    "freeze_commit_must_be_ancestor_of_HEAD": True,
    "all_bound_blobs_must_match_freeze_commit_and_worktree": True,
    "freeze_result_must_be_canonical_and_preserve_unimplemented_boundary": True,
    "raw_campaign_bundle_must_not_be_loaded": True,
}

EXPECTED_QUARTER = {
    "declared_polynomial_class": "P3_on_each_accepted_normalized_time_interval",
    "direct_formula": "b_k=sum_n>=k[a_n*C(n,k)*j^(n-k)/4^n]",
    "quarter_interval_indices": [0, 1, 2, 3],
    "quarter_restriction_determinant": "1/4096",
    "all_four_maps_must_be_injective": True,
    "direct_quarter_maps_must_equal_two_composed_half_maps": True,
    "complete_basis_and_probe_controls_required": True,
    "continuous_D01_subinterval_count": 2,
    "continuous_D12_subinterval_count": 4,
    "finite_point_sampling_used_as_continuum_proof": False,
}

EXPECTED_INTERVAL = {
    "minimum_observed_order": "3/2",
    "squared_multiplier": 8,
    "conservative_test": "8*U12^2<=L01^2",
    "universal_implication_chain": "8*d12^2<=8*U12^2<=L01^2<=d01^2",
    "nonnegative_outward_intervals_required": True,
    "equality_passes": True,
    "one_rational_unit_below_fails": True,
    "exact_zero_passes_without_order_claim": True,
    "enclosure_dominated_passes_without_order_claim": True,
    "enclosure_dominated_fine_upper_must_not_exceed_outer_upper": True,
    "resolved_failure_is_retryable_numerical_admission": True,
    "inconclusive_is_retryable_numerical_admission": True,
    "absolute_magnitude_threshold_used": False,
    "physical_signal_normalization_used": False,
}

EXPECTED_CHANNEL = {
    "channel_order": list(TDG6_BINDER_CHANNELS),
    "channel_count": 18,
    "complete_admission_is_all_of": True,
    "one_failed_channel_vetoes_commit": True,
    "per_channel_debit": "D12_certified_upper",
    "debits_accumulate_componentwise_by_nonnegative_addition": True,
    "cancellation_permitted": False,
    "Richardson_division_used_for_admission_debit": False,
    "debit_is_rigorous_global_PDE_error_bound": False,
    "future_constraint_to_Raychaudhuri_stability_map_required": True,
}

EXPECTED_RUNTIME = {
    "TDG5_prepared_pair_has_state_monitor_and_causal_shadow_boundary": True,
    "TDG5_prepare_uses_two_transaction_clones_and_three_shadow_attempts": True,
    "TDG5_commit_revalidates_time_state_step_serial_monitor_and_causal_boundary": True,
    "TDG5_commit_adopts_only_fine_monitor_and_causal_state": True,
    "TDG5_shadow_attempt_currently_uses_no_tracer_preaccept": True,
    "PROTO7_exposes_tracer_preview_commit_and_rollback_inputs": True,
    "PROTO13_checkpoint_exposes_state_monitor_causal_and_tracer_extension_surface": True,
    "TDG6_temporal_debit_and_retry_fields_currently_absent": True,
    "production_compositor_extension_is_source_shape_feasible": True,
    "production_compositor_implemented": False,
    "production_trajectory_authorized": False,
}

EXPECTED_AUTHORIZATION = {
    "authorized_object": "implementation_and_synthetic_qualification_of_the_frozen_TDG6_three_level_GR0_production_compositor",
    "same_bitwise_initial_state_for_all_three_levels_required": True,
    "same_method_source_projector_spatial_operator_and_macro_interval_required": True,
    "coarse_and_medium_paths_are_evidence_only": True,
    "only_complete_four_quarter_step_path_is_committable": True,
    "fine_state_monitor_causal_tracer_and_debit_ledgers_commit_atomically": True,
    "all_shadow_paths_require_independent_monitor_causal_and_tracer_state": True,
    "unchanged_source_health_CFL_boundary_scale_and_transaction_gates_required": True,
    "complete_18_channel_D01_and_D12_interval_evidence_required": True,
    "temporal_retry_counter_separate_from_source_and_CFL_counters": True,
    "rejection_evidence_must_be_durable_before_retry": True,
    "checkpoint_must_round_trip_debit_retry_and_rejection_ledgers": True,
    "unexpected_exception_is_invalid_not_physical": True,
    "threshold_or_protocol_mutation_authorized": False,
    "production_trajectory_authorized": False,
    "historical_PROTO13_stop_preserved": True,
    "historical_CAL10_result_reclassified": False,
}

EXPECTED_LIMITATIONS = {
    "piecewise_cubic_extension_is_not_exact_PDE_history": True,
    "three_level_order_pass_does_not_prove_trajectory_asymptotics": True,
    "local_debit_does_not_prove_a_global_error_bound": True,
    "source_shape_feasibility_is_not_runtime_implementation": True,
    "tracer_and_checkpoint_extensions_remain_to_be_implemented": True,
    "spatial_ladders_and_cross_method_agreement_remain_required": True,
    "complete_DEF1_error_ledger_remains_required": True,
    "no_GR0_SGBL_or_FGCQR_outcome_is_classified": True,
}

EXPECTED_PROOF = {
    "freeze_commit_and_all_bound_blobs_must_match": True,
    "freeze_result_must_be_canonical_and_preserve_the_unimplemented_boundary": True,
    "quarter_maps_must_be_derived_directly_without_importing_the_freeze_design_module": True,
    "quarter_maps_must_match_two_independently_composed_half_maps": True,
    "quarter_determinants_and_basis_controls_must_be_exact": True,
    "interval_implication_and_squared_multiplier_must_be_rederived_exactly": True,
    "channel_all_of_and_nonnegative_debit_semantics_must_be_independently_controlled": True,
    "runtime_feasibility_must_be_checked_against_immutable_source_syntax": True,
    "tracer_and_checkpoint_missing-production-fields_must_remain_explicit": True,
    "freeze_crosscheck_may_not_replace_independent_derivation": True,
    "no_history_checkpoint_state_or_candidate_data_may_be_loaded": True,
    "lineage_theorem_runtime_authorization_nonclaim_and_claim_mutations_must_fail_closed": True,
    "canonical_hash_bound_result_required": True,
}

EXPECTED_SUCCESSOR = {
    "TDG5_runtime_refinement_pair_implemented": True,
    "TDG6_threshold_and_admission_design_frozen": True,
    "runtime_thresholds_frozen": True,
    "replacement_temporal_admission_defined": True,
    "TDG6_independent_binder_completed": True,
    "production_compositor_implementation_authorized": True,
    "production_compositor_implemented": False,
    "PROTO14_frozen": False,
    "fresh_GR0_dynamic_calibration_completed": False,
    "GR0_case_eligible": False,
    "SGBL_execution_authorized": False,
    "FGCQR_holdout_execution_authorized": False,
    "DEF1_execution_authorized": False,
}

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
    "production_compositor_implemented": False,
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
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _bytes_sha(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def _canonical(value: Mapping[str, Any]) -> str:
    return json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n"


def _canonical_bytes(value: Mapping[str, Any]) -> bytes:
    return _canonical(value).encode("utf-8")


def _atomic_write(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def _git(*args: str, check: bool = True) -> subprocess.CompletedProcess[bytes]:
    return subprocess.run(
        ["git", *args],
        cwd=REPOSITORY,
        check=check,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )


def _strict_equal(name: str, observed: object, expected: object) -> None:
    if observed != expected:
        raise ValueError(f"{name} differs from the frozen PREF21 contract")


def validate_config_data(config: Mapping[str, Any]) -> None:
    top_keys = {
        "schema_version",
        "artifact_id",
        "project_version",
        "metric_signature",
        "riemann_convention",
        *TOP_LEVEL_PATHS,
        "scope",
        "immutable_lineage",
        "quarter_restriction",
        "interval_order",
        "channel_debit",
        "runtime_feasibility",
        "implementation_authorization",
        "limitation_audit",
        "proof_contract",
        "successor_boundary",
        "claims",
    }
    _strict_equal("top-level keys", set(config), top_keys)
    _strict_equal("schema_version", config.get("schema_version"), 1)
    _strict_equal("artifact_id", config.get("artifact_id"), ARTIFACT_ID)
    _strict_equal("project_version", config.get("project_version"), "0.11.0")
    _strict_equal("metric_signature", config.get("metric_signature"), "-+++")
    _strict_equal(
        "riemann_convention",
        config.get("riemann_convention"),
        "plus_partial_mu_gamma_nu",
    )
    for key, expected in TOP_LEVEL_PATHS.items():
        _strict_equal(key, config.get(key), expected)
    for name, expected in (
        ("scope", EXPECTED_SCOPE),
        ("immutable_lineage", EXPECTED_LINEAGE),
        ("quarter_restriction", EXPECTED_QUARTER),
        ("interval_order", EXPECTED_INTERVAL),
        ("channel_debit", EXPECTED_CHANNEL),
        ("runtime_feasibility", EXPECTED_RUNTIME),
        ("implementation_authorization", EXPECTED_AUTHORIZATION),
        ("limitation_audit", EXPECTED_LIMITATIONS),
        ("proof_contract", EXPECTED_PROOF),
        ("successor_boundary", EXPECTED_SUCCESSOR),
        ("claims", EXPECTED_CLAIMS),
    ):
        _strict_equal(name, config.get(name), expected)


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    value = tomllib.loads(path.read_text(encoding="utf-8"))
    validate_config_data(value)
    return value


def _verify_lineage(config: Mapping[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    lineage = config["immutable_lineage"]
    commit = lineage["freeze_commit"]
    if len(commit) != 40 or any(character not in "0123456789abcdef" for character in commit):
        raise ValueError("freeze_commit must be a full lowercase SHA-1")
    _git("cat-file", "-e", f"{commit}^{{commit}}")
    if _git("merge-base", "--is-ancestor", commit, "HEAD", check=False).returncode != 0:
        raise ValueError("freeze commit is not an ancestor of HEAD")

    binding = {
        "freeze_config": "freeze_config_sha256",
        "freeze_result": "freeze_result_sha256",
        "freeze_reproducer": "freeze_reproducer_sha256",
        "freeze_module": "freeze_module_sha256",
        "freeze_document": "freeze_document_sha256",
        "TDG5_runtime": "TDG5_runtime_sha256",
        "PROTO7_runtime": "PROTO7_runtime_sha256",
        "PROTO7_runner": "PROTO7_runner_sha256",
        "PROTO13_runner": "PROTO13_runner_sha256",
    }
    verified: dict[str, str] = {}
    for path_key, hash_key in binding.items():
        relative = str(config[path_key])
        expected = str(lineage[hash_key])
        blob = _git("show", f"{commit}:{relative}").stdout
        worktree = (REPOSITORY / relative).read_bytes()
        if _bytes_sha(blob) != expected or _bytes_sha(worktree) != expected or blob != worktree:
            raise ValueError(f"immutable PREF21 input differs: {relative}")
        verified[relative] = expected

    freeze_path = REPOSITORY / str(config["freeze_result"])
    freeze_raw = freeze_path.read_bytes()
    freeze = json.loads(freeze_raw)
    if not isinstance(freeze, dict) or freeze_raw != _canonical_bytes(freeze):
        raise ValueError("TDG6 freeze result is not canonical JSON")
    if (
        freeze.get("artifact_id") != "FGC-1-TDG6-FRZ1"
        or freeze.get("gate_status") != "PASS_THRESHOLD_AND_ADMISSION_DESIGN_FROZEN"
    ):
        raise ValueError("TDG6 freeze identity differs")
    payload = freeze.get("artifact_payload")
    if not isinstance(payload, Mapping):
        raise ValueError("TDG6 freeze payload is missing")
    boundary = payload.get("claim_boundary")
    successor = payload.get("successor_boundary")
    if not isinstance(boundary, Mapping) or not isinstance(successor, Mapping):
        raise ValueError("TDG6 freeze boundary is missing")
    if (
        boundary.get("TDG6_independent_binder_completed") is not False
        or boundary.get("production_compositor_implementation_authorized") is not False
        or boundary.get("production_compositor_implemented") is not False
        or boundary.get("actual_PROTO13_history_arrays_consumed") is not False
        or boundary.get("historical_CAL10_result_reclassified") is not False
        or successor.get("PROTO14_frozen") is not False
        or successor.get("GR0_case_eligible") is not False
        or successor.get("FGCQR_holdout_execution_authorized") is not False
    ):
        raise ValueError("TDG6 freeze no longer preserves its unimplemented boundary")
    return freeze, {
        "freeze_commit": commit,
        "freeze_commit_is_ancestor_of_HEAD": True,
        "bound_source_sha256": verified,
        "freeze_result_verified_canonical_and_unimplemented": True,
        "raw_campaign_bundle_loaded": False,
    }


def _freeze_crosscheck(
    freeze: Mapping[str, Any], binder: Mapping[str, Any]
) -> dict[str, bool]:
    payload = freeze["artifact_payload"]
    route = payload["selected_route"]
    threshold = payload["threshold_contract"]
    error = payload["error_ledger"]
    quarter = binder["independent_quarter_restriction"]
    interval = binder["independent_interval_order_theorem"]
    channel = binder["independent_channel_debit_reduction"]
    runtime = binder["immutable_runtime_extension_audit"]
    controls = {
        "three_level_counts_match": route["level_proposal_counts"] == [1, 2, 4],
        "complete_channel_order_matches": route["complete_state_channel_order"]
        == list(TDG6_BINDER_CHANNELS),
        "quarter_maps_supply_the_frozen_fine_partition": quarter[
            "four_quarter_maps_are_exact_and_injective"
        ],
        "minimum_order_matches": threshold[
            "minimum_observed_temporal_refinement_order"
        ]
        == interval["minimum_observed_order"],
        "squared_test_matches": (
            threshold["resolved_order_test"]
            == "8*D12_upper_squared<=D01_lower_squared"
            and interval["conservative_test"] == "8*U12^2<=L01^2"
        ),
        "equality_semantics_match": threshold["equality_at_three_halves_passes"]
        == interval["exact_boundary_passes"],
        "no_absolute_threshold_matches": threshold[
            "absolute_state_error_threshold_defined"
        ]
        is False
        and interval["no_absolute_magnitude_threshold_is_used"],
        "debit_semantics_match": error["per_channel_finest_pair_debit"]
        == "D12_certified_upper"
        and channel["finest_upper_bound_is_retained_per_channel"],
        "noncancellation_matches": error[
            "accepted_macro_step_debits_accumulate_by_nonnegative_addition_without_cancellation"
        ]
        and channel[
            "accepted_debits_accumulate_by_componentwise_nonnegative_addition"
        ],
        "production_remains_unimplemented": runtime[
            "production_compositor_implemented"
        ]
        is False,
    }
    if set(controls.values()) != {True}:
        raise ValueError(f"independent binder differs from TDG6 freeze: {controls}")
    return controls


def _expect_config_rejection(config: Mapping[str, Any], mutate: Any) -> bool:
    attacked = copy.deepcopy(config)
    mutate(attacked)
    try:
        validate_config_data(attacked)
    except (TypeError, ValueError):
        return True
    return False


def _mutation_controls(
    config: Mapping[str, Any],
    *,
    tdg5_source: str,
    proto7_source: str,
    proto13_source: str,
) -> dict[str, bool]:
    runtime_preview = runtime_extension_certificate(
        tdg5_runtime_source=tdg5_source,
        proto7_runner_source=proto7_source.replace("preview_advance", "preview_removed", 1),
        proto13_runner_source=proto13_source,
    )
    runtime_commit = runtime_extension_certificate(
        tdg5_runtime_source=tdg5_source.replace(
            "prepared.fine_monitor_state", "prepared.original_monitor_state", 1
        ),
        proto7_runner_source=proto7_source,
        proto13_runner_source=proto13_source,
    )
    runtime_checkpoint = runtime_extension_certificate(
        tdg5_runtime_source=tdg5_source,
        proto7_runner_source=proto7_source,
        proto13_runner_source=proto13_source.replace(
            "_write_checkpoint", "_write_removed"
        ),
    )
    controls = {
        "freeze_commit_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["immutable_lineage"].__setitem__(
                "freeze_commit", "0" * 40
            ),
        ),
        "freeze_result_hash_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["immutable_lineage"].__setitem__(
                "freeze_result_sha256", "0" * 64
            ),
        ),
        "quarter_determinant_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["quarter_restriction"].__setitem__(
                "quarter_restriction_determinant", "1/1024"
            ),
        ),
        "quarter_composition_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["quarter_restriction"].__setitem__(
                "direct_quarter_maps_must_equal_two_composed_half_maps", False
            ),
        ),
        "squared_multiplier_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["interval_order"].__setitem__(
                "squared_multiplier", 7
            ),
        ),
        "equality_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["interval_order"].__setitem__(
                "equality_passes", False
            ),
        ),
        "physical_normalization_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["interval_order"].__setitem__(
                "physical_signal_normalization_used", True
            ),
        ),
        "channel_drop_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["channel_debit"].__setitem__(
                "channel_order", value["channel_debit"]["channel_order"][:-1]
            ),
        ),
        "all_of_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["channel_debit"].__setitem__(
                "complete_admission_is_all_of", False
            ),
        ),
        "fine_commit_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["implementation_authorization"].__setitem__(
                "only_complete_four_quarter_step_path_is_committable", False
            ),
        ),
        "tracer_obligation_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["implementation_authorization"].__setitem__(
                "all_shadow_paths_require_independent_monitor_causal_and_tracer_state",
                False,
            ),
        ),
        "checkpoint_obligation_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["implementation_authorization"].__setitem__(
                "checkpoint_must_round_trip_debit_retry_and_rejection_ledgers", False
            ),
        ),
        "trajectory_authorization_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["implementation_authorization"].__setitem__(
                "production_trajectory_authorized", True
            ),
        ),
        "historical_reclassification_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["implementation_authorization"].__setitem__(
                "historical_CAL10_result_reclassified", True
            ),
        ),
        "candidate_promotion_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["claims"].__setitem__(
                "FGCQR_holdout_execution_authorized", True
            ),
        ),
        "invalid_quarter_index_rejected": _invalid_quarter_index_rejected(),
        "tracer_preview_source_mutation_rejected": not runtime_preview[
            "production_compositor_extension_is_source_shape_feasible"
        ],
        "fine_monitor_source_mutation_rejected": not runtime_commit[
            "production_compositor_extension_is_source_shape_feasible"
        ],
        "checkpoint_writer_source_mutation_rejected": not runtime_checkpoint[
            "production_compositor_extension_is_source_shape_feasible"
        ],
    }
    if set(controls.values()) != {True}:
        raise ValueError(f"PREF21 mutation control failed: {controls}")
    return controls


def _invalid_quarter_index_rejected() -> bool:
    try:
        direct_restriction_matrix(divisor=4, interval_index=4)
    except ValueError:
        return True
    return False


def record(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    config = load_config(config_path)
    freeze, lineage = _verify_lineage(config)
    tdg5_source = (REPOSITORY / str(config["TDG5_runtime"])).read_text(
        encoding="utf-8"
    )
    proto7_source = (REPOSITORY / str(config["PROTO7_runner"])).read_text(
        encoding="utf-8"
    )
    proto13_source = (REPOSITORY / str(config["PROTO13_runner"])).read_text(
        encoding="utf-8"
    )
    binder = tdg6_independent_binder_certificate(
        tdg5_runtime_source=tdg5_source,
        proto7_runner_source=proto7_source,
        proto13_runner_source=proto13_source,
    )
    if (
        binder["all_independent_binder_controls_pass"] is not True
        or binder["production_compositor_implementation_authorized"] is not True
        or binder["production_compositor_implemented"] is not False
    ):
        raise ValueError("PREF21 independent binder did not authorize implementation")
    crosscheck = _freeze_crosscheck(freeze, binder)
    mutations = _mutation_controls(
        config,
        tdg5_source=tdg5_source,
        proto7_source=proto7_source,
        proto13_source=proto13_source,
    )
    return {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": config["project_version"],
        "metric_signature": config["metric_signature"],
        "riemann_convention": config["riemann_convention"],
        "classification": "independent_temporal_admission_binder_and_runtime_implementation_authorization",
        "gate_status": "PASS_PRODUCTION_COMPOSITOR_IMPLEMENTATION_AUTHORIZED",
        "source_config_sha256": _sha(config_path),
        "derivation_document": _rel(DERIVATION_DOCUMENT),
        "derivation_document_sha256": _sha(DERIVATION_DOCUMENT),
        "generated_by": _rel(Path(__file__)),
        "implementation_sha256": {
            _rel(config_path): _sha(config_path),
            _rel(DERIVATION_DOCUMENT): _sha(DERIVATION_DOCUMENT),
            _rel(Path(__file__)): _sha(Path(__file__)),
            _rel(BINDER_MODULE): _sha(BINDER_MODULE),
        },
        "artifact_payload": {
            "immutable_lineage": lineage,
            "independent_binder": binder,
            "freeze_crosscheck": crosscheck,
            "mutation_controls": mutations,
            "implementation_authorization": dict(config["implementation_authorization"]),
            "limitation_audit": dict(config["limitation_audit"]),
            "proof_contract": dict(config["proof_contract"]),
            "successor_boundary": dict(config["successor_boundary"]),
            "claim_boundary": dict(config["claims"]),
        },
        "nonclaims": [
            "PREF21 reads no PROTO13 history or checkpoint and advances no state.",
            "The quarter restriction and interval theorem apply to the declared numerical extension, not the exact PDE history.",
            "Source-shape feasibility is not a production compositor implementation or runtime qualification.",
            "The temporal debit is not a rigorous global PDE error bound and still requires the future observable stability map.",
            "CAL10 remains historical and unreclassified; no GR-0 case is eligible.",
            "No SGB-L, FGC-QR, collapse, defocusing, transition, retained-EFT, dark-sector, child-domain, or varying-local-c claim follows.",
        ],
    }


def _validate_stored_record(
    value: Mapping[str, Any], config_path: Path, output_path: Path
) -> None:
    expected = record(config_path)
    if value != expected:
        raise ValueError("stored PREF21 record differs from independent reproduction")
    if value.get("artifact_id") != ARTIFACT_ID:
        raise ValueError("stored PREF21 artifact identity differs")
    if value.get("gate_status") != "PASS_PRODUCTION_COMPOSITOR_IMPLEMENTATION_AUTHORIZED":
        raise ValueError("stored PREF21 gate status differs")
    implementation = value.get("implementation_sha256")
    if not isinstance(implementation, Mapping) or set(implementation) != {
        _rel(config_path),
        _rel(DERIVATION_DOCUMENT),
        _rel(Path(__file__)),
        _rel(BINDER_MODULE),
    }:
        raise ValueError("stored PREF21 implementation hash ledger differs")
    payload = value.get("artifact_payload")
    if not isinstance(payload, Mapping):
        raise ValueError("stored PREF21 payload is missing")
    binder = payload.get("independent_binder")
    mutations = payload.get("mutation_controls")
    if (
        not isinstance(binder, Mapping)
        or binder.get("all_independent_binder_controls_pass") is not True
        or binder.get("production_compositor_implementation_authorized") is not True
        or binder.get("production_compositor_implemented") is not False
        or not isinstance(mutations, Mapping)
        or set(mutations.values()) != {True}
    ):
        raise ValueError("stored PREF21 theorem or mutation ledger differs")
    if output_path != DEFAULT_OUTPUT and output_path.exists():
        # The caller owns alternate outputs; this check merely makes the
        # otherwise-unused argument part of the explicit validation surface.
        output_path.resolve()


def verify_canonical(
    config_path: Path = DEFAULT_CONFIG, output_path: Path = DEFAULT_OUTPUT
) -> dict[str, Any]:
    raw = output_path.read_bytes()
    stored = json.loads(raw)
    if not isinstance(stored, dict) or raw != _canonical_bytes(stored):
        raise ValueError("PREF21 stored result is not canonical JSON")
    _validate_stored_record(stored, config_path, output_path)
    return stored


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    if args.write and args.verify:
        parser.error("--write and --verify are mutually exclusive")
    config_path = args.config.resolve()
    output_path = args.output.resolve()
    if args.write:
        value = record(config_path)
        _atomic_write(output_path, _canonical_bytes(value))
        print(f"wrote {_rel(output_path)}")
        return
    if args.verify:
        value = verify_canonical(config_path, output_path)
        print(
            json.dumps(
                {
                    "artifact_id": value["artifact_id"],
                    "gate_status": value["gate_status"],
                    "verified": True,
                },
                sort_keys=True,
            )
        )
        return
    print(_canonical(record(config_path)), end="")


if __name__ == "__main__":
    main()
