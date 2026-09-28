#!/usr/bin/env python3
"""Reproduce the prospective FGC-1-TDG6-FRZ1 threshold-design freeze."""

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
from typing import Any, Callable, Mapping


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.tdg6_temporal_admission_design import (  # noqa: E402
    TDG6_COMPLETE_STATE_CHANNELS,
    TDG6_LEVEL_STEP_COUNTS,
    TDG6_MAXIMUM_TEMPORAL_RETRIES_PER_STEP,
    TDG6_MINIMUM_MACRO_STEP,
    TDG6_MINIMUM_OBSERVED_ORDER,
    TDG6_RETRY_FACTOR,
    tdg6_design_preflight,
)


ARTIFACT_ID = "FGC-1-TDG6-FRZ1"
PROJECT_VERSION = "0.11.0"
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-tdg6-frz1.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-tdg6-frz1.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-tdg6-frz1.md"
DESIGN_MODULE = (
    REPOSITORY
    / "src/recursive_horizons/fgc/evolution/tdg6_temporal_admission_design.py"
)

SOURCE_PATHS = {
    "source_config": REPOSITORY / "configs/fgc/fgc-1-tdg5-imp1.toml",
    "source_result": REPOSITORY / "results/fgc-1-tdg5-imp1.json",
    "source_reproducer": REPOSITORY / "scripts/reproduce_fgc_tdg5_imp1.py",
    "source_runtime_module": (
        REPOSITORY
        / "src/recursive_horizons/fgc/evolution/tdg5_stage_complete_refinement_runtime.py"
    ),
    "source_document": REPOSITORY / "docs/fgc-tdg5-imp1.md",
    "protocol_v13": REPOSITORY / "configs/fgc/fgc-2-sf1-protocol-v13.toml",
    "PROTO13_runner": REPOSITORY / "scripts/run_fgc_gr0_calibration_v13.py",
    "PROTO7_runtime": (
        REPOSITORY / "src/recursive_horizons/fgc/evolution/proto7_runtime.py"
    ),
    "PROTO7_runner": REPOSITORY / "scripts/run_fgc_gr0_calibration_v7.py",
}

EXPECTED_LINEAGE = {
    "checkpoint_commit": "986ceb1afb7fa7574c3cda7c623cc94df69471ce",
    "source_config_sha256": "96e31b0178c91dc86563cc745d16fa1b34c76f433829d410a772a3e567ea523a",
    "source_result_sha256": "cc95b04f5cea1c1f36443bea78d3851cba4e763a09551fb208f795cdbe232921",
    "source_reproducer_sha256": "f3291ff4e8508b491df8bbb1e3438a06c69d89200e0d19f6205a325ad6637694",
    "source_runtime_module_sha256": "75fc11943f459ca7549f75a812fe961e55ae04bbdb13450c1fdd509b5b53a080",
    "source_document_sha256": "60e42dde3766b3f2d071ac379f7adb6fefb1bcdb203c6273d12e9a6992a5f31f",
    "protocol_v13_sha256": "fc9398bad6ea4448f659188641e6da6bd7c2e295a443fe73bf73eb0308480692",
    "PROTO13_runner_sha256": "bdc31c8ab4205ff86a134e9744aaaf8d085344aa1ba6a19cebb3fdbad94286d7",
    "PROTO7_runtime_sha256": "ffa75767c1fdbae42a09c8ef50276183861fc9255ad0061e5a1ee3673f142734",
    "PROTO7_runner_sha256": "04324e7066dc5a3c4ea3f6200612cc6c1a51fae5adba0c28ecb10ebec2f5d2e7",
    "checkpoint_commit_must_be_ancestor_of_HEAD": True,
    "all_bound_blobs_must_match_checkpoint_commit_and_worktree": True,
    "source_result_must_be_canonical_and_authorize_threshold_design_only": True,
    "raw_campaign_bundle_must_not_be_loaded": True,
}

EXPECTED_CLAIMS = {
    "TDG4_sampling_identifiability_theorem_completed": True,
    "current_temporal_spectral_rule_retired_as_continuum_admission": True,
    "current_temporal_spectral_rule_retained_as_sampled_diagnostic": True,
    "TDG5_runtime_refinement_pair_implemented": True,
    "TDG5_runtime_refinement_pair_synthetic_validation_passed": True,
    "TDG6_threshold_and_admission_design_frozen": True,
    "runtime_thresholds_frozen": True,
    "replacement_temporal_admission_defined": True,
    "TDG6_independent_binder_completed": False,
    "production_compositor_implementation_authorized": False,
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
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _bytes_sha(value: bytes) -> str:
    return sha256(value).hexdigest()


def _canonical(value: Mapping[str, Any]) -> str:
    return json.dumps(
        value,
        indent=2,
        sort_keys=True,
        ensure_ascii=True,
        allow_nan=False,
    ) + "\n"


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


def _load_config(path: Path) -> dict[str, Any]:
    with path.open("rb") as handle:
        return tomllib.load(handle)


def _validate_source_result(path: Path) -> dict[str, Any]:
    payload = path.read_bytes()
    result = json.loads(payload)
    if _canonical(result).encode("utf-8") != payload:
        raise ValueError("TDG5-IMP1 source result is not canonical")
    if (
        result.get("artifact_id") != "FGC-1-TDG5-IMP1"
        or result.get("classification")
        != "implemented_and_synthetically_validated_TDG5_atomic_runtime_pair_threshold_design_only_authorized"
    ):
        raise ValueError("TDG5-IMP1 source result has the wrong boundary")
    status = result.get("gate_status", {})
    if (
        status.get("runtime_refinement_pair_implemented") is not True
        or status.get("runtime_refinement_pair_synthetic_validation_passed")
        is not True
        or status.get("runtime_threshold_design_authorized") is not True
        or status.get("runtime_thresholds_frozen") is not False
        or status.get("replacement_temporal_admission_defined") is not False
        or status.get("GR0_case_eligible") is not False
        or status.get("FGCQR_holdout_execution_authorized") is not False
    ):
        raise ValueError("TDG5-IMP1 gate status differs")
    claim = result.get("artifact_payload", {}).get("claim_boundary", {})
    if (
        claim.get("runtime_refinement_pair_implemented") is not True
        or claim.get("runtime_refinement_pair_synthetic_validation_passed")
        is not True
        or claim.get("runtime_threshold_design_authorized") is not True
        or claim.get("runtime_thresholds_frozen") is not False
        or claim.get("replacement_temporal_admission_defined") is not False
        or claim.get("production_compositor_integrated") is not False
        or claim.get("GR0_case_eligible") is not False
        or claim.get("candidate_branch_opened") is not False
        or claim.get("physical_question_answered") is not False
    ):
        raise ValueError("TDG5-IMP1 source result was promoted beyond its scope")
    return result


def _validate_config(value: Mapping[str, Any]) -> None:
    expected_top = {
        "schema_version",
        "artifact_id",
        "project_version",
        "metric_signature",
        "riemann_convention",
        *SOURCE_PATHS,
        "design_module",
        "scope",
        "immutable_lineage",
        "selected_route",
        "thresholds",
        "error_ledger",
        "production_compositor",
        "exact_controls",
        "alternatives",
        "proof_contract",
        "successor_boundary",
        "claims",
    }
    _strict_keys("TDG6 config", value, expected_top)
    if (
        value["schema_version"] != 1
        or value["artifact_id"] != ARTIFACT_ID
        or value["project_version"] != PROJECT_VERSION
        or value["metric_signature"] != "-+++"
        or value["riemann_convention"] != "plus_partial_mu_gamma_nu"
        or value["design_module"] != _rel(DESIGN_MODULE)
    ):
        raise ValueError("TDG6 identity differs")
    for key, path in SOURCE_PATHS.items():
        if value[key] != _rel(path):
            raise ValueError(f"TDG6 source path differs: {key}")

    scope = value["scope"]
    if (
        scope.get("target") != "FGC-2-SF1-TDG6"
        or scope.get("calibration_branch") != "GR-0"
        or scope.get("source_compact_result_read") is not True
        or scope.get("actual_PROTO13_history_arrays_consumed") is not False
        or scope.get("campaign_checkpoint_loaded") is not False
        or scope.get("campaign_checkpoint_resumed") is not False
        or scope.get("production_state_advanced") is not False
        or scope.get("synthetic_or_exact_design_controls_only") is not True
        or scope.get("SGBL_trajectory_read") is not False
        or scope.get("FGCQR_trajectory_read") is not False
        or scope.get("physical_or_candidate_question_answered") is not False
    ):
        raise ValueError("TDG6 scope differs")

    if value["immutable_lineage"] != EXPECTED_LINEAGE:
        raise ValueError("TDG6 immutable lineage differs")

    route = value["selected_route"]
    if (
        route.get("name")
        != "three_level_same_grid_refinement_with_order_gate_and_public_debit"
        or tuple(route.get("level_proposal_counts", ())) != TDG6_LEVEL_STEP_COUNTS
        or tuple(route.get("complete_state_channel_order", ()))
        != TDG6_COMPLETE_STATE_CHANNELS
        or route.get("all_levels_start_from_one_bitwise_identical_accepted_state")
        is not True
        or route.get("outer_and_medium_levels_are_evidence_only") is not True
        or route.get("only_the_complete_four_quarter_step_fine_level_is_committable")
        is not True
        or route.get("all_four_fine_steps_commit_atomically_or_none_commits")
        is not True
        or route.get("old_64_sample_temporal_spectrum_is_admission_veto")
        is not False
        or route.get("old_64_sample_temporal_spectrum_remains_public_diagnostic")
        is not True
    ):
        raise ValueError("TDG6 selected route differs")

    thresholds = value["thresholds"]
    if (
        thresholds.get("minimum_observed_temporal_refinement_order")
        != str(TDG6_MINIMUM_OBSERVED_ORDER)
        or thresholds.get("resolved_order_test")
        != "8*D12_upper_squared<=D01_lower_squared"
        or thresholds.get("equality_at_three_halves_passes") is not True
        or thresholds.get("enclosure_dominated_class_makes_no_order_claim")
        is not True
        or thresholds.get("temporal_retry_factor") != str(TDG6_RETRY_FACTOR)
        or thresholds.get("maximum_temporal_retries_per_macro_step")
        != TDG6_MAXIMUM_TEMPORAL_RETRIES_PER_STEP
        or thresholds.get("minimum_macro_step") != str(TDG6_MINIMUM_MACRO_STEP)
        or thresholds.get("absolute_state_error_threshold_required") is not False
        or thresholds.get("absolute_state_error_threshold_defined") is not False
        or thresholds.get("observed_trapped_or_Raychaudhuri_margin_used_to_choose_threshold")
        is not False
        or thresholds.get("historical_CAL10_values_used_to_choose_threshold")
        is not False
    ):
        raise ValueError("TDG6 thresholds differ")

    ledger = value["error_ledger"]
    if (
        ledger.get("per_channel_finest_pair_debit") != "D12_certified_upper"
        or ledger.get("debit_divided_by_formal_method_order_denominator_for_admission")
        is not False
        or ledger.get("accepted_macro_step_debits_accumulate_by_nonnegative_addition_without_cancellation")
        is not True
        or ledger.get("common_event_temporal_debit_vector_has_18_channels")
        is not True
        or ledger.get("future_trapped_and_DEF1_sign_margins_must_pay_the_temporal_debit")
        is not True
        or ledger.get("future_constraint_to_Raychaudhuri_stability_map_remains_required")
        is not True
        or ledger.get("debit_is_rigorous_global_PDE_error_bound") is not False
        or ledger.get("debit_is_proof_of_trajectory_asymptotics") is not False
    ):
        raise ValueError("TDG6 error ledger differs")

    compositor = value["production_compositor"]
    if (
        compositor.get("macro_step_is_selected_only_from_inherited_CFL_and_target_event_distance_before_refinement")
        is not True
        or compositor.get("coarse_medium_and_fine_transactions_have_independent_monitor_causal_and_tracer_shadows")
        is not True
        or compositor.get("only_fine_tracer_monitor_causal_state_and_fields_commit")
        is not True
        or compositor.get("commit_revalidates_state_time_step_serial_monitor_causal_tracer_and_debit_ledgers")
        is not True
        or compositor.get("temporal_retry_counter_is_separate_from_source_and_CFL_counters")
        is not True
        or compositor.get("unexpected_exception_is_invalid_not_physical") is not True
    ):
        raise ValueError("TDG6 production compositor contract differs")

    controls = value["exact_controls"]
    if (
        controls.get("exact_three_halves_squared_boundary") != "pass"
        or controls.get("one_rational_unit_below_squared_boundary") != "fail"
        or controls.get("arbitrarily_small_nonconvergent_pair")
        != "fail_and_retry"
        or controls.get("large_convergent_pair")
        != "pass_but_retain_large_debit"
        or controls.get("RK4_total_shadow_stage_records") != 35
        or controls.get("RK4_committed_fine_stage_records") != 20
        or controls.get("SSPRK3_total_shadow_stage_records") != 28
        or controls.get("SSPRK3_committed_fine_stage_records") != 16
    ):
        raise ValueError("TDG6 exact controls differ")

    alternatives = value["alternatives"]
    if (
        alternatives.get("fixed_absolute_or_relative_state_tolerance_selected")
        is not False
        or alternatives.get("single_coarse_fine_pair_selected") is not False
        or alternatives.get("more_finite_history_samples_selected") is not False
        or alternatives.get("full_nonlinear_a_posteriori_PDE_stability_bound_selected")
        is not False
        or alternatives.get("three_level_order_plus_debit_route_selected")
        is not True
    ):
        raise ValueError("TDG6 alternative decision differs")

    proof = value["proof_contract"]
    if not all(proof.values()):
        raise ValueError("TDG6 proof contract must be all-of")

    successor = value["successor_boundary"]
    if (
        successor.get("TDG6_threshold_and_admission_design_frozen") is not True
        or successor.get("runtime_thresholds_frozen") is not True
        or successor.get("replacement_temporal_admission_defined") is not True
        or successor.get("TDG6_independent_binder_completed") is not False
        or successor.get("production_compositor_implementation_authorized")
        is not False
        or successor.get("production_compositor_implemented") is not False
        or successor.get("PROTO14_frozen") is not False
        or successor.get("GR0_case_eligible") is not False
        or successor.get("SGBL_execution_authorized") is not False
        or successor.get("FGCQR_holdout_execution_authorized") is not False
        or successor.get("DEF1_execution_authorized") is not False
    ):
        raise ValueError("TDG6 successor boundary differs")
    if value["claims"] != EXPECTED_CLAIMS:
        raise ValueError("TDG6 claims differ")


def _validate_lineage(config: Mapping[str, Any]) -> dict[str, Any]:
    commit = config["immutable_lineage"]["checkpoint_commit"]
    ancestry = _git("merge-base", "--is-ancestor", commit, "HEAD", check=False)
    if ancestry.returncode != 0:
        raise ValueError("TDG6 checkpoint is not an ancestor of HEAD")
    observed: dict[str, str] = {}
    for key, path in SOURCE_PATHS.items():
        relative = _rel(path)
        expected = config["immutable_lineage"][f"{key}_sha256"]
        committed = _git("show", f"{commit}:{relative}").stdout
        if _bytes_sha(committed) != expected:
            raise ValueError(f"TDG6 committed source differs: {key}")
        if not path.is_file() or _sha(path) != expected:
            raise ValueError(f"TDG6 worktree source differs: {key}")
        observed[relative] = expected
    _validate_source_result(SOURCE_PATHS["source_result"])
    return {
        "checkpoint_commit": commit,
        "checkpoint_commit_is_ancestor_of_HEAD": True,
        "bound_source_sha256": observed,
        "raw_campaign_bundle_loaded": False,
    }


def _expect_rejection(
    config: Mapping[str, Any], mutate: Callable[[dict[str, Any]], None]
) -> bool:
    changed = deepcopy(config)
    mutate(changed)
    try:
        _validate_config(changed)
    except (KeyError, TypeError, ValueError):
        return True
    return False


def _mutation_controls(config: Mapping[str, Any]) -> dict[str, bool]:
    attacks = {
        "minimum_order_mutation_rejected": lambda value: value["thresholds"].__setitem__(
            "minimum_observed_temporal_refinement_order", "149/100"
        ),
        "squared_boundary_mutation_rejected": lambda value: value["thresholds"].__setitem__(
            "resolved_order_test", "7*D12_upper_squared<=D01_lower_squared"
        ),
        "absolute_tolerance_injection_rejected": lambda value: value["thresholds"].__setitem__(
            "absolute_state_error_threshold_defined", True
        ),
        "outcome_scaling_injection_rejected": lambda value: value["thresholds"].__setitem__(
            "observed_trapped_or_Raychaudhuri_margin_used_to_choose_threshold",
            True,
        ),
        "two_level_path_mutation_rejected": lambda value: value["selected_route"].__setitem__(
            "level_proposal_counts", [1, 2]
        ),
        "medium_commit_mutation_rejected": lambda value: value["selected_route"].__setitem__(
            "outer_and_medium_levels_are_evidence_only", False
        ),
        "retry_count_mutation_rejected": lambda value: value["thresholds"].__setitem__(
            "maximum_temporal_retries_per_macro_step", 31
        ),
        "lineage_hash_mutation_rejected": lambda value: value["immutable_lineage"].__setitem__(
            "source_result_sha256", "0" * 64
        ),
        "PROTO14_promotion_rejected": lambda value: value["claims"].__setitem__(
            "PROTO14_frozen", True
        ),
        "candidate_promotion_rejected": lambda value: value["claims"].__setitem__(
            "FGCQR_holdout_execution_authorized", True
        ),
    }
    result = {
        name: _expect_rejection(config, mutation)
        for name, mutation in attacks.items()
    }
    if not all(result.values()):
        raise AssertionError("TDG6 mutation control failed")
    return result


def record(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    config = _load_config(config_path)
    _validate_config(config)
    lineage = _validate_lineage(config)
    preflight = tdg6_design_preflight()
    if preflight.get("controls_passed") is not True:
        raise AssertionError("TDG6 exact preflight did not pass")
    mutations = _mutation_controls(config)

    implementation = {
        _rel(config_path): _sha(config_path),
        _rel(OWNER_DOCUMENT): _sha(OWNER_DOCUMENT),
        _rel(Path(__file__)): _sha(Path(__file__)),
        _rel(DESIGN_MODULE): _sha(DESIGN_MODULE),
    }
    payload = {
        "immutable_lineage": lineage,
        "selected_route": dict(config["selected_route"]),
        "threshold_contract": dict(config["thresholds"]),
        "error_ledger": dict(config["error_ledger"]),
        "production_compositor_contract": dict(config["production_compositor"]),
        "alternative_decision": dict(config["alternatives"]),
        "exact_design_preflight": preflight,
        "mutation_controls": mutations,
        "proof_contract": dict(config["proof_contract"]),
        "successor_boundary": dict(config["successor_boundary"]),
        "claim_boundary": {
            "actual_PROTO13_history_arrays_consumed": False,
            "campaign_checkpoint_loaded": False,
            "campaign_checkpoint_resumed": False,
            "production_state_advanced": False,
            "historical_CAL10_result_reclassified": False,
            "historical_PROTO13_stop_preserved": True,
            "runtime_thresholds_frozen": True,
            "replacement_temporal_admission_defined": True,
            "absolute_state_error_threshold_defined": False,
            "physical_signal_used_for_threshold_selection": False,
            "TDG6_independent_binder_completed": False,
            "production_compositor_implementation_authorized": False,
            "production_compositor_implemented": False,
            "PROTO14_frozen": False,
            "GR0_case_eligible": False,
            "candidate_branch_opened": False,
            "physical_question_answered": False,
            "rigorous_global_PDE_error_bound_proved": False,
        },
    }
    return {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "classification": "prospective_three_level_temporal_admission_threshold_freeze",
        "gate_status": "PASS_THRESHOLD_AND_ADMISSION_DESIGN_FROZEN",
        "metric_signature": "-+++",
        "riemann_convention": "plus_partial_mu_gamma_nu",
        "source_config_sha256": _sha(config_path),
        "derivation_document": _rel(OWNER_DOCUMENT),
        "derivation_document_sha256": _sha(OWNER_DOCUMENT),
        "implementation_sha256": implementation,
        "generated_by": _rel(Path(__file__)),
        "artifact_payload": payload,
        "nonclaims": [
            "The prospective threshold freeze reads no PROTO13 history or checkpoint and does not resume CAL10.",
            "The three-level piecewise-cubic discriminator is a numerical class, not the exact PDE history.",
            "The finest-pair debit is not a rigorous global PDE error bound or a replacement for the DEF1 stability map.",
            "No production compositor or successor trajectory is implemented or authorized by this freeze.",
            "CAL10 remains historical and unreclassified; no GR-0 case is eligible.",
            "No SGB-L, FGC-QR, collapse, defocusing, transition, retained-EFT, dark-sector, child-domain, or varying-local-c claim follows.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--verify", action="store_true")
    arguments = parser.parse_args()
    value = record(arguments.config)
    payload = _canonical(value).encode("utf-8")
    if arguments.verify:
        if not arguments.output.is_file() or arguments.output.read_bytes() != payload:
            raise SystemExit("FGC-1-TDG6-FRZ1 canonical result differs")
    else:
        _atomic_write(arguments.output, payload)
    print(_canonical(value), end="")


if __name__ == "__main__":
    main()
