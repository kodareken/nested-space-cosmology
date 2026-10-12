#!/usr/bin/env python3
"""Bind the post-freeze TDG4 finite-sampling theorem and its exact scope."""

from __future__ import annotations

import argparse
from copy import deepcopy
from fractions import Fraction
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

from recursive_horizons.fgc.evolution.tdg4_sampling_identifiability import (  # noqa: E402
    sampling_identifiability_preflight,
)
from recursive_horizons.fgc.evolution.tdg4_sampling_identifiability_theorem import (  # noqa: E402
    TDG4_THEOREM_FUNCTION_SPACE,
    TDG4_THEOREM_INSUFFICIENT_PREMISES,
    TDG4_THEOREM_SAMPLE_COUNT,
    TDG4_THEOREM_SUFFICIENT_ROUTES,
    TDG4_THEOREM_TOP_BINS,
    factorization_control,
    sampling_identifiability_theorem_certificate,
    smooth_kernel_witness_theorem,
    uniform_alias_crosscheck,
)


ARTIFACT_ID = "FGC-1-TDG4-PREF19"
PROJECT_VERSION = "0.11.0"
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-tdg4-pref19.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-tdg4-pref19.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-tdg4-pref19.md"
FREEZE_RESULT = REPOSITORY / "results/fgc-1-tdg4-frz1.json"
FREEZE_CONFIG = REPOSITORY / "configs/fgc/fgc-1-tdg4-frz1.toml"
FREEZE_REPRODUCER = REPOSITORY / "scripts/reproduce_fgc_tdg4_frz1.py"
FREEZE_MODULE = (
    REPOSITORY
    / "src/recursive_horizons/fgc/evolution/tdg4_sampling_identifiability.py"
)
FREEZE_DOCUMENT = REPOSITORY / "docs/fgc-tdg4-frz1.md"
PROTOCOL_V4 = REPOSITORY / "configs/fgc/fgc-2-sf1-protocol-v4.toml"
PROTOCOL_V13 = REPOSITORY / "configs/fgc/fgc-2-sf1-protocol-v13.toml"
CAL10_CONFIG = REPOSITORY / "configs/fgc/fgc-1-cal10-pref15.toml"
CAL10_RESULT = REPOSITORY / "results/fgc-1-cal10-pref15.json"
THEOREM_MODULE = (
    REPOSITORY
    / "src/recursive_horizons/fgc/evolution/tdg4_sampling_identifiability_theorem.py"
)
IMPLEMENTATION = (Path(__file__).resolve(), THEOREM_MODULE)


EXPECTED_CLAIMS = {
    "TDG3_actual_histories_diagnosed": True,
    "TDG3_native_grid_outcome_is_mixed": True,
    "TDG4_sampling_identifiability_theorem_design_frozen": True,
    "TDG4_exact_no_history_controls_pass": True,
    "TDG4_sampling_identifiability_theorem_completed": True,
    "finite_samples_alone_proved_insufficient_for_continuum_top_band_bound": True,
    "current_temporal_spectral_rule_retired_as_continuum_admission": True,
    "current_temporal_spectral_rule_retained_as_sampled_diagnostic": True,
    "replacement_temporal_gate_design_authorized": True,
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
        "target": "FGC-2-SF1-TDG4",
        "role": "post_freeze_sampling_identifiability_theorem_binder",
        "calibration_branch": "GR-0",
        "actual_terminal_history_arrays_consumed": False,
        "campaign_checkpoint_loaded": False,
        "campaign_checkpoint_resumed": False,
        "state_advanced": False,
        "SGBL_trajectory_read": False,
        "FGCQR_trajectory_read": False,
        "physical_or_candidate_question_answered": False,
    }


def _expected_lineage() -> dict[str, Any]:
    return {
        "freeze_commit": "8d0df53b32f6668fd29dc085c350f77b3e0c0252",
        "freeze_config_sha256": "b9990512f7cf70b00b05df587454c513c66651a7a740bb900c4e6b587d456cf5",
        "freeze_result_sha256": "d0f6db977a0c7a62db176b12f22479a157ebeb2d7e2483273c43917df8511eb4",
        "freeze_reproducer_sha256": "65582bd00988e551ba0b42f716de01b272b3ac20de067b8b395e594ab3ca1cde",
        "freeze_module_sha256": "31d816df99682114c8dfebff545743bc9624bccad54aff105b094f635f7f0b8d",
        "freeze_document_sha256": "b1d6fa73cd6b5bdfd0b35faa726223dd0e0d6b041ce98fa5f1af4ffd3d50b15c",
        "protocol_v4_sha256": "02a25337e040b97fdd115bc7da4aa87fd09755366599bcade064644ffb3cb87a",
        "protocol_v13_sha256": "fc9398bad6ea4448f659188641e6da6bd7c2e295a443fe73bf73eb0308480692",
        "CAL10_config_sha256": "3661ba8cb6b933ab47daf8897c577b4631544605175dbe15d0c084542cf6c9eb",
        "CAL10_result_sha256": "1ebbdc7222baee8134a6d466589770ee445aee594f10630fc09aa6264cb4ac18",
        "freeze_commit_must_be_ancestor_of_HEAD": True,
        "all_bound_blobs_must_match_freeze_commit_and_worktree": True,
        "freeze_result_must_be_canonical_and_prospective": True,
        "raw_campaign_bundle_must_not_be_loaded": True,
    }


def _expected_theorem() -> dict[str, Any]:
    return {
        "sample_count": TDG4_THEOREM_SAMPLE_COUNT,
        "top_bins": list(TDG4_THEOREM_TOP_BINS),
        "function_space": TDG4_THEOREM_FUNCTION_SPACE,
        "exact_identifiability_iff": "T annihilates ker(S)",
        "equivalent_factorization": (
            "there exists A on range(S) with T=A composed_with S"
        ),
        "smooth_between_sample_kernel_witness_exists_for_every_declared_bin": True,
        "uniform_64_node_alias_crosscheck_exact_for_every_declared_bin": True,
        "finite_samples_alone_bound_continuum_top_band_field_power": False,
        "finite_samples_alone_bound_continuum_top_band_derivative_power": False,
        "both_powers_are_unbounded_on_each_unrestricted_sample_fibre": True,
        "actual_unsampled_power_of_PROTO13_histories_inferred": False,
        "PDE_solution_manifold_restriction_proved": False,
    }


def _expected_retirement() -> dict[str, Any]:
    return {
        "retired_object": (
            "PROTO4_inherited_finite_sample_temporal_spectral_rule_as_a_continuum_admission"
        ),
        "current_rule_retained_as_sampled_diagnostic": True,
        "historical_PROTO13_stop_preserved": True,
        "historical_CAL10_result_reclassified": False,
        "historical_thresholds_changed": False,
        "actual_history_classifications_changed": False,
        "physical_or_candidate_result_inferred": False,
        "replacement_temporal_admission_defined": False,
        "replacement_gate_acceptance_inferred": False,
    }


def _expected_successor() -> dict[str, Any]:
    return {
        "TDG3_actual_histories_diagnosed": True,
        "TDG3_native_grid_outcome_is_mixed": True,
        "TDG4_sampling_identifiability_theorem_design_frozen": True,
        "TDG4_sampling_identifiability_theorem_completed": True,
        "finite_samples_alone_proved_insufficient_for_continuum_top_band_bound": True,
        "current_temporal_spectral_rule_retired_as_continuum_admission": True,
        "current_temporal_spectral_rule_retained_as_sampled_diagnostic": True,
        "replacement_temporal_gate_design_authorized": True,
        "replacement_temporal_admission_defined": False,
        "PROTO14_frozen": False,
        "fresh_GR0_dynamic_calibration_completed": False,
        "GR0_case_eligible": False,
        "SGBL_execution_authorized": False,
        "FGCQR_holdout_execution_authorized": False,
        "DEF1_execution_authorized": False,
    }


def validate_config_data(config: Mapping[str, Any]) -> None:
    _strict_keys(
        "PREF19 config",
        config,
        {
            "schema_version",
            "artifact_id",
            "project_version",
            "metric_signature",
            "riemann_convention",
            "freeze_result",
            "freeze_config",
            "protocol_v4",
            "protocol_v13",
            "CAL10_config",
            "CAL10_result",
            "scope",
            "immutable_lineage",
            "theorem_result",
            "retirement_scope",
            "assumption_audit",
            "proof_contract",
            "successor_boundary",
            "claims",
        },
    )
    expected_paths = {
        "freeze_result": _rel(FREEZE_RESULT),
        "freeze_config": _rel(FREEZE_CONFIG),
        "protocol_v4": _rel(PROTOCOL_V4),
        "protocol_v13": _rel(PROTOCOL_V13),
        "CAL10_config": _rel(CAL10_CONFIG),
        "CAL10_result": _rel(CAL10_RESULT),
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
        or config.get("theorem_result") != _expected_theorem()
        or config.get("retirement_scope") != _expected_retirement()
    ):
        raise ValueError("PREF19 identity, lineage, theorem, or scope differs")
    audit = config.get("assumption_audit", {})
    if audit != {
        "sufficient_routes": list(TDG4_THEOREM_SUFFICIENT_ROUTES),
        "insufficient_premises": list(TDG4_THEOREM_INSUFFICIENT_PREMISES),
        "replacement_design_must_choose_and_quantify_at_least_one_sufficient_route": True,
        "replacement_bound_must_be_independent_of_the_observed_target_margin": True,
    }:
        raise ValueError("PREF19 assumption audit differs")
    if set(config.get("proof_contract", {}).values()) != {True}:
        raise ValueError("PREF19 proof contract must be all-of")
    if config.get("successor_boundary") != _expected_successor():
        raise ValueError("PREF19 successor boundary differs")
    if config.get("claims") != EXPECTED_CLAIMS:
        raise ValueError("PREF19 claim boundary differs")


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    with path.open("rb") as handle:
        config = tomllib.load(handle)
    validate_config_data(config)
    return config


def _load_toml(path: Path) -> dict[str, Any]:
    with path.open("rb") as handle:
        return tomllib.load(handle)


def _verify_lineage(config: Mapping[str, Any]) -> dict[str, Any]:
    lineage = config["immutable_lineage"]
    commit = lineage["freeze_commit"]
    _git("cat-file", "-e", f"{commit}^{{commit}}")
    if _git("merge-base", "--is-ancestor", commit, "HEAD", check=False).returncode:
        raise ValueError("PREF19 freeze commit is not an ancestor of HEAD")
    tracked = {
        _rel(FREEZE_CONFIG): lineage["freeze_config_sha256"],
        _rel(FREEZE_RESULT): lineage["freeze_result_sha256"],
        _rel(FREEZE_REPRODUCER): lineage["freeze_reproducer_sha256"],
        _rel(FREEZE_MODULE): lineage["freeze_module_sha256"],
        _rel(FREEZE_DOCUMENT): lineage["freeze_document_sha256"],
        _rel(PROTOCOL_V4): lineage["protocol_v4_sha256"],
        _rel(PROTOCOL_V13): lineage["protocol_v13_sha256"],
        _rel(CAL10_CONFIG): lineage["CAL10_config_sha256"],
        _rel(CAL10_RESULT): lineage["CAL10_result_sha256"],
    }
    for relative, expected_hash in tracked.items():
        committed = _git("show", f"{commit}:{relative}").stdout
        current = (REPOSITORY / relative).read_bytes()
        if (
            committed != current
            or _bytes_sha(committed) != expected_hash
            or _bytes_sha(current) != expected_hash
        ):
            raise ValueError(f"PREF19 bound blob differs: {relative}")

    freeze_raw = FREEZE_RESULT.read_bytes()
    freeze = json.loads(freeze_raw)
    freeze_status = freeze.get("gate_status", {})
    if (
        freeze_raw != _canonical_bytes(freeze)
        or freeze.get("artifact_id") != "FGC-1-TDG4-FRZ1"
        or freeze.get("classification")
        != "prospective_TDG4_sampling_identifiability_theorem_freeze"
        or freeze_status.get(
            "TDG4_sampling_identifiability_theorem_design_frozen"
        )
        is not True
        or freeze_status.get("TDG4_exact_no_history_controls_pass") is not True
        or freeze_status.get("TDG4_sampling_identifiability_theorem_completed")
        is not False
        or freeze_status.get(
            "current_temporal_spectral_rule_retired_as_continuum_admission"
        )
        is not False
        or freeze_status.get("replacement_temporal_gate_design_authorized")
        is not False
        or freeze.get("artifact_payload", {})
        .get("claim_boundary", {})
        .get("actual_terminal_history_arrays_consumed")
        is not False
    ):
        raise ValueError("PREF19 prospective freeze record differs")

    protocol_v4 = _load_toml(PROTOCOL_V4)
    spectrum = protocol_v4.get("numerical_admission", {}).get("spectrum", {})
    if (
        protocol_v4.get("artifact_id") != "FGC-2-SF1-PROTO4"
        or spectrum.get("minimum_temporal_samples") != 64
        or spectrum.get("maximum_top_band_field_power_fraction")
        != "1/1048576"
        or spectrum.get("maximum_top_band_derivative_weighted_power_fraction")
        != "1/1024"
        or spectrum.get("maximum_nested_refinement_tail_ratio") != "1/4"
    ):
        raise ValueError("PREF19 PROTO4 temporal rule differs")

    protocol_v13 = _load_toml(PROTOCOL_V13)
    amendment = protocol_v13.get("amendment", {})
    inheritance = protocol_v13.get("inheritance", {})
    if (
        protocol_v13.get("artifact_id") != "FGC-2-SF1-PROTO13"
        or amendment.get("absolute_spectral_budgets_and_direct_ratio_ceiling_unchanged")
        is not True
        or inheritance.get(
            "PROTO12_CFL_dissipation_retries_health_boundary_scale_and_transaction_rules_inherited_unchanged"
        )
        is not True
    ):
        raise ValueError("PREF19 PROTO13 inheritance differs")

    cal10_config = _load_toml(CAL10_CONFIG)
    temporal = cal10_config.get("temporal_stop", {})
    if (
        cal10_config.get("artifact_id") != "FGC-1-CAL10-PREF15"
        or temporal.get("minimum_sample_count") != 64
        or temporal.get("observed_sample_count") != 64
        or temporal.get("both_method_admissions_failed") is not True
        or temporal.get("both_individual_budget_layers_failed") is not True
        or temporal.get("both_nested_ratio_layers_failed") is not True
        or temporal.get("numerical_or_physical_cause_derived") is not False
    ):
        raise ValueError("PREF19 CAL10 contract differs")

    cal10_raw = CAL10_RESULT.read_bytes()
    cal10 = json.loads(cal10_raw)
    diagnosis = cal10.get("artifact_payload", {}).get("campaign_diagnosis", {})
    boundary = cal10.get("artifact_payload", {}).get("claim_boundary", {})
    if (
        cal10_raw != _canonical_bytes(cal10)
        or cal10.get("artifact_id") != "FGC-1-CAL10-PREF15"
        or cal10.get("gate_status", {}).get("PROTO13_temporal_admission_failed")
        is not True
        or cal10.get("gate_status", {}).get("PROTO13_temporal_failure_cause_derived")
        is not False
        or diagnosis.get("temporal_minimum_samples_reached") is not True
        or diagnosis.get("both_temporal_method_admissions_failed") is not True
        or diagnosis.get("temporal_failure_cause_derived") is not False
        or any(
            method.get("sample_count") != 64
            for method in diagnosis.get("temporal_methods", {}).values()
        )
        or boundary.get("temporal_stop_is_not_a_candidate_action_result")
        is not True
        or boundary.get("temporal_stop_is_not_a_physical_obstruction") is not True
    ):
        raise ValueError("PREF19 CAL10 compact evidence differs")
    return {
        "freeze_commit": commit,
        "freeze_commit_is_ancestor_of_HEAD": True,
        "bound_blob_sha256": tracked,
        "freeze_result_verified_prospective": True,
        "PROTO4_temporal_rule": {
            "minimum_temporal_samples": spectrum["minimum_temporal_samples"],
            "maximum_top_band_field_power_fraction": spectrum[
                "maximum_top_band_field_power_fraction"
            ],
            "maximum_top_band_derivative_weighted_power_fraction": spectrum[
                "maximum_top_band_derivative_weighted_power_fraction"
            ],
            "maximum_nested_refinement_tail_ratio": spectrum[
                "maximum_nested_refinement_tail_ratio"
            ],
        },
        "PROTO13_absolute_spectral_contract_inherited_unchanged": True,
        "CAL10_observed_sample_count": temporal["observed_sample_count"],
        "CAL10_both_method_admissions_failed": True,
        "CAL10_failure_cause_was_underived": True,
        "actual_terminal_history_arrays_consumed": False,
        "campaign_checkpoint_loaded": False,
        "campaign_checkpoint_resumed": False,
        "state_advanced": False,
    }


def _freeze_crosscheck(theorem: Mapping[str, Any]) -> dict[str, Any]:
    preflight = sampling_identifiability_preflight()
    if (
        preflight.get("all_exact_controls_pass") is not True
        or preflight.get("actual_terminal_histories_consumed") is not False
        or preflight.get("campaign_checkpoint_loaded") is not False
        or preflight.get("state_advanced") is not False
    ):
        raise ValueError("PREF19 independent freeze preflight differs")
    for frequency in TDG4_THEOREM_TOP_BINS:
        key = str(frequency)
        frozen_witness = preflight["generic_smooth_kernel_witnesses"][key]
        theorem_witness = theorem["smooth_kernel_witnesses"][key]
        frozen_alias = preflight["exact_uniform_alias_witnesses"][key]
        theorem_alias = theorem["uniform_alias_crosschecks"][key]
        if (
            frozen_witness.get("sample_gap")
            != theorem_witness.get("selected_gap")
            or frozen_witness.get("support_interval")
            != theorem_witness.get("support_interval")
            or frozen_witness.get("target_coefficient_nonzero_for_nonzero_amplitude")
            is not theorem_witness.get(
                "target_coefficient_nonzero_for_every_nonzero_amplitude"
            )
            or frozen_alias.get("target_complex_coefficient")
            != theorem_alias.get("target_complex_coefficient")
            or frozen_alias.get("partner_complex_coefficient")
            != theorem_alias.get("partner_complex_coefficient")
            or frozen_alias.get("unit_amplitude_positive_bin_field_power")
            != theorem_alias.get("positive_bin_field_power")
        ):
            raise ValueError(f"PREF19 freeze crosscheck differs at bin {frequency}")
    return {
        "freeze_exact_preflight_reexecuted": True,
        "factorization_controls_match_frozen_design": True,
        "smooth_support_intervals_match_frozen_design": True,
        "uniform_alias_coefficients_and_power_match_frozen_design": True,
        "all_five_declared_top_bins_crosschecked": True,
        "actual_terminal_histories_consumed": False,
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
    folded_nodes_rejected = False
    invalid_bin_rejected = False
    matrix_domain_rejected = False
    try:
        smooth_kernel_witness_theorem(
            (0, Fraction(1, 2), Fraction(1, 2), 1), target_bin=28
        )
    except (TypeError, ValueError):
        folded_nodes_rejected = True
    try:
        uniform_alias_crosscheck(sample_count=64, target_bin=63)
    except (TypeError, ValueError):
        invalid_bin_rejected = True
    try:
        factorization_control(((1, 0),), ((1, 0, 0),))
    except (TypeError, ValueError):
        matrix_domain_rejected = True
    controls = {
        "folded_node_mutation_rejected": folded_nodes_rejected,
        "invalid_alias_bin_mutation_rejected": invalid_bin_rejected,
        "matrix_domain_mutation_rejected": matrix_domain_rejected,
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
        "protocol_hash_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["immutable_lineage"].__setitem__(
                "protocol_v4_sha256", "0" * 64
            ),
        ),
        "sample_count_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["theorem_result"].__setitem__(
                "sample_count", 63
            ),
        ),
        "top_bin_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["theorem_result"].__setitem__(
                "top_bins", [27, 28, 29, 30, 31]
            ),
        ),
        "actual_history_inference_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["theorem_result"].__setitem__(
                "actual_unsampled_power_of_PROTO13_histories_inferred", True
            ),
        ),
        "PDE_manifold_promotion_rejected": _expect_config_rejection(
            config,
            lambda value: value["theorem_result"].__setitem__(
                "PDE_solution_manifold_restriction_proved", True
            ),
        ),
        "historical_reclassification_rejected": _expect_config_rejection(
            config,
            lambda value: value["retirement_scope"].__setitem__(
                "historical_CAL10_result_reclassified", True
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
    }
    if set(controls.values()) != {True}:
        raise ValueError(f"PREF19 mutation control failed: {controls}")
    return controls


def record(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    config = load_config(config_path)
    lineage = _verify_lineage(config)
    theorem = sampling_identifiability_theorem_certificate()
    if (
        theorem.get("sampling_identifiability_theorem_completed") is not True
        or theorem.get("sample_count") != TDG4_THEOREM_SAMPLE_COUNT
        or tuple(theorem.get("target_bins", ())) != TDG4_THEOREM_TOP_BINS
        or theorem.get("function_space") != TDG4_THEOREM_FUNCTION_SPACE
        or theorem.get("finite_samples_alone_bound_field_power") is not False
        or theorem.get("finite_samples_alone_bound_derivative_power") is not False
        or theorem.get("actual_unsampled_power_of_the_PROTO13_histories_inferred")
        is not False
        or theorem.get("PDE_solution_manifold_restriction_proved") is not False
        or theorem.get("actual_terminal_histories_consumed") is not False
        or theorem.get("campaign_checkpoint_loaded") is not False
        or theorem.get("state_advanced") is not False
    ):
        raise ValueError("PREF19 theorem certificate differs")
    crosscheck = _freeze_crosscheck(theorem)
    mutations = _mutation_controls(config)
    result = {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "metric_signature": "-+++",
        "riemann_convention": "plus_partial_mu_gamma_nu",
        "classification": (
            "completed_TDG4_sampling_identifiability_theorem_"
            "current_sampled_rule_retired_as_continuum_admission"
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
            "theorem_certificate": theorem,
            "freeze_crosscheck": crosscheck,
            "retirement_scope": deepcopy(config["retirement_scope"]),
            "assumption_audit": deepcopy(config["assumption_audit"]),
            "mutation_controls": mutations,
            "proof_contract": dict(config["proof_contract"]),
            "successor_boundary": dict(config["successor_boundary"]),
            "claim_boundary": {
                "actual_terminal_history_arrays_consumed": False,
                "campaign_checkpoint_loaded": False,
                "state_advanced": False,
                "actual_unsampled_PROTO13_history_power_inferred": False,
                "PDE_solution_manifold_restriction_proved": False,
                "historical_PROTO13_stop_preserved": True,
                "historical_CAL10_result_reclassified": False,
                "current_rule_retired_only_as_continuum_admission": True,
                "current_rule_retained_as_sampled_diagnostic": True,
                "replacement_temporal_gate_design_authorized": True,
                "replacement_temporal_admission_defined": False,
                "new_trajectory_authorized": False,
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
    theorem = payload.get("theorem_certificate", {})
    boundary = payload.get("claim_boundary", {})
    if (
        value.get("schema_version") != 1
        or value.get("artifact_id") != ARTIFACT_ID
        or value.get("project_version") != PROJECT_VERSION
        or value.get("classification")
        != "completed_TDG4_sampling_identifiability_theorem_current_sampled_rule_retired_as_continuum_admission"
        or value.get("generated_by") != "scripts/reproduce_fgc_tdg4_pref19.py"
        or value.get("gate_status") != EXPECTED_CLAIMS
        or value.get("nonclaims")
        != {
            key: False
            for key, item in sorted(EXPECTED_CLAIMS.items())
            if item is False
        }
        or value.get("source_config_sha256")
        != {_rel(DEFAULT_CONFIG): _sha(DEFAULT_CONFIG)}
        or payload.get("retirement_scope") != config["retirement_scope"]
        or payload.get("assumption_audit") != config["assumption_audit"]
        or payload.get("successor_boundary") != config["successor_boundary"]
        or theorem.get("sampling_identifiability_theorem_completed") is not True
        or theorem.get("actual_unsampled_power_of_the_PROTO13_histories_inferred")
        is not False
        or theorem.get("PDE_solution_manifold_restriction_proved") is not False
        or boundary.get("current_rule_retired_only_as_continuum_admission")
        is not True
        or boundary.get("historical_PROTO13_stop_preserved") is not True
        or boundary.get("historical_CAL10_result_reclassified") is not False
        or boundary.get("replacement_temporal_gate_design_authorized") is not True
        or boundary.get("replacement_temporal_admission_defined") is not False
        or boundary.get("new_trajectory_authorized") is not False
        or boundary.get("candidate_branch_opened") is not False
    ):
        raise ValueError("PREF19 stored result boundary differs")
    implementation = value.get("implementation_sha256", {})
    if (
        not isinstance(implementation, Mapping)
        or set(implementation) != {_rel(path) for path in IMPLEMENTATION}
        or any(
            _sha(REPOSITORY / relative) != expected
            for relative, expected in implementation.items()
        )
    ):
        raise ValueError("PREF19 implementation hash ledger differs")
    document = REPOSITORY / str(value.get("derivation_document", ""))
    if (
        value.get("derivation_document") != "docs/fgc-tdg4-pref19.md"
        or not document.is_file()
        or _sha(document) != value.get("derivation_document_sha256")
    ):
        raise ValueError("PREF19 derivation document differs")
    if (
        not isinstance(payload.get("mutation_controls"), Mapping)
        or set(payload["mutation_controls"].values()) != {True}
        or not isinstance(payload.get("proof_contract"), Mapping)
        or set(payload["proof_contract"].values()) != {True}
    ):
        raise ValueError("PREF19 proof or mutation ledger differs")


def verify_canonical(
    config_path: Path = DEFAULT_CONFIG,
    output_path: Path = DEFAULT_OUTPUT,
) -> dict[str, Any]:
    config = load_config(config_path)
    raw = output_path.read_bytes()
    stored = json.loads(raw)
    if not isinstance(stored, dict) or raw != _canonical_bytes(stored):
        raise ValueError("PREF19 stored result is not canonical JSON")
    _validate_stored_record(stored, config)
    observed = record(config_path)
    if stored != observed:
        raise ValueError("PREF19 stored result differs from reproduction")
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
        print(
            f"PASS {ARTIFACT_ID}: theorem="
            f"{value['gate_status']['TDG4_sampling_identifiability_theorem_completed']} "
            "replacement_defined="
            f"{value['gate_status']['replacement_temporal_admission_defined']}"
        )
    else:
        print(_canonical(record(args.config)), end="")


if __name__ == "__main__":
    main()
