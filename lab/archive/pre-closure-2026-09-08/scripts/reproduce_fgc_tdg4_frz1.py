#!/usr/bin/env python3
"""Reproduce the prospective FGC-1-TDG4-FRZ1 theorem freeze."""

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
    TDG4_FUNCTION_SPACE,
    TDG4_INSUFFICIENT_PREMISES,
    TDG4_SAMPLE_COUNT,
    TDG4_SUFFICIENT_ASSUMPTION_ROUTES,
    TDG4_TOP_BINS,
    finite_dimensional_sampling_assessment,
    generic_smooth_kernel_witness,
    lipschitz_coefficient_error_bound,
    sampling_identifiability_preflight,
    uniform_alias_witness,
)


ARTIFACT_ID = "FGC-1-TDG4-FRZ1"
PROJECT_VERSION = "0.11.0"
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-tdg4-frz1.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-tdg4-frz1.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-tdg4-frz1.md"
SOURCE_RESULT = REPOSITORY / "results/fgc-1-tdg3-pref18.json"
SOURCE_CONFIG = REPOSITORY / "configs/fgc/fgc-1-tdg3-pref18.toml"
SOURCE_REPRODUCER = REPOSITORY / "scripts/reproduce_fgc_tdg3_pref18.py"
SOURCE_DISCRIMINATOR = (
    REPOSITORY
    / "src/recursive_horizons/fgc/evolution/tdg3_native_tail_diagnosis.py"
)
SOURCE_DOCUMENT = REPOSITORY / "docs/fgc-tdg3-pref18.md"
THEOREM_MODULE = (
    REPOSITORY
    / "src/recursive_horizons/fgc/evolution/tdg4_sampling_identifiability.py"
)
IMPLEMENTATION = (Path(__file__).resolve(), THEOREM_MODULE)


EXPECTED_CLAIMS = {
    "TDG3_actual_histories_diagnosed": True,
    "TDG3_native_grid_outcome_is_mixed": True,
    "TDG4_sampling_identifiability_theorem_design_frozen": True,
    "TDG4_exact_no_history_controls_pass": True,
    "TDG4_sampling_identifiability_theorem_completed": False,
    "finite_samples_alone_proved_insufficient_for_continuum_top_band_bound": False,
    "current_temporal_spectral_rule_retired_as_continuum_admission": False,
    "replacement_temporal_gate_design_authorized": False,
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


def _expected_lineage() -> dict[str, Any]:
    return {
        "checkpoint_commit": "63f5278d82497380eba57a32331472541b22f5a6",
        "source_result_sha256": "1e12db5903001aeba31338a1d10df4ed96ebab4771509c07d25fd4524ce4f181",
        "source_config_sha256": "c84c8ee99633effb7bb72a5d9602142c8e8c542edfc87f11eb27064050500f78",
        "source_reproducer_sha256": "57342031fef9d7acd91676310ab553c407385da1cafc74c05a54831a7244e3e7",
        "source_discriminator_module_sha256": "9b53ea697cbc8f9026db197849d0bf0b355b92dda0dc22e6bb58d4a4369fb984",
        "source_document_sha256": "d78aff2c2bed81807a6999203701a762e02d95439266d0b1b8ea3e5c4a2b66d0",
        "checkpoint_commit_must_be_ancestor_of_HEAD": True,
        "tracked_source_blobs_must_match_commit_and_worktree": True,
        "source_compact_result_must_verify_without_raw_history": True,
        "raw_campaign_bundle_must_not_be_loaded": True,
    }


def _expected_scope() -> dict[str, Any]:
    return {
        "target": "FGC-2-SF1-TDG4",
        "role": "prospective_finite_sampling_identifiability_theorem",
        "calibration_branch": "GR-0",
        "sampled_history_shape": TDG4_SAMPLE_COUNT,
        "source_compact_result_read": True,
        "actual_terminal_history_arrays_consumed": False,
        "campaign_checkpoint_loaded": False,
        "campaign_checkpoint_resumed": False,
        "state_advanced": False,
        "SGBL_trajectory_read": False,
        "FGCQR_trajectory_read": False,
        "replacement_temporal_admission_defined": False,
        "physical_or_candidate_question_answered": False,
    }


def _expected_successor() -> dict[str, Any]:
    return {
        "TDG3_actual_histories_diagnosed": True,
        "TDG4_sampling_identifiability_theorem_design_frozen": True,
        "TDG4_theorem_execution_authorized_after_freeze_commit": True,
        "TDG4_sampling_identifiability_theorem_completed": False,
        "current_temporal_spectral_rule_retired_as_continuum_admission": False,
        "replacement_temporal_gate_design_authorized": False,
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
        "TDG4 freeze config",
        config,
        {
            "schema_version",
            "artifact_id",
            "project_version",
            "metric_signature",
            "riemann_convention",
            "source_result",
            "scope",
            "immutable_lineage",
            "frozen_theorem",
            "assumption_audit",
            "proof_contract",
            "successor_boundary",
            "claims",
        },
    )
    if (
        config.get("schema_version") != 1
        or config.get("artifact_id") != ARTIFACT_ID
        or config.get("project_version") != PROJECT_VERSION
        or config.get("metric_signature") != "-+++"
        or config.get("riemann_convention") != "plus_partial_mu_gamma_nu"
        or config.get("source_result") != _rel(SOURCE_RESULT)
    ):
        raise ValueError("TDG4 freeze identity differs")
    if config.get("scope") != _expected_scope():
        raise ValueError("TDG4 freeze scope differs")
    if config.get("immutable_lineage") != _expected_lineage():
        raise ValueError("TDG4 freeze lineage differs")
    theorem = config.get("frozen_theorem", {})
    _strict_keys(
        "TDG4 frozen theorem",
        theorem,
        {
            "sample_count",
            "phase_interval",
            "top_bins",
            "function_space",
            "sampling_operator",
            "target_operator",
            "target_field_power",
            "target_phase_derivative_power",
            "theorem_question",
            "actual_history_values_required",
            "replacement_admission_authorized",
            "kernel_criterion",
            "smooth_witness",
            "uniform_alias_control",
            "quantitative_regular_bound",
        },
    )
    if (
        theorem.get("sample_count") != TDG4_SAMPLE_COUNT
        or theorem.get("phase_interval") != [0, 1]
        or tuple(theorem.get("top_bins", ())) != TDG4_TOP_BINS
        or theorem.get("function_space") != TDG4_FUNCTION_SPACE
        or theorem.get("actual_history_values_required") is not False
        or theorem.get("replacement_admission_authorized") is not False
    ):
        raise ValueError("TDG4 theorem shape differs")
    kernel = theorem.get("kernel_criterion", {})
    smooth = theorem.get("smooth_witness", {})
    alias = theorem.get("uniform_alias_control", {})
    regular = theorem.get("quantitative_regular_bound", {})
    if (
        kernel.get("exact_identifiability_iff") != "T annihilates ker(S)"
        or kernel.get("amplitude_is_unrestricted_without_an_external_norm_bound")
        is not True
        or smooth.get("sample_values_and_all_sample_node_derivatives_vanish")
        is not True
        or smooth.get("every_declared_top_bin_must_receive_a_nonzero_coefficient")
        is not True
        or alias.get("carrier_frequency") != 63
        or alias.get("target_complex_coefficient") != "1/4"
        or alias.get("unit_amplitude_positive_bin_field_power") != "1/8"
        or regular.get("sample_values_do_not_supply_L") is not True
    ):
        raise ValueError("TDG4 theorem proof route differs")
    audit = config.get("assumption_audit", {})
    if (
        tuple(audit.get("sufficient_routes", ()))
        != TDG4_SUFFICIENT_ASSUMPTION_ROUTES
        or tuple(audit.get("insufficient_premises", ()))
        != TDG4_INSUFFICIENT_PREMISES
        or audit.get("every_sufficient_route_must_be_independent_of_the_observed_target_margin")
        is not True
        or audit.get("mere_additional_finite_sampling_may_not_be_promoted_to_a_continuum_bound")
        is not True
    ):
        raise ValueError("TDG4 assumption audit differs")
    if set(config.get("proof_contract", {}).values()) != {True}:
        raise ValueError("TDG4 proof contract must be all-of")
    if config.get("successor_boundary") != _expected_successor():
        raise ValueError("TDG4 successor boundary differs")
    if config.get("claims") != EXPECTED_CLAIMS:
        raise ValueError("TDG4 claims differ")


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    with path.open("rb") as handle:
        config = tomllib.load(handle)
    validate_config_data(config)
    return config


def _verify_source_lineage(config: Mapping[str, Any]) -> dict[str, Any]:
    lineage = config["immutable_lineage"]
    commit = lineage["checkpoint_commit"]
    _git("cat-file", "-e", f"{commit}^{{commit}}")
    if _git("merge-base", "--is-ancestor", commit, "HEAD", check=False).returncode:
        raise ValueError("TDG4 source commit is not an ancestor of HEAD")
    tracked = {
        _rel(SOURCE_RESULT): lineage["source_result_sha256"],
        _rel(SOURCE_CONFIG): lineage["source_config_sha256"],
        _rel(SOURCE_REPRODUCER): lineage["source_reproducer_sha256"],
        _rel(SOURCE_DISCRIMINATOR): lineage[
            "source_discriminator_module_sha256"
        ],
        _rel(SOURCE_DOCUMENT): lineage["source_document_sha256"],
    }
    for relative, expected_hash in tracked.items():
        committed = _git("show", f"{commit}:{relative}").stdout
        current = (REPOSITORY / relative).read_bytes()
        if (
            committed != current
            or _bytes_sha(committed) != expected_hash
            or _bytes_sha(current) != expected_hash
        ):
            raise ValueError(f"TDG4 source tracked blob differs: {relative}")
    raw = SOURCE_RESULT.read_bytes()
    source = json.loads(raw)
    compact = source.get("artifact_payload", {}).get("compact_findings", {})
    if (
        raw != _canonical_bytes(source)
        or source.get("artifact_id") != "FGC-1-TDG3-PREF18"
        or source.get("classification")
        != "completed_TDG3_mixed_native_grid_history_diagnosis_sampling_identifiability_theorem_required"
        or source.get("gate_status", {}).get("TDG3_actual_histories_diagnosed")
        is not True
        or source.get("gate_status", {}).get("replacement_temporal_admission_defined")
        is not False
        or compact.get("total_signal_ladders") != 576
        or compact.get("combined_power_classifications") != 1152
        or compact.get("native_estimator_disagreement_count") != 33
        or compact.get("TDG2_finest_zero_survive_native_zero") != 92
        or compact.get("replacement_temporal_admission_earned") is not False
    ):
        raise ValueError("TDG4 source compact result differs")
    return {
        "checkpoint_commit": commit,
        "checkpoint_commit_is_ancestor_of_HEAD": True,
        "tracked_blob_sha256": tracked,
        "source_compact_result_verified": True,
        "source_compact_result_read": True,
        "actual_terminal_history_arrays_consumed": False,
        "campaign_checkpoint_loaded": False,
        "campaign_checkpoint_resumed": False,
        "state_advanced": False,
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


def _mutation_controls(
    config: Mapping[str, Any], preflight: Mapping[str, Any]
) -> dict[str, bool]:
    folded_nodes_rejected = False
    invalid_bin_rejected = False
    invalid_alias_rejected = False
    matrix_shape_rejected = False
    negative_bound_rejected = False
    try:
        generic_smooth_kernel_witness((0, Fraction(1, 2), Fraction(1, 2), 1), target_bin=29)
    except (TypeError, ValueError):
        folded_nodes_rejected = True
    try:
        generic_smooth_kernel_witness((0, 1), target_bin=0)
    except (TypeError, ValueError):
        invalid_bin_rejected = True
    try:
        uniform_alias_witness(sample_count=64, target_bin=63)
    except (TypeError, ValueError):
        invalid_alias_rejected = True
    try:
        finite_dimensional_sampling_assessment(((1, 0),), ((1, 0, 0),))
    except (TypeError, ValueError):
        matrix_shape_rejected = True
    try:
        lipschitz_coefficient_error_bound((0, 1), derivative_bound=-1)
    except (TypeError, ValueError):
        negative_bound_rejected = True
    controls = {
        "folded_node_mutation_rejected": folded_nodes_rejected,
        "invalid_target_bin_rejected": invalid_bin_rejected,
        "invalid_uniform_alias_rejected": invalid_alias_rejected,
        "matrix_domain_mutation_rejected": matrix_shape_rejected,
        "negative_regular_bound_rejected": negative_bound_rejected,
        "sample_count_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["frozen_theorem"].__setitem__(
                "sample_count", 63
            ),
        ),
        "top_bin_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["frozen_theorem"].__setitem__(
                "top_bins", [27, 28, 29, 30, 31]
            ),
        ),
        "function_space_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["frozen_theorem"].__setitem__(
                "function_space", "piecewise_linear_only"
            ),
        ),
        "source_hash_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["immutable_lineage"].__setitem__(
                "source_result_sha256", "0" * 64
            ),
        ),
        "sufficient_route_removal_rejected": _expect_config_rejection(
            config,
            lambda value: value["assumption_audit"].__setitem__(
                "sufficient_routes",
                ["finite_dimensional_class_with_injective_stable_sampling"],
            ),
        ),
        "raw_history_scope_promotion_rejected": _expect_config_rejection(
            config,
            lambda value: value["scope"].__setitem__(
                "actual_terminal_history_arrays_consumed", True
            ),
        ),
        "theorem_completion_promotion_rejected": _expect_config_rejection(
            config,
            lambda value: value["claims"].__setitem__(
                "TDG4_sampling_identifiability_theorem_completed", True
            ),
        ),
        "replacement_design_promotion_rejected": _expect_config_rejection(
            config,
            lambda value: value["claims"].__setitem__(
                "replacement_temporal_gate_design_authorized", True
            ),
        ),
        "GR0_eligibility_promotion_rejected": _expect_config_rejection(
            config,
            lambda value: value["claims"].__setitem__(
                "GR0_case_eligible", True
            ),
        ),
        "exact_preflight_complete": preflight.get("all_exact_controls_pass")
        is True,
    }
    if set(controls.values()) != {True}:
        raise ValueError(f"TDG4 mutation control failed: {controls}")
    return controls


def record(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    config = load_config(config_path)
    lineage = _verify_source_lineage(config)
    preflight = sampling_identifiability_preflight()
    if (
        preflight.get("all_exact_controls_pass") is not True
        or preflight.get("actual_terminal_histories_consumed") is not False
        or preflight.get("campaign_checkpoint_loaded") is not False
        or preflight.get("state_advanced") is not False
        or preflight.get("replacement_temporal_admission_defined") is not False
        or preflight.get("PROTO14_frozen") is not False
        or tuple(preflight.get("top_bins", ())) != TDG4_TOP_BINS
    ):
        raise ValueError("TDG4 exact preflight differs")
    mutations = _mutation_controls(config, preflight)
    result = {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "metric_signature": "-+++",
        "riemann_convention": "plus_partial_mu_gamma_nu",
        "classification": "prospective_TDG4_sampling_identifiability_theorem_freeze",
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
            "frozen_theorem": deepcopy(config["frozen_theorem"]),
            "assumption_audit": deepcopy(config["assumption_audit"]),
            "exact_no_history_preflight": preflight,
            "mutation_controls": mutations,
            "proof_contract": dict(config["proof_contract"]),
            "successor_boundary": dict(config["successor_boundary"]),
            "claim_boundary": {
                "actual_terminal_history_arrays_consumed": False,
                "campaign_checkpoint_loaded": False,
                "state_advanced": False,
                "sampling_identifiability_theorem_completed": False,
                "finite_samples_alone_proved_insufficient": False,
                "current_temporal_spectral_rule_retired": False,
                "replacement_temporal_gate_design_authorized": False,
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
    if (
        value.get("schema_version") != 1
        or value.get("artifact_id") != ARTIFACT_ID
        or value.get("project_version") != PROJECT_VERSION
        or value.get("classification")
        != "prospective_TDG4_sampling_identifiability_theorem_freeze"
        or value.get("generated_by") != "scripts/reproduce_fgc_tdg4_frz1.py"
        or value.get("gate_status") != EXPECTED_CLAIMS
        or value.get("nonclaims")
        != {
            key: False
            for key, item in sorted(EXPECTED_CLAIMS.items())
            if item is False
        }
        or value.get("source_config_sha256")
        != {_rel(DEFAULT_CONFIG): _sha(DEFAULT_CONFIG)}
        or payload.get("frozen_theorem") != config["frozen_theorem"]
        or payload.get("assumption_audit") != config["assumption_audit"]
        or payload.get("successor_boundary") != config["successor_boundary"]
    ):
        raise ValueError("TDG4 stored result boundary differs")
    implementation = value.get("implementation_sha256", {})
    if (
        not isinstance(implementation, Mapping)
        or set(implementation) != {_rel(path) for path in IMPLEMENTATION}
        or any(
            _sha(REPOSITORY / relative) != expected
            for relative, expected in implementation.items()
        )
    ):
        raise ValueError("TDG4 implementation hash ledger differs")
    document = REPOSITORY / str(value.get("derivation_document", ""))
    if (
        value.get("derivation_document") != "docs/fgc-tdg4-frz1.md"
        or not document.is_file()
        or _sha(document) != value.get("derivation_document_sha256")
    ):
        raise ValueError("TDG4 derivation document differs")
    preflight = payload.get("exact_no_history_preflight", {})
    if (
        preflight.get("all_exact_controls_pass") is not True
        or set(preflight.get("generic_smooth_kernel_witnesses", {}))
        != {str(value) for value in TDG4_TOP_BINS}
        or set(preflight.get("exact_uniform_alias_witnesses", {}))
        != {str(value) for value in TDG4_TOP_BINS}
        or preflight.get("actual_terminal_histories_consumed") is not False
    ):
        raise ValueError("TDG4 exact preflight evidence differs")
    mutations = payload.get("mutation_controls", {})
    if not isinstance(mutations, Mapping) or set(mutations.values()) != {True}:
        raise ValueError("TDG4 mutation controls differ")


def verify_canonical(
    config_path: Path = DEFAULT_CONFIG,
    output_path: Path = DEFAULT_OUTPUT,
) -> dict[str, Any]:
    config = load_config(config_path)
    raw = output_path.read_bytes()
    stored = json.loads(raw)
    if not isinstance(stored, dict) or raw != _canonical_bytes(stored):
        raise ValueError("TDG4 stored result is not canonical JSON")
    _validate_stored_record(stored, config)
    observed = record(config_path)
    if stored != observed:
        raise ValueError("TDG4 stored result differs from reproduction")
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
            f"PASS {ARTIFACT_ID}: exact_preflight="
            f"{value['artifact_payload']['exact_no_history_preflight']['all_exact_controls_pass']}"
        )
    else:
        print(_canonical(record(args.config)), end="")


if __name__ == "__main__":
    main()
