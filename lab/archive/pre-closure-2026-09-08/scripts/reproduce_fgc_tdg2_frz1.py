#!/usr/bin/env python3
"""Reproduce the prospective TDG2 absolute-tail discriminator freeze."""

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

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.tdg2_absolute_tail_diagnosis import (  # noqa: E402
    TDG2_ABSOLUTE_AMPLITUDE_FLOOR,
    TDG2_CLASS_NAMES,
    TDG2_DECIMAL_PRECISION_HIGH,
    TDG2_DECIMAL_PRECISION_LOW,
    TDG2_FIELD_NAMES,
    TDG2_MEASURE_NAMES,
    TDG2_METHOD_POINT_COUNTS,
    TDG2_MINIMUM_HISTORY_ORDER,
    TDG2_MINIMUM_ZERO_POWER_ORDER,
    TDG2_SAMPLE_COUNT,
    TDG2_SYNTHETIC_CONTROL_NAMES,
    TDG2_TAPER_FRACTION,
    TDG2_TRACER_COUNT,
    TDG2_UNRESOLVED_TOP_FRACTION,
    absolute_tail_observation,
    synthetic_absolute_tail_controls,
)


ARTIFACT_ID = "FGC-1-TDG2-FRZ1"
PROJECT_VERSION = "0.11.0"
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-tdg2-frz1.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-tdg2-frz1.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-tdg2-frz1.md"
SOURCE_RESULT = REPOSITORY / "results/fgc-1-tdg1-pref16.json"
SOURCE_CONFIG = REPOSITORY / "configs/fgc/fgc-1-tdg1-pref16.toml"
SOURCE_REPRODUCER = REPOSITORY / "scripts/reproduce_fgc_tdg1_pref16.py"
SOURCE_DIAGNOSIS_MODULE = (
    REPOSITORY
    / "src/recursive_horizons/fgc/evolution/tdg1_temporal_diagnosis.py"
)
DISCRIMINATOR_MODULE = (
    REPOSITORY
    / "src/recursive_horizons/fgc/evolution/tdg2_absolute_tail_diagnosis.py"
)
IMPLEMENTATION = (Path(__file__).resolve(), DISCRIMINATOR_MODULE)


EXPECTED_CLAIMS = {
    "TDG1_diagnostic_contract_frozen": True,
    "TDG1_terminal_histories_diagnosed": True,
    "normalized_raw_temporal_tail_ratio_is_valid_general_spatial_convergence_test": False,
    "TDG2_absolute_tail_discriminator_frozen": True,
    "TDG2_actual_histories_diagnosed": False,
    "individual_temporal_budget_failure_cause_fully_derived": False,
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


def _fraction(value: str) -> float:
    numerator, separator, denominator = value.partition("/")
    return float(int(numerator) / int(denominator)) if separator else float(value)


def _expected_lineage() -> dict[str, Any]:
    return {
        "checkpoint_commit": "873e13b5ce74624e2cbc096694eb1180b4dbc0a3",
        "source_result_sha256": "1713a05b6442461962253ec9abb84f50dd6c376e91e7e13447c4682256b0c3a6",
        "source_config_sha256": "d598b9f59fc50dab23f05150cd6471ec66e28f36cddfe60438c8402bdf70a519",
        "source_reproducer_sha256": "2fd21140147b8ddf326cc61fe4ecc29832cbad775a0af03a3854e40f4e209c00",
        "source_diagnosis_module_sha256": "3aaf100a8c9abb691978ac3ddf1f58210ed59107dc6018cf3f8b70c99de8ba44",
        "source_checkpoint_sha256": "0432b266081b1ec1af3ac2a9efb99bcdd67dc8afef837241ec27fe5d37432f97",
        "checkpoint_commit_must_be_ancestor_of_HEAD": True,
        "tracked_source_blobs_must_match_checkpoint_commit_and_worktree": True,
        "source_checkpoint_is_read_only": True,
        "source_checkpoint_may_not_resume": True,
    }


def _expected_frozen_discriminator() -> dict[str, Any]:
    return {
        "sample_count": TDG2_SAMPLE_COUNT,
        "tracer_count": TDG2_TRACER_COUNT,
        "field_names": list(TDG2_FIELD_NAMES),
        "primary_method": "RK4",
        "primary_point_counts": list(TDG2_METHOD_POINT_COUNTS["RK4"]),
        "comparator_method": "SSPRK3",
        "comparator_point_counts": list(TDG2_METHOD_POINT_COUNTS["SSPRK3"]),
        "taper_fraction": "1/8",
        "unresolved_top_fraction": "1/8",
        "top_band_first_bin": 28,
        "nyquist_bin": 32,
        "top_band_bin_count": 5,
        "declared_absolute_DFT_amplitude_floor": "1e-30",
        "minimum_history_difference_order": "3/2",
        "minimum_zero_power_order": "3",
        "power_measures": list(TDG2_MEASURE_NAMES),
        "normalization_by_total_power_forbidden": True,
        "common_proper_time_intersection_required": True,
        "exactly_64_uniform_aligned_samples_required": True,
        "classification_is_per_method_per_tracer_per_field_per_measure": True,
    }


def _expected_synthetic() -> dict[str, Any]:
    expected = {
        "exact_zero": "instrument_floor_or_interpolation_enclosure_dominated",
        "subfloor_nonzero": "instrument_floor_or_interpolation_enclosure_dominated",
        "power_contracts_to_zero": "contracts_to_zero_at_required_power_order",
        "nonzero_resolution_independent": "nonzero_resolution_independent_within_enclosure",
        "nonzero_convergent": "converges_to_nonzero_resolved_power",
        "unresolved_nonmonotone": "unresolved_absolute_tail_behavior",
        "interpolation_enclosure_dominated": "instrument_floor_or_interpolation_enclosure_dominated",
    }
    return {
        "control_names": list(TDG2_SYNTHETIC_CONTROL_NAMES),
        **{f"{name}_expected": value for name, value in expected.items()},
        "both_power_measures_must_match_every_expected_class": True,
        "controls_consume_no_campaign_histories": True,
    }


def validate_config_data(config: Mapping[str, Any]) -> None:
    _strict_keys(
        "TDG2 freeze config",
        config,
        {
            "schema_version",
            "artifact_id",
            "project_version",
            "metric_signature",
            "riemann_convention",
            "source_result",
            "source_config",
            "source_checkpoint",
            "scope",
            "immutable_lineage",
            "frozen_discriminator",
            "arithmetic_contract",
            "interpolation_contract",
            "classification_contract",
            "synthetic_controls",
            "outcome_classes",
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
        or config.get("source_config") != _rel(SOURCE_CONFIG)
        or config.get("source_checkpoint")
        != "runs/fgc-2-sf1/proto13/calibration/latest-checkpoint.npz"
    ):
        raise ValueError("TDG2 freeze top-level contract differs")
    if config.get("scope") != {
        "target": "FGC-2-SF1-TDG2",
        "role": "prospective_no_history_absolute_tail_discriminator_freeze",
        "calibration_branch": "GR-0",
        "source_terminal_event": 63,
        "source_terminal_coordinate_time": "63/16",
        "TDG1_terminal_diagnosis_already_bound": True,
        "source_checkpoint_history_arrays_loaded": False,
        "source_checkpoint_resumed": False,
        "state_advanced": False,
        "SGBL_trajectory_read": False,
        "FGCQR_trajectory_read": False,
        "replacement_temporal_admission_defined": False,
        "physical_or_candidate_question_answered": False,
    }:
        raise ValueError("TDG2 freeze scope differs")
    if config.get("immutable_lineage") != _expected_lineage():
        raise ValueError("TDG2 immutable lineage differs")
    frozen = config.get("frozen_discriminator", {})
    if frozen != _expected_frozen_discriminator():
        raise ValueError("TDG2 discriminator constants differ")
    if (
        _fraction(frozen["taper_fraction"]) != TDG2_TAPER_FRACTION
        or _fraction(frozen["unresolved_top_fraction"])
        != TDG2_UNRESOLVED_TOP_FRACTION
        or _fraction(frozen["declared_absolute_DFT_amplitude_floor"])
        != TDG2_ABSOLUTE_AMPLITUDE_FLOOR
        or _fraction(frozen["minimum_history_difference_order"])
        != TDG2_MINIMUM_HISTORY_ORDER
        or _fraction(frozen["minimum_zero_power_order"])
        != TDG2_MINIMUM_ZERO_POWER_ORDER
    ):
        raise ValueError("TDG2 numerical constants differ")
    if config.get("arithmetic_contract") != {
        "low_decimal_precision": TDG2_DECIMAL_PRECISION_LOW,
        "high_decimal_precision": TDG2_DECIMAL_PRECISION_HIGH,
        "independent_DFT_uses_only_top_band_bins": True,
        "binary64_rFFT_is_audited_against_decimal_DFT": True,
        "power_interval_contains_binary64_decimal_disagreement": True,
        "power_interval_uses_outward_binary64_rounding": True,
        "coefficient_interval_contains_binary64_decimal_disagreement": True,
        "two_decimal_precisions_must_agree_within_public_debit": True,
        "arithmetic_disagreement_may_not_be_silently_dropped": True,
    }:
        raise ValueError("TDG2 arithmetic contract differs")
    if config.get("interpolation_contract") != {
        "alignment": "piecewise_linear_on_common_proper_time_intersection",
        "round_trip_infinity_debit_is_public": True,
        "round_trip_debit_is_used_only_as_conditional_pointwise_perturbation_radius": True,
        "round_trip_debit_is_not_claimed_as_continuum_interpolation_theorem": True,
        "coefficient_perturbation_bound": "E_k <= ||window||_1 delta",
        "field_power_perturbation_bound": "sum m_k (2 |F_k| E_k + E_k^2) / N^2",
        "derivative_power_perturbation_bound": "sum omega_k^2 m_k (2 |F_k| E_k + E_k^2) / N^2",
        "arithmetic_and_interpolation_debits_remain_separate": True,
    }:
        raise ValueError("TDG2 interpolation contract differs")
    if config.get("classification_contract") != {
        "class_names": list(TDG2_CLASS_NAMES),
        "zero_class_requires_both_power_ratios_to_have_conservative_order_at_least_3": True,
        "resolution_independent_nonzero_class_requires_positive_three_interval_intersection": True,
        "nonzero_convergent_class_requires_same_direction_difference_contraction_at_order_at_least_3_over_2": True,
        "nonzero_convergent_decreasing_class_requires_positive_geometric_continuum_lower_bound": True,
        "floor_class_requires_all_three_top_band_coefficient_intervals_below_1e_minus_30": True,
        "finest_interval_touching_zero_is_enclosure_dominated": True,
        "nonmonotone_or_below_order_evidence_is_unresolved": True,
        "field_and_derivative_power_are_classified_separately": True,
        "no_class_is_a_pass_fail_temporal_admission": True,
    }:
        raise ValueError("TDG2 classification contract differs")
    if config.get("synthetic_controls") != _expected_synthetic():
        raise ValueError("TDG2 synthetic controls differ")
    if config.get("outcome_classes") != {
        "synthetic_discriminator_validated": "freeze_contract_reproduces_and_actual_history_diagnosis_may_execute",
        "synthetic_discriminator_not_validated": "implementation_or_contract_failure_no_actual_history_read",
        "actual_history_diagnosis_pending": "no_terminal_history_conclusion_at_freeze",
        "actual_history_result_must_preserve_all_per_signal_classes": True,
        "actual_history_result_may_not_retroactively_pass_PROTO13": True,
        "actual_history_result_may_not_directly_define_replacement_gate": True,
    }:
        raise ValueError("TDG2 outcome classes differ")
    if set(config.get("proof_contract", {}).values()) != {True}:
        raise ValueError("TDG2 proof contract must be all-of")
    if config.get("successor_boundary") != {
        "TDG1_terminal_histories_diagnosed": True,
        "TDG2_absolute_tail_discriminator_frozen": True,
        "TDG2_actual_history_execution_authorized": True,
        "TDG2_actual_history_execution_completed": False,
        "replacement_temporal_gate_design_authorized": False,
        "replacement_temporal_admission_defined": False,
        "PROTO14_frozen": False,
        "fresh_GR0_dynamic_calibration_completed": False,
        "GR0_case_eligible": False,
        "SGBL_execution_authorized": False,
        "FGCQR_holdout_execution_authorized": False,
        "DEF1_execution_authorized": False,
    }:
        raise ValueError("TDG2 successor boundary differs")
    if config.get("claims") != EXPECTED_CLAIMS:
        raise ValueError("TDG2 claims differ")


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    with path.open("rb") as handle:
        config = tomllib.load(handle)
    validate_config_data(config)
    return config


def _verify_optional_raw_hash(path: Path, expected_hash: str) -> None:
    """Verify the ignored checkpoint when present without loading its arrays."""

    if not path.exists():
        return
    if not path.is_file() or _sha(path) != expected_hash:
        raise ValueError("TDG2 source checkpoint hash differs")


def _verify_lineage(config: Mapping[str, Any]) -> dict[str, Any]:
    lineage = config["immutable_lineage"]
    commit = lineage["checkpoint_commit"]
    _git("cat-file", "-e", f"{commit}^{{commit}}")
    if _git("merge-base", "--is-ancestor", commit, "HEAD", check=False).returncode:
        raise ValueError("TDG2 checkpoint commit is not an ancestor of HEAD")
    tracked = {
        _rel(SOURCE_RESULT): lineage["source_result_sha256"],
        _rel(SOURCE_CONFIG): lineage["source_config_sha256"],
        _rel(SOURCE_REPRODUCER): lineage["source_reproducer_sha256"],
        _rel(SOURCE_DIAGNOSIS_MODULE): lineage["source_diagnosis_module_sha256"],
    }
    blobs: dict[str, str] = {}
    for relative, expected_hash in tracked.items():
        committed = _git("show", f"{commit}:{relative}").stdout
        current = (REPOSITORY / relative).read_bytes()
        if committed != current or _bytes_sha(committed) != expected_hash:
            raise ValueError(f"TDG2 immutable source differs: {relative}")
        blobs[relative] = expected_hash
    source_result_bytes = SOURCE_RESULT.read_bytes()
    source_result = json.loads(source_result_bytes)
    if (
        source_result_bytes != _canonical_bytes(source_result)
        or source_result.get("artifact_id") != "FGC-1-TDG1-PREF16"
        or source_result.get("classification")
        != "completed_TDG1_history_diagnosis_TDG2_required"
        or source_result.get("gate_status", {}).get("TDG2_frozen") is not False
        or source_result.get("artifact_payload", {})
        .get("successor_boundary", {})
        .get("TDG2_absolute_tail_discriminator_design_may_begin")
        is not True
    ):
        raise ValueError("TDG2 source PREF16 boundary differs")
    checkpoint = REPOSITORY / config["source_checkpoint"]
    _verify_optional_raw_hash(checkpoint, lineage["source_checkpoint_sha256"])
    return {
        "checkpoint_commit": commit,
        "checkpoint_commit_is_ancestor_of_HEAD": True,
        "tracked_blob_sha256": blobs,
        "source_checkpoint_sha256": lineage["source_checkpoint_sha256"],
        "source_checkpoint_present_and_hash_verified": checkpoint.is_file(),
        "source_checkpoint_history_arrays_loaded": False,
        "source_checkpoint_resumed": False,
    }


def _expect_config_rejection(
    config: Mapping[str, Any],
    mutate: Callable[[dict[str, Any]], None],
) -> bool:
    attacked = deepcopy(config)
    mutate(attacked)
    try:
        validate_config_data(attacked)
    except (TypeError, ValueError):
        return True
    return False


def _one_bit_evidence_control() -> bool:
    proper = np.linspace(0.0, 63.0 / 16.0, TDG2_SAMPLE_COUNT, dtype=np.float64)
    index = np.arange(TDG2_SAMPLE_COUNT, dtype=np.float64)
    values = 1.0e-3 * np.sin(2.0 * np.pi * 29.0 * index / TDG2_SAMPLE_COUNT)
    original = absolute_tail_observation(
        proper,
        values,
        round_trip_interpolation_infinity=0.0,
    )
    bits = values.view(np.uint64).copy()
    bits[17] ^= np.uint64(1 << 44)
    mutated = absolute_tail_observation(
        proper,
        bits.view(np.float64),
        round_trip_interpolation_infinity=0.0,
    )
    return (
        original.field_tail_power.observed
        != mutated.field_tail_power.observed
        and original.derivative_tail_power.observed
        != mutated.derivative_tail_power.observed
    )


def _mutation_controls(config: Mapping[str, Any]) -> dict[str, bool]:
    controls = {
        "one_bit_signal_evidence_mutation_detected": _one_bit_evidence_control(),
        "zero_power_order_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["frozen_discriminator"].__setitem__(
                "minimum_zero_power_order", "2"
            ),
        ),
        "decimal_precision_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["arithmetic_contract"].__setitem__(
                "high_decimal_precision", 80
            ),
        ),
        "method_ladder_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["frozen_discriminator"].__setitem__(
                "primary_point_counts", [1025, 2049, 4097]
            ),
        ),
        "interpolation_theorem_promotion_rejected": _expect_config_rejection(
            config,
            lambda value: value["interpolation_contract"].__setitem__(
                "round_trip_debit_is_not_claimed_as_continuum_interpolation_theorem",
                False,
            ),
        ),
        "source_hash_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["immutable_lineage"].__setitem__(
                "source_result_sha256", "0" * 64
            ),
        ),
        "raw_history_read_promotion_rejected": _expect_config_rejection(
            config,
            lambda value: value["scope"].__setitem__(
                "source_checkpoint_history_arrays_loaded", True
            ),
        ),
        "replacement_gate_promotion_rejected": _expect_config_rejection(
            config,
            lambda value: value["claims"].__setitem__(
                "replacement_temporal_admission_defined", True
            ),
        ),
        "GR0_eligibility_promotion_rejected": _expect_config_rejection(
            config,
            lambda value: value["claims"].__setitem__("GR0_case_eligible", True),
        ),
    }
    if not all(controls.values()):
        raise RuntimeError("TDG2 mutation controls did not all fail closed")
    return controls


def _compact_synthetic(value: Mapping[str, Any]) -> dict[str, Any]:
    controls: dict[str, Any] = {}
    for name, record in value["controls"].items():
        diagnosis = record["diagnosis"]
        controls[name] = {
            "expected_classification": record["expected_classification"],
            "observed_classifications": record["observed_classifications"],
            "classification_matched": record["classification_matched"],
            "round_trip_interpolation_debits": [
                item["round_trip_interpolation_infinity"]
                for item in diagnosis["observations"]
            ],
            "field_observed_powers": list(
                diagnosis["classifications"]["field_tail_power"]
                ["observed_powers"]
            ),
            "derivative_observed_powers": list(
                diagnosis["classifications"]["derivative_tail_power"]
                ["observed_powers"]
            ),
            "field_zero_power_order_lower_bounds": list(
                diagnosis["classifications"]["field_tail_power"]
                ["zero_power_order_lower_bounds"]
            ),
            "field_difference_order_lower_bound": diagnosis["classifications"]
            ["field_tail_power"]["difference_order_lower_bound"],
            "classification_is_not_a_replacement_admission": diagnosis[
                "classification_is_not_a_replacement_admission"
            ],
        }
    return {
        "controls": controls,
        "all_control_classes_separated": value["all_control_classes_separated"],
        "actual_terminal_histories_consumed": value[
            "actual_terminal_histories_consumed"
        ],
        "replacement_temporal_admission_defined": value[
            "replacement_temporal_admission_defined"
        ],
        "mathematical_consequence": value["mathematical_consequence"],
    }


def record(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    config = load_config(config_path)
    lineage = _verify_lineage(config)
    synthetic = synthetic_absolute_tail_controls()
    if synthetic["all_control_classes_separated"] is not True:
        raise RuntimeError("TDG2 synthetic discriminator did not reproduce")
    mutations = _mutation_controls(config)
    payload = {
        "immutable_lineage": lineage,
        "frozen_discriminator": {
            **config["frozen_discriminator"],
            "arithmetic_contract": config["arithmetic_contract"],
            "interpolation_contract": config["interpolation_contract"],
            "classification_contract": config["classification_contract"],
            "discriminator_module_sha256": _sha(DISCRIMINATOR_MODULE),
        },
        "synthetic_preflight": _compact_synthetic(synthetic),
        "mutation_controls": mutations,
        "claim_boundary": {
            "actual_terminal_histories_diagnosed": False,
            "historical_PROTO13_result_reclassified": False,
            "replacement_temporal_admission_defined": False,
            "new_trajectory_authorized": False,
            "candidate_branch_opened": False,
            "physical_or_candidate_question_answered": False,
        },
        "successor_boundary": config["successor_boundary"],
    }
    return {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "metric_signature": "-+++",
        "riemann_convention": "plus_partial_mu_gamma_nu",
        "classification": "prospective_TDG2_absolute_tail_discriminator_freeze",
        "artifact_payload": payload,
        "gate_status": config["claims"],
        "nonclaims": {
            key: value for key, value in config["claims"].items() if value is False
        },
        "source_config_sha256": {_rel(config_path): _sha(config_path)},
        "implementation_sha256": {_rel(path): _sha(path) for path in IMPLEMENTATION},
        "derivation_document": _rel(OWNER_DOCUMENT),
        "derivation_document_sha256": (
            _sha(OWNER_DOCUMENT) if OWNER_DOCUMENT.is_file() else None
        ),
        "generated_by": _rel(Path(__file__)),
    }


def _validate_stored_record(
    value: Mapping[str, Any], config: Mapping[str, Any]
) -> None:
    if (
        value.get("artifact_id") != ARTIFACT_ID
        or value.get("classification")
        != "prospective_TDG2_absolute_tail_discriminator_freeze"
        or value.get("gate_status") != EXPECTED_CLAIMS
        or value.get("artifact_payload", {})
        .get("synthetic_preflight", {})
        .get("all_control_classes_separated")
        is not True
        or value.get("artifact_payload", {})
        .get("immutable_lineage", {})
        .get("source_checkpoint_history_arrays_loaded")
        is not False
        or value.get("artifact_payload", {})
        .get("claim_boundary", {})
        .get("new_trajectory_authorized")
        is not False
        or value.get("source_config_sha256")
        != {_rel(DEFAULT_CONFIG): _sha(DEFAULT_CONFIG)}
    ):
        raise ValueError("TDG2 stored freeze claim boundary differs")
    if config["claims"] != EXPECTED_CLAIMS:
        raise ValueError("TDG2 stored freeze config claims differ")


def verify_canonical(
    config_path: Path = DEFAULT_CONFIG,
    output_path: Path = DEFAULT_OUTPUT,
) -> dict[str, Any]:
    config = load_config(config_path)
    observed = record(config_path)
    raw = output_path.read_bytes()
    stored = json.loads(raw)
    if not isinstance(stored, dict) or raw != _canonical_bytes(stored):
        raise ValueError("TDG2 stored freeze is not canonical JSON")
    _validate_stored_record(stored, config)
    if stored != observed:
        raise ValueError("TDG2 stored freeze differs from reproduction")
    return observed


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    if args.write and args.verify:
        raise SystemExit("choose either --write or --verify")
    value = record(args.config)
    if args.write:
        _atomic_write(args.output, _canonical_bytes(value))
        print(f"wrote {_rel(args.output)}")
    elif args.verify:
        verify_canonical(args.config, args.output)
        print(f"PASS {ARTIFACT_ID}")
    else:
        print(_canonical(value), end="")


if __name__ == "__main__":
    main()
