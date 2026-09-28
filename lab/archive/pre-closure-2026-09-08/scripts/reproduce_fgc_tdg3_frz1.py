#!/usr/bin/env python3
"""Reproduce the prospective FGC-1-TDG3-FRZ1 native-grid freeze."""

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

from recursive_horizons.fgc.evolution.tdg3_native_tail_diagnosis import (  # noqa: E402
    TDG3_COMBINED_CLASS_NAMES,
    TDG3_ESTIMATOR_NAMES,
    TDG3_FIELD_NAMES,
    TDG3_MEASURE_NAMES,
    TDG3_NORMALIZED_COEFFICIENT_FLOOR,
    TDG3_POWER_CLASS_NAMES,
    TDG3_SAMPLE_COUNT,
    TDG3_SYNTHETIC_CONTROL_NAMES,
    TDG3_TOP_BINS,
    native_tail_observation,
    synthetic_native_tail_controls,
)


ARTIFACT_ID = "FGC-1-TDG3-FRZ1"
PROJECT_VERSION = "0.11.0"
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-tdg3-frz1.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-tdg3-frz1.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-tdg3-frz1.md"
SOURCE_RESULT = REPOSITORY / "results/fgc-1-tdg2-pref17.json"
SOURCE_CONFIG = REPOSITORY / "configs/fgc/fgc-1-tdg2-pref17.toml"
SOURCE_REPRODUCER = REPOSITORY / "scripts/reproduce_fgc_tdg2_pref17.py"
SOURCE_DISCRIMINATOR = (
    REPOSITORY
    / "src/recursive_horizons/fgc/evolution/tdg2_absolute_tail_diagnosis.py"
)
SOURCE_DOCUMENT = REPOSITORY / "docs/fgc-tdg2-pref17.md"
SOURCE_CHECKPOINT = (
    REPOSITORY / "runs/fgc-2-sf1/proto13/calibration/latest-checkpoint.npz"
)
DISCRIMINATOR_MODULE = (
    REPOSITORY
    / "src/recursive_horizons/fgc/evolution/tdg3_native_tail_diagnosis.py"
)
IMPLEMENTATION = (Path(__file__).resolve(), DISCRIMINATOR_MODULE)


EXPECTED_CLAIMS = {
    "TDG2_actual_histories_diagnosed": True,
    "TDG2_absolute_tail_outcome_is_mixed": True,
    "TDG3_interpolation_ownership_discriminator_frozen": True,
    "TDG3_actual_histories_diagnosed": False,
    "common_grid_interpolation_removed_from_TDG3_measurement": True,
    "piecewise_linear_or_quadrature_surrogate_is_a_continuum_enclosure": False,
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


def _expected_lineage() -> dict[str, Any]:
    return {
        "checkpoint_commit": "8585dbef39a1e706390a17a002cc403445e18cf2",
        "source_result_sha256": "1a0d9ee99b2eb2f15f10b0625c1cbfe63db15d26692ff5f7b1975d370db593a3",
        "source_config_sha256": "3f2ec6acf07779c4ae0d5b47fe235f66724d4f3b86975eb391dfae02d325c000",
        "source_reproducer_sha256": "f0a939c5f46a0e4fb57b9118e6107788c5f585555802a07a5e8324a63111e776",
        "source_discriminator_module_sha256": "d8332ab5b6a3d3ba9db3a129e71a4147c88a3da5c89624914f2f9c757ac74595",
        "source_document_sha256": "8495b96c84bb5cba4af869701c732fd40793b0cc9b951990f1a71a64b129cd01",
        "source_checkpoint_sha256": "0432b266081b1ec1af3ac2a9efb99bcdd67dc8afef837241ec27fe5d37432f97",
        "checkpoint_commit_must_be_ancestor_of_HEAD": True,
        "tracked_source_blobs_must_match_commit_and_worktree": True,
        "source_compact_result_must_be_verified_without_raw_reproduction": True,
        "source_checkpoint_hash_must_match_if_present": True,
    }


def validate_config_data(config: Mapping[str, Any]) -> None:
    _strict_keys(
        "TDG3 freeze config",
        config,
        {
            "schema_version",
            "artifact_id",
            "project_version",
            "metric_signature",
            "riemann_convention",
            "source_result",
            "source_checkpoint",
            "scope",
            "immutable_lineage",
            "frozen_discriminator",
            "synthetic_controls",
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
        or config.get("source_checkpoint") != _rel(SOURCE_CHECKPOINT)
    ):
        raise ValueError("TDG3 freeze identity differs")
    scope = config.get("scope", {})
    if scope != {
        "target": "FGC-2-SF1-TDG3",
        "role": "prospective_native_grid_interpolation_ownership_discriminator",
        "calibration_branch": "GR-0",
        "source_terminal_event": 63,
        "source_terminal_coordinate_time": "63/16",
        "actual_terminal_histories_consumed": False,
        "source_checkpoint_history_arrays_loaded": False,
        "source_checkpoint_resumed": False,
        "state_advanced": False,
        "common_grid_resampling_performed": False,
        "SGBL_trajectory_read": False,
        "FGCQR_trajectory_read": False,
        "replacement_temporal_admission_defined": False,
        "physical_or_candidate_question_answered": False,
    }:
        raise ValueError("TDG3 freeze scope differs")
    if config.get("immutable_lineage") != _expected_lineage():
        raise ValueError("TDG3 freeze lineage differs")
    frozen = config.get("frozen_discriminator", {})
    expected_frozen_keys = {
        "sample_count",
        "tracer_count",
        "field_names",
        "primary_point_counts",
        "comparator_point_counts",
        "native_phase_interval",
        "phase_map",
        "common_proper_time_overlap_used",
        "common_grid_resampling_allowed",
        "taper_fraction",
        "top_bins",
        "power_measures",
        "estimator_names",
        "normalized_coefficient_floor",
        "minimum_history_difference_order",
        "minimum_zero_power_order",
        "replacement_admission_authorized",
        "primary_estimator",
        "comparator_estimator",
        "arithmetic_contract",
        "classification_contract",
    }
    _strict_keys("TDG3 frozen discriminator", frozen, expected_frozen_keys)
    if (
        frozen["sample_count"] != TDG3_SAMPLE_COUNT
        or frozen["tracer_count"] != 48
        or tuple(frozen["field_names"]) != tuple(TDG3_FIELD_NAMES)
        or frozen["primary_point_counts"] != [2049, 4097, 8193]
        or frozen["comparator_point_counts"] != [4097, 8193, 16385]
        or frozen["native_phase_interval"] != [0, 1]
        or frozen["phase_map"]
        != "(tau-tau_0)/(tau_63-tau_0)_per_native_history"
        or frozen["common_proper_time_overlap_used"] is not False
        or frozen["common_grid_resampling_allowed"] is not False
        or frozen["taper_fraction"] != "1/8"
        or tuple(frozen["top_bins"]) != TDG3_TOP_BINS
        or tuple(frozen["power_measures"]) != TDG3_MEASURE_NAMES
        or tuple(frozen["estimator_names"]) != TDG3_ESTIMATOR_NAMES
        or float(frozen["normalized_coefficient_floor"])
        != TDG3_NORMALIZED_COEFFICIENT_FLOOR
        or frozen["minimum_history_difference_order"] != "3/2"
        or frozen["minimum_zero_power_order"] != "3"
        or frozen["replacement_admission_authorized"] is not False
    ):
        raise ValueError("TDG3 frozen measurement differs")
    primary = frozen["primary_estimator"]
    comparator = frozen["comparator_estimator"]
    arithmetic = frozen["arithmetic_contract"]
    classification = frozen["classification_contract"]
    if (
        primary.get("name") != TDG3_ESTIMATOR_NAMES[0]
        or primary.get("resampled_values_allowed") is not False
        or primary.get("continuum_enclosure_claimed") is not False
        or comparator.get("name") != TDG3_ESTIMATOR_NAMES[1]
        or comparator.get("resampled_values_allowed") is not False
        or comparator.get("continuum_enclosure_claimed") is not False
        or arithmetic.get("low_decimal_precision") != 72
        or arithmetic.get("high_decimal_precision") != 96
        or arithmetic.get("arithmetic_debit_is_the_only_interval_debit")
        is not True
        or arithmetic.get("inter_estimator_spread_is_public_but_not_an_error_enclosure")
        is not True
        or tuple(classification.get("power_class_names", ()))
        != TDG3_POWER_CLASS_NAMES
        or tuple(classification.get("combined_class_names", ()))
        != TDG3_COMBINED_CLASS_NAMES
        or classification.get("no_class_is_a_pass_fail_temporal_admission")
        is not True
        or classification.get("no_class_retroactively_reclassifies_PROTO13_or_TDG2")
        is not True
    ):
        raise ValueError("TDG3 estimator or classification contract differs")
    synthetic = config.get("synthetic_controls", {})
    if (
        tuple(synthetic.get("control_names", ()))
        != TDG3_SYNTHETIC_CONTROL_NAMES
        or set(
            value
            for key, value in synthetic.items()
            if key != "control_names"
        )
        != {True, False}
        or synthetic.get("every_control_must_match_both_power_measures")
        is not True
        or synthetic.get("actual_terminal_histories_consumed") is not False
        or synthetic.get("common_grid_resampling_performed") is not False
        or synthetic.get("continuum_surrogate_enclosure_claimed") is not False
        or synthetic.get("replacement_temporal_admission_defined") is not False
    ):
        raise ValueError("TDG3 synthetic contract differs")
    if set(config.get("proof_contract", {}).values()) != {True}:
        raise ValueError("TDG3 proof contract must be all-of")
    if config.get("successor_boundary") != {
        "TDG2_actual_histories_diagnosed": True,
        "TDG3_interpolation_ownership_discriminator_frozen": True,
        "TDG3_actual_history_execution_authorized_after_freeze_commit": True,
        "TDG3_actual_history_execution_completed": False,
        "replacement_temporal_gate_design_authorized": False,
        "replacement_temporal_admission_defined": False,
        "PROTO14_frozen": False,
        "fresh_GR0_dynamic_calibration_completed": False,
        "GR0_case_eligible": False,
        "SGBL_execution_authorized": False,
        "FGCQR_holdout_execution_authorized": False,
        "DEF1_execution_authorized": False,
    }:
        raise ValueError("TDG3 successor boundary differs")
    if config.get("claims") != EXPECTED_CLAIMS:
        raise ValueError("TDG3 claims differ")


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
        raise ValueError("TDG3 checkpoint commit is not an ancestor of HEAD")
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
            raise ValueError(f"TDG3 source tracked blob differs: {relative}")
    raw = SOURCE_RESULT.read_bytes()
    source = json.loads(raw)
    compact = source.get("artifact_payload", {}).get("compact_findings", {})
    if (
        raw != _canonical_bytes(source)
        or source.get("artifact_id") != "FGC-1-TDG2-PREF17"
        or source.get("classification")
        != "completed_TDG2_mixed_absolute_tail_diagnosis_TDG3_required"
        or source.get("gate_status", {}).get("TDG2_actual_histories_diagnosed")
        is not True
        or source.get("gate_status", {}).get("replacement_temporal_admission_defined")
        is not False
        or compact.get("total_signal_ladders") != 576
        or compact.get("power_classifications") != 1152
        or compact.get("enclosure_dominated_finest_failures_field") != 101
        or compact.get("enclosure_dominated_finest_failures_derivative") != 102
    ):
        raise ValueError("TDG3 source compact result differs")
    checkpoint_present = SOURCE_CHECKPOINT.is_file()
    if checkpoint_present and _sha(SOURCE_CHECKPOINT) != lineage["source_checkpoint_sha256"]:
        raise ValueError("TDG3 optional checkpoint hash differs")
    return {
        "checkpoint_commit": commit,
        "checkpoint_commit_is_ancestor_of_HEAD": True,
        "tracked_blob_sha256": tracked,
        "source_compact_result_verified_without_raw_reproduction": True,
        "source_diagnostic_record_read": True,
        "source_checkpoint_present": checkpoint_present,
        "source_checkpoint_sha256": lineage["source_checkpoint_sha256"],
        "source_checkpoint_history_arrays_loaded": False,
        "source_checkpoint_resumed": False,
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
    config: Mapping[str, Any], synthetic: Mapping[str, Any]
) -> dict[str, bool]:
    coordinate = np.linspace(0.0, 63.0 / 16.0, TDG3_SAMPLE_COUNT)
    phase = (coordinate - coordinate[0]) / (coordinate[-1] - coordinate[0])
    signal = 1.0e-3 * np.sin(2.0 * np.pi * 29.0 * phase)
    baseline = native_tail_observation(
        coordinate, signal, estimator=TDG3_ESTIMATOR_NAMES[0]
    )
    bits = signal.view(np.uint64).copy()
    bits[17] ^= np.uint64(1 << 44)
    mutated = native_tail_observation(
        coordinate,
        bits.view(np.float64),
        estimator=TDG3_ESTIMATOR_NAMES[0],
    )
    unknown_estimator_rejected = False
    try:
        native_tail_observation(coordinate, signal, estimator="cubic_spline")
    except (TypeError, ValueError):
        unknown_estimator_rejected = True
    controls = {
        "one_bit_signal_mutation_detected": baseline != mutated,
        "unknown_estimator_rejected": unknown_estimator_rejected,
        "estimator_config_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["frozen_discriminator"].__setitem__(
                "estimator_names", ["piecewise_linear_exact_integral"]
            ),
        ),
        "top_bin_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["frozen_discriminator"].__setitem__(
                "top_bins", [27, 28, 29, 30, 31]
            ),
        ),
        "minimum_order_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["frozen_discriminator"].__setitem__(
                "minimum_zero_power_order", "2"
            ),
        ),
        "decimal_precision_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["frozen_discriminator"][
                "arithmetic_contract"
            ].__setitem__("high_decimal_precision", 80),
        ),
        "source_hash_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["immutable_lineage"].__setitem__(
                "source_result_sha256", "0" * 64
            ),
        ),
        "resampling_promotion_rejected": _expect_config_rejection(
            config,
            lambda value: value["scope"].__setitem__(
                "common_grid_resampling_performed", True
            ),
        ),
        "continuum_enclosure_promotion_rejected": _expect_config_rejection(
            config,
            lambda value: value["claims"].__setitem__(
                "piecewise_linear_or_quadrature_surrogate_is_a_continuum_enclosure",
                True,
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
        "synthetic_preflight_complete": synthetic.get(
            "all_control_classes_separated"
        )
        is True,
    }
    if set(controls.values()) != {True}:
        raise ValueError(f"TDG3 mutation control failed: {controls}")
    return controls


def record(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    config = load_config(config_path)
    lineage = _verify_source_lineage(config)
    synthetic = synthetic_native_tail_controls()
    if (
        synthetic.get("all_control_classes_separated") is not True
        or synthetic.get("actual_terminal_histories_consumed") is not False
        or synthetic.get("common_grid_resampling_performed") is not False
        or synthetic.get("continuum_surrogate_enclosure_claimed") is not False
        or synthetic.get("replacement_temporal_admission_defined") is not False
        or synthetic.get("common_grid_interpolation_owner_control", {}).get(
            "TDG2_identifies_interpolation_owner"
        )
        is not True
        or synthetic.get("common_grid_interpolation_owner_control", {}).get(
            "TDG3_native_estimators_agree_nonzero"
        )
        is not True
    ):
        raise ValueError("TDG3 synthetic preflight differs")
    mutations = _mutation_controls(config, synthetic)
    result = {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "metric_signature": "-+++",
        "riemann_convention": "plus_partial_mu_gamma_nu",
        "classification": "prospective_TDG3_native_grid_discriminator_freeze",
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
            "frozen_discriminator": deepcopy(config["frozen_discriminator"]),
            "synthetic_preflight": synthetic,
            "mutation_controls": mutations,
            "proof_contract": dict(config["proof_contract"]),
            "successor_boundary": dict(config["successor_boundary"]),
            "claim_boundary": {
                "actual_terminal_histories_diagnosed": False,
                "common_grid_resampling_performed": False,
                "continuum_surrogate_enclosure_claimed": False,
                "historical_PROTO13_or_TDG2_result_reclassified": False,
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
        != "prospective_TDG3_native_grid_discriminator_freeze"
        or value.get("generated_by") != "scripts/reproduce_fgc_tdg3_frz1.py"
        or value.get("gate_status") != EXPECTED_CLAIMS
        or value.get("nonclaims")
        != {
            key: False
            for key, item in sorted(EXPECTED_CLAIMS.items())
            if item is False
        }
        or value.get("source_config_sha256")
        != {_rel(DEFAULT_CONFIG): _sha(DEFAULT_CONFIG)}
        or payload.get("frozen_discriminator") != config["frozen_discriminator"]
        or payload.get("successor_boundary") != config["successor_boundary"]
    ):
        raise ValueError("TDG3 stored result boundary differs")
    implementation = value.get("implementation_sha256", {})
    if (
        not isinstance(implementation, Mapping)
        or set(implementation) != {_rel(path) for path in IMPLEMENTATION}
        or any(
            _sha(REPOSITORY / relative) != expected
            for relative, expected in implementation.items()
        )
    ):
        raise ValueError("TDG3 implementation hash ledger differs")
    document = REPOSITORY / str(value.get("derivation_document", ""))
    if (
        value.get("derivation_document") != "docs/fgc-tdg3-frz1.md"
        or not document.is_file()
        or _sha(document) != value.get("derivation_document_sha256")
    ):
        raise ValueError("TDG3 derivation document differs")
    synthetic = payload.get("synthetic_preflight", {})
    if (
        synthetic.get("all_control_classes_separated") is not True
        or set(synthetic.get("controls", {})) != set(TDG3_SYNTHETIC_CONTROL_NAMES)
        or synthetic.get("actual_terminal_histories_consumed") is not False
    ):
        raise ValueError("TDG3 synthetic evidence differs")
    mutations = payload.get("mutation_controls", {})
    if not isinstance(mutations, Mapping) or set(mutations.values()) != {True}:
        raise ValueError("TDG3 mutation controls differ")


def verify_canonical(
    config_path: Path = DEFAULT_CONFIG,
    output_path: Path = DEFAULT_OUTPUT,
) -> dict[str, Any]:
    config = load_config(config_path)
    raw = output_path.read_bytes()
    stored = json.loads(raw)
    if not isinstance(stored, dict) or raw != _canonical_bytes(stored):
        raise ValueError("TDG3 stored result is not canonical JSON")
    _validate_stored_record(stored, config)
    observed = record(config_path)
    if stored != observed:
        raise ValueError("TDG3 stored result differs from reproduction")
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
            f"PASS {ARTIFACT_ID}: synthetic="
            f"{value['artifact_payload']['synthetic_preflight']['all_control_classes_separated']}"
        )
    else:
        print(_canonical(record(args.config)), end="")


if __name__ == "__main__":
    main()
