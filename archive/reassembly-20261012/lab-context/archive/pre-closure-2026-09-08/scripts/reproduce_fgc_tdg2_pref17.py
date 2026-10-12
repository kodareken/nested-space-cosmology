#!/usr/bin/env python3
"""Bind TDG2's read-only diagnosis of the immutable PROTO13 histories."""

from __future__ import annotations

import argparse
from collections import Counter
from copy import deepcopy
from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import tomllib
from typing import Any, Callable, Mapping, Sequence

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts import reproduce_fgc_cal10_pref15 as cal10  # noqa: E402
from recursive_horizons.fgc.evolution.tdg1_temporal_diagnosis import (  # noqa: E402
    TDG1_FIELD_NAMES,
    temporal_signal_views,
)
from recursive_horizons.fgc.evolution.tdg2_absolute_tail_diagnosis import (  # noqa: E402
    TDG2_CLASS_NAMES,
    TDG2_MEASURE_NAMES,
    TDG2_METHOD_POINT_COUNTS,
    diagnose_absolute_tail_ladder,
)


ARTIFACT_ID = "FGC-1-TDG2-PREF17"
PROJECT_VERSION = "0.11.0"
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-tdg2-pref17.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-tdg2-pref17.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-tdg2-pref17.md"
FREEZE_RESULT = REPOSITORY / "results/fgc-1-tdg2-frz1.json"
FREEZE_CONFIG = REPOSITORY / "configs/fgc/fgc-1-tdg2-frz1.toml"
FREEZE_REPRODUCER = REPOSITORY / "scripts/reproduce_fgc_tdg2_frz1.py"
DISCRIMINATOR_MODULE = (
    REPOSITORY
    / "src/recursive_horizons/fgc/evolution/tdg2_absolute_tail_diagnosis.py"
)
FREEZE_DOCUMENT = REPOSITORY / "docs/fgc-tdg2-frz1.md"
LEGACY_DIAGNOSIS_MODULE = (
    REPOSITORY
    / "src/recursive_horizons/fgc/evolution/tdg1_temporal_diagnosis.py"
)
CAL10_REPRODUCER = REPOSITORY / "scripts/reproduce_fgc_cal10_pref15.py"
IMPLEMENTATION = (Path(__file__).resolve(),)
METHOD_KEYS = {
    "RK4": ("RK4-2049", "RK4-4097", "RK4-8193"),
    "SSPRK3": ("SSPRK3-4097", "SSPRK3-8193", "SSPRK3-16385"),
}


EXPECTED_CLAIMS = {
    "TDG1_terminal_histories_diagnosed": True,
    "normalized_raw_temporal_tail_ratio_is_valid_general_spatial_convergence_test": False,
    "TDG2_absolute_tail_discriminator_frozen": True,
    "TDG2_actual_histories_diagnosed": True,
    "TDG2_absolute_tail_outcome_is_mixed": True,
    "individual_temporal_budget_failure_cause_fully_derived": False,
    "replacement_temporal_admission_defined": False,
    "TDG3_frozen": False,
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


def _as_float(value: str) -> float:
    numerator, separator, denominator = value.partition("/")
    return float(int(numerator) / int(denominator)) if separator else float(value)


def _expected_lineage() -> dict[str, Any]:
    return {
        "freeze_commit": "a1b6b99c625ee04e4c8d931545c471e56198d540",
        "freeze_result_sha256": "4e35f8cf3b23970e4031e6fbc7ec7322ac802dfb4cca7982d968fbd71e3b80bb",
        "freeze_config_sha256": "9b3472b54819109ba485f654ff5915c8d96c6fc2d9d14496d6db97985568a179",
        "freeze_reproducer_sha256": "2544cca6f4b4a0a541aade3d0dd4369636e5c163a798d28b7cad180e67adb46b",
        "discriminator_module_sha256": "d8332ab5b6a3d3ba9db3a129e71a4147c88a3da5c89624914f2f9c757ac74595",
        "freeze_document_sha256": "0c0030ec93202641331836dee1d3263c9ebe2440cc378805254b4e25ef53839d",
        "legacy_diagnosis_module_sha256": "3aaf100a8c9abb691978ac3ddf1f58210ed59107dc6018cf3f8b70c99de8ba44",
        "CAL10_reproducer_sha256": "71e1e4d211e46d9679e9772eefb0ab90440a648c009efbe570716b5dcf2659e6",
        "manifest_sha256": "7fbe3835da9d23ed868bf5862250a9c9cad31bfb948091fc1290454867f084f6",
        "event_log_sha256": "cc26095e48de376f0e7924a8b5f2c5126e3a0f21f6d6036942cfc4031d248bdd",
        "checkpoint_sha256": "0432b266081b1ec1af3ac2a9efb99bcdd67dc8afef837241ec27fe5d37432f97",
        "campaign_result_sha256": "83475a061e2dd4e354e357a3c0cca66e44c60312f9c420585dd989598779aa0b",
        "freeze_commit_must_be_ancestor_of_HEAD": True,
        "tracked_freeze_blobs_must_match_commit_and_worktree": True,
        "raw_bundle_must_be_complete_and_hash_exact_when_present": True,
        "clean_clone_may_verify_only_the_compact_result": True,
    }


def _expected_observed() -> dict[str, Any]:
    return {
        "method_count": 2,
        "signal_ladders_per_method": 288,
        "total_signal_ladders": 576,
        "resolution_observations": 1728,
        "power_classifications": 1152,
        "legacy_finest_individual_failures": 152,
        "legacy_finest_above_floor_failures": 103,
        "legacy_persistent_all_resolution_failures": 148,
        "all_three_below_DFT_floor_ladders": 96,
        "nonzero_interpolation_debit_ladders": 474,
        "zero_contracting_finest_failures_field": 46,
        "nonzero_resolution_independent_finest_failures_field": 5,
        "enclosure_dominated_finest_failures_field": 101,
        "zero_contracting_finest_failures_derivative": 46,
        "nonzero_resolution_independent_finest_failures_derivative": 4,
        "enclosure_dominated_finest_failures_derivative": 102,
        "nonzero_convergent_classifications": 0,
        "unresolved_classifications": 0,
        "replacement_temporal_admission_earned": False,
        "RK4": {
            "point_counts": [2049, 4097, 8193],
            "field_zero_contracting": 35,
            "field_nonzero_resolution_independent": 94,
            "field_enclosure_dominated": 159,
            "derivative_zero_contracting": 35,
            "derivative_nonzero_resolution_independent": 94,
            "derivative_enclosure_dominated": 159,
            "finest_individual_failures": 93,
            "finest_above_floor_failures": 66,
            "persistent_all_resolution_failures": 89,
            "finest_failure_zero_contracting": 35,
            "finest_failure_nonzero_resolution_independent": 4,
            "finest_failure_enclosure_dominated": 54,
            "finest_failure_all_three_DFT_floor": 4,
            "finest_failure_interpolation_interval_reaches_zero": 47,
            "finest_failure_cross_resolution_intervals_overlap": 3,
        },
        "SSPRK3": {
            "point_counts": [4097, 8193, 16385],
            "field_zero_contracting": 11,
            "field_nonzero_resolution_independent": 91,
            "field_enclosure_dominated": 186,
            "derivative_zero_contracting": 11,
            "derivative_nonzero_resolution_independent": 90,
            "derivative_enclosure_dominated": 187,
            "finest_individual_failures": 59,
            "finest_above_floor_failures": 37,
            "persistent_all_resolution_failures": 59,
            "field_finest_failure_zero_contracting": 11,
            "field_finest_failure_nonzero_resolution_independent": 1,
            "field_finest_failure_enclosure_dominated": 47,
            "derivative_finest_failure_zero_contracting": 11,
            "derivative_finest_failure_nonzero_resolution_independent": 0,
            "derivative_finest_failure_enclosure_dominated": 48,
            "field_finest_failure_all_three_DFT_floor": 10,
            "field_finest_failure_interpolation_interval_reaches_zero": 37,
            "derivative_finest_failure_all_three_DFT_floor": 10,
            "derivative_finest_failure_interpolation_interval_reaches_zero": 38,
        },
        "cross_method": {
            "field_method_agreement_count": 259,
            "derivative_method_agreement_count": 258,
            "signal_count": 288,
            "field_both_methods_zero_contracting": 10,
            "field_both_methods_nonzero_resolution_independent": 91,
            "field_both_methods_enclosure_dominated": 158,
            "derivative_both_methods_zero_contracting": 10,
            "derivative_both_methods_nonzero_resolution_independent": 90,
            "derivative_both_methods_enclosure_dominated": 158,
        },
    }


def validate_config_data(config: Mapping[str, Any]) -> None:
    _strict_keys(
        "PREF17 config",
        config,
        {
            "schema_version",
            "artifact_id",
            "project_version",
            "metric_signature",
            "riemann_convention",
            "freeze_result",
            "freeze_config",
            "campaign_manifest",
            "campaign_event_log",
            "campaign_checkpoint",
            "campaign_result",
            "scope",
            "immutable_lineage",
            "observed_diagnosis",
            "conclusions",
            "proof_contract",
            "successor_boundary",
            "claims",
        },
    )
    expected_paths = {
        "freeze_result": _rel(FREEZE_RESULT),
        "freeze_config": _rel(FREEZE_CONFIG),
        "campaign_manifest": "runs/fgc-2-sf1/proto13/calibration/manifest.json",
        "campaign_event_log": "runs/fgc-2-sf1/proto13/calibration/events.jsonl",
        "campaign_checkpoint": "runs/fgc-2-sf1/proto13/calibration/latest-checkpoint.npz",
        "campaign_result": "runs/fgc-2-sf1/proto13/calibration/campaign-result.json",
    }
    if (
        config.get("schema_version") != 1
        or config.get("artifact_id") != ARTIFACT_ID
        or config.get("project_version") != PROJECT_VERSION
        or config.get("metric_signature") != "-+++"
        or config.get("riemann_convention") != "plus_partial_mu_gamma_nu"
        or any(config.get(key) != value for key, value in expected_paths.items())
    ):
        raise ValueError("PREF17 top-level contract differs")
    if config.get("scope") != {
        "target": "FGC-2-SF1-TDG2",
        "role": "post_freeze_read_only_absolute_tail_history_diagnosis",
        "calibration_branch": "GR-0",
        "source_terminal_event": 63,
        "source_terminal_coordinate_time": "63/16",
        "terminal_histories_read": True,
        "terminal_checkpoint_resumed": False,
        "state_advanced": False,
        "SGBL_trajectory_read": False,
        "FGCQR_trajectory_read": False,
        "replacement_temporal_admission_defined": False,
        "physical_or_candidate_question_answered": False,
    }:
        raise ValueError("PREF17 scope differs")
    if _as_float(config["scope"]["source_terminal_coordinate_time"]) != 63 / 16:
        raise ValueError("PREF17 source time differs")
    if config.get("immutable_lineage") != _expected_lineage():
        raise ValueError("PREF17 immutable lineage differs")
    if config.get("observed_diagnosis") != _expected_observed():
        raise ValueError("PREF17 observed diagnosis differs")
    if config.get("conclusions") != {
        "historical_PROTO13_failure_preserved": True,
        "TDG2_absolute_tail_diagnosis_completed": True,
        "absolute_tail_outcome_is_mixed": True,
        "some_legacy_failures_contract_toward_zero": True,
        "small_nonzero_subset_is_resolution_independent_within_conditional_enclosure": True,
        "majority_of_finest_legacy_failures_remain_enclosure_dominated": True,
        "interpolation_debit_is_the_dominant_unresolved_owner": True,
        "binary64_DFT_arithmetic_is_not_the_owner_of_any_finest_legacy_failure": True,
        "replacement_temporal_admission_not_yet_justified": True,
        "TDG3_interpolation_ownership_discriminator_design_may_begin": True,
        "PROTO14_design_may_not_yet_begin": True,
    }:
        raise ValueError("PREF17 conclusions differ")
    if set(config.get("proof_contract", {}).values()) != {True}:
        raise ValueError("PREF17 proof contract must be all-of")
    if config.get("successor_boundary") != {
        "TDG2_absolute_tail_discriminator_frozen": True,
        "TDG2_actual_histories_diagnosed": True,
        "TDG3_interpolation_ownership_discriminator_design_may_begin": True,
        "TDG3_frozen": False,
        "replacement_temporal_gate_design_authorized": False,
        "replacement_temporal_admission_defined": False,
        "PROTO14_frozen": False,
        "fresh_GR0_dynamic_calibration_completed": False,
        "GR0_case_eligible": False,
        "SGBL_execution_authorized": False,
        "FGCQR_holdout_execution_authorized": False,
        "DEF1_execution_authorized": False,
    }:
        raise ValueError("PREF17 successor boundary differs")
    if config.get("claims") != EXPECTED_CLAIMS:
        raise ValueError("PREF17 claims differ")


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    with path.open("rb") as handle:
        config = tomllib.load(handle)
    validate_config_data(config)
    return config


def _raw_paths(config: Mapping[str, Any]) -> tuple[Path, Path, Path, Path]:
    return tuple(
        REPOSITORY / config[name]
        for name in (
            "campaign_manifest",
            "campaign_event_log",
            "campaign_checkpoint",
            "campaign_result",
        )
    )


def _raw_bundle_is_complete_or_absent(paths: Sequence[Path]) -> bool:
    present = [path.is_file() for path in paths]
    if any(present) and not all(present):
        raise ValueError("PREF17 raw campaign bundle is partial")
    return all(present)


def _verify_freeze_lineage(config: Mapping[str, Any]) -> dict[str, Any]:
    lineage = config["immutable_lineage"]
    commit = lineage["freeze_commit"]
    _git("cat-file", "-e", f"{commit}^{{commit}}")
    if _git("merge-base", "--is-ancestor", commit, "HEAD", check=False).returncode:
        raise ValueError("PREF17 freeze commit is not an ancestor of HEAD")
    tracked = {
        _rel(FREEZE_RESULT): lineage["freeze_result_sha256"],
        _rel(FREEZE_CONFIG): lineage["freeze_config_sha256"],
        _rel(FREEZE_REPRODUCER): lineage["freeze_reproducer_sha256"],
        _rel(DISCRIMINATOR_MODULE): lineage["discriminator_module_sha256"],
        _rel(FREEZE_DOCUMENT): lineage["freeze_document_sha256"],
        _rel(LEGACY_DIAGNOSIS_MODULE): lineage["legacy_diagnosis_module_sha256"],
        _rel(CAL10_REPRODUCER): lineage["CAL10_reproducer_sha256"],
    }
    for relative, expected_hash in tracked.items():
        committed = _git("show", f"{commit}:{relative}").stdout
        current = (REPOSITORY / relative).read_bytes()
        if (
            committed != current
            or _bytes_sha(committed) != expected_hash
            or _bytes_sha(current) != expected_hash
        ):
            raise ValueError(f"PREF17 frozen tracked blob differs: {relative}")
    freeze_raw = FREEZE_RESULT.read_bytes()
    freeze = json.loads(freeze_raw)
    if (
        freeze_raw != _canonical_bytes(freeze)
        or freeze.get("artifact_id") != "FGC-1-TDG2-FRZ1"
        or freeze.get("gate_status", {}).get("TDG2_actual_histories_diagnosed")
        is not False
        or freeze.get("artifact_payload", {})
        .get("successor_boundary", {})
        .get("TDG2_actual_history_execution_authorized")
        is not True
    ):
        raise ValueError("PREF17 source freeze boundary differs")
    return {
        "freeze_commit": commit,
        "freeze_commit_is_ancestor_of_HEAD": True,
        "tracked_blob_sha256": tracked,
        "source_freeze_authorized_actual_history_execution": True,
    }


def _history_inputs(
    histories: Mapping[str, Mapping[str, np.ndarray]],
) -> tuple[dict[str, list[np.ndarray]], dict[str, list[np.ndarray]]]:
    proper = {
        method: [histories[key]["event_proper_times"] for key in keys]
        for method, keys in METHOD_KEYS.items()
    }
    fields = {
        method: [histories[key]["event_fields"] for key in keys]
        for method, keys in METHOD_KEYS.items()
    }
    return proper, fields


def _classification_evidence(
    diagnosis: Mapping[str, Any], measure: str
) -> dict[str, Any]:
    record = diagnosis["classifications"][measure]
    return {
        "classification": record["classification"],
        "observed_powers": list(record["observed_powers"]),
        "lower_bounds": list(record["lower_bounds"]),
        "upper_bounds": list(record["upper_bounds"]),
        "zero_power_order_lower_bounds": list(
            record["zero_power_order_lower_bounds"]
        ),
        "difference_order_lower_bound": record["difference_order_lower_bound"],
        "conservative_continuum_lower_bound": record[
            "conservative_continuum_lower_bound"
        ],
        "every_resolution_below_declared_amplitude_floor": record[
            "every_resolution_below_declared_amplitude_floor"
        ],
        "enclosure_controls_classification": record[
            "enclosure_controls_classification"
        ],
        "replacement_admission_authorized": record[
            "replacement_admission_authorized"
        ],
    }


def _enclosure_cause(diagnosis: Mapping[str, Any], measure: str) -> str | None:
    classification = diagnosis["classifications"][measure]
    if (
        classification["classification"]
        != "instrument_floor_or_interpolation_enclosure_dominated"
    ):
        return None
    observations = diagnosis["observations"]
    if all(
        record["top_band_is_below_declared_amplitude_floor"]
        for record in observations
    ):
        return "all_three_DFT_floor"
    if classification["lower_bounds"][2] == 0.0:
        power = observations[2][measure]
        if power["interpolation_debit"] > power["arithmetic_debit"]:
            return "finest_interpolation_interval_reaches_zero"
        return "finest_arithmetic_interval_reaches_zero"
    return "cross_resolution_intervals_overlap"


def _signal_record(
    *,
    method: str,
    tracer: int,
    field_index: int,
    proper: Sequence[np.ndarray],
    fields: Sequence[np.ndarray],
) -> dict[str, Any]:
    field_name = TDG1_FIELD_NAMES[field_index]
    diagnosis = diagnose_absolute_tail_ladder(
        [record[:, tracer] for record in proper],
        [record[:, tracer, field_index] for record in fields],
        method=method,
        point_counts=TDG2_METHOD_POINT_COUNTS[method],
    )
    legacy_failed: list[bool] = []
    legacy_floor: list[bool] = []
    for times, values in zip(proper, fields, strict=True):
        view = temporal_signal_views(
            times[:, tracer], values[:, tracer, field_index]
        )["views"]["inherited_compact"]
        legacy_failed.append(not view["budget"]["individual_admission_passed"])
        legacy_floor.append(bool(view["below_declared_absolute_amplitude_floor"]))
    observations = []
    for count, record in zip(
        TDG2_METHOD_POINT_COUNTS[method], diagnosis["observations"], strict=True
    ):
        observations.append(
            {
                "point_count": count,
                "round_trip_interpolation_infinity": record[
                    "round_trip_interpolation_infinity"
                ],
                "top_band_is_below_declared_amplitude_floor": record[
                    "top_band_is_below_declared_amplitude_floor"
                ],
                "peak_top_band_coefficient_amplitude": record[
                    "peak_top_band_coefficient_amplitude"
                ],
                "binary64_peak_top_band_coefficient_amplitude": record[
                    "binary64_peak_top_band_coefficient_amplitude"
                ],
                "coefficient_arithmetic_debit": record[
                    "coefficient_arithmetic_debit"
                ],
                "coefficient_interpolation_debit": record[
                    "coefficient_interpolation_debit"
                ],
                "field_tail_power": record["field_tail_power"],
                "derivative_tail_power": record["derivative_tail_power"],
            }
        )
    return {
        "method": method,
        "tracer": tracer,
        "field": field_name,
        "legacy_individual_failure_pattern": "".join(
            "F" if value else "P" for value in legacy_failed
        ),
        "legacy_below_amplitude_floor_pattern": list(legacy_floor),
        "field_tail_power": _classification_evidence(
            diagnosis, "field_tail_power"
        ),
        "derivative_tail_power": _classification_evidence(
            diagnosis, "derivative_tail_power"
        ),
        "field_enclosure_cause": _enclosure_cause(
            diagnosis, "field_tail_power"
        ),
        "derivative_enclosure_cause": _enclosure_cause(
            diagnosis, "derivative_tail_power"
        ),
        "observations": observations,
    }


def _class_count(records: Sequence[Mapping[str, Any]], measure: str) -> Counter[str]:
    return Counter(record[measure]["classification"] for record in records)


def _method_findings(
    method: str, records: Sequence[Mapping[str, Any]]
) -> dict[str, Any]:
    if len(records) != 288:
        raise ValueError(f"PREF17 {method} signal count differs")
    field = _class_count(records, "field_tail_power")
    derivative = _class_count(records, "derivative_tail_power")
    finest_failures = [
        record
        for record in records
        if record["legacy_individual_failure_pattern"][-1] == "F"
    ]
    above_floor = [
        record
        for record in finest_failures
        if not record["legacy_below_amplitude_floor_pattern"][-1]
    ]
    persistent = [
        record
        for record in records
        if record["legacy_individual_failure_pattern"] == "FFF"
    ]
    all_floor = sum(
        all(record["legacy_below_amplitude_floor_pattern"])
        for record in records
    )
    nonzero_interpolation = sum(
        any(
            observation["round_trip_interpolation_infinity"] > 0.0
            for observation in record["observations"]
        )
        for record in records
    )
    by_field = {
        field_name: {
            measure: dict(
                sorted(
                    _class_count(
                        [item for item in records if item["field"] == field_name],
                        measure,
                    ).items()
                )
            )
            for measure in TDG2_MEASURE_NAMES
        }
        for field_name in TDG1_FIELD_NAMES
    }
    pattern_by_field_class = {
        measure: {
            " | ".join(pair): count
            for pair, count in sorted(
                Counter(
                    (
                        record["legacy_individual_failure_pattern"],
                        record[measure]["classification"],
                    )
                    for record in records
                ).items()
            )
        }
        for measure in TDG2_MEASURE_NAMES
    }
    enclosure_causes = {
        measure: dict(
            sorted(
                Counter(
                    record[
                        "field_enclosure_cause"
                        if measure == "field_tail_power"
                        else "derivative_enclosure_cause"
                    ]
                    for record in records
                    if record[measure]["classification"]
                    == "instrument_floor_or_interpolation_enclosure_dominated"
                ).items()
            )
        )
        for measure in TDG2_MEASURE_NAMES
    }
    finest_enclosure_causes = {
        measure: dict(
            sorted(
                Counter(
                    record[
                        "field_enclosure_cause"
                        if measure == "field_tail_power"
                        else "derivative_enclosure_cause"
                    ]
                    for record in finest_failures
                    if record[measure]["classification"]
                    == "instrument_floor_or_interpolation_enclosure_dominated"
                ).items()
            )
        )
        for measure in TDG2_MEASURE_NAMES
    }
    return {
        "method": method,
        "point_counts": list(TDG2_METHOD_POINT_COUNTS[method]),
        "signal_count": len(records),
        "classification_counts": {
            "field_tail_power": dict(sorted(field.items())),
            "derivative_tail_power": dict(sorted(derivative.items())),
        },
        "legacy_finest_individual_failure_count": len(finest_failures),
        "legacy_finest_above_floor_failure_count": len(above_floor),
        "legacy_persistent_all_resolution_failure_count": len(persistent),
        "all_three_below_DFT_floor_count": all_floor,
        "nonzero_interpolation_debit_ladder_count": nonzero_interpolation,
        "finest_failure_classification_counts": {
            measure: dict(sorted(_class_count(finest_failures, measure).items()))
            for measure in TDG2_MEASURE_NAMES
        },
        "persistent_failure_classification_counts": {
            measure: dict(sorted(_class_count(persistent, measure).items()))
            for measure in TDG2_MEASURE_NAMES
        },
        "above_floor_finest_failure_classification_counts": {
            measure: dict(sorted(_class_count(above_floor, measure).items()))
            for measure in TDG2_MEASURE_NAMES
        },
        "pattern_by_class": pattern_by_field_class,
        "classification_counts_by_field": by_field,
        "enclosure_cause_counts": enclosure_causes,
        "finest_failure_enclosure_cause_counts": finest_enclosure_causes,
    }


def _cross_method_findings(
    records_by_method: Mapping[str, Sequence[Mapping[str, Any]]],
) -> dict[str, Any]:
    indexed = {
        method: {
            (record["tracer"], record["field"]): record for record in records
        }
        for method, records in records_by_method.items()
    }
    if set(indexed["RK4"]) != set(indexed["SSPRK3"]):
        raise ValueError("PREF17 cross-method signal keys differ")
    output: dict[str, Any] = {"signal_count": len(indexed["RK4"])}
    for measure in TDG2_MEASURE_NAMES:
        pairs = Counter(
            (
                indexed["RK4"][key][measure]["classification"],
                indexed["SSPRK3"][key][measure]["classification"],
            )
            for key in indexed["RK4"]
        )
        both_finest_fail = Counter()
        either_finest_fail = Counter()
        for key in indexed["RK4"]:
            rk = indexed["RK4"][key]
            ss = indexed["SSPRK3"][key]
            pair = (rk[measure]["classification"], ss[measure]["classification"])
            rk_fail = rk["legacy_individual_failure_pattern"][-1] == "F"
            ss_fail = ss["legacy_individual_failure_pattern"][-1] == "F"
            if rk_fail and ss_fail:
                both_finest_fail[pair] += 1
            if rk_fail or ss_fail:
                either_finest_fail[pair] += 1
        prefix = "field" if measure == "field_tail_power" else "derivative"
        output[f"{prefix}_method_agreement_count"] = sum(
            count for (left, right), count in pairs.items() if left == right
        )
        output[f"{prefix}_classification_pairs"] = {
            " | ".join(pair): count for pair, count in sorted(pairs.items())
        }
        output[f"{prefix}_both_finest_fail_pairs"] = {
            " | ".join(pair): count
            for pair, count in sorted(both_finest_fail.items())
        }
        output[f"{prefix}_either_finest_fail_pairs"] = {
            " | ".join(pair): count
            for pair, count in sorted(either_finest_fail.items())
        }
    return output


def diagnose_actual_histories(
    histories: Mapping[str, Mapping[str, np.ndarray]],
) -> dict[str, Any]:
    proper, fields = _history_inputs(histories)
    records_by_method: dict[str, list[dict[str, Any]]] = {}
    methods: dict[str, Any] = {}
    for method in ("RK4", "SSPRK3"):
        records = [
            _signal_record(
                method=method,
                tracer=tracer,
                field_index=field_index,
                proper=proper[method],
                fields=fields[method],
            )
            for tracer in range(48)
            for field_index in range(len(TDG1_FIELD_NAMES))
        ]
        records_by_method[method] = records
        methods[method] = _method_findings(method, records)
    return {
        "methods": methods,
        "cross_method": _cross_method_findings(records_by_method),
        "signal_records": records_by_method,
        "frozen_class_names": list(TDG2_CLASS_NAMES),
        "frozen_measure_names": list(TDG2_MEASURE_NAMES),
        "diagnosis_is_not_a_replacement_admission": True,
    }


def _compact_findings(diagnosis: Mapping[str, Any]) -> dict[str, Any]:
    methods = diagnosis["methods"]
    rk = methods["RK4"]
    ss = methods["SSPRK3"]
    zero = "contracts_to_zero_at_required_power_order"
    nonzero = "nonzero_resolution_independent_within_enclosure"
    enclosure = "instrument_floor_or_interpolation_enclosure_dominated"

    def count(record: Mapping[str, Any], measure: str, name: str) -> int:
        return int(record["classification_counts"][measure].get(name, 0))

    def finest(record: Mapping[str, Any], measure: str, name: str) -> int:
        return int(record["finest_failure_classification_counts"][measure].get(name, 0))

    cross = diagnosis["cross_method"]
    compact = {
        "method_count": 2,
        "signal_ladders_per_method": 288,
        "total_signal_ladders": 576,
        "resolution_observations": 1728,
        "power_classifications": 1152,
        "legacy_finest_individual_failures": (
            rk["legacy_finest_individual_failure_count"]
            + ss["legacy_finest_individual_failure_count"]
        ),
        "legacy_finest_above_floor_failures": (
            rk["legacy_finest_above_floor_failure_count"]
            + ss["legacy_finest_above_floor_failure_count"]
        ),
        "legacy_persistent_all_resolution_failures": (
            rk["legacy_persistent_all_resolution_failure_count"]
            + ss["legacy_persistent_all_resolution_failure_count"]
        ),
        "all_three_below_DFT_floor_ladders": (
            rk["all_three_below_DFT_floor_count"]
            + ss["all_three_below_DFT_floor_count"]
        ),
        "nonzero_interpolation_debit_ladders": (
            rk["nonzero_interpolation_debit_ladder_count"]
            + ss["nonzero_interpolation_debit_ladder_count"]
        ),
        "zero_contracting_finest_failures_field": finest(
            rk, "field_tail_power", zero
        )
        + finest(ss, "field_tail_power", zero),
        "nonzero_resolution_independent_finest_failures_field": finest(
            rk, "field_tail_power", nonzero
        )
        + finest(ss, "field_tail_power", nonzero),
        "enclosure_dominated_finest_failures_field": finest(
            rk, "field_tail_power", enclosure
        )
        + finest(ss, "field_tail_power", enclosure),
        "zero_contracting_finest_failures_derivative": finest(
            rk, "derivative_tail_power", zero
        )
        + finest(ss, "derivative_tail_power", zero),
        "nonzero_resolution_independent_finest_failures_derivative": finest(
            rk, "derivative_tail_power", nonzero
        )
        + finest(ss, "derivative_tail_power", nonzero),
        "enclosure_dominated_finest_failures_derivative": finest(
            rk, "derivative_tail_power", enclosure
        )
        + finest(ss, "derivative_tail_power", enclosure),
        "nonzero_convergent_classifications": sum(
            count(record, measure, "converges_to_nonzero_resolved_power")
            for record in (rk, ss)
            for measure in TDG2_MEASURE_NAMES
        ),
        "unresolved_classifications": sum(
            count(record, measure, "unresolved_absolute_tail_behavior")
            for record in (rk, ss)
            for measure in TDG2_MEASURE_NAMES
        ),
        "replacement_temporal_admission_earned": False,
        "RK4": {
            "point_counts": rk["point_counts"],
            "field_zero_contracting": count(rk, "field_tail_power", zero),
            "field_nonzero_resolution_independent": count(
                rk, "field_tail_power", nonzero
            ),
            "field_enclosure_dominated": count(rk, "field_tail_power", enclosure),
            "derivative_zero_contracting": count(
                rk, "derivative_tail_power", zero
            ),
            "derivative_nonzero_resolution_independent": count(
                rk, "derivative_tail_power", nonzero
            ),
            "derivative_enclosure_dominated": count(
                rk, "derivative_tail_power", enclosure
            ),
            "finest_individual_failures": rk[
                "legacy_finest_individual_failure_count"
            ],
            "finest_above_floor_failures": rk[
                "legacy_finest_above_floor_failure_count"
            ],
            "persistent_all_resolution_failures": rk[
                "legacy_persistent_all_resolution_failure_count"
            ],
            "finest_failure_zero_contracting": finest(
                rk, "field_tail_power", zero
            ),
            "finest_failure_nonzero_resolution_independent": finest(
                rk, "field_tail_power", nonzero
            ),
            "finest_failure_enclosure_dominated": finest(
                rk, "field_tail_power", enclosure
            ),
            "finest_failure_all_three_DFT_floor": rk[
                "finest_failure_enclosure_cause_counts"
            ]["field_tail_power"].get("all_three_DFT_floor", 0),
            "finest_failure_interpolation_interval_reaches_zero": rk[
                "finest_failure_enclosure_cause_counts"
            ]["field_tail_power"].get(
                "finest_interpolation_interval_reaches_zero", 0
            ),
            "finest_failure_cross_resolution_intervals_overlap": rk[
                "finest_failure_enclosure_cause_counts"
            ]["field_tail_power"].get("cross_resolution_intervals_overlap", 0),
        },
        "SSPRK3": {
            "point_counts": ss["point_counts"],
            "field_zero_contracting": count(ss, "field_tail_power", zero),
            "field_nonzero_resolution_independent": count(
                ss, "field_tail_power", nonzero
            ),
            "field_enclosure_dominated": count(ss, "field_tail_power", enclosure),
            "derivative_zero_contracting": count(
                ss, "derivative_tail_power", zero
            ),
            "derivative_nonzero_resolution_independent": count(
                ss, "derivative_tail_power", nonzero
            ),
            "derivative_enclosure_dominated": count(
                ss, "derivative_tail_power", enclosure
            ),
            "finest_individual_failures": ss[
                "legacy_finest_individual_failure_count"
            ],
            "finest_above_floor_failures": ss[
                "legacy_finest_above_floor_failure_count"
            ],
            "persistent_all_resolution_failures": ss[
                "legacy_persistent_all_resolution_failure_count"
            ],
            "field_finest_failure_zero_contracting": finest(
                ss, "field_tail_power", zero
            ),
            "field_finest_failure_nonzero_resolution_independent": finest(
                ss, "field_tail_power", nonzero
            ),
            "field_finest_failure_enclosure_dominated": finest(
                ss, "field_tail_power", enclosure
            ),
            "derivative_finest_failure_zero_contracting": finest(
                ss, "derivative_tail_power", zero
            ),
            "derivative_finest_failure_nonzero_resolution_independent": finest(
                ss, "derivative_tail_power", nonzero
            ),
            "derivative_finest_failure_enclosure_dominated": finest(
                ss, "derivative_tail_power", enclosure
            ),
            "field_finest_failure_all_three_DFT_floor": ss[
                "finest_failure_enclosure_cause_counts"
            ]["field_tail_power"].get("all_three_DFT_floor", 0),
            "field_finest_failure_interpolation_interval_reaches_zero": ss[
                "finest_failure_enclosure_cause_counts"
            ]["field_tail_power"].get(
                "finest_interpolation_interval_reaches_zero", 0
            ),
            "derivative_finest_failure_all_three_DFT_floor": ss[
                "finest_failure_enclosure_cause_counts"
            ]["derivative_tail_power"].get("all_three_DFT_floor", 0),
            "derivative_finest_failure_interpolation_interval_reaches_zero": ss[
                "finest_failure_enclosure_cause_counts"
            ]["derivative_tail_power"].get(
                "finest_interpolation_interval_reaches_zero", 0
            ),
        },
        "cross_method": {
            "field_method_agreement_count": cross[
                "field_method_agreement_count"
            ],
            "derivative_method_agreement_count": cross[
                "derivative_method_agreement_count"
            ],
            "signal_count": cross["signal_count"],
            "field_both_methods_zero_contracting": cross[
                "field_classification_pairs"
            ].get(f"{zero} | {zero}", 0),
            "field_both_methods_nonzero_resolution_independent": cross[
                "field_classification_pairs"
            ].get(f"{nonzero} | {nonzero}", 0),
            "field_both_methods_enclosure_dominated": cross[
                "field_classification_pairs"
            ].get(f"{enclosure} | {enclosure}", 0),
            "derivative_both_methods_zero_contracting": cross[
                "derivative_classification_pairs"
            ].get(f"{zero} | {zero}", 0),
            "derivative_both_methods_nonzero_resolution_independent": cross[
                "derivative_classification_pairs"
            ].get(f"{nonzero} | {nonzero}", 0),
            "derivative_both_methods_enclosure_dominated": cross[
                "derivative_classification_pairs"
            ].get(f"{enclosure} | {enclosure}", 0),
        },
    }
    expected = _expected_observed()
    if compact != expected:
        mismatches = {
            key: {"observed": compact.get(key), "expected": expected.get(key)}
            for key in sorted(set(compact) | set(expected))
            if compact.get(key) != expected.get(key)
        }
        raise ValueError(
            "PREF17 compact diagnosis differs from frozen outcome: "
            + json.dumps(mismatches, sort_keys=True, allow_nan=False)
        )
    return compact


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
    config: Mapping[str, Any],
    histories: Mapping[str, Mapping[str, np.ndarray]],
    baseline: Mapping[str, Any],
) -> dict[str, bool]:
    proper, fields = _history_inputs(histories)
    target = [record.copy() for record in fields["RK4"]]
    index = np.unravel_index(int(np.argmax(np.abs(target[0]))), target[0].shape)
    target[0].view(np.uint64)[index] ^= np.uint64(1 << 44)
    tracer = int(index[1])
    field_index = int(index[2])
    mutated = _signal_record(
        method="RK4",
        tracer=tracer,
        field_index=field_index,
        proper=proper["RK4"],
        fields=target,
    )
    original = next(
        record
        for record in baseline["signal_records"]["RK4"]
        if record["tracer"] == tracer
        and record["field"] == TDG1_FIELD_NAMES[field_index]
    )
    controls = {
        "one_bit_history_mutation_detected": mutated != original,
        "method_ladder_mutation_rejected": False,
        "freeze_hash_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["immutable_lineage"].__setitem__(
                "freeze_result_sha256", "0" * 64
            ),
        ),
        "observed_count_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["observed_diagnosis"].__setitem__(
                "legacy_finest_individual_failures", 151
            ),
        ),
        "interpolation_theorem_promotion_rejected": _expect_config_rejection(
            config,
            lambda value: value["conclusions"].__setitem__(
                "interpolation_debit_is_the_dominant_unresolved_owner", False
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
    try:
        diagnose_absolute_tail_ladder(
            [record[:, tracer] for record in proper["RK4"]],
            [record[:, tracer, field_index] for record in fields["RK4"]],
            method="RK4",
            point_counts=(2049, 4097, 16385),
        )
    except (TypeError, ValueError):
        controls["method_ladder_mutation_rejected"] = True
    if set(controls.values()) != {True}:
        raise ValueError(f"PREF17 mutation control failed: {controls}")
    return controls


def record(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    config = load_config(config_path)
    freeze_lineage = _verify_freeze_lineage(config)
    raw_paths = _raw_paths(config)
    if not _raw_bundle_is_complete_or_absent(raw_paths):
        raise FileNotFoundError("PREF17 raw campaign bundle is absent")
    lineage = config["immutable_lineage"]
    labels = ("manifest", "event_log", "checkpoint", "campaign_result")
    expected_hashes = {
        "manifest": lineage["manifest_sha256"],
        "event_log": lineage["event_log_sha256"],
        "checkpoint": lineage["checkpoint_sha256"],
        "campaign_result": lineage["campaign_result_sha256"],
    }
    raw_hashes = {
        label: _sha(path) for label, path in zip(labels, raw_paths, strict=True)
    }
    if raw_hashes != expected_hashes:
        raise ValueError("PREF17 raw campaign hash differs")

    cal10_record = cal10.verify_canonical()
    _, event_path, checkpoint_path, result_path = raw_paths
    events, event_raw = cal10._load_events(event_path)
    campaign_result = cal10._load_json(result_path)
    metadata, states, grids, histories = cal10._restore_terminal_checkpoint(
        checkpoint_path,
        event_log_raw=event_raw,
        campaign_result=campaign_result,
    )
    if (
        len(events) != 41
        or metadata.get("terminal") is not True
        or metadata.get("completed_common_event_index") != 63
        or set(states) != set(key for keys in METHOD_KEYS.values() for key in keys)
        or set(grids) != set(states)
        or set(histories) != set(states)
        or cal10_record.get("gate_status", {}).get(
            "PROTO13_terminal_checkpoint_recomputed"
        )
        is not True
    ):
        raise ValueError("PREF17 CAL10 terminal restoration differs")
    diagnosis = diagnose_actual_histories(histories)
    compact = _compact_findings(diagnosis)
    mutations = _mutation_controls(config, histories, diagnosis)
    result = {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "metric_signature": "-+++",
        "riemann_convention": "plus_partial_mu_gamma_nu",
        "classification": "completed_TDG2_mixed_absolute_tail_diagnosis_TDG3_required",
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
            "immutable_freeze_lineage": freeze_lineage,
            "immutable_raw_bundle": {
                "raw_sha256": raw_hashes,
                "terminal_event_index": 63,
                "terminal_coordinate_time": 63 / 16,
                "history_arrays_loaded_read_only": True,
                "terminal_checkpoint_resumed": False,
                "state_advanced": False,
                "runtime_environment": campaign_result["runtime_environment"],
                "numerical_runtime_contract": campaign_result[
                    "numerical_runtime_contract"
                ],
            },
            "actual_history_diagnosis": diagnosis,
            "compact_findings": compact,
            "mutation_controls": mutations,
            "conclusions": dict(config["conclusions"]),
            "successor_boundary": dict(config["successor_boundary"]),
            "claim_boundary": {
                "historical_PROTO13_failure_reclassified": False,
                "absolute_tail_outcome_is_mixed": True,
                "individual_budget_obstruction_resolved": False,
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
    return cal10._serial(result)


def _validate_stored_record(
    value: Mapping[str, Any], config: Mapping[str, Any]
) -> None:
    payload = value.get("artifact_payload", {})
    if (
        value.get("schema_version") != 1
        or value.get("artifact_id") != ARTIFACT_ID
        or value.get("project_version") != PROJECT_VERSION
        or value.get("classification")
        != "completed_TDG2_mixed_absolute_tail_diagnosis_TDG3_required"
        or value.get("generated_by") != "scripts/reproduce_fgc_tdg2_pref17.py"
        or value.get("gate_status") != EXPECTED_CLAIMS
        or value.get("nonclaims")
        != {
            key: False
            for key, item in sorted(EXPECTED_CLAIMS.items())
            if item is False
        }
        or value.get("source_config_sha256")
        != {_rel(DEFAULT_CONFIG): _sha(DEFAULT_CONFIG)}
        or payload.get("compact_findings") != _expected_observed()
        or payload.get("conclusions") != config["conclusions"]
        or payload.get("successor_boundary") != config["successor_boundary"]
        or payload.get("claim_boundary")
        != {
            "historical_PROTO13_failure_reclassified": False,
            "absolute_tail_outcome_is_mixed": True,
            "individual_budget_obstruction_resolved": False,
            "replacement_temporal_admission_defined": False,
            "new_trajectory_authorized": False,
            "GR0_case_eligible": False,
            "candidate_branch_opened": False,
            "physical_question_answered": False,
        }
    ):
        raise ValueError("PREF17 stored result boundary differs")
    implementation = value.get("implementation_sha256", {})
    if (
        not isinstance(implementation, Mapping)
        or set(implementation) != {_rel(path) for path in IMPLEMENTATION}
        or any(
            _sha(REPOSITORY / relative) != expected
            for relative, expected in implementation.items()
        )
    ):
        raise ValueError("PREF17 implementation hash ledger differs")
    document = REPOSITORY / str(value.get("derivation_document", ""))
    if (
        value.get("derivation_document") != "docs/fgc-tdg2-pref17.md"
        or not document.is_file()
        or _sha(document) != value.get("derivation_document_sha256")
    ):
        raise ValueError("PREF17 derivation document differs")
    diagnosis = payload.get("actual_history_diagnosis", {})
    if (
        sum(len(records) for records in diagnosis.get("signal_records", {}).values())
        != 576
        or diagnosis.get("diagnosis_is_not_a_replacement_admission") is not True
    ):
        raise ValueError("PREF17 complete signal evidence differs")
    mutations = payload.get("mutation_controls", {})
    if not isinstance(mutations, Mapping) or set(mutations.values()) != {True}:
        raise ValueError("PREF17 mutation controls differ")


def verify_canonical(
    config_path: Path = DEFAULT_CONFIG,
    output_path: Path = DEFAULT_OUTPUT,
) -> dict[str, Any]:
    config = load_config(config_path)
    raw = output_path.read_bytes()
    stored = json.loads(raw)
    if not isinstance(stored, dict) or raw != _canonical_bytes(stored):
        raise ValueError("PREF17 stored result is not canonical JSON")
    _validate_stored_record(stored, config)
    paths = _raw_paths(config)
    if _raw_bundle_is_complete_or_absent(paths):
        observed = record(config_path)
        if stored != observed:
            raise ValueError("PREF17 stored result differs from raw reproduction")
    return stored


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.write and args.check:
        raise SystemExit("choose either --write or --check")
    if args.write:
        value = record(args.config)
        _atomic_write(args.output, _canonical_bytes(value))
        print(f"wrote {_rel(args.output)}")
    elif args.check:
        value = verify_canonical(args.config, args.output)
        print(
            f"{ARTIFACT_ID}: histories={value['gate_status']['TDG2_actual_histories_diagnosed']} "
            f"replacement={value['gate_status']['replacement_temporal_admission_defined']}"
        )
    else:
        print(_canonical(record(args.config)), end="")


if __name__ == "__main__":
    main()
