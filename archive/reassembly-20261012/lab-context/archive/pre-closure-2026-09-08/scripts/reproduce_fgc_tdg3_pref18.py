#!/usr/bin/env python3
"""Bind TDG3's read-only native-grid diagnosis of the PROTO13 histories."""

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
from scripts import reproduce_fgc_tdg2_pref17 as pref17  # noqa: E402
from scripts import reproduce_fgc_tdg3_frz1 as freeze  # noqa: E402
from recursive_horizons.fgc.evolution.tdg3_native_tail_diagnosis import (  # noqa: E402
    TDG3_COMBINED_CLASS_NAMES,
    TDG3_ESTIMATOR_NAMES,
    TDG3_FIELD_NAMES,
    TDG3_MEASURE_NAMES,
    TDG3_METHOD_POINT_COUNTS,
    TDG3_POWER_CLASS_NAMES,
    diagnose_native_tail_ladder,
)


ARTIFACT_ID = "FGC-1-TDG3-PREF18"
PROJECT_VERSION = "0.11.0"
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-tdg3-pref18.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-tdg3-pref18.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-tdg3-pref18.md"
FREEZE_RESULT = REPOSITORY / "results/fgc-1-tdg3-frz1.json"
FREEZE_CONFIG = REPOSITORY / "configs/fgc/fgc-1-tdg3-frz1.toml"
FREEZE_REPRODUCER = REPOSITORY / "scripts/reproduce_fgc_tdg3_frz1.py"
DISCRIMINATOR_MODULE = (
    REPOSITORY
    / "src/recursive_horizons/fgc/evolution/tdg3_native_tail_diagnosis.py"
)
FREEZE_DOCUMENT = REPOSITORY / "docs/fgc-tdg3-frz1.md"
SOURCE_RESULT = REPOSITORY / "results/fgc-1-tdg2-pref17.json"
SOURCE_CONFIG = REPOSITORY / "configs/fgc/fgc-1-tdg2-pref17.toml"
SOURCE_REPRODUCER = REPOSITORY / "scripts/reproduce_fgc_tdg2_pref17.py"
SOURCE_DOCUMENT = REPOSITORY / "docs/fgc-tdg2-pref17.md"
CAL10_REPRODUCER = REPOSITORY / "scripts/reproduce_fgc_cal10_pref15.py"
IMPLEMENTATION = (Path(__file__).resolve(),)
METHOD_KEYS = pref17.METHOD_KEYS

TDG2_MEASURE_BY_TDG3 = {
    "phase_field_tail_power": "field_tail_power",
    "phase_derivative_tail_power": "derivative_tail_power",
}
TDG2_CAUSE_KEY_BY_TDG3 = {
    "phase_field_tail_power": "field_enclosure_cause",
    "phase_derivative_tail_power": "derivative_enclosure_cause",
}
TDG2_ZERO = "contracts_to_zero_at_required_power_order"
TDG2_NONZERO = "nonzero_resolution_independent_within_enclosure"
TDG2_ENCLOSURE = "instrument_floor_or_interpolation_enclosure_dominated"
TDG3_ZERO = "native_estimators_agree_zero_contracting"
TDG3_NONZERO = "native_estimators_agree_nonzero"
TDG3_FLOOR = "native_estimators_agree_floor_dominated"
TDG3_UNRESOLVED = "native_estimators_agree_unresolved"
TDG3_DISAGREE = "native_estimators_disagree"


EXPECTED_CLAIMS = {
    "TDG2_actual_histories_diagnosed": True,
    "TDG2_absolute_tail_outcome_is_mixed": True,
    "TDG3_interpolation_ownership_discriminator_frozen": True,
    "TDG3_actual_histories_diagnosed": True,
    "TDG3_native_grid_outcome_is_mixed": True,
    "common_grid_interpolation_removed_from_TDG3_measurement": True,
    "common_grid_interpolation_is_the_sole_temporal_failure_owner": False,
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


def _as_float(value: str) -> float:
    numerator, separator, denominator = value.partition("/")
    return float(int(numerator) / int(denominator)) if separator else float(value)


def _expected_lineage() -> dict[str, Any]:
    return {
        "freeze_commit": "5c8d3cc53155ca42db38ddc2005ed60835743309",
        "freeze_result_sha256": "8076210bcaa33f58c5c053c8731df10809427564d3a4594c9ea734acd3cdc770",
        "freeze_config_sha256": "de72fdab1d53b713e321ba17f5177dcddb9de825fef61cc7e0fa31974e455d78",
        "freeze_reproducer_sha256": "18ad0f47d2a8de0f15f4e967a81a99537c09c83f3378487fc6a74e8fc50633c5",
        "discriminator_module_sha256": "9b53ea697cbc8f9026db197849d0bf0b355b92dda0dc22e6bb58d4a4369fb984",
        "freeze_document_sha256": "a8d63e5d735fd729b56404b32eebc2c79f5dba3dfa6097acc3b1cb11db03f4dd",
        "source_result_sha256": "1a0d9ee99b2eb2f15f10b0625c1cbfe63db15d26692ff5f7b1975d370db593a3",
        "source_config_sha256": "3f2ec6acf07779c4ae0d5b47fe235f66724d4f3b86975eb391dfae02d325c000",
        "source_reproducer_sha256": "f0a939c5f46a0e4fb57b9118e6107788c5f585555802a07a5e8324a63111e776",
        "source_document_sha256": "8495b96c84bb5cba4af869701c732fd40793b0cc9b951990f1a71a64b129cd01",
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
        "combined_power_classifications": 1152,
        "estimator_power_classifications": 2304,
        "native_estimator_agreement_count": 1119,
        "native_estimator_disagreement_count": 33,
        "native_estimators_agree_zero_contracting": 156,
        "native_estimators_agree_nonzero": 167,
        "native_estimators_agree_floor_dominated": 242,
        "native_estimators_agree_unresolved": 554,
        "TDG2_finest_zero_classifications": 92,
        "TDG2_finest_zero_survive_native_zero": 92,
        "TDG2_finest_nonzero_classifications": 9,
        "TDG2_finest_nonzero_survive_native_nonzero": 0,
        "TDG2_finest_nonzero_become_native_zero": 6,
        "TDG2_finest_nonzero_become_native_unresolved": 3,
        "TDG2_finest_enclosure_classifications": 203,
        "TDG2_finest_enclosure_become_native_zero": 56,
        "TDG2_finest_enclosure_become_native_nonzero": 15,
        "TDG2_finest_enclosure_become_native_floor": 28,
        "TDG2_finest_enclosure_remain_native_unresolved": 101,
        "TDG2_finest_enclosure_native_disagreement": 3,
        "TDG2_finest_interpolation_owned_classifications": 169,
        "TDG2_finest_interpolation_become_native_zero": 50,
        "TDG2_finest_interpolation_become_native_nonzero": 15,
        "TDG2_finest_interpolation_remain_native_unresolved": 101,
        "TDG2_finest_interpolation_native_disagreement": 3,
        "replacement_temporal_admission_earned": False,
        "RK4": {
            "point_counts": [2049, 4097, 8193],
            "field_zero_contracting": 59,
            "field_nonzero": 34,
            "field_floor_dominated": 40,
            "field_unresolved": 147,
            "field_estimator_disagreement": 8,
            "derivative_zero_contracting": 59,
            "derivative_nonzero": 35,
            "derivative_floor_dominated": 40,
            "derivative_unresolved": 145,
            "derivative_estimator_disagreement": 9,
        },
        "SSPRK3": {
            "point_counts": [4097, 8193, 16385],
            "field_zero_contracting": 19,
            "field_nonzero": 47,
            "field_floor_dominated": 81,
            "field_unresolved": 132,
            "field_estimator_disagreement": 9,
            "derivative_zero_contracting": 19,
            "derivative_nonzero": 51,
            "derivative_floor_dominated": 81,
            "derivative_unresolved": 130,
            "derivative_estimator_disagreement": 7,
        },
        "cross_method": {
            "signal_count": 288,
            "field_same_combined_class": 170,
            "derivative_same_combined_class": 167,
            "field_both_zero_contracting": 19,
            "field_both_nonzero": 8,
            "field_both_floor_dominated": 40,
            "field_both_unresolved": 102,
            "field_both_estimator_disagreement": 1,
            "derivative_both_zero_contracting": 19,
            "derivative_both_nonzero": 10,
            "derivative_both_floor_dominated": 40,
            "derivative_both_unresolved": 98,
            "derivative_both_estimator_disagreement": 0,
        },
    }


def validate_config_data(config: Mapping[str, Any]) -> None:
    _strict_keys(
        "PREF18 config",
        config,
        {
            "schema_version",
            "artifact_id",
            "project_version",
            "metric_signature",
            "riemann_convention",
            "freeze_result",
            "freeze_config",
            "source_result",
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
        "source_result": _rel(SOURCE_RESULT),
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
        raise ValueError("PREF18 top-level contract differs")
    if config.get("scope") != {
        "target": "FGC-2-SF1-TDG3",
        "role": "post_freeze_read_only_native_grid_history_diagnosis",
        "calibration_branch": "GR-0",
        "source_terminal_event": 63,
        "source_terminal_coordinate_time": "63/16",
        "terminal_histories_read": True,
        "terminal_checkpoint_resumed": False,
        "state_advanced": False,
        "common_grid_resampling_performed": False,
        "SGBL_trajectory_read": False,
        "FGCQR_trajectory_read": False,
        "replacement_temporal_admission_defined": False,
        "physical_or_candidate_question_answered": False,
    }:
        raise ValueError("PREF18 scope differs")
    if _as_float(config["scope"]["source_terminal_coordinate_time"]) != 63 / 16:
        raise ValueError("PREF18 source time differs")
    if config.get("immutable_lineage") != _expected_lineage():
        raise ValueError("PREF18 immutable lineage differs")
    if config.get("observed_diagnosis") != _expected_observed():
        raise ValueError("PREF18 observed diagnosis differs")
    if config.get("conclusions") != {
        "historical_PROTO13_failure_preserved": True,
        "TDG3_native_grid_history_diagnosis_completed": True,
        "common_grid_interpolation_was_a_real_but_not_sole_owner": True,
        "all_TDG2_zero_contracting_finest_failures_survive_both_native_estimators": True,
        "TDG2_tentative_nonzero_finest_subset_is_not_confirmed": True,
        "majority_of_TDG2_enclosure_dominated_finest_classifications_remain_unresolved_or_disagree": True,
        "native_estimator_agreement_does_not_supply_a_continuum_enclosure": True,
        "replacement_temporal_admission_not_yet_justified": True,
        "TDG4_sampling_identifiability_theorem_design_may_begin": True,
        "PROTO14_design_may_not_yet_begin": True,
    }:
        raise ValueError("PREF18 conclusions differ")
    if set(config.get("proof_contract", {}).values()) != {True}:
        raise ValueError("PREF18 proof contract must be all-of")
    if config.get("successor_boundary") != {
        "TDG3_interpolation_ownership_discriminator_frozen": True,
        "TDG3_actual_histories_diagnosed": True,
        "TDG4_sampling_identifiability_theorem_design_may_begin": True,
        "replacement_temporal_gate_design_authorized": False,
        "replacement_temporal_admission_defined": False,
        "PROTO14_frozen": False,
        "fresh_GR0_dynamic_calibration_completed": False,
        "GR0_case_eligible": False,
        "SGBL_execution_authorized": False,
        "FGCQR_holdout_execution_authorized": False,
        "DEF1_execution_authorized": False,
    }:
        raise ValueError("PREF18 successor boundary differs")
    if config.get("claims") != EXPECTED_CLAIMS:
        raise ValueError("PREF18 claims differ")


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
        raise ValueError("PREF18 raw campaign bundle is partial")
    return all(present)


def _verify_freeze_lineage(
    config: Mapping[str, Any],
) -> tuple[dict[str, Any], dict[str, Any]]:
    lineage = config["immutable_lineage"]
    commit = lineage["freeze_commit"]
    _git("cat-file", "-e", f"{commit}^{{commit}}")
    if _git("merge-base", "--is-ancestor", commit, "HEAD", check=False).returncode:
        raise ValueError("PREF18 freeze commit is not an ancestor of HEAD")
    tracked = {
        _rel(FREEZE_RESULT): lineage["freeze_result_sha256"],
        _rel(FREEZE_CONFIG): lineage["freeze_config_sha256"],
        _rel(FREEZE_REPRODUCER): lineage["freeze_reproducer_sha256"],
        _rel(DISCRIMINATOR_MODULE): lineage["discriminator_module_sha256"],
        _rel(FREEZE_DOCUMENT): lineage["freeze_document_sha256"],
        _rel(SOURCE_RESULT): lineage["source_result_sha256"],
        _rel(SOURCE_CONFIG): lineage["source_config_sha256"],
        _rel(SOURCE_REPRODUCER): lineage["source_reproducer_sha256"],
        _rel(SOURCE_DOCUMENT): lineage["source_document_sha256"],
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
            raise ValueError(f"PREF18 frozen tracked blob differs: {relative}")
    freeze_record = freeze.verify_canonical()
    source_record = pref17.verify_canonical()
    if (
        freeze_record.get("artifact_id") != "FGC-1-TDG3-FRZ1"
        or freeze_record.get("gate_status", {}).get(
            "TDG3_interpolation_ownership_discriminator_frozen"
        )
        is not True
        or freeze_record.get("gate_status", {}).get("TDG3_actual_histories_diagnosed")
        is not False
        or source_record.get("artifact_id") != "FGC-1-TDG2-PREF17"
        or source_record.get("gate_status", {}).get("TDG2_actual_histories_diagnosed")
        is not True
    ):
        raise ValueError("PREF18 source boundary differs")
    return (
        {
            "freeze_commit": commit,
            "freeze_commit_is_ancestor_of_HEAD": True,
            "tracked_blob_sha256": tracked,
            "freeze_authorized_actual_history_execution": True,
            "source_TDG2_actual_histories_diagnosed": True,
            "freeze_loaded_actual_history": False,
        },
        source_record,
    )


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


def _classification_evidence(record: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "classification": record["classification"],
        "observed_powers": list(record["observed_powers"]),
        "lower_bounds": list(record["lower_bounds"]),
        "upper_bounds": list(record["upper_bounds"]),
        "zero_power_order_lower_bounds": list(
            record["zero_power_order_lower_bounds"]
        ),
        "difference_order_lower_bound": record["difference_order_lower_bound"],
        "conservative_native_surrogate_lower_bound": record[
            "conservative_native_surrogate_lower_bound"
        ],
        "every_resolution_below_declared_amplitude_floor": record[
            "every_resolution_below_declared_amplitude_floor"
        ],
        "replacement_admission_authorized": record[
            "replacement_admission_authorized"
        ],
    }


def _observation_evidence(
    observations: Sequence[Mapping[str, Any]], point_counts: Sequence[int]
) -> list[dict[str, Any]]:
    output = []
    for point_count, record in zip(point_counts, observations, strict=True):
        output.append(
            {
                "point_count": point_count,
                "proper_time_start": record["proper_time_start"],
                "proper_time_stop": record["proper_time_stop"],
                "proper_time_span": record["proper_time_span"],
                "minimum_normalized_spacing": record[
                    "minimum_normalized_spacing"
                ],
                "maximum_normalized_spacing": record[
                    "maximum_normalized_spacing"
                ],
                "normalized_spacing_ratio": record["normalized_spacing_ratio"],
                "peak_coefficient_amplitude": record[
                    "peak_coefficient_amplitude"
                ],
                "binary64_peak_coefficient_amplitude": record[
                    "binary64_peak_coefficient_amplitude"
                ],
                "coefficient_arithmetic_debit": record[
                    "coefficient_arithmetic_debit"
                ],
                "declared_normalized_coefficient_floor": record[
                    "declared_normalized_coefficient_floor"
                ],
                "top_band_is_below_declared_amplitude_floor": record[
                    "top_band_is_below_declared_amplitude_floor"
                ],
                "phase_field_tail_power": dict(record["phase_field_tail_power"]),
                "phase_derivative_tail_power": dict(
                    record["phase_derivative_tail_power"]
                ),
                "common_grid_resampling_performed": record[
                    "common_grid_resampling_performed"
                ],
                "continuum_surrogate_enclosure_claimed": record[
                    "continuum_surrogate_enclosure_claimed"
                ],
            }
        )
    return output


def _source_evidence(record: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "legacy_individual_failure_pattern": record[
            "legacy_individual_failure_pattern"
        ],
        "legacy_below_amplitude_floor_pattern": list(
            record["legacy_below_amplitude_floor_pattern"]
        ),
        "field_tail_classification": record["field_tail_power"]["classification"],
        "derivative_tail_classification": record["derivative_tail_power"][
            "classification"
        ],
        "field_enclosure_cause": record["field_enclosure_cause"],
        "derivative_enclosure_cause": record["derivative_enclosure_cause"],
    }


def _signal_record(
    *,
    method: str,
    tracer: int,
    field_index: int,
    proper: Sequence[np.ndarray],
    fields: Sequence[np.ndarray],
    source_record: Mapping[str, Any],
) -> dict[str, Any]:
    point_counts = TDG3_METHOD_POINT_COUNTS[method]
    diagnosis = diagnose_native_tail_ladder(
        [record[:, tracer] for record in proper],
        [record[:, tracer, field_index] for record in fields],
        method=method,
        point_counts=point_counts,
    )
    estimators: dict[str, Any] = {}
    for estimator in TDG3_ESTIMATOR_NAMES:
        estimators[estimator] = {
            "classifications": {
                measure: _classification_evidence(
                    diagnosis["classifications"][estimator][measure]
                )
                for measure in TDG3_MEASURE_NAMES
            },
            "observations": _observation_evidence(
                diagnosis["observations"][estimator], point_counts
            ),
        }
    spreads = {}
    for measure in TDG3_MEASURE_NAMES:
        left = estimators[TDG3_ESTIMATOR_NAMES[0]]["classifications"][measure][
            "observed_powers"
        ]
        right = estimators[TDG3_ESTIMATOR_NAMES[1]]["classifications"][measure][
            "observed_powers"
        ]
        spreads[measure] = [abs(a - b) for a, b in zip(left, right, strict=True)]
    return {
        "method": method,
        "tracer": tracer,
        "field": TDG3_FIELD_NAMES[field_index],
        "TDG2_source": _source_evidence(source_record),
        "estimators": estimators,
        "combined_classifications": dict(diagnosis["combined_classifications"]),
        "inter_estimator_absolute_power_spread": spreads,
        "common_grid_resampling_performed": False,
        "continuum_surrogate_enclosure_claimed": False,
        "classification_is_not_a_replacement_admission": True,
    }


def _class_count(
    records: Sequence[Mapping[str, Any]], measure: str
) -> Counter[str]:
    return Counter(record["combined_classifications"][measure] for record in records)


def _estimator_class_count(
    records: Sequence[Mapping[str, Any]], estimator: str, measure: str
) -> Counter[str]:
    return Counter(
        record["estimators"][estimator]["classifications"][measure][
            "classification"
        ]
        for record in records
    )


def _tdg2_class(record: Mapping[str, Any], measure: str) -> str:
    key = TDG2_MEASURE_BY_TDG3[measure]
    return record["TDG2_source"][f"{key.removesuffix('_power')}_classification"]


def _method_findings(
    method: str, records: Sequence[Mapping[str, Any]]
) -> dict[str, Any]:
    if len(records) != 288:
        raise ValueError(f"PREF18 {method} signal count differs")
    combined_counts = {
        measure: dict(sorted(_class_count(records, measure).items()))
        for measure in TDG3_MEASURE_NAMES
    }
    estimator_counts = {
        estimator: {
            measure: dict(
                sorted(_estimator_class_count(records, estimator, measure).items())
            )
            for measure in TDG3_MEASURE_NAMES
        }
        for estimator in TDG3_ESTIMATOR_NAMES
    }
    estimator_pairs = {}
    tdg2_cross_tabs = {}
    finest_cross_tabs = {}
    finest_cause_cross_tabs = {}
    for measure in TDG3_MEASURE_NAMES:
        estimator_pairs[measure] = {
            " | ".join(pair): count
            for pair, count in sorted(
                Counter(
                    (
                        record["estimators"][TDG3_ESTIMATOR_NAMES[0]][
                            "classifications"
                        ][measure]["classification"],
                        record["estimators"][TDG3_ESTIMATOR_NAMES[1]][
                            "classifications"
                        ][measure]["classification"],
                    )
                    for record in records
                ).items()
            )
        }
        tdg2_cross_tabs[measure] = {
            " | ".join(pair): count
            for pair, count in sorted(
                Counter(
                    (
                        _tdg2_class(record, measure),
                        record["combined_classifications"][measure],
                    )
                    for record in records
                ).items()
            )
        }
        finest = [
            record
            for record in records
            if record["TDG2_source"]["legacy_individual_failure_pattern"][-1]
            == "F"
        ]
        finest_cross_tabs[measure] = {
            " | ".join(pair): count
            for pair, count in sorted(
                Counter(
                    (
                        _tdg2_class(record, measure),
                        record["combined_classifications"][measure],
                    )
                    for record in finest
                ).items()
            )
        }
        cause_key = TDG2_CAUSE_KEY_BY_TDG3[measure]
        finest_cause_cross_tabs[measure] = {
            " | ".join(pair): count
            for pair, count in sorted(
                Counter(
                    (
                        record["TDG2_source"][cause_key],
                        record["combined_classifications"][measure],
                    )
                    for record in finest
                    if record["TDG2_source"][cause_key] is not None
                ).items()
            )
        }
    by_field = {
        field_name: {
            measure: dict(
                sorted(
                    _class_count(
                        [record for record in records if record["field"] == field_name],
                        measure,
                    ).items()
                )
            )
            for measure in TDG3_MEASURE_NAMES
        }
        for field_name in TDG3_FIELD_NAMES
    }
    return {
        "method": method,
        "point_counts": list(TDG3_METHOD_POINT_COUNTS[method]),
        "signal_count": len(records),
        "combined_classification_counts": combined_counts,
        "estimator_classification_counts": estimator_counts,
        "estimator_classification_pairs": estimator_pairs,
        "TDG2_to_TDG3_cross_tabs": tdg2_cross_tabs,
        "TDG2_finest_failure_cross_tabs": finest_cross_tabs,
        "TDG2_finest_enclosure_cause_cross_tabs": finest_cause_cross_tabs,
        "combined_classification_counts_by_field": by_field,
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
        raise ValueError("PREF18 cross-method signal keys differ")
    output: dict[str, Any] = {"signal_count": len(indexed["RK4"])}
    for measure in TDG3_MEASURE_NAMES:
        pairs = Counter(
            (
                indexed["RK4"][key]["combined_classifications"][measure],
                indexed["SSPRK3"][key]["combined_classifications"][measure],
            )
            for key in indexed["RK4"]
        )
        prefix = "field" if measure == "phase_field_tail_power" else "derivative"
        output[f"{prefix}_same_combined_class_count"] = sum(
            count for (left, right), count in pairs.items() if left == right
        )
        output[f"{prefix}_combined_classification_pairs"] = {
            " | ".join(pair): count for pair, count in sorted(pairs.items())
        }
    return output


def diagnose_actual_histories(
    histories: Mapping[str, Mapping[str, np.ndarray]],
    source_diagnosis: Mapping[str, Any],
) -> dict[str, Any]:
    proper, fields = _history_inputs(histories)
    source_index = {
        method: {
            (record["tracer"], record["field"]): record for record in records
        }
        for method, records in source_diagnosis["signal_records"].items()
    }
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
                source_record=source_index[method][
                    (tracer, TDG3_FIELD_NAMES[field_index])
                ],
            )
            for tracer in range(48)
            for field_index in range(len(TDG3_FIELD_NAMES))
        ]
        records_by_method[method] = records
        methods[method] = _method_findings(method, records)
    return {
        "methods": methods,
        "cross_method": _cross_method_findings(records_by_method),
        "signal_records": records_by_method,
        "frozen_estimator_names": list(TDG3_ESTIMATOR_NAMES),
        "frozen_power_class_names": list(TDG3_POWER_CLASS_NAMES),
        "frozen_combined_class_names": list(TDG3_COMBINED_CLASS_NAMES),
        "frozen_measure_names": list(TDG3_MEASURE_NAMES),
        "common_grid_resampling_performed": False,
        "continuum_surrogate_enclosure_claimed": False,
        "diagnosis_is_not_a_replacement_admission": True,
    }


def _combined_count(
    diagnosis: Mapping[str, Any], class_name: str
) -> int:
    return sum(
        int(method["combined_classification_counts"][measure].get(class_name, 0))
        for method in diagnosis["methods"].values()
        for measure in TDG3_MEASURE_NAMES
    )


def _finest_records(
    diagnosis: Mapping[str, Any], source_class: str
) -> list[tuple[Mapping[str, Any], str]]:
    return [
        (record, measure)
        for records in diagnosis["signal_records"].values()
        for record in records
        for measure in TDG3_MEASURE_NAMES
        if record["TDG2_source"]["legacy_individual_failure_pattern"][-1] == "F"
        and _tdg2_class(record, measure) == source_class
    ]


def _count_combined(
    records: Sequence[tuple[Mapping[str, Any], str]], class_name: str
) -> int:
    return sum(
        record["combined_classifications"][measure] == class_name
        for record, measure in records
    )


def _compact_method(
    diagnosis: Mapping[str, Any], method: str
) -> dict[str, Any]:
    record = diagnosis["methods"][method]
    field = record["combined_classification_counts"]["phase_field_tail_power"]
    derivative = record["combined_classification_counts"][
        "phase_derivative_tail_power"
    ]
    return {
        "point_counts": record["point_counts"],
        "field_zero_contracting": int(field.get(TDG3_ZERO, 0)),
        "field_nonzero": int(field.get(TDG3_NONZERO, 0)),
        "field_floor_dominated": int(field.get(TDG3_FLOOR, 0)),
        "field_unresolved": int(field.get(TDG3_UNRESOLVED, 0)),
        "field_estimator_disagreement": int(field.get(TDG3_DISAGREE, 0)),
        "derivative_zero_contracting": int(derivative.get(TDG3_ZERO, 0)),
        "derivative_nonzero": int(derivative.get(TDG3_NONZERO, 0)),
        "derivative_floor_dominated": int(derivative.get(TDG3_FLOOR, 0)),
        "derivative_unresolved": int(derivative.get(TDG3_UNRESOLVED, 0)),
        "derivative_estimator_disagreement": int(
            derivative.get(TDG3_DISAGREE, 0)
        ),
    }


def _compact_cross_method(diagnosis: Mapping[str, Any]) -> dict[str, Any]:
    cross = diagnosis["cross_method"]

    def pair_count(measure_prefix: str, class_name: str) -> int:
        return int(
            cross[f"{measure_prefix}_combined_classification_pairs"].get(
                f"{class_name} | {class_name}", 0
            )
        )

    return {
        "signal_count": cross["signal_count"],
        "field_same_combined_class": cross["field_same_combined_class_count"],
        "derivative_same_combined_class": cross[
            "derivative_same_combined_class_count"
        ],
        "field_both_zero_contracting": pair_count("field", TDG3_ZERO),
        "field_both_nonzero": pair_count("field", TDG3_NONZERO),
        "field_both_floor_dominated": pair_count("field", TDG3_FLOOR),
        "field_both_unresolved": pair_count("field", TDG3_UNRESOLVED),
        "field_both_estimator_disagreement": pair_count("field", TDG3_DISAGREE),
        "derivative_both_zero_contracting": pair_count("derivative", TDG3_ZERO),
        "derivative_both_nonzero": pair_count("derivative", TDG3_NONZERO),
        "derivative_both_floor_dominated": pair_count("derivative", TDG3_FLOOR),
        "derivative_both_unresolved": pair_count("derivative", TDG3_UNRESOLVED),
        "derivative_both_estimator_disagreement": pair_count(
            "derivative", TDG3_DISAGREE
        ),
    }


def _compact_findings(diagnosis: Mapping[str, Any]) -> dict[str, Any]:
    zero = _finest_records(diagnosis, TDG2_ZERO)
    nonzero = _finest_records(diagnosis, TDG2_NONZERO)
    enclosure = _finest_records(diagnosis, TDG2_ENCLOSURE)
    interpolation = [
        (record, measure)
        for record, measure in enclosure
        if record["TDG2_source"][TDG2_CAUSE_KEY_BY_TDG3[measure]]
        == "finest_interpolation_interval_reaches_zero"
    ]
    compact = {
        "method_count": 2,
        "signal_ladders_per_method": 288,
        "total_signal_ladders": sum(
            len(records) for records in diagnosis["signal_records"].values()
        ),
        "combined_power_classifications": sum(
            len(records) * len(TDG3_MEASURE_NAMES)
            for records in diagnosis["signal_records"].values()
        ),
        "estimator_power_classifications": sum(
            len(records) * len(TDG3_MEASURE_NAMES) * len(TDG3_ESTIMATOR_NAMES)
            for records in diagnosis["signal_records"].values()
        ),
        "native_estimator_agreement_count": 1152
        - _combined_count(diagnosis, TDG3_DISAGREE),
        "native_estimator_disagreement_count": _combined_count(
            diagnosis, TDG3_DISAGREE
        ),
        "native_estimators_agree_zero_contracting": _combined_count(
            diagnosis, TDG3_ZERO
        ),
        "native_estimators_agree_nonzero": _combined_count(
            diagnosis, TDG3_NONZERO
        ),
        "native_estimators_agree_floor_dominated": _combined_count(
            diagnosis, TDG3_FLOOR
        ),
        "native_estimators_agree_unresolved": _combined_count(
            diagnosis, TDG3_UNRESOLVED
        ),
        "TDG2_finest_zero_classifications": len(zero),
        "TDG2_finest_zero_survive_native_zero": _count_combined(zero, TDG3_ZERO),
        "TDG2_finest_nonzero_classifications": len(nonzero),
        "TDG2_finest_nonzero_survive_native_nonzero": _count_combined(
            nonzero, TDG3_NONZERO
        ),
        "TDG2_finest_nonzero_become_native_zero": _count_combined(
            nonzero, TDG3_ZERO
        ),
        "TDG2_finest_nonzero_become_native_unresolved": _count_combined(
            nonzero, TDG3_UNRESOLVED
        ),
        "TDG2_finest_enclosure_classifications": len(enclosure),
        "TDG2_finest_enclosure_become_native_zero": _count_combined(
            enclosure, TDG3_ZERO
        ),
        "TDG2_finest_enclosure_become_native_nonzero": _count_combined(
            enclosure, TDG3_NONZERO
        ),
        "TDG2_finest_enclosure_become_native_floor": _count_combined(
            enclosure, TDG3_FLOOR
        ),
        "TDG2_finest_enclosure_remain_native_unresolved": _count_combined(
            enclosure, TDG3_UNRESOLVED
        ),
        "TDG2_finest_enclosure_native_disagreement": _count_combined(
            enclosure, TDG3_DISAGREE
        ),
        "TDG2_finest_interpolation_owned_classifications": len(interpolation),
        "TDG2_finest_interpolation_become_native_zero": _count_combined(
            interpolation, TDG3_ZERO
        ),
        "TDG2_finest_interpolation_become_native_nonzero": _count_combined(
            interpolation, TDG3_NONZERO
        ),
        "TDG2_finest_interpolation_remain_native_unresolved": _count_combined(
            interpolation, TDG3_UNRESOLVED
        ),
        "TDG2_finest_interpolation_native_disagreement": _count_combined(
            interpolation, TDG3_DISAGREE
        ),
        "replacement_temporal_admission_earned": False,
        "RK4": _compact_method(diagnosis, "RK4"),
        "SSPRK3": _compact_method(diagnosis, "SSPRK3"),
        "cross_method": _compact_cross_method(diagnosis),
    }
    expected = _expected_observed()
    if compact != expected:
        mismatches = {
            key: {"observed": compact.get(key), "expected": expected.get(key)}
            for key in sorted(set(compact) | set(expected))
            if compact.get(key) != expected.get(key)
        }
        raise ValueError(
            "PREF18 compact diagnosis differs from observed outcome: "
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
    source_diagnosis: Mapping[str, Any],
    baseline: Mapping[str, Any],
) -> dict[str, bool]:
    proper, fields = _history_inputs(histories)
    source_index = {
        method: {
            (record["tracer"], record["field"]): record for record in records
        }
        for method, records in source_diagnosis["signal_records"].items()
    }
    target = [record.copy() for record in fields["RK4"]]
    interior = target[0][8:56]
    inner_index = np.unravel_index(
        int(np.argmax(np.abs(interior))), interior.shape
    )
    index = (int(inner_index[0]) + 8, int(inner_index[1]), int(inner_index[2]))
    target[0].view(np.uint64)[index] ^= np.uint64(1 << 44)
    tracer = int(index[1])
    field_index = int(index[2])
    mutated = _signal_record(
        method="RK4",
        tracer=tracer,
        field_index=field_index,
        proper=proper["RK4"],
        fields=target,
        source_record=source_index["RK4"][(tracer, TDG3_FIELD_NAMES[field_index])],
    )
    original = next(
        record
        for record in baseline["signal_records"]["RK4"]
        if record["tracer"] == tracer
        and record["field"] == TDG3_FIELD_NAMES[field_index]
    )
    attacked_diagnosis = deepcopy(baseline)
    original_combined = attacked_diagnosis["signal_records"]["RK4"][0][
        "combined_classifications"
    ]["phase_field_tail_power"]
    attacked_diagnosis["signal_records"]["RK4"][0][
        "combined_classifications"
    ]["phase_field_tail_power"] = (
        TDG3_DISAGREE if original_combined != TDG3_DISAGREE else TDG3_ZERO
    )
    attacked_diagnosis["methods"]["RK4"] = _method_findings(
        "RK4", attacked_diagnosis["signal_records"]["RK4"]
    )
    attacked_diagnosis["cross_method"] = _cross_method_findings(
        attacked_diagnosis["signal_records"]
    )
    controls = {
        "one_bit_history_mutation_detected": mutated != original,
        "method_ladder_mutation_rejected": False,
        "combined_class_mutation_detected": False,
        "freeze_hash_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["immutable_lineage"].__setitem__(
                "freeze_result_sha256", "0" * 64
            ),
        ),
        "observed_count_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["observed_diagnosis"].__setitem__(
                "native_estimator_disagreement_count", 32
            ),
        ),
        "no_resampling_claim_mutation_rejected": _expect_config_rejection(
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
    }
    try:
        diagnose_native_tail_ladder(
            [record[:, tracer] for record in proper["RK4"]],
            [record[:, tracer, field_index] for record in fields["RK4"]],
            method="RK4",
            point_counts=(2049, 4097, 16385),
        )
    except (TypeError, ValueError):
        controls["method_ladder_mutation_rejected"] = True
    try:
        _compact_findings(attacked_diagnosis)
    except (TypeError, ValueError):
        controls["combined_class_mutation_detected"] = True
    if set(controls.values()) != {True}:
        raise ValueError(f"PREF18 mutation control failed: {controls}")
    return controls


def record(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    config = load_config(config_path)
    freeze_lineage, source_record = _verify_freeze_lineage(config)
    raw_paths = _raw_paths(config)
    if not _raw_bundle_is_complete_or_absent(raw_paths):
        raise FileNotFoundError("PREF18 raw campaign bundle is absent")
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
        raise ValueError("PREF18 raw campaign hash differs")

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
    ):
        raise ValueError("PREF18 CAL10 terminal restoration differs")
    source_diagnosis = source_record["artifact_payload"]["actual_history_diagnosis"]
    diagnosis = diagnose_actual_histories(histories, source_diagnosis)
    compact = _compact_findings(diagnosis)
    mutations = _mutation_controls(
        config, histories, source_diagnosis, diagnosis
    )
    result = {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "metric_signature": "-+++",
        "riemann_convention": "plus_partial_mu_gamma_nu",
        "classification": (
            "completed_TDG3_mixed_native_grid_history_diagnosis_"
            "sampling_identifiability_theorem_required"
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
            "immutable_freeze_lineage": freeze_lineage,
            "immutable_raw_bundle": {
                "raw_sha256": raw_hashes,
                "terminal_event_index": 63,
                "terminal_coordinate_time": 63 / 16,
                "history_arrays_loaded_read_only": True,
                "terminal_checkpoint_resumed": False,
                "state_advanced": False,
                "common_grid_resampling_performed": False,
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
                "historical_PROTO13_or_TDG2_result_reclassified": False,
                "TDG3_native_grid_outcome_is_mixed": True,
                "common_grid_interpolation_is_the_sole_failure_owner": False,
                "continuum_surrogate_enclosure_claimed": False,
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
    expected_classification = (
        "completed_TDG3_mixed_native_grid_history_diagnosis_"
        "sampling_identifiability_theorem_required"
    )
    expected_claim_boundary = {
        "historical_PROTO13_or_TDG2_result_reclassified": False,
        "TDG3_native_grid_outcome_is_mixed": True,
        "common_grid_interpolation_is_the_sole_failure_owner": False,
        "continuum_surrogate_enclosure_claimed": False,
        "individual_budget_obstruction_resolved": False,
        "replacement_temporal_admission_defined": False,
        "new_trajectory_authorized": False,
        "GR0_case_eligible": False,
        "candidate_branch_opened": False,
        "physical_question_answered": False,
    }
    if (
        value.get("schema_version") != 1
        or value.get("artifact_id") != ARTIFACT_ID
        or value.get("project_version") != PROJECT_VERSION
        or value.get("metric_signature") != "-+++"
        or value.get("riemann_convention") != "plus_partial_mu_gamma_nu"
        or value.get("classification") != expected_classification
        or value.get("generated_by") != "scripts/reproduce_fgc_tdg3_pref18.py"
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
        or payload.get("claim_boundary") != expected_claim_boundary
    ):
        raise ValueError("PREF18 stored result boundary differs")
    implementation = value.get("implementation_sha256", {})
    if (
        not isinstance(implementation, Mapping)
        or set(implementation) != {_rel(path) for path in IMPLEMENTATION}
        or any(
            _sha(REPOSITORY / relative) != expected
            for relative, expected in implementation.items()
        )
    ):
        raise ValueError("PREF18 implementation hash ledger differs")
    document = REPOSITORY / str(value.get("derivation_document", ""))
    if (
        value.get("derivation_document") != "docs/fgc-tdg3-pref18.md"
        or not document.is_file()
        or _sha(document) != value.get("derivation_document_sha256")
    ):
        raise ValueError("PREF18 derivation document differs")
    diagnosis = payload.get("actual_history_diagnosis", {})
    signal_records = diagnosis.get("signal_records", {})
    if (
        sum(len(records) for records in signal_records.values()) != 576
        or diagnosis.get("common_grid_resampling_performed") is not False
        or diagnosis.get("continuum_surrogate_enclosure_claimed") is not False
        or diagnosis.get("diagnosis_is_not_a_replacement_admission") is not True
        or diagnosis.get("frozen_estimator_names") != list(TDG3_ESTIMATOR_NAMES)
        or diagnosis.get("frozen_power_class_names") != list(TDG3_POWER_CLASS_NAMES)
        or diagnosis.get("frozen_combined_class_names")
        != list(TDG3_COMBINED_CLASS_NAMES)
        or diagnosis.get("frozen_measure_names") != list(TDG3_MEASURE_NAMES)
    ):
        raise ValueError("PREF18 complete signal evidence differs")
    for records in signal_records.values():
        for item in records:
            if (
                item.get("common_grid_resampling_performed") is not False
                or item.get("continuum_surrogate_enclosure_claimed") is not False
                or item.get("classification_is_not_a_replacement_admission")
                is not True
            ):
                raise ValueError("PREF18 signal claim boundary differs")
            for estimator in TDG3_ESTIMATOR_NAMES:
                observations = item.get("estimators", {}).get(estimator, {}).get(
                    "observations", []
                )
                if len(observations) != 3 or any(
                    record.get("common_grid_resampling_performed") is not False
                    or record.get("continuum_surrogate_enclosure_claimed") is not False
                    for record in observations
                ):
                    raise ValueError("PREF18 native observation boundary differs")
    mutations = payload.get("mutation_controls", {})
    if not isinstance(mutations, Mapping) or set(mutations.values()) != {True}:
        raise ValueError("PREF18 mutation controls differ")


def verify_canonical(
    config_path: Path = DEFAULT_CONFIG,
    output_path: Path = DEFAULT_OUTPUT,
) -> dict[str, Any]:
    config = load_config(config_path)
    raw = output_path.read_bytes()
    stored = json.loads(raw)
    if not isinstance(stored, dict) or raw != _canonical_bytes(stored):
        raise ValueError("PREF18 stored result is not canonical JSON")
    _validate_stored_record(stored, config)
    paths = _raw_paths(config)
    if _raw_bundle_is_complete_or_absent(paths):
        observed = record(config_path)
        if stored != observed:
            raise ValueError("PREF18 stored result differs from raw reproduction")
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
            f"{ARTIFACT_ID}: histories="
            f"{value['gate_status']['TDG3_actual_histories_diagnosed']} "
            f"replacement="
            f"{value['gate_status']['replacement_temporal_admission_defined']}"
        )
    else:
        print(_canonical(record(args.config)), end="")


if __name__ == "__main__":
    main()
