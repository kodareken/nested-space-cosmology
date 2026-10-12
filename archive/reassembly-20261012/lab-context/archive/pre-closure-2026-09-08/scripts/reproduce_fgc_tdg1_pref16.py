#!/usr/bin/env python3
"""Bind TDG1's read-only diagnosis of the immutable PROTO13 histories."""

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
from typing import Any, Callable, Mapping, Sequence

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts import reproduce_fgc_cal10_pref15 as cal10  # noqa: E402
from recursive_horizons.fgc.evolution.tdg1_temporal_diagnosis import (  # noqa: E402
    TDG1_FIELD_NAMES,
    diagnose_terminal_method_histories,
    diagnose_terminal_temporal_histories,
)


ARTIFACT_ID = "FGC-1-TDG1-PREF16"
PROJECT_VERSION = "0.11.0"
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-tdg1-pref16.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-tdg1-pref16.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-tdg1-pref16.md"
FREEZE_RESULT = REPOSITORY / "results/fgc-1-tdg1-frz1.json"
FREEZE_CONFIG = REPOSITORY / "configs/fgc/fgc-1-tdg1-frz1.toml"
DIAGNOSIS_MODULE = (
    REPOSITORY
    / "src/recursive_horizons/fgc/evolution/tdg1_temporal_diagnosis.py"
)
FREEZE_REPRODUCER = REPOSITORY / "scripts/reproduce_fgc_tdg1_frz1.py"
CAL10_REPRODUCER = REPOSITORY / "scripts/reproduce_fgc_cal10_pref15.py"
IMPLEMENTATION = (
    Path(__file__).resolve(),
    DIAGNOSIS_MODULE,
    CAL10_REPRODUCER,
)
METHOD_KEYS = {
    "RK4": ("RK4-2049", "RK4-4097", "RK4-8193"),
    "SSPRK3": ("SSPRK3-4097", "SSPRK3-8193", "SSPRK3-16385"),
}


EXPECTED_CLAIMS = {
    "TDG1_diagnostic_contract_frozen": True,
    "TDG1_synthetic_counterexample_confirmed": True,
    "TDG1_terminal_histories_diagnosed": True,
    "normalized_raw_temporal_tail_ratio_is_valid_general_spatial_convergence_test": False,
    "individual_temporal_budget_failure_cause_fully_derived": False,
    "replacement_temporal_admission_defined": False,
    "TDG2_frozen": False,
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
        "freeze_commit": "e8310bdbc893f76c67daf8fc993bdcd35259184d",
        "freeze_result_sha256": "11beff809d1c1af277e119dfc9fd54b571c2c7bd70de3da96c47e9945b28b7f3",
        "freeze_config_sha256": "4a1f41a40fedb4aecc2adafd79d90a5b96c81d9aaf5001bac1bcfd243d92e4e3",
        "diagnosis_module_sha256": "3aaf100a8c9abb691978ac3ddf1f58210ed59107dc6018cf3f8b70c99de8ba44",
        "freeze_reproducer_sha256": "00256721f72d696360407aab161f057ef0dbbc7a4795f87f8b63edd6a1ab8b2b",
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
        "resolution_count": 6,
        "signals_per_resolution": 288,
        "total_signal_resolution_records": 1728,
        "inherited_individual_failure_count": 544,
        "failed_below_declared_absolute_amplitude_floor": 149,
        "failed_above_declared_absolute_amplitude_floor": 395,
        "normalized_raw_tail_ratio_gate_structurally_invalid": True,
        "individual_budget_failure_cause_fully_derived": False,
        "history_difference_is_replacement_admission": False,
        "RK4": {
            "point_counts": [2049, 4097, 8193],
            "individual_failure_counts": [110, 106, 93],
            "subfloor_failure_counts": [3, 16, 27],
            "above_floor_failure_counts": [107, 90, 66],
            "all_view_persistent_failure_counts": [89, 77, 48],
            "history_difference_exactly_identical": 15,
            "history_difference_fine_pair_exact": 8,
            "history_difference_raw_order_pass": 237,
            "history_difference_raw_order_below_required": 28,
            "history_difference_nonmonotone": 0,
            "history_difference_interpolation_limited": 83,
        },
        "SSPRK3": {
            "point_counts": [4097, 8193, 16385],
            "individual_failure_counts": [96, 80, 59],
            "subfloor_failure_counts": [44, 37, 22],
            "above_floor_failure_counts": [52, 43, 37],
            "all_view_persistent_failure_counts": [44, 26, 22],
            "history_difference_exactly_identical": 33,
            "history_difference_fine_pair_exact": 1,
            "history_difference_raw_order_pass": 233,
            "history_difference_raw_order_below_required": 21,
            "history_difference_nonmonotone": 0,
            "history_difference_interpolation_limited": 51,
        },
    }


def validate_config_data(config: Mapping[str, Any]) -> None:
    _strict_keys(
        "PREF16 config",
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
        "freeze_result": "results/fgc-1-tdg1-frz1.json",
        "freeze_config": "configs/fgc/fgc-1-tdg1-frz1.toml",
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
        raise ValueError("PREF16 top-level contract differs")
    if config.get("scope") != {
        "target": "FGC-2-SF1-TDG1",
        "role": "post_freeze_read_only_terminal_history_diagnosis",
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
        raise ValueError("PREF16 scope differs")
    if _as_float(config["scope"]["source_terminal_coordinate_time"]) != 63 / 16:
        raise ValueError("PREF16 terminal time differs")
    if config.get("immutable_lineage") != _expected_lineage():
        raise ValueError("PREF16 immutable lineage differs")
    if config.get("observed_diagnosis") != _expected_observed():
        raise ValueError("PREF16 observed diagnosis differs")
    if config.get("conclusions") != {
        "historical_PROTO13_failure_preserved": True,
        "normalized_ratio_cause_derived_structurally": True,
        "individual_budget_obstruction_remains": True,
        "individual_failures_contract_with_resolution_but_do_not_vanish": True,
        "raw_history_differences_are_mostly_convergent_but_not_complete": True,
        "interpolation_limits_a_nonzero_subset": True,
        "absolute_tail_convergence_not_yet_diagnosed": True,
        "TDG2_absolute_tail_discriminator_design_may_begin": True,
        "replacement_gate_design_may_not_yet_begin": True,
        "PROTO14_design_may_not_yet_begin": True,
    }:
        raise ValueError("PREF16 conclusions differ")
    if set(config.get("proof_contract", {}).values()) != {True}:
        raise ValueError("PREF16 proof contract must be all-of")
    if config.get("successor_boundary") != {
        "TDG1_terminal_histories_diagnosed": True,
        "TDG2_absolute_tail_discriminator_design_may_begin": True,
        "TDG2_frozen": False,
        "replacement_temporal_gate_design_authorized": False,
        "replacement_temporal_admission_defined": False,
        "PROTO14_frozen": False,
        "fresh_GR0_dynamic_calibration_completed": False,
        "GR0_case_eligible": False,
        "SGBL_execution_authorized": False,
        "FGCQR_holdout_execution_authorized": False,
        "DEF1_execution_authorized": False,
    }:
        raise ValueError("PREF16 successor boundary differs")
    if config.get("claims") != EXPECTED_CLAIMS:
        raise ValueError("PREF16 claims differ")


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    with path.open("rb") as handle:
        config = tomllib.load(handle)
    validate_config_data(config)
    return config


def _raw_paths(config: Mapping[str, Any]) -> tuple[Path, Path, Path, Path]:
    return tuple(
        REPOSITORY / str(config[key])
        for key in (
            "campaign_manifest",
            "campaign_event_log",
            "campaign_checkpoint",
            "campaign_result",
        )
    )  # type: ignore[return-value]


def _raw_bundle_is_complete_or_absent(paths: Sequence[Path]) -> bool:
    present = [path.is_file() for path in paths]
    if not any(present):
        return False
    if not all(present):
        raise ValueError("PREF16 raw campaign bundle is partial")
    return True


def _verify_freeze_lineage(config: Mapping[str, Any]) -> dict[str, Any]:
    lineage = config["immutable_lineage"]
    commit = lineage["freeze_commit"]
    _git("cat-file", "-e", f"{commit}^{{commit}}")
    if _git("merge-base", "--is-ancestor", commit, "HEAD", check=False).returncode:
        raise ValueError("PREF16 freeze commit is not an ancestor of HEAD")
    tracked = {
        config["freeze_result"]: lineage["freeze_result_sha256"],
        config["freeze_config"]: lineage["freeze_config_sha256"],
        _rel(DIAGNOSIS_MODULE): lineage["diagnosis_module_sha256"],
        _rel(FREEZE_REPRODUCER): lineage["freeze_reproducer_sha256"],
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
            raise ValueError(f"PREF16 frozen tracked blob differs: {relative}")
    freeze = json.loads(FREEZE_RESULT.read_bytes())
    if (
        freeze.get("artifact_id") != "FGC-1-TDG1-FRZ1"
        or freeze.get("gate_status", {}).get(
            "TDG1_terminal_histories_diagnosed"
        )
        is not False
        or freeze.get("artifact_payload", {})
        .get("successor_boundary", {})
        .get("TDG1_actual_history_execution_authorized")
        is not True
    ):
        raise ValueError("PREF16 source freeze boundary differs")
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


def _compact_findings(diagnosis: Mapping[str, Any]) -> dict[str, Any]:
    methods = diagnosis["methods"]
    record: dict[str, Any] = {
        "method_count": len(methods),
        "resolution_count": sum(
            len(method["individual_signal_views_by_resolution"])
            for method in methods.values()
        ),
        "signals_per_resolution": 288,
        "total_signal_resolution_records": 1728,
        "inherited_individual_failure_count": 0,
        "failed_below_declared_absolute_amplitude_floor": diagnosis[
            "aggregate"
        ]["subfloor_inherited_individual_failure_count"],
        "failed_above_declared_absolute_amplitude_floor": diagnosis[
            "aggregate"
        ]["above_floor_inherited_individual_failure_count"],
        "normalized_raw_tail_ratio_gate_structurally_invalid": diagnosis[
            "aggregate"
        ]["inherited_raw_tail_ratio_gate_structurally_invalid"],
        "individual_budget_failure_cause_fully_derived": diagnosis[
            "aggregate"
        ]["individual_budget_failure_cause_fully_derived"],
        "history_difference_is_replacement_admission": False,
    }
    for method, expected_counts in (
        ("RK4", (2049, 4097, 8193)),
        ("SSPRK3", (4097, 8193, 16385)),
    ):
        item = methods[method]
        by_resolution = item["individual_signal_views_by_resolution"]
        if any(by_resolution[str(count)]["signal_count"] != 288 for count in expected_counts):
            raise ValueError(f"PREF16 {method} signal count differs")
        failures = [
            by_resolution[str(count)]["inherited_individual_failure_count"]
            for count in expected_counts
        ]
        record["inherited_individual_failure_count"] += sum(failures)
        classes = item["history_difference_classification_counts"]
        if sum(classes.values()) != 288:
            raise ValueError(f"PREF16 {method} history-difference count differs")
        record[method] = {
            "point_counts": list(expected_counts),
            "individual_failure_counts": failures,
            "subfloor_failure_counts": [
                by_resolution[str(count)][
                    "failed_below_declared_absolute_amplitude_floor"
                ]
                for count in expected_counts
            ],
            "above_floor_failure_counts": [
                by_resolution[str(count)][
                    "failed_above_declared_absolute_amplitude_floor"
                ]
                for count in expected_counts
            ],
            "all_view_persistent_failure_counts": [
                by_resolution[str(count)][
                    "failed_persisting_under_all_fixed_views"
                ]
                for count in expected_counts
            ],
            "history_difference_exactly_identical": classes[
                "exactly_resolution_identical"
            ],
            "history_difference_fine_pair_exact": classes[
                "fine_pair_exact_after_nonzero_coarse_difference"
            ],
            "history_difference_raw_order_pass": classes["raw_order_pass"],
            "history_difference_raw_order_below_required": classes[
                "raw_order_below_required"
            ],
            "history_difference_nonmonotone": classes[
                "raw_difference_nonmonotone"
            ],
            "history_difference_interpolation_limited": item[
                "history_difference_interpolation_limited_count"
            ],
        }
    if record != _expected_observed():
        raise ValueError("PREF16 compact diagnosis differs from frozen outcome")
    return record


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
    mutated_rk4 = [array.copy() for array in fields["RK4"]]
    target = mutated_rk4[0]
    index = np.unravel_index(int(np.argmax(np.abs(target))), target.shape)
    # Flip one high mantissa bit.  Lower mantissa flips are intentionally
    # swallowed by the published diagnostic's arithmetic enclosures, whereas
    # this bit remains finite and must alter the bound terminal diagnosis.
    target.view(np.uint64)[index] ^= np.uint64(1 << 44)
    changed = diagnose_terminal_method_histories(
        method="RK4",
        point_counts=(2049, 4097, 8193),
        proper_times_by_resolution=proper["RK4"],
        field_histories_by_resolution=mutated_rk4,
    )
    controls = {
        "one_bit_history_mutation_detected": changed != baseline["methods"]["RK4"],
        "method_ladder_mutation_rejected": False,
        "method_order_mutation_rejected": False,
        "freeze_hash_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["immutable_lineage"].__setitem__(
                "freeze_result_sha256", "0" * 64
            ),
        ),
        "observed_count_mutation_rejected": _expect_config_rejection(
            config,
            lambda value: value["observed_diagnosis"].__setitem__(
                "inherited_individual_failure_count", 543
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
            lambda value: value["claims"].__setitem__(
                "GR0_case_eligible", True
            ),
        ),
    }
    try:
        diagnose_terminal_method_histories(
            method="RK4",
            point_counts=(2049, 4097, 16385),
            proper_times_by_resolution=proper["RK4"],
            field_histories_by_resolution=fields["RK4"],
        )
    except (TypeError, ValueError):
        controls["method_ladder_mutation_rejected"] = True
    reversed_proper = {"SSPRK3": proper["SSPRK3"], "RK4": proper["RK4"]}
    reversed_fields = {"SSPRK3": fields["SSPRK3"], "RK4": fields["RK4"]}
    try:
        diagnose_terminal_temporal_histories(
            proper_times_by_method=reversed_proper,
            field_histories_by_method=reversed_fields,
        )
    except (TypeError, ValueError):
        controls["method_order_mutation_rejected"] = True
    if set(controls.values()) != {True}:
        raise ValueError(f"PREF16 mutation control failed: {controls}")
    return controls


def record(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    config = load_config(config_path)
    freeze_lineage = _verify_freeze_lineage(config)
    raw_paths = _raw_paths(config)
    if not _raw_bundle_is_complete_or_absent(raw_paths):
        raise FileNotFoundError("PREF16 raw campaign bundle is absent")
    lineage = config["immutable_lineage"]
    labels = ("manifest", "event_log", "checkpoint", "campaign_result")
    expected_hashes = {
        "manifest": lineage["manifest_sha256"],
        "event_log": lineage["event_log_sha256"],
        "checkpoint": lineage["checkpoint_sha256"],
        "campaign_result": lineage["campaign_result_sha256"],
    }
    raw_hashes = {
        label: _sha(path)
        for label, path in zip(labels, raw_paths, strict=True)
    }
    if raw_hashes != expected_hashes:
        raise ValueError("PREF16 raw campaign hash differs")

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
        raise ValueError("PREF16 CAL10 terminal restoration differs")
    proper, fields = _history_inputs(histories)
    diagnosis = diagnose_terminal_temporal_histories(
        proper_times_by_method=proper,
        field_histories_by_method=fields,
    )
    compact = _compact_findings(diagnosis)
    mutations = _mutation_controls(config, histories, diagnosis)
    result = {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "metric_signature": "-+++",
        "riemann_convention": "plus_partial_mu_gamma_nu",
        "classification": "completed_TDG1_history_diagnosis_TDG2_required",
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
                "normalized_ratio_instrument_defect_confirmed": True,
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
    compact = payload.get("compact_findings", {})
    boundary = payload.get("claim_boundary", {})
    successor = payload.get("successor_boundary", {})
    if (
        value.get("schema_version") != 1
        or value.get("artifact_id") != ARTIFACT_ID
        or value.get("project_version") != PROJECT_VERSION
        or value.get("classification")
        != "completed_TDG1_history_diagnosis_TDG2_required"
        or value.get("generated_by") != "scripts/reproduce_fgc_tdg1_pref16.py"
        or value.get("gate_status") != EXPECTED_CLAIMS
        or value.get("nonclaims")
        != {
            key: False
            for key, item in sorted(EXPECTED_CLAIMS.items())
            if item is False
        }
        or value.get("source_config_sha256")
        != {_rel(DEFAULT_CONFIG): _sha(DEFAULT_CONFIG)}
        or compact != _expected_observed()
        or payload.get("conclusions") != config["conclusions"]
        or successor != config["successor_boundary"]
        or boundary
        != {
            "historical_PROTO13_failure_reclassified": False,
            "normalized_ratio_instrument_defect_confirmed": True,
            "individual_budget_obstruction_resolved": False,
            "replacement_temporal_admission_defined": False,
            "new_trajectory_authorized": False,
            "GR0_case_eligible": False,
            "candidate_branch_opened": False,
            "physical_question_answered": False,
        }
    ):
        raise ValueError("PREF16 stored result boundary differs")
    implementation = value.get("implementation_sha256", {})
    if (
        not isinstance(implementation, Mapping)
        or set(implementation) != {_rel(path) for path in IMPLEMENTATION}
        or any(
            _sha(REPOSITORY / relative) != expected
            for relative, expected in implementation.items()
        )
    ):
        raise ValueError("PREF16 implementation hash ledger differs")
    document = REPOSITORY / str(value.get("derivation_document", ""))
    if (
        value.get("derivation_document") != "docs/fgc-tdg1-pref16.md"
        or not document.is_file()
        or _sha(document) != value.get("derivation_document_sha256")
    ):
        raise ValueError("PREF16 derivation document differs")
    mutations = payload.get("mutation_controls", {})
    if not isinstance(mutations, Mapping) or set(mutations.values()) != {True}:
        raise ValueError("PREF16 mutation controls differ")


def verify_canonical(
    config_path: Path = DEFAULT_CONFIG,
    output_path: Path = DEFAULT_OUTPUT,
) -> dict[str, Any]:
    config = load_config(config_path)
    raw = output_path.read_bytes()
    stored = json.loads(raw)
    if not isinstance(stored, dict) or raw != _canonical_bytes(stored):
        raise ValueError("PREF16 stored result is not canonical JSON")
    _validate_stored_record(stored, config)
    paths = _raw_paths(config)
    if _raw_bundle_is_complete_or_absent(paths):
        observed = record(config_path)
        if stored != observed:
            raise ValueError("PREF16 stored result differs from raw reproduction")
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
            f"{ARTIFACT_ID}: histories={value['gate_status']['TDG1_terminal_histories_diagnosed']} "
            f"replacement={value['gate_status']['replacement_temporal_admission_defined']}"
        )
    else:
        print(_canonical(record(args.config)), end="")


if __name__ == "__main__":
    main()
