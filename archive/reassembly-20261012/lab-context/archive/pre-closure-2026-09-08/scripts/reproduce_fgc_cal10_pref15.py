#!/usr/bin/env python3
"""Reproduce PREF15's independent terminal PROTO13 campaign certificate."""

from __future__ import annotations

import argparse
from copy import deepcopy
from dataclasses import asdict, is_dataclass
from fractions import Fraction
from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import tomllib
from typing import Any, Mapping, Sequence

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.cal10_proto13_campaign_diagnosis import (  # noqa: E402
    EXPECTED_FIELDS,
    EXPECTED_METHODS,
    EXPECTED_POINT_COUNTS,
    TEMPORAL_MAXIMUM_NESTED_TAIL_RATIO,
    diagnose_proto13_campaign,
)
from recursive_horizons.fgc.evolution.health_monitor import (  # noqa: E402
    compact_vacuum_buffer_window,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    EvolutionState,
    UniformRadialGrid,
    array_content_sha256,
)
from recursive_horizons.fgc.evolution.proto4_admission import (  # noqa: E402
    SpectralThresholds,
    windowed_spectral_power_budget,
)
from recursive_horizons.fgc.evolution.proto5_runtime import (  # noqa: E402
    proto5_temporal_spectral_admission,
)
from recursive_horizons.fgc.evolution.proto13_runtime import (  # noqa: E402
    proto13_gr0_common_event,
    proto13_method_owned_trapped_assessment,
)


Q = Fraction
ARTIFACT_ID = "FGC-1-CAL10-PREF15"
PROJECT_VERSION = "0.11.0"
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-cal10-pref15.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-cal10-pref15.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-cal10-pref15.md"
PRO13_FREEZE_RESULT = REPOSITORY / "results/fgc-1-pro13-frz1.json"
DIAGNOSIS_MODULE = (
    REPOSITORY
    / "src/recursive_horizons/fgc/evolution/cal10_proto13_campaign_diagnosis.py"
)
IMPLEMENTATION = (
    Path(__file__).resolve(),
    DIAGNOSIS_MODULE,
    REPOSITORY / "src/recursive_horizons/fgc/evolution/proto13_runtime.py",
    REPOSITORY / "src/recursive_horizons/fgc/evolution/proto5_runtime.py",
)


EXPECTED_CLAIMS = {
    "PROTO13_campaign_terminated_normally": True,
    "PROTO13_terminal_checkpoint_recomputed": True,
    "PROTO13_temporal_gate_reached": True,
    "PROTO13_terminal_common_spatial_and_constraint_admissions_passed": True,
    "PROTO13_temporal_admission_failed": True,
    "PROTO13_temporal_failure_cause_derived": False,
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


def _canonical_line(value: Mapping[str, Any]) -> bytes:
    return (
        json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        )
        + "\n"
    ).encode("utf-8")


def _serial(value: Any) -> Any:
    if is_dataclass(value):
        return _serial(asdict(value))
    if isinstance(value, Mapping):
        return {str(key): _serial(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_serial(item) for item in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


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


def _load_json(path: Path, *, canonical: bool = True) -> dict[str, Any]:
    raw = path.read_bytes()
    value = json.loads(raw)
    if not isinstance(value, dict):
        raise ValueError(f"{_rel(path)} must contain a JSON object")
    if canonical and raw != _canonical_bytes(value):
        raise ValueError(f"{_rel(path)} must contain canonical JSON")
    return value


def _load_events(path: Path) -> tuple[list[dict[str, Any]], bytes]:
    raw = path.read_bytes()
    if not raw.endswith(b"\n"):
        raise ValueError("PROTO13 event log must end in a newline")
    records: list[dict[str, Any]] = []
    rebuilt = bytearray()
    for line in raw.splitlines():
        value = json.loads(line)
        if not isinstance(value, dict):
            raise ValueError("PROTO13 event records must be objects")
        records.append(value)
        rebuilt.extend(_canonical_line(value))
    if bytes(rebuilt) != raw:
        raise ValueError("PROTO13 event log is not canonical JSONL")
    return records, raw


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
    return float(Q(value))


def validate_config_data(config: Mapping[str, Any]) -> None:
    _strict_keys(
        "PREF15 config",
        config,
        {
            "schema_version",
            "artifact_id",
            "project_version",
            "metric_signature",
            "riemann_convention",
            "run_plan_config",
            "runtime_authorization_result",
            "campaign_manifest",
            "campaign_event_log",
            "campaign_checkpoint",
            "campaign_result",
            "scope",
            "immutable_campaign",
            "method_owned_ladders",
            "temporal_stop",
            "proof_contract",
            "successor_boundary",
            "claims",
        },
    )
    expected_paths = {
        "run_plan_config": "configs/fgc/fgc-1-pro13-run1.toml",
        "runtime_authorization_result": "results/fgc-1-hlt11-mon11.json",
        "campaign_manifest": "runs/fgc-2-sf1/proto13/calibration/manifest.json",
        "campaign_event_log": "runs/fgc-2-sf1/proto13/calibration/events.jsonl",
        "campaign_checkpoint": (
            "runs/fgc-2-sf1/proto13/calibration/latest-checkpoint.npz"
        ),
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
        raise ValueError("PREF15 top-level configuration contract differs")

    scope = config.get("scope", {})
    immutable = config.get("immutable_campaign", {})
    ladders = config.get("method_owned_ladders", {})
    temporal = config.get("temporal_stop", {})
    proof = config.get("proof_contract", {})
    successor = config.get("successor_boundary", {})
    claims = config.get("claims", {})
    if (
        scope
        != {
            "target_protocol": "FGC-2-SF1-PROTO13",
            "role": "post_run_terminal_checkpoint_recomputed_GR0_temporal_gate_result",
            "calibration_branch": "GR-0",
            "amplitude": "3",
            "GR0_numerical_trajectory_read": True,
            "SGBL_trajectory_read": False,
            "FGCQR_trajectory_read": False,
            "mechanism_question_answered": False,
            "physical_obstruction_inferred_from_temporal_stop": False,
        }
        or immutable.get("authorization_commit")
        != "682e6cd91419853a4e5bf132c163326f9cc6a029"
        or immutable.get("campaign_id")
        != "36c0e18066cf1646e63ea4bb898f202dc05443dd190e9c053d69940f22665cb8"
        or immutable.get("event_log_line_count") != 41
        or immutable.get("event_log_byte_count") != 2474120
        or immutable.get("terminal_classification")
        != "calibration_failed_no_eligible_GR0_case"
        or immutable.get("terminal_event_index") != 63
        or _as_float(immutable.get("terminal_coordinate_time", "")) != 63 / 16
        or immutable.get("terminal_stop_reason")
        != "causal_past_temporal_spectral_admission"
        or immutable.get("resume_recovery_action")
        != "event_log_already_matches_checkpoint"
        or ladders
        != {
            "primary_method": "RK4",
            "primary_point_counts": [2049, 4097, 8193],
            "comparator_method": "SSPRK3",
            "comparator_point_counts": [4097, 8193, 16385],
            "common_physical_point_count": 2049,
            "terminal_common_spatial_and_constraint_admissions_passed": True,
            "terminal_trapped_sign_passed": False,
        }
        or temporal
        != {
            "minimum_sample_count": 64,
            "observed_sample_count": 64,
            "tracer_count": 48,
            "field_count": 6,
            "maximum_nested_tail_ratio": "1/4",
            "both_method_admissions_failed": True,
            "both_individual_budget_layers_failed": True,
            "both_nested_ratio_layers_failed": True,
            "RK4_failed_field_ratio_count": 510,
            "RK4_failed_derivative_ratio_count": 511,
            "SSPRK3_failed_field_ratio_count": 453,
            "SSPRK3_failed_derivative_ratio_count": 453,
            "threshold_may_not_be_changed_in_this_artifact": True,
            "numerical_or_physical_cause_derived": False,
        }
        or claims != EXPECTED_CLAIMS
        or set(proof.values()) != {True}
        or successor
        != {
            "temporal_gate_diagnosis_may_begin": True,
            "successor_protocol_frozen": False,
            "terminal_checkpoint_may_not_resume": True,
            "fresh_GR0_dynamic_calibration_completed": False,
            "GR0_case_eligible": False,
            "SGBL_execution_authorized": False,
            "FGCQR_holdout_execution_authorized": False,
            "DEF1_execution_authorized": False,
            "classical_spherical_diagnostic_authorized": False,
        }
    ):
        raise ValueError("PREF15 configuration contract differs")
    for key in (
        "manifest_sha256",
        "event_log_sha256",
        "checkpoint_sha256",
        "campaign_result_sha256",
        "resume_checkpoint_event_log_sha256",
    ):
        value = immutable.get(key)
        if not isinstance(value, str) or len(value) != 64:
            raise ValueError(f"PREF15 immutable {key} differs")


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
        raise ValueError("PREF15 raw campaign bundle is partial")
    return True


def _verify_authorization_commit(
    config: Mapping[str, Any], manifest: Mapping[str, Any]
) -> dict[str, Any]:
    commit = config["immutable_campaign"]["authorization_commit"]
    _git("cat-file", "-e", f"{commit}^{{commit}}")
    if _git("merge-base", "--is-ancestor", commit, "HEAD", check=False).returncode != 0:
        raise ValueError("PREF15 authorization commit is not an ancestor of HEAD")
    if manifest.get("authorization_checkpoint_commit") != commit:
        raise ValueError("PREF15 manifest authorization commit differs")
    paths = (
        config["run_plan_config"],
        config["runtime_authorization_result"],
        _rel(PRO13_FREEZE_RESULT),
        "scripts/run_fgc_gr0_calibration_v13.py",
        "src/recursive_horizons/fgc/evolution/proto13_runtime.py",
    )
    blob_hashes: dict[str, str] = {}
    for relative in paths:
        blob = _git("show", f"{commit}:{relative}").stdout
        blob_hashes[relative] = _bytes_sha(blob)
    if (
        blob_hashes[config["run_plan_config"]] != manifest.get("plan_sha256")
        or blob_hashes[config["runtime_authorization_result"]]
        != manifest.get("authorization_sha256")
        or blob_hashes["scripts/run_fgc_gr0_calibration_v13.py"]
        != manifest.get("implementation_sha256", {}).get(
            "scripts/run_fgc_gr0_calibration_v13.py"
        )
        or blob_hashes["src/recursive_horizons/fgc/evolution/proto13_runtime.py"]
        != manifest.get("implementation_sha256", {}).get(
            "src/recursive_horizons/fgc/evolution/proto13_runtime.py"
        )
    ):
        raise ValueError("PREF15 immutable authorization blob binding differs")
    _verify_immutable_pro13_freeze(config)
    return {
        "authorization_commit": commit,
        "authorization_commit_is_ancestor_of_HEAD": True,
        "historical_PRO13_freeze_verified_at_authorization_commit": True,
        "authorization_blob_sha256": blob_hashes,
    }


def _verify_immutable_pro13_freeze(
    config: Mapping[str, Any],
) -> dict[str, Any]:
    """Verify PRO13's historical namespace observation after its successor ran.

    Namespace absence was observed before launch and cannot truthfully be
    re-observed after the authorized campaign created the calibration root.
    The immutable authorization commit must instead contain the exact tracked
    PRO13 record that is still present now.
    """

    commit = config["immutable_campaign"]["authorization_commit"]
    relative = _rel(PRO13_FREEZE_RESULT)
    committed = _git("show", f"{commit}:{relative}").stdout
    current = PRO13_FREEZE_RESULT.read_bytes()
    if committed != current:
        raise ValueError("PREF15 historical PRO13 freeze differs from authorization commit")
    record = json.loads(current)
    namespace = record.get("artifact_payload", {}).get(
        "namespace_precondition", {}
    )
    status = record.get("gate_status", {})
    if (
        current != _canonical_bytes(record)
        or record.get("artifact_id") != "FGC-1-PRO13-FRZ1"
        or namespace
        != {
            "both_new_namespaces_absent": True,
            "freeze_created_no_namespace": True,
            "records": [
                {
                    "absent_before_freeze": True,
                    "created_by_freeze": False,
                    "path": "runs/fgc-2-sf1/proto13/calibration",
                },
                {
                    "absent_before_freeze": True,
                    "created_by_freeze": False,
                    "path": "runs/fgc-2-sf1/proto13/holdout",
                },
            ],
        }
        or status.get("PROTO13_frozen") is not True
        or status.get("PROTO13_successor_runtime_implemented") is not False
        or status.get("FGCQR_holdout_execution_authorized") is not False
    ):
        raise ValueError("PREF15 historical PRO13 freeze boundary differs")
    return record


def verify_immutable_pro13_freeze(
    config_path: Path = DEFAULT_CONFIG,
) -> dict[str, Any]:
    """Public successor-time verifier for the immutable PRO13 freeze."""

    return _verify_immutable_pro13_freeze(load_config(config_path))


def _restore_terminal_checkpoint(
    path: Path,
    *,
    event_log_raw: bytes,
    campaign_result: Mapping[str, Any],
) -> tuple[
    dict[str, Any],
    dict[str, EvolutionState],
    dict[str, UniformRadialGrid],
    dict[str, dict[str, np.ndarray]],
]:
    expected_members = {
        f"{method}-{count}"
        for method, counts in EXPECTED_POINT_COUNTS.items()
        for count in counts
    }
    states: dict[str, EvolutionState] = {}
    grids: dict[str, UniformRadialGrid] = {}
    histories: dict[str, dict[str, np.ndarray]] = {}
    with np.load(path, allow_pickle=False) as archive:
        metadata = json.loads(bytes(archive["metadata_utf8"].tolist()).decode("utf-8"))
        embedded_log = bytes(archive["event_log_utf8"].tolist())
        if (
            metadata.get("terminal") is not True
            or metadata.get("completed_common_event_index") != 63
            or metadata.get("consecutive_qualified_events") != 0
            or metadata.get("terminal_result") != campaign_result
            or embedded_log != event_log_raw
            or metadata.get("event_log")
            != {
                "sha256": _bytes_sha(event_log_raw),
                "line_count": 41,
                "byte_count": len(event_log_raw),
            }
            or set(metadata.get("members", {})) != expected_members
        ):
            raise ValueError("PREF15 terminal checkpoint identity differs")
        result_hashes = campaign_result["amplitude_records"][0][
            "member_final_state_hashes"
        ]
        for key in sorted(expected_members):
            stored = metadata["members"][key]
            prefix = key.replace("-", "_")
            state = EvolutionState(
                archive[f"{prefix}_u"].copy(),
                archive[f"{prefix}_p"].copy(),
                archive[f"{prefix}_q"].copy(),
            )
            observed_hash = array_content_sha256(state.u, state.p, state.q)
            point_count = int(key.rsplit("-", 1)[1])
            if (
                state.shape != (point_count, 6)
                or stored.get("state_sha256") != observed_hash
                or result_hashes.get(key) != observed_hash
                or stored.get("time") != 63 / 16
                or stored.get("runtime_monitor_state", {}).get(
                    "accepted_stage_count"
                )
                != stored.get("transaction_serial")
                or stored.get("source_retry_count") != 0
            ):
                raise ValueError(f"PREF15 terminal member differs: {key}")
            event_times = archive[f"{prefix}_event_proper_times"].copy()
            event_fields = archive[f"{prefix}_event_fields"].copy()
            tracer_positions = archive[f"{prefix}_tracer_positions"].copy()
            tracer_proper_times = archive[f"{prefix}_tracer_proper_times"].copy()
            if (
                event_times.shape != (64, 48)
                or event_fields.shape != (64, 48, 6)
                or tracer_positions.shape != (48,)
                or tracer_proper_times.shape != (48,)
                or np.any(np.diff(event_times, axis=0) <= 0.0)
                or not all(
                    np.all(np.isfinite(value))
                    for value in (
                        event_times,
                        event_fields,
                        tracer_positions,
                        tracer_proper_times,
                    )
                )
            ):
                raise ValueError(f"PREF15 terminal tracer history differs: {key}")
            states[key] = state
            grids[key] = UniformRadialGrid(0.0, 128.0, point_count)
            histories[key] = {
                "event_proper_times": event_times,
                "event_fields": event_fields,
                "tracer_positions": tracer_positions,
                "tracer_proper_times": tracer_proper_times,
            }
    return metadata, states, grids, histories


def _require_state_hash(state: EvolutionState, expected: str) -> None:
    if array_content_sha256(state.u, state.p, state.q) != expected:
        raise ValueError("PREF15 terminal state-content hash differs")


def _method_keys(method: str) -> tuple[str, str, str]:
    return tuple(f"{method}-{count}" for count in EXPECTED_POINT_COUNTS[method])  # type: ignore[return-value]


def _recompute_terminal_assessment(
    metadata: Mapping[str, Any],
    states: Mapping[str, EvolutionState],
    grids: Mapping[str, UniformRadialGrid],
    histories: Mapping[str, Mapping[str, np.ndarray]],
) -> dict[str, Any]:
    common: dict[str, Any] = {}
    temporal: dict[str, Any] = {}
    for method in EXPECTED_METHODS:
        keys = _method_keys(method)
        common[method] = proto13_gr0_common_event(
            [states[key] for key in keys],
            [grids[key] for key in keys],
            accepted_stage_counts=[
                metadata["members"][key]["runtime_monitor_state"][
                    "accepted_stage_count"
                ]
                for key in keys
            ],
            method=method,
            coordinate_time=63 / 16,
            cutoff=16.0,
            measurement_radius_maximum=24.0,
            taper_fraction=1.0 / 8.0,
            fixed_outer_rows=4,
        )
        temporal[method] = proto5_temporal_spectral_admission(
            point_counts=EXPECTED_POINT_COUNTS[method],
            proper_times_by_resolution=[
                histories[key]["event_proper_times"] for key in keys
            ],
            field_histories_by_resolution=[
                histories[key]["event_fields"] for key in keys
            ],
            field_names=EXPECTED_FIELDS,
            cutoff=16.0,
            taper_fraction=1.0 / 8.0,
        )
    trapped = proto13_method_owned_trapped_assessment(
        primary_states=[states[key] for key in _method_keys("RK4")],
        primary_grids=[grids[key] for key in _method_keys("RK4")],
        comparator_states=[states[key] for key in _method_keys("SSPRK3")],
        comparator_grids=[grids[key] for key in _method_keys("SSPRK3")],
        primary_common_event=common["RK4"],
        comparator_common_event=common["SSPRK3"],
        measurement_radius_maximum=24.0,
        minimum_observed_order=1.5,
        positive_margin_factor=4.0,
    )
    return {
        "primary_common_event": _serial(common["RK4"]),
        "comparator_common_event": _serial(common["SSPRK3"]),
        "temporal_spectral_admission": {
            "available": True,
            "sample_count": 64,
            "methods": {method: _serial(temporal[method]) for method in EXPECTED_METHODS},
            "admission_passed": False,
        },
        "trapped_assessment": _serial(trapped),
    }


def _individual_budget_summary(
    histories: Mapping[str, Mapping[str, np.ndarray]],
) -> dict[str, Any]:
    limits = SpectralThresholds()
    result: dict[str, Any] = {}
    for method in EXPECTED_METHODS:
        resolutions: dict[str, Any] = {}
        for key in _method_keys(method):
            times = histories[key]["event_proper_times"]
            fields = histories[key]["event_fields"]
            failed_by_field = {field: 0 for field in EXPECTED_FIELDS}
            causes = {"field_tail": 0, "derivative_tail": 0, "RMS_scale": 0}
            maximum = {
                "field_tail_fraction": 0.0,
                "derivative_tail_fraction": 0.0,
                "RMS_scale_over_cutoff": 0.0,
            }
            smallest_failed_total_power: float | None = None
            largest_failed_total_power = 0.0
            failed_signal_count = 0
            for tracer in range(48):
                proper = times[:, tracer]
                uniform = np.linspace(proper[0], proper[-1], proper.size)
                window = compact_vacuum_buffer_window(
                    uniform,
                    support_minimum=float(uniform[0]),
                    support_maximum=float(uniform[-1]),
                    taper_width=float((uniform[-1] - uniform[0]) / 8.0),
                )
                for field_index, field in enumerate(EXPECTED_FIELDS):
                    values = np.interp(
                        uniform, proper, fields[:, tracer, field_index]
                    )
                    budget = windowed_spectral_power_budget(
                        values, uniform, cutoff=16.0, window=window
                    )
                    maximum["field_tail_fraction"] = max(
                        maximum["field_tail_fraction"],
                        budget.top_band_field_power_fraction,
                    )
                    maximum["derivative_tail_fraction"] = max(
                        maximum["derivative_tail_fraction"],
                        budget.top_band_derivative_weighted_power_fraction,
                    )
                    maximum["RMS_scale_over_cutoff"] = max(
                        maximum["RMS_scale_over_cutoff"],
                        budget.rms_scale_over_cutoff,
                    )
                    local = False
                    if (
                        budget.top_band_field_power_fraction
                        >= limits.maximum_top_band_field_power_fraction
                    ):
                        causes["field_tail"] += 1
                        local = True
                    if (
                        budget.top_band_derivative_weighted_power_fraction
                        >= limits.maximum_top_band_derivative_power_fraction
                    ):
                        causes["derivative_tail"] += 1
                        local = True
                    if budget.rms_scale_over_cutoff >= limits.maximum_rms_scale_over_cutoff:
                        causes["RMS_scale"] += 1
                        local = True
                    if local:
                        failed_signal_count += 1
                        failed_by_field[field] += 1
                        largest_failed_total_power = max(
                            largest_failed_total_power, budget.total_field_power
                        )
                        if smallest_failed_total_power is None:
                            smallest_failed_total_power = budget.total_field_power
                        else:
                            smallest_failed_total_power = min(
                                smallest_failed_total_power, budget.total_field_power
                            )
            grid_count = int(key.rsplit("-", 1)[1])
            resolutions[str(grid_count)] = {
                "point_count": grid_count,
                "signal_count": 48 * len(EXPECTED_FIELDS),
                "failed_signal_count": failed_signal_count,
                "failed_by_field": failed_by_field,
                "failure_cause_counts": causes,
                "maximum_observed": maximum,
                "smallest_failed_total_field_power": smallest_failed_total_power,
                "largest_failed_total_field_power": largest_failed_total_power,
            }
        result[method] = resolutions
    return result


def _mutation_controls(
    config: Mapping[str, Any],
    events: Sequence[Mapping[str, Any]],
    campaign_result: Mapping[str, Any],
    expected_terminal: Mapping[str, Any],
    metadata: Mapping[str, Any],
    states: Mapping[str, EvolutionState],
    histories: Mapping[str, Mapping[str, np.ndarray]],
) -> dict[str, bool]:
    controls: dict[str, bool] = {}

    changed_threshold = deepcopy(config)
    changed_threshold["temporal_stop"]["maximum_nested_tail_ratio"] = "1/3"
    try:
        validate_config_data(changed_threshold)
    except (TypeError, ValueError):
        controls["threshold"] = True
    else:
        controls["threshold"] = False

    changed_claim = deepcopy(config)
    changed_claim["claims"]["FGCQR_holdout_execution_authorized"] = True
    try:
        validate_config_data(changed_claim)
    except (TypeError, ValueError):
        controls["claim"] = True
    else:
        controls["claim"] = False

    changed_event = deepcopy(list(events))
    changed_event[-1]["assessment"]["temporal_spectral_admission"][
        "admission_passed"
    ] = True
    try:
        diagnose_proto13_campaign(changed_event, campaign_result)
    except (TypeError, ValueError):
        controls["terminal_event"] = True
    else:
        controls["terminal_event"] = False

    changed_method = deepcopy(list(events))
    changed_method[-1]["assessment"]["temporal_spectral_admission"]["methods"][
        "SSPRK3"
    ]["point_counts"][-1] = 8193
    try:
        diagnose_proto13_campaign(changed_method, campaign_result)
    except (TypeError, ValueError):
        controls["method_ladder"] = True
    else:
        controls["method_ladder"] = False

    changed_result = deepcopy(campaign_result)
    changed_result["FGCQR_outcome_read"] = True
    try:
        diagnose_proto13_campaign(events, changed_result)
    except (TypeError, ValueError):
        controls["campaign_claim"] = True
    else:
        controls["campaign_claim"] = False

    mutated_u = states["RK4-2049"].u.copy()
    mutated_u.view(np.uint64).flat[0] ^= np.uint64(1)
    mutated_state = EvolutionState(
        mutated_u,
        states["RK4-2049"].p.copy(),
        states["RK4-2049"].q.copy(),
    )
    try:
        _require_state_hash(
            mutated_state, metadata["members"]["RK4-2049"]["state_sha256"]
        )
    except ValueError:
        controls["one_bit_state"] = True
    else:
        controls["one_bit_state"] = False

    changed_histories = {
        key: {name: value.copy() for name, value in record.items()}
        for key, record in histories.items()
    }
    history_bits = changed_histories["RK4-2049"]["event_fields"].view(np.uint64)
    middle = changed_histories["RK4-2049"]["event_fields"][32]
    tracer_index, field_index = np.unravel_index(
        int(np.argmax(np.abs(middle))), middle.shape
    )
    history_bits[32, tracer_index, field_index] ^= np.uint64(1 << 48)
    changed_temporal = proto5_temporal_spectral_admission(
        point_counts=EXPECTED_POINT_COUNTS["RK4"],
        proper_times_by_resolution=[
            changed_histories[key]["event_proper_times"]
            for key in _method_keys("RK4")
        ],
        field_histories_by_resolution=[
            changed_histories[key]["event_fields"] for key in _method_keys("RK4")
        ],
        field_names=EXPECTED_FIELDS,
        cutoff=16.0,
        taper_fraction=1.0 / 8.0,
    )
    controls["history_recomputation"] = (
        _serial(changed_temporal)
        != expected_terminal["temporal_spectral_admission"]["methods"]["RK4"]
    )
    return controls


def record(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    config = load_config(config_path)
    manifest_path, event_path, checkpoint_path, result_path = _raw_paths(config)
    if not _raw_bundle_is_complete_or_absent(
        (manifest_path, event_path, checkpoint_path, result_path)
    ):
        raise FileNotFoundError("PREF15 raw campaign bundle is absent")
    immutable = config["immutable_campaign"]
    observed_hashes = {
        "manifest": _sha(manifest_path),
        "event_log": _sha(event_path),
        "checkpoint": _sha(checkpoint_path),
        "campaign_result": _sha(result_path),
    }
    expected_hashes = {
        "manifest": immutable["manifest_sha256"],
        "event_log": immutable["event_log_sha256"],
        "checkpoint": immutable["checkpoint_sha256"],
        "campaign_result": immutable["campaign_result_sha256"],
    }
    if observed_hashes != expected_hashes:
        raise ValueError("PREF15 raw campaign hash differs")

    manifest = _load_json(manifest_path)
    events, event_raw = _load_events(event_path)
    campaign_result = _load_json(result_path)
    authorization = _verify_authorization_commit(config, manifest)
    if (
        manifest.get("runner_id") != "FGC-1-CAL10-RUN1-RUNNER"
        or manifest.get("campaign_id") != immutable["campaign_id"]
        or manifest.get("outcome_fields_present_at_creation") is not False
        or manifest.get("restart_coordinate_time") != 23 / 16
        or manifest.get("numerical_runtime_contract")
        != campaign_result.get("numerical_runtime_contract")
        or manifest.get("runtime_environment")
        != campaign_result.get("runtime_environment")
    ):
        raise ValueError("PREF15 manifest/result runtime boundary differs")
    metadata, states, grids, histories = _restore_terminal_checkpoint(
        checkpoint_path,
        event_log_raw=event_raw,
        campaign_result=campaign_result,
    )
    recomputed = _recompute_terminal_assessment(metadata, states, grids, histories)
    terminal_assessment = events[-1]["assessment"]
    for key in (
        "primary_common_event",
        "comparator_common_event",
        "temporal_spectral_admission",
        "trapped_assessment",
    ):
        if recomputed[key] != terminal_assessment[key]:
            raise ValueError(f"PREF15 terminal {key} recomputation differs")
    diagnosis = diagnose_proto13_campaign(events, campaign_result)
    if (
        diagnosis["temporal_methods"]["RK4"]["ratio_summaries"][
            "field_power_tail_ratios"
        ]["failed_comparison_count"]
        != config["temporal_stop"]["RK4_failed_field_ratio_count"]
        or diagnosis["temporal_methods"]["RK4"]["ratio_summaries"][
            "derivative_power_tail_ratios"
        ]["failed_comparison_count"]
        != config["temporal_stop"]["RK4_failed_derivative_ratio_count"]
        or diagnosis["temporal_methods"]["SSPRK3"]["ratio_summaries"][
            "field_power_tail_ratios"
        ]["failed_comparison_count"]
        != config["temporal_stop"]["SSPRK3_failed_field_ratio_count"]
        or diagnosis["temporal_methods"]["SSPRK3"]["ratio_summaries"][
            "derivative_power_tail_ratios"
        ]["failed_comparison_count"]
        != config["temporal_stop"]["SSPRK3_failed_derivative_ratio_count"]
    ):
        raise ValueError("PREF15 temporal failure counts differ")

    member_records = {}
    for key in sorted(states):
        stored = metadata["members"][key]
        member_records[key] = {
            "point_count": grids[key].point_count,
            "coordinate_time": stored["time"],
            "step_index": stored["step_index"],
            "transaction_serial": stored["transaction_serial"],
            "accepted_stage_count": stored["runtime_monitor_state"][
                "accepted_stage_count"
            ],
            "source_retry_count": stored["source_retry_count"],
            "CFL_retry_count": stored["CFL_retry_count"],
            "state_sha256": stored["state_sha256"],
            "event_history_shape": list(histories[key]["event_fields"].shape),
            "proper_time_history_shape": list(
                histories[key]["event_proper_times"].shape
            ),
            "accumulated_characteristic_distance": stored["causal_state"][
                "accumulated_characteristic_distance"
            ],
        }

    mutations = _mutation_controls(
        config,
        events,
        campaign_result,
        recomputed,
        metadata,
        states,
        histories,
    )
    if set(mutations.values()) != {True}:
        raise ValueError(f"PREF15 mutation control did not fail closed: {mutations}")
    result = {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "classification": "completed_PROTO13_GR0_temporal_gate_stop",
        "generated_by": "scripts/reproduce_fgc_cal10_pref15.py",
        "metric_signature": "-+++",
        "riemann_convention": "plus_partial_mu_gamma_nu",
        "gate_status": dict(EXPECTED_CLAIMS),
        "nonclaims": {
            key: False for key, value in sorted(EXPECTED_CLAIMS.items()) if value is False
        },
        "source_config_sha256": {
            _rel(config_path): _sha(config_path),
            config["run_plan_config"]: authorization[
                "authorization_blob_sha256"
            ][config["run_plan_config"]],
            config["runtime_authorization_result"]: authorization[
                "authorization_blob_sha256"
            ][config["runtime_authorization_result"]],
        },
        "implementation_sha256": {
            _rel(path): _sha(path) for path in IMPLEMENTATION
        },
        "artifact_payload": {
            "immutable_campaign": {
                **authorization,
                "campaign_id": manifest["campaign_id"],
                "raw_sha256": observed_hashes,
                "event_log_line_count": len(events),
                "event_log_byte_count": len(event_raw),
                "manifest_runtime_environment": manifest[
                    "runtime_environment"
                ],
            },
            "terminal_checkpoint": {
                "terminal": True,
                "completed_common_event_index": metadata[
                    "completed_common_event_index"
                ],
                "consecutive_qualified_events": metadata[
                    "consecutive_qualified_events"
                ],
                "event_log_identity": metadata["event_log"],
                "terminal_result_matches_external": True,
                "event_log_matches_external": True,
                "members": member_records,
            },
            "campaign_diagnosis": diagnosis,
            "independent_terminal_recomputation": {
                "primary_common_event_matches_terminal_event": True,
                "comparator_common_event_matches_terminal_event": True,
                "trapped_assessment_matches_terminal_event": True,
                "temporal_admission_matches_checkpoint_histories": True,
                "all_serialized_temporal_ratios_match_recomputation": True,
                "terminal_primary_admission_passed": recomputed[
                    "primary_common_event"
                ]["admission_passed"],
                "terminal_comparator_admission_passed": recomputed[
                    "comparator_common_event"
                ]["admission_passed"],
                "terminal_trapped_sign_passed": recomputed[
                    "trapped_assessment"
                ]["trapped_sign_passed"],
                "terminal_temporal_admission_passed": recomputed[
                    "temporal_spectral_admission"
                ]["admission_passed"],
            },
            "individual_temporal_budget_diagnosis": _individual_budget_summary(
                histories
            ),
            "resume_provenance": campaign_result["resume_recovery"],
            "mutation_controls": mutations,
            "successor_boundary": dict(config["successor_boundary"]),
            "claim_boundary": {
                "campaign_terminated_normally_but_calibration_not_completed": True,
                "temporal_stop_is_not_a_candidate_action_result": True,
                "temporal_stop_is_not_a_physical_obstruction": True,
                "terminal_checkpoint_cannot_be_resumed": True,
                "threshold_change_or_post_outcome_reclassification_authorized": False,
                "SGBL_or_FGCQR_outcome_read": False,
                "mechanism_question_answered": False,
            },
        },
        "derivation_document": _rel(OWNER_DOCUMENT),
        "derivation_document_sha256": _sha(OWNER_DOCUMENT),
    }
    return _serial(result)


def _validate_stored_record(record: Mapping[str, Any], config: Mapping[str, Any]) -> None:
    payload = record.get("artifact_payload", {})
    campaign = payload.get("immutable_campaign", {})
    diagnosis = payload.get("campaign_diagnosis", {})
    recomputed = payload.get("independent_terminal_recomputation", {})
    successor = payload.get("successor_boundary", {})
    boundary = payload.get("claim_boundary", {})
    raw = campaign.get("raw_sha256", {})
    immutable = config["immutable_campaign"]
    if (
        record.get("schema_version") != 1
        or record.get("artifact_id") != ARTIFACT_ID
        or record.get("project_version") != PROJECT_VERSION
        or record.get("classification")
        != "completed_PROTO13_GR0_temporal_gate_stop"
        or record.get("generated_by") != "scripts/reproduce_fgc_cal10_pref15.py"
        or record.get("gate_status") != EXPECTED_CLAIMS
        or record.get("nonclaims")
        != {
            key: False
            for key, value in sorted(EXPECTED_CLAIMS.items())
            if value is False
        }
        or record.get("derivation_document") != _rel(OWNER_DOCUMENT)
        or record.get("derivation_document_sha256") != _sha(OWNER_DOCUMENT)
        or raw
        != {
            "manifest": immutable["manifest_sha256"],
            "event_log": immutable["event_log_sha256"],
            "checkpoint": immutable["checkpoint_sha256"],
            "campaign_result": immutable["campaign_result_sha256"],
        }
        or campaign.get("authorization_commit")
        != immutable["authorization_commit"]
        or diagnosis.get("terminal_event_index") != 63
        or diagnosis.get("terminal_stop", {}).get("reason")
        != "causal_past_temporal_spectral_admission"
        or diagnosis.get("both_temporal_method_admissions_failed") is not True
        or diagnosis.get("temporal_failure_cause_derived") is not False
        or recomputed.get("temporal_admission_matches_checkpoint_histories")
        is not True
        or recomputed.get("terminal_temporal_admission_passed") is not False
        or set(payload.get("mutation_controls", {}).values()) != {True}
        or successor != config["successor_boundary"]
        or boundary
        != {
            "campaign_terminated_normally_but_calibration_not_completed": True,
            "temporal_stop_is_not_a_candidate_action_result": True,
            "temporal_stop_is_not_a_physical_obstruction": True,
            "terminal_checkpoint_cannot_be_resumed": True,
            "threshold_change_or_post_outcome_reclassification_authorized": False,
            "SGBL_or_FGCQR_outcome_read": False,
            "mechanism_question_answered": False,
        }
    ):
        raise ValueError("PREF15 stored certificate boundary differs")


def verify_canonical(
    config_path: Path = DEFAULT_CONFIG,
    output_path: Path = DEFAULT_OUTPUT,
) -> dict[str, Any]:
    config = load_config(config_path)
    raw_complete = _raw_bundle_is_complete_or_absent(_raw_paths(config))
    stored = _load_json(output_path)
    _validate_stored_record(stored, config)
    if raw_complete:
        reproduced = record(config_path)
        if reproduced != stored:
            raise ValueError("PREF15 canonical result does not reproduce")
    return stored


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--check-pro13-freeze", action="store_true")
    args = parser.parse_args()
    if args.check_pro13_freeze:
        result = verify_immutable_pro13_freeze(args.config.resolve())
        print(
            "FGC-1-PRO13-FRZ1: historical namespace evidence verified at "
            "immutable authorization commit"
        )
        return
    if args.check:
        result = verify_canonical(args.config.resolve(), args.output.resolve())
    else:
        result = record(args.config.resolve())
        _atomic_write(args.output.resolve(), _canonical_bytes(result))
    print(
        f"{ARTIFACT_ID}: terminal={result['gate_status']['PROTO13_campaign_terminated_normally']} "
        f"eligible={result['gate_status']['GR0_case_eligible']} "
        f"temporal_cause={result['gate_status']['PROTO13_temporal_failure_cause_derived']}"
    )


if __name__ == "__main__":
    main()
