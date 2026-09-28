#!/usr/bin/env python3
"""Reproduce the CAL4 diagnosis of the PROTO7 common-event stop.

Fast verification binds the immutable PROTO7 campaign, replay cache, stored
certificate, source, and owner document.  ``--replay`` reconstructs the
amplitude-5/2 slice from the immutable t=0 inputs and compares it with both the
public event and the stored diagnosis.  No candidate branch is constructible.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import tempfile
import tomllib
from typing import Any, Mapping, Sequence

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
import sys

sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts import run_fgc_gr0_calibration as inherited  # noqa: E402
from scripts import run_fgc_gr0_calibration_v7 as runner  # noqa: E402
from recursive_horizons.fgc.evolution.cal4_common_event_diagnosis import (  # noqa: E402
    finest_pair_nested_spectral_admission,
    owned_constraint_admission,
    owned_semidiscrete_constraint_snapshot,
)
from recursive_horizons.fgc.evolution.health_monitor import (  # noqa: E402
    compact_vacuum_buffer_window,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    EvolutionState,
    array_content_sha256,
)
from recursive_horizons.fgc.evolution.proto4_admission import (  # noqa: E402
    PROTO4_CONSTRAINT_ORDER,
    SpectralThresholds,
    spectral_field_budgets,
)
from recursive_horizons.fgc.evolution.proto5_runtime import (  # noqa: E402
    proper_radial_profile,
)


ARTIFACT_ID = "FGC-1-CAL4-PREF6"
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-cal4-pref6.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-cal4-pref6.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-cal4-pref6.md"
DIAGNOSIS_MODULE = (
    REPOSITORY
    / "src/recursive_horizons/fgc/evolution/cal4_common_event_diagnosis.py"
)
IMPLEMENTATION = (
    Path(__file__).resolve(),
    DIAGNOSIS_MODULE,
    REPOSITORY / "scripts/run_fgc_gr0_calibration_v7.py",
)
Q = Fraction


def _rel(path: Path) -> str:
    return path.resolve().relative_to(REPOSITORY).as_posix()


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _canonical(value: Mapping[str, Any]) -> str:
    return json.dumps(
        value,
        sort_keys=True,
        indent=2,
        ensure_ascii=True,
        allow_nan=False,
    ) + "\n"


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


def _load_json(path: Path, artifact_id: str | None = None) -> dict[str, Any]:
    raw = path.read_bytes()
    data = json.loads(raw)
    if not isinstance(data, dict):
        raise ValueError(f"{_rel(path)} must contain one JSON object")
    if raw != _canonical(data).encode("utf-8"):
        raise ValueError(f"{_rel(path)} is not canonical JSON")
    if artifact_id is not None and data.get("artifact_id") != artifact_id:
        raise ValueError(f"{_rel(path)} artifact identifier differs")
    return data


def _load_events(path: Path) -> list[dict[str, Any]]:
    raw = path.read_bytes()
    if not raw.endswith(b"\n"):
        raise ValueError("campaign event log must end in a newline")
    records: list[dict[str, Any]] = []
    rebuilt = bytearray()
    for line in raw.splitlines():
        value = json.loads(line)
        if not isinstance(value, dict):
            raise ValueError("campaign event log entries must be JSON objects")
        records.append(value)
        rebuilt.extend(_canonical_line(value))
    if bytes(rebuilt) != raw:
        raise ValueError("campaign event log is not canonical JSONL")
    return records


def _atomic_write(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    try:
        with open(descriptor, "wb", closefd=True) as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        Path(temporary).replace(path)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise


def _atomic_npz(path: Path, arrays: Mapping[str, np.ndarray]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    try:
        with os.fdopen(descriptor, "wb") as handle:
            np.savez_compressed(handle, **arrays)
            handle.flush()
            os.fsync(handle.fileno())
        Path(temporary).replace(path)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    with path.open("rb") as handle:
        config = tomllib.load(handle)
    if config.get("artifact_id") != ARTIFACT_ID or config.get("schema_version") != 1:
        raise ValueError("CAL4 config identity differs")
    scope = config.get("scope", {})
    immutable = config.get("immutable_campaign", {})
    ownership = config.get("constraint_ownership_diagnosis", {})
    spectral = config.get("spectral_ownership_diagnosis", {})
    prospective = config.get("prospective_repair", {})
    claims = config.get("claims", {})
    if (
        scope.get("target_protocol") != "FGC-2-SF1-PROTO7"
        or scope.get("GR0_calibration_trajectory_read") is not True
        or scope.get("SGBL_trajectory_read") is not False
        or scope.get("FGCQR_trajectory_read") is not False
        or scope.get("first_nonzero_common_event_completed") is not True
        or scope.get("mechanism_question_answered") is not False
        or immutable.get("authorization_commit")
        != "b8f1921985191b9ad488ecfd2f183c19936996b1"
        or immutable.get("amplitudes_in_order") != ["5/2", "3"]
        or immutable.get("event_log_line_count") != 8
        or immutable.get("first_nonzero_common_event_time") != "1/16"
        or ownership.get("fixed_outer_rows") != 4
        or ownership.get("roundoff_operation_budget_per_diagnostic_evaluation")
        != 4096
        or ownership.get("accumulated_roundoff_changes_order_classification_only")
        is not True
        or spectral.get("finest_pair_individual_budgets_required") is not True
        or spectral.get("all_adjacent_nested_tail_ratios_required") is not True
        or spectral.get("no_spectral_threshold_changes") is not True
        or prospective.get("new_protocol_required") != "FGC-2-SF1-PROTO8"
        or prospective.get("physical_inputs_unchanged") is not True
        or prospective.get("constraint_magnitude_guards_unchanged") is not True
        or prospective.get("spectral_thresholds_unchanged") is not True
        or claims.get("PROTO7_common_event_ownership_contract_obstructed")
        is not True
        or claims.get("FGCQR_holdout_execution_authorized") is not False
        or claims.get("FGCQR_mechanism_rejected") is not False
        or claims.get("general_gradient_route_rejected") is not False
    ):
        raise ValueError("CAL4 config violates its frozen semantic boundary")
    for key in (
        "runtime_authorization_result",
        "run_plan_config",
        "campaign_manifest",
        "campaign_event_log",
        "campaign_checkpoint",
        "campaign_result",
        "deterministic_replay_cache",
        "deterministic_replay_state_cache",
    ):
        if not isinstance(config.get(key), str):
            raise ValueError(f"CAL4 config path is absent: {key}")
    return config


def _campaign_paths(config: Mapping[str, Any]) -> dict[str, Path]:
    return {
        key: REPOSITORY / str(config[key])
        for key in (
            "runtime_authorization_result",
            "run_plan_config",
            "campaign_manifest",
            "campaign_event_log",
            "campaign_checkpoint",
            "campaign_result",
            "deterministic_replay_cache",
            "deterministic_replay_state_cache",
        )
    }


def _validate_campaign(
    config: Mapping[str, Any],
) -> tuple[
    dict[str, Path],
    dict[str, Any],
    list[dict[str, Any]],
    dict[str, Any],
    dict[str, Any],
]:
    paths = _campaign_paths(config)
    for key in (
        "runtime_authorization_result",
        "run_plan_config",
        "campaign_manifest",
        "campaign_event_log",
        "campaign_checkpoint",
        "campaign_result",
    ):
        if not paths[key].is_file():
            raise FileNotFoundError(paths[key])
    immutable = config["immutable_campaign"]
    hash_fields = {
        "campaign_manifest": "manifest_sha256",
        "campaign_event_log": "event_log_sha256",
        "campaign_checkpoint": "checkpoint_sha256",
        "campaign_result": "campaign_result_sha256",
    }
    for path_key, config_key in hash_fields.items():
        if _sha(paths[path_key]) != immutable[config_key]:
            raise ValueError(f"immutable CAL4 campaign hash differs: {path_key}")

    manifest = _load_json(paths["campaign_manifest"])
    result = _load_json(paths["campaign_result"])
    events = _load_events(paths["campaign_event_log"])
    authorization = _load_json(paths["runtime_authorization_result"], "FGC-1-HLT5-MON5")
    if (
        manifest.get("campaign_id") != immutable["campaign_id"]
        or manifest.get("git_commit") != immutable["authorization_commit"]
        or manifest.get("runner_id") != "FGC-1-CAL4-RUN1-RUNNER"
        or manifest.get("plan_sha256") != _sha(paths["run_plan_config"])
        or manifest.get("authorization_sha256")
        != _sha(paths["runtime_authorization_result"])
        or result.get("campaign_id") != immutable["campaign_id"]
        or result.get("classification") != immutable["terminal_classification"]
        or result.get("holdout_execution_authorized") is not False
        or result.get("mechanism_question_answered") is not False
        or result.get("FGCQR_outcome_read") is not False
        or result.get("SGBL_outcome_read") is not False
        or len(events) != immutable["event_log_line_count"]
    ):
        raise ValueError("immutable PROTO7 campaign identity or boundary differs")
    expected_sequence = [
        ("5/2", "common_event", 0),
        ("5/2", "rejected_unaccepted_source_only_proposal", None),
        ("5/2", "rejected_unaccepted_source_only_proposal", None),
        ("5/2", "common_event", 1),
        ("3", "common_event", 0),
        ("3", "rejected_unaccepted_source_only_proposal", None),
        ("3", "rejected_unaccepted_source_only_proposal", None),
        ("3", "common_event", 1),
    ]
    observed = [
        (item.get("amplitude"), item.get("event_type"), item.get("event_index"))
        for item in events
    ]
    if observed != expected_sequence:
        raise ValueError("PROTO7 campaign event sequence differs")
    for amplitude in ("5/2", "3"):
        t0 = next(
            item
            for item in events
            if item.get("amplitude") == amplitude and item.get("event_index") == 0
        )
        t1 = next(
            item
            for item in events
            if item.get("amplitude") == amplitude and item.get("event_index") == 1
        )
        retries = [
            item
            for item in events
            if item.get("amplitude") == amplitude
            and item.get("event_type")
            == "rejected_unaccepted_source_only_proposal"
        ]
        if (
            t0["assessment"]["primary_common_event"]["admission_passed"] is not True
            or t0["assessment"]["comparator_common_event"]["admission_passed"]
            is not True
            or t1["assessment"]["coordinate_time"] != 0.0625
            or len(retries) != 2
            or any(item.get("member") != "RK4-4097" for item in retries)
            or any(item.get("fields_bitwise_preserved") is not True for item in retries)
            or any(item.get("time_advanced") is not False for item in retries)
            or any(item.get("accepted_stage_count_advanced") is not False for item in retries)
            or any(item.get("causal_debit_advanced") is not False for item in retries)
            or any(item.get("external_transaction_state_preserved") is not True for item in retries)
            or any(item.get("retryable_source_only") is not True for item in retries)
            or any(item.get("non_source_failure_vetoed_retry") is not False for item in retries)
            or any(item.get("complete_failure_set") != ["newton_residual_limit"] for item in retries)
        ):
            raise ValueError("PROTO7 common-event or retry evidence differs")
    records = result.get("amplitude_records", [])
    if [item.get("amplitude") for item in records] != ["5/2", "3"]:
        raise ValueError("PROTO7 amplitude ordering differs")
    for item in records:
        stop = item.get("stop", {})
        if (
            item.get("last_completed_common_event_index") != 1
            or item.get("last_completed_coordinate_time") != 0.0625
            or stop.get("classification") != immutable["common_stop_classification"]
            or stop.get("reason") != immutable["common_stop_reason"]
        ):
            raise ValueError("PROTO7 terminal common-event record differs")
    return paths, manifest, events, result, authorization


def _spectral_thresholds(config: Mapping[str, Any]) -> SpectralThresholds:
    values = config["spectral_ownership_diagnosis"]
    return SpectralThresholds(
        unresolved_top_fraction=float(Q(values["unresolved_top_fraction"])),
        maximum_top_band_field_power_fraction=float(
            Q(values["maximum_top_band_field_power_fraction"])
        ),
        maximum_top_band_derivative_power_fraction=float(
            Q(values["maximum_top_band_derivative_power_fraction"])
        ),
        maximum_rms_scale_over_cutoff=float(
            Q(values["maximum_rms_scale_over_cutoff"])
        ),
        maximum_nested_tail_ratio=float(Q(values["maximum_nested_tail_ratio"])),
    )


def _method_diagnosis(
    config: Mapping[str, Any],
    plan: Mapping[str, Any],
    members: Mapping[str, runner.Proto7RunMember],
    *,
    method: str,
) -> dict[str, Any]:
    method_config = plan["primary_method" if method == "RK4" else "comparator_method"]
    order = int(method_config["spatial_order"])
    counts = tuple(int(item) for item in plan["numerics"]["resolutions"])
    selected = tuple(members[f"{method}-{count}"] for count in counts)
    snapshots = tuple(
        owned_semidiscrete_constraint_snapshot(
            member.state,
            member.initial.grid,
            coordinate_time=member.time,
            diagnostic_spatial_order=order,
            fixed_outer_rows=int(config["constraint_ownership_diagnosis"]["fixed_outer_rows"]),
            accepted_stage_count=member.transaction.state.accepted_stage_count,
            planck_mass=float(Q(plan["physical_inputs"]["planck_mass"])),
            scalar_mass=float(Q(plan["physical_inputs"]["scalar_mass"])),
            quartic_coupling=float(Q(plan["physical_inputs"]["quartic_coupling"])),
            length_unit=float(Q(plan["physical_inputs"]["length_unit_L0"])),
        )
        for member in selected
    )
    constraint = owned_constraint_admission(
        snapshots,
        method=method,
        coarsest_guard_maximum=float(Q(method_config["coarsest_constraint_guard"])),
        finest_guard_maximum=float(Q(method_config["finest_constraint_guard"])),
        minimum_finest_pair_order=float(
            Q(config["constraint_ownership_diagnosis"]["minimum_finest_pair_order"])
        ),
    )
    budgets = []
    budget_evidence: dict[str, dict[str, Any]] = {}
    threshold = _spectral_thresholds(config)
    cutoff = float(Q(plan["physical_inputs"]["cutoff_Lambda"]))
    radius_max = float(Q(plan["physical_inputs"]["measurement_radius_maximum"]))
    taper_fraction = float(Q(plan["numerics"]["proper_spectral_taper_fraction"]))
    for member in selected:
        profile = proper_radial_profile(
            member.state.u,
            member.state.q,
            member.initial.grid.coordinates,
            cutoff=cutoff,
            measurement_radius_maximum=radius_max,
        )
        proper = profile.proper_coordinates
        window = compact_vacuum_buffer_window(
            proper,
            support_minimum=float(proper[0]),
            support_maximum=float(proper[-1]),
            taper_width=float((proper[-1] - proper[0]) * taper_fraction),
        )
        record = spectral_field_budgets(
            profile.deviations,
            proper,
            cutoff=cutoff,
            window=window,
            thresholds=threshold,
        )
        budgets.append(record)
        budget_evidence[str(member.point_count)] = {
            "maximum_round_trip_interpolation_infinity": float(
                np.max(profile.round_trip_interpolation_infinity, initial=0.0)
            ),
            "fields": {name: asdict(value) for name, value in record.items()},
        }
    spectral = finest_pair_nested_spectral_admission(
        counts,
        budgets,
        thresholds=threshold,
    )

    reduction_names = PROTO4_CONSTRAINT_ORDER[4:]
    old_full_failures = []
    for name in reduction_names:
        full_values = constraint.full_domain_component_infinity[name]
        old_enclosure = 4096 * np.finfo(np.float64).eps
        effective = tuple(max(value - old_enclosure, 0.0) for value in full_values)
        medium, fine = effective[-2:]
        order_value = (
            None
            if medium == 0.0 or fine == 0.0
            else float(np.log2(medium / fine))
        )
        if order_value is not None and order_value < 1.5:
            old_full_failures.append(name)
    boundary_localized = all(
        all(
            radius > last_owned
            for radius, last_owned in zip(
                constraint.full_component_maximum_radii[name],
                constraint.last_owned_radii,
                strict=True,
            )
        )
        for name in old_full_failures
    )
    reduction_owned_enclosed = all(
        all(
            value <= enclosure
            for value, enclosure in zip(
                constraint.owned_domain_component_infinity[name],
                constraint.accumulated_roundoff_enclosures,
                strict=True,
            )
        )
        for name in reduction_names
    )
    return {
        "method": method,
        "member_state_sha256": {
            member.key: array_content_sha256(
                member.state.u, member.state.p, member.state.q
            )
            for member in selected
        },
        "constraint": asdict(constraint),
        "old_full_domain_reduction_order_failures": old_full_failures,
        "old_full_domain_reduction_failures_boundary_localized": boundary_localized,
        "owned_reduction_residuals_within_accumulated_roundoff_enclosure": reduction_owned_enclosed,
        "spectral": asdict(spectral),
        "spectral_raw_budgets": budget_evidence,
        "prospective_common_event_admission_passed": (
            constraint.admission_passed and spectral.admission_passed
        ),
    }


def _diagnose_members(
    config: Mapping[str, Any],
    plan: Mapping[str, Any],
    members: Mapping[str, runner.Proto7RunMember],
) -> dict[str, Any]:
    methods = {
        name: _method_diagnosis(config, plan, members, method=name)
        for name in ("RK4", "SSPRK3")
    }
    return {
        "methods": methods,
        "both_methods_prospectively_admitted": all(
            value["prospective_common_event_admission_passed"]
            for value in methods.values()
        ),
    }


def _restore_amplitude_3(
    config: Mapping[str, Any],
    paths: Mapping[str, Path],
    events: Sequence[Mapping[str, Any]],
    campaign_result: Mapping[str, Any],
) -> dict[str, Any]:
    plan, authorization = runner._validate_authorization(
        paths["run_plan_config"],
        paths["runtime_authorization_result"],
        require_fresh_namespace=False,
    )
    metadata, members, checkpoint_log = runner._restore_checkpoint(
        paths["campaign_checkpoint"],
        plan=plan,
        authorization=authorization,
    )
    if (
        metadata.get("terminal") is not True
        or metadata.get("amplitude") != "3"
        or metadata.get("completed_common_event_index") != 1
        or checkpoint_log != paths["campaign_event_log"].read_bytes()
    ):
        raise ValueError("terminal checkpoint does not own the amplitude-3 t1 slice")
    public = next(
        item["assessment"]
        for item in events
        if item.get("amplitude") == "3" and item.get("event_index") == 1
    )
    reproduced = inherited._serial(inherited._event_assessment(plan, members, 0.0625))
    if reproduced != public:
        raise ValueError("amplitude-3 checkpoint assessment differs from its event")
    expected_hashes = campaign_result["amplitude_records"][1]["member_final_state_hashes"]
    actual_hashes = {
        key: array_content_sha256(member.state.u, member.state.p, member.state.q)
        for key, member in members.items()
    }
    if actual_hashes != expected_hashes:
        raise ValueError("amplitude-3 checkpoint member hashes differ")
    diagnosis = _diagnose_members(config, plan, members)
    diagnosis["public_event_assessment_reproduced"] = True
    diagnosis["checkpoint_member_state_sha256"] = actual_hashes
    return diagnosis


def replay_amplitude_5over2(
    config_path: Path = DEFAULT_CONFIG,
    *,
    write_cache: bool = False,
) -> dict[str, Any]:
    config = load_config(config_path)
    paths, manifest, events, campaign_result, _authorization_record = _validate_campaign(
        config
    )
    plan, authorization = runner._validate_authorization(
        paths["run_plan_config"],
        paths["runtime_authorization_result"],
        require_fresh_namespace=False,
    )
    members = runner._build_members(plan, authorization, "5/2")
    retry_records: list[dict[str, Any]] = []
    target = float(Q(config["immutable_campaign"]["first_nonzero_common_event_time"]))
    for key in members:
        members[key].advance_to(
            target,
            retry_factor=float(Q(plan["numerics"]["retry_factor"])),
            maximum_CFL_retries=int(plan["numerics"]["maximum_CFL_retries_per_step"]),
            maximum_source_retries=int(
                plan["numerics"]["maximum_source_retries_per_step"]
            ),
            minimum_step_size=float(Q(plan["numerics"]["minimum_step_size"])),
            rejected_trial_sink=lambda record: retry_records.append(dict(record)),
        )
        members[key].append_common_event()
    reproduced = inherited._serial(inherited._event_assessment(plan, members, target))
    public = next(
        item["assessment"]
        for item in events
        if item.get("amplitude") == "5/2" and item.get("event_index") == 1
    )
    public_retries = [
        item
        for item in events
        if item.get("amplitude") == "5/2"
        and item.get("event_type") == "rejected_unaccepted_source_only_proposal"
    ]
    expected_hashes = campaign_result["amplitude_records"][0]["member_final_state_hashes"]
    actual_hashes = {
        key: array_content_sha256(member.state.u, member.state.p, member.state.q)
        for key, member in members.items()
    }
    assessment_matches = reproduced == public
    serialized_retries = inherited._serial(retry_records)
    retries_match = serialized_retries == public_retries
    hashes_match = actual_hashes == expected_hashes
    state_arrays = {
        f"{key.replace('-', '_')}_{field}": np.asarray(getattr(member.state, field))
        for key, member in members.items()
        for field in ("u", "p", "q")
    }
    if write_cache:
        _atomic_npz(paths["deterministic_replay_state_cache"], state_arrays)
    state_cache_sha256 = (
        _sha(paths["deterministic_replay_state_cache"])
        if paths["deterministic_replay_state_cache"].is_file()
        else None
    )
    payload = {
        "schema_version": 1,
        "cache_role": "deterministic_CAL4_amplitude_5over2_replay",
        "campaign_id": manifest["campaign_id"],
        "campaign_event_log_sha256": _sha(paths["campaign_event_log"]),
        "authorization_commit": manifest["git_commit"],
        "amplitude": "5/2",
        "coordinate_time": target,
        "public_event_assessment_reproduced": assessment_matches,
        "public_retry_records_reproduced": retries_match,
        "public_member_state_hashes_reproduced": hashes_match,
        "public_retry_records_sha256": hashlib.sha256(
            _canonical({"records": public_retries}).encode("utf-8")
        ).hexdigest(),
        "replayed_retry_records_sha256": hashlib.sha256(
            _canonical({"records": serialized_retries}).encode("utf-8")
        ).hexdigest(),
        "state_cache_sha256": state_cache_sha256,
        "member_state_sha256": actual_hashes,
        "source_retry_count": len(retry_records),
        "diagnosis": _diagnose_members(config, plan, members),
    }
    if write_cache:
        _atomic_write(paths["deterministic_replay_cache"], _canonical(payload).encode())
    if not assessment_matches:
        raise ValueError("amplitude-5/2 deterministic replay differs from its event")
    if not retries_match:
        raise ValueError("amplitude-5/2 deterministic retries differ from event log")
    if not hashes_match:
        raise ValueError("amplitude-5/2 deterministic state hashes differ")
    return payload


def _load_replay_cache(
    config: Mapping[str, Any],
    paths: Mapping[str, Path],
    events: Sequence[Mapping[str, Any]],
    campaign_result: Mapping[str, Any],
) -> dict[str, Any]:
    path = paths["deterministic_replay_cache"]
    if not path.is_file():
        raise FileNotFoundError(
            f"{path} is absent; run {Path(__file__).name} --replay-only first"
        )
    cache = _load_json(path)
    state_path = paths["deterministic_replay_state_cache"]
    if not state_path.is_file():
        raise FileNotFoundError(state_path)
    expected_hashes = campaign_result["amplitude_records"][0]["member_final_state_hashes"]
    public_retry_count = sum(
        item.get("amplitude") == "5/2"
        and item.get("event_type") == "rejected_unaccepted_source_only_proposal"
        for item in events
    )
    if (
        cache.get("cache_role") != "deterministic_CAL4_amplitude_5over2_replay"
        or cache.get("campaign_id") != config["immutable_campaign"]["campaign_id"]
        or cache.get("campaign_event_log_sha256")
        != config["immutable_campaign"]["event_log_sha256"]
        or cache.get("authorization_commit")
        != config["immutable_campaign"]["authorization_commit"]
        or cache.get("amplitude") != "5/2"
        or cache.get("coordinate_time") != 0.0625
        or cache.get("public_event_assessment_reproduced") is not True
        or cache.get("public_retry_records_reproduced") is not True
        or cache.get("public_member_state_hashes_reproduced") is not True
        or cache.get("member_state_sha256") != expected_hashes
        or cache.get("source_retry_count") != public_retry_count
        or cache.get("public_retry_records_sha256")
        != cache.get("replayed_retry_records_sha256")
        or cache.get("state_cache_sha256") != _sha(state_path)
    ):
        raise ValueError("CAL4 deterministic replay cache differs")
    with np.load(state_path, allow_pickle=False) as archive:
        for member, expected in expected_hashes.items():
            prefix = member.replace("-", "_")
            observed = array_content_sha256(
                archive[f"{prefix}_u"],
                archive[f"{prefix}_p"],
                archive[f"{prefix}_q"],
            )
            if observed != expected:
                raise ValueError("CAL4 deterministic replay state cache differs")
    return cache


def record(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    config = load_config(config_path)
    paths, manifest, events, campaign_result, authorization = _validate_campaign(config)
    replay = _load_replay_cache(config, paths, events, campaign_result)
    amplitude_3 = _restore_amplitude_3(
        config, paths, events, campaign_result
    )
    amplitude_5 = replay["diagnosis"]
    if amplitude_5.get("both_methods_prospectively_admitted") is not True:
        raise ValueError("amplitude-5/2 does not clear the prospective CAL4 repair")
    if amplitude_3.get("both_methods_prospectively_admitted") is not False:
        raise ValueError("prospective CAL4 repair unexpectedly rescues amplitude 3")
    for diagnosis in (amplitude_5, amplitude_3):
        for method in ("RK4", "SSPRK3"):
            item = diagnosis["methods"][method]
            if (
                item["old_full_domain_reduction_order_failures"]
                != ["reduction_alpha", "reduction_shift", "reduction_lambda", "reduction_R"]
                or item[
                    "old_full_domain_reduction_failures_boundary_localized"
                ]
                is not True
                or item[
                    "owned_reduction_residuals_within_accumulated_roundoff_enclosure"
                ]
                is not True
            ):
                raise ValueError("CAL4 reduction-ownership diagnosis differs")

    immutable = config["immutable_campaign"]
    payload = {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": config["project_version"],
        "metric_signature": config["metric_signature"],
        "riemann_convention": config["riemann_convention"],
        "artifact_payload": {
            "immutable_campaign": {
                "authorization_commit": immutable["authorization_commit"],
                "campaign_id": immutable["campaign_id"],
                "manifest_sha256": immutable["manifest_sha256"],
                "event_log_sha256": immutable["event_log_sha256"],
                "checkpoint_sha256": immutable["checkpoint_sha256"],
                "campaign_result_sha256": immutable["campaign_result_sha256"],
                "event_log_line_count": len(events),
                "elapsed_wall_seconds": campaign_result["elapsed_wall_seconds"],
                "terminal_classification": campaign_result["classification"],
                "first_nonzero_common_event_completed_by_both_amplitudes": True,
                "amplitude_stop_records": campaign_result["amplitude_records"],
                "PROTO7_source_only_retry_counts": {
                    item["amplitude"]: item["member_source_retry_counts"]
                    for item in campaign_result["amplitude_records"]
                },
                "manifest_implementation_sha256": manifest["implementation_sha256"],
                "HLT5_gate_status": authorization["gate_status"],
            },
            "common_event_ownership_diagnosis": {
                "amplitude_5over2": amplitude_5,
                "amplitude_3": amplitude_3,
                "deterministic_replay_cache_sha256": _sha(
                    paths["deterministic_replay_cache"]
                ),
                "deterministic_replay_state_cache_sha256": _sha(
                    paths["deterministic_replay_state_cache"]
                ),
                "unchanged_thresholds": {
                    "constraint_minimum_finest_pair_order": "3/2",
                    "spectral_maximum_top_band_field_power_fraction": "1/1048576",
                    "spectral_maximum_top_band_derivative_power_fraction": "1/1024",
                    "spectral_maximum_rms_scale_over_cutoff": "1/4",
                    "spectral_maximum_nested_tail_ratio": "1/4",
                },
            },
            "decision": {
                "PROTO7_common_event_ownership_contract_obstructed": True,
                "new_protocol_required": "FGC-2-SF1-PROTO8",
                "reason": (
                    "the_old_gate_mixes_projector_owned_outer_reduction_error_and_"
                    "coarse_grid_spectral_truncation_into_a_veto_even_though_the_"
                    "evolution_owned_constraints_and_finest_pair_spectra_pass_for_"
                    "amplitude_5over2_under_every_unchanged_threshold"
                ),
                "prospective_repair": dict(config["prospective_repair"]),
            },
            "epistemic_boundary": {
                "PROTO7_transaction_repair_succeeded_through_t1": True,
                "GR0_amplitude_selected": False,
                "GR0_calibration_completed": False,
                "SGBL_or_FGCQR_trajectory_read": False,
                "candidate_or_general_gradient_route_rejected": False,
                "dynamic_trapped_sphere_classified": False,
                "mechanism_question_answered": False,
            },
        },
        "gate_status": dict(config["claims"]),
        "implementation_sha256": {_rel(path): _sha(path) for path in IMPLEMENTATION},
        "derivation_document": _rel(OWNER_DOCUMENT),
        "derivation_document_sha256": _sha(OWNER_DOCUMENT),
        "source_config_sha256": {_rel(config_path): _sha(config_path)},
        "nonclaims": {
            "fresh_GR0_calibration_completed": False,
            "trapped_sphere_observed": False,
            "SGBL_control_completed": False,
            "FGCQR_holdout_evolved": False,
            "regulator_activation_observed": False,
            "positive_Raychaudhuri_margin_observed": False,
            "FGCQR_mechanism_rejected": False,
            "general_gradient_route_rejected": False,
            "retained_EFT_validity": False,
            "physical_transition": False,
        },
    }
    return payload


def _validate_stored(result: Mapping[str, Any], config: Mapping[str, Any]) -> None:
    paths, _manifest, events, campaign_result, _authorization = _validate_campaign(
        config
    )
    replay = _load_replay_cache(config, paths, events, campaign_result)
    artifact = result.get("artifact_payload", {})
    diagnosis = artifact.get("common_event_ownership_diagnosis", {})
    decision = artifact.get("decision", {})
    boundary = artifact.get("epistemic_boundary", {})
    if (
        result.get("artifact_id") != ARTIFACT_ID
        or diagnosis.get("deterministic_replay_cache_sha256")
        != _sha(paths["deterministic_replay_cache"])
        or diagnosis.get("deterministic_replay_state_cache_sha256")
        != _sha(paths["deterministic_replay_state_cache"])
        or diagnosis.get("amplitude_5over2") != replay.get("diagnosis")
        or diagnosis.get("amplitude_3", {}).get(
            "both_methods_prospectively_admitted"
        )
        is not False
        or diagnosis.get("unchanged_thresholds")
        != {
            "constraint_minimum_finest_pair_order": "3/2",
            "spectral_maximum_top_band_field_power_fraction": "1/1048576",
            "spectral_maximum_top_band_derivative_power_fraction": "1/1024",
            "spectral_maximum_rms_scale_over_cutoff": "1/4",
            "spectral_maximum_nested_tail_ratio": "1/4",
        }
        or decision.get("PROTO7_common_event_ownership_contract_obstructed")
        is not True
        or decision.get("new_protocol_required") != "FGC-2-SF1-PROTO8"
        or decision.get("prospective_repair") != config["prospective_repair"]
        or boundary.get("GR0_calibration_completed") is not False
        or boundary.get("SGBL_or_FGCQR_trajectory_read") is not False
        or boundary.get("candidate_or_general_gradient_route_rejected") is not False
        or boundary.get("mechanism_question_answered") is not False
        or result.get("gate_status") != config["claims"]
        or not result.get("nonclaims")
        or any(value is not False for value in result.get("nonclaims", {}).values())
    ):
        raise ValueError("stored CAL4 result violates its semantic boundary")
    if (
        result.get("implementation_sha256")
        != {_rel(path): _sha(path) for path in IMPLEMENTATION}
        or result.get("derivation_document") != _rel(OWNER_DOCUMENT)
        or result.get("derivation_document_sha256") != _sha(OWNER_DOCUMENT)
        or result.get("source_config_sha256")
        != {_rel(DEFAULT_CONFIG): _sha(DEFAULT_CONFIG)}
    ):
        raise ValueError("stored CAL4 implementation or derivation hash differs")


def load_canonical_result(path: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    return _load_json(path, ARTIFACT_ID)


def verify_canonical(
    path: Path = DEFAULT_OUTPUT,
    config_path: Path = DEFAULT_CONFIG,
    *,
    replay: bool = False,
) -> None:
    config = load_config(config_path)
    actual = load_canonical_result(path)
    _validate_stored(actual, config)
    if replay:
        reconstructed = replay_amplitude_5over2(config_path, write_cache=False)
        stored = actual["artifact_payload"]["common_event_ownership_diagnosis"]
        if reconstructed["diagnosis"] != stored["amplitude_5over2"]:
            raise ValueError("stored CAL4 result differs from deterministic replay")
        if actual != record(config_path):
            raise ValueError("stored CAL4 result differs from full reconstruction")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--replay", action="store_true")
    parser.add_argument("--replay-only", action="store_true")
    args = parser.parse_args()
    config = args.config.resolve()
    output = args.output.resolve()
    if args.replay_only:
        payload = replay_amplitude_5over2(config, write_cache=True)
        cache = REPOSITORY / load_config(config)["deterministic_replay_cache"]
        print(
            f"wrote {_rel(cache)} "
            f"(prospective_pass={str(payload['diagnosis']['both_methods_prospectively_admitted']).lower()})"
        )
        return
    if args.check:
        verify_canonical(output, config, replay=args.replay)
        print(f"verified {_rel(output)} (full_replay={str(args.replay).lower()})")
        return
    if args.replay or not (
        REPOSITORY / load_config(config)["deterministic_replay_cache"]
    ).is_file():
        replay_amplitude_5over2(config, write_cache=True)
    payload = record(config)
    _atomic_write(output, _canonical(payload).encode("utf-8"))
    print(
        f"wrote {_rel(output)} "
        "(PROTO7 common-event ownership obstructed=true; mechanism answered=false)"
    )


if __name__ == "__main__":
    main()
