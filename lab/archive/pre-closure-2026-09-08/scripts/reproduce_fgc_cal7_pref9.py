#!/usr/bin/env python3
"""Reproduce CAL7/PREF9's PROTO10 centre-roundoff diagnosis."""

from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import tomllib
from typing import Any, Mapping

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
import sys

sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.cal7_center_roundoff_diagnosis import (  # noqa: E402
    EXPECTED_AMPLITUDES,
    EXPECTED_POINT_COUNTS,
    diagnose_campaign_pair,
    initial_source_floor_record,
    well_balanced_projected_state,
)
from recursive_horizons.fgc.evolution.calibration_runtime import (  # noqa: E402
    project_gr0_semidiscrete_state,
)
from recursive_horizons.fgc.evolution.gr0_calibration import (  # noqa: E402
    construct_gr0_grid_initial_data,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    array_content_sha256,
)
from recursive_horizons.fgc.initial_data_preflight import (  # noqa: E402
    PulseParameters,
)


ARTIFACT_ID = "FGC-1-CAL7-PREF9"
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-cal7-pref9.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-cal7-pref9.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-cal7-pref9.md"
DIAGNOSIS_MODULE = (
    REPOSITORY
    / "src/recursive_horizons/fgc/evolution/cal7_center_roundoff_diagnosis.py"
)
IMPLEMENTATION = (Path(__file__).resolve(), DIAGNOSIS_MODULE)
Q = Fraction


EXPECTED_CLAIMS = {
    "PROTO10_campaign_terminated_normally": True,
    "PROTO10_guarded_common_event_contract_exercised": True,
    "PROTO10_GR0_case_eligible": False,
    "PROTO10_center_adjacent_binary64_source_floor_localized": True,
    "PROTO11_well_balanced_reference_map_required": True,
    "PROTO11_frozen": False,
    "fresh_GR0_dynamic_calibration_completed": False,
    "classical_spherical_diagnostic_authorized": False,
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
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _bytes_sha(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


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


def _load_json(path: Path, artifact_id: str | None = None) -> dict[str, Any]:
    raw = path.read_bytes()
    value = json.loads(raw)
    if not isinstance(value, dict) or raw != _canonical(value).encode("utf-8"):
        raise ValueError(f"{_rel(path)} must contain canonical JSON")
    if artifact_id is not None and value.get("artifact_id") != artifact_id:
        raise ValueError(f"{_rel(path)} artifact identifier differs")
    return value


def _load_events(path: Path) -> tuple[list[dict[str, Any]], bytes]:
    raw = path.read_bytes()
    if not raw.endswith(b"\n"):
        raise ValueError("campaign event log must end in a newline")
    records: list[dict[str, Any]] = []
    rebuilt = bytearray()
    for line in raw.splitlines():
        value = json.loads(line)
        if not isinstance(value, dict):
            raise ValueError("campaign event records must be JSON objects")
        records.append(value)
        rebuilt.extend(_canonical_line(value))
    if bytes(rebuilt) != raw:
        raise ValueError("campaign event log is not canonical JSONL")
    return records, raw


def _git(*args: str, check: bool = True) -> subprocess.CompletedProcess[bytes]:
    return subprocess.run(
        ["git", *args],
        cwd=REPOSITORY,
        check=check,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )


def _validate_commit_lineage(commit: str, implementation: Mapping[str, str]) -> None:
    _git("cat-file", "-e", f"{commit}^{{commit}}")
    if _git("merge-base", "--is-ancestor", commit, "HEAD", check=False).returncode != 0:
        raise ValueError("campaign authorization commit is not an ancestor of HEAD")
    for relative, expected in implementation.items():
        blob = _git("show", f"{commit}:{relative}").stdout
        if _bytes_sha(blob) != expected:
            raise ValueError(f"immutable campaign implementation differs: {relative}")


def _strict_keys(name: str, value: Mapping[str, Any], expected: set[str]) -> None:
    if set(value) != expected:
        raise ValueError(f"{name} keys differ")


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    with path.open("rb") as handle:
        config = tomllib.load(handle)
    _strict_keys(
        "CAL7 config",
        config,
        {
            "schema_version",
            "artifact_id",
            "project_version",
            "metric_signature",
            "riemann_convention",
            "runtime_authorization_result",
            "run_plan_config",
            "campaign_manifest",
            "campaign_event_log",
            "campaign_checkpoint",
            "campaign_result",
            "scope",
            "immutable_campaign",
            "center_roundoff_diagnosis",
            "prospective_revision",
            "proof_contract",
            "claims",
        },
    )
    scope = config["scope"]
    immutable = config["immutable_campaign"]
    diagnosis = config["center_roundoff_diagnosis"]
    revision = config["prospective_revision"]
    proof = config["proof_contract"]
    if (
        config["schema_version"] != 1
        or config["artifact_id"] != ARTIFACT_ID
        or scope.get("target_protocol") != "FGC-2-SF1-PROTO10"
        or scope.get("role")
        != "post_calibration_center_adjacent_binary64_source_floor_diagnosis"
        or scope.get("GR0_calibration_trajectory_read") is not True
        or scope.get("SGBL_trajectory_read") is not False
        or scope.get("FGCQR_trajectory_read") is not False
        or scope.get("mechanism_question_answered") is not False
        or immutable.get("terminal_classification")
        != "calibration_failed_no_eligible_GR0_case"
        or immutable.get("amplitudes_in_order") != list(EXPECTED_AMPLITUDES)
        or immutable.get("event_log_line_count") != 264
        or immutable.get("common_event_count") != 2
        or immutable.get("source_only_rejection_count") != 262
        or immutable.get("per_amplitude_source_only_rejection_count") != 131
        or immutable.get("terminal_member") != "RK4-8193"
        or immutable.get("terminal_reason") != "newton_residual_limit"
        or immutable.get("terminal_exhaustion_kind") != "minimum_step_size"
        or immutable.get("last_completed_common_event_index") != 0
        or immutable.get("target_common_event_index") != 1
        or immutable.get("accepted_step_index") != 15
        or immutable.get("cumulative_8193_source_retry_count") != 129
        or immutable.get("source_residual_infinity_max") != "1/1000000000000"
        or immutable.get("minimum_step_size") != "1/1073741824"
        or diagnosis.get("point_counts") != list(EXPECTED_POINT_COUNTS)
        or diagnosis.get("primary_method") != "RK4"
        or diagnosis.get("primary_spatial_order") != 4
        or diagnosis.get("maximum_grid_index") != 1
        or diagnosis.get("maximum_field") != "alpha"
        or diagnosis.get("maximum_radii") != ["1/16", "1/32", "1/64"]
        or diagnosis.get("expected_refinement_ratios") != ["4", "4"]
        or diagnosis.get("roundoff_coefficient_target") != "83/96"
        or diagnosis.get("well_balanced_reference_u")
        != ["1", "0", "1", "r", "0", "0"]
        or diagnosis.get("well_balanced_reference_q")
        != ["0", "0", "0", "1", "0", "0"]
        or diagnosis.get("well_balanced_control_is_not_a_runtime_result") is not True
        or revision.get("new_protocol_required") != "FGC-2-SF1-PROTO11"
        or any(value is not True for key, value in revision.items() if key != "new_protocol_required")
        or any(value is not True for value in proof.values())
        or config["claims"] != EXPECTED_CLAIMS
    ):
        raise ValueError("CAL7/PREF9 configuration contract differs")
    return config


def _validate_campaign(
    config: Mapping[str, Any],
) -> tuple[dict[str, Any], list[dict[str, Any]], bytes, dict[str, Any]]:
    immutable = config["immutable_campaign"]
    manifest_path = REPOSITORY / config["campaign_manifest"]
    event_path = REPOSITORY / config["campaign_event_log"]
    checkpoint_path = REPOSITORY / config["campaign_checkpoint"]
    result_path = REPOSITORY / config["campaign_result"]
    for path, expected in (
        (manifest_path, immutable["manifest_sha256"]),
        (event_path, immutable["event_log_sha256"]),
        (checkpoint_path, immutable["checkpoint_sha256"]),
        (result_path, immutable["campaign_result_sha256"]),
    ):
        if _sha(path) != expected:
            raise ValueError(f"immutable campaign hash differs: {_rel(path)}")
    manifest = _load_json(manifest_path)
    result = _load_json(result_path)
    events, event_bytes = _load_events(event_path)
    if (
        manifest.get("campaign_id") != immutable["campaign_id"]
        or manifest.get("git_commit") != immutable["authorization_commit"]
        or manifest.get("runner_id") != "FGC-1-CAL7-RUN1-RUNNER"
        or manifest.get("ordered_amplitudes") != list(EXPECTED_AMPLITUDES)
        or manifest.get("outcome_fields_present_at_creation") is not False
        or result.get("campaign_id") != immutable["campaign_id"]
        or result.get("runner_id") != "FGC-1-CAL7-RUN1-RUNNER"
        or result.get("classification") != immutable["terminal_classification"]
        or result.get("event_log_sha256") != immutable["event_log_sha256"]
        or result.get("selected_amplitude") is not None
        or result.get("holdout_execution_authorized") is not False
        or result.get("mechanism_question_answered") is not False
        or result.get("SGBL_outcome_read") is not False
        or result.get("FGCQR_outcome_read") is not False
        or len(events) != immutable["event_log_line_count"]
    ):
        raise ValueError("immutable campaign manifest/result contract differs")
    authorization_path = REPOSITORY / config["runtime_authorization_result"]
    plan_path = REPOSITORY / config["run_plan_config"]
    if (
        _sha(authorization_path) != manifest.get("authorization_sha256")
        or _sha(plan_path) != manifest.get("plan_sha256")
        or _rel(authorization_path) != manifest.get("authorization_path")
        or _rel(plan_path) != manifest.get("plan_path")
    ):
        raise ValueError("campaign authorization or plan binding differs")
    _validate_commit_lineage(
        immutable["authorization_commit"],
        manifest.get("implementation_sha256", {}),
    )
    with np.load(checkpoint_path, allow_pickle=False) as checkpoint:
        checkpoint_event_bytes = bytes(checkpoint["event_log_utf8"])
        metadata = json.loads(bytes(checkpoint["metadata_utf8"]).decode("utf-8"))
        if (
            checkpoint_event_bytes != event_bytes
            or metadata.get("terminal") is not True
            or metadata.get("campaign_id") != immutable["campaign_id"]
            or metadata.get("completed_common_event_index") != 0
            or metadata.get("event_log", {}).get("sha256") != immutable["event_log_sha256"]
            or metadata.get("event_log", {}).get("line_count") != len(events)
            or metadata.get("completed_amplitude_records") != result.get("amplitude_records")
        ):
            raise ValueError("terminal checkpoint metadata differs")
        final_record = result["amplitude_records"][-1]
        for member, expected_hash in final_record["member_final_state_hashes"].items():
            method, count = member.split("-")
            observed = array_content_sha256(
                checkpoint[f"{method}_{count}_u"],
                checkpoint[f"{method}_{count}_p"],
                checkpoint[f"{method}_{count}_q"],
            )
            if observed != expected_hash:
                raise ValueError(f"terminal checkpoint state differs: {member}")
    return manifest, events, event_bytes, result


def _initial_source_records(
    config: Mapping[str, Any],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], dict[str, Any]]:
    with (REPOSITORY / config["run_plan_config"]).open("rb") as handle:
        plan = tomllib.load(handle)
    authorization = _load_json(
        REPOSITORY / config["runtime_authorization_result"],
        "FGC-1-HLT8-MON8",
    )
    frozen_inputs = {
        (item["amplitude"], item["method"], item["point_count"]): item
        for item in authorization["artifact_payload"]["frozen_run_inputs"]
    }
    source_checks = {
        (item["amplitude"], item["method"], item["point_count"]): item
        for item in authorization["artifact_payload"]["accepted_state_source_prechecks"]
    }
    physical = plan["physical_inputs"]
    thresholds = plan["universal_thresholds"]
    raw_records: list[dict[str, Any]] = []
    balanced_records: list[dict[str, Any]] = []
    balanced_hashes: list[dict[str, Any]] = []
    for amplitude in EXPECTED_AMPLITUDES:
        parameters = PulseParameters(
            chi_amplitude=float(Q(amplitude)),
            center=float(Q(physical["chi_center"])),
            half_width=float(Q(physical["chi_half_width"])),
            phi_amplitude=float(Q(physical["phi_seed_amplitude"])),
            planck_mass=float(Q(physical["planck_mass"])),
            scalar_mass=float(Q(physical["scalar_mass"])),
            quartic_coupling=float(Q(physical["quartic_coupling"])),
        )
        for method in (plan["primary_method"], plan["comparator_method"]):
            label = method["method_label"]
            order = method["spatial_order"]
            for point_count in EXPECTED_POINT_COUNTS:
                initial = construct_gr0_grid_initial_data(
                    parameters,
                    point_count=point_count,
                    outer_radius=float(Q(physical["outer_radius"])),
                    constraint_method=method["constraint_solve_method"],
                    diagnostic_spatial_order=order,
                )
                raw_state = project_gr0_semidiscrete_state(initial, spatial_order=order)
                key = (amplitude, label, point_count)
                if array_content_sha256(raw_state.u, raw_state.p, raw_state.q) != frozen_inputs[key]["projected_state_sha256"]:
                    raise ValueError(f"reconstructed PROTO10 input differs: {key}")
                raw = initial_source_floor_record(
                    initial,
                    amplitude=amplitude,
                    method=label,
                    spatial_order=order,
                    representation="PROTO10_raw",
                    residual_tolerance=float(Q(thresholds["source_residual_infinity_max"])),
                    condition_number_maximum=float(Q(thresholds["kinetic_condition_number_max"])),
                )
                if raw.residual_infinity != source_checks[key]["source_residual_infinity"]:
                    raise ValueError(f"reconstructed PROTO10 source residual differs: {key}")
                balanced = initial_source_floor_record(
                    initial,
                    amplitude=amplitude,
                    method=label,
                    spatial_order=order,
                    representation="well_balanced_control",
                    residual_tolerance=float(Q(thresholds["source_residual_infinity_max"])),
                    condition_number_maximum=float(Q(thresholds["kinetic_condition_number_max"])),
                )
                balanced_state = well_balanced_projected_state(initial, spatial_order=order)
                raw_records.append(raw.as_record())
                balanced_records.append(balanced.as_record())
                balanced_hashes.append(
                    {
                        "amplitude": amplitude,
                        "method": label,
                        "point_count": point_count,
                        "projected_state_sha256": array_content_sha256(
                            balanced_state.u,
                            balanced_state.p,
                            balanced_state.q,
                        ),
                        "differs_from_PROTO10_input": array_content_sha256(
                            balanced_state.u,
                            balanced_state.p,
                            balanced_state.q,
                        )
                        != frozen_inputs[key]["projected_state_sha256"],
                    }
                )
    raw_primary = [
        item
        for item in raw_records
        if item["amplitude"] == "5/2" and item["method"] == "RK4"
    ]
    raw_other = [
        item
        for item in raw_records
        if item["amplitude"] == "3" and item["method"] == "RK4"
    ]
    if [item["point_count"] for item in raw_primary] != list(EXPECTED_POINT_COUNTS):
        raise ValueError("primary raw source-floor order differs")
    if any(
        item["maximum_grid_index"] != 1
        or item["maximum_field"] != "alpha"
        for item in raw_primary
    ):
        raise ValueError("primary raw source floor is not the first alpha row")
    raw_residuals = [item["residual_infinity"] for item in raw_primary]
    configured_raw = [float(value) for value in config["center_roundoff_diagnosis"]["raw_initial_residuals"]]
    if raw_residuals != configured_raw or [
        item["residual_infinity"] for item in raw_other
    ] != raw_residuals:
        raise ValueError("configured or cross-amplitude raw source floor differs")
    ratios = [raw_residuals[index + 1] / raw_residuals[index] for index in range(2)]
    coefficient_target = float(Q(config["center_roundoff_diagnosis"]["roundoff_coefficient_target"]))
    coefficient_errors = [
        abs(item["roundoff_coefficient"] - coefficient_target)
        for item in raw_primary
    ]
    tolerance = float(Q(thresholds["source_residual_infinity_max"]))
    if (
        any(abs(value - 4.0) > 5.0e-13 for value in ratios)
        or max(coefficient_errors) > 1.0e-12
        or any(item["residual_infinity"] >= tolerance for item in balanced_records)
        or not all(
            item["differs_from_PROTO10_input"]
            for item in balanced_hashes
            if item["method"] == "RK4"
        )
    ):
        raise ValueError("centre-roundoff or well-balanced control gate differs")
    comparison = {
        "primary_raw_refinement_ratios": ratios,
        "primary_roundoff_coefficient_target": coefficient_target,
        "primary_roundoff_coefficient_maximum_absolute_error": max(coefficient_errors),
        "binary64_epsilon": np.finfo(np.float64).eps,
        "all_well_balanced_source_residuals_pass_unchanged_raw_tolerance": True,
        "maximum_well_balanced_source_residual": max(
            item["residual_infinity"] for item in balanced_records
        ),
        "fine_primary_residual_reduction_factor": (
            raw_primary[-1]["residual_infinity"]
            / next(
                item["residual_infinity"]
                for item in balanced_records
                if item["amplitude"] == "5/2"
                and item["method"] == "RK4"
                and item["point_count"] == 8193
            )
        ),
        "well_balanced_control_changes_every_primary_projected_input_hash": True,
        "comparator_hash_may_remain_identical_when_its_stencil_is_already_well_balanced": True,
        "well_balanced_control_is_not_a_runtime_or_campaign_result": True,
    }
    return raw_records, balanced_records, {
        "comparison": comparison,
        "well_balanced_input_hashes": balanced_hashes,
    }


def reproduce(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    config = load_config(config_path)
    manifest, events, event_bytes, campaign_result = _validate_campaign(config)
    campaign_diagnosis = diagnose_campaign_pair(events, campaign_result)
    diagnosis_limits = config["center_roundoff_diagnosis"]
    if (
        campaign_diagnosis["terminal_time_absolute_difference"]
        > float(Q(diagnosis_limits["terminal_time_absolute_difference_max"]))
        or campaign_diagnosis["terminal_residual_relative_difference"]
        > float(Q(diagnosis_limits["terminal_residual_relative_difference_max"]))
        or campaign_diagnosis["paired_maximum_time_absolute_difference"]
        > float(Q(diagnosis_limits["paired_trace_time_absolute_difference_max"]))
    ):
        raise ValueError("paired PROTO10 obstruction is not amplitude-independent enough")
    raw_records, balanced_records, balanced_summary = _initial_source_records(config)
    runtime_authorization = REPOSITORY / config["runtime_authorization_result"]
    run_plan = REPOSITORY / config["run_plan_config"]
    record = {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": config["project_version"],
        "classification": (
            "post_calibration_PROTO10_center_adjacent_binary64_source_floor_diagnosis"
        ),
        "generated_by": "scripts/reproduce_fgc_cal7_pref9.py",
        "derivation_document": _rel(OWNER_DOCUMENT),
        "derivation_document_sha256": _sha(OWNER_DOCUMENT),
        "source_config_sha256": {_rel(config_path): _sha(config_path)},
        "predecessor_sha256": {
            _rel(runtime_authorization): _sha(runtime_authorization),
            _rel(run_plan): _sha(run_plan),
        },
        "implementation_sha256": {
            _rel(path): _sha(path) for path in IMPLEMENTATION
        },
        "scope_bindings": dict(config["scope"]),
        "gate_status": dict(EXPECTED_CLAIMS),
        "artifact_payload": {
            "immutable_campaign": {
                "authorization_commit": config["immutable_campaign"]["authorization_commit"],
                "campaign_id": manifest["campaign_id"],
                "manifest_sha256": config["immutable_campaign"]["manifest_sha256"],
                "event_log_sha256": _bytes_sha(event_bytes),
                "checkpoint_sha256": config["immutable_campaign"]["checkpoint_sha256"],
                "campaign_result_sha256": config["immutable_campaign"]["campaign_result_sha256"],
                "elapsed_wall_seconds": campaign_result["elapsed_wall_seconds"],
                "classification": campaign_result["classification"],
            },
            "campaign_diagnosis": campaign_diagnosis,
            "raw_initial_source_floor": raw_records,
            "well_balanced_initial_control": balanced_records,
            "well_balanced_control_summary": balanced_summary,
            "prospective_revision": dict(config["prospective_revision"]),
            "epistemic_boundary": {
                "observed_result": "repeatable centre-adjacent binary64 source-floor obstruction in the frozen PROTO10 GR-0 calibration",
                "localized_numerical_cause": "full-state SBP differentiation leaves a first-positive-row residual proportional to binary64_epsilon/r_squared",
                "not_observed": [
                    "an evolved common event after t=0",
                    "a trapped sphere",
                    "SGB-L or FGC-QR dynamics",
                    "regulator activation or metric-null defocusing",
                    "a continuum branch loss",
                    "a retained-EFT or physical transition",
                ],
                "repair_status": "well-balanced t0 control passes, but PROTO11 is not frozen or run",
            },
        },
        "nonclaims": {
            key: value
            for key, value in EXPECTED_CLAIMS.items()
            if value is False
        },
    }
    return record


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    record = reproduce(args.config.resolve())
    payload = _canonical(record).encode("utf-8")
    output = args.output.resolve()
    if args.check:
        if not output.is_file() or output.read_bytes() != payload:
            raise SystemExit(f"stored result differs: {output}")
        print(f"verified {output}")
        return
    _atomic_write(output, payload)
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
