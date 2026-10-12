#!/usr/bin/env python3
"""Reproduce RSP1/PREF12's checkpoint-recomputed post-run certificate."""

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
from typing import Any, Mapping, Sequence


REPOSITORY = Path(__file__).resolve().parents[1]
import sys

sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts import run_fgc_gr0_calibration as inherited  # noqa: E402
from scripts import run_fgc_rsp1_resolution_study as runner  # noqa: E402
from recursive_horizons.fgc.evolution.rsp1_campaign_diagnosis import (  # noqa: E402
    EXPECTED_CLASSIFICATION,
    EXPECTED_MEMBER_KEYS,
    EXPECTED_METHODS,
    EXPECTED_POINT_COUNTS,
    diagnose_rsp1_campaign,
)


ARTIFACT_ID = "FGC-1-RSP1-PREF12"
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-rsp1-pref12.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-rsp1-pref12.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-rsp1-pref12.md"
DIAGNOSIS_MODULE = (
    REPOSITORY
    / "src/recursive_horizons/fgc/evolution/rsp1_campaign_diagnosis.py"
)
IMPLEMENTATION = (
    Path(__file__).resolve(),
    DIAGNOSIS_MODULE,
    REPOSITORY / "scripts/run_fgc_rsp1_resolution_study.py",
    REPOSITORY
    / "src/recursive_horizons/fgc/evolution/rsp1_resolution_runtime.py",
)
Q = Fraction


EXPECTED_CLAIMS = {
    "RSP1_study_terminated_normally": True,
    "RSP1_all_six_members_reached_endpoint": True,
    "RSP1_terminal_checkpoint_recomputed": True,
    "RSP1_source_retry_total_zero": True,
    "amplitude_three_spectral_veto_cleared_for_successor_design": True,
    "RSP1_generic_all_field_raw_spectral_admission_passed": False,
    "PROTO12_design_may_begin": True,
    "PROTO12_frozen": False,
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
        raise ValueError("RSP1 event log must end in a newline")
    records: list[dict[str, Any]] = []
    rebuilt = bytearray()
    for line in raw.splitlines():
        value = json.loads(line)
        if not isinstance(value, dict):
            raise ValueError("RSP1 event records must be JSON objects")
        records.append(value)
        rebuilt.extend(_canonical_line(value))
    if bytes(rebuilt) != raw:
        raise ValueError("RSP1 event log is not canonical JSONL")
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
    if _git("merge-base", "--is-ancestor", commit, "HEAD", check=False).returncode:
        raise ValueError("RSP1 authorization commit is not an ancestor of HEAD")
    if not implementation:
        raise ValueError("RSP1 implementation lineage is absent")
    for relative, expected in implementation.items():
        blob = _git("show", f"{commit}:{relative}").stdout
        if _bytes_sha(blob) != expected:
            raise ValueError(f"immutable RSP1 implementation differs: {relative}")


def _strict_keys(name: str, value: Mapping[str, Any], expected: set[str]) -> None:
    if set(value) != expected:
        raise ValueError(f"{name} keys differ")


def _as_float(value: str) -> float:
    return float(Q(value))


def _expected_terminal_members(config: Mapping[str, Any]) -> dict[str, Any]:
    records = config["terminal_members"]
    expected_tables = {
        "RK4_4097",
        "RK4_8193",
        "RK4_16385",
        "SSPRK3_4097",
        "SSPRK3_8193",
        "SSPRK3_16385",
    }
    if set(records) != expected_tables:
        raise ValueError("RSP1 terminal-member tables differ")
    answer: dict[str, Any] = {}
    for table in sorted(records):
        item = records[table]
        _strict_keys(
            f"RSP1 terminal member {table}",
            item,
            {
                "key",
                "step_index",
                "transaction_serial",
                "accepted_stage_count",
                "CFL_retry_count",
                "source_retry_count",
                "state_sha256",
            },
        )
        key = item["key"]
        if key not in EXPECTED_MEMBER_KEYS or key in answer:
            raise ValueError("RSP1 terminal-member identity differs")
        if (
            item["step_index"] <= 0
            or item["transaction_serial"] <= 0
            or item["accepted_stage_count"] <= 0
            or item["CFL_retry_count"] < 0
            or item["source_retry_count"] != 0
            or not isinstance(item["state_sha256"], str)
            or len(item["state_sha256"]) != 64
        ):
            raise ValueError(f"RSP1 terminal-member contract differs: {key}")
        answer[key] = dict(item)
    return answer


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    with path.open("rb") as handle:
        config = tomllib.load(handle)
    _strict_keys(
        "RSP1/PREF12 config",
        config,
        {
            "schema_version",
            "artifact_id",
            "project_version",
            "metric_signature",
            "riemann_convention",
            "runtime_authorization_result",
            "run_plan_config",
            "study_manifest",
            "study_event_log",
            "study_checkpoint",
            "study_result",
            "scope",
            "immutable_study",
            "initial_target",
            "evolved_target",
            "terminal_members",
            "proof_contract",
            "successor_boundary",
            "claims",
        },
    )
    scope = config["scope"]
    immutable = config["immutable_study"]
    initial = config["initial_target"]
    evolved = config["evolved_target"]
    successor = config["successor_boundary"]
    if (
        config["schema_version"] != 1
        or config["artifact_id"] != ARTIFACT_ID
        or scope
        != {
            "role": "post_run_amplitude_three_direct_resolution_spectrum_reduction",
            "branch": "GR-0",
            "amplitude": "3",
            "RSP1_trajectory_read": True,
            "SGBL_trajectory_read": False,
            "FGCQR_trajectory_read": False,
            "collapse_or_trapped_outcome_classified": False,
            "mechanism_question_answered": False,
        }
        or immutable.get("authorization_commit")
        != "7f448bf6e5f0736e025c484deb664ba7887269da"
        or immutable.get("study_id")
        != "0de794592ecafb7b93aff5b316d5b393a4007008c877535baf826a3c8d45f53f"
        or immutable.get("terminal_classification") != EXPECTED_CLASSIFICATION
        or immutable.get("event_log_line_count") != 2
        or immutable.get("target_coordinate_time") != "1/16"
        or immutable.get("methods_in_order") != list(EXPECTED_METHODS)
        or immutable.get("point_counts_in_order")
        != list(EXPECTED_POINT_COUNTS)
        or immutable.get("source_retry_total") != 0
        or initial.get("both_pairs_fail_strict_one_quarter_ceiling") is not True
        or initial.get("endpoint_outcome_not_read_at_authorization") is not True
        or evolved.get("field") != "phi_over_Lambda"
        or evolved.get("quantity")
        != "top_band_derivative_weighted_power_fraction"
        or evolved.get("direct_ratio_ceiling") != "1/4"
        or evolved.get("both_pairs_pass_strict_one_quarter_ceiling") is not True
        or evolved.get("both_constraint_admissions_pass") is not True
        or evolved.get("both_complete_profile_contractions_pass") is not True
        or evolved.get("both_target_absolute_budget_guards_pass") is not True
        or evolved.get("generic_all_field_raw_spectral_admission_passes")
        is not False
        or evolved.get("generic_every_nested_tail_ratio_passes") is not False
        or any(value is not True for value in config["proof_contract"].values())
        or successor.get("next_gate") != "FGC-1-PRO12-FRZ1"
        or successor.get("PROTO12_design_may_begin") is not True
        or any(
            value is not True
            for key, value in successor.items()
            if key != "next_gate"
        )
        or config["claims"] != EXPECTED_CLAIMS
    ):
        raise ValueError("RSP1/PREF12 configuration contract differs")
    for name in (
        "manifest_sha256",
        "event_log_sha256",
        "checkpoint_sha256",
        "study_result_sha256",
    ):
        if not isinstance(immutable.get(name), str) or len(immutable[name]) != 64:
            raise ValueError(f"RSP1 immutable hash differs: {name}")
    for method in EXPECTED_METHODS:
        initial_ratios = [_as_float(item) for item in initial[f"{method}_ratios"]]
        evolved_ratios = [_as_float(item) for item in evolved[f"{method}_ratios"]]
        if (
            len(initial_ratios) != 2
            or not all(item >= 0.25 for item in initial_ratios)
            or len(evolved_ratios) != 2
            or not all(item < 0.25 for item in evolved_ratios)
        ):
            raise ValueError(f"RSP1 configured target ratios differ: {method}")
    _expected_terminal_members(config)
    return config


def _validate_initial_event(
    authorization: Mapping[str, Any], initial_event: Mapping[str, Any]
) -> None:
    assessment = initial_event["assessment"]
    frozen_events = {
        item["method"]: item
        for item in authorization["artifact_payload"]["initial_premise_events"]
    }
    if set(frozen_events) != set(EXPECTED_METHODS):
        raise ValueError("RSP1 authorization initial events differ")
    if (
        assessment["primary_common_event"]
        != frozen_events["RK4"]["complete_common_event"]
        or assessment["comparator_common_event"]
        != frozen_events["SSPRK3"]["complete_common_event"]
    ):
        raise ValueError("RSP1 initial event differs from the authorization")
    frozen_inputs = {
        f"{item['method']}-{item['point_count']}": item["projected_state_sha256"]
        for item in authorization["artifact_payload"]["frozen_run_inputs"]
    }
    members = assessment["member_diagnostics"]
    if set(frozen_inputs) != set(EXPECTED_MEMBER_KEYS) or set(members) != set(
        EXPECTED_MEMBER_KEYS
    ):
        raise ValueError("RSP1 t0 member set differs")
    for key in EXPECTED_MEMBER_KEYS:
        item = members[key]
        if (
            item["step_index"] != 0
            or item["transaction_serial"] != 0
            or item["accepted_stage_count"] != 0
            or item["CFL_retry_count"] != 0
            or item["source_retry_count"] != 0
            or item["state_sha256"] != frozen_inputs[key]
        ):
            raise ValueError(f"RSP1 t0 member differs: {key}")


def _validate_study(
    config: Mapping[str, Any],
) -> tuple[
    dict[str, Any],
    list[dict[str, Any]],
    bytes,
    dict[str, Any],
    dict[str, Any],
    dict[str, Any],
]:
    immutable = config["immutable_study"]
    manifest_path = REPOSITORY / config["study_manifest"]
    event_path = REPOSITORY / config["study_event_log"]
    checkpoint_path = REPOSITORY / config["study_checkpoint"]
    result_path = REPOSITORY / config["study_result"]
    for path, expected in (
        (manifest_path, immutable["manifest_sha256"]),
        (event_path, immutable["event_log_sha256"]),
        (checkpoint_path, immutable["checkpoint_sha256"]),
        (result_path, immutable["study_result_sha256"]),
    ):
        if not path.is_file() or _sha(path) != expected:
            raise ValueError(f"immutable RSP1 study hash differs: {_rel(path)}")

    manifest = _load_json(manifest_path)
    result = _load_json(result_path)
    events, event_bytes = _load_events(event_path)
    authorization_path = REPOSITORY / config["runtime_authorization_result"]
    plan_path = REPOSITORY / config["run_plan_config"]
    authorization = _load_json(authorization_path, "FGC-1-RSP1-FRZ1")
    if (
        manifest.get("schema_version") != 1
        or manifest.get("runner_id") != runner.RUNNER_ID
        or manifest.get("study_id") != immutable["study_id"]
        or manifest.get("git_commit") != immutable["authorization_commit"]
        or manifest.get("plan_path") != _rel(plan_path)
        or manifest.get("plan_sha256") != _sha(plan_path)
        or manifest.get("authorization_path") != _rel(authorization_path)
        or manifest.get("authorization_sha256") != _sha(authorization_path)
        or manifest.get("point_counts") != list(EXPECTED_POINT_COUNTS)
        or manifest.get("methods") != list(EXPECTED_METHODS)
        or manifest.get("amplitude") != "3"
        or manifest.get("outcome_fields_present_at_creation") is not False
        or result.get("event_log_sha256") != immutable["event_log_sha256"]
        or len(events) != immutable["event_log_line_count"]
    ):
        raise ValueError("RSP1 manifest/result contract differs")
    frozen_plan = authorization["artifact_payload"]["frozen_run_plan"]
    if (
        manifest.get("input_manifest_sha256")
        != frozen_plan["input_manifest_sha256"]
        or frozen_plan["study_id"] != immutable["study_id"]
        or authorization["gate_status"]["RSP1_execution_authorized"] is not True
        or authorization["gate_status"]["amplitude_three_spectral_veto_cleared"]
        is not False
    ):
        raise ValueError("RSP1 authorization binding differs")
    _validate_commit_lineage(
        immutable["authorization_commit"],
        manifest.get("implementation_sha256", {}),
    )
    _validate_initial_event(authorization, events[0])

    plan, reproduced_authorization = runner._validate_authorization(
        plan_path,
        authorization_path,
        require_fresh_namespace=False,
    )
    if reproduced_authorization != authorization:
        raise ValueError("RSP1 authorization replay differs")
    metadata, members, checkpoint_event_bytes = runner._restore_checkpoint(
        checkpoint_path,
        plan=plan,
        authorization=authorization,
    )
    if (
        checkpoint_event_bytes != event_bytes
        or metadata.get("terminal") is not True
        or metadata.get("campaign_id") != immutable["study_id"]
        or metadata.get("study_id") != immutable["study_id"]
        or metadata.get("amplitude") != "3"
        or metadata.get("completed_common_event_index") != 1
        or metadata.get("event_log")
        != inherited._event_log_identity(event_bytes)
        or metadata.get("terminal_result")
        != {
            "classification": EXPECTED_CLASSIFICATION,
            "amplitude_three_spectral_veto_cleared_for_successor_design": True,
        }
        or set(metadata.get("members", {})) != set(EXPECTED_MEMBER_KEYS)
    ):
        raise ValueError("RSP1 terminal checkpoint metadata differs")

    recomputed_endpoint = inherited._serial(
        runner._event_assessment(
            plan,
            members,
            _as_float(immutable["target_coordinate_time"]),
        )
    )
    if (
        recomputed_endpoint != result.get("endpoint_assessment")
        or recomputed_endpoint != events[-1].get("assessment")
    ):
        raise ValueError("RSP1 checkpoint endpoint recomputation differs")
    return (
        manifest,
        events,
        event_bytes,
        result,
        authorization,
        recomputed_endpoint,
    )


def _validate_configured_diagnosis(
    config: Mapping[str, Any], diagnosis: Mapping[str, Any]
) -> None:
    initial = config["initial_target"]
    evolved = config["evolved_target"]
    configured_members = _expected_terminal_members(config)
    observed_members = diagnosis["members"]
    if (
        diagnosis["classification"] != EXPECTED_CLASSIFICATION
        or diagnosis["event_log_line_count"]
        != config["immutable_study"]["event_log_line_count"]
        or diagnosis["total_source_retry_count"] != 0
        or diagnosis["all_six_members_reached_endpoint"] is not True
        or diagnosis["amplitude_three_spectral_veto_cleared_for_successor_design"]
        is not True
        or diagnosis["generic_all_field_raw_spectral_admission_passed"]
        is not False
        or diagnosis["GR0_case_eligible"] is not False
        or diagnosis["candidate_or_mechanism_outcome_read"] is not False
    ):
        raise ValueError("configured RSP1 diagnosis differs")
    for method in EXPECTED_METHODS:
        initial_method = diagnosis["initial_premise"]["methods"][method]
        evolved_method = diagnosis["evolved_endpoint"]["methods"][method]
        if (
            initial_method["target_derivative_tail_ratios"]
            != [_as_float(item) for item in initial[f"{method}_ratios"]]
            or evolved_method["target_derivative_tail_ratios"]
            != [_as_float(item) for item in evolved[f"{method}_ratios"]]
            or evolved_method["minimum_finite_finest_pair_order"]
            != _as_float(evolved[f"{method}_minimum_finite_constraint_order"])
            or not all(evolved_method["target_direct_pair_passed"])
            or evolved_method["generic_raw_all_field_spectral_admission_passed"]
            is not False
        ):
            raise ValueError(f"configured RSP1 method diagnosis differs: {method}")
    if set(observed_members) != set(configured_members):
        raise ValueError("configured RSP1 member diagnosis differs")
    for key, expected in configured_members.items():
        observed = observed_members[key]
        for field in (
            "step_index",
            "transaction_serial",
            "accepted_stage_count",
            "CFL_retry_count",
            "source_retry_count",
            "state_sha256",
        ):
            if observed[field] != expected[field]:
                raise ValueError(f"configured RSP1 member differs: {key} {field}")


def reproduce(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    config = load_config(config_path)
    (
        manifest,
        events,
        event_bytes,
        study_result,
        authorization,
        recomputed_endpoint,
    ) = _validate_study(config)
    diagnosis = diagnose_rsp1_campaign(
        events,
        study_result,
        study_id=config["immutable_study"]["study_id"],
        direct_tail_ceiling=_as_float(
            config["evolved_target"]["direct_ratio_ceiling"]
        ),
    )
    _validate_configured_diagnosis(config, diagnosis)
    authorization_path = REPOSITORY / config["runtime_authorization_result"]
    run_plan = REPOSITORY / config["run_plan_config"]
    record = {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": config["project_version"],
        "classification": (
            "post_run_amplitude_three_direct_spectrum_clearance_for_successor_design"
        ),
        "generated_by": "scripts/reproduce_fgc_rsp1_pref12.py",
        "derivation_document": _rel(OWNER_DOCUMENT),
        "derivation_document_sha256": _sha(OWNER_DOCUMENT),
        "source_config_sha256": {_rel(config_path): _sha(config_path)},
        "predecessor_sha256": {
            _rel(authorization_path): _sha(authorization_path),
            _rel(run_plan): _sha(run_plan),
        },
        "implementation_sha256": {_rel(path): _sha(path) for path in IMPLEMENTATION},
        "scope_bindings": dict(config["scope"]),
        "gate_status": dict(EXPECTED_CLAIMS),
        "artifact_payload": {
            "immutable_study": {
                "authorization_commit": config["immutable_study"][
                    "authorization_commit"
                ],
                "study_id": manifest["study_id"],
                "manifest_sha256": config["immutable_study"]["manifest_sha256"],
                "event_log_sha256": _bytes_sha(event_bytes),
                "checkpoint_sha256": config["immutable_study"][
                    "checkpoint_sha256"
                ],
                "study_result_sha256": config["immutable_study"][
                    "study_result_sha256"
                ],
                "elapsed_wall_seconds": study_result["elapsed_wall_seconds"],
                "classification": study_result["classification"],
            },
            "checkpoint_recomputation": {
                "all_six_terminal_states_restored": True,
                "checkpoint_event_log_equals_raw_event_log": True,
                "recomputed_endpoint_equals_terminal_event": True,
                "recomputed_endpoint_equals_study_result": True,
                "recomputed_endpoint_sha256": _bytes_sha(
                    _canonical(recomputed_endpoint).encode("utf-8")
                ),
                "authorization_reproduced": True,
                "initial_event_equals_frozen_authorization": True,
            },
            "study_diagnosis": diagnosis,
            "successor_boundary": dict(config["successor_boundary"]),
            "epistemic_boundary": {
                "observed_results": [
                    "all six frozen amplitude-3 GR-0 members reached t=1/16 with zero source retries under SRC3",
                    "both direct adjacent phi/Lambda derivative-tail ratios are strictly below one quarter for RK4 and SSPRK3",
                    "both constraint admissions, complete profile contractions, and target absolute-budget guards pass",
                    "the terminal checkpoint independently recomputes the exact endpoint serialized by the event log and result",
                ],
                "still_public_limits": [
                    "the generic all-field raw nested-tail aggregator is false for both methods and is not relabelled by RSP1",
                    "RSP1 is one amplitude, one endpoint, spherical GR-0, and a targeted numerical resolution question",
                    "CFL trial retries are numerical step selection and are serialized separately from the zero source-retry count",
                ],
                "not_observed": [
                    "a completed fresh GR-0 calibration or eligible calibration case",
                    "a collapse, trapped sphere, regulator activation, or metric-null defocusing result",
                    "SGB-L or FGC-QR dynamics",
                    "a retained-EFT, singularity-resolution, child-domain, dark-sector, or varying-local-c claim",
                ],
                "successor_status": (
                    "PROTO12 design may begin, but its complete fresh-run contract must be frozen before any trajectory is read"
                ),
            },
            "authorization_identity": {
                "artifact_id": authorization["artifact_id"],
                "study_id": authorization["artifact_payload"]["frozen_run_plan"][
                    "study_id"
                ],
            },
        },
        "nonclaims": {
            key: value for key, value in EXPECTED_CLAIMS.items() if value is False
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
