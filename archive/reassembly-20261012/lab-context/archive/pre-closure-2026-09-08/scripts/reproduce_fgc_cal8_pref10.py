#!/usr/bin/env python3
"""Reproduce CAL8/PREF10's outcome-neutral PROTO11 campaign diagnosis."""

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

from recursive_horizons.fgc.evolution.cal8_proto11_campaign_diagnosis import (  # noqa: E402
    EXPECTED_AMPLITUDES,
    diagnose_proto11_campaign,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    array_content_sha256,
)


ARTIFACT_ID = "FGC-1-CAL8-PREF10"
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-cal8-pref10.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-cal8-pref10.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-cal8-pref10.md"
DIAGNOSIS_MODULE = (
    REPOSITORY
    / "src/recursive_horizons/fgc/evolution/cal8_proto11_campaign_diagnosis.py"
)
IMPLEMENTATION = (Path(__file__).resolve(), DIAGNOSIS_MODULE)
Q = Fraction


EXPECTED_CLAIMS = {
    "PROTO11_campaign_terminated_normally": True,
    "PROTO11_both_amplitudes_reached_evolved_common_event": True,
    "PROTO11_GR0_case_eligible": False,
    "PROTO11_evolved_affine_source_arithmetic_obstruction_localized": True,
    "PROTO11_direct_coarse_phi_derivative_veto_localized": True,
    "SRC2_affine_source_arithmetic_preflight_required": True,
    "PROTO12_frozen": False,
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
    if _git("merge-base", "--is-ancestor", commit, "HEAD", check=False).returncode:
        raise ValueError("campaign authorization commit is not an ancestor of HEAD")
    if not implementation:
        raise ValueError("campaign implementation lineage is absent")
    for relative, expected in implementation.items():
        blob = _git("show", f"{commit}:{relative}").stdout
        if _bytes_sha(blob) != expected:
            raise ValueError(f"immutable campaign implementation differs: {relative}")


def _strict_keys(name: str, value: Mapping[str, Any], expected: set[str]) -> None:
    if set(value) != expected:
        raise ValueError(f"{name} keys differ")


def _as_float(value: str) -> float:
    return float(Q(value))


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    with path.open("rb") as handle:
        config = tomllib.load(handle)
    _strict_keys(
        "CAL8 config",
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
            "source_obstruction",
            "spectral_obstruction",
            "prospective_revision",
            "proof_contract",
            "claims",
        },
    )
    scope = config["scope"]
    immutable = config["immutable_campaign"]
    source = config["source_obstruction"]
    spectral = config["spectral_obstruction"]
    revision = config["prospective_revision"]
    proof = config["proof_contract"]
    if (
        config["schema_version"] != 1
        or config["artifact_id"] != ARTIFACT_ID
        or scope.get("target_protocol") != "FGC-2-SF1-PROTO11"
        or scope.get("role")
        != "post_calibration_evolved_source_arithmetic_and_direct_spectral_obstruction_diagnosis"
        or scope.get("calibration_branch") != "GR-0"
        or scope.get("GR0_calibration_trajectory_read") is not True
        or scope.get("SGBL_trajectory_read") is not False
        or scope.get("FGCQR_trajectory_read") is not False
        or scope.get("collapse_or_trapped_outcome_classified") is not False
        or scope.get("mechanism_question_answered") is not False
        or immutable.get("authorization_commit")
        != "199c86e55a7504e780ccce5a840648528482f819"
        or immutable.get("campaign_id")
        != "63ea00c3a781d3a024c69d9841cadc74bf994cd7ed3750bc12cd29bbf24c6548"
        or immutable.get("terminal_classification")
        != "calibration_failed_no_eligible_GR0_case"
        or immutable.get("amplitudes_in_order") != list(EXPECTED_AMPLITUDES)
        or immutable.get("event_log_line_count") != 168
        or immutable.get("common_event_count") != 4
        or immutable.get("source_only_rejection_count") != 164
        or immutable.get("selected_amplitude_absent") is not True
        or immutable.get("holdout_execution_authorized") is not False
        or source.get("amplitude") != "5/2"
        or source.get("first_evolved_common_event_index") != 1
        or source.get("first_evolved_common_event_time") != "1/16"
        or source.get("first_evolved_common_event_fully_admitted") is not True
        or source.get("target_common_event_index") != 2
        or source.get("target_common_event_time") != "1/8"
        or source.get("terminal_member") != "RK4-8193"
        or source.get("terminal_reason") != "newton_residual_limit"
        or source.get("terminal_classification")
        != "scientific_source_retry_exhausted"
        or source.get("terminal_exhaustion_kind") != "minimum_step_size"
        or source.get("source_only_rejection_count") != 164
        or source.get("accepted_time_group_count") != 16
        or source.get("accepted_increment_time_ulp_budget") != 2
        or source.get("last_accepted_step_index") != 97
        or source.get("last_accepted_transaction_serial") != 485
        or source.get("last_retry_count_for_accepted_step") != 20
        or source.get("cumulative_source_retry_count") != 164
        or source.get("failed_retry_counts_by_accepted_time")
        != [1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 15, 17, 18, 19, 20]
        or spectral.get("amplitude") != "3"
        or spectral.get("evolved_common_event_index") != 1
        or spectral.get("evolved_common_event_time") != "1/16"
        or spectral.get("source_only_rejection_count") != 0
        or spectral.get("both_constraint_admissions_passed") is not True
        or spectral.get("both_finest_pair_individual_budgets_passed") is not True
        or spectral.get("both_whole_profile_contractions_passed") is not True
        or spectral.get("failed_direct_coarse_to_medium_quantity")
        != "phi_over_Lambda_derivative_tail"
        or spectral.get("direct_tail_ceiling") != "1/4"
        or spectral.get(
            "all_other_coarse_to_medium_direct_field_and_derivative_tails_passed"
        )
        is not True
        or spectral.get("trapped_sign_passed") is not False
        or spectral.get("temporal_spectrum_available") is not False
        or revision.get("next_gate") != "FGC-1-SRC2-PREF11"
        or any(value is not True for key, value in revision.items() if key not in {"next_gate", "purpose"})
        or revision.get("purpose")
        != "capture_and_attack_the_evolved_affine_source_arithmetic_wall_before_any_successor_protocol"
        or any(value is not True for value in proof.values())
        or config["claims"] != EXPECTED_CLAIMS
    ):
        raise ValueError("CAL8/PREF10 configuration contract differs")
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
        or manifest.get("runner_id") != "FGC-1-CAL8-RUN1-RUNNER"
        or manifest.get("ordered_amplitudes") != list(EXPECTED_AMPLITUDES)
        or manifest.get("outcome_fields_present_at_creation") is not False
        or result.get("campaign_id") != immutable["campaign_id"]
        or result.get("runner_id") != "FGC-1-CAL8-RUN1-RUNNER"
        or result.get("classification") != immutable["terminal_classification"]
        or result.get("event_log_sha256") != immutable["event_log_sha256"]
        or result.get("selected_amplitude") is not None
        or result.get("holdout_execution_authorized") is not False
        or result.get("mechanism_question_answered") is not False
        or result.get("SGBL_outcome_read") is not False
        or result.get("FGCQR_outcome_read") is not False
        or result.get("retained_EFT_evolution_authorized") is not False
        or len(events) != immutable["event_log_line_count"]
    ):
        raise ValueError("immutable campaign manifest/result contract differs")
    authorization_path = REPOSITORY / config["runtime_authorization_result"]
    plan_path = REPOSITORY / config["run_plan_config"]
    authorization = _load_json(authorization_path, "FGC-1-HLT9-MON9")
    if (
        _sha(authorization_path) != manifest.get("authorization_sha256")
        or _sha(plan_path) != manifest.get("plan_sha256")
        or _rel(authorization_path) != manifest.get("authorization_path")
        or _rel(plan_path) != manifest.get("plan_path")
        or authorization.get("gate_status", {}).get(
            "PROTO11_fresh_GR0_dynamic_calibration_authorized"
        )
        is not True
        or authorization.get("gate_status", {}).get(
            "FGCQR_holdout_execution_authorized"
        )
        is not False
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
            or metadata.get("amplitude") != "3"
            or metadata.get("amplitude_index") != 1
            or metadata.get("completed_common_event_index") != 1
            or metadata.get("event_log", {}).get("sha256")
            != immutable["event_log_sha256"]
            or metadata.get("event_log", {}).get("line_count") != len(events)
            or metadata.get("completed_amplitude_records")
            != result.get("amplitude_records")
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


def _validate_configured_diagnosis(
    config: Mapping[str, Any], diagnosis: Mapping[str, Any]
) -> None:
    immutable = config["immutable_campaign"]
    source_config = config["source_obstruction"]
    spectral_config = config["spectral_obstruction"]
    source = diagnosis["amplitude_five_halves_source_obstruction"]
    spectral = diagnosis["amplitude_three_direct_spectral_obstruction"]
    groups = source["groups"]
    retry_counts = [item["failed_retry_count"] for item in groups]
    methods = spectral["methods"]
    expected_phi = {
        "RK4": [
            _as_float(source_value)
            for source_value in (
                spectral_config["RK4_coarse_to_medium_phi_derivative_tail_ratio"],
                spectral_config["RK4_medium_to_fine_phi_derivative_tail_ratio"],
            )
        ],
        "SSPRK3": [
            _as_float(source_value)
            for source_value in (
                spectral_config[
                    "SSPRK3_coarse_to_medium_phi_derivative_tail_ratio"
                ],
                spectral_config[
                    "SSPRK3_medium_to_fine_phi_derivative_tail_ratio"
                ],
            )
        ],
    }
    if (
        diagnosis["event_log_line_count"] != immutable["event_log_line_count"]
        or diagnosis["common_event_count"] != immutable["common_event_count"]
        or diagnosis["source_only_rejection_count"]
        != immutable["source_only_rejection_count"]
        or diagnosis["both_amplitudes_reached_one_evolved_common_event"] is not True
        or diagnosis["no_GR0_amplitude_eligible"] is not True
        or diagnosis["candidate_or_mechanism_outcome_read"] is not False
        or source["rejected_source_only_proposal_count"]
        != source_config["source_only_rejection_count"]
        or source["accepted_time_group_count"]
        != source_config["accepted_time_group_count"]
        or retry_counts != source_config["failed_retry_counts_by_accepted_time"]
        or source["last_accepted_time"] != _as_float(source_config["last_accepted_time"])
        or source["last_failed_stage_time"]
        != _as_float(source_config["last_failed_stage_time"])
        or source["last_legal_failed_step_size"]
        != _as_float(source_config["last_legal_failed_step_size"])
        or source["next_forbidden_step_size"]
        != _as_float(source_config["next_forbidden_step_size"])
        or source["frozen_minimum_step_size"]
        != _as_float(source_config["minimum_step_size"])
        or source["raw_source_residual_limit"]
        != _as_float(source_config["raw_source_residual_maximum"])
        or source["last_accepted_step_index"]
        != source_config["last_accepted_step_index"]
        or source["last_accepted_transaction_serial"]
        != source_config["last_accepted_transaction_serial"]
        or spectral["sole_direct_coarse_to_medium_veto"]
        != "phi_over_Lambda derivative tail"
        or spectral["direct_tail_ceiling"]
        != _as_float(spectral_config["direct_tail_ceiling"])
        or any(
            methods[method]["phi_derivative_tail_ratios"] != ratios
            for method, ratios in expected_phi.items()
        )
    ):
        raise ValueError("configured CAL8 campaign diagnosis differs")


def reproduce(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    config = load_config(config_path)
    manifest, events, event_bytes, campaign_result = _validate_campaign(config)
    source_config = config["source_obstruction"]
    spectral_config = config["spectral_obstruction"]
    diagnosis = diagnose_proto11_campaign(
        events,
        campaign_result,
        raw_source_limit=_as_float(source_config["raw_source_residual_maximum"]),
        minimum_step_size=_as_float(source_config["minimum_step_size"]),
        direct_tail_ceiling=_as_float(spectral_config["direct_tail_ceiling"]),
        accepted_increment_time_ulp_budget=source_config[
            "accepted_increment_time_ulp_budget"
        ],
    )
    _validate_configured_diagnosis(config, diagnosis)
    runtime_authorization = REPOSITORY / config["runtime_authorization_result"]
    run_plan = REPOSITORY / config["run_plan_config"]
    record = {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": config["project_version"],
        "classification": (
            "post_calibration_PROTO11_evolved_source_arithmetic_and_direct_spectral_obstruction_diagnosis"
        ),
        "generated_by": "scripts/reproduce_fgc_cal8_pref10.py",
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
                "authorization_commit": config["immutable_campaign"][
                    "authorization_commit"
                ],
                "campaign_id": manifest["campaign_id"],
                "manifest_sha256": config["immutable_campaign"]["manifest_sha256"],
                "event_log_sha256": _bytes_sha(event_bytes),
                "checkpoint_sha256": config["immutable_campaign"][
                    "checkpoint_sha256"
                ],
                "campaign_result_sha256": config["immutable_campaign"][
                    "campaign_result_sha256"
                ],
                "elapsed_wall_seconds": campaign_result["elapsed_wall_seconds"],
                "classification": campaign_result["classification"],
            },
            "campaign_diagnosis": diagnosis,
            "prospective_revision": dict(config["prospective_revision"]),
            "epistemic_boundary": {
                "observed_results": [
                    "amplitude 5/2 passed one fully admitted evolved common event before a repeatable RK4-8193 affine source residual wall exhausted the frozen binary-halving budget",
                    "amplitude 3 reached the same evolved event without source retry and failed only the unchanged direct coarse-to-medium phi derivative-tail veto",
                ],
                "localized_numerical_questions": [
                    "whether stable affine source arithmetic can satisfy the unchanged complete raw 1e-12 residual gate",
                    "whether the amplitude-3 coarse phi derivative tail is pre-asymptotic or a persistent continuum-resolution obstruction",
                ],
                "not_observed": [
                    "an eligible GR-0 calibration case",
                    "a qualified trapped common event",
                    "SGB-L or FGC-QR dynamics",
                    "regulator activation or metric-null defocusing",
                    "a continuum branch or hyperbolicity loss",
                    "a retained-EFT or physical transition",
                ],
                "successor_status": (
                    "SRC2 replay and arithmetic preflight is required; PROTO12 is not frozen"
                ),
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
