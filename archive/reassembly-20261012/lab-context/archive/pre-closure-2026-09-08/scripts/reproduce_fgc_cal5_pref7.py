#!/usr/bin/env python3
"""Reproduce the CAL5 diagnosis of the PROTO8 GR-0 calibration stops.

The certificate binds the immutable HLT6-authorized campaign, verifies its
terminal checkpoint and public event sequence, and reduces the failed
common-event gates without advancing a trajectory.  It authorizes no
candidate branch; its only successor decision is that a separately frozen
one-level resolution-ladder experiment is warranted under unchanged rules.
"""

from __future__ import annotations

import argparse
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

from recursive_horizons.fgc.evolution.cal5_resolution_diagnosis import (  # noqa: E402
    diagnose_proto8_campaign_events,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    array_content_sha256,
)


ARTIFACT_ID = "FGC-1-CAL5-PREF7"
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-cal5-pref7.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-cal5-pref7.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-cal5-pref7.md"
DIAGNOSIS_MODULE = (
    REPOSITORY
    / "src/recursive_horizons/fgc/evolution/cal5_resolution_diagnosis.py"
)
IMPLEMENTATION = (Path(__file__).resolve(), DIAGNOSIS_MODULE)

EXPECTED_CLAIMS = {
    "PROTO8_campaign_terminated_normally": True,
    "PROTO8_common_event_ownership_repair_exercised": True,
    "PROTO8_GR0_case_eligible": False,
    "PROTO9_resolution_ladder_revision_required": True,
    "PROTO9_frozen": False,
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
            raise ValueError("campaign event records must be JSON objects")
        records.append(value)
        rebuilt.extend(_canonical_line(value))
    if bytes(rebuilt) != raw:
        raise ValueError("campaign event log is not canonical JSONL")
    return records


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
            raise ValueError(f"immutable implementation blob differs: {relative}")


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    with path.open("rb") as handle:
        config = tomllib.load(handle)
    scope = config.get("scope", {})
    immutable = config.get("immutable_campaign", {})
    diagnosis = config.get("resolution_diagnosis", {})
    revision = config.get("prospective_revision", {})
    proof = config.get("proof_contract", {})
    if (
        config.get("schema_version") != 1
        or config.get("artifact_id") != ARTIFACT_ID
        or scope.get("target_protocol") != "FGC-2-SF1-PROTO8"
        or scope.get("role")
        != "post_calibration_resolution_ladder_adequacy_diagnosis"
        or scope.get("GR0_calibration_trajectory_read") is not True
        or scope.get("SGBL_trajectory_read") is not False
        or scope.get("FGCQR_trajectory_read") is not False
        or scope.get("mechanism_question_answered") is not False
        or immutable.get("terminal_classification")
        != "calibration_failed_no_eligible_GR0_case"
        or immutable.get("amplitudes_in_order") != ["5/2", "3"]
        or immutable.get("event_log_line_count") != 9
        or immutable.get("source_only_retry_count") != 4
        or immutable.get("source_only_retry_member") != "RK4-4097"
        or immutable.get("amplitude_5over2_stop_time") != "1/8"
        or immutable.get("amplitude_3_stop_time") != "1/16"
        or diagnosis.get("current_point_counts") != [1025, 2049, 4097]
        or diagnosis.get("smallest_prospective_point_counts")
        != [2049, 4097, 8193]
        or diagnosis.get("minimum_finest_pair_order") != "3/2"
        or diagnosis.get("maximum_nested_tail_ratio") != "1/4"
        or diagnosis.get("conditional_Richardson_projection_is_planning_only")
        is not True
        or diagnosis.get(
            "new_8193_constraint_and_spectrum_data_required_before_admission"
        )
        is not True
        or diagnosis.get("no_threshold_change") is not True
        or revision.get("new_protocol_required") != "FGC-2-SF1-PROTO9"
        or revision.get("only_resolution_ladder_may_change") is not True
        or revision.get("point_counts") != [2049, 4097, 8193]
        or revision.get("fresh_output_namespace_required") is not True
        or any(
            revision.get(key) is not True
            for key in (
                "physical_inputs_unchanged",
                "amplitude_order_unchanged",
                "evolution_equations_unchanged",
                "source_solver_and_raw_residual_tolerance_unchanged",
                "methods_CFL_and_retry_rules_unchanged",
                "constraint_ownership_magnitude_and_order_rules_unchanged",
                "spectral_ownership_budgets_and_tail_rules_unchanged",
                "trapped_sign_boundary_and_physical_rules_unchanged",
            )
        )
        or any(value is not True for value in proof.values())
        or config.get("claims") != EXPECTED_CLAIMS
    ):
        raise ValueError("CAL5/PREF7 config violates its frozen semantic boundary")
    commit = immutable.get("authorization_commit")
    if not isinstance(commit, str) or len(commit) != 40:
        raise ValueError("immutable_campaign.authorization_commit must be a Git SHA-1")
    for key in (
        "campaign_id",
        "manifest_sha256",
        "event_log_sha256",
        "checkpoint_sha256",
        "campaign_result_sha256",
    ):
        value = immutable.get(key)
        if not isinstance(value, str) or len(value) != 64:
            raise ValueError(f"immutable_campaign.{key} must be a SHA-256 value")
    return config


def _paths(config: Mapping[str, Any]) -> dict[str, Path]:
    names = {
        "authorization": "runtime_authorization_result",
        "plan": "run_plan_config",
        "manifest": "campaign_manifest",
        "events": "campaign_event_log",
        "checkpoint": "campaign_checkpoint",
        "result": "campaign_result",
    }
    paths = {key: (REPOSITORY / config[field]).resolve() for key, field in names.items()}
    for key, path in paths.items():
        try:
            path.relative_to(REPOSITORY)
        except ValueError as exc:
            raise ValueError(f"CAL5 {key} path escapes the repository") from exc
        if not path.is_file():
            raise ValueError(f"CAL5 {key} path is missing: {_rel(path)}")
    return paths


def _validate_checkpoint(
    path: Path,
    *,
    event_bytes: bytes,
    campaign_id: str,
    amplitude_records: list[Mapping[str, Any]],
) -> dict[str, Any]:
    with np.load(path, allow_pickle=False) as archive:
        metadata = json.loads(bytes(archive["metadata_utf8"]).decode("utf-8"))
        embedded_events = bytes(archive["event_log_utf8"])
        if embedded_events != event_bytes:
            raise ValueError("terminal checkpoint event log differs")
        if (
            metadata.get("terminal") is not True
            or metadata.get("campaign_id") != campaign_id
            or metadata.get("amplitude") != "3"
            or metadata.get("amplitude_index") != 1
            or metadata.get("completed_common_event_index") != 1
            or metadata.get("completed_amplitude_records") != amplitude_records
            or metadata.get("event_log", {}).get("line_count") != 9
            or metadata.get("event_log", {}).get("sha256") != _bytes_sha(event_bytes)
        ):
            raise ValueError("terminal checkpoint metadata differs")
        expected = amplitude_records[1]["member_final_state_hashes"]
        for member, digest in expected.items():
            prefix = member.replace("-", "_")
            observed = array_content_sha256(
                archive[f"{prefix}_u"],
                archive[f"{prefix}_p"],
                archive[f"{prefix}_q"],
            )
            if observed != digest:
                raise ValueError(f"terminal checkpoint state differs: {member}")
    return metadata


def _validate_campaign(
    config: Mapping[str, Any],
) -> tuple[
    dict[str, Path],
    dict[str, Any],
    list[dict[str, Any]],
    dict[str, Any],
    dict[str, Any],
    dict[str, Any],
]:
    paths = _paths(config)
    immutable = config["immutable_campaign"]
    expected_hashes = {
        "manifest": immutable["manifest_sha256"],
        "events": immutable["event_log_sha256"],
        "checkpoint": immutable["checkpoint_sha256"],
        "result": immutable["campaign_result_sha256"],
    }
    for key, expected in expected_hashes.items():
        if _sha(paths[key]) != expected:
            raise ValueError(f"immutable PROTO8 {key} hash differs")

    authorization = _load_json(paths["authorization"], "FGC-1-HLT6-MON6")
    manifest = _load_json(paths["manifest"])
    campaign_result = _load_json(paths["result"])
    events = _load_events(paths["events"])
    event_bytes = paths["events"].read_bytes()
    if (
        manifest.get("runner_id") != "FGC-1-CAL5-RUN1-RUNNER"
        or manifest.get("campaign_id") != immutable["campaign_id"]
        or manifest.get("git_commit") != immutable["authorization_commit"]
        or manifest.get("ordered_amplitudes") != ["5/2", "3"]
        or manifest.get("plan_path") != config["run_plan_config"]
        or manifest.get("plan_sha256") != _sha(paths["plan"])
        or manifest.get("authorization_path")
        != config["runtime_authorization_result"]
        or manifest.get("authorization_sha256") != _sha(paths["authorization"])
        or manifest.get("outcome_fields_present_at_creation") is not False
        or manifest.get("implementation_sha256")
        != authorization.get("implementation_sha256")
    ):
        raise ValueError("PROTO8 manifest lineage differs")
    gates = authorization.get("gate_status", {})
    if (
        gates.get("PROTO8_fresh_GR0_dynamic_calibration_authorized") is not True
        or gates.get("FGCQR_holdout_execution_authorized") is not False
        or gates.get("retained_EFT_evolution_authorized") is not False
    ):
        raise ValueError("HLT6 authorization boundary differs")
    _validate_commit_lineage(
        immutable["authorization_commit"], manifest["implementation_sha256"]
    )
    if (
        campaign_result.get("runner_id") != manifest["runner_id"]
        or campaign_result.get("campaign_id") != immutable["campaign_id"]
        or campaign_result.get("classification")
        != immutable["terminal_classification"]
        or campaign_result.get("selected_amplitude") is not None
        or campaign_result.get("event_log_sha256") != immutable["event_log_sha256"]
        or campaign_result.get("PROTO8_common_event_ownership_enabled") is not True
        or campaign_result.get("SGBL_outcome_read") is not False
        or campaign_result.get("FGCQR_outcome_read") is not False
        or campaign_result.get("holdout_execution_authorized") is not False
        or campaign_result.get("mechanism_question_answered") is not False
        or campaign_result.get("retained_EFT_evolution_authorized") is not False
        or len(campaign_result.get("amplitude_records", ())) != 2
    ):
        raise ValueError("PROTO8 terminal result boundary differs")
    metadata = _validate_checkpoint(
        paths["checkpoint"],
        event_bytes=event_bytes,
        campaign_id=immutable["campaign_id"],
        amplitude_records=campaign_result["amplitude_records"],
    )
    diagnosis = diagnose_proto8_campaign_events(events)
    stops = campaign_result["amplitude_records"]
    if (
        stops[0]["amplitude"] != "5/2"
        or stops[0]["stop"]
        != {
            "classification": "common_event_constraint_or_spectral_stop",
            "coordinate_time": 0.125,
            "event_index": 2,
            "reason": "common_event_constraint_or_spatial_spectral_admission",
        }
        or stops[1]["amplitude"] != "3"
        or stops[1]["stop"]
        != {
            "classification": "common_event_constraint_or_spectral_stop",
            "coordinate_time": 0.0625,
            "event_index": 1,
            "reason": "common_event_constraint_or_spatial_spectral_admission",
        }
    ):
        raise ValueError("PROTO8 amplitude stop records differ")
    return paths, manifest, events, campaign_result, metadata, diagnosis


def record(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    config = load_config(config_path)
    paths, manifest, events, campaign_result, metadata, diagnosis = _validate_campaign(
        config
    )
    immutable = config["immutable_campaign"]
    return {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": config["project_version"],
        "metric_signature": config["metric_signature"],
        "riemann_convention": config["riemann_convention"],
        "classification": (
            "post_PROTO8_GR0_resolution_ladder_adequacy_diagnosis"
        ),
        "generated_by": _rel(Path(__file__)),
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
                "selected_amplitude": campaign_result["selected_amplitude"],
                "amplitude_stop_records": campaign_result["amplitude_records"],
                "manifest_implementation_sha256": manifest[
                    "implementation_sha256"
                ],
                "terminal_checkpoint_amplitude": metadata["amplitude"],
                "terminal_checkpoint_event_index": metadata[
                    "completed_common_event_index"
                ],
            },
            "resolution_ladder_diagnosis": diagnosis,
            "decision": {
                "PROTO8_ownership_repair_exercised_by_evolved_data": True,
                "PROTO8_GR0_case_eligible": False,
                "current_failure_is_a_candidate_or_gradient_obstruction": False,
                "smallest_discriminating_successor": (
                    "shift_the_entire_nested_grid_ladder_one_level_while_"
                    "retaining_every_threshold_and_physical_input"
                ),
                "new_protocol_required": "FGC-2-SF1-PROTO9",
                "prospective_revision": dict(config["prospective_revision"]),
            },
            "epistemic_boundary": {
                "conditional_Richardson_forecasts_are_observations": False,
                "new_8193_data_exist": False,
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
            "new_8193_resolution_admitted": False,
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


def load_canonical_result(path: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    return _load_json(path, ARTIFACT_ID)


def verify_canonical(
    path: Path = DEFAULT_OUTPUT,
    config_path: Path = DEFAULT_CONFIG,
) -> None:
    actual = load_canonical_result(path)
    expected = record(config_path)
    if actual != expected:
        raise ValueError("stored CAL5/PREF7 result differs from reconstruction")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    config = args.config.resolve()
    output = args.output.resolve()
    if args.check:
        verify_canonical(output, config)
        print(f"verified {_rel(output)} (PROTO9 frozen=false; mechanism answered=false)")
        return
    payload = record(config)
    _atomic_write(output, _canonical(payload).encode("utf-8"))
    print(
        f"wrote {_rel(output)} "
        "(PROTO9 resolution-ladder revision required=true; mechanism answered=false)"
    )


if __name__ == "__main__":
    main()
