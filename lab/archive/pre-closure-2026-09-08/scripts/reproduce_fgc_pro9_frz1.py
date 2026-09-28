#!/usr/bin/env python3
"""Reproduce the immutable, outcome-neutral PROTO9 resolution freeze."""

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


REPOSITORY = Path(__file__).resolve().parents[1]
import sys

sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.protocol_v8 import (  # noqa: E402
    validate_sf1_protocol_v8,
)
from recursive_horizons.fgc.evolution.protocol_v9 import (  # noqa: E402
    PROTO9_CLAIMS,
    validate_sf1_protocol_v9,
)


ARTIFACT_ID = "FGC-1-PRO9-FRZ1"
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-pro9-frz1.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-pro9-frz1.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-pro9-frz1.md"
PROTOCOL_MODULE = REPOSITORY / "src/recursive_horizons/fgc/evolution/protocol_v9.py"
IMPLEMENTATION = (Path(__file__).resolve(), PROTOCOL_MODULE)


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


def _load_canonical_json_bytes(payload: bytes, artifact_id: str) -> dict[str, Any]:
    value = json.loads(payload)
    if not isinstance(value, dict) or value.get("artifact_id") != artifact_id:
        raise ValueError(f"{artifact_id} canonical object differs")
    if payload != _canonical(value).encode("utf-8"):
        raise ValueError(f"{artifact_id} JSON is not canonical")
    return value


def _load_canonical_json(path: Path, artifact_id: str) -> dict[str, Any]:
    return _load_canonical_json_bytes(path.read_bytes(), artifact_id)


def _git(*args: str, check: bool = True) -> subprocess.CompletedProcess[bytes]:
    return subprocess.run(
        ["git", *args],
        cwd=REPOSITORY,
        check=check,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )


def _checkpoint_blob(commit: str, relative: str, expected_sha: str) -> bytes:
    payload = _git("show", f"{commit}:{relative}").stdout
    if _bytes_sha(payload) != expected_sha:
        raise ValueError(f"immutable checkpoint blob differs: {relative}")
    return payload


def load_protocol(path: Path) -> dict[str, Any]:
    with path.open("rb") as handle:
        protocol = tomllib.load(handle)
    validate_sf1_protocol_v9(protocol)
    return protocol


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    with path.open("rb") as handle:
        config = tomllib.load(handle)
    scope = config.get("scope", {})
    lineage = config.get("immutable_lineage", {})
    proof = config.get("proof_contract", {})
    if (
        config.get("schema_version") != 1
        or config.get("artifact_id") != ARTIFACT_ID
        or scope.get("target_protocol") != "FGC-2-SF1-PROTO9"
        or scope.get("predecessor_protocol") != "FGC-2-SF1-PROTO8"
        or scope.get("diagnosis_artifact") != "FGC-1-CAL5-PREF7"
        or scope.get("freeze_role")
        != "premise_only_post_calibration_resolution_ladder_protocol_and_lineage_certificate"
        or scope.get("PROTO8_GR0_calibration_trajectory_disclosed") is not True
        or scope.get("FGCQR_evolution_outcomes_inspected") is not False
        or scope.get("SGBL_evolution_outcomes_inspected") is not False
        or scope.get("mechanism_question_answered") is not False
        or any(value is not True for value in proof.values())
        or config.get("claims") != PROTO9_CLAIMS
    ):
        raise ValueError("PRO9-FRZ1 config violates its fail-closed scope")
    commit = lineage.get("checkpoint_commit")
    if not isinstance(commit, str) or len(commit) != 40:
        raise ValueError("PRO9 checkpoint commit must be a Git SHA-1")
    for key in (
        "predecessor_protocol_sha256",
        "diagnosis_config_sha256",
        "diagnosis_result_sha256",
    ):
        value = lineage.get(key)
        if not isinstance(value, str) or len(value) != 64:
            raise ValueError(f"immutable_lineage.{key} must be SHA-256")
    return config


def _validate_lineage(config: Mapping[str, Any]) -> dict[str, Any]:
    lineage = config["immutable_lineage"]
    commit = lineage["checkpoint_commit"]
    _git("cat-file", "-e", f"{commit}^{{commit}}")
    if _git("merge-base", "--is-ancestor", commit, "HEAD", check=False).returncode != 0:
        raise ValueError("PRO9 checkpoint is not an ancestor of HEAD")

    predecessor_bytes = _checkpoint_blob(
        commit,
        lineage["predecessor_protocol_config"],
        lineage["predecessor_protocol_sha256"],
    )
    predecessor = tomllib.loads(predecessor_bytes.decode("utf-8"))
    predecessor_validation = validate_sf1_protocol_v8(predecessor)
    diagnosis_config_bytes = _checkpoint_blob(
        commit,
        lineage["diagnosis_config"],
        lineage["diagnosis_config_sha256"],
    )
    diagnosis_config = tomllib.loads(diagnosis_config_bytes.decode("utf-8"))
    if (
        diagnosis_config.get("artifact_id") != "FGC-1-CAL5-PREF7"
        or diagnosis_config.get("claims", {}).get(
            "PROTO9_resolution_ladder_revision_required"
        )
        is not True
        or diagnosis_config.get("claims", {}).get("PROTO9_frozen") is not False
    ):
        raise ValueError("immutable CAL5 config boundary differs")
    diagnosis_result_bytes = _checkpoint_blob(
        commit,
        lineage["diagnosis_result"],
        lineage["diagnosis_result_sha256"],
    )
    diagnosis = _load_canonical_json_bytes(
        diagnosis_result_bytes, "FGC-1-CAL5-PREF7"
    )
    payload = diagnosis.get("artifact_payload", {})
    decision = payload.get("decision", {})
    boundary = payload.get("epistemic_boundary", {})
    ladder = payload.get("resolution_ladder_diagnosis", {})
    if (
        diagnosis.get("classification")
        != "post_PROTO8_GR0_resolution_ladder_adequacy_diagnosis"
        or decision.get("new_protocol_required") != "FGC-2-SF1-PROTO9"
        or decision.get("PROTO8_GR0_case_eligible") is not False
        or decision.get("current_failure_is_a_candidate_or_gradient_obstruction")
        is not False
        or ladder.get("current_point_counts") != [1025, 2049, 4097]
        or ladder.get("smallest_prospective_point_counts")
        != [2049, 4097, 8193]
        or ladder.get("threshold_changes_required") is not False
        or ladder.get("new_8193_data_required_before_admission") is not True
        or any(value is not False for value in boundary.values())
        or any(value is not False for value in diagnosis.get("nonclaims", {}).values())
    ):
        raise ValueError("immutable CAL5 diagnosis or nonclaim boundary differs")
    return {
        "checkpoint_commit": commit,
        "checkpoint_is_ancestor_of_HEAD": True,
        "predecessor_protocol": predecessor_validation,
        "predecessor_protocol_sha256": lineage["predecessor_protocol_sha256"],
        "diagnosis": {
            "artifact_id": diagnosis["artifact_id"],
            "classification": diagnosis["classification"],
            "campaign_id": payload["immutable_campaign"]["campaign_id"],
            "terminal_classification": payload["immutable_campaign"][
                "terminal_classification"
            ],
            "selected_amplitude": payload["immutable_campaign"][
                "selected_amplitude"
            ],
            "current_point_counts": ladder["current_point_counts"],
            "prospective_point_counts": ladder[
                "smallest_prospective_point_counts"
            ],
            "new_8193_data_required": ladder[
                "new_8193_data_required_before_admission"
            ],
            "config_sha256": lineage["diagnosis_config_sha256"],
            "result_sha256": lineage["diagnosis_result_sha256"],
            "claims_all_candidate_and_physical_false": True,
        },
    }


def record(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    config = load_config(path)
    protocol_path = (REPOSITORY / config["protocol_config"]).resolve()
    protocol = load_protocol(protocol_path)
    validation = validate_sf1_protocol_v9(protocol)
    lineage = _validate_lineage(config)
    roots = (
        REPOSITORY / validation["calibration_output_root"],
        REPOSITORY / validation["holdout_output_root"],
    )
    namespace = [
        {"path": _rel(root), "exists": root.exists()} for root in roots
    ]
    if any(item["exists"] for item in namespace):
        raise ValueError("PROTO9 output namespace is not fresh at freeze time")
    return {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": config["project_version"],
        "metric_signature": config["metric_signature"],
        "riemann_convention": config["riemann_convention"],
        "classification": "outcome_neutral_PROTO9_resolution_ladder_protocol_freeze",
        "generated_by": _rel(Path(__file__)),
        "artifact_payload": {
            "immutable_lineage": lineage,
            "protocol_validation": validation,
            "resolution_revision": {
                "old_point_counts": [1025, 2049, 4097],
                "point_counts": [2049, 4097, 8193],
                "only_resolution_ladder_changes": True,
                "constraint_rules_unchanged": True,
                "spectral_rules_unchanged": True,
                "all_other_physical_and_numerical_rules_unchanged": True,
                "conditional_projection_is_admission": False,
                "new_8193_data_exist": False,
            },
            "namespace_precondition": {
                "records": namespace,
                "both_fresh_namespaces_absent": True,
                "freeze_created_no_namespace": True,
            },
            "decision": (
                "PROTO9_is_frozen_but_HLT7_fresh_calibration_and_every_"
                "candidate_outcome_remain_closed"
            ),
            "epistemic_boundary": {
                "PROTO9_trajectory_read": False,
                "PROTO9_trajectory_advanced": False,
                "new_8193_state_constructed": False,
                "GR0_calibration_completed": False,
                "SGBL_or_FGCQR_trajectory_read": False,
                "mechanism_question_answered": False,
            },
        },
        "gate_status": dict(config["claims"]),
        "implementation_sha256": {_rel(item): _sha(item) for item in IMPLEMENTATION},
        "derivation_document": _rel(OWNER_DOCUMENT),
        "derivation_document_sha256": _sha(OWNER_DOCUMENT),
        "source_config_sha256": {
            _rel(path): _sha(path),
            _rel(protocol_path): _sha(protocol_path),
        },
        "nonclaims": {
            "PROTO9_runtime_implemented": False,
            "fresh_GR0_calibration_authorized": False,
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
    return _load_canonical_json(path, ARTIFACT_ID)


def verify_canonical(
    path: Path = DEFAULT_OUTPUT,
    config_path: Path = DEFAULT_CONFIG,
) -> None:
    if load_canonical_result(path) != record(config_path):
        raise ValueError("stored PRO9-FRZ1 result differs from reconstruction")


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
        print(f"verified {_rel(output)} (PROTO9 frozen=true; runtime=false)")
        return
    payload = record(config)
    _atomic_write(output, _canonical(payload).encode("utf-8"))
    print(
        f"wrote {_rel(output)} "
        "(PROTO9 frozen=true; HLT7=false; trajectory_read=false)"
    )


if __name__ == "__main__":
    main()
