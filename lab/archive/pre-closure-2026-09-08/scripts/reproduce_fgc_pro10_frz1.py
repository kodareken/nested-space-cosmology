#!/usr/bin/env python3
"""Reproduce the immutable, outcome-neutral PROTO10 spectral freeze."""

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

from recursive_horizons.fgc.evolution.protocol_v9 import (  # noqa: E402
    validate_sf1_protocol_v9,
)
from recursive_horizons.fgc.evolution.protocol_v10 import (  # noqa: E402
    PROTO10_AMENDMENT,
    PROTO10_CLAIMS,
    validate_sf1_protocol_v10,
)


ARTIFACT_ID = "FGC-1-PRO10-FRZ1"
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-pro10-frz1.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-pro10-frz1.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-pro10-frz1.md"
PROTOCOL_MODULE = REPOSITORY / "src/recursive_horizons/fgc/evolution/protocol_v10.py"
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
    validate_sf1_protocol_v10(protocol)
    return protocol


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    with path.open("rb") as handle:
        config = tomllib.load(handle)
    if set(config) != {
        "schema_version",
        "artifact_id",
        "project_version",
        "metric_signature",
        "riemann_convention",
        "protocol_config",
        "scope",
        "immutable_lineage",
        "proof_contract",
        "claims",
    }:
        raise ValueError("PRO10-FRZ1 config root differs")
    scope = config["scope"]
    lineage = config["immutable_lineage"]
    proof = config["proof_contract"]
    if (
        config["schema_version"] != 1
        or config["artifact_id"] != ARTIFACT_ID
        or config["project_version"] != "0.11.0"
        or config["metric_signature"] != "-+++"
        or config["riemann_convention"] != "plus_partial_mu_gamma_nu"
        or config["protocol_config"]
        != "configs/fgc/fgc-2-sf1-protocol-v10.toml"
        or scope
        != {
            "target_protocol": "FGC-2-SF1-PROTO10",
            "predecessor_protocol": "FGC-2-SF1-PROTO9",
            "diagnosis_artifact": "FGC-1-CAL6-PREF8",
            "calibration_branch": "GR-0",
            "holdout_branch": "FGC-QR",
            "freeze_role": (
                "premise_only_finest_pair_spectral_interpretation_protocol_and_lineage_certificate"
            ),
            "PROTO9_GR0_trajectory_inspected": False,
            "FGCQR_evolution_outcomes_inspected": False,
            "SGBL_evolution_outcomes_inspected": False,
            "mechanism_question_answered": False,
        }
        or any(value is not True for value in proof.values())
        or config["claims"] != PROTO10_CLAIMS
    ):
        raise ValueError("PRO10-FRZ1 config violates its fail-closed scope")
    if set(lineage) != {
        "checkpoint_commit",
        "predecessor_protocol_config",
        "predecessor_protocol_sha256",
        "diagnosis_config",
        "diagnosis_config_sha256",
        "diagnosis_result",
        "diagnosis_result_sha256",
    }:
        raise ValueError("PRO10 immutable lineage keys differ")
    commit = lineage["checkpoint_commit"]
    if not isinstance(commit, str) or len(commit) != 40:
        raise ValueError("PRO10 checkpoint commit must be a Git SHA-1")
    for key in (
        "predecessor_protocol_sha256",
        "diagnosis_config_sha256",
        "diagnosis_result_sha256",
    ):
        value = lineage[key]
        if not isinstance(value, str) or len(value) != 64:
            raise ValueError(f"immutable_lineage.{key} must be SHA-256")
    amendment_projection = {
        "checkpoint_commit": PROTO10_AMENDMENT["diagnosis_checkpoint_commit"],
        "predecessor_protocol_config": PROTO10_AMENDMENT[
            "predecessor_protocol_config"
        ],
        "predecessor_protocol_sha256": PROTO10_AMENDMENT[
            "predecessor_protocol_sha256"
        ],
        "diagnosis_config": PROTO10_AMENDMENT["diagnosis_config"],
        "diagnosis_config_sha256": PROTO10_AMENDMENT["diagnosis_config_sha256"],
        "diagnosis_result": PROTO10_AMENDMENT["diagnosis_result"],
        "diagnosis_result_sha256": PROTO10_AMENDMENT["diagnosis_result_sha256"],
    }
    if dict(lineage) != amendment_projection:
        raise ValueError("PRO10 lineage differs from the frozen protocol amendment")
    return config


def _validate_lineage(config: Mapping[str, Any]) -> dict[str, Any]:
    lineage = config["immutable_lineage"]
    commit = lineage["checkpoint_commit"]
    _git("cat-file", "-e", f"{commit}^{{commit}}")
    if _git("merge-base", "--is-ancestor", commit, "HEAD", check=False).returncode != 0:
        raise ValueError("PRO10 checkpoint is not an ancestor of HEAD")

    predecessor_bytes = _checkpoint_blob(
        commit,
        lineage["predecessor_protocol_config"],
        lineage["predecessor_protocol_sha256"],
    )
    predecessor = tomllib.loads(predecessor_bytes.decode("utf-8"))
    predecessor_validation = validate_sf1_protocol_v9(predecessor)

    diagnosis_config_bytes = _checkpoint_blob(
        commit,
        lineage["diagnosis_config"],
        lineage["diagnosis_config_sha256"],
    )
    diagnosis_config = tomllib.loads(diagnosis_config_bytes.decode("utf-8"))
    prospective = diagnosis_config.get("prospective_revision", {})
    diagnosis_claims = diagnosis_config.get("claims", {})
    if (
        diagnosis_config.get("artifact_id") != "FGC-1-CAL6-PREF8"
        or diagnosis_config.get("scope", {}).get("target_protocol")
        != "FGC-2-SF1-PROTO9"
        or diagnosis_config.get("scope", {}).get("PROTO9_trajectory_read")
        is not False
        or prospective.get("new_protocol_required") != "FGC-2-SF1-PROTO10"
        or prospective.get("only_finest_pair_nested_tail_interpretation_may_change")
        is not True
        or prospective.get("resolution_ladder_unchanged") is not True
        or prospective.get("coarse_to_medium_direct_ratio_veto_unchanged")
        is not True
        or diagnosis_claims.get("PROTO10_spectral_contract_revision_required")
        is not True
        or diagnosis_claims.get("PROTO10_frozen") is not False
        or diagnosis_claims.get("PROTO10_fresh_GR0_dynamic_calibration_authorized")
        is not False
    ):
        raise ValueError("immutable CAL6/PREF8 config boundary differs")

    diagnosis_result_bytes = _checkpoint_blob(
        commit,
        lineage["diagnosis_result"],
        lineage["diagnosis_result_sha256"],
    )
    diagnosis = _load_canonical_json_bytes(
        diagnosis_result_bytes, "FGC-1-CAL6-PREF8"
    )
    payload = diagnosis.get("artifact_payload", {})
    aggregate = payload.get("aggregate", {})
    decision = payload.get("decision", {})
    boundary = payload.get("epistemic_boundary", {})
    if (
        diagnosis.get("classification")
        != "pre_trajectory_proper_spectral_ratio_conditioning_diagnosis"
        or diagnosis.get("source_config_sha256", {}).get(lineage["diagnosis_config"])
        != lineage["diagnosis_config_sha256"]
        or aggregate.get("case_count") != 4
        or aggregate.get("raw_failed_finest_pair_metric_count") != 16
        or aggregate.get("prospective_diagnostically_saturated_metric_count") != 16
        or aggregate.get("prospective_directly_resolved_finest_pair_metric_count")
        != 32
        or aggregate.get("all_raw_PROTO9_admissions_failed") is not True
        or aggregate.get("all_prospective_resolved_or_saturated_admissions_passed")
        is not True
        or decision.get("PROTO9_ratio_only_failure_is_candidate_or_gradient_evidence")
        is not False
        or decision.get("HLT7_finest_pair_ratio_is_conditioned_below_the_diagnostic_map_scale")
        is not True
        or decision.get("new_protocol_required") != "FGC-2-SF1-PROTO10"
        or decision.get("current_campaign_authorized") is not False
        or any(value is not False for value in boundary.values())
        or any(value is not False for value in diagnosis.get("nonclaims", {}).values())
    ):
        raise ValueError("immutable CAL6/PREF8 diagnosis or nonclaim boundary differs")
    return {
        "checkpoint_commit": commit,
        "checkpoint_is_ancestor_of_HEAD": True,
        "predecessor_protocol": predecessor_validation,
        "predecessor_protocol_sha256": lineage["predecessor_protocol_sha256"],
        "diagnosis": {
            "artifact_id": diagnosis["artifact_id"],
            "classification": diagnosis["classification"],
            "config_sha256": lineage["diagnosis_config_sha256"],
            "result_sha256": lineage["diagnosis_result_sha256"],
            "case_count": aggregate["case_count"],
            "raw_failed_finest_pair_metric_count": aggregate[
                "raw_failed_finest_pair_metric_count"
            ],
            "diagnostically_saturated_metric_count": aggregate[
                "prospective_diagnostically_saturated_metric_count"
            ],
            "directly_resolved_finest_pair_metric_count": aggregate[
                "prospective_directly_resolved_finest_pair_metric_count"
            ],
            "maximum_top_band_erasure_over_round_trip": aggregate[
                "maximum_top_band_erasure_over_round_trip"
            ],
            "maximum_complete_profile_finest_pair_difference_ratio": aggregate[
                "maximum_complete_profile_finest_pair_difference_ratio"
            ],
            "maximum_round_trip_adjacent_ratio": aggregate[
                "maximum_round_trip_adjacent_ratio"
            ],
            "all_prospective_design_admissions_passed": True,
            "current_campaign_authorized": False,
            "claims_candidate_and_physical_remain_false": True,
        },
    }


def record(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    config = load_config(path)
    protocol_path = (REPOSITORY / config["protocol_config"]).resolve()
    protocol = load_protocol(protocol_path)
    validation = validate_sf1_protocol_v10(protocol)
    lineage = _validate_lineage(config)
    roots = (
        REPOSITORY / validation["calibration_output_root"],
        REPOSITORY / validation["holdout_output_root"],
    )
    namespace = [{"path": _rel(root), "exists": root.exists()} for root in roots]
    if any(item["exists"] for item in namespace):
        raise ValueError("PROTO10 output namespace is not fresh at freeze time")
    diagnosis_summary = lineage["diagnosis"]
    return {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": config["project_version"],
        "metric_signature": config["metric_signature"],
        "riemann_convention": config["riemann_convention"],
        "classification": (
            "outcome_neutral_PROTO10_finest_pair_spectral_interpretation_freeze"
        ),
        "generated_by": _rel(Path(__file__)),
        "artifact_payload": {
            "immutable_lineage": lineage,
            "protocol_validation": validation,
            "spectral_revision": {
                "point_counts": [2049, 4097, 8193],
                "resolution_ladder_unchanged": True,
                "maximum_nested_tail_ratio": "1/4",
                "raw_ratio_ceiling_unchanged": True,
                "coarse_to_medium_direct_ratio_veto_unchanged": True,
                "medium_and_fine_absolute_budgets_unchanged": True,
                "finest_pair_direct_route_remains_primary": True,
                "finest_pair_guarded_saturation_route_added": True,
                "tail_erasure_map_witness_required": True,
                "round_trip_contraction_required": True,
                "complete_profile_infinity_and_RMS_contraction_required": True,
                "round_trip_is_a_continuum_error_bound": False,
                "diagnostic_saturation_is_physical_resolution_evidence": False,
                "new_epsilon_floor_or_numeric_tolerance_added": False,
                "immutable_design_witness": {
                    "raw_failed_finest_pair_metric_count": diagnosis_summary[
                        "raw_failed_finest_pair_metric_count"
                    ],
                    "diagnostically_saturated_metric_count": diagnosis_summary[
                        "diagnostically_saturated_metric_count"
                    ],
                    "directly_resolved_finest_pair_metric_count": diagnosis_summary[
                        "directly_resolved_finest_pair_metric_count"
                    ],
                    "maximum_top_band_erasure_over_round_trip": diagnosis_summary[
                        "maximum_top_band_erasure_over_round_trip"
                    ],
                    "maximum_complete_profile_finest_pair_difference_ratio": diagnosis_summary[
                        "maximum_complete_profile_finest_pair_difference_ratio"
                    ],
                    "maximum_round_trip_adjacent_ratio": diagnosis_summary[
                        "maximum_round_trip_adjacent_ratio"
                    ],
                    "is_calibration_or_mechanism_evidence": False,
                },
            },
            "namespace_precondition": {
                "records": namespace,
                "both_fresh_namespaces_absent": True,
                "freeze_created_no_namespace": True,
            },
            "decision": (
                "PROTO10_is_frozen_but_HLT8_fresh_calibration_and_every_"
                "candidate_outcome_remain_closed"
            ),
            "epistemic_boundary": {
                "PROTO10_trajectory_read": False,
                "PROTO10_trajectory_advanced": False,
                "fresh_state_constructed": False,
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
            "PROTO10_runtime_implemented": False,
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
        raise ValueError("stored PRO10-FRZ1 result differs from reconstruction")


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
        print(f"verified {_rel(output)} (PROTO10 frozen=true; runtime=false)")
        return
    payload = record(config)
    _atomic_write(output, _canonical(payload).encode("utf-8"))
    print(
        f"wrote {_rel(output)} "
        "(PROTO10 frozen=true; HLT8=false; trajectory_read=false)"
    )


if __name__ == "__main__":
    main()
