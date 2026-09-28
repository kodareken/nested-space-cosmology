#!/usr/bin/env python3
"""Reproduce the immutable, outcome-neutral PROTO11 numerical-map freeze."""

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

from recursive_horizons.fgc.evolution.protocol_v10 import (  # noqa: E402
    validate_sf1_protocol_v10,
)
from recursive_horizons.fgc.evolution.protocol_v11 import (  # noqa: E402
    PROTO11_AMENDMENT,
    PROTO11_CLAIMS,
    validate_sf1_protocol_v11,
)


ARTIFACT_ID = "FGC-1-PRO11-FRZ1"
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-pro11-frz1.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-pro11-frz1.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-pro11-frz1.md"
PROTOCOL_MODULE = REPOSITORY / "src/recursive_horizons/fgc/evolution/protocol_v11.py"
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
    validate_sf1_protocol_v11(protocol)
    return protocol


EXPECTED_SCOPE = {
    "target_protocol": "FGC-2-SF1-PROTO11",
    "predecessor_protocol": "FGC-2-SF1-PROTO10",
    "diagnosis_artifact": "FGC-1-CAL7-PREF9",
    "calibration_branch": "GR-0",
    "holdout_branch": "FGC-QR",
    "freeze_role": (
        "premise_only_well_balanced_reference_map_protocol_and_lineage_certificate"
    ),
    "PROTO10_GR0_calibration_trajectory_inspected": True,
    "PROTO10_candidate_or_holdout_outcome_observed": False,
    "collapse_or_trapped_outcome_classified": False,
    "FGCQR_evolution_outcomes_inspected": False,
    "SGBL_evolution_outcomes_inspected": False,
    "mechanism_question_answered": False,
}

EXPECTED_PROOF_KEYS = {
    "PROTO11_overlay_must_validate_exactly",
    "checkpoint_must_be_ancestor_of_HEAD",
    "predecessor_and_diagnosis_must_be_read_from_checkpoint",
    "immutable_blob_hashes_must_match_PROTO11_amendment",
    "immutable_PROTO10_must_validate_unchanged",
    "immutable_CAL7_must_be_canonical_and_hash_bound",
    "CAL7_must_remain_a_numerical_premise_diagnosis_not_a_mechanism_result",
    "only_semidiscrete_reference_state_differentiation_may_change",
    "reference_state_and_field_order_must_be_exact_and_fixed",
    "initial_q_q_r_and_reduction_rules_must_match_the_frozen_algebra",
    "interior_q_reprojection_damping_and_constraint_cleaning_are_forbidden",
    "raw_source_gate_and_complete_unredefined_source_must_remain_unchanged",
    "all_physical_inputs_equations_solvers_methods_CFL_retries_stops_spectra_observables_and_claims_must_remain_unchanged",
    "projected_input_hashes_must_be_rebuilt_under_PROTO11",
    "exact_Minkowski_and_nontrivial_REF1_controls_are_required_before_runtime",
    "new_output_namespaces_must_be_declared_but_not_accessed",
    "canonical_hash_bound_result_required",
}


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
        raise ValueError("PRO11-FRZ1 config root differs")
    lineage = config["immutable_lineage"]
    proof = config["proof_contract"]
    if (
        config["schema_version"] != 1
        or config["artifact_id"] != ARTIFACT_ID
        or config["project_version"] != "0.11.0"
        or config["metric_signature"] != "-+++"
        or config["riemann_convention"] != "plus_partial_mu_gamma_nu"
        or config["protocol_config"]
        != "configs/fgc/fgc-2-sf1-protocol-v11.toml"
        or config["scope"] != EXPECTED_SCOPE
        or set(proof) != EXPECTED_PROOF_KEYS
        or any(value is not True for value in proof.values())
        or config["claims"] != PROTO11_CLAIMS
    ):
        raise ValueError("PRO11-FRZ1 config violates its fail-closed scope")
    expected_lineage = {
        "checkpoint_commit": PROTO11_AMENDMENT["diagnosis_checkpoint_commit"],
        "predecessor_protocol_config": PROTO11_AMENDMENT[
            "predecessor_protocol_config"
        ],
        "predecessor_protocol_sha256": PROTO11_AMENDMENT[
            "predecessor_protocol_sha256"
        ],
        "diagnosis_config": PROTO11_AMENDMENT["diagnosis_config"],
        "diagnosis_config_sha256": PROTO11_AMENDMENT["diagnosis_config_sha256"],
        "diagnosis_result": PROTO11_AMENDMENT["diagnosis_result"],
        "diagnosis_result_sha256": PROTO11_AMENDMENT["diagnosis_result_sha256"],
    }
    if dict(lineage) != expected_lineage:
        raise ValueError("PRO11 lineage differs from the frozen protocol amendment")
    return config


def _validate_lineage(config: Mapping[str, Any]) -> dict[str, Any]:
    lineage = config["immutable_lineage"]
    commit = lineage["checkpoint_commit"]
    _git("cat-file", "-e", f"{commit}^{{commit}}")
    if _git("merge-base", "--is-ancestor", commit, "HEAD", check=False).returncode != 0:
        raise ValueError("PRO11 checkpoint is not an ancestor of HEAD")

    predecessor_bytes = _checkpoint_blob(
        commit,
        lineage["predecessor_protocol_config"],
        lineage["predecessor_protocol_sha256"],
    )
    predecessor = tomllib.loads(predecessor_bytes.decode("utf-8"))
    predecessor_validation = validate_sf1_protocol_v10(predecessor)

    diagnosis_config_bytes = _checkpoint_blob(
        commit,
        lineage["diagnosis_config"],
        lineage["diagnosis_config_sha256"],
    )
    diagnosis_config = tomllib.loads(diagnosis_config_bytes.decode("utf-8"))
    prospective = diagnosis_config.get("prospective_revision", {})
    diagnosis_claims = diagnosis_config.get("claims", {})
    if (
        diagnosis_config.get("artifact_id") != "FGC-1-CAL7-PREF9"
        or diagnosis_config.get("scope", {}).get("target_protocol")
        != "FGC-2-SF1-PROTO10"
        or diagnosis_config.get("scope", {}).get("GR0_calibration_trajectory_read")
        is not True
        or diagnosis_config.get("scope", {}).get("SGBL_trajectory_read")
        is not False
        or diagnosis_config.get("scope", {}).get("FGCQR_trajectory_read")
        is not False
        or prospective.get("new_protocol_required") != "FGC-2-SF1-PROTO11"
        or prospective.get("only_semidiscrete_reference_state_differentiation_may_change")
        is not True
        or prospective.get("raw_source_residual_and_tolerance_remain_public_and_unchanged")
        is not True
        or prospective.get("fresh_campaign_required") is not True
        or diagnosis_claims.get("PROTO11_well_balanced_reference_map_required")
        is not True
        or diagnosis_claims.get("PROTO11_frozen") is not False
        or diagnosis_claims.get("PROTO10_GR0_case_eligible") is not False
    ):
        raise ValueError("immutable CAL7/PREF9 config boundary differs")

    diagnosis_result_bytes = _checkpoint_blob(
        commit,
        lineage["diagnosis_result"],
        lineage["diagnosis_result_sha256"],
    )
    diagnosis = _load_canonical_json_bytes(
        diagnosis_result_bytes, "FGC-1-CAL7-PREF9"
    )
    payload = diagnosis.get("artifact_payload", {})
    campaign = payload.get("campaign_diagnosis", {})
    comparison = payload.get("well_balanced_control_summary", {}).get(
        "comparison", {}
    )
    raw = payload.get("raw_initial_source_floor", [])
    controlled = payload.get("well_balanced_initial_control", [])
    boundary = payload.get("epistemic_boundary", {})
    if (
        diagnosis.get("classification")
        != "post_calibration_PROTO10_center_adjacent_binary64_source_floor_diagnosis"
        or diagnosis.get("source_config_sha256", {}).get(lineage["diagnosis_config"])
        != lineage["diagnosis_config_sha256"]
        or campaign.get("event_count") != 264
        or campaign.get("common_event_count") != 2
        or campaign.get("source_only_rejection_count") != 262
        or campaign.get("no_non_source_failure_observed") is not True
        or campaign.get("all_rejections_exactly_rolled_back") is not True
        or len(raw) != 12
        or len(controlled) != 12
        or comparison.get("all_well_balanced_source_residuals_pass_unchanged_raw_tolerance")
        is not True
        or comparison.get("maximum_well_balanced_source_residual")
        != 7.579122514774399e-14
        or comparison.get("fine_primary_residual_reduction_factor")
        != 10.375000000000075
        or comparison.get("well_balanced_control_is_not_a_runtime_or_campaign_result")
        is not True
        or boundary.get("repair_status")
        != "well-balanced t0 control passes, but PROTO11 is not frozen or run"
        or diagnosis.get("gate_status", {}).get("PROTO11_frozen") is not False
        or diagnosis.get("gate_status", {}).get("FGCQR_holdout_execution_authorized")
        is not False
        or any(value is not False for value in diagnosis.get("nonclaims", {}).values())
    ):
        raise ValueError("immutable CAL7/PREF9 diagnosis or nonclaim boundary differs")
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
            "campaign_event_count": campaign["event_count"],
            "source_only_rejection_count": campaign["source_only_rejection_count"],
            "raw_initial_record_count": len(raw),
            "well_balanced_control_record_count": len(controlled),
            "maximum_well_balanced_source_residual": comparison[
                "maximum_well_balanced_source_residual"
            ],
            "fine_primary_residual_reduction_factor": comparison[
                "fine_primary_residual_reduction_factor"
            ],
            "is_runtime_calibration_or_mechanism_evidence": False,
        },
    }


def record(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    config = load_config(path)
    protocol_path = (REPOSITORY / config["protocol_config"]).resolve()
    protocol = load_protocol(protocol_path)
    validation = validate_sf1_protocol_v11(protocol)
    lineage = _validate_lineage(config)
    roots = (
        REPOSITORY / validation["calibration_output_root"],
        REPOSITORY / validation["holdout_output_root"],
    )
    namespace = [{"path": _rel(root), "exists": root.exists()} for root in roots]
    if any(item["exists"] for item in namespace):
        raise ValueError("PROTO11 output namespace is not fresh at freeze time")
    diagnosis = lineage["diagnosis"]
    return {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": config["project_version"],
        "metric_signature": config["metric_signature"],
        "riemann_convention": config["riemann_convention"],
        "classification": (
            "outcome_neutral_PROTO11_well_balanced_reference_map_freeze"
        ),
        "generated_by": _rel(Path(__file__)),
        "scope_bindings": dict(config["scope"]),
        "artifact_payload": {
            "immutable_lineage": lineage,
            "protocol_validation": validation,
            "numerical_map_revision": {
                "field_order": ["alpha", "shift", "lambda", "R", "phi", "chi"],
                "reference_u": ["1", "0", "1", "r", "0", "0"],
                "reference_q": ["0", "0", "0", "1", "0", "0"],
                "initial_q_rule": "D_h(u-u_ref)+q_ref",
                "source_q_r_rule": "D_h(q-q_ref)",
                "reduction_constraint_rule": "q-(D_h(u-u_ref)+q_ref)",
                "p_r_du_dt_and_dq_dt_unchanged": True,
                "q_remains_independently_evolved": True,
                "interior_q_reprojection_added": False,
                "damping_or_constraint_cleaning_added": False,
                "continuum_equation_or_source_row_changed": False,
                "raw_source_residual_maximum": "1/1000000000000",
                "new_epsilon_floor_or_numeric_tolerance_added": False,
                "projected_input_hashes_must_be_rebuilt": True,
                "exact_Minkowski_and_generic_REF1_controls_required": True,
                "immutable_design_witness": {
                    "campaign_event_count": diagnosis["campaign_event_count"],
                    "source_only_rejection_count": diagnosis[
                        "source_only_rejection_count"
                    ],
                    "raw_initial_record_count": diagnosis[
                        "raw_initial_record_count"
                    ],
                    "well_balanced_control_record_count": diagnosis[
                        "well_balanced_control_record_count"
                    ],
                    "maximum_well_balanced_source_residual": diagnosis[
                        "maximum_well_balanced_source_residual"
                    ],
                    "fine_primary_residual_reduction_factor": diagnosis[
                        "fine_primary_residual_reduction_factor"
                    ],
                    "is_runtime_calibration_or_mechanism_evidence": False,
                },
            },
            "namespace_precondition": {
                "records": namespace,
                "both_fresh_namespaces_absent": True,
                "freeze_created_no_namespace": True,
            },
            "decision": (
                "PROTO11_is_frozen_but_HLT9_runtime_fresh_calibration_and_every_"
                "candidate_outcome_remain_closed"
            ),
            "epistemic_boundary": {
                "PROTO11_runtime_implemented": False,
                "PROTO11_trajectory_read": False,
                "PROTO11_trajectory_advanced": False,
                "fresh_state_constructed_under_PROTO11": False,
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
            "PROTO11_runtime_implemented": False,
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
        raise ValueError("stored PRO11-FRZ1 result differs from reconstruction")


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
        print(f"verified {_rel(output)} (PROTO11 frozen=true; runtime=false)")
        return
    payload = record(config)
    _atomic_write(output, _canonical(payload).encode("utf-8"))
    print(
        f"wrote {_rel(output)} "
        "(PROTO11 frozen=true; HLT9=false; trajectory_read=false)"
    )


if __name__ == "__main__":
    main()
