#!/usr/bin/env python3
"""Regenerate the immutable, outcome-neutral PROTO8 freeze certificate."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
import tomllib
from typing import Any, Mapping


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-pro8-frz1.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-pro8-frz1.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-pro8-frz1.md"

from recursive_horizons.fgc.evolution.protocol_v7 import (  # noqa: E402
    SF1_PROTOCOL_V7_ARTIFACT_ID,
    validate_sf1_protocol_v7,
)
from recursive_horizons.fgc.evolution.protocol_v8 import (  # noqa: E402
    PROTO8_AMENDMENT,
    PROTO8_CLAIMS,
    SF1_PROTOCOL_V8_ARTIFACT_ID,
    validate_sf1_protocol_v8,
)


ARTIFACT_ID = "FGC-1-PRO8-FRZ1"
PROJECT_VERSION = "0.11.0"
EXPECTED_SCOPE = {
    "target_protocol": SF1_PROTOCOL_V8_ARTIFACT_ID,
    "predecessor_protocol": SF1_PROTOCOL_V7_ARTIFACT_ID,
    "diagnosis_artifact": "FGC-1-CAL4-PREF6",
    "calibration_branch": "GR-0",
    "holdout_branch": "FGC-QR",
    "freeze_role": "premise_only_post_calibration_common_event_diagnostic_ownership_protocol_and_lineage_certificate",
    "PROTO7_GR0_calibration_trajectory_disclosed": True,
    "first_nonzero_common_event_completed": True,
    "FGCQR_evolution_outcomes_inspected": False,
    "SGBL_evolution_outcomes_inspected": False,
    "mechanism_question_answered": False,
}
EXPECTED_LINEAGE = {
    "checkpoint_commit": "219ba9493899029bc51bd313e3c5f689f3809c2c",
    "predecessor_protocol_config": "configs/fgc/fgc-2-sf1-protocol-v7.toml",
    "predecessor_protocol_sha256": "ea0563b205bf5c6c90d8f728afc56a13f624953fffc67c6d6b9a97c03c35f3fb",
    "diagnosis_config": "configs/fgc/fgc-1-cal4-pref6.toml",
    "diagnosis_config_sha256": "687322014f16a5c48f6b3947f4cfd18a0f60f603f3c4a3abda842b19e3c8f1e2",
    "diagnosis_result": "results/fgc-1-cal4-pref6.json",
    "diagnosis_result_sha256": "e03d647e9d068b35be79ba817a3ab6c36b93f29c29de5e395a2106030d958815",
}
EXPECTED_PROOF = {
    "PROTO8_overlay_must_validate_exactly",
    "checkpoint_must_be_ancestor_of_HEAD",
    "predecessor_and_diagnosis_must_be_read_from_checkpoint",
    "immutable_blob_hashes_must_match_PROTO8_amendment",
    "immutable_PROTO7_must_validate_unchanged",
    "immutable_CAL4_must_be_canonical_and_nonmechanism",
    "outer_constraint_ownership_must_match_projector_and_method_stencil",
    "accumulated_roundoff_must_change_order_classification_only",
    "raw_constraint_magnitude_guards_must_remain_unchanged",
    "coarsest_spectrum_must_remain_a_public_convergence_witness",
    "medium_and_finest_spectra_must_pass_every_unchanged_absolute_budget",
    "every_adjacent_nested_tail_must_still_contract",
    "all_physical_inputs_equations_solvers_thresholds_and_stops_must_remain_unchanged",
    "new_output_namespaces_must_be_declared_but_not_accessed",
    "canonical_hash_bound_result_required",
}
IMPLEMENTATION = (
    REPOSITORY / "src/recursive_horizons/fgc/evolution/protocol_v7.py",
    REPOSITORY / "src/recursive_horizons/fgc/evolution/protocol_v8.py",
    Path(__file__).resolve(),
)


def _canonical(value: object) -> str:
    return (
        json.dumps(value, sort_keys=True, indent=2, ensure_ascii=True, allow_nan=False)
        + "\n"
    )


def _unique_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    answer: dict[str, Any] = {}
    for key, value in pairs:
        if key in answer:
            raise ValueError(f"duplicate JSON key: {key}")
        answer[key] = value
    return answer


def _load_toml_bytes(source: bytes, name: str) -> dict[str, Any]:
    try:
        value = tomllib.loads(source.decode("utf-8"))
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as exc:
        raise ValueError(f"{name} is not valid UTF-8 TOML") from exc
    if not isinstance(value, dict):
        raise ValueError(f"{name} root must be a table")
    return value


def _load_json_bytes(source: bytes, name: str) -> dict[str, Any]:
    try:
        value = json.loads(source.decode("utf-8"), object_pairs_hook=_unique_pairs)
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"{name} is not valid unique-key JSON") from exc
    if not isinstance(value, dict) or source.decode("utf-8") != _canonical(value):
        raise ValueError(f"{name} is not canonical sorted JSON")
    return value


def _repo_file(name: str, value: object) -> Path:
    if not isinstance(value, str) or not value:
        raise TypeError(f"{name} must be a repository-relative path")
    path = (REPOSITORY / value).resolve()
    try:
        relative = path.relative_to(REPOSITORY).as_posix()
    except ValueError as exc:
        raise ValueError(f"{name} escapes the repository") from exc
    if relative != value or not path.is_file():
        raise ValueError(f"{name} is absent or noncanonical")
    return path


def _git_blob(commit: str, relative: str, name: str) -> bytes:
    completed = subprocess.run(
        ["git", "show", f"{commit}:{relative}"],
        cwd=REPOSITORY,
        check=False,
        capture_output=True,
    )
    if completed.returncode != 0:
        raise ValueError(f"cannot read immutable {name} blob from {commit}")
    return completed.stdout


def load_canonical_result(path: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    return _load_json_bytes(path.read_bytes(), "PRO8-FRZ1 result")


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    source = path.read_bytes()
    raw = _load_toml_bytes(source, "PRO8-FRZ1 config")
    expected_root = {
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
    }
    if set(raw) != expected_root:
        raise ValueError("PRO8-FRZ1 root keys differ")
    if (
        raw["schema_version"] != 1
        or raw["artifact_id"] != ARTIFACT_ID
        or raw["project_version"] != PROJECT_VERSION
        or raw["metric_signature"] != "-+++"
        or raw["riemann_convention"] != "plus_partial_mu_gamma_nu"
    ):
        raise ValueError("PRO8-FRZ1 identity or convention differs")
    if raw["scope"] != EXPECTED_SCOPE:
        raise ValueError("PRO8-FRZ1 scope differs")
    if raw["immutable_lineage"] != EXPECTED_LINEAGE:
        raise ValueError("PRO8-FRZ1 immutable lineage differs")
    if set(raw["proof_contract"]) != EXPECTED_PROOF or any(
        value is not True for value in raw["proof_contract"].values()
    ):
        raise ValueError("PRO8-FRZ1 proof contract differs")
    if raw["claims"] != PROTO8_CLAIMS:
        raise ValueError("PRO8-FRZ1 claims must remain fail-closed")

    protocol_path = _repo_file("protocol_config", raw["protocol_config"])
    protocol_source = protocol_path.read_bytes()
    protocol_raw = _load_toml_bytes(protocol_source, "PROTO8")
    protocol_certificate = validate_sf1_protocol_v8(protocol_raw)
    if protocol_raw["amendment"] != PROTO8_AMENDMENT:
        raise ValueError("PROTO8 amendment differs from its validator")
    amendment_lineage = {
        "checkpoint_commit": protocol_raw["amendment"]["diagnosis_checkpoint_commit"],
        "predecessor_protocol_config": protocol_raw["amendment"][
            "predecessor_protocol_config"
        ],
        "predecessor_protocol_sha256": protocol_raw["amendment"][
            "predecessor_protocol_sha256"
        ],
        "diagnosis_config": protocol_raw["amendment"]["diagnosis_config"],
        "diagnosis_config_sha256": protocol_raw["amendment"][
            "diagnosis_config_sha256"
        ],
        "diagnosis_result": protocol_raw["amendment"]["diagnosis_result"],
        "diagnosis_result_sha256": protocol_raw["amendment"][
            "diagnosis_result_sha256"
        ],
    }
    if amendment_lineage != EXPECTED_LINEAGE:
        raise ValueError("PROTO8 amendment and PRO8-FRZ1 lineage differ")
    return {
        "raw": raw,
        "source_sha256": sha256(source).hexdigest(),
        "protocol_path": protocol_path,
        "protocol_source": protocol_source,
        "protocol_certificate": protocol_certificate,
    }


def _validate_immutable_cal4(diagnosis: Mapping[str, Any], lineage: Mapping[str, Any]) -> dict[str, Any]:
    gates = diagnosis.get("gate_status", {})
    payload = diagnosis.get("artifact_payload", {})
    campaign = payload.get("immutable_campaign", {})
    common = payload.get("common_event_ownership_diagnosis", {})
    amp5 = common.get("amplitude_5over2", {})
    amp3 = common.get("amplitude_3", {})
    decision = payload.get("decision", {})
    boundary = payload.get("epistemic_boundary", {})
    amp5_methods = amp5.get("methods", {})
    amp3_methods = amp3.get("methods", {})
    if (
        diagnosis.get("artifact_id") != "FGC-1-CAL4-PREF6"
        or gates.get("PROTO7_common_event_ownership_contract_obstructed") is not True
        or gates.get("PROTO8_common_event_admission_revision_required") is not True
        or gates.get("fresh_GR0_dynamic_calibration_completed") is not False
        or gates.get("FGCQR_holdout_execution_authorized") is not False
        or campaign.get("first_nonzero_common_event_completed_by_both_amplitudes")
        is not True
        or campaign.get("terminal_classification")
        != "calibration_failed_no_eligible_GR0_case"
        or amp5.get("both_methods_prospectively_admitted") is not True
        or amp3.get("both_methods_prospectively_admitted") is not False
        or set(amp5_methods) != {"RK4", "SSPRK3"}
        or set(amp3_methods) != {"RK4", "SSPRK3"}
        or any(
            item.get("prospective_common_event_admission_passed") is not True
            or item.get("old_full_domain_reduction_failures_boundary_localized")
            is not True
            or item.get(
                "owned_reduction_residuals_within_accumulated_roundoff_enclosure"
            )
            is not True
            or item.get("spectral", {}).get(
                "finest_pair_individual_budgets_passed"
            )
            is not True
            or item.get("spectral", {}).get("every_nested_tail_ratio_passed")
            is not True
            for item in amp5_methods.values()
        )
        or decision.get("new_protocol_required") != SF1_PROTOCOL_V8_ARTIFACT_ID
        or boundary.get("PROTO7_transaction_repair_succeeded_through_t1") is not True
        or boundary.get("GR0_calibration_completed") is not False
        or boundary.get("SGBL_or_FGCQR_trajectory_read") is not False
        or boundary.get("mechanism_question_answered") is not False
        or boundary.get("candidate_or_general_gradient_route_rejected") is not False
        or any(value is not False for value in diagnosis.get("nonclaims", {}).values())
        or not diagnosis.get("nonclaims")
    ):
        raise ValueError("immutable CAL4 is not the bounded nonmechanism diagnosis")
    if (
        diagnosis.get("source_config_sha256", {}).get(lineage["diagnosis_config"])
        != lineage["diagnosis_config_sha256"]
    ):
        raise ValueError("immutable CAL4 does not bind its diagnosis config")
    return {
        "artifact_id": diagnosis["artifact_id"],
        "config_blob_sha256": lineage["diagnosis_config_sha256"],
        "result_blob_sha256": lineage["diagnosis_result_sha256"],
        "PROTO7_common_event_ownership_obstructed": True,
        "PROTO8_common_event_admission_revision_required": True,
        "first_nonzero_common_event_completed_by_both_amplitudes": True,
        "amplitude_5over2_prospectively_admitted": True,
        "amplitude_3_prospectively_admitted": False,
        "FGCQR_outcome_read": False,
        "mechanism_question_answered": False,
    }


def _immutable_lineage(configuration: Mapping[str, Any]) -> dict[str, Any]:
    lineage = configuration["raw"]["immutable_lineage"]
    commit = lineage["checkpoint_commit"]
    if subprocess.run(
        ["git", "merge-base", "--is-ancestor", commit, "HEAD"],
        cwd=REPOSITORY,
        check=False,
        capture_output=True,
    ).returncode != 0:
        raise ValueError("PROTO8 diagnosis checkpoint is not an ancestor of HEAD")

    predecessor_blob = _git_blob(
        commit, lineage["predecessor_protocol_config"], "PROTO7 protocol"
    )
    if sha256(predecessor_blob).hexdigest() != lineage["predecessor_protocol_sha256"]:
        raise ValueError("immutable PROTO7 blob hash differs")
    predecessor = validate_sf1_protocol_v7(
        _load_toml_bytes(predecessor_blob, "immutable PROTO7")
    )
    if (
        predecessor["artifact_id"] != SF1_PROTOCOL_V7_ARTIFACT_ID
        or predecessor["protocol_version"] != 7
        or predecessor["frozen"] is not True
        or predecessor["outcome_neutral_contract_validated"] is not True
        or any(value is not False for value in predecessor["claims"].values())
    ):
        raise ValueError("immutable PROTO7 is not frozen and fail-closed")

    diagnosis_config_blob = _git_blob(
        commit, lineage["diagnosis_config"], "CAL4 config"
    )
    diagnosis_blob = _git_blob(commit, lineage["diagnosis_result"], "CAL4 diagnosis")
    if sha256(diagnosis_config_blob).hexdigest() != lineage["diagnosis_config_sha256"]:
        raise ValueError("immutable CAL4 config hash differs")
    if sha256(diagnosis_blob).hexdigest() != lineage["diagnosis_result_sha256"]:
        raise ValueError("immutable CAL4 result hash differs")
    diagnosis = _load_json_bytes(diagnosis_blob, "immutable CAL4")
    diagnosis_certificate = _validate_immutable_cal4(diagnosis, lineage)

    for relative, expected, name in (
        (
            lineage["predecessor_protocol_config"],
            lineage["predecessor_protocol_sha256"],
            "current PROTO7",
        ),
        (
            lineage["diagnosis_config"],
            lineage["diagnosis_config_sha256"],
            "current CAL4 config",
        ),
        (
            lineage["diagnosis_result"],
            lineage["diagnosis_result_sha256"],
            "current CAL4 result",
        ),
    ):
        if sha256((REPOSITORY / relative).read_bytes()).hexdigest() != expected:
            raise ValueError(f"{name} differs from its immutable checkpoint blob")

    return {
        "checkpoint_commit": commit,
        "checkpoint_is_ancestor_of_HEAD": True,
        "predecessor_protocol": {
            "artifact_id": predecessor["artifact_id"],
            "protocol_version": predecessor["protocol_version"],
            "blob_sha256": lineage["predecessor_protocol_sha256"],
            "semantic_holdout_contract_sha256": predecessor[
                "semantic_holdout_contract_sha256"
            ],
            "claims_all_false": all(
                value is False for value in predecessor["claims"].values()
            ),
        },
        "diagnosis": diagnosis_certificate,
    }


def record(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    configuration = load_config(path)
    lineage = _immutable_lineage(configuration)
    protocol = configuration["protocol_certificate"]
    raw = configuration["raw"]
    protocol_relative = configuration["protocol_path"].relative_to(REPOSITORY).as_posix()
    gate_status = {
        "PROTO8_outcome_neutral_protocol_frozen": True,
        "PROTO8_immutable_lineage_verified": True,
        "PROTO8_evolution_owned_constraint_contract_frozen": True,
        "PROTO8_accumulated_roundoff_order_contract_frozen": True,
        "PROTO8_finest_pair_nested_spectral_contract_frozen": True,
        "PROTO8_physical_and_numerical_thresholds_inherited_unchanged": True,
        **raw["claims"],
    }
    return {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "classification": "immutable_outcome_neutral_PROTO8_common_event_diagnostic_ownership_protocol_freeze_not_run_authorization",
        "generated_by": "scripts/reproduce_fgc_pro8_frz1.py",
        "scope_bindings": dict(raw["scope"]),
        "source_config_sha256": {
            path.relative_to(REPOSITORY).as_posix(): configuration["source_sha256"],
            protocol_relative: sha256(configuration["protocol_source"]).hexdigest(),
        },
        "predecessor_sha256": {
            raw["immutable_lineage"]["predecessor_protocol_config"]: raw[
                "immutable_lineage"
            ]["predecessor_protocol_sha256"],
            raw["immutable_lineage"]["diagnosis_config"]: raw["immutable_lineage"][
                "diagnosis_config_sha256"
            ],
            raw["immutable_lineage"]["diagnosis_result"]: raw["immutable_lineage"][
                "diagnosis_result_sha256"
            ],
        },
        "implementation_sha256": {
            item.relative_to(REPOSITORY).as_posix(): sha256(item.read_bytes()).hexdigest()
            for item in IMPLEMENTATION
        },
        "derivation_document": OWNER_DOCUMENT.relative_to(REPOSITORY).as_posix(),
        "derivation_document_sha256": sha256(OWNER_DOCUMENT.read_bytes()).hexdigest(),
        "certificate_contract": sorted(raw["proof_contract"]),
        "gate_status": gate_status,
        "artifact_payload": {
            "protocol_certificate": protocol,
            "immutable_lineage": lineage,
            "frozen_repair": {
                "projector_fixed_outer_rows": 4,
                "RK4_excluded_outer_rows": 7,
                "SSPRK3_excluded_outer_rows": 5,
                "full_and_owned_constraint_residuals_are_public": True,
                "owned_raw_constraint_norms_decide_unchanged_magnitude_guards": True,
                "accumulated_roundoff_changes_order_classification_only": True,
                "accumulated_roundoff_formula": "4096*binary64_epsilon*(1+accepted_stage_count)",
                "minimum_finest_pair_constraint_order": "3/2",
                "coarsest_spectrum_is_public_convergence_witness": True,
                "medium_and_finest_absolute_spectral_budgets_required": True,
                "all_adjacent_nested_tail_ratios_required": True,
                "maximum_nested_tail_ratio": "1/4",
                "all_physical_and_numerical_thresholds_are_unchanged": True,
                "eligible_amplitudes_in_order": ["5/2", "3"],
            },
            "epistemic_boundary": {
                "PROTO7_GR0_calibration_trajectory_was_seen": True,
                "first_nonzero_common_event_completed": True,
                "PROTO8_trajectory_read": False,
                "FGCQR_outcome_read": False,
                "SGBL_outcome_read": False,
                "mechanism_question_answered": False,
                "physical_experiment_changed": False,
                "numerical_threshold_changed": False,
            },
            "decision": "PROTO8_is_frozen_but_HLT6_fresh_calibration_and_every_candidate_outcome_remain_closed",
        },
        "nonclaims": {
            "fresh_GR0_calibration_completed": False,
            "SGBL_comparison_completed": False,
            "FGCQR_holdout_evolved": False,
            "trapped_sphere_observed": False,
            "regulator_activation_observed": False,
            "positive_Raychaudhuri_margin_observed": False,
            "FGCQR_mechanism_rejected": False,
            "general_gradient_route_rejected": False,
            "retained_EFT_validity": False,
            "physical_transition": False,
        },
    }


def write(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(_canonical(payload), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    payload = record(args.config.resolve())
    write(args.output.resolve(), payload)
    print(
        f"wrote {args.output} "
        "(PROTO8 frozen=true; fresh calibration authorized=false; trajectory_read=false)"
    )


if __name__ == "__main__":
    main()
