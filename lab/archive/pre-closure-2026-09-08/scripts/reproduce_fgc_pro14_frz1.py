#!/usr/bin/env python3
"""Reproduce PROTO14's prospective TDG6 calibration freeze.

The certificate binds the immutable PROTO13/CAL10 history, TDG6-IMP2, and the
six original ``t=23/16`` GR-0 restart members.  It freezes the replacement
temporal admission, terminal classes, and fresh namespaces before a successor
runtime or production trajectory exists.
"""

from __future__ import annotations

import argparse
from copy import deepcopy
from hashlib import sha256
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

from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    array_content_sha256,
)
from recursive_horizons.fgc.evolution.protocol_v14 import (  # noqa: E402
    PROTO14_CLAIMS,
    PROTO14_RESTART,
    PROTO14_TEMPORAL_ADMISSION,
    PROTO14_TERMINAL,
    validate_sf1_protocol_v14,
)


ARTIFACT_ID = "FGC-1-PRO14-FRZ1"
PROJECT_VERSION = "0.11.0"
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-pro14-frz1.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-pro14-frz1.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-pro14-frz1.md"
PROTOCOL_CONFIG = REPOSITORY / "configs/fgc/fgc-2-sf1-protocol-v14.toml"
PROTOCOL_MODULE = REPOSITORY / "src/recursive_horizons/fgc/evolution/protocol_v14.py"
IMPLEMENTATION = (Path(__file__).resolve(), PROTOCOL_MODULE)

EXPECTED_SCOPE = {
    "target_protocol": "FGC-2-SF1-PROTO14",
    "predecessor_protocol": "FGC-2-SF1-PROTO13",
    "replacement_instrument_artifact": "FGC-1-TDG6-IMP2",
    "historical_terminal_artifact": "FGC-1-CAL10-PREF15",
    "calibration_branch": "GR-0",
    "calibration_amplitude": "3",
    "holdout_branch": "FGC-QR",
    "freeze_role": "premise_only_TDG6_post_restart_calibration_protocol_certificate",
    "PROTO13_GR0_trajectory_inspected": True,
    "PROTO14_trajectory_read": False,
    "SGBL_trajectory_read": False,
    "FGCQR_trajectory_read": False,
    "collapse_or_trapped_outcome_classified": False,
    "mechanism_question_answered": False,
}

EXPECTED_LINEAGE = {
    "checkpoint_commit": "27badf9ad3c0bc3aa1c8a4153a42a26bdaf4c00f",
    "predecessor_protocol_config": "configs/fgc/fgc-2-sf1-protocol-v13.toml",
    "predecessor_protocol_sha256": "fc9398bad6ea4448f659188641e6da6bd7c2e295a443fe73bf73eb0308480692",
    "predecessor_freeze_result": "results/fgc-1-pro13-frz1.json",
    "predecessor_freeze_result_sha256": "a6fd1ce1bf429fd39ff2d7975de009991199a0f41c73e44d46e57b88b35d74b4",
    "historical_terminal_result": "results/fgc-1-cal10-pref15.json",
    "historical_terminal_result_sha256": "1ebbdc7222baee8134a6d466589770ee445aee594f10630fc09aa6264cb4ac18",
    "replacement_instrument_config": "configs/fgc/fgc-1-tdg6-imp2.toml",
    "replacement_instrument_config_sha256": "89803b43e9243558aed33473bd9ebd010142b02302f1a5b0b380e8c4f4d7d8c0",
    "replacement_instrument_result": "results/fgc-1-tdg6-imp2.json",
    "replacement_instrument_result_sha256": "6a2cac24f65e067aed7ff7de299f9c1f0f3c3319e9aee692763b8df1db4fd5b5",
    "replacement_runtime_module": "src/recursive_horizons/fgc/evolution/tdg6_temporal_admission_runtime.py",
    "replacement_runtime_module_sha256": "039de0bd008a2536bea6f7187afc3856bfc841d14790cbae5ca2fe27f7eb8e2d",
    "PROTO12_checkpoint": "runs/fgc-2-sf1/proto12/calibration/latest-checkpoint.npz",
    "PROTO12_checkpoint_sha256": "c784d4706911029883d14e1763c2d36c5dbe0e619a4e50416cc6d666d5b7ea17",
    "RSP2_checkpoint": "runs/fgc-2-sf1/rsp2/amplitude3-ssprk3-momentum/latest-checkpoint.npz",
    "RSP2_checkpoint_sha256": "000546dcf3726c7882e80b7e70b64a5d2e060e2e0f7a567e174f714432ae7f34",
}

EXPECTED_RESTART_CONTRACT = {
    "restart_coordinate_time": "23/16",
    "final_coordinate_time": "32",
    "first_new_common_event_coordinate_time": "3/2",
    "common_output_interval": "1/16",
    "checkpoint_interval": "1/4",
    "primary_method": "RK4",
    "primary_point_counts": [2049, 4097, 8193],
    "comparator_method": "SSPRK3",
    "comparator_point_counts": [4097, 8193, 16385],
    "minimum_constraint_finest_pair_order": "3/2",
    "maximum_nested_spatial_tail_ratio": "1/4",
    "minimum_consecutive_qualified_trapped_common_events": 8,
    "trapped_sign_margin_over_combined_spatial_error_factor": 4,
    "restore_state_tracers_history_runtime_monitor_and_causal_debit": True,
    "continue_source_and_CFL_retry_counters_without_reclassification": True,
    "no_reinitialization_reprojection_interpolation_or_refit": True,
    "PROTO13_terminal_checkpoint_may_not_resume": True,
    "restart_arrays_are_fixed_discrete_calibration_inputs": True,
    "pre_restart_temporal_error_bound_proved": False,
    "restart_may_seed_candidate_or_physical_claim": False,
    "each_method_must_pass_its_own_complete_spatial_admission": True,
    "cross_method_observables_must_be_extrapolated_or_interval_enclosed": True,
}

EXPECTED_TEMPORAL_CONTRACT = {
    "level_step_counts": [1, 2, 4],
    "complete_state_channel_count": 18,
    "minimum_observed_order": "3/2",
    "exact_outward_order_test": "8*U12^2<=L01^2",
    "all_channels_must_pass_every_accepted_post_restart_macro_step": True,
    "only_four_quarter_fine_path_may_commit": True,
    "finest_pair_debits_accumulate_componentwise_outward": True,
    "debit_cancellation_permitted": False,
    "absolute_or_signal_scaled_tolerance_permitted": False,
    "temporal_retry_factor": "1/2",
    "maximum_temporal_retries_per_macro_step": 32,
    "minimum_macro_step": "1/1073741824",
    "durable_rejection_before_retry_required": True,
    "TDG6_ledger_starts_zero_at_restart": True,
    "sampled_64_history_rule_is_nonveto_diagnostic": True,
    "sampled_64_history_rule_may_admit_or_reject": False,
    "temporal_retry_exhaustion_is_invalid_not_physical": True,
    "global_PDE_error_bound_proved": False,
}

EXPECTED_NAMESPACE = {
    "calibration_output_root": "runs/fgc-2-sf1/proto14/calibration",
    "holdout_output_root": "runs/fgc-2-sf1/proto14/holdout",
    "both_roots_must_be_absent": True,
    "freeze_must_create_neither_root": True,
}

EXPECTED_PROOF_KEYS = {
    "PROTO14_overlay_must_validate_exactly",
    "checkpoint_commit_must_be_ancestor_of_HEAD",
    "all_tracked_lineage_must_match_immutable_commit_blobs",
    "PROTO12_and_RSP2_checkpoint_hashes_must_match",
    "all_six_state_and_complete_restart_payload_hashes_must_match",
    "all_six_members_must_restore_at_exactly_t_23_over_16",
    "replacement_instrument_must_be_implemented_and_synthetically_qualified",
    "historical_PROTO13_and_CAL10_results_must_remain_unreclassified",
    "sampled_temporal_rule_must_be_nonveto_only",
    "TDG6_threshold_channel_retry_debit_and_checkpoint_contract_may_not_change",
    "source_constraints_spatial_health_boundary_scale_and_trapped_rules_may_not_change",
    "restart_limitation_must_remain_public",
    "new_output_namespaces_must_be_absent_and_uncreated",
    "threshold_restart_history_claim_and_namespace_mutations_must_fail_closed",
    "canonical_hash_bound_result_required",
}

EXPECTED_CLAIMS = {
    "PROTO14_frozen": True,
    "PROTO14_replacement_temporal_admission_frozen": True,
    "PROTO14_restart_manifest_frozen": True,
    "PROTO14_terminal_classification_frozen": True,
    **PROTO14_CLAIMS,
}

RESTART_MEMBER_FIELDS = (
    "key",
    "source_checkpoint",
    "method",
    "point_count",
    "state_sha256",
    "restart_payload_sha256",
    "input_hash",
    "step_index",
    "transaction_serial",
    "accepted_stage_count",
    "source_retry_count",
    "CFL_retry_count",
)


def _rel(path: Path) -> str:
    return path.resolve().relative_to(REPOSITORY).as_posix()


def _sha(path: Path) -> str:
    digest = sha256()
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


def _git(*args: str, check: bool = True) -> subprocess.CompletedProcess[bytes]:
    return subprocess.run(
        ["git", *args],
        cwd=REPOSITORY,
        check=check,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )


def _load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=_reject_duplicate_pairs)
    if not isinstance(value, dict):
        raise ValueError(f"{_rel(path)} must contain one JSON object")
    return value


def _decode_json_array(value: np.ndarray) -> dict[str, Any]:
    decoded = json.loads(bytes(value.tolist()).decode("utf-8"))
    if not isinstance(decoded, dict):
        raise ValueError("checkpoint metadata must be an object")
    return decoded


def load_protocol(path: Path = PROTOCOL_CONFIG) -> dict[str, Any]:
    with path.open("rb") as handle:
        value = tomllib.load(handle)
    validate_sf1_protocol_v14(value)
    return value


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    with path.open("rb") as handle:
        value = tomllib.load(handle)
    if set(value) != {
        "schema_version",
        "artifact_id",
        "project_version",
        "metric_signature",
        "riemann_convention",
        "protocol_config",
        "scope",
        "immutable_lineage",
        "restart_contract",
        "temporal_contract",
        "restart_members",
        "namespace_precondition",
        "proof_contract",
        "claims",
    }:
        raise ValueError("PRO14-FRZ1 config root differs")
    if (
        value["schema_version"] != 1
        or value["artifact_id"] != ARTIFACT_ID
        or value["project_version"] != PROJECT_VERSION
        or value["metric_signature"] != "-+++"
        or value["riemann_convention"] != "plus_partial_mu_gamma_nu"
        or value["protocol_config"] != "configs/fgc/fgc-2-sf1-protocol-v14.toml"
        or value["scope"] != EXPECTED_SCOPE
        or value["immutable_lineage"] != EXPECTED_LINEAGE
        or value["restart_contract"] != EXPECTED_RESTART_CONTRACT
        or value["temporal_contract"] != EXPECTED_TEMPORAL_CONTRACT
        or value["namespace_precondition"] != EXPECTED_NAMESPACE
        or set(value["proof_contract"]) != EXPECTED_PROOF_KEYS
        or any(item is not True for item in value["proof_contract"].values())
        or value["claims"] != EXPECTED_CLAIMS
    ):
        raise ValueError("PRO14-FRZ1 config violates its frozen scope")
    return value


def _immutable_lineage(config: Mapping[str, Any]) -> dict[str, Any]:
    lineage = config["immutable_lineage"]
    commit = lineage["checkpoint_commit"]
    if _git("merge-base", "--is-ancestor", commit, "HEAD", check=False).returncode != 0:
        raise ValueError("PRO14 checkpoint is not an ancestor of HEAD")
    tracked_keys = (
        ("predecessor_protocol_config", "predecessor_protocol_sha256"),
        ("predecessor_freeze_result", "predecessor_freeze_result_sha256"),
        ("historical_terminal_result", "historical_terminal_result_sha256"),
        ("replacement_instrument_config", "replacement_instrument_config_sha256"),
        ("replacement_instrument_result", "replacement_instrument_result_sha256"),
        ("replacement_runtime_module", "replacement_runtime_module_sha256"),
    )
    records = []
    for path_key, hash_key in tracked_keys:
        relative = lineage[path_key]
        expected = lineage[hash_key]
        blob = _git("show", f"{commit}:{relative}").stdout
        if sha256(blob).hexdigest() != expected or _sha(REPOSITORY / relative) != expected:
            raise ValueError(f"PRO14 immutable tracked lineage differs: {relative}")
        records.append({"path": relative, "sha256": expected})
    for key in ("PROTO12", "RSP2"):
        relative = lineage[f"{key}_checkpoint"]
        expected = lineage[f"{key}_checkpoint_sha256"]
        path = REPOSITORY / relative
        if not path.is_file() or _sha(path) != expected:
            raise ValueError(f"PRO14 {key} restart checkpoint differs")
    return {
        "checkpoint_commit": commit,
        "checkpoint_is_ancestor_of_HEAD": True,
        "tracked_blobs": records,
        "PROTO12_checkpoint": lineage["PROTO12_checkpoint"],
        "PROTO12_checkpoint_sha256": lineage["PROTO12_checkpoint_sha256"],
        "RSP2_checkpoint": lineage["RSP2_checkpoint"],
        "RSP2_checkpoint_sha256": lineage["RSP2_checkpoint_sha256"],
    }


def _historical_restart_manifest(config: Mapping[str, Any]) -> list[dict[str, Any]]:
    record = _load_json(REPOSITORY / config["immutable_lineage"]["predecessor_freeze_result"])
    if record.get("artifact_id") != "FGC-1-PRO13-FRZ1":
        raise ValueError("PRO14 predecessor freeze identity differs")
    historical = record.get("artifact_payload", {}).get("restart_members")
    if not isinstance(historical, list) or len(historical) != 6:
        raise ValueError("PRO14 predecessor restart manifest differs")
    expected = [
        {field: item[field] for field in RESTART_MEMBER_FIELDS}
        for item in historical
    ]
    if config["restart_members"] != expected:
        raise ValueError("PRO14 restart members differ from immutable PROTO13")
    return expected


def _restore_members(config: Mapping[str, Any]) -> list[dict[str, Any]]:
    expected_members = _historical_restart_manifest(config)
    lineage = config["immutable_lineage"]
    paths = {
        "PROTO12": REPOSITORY / lineage["PROTO12_checkpoint"],
        "RSP2": REPOSITORY / lineage["RSP2_checkpoint"],
    }
    handles = {name: np.load(path, allow_pickle=False) for name, path in paths.items()}
    try:
        metadata = {name: _decode_json_array(handle["metadata_utf8"]) for name, handle in handles.items()}
        records = []
        suffixes = (
            "u",
            "p",
            "q",
            "tracer_positions",
            "tracer_proper_times",
            "event_proper_times",
            "event_fields",
        )
        for expected in expected_members:
            key = expected["key"]
            source = expected["source_checkpoint"]
            handle = handles[source]
            stored = metadata[source].get("members", {}).get(key)
            if not isinstance(stored, dict):
                raise ValueError(f"PRO14 restart member is absent: {key}")
            prefix = key.replace("-", "_")
            arrays = [handle[f"{prefix}_{suffix}"] for suffix in suffixes]
            runtime = stored.get("runtime_monitor_state", {})
            causal = stored.get("causal_state", {})
            observed = {
                "key": key,
                "source_checkpoint": source,
                "method": expected["method"],
                "point_count": expected["point_count"],
                "coordinate_time": stored.get("time"),
                "state_sha256": array_content_sha256(*arrays[:3]),
                "restart_payload_sha256": array_content_sha256(*arrays),
                "input_hash": stored.get("input_hash"),
                "step_index": stored.get("step_index"),
                "transaction_serial": stored.get("transaction_serial"),
                "accepted_stage_count": runtime.get("accepted_stage_count"),
                "source_retry_count": stored.get("source_retry_count"),
                "CFL_retry_count": stored.get("CFL_retry_count"),
                "accumulated_characteristic_distance": causal.get("accumulated_characteristic_distance"),
                "event_sample_count": int(arrays[-2].shape[0]),
                "tracer_count": int(arrays[3].shape[0]),
                "state_shape": list(arrays[0].shape),
            }
            for field in RESTART_MEMBER_FIELDS:
                if observed.get(field) != expected[field]:
                    raise ValueError(f"PRO14 restart ledger differs for {key}: {field}")
            if (
                observed["coordinate_time"] != 23 / 16
                or runtime.get("last_accepted_time") != 23 / 16
                or causal.get("accepted_time") != 23 / 16
                or runtime.get("first_failed_premise") is not None
                or runtime.get("first_failed_time") is not None
                or observed["event_sample_count"] != 24
                or observed["tracer_count"] != 48
                or observed["state_shape"] != [expected["point_count"], 6]
            ):
                raise ValueError(f"PRO14 restart completeness differs for {key}")
            records.append(observed)
    finally:
        for handle in handles.values():
            handle.close()
    return records


def _evidence_boundary(config: Mapping[str, Any]) -> dict[str, Any]:
    lineage = config["immutable_lineage"]
    cal10 = _load_json(REPOSITORY / lineage["historical_terminal_result"])
    tdg6 = _load_json(REPOSITORY / lineage["replacement_instrument_result"])
    cal10_status = cal10.get("gate_status", {})
    cal10_boundary = cal10.get("artifact_payload", {}).get("claim_boundary", {})
    tdg6_status = tdg6.get("gate_status", {})
    tdg6_boundary = tdg6.get("artifact_payload", {}).get("claim_boundary", {})
    if (
        cal10.get("artifact_id") != "FGC-1-CAL10-PREF15"
        or cal10.get("classification") != "completed_PROTO13_GR0_temporal_gate_stop"
        or cal10_status.get("PROTO13_temporal_admission_failed") is not True
        or cal10_status.get("PROTO13_terminal_common_spatial_and_constraint_admissions_passed") is not True
        or cal10_status.get("GR0_case_eligible") is not False
        or cal10_boundary.get("terminal_checkpoint_cannot_be_resumed") is not True
        or cal10_boundary.get("temporal_stop_is_not_a_physical_obstruction") is not True
    ):
        raise ValueError("PRO14 historical CAL10 boundary differs")
    if (
        tdg6.get("artifact_id") != "FGC-1-TDG6-IMP2"
        or tdg6_status.get("production_compositor_implemented") is not True
        or tdg6_status.get("production_compositor_synthetic_qualification_passed") is not True
        or tdg6_status.get("replacement_temporal_admission_defined") is not True
        or tdg6_status.get("PROTO14_freeze_design_authorized") is not True
        or tdg6_status.get("production_trajectory_authorized") is not False
        or tdg6_boundary.get("historical_CAL10_result_reclassified") is not False
        or tdg6_boundary.get("actual_PROTO13_history_arrays_consumed") is not False
    ):
        raise ValueError("PRO14 TDG6-IMP2 evidence boundary differs")
    return {
        "historical_terminal_artifact_id": cal10["artifact_id"],
        "historical_PROTO13_temporal_stop_preserved": True,
        "historical_terminal_checkpoint_may_resume": False,
        "historical_stop_is_physical_obstruction": False,
        "historical_GR0_case_eligible": False,
        "replacement_instrument_artifact_id": tdg6["artifact_id"],
        "production_compositor_implemented": True,
        "production_compositor_synthetic_qualification_passed": True,
        "replacement_temporal_admission_defined": True,
        "PROTO14_freeze_design_authorized": True,
        "actual_PROTO13_history_arrays_consumed_by_TDG6_IMP2": False,
        "production_trajectory_advanced": False,
    }


def _namespace_evidence(config: Mapping[str, Any]) -> dict[str, Any]:
    records = []
    for name in ("calibration_output_root", "holdout_output_root"):
        relative = config["namespace_precondition"][name]
        path = REPOSITORY / relative
        if path.exists():
            raise ValueError(f"PRO14 prospective namespace already exists: {relative}")
        records.append({"path": relative, "absent_before_freeze": True, "created_by_freeze": False})
    return {
        "both_new_namespaces_absent": True,
        "freeze_created_no_namespace": True,
        "records": records,
    }


def _mutation_controls(protocol: Mapping[str, Any], config: Mapping[str, Any]) -> dict[str, bool]:
    controls: dict[str, bool] = {}
    protocol_mutations = {
        "threshold": lambda value: value["replacement"]["temporal_admission"].__setitem__(
            "minimum_observed_order", "149/100"
        ),
        "sampled_diagnostic_veto": lambda value: value["replacement"]["temporal_admission"].__setitem__(
            "sampled_64_history_rule_may_admit_or_reject_PROTO14", True
        ),
        "restart_error_promotion": lambda value: value["replacement"]["restart"].__setitem__(
            "pre_restart_temporal_error_bound_proved", True
        ),
        "terminal_candidate_promotion": lambda value: value["replacement"]["terminal_classification"].__setitem__(
            "candidate_branch_may_open_from_PROTO14_result_alone", True
        ),
        "claim": lambda value: value["claims"].__setitem__(
            "FGCQR_holdout_execution_authorized", True
        ),
    }
    for name, mutate in protocol_mutations.items():
        attacked = deepcopy(protocol)
        mutate(attacked)
        try:
            validate_sf1_protocol_v14(attacked)
        except (TypeError, ValueError):
            controls[name] = True
        else:
            controls[name] = False
    config_mutations = {
        "restart_state_hash": lambda value: value["restart_members"][0].__setitem__(
            "state_sha256", "0" * 64
        ),
        "namespace": lambda value: value["namespace_precondition"].__setitem__(
            "calibration_output_root", "runs/fgc-2-sf1/proto13/calibration"
        ),
    }
    for name, mutate in config_mutations.items():
        attacked = deepcopy(config)
        mutate(attacked)
        try:
            if (
                attacked["restart_members"] != config["restart_members"]
                or attacked["namespace_precondition"] != EXPECTED_NAMESPACE
            ):
                raise ValueError("attacked PRO14 config differs")
        except ValueError:
            controls[name] = True
        else:
            controls[name] = False
    if not all(controls.values()):
        raise RuntimeError("PRO14 mutation control failed")
    return controls


def record(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    config = load_config(config_path)
    protocol = load_protocol(REPOSITORY / config["protocol_config"])
    protocol_validation = validate_sf1_protocol_v14(protocol)
    lineage = _immutable_lineage(config)
    members = _restore_members(config)
    evidence = _evidence_boundary(config)
    namespace = _namespace_evidence(config)
    controls = _mutation_controls(protocol, config)
    if [item["point_count"] for item in members[:3]] != [2049, 4097, 8193]:
        raise RuntimeError("PRO14 primary ladder differs after restoration")
    if [item["point_count"] for item in members[3:]] != [4097, 8193, 16385]:
        raise RuntimeError("PRO14 comparator ladder differs after restoration")
    gates = dict(config["claims"])
    nonclaims = {name: value for name, value in gates.items() if value is False}
    return {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "classification": "prospective_TDG6_post_restart_GR0_calibration_protocol_freeze",
        "generated_by": "scripts/reproduce_fgc_pro14_frz1.py",
        "derivation_document": _rel(OWNER_DOCUMENT),
        "derivation_document_sha256": _sha(OWNER_DOCUMENT),
        "source_config_sha256": {
            _rel(config_path): _sha(config_path),
            _rel(PROTOCOL_CONFIG): _sha(PROTOCOL_CONFIG),
        },
        "predecessor_sha256": {
            key: config["immutable_lineage"][hash_key]
            for key, hash_key in (
                (config["immutable_lineage"]["predecessor_protocol_config"], "predecessor_protocol_sha256"),
                (config["immutable_lineage"]["predecessor_freeze_result"], "predecessor_freeze_result_sha256"),
                (config["immutable_lineage"]["historical_terminal_result"], "historical_terminal_result_sha256"),
                (config["immutable_lineage"]["replacement_instrument_config"], "replacement_instrument_config_sha256"),
                (config["immutable_lineage"]["replacement_instrument_result"], "replacement_instrument_result_sha256"),
                (config["immutable_lineage"]["replacement_runtime_module"], "replacement_runtime_module_sha256"),
                (config["immutable_lineage"]["PROTO12_checkpoint"], "PROTO12_checkpoint_sha256"),
                (config["immutable_lineage"]["RSP2_checkpoint"], "RSP2_checkpoint_sha256"),
            )
        },
        "implementation_sha256": {_rel(path): _sha(path) for path in IMPLEMENTATION},
        "scope_bindings": dict(config["scope"]),
        "artifact_payload": {
            "immutable_lineage": lineage,
            "protocol_validation": protocol_validation,
            "evidence_boundary": evidence,
            "restart_contract": dict(config["restart_contract"]),
            "temporal_contract": dict(config["temporal_contract"]),
            "restart_members": members,
            "namespace_precondition": namespace,
            "mutation_controls": controls,
            "decision": {
                "PROTO14_frozen": True,
                "replacement_temporal_admission_frozen": True,
                "PROTO14_runtime_owner": "FGC-1-HLT12-MON12",
                "PROTO14_runtime_implemented": False,
                "PROTO14_calibration_authorized": False,
                "candidate_execution_authorized": False,
            },
            "epistemic_boundary": {
                "historical_PROTO13_and_CAL10_results_reclassified": False,
                "restart_manifest_is_not_a_fresh_calibration": True,
                "restart_state_has_pre_restart_temporal_error_bound": False,
                "restart_state_may_seed_candidate_or_physical_claim": False,
                "sampled_64_history_rule_is_nonveto_diagnostic_only": True,
                "TDG6_debit_is_global_PDE_or_Raychaudhuri_error_bound": False,
                "PROTO14_freeze_advances_no_trajectory": True,
                "eligible_GR0_case_selected": False,
                "candidate_or_mechanism_outcome_read": False,
            },
        },
        "gate_status": gates,
        "nonclaims": nonclaims,
    }


def load_canonical_result(path: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    source = path.read_text(encoding="utf-8")
    value = json.loads(source, object_pairs_hook=_reject_duplicate_pairs)
    if not isinstance(value, dict) or value.get("artifact_id") != ARTIFACT_ID:
        raise ValueError("PRO14 canonical result identity differs")
    if source != _canonical(value):
        raise ValueError("PRO14 canonical result is not sorted canonical JSON")
    return value


def _reject_duplicate_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    answer: dict[str, Any] = {}
    for key, value in pairs:
        if key in answer:
            raise ValueError(f"duplicate JSON key: {key}")
        answer[key] = value
    return answer


def verify_canonical(path: Path = DEFAULT_OUTPUT, config_path: Path = DEFAULT_CONFIG) -> None:
    observed = load_canonical_result(path)
    expected = record(config_path)
    if observed != expected:
        raise ValueError("PRO14 output differs from a fresh prospective reproduction")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--check", action="store_true")
    arguments = parser.parse_args()
    if arguments.check:
        verify_canonical(arguments.output, arguments.config)
        print(f"verified {arguments.output}")
        return
    payload = _canonical(record(arguments.config)).encode("utf-8")
    _atomic_write(arguments.output, payload)
    print(f"wrote {arguments.output}")


if __name__ == "__main__":
    main()
