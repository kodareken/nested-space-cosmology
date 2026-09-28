#!/usr/bin/env python3
"""Reproduce PROTO13's prospective method-owned restart freeze.

The certificate consumes the immutable PREF14 checkpoint result and the two
ignored restart checkpoints.  It restores six GR-0 members at ``t=23/16``,
binds their full state/tracer/history payloads, freezes the method-owned
resolution ladders and new namespaces, and advances no trajectory.
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
from recursive_horizons.fgc.evolution.protocol_v13 import (  # noqa: E402
    PROTO13_AMENDMENT,
    PROTO13_CLAIMS,
    PROTO13_METHOD_OWNED_LADDERS,
    PROTO13_RESTART,
    validate_sf1_protocol_v13,
)


ARTIFACT_ID = "FGC-1-PRO13-FRZ1"
PROJECT_VERSION = "0.11.0"
DEFAULT_CONFIG = REPOSITORY / "configs/fgc/fgc-1-pro13-frz1.toml"
DEFAULT_OUTPUT = REPOSITORY / "results/fgc-1-pro13-frz1.json"
OWNER_DOCUMENT = REPOSITORY / "docs/fgc-pro13-frz1.md"
PROTOCOL_CONFIG = REPOSITORY / "configs/fgc/fgc-2-sf1-protocol-v13.toml"
PROTOCOL_MODULE = REPOSITORY / "src/recursive_horizons/fgc/evolution/protocol_v13.py"
IMPLEMENTATION = (Path(__file__).resolve(), PROTOCOL_MODULE)

EXPECTED_SCOPE = {
    "target_protocol": "FGC-2-SF1-PROTO13",
    "predecessor_protocol": "FGC-2-SF1-PROTO12",
    "resolution_evidence_artifact": "FGC-1-RSP2-PREF14",
    "calibration_branch": "GR-0",
    "calibration_amplitude": "3",
    "holdout_branch": "FGC-QR",
    "freeze_role": "premise_only_method_owned_ladder_restart_protocol_certificate",
    "PROTO12_and_RSP2_GR0_trajectories_inspected": True,
    "PROTO13_trajectory_read": False,
    "SGBL_trajectory_read": False,
    "FGCQR_trajectory_read": False,
    "collapse_or_trapped_outcome_classified": False,
    "mechanism_question_answered": False,
}

EXPECTED_LINEAGE = {
    "checkpoint_commit": "ef5c78eb5c981883ea47dacb7a968bb660eff4dd",
    "predecessor_protocol_config": "configs/fgc/fgc-2-sf1-protocol-v12.toml",
    "predecessor_protocol_sha256": "2c33ef4a703e918683c886eb216bfc77af5a48d477c1e1556e0b5176d44fbe7f",
    "resolution_evidence_config": "configs/fgc/fgc-1-rsp2-pref14.toml",
    "resolution_evidence_config_sha256": "d9d9a6ebae797f5f59646b7ae7d998740516c775038e2e2ff356577b35373b42",
    "resolution_evidence_result": "results/fgc-1-rsp2-pref14.json",
    "resolution_evidence_result_sha256": "dd09cbfff2f3f22729adc6580fd577ab183a68960fedaefcc4e3a90d0bc1b0c1",
    "calibration_diagnosis_result": "results/fgc-1-cal9-pref13.json",
    "calibration_diagnosis_result_sha256": "69a2aeff732081d1d1be2e6fd49d020335e2fa13f737d5e21484d080fb1f5fd2",
    "PROTO12_checkpoint": "runs/fgc-2-sf1/proto12/calibration/latest-checkpoint.npz",
    "PROTO12_checkpoint_sha256": "c784d4706911029883d14e1763c2d36c5dbe0e619a4e50416cc6d666d5b7ea17",
    "RSP2_checkpoint": "runs/fgc-2-sf1/rsp2/amplitude3-ssprk3-momentum/latest-checkpoint.npz",
    "RSP2_checkpoint_sha256": "000546dcf3726c7882e80b7e70b64a5d2e060e2e0f7a567e174f714432ae7f34",
}

EXPECTED_RESTART = {
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
    "maximum_nested_tail_ratio": "1/4",
    "minimum_consecutive_qualified_trapped_common_events": 8,
    "trapped_sign_margin_over_combined_error_factor": 4,
    "restore_state_tracers_history_runtime_monitor_and_causal_debit": True,
    "continue_retry_counters_without_reclassification": True,
    "no_reinitialization_reprojection_interpolation_or_refit": True,
    "each_method_must_pass_its_own_complete_admission": True,
    "cross_method_observables_must_be_extrapolated_or_interval_enclosed": True,
}

EXPECTED_MEMBERS = [
    {
        "key": "RK4-2049",
        "source_checkpoint": "PROTO12",
        "method": "RK4",
        "point_count": 2049,
        "state_sha256": "d119a3b0ab450b0f0916b801254bbbfa964c8e5e95223cb8e51d9034f62c665e",
        "restart_payload_sha256": "4bb04ae97c09f66881084b5816d14ce300c54152ae11c182f53b64d5e887c072",
        "input_hash": "3bfc9a137f20353e9ead2f742f7780b049a11ab576f6e4de98fe4e5d0ad5c0d3",
        "step_index": 342,
        "transaction_serial": 1710,
        "accepted_stage_count": 1710,
        "source_retry_count": 0,
        "CFL_retry_count": 217,
    },
    {
        "key": "RK4-4097",
        "source_checkpoint": "PROTO12",
        "method": "RK4",
        "point_count": 4097,
        "state_sha256": "9d27ca221fbac2a5b3b86a8d9aef49a8b8ac451b255b80529563b8cb7f8b73a9",
        "restart_payload_sha256": "cc6d3b6fbb3a5b6abd93f69f7f59ba67c38e9b7d332fac991a88847970b57aac",
        "input_hash": "ffbde443e99535e74900b1e04c1e60cc27da105dd5d605286e77a6c77620c0bb",
        "step_index": 673,
        "transaction_serial": 3365,
        "accepted_stage_count": 3365,
        "source_retry_count": 0,
        "CFL_retry_count": 443,
    },
    {
        "key": "RK4-8193",
        "source_checkpoint": "PROTO12",
        "method": "RK4",
        "point_count": 8193,
        "state_sha256": "d5b4cfe25199862fff5be7d524cc1c3a6030cc634d268cedb626c3ced7220f4a",
        "restart_payload_sha256": "b95ea0be0779f712fb0c354f9c10f28b9ccc6c737fb403d5710dbd07e0e506f9",
        "input_hash": "076bebfc71482793017d2d306001b3416f5ad7a4088f35351d4d26322a09d12c",
        "step_index": 976,
        "transaction_serial": 4880,
        "accepted_stage_count": 4880,
        "source_retry_count": 0,
        "CFL_retry_count": 161,
    },
    {
        "key": "SSPRK3-4097",
        "source_checkpoint": "PROTO12",
        "method": "SSPRK3",
        "point_count": 4097,
        "state_sha256": "ab96602c0f05633375710efee2be01be84e24c32a31d07d6ef23390788a72162",
        "restart_payload_sha256": "dc056db806e5d983852135138dd48dbb9b995a8585f56900258bece7e16e7016",
        "input_hash": "3cc3dc3a262f83cc19f07c49d477f5e43309f0c9e0f16bdeaf6d0c7dc02bd48f",
        "step_index": 631,
        "transaction_serial": 2524,
        "accepted_stage_count": 2524,
        "source_retry_count": 0,
        "CFL_retry_count": 358,
    },
    {
        "key": "SSPRK3-8193",
        "source_checkpoint": "PROTO12",
        "method": "SSPRK3",
        "point_count": 8193,
        "state_sha256": "c46abc3f0f72fa2a13c1fee7bb329b94698cc2c26003406ae3eb1592f11b648c",
        "restart_payload_sha256": "b262bf4bd180717d42a7fa4d502f006bac0950669bc9e86c474ad6e4a21ce7ad",
        "input_hash": "d4bd39d36dfe749cc6ec04e5f112fc5544e3defdce67315ee871d2d8cfc75b5f",
        "step_index": 987,
        "transaction_serial": 3948,
        "accepted_stage_count": 3948,
        "source_retry_count": 0,
        "CFL_retry_count": 187,
    },
    {
        "key": "SSPRK3-16385",
        "source_checkpoint": "RSP2",
        "method": "SSPRK3",
        "point_count": 16385,
        "state_sha256": "91e64d0864390cb4ac63cfdbdbbcb8ca5ad6c2f47cdd4d6741e05ef560adad0d",
        "restart_payload_sha256": "62784c8cc63eedc210366f146f28d535343b966e1f988a4298894e78ec560884",
        "input_hash": "2946239e5d23b4e7a3cd74419a0b9c214eb89195f72bc84d4b89f39b6ef01963",
        "step_index": 1771,
        "transaction_serial": 7084,
        "accepted_stage_count": 7084,
        "source_retry_count": 0,
        "CFL_retry_count": 0,
    },
]

EXPECTED_NAMESPACE = {
    "calibration_output_root": "runs/fgc-2-sf1/proto13/calibration",
    "holdout_output_root": "runs/fgc-2-sf1/proto13/holdout",
    "both_roots_must_be_absent": True,
    "freeze_must_create_neither_root": True,
}

EXPECTED_PROOF_KEYS = {
    "PROTO13_overlay_must_validate_exactly",
    "checkpoint_commit_must_be_ancestor_of_HEAD",
    "all_tracked_lineage_must_match_immutable_commit_blobs",
    "PROTO12_and_RSP2_checkpoint_hashes_must_match",
    "all_six_state_and_complete_restart_payload_hashes_must_match",
    "all_six_members_must_restore_at_exactly_t_23_over_16",
    "stage_step_transaction_retry_and_input_ledgers_must_match",
    "method_owned_ladders_must_be_exact_and_disjointly_owned",
    "RSP2_old_miss_and_finer_pass_must_remain_public",
    "source_constraints_spatial_health_boundary_and_scale_rules_may_not_change",
    "cross_method_raw_finest_grid_collocation_is_forbidden",
    "new_output_namespaces_must_be_absent_and_uncreated",
    "threshold_ladder_restart_claim_and_namespace_mutations_must_fail_closed",
    "canonical_hash_bound_result_required",
}

EXPECTED_CLAIMS = {
    "PROTO13_frozen": True,
    "PROTO13_method_owned_ladders_frozen": True,
    "PROTO13_restart_manifest_frozen": True,
    **PROTO13_CLAIMS,
}


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


def _decode_json_array(value: np.ndarray) -> dict[str, Any]:
    decoded = json.loads(bytes(value.tolist()).decode("utf-8"))
    if not isinstance(decoded, dict):
        raise ValueError("checkpoint metadata must be an object")
    return decoded


def load_protocol(path: Path = PROTOCOL_CONFIG) -> dict[str, Any]:
    with path.open("rb") as handle:
        value = tomllib.load(handle)
    validate_sf1_protocol_v13(value)
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
        "restart_members",
        "namespace_precondition",
        "proof_contract",
        "claims",
    }:
        raise ValueError("PRO13-FRZ1 config root differs")
    if (
        value["schema_version"] != 1
        or value["artifact_id"] != ARTIFACT_ID
        or value["project_version"] != PROJECT_VERSION
        or value["metric_signature"] != "-+++"
        or value["riemann_convention"] != "plus_partial_mu_gamma_nu"
        or value["protocol_config"]
        != "configs/fgc/fgc-2-sf1-protocol-v13.toml"
        or value["scope"] != EXPECTED_SCOPE
        or value["immutable_lineage"] != EXPECTED_LINEAGE
        or value["restart_contract"] != EXPECTED_RESTART
        or value["restart_members"] != EXPECTED_MEMBERS
        or value["namespace_precondition"] != EXPECTED_NAMESPACE
        or set(value["proof_contract"]) != EXPECTED_PROOF_KEYS
        or any(item is not True for item in value["proof_contract"].values())
        or value["claims"] != EXPECTED_CLAIMS
    ):
        raise ValueError("PRO13-FRZ1 config violates its frozen scope")
    return value


def _immutable_lineage(config: Mapping[str, Any]) -> dict[str, Any]:
    lineage = config["immutable_lineage"]
    commit = lineage["checkpoint_commit"]
    ancestor = _git("merge-base", "--is-ancestor", commit, "HEAD", check=False)
    if ancestor.returncode != 0:
        raise ValueError("PRO13 checkpoint is not an ancestor of HEAD")
    tracked = (
        (
            lineage["predecessor_protocol_config"],
            lineage["predecessor_protocol_sha256"],
        ),
        (
            lineage["resolution_evidence_config"],
            lineage["resolution_evidence_config_sha256"],
        ),
        (
            lineage["resolution_evidence_result"],
            lineage["resolution_evidence_result_sha256"],
        ),
        (
            lineage["calibration_diagnosis_result"],
            lineage["calibration_diagnosis_result_sha256"],
        ),
    )
    records = []
    for relative, expected in tracked:
        blob = _git("show", f"{commit}:{relative}").stdout
        observed = sha256(blob).hexdigest()
        if observed != expected or _sha(REPOSITORY / relative) != expected:
            raise ValueError(f"PRO13 immutable tracked lineage differs: {relative}")
        records.append({"path": relative, "sha256": expected})
    for key in ("PROTO12", "RSP2"):
        relative = lineage[f"{key}_checkpoint"]
        expected = lineage[f"{key}_checkpoint_sha256"]
        path = REPOSITORY / relative
        if not path.is_file() or _sha(path) != expected:
            raise ValueError(f"PRO13 {key} restart checkpoint differs")
    return {
        "checkpoint_commit": commit,
        "checkpoint_is_ancestor_of_HEAD": True,
        "tracked_blobs": records,
        "PROTO12_checkpoint": lineage["PROTO12_checkpoint"],
        "PROTO12_checkpoint_sha256": lineage["PROTO12_checkpoint_sha256"],
        "RSP2_checkpoint": lineage["RSP2_checkpoint"],
        "RSP2_checkpoint_sha256": lineage["RSP2_checkpoint_sha256"],
    }


def _restore_members(config: Mapping[str, Any]) -> list[dict[str, Any]]:
    lineage = config["immutable_lineage"]
    paths = {
        "PROTO12": REPOSITORY / lineage["PROTO12_checkpoint"],
        "RSP2": REPOSITORY / lineage["RSP2_checkpoint"],
    }
    handles = {
        name: np.load(path, allow_pickle=False) for name, path in paths.items()
    }
    try:
        metadata = {
            name: _decode_json_array(handle["metadata_utf8"])
            for name, handle in handles.items()
        }
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
        for expected in config["restart_members"]:
            key = expected["key"]
            source = expected["source_checkpoint"]
            handle = handles[source]
            member = metadata[source].get("members", {}).get(key)
            if not isinstance(member, dict):
                raise ValueError(f"PRO13 restart member is absent: {key}")
            prefix = key.replace("-", "_")
            arrays = [handle[f"{prefix}_{suffix}"] for suffix in suffixes]
            state_hash = array_content_sha256(*arrays[:3])
            restart_hash = array_content_sha256(*arrays)
            runtime = member.get("runtime_monitor_state", {})
            causal = member.get("causal_state", {})
            observed = {
                "key": key,
                "source_checkpoint": source,
                "method": expected["method"],
                "point_count": expected["point_count"],
                "coordinate_time": member.get("time"),
                "state_sha256": state_hash,
                "restart_payload_sha256": restart_hash,
                "input_hash": member.get("input_hash"),
                "step_index": member.get("step_index"),
                "transaction_serial": member.get("transaction_serial"),
                "accepted_stage_count": runtime.get("accepted_stage_count"),
                "source_retry_count": member.get("source_retry_count"),
                "CFL_retry_count": member.get("CFL_retry_count"),
                "accumulated_characteristic_distance": causal.get(
                    "accumulated_characteristic_distance"
                ),
                "previous_speed_upper": causal.get("previous_speed_upper"),
                "event_sample_count": int(arrays[-2].shape[0]),
                "tracer_count": int(arrays[3].shape[0]),
                "state_shape": list(arrays[0].shape),
            }
            for name, value in expected.items():
                if name in {"method", "point_count", "source_checkpoint", "key"}:
                    continue
                if observed.get(name) != value:
                    raise ValueError(f"PRO13 restart ledger differs for {key}: {name}")
            if (
                observed["coordinate_time"] != 1.4375
                or runtime.get("last_accepted_time") != 1.4375
                or causal.get("accepted_time") != 1.4375
                or runtime.get("first_failed_premise") is not None
                or runtime.get("first_failed_time") is not None
                or observed["event_sample_count"] != 24
                or observed["tracer_count"] != 48
                or observed["state_shape"] != [expected["point_count"], 6]
            ):
                raise ValueError(f"PRO13 restart completeness differs for {key}")
            records.append(observed)
    finally:
        for handle in handles.values():
            handle.close()
    if [item["key"] for item in records] != [item["key"] for item in EXPECTED_MEMBERS]:
        raise ValueError("PRO13 restart member order differs")
    return records


def _resolution_evidence(config: Mapping[str, Any]) -> dict[str, Any]:
    path = REPOSITORY / config["immutable_lineage"]["resolution_evidence_result"]
    record = json.loads(path.read_text(encoding="utf-8"))
    status = record.get("gate_status", {})
    reduction = record.get("artifact_payload", {}).get(
        "independent_constraint_reduction", {}
    )
    if (
        record.get("artifact_id") != "FGC-1-RSP2-PREF14"
        or status.get("RSP2_target_order_cleared") is not True
        or status.get("RSP2_complete_constraint_admission_passed") is not True
        or status.get("PROTO13_design_may_begin") is not True
        or status.get("PROTO13_frozen") is not False
        or reduction.get("point_counts") != [4097, 8193, 16385]
        or reduction.get("target_adjacent_pair_orders")
        != [1.4990947384363291, 1.9817436463819333]
        or reduction.get("minimum_finite_finest_pair_order")
        != 1.8428133958883794
    ):
        raise ValueError("PRO13 RSP2 evidence differs")
    return {
        "artifact_id": record["artifact_id"],
        "preserved_coarse_pair_order": 1.4990947384363291,
        "finer_pair_order": 1.9817436463819333,
        "minimum_finite_complete_constraint_order": 1.8428133958883794,
        "target_and_complete_admission_passed": True,
        "PROTO13_design_may_begin": True,
        "is_fresh_PROTO13_calibration": False,
        "is_candidate_or_physical_evidence": False,
    }


def _namespace_evidence(config: Mapping[str, Any]) -> dict[str, Any]:
    records = []
    for name in ("calibration_output_root", "holdout_output_root"):
        relative = config["namespace_precondition"][name]
        path = REPOSITORY / relative
        if path.exists():
            raise ValueError(f"PRO13 prospective namespace already exists: {relative}")
        records.append(
            {
                "path": relative,
                "absent_before_freeze": True,
                "created_by_freeze": False,
            }
        )
    return {
        "both_new_namespaces_absent": True,
        "freeze_created_no_namespace": True,
        "records": records,
    }


def _mutation_controls(protocol: Mapping[str, Any], config: Mapping[str, Any]) -> dict[str, bool]:
    controls: dict[str, bool] = {}
    for name, mutate in {
        "primary_ladder": lambda value: value["replacement"]["method_owned_ladders"].__setitem__(
            "primary_point_counts", [4097, 8193, 16385]
        ),
        "comparator_ladder": lambda value: value["replacement"]["method_owned_ladders"].__setitem__(
            "comparator_point_counts", [2049, 4097, 8193]
        ),
        "threshold": lambda value: value["replacement"]["method_owned_ladders"].__setitem__(
            "minimum_constraint_finest_pair_order", "149/100"
        ),
        "restart_time": lambda value: value["replacement"]["restart"].__setitem__(
            "restart_coordinate_time", "3/2"
        ),
        "claim": lambda value: value["claims"].__setitem__(
            "FGCQR_holdout_execution_authorized", True
        ),
    }.items():
        attacked = deepcopy(protocol)
        mutate(attacked)
        try:
            validate_sf1_protocol_v13(attacked)
        except (TypeError, ValueError):
            controls[name] = True
        else:
            controls[name] = False
    for name, mutate in {
        "restart_state_hash": lambda value: value["restart_members"][0].__setitem__(
            "state_sha256", "0" * 64
        ),
        "namespace": lambda value: value["namespace_precondition"].__setitem__(
            "calibration_output_root", "runs/fgc-2-sf1/proto12/calibration"
        ),
    }.items():
        attacked = deepcopy(config)
        mutate(attacked)
        try:
            # Validate the pure serialized contract without touching namespaces.
            if (
                attacked["restart_members"] != EXPECTED_MEMBERS
                or attacked["namespace_precondition"] != EXPECTED_NAMESPACE
            ):
                raise ValueError("attacked PRO13 config differs")
        except ValueError:
            controls[name] = True
        else:
            controls[name] = False
    if not all(controls.values()):
        raise RuntimeError("PRO13 mutation control failed")
    return controls


def record(config_path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    config = load_config(config_path)
    protocol = load_protocol(REPOSITORY / config["protocol_config"])
    protocol_validation = validate_sf1_protocol_v13(protocol)
    lineage = _immutable_lineage(config)
    members = _restore_members(config)
    resolution = _resolution_evidence(config)
    namespace = _namespace_evidence(config)
    controls = _mutation_controls(protocol, config)
    if [item["point_count"] for item in members[:3]] != [2049, 4097, 8193]:
        raise RuntimeError("PRO13 primary ladder differs after restoration")
    if [item["point_count"] for item in members[3:]] != [4097, 8193, 16385]:
        raise RuntimeError("PROTO13 comparator ladder differs after restoration")
    gates = dict(config["claims"])
    nonclaims = {name: value for name, value in gates.items() if value is False}
    return {
        "schema_version": 1,
        "artifact_id": ARTIFACT_ID,
        "project_version": PROJECT_VERSION,
        "classification": "prospective_method_owned_ladder_restart_protocol_freeze",
        "generated_by": "scripts/reproduce_fgc_pro13_frz1.py",
        "derivation_document": _rel(OWNER_DOCUMENT),
        "derivation_document_sha256": _sha(OWNER_DOCUMENT),
        "source_config_sha256": {
            _rel(config_path): _sha(config_path),
            _rel(PROTOCOL_CONFIG): _sha(PROTOCOL_CONFIG),
        },
        "predecessor_sha256": {
            config["immutable_lineage"]["resolution_evidence_result"]: config[
                "immutable_lineage"
            ]["resolution_evidence_result_sha256"],
            config["immutable_lineage"]["calibration_diagnosis_result"]: config[
                "immutable_lineage"
            ]["calibration_diagnosis_result_sha256"],
            config["immutable_lineage"]["PROTO12_checkpoint"]: config[
                "immutable_lineage"
            ]["PROTO12_checkpoint_sha256"],
            config["immutable_lineage"]["RSP2_checkpoint"]: config[
                "immutable_lineage"
            ]["RSP2_checkpoint_sha256"],
        },
        "implementation_sha256": {_rel(path): _sha(path) for path in IMPLEMENTATION},
        "scope_bindings": dict(config["scope"]),
        "artifact_payload": {
            "immutable_lineage": lineage,
            "protocol_validation": protocol_validation,
            "resolution_evidence": resolution,
            "restart_contract": dict(config["restart_contract"]),
            "restart_members": members,
            "namespace_precondition": namespace,
            "mutation_controls": controls,
            "decision": {
                "PROTO13_frozen": True,
                "PROTO13_runtime_owner": "FGC-1-HLT11-MON11",
                "PROTO13_runtime_implemented": False,
                "PROTO13_calibration_authorized": False,
                "candidate_execution_authorized": False,
            },
            "epistemic_boundary": {
                "RSP2_pass_is_numerical_premise_evidence_only": True,
                "restart_manifest_is_not_a_fresh_calibration": True,
                "PROTO13_freeze_advances_no_trajectory": True,
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
        raise ValueError("PRO13 canonical result identity differs")
    if source != _canonical(value):
        raise ValueError("PRO13 canonical result is not sorted canonical JSON")
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
        raise ValueError("PRO13 output differs from a fresh prospective reproduction")


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
