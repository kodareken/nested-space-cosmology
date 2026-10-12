"""AUTH1's compact-to-complete production GenesisSpec adapter.

This module deliberately separates two operations.  ``build_calibration_genesis``
consumes only the compact PREF26 evidence already committed into an AUTH1
record; it never opens a raw checkpoint.  ``reimport_selected_sources`` is a
separate, explicitly requested same-invocation check which may reopen the two
declared predecessor containers before a launch authority is emitted.
"""
from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

from .proto17_pure_construction import Proto17ConstructionError, build_genesis, digest
from .protocol_v17 import MEMBER_KEYS
from .proto18_production_inputs import (
    Proto18ProductionInputError, VerifiedProductionCorpus, verify_production_corpus,
)


class Proto18Auth1InputError(ValueError):
    """Compact PREF26 evidence cannot lawfully seed a calibration GenesisSpec."""


_REQUIRED_EVIDENCE = frozenset({
    "target_protocol", "authority_identity", "restart_coordinate_time",
    "active_event_target_time", "common_event_index", "member_keys_in_canonical_order",
    "member_descriptors", "member_descriptor_set_sha256", "physical_input_manifest",
    "physical_input_manifest_sha256", "source_archive_manifest_sha256",
    "historical_replay_sha256", "final_PROTO17_GenesisSpec_construction_deferred_to_AUTH1",
    "future_runtime_adapter_and_runner_source_hashes_present", "AUTH1_input_evidence_sha256",
})
_SOURCE_HASH_FIELDS = (
    "protocol_config_sha256", "protocol_freeze_result_sha256", "hlt13_result_sha256",
    "run_plan_sha256", "runtime_module_sha256", "adapter_module_sha256", "runner_sha256",
)


@dataclass(frozen=True, slots=True)
class Auth1GenesisInputs:
    genesis_spec: Mapping[str, Any]
    genesis_checkpoint: Mapping[str, Any]
    auth1_input_evidence_sha256: str


def _sha(value: object, label: str) -> str:
    if not isinstance(value, str) or len(value) != 64:
        raise Proto18Auth1InputError(f"{label} must be a SHA-256 digest")
    try:
        int(value, 16)
    except ValueError as exc:
        raise Proto18Auth1InputError(f"{label} must be hexadecimal") from exc
    if value.lower() != value:
        raise Proto18Auth1InputError(f"{label} must be lowercase")
    return value


def validate_pref26_auth1_evidence(value: Mapping[str, Any]) -> dict[str, Any]:
    """Validate compact PREF26 evidence without reading any raw source path."""
    if not isinstance(value, Mapping) or set(value) != _REQUIRED_EVIDENCE:
        raise Proto18Auth1InputError("PREF26 AUTH1 evidence fields differ")
    evidence = deepcopy(dict(value))
    observed = evidence.pop("AUTH1_input_evidence_sha256")
    if _sha(observed, "AUTH1 input evidence") != digest(evidence):
        raise Proto18Auth1InputError("PREF26 AUTH1 evidence digest differs")
    evidence["AUTH1_input_evidence_sha256"] = observed
    if (evidence["target_protocol"] != "FGC-2-SF1-PROTO17"
            or evidence["authority_identity"] != "FGC-1-PRO18-AUTH1"
            or evidence["common_event_index"] != 23
            or evidence["member_keys_in_canonical_order"] != list(MEMBER_KEYS)
            or evidence["final_PROTO17_GenesisSpec_construction_deferred_to_AUTH1"] is not True
            or evidence["future_runtime_adapter_and_runner_source_hashes_present"] is not False):
        raise Proto18Auth1InputError("PREF26 AUTH1 evidence identity differs")
    for name in (
        "member_descriptor_set_sha256", "physical_input_manifest_sha256",
        "source_archive_manifest_sha256", "historical_replay_sha256",
    ):
        _sha(evidence[name], name)
    if digest(evidence["member_descriptors"]) != evidence["member_descriptor_set_sha256"]:
        raise Proto18Auth1InputError("PREF26 member descriptor digest differs")
    if digest(evidence["physical_input_manifest"]) != evidence["physical_input_manifest_sha256"]:
        raise Proto18Auth1InputError("PREF26 physical input manifest digest differs")
    return evidence


def source_identities_from_pref26_result(
    pref26_result: Mapping[str, Any],
    auth1_input_evidence: Mapping[str, Any],
) -> dict[str, dict[str, Any]]:
    """Extract the two immutable raw-source identities bound by PREF26.

    HLT15 calls this only after the complete PREF26 result blob has been
    authenticated through AUTH1's external Git authority.  The returned
    records are the exact pre-decompression identities accepted by
    :func:`verify_production_corpus`; no raw source path is opened here.
    """
    evidence = validate_pref26_auth1_evidence(auth1_input_evidence)
    if (not isinstance(pref26_result, Mapping)
            or pref26_result.get("artifact_id") != "FGC-1-PRO18-PREF26"):
        raise Proto18Auth1InputError("PREF26 result identity differs")
    payload = pref26_result.get("artifact_payload")
    manifest = payload.get("source_archive_manifest") if isinstance(payload, Mapping) else None
    if not isinstance(manifest, Mapping) or set(manifest) != {
        "sources", "resource_totals", "source_archive_manifest_sha256",
    }:
        raise Proto18Auth1InputError("PREF26 source archive manifest differs")
    bare = dict(manifest)
    observed = bare.pop("source_archive_manifest_sha256")
    if (_sha(observed, "PREF26 source archive manifest") != digest(bare)
            or observed != evidence["source_archive_manifest_sha256"]):
        raise Proto18Auth1InputError("PREF26 source archive manifest digest differs")
    sources = manifest["sources"]
    if (not isinstance(sources, list) or len(sources) != 2
            or [item.get("source_id") if isinstance(item, Mapping) else None for item in sources]
            != ["PROTO12", "RSP2"]):
        raise Proto18Auth1InputError("PREF26 source identity order differs")
    required = {
        "source_id", "checkpoint", "checkpoint_sha256", "checkpoint_byte_count",
        "event_log", "event_log_sha256", "event_log_byte_count",
    }
    identities: dict[str, dict[str, Any]] = {}
    for item in sources:
        if (not isinstance(item, Mapping) or not required <= set(item)
                or item.get("embedded_event_log_bytes_equal_external") is not True
                or item.get("embedded_event_log_sha256") != item.get("event_log_sha256")):
            raise Proto18Auth1InputError("PREF26 source identity record differs")
        identity = {name: item[name] for name in required}
        source = identity["source_id"]
        if (not isinstance(source, str) or source in identities
                or any(not isinstance(identity[name], str) or not identity[name]
                       for name in ("checkpoint", "event_log"))
                or _sha(identity["checkpoint_sha256"], f"{source} checkpoint")
                   != identity["checkpoint_sha256"]
                or _sha(identity["event_log_sha256"], f"{source} event log")
                   != identity["event_log_sha256"]
                or any(isinstance(identity[name], bool) or not isinstance(identity[name], int)
                       or identity[name] <= 0
                       for name in ("checkpoint_byte_count", "event_log_byte_count"))):
            raise Proto18Auth1InputError("PREF26 source byte identity differs")
        identities[source] = identity
    return identities


def build_calibration_genesis(
    auth1_input_evidence: Mapping[str, Any], *, campaign_id: str, namespace: str,
    source_hashes: Mapping[str, str], numerical_environment: Mapping[str, Any],
) -> Auth1GenesisInputs:
    """Create a complete calibration-only GenesisSpec from compact evidence.

    No path argument intentionally exists: AUTH1 construction cannot accidentally
    reopen raw containers.  The caller supplies prospective authority hashes for
    the protocol/runtime/adapter/runner that AUTH1 will itself pin.
    """
    evidence = validate_pref26_auth1_evidence(auth1_input_evidence)
    if not isinstance(campaign_id, str) or not campaign_id or ":" in campaign_id:
        raise Proto18Auth1InputError("calibration campaign identity differs")
    if not isinstance(namespace, str) or namespace != "runs/fgc-2-sf1/proto17/calibration":
        raise Proto18Auth1InputError("AUTH1 permits only the calibration namespace")
    if not isinstance(source_hashes, Mapping) or set(source_hashes) != set(_SOURCE_HASH_FIELDS):
        raise Proto18Auth1InputError("AUTH1 source-hash tuple differs")
    hashes = {name: _sha(source_hashes[name], name) for name in _SOURCE_HASH_FIELDS}
    if not isinstance(numerical_environment, Mapping) or not numerical_environment:
        raise Proto18Auth1InputError("numerical environment is absent")
    spec: dict[str, Any] = {
        "genesis_algorithm_id": "FGC-2-SF1-PROTO17-pure-construction-v1",
        "canonical_json_encoding": "utf-8;json-sort-keys;separators=(',',':');ensure-ascii=true;allow-nan=false",
        "state_object_schema_id": "FGC-2-SF1-PROTO17-state-object-v1",
        "tdg6_ledger_schema_id": "FGC-2-SF1-PROTO17-tdg6-ledger-v1",
        "cursor_schema_id": "FGC-2-SF1-PROTO17-cursor-v1",
        "checkpoint_schema_id": "FGC-2-SF1-PROTO17-checkpoint-v1",
        "protocol_artifact_id": evidence["target_protocol"],
        **hashes,
        "numerical_environment": deepcopy(dict(numerical_environment)),
        "campaign_id": campaign_id,
        "namespace": namespace,
        "journal_root_parent_sha256": "0" * 64,
        "common_event_index": evidence["common_event_index"],
        "restart_coordinate_time": deepcopy(evidence["restart_coordinate_time"]),
        "active_event_target_time": deepcopy(evidence["active_event_target_time"]),
        "member_keys_in_canonical_order": deepcopy(evidence["member_keys_in_canonical_order"]),
        "physical_input_manifest": deepcopy(evidence["physical_input_manifest"]),
        "physical_input_manifest_sha256": evidence["physical_input_manifest_sha256"],
        "member_descriptors": deepcopy(evidence["member_descriptors"]),
        # Derived fields are supplied after an initial pure derivation.
        "member_descriptor_set_sha256": evidence["member_descriptor_set_sha256"],
        "genesis_checkpoint_sha256": "0" * 64,
        "cursor_set_sha256": "0" * 64,
        "ledger_set_sha256": "0" * 64,
        "accepted_state_set_sha256": "0" * 64,
    }
    try:
        provisional = build_genesis(spec, _validate_derived=False)
        derived = provisional["checkpoint"]
        spec.update({
            "genesis_checkpoint_sha256": derived["checkpoint_sha256"],
            "cursor_set_sha256": derived["cursor_set_sha256"],
            "ledger_set_sha256": derived["ledger_set_sha256"],
            "accepted_state_set_sha256": derived["accepted_state_set_sha256"],
        })
        complete = build_genesis(spec)
    except (Proto17ConstructionError, KeyError, TypeError, ValueError) as exc:
        raise Proto18Auth1InputError("compact evidence cannot derive a legal GenesisSpec") from exc
    return Auth1GenesisInputs(
        genesis_spec=complete["genesis_spec"], genesis_checkpoint=complete["checkpoint"],
        auth1_input_evidence_sha256=evidence["AUTH1_input_evidence_sha256"],
    )


def reimport_selected_sources(
    repository_root: Path, pro18_config: Mapping[str, Any], *,
    expected_source_identities: Mapping[str, Any] | None = None,
) -> VerifiedProductionCorpus:
    """Explicit raw reimport used only in AUTH1's same-invocation gate."""
    try:
        return verify_production_corpus(
            repository_root, pro18_config,
            source_identities=expected_source_identities,
        )
    except Proto18ProductionInputError as exc:
        raise Proto18Auth1InputError("AUTH1 raw reimport failed") from exc


def reimport_matches_pref26_evidence(
    repository_root: Path, pro18_config: Mapping[str, Any],
    auth1_input_evidence: Mapping[str, Any], *,
    expected_source_identities: Mapping[str, Any] | None = None,
) -> VerifiedProductionCorpus:
    """HLT15 raw reimport and crosscheck against AUTH1's compact evidence.

    This is intentionally distinct from Genesis construction.  It proves that
    the same-invocation raw inputs still have the selected state/array identity
    already bound by PREF26 before a future root is observed.
    """
    evidence = validate_pref26_auth1_evidence(auth1_input_evidence)
    corpus = reimport_selected_sources(
        repository_root, pro18_config, expected_source_identities=expected_source_identities,
    )
    descriptors = evidence["member_descriptors"]
    if not isinstance(descriptors, list) or len(descriptors) != len(MEMBER_KEYS):
        raise Proto18Auth1InputError("compact descriptors are absent during reimport")
    for descriptor in descriptors:
        key = descriptor.get("member_key") if isinstance(descriptor, Mapping) else None
        member = corpus.members.get(key)
        state = descriptor.get("state_object") if isinstance(descriptor, Mapping) else None
        if (member is None or not isinstance(state, Mapping)
                or descriptor.get("method") != member.selector.method
                or descriptor.get("point_count") != member.selector.point_count
                or descriptor.get("step_index") != member.selector.step_index
                or descriptor.get("transaction_serial") != member.selector.transaction_serial
                or state.get("input_hash") != member.selector.input_sha256
                or state.get("source_retry_count") != member.selector.source_retry_count
                or state.get("CFL_retry_count") != member.selector.cfl_retry_count
                or member.physical_state_sha256 != state.get("physical_state_sha256")
                or member.restart_payload_sha256 != state.get("restart_payload_sha256")):
            raise Proto18Auth1InputError("HLT15 raw reimport differs from compact state identity")
        arrays = state.get("physical_arrays")
        if not isinstance(arrays, list) or {
            item.get("logical_name"): item.get("little_endian_c_bytes_sha256")
            for item in arrays if isinstance(item, Mapping)
        } != dict(member.array_sha256):
            raise Proto18Auth1InputError("HLT15 raw reimport differs from compact array identity")
    return corpus


__all__ = [
    "Auth1GenesisInputs", "Proto18Auth1InputError", "build_calibration_genesis",
    "reimport_selected_sources", "reimport_matches_pref26_evidence", "validate_pref26_auth1_evidence",
    "source_identities_from_pref26_result",
]
