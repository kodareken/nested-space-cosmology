"""Compose the read-only PREF26 source and legacy-history evidence.

This module owns no output namespace and performs no evolution.  It binds the
exact source identities frozen by PREF26, delegates byte-safe archive parsing
to :mod:`proto18_production_inputs`, delegates legacy grammar checks to
:mod:`proto18_historical_replay`, and returns only compact AUTH1 input
evidence.  The final PROTO17 GenesisSpec remains an AUTH1 responsibility
because its production runtime, adapter, and runner hashes do not exist here.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

from .proto17_pure_construction import digest, time_identity, validate_ledger
from .proto18_historical_replay import (
    Proto18HistoricalReplayError,
    replay_proto12_history,
    replay_rsp2_history,
)
from .proto18_production_inputs import (
    ARRAY_SUFFIXES,
    Proto18ProductionArchiveError,
    Proto18ProductionInputError,
    Proto18ProductionPathError,
    Proto18ProductionResourceError,
    Proto18ProductionSelectorError,
    Proto18ProductionSourceIdentityError,
    VerifiedProductionCorpus,
    verify_production_corpus,
)


class Proto18Pref26BinderError(ValueError):
    """A frozen PREF26 premise failed before any authority or namespace step."""

    def __init__(self, stop_id: str, detail: str) -> None:
        super().__init__(f"{stop_id}: {detail}")
        self.stop_id = stop_id
        self.detail = detail


@dataclass(frozen=True, slots=True)
class Pref26Evidence:
    source_archive_manifest: Mapping[str, Any]
    historical_replay: Mapping[str, Any]
    auth1_input_evidence: Mapping[str, Any]
    duty_status: Mapping[str, Any]


_EXPECTED_SOURCES = {
    "PROTO12": {
        "checkpoint": "runs/fgc-2-sf1/proto12/calibration/latest-checkpoint.npz",
        "checkpoint_sha256": "c784d4706911029883d14e1763c2d36c5dbe0e619a4e50416cc6d666d5b7ea17",
        "checkpoint_byte_count": 4_481_907,
        "checkpoint_npz_member_count": 44,
        "checkpoint_total_uncompressed_byte_count": 25_630_600,
        "checkpoint_largest_member_uncompressed_byte_count": 21_096_578,
        "event_log": "runs/fgc-2-sf1/proto12/calibration/events.jsonl",
        "event_log_sha256": "c2d35d07db2f6c983a449cdcf2db469360565cf392d78348dcafbd1d81234ba7",
        "event_log_byte_count": 21_096_450,
        "event_log_line_count": 49,
        "legacy_terminal_classification": "calibration_failed_no_eligible_GR0_case",
        "selector_keys": ["RK4-2049", "RK4-4097", "RK4-8193", "SSPRK3-4097", "SSPRK3-8193"],
    },
    "RSP2": {
        "checkpoint": "runs/fgc-2-sf1/rsp2/amplitude3-ssprk3-momentum/latest-checkpoint.npz",
        "checkpoint_sha256": "000546dcf3726c7882e80b7e70b64a5d2e060e2e0f7a567e174f714432ae7f34",
        "checkpoint_byte_count": 1_460_326,
        "checkpoint_npz_member_count": 9,
        "checkpoint_total_uncompressed_byte_count": 2_452_515,
        "checkpoint_largest_member_uncompressed_byte_count": 786_608,
        "event_log": "runs/fgc-2-sf1/rsp2/amplitude3-ssprk3-momentum/events.jsonl",
        "event_log_sha256": "ca9168c36ae0592fc69d20bfb32f1b501100d162b50d500a2f6b9529cdac2937",
        "event_log_byte_count": 15_060,
        "event_log_line_count": 25,
        "legacy_terminal_classification": "completed_target_and_complete_constraint_pass",
        "selector_keys": ["SSPRK3-16385"],
    },
}
_EXPECTED_RESOURCE_BOUNDS = {
    "maximum_total_source_container_bytes": 5_942_233,
    "maximum_total_legacy_event_log_bytes": 21_111_510,
    "maximum_total_source_bytes": 27_053_743,
    "maximum_total_NPZ_members": 53,
    "maximum_total_uncompressed_NPZ_bytes": 28_083_115,
    "maximum_single_uncompressed_NPZ_member_bytes": 21_096_578,
    "maximum_legacy_event_log_lines": 74,
    "network_access_permitted": False,
    "write_access_to_source_containers_or_logs": False,
    "write_access_to_future_output_roots": False,
}
_EXPECTED_STOPS = {
    "predecessor_lineage_drift": "PREF26_PREDECESSOR_LINEAGE_DRIFT",
    "source_path_missing_or_unsafe": "PREF26_DECLARED_SOURCE_PATH_MISSING_OR_UNSAFE",
    "source_byte_or_hash_drift": "PREF26_SOURCE_BYTE_OR_HASH_DRIFT",
    "resource_bound_exceeded": "PREF26_SOURCE_RESOURCE_BOUND_EXCEEDED",
    "NPZ_schema_or_selector_drift": "PREF26_NPZ_SCHEMA_OR_SELECTOR_DRIFT",
    "legacy_history_malformed_or_semantically_inconsistent": "PREF26_LEGACY_HISTORY_SEMANTIC_REPLAY_FAILED",
    "legacy_history_write_attempt": "PREF26_LEGACY_HISTORY_WRITE_FORBIDDEN",
    "future_output_root_access_attempt": "PREF26_FUTURE_OUTPUT_ROOT_ACCESS_FORBIDDEN",
    "authority_or_launch_promotion": "PREF26_AUTHORITY_OR_LAUNCH_PROMOTION_FORBIDDEN",
    "invalid_or_nonconverged_binder": "PREF26_INVALID_OR_NONCONVERGED_BINDER",
}


def _fail(config: Mapping[str, Any], owner: str, detail: str) -> None:
    stops = config.get("typed_stops")
    stop = stops.get(owner) if isinstance(stops, Mapping) else None
    if not isinstance(stop, str) or not stop:
        stop = "PREF26_INVALID_OR_NONCONVERGED_BINDER"
    raise Proto18Pref26BinderError(stop, detail)


def _source_contract(config: Mapping[str, Any]) -> Mapping[str, Mapping[str, Any]]:
    contract = config.get("source_contract")
    if not isinstance(contract, Mapping):
        _fail(config, "invalid_or_nonconverged_binder", "source contract is absent")
    sources = contract.get("source")
    expected_header = {
        "container_format": "NPZ", "source_container_count": 2,
        "legacy_event_log_count": 2, "logical_selector_count": 6,
        "selector_authority": "FGC-1-PRO18-FRZ1 selectors.member",
        "source_container_grammar": (
            "exactly_two_shared_NPZ_containers_plus_matching_legacy_events_jsonl; "
            "six selectors may share a container but never a selector identity"
        ),
        "canonical_selector_order": [
            "RK4-2049", "RK4-4097", "RK4-8193", "SSPRK3-4097",
            "SSPRK3-8193", "SSPRK3-16385",
        ],
        "restart_coordinate_time": "23/16", "active_target_time": "3/2",
    }
    if ({key: value for key, value in contract.items() if key != "source"} != expected_header
            or not isinstance(sources, list) or len(sources) != 2):
        _fail(config, "invalid_or_nonconverged_binder", "source contract grammar differs")
    parsed: dict[str, Mapping[str, Any]] = {}
    for item in sources:
        if not isinstance(item, Mapping) or item.get("source_id") not in _EXPECTED_SOURCES:
            _fail(config, "invalid_or_nonconverged_binder", "source record identity differs")
        source = str(item["source_id"])
        expected = {"source_id": source, **_EXPECTED_SOURCES[source]}
        if dict(item) != expected or source in parsed:
            _fail(config, "invalid_or_nonconverged_binder", f"{source} frozen source record differs")
        parsed[source] = item
    if tuple(parsed) != ("PROTO12", "RSP2"):
        _fail(config, "invalid_or_nonconverged_binder", "source order differs")
    return parsed


def validate_pref26_contract(config: Mapping[str, Any]) -> None:
    """Validate the operation-critical PREF26 contract before raw access."""
    if not isinstance(config, Mapping):
        raise Proto18Pref26BinderError(
            "PREF26_INVALID_OR_NONCONVERGED_BINDER", "configuration is not a mapping",
        )
    if (config.get("schema_version"), config.get("artifact_id"), config.get("project_version"),
            config.get("target_protocol")) != (1, "FGC-1-PRO18-PREF26", "0.11.0", "FGC-2-SF1-PROTO17"):
        _fail(config, "invalid_or_nonconverged_binder", "artifact identity differs")
    scope = config.get("scope")
    if (not isinstance(scope, Mapping)
            or scope.get("PREF26_read_only_source_execution_authorized") is not True
            or scope.get("raw_source_checkpoint_opening_permitted_only_for_declared_read_only_verification") is not True
            or scope.get("historical_event_log_opening_permitted_only_for_declared_read_only_semantic_replay") is not True
            or any(scope.get(key) is not False for key in (
                "future_PROTO17_output_root_observation", "future_namespace_created",
                "pretrajectory_operation", "trajectory_read", "mechanism_question_answered",
            ))):
        _fail(config, "authority_or_launch_promotion", "read-only scope differs")
    _source_contract(config)
    if config.get("resource_bounds") != _EXPECTED_RESOURCE_BOUNDS:
        _fail(config, "resource_bound_exceeded", "resource contract differs")
    if config.get("typed_stops") != _EXPECTED_STOPS:
        _fail(config, "invalid_or_nonconverged_binder", "typed-stop vocabulary differs")
    replay = config.get("semantic_replay")
    if (not isinstance(replay, Mapping)
            or replay.get("PROTO14_must_not_be_used_as_source") is not True
            or replay.get("full_prior_PDE_step_reexecution_claimed") is not False
            or any(replay.get(key) is not True for key in (
                "source_specific_legacy_grammar_required", "append_only_order_required",
                "terminal_identity_required", "checkpoint_embedded_history_linkage_required",
                "semantic_replay_must_not_normalize_or_rewrite_history",
            ))):
        _fail(config, "legacy_history_malformed_or_semantically_inconsistent", "semantic replay contract differs")


def _enforce_raw_contract(
    config: Mapping[str, Any], corpus: VerifiedProductionCorpus,
) -> dict[str, Any]:
    expected = _source_contract(config)
    records: list[dict[str, Any]] = []
    total_archive = total_log = total_members = total_uncompressed = total_lines = 0
    largest = 0
    for source in ("PROTO12", "RSP2"):
        archive = corpus.archives[source]
        frozen = expected[source]
        actual = {
            "source_id": source,
            "checkpoint": archive.spec.archive_path,
            "checkpoint_sha256": archive.archive_snapshot.sha256,
            "checkpoint_byte_count": archive.archive_snapshot.size,
            "checkpoint_npz_member_count": len(archive.archive_member_names),
            "checkpoint_total_uncompressed_byte_count": archive.total_uncompressed_bytes,
            "checkpoint_largest_member_uncompressed_byte_count": archive.largest_uncompressed_member_bytes,
            "event_log": archive.spec.event_log_path,
            "event_log_sha256": archive.event_snapshot.sha256,
            "event_log_byte_count": archive.event_snapshot.size,
            "event_log_line_count": len(archive.event_snapshot.payload[:-1].split(b"\n")) if archive.event_snapshot.payload.endswith(b"\n") else -1,
            "legacy_terminal_classification": frozen["legacy_terminal_classification"],
            "selector_keys": [item.key for item in archive.spec.selectors],
        }
        if actual != dict(frozen):
            owner = "source_byte_or_hash_drift"
            if (actual["checkpoint_npz_member_count"], actual["checkpoint_total_uncompressed_byte_count"],
                    actual["checkpoint_largest_member_uncompressed_byte_count"], actual["event_log_line_count"]) != (
                        frozen["checkpoint_npz_member_count"], frozen["checkpoint_total_uncompressed_byte_count"],
                        frozen["checkpoint_largest_member_uncompressed_byte_count"], frozen["event_log_line_count"],
                    ):
                owner = "resource_bound_exceeded"
            _fail(config, owner, f"{source} raw source identity differs")
        total_archive += actual["checkpoint_byte_count"]
        total_log += actual["event_log_byte_count"]
        total_members += actual["checkpoint_npz_member_count"]
        total_uncompressed += actual["checkpoint_total_uncompressed_byte_count"]
        total_lines += actual["event_log_line_count"]
        largest = max(largest, actual["checkpoint_largest_member_uncompressed_byte_count"])
        records.append({
            **actual,
            "metadata_raw_sha256": archive.metadata_sha256,
            "embedded_event_log_sha256": archive.embedded_event_log_sha256,
            "embedded_event_log_bytes_equal_external": (
                archive.embedded_event_log_payload == archive.event_snapshot.payload
            ),
            "archive_member_names": list(archive.archive_member_names),
            "archive_member_set_sha256": digest(list(archive.archive_member_names)),
        })
    totals = {
        "source_container_bytes": total_archive,
        "legacy_event_log_bytes": total_log,
        "total_source_bytes": total_archive + total_log,
        "NPZ_members": total_members,
        "uncompressed_NPZ_bytes": total_uncompressed,
        "largest_uncompressed_NPZ_member_bytes": largest,
        "legacy_event_log_lines": total_lines,
    }
    bounds = config["resource_bounds"]
    comparisons = {
        "source_container_bytes": bounds["maximum_total_source_container_bytes"],
        "legacy_event_log_bytes": bounds["maximum_total_legacy_event_log_bytes"],
        "total_source_bytes": bounds["maximum_total_source_bytes"],
        "NPZ_members": bounds["maximum_total_NPZ_members"],
        "uncompressed_NPZ_bytes": bounds["maximum_total_uncompressed_NPZ_bytes"],
        "largest_uncompressed_NPZ_member_bytes": bounds["maximum_single_uncompressed_NPZ_member_bytes"],
        "legacy_event_log_lines": bounds["maximum_legacy_event_log_lines"],
    }
    if totals != comparisons:
        _fail(config, "resource_bound_exceeded", "aggregate source resource identity differs")
    manifest = {"sources": records, "resource_totals": totals}
    return {**manifest, "source_archive_manifest_sha256": digest(manifest)}


def _descriptor(member: Any) -> dict[str, Any]:
    selector = member.selector
    source_bundle_key = f"{selector.source}:{selector.key}"
    shapes = {
        "u": [selector.point_count, 6], "p": [selector.point_count, 6],
        "q": [selector.point_count, 6], "tracer_positions": [48],
        "tracer_proper_times": [48], "event_proper_times": [24, 48],
        "event_fields": [24, 48, 6],
    }
    arrays = [
        {
            "logical_name": name, "source_bundle_key": source_bundle_key,
            "npz_storage_key": f"{selector.array_prefix}_{name}",
            "dtype": "<f8", "layout": "C", "shape": shapes[name],
            "little_endian_c_bytes_sha256": member.array_sha256[name],
        }
        for name in ARRAY_SUFFIXES
    ]
    accepted = time_identity(23 / 16)
    state = {
        "schema_id": "FGC-2-SF1-PROTO17-state-object-v1",
        "member_key": selector.key, "method": selector.method,
        "point_count": selector.point_count, "coordinate_time": accepted,
        "accepted_boundary_time": accepted, "step_index": selector.step_index,
        "transaction_serial": selector.transaction_serial,
        "input_hash": selector.input_sha256,
        "source_retry_count": selector.source_retry_count,
        "CFL_retry_count": selector.cfl_retry_count,
        "runtime_monitor_state": dict(member.metadata["runtime_monitor_state"]),
        "causal_state": dict(member.metadata["causal_state"]),
        "tracer_state": {"tracer_count": selector.tracer_count},
        "event_history_state": {"sample_count": selector.event_sample_count},
        "physical_arrays": arrays,
        "physical_state_sha256": member.physical_state_sha256,
        "restart_payload_sha256": member.restart_payload_sha256,
    }
    ledger = validate_ledger({
        "last_accepted_time_hex": accepted["binary64_hex"],
        "accepted_macro_step_count": 0,
        "cumulative_temporal_retry_count": 0,
        "current_macro_step_temporal_retry_count": 0,
        "last_accepted_macro_step_temporal_retry_count": 0,
        "accumulated_debit_vector_hex": ["0x0.0p+0"] * 18,
        "serialized_temporal_rejections": [],
    }, accepted_time=accepted)
    return {
        "member_key": selector.key, "method": selector.method,
        "point_count": selector.point_count, "source_checkpoint": selector.source,
        "source_bundle_key": source_bundle_key, "accepted_boundary_time": accepted,
        "step_index": selector.step_index,
        "transaction_serial": selector.transaction_serial,
        "state_object": state, "state_object_canonical_sha256": digest(state),
        "initial_TDG6_ledger": ledger,
        "initial_TDG6_ledger_sha256": digest(ledger),
    }


def _auth1_input_evidence(
    corpus: VerifiedProductionCorpus, source_manifest: Mapping[str, Any],
    replay: Mapping[str, Any], config: Mapping[str, Any],
) -> dict[str, Any]:
    descriptors = [_descriptor(corpus.members[key]) for key in corpus.members]
    bundles: list[dict[str, Any]] = []
    for descriptor in descriptors:
        source = descriptor["source_checkpoint"]
        member = corpus.members[descriptor["member_key"]]
        archive_path = corpus.archives[source].spec.archive_path
        prefix = "runs/fgc-2-sf1/"
        if not archive_path.startswith(prefix):
            _fail(config, "source_path_missing_or_unsafe", "source archive root differs")
        bundles.append({
            "source_bundle_key": descriptor["source_bundle_key"],
            "source_checkpoint": source,
            "relative_input_path": archive_path[len(prefix):],
            "raw_file_sha256": member.source_archive_sha256,
            "required_npz_member_keys": [
                item["npz_storage_key"]
                for item in descriptor["state_object"]["physical_arrays"]
            ] + ["metadata_utf8"],
        })
    physical_manifest = {
        "format": "FGC-2-SF1-PROTO17-physical-input-manifest-v1",
        "relative_input_root": "runs/fgc-2-sf1",
        "bundle_members": bundles,
    }
    bare = {
        "target_protocol": "FGC-2-SF1-PROTO17",
        "authority_identity": "FGC-1-PRO18-AUTH1",
        "restart_coordinate_time": time_identity(23 / 16),
        "active_event_target_time": time_identity(3 / 2),
        "common_event_index": 23,
        "member_keys_in_canonical_order": list(corpus.members),
        "member_descriptors": descriptors,
        "member_descriptor_set_sha256": digest(descriptors),
        "physical_input_manifest": physical_manifest,
        "physical_input_manifest_sha256": digest(physical_manifest),
        "source_archive_manifest_sha256": source_manifest["source_archive_manifest_sha256"],
        "historical_replay_sha256": digest(replay),
        "final_PROTO17_GenesisSpec_construction_deferred_to_AUTH1": True,
        "future_runtime_adapter_and_runner_source_hashes_present": False,
    }
    return {**bare, "AUTH1_input_evidence_sha256": digest(bare)}


def bind_pref26(
    repository_root: Path, pref26_config: Mapping[str, Any], pro18_config: Mapping[str, Any],
) -> Pref26Evidence:
    """Run the complete read-only two-container PREF26 evidence operation."""
    validate_pref26_contract(pref26_config)
    source_identities = _source_contract(pref26_config)
    try:
        corpus = verify_production_corpus(
            repository_root,
            pro18_config,
            source_identities=source_identities,
        )
    except Proto18ProductionPathError as error:
        _fail(pref26_config, "source_path_missing_or_unsafe", str(error))
    except Proto18ProductionSourceIdentityError as error:
        _fail(pref26_config, "source_byte_or_hash_drift", str(error))
    except Proto18ProductionResourceError as error:
        _fail(pref26_config, "resource_bound_exceeded", str(error))
    except Proto18ProductionSelectorError as error:
        _fail(pref26_config, "NPZ_schema_or_selector_drift", str(error))
    except Proto18ProductionArchiveError as error:
        _fail(pref26_config, "NPZ_schema_or_selector_drift", str(error))
    except Proto18ProductionInputError as error:
        _fail(pref26_config, "invalid_or_nonconverged_binder", str(error))
    source_manifest = _enforce_raw_contract(pref26_config, corpus)
    hashes = {key: member.physical_state_sha256 for key, member in corpus.members.items()}
    try:
        replay = {
            "PROTO12": replay_proto12_history(
                corpus.archives["PROTO12"].event_snapshot.payload,
                corpus.archives["PROTO12"].embedded_event_log_payload,
                corpus.archives["PROTO12"].metadata,
                {key: hashes[key] for key in hashes if key != "SSPRK3-16385"},
            ),
            "RSP2": replay_rsp2_history(
                corpus.archives["RSP2"].event_snapshot.payload,
                corpus.archives["RSP2"].embedded_event_log_payload,
                corpus.archives["RSP2"].metadata,
                hashes["SSPRK3-16385"],
            ),
        }
    except Proto18HistoricalReplayError as error:
        _fail(
            pref26_config,
            "legacy_history_malformed_or_semantically_inconsistent",
            str(error),
        )
    replay = {**replay, "historical_replay_sha256": digest(replay)}
    auth1_inputs = _auth1_input_evidence(
        corpus, source_manifest, replay, pref26_config,
    )
    duties = {
        "raw_bundle_byte_and_hash_recomputation": {"required": True, "completed": True},
        "external_Git_authority_blob_and_live_import_validation": {"required": True, "completed": False},
        "semantic_historical_journal_payload_replay": {"required": True, "completed": True},
        "namespace_reuse_and_foreign_store_rejection": {"required": True, "completed": False},
    }
    return Pref26Evidence(source_manifest, replay, auth1_inputs, duties)


__all__ = [
    "Pref26Evidence", "Proto18Pref26BinderError", "bind_pref26",
    "validate_pref26_contract",
]
