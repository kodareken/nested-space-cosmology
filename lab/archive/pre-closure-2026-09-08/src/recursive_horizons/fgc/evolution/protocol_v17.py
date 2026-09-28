"""Pure schema constants for the prospective PROTO17 reference oracle.

This module deliberately owns no configuration file, runtime store, or
authorization decision.  It fixes the input/output vocabulary used by the
side-effect-free construction oracle only.
"""

from __future__ import annotations

from typing import Any, Mapping


SF1_PROTOCOL_V17_ARTIFACT_ID = "FGC-2-SF1-PROTO17"
ROOT_SHA256 = "0" * 64
MEMBER_KEYS = (
    "RK4-2049", "RK4-4097", "RK4-8193",
    "SSPRK3-4097", "SSPRK3-8193", "SSPRK3-16385",
)
CADENCE_RATIONAL = "1/16"

LEDGER_FIELDS = (
    "last_accepted_time_hex",
    "accepted_macro_step_count",
    "cumulative_temporal_retry_count",
    "current_macro_step_temporal_retry_count",
    "last_accepted_macro_step_temporal_retry_count",
    "accumulated_debit_vector_hex",
    "serialized_temporal_rejections",
)
CURSOR_FIELDS = (
    "protocol_artifact_id", "campaign_id", "member_key", "method",
    "point_count", "committed_common_event_index", "accepted_boundary_time",
    "accepted_state_sha256", "previous_step_index", "previous_transaction_serial",
    "TDG6_ledger_sha256", "cursor_generation", "event_target_time",
    "attempt_serial", "attempt_id", "mode", "retry_successor_payload_or_none",
    "journal_tip_sha256", "cursor_chain_parent_sha256", "cursor_chain_sha256",
)
CHECKPOINT_FIELDS = (
    "protocol_artifact_id", "campaign_id", "campaign_generation",
    "committed_common_event_index", "active_event_target_time", "disposition",
    "terminal_member_key", "terminal_record_sha256", "member_keys", "cursors",
    "ledgers", "states", "cursor_set_sha256", "ledger_set_sha256",
    "accepted_state_set_sha256", "journal_tip_sha256", "terminal_lock",
    "attempt_high_water", "cursor_generation_high_water",
    "parent_checkpoint_sha256", "terminal_closure", "checkpoint_sha256",
)
RECEIPT_FIELDS = (
    "protocol_artifact_id", "campaign_id", "previous_campaign_generation",
    "previous_checkpoint_sha256", "previous_committed_common_event_index",
    "completed_common_event_index", "completed_event_time", "next_event_target_time",
    "member_keys", "prior_cursor_set_sha256", "prior_ledger_set_sha256",
    "prior_accepted_state_set_sha256",
)

# These are the visible, non-derived values needed to reconstruct generation
# zero.  The materialized spec adds checksums only after the checkpoint exists.
GENESIS_INPUT_FIELDS = (
    "genesis_algorithm_id", "canonical_json_encoding", "state_object_schema_id",
    "tdg6_ledger_schema_id", "cursor_schema_id", "checkpoint_schema_id",
    "protocol_artifact_id", "protocol_config_sha256", "protocol_freeze_result_sha256",
    "hlt13_result_sha256", "run_plan_sha256", "runtime_module_sha256",
    "adapter_module_sha256", "runner_sha256", "numerical_environment",
    "campaign_id", "namespace", "journal_root_parent_sha256", "common_event_index",
    "restart_coordinate_time", "active_event_target_time", "member_keys_in_canonical_order",
    "physical_input_manifest", "physical_input_manifest_sha256", "member_descriptors",
    "member_descriptor_set_sha256", "genesis_checkpoint_sha256", "cursor_set_sha256",
    "ledger_set_sha256", "accepted_state_set_sha256",
)
MEMBER_DESCRIPTOR_FIELDS = (
    "member_key", "method", "point_count", "source_checkpoint", "source_bundle_key",
    "accepted_boundary_time", "step_index", "transaction_serial", "state_object",
    "state_object_canonical_sha256", "initial_TDG6_ledger", "initial_TDG6_ledger_sha256",
)
STATE_OBJECT_FIELDS = (
    "schema_id", "member_key", "method", "point_count", "coordinate_time",
    "accepted_boundary_time", "step_index", "transaction_serial", "input_hash", "source_retry_count",
    "CFL_retry_count", "runtime_monitor_state", "causal_state", "tracer_state",
    "event_history_state", "physical_arrays", "physical_state_sha256", "restart_payload_sha256",
)
ARRAY_DESCRIPTOR_FIELDS = (
    "logical_name", "source_bundle_key", "npz_storage_key", "dtype", "layout",
    "shape", "little_endian_c_bytes_sha256",
)
STATE_ARRAY_NAMES = (
    "u", "p", "q", "tracer_positions", "tracer_proper_times",
    "event_proper_times", "event_fields",
)
DERIVED_GENESIS_FIELDS = (
    "member_descriptor_set_sha256", "genesis_checkpoint_sha256",
    "cursor_set_sha256", "ledger_set_sha256", "accepted_state_set_sha256",
)


def validate_proto17_schema(value: Mapping[str, Any]) -> dict[str, Any]:
    """Validate a minimal, non-authorizing PROTO17 schema declaration.

    This deliberately accepts no location, runtime, or outcome fields.  Those
    belong to a later freeze/integration surface, not the pure oracle.
    """
    if not isinstance(value, Mapping):
        raise TypeError("PROTO17 schema must be a mapping")
    expected = {
        "schema_version", "artifact_id", "protocol_version", "frozen",
        "member_keys", "cadence", "genesis_input_fields",
        "member_descriptor_fields", "receipt_fields", "claims",
    }
    if set(value) != expected:
        raise ValueError("PROTO17 schema keys differ")
    if (
        value["schema_version"] != 1
        or value["artifact_id"] != SF1_PROTOCOL_V17_ARTIFACT_ID
        or value["protocol_version"] != 17
        or value["frozen"] is not True
        or tuple(value["member_keys"]) != MEMBER_KEYS
        or value["cadence"] != CADENCE_RATIONAL
        or tuple(value["genesis_input_fields"]) != GENESIS_INPUT_FIELDS
        or tuple(value["member_descriptor_fields"]) != MEMBER_DESCRIPTOR_FIELDS
        or tuple(value["receipt_fields"]) != RECEIPT_FIELDS
    ):
        raise ValueError("PROTO17 schema identity differs")
    claims = value["claims"]
    if not isinstance(claims, Mapping) or set(claims) != {
        "PROTO17_pure_reference_oracle_implemented",
        "PROTO17_runtime_implemented",
        "PROTO17_namespace_authorized",
        "PROTO17_trajectory_read",
    }:
        raise ValueError("PROTO17 pure-oracle claims differ")
    if claims["PROTO17_pure_reference_oracle_implemented"] is not True or any(
        claims[key] is not False for key in claims if key != "PROTO17_pure_reference_oracle_implemented"
    ):
        raise ValueError("PROTO17 pure-oracle decision boundary differs")
    return dict(value)


__all__ = [
    "CADENCE_RATIONAL", "CHECKPOINT_FIELDS", "CURSOR_FIELDS",
    "DERIVED_GENESIS_FIELDS", "GENESIS_INPUT_FIELDS", "LEDGER_FIELDS",
    "MEMBER_DESCRIPTOR_FIELDS", "MEMBER_KEYS", "RECEIPT_FIELDS", "ROOT_SHA256",
    "STATE_ARRAY_NAMES", "STATE_OBJECT_FIELDS", "ARRAY_DESCRIPTOR_FIELDS",
    "SF1_PROTOCOL_V17_ARTIFACT_ID", "validate_proto17_schema",
]
