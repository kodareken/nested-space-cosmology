"""Static completeness audit for the sealed PROTO16 transition contract.

This module does not invent a runtime completion.  It reports what the frozen
schema fixes and what it does not fix, so a later implementation cannot turn
omitted construction rules into implicit authority.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
import json
from typing import Any, Mapping


MEMBERS = (
    "RK4-2049",
    "RK4-4097",
    "RK4-8193",
    "SSPRK3-4097",
    "SSPRK3-8193",
    "SSPRK3-16385",
)
RECEIPT_FIELDS = (
    "protocol_artifact_id",
    "campaign_id",
    "previous_campaign_generation",
    "previous_checkpoint_sha256",
    "previous_committed_common_event_index",
    "completed_common_event_index",
    "completed_event_time",
    "next_event_target_time",
    "member_keys",
    "prior_cursor_set_sha256",
    "prior_ledger_set_sha256",
    "prior_accepted_state_set_sha256",
)
SUCCESSOR_PRESENT = (
    "new_checkpoint_binds_derived_successor_cursor_set_after_receipt_hash_exists",
    "new_cursor_generation_and_attempt_serial_are_prior_plus_one",
    "new_attempt_id_uses_completed_common_event_index",
    "new_cursor_parent_is_prior_cursor_hash",
    "new_cursor_journal_tip_is_commit_record_hash",
)
SUCCESSOR_MISSING = (
    "complete_byte_for_byte_successor_cursor_copy_set",
    "attempt_high_water_map_rule",
    "cursor_high_water_map_rule",
    "exact_checkpoint_payload_construction",
    "exact_successor_cursor_set_byte_equality",
)
GENESIS_PRESENT = (
    "genesis_spec_required_fields",
    "member_descriptor_required_fields",
    "genesis_checkpoint_is_derivable_without_namespace_generated_value",
    "genesis_spec_and_checkpoint_must_not_embed_genesis_spec_sha256",
)
GENESIS_MISSING = (
    "complete_initial_TDG6_ledger_mapping_bytes",
    "canonical_state_object_descriptor_bytes",
    "cursor_origin_construction_rule",
    "pure_generation_zero_checkpoint_algorithm",
    "derived_checkpoint_byte_equality_rule",
)


def canonical(value: Any) -> bytes:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False
    ).encode()


def digest(value: Any) -> str:
    return sha256(canonical(value)).hexdigest()


def time_identity(value: Fraction) -> dict[str, str]:
    return {
        "rational": f"{value.numerator}/{value.denominator}",
        "binary64_hex": float(value).hex(),
    }


def receipt_fixture() -> dict[str, Any]:
    """Construct only the acyclic outer journal record shape the schema fixes."""
    payload = {
        "protocol_artifact_id": "FGC-2-SF1-PROTO16",
        "campaign_id": "fixture-campaign",
        "previous_campaign_generation": 4,
        "previous_checkpoint_sha256": "c" * 64,
        "previous_committed_common_event_index": 23,
        "completed_common_event_index": 24,
        "completed_event_time": time_identity(Fraction(3, 2)),
        "next_event_target_time": time_identity(Fraction(25, 16)),
        "member_keys": list(MEMBERS),
        "prior_cursor_set_sha256": "a" * 64,
        "prior_ledger_set_sha256": "b" * 64,
        "prior_accepted_state_set_sha256": "d" * 64,
        "sequence": 1,
    }
    record = {
        "record_kind": "COMMON_EVENT_COMMIT",
        "journal_parent_sha256": "e" * 64,
        "payload": payload,
    }
    record["record_sha256"] = digest(record)
    return record


def audit(protocol: Mapping[str, Any]) -> dict[str, Any]:
    """Audit the actual sealed TOML mapping rather than hardcoded assumptions."""
    event = protocol["event_commit"]
    genesis = protocol["trusted_genesis"]
    if event["receipt_payload_required_fields"] != list(RECEIPT_FIELDS):
        raise ValueError("sealed receipt field list differs")
    if any(event.get(key) is not True for key in SUCCESSOR_PRESENT):
        raise ValueError("sealed successor rule differs")
    if any(key in event for key in SUCCESSOR_MISSING):
        raise ValueError("a claimed missing successor rule is present")
    if any(key not in genesis for key in GENESIS_PRESENT):
        raise ValueError("sealed GenesisSpec rule differs")
    if any(key in genesis for key in GENESIS_MISSING):
        raise ValueError("a claimed missing GenesisSpec rule is present")
    record = receipt_fixture()
    payload = record["payload"]
    payload_without_sequence = tuple(key for key in payload if key != "sequence")
    acyclic = (
        payload_without_sequence == RECEIPT_FIELDS
        and not any("successor" in key or key.startswith("new_") for key in payload)
        and record["record_sha256"]
        == digest(
            {key: value for key, value in record.items() if key != "record_sha256"}
        )
    )
    return {
        "receipt": {
            "record_sha256": record["record_sha256"],
            "payload_sequence": payload["sequence"],
            "receipt_acyclic": acyclic,
            "present_payload_fields": list(RECEIPT_FIELDS),
            "sealed_payload_fields_match": True,
        },
        "successor_schema": {
            "present_frozen_rules": list(SUCCESSOR_PRESENT),
            "missing_exact_rules": list(SUCCESSOR_MISSING),
            "exact_derivation_frozen": False,
        },
        "genesis_schema": {
            "present_frozen_rules": list(GENESIS_PRESENT),
            "missing_exact_rules": list(GENESIS_MISSING),
            "exact_derivation_frozen": False,
        },
        "conclusion": {
            "PROTO16_receipt_hash_cycle_removed": acyclic,
            "PROTO16_exact_successor_and_genesis_derivation_not_frozen": True,
            "PROTO17_successor_freeze_design_authorized": True,
            "HLT14_implementation_and_synthetic_qualification_authorized": False,
        },
    }
