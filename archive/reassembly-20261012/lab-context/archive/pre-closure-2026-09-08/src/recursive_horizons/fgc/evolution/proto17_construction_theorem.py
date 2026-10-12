"""Independent structural rederivation for the PROTO17 construction boundary.

This module deliberately does not import the PROTO17 protocol or oracle.  It
parses supplied TOML *bytes*, checks the narrow frozen vocabulary independently,
then derives a synthetic generation-zero checkpoint and one successful common
event edge.  It owns no filesystem, namespace, runtime, or numerical state.
"""

from __future__ import annotations

from copy import deepcopy
from fractions import Fraction
from hashlib import sha256
import json
import math
import tomllib
from typing import Any, Mapping


ROOT = "0" * 64
ARTIFACT = "FGC-2-SF1-PROTO17"
MEMBERS = (
    "RK4-2049", "RK4-4097", "RK4-8193",
    "SSPRK3-4097", "SSPRK3-8193", "SSPRK3-16385",
)
CADENCE = Fraction(1, 16)
LEDGER_FIELDS = (
    "last_accepted_time_hex", "accepted_macro_step_count",
    "cumulative_temporal_retry_count", "current_macro_step_temporal_retry_count",
    "last_accepted_macro_step_temporal_retry_count", "accumulated_debit_vector_hex",
    "serialized_temporal_rejections",
)
STATE_FIELDS = (
    "schema_id", "member_key", "method", "point_count", "coordinate_time",
    "accepted_boundary_time", "step_index", "transaction_serial", "input_hash",
    "source_retry_count", "CFL_retry_count", "runtime_monitor_state", "causal_state",
    "tracer_state", "event_history_state", "physical_arrays", "physical_state_sha256",
    "restart_payload_sha256",
)
ARRAY_FIELDS = (
    "logical_name", "source_bundle_key", "npz_storage_key", "dtype", "layout",
    "shape", "little_endian_c_bytes_sha256",
)
ARRAY_NAMES = (
    "u", "p", "q", "tracer_positions", "tracer_proper_times",
    "event_proper_times", "event_fields",
)
STATE_OBJECT_SCHEMA_ID = "FGC-2-SF1-PROTO17-state-object-v1"
GENESIS_FIELDS = (
    "genesis_algorithm_id", "canonical_json_encoding", "state_object_schema_id",
    "tdg6_ledger_schema_id", "cursor_schema_id", "checkpoint_schema_id",
    "protocol_artifact_id", "protocol_config_sha256", "protocol_freeze_result_sha256",
    "hlt13_result_sha256", "run_plan_sha256", "runtime_module_sha256",
    "adapter_module_sha256", "runner_sha256", "numerical_environment", "campaign_id",
    "namespace", "journal_root_parent_sha256", "common_event_index", "restart_coordinate_time",
    "active_event_target_time", "member_keys_in_canonical_order", "physical_input_manifest",
    "physical_input_manifest_sha256", "member_descriptors", "member_descriptor_set_sha256",
    "genesis_checkpoint_sha256", "cursor_set_sha256", "ledger_set_sha256",
    "accepted_state_set_sha256",
)
MEMBER_DESCRIPTOR_FIELDS = (
    "member_key", "method", "point_count", "source_checkpoint", "source_bundle_key",
    "accepted_boundary_time", "step_index", "transaction_serial", "state_object",
    "state_object_canonical_sha256", "initial_TDG6_ledger", "initial_TDG6_ledger_sha256",
)
CHECKPOINT_FIELDS = (
    "protocol_artifact_id", "campaign_id", "campaign_generation", "committed_common_event_index",
    "active_event_target_time", "disposition", "terminal_member_key", "terminal_record_sha256",
    "member_keys", "cursors", "ledgers", "states", "cursor_set_sha256", "ledger_set_sha256",
    "accepted_state_set_sha256", "journal_tip_sha256", "terminal_lock", "attempt_high_water",
    "cursor_generation_high_water", "parent_checkpoint_sha256", "terminal_closure", "checkpoint_sha256",
)
SEALED_SYNTHETIC_HASHES = {
    "genesis_checkpoint_sha256": "01e8214dcbbef025fe7883670043659d3ec66bb47e9b3e9dd4b8089b8ba5f8de",
    "common_event_record_sha256": "17f925e8489341392ac5df4fba1cc47b66a27153e64eca4455c3f701f37cbc17",
    "successor_checkpoint_sha256": "c4c31f53a6a97632addafd54d6b3ec08bdcffe4aa4321269f8c8bca957e5380b",
}


class Proto17TheoremError(ValueError):
    """A claimed PROTO17 construction premise is structurally false."""


def canonical(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
                      allow_nan=False).encode("utf-8")


def digest(value: object) -> str:
    return sha256(canonical(value)).hexdigest()


def _sha(value: object, name: str) -> str:
    if not isinstance(value, str) or len(value) != 64 or value.lower() != value:
        raise Proto17TheoremError(f"{name} is not a lowercase SHA-256 digest")
    try:
        int(value, 16)
    except ValueError as error:
        raise Proto17TheoremError(f"{name} is not hexadecimal") from error
    return value


def _integer(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise Proto17TheoremError(f"{name} is not a nonnegative integer")
    return value


def _time(value: object) -> dict[str, str]:
    if not isinstance(value, Mapping) or set(value) != {"rational", "binary64_hex"}:
        raise Proto17TheoremError("time identity fields differ")
    if not isinstance(value["rational"], str) or not isinstance(value["binary64_hex"], str):
        raise Proto17TheoremError("time identity values differ")
    try:
        rational = Fraction(value["rational"])
        floating = float.fromhex(value["binary64_hex"])
    except (ValueError, ZeroDivisionError) as error:
        raise Proto17TheoremError("time identity is malformed") from error
    if not math.isfinite(floating) or str(rational) != value["rational"] or float(rational).hex() != floating.hex():
        raise Proto17TheoremError("time rational/binary64 identity differs")
    return {"rational": str(rational), "binary64_hex": floating.hex()}


def _time_of(rational: Fraction) -> dict[str, str]:
    return {"rational": str(rational), "binary64_hex": float(rational).hex()}


def _method_grid(key: str) -> tuple[str, int]:
    if key not in MEMBERS:
        raise Proto17TheoremError("member key differs from frozen six")
    method, points = key.split("-", 1)
    return method, int(points)


def _array_shapes(points: int) -> dict[str, list[int]]:
    return {
        "u": [points, 6],
        "p": [points, 6],
        "q": [points, 6],
        "tracer_positions": [48],
        "tracer_proper_times": [48],
        "event_proper_times": [24, 48],
        "event_fields": [24, 48, 6],
    }


def _validate_visible_state_object(
    value: Mapping[str, object], *, member_key: str, source_bundle_key: str
) -> dict[str, object]:
    """Validate the complete frozen state-object vocabulary without opening bytes."""
    if not isinstance(value, Mapping) or set(value) != set(STATE_FIELDS):
        raise Proto17TheoremError("state-object fields differ")
    state = dict(value)
    method, points = _method_grid(member_key)
    if (state["schema_id"], state["member_key"], state["method"], state["point_count"]) != (
        STATE_OBJECT_SCHEMA_ID, member_key, method, points,
    ):
        raise Proto17TheoremError("state-object identity differs")
    state["coordinate_time"] = _time(state["coordinate_time"])
    state["accepted_boundary_time"] = _time(state["accepted_boundary_time"])
    for field in ("step_index", "transaction_serial", "source_retry_count", "CFL_retry_count"):
        _integer(state[field], f"state-object {field}")
    for field in ("input_hash", "physical_state_sha256", "restart_payload_sha256"):
        _sha(state[field], f"state-object {field}")
    for field in ("runtime_monitor_state", "causal_state", "tracer_state", "event_history_state"):
        if not isinstance(state[field], Mapping):
            raise Proto17TheoremError(f"state-object {field} differs")
    arrays = state["physical_arrays"]
    if not isinstance(arrays, list) or len(arrays) != len(ARRAY_NAMES):
        raise Proto17TheoremError("state-object physical-array set differs")
    expected_shapes = _array_shapes(points)
    normalized_arrays: list[dict[str, object]] = []
    for name, supplied in zip(ARRAY_NAMES, arrays, strict=True):
        if not isinstance(supplied, Mapping) or set(supplied) != set(ARRAY_FIELDS):
            raise Proto17TheoremError("physical-array fields differ")
        item = dict(supplied)
        if (item["logical_name"], item["source_bundle_key"], item["dtype"], item["layout"], item["shape"]) != (
            name, source_bundle_key, "<f8", "C", expected_shapes[name],
        ):
            raise Proto17TheoremError("physical-array visible identity differs")
        if not isinstance(item["npz_storage_key"], str) or not item["npz_storage_key"]:
            raise Proto17TheoremError("physical-array storage key differs")
        _sha(item["little_endian_c_bytes_sha256"], "physical-array byte digest")
        normalized_arrays.append(item)
    state["physical_arrays"] = normalized_arrays
    return state


def parse_frozen_contract(protocol_bytes: bytes, freeze_bytes: bytes) -> dict[str, object]:
    """Independently parse only the theorem-relevant frozen configuration bytes."""
    try:
        protocol = tomllib.loads(protocol_bytes.decode("utf-8"))
        freeze = tomllib.loads(freeze_bytes.decode("utf-8"))
    except (UnicodeDecodeError, tomllib.TOMLDecodeError) as error:
        raise Proto17TheoremError("PROTO17 TOML bytes are malformed") from error
    if protocol.get("artifact_id") != ARTIFACT or protocol.get("protocol_version") != 17 or protocol.get("frozen") is not True:
        raise Proto17TheoremError("PROTO17 protocol identity differs")
    if freeze.get("artifact_id") != "FGC-1-PRO17-FRZ1" or freeze.get("claims") != protocol.get("claims"):
        raise Proto17TheoremError("PROTO17 freeze identity/claims differ")
    restart = protocol.get("restart_inputs")
    successor = protocol.get("common_event_successor")
    genesis = protocol.get("trusted_genesis")
    if not isinstance(restart, Mapping) or not isinstance(successor, Mapping) or not isinstance(genesis, Mapping):
        raise Proto17TheoremError("PROTO17 theorem sections differ")
    if (restart.get("restart_coordinate_time"), restart.get("genesis_active_event_target_time"),
            tuple(restart.get("member_keys_in_canonical_order", ()))) != ("23/16", "3/2", MEMBERS):
        raise Proto17TheoremError("PROTO17 restart contract differs")
    required_successor = (
        "receipt_is_predecessor_only", "receipt_self_successor_and_checkpoint_hashes_forbidden",
        "outer_record_sha256_must_exist_before_successor_construction",
        "successor_cursor_preserves_state_ledger_time_method_grid_and_accepted_serials",
        "successor_cursor_generation_and_attempt_serial_are_prior_plus_one",
        "successor_cursor_parent_is_prior_cursor_hash", "successor_cursor_journal_tip_is_outer_record_hash",
        "successor_high_water_maps_equal_new_cursor_values",
        "successor_checkpoint_preserves_complete_state_and_ledger_sets",
        "only_one_lone_durable_common_event_suffix_is_reconstructible",
        "all_other_uncheckpointed_suffixes_are_invalid",
        "receipt_sequence_is_hash_contiguous_supplied_prior_journal_length_plus_one",
        "prior_journal_record_payload_semantics_and_replay_are_deferred_to_HLT14_runtime",
    )
    if any(successor.get(item) is not True for item in required_successor):
        raise Proto17TheoremError("PROTO17 successor theorem contract differs")
    if (tuple(genesis.get("genesis_spec_required_fields", ())) != GENESIS_FIELDS
            or tuple(genesis.get("member_descriptor_required_fields", ())) != MEMBER_DESCRIPTOR_FIELDS
            or tuple(genesis.get("state_object_required_fields", ())) != STATE_FIELDS
            or tuple(genesis.get("physical_array_required_fields", ())) != ARRAY_FIELDS
            or tuple(genesis.get("tdg6_ledger_required_fields", ())) != LEDGER_FIELDS):
        raise Proto17TheoremError("PROTO17 visible state/ledger vocabulary differs")
    if genesis.get("physical_array_raw_bytes_are_not_opened_or_recomputed_by_this_freeze") is not True:
        raise Proto17TheoremError("PROTO17 raw-byte boundary differs")
    return {"protocol_sha256": sha256(protocol_bytes).hexdigest(), "freeze_sha256": sha256(freeze_bytes).hexdigest(),
            "historical_record_payload_semantics_deferred": True, "raw_bytes_deferred": True}


def _ledger(time: Mapping[str, object], *, zero: bool = True) -> dict[str, object]:
    identity = _time(time)
    answer = {"last_accepted_time_hex": identity["binary64_hex"], "accepted_macro_step_count": 0,
              "cumulative_temporal_retry_count": 0, "current_macro_step_temporal_retry_count": 0,
              "last_accepted_macro_step_temporal_retry_count": 0,
              "accumulated_debit_vector_hex": ["0x0.0p+0"] * 18,
              "serialized_temporal_rejections": []}
    if not zero:
        answer["accepted_macro_step_count"] = 1
    return answer


def _validate_ledger(value: Mapping[str, object], accepted: Mapping[str, object]) -> dict[str, object]:
    if not isinstance(value, Mapping) or set(value) != set(LEDGER_FIELDS):
        raise Proto17TheoremError("ledger fields differ")
    result = dict(value)
    if result["last_accepted_time_hex"] != _time(accepted)["binary64_hex"]:
        raise Proto17TheoremError("ledger accepted time differs")
    for field in LEDGER_FIELDS[1:5]:
        _integer(result[field], f"ledger {field}")
    if (result["current_macro_step_temporal_retry_count"] > result["cumulative_temporal_retry_count"]
            or result["last_accepted_macro_step_temporal_retry_count"] > result["cumulative_temporal_retry_count"]):
        raise Proto17TheoremError("ledger retry counters differ")
    debit = result["accumulated_debit_vector_hex"]
    if not isinstance(debit, list) or len(debit) != 18:
        raise Proto17TheoremError("ledger debit vector differs")
    for item in debit:
        try:
            numeric = float.fromhex(item)
        except (TypeError, ValueError) as error:
            raise Proto17TheoremError("ledger debit is malformed") from error
        if not math.isfinite(numeric) or numeric < 0 or numeric.hex() != item:
            raise Proto17TheoremError("ledger debit differs")
    rejected = result["serialized_temporal_rejections"]
    if not isinstance(rejected, list) or len(rejected) != result["cumulative_temporal_retry_count"]:
        raise Proto17TheoremError("ledger rejection count differs")
    return result


def _cursor(payload: Mapping[str, object]) -> dict[str, object]:
    required = {"protocol_artifact_id", "campaign_id", "member_key", "method", "point_count",
                "committed_common_event_index", "accepted_boundary_time", "accepted_state_sha256",
                "previous_step_index", "previous_transaction_serial", "TDG6_ledger_sha256",
                "cursor_generation", "event_target_time", "attempt_serial", "attempt_id", "mode",
                "retry_successor_payload_or_none", "journal_tip_sha256", "cursor_chain_parent_sha256"}
    result = dict(payload)
    supplied = result.pop("cursor_chain_sha256", None)
    if set(result) != required:
        raise Proto17TheoremError("cursor fields differ")
    key = result["member_key"]; method, points = _method_grid(str(key))
    if result["protocol_artifact_id"] != ARTIFACT or result["method"] != method or result["point_count"] != points:
        raise Proto17TheoremError("cursor identity differs")
    if result["mode"] != "FRESH_READY" or result["retry_successor_payload_or_none"] is not None:
        raise Proto17TheoremError("cursor is not fresh ready")
    for field in ("committed_common_event_index", "previous_step_index", "previous_transaction_serial", "cursor_generation", "attempt_serial"):
        _integer(result[field], field)
    result["accepted_boundary_time"] = _time(result["accepted_boundary_time"])
    result["event_target_time"] = _time(result["event_target_time"])
    for field in ("accepted_state_sha256", "TDG6_ledger_sha256", "journal_tip_sha256", "cursor_chain_parent_sha256"):
        _sha(result[field], field)
    if result["attempt_id"] != f"{result['campaign_id']}:{key}:{result['committed_common_event_index']}:{result['attempt_serial']}":
        raise Proto17TheoremError("cursor attempt identity differs")
    result["cursor_chain_sha256"] = digest(result)
    if supplied is not None and _sha(supplied, "cursor chain") != result["cursor_chain_sha256"]:
        raise Proto17TheoremError("cursor hash differs")
    return result


def _sets(cursors: Mapping[str, Mapping[str, object]], ledgers: Mapping[str, Mapping[str, object]], states: Mapping[str, str]) -> tuple[str, str, str]:
    return (digest([{ "member_key": key, "cursor_chain_sha256": cursors[key]["cursor_chain_sha256"],
                     "TDG6_ledger_sha256": cursors[key]["TDG6_ledger_sha256"],
                     "accepted_state_sha256": states[key], "accepted_boundary_time": cursors[key]["accepted_boundary_time"],
                     "cursor_generation": cursors[key]["cursor_generation"], "committed_common_event_index": cursors[key]["committed_common_event_index"]} for key in MEMBERS]),
            digest({key: digest(ledgers[key]) for key in MEMBERS}), digest({key: states[key] for key in MEMBERS}))


def _checkpoint(*, generation: int, parent: str, event: int, target: Mapping[str, object], cursors: Mapping[str, Mapping[str, object]], ledgers: Mapping[str, Mapping[str, object]], states: Mapping[str, str], tip: str, high: Mapping[str, int]) -> dict[str, object]:
    _integer(generation, "generation"); _integer(event, "event"); _sha(parent, "checkpoint parent"); _sha(tip, "journal tip")
    if (generation == 0) != (parent == ROOT) or (generation == 0) != (tip == ROOT):
        raise Proto17TheoremError("generation/root relations differ")
    if set(cursors) != set(MEMBERS) or set(ledgers) != set(MEMBERS) or set(states) != set(MEMBERS) or set(high) != set(MEMBERS):
        raise Proto17TheoremError("checkpoint six-member sets differ")
    normalized_cursors: dict[str, dict[str, object]] = {}; normalized_ledgers: dict[str, dict[str, object]] = {}
    identity_target = _time(target)
    for key in MEMBERS:
        current = _cursor(cursors[key]); ledger = _validate_ledger(ledgers[key], current["accepted_boundary_time"])
        if current["committed_common_event_index"] != event or current["event_target_time"] != identity_target or current["accepted_state_sha256"] != _sha(states[key], "state address") or current["TDG6_ledger_sha256"] != digest(ledger) or current["attempt_serial"] != _integer(high[key], "high water") or current["cursor_generation"] != high[key]:
            raise Proto17TheoremError("checkpoint cursor/ledger/state binding differs")
        normalized_cursors[key] = current; normalized_ledgers[key] = ledger
    cursor_set, ledger_set, state_set = _sets(normalized_cursors, normalized_ledgers, states)
    bare = {"protocol_artifact_id": ARTIFACT, "campaign_id": normalized_cursors[MEMBERS[0]]["campaign_id"],
            "campaign_generation": generation, "committed_common_event_index": event, "active_event_target_time": identity_target,
            "disposition": "nonterminal", "terminal_member_key": None, "terminal_record_sha256": None,
            "member_keys": list(MEMBERS), "cursors": normalized_cursors, "ledgers": normalized_ledgers,
            "states": dict(states), "cursor_set_sha256": cursor_set, "ledger_set_sha256": ledger_set,
            "accepted_state_set_sha256": state_set, "journal_tip_sha256": tip, "terminal_lock": False,
            "attempt_high_water": dict(high), "cursor_generation_high_water": dict(high),
            "parent_checkpoint_sha256": parent, "terminal_closure": None}
    return {**bare, "checkpoint_sha256": digest(bare)}


def _record(kind: str, parent: str, payload: Mapping[str, object]) -> dict[str, object]:
    bare = {"record_kind": kind, "journal_parent_sha256": parent, "payload": dict(payload)}
    return {**bare, "record_sha256": digest(bare)}


def _validate_prior(records: list[Mapping[str, object]], tip: str) -> int:
    parent = ROOT
    allowed = {"TDG6_REJECTION", "CURSOR_TRANSITION", "ACCEPTED_FINE_COMMIT", "TERMINAL_INVALID", "COMMON_EVENT_ROLLBACK", "COMMON_EVENT_COMMIT"}
    for sequence, record in enumerate(records, 1):
        if not isinstance(record, Mapping) or set(record) != {"record_kind", "journal_parent_sha256", "payload", "record_sha256"} or record["record_kind"] not in allowed:
            raise Proto17TheoremError("prior journal envelope differs")
        if not isinstance(record["payload"], Mapping) or record["payload"].get("sequence") != sequence or record["journal_parent_sha256"] != parent:
            raise Proto17TheoremError("prior journal contiguity differs")
        bare = dict(record); observed = _sha(bare.pop("record_sha256"), "prior journal hash")
        if digest(bare) != observed:
            raise Proto17TheoremError("prior journal hash differs")
        parent = observed
    if parent != tip:
        raise Proto17TheoremError("prior journal does not end at checkpoint tip")
    return len(records) + 1


def generation_zero_fixture() -> dict[str, object]:
    """Return the exact checkpoint derived from the complete visible fixture."""
    return deepcopy(sealed_synthetic_genesis()["checkpoint"])


def sealed_synthetic_genesis() -> dict[str, object]:
    """Independently reproduce the complete sealed PROTO17 synthetic GenesisSpec.

    The fixture has no physical bytes: array identities are declared external
    inputs exactly as the freeze states.  Its shapes and descriptors are still
    complete enough to derive every generation-zero content address.
    """
    restart, target = _time_of(Fraction(23, 16)), _time_of(Fraction(3, 2))
    descriptors: list[dict[str, object]] = []; bundles: list[dict[str, object]] = []
    for key in MEMBERS:
        method, points = _method_grid(key)
        shapes = _array_shapes(points)
        arrays = [{"logical_name": name, "source_bundle_key": f"fixture-{key}",
                   "npz_storage_key": f"arr_{name}", "dtype": "<f8", "layout": "C",
                   "shape": shapes[name], "little_endian_c_bytes_sha256": digest({"fixture": key, "array": name})}
                  for name in ARRAY_NAMES]
        state = {"schema_id": "FGC-2-SF1-PROTO17-state-object-v1", "member_key": key,
                 "method": method, "point_count": points, "coordinate_time": restart,
                 "accepted_boundary_time": restart, "step_index": 1, "transaction_serial": 2,
                 "input_hash": digest({"fixture": key, "input": True}), "source_retry_count": 0,
                 "CFL_retry_count": 0, "runtime_monitor_state": {"accepted_stage_count": 1},
                 "causal_state": {"remaining_buffer": "1"}, "tracer_state": {"tracer_count": 48},
                 "event_history_state": {"sample_count": 24}, "physical_arrays": arrays,
                 "physical_state_sha256": digest({"fixture": key, "kind": "u-p-q"}),
                 "restart_payload_sha256": digest({"fixture": key, "restart": True})}
        state = _validate_visible_state_object(state, member_key=key, source_bundle_key=f"fixture-{key}")
        ledger = _ledger(restart)
        descriptor = {"member_key": key, "method": method, "point_count": points,
                      "source_checkpoint": "RSP2" if key == "SSPRK3-16385" else "PROTO12",
                      "source_bundle_key": f"fixture-{key}", "accepted_boundary_time": restart,
                      "step_index": 1, "transaction_serial": 2, "state_object": state,
                      "state_object_canonical_sha256": digest(state), "initial_TDG6_ledger": ledger,
                      "initial_TDG6_ledger_sha256": digest(ledger)}
        if set(descriptor) != set(MEMBER_DESCRIPTOR_FIELDS):
            raise AssertionError("sealed independent descriptor schema drift")
        descriptors.append(descriptor)
        bundles.append({"source_bundle_key": f"fixture-{key}", "source_checkpoint": descriptor["source_checkpoint"],
                        "relative_input_path": f"fixture/{key}.npz", "raw_file_sha256": digest({"fixture": key, "raw": True}),
                        "required_npz_member_keys": [item["npz_storage_key"] for item in arrays] + ["metadata_utf8"]})
    manifest = {"format": "FGC-2-SF1-PROTO17-physical-input-manifest-v1", "relative_input_root": "synthetic/no-read", "bundle_members": bundles}
    source = {"genesis_algorithm_id": "FGC-2-SF1-PROTO17-pure-construction-v1",
              "canonical_json_encoding": "utf-8;json-sort-keys;separators=(',',':');ensure-ascii=true;allow-nan=false",
              "state_object_schema_id": "FGC-2-SF1-PROTO17-state-object-v1", "tdg6_ledger_schema_id": "FGC-2-SF1-PROTO17-tdg6-ledger-v1",
              "cursor_schema_id": "FGC-2-SF1-PROTO17-cursor-v1", "checkpoint_schema_id": "FGC-2-SF1-PROTO17-checkpoint-v1",
              "protocol_artifact_id": ARTIFACT, "protocol_config_sha256": digest("protocol"), "protocol_freeze_result_sha256": digest("freeze"),
              "hlt13_result_sha256": digest("hlt13"), "run_plan_sha256": digest("plan"), "runtime_module_sha256": digest("runtime"),
              "adapter_module_sha256": digest("adapter"), "runner_sha256": digest("runner"), "numerical_environment": {"python": "fixture"},
              "campaign_id": "fixture-campaign", "namespace": "synthetic/no-write", "journal_root_parent_sha256": ROOT,
              "common_event_index": 23, "restart_coordinate_time": restart, "active_event_target_time": target,
              "member_keys_in_canonical_order": list(MEMBERS), "physical_input_manifest": manifest,
              "physical_input_manifest_sha256": digest(manifest), "member_descriptors": descriptors}
    states = {item["member_key"]: item["state_object_canonical_sha256"] for item in descriptors}; ledgers = {item["member_key"]: item["initial_TDG6_ledger"] for item in descriptors}; cursors = {}
    for item in descriptors:
        key, method, points = item["member_key"], item["method"], item["point_count"]
        cursors[key] = _cursor({"protocol_artifact_id": ARTIFACT, "campaign_id": "fixture-campaign", "member_key": key, "method": method, "point_count": points, "committed_common_event_index": 23, "accepted_boundary_time": restart, "accepted_state_sha256": states[key], "previous_step_index": 1, "previous_transaction_serial": 2, "TDG6_ledger_sha256": digest(ledgers[key]), "cursor_generation": 0, "event_target_time": target, "attempt_serial": 0, "attempt_id": f"fixture-campaign:{key}:23:0", "mode": "FRESH_READY", "retry_successor_payload_or_none": None, "journal_tip_sha256": ROOT, "cursor_chain_parent_sha256": ROOT})
    high = {key: 0 for key in MEMBERS}; checkpoint = _checkpoint(generation=0, parent=ROOT, event=23, target=target, cursors=cursors, ledgers=ledgers, states=states, tip=ROOT, high=high)
    materialized = {**source, "member_descriptor_set_sha256": digest(descriptors), "genesis_checkpoint_sha256": checkpoint["checkpoint_sha256"], "cursor_set_sha256": checkpoint["cursor_set_sha256"], "ledger_set_sha256": checkpoint["ledger_set_sha256"], "accepted_state_set_sha256": checkpoint["accepted_state_set_sha256"]}
    if set(materialized) != set(GENESIS_FIELDS):
        raise AssertionError("sealed independent GenesisSpec schema drift")
    store_plan = {"state_objects": {states[key]: descriptors[index]["state_object"] for index, key in enumerate(MEMBERS)}, "state_paths": {key: f"states/{states[key]}.json" for key in MEMBERS}, "journal_records": [], "checkpoint_path": f"checkpoints/00000000000000000000-{checkpoint['checkpoint_sha256']}.json", "checkpoint": checkpoint}
    return {"genesis_spec": materialized, "checkpoint": checkpoint, "store_plan": store_plan}


def validate_sealed_synthetic_fixture(value: Mapping[str, object]) -> dict[str, object]:
    """Validate the complete independent fixture by fresh reconstruction.

    This is intentionally a qualification-only exact fixture validator, not an
    authority for production GenesisSpec values.
    """
    expected = sealed_synthetic_genesis()
    if not isinstance(value, Mapping) or dict(value) != expected:
        raise Proto17TheoremError("complete independently derived synthetic GenesisSpec differs")
    return expected


def synchronized_predecessor() -> tuple[dict[str, object], list[dict[str, object]]]:
    origin = generation_zero_fixture(); target = origin["active_event_target_time"]; cursors = deepcopy(origin["cursors"]); ledgers = deepcopy(origin["ledgers"]); records = []; tip = ROOT
    for sequence, key in enumerate(MEMBERS, 1):
        ledgers[key] = _ledger(target, zero=False); old = cursors[key]; update = {**old, "accepted_boundary_time": target, "TDG6_ledger_sha256": digest(ledgers[key]), "cursor_generation": 1, "attempt_serial": 1, "attempt_id": f"{origin['campaign_id']}:{key}:23:1", "journal_tip_sha256": tip, "cursor_chain_parent_sha256": old["cursor_chain_sha256"]}; update.pop("cursor_chain_sha256")
        cursors[key] = _cursor(update); record = _record("ACCEPTED_FINE_COMMIT", tip, {"sequence": sequence, "member_key": key, "synthetic": True}); records.append(record); tip = record["record_sha256"]
    high = {key: 1 for key in MEMBERS}
    return _checkpoint(generation=1, parent=origin["checkpoint_sha256"], event=23, target=target, cursors=cursors, ledgers=ledgers, states=origin["states"], tip=tip, high=high), records


def sealed_synchronized_predecessor() -> tuple[dict[str, object], list[dict[str, object]]]:
    """Advance the independent sealed fixture using PROTO15-compatible records."""
    origin = sealed_synthetic_genesis()["checkpoint"]; target = origin["active_event_target_time"]
    cursors, ledgers, records, tip = deepcopy(origin["cursors"]), deepcopy(origin["ledgers"]), [], ROOT
    for sequence, key in enumerate(MEMBERS, 1):
        old = cursors[key]; ledgers[key] = _ledger(target, zero=False)
        update = {**old, "accepted_boundary_time": target, "TDG6_ledger_sha256": digest(ledgers[key]),
                  "cursor_generation": 1, "attempt_serial": 1, "attempt_id": f"fixture-campaign:{key}:23:1",
                  "journal_tip_sha256": tip, "cursor_chain_parent_sha256": old["cursor_chain_sha256"]}
        update.pop("cursor_chain_sha256"); new = _cursor(update)
        payload = {"member_key": key, "predecessor_cursor_chain_sha256": old["cursor_chain_sha256"],
                   "successor_cursor_chain_sha256": new["cursor_chain_sha256"], "accepted_state_sha256": origin["states"][key],
                   "TDG6_ledger_sha256": new["TDG6_ledger_sha256"], "executed_plan": {"synthetic": True},
                   "cursor": new, "sequence": sequence}
        record = _record("ACCEPTED_FINE_COMMIT", tip, payload); records.append(record); tip = record["record_sha256"]; cursors[key] = new
    high = {key: 1 for key in MEMBERS}
    return _checkpoint(generation=1, parent=origin["checkpoint_sha256"], event=23, target=target, cursors=cursors, ledgers=ledgers, states=origin["states"], tip=tip, high=high), records


def construct_successor(previous: Mapping[str, object], prior_records: list[Mapping[str, object]]) -> dict[str, object]:
    prior = _checkpoint(generation=previous["campaign_generation"], parent=previous["parent_checkpoint_sha256"], event=previous["committed_common_event_index"], target=previous["active_event_target_time"], cursors=previous["cursors"], ledgers=previous["ledgers"], states=previous["states"], tip=previous["journal_tip_sha256"], high=previous["attempt_high_water"])
    if prior != previous:
        raise Proto17TheoremError("predecessor checkpoint differs from independent construction")
    sequence = _validate_prior(prior_records, prior["journal_tip_sha256"]); target = prior["active_event_target_time"]
    if any(prior["cursors"][key]["accepted_boundary_time"] != target or prior["ledgers"][key]["current_macro_step_temporal_retry_count"] != 0 for key in MEMBERS):
        raise Proto17TheoremError("predecessor is not synchronously fresh")
    completed = prior["committed_common_event_index"] + 1; next_target = _time_of(Fraction(target["rational"]) + CADENCE)
    payload = {"protocol_artifact_id": ARTIFACT, "campaign_id": prior["campaign_id"], "previous_campaign_generation": prior["campaign_generation"], "previous_checkpoint_sha256": prior["checkpoint_sha256"], "previous_committed_common_event_index": prior["committed_common_event_index"], "completed_common_event_index": completed, "completed_event_time": target, "next_event_target_time": next_target, "member_keys": list(MEMBERS), "prior_cursor_set_sha256": prior["cursor_set_sha256"], "prior_ledger_set_sha256": prior["ledger_set_sha256"], "prior_accepted_state_set_sha256": prior["accepted_state_set_sha256"], "sequence": sequence}
    receipt = _record("COMMON_EVENT_COMMIT", prior["journal_tip_sha256"], payload)
    cursors = {}; high = {}
    for key in MEMBERS:
        old = prior["cursors"][key]; serial = old["attempt_serial"] + 1; update = {**old, "committed_common_event_index": completed, "cursor_generation": old["cursor_generation"] + 1, "event_target_time": next_target, "attempt_serial": serial, "attempt_id": f"{prior['campaign_id']}:{key}:{completed}:{serial}", "journal_tip_sha256": receipt["record_sha256"], "cursor_chain_parent_sha256": old["cursor_chain_sha256"]}; update.pop("cursor_chain_sha256"); cursors[key] = _cursor(update); high[key] = serial
    successor = _checkpoint(generation=prior["campaign_generation"] + 1, parent=prior["checkpoint_sha256"], event=completed, target=next_target, cursors=cursors, ledgers=prior["ledgers"], states=prior["states"], tip=receipt["record_sha256"], high=high)
    return {"receipt": receipt, "checkpoint": successor}


def recover_only_lone_successor(previous: Mapping[str, object], prior_records: list[Mapping[str, object]], suffix: list[Mapping[str, object]]) -> dict[str, object]:
    if not isinstance(suffix, list) or len(suffix) != 1:
        raise Proto17TheoremError("recovery suffix is not exactly one record")
    expected = construct_successor(previous, prior_records)
    if dict(suffix[0]) != expected["receipt"]:
        raise Proto17TheoremError("recovery suffix differs from exact receipt")
    return expected


def qualification(protocol_bytes: bytes, freeze_bytes: bytes) -> dict[str, object]:
    """Run independent, nonphysical structural controls for a future binder."""
    contract = parse_frozen_contract(protocol_bytes, freeze_bytes)
    sealed = sealed_synthetic_genesis(); genesis = sealed["checkpoint"]
    predecessor, history = sealed_synchronized_predecessor(); edge = construct_successor(predecessor, history)
    mutations = {
        "nonzero_generation_root_journal": lambda r, c: c.__setitem__("journal_tip_sha256", ROOT),
        "receipt_cycle_injection": lambda r, c: r["payload"].__setitem__("successor_checkpoint_sha256", c["checkpoint_sha256"]),
        "receipt_sequence": lambda r, c: r["payload"].__setitem__("sequence", 99),
        "successor_untouched_method": lambda r, c: c["cursors"]["RK4-2049"].__setitem__("method", "SSPRK3"),
        "successor_high_water": lambda r, c: c["attempt_high_water"].__setitem__("RK4-2049", 98),
        "successor_checkpoint_hash": lambda r, c: c.__setitem__("checkpoint_sha256", ROOT),
    }
    results: dict[str, bool] = {}
    for name, mutate in mutations.items():
        receipt, successor = deepcopy(edge["receipt"]), deepcopy(edge["checkpoint"]); mutate(receipt, successor)
        try:
            expected = construct_successor(predecessor, history)
            if receipt != expected["receipt"] or successor != expected["checkpoint"]:
                raise Proto17TheoremError("altered successor differs")
        except Proto17TheoremError:
            results[name] = True
        else:
            results[name] = False
    fixture_mutations = {
        "physical_state_digest": lambda x: x["genesis_spec"]["member_descriptors"][0]["state_object"].__setitem__("physical_state_sha256", ROOT),
        "state_object_digest": lambda x: x["genesis_spec"]["member_descriptors"][0].__setitem__("state_object_canonical_sha256", ROOT),
        "manifest_digest": lambda x: x["genesis_spec"].__setitem__("physical_input_manifest_sha256", ROOT),
        "array_storage_key": lambda x: x["genesis_spec"]["member_descriptors"][0]["state_object"]["physical_arrays"][0].__setitem__("npz_storage_key", "wrong"),
        "input_hash": lambda x: x["genesis_spec"]["member_descriptors"][0]["state_object"].__setitem__("input_hash", ROOT),
        "zero_ledger": lambda x: x["genesis_spec"]["member_descriptors"][0]["initial_TDG6_ledger"].__setitem__("accepted_macro_step_count", 1),
        "frozen_restart": lambda x: x["genesis_spec"].__setitem__("restart_coordinate_time", _time_of(Fraction(3, 2))),
        "foreign_protocol": lambda x: x["genesis_spec"].__setitem__("protocol_artifact_id", "foreign"),
    }
    for name, mutate in fixture_mutations.items():
        altered = deepcopy(sealed); mutate(altered)
        try:
            validate_sealed_synthetic_fixture(altered)
        except Proto17TheoremError:
            results[name] = True
        else:
            results[name] = False
    recovered = recover_only_lone_successor(predecessor, history, [edge["receipt"]])
    hashes = {"genesis_checkpoint_sha256": genesis["checkpoint_sha256"], "common_event_record_sha256": edge["receipt"]["record_sha256"], "successor_checkpoint_sha256": edge["checkpoint"]["checkpoint_sha256"]}
    exact_names = {"physical_state_digest", "state_object_digest", "manifest_digest", "array_storage_key", "input_hash", "zero_ledger", "frozen_restart", "foreign_protocol", "receipt_cycle_injection", "receipt_sequence", "successor_untouched_method", "successor_high_water", "successor_checkpoint_hash", "nonzero_generation_root_journal"}
    return {**contract, **hashes, "sealed_hashes_reproduced": hashes == SEALED_SYNTHETIC_HASHES,
            "store_plan_is_complete": set(sealed["store_plan"]) == {"state_objects", "state_paths", "journal_records", "checkpoint_path", "checkpoint"},
            "lone_suffix_reconstructs_exactly": recovered == edge, "mutations_rejected": results,
            "qualified_mutation_names": sorted(exact_names), "qualified_mutation_count": len(exact_names),
            "all_mutations_rejected": set(results) == exact_names and all(results.values()), "structural_only": True, "runtime_or_physical_claim": False}


__all__ = ["Proto17TheoremError", "canonical", "construct_successor", "generation_zero_fixture", "parse_frozen_contract", "qualification", "recover_only_lone_successor", "sealed_synthetic_genesis", "sealed_synchronized_predecessor", "synchronized_predecessor", "validate_sealed_synthetic_fixture"]
