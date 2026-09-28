"""Side-effect-free reference construction for the prospective PROTO17 edge.

This is a semantic oracle, not a campaign runtime.  It reads only supplied
mappings, produces only mappings, and never touches a namespace, filesystem,
or numerical trajectory.  Its purpose is to make the generation-zero and
six-member successful-event constructions exact before a runtime is allowed
to implement them.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
import json
import math
from typing import Any, Mapping

from .protocol_v17 import (
    CADENCE_RATIONAL,
    ARRAY_DESCRIPTOR_FIELDS,
    CHECKPOINT_FIELDS,
    CURSOR_FIELDS,
    DERIVED_GENESIS_FIELDS,
    GENESIS_INPUT_FIELDS,
    LEDGER_FIELDS,
    MEMBER_DESCRIPTOR_FIELDS,
    MEMBER_KEYS,
    RECEIPT_FIELDS,
    ROOT_SHA256,
    SF1_PROTOCOL_V17_ARTIFACT_ID,
    STATE_ARRAY_NAMES,
    STATE_OBJECT_FIELDS,
)


_KNOWN_JOURNAL_KINDS = frozenset({
    "TDG6_REJECTION", "CURSOR_TRANSITION", "ACCEPTED_FINE_COMMIT",
    "TERMINAL_INVALID", "COMMON_EVENT_ROLLBACK", "COMMON_EVENT_COMMIT",
})


class Proto17ConstructionError(ValueError):
    """The supplied pure construction input is not an exact legal predecessor."""


class Proto17RecoveryInvalidStop(Proto17ConstructionError):
    """A journal suffix is not the one exact recoverable successful edge."""


def canonical(value: object) -> bytes:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=True,
        allow_nan=False,
    ).encode("utf-8")


def digest(value: object) -> str:
    return sha256(canonical(value)).hexdigest()


def decode_canonical(encoded: bytes, *, context: str = "canonical JSON") -> object:
    """Decode canonical JSON while rejecting duplicate keys and alternate bytes."""
    def duplicates(pairs: list[tuple[str, object]]) -> dict[str, object]:
        answer: dict[str, object] = {}
        for key, item in pairs:
            if key in answer:
                raise ValueError(f"duplicate key {key}")
            answer[key] = item
        return answer

    try:
        value = json.loads(encoded.decode("utf-8"), object_pairs_hook=duplicates)
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as cause:
        raise Proto17ConstructionError(f"{context} is malformed") from cause
    if canonical(value) != encoded:
        raise Proto17ConstructionError(f"{context} is noncanonical")
    return value


def _sha(name: str, value: object) -> str:
    if not isinstance(value, str) or len(value) != 64 or value.lower() != value:
        raise Proto17ConstructionError(f"{name} must be lowercase SHA-256 text")
    try:
        int(value, 16)
    except ValueError as cause:
        raise Proto17ConstructionError(f"{name} must be hexadecimal") from cause
    return value


def _nonnegative(name: str, value: object) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise Proto17ConstructionError(f"{name} must be a nonnegative integer")
    return value


def _finite_hex(name: str, value: object) -> float:
    if not isinstance(value, str):
        raise Proto17ConstructionError(f"{name} must be binary64 hexadecimal text")
    try:
        answer = float.fromhex(value)
    except ValueError as cause:
        raise Proto17ConstructionError(f"{name} is malformed") from cause
    if not math.isfinite(answer) or answer.hex() != value:
        raise Proto17ConstructionError(f"{name} is not canonical finite binary64 text")
    return answer


def time_identity(value: object) -> dict[str, str]:
    """Validate/construct the PROTO15 rational-plus-binary64 identity."""
    if isinstance(value, Mapping):
        if set(value) != {"rational", "binary64_hex"}:
            raise Proto17ConstructionError("time identity fields differ")
        rational, binary = value["rational"], value["binary64_hex"]
        if not isinstance(rational, str) or not isinstance(binary, str):
            raise Proto17ConstructionError("time identity values must be text")
        try:
            fraction = Fraction(rational)
        except (ValueError, ZeroDivisionError) as cause:
            raise Proto17ConstructionError("time rational is malformed") from cause
        observed = _finite_hex("time binary64 identity", binary)
        if str(fraction) != rational or float(fraction).hex() != observed.hex():
            raise Proto17ConstructionError("time rational and binary64 identity differ")
        return {"rational": str(fraction), "binary64_hex": binary}
    if isinstance(value, bool) or not isinstance(value, (int, float, Fraction)):
        raise Proto17ConstructionError("time identity must be numeric or a mapping")
    numeric = float(value)
    if not math.isfinite(numeric):
        raise Proto17ConstructionError("time identity must be finite")
    fraction = value if isinstance(value, Fraction) else Fraction.from_float(numeric)
    if float(fraction).hex() != numeric.hex():
        raise Proto17ConstructionError("time does not round trip through binary64")
    return {"rational": str(fraction), "binary64_hex": numeric.hex()}


def _next_time(current: Mapping[str, object]) -> dict[str, str]:
    current_identity = time_identity(current)
    return time_identity(Fraction(current_identity["rational"]) + Fraction(CADENCE_RATIONAL))


def _method_grid(key: str) -> tuple[str, int]:
    if key not in MEMBER_KEYS:
        raise Proto17ConstructionError("member key differs from canonical six")
    method, points = key.split("-", 1)
    return method, int(points)


def _safe_relative_path(value: object) -> str:
    if not isinstance(value, str) or not value or value.startswith("/"):
        raise Proto17ConstructionError("path must be a safe relative path")
    parts = value.split("/")
    if any(part in {"", ".", ".."} for part in parts):
        raise Proto17ConstructionError("path must not contain empty, dot, or parent segments")
    return value


def validate_ledger(value: Mapping[str, object], *, accepted_time: Mapping[str, object] | None = None) -> dict[str, object]:
    if not isinstance(value, Mapping) or tuple(value.keys()) != LEDGER_FIELDS:
        # Mapping equality deliberately ignores source order, but the pure
        # oracle returns canonical order.  Reject extra/missing keys here.
        if not isinstance(value, Mapping) or set(value) != set(LEDGER_FIELDS):
            raise Proto17ConstructionError("complete TDG6 ledger fields differ")
    ledger = dict(value)
    last = _finite_hex("ledger last accepted time", ledger["last_accepted_time_hex"])
    for name in (
        "accepted_macro_step_count", "cumulative_temporal_retry_count",
        "current_macro_step_temporal_retry_count",
        "last_accepted_macro_step_temporal_retry_count",
    ):
        _nonnegative(f"ledger {name}", ledger[name])
    debit = ledger["accumulated_debit_vector_hex"]
    if not isinstance(debit, list) or len(debit) != 18:
        raise Proto17ConstructionError("ledger debit vector must contain all eighteen channels")
    for index, item in enumerate(debit):
        if _finite_hex(f"ledger debit {index}", item) < 0.0:
            raise Proto17ConstructionError("ledger debit must be nonnegative")
    rejects = ledger["serialized_temporal_rejections"]
    if not isinstance(rejects, list) or any(not isinstance(item, str) for item in rejects):
        raise Proto17ConstructionError("ledger rejection record differs")
    if (ledger["current_macro_step_temporal_retry_count"] > ledger["cumulative_temporal_retry_count"]
            or ledger["last_accepted_macro_step_temporal_retry_count"] > ledger["cumulative_temporal_retry_count"]
            or len(rejects) != ledger["cumulative_temporal_retry_count"]):
        raise Proto17ConstructionError("ledger retry counts and rejection history differ")
    for item in rejects:
        decoded = decode_canonical(item.encode("utf-8"), context="ledger rejection record")
        if not isinstance(decoded, Mapping) or decoded.get("event_type") != "rejected_TDG6_temporal_admission":
            raise Proto17ConstructionError("ledger rejection event identity differs")
    if accepted_time is not None:
        accepted = time_identity(accepted_time)
        if last.hex() != accepted["binary64_hex"]:
            raise Proto17ConstructionError("ledger last accepted time differs from cursor")
    bare = {
        "last_accepted_time_hex": last.hex(),
        "accepted_macro_step_count": ledger["accepted_macro_step_count"],
        "cumulative_temporal_retry_count": ledger["cumulative_temporal_retry_count"],
        "current_macro_step_temporal_retry_count": ledger["current_macro_step_temporal_retry_count"],
        "last_accepted_macro_step_temporal_retry_count": ledger["last_accepted_macro_step_temporal_retry_count"],
        "accumulated_debit_vector_hex": list(debit),
        "serialized_temporal_rejections": list(rejects),
    }
    return bare


def cursor(payload: Mapping[str, object]) -> dict[str, object]:
    """Canonicalize and validate one complete PROTO15-compatible cursor."""
    if not isinstance(payload, Mapping):
        raise Proto17ConstructionError("cursor must be a mapping")
    value = dict(payload)
    observed = value.pop("cursor_chain_sha256", None)
    if set(value) != set(CURSOR_FIELDS) - {"cursor_chain_sha256"}:
        raise Proto17ConstructionError("cursor fields differ")
    key = value["member_key"]
    method, points = _method_grid(str(key))
    if value["method"] != method or value["point_count"] != points:
        raise Proto17ConstructionError("cursor method/grid differs from member key")
    if value["protocol_artifact_id"] != SF1_PROTOCOL_V17_ARTIFACT_ID:
        raise Proto17ConstructionError("cursor protocol identity differs")
    if not isinstance(value["campaign_id"], str) or not value["campaign_id"] or ":" in value["campaign_id"]:
        raise Proto17ConstructionError("cursor campaign identity differs")
    committed = _nonnegative("cursor event", value["committed_common_event_index"])
    generation = _nonnegative("cursor generation", value["cursor_generation"])
    attempt = _nonnegative("cursor attempt", value["attempt_serial"])
    for name in ("previous_step_index", "previous_transaction_serial"):
        _nonnegative(name, value[name])
    accepted = time_identity(value["accepted_boundary_time"])
    target = time_identity(value["event_target_time"])
    for name in (
        "accepted_state_sha256", "TDG6_ledger_sha256", "journal_tip_sha256",
        "cursor_chain_parent_sha256",
    ):
        _sha(name, value[name])
    if value["mode"] != "FRESH_READY" or value["retry_successor_payload_or_none"] is not None:
        raise Proto17ConstructionError("PROTO17 successful-edge cursor must be fresh ready")
    expected_id = f"{value['campaign_id']}:{key}:{committed}:{attempt}"
    if value["attempt_id"] != expected_id:
        raise Proto17ConstructionError("cursor attempt identity differs")
    value["accepted_boundary_time"] = accepted
    value["event_target_time"] = target
    value["cursor_chain_sha256"] = digest(value)
    if observed is not None and _sha("cursor chain", observed) != value["cursor_chain_sha256"]:
        raise Proto17ConstructionError("cursor digest differs")
    return value


def _sets(cursors: Mapping[str, Mapping[str, object]], ledgers: Mapping[str, Mapping[str, object]], states: Mapping[str, str]) -> tuple[str, str, str]:
    if set(cursors) != set(MEMBER_KEYS) or set(ledgers) != set(MEMBER_KEYS) or set(states) != set(MEMBER_KEYS):
        raise Proto17ConstructionError("set construction requires canonical six members")
    ordered = [
        {
            "member_key": key,
            "cursor_chain_sha256": cursors[key]["cursor_chain_sha256"],
            "TDG6_ledger_sha256": cursors[key]["TDG6_ledger_sha256"],
            "accepted_state_sha256": states[key],
            "accepted_boundary_time": cursors[key]["accepted_boundary_time"],
            "cursor_generation": cursors[key]["cursor_generation"],
            "committed_common_event_index": cursors[key]["committed_common_event_index"],
        }
        for key in MEMBER_KEYS
    ]
    return (
        digest(ordered),
        digest({key: digest(ledgers[key]) for key in MEMBER_KEYS}),
        digest({key: states[key] for key in MEMBER_KEYS}),
    )


def checkpoint(
    *, protocol_artifact_id: str, campaign_id: str, campaign_generation: int,
    committed_common_event_index: int, active_event_target_time: Mapping[str, object],
    cursors: Mapping[str, Mapping[str, object]], ledgers: Mapping[str, Mapping[str, object]],
    states: Mapping[str, str], journal_tip_sha256: str, attempt_high_water: Mapping[str, int],
    cursor_generation_high_water: Mapping[str, int], parent_checkpoint_sha256: str,
) -> dict[str, object]:
    """Construct an exact nonterminal checkpoint from fully supplied state."""
    _nonnegative("campaign generation", campaign_generation)
    _nonnegative("common event index", committed_common_event_index)
    target = time_identity(active_event_target_time)
    if protocol_artifact_id != SF1_PROTOCOL_V17_ARTIFACT_ID:
        raise Proto17ConstructionError("checkpoint protocol identity differs")
    _sha("journal tip", journal_tip_sha256)
    _sha("parent checkpoint", parent_checkpoint_sha256)
    if (campaign_generation == 0) != (parent_checkpoint_sha256 == ROOT_SHA256):
        raise Proto17ConstructionError("checkpoint generation/root-parent relation differs")
    if (campaign_generation == 0) != (journal_tip_sha256 == ROOT_SHA256):
        raise Proto17ConstructionError("checkpoint generation/journal-root relation differs")
    if set(attempt_high_water) != set(MEMBER_KEYS) or set(cursor_generation_high_water) != set(MEMBER_KEYS):
        raise Proto17ConstructionError("checkpoint high-water mappings differ")
    normalized_cursors: dict[str, dict[str, object]] = {}
    normalized_ledgers: dict[str, dict[str, object]] = {}
    normalized_states: dict[str, str] = {}
    for key in MEMBER_KEYS:
        cur = cursor(cursors[key])
        state = _sha("checkpoint state", states[key])
        ledger = validate_ledger(ledgers[key], accepted_time=cur["accepted_boundary_time"])
        if (
            cur["protocol_artifact_id"] != protocol_artifact_id
            or cur["campaign_id"] != campaign_id
            or cur["committed_common_event_index"] != committed_common_event_index
            or cur["event_target_time"] != target
            or cur["accepted_state_sha256"] != state
            or cur["TDG6_ledger_sha256"] != digest(ledger)
            or cur["attempt_serial"] != _nonnegative("attempt high-water", attempt_high_water[key])
            or cur["cursor_generation"] != _nonnegative("cursor high-water", cursor_generation_high_water[key])
        ):
            raise Proto17ConstructionError("checkpoint member cross-binding differs")
        normalized_cursors[key] = cur
        normalized_ledgers[key] = ledger
        normalized_states[key] = state
    cursor_set, ledger_set, state_set = _sets(normalized_cursors, normalized_ledgers, normalized_states)
    bare: dict[str, object] = {
        "protocol_artifact_id": protocol_artifact_id,
        "campaign_id": campaign_id,
        "campaign_generation": campaign_generation,
        "committed_common_event_index": committed_common_event_index,
        "active_event_target_time": target,
        "disposition": "nonterminal",
        "terminal_member_key": None,
        "terminal_record_sha256": None,
        "member_keys": list(MEMBER_KEYS),
        "cursors": normalized_cursors,
        "ledgers": normalized_ledgers,
        "states": normalized_states,
        "cursor_set_sha256": cursor_set,
        "ledger_set_sha256": ledger_set,
        "accepted_state_set_sha256": state_set,
        "journal_tip_sha256": journal_tip_sha256,
        "terminal_lock": False,
        "attempt_high_water": {key: attempt_high_water[key] for key in MEMBER_KEYS},
        "cursor_generation_high_water": {key: cursor_generation_high_water[key] for key in MEMBER_KEYS},
        "parent_checkpoint_sha256": parent_checkpoint_sha256,
        "terminal_closure": None,
    }
    bare["checkpoint_sha256"] = digest(bare)
    return bare


def validate_checkpoint(value: Mapping[str, object]) -> dict[str, object]:
    if not isinstance(value, Mapping) or set(value) != set(CHECKPOINT_FIELDS):
        raise Proto17ConstructionError("checkpoint fields differ")
    observed = _sha("checkpoint", value["checkpoint_sha256"])
    bare = dict(value)
    bare.pop("checkpoint_sha256")
    if digest(bare) != observed:
        raise Proto17ConstructionError("checkpoint digest differs")
    if value["member_keys"] != list(MEMBER_KEYS):
        raise Proto17ConstructionError("checkpoint member order differs")
    if value["disposition"] != "nonterminal" or value["terminal_lock"] is not False or any(
        value[name] is not None for name in ("terminal_member_key", "terminal_record_sha256", "terminal_closure")
    ):
        raise Proto17ConstructionError("pure successful-edge checkpoint terminal fields differ")
    return checkpoint(
        protocol_artifact_id=value["protocol_artifact_id"], campaign_id=value["campaign_id"],
        campaign_generation=_nonnegative("campaign generation", value["campaign_generation"]),
        committed_common_event_index=_nonnegative("common event", value["committed_common_event_index"]),
        active_event_target_time=value["active_event_target_time"], cursors=value["cursors"],
        ledgers=value["ledgers"], states=value["states"], journal_tip_sha256=value["journal_tip_sha256"],
        attempt_high_water=value["attempt_high_water"], cursor_generation_high_water=value["cursor_generation_high_water"],
        parent_checkpoint_sha256=value["parent_checkpoint_sha256"],
    )


def _validate_descriptor(item: Mapping[str, object], *, event: int, restart_time: Mapping[str, object]) -> dict[str, object]:
    if not isinstance(item, Mapping) or set(item) != set(MEMBER_DESCRIPTOR_FIELDS):
        raise Proto17ConstructionError("genesis member descriptor fields differ")
    value = dict(item)
    key = value["member_key"]
    method, points = _method_grid(str(key))
    if value["method"] != method or value["point_count"] != points:
        raise Proto17ConstructionError("genesis descriptor method/grid differs")
    state = value["state_object"]
    if not isinstance(state, Mapping) or set(state) != set(STATE_OBJECT_FIELDS):
        raise Proto17ConstructionError("genesis descriptor omits complete visible state object")
    state = dict(state)
    if state["schema_id"] != "FGC-2-SF1-PROTO17-state-object-v1":
        raise Proto17ConstructionError("state-object schema identity differs")
    physical = _sha("physical state", state["physical_state_sha256"])
    accepted = time_identity(value["accepted_boundary_time"])
    if accepted != restart_time:
        raise Proto17ConstructionError("genesis descriptor accepted time differs from restart time")
    if (
        state["member_key"] != key or state["method"] != method or state["point_count"] != points
        or time_identity(state["coordinate_time"]) != accepted
        or time_identity(state["accepted_boundary_time"]) != accepted
        or state["step_index"] != value["step_index"]
        or state["transaction_serial"] != value["transaction_serial"]
    ):
        raise Proto17ConstructionError("state object identity/time/counter binding differs")
    _nonnegative("state object step", state["step_index"])
    _nonnegative("state object transaction", state["transaction_serial"])
    _sha("state object input hash", state["input_hash"])
    arrays = state["physical_arrays"]
    if not isinstance(arrays, list) or len(arrays) != len(STATE_ARRAY_NAMES):
        raise Proto17ConstructionError("state object must bind exactly seven array descriptors")
    normalized_arrays: list[dict[str, object]] = []
    for expected_name, array in zip(STATE_ARRAY_NAMES, arrays, strict=True):
        if not isinstance(array, Mapping) or set(array) != set(ARRAY_DESCRIPTOR_FIELDS):
            raise Proto17ConstructionError("array descriptor fields differ")
        descriptor = dict(array)
        if (descriptor["logical_name"] != expected_name
                or not isinstance(descriptor["source_bundle_key"], str) or not descriptor["source_bundle_key"]
                or not isinstance(descriptor["npz_storage_key"], str) or not descriptor["npz_storage_key"]):
            raise Proto17ConstructionError("array descriptor identity differs")
        expected_shape = {
            "u": [points, 6], "p": [points, 6], "q": [points, 6],
            "tracer_positions": [48], "tracer_proper_times": [48],
            "event_proper_times": [24, 48], "event_fields": [24, 48, 6],
        }[expected_name]
        if (descriptor["dtype"] != "<f8" or descriptor["layout"] != "C"
                or not isinstance(descriptor["shape"], list) or not descriptor["shape"]):
            raise Proto17ConstructionError("array descriptor layout differs")
        if any(isinstance(dimension, bool) or not isinstance(dimension, int) or dimension < 1 for dimension in descriptor["shape"]):
            raise Proto17ConstructionError("array descriptor shape differs")
        if descriptor["shape"] != expected_shape:
            raise Proto17ConstructionError("array descriptor shape is not the frozen member shape")
        _sha("array little-endian C bytes", descriptor["little_endian_c_bytes_sha256"])
        normalized_arrays.append(descriptor)
    if len({array["npz_storage_key"] for array in normalized_arrays}) != len(STATE_ARRAY_NAMES):
        raise Proto17ConstructionError("physical array NPZ storage keys must be unique")
    for name in ("runtime_monitor_state", "causal_state", "tracer_state", "event_history_state"):
        if not isinstance(state[name], Mapping) or not state[name]:
            raise Proto17ConstructionError(f"state object {name} differs")
    _nonnegative("state object source retry count", state["source_retry_count"])
    _nonnegative("state object CFL retry count", state["CFL_retry_count"])
    normalized_state = {
        "schema_id": state["schema_id"],
        "member_key": key, "method": method, "point_count": points,
        "coordinate_time": accepted, "accepted_boundary_time": accepted,
        "step_index": state["step_index"], "transaction_serial": state["transaction_serial"],
        "input_hash": state["input_hash"],
        "source_retry_count": state["source_retry_count"], "CFL_retry_count": state["CFL_retry_count"],
        "physical_state_sha256": physical,
        "physical_arrays": normalized_arrays,
        "runtime_monitor_state": dict(state["runtime_monitor_state"]),
        "causal_state": dict(state["causal_state"]),
        "tracer_state": dict(state["tracer_state"]),
        "event_history_state": dict(state["event_history_state"]),
        "restart_payload_sha256": _sha("state restart payload", state["restart_payload_sha256"]),
    }
    state_hash = digest(normalized_state)
    if _sha("state-object canonical", value["state_object_canonical_sha256"]) != state_hash:
        raise Proto17ConstructionError("state-object canonical digest differs")
    ledger = validate_ledger(value["initial_TDG6_ledger"], accepted_time=accepted)
    if ledger["current_macro_step_temporal_retry_count"] != 0:
        raise Proto17ConstructionError("genesis descriptor has pending TDG6 retry")
    if (ledger["accepted_macro_step_count"] != 0
            or ledger["cumulative_temporal_retry_count"] != 0
            or ledger["last_accepted_macro_step_temporal_retry_count"] != 0
            or ledger["serialized_temporal_rejections"]
            or any(item != "0x0.0p+0" for item in ledger["accumulated_debit_vector_hex"])):
        raise Proto17ConstructionError("genesis TDG6 ledger is not the frozen zero initialization")
    if _sha("descriptor ledger", value["initial_TDG6_ledger_sha256"]) != digest(ledger):
        raise Proto17ConstructionError("descriptor ledger digest differs")
    expected_source = "RSP2" if key == "SSPRK3-16385" else "PROTO12"
    if value["source_checkpoint"] != expected_source:
        raise Proto17ConstructionError("source checkpoint identity differs")
    _nonnegative("descriptor step", value["step_index"])
    _nonnegative("descriptor transaction", value["transaction_serial"])
    if not isinstance(value["source_bundle_key"], str) or not value["source_bundle_key"]:
        raise Proto17ConstructionError("source bundle key differs")
    if any(array["source_bundle_key"] != value["source_bundle_key"] for array in normalized_arrays):
        raise Proto17ConstructionError("physical arrays do not bind descriptor source bundle")
    bare = {
        **value,
        "accepted_boundary_time": accepted,
        "state_object": normalized_state,
        "initial_TDG6_ledger": ledger,
    }
    return bare


def build_genesis(spec: Mapping[str, object], *, _validate_derived: bool = True) -> dict[str, object]:
    """Derive generation zero from visible state objects and full TDG6 ledgers."""
    if not isinstance(spec, Mapping) or set(spec) != set(GENESIS_INPUT_FIELDS):
        raise Proto17ConstructionError("genesis input fields differ")
    supplied = dict(spec)
    expected_derived = {name: supplied[name] for name in DERIVED_GENESIS_FIELDS}
    source = {name: supplied[name] for name in GENESIS_INPUT_FIELDS if name not in DERIVED_GENESIS_FIELDS}
    if (
        source["genesis_algorithm_id"] != "FGC-2-SF1-PROTO17-pure-construction-v1"
        or source["canonical_json_encoding"] != "utf-8;json-sort-keys;separators=(',',':');ensure-ascii=true;allow-nan=false"
        or source["state_object_schema_id"] != "FGC-2-SF1-PROTO17-state-object-v1"
        or source["tdg6_ledger_schema_id"] != "FGC-2-SF1-PROTO17-tdg6-ledger-v1"
        or source["cursor_schema_id"] != "FGC-2-SF1-PROTO17-cursor-v1"
        or source["checkpoint_schema_id"] != "FGC-2-SF1-PROTO17-checkpoint-v1"
    ):
        raise Proto17ConstructionError("genesis construction/schema identity differs")
    if source["journal_root_parent_sha256"] != ROOT_SHA256:
        raise Proto17ConstructionError("genesis journal root differs")
    if source["protocol_artifact_id"] != SF1_PROTOCOL_V17_ARTIFACT_ID:
        raise Proto17ConstructionError("genesis protocol identity differs")
    if not isinstance(source["campaign_id"], str) or not source["campaign_id"] or ":" in source["campaign_id"]:
        raise Proto17ConstructionError("genesis campaign identity differs")
    for name in (
        "protocol_config_sha256", "protocol_freeze_result_sha256", "hlt13_result_sha256",
        "run_plan_sha256", "runtime_module_sha256", "adapter_module_sha256",
        "runner_sha256",
    ):
        _sha(name, source[name])
    manifest = source["physical_input_manifest"]
    if not isinstance(manifest, Mapping) or set(manifest) != {
        "format", "relative_input_root", "bundle_members",
    }:
        raise Proto17ConstructionError("physical input manifest fields differ")
    if (manifest["format"] != "FGC-2-SF1-PROTO17-physical-input-manifest-v1"
            or not isinstance(manifest["relative_input_root"], str)
            or not manifest["relative_input_root"]
            or not isinstance(manifest["bundle_members"], list)
            or len(manifest["bundle_members"]) != len(MEMBER_KEYS)):
        raise Proto17ConstructionError("physical input manifest values differ")
    _safe_relative_path(manifest["relative_input_root"])
    bundles: dict[str, Mapping[str, object]] = {}
    for entry in manifest["bundle_members"]:
        required = {
            "source_bundle_key", "source_checkpoint", "relative_input_path",
            "raw_file_sha256", "required_npz_member_keys",
        }
        if (not isinstance(entry, Mapping) or set(entry) != required
                or not isinstance(entry["source_bundle_key"], str) or not entry["source_bundle_key"]
                or entry["source_checkpoint"] not in {"PROTO12", "RSP2"}
                or not isinstance(entry["relative_input_path"], str)
                or not isinstance(entry["required_npz_member_keys"], list)
                or len(entry["required_npz_member_keys"]) != len(STATE_ARRAY_NAMES) + 1):
            raise Proto17ConstructionError("physical input bundle member differs")
        _safe_relative_path(entry["relative_input_path"])
        _sha("physical input raw file", entry["raw_file_sha256"])
        if entry["source_bundle_key"] in bundles:
            raise Proto17ConstructionError("physical input bundle key is duplicated")
        bundles[entry["source_bundle_key"]] = entry
    # The manifest digest is over the bare manifest.  It is a sibling field,
    # never an item inside the manifest it authenticates.
    if _sha("physical input manifest", source["physical_input_manifest_sha256"]) != digest(dict(manifest)):
        raise Proto17ConstructionError("physical input manifest digest differs")
    if not isinstance(source["numerical_environment"], Mapping) or not source["numerical_environment"]:
        raise Proto17ConstructionError("genesis numerical environment differs")
    if not isinstance(source["namespace"], str) or not source["namespace"]:
        raise Proto17ConstructionError("genesis namespace differs")
    event = _nonnegative("genesis event", source["common_event_index"])
    restart_time = time_identity(source["restart_coordinate_time"])
    target = time_identity(source["active_event_target_time"])
    if source["member_keys_in_canonical_order"] != list(MEMBER_KEYS):
        raise Proto17ConstructionError("genesis canonical member key list differs")
    if (event != 23 or restart_time != time_identity(Fraction(23, 16))
            or target != time_identity(Fraction(3, 2))):
        raise Proto17ConstructionError("genesis event/restart/target identity differs from frozen PROTO17")
    descriptors = source["member_descriptors"]
    if not isinstance(descriptors, list) or len(descriptors) != len(MEMBER_KEYS):
        raise Proto17ConstructionError("genesis must visibly contain six descriptors")
    parsed = [_validate_descriptor(item, event=event, restart_time=restart_time) for item in descriptors]
    if [item["member_key"] for item in parsed] != list(MEMBER_KEYS):
        raise Proto17ConstructionError("genesis descriptor order differs")
    for item in parsed:
        bundle = bundles.get(item["source_bundle_key"])
        if bundle is None or bundle["source_checkpoint"] != item["source_checkpoint"]:
            raise Proto17ConstructionError("descriptor source bundle does not bind physical manifest")
        if bundle["required_npz_member_keys"] != [
            array["npz_storage_key"] for array in item["state_object"]["physical_arrays"]
        ] + ["metadata_utf8"]:
            raise Proto17ConstructionError("physical manifest npz keys differ from state-object descriptors")
    if set(bundles) != {item["source_bundle_key"] for item in parsed} or len({item["source_bundle_key"] for item in parsed}) != len(MEMBER_KEYS):
        raise Proto17ConstructionError("physical manifest and descriptors do not form one six-member bijection")
    states = {key: parsed[index]["state_object_canonical_sha256"] for index, key in enumerate(MEMBER_KEYS)}
    ledgers = {key: parsed[index]["initial_TDG6_ledger"] for index, key in enumerate(MEMBER_KEYS)}
    cursors: dict[str, dict[str, object]] = {}
    for index, key in enumerate(MEMBER_KEYS):
        descriptor = parsed[index]
        method, points = _method_grid(key)
        cursors[key] = cursor({
            "protocol_artifact_id": source["protocol_artifact_id"], "campaign_id": source["campaign_id"],
            "member_key": key, "method": method, "point_count": points,
            "committed_common_event_index": event, "accepted_boundary_time": restart_time,
            "accepted_state_sha256": states[key], "previous_step_index": descriptor["step_index"],
            "previous_transaction_serial": descriptor["transaction_serial"],
            "TDG6_ledger_sha256": digest(ledgers[key]), "cursor_generation": 0,
            "event_target_time": target, "attempt_serial": 0,
            "attempt_id": f"{source['campaign_id']}:{key}:{event}:0", "mode": "FRESH_READY",
            "retry_successor_payload_or_none": None, "journal_tip_sha256": ROOT_SHA256,
            "cursor_chain_parent_sha256": ROOT_SHA256,
        })
    high = {key: 0 for key in MEMBER_KEYS}
    genesis_checkpoint = checkpoint(
        protocol_artifact_id=source["protocol_artifact_id"], campaign_id=source["campaign_id"],
        campaign_generation=0, committed_common_event_index=event, active_event_target_time=target,
        cursors=cursors, ledgers=ledgers, states=states, journal_tip_sha256=ROOT_SHA256,
        attempt_high_water=high, cursor_generation_high_water=high,
        parent_checkpoint_sha256=ROOT_SHA256,
    )
    descriptor_hash = digest(parsed)
    materialized = {
        **source,
        "member_descriptors": parsed,
        "member_descriptor_set_sha256": descriptor_hash,
        "genesis_checkpoint_sha256": genesis_checkpoint["checkpoint_sha256"],
        "cursor_set_sha256": genesis_checkpoint["cursor_set_sha256"],
        "ledger_set_sha256": genesis_checkpoint["ledger_set_sha256"],
        "accepted_state_set_sha256": genesis_checkpoint["accepted_state_set_sha256"],
    }
    if set(materialized) != set(GENESIS_INPUT_FIELDS):
        raise AssertionError("internal genesis materialization drift")
    if _validate_derived and expected_derived != {
        name: materialized[name] for name in DERIVED_GENESIS_FIELDS
    }:
        raise Proto17ConstructionError("genesis derived checksum fields differ")
    if len(set(states.values())) != len(MEMBER_KEYS):
        raise Proto17ConstructionError("generation-zero state-object addresses must be unique")
    store_plan = {
        "state_objects": {states[key]: parsed[index]["state_object"] for index, key in enumerate(MEMBER_KEYS)},
        "state_paths": {key: f"states/{states[key]}.json" for key in MEMBER_KEYS},
        "journal_records": [],
        "checkpoint_path": (
            f"checkpoints/{genesis_checkpoint['campaign_generation']:020d}-"
            f"{genesis_checkpoint['checkpoint_sha256']}.json"
        ),
        "checkpoint": genesis_checkpoint,
    }
    return {"genesis_spec": materialized, "checkpoint": genesis_checkpoint, "store_plan": store_plan}


def _next_sequence(previous: Mapping[str, object], prior_journal_records: list[Mapping[str, object]]) -> int:
    """Validate hash-contiguous structural journal evidence and derive sequence.

    Historical record-payload replay is deliberately deferred to the later
    runtime/recovery authority.  This pure oracle consumes only a predecessor
    checkpoint/history pair that that authority has already authenticated.
    """
    prior = validate_checkpoint(previous)
    if not isinstance(prior_journal_records, list):
        raise Proto17ConstructionError("prior journal evidence must be a list")
    parent = ROOT_SHA256
    for expected, record in enumerate(prior_journal_records, 1):
        if not isinstance(record, Mapping) or set(record) != {
            "record_kind", "journal_parent_sha256", "payload", "record_sha256",
        }:
            raise Proto17ConstructionError("prior journal envelope fields differ")
        if record["record_kind"] not in _KNOWN_JOURNAL_KINDS:
            raise Proto17ConstructionError("prior journal kind is not in the frozen lifecycle vocabulary")
        if not isinstance(record["payload"], Mapping) or record["payload"].get("sequence") != expected:
            raise Proto17ConstructionError("prior journal sequence differs")
        bare = dict(record)
        observed = _sha("prior journal record", bare.pop("record_sha256"))
        if digest(bare) != observed or record["journal_parent_sha256"] != parent:
            raise Proto17ConstructionError("prior journal hash chain differs")
        parent = observed
    if parent != prior["journal_tip_sha256"]:
        raise Proto17ConstructionError("prior journal evidence does not end at checkpoint tip")
    return len(prior_journal_records) + 1


def _receipt_from_predecessor(previous: Mapping[str, object], *, prior_journal_records: list[Mapping[str, object]]) -> dict[str, object]:
    prior = validate_checkpoint(previous)
    sequence = _next_sequence(prior, prior_journal_records)
    target = prior["active_event_target_time"]
    for key in MEMBER_KEYS:
        cur = prior["cursors"][key]
        ledger = prior["ledgers"][key]
        if (
            cur["mode"] != "FRESH_READY"
            or cur["retry_successor_payload_or_none"] is not None
            or cur["accepted_boundary_time"] != target
            or ledger["current_macro_step_temporal_retry_count"] != 0
        ):
            raise Proto17ConstructionError("common event predecessor is not synchronously fresh")
    event = prior["committed_common_event_index"]
    payload: dict[str, object] = {
        "protocol_artifact_id": prior["protocol_artifact_id"], "campaign_id": prior["campaign_id"],
        "previous_campaign_generation": prior["campaign_generation"],
        "previous_checkpoint_sha256": prior["checkpoint_sha256"],
        "previous_committed_common_event_index": event,
        "completed_common_event_index": event + 1,
        "completed_event_time": target, "next_event_target_time": _next_time(target),
        "member_keys": list(MEMBER_KEYS), "prior_cursor_set_sha256": prior["cursor_set_sha256"],
        "prior_ledger_set_sha256": prior["ledger_set_sha256"],
        "prior_accepted_state_set_sha256": prior["accepted_state_set_sha256"],
        "sequence": sequence,
    }
    if set(payload) != set(RECEIPT_FIELDS) | {"sequence"}:
        raise AssertionError("internal receipt fields drift")
    bare = {
        "record_kind": "COMMON_EVENT_COMMIT",
        "journal_parent_sha256": prior["journal_tip_sha256"],
        "payload": payload,
    }
    return {**bare, "record_sha256": digest(bare)}


def construct_common_event(
    previous: Mapping[str, object], *, prior_journal_records: list[Mapping[str, object]],
) -> dict[str, object]:
    """Construct the sole legal successful event receipt and successor checkpoint."""
    prior = validate_checkpoint(previous)
    receipt = _receipt_from_predecessor(prior, prior_journal_records=prior_journal_records)
    payload = receipt["payload"]
    record_sha = receipt["record_sha256"]
    event = payload["completed_common_event_index"]
    target = payload["next_event_target_time"]
    cursors: dict[str, dict[str, object]] = {}
    attempts: dict[str, int] = {}
    generations: dict[str, int] = {}
    for key in MEMBER_KEYS:
        old = prior["cursors"][key]
        next_attempt = prior["attempt_high_water"][key] + 1
        next_generation = prior["cursor_generation_high_water"][key] + 1
        updated = {**old,
            "committed_common_event_index": event,
            "cursor_generation": next_generation,
            "event_target_time": target,
            "attempt_serial": next_attempt,
            "attempt_id": f"{prior['campaign_id']}:{key}:{event}:{next_attempt}",
            "journal_tip_sha256": record_sha,
            "cursor_chain_parent_sha256": old["cursor_chain_sha256"],
        }
        updated.pop("cursor_chain_sha256", None)
        cursors[key] = cursor(updated)
        attempts[key] = next_attempt
        generations[key] = next_generation
    successor = checkpoint(
        protocol_artifact_id=prior["protocol_artifact_id"], campaign_id=prior["campaign_id"],
        campaign_generation=prior["campaign_generation"] + 1,
        committed_common_event_index=event, active_event_target_time=target,
        cursors=cursors, ledgers=prior["ledgers"], states=prior["states"],
        journal_tip_sha256=record_sha, attempt_high_water=attempts,
        cursor_generation_high_water=generations,
        parent_checkpoint_sha256=prior["checkpoint_sha256"],
    )
    return {"receipt": receipt, "checkpoint": successor}


def recover_lone_common_event(
    previous: Mapping[str, object], suffix: list[Mapping[str, object]], *,
    prior_journal_records: list[Mapping[str, object]],
) -> dict[str, object]:
    """Recover only the exact one-record successful suffix; reject every other tail."""
    if not isinstance(suffix, list) or len(suffix) != 1:
        raise Proto17RecoveryInvalidStop("successful suffix must contain exactly one record")
    expected = construct_common_event(previous, prior_journal_records=prior_journal_records)
    observed = suffix[0]
    if not isinstance(observed, Mapping) or dict(observed) != expected["receipt"]:
        raise Proto17RecoveryInvalidStop("successful suffix differs from exact predecessor construction")
    return expected


def validate_common_event_successor(
    previous: Mapping[str, object], *, prior_journal_records: list[Mapping[str, object]],
    receipt: Mapping[str, object], successor_checkpoint: Mapping[str, object],
) -> dict[str, object]:
    """Verify a supplied successful edge against the one pure construction.

    Exact mapping equality intentionally enforces the copy whitelist, every
    high-water value, receipt acyclicity, and every successor-checkpoint field.
    """
    expected = construct_common_event(previous, prior_journal_records=prior_journal_records)
    if not isinstance(receipt, Mapping) or dict(receipt) != expected["receipt"]:
        raise Proto17ConstructionError("common-event receipt differs from exact construction")
    if not isinstance(successor_checkpoint, Mapping) or dict(successor_checkpoint) != expected["checkpoint"]:
        raise Proto17ConstructionError("common-event successor checkpoint differs from exact construction")
    return expected


def synthetic_genesis_fixture(*, synchronized: bool = False) -> dict[str, object]:
    """Return a deterministic, nonphysical fixture for integration qualification.

    The fixture is intentionally synthetic: its hashes bind descriptors only,
    not numerical arrays or a run namespace.
    """
    target = time_identity(Fraction(3, 2))
    restart = target if synchronized else time_identity(Fraction(23, 16))
    descriptors: list[dict[str, object]] = []
    bundles: list[dict[str, object]] = []
    for key in MEMBER_KEYS:
        method, points = _method_grid(key)
        arrays = [
            {
                "logical_name": name, "source_bundle_key": f"fixture-{key}",
                "npz_storage_key": f"arr_{name}", "dtype": "<f8", "layout": "C",
                "shape": {
                    "u": [points, 6], "p": [points, 6], "q": [points, 6],
                    "tracer_positions": [48], "tracer_proper_times": [48],
                    "event_proper_times": [24, 48], "event_fields": [24, 48, 6],
                }[name],
                "little_endian_c_bytes_sha256": digest({"fixture": key, "array": name}),
            }
            for name in STATE_ARRAY_NAMES
        ]
        physical = digest({"fixture": key, "kind": "u-p-q"})
        state = {
            "schema_id": "FGC-2-SF1-PROTO17-state-object-v1",
            "member_key": key, "method": method, "point_count": points,
            "coordinate_time": restart, "accepted_boundary_time": restart,
            "step_index": 1, "transaction_serial": 2, "input_hash": digest({"fixture": key, "input": True}),
            "source_retry_count": 0, "CFL_retry_count": 0,
            "physical_state_sha256": physical,
            "physical_arrays": arrays,
            "runtime_monitor_state": {"accepted_stage_count": 1},
            "causal_state": {"remaining_buffer": "1"},
            "tracer_state": {"tracer_count": 48},
            "event_history_state": {"sample_count": 24},
            "restart_payload_sha256": digest({"fixture": key, "restart": True}),
        }
        ledger = {
            "last_accepted_time_hex": restart["binary64_hex"],
            "accepted_macro_step_count": 0,
            "cumulative_temporal_retry_count": 0,
            "current_macro_step_temporal_retry_count": 0,
            "last_accepted_macro_step_temporal_retry_count": 0,
            "accumulated_debit_vector_hex": ["0x0.0p+0"] * 18,
            "serialized_temporal_rejections": [],
        }
        descriptors.append({
            "member_key": key, "method": method, "point_count": points,
            "source_checkpoint": "RSP2" if key == "SSPRK3-16385" else "PROTO12",
            "state_object": state, "state_object_canonical_sha256": digest(state),
            "accepted_boundary_time": restart, "step_index": 1, "transaction_serial": 2,
            "initial_TDG6_ledger": ledger, "initial_TDG6_ledger_sha256": digest(ledger),
            "source_bundle_key": f"fixture-{key}",
        })
        bundles.append({
            "source_bundle_key": f"fixture-{key}",
            "source_checkpoint": "RSP2" if key == "SSPRK3-16385" else "PROTO12",
            "relative_input_path": f"fixture/{key}.npz",
            "raw_file_sha256": digest({"fixture": key, "raw": True}),
            "required_npz_member_keys": [item["npz_storage_key"] for item in arrays] + ["metadata_utf8"],
        })
    manifest = {
        "format": "FGC-2-SF1-PROTO17-physical-input-manifest-v1",
        "relative_input_root": "synthetic/no-read", "bundle_members": bundles,
    }
    bare = {
        "genesis_algorithm_id": "FGC-2-SF1-PROTO17-pure-construction-v1",
        "canonical_json_encoding": "utf-8;json-sort-keys;separators=(',',':');ensure-ascii=true;allow-nan=false",
        "state_object_schema_id": "FGC-2-SF1-PROTO17-state-object-v1",
        "tdg6_ledger_schema_id": "FGC-2-SF1-PROTO17-tdg6-ledger-v1",
        "cursor_schema_id": "FGC-2-SF1-PROTO17-cursor-v1",
        "checkpoint_schema_id": "FGC-2-SF1-PROTO17-checkpoint-v1",
        "protocol_artifact_id": "FGC-2-SF1-PROTO17", "protocol_config_sha256": digest("protocol"),
        "protocol_freeze_result_sha256": digest("freeze"), "hlt13_result_sha256": digest("hlt13"),
        "run_plan_sha256": digest("plan"), "runtime_module_sha256": digest("runtime"),
        "adapter_module_sha256": digest("adapter"), "runner_sha256": digest("runner"),
        "numerical_environment": {"python": "fixture"}, "campaign_id": "fixture-campaign",
        "namespace": "synthetic/no-write", "journal_root_parent_sha256": ROOT_SHA256,
        "common_event_index": 23, "restart_coordinate_time": restart,
        "active_event_target_time": target, "member_keys_in_canonical_order": list(MEMBER_KEYS),
        "member_descriptors": descriptors, "physical_input_manifest": manifest,
        "physical_input_manifest_sha256": digest(manifest),
        "member_descriptor_set_sha256": "0" * 64, "genesis_checkpoint_sha256": "0" * 64,
        "cursor_set_sha256": "0" * 64, "ledger_set_sha256": "0" * 64,
        "accepted_state_set_sha256": "0" * 64,
    }
    return build_genesis(bare, _validate_derived=False)["genesis_spec"]


def derive_generation_zero(spec_or_inputs: Mapping[str, object]) -> dict[str, object]:
    """Stable integration alias for :func:`build_genesis`."""
    return build_genesis(spec_or_inputs)


def derive_common_event_successor(
    checkpoint_value: Mapping[str, object], prior_journal_records: list[Mapping[str, object]],
) -> dict[str, object]:
    """Stable integration alias that derives sequence from complete journal evidence."""
    return construct_common_event(checkpoint_value, prior_journal_records=prior_journal_records)


def synthetic_common_event_fixture() -> dict[str, object]:
    """Create a synthetic, journal-backed synchronized predecessor fixture.

    Generation zero remains the actual frozen 23/16 -> 3/2 handoff.  Six
    synthetic accepted-fine records then advance its six cursors to the target
    before the common-event edge is exercised; no impossible target-at-genesis
    shortcut is used.
    """
    origin = build_genesis(synthetic_genesis_fixture())["checkpoint"]
    target = origin["active_event_target_time"]
    cursors = {key: dict(origin["cursors"][key]) for key in MEMBER_KEYS}
    ledgers = {key: dict(origin["ledgers"][key]) for key in MEMBER_KEYS}
    records: list[dict[str, object]] = []
    tip = ROOT_SHA256
    for sequence, key in enumerate(MEMBER_KEYS, 1):
        old = cursors[key]
        ledger = dict(ledgers[key])
        ledger["last_accepted_time_hex"] = target["binary64_hex"]
        ledger["accepted_macro_step_count"] = 1
        ledgers[key] = validate_ledger(ledger, accepted_time=target)
        updated = {**old,
            "accepted_boundary_time": target, "TDG6_ledger_sha256": digest(ledgers[key]),
            "cursor_generation": 1, "attempt_serial": 1,
            "attempt_id": f"{origin['campaign_id']}:{key}:{origin['committed_common_event_index']}:1",
            "journal_tip_sha256": tip, "cursor_chain_parent_sha256": old["cursor_chain_sha256"],
        }
        updated.pop("cursor_chain_sha256", None)
        new = cursor(updated)
        payload = {
            "member_key": key, "predecessor_cursor_chain_sha256": old["cursor_chain_sha256"],
            "successor_cursor_chain_sha256": new["cursor_chain_sha256"],
            "accepted_state_sha256": origin["states"][key],
            "TDG6_ledger_sha256": new["TDG6_ledger_sha256"],
            "executed_plan": {"synthetic": True}, "cursor": new, "sequence": sequence,
        }
        envelope = {"record_kind": "ACCEPTED_FINE_COMMIT", "journal_parent_sha256": tip, "payload": payload}
        record = {**envelope, "record_sha256": digest(envelope)}
        records.append(record)
        tip = record["record_sha256"]
        cursors[key] = new
    high = {key: 1 for key in MEMBER_KEYS}
    predecessor = checkpoint(
        protocol_artifact_id=origin["protocol_artifact_id"], campaign_id=origin["campaign_id"],
        campaign_generation=1, committed_common_event_index=origin["committed_common_event_index"],
        active_event_target_time=target, cursors=cursors, ledgers=ledgers, states=origin["states"],
        journal_tip_sha256=tip, attempt_high_water=high, cursor_generation_high_water=high,
        parent_checkpoint_sha256=origin["checkpoint_sha256"],
    )
    return {"predecessor": predecessor, "prior_journal_records": records}


def synthetic_common_event_predecessor() -> dict[str, object]:
    """Compatibility helper returning the journal-backed predecessor checkpoint."""
    return synthetic_common_event_fixture()["predecessor"]


def run_pure_construction_qualification() -> dict[str, object]:
    """Deterministic synthetic proof summary used by a future freeze reproducer."""
    fixture = synthetic_genesis_fixture()
    genesis = build_genesis(fixture)
    common_fixture = synthetic_common_event_fixture()
    predecessor = common_fixture["predecessor"]
    prior_journal = common_fixture["prior_journal_records"]
    transition = construct_common_event(predecessor, prior_journal_records=prior_journal)
    recovered = recover_lone_common_event(
        predecessor, [transition["receipt"]], prior_journal_records=prior_journal,
    )
    mutations: dict[str, bool] = {}
    for name, altered in {
        "physical_state_digest": (lambda x: x["member_descriptors"][0]["state_object"].__setitem__("physical_state_sha256", "0" * 64)),
        "state_object_digest": (lambda x: x["member_descriptors"][0].__setitem__("state_object_canonical_sha256", "0" * 64)),
        "manifest_digest": (lambda x: x.__setitem__("physical_input_manifest_sha256", "0" * 64)),
        "array_storage_key": (lambda x: x["member_descriptors"][0]["state_object"]["physical_arrays"][0].__setitem__("npz_storage_key", "wrong")),
        "input_hash": (lambda x: x["member_descriptors"][0]["state_object"].__setitem__("input_hash", "0" * 64)),
        "zero_ledger": (lambda x: x["member_descriptors"][0]["initial_TDG6_ledger"].__setitem__("accepted_macro_step_count", 1)),
        "frozen_restart": (lambda x: x.__setitem__("restart_coordinate_time", time_identity(Fraction(3, 2)))),
        "foreign_protocol": (lambda x: x.__setitem__("protocol_artifact_id", "foreign")),
    }.items():
        clone = json.loads(canonical(fixture))
        altered(clone)
        try:
            build_genesis(clone)
        except Proto17ConstructionError:
            mutations[name] = True
        else:
            mutations[name] = False
    for name, altered in {
        "receipt_cycle_injection": lambda r, c: r["payload"].__setitem__("successor_checkpoint_sha256", c["checkpoint_sha256"]),
        "receipt_sequence": lambda r, _c: r["payload"].__setitem__("sequence", 2),
        "successor_untouched_method": lambda _r, c: c["cursors"]["RK4-2049"].__setitem__("method", "SSPRK3"),
        "successor_high_water": lambda _r, c: c["attempt_high_water"].__setitem__("RK4-2049", 99),
        "successor_checkpoint_hash": lambda _r, c: c.__setitem__("checkpoint_sha256", "0" * 64),
    }.items():
        receipt = json.loads(canonical(transition["receipt"]))
        successor = json.loads(canonical(transition["checkpoint"]))
        altered(receipt, successor)
        try:
            validate_common_event_successor(
                predecessor, prior_journal_records=prior_journal, receipt=receipt,
                successor_checkpoint=successor,
            )
        except Proto17ConstructionError:
            mutations[name] = True
        else:
            mutations[name] = False
    try:
        checkpoint(
            protocol_artifact_id=predecessor["protocol_artifact_id"], campaign_id=predecessor["campaign_id"],
            campaign_generation=1, committed_common_event_index=predecessor["committed_common_event_index"],
            active_event_target_time=predecessor["active_event_target_time"], cursors=predecessor["cursors"],
            ledgers=predecessor["ledgers"], states=predecessor["states"], journal_tip_sha256=ROOT_SHA256,
            attempt_high_water=predecessor["attempt_high_water"],
            cursor_generation_high_water=predecessor["cursor_generation_high_water"],
            parent_checkpoint_sha256=predecessor["parent_checkpoint_sha256"],
        )
    except Proto17ConstructionError:
        mutations["nonzero_generation_root_journal"] = True
    else:
        mutations["nonzero_generation_root_journal"] = False
    return {
        "genesis_checkpoint_sha256": genesis["checkpoint"]["checkpoint_sha256"],
        "common_event_record_sha256": transition["receipt"]["record_sha256"],
        "successor_checkpoint_sha256": transition["checkpoint"]["checkpoint_sha256"],
        "lone_suffix_reconstructs_byte_identically": recovered == transition,
        "mutations_rejected": mutations,
        "all_qualified_mutations_rejected": all(mutations.values()),
        "qualified_mutation_count": len(mutations),
        "qualified_mutation_names": sorted(mutations),
        "side_effect_free_reference_oracle_only": True,
        "historical_record_payload_semantics_deferred_to_runtime": True,
    }


__all__ = [
    "Proto17ConstructionError", "Proto17RecoveryInvalidStop", "build_genesis",
    "canonical", "checkpoint", "construct_common_event", "cursor", "decode_canonical",
    "derive_common_event_successor", "derive_generation_zero", "digest",
    "recover_lone_common_event", "run_pure_construction_qualification",
    "synthetic_common_event_fixture", "synthetic_common_event_predecessor",
    "synthetic_genesis_fixture", "time_identity", "validate_checkpoint",
    "validate_common_event_successor",
    "validate_ledger",
]
