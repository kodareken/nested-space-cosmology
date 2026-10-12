"""Synthetic-only durable runtime for the sealed PROTO17 construction.

This module intentionally accepts already verified GenesisSpec/state objects.
It owns temporary-store materialization, semantic accepted-record replay, the
new common-event receipt/checkpoint edge, and recovery.  Git authority and raw
NPZ verification are deliberately outside this bounded layer.
"""

from __future__ import annotations

from copy import deepcopy
from fractions import Fraction
from hashlib import sha256
import json
import os
from pathlib import Path
import shutil
import tempfile
from typing import Callable, Mapping, Any

from .proto17_pure_construction import (
    Proto17ConstructionError,
    build_genesis,
    canonical,
    checkpoint,
    construct_common_event,
    digest,
    time_identity,
    validate_checkpoint,
    validate_ledger,
)
from .protocol_v17 import MEMBER_KEYS, ROOT_SHA256


FaultHook = Callable[[str, str], None]
_PRODUCTION_ROOTS = (
    Path("runs/fgc-2-sf1/proto17/calibration"),
    Path("runs/fgc-2-sf1/proto17/holdout"),
)


class Proto17HLT14Error(ValueError):
    """A temporary trusted-genesis store violates the sealed construction."""


class Proto17HLT14RecoveryStop(Proto17HLT14Error):
    """Recovery found an unrecoverable suffix or semantic replay mismatch."""


def _sha(value: object, name: str) -> str:
    if not isinstance(value, str) or len(value) != 64 or value.lower() != value:
        raise Proto17HLT14Error(f"{name} is not a lowercase SHA-256 digest")
    try:
        int(value, 16)
    except ValueError as error:
        raise Proto17HLT14Error(f"{name} is not hexadecimal") from error
    return value


def _safe_root(root: Path) -> Path:
    answer = Path(root).resolve()
    for production in _PRODUCTION_ROOTS:
        if answer == production.resolve():
            raise Proto17HLT14Error("HLT14 refuses a production PROTO17 namespace")
    temporary = Path(tempfile.gettempdir()).resolve()
    try:
        answer.relative_to(temporary)
    except ValueError as error:
        raise Proto17HLT14Error("HLT14 accepts only caller-supplied temporary roots") from error
    return answer


def _hook(hook: FaultHook | None, phase: str, target: str) -> None:
    if hook is not None:
        hook(phase, target)


def _atomic(path: Path, payload: bytes, hook: FaultHook | None, target: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp-{os.getpid()}")
    try:
        _hook(hook, "before_write", target)
        with temporary.open("xb") as handle:
            handle.write(payload)
            _hook(hook, "after_write", target)
            _hook(hook, "before_flush", target)
            handle.flush()
            _hook(hook, "after_flush", target)
            _hook(hook, "before_fsync", target)
            os.fsync(handle.fileno())
            _hook(hook, "after_fsync", target)
        _hook(hook, "before_replace", target)
        os.replace(temporary, path)
        _hook(hook, "after_replace", target)
        descriptor = os.open(path.parent, os.O_RDONLY)
        try:
            _hook(hook, "before_dir_fsync", target)
            os.fsync(descriptor)
            _hook(hook, "after_dir_fsync", target)
        finally:
            os.close(descriptor)
    finally:
        temporary.unlink(missing_ok=True)


def _decode_canonical(payload: bytes, context: str) -> dict[str, Any]:
    def no_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        answer: dict[str, Any] = {}
        for key, value in pairs:
            if key in answer:
                raise ValueError(key)
            answer[key] = value
        return answer
    try:
        value = json.loads(payload.decode("utf-8"), object_pairs_hook=no_duplicates)
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as error:
        raise Proto17HLT14Error(f"{context} is malformed") from error
    if not isinstance(value, dict) or canonical(value) != payload:
        raise Proto17HLT14Error(f"{context} is noncanonical")
    return value


def _frame(record: Mapping[str, object]) -> bytes:
    encoded = canonical(dict(record))
    return f"{len(encoded):08x} ".encode("ascii") + encoded + b"\n"


def _unframe(payload: bytes) -> dict[str, Any]:
    try:
        prefix, rest = payload.split(b" ", 1)
        encoded, newline = rest[:-1], rest[-1:]
        if newline != b"\n" or len(prefix) != 8 or prefix.decode("ascii") != f"{len(encoded):08x}":
            raise ValueError
    except (ValueError, UnicodeDecodeError) as error:
        raise Proto17HLT14Error("journal frame is malformed") from error
    record = _decode_canonical(encoded, "journal record")
    if set(record) != {"record_kind", "journal_parent_sha256", "payload", "record_sha256"}:
        raise Proto17HLT14Error("journal record fields differ")
    bare = dict(record)
    observed = _sha(bare.pop("record_sha256"), "journal record")
    if digest(bare) != observed:
        raise Proto17HLT14Error("journal record digest differs")
    return record


class Proto17HLT14Runtime:
    """Temporary-root implementation of the pure PROTO17 lifecycle edge."""

    def __init__(self, genesis_spec: Mapping[str, object]) -> None:
        try:
            derived = build_genesis(genesis_spec)
        except (Proto17ConstructionError, TypeError, ValueError) as error:
            raise Proto17HLT14Error("verified GenesisSpec does not derive") from error
        self.spec = deepcopy(derived["genesis_spec"])
        self.genesis = deepcopy(derived["checkpoint"])
        self.store_plan = deepcopy(derived["store_plan"])

    @staticmethod
    def _state_path(root: Path, digest_value: str) -> Path:
        return root / "states" / f"{digest_value}.json"

    @staticmethod
    def _checkpoint_path(root: Path, value: Mapping[str, object]) -> Path:
        return root / "checkpoints" / f"{int(value['campaign_generation']):020d}-{value['checkpoint_sha256']}.json"

    @staticmethod
    def _record_path(root: Path, record: Mapping[str, object]) -> Path:
        return root / "journal" / f"{int(record['payload']['sequence']):020d}-{record['record_sha256']}.journal"

    def materialize_generation_zero(self, root: Path, *, fault_hook: FaultHook | None = None) -> dict[str, object]:
        root = _safe_root(root)
        if root.exists():
            raise Proto17HLT14Error("temporary trusted-genesis root must be absent")
        root.parent.mkdir(parents=True, exist_ok=True)
        staging = Path(tempfile.mkdtemp(prefix=f".{root.name}.hlt14-stage-", dir=root.parent))
        adopted = False
        try:
            for digest_value, state in self.store_plan["state_objects"].items():
                _sha(digest_value, "state object address")
                payload = canonical(state)
                if sha256(payload).hexdigest() != digest_value:
                    raise Proto17HLT14Error("derived state object content address differs")
                _atomic(self._state_path(staging, digest_value), payload, fault_hook, "state_object")
            checkpoint_value = deepcopy(self.genesis)
            self._write_checkpoint(staging, checkpoint_value, fault_hook, target="genesis_checkpoint")
            descriptor = os.open(staging, os.O_RDONLY)
            try:
                _hook(fault_hook, "before_namespace_dir_fsync", "namespace")
                os.fsync(descriptor)
                _hook(fault_hook, "after_namespace_dir_fsync", "namespace")
            finally:
                os.close(descriptor)
            _hook(fault_hook, "before_namespace_rename", "namespace")
            if root.exists():
                raise Proto17HLT14Error("temporary trusted-genesis target appeared during staging")
            os.replace(staging, root)
            adopted = True
            _hook(fault_hook, "after_namespace_rename", "namespace")
            parent = os.open(root.parent, os.O_RDONLY)
            try:
                _hook(fault_hook, "before_parent_dir_fsync", "namespace")
                os.fsync(parent)
                _hook(fault_hook, "after_parent_dir_fsync", "namespace")
            finally:
                os.close(parent)
            return checkpoint_value
        except BaseException:
            if adopted and root.exists():
                shutil.rmtree(root)
            if staging.exists():
                shutil.rmtree(staging)
            raise

    def _write_checkpoint(self, root: Path, value: Mapping[str, object], hook: FaultHook | None, *, target: str = "checkpoint") -> None:
        try:
            normalized = validate_checkpoint(value)
        except (Proto17ConstructionError, ValueError, TypeError) as error:
            raise Proto17HLT14Error("checkpoint is not a legal PROTO17 construction") from error
        if normalized != value:
            raise Proto17HLT14Error("checkpoint is not canonical construction bytes")
        path = self._checkpoint_path(root, normalized)
        if path.exists():
            if path.read_bytes() != canonical(normalized):
                raise Proto17HLT14Error("same checkpoint generation has forked")
            return
        if list((root / "checkpoints").glob(f"{int(normalized['campaign_generation']):020d}-*.json")):
            raise Proto17HLT14Error("same checkpoint generation has forked")
        _atomic(path, canonical(normalized), hook, target)

    def _records(self, root: Path) -> list[dict[str, Any]]:
        paths = sorted((root / "journal").glob("*.journal"))
        parent = ROOT_SHA256
        records: list[dict[str, Any]] = []
        for sequence, path in enumerate(paths, 1):
            record = _unframe(path.read_bytes())
            if (record["payload"].get("sequence") != sequence
                    or record["journal_parent_sha256"] != parent
                    or path.name != f"{sequence:020d}-{record['record_sha256']}.journal"):
                raise Proto17HLT14Error("journal chain differs")
            parent = record["record_sha256"]
            records.append(record)
        return records

    def _load_checkpoints(self, root: Path) -> list[dict[str, Any]]:
        paths = sorted((root / "checkpoints").glob("*.json"))
        if not paths:
            raise Proto17HLT14Error("trusted-genesis store has no checkpoint")
        checkpoints: list[dict[str, Any]] = []
        for generation, path in enumerate(paths):
            value = _decode_canonical(path.read_bytes(), "checkpoint")
            if value.get("campaign_generation") != generation or path.name != f"{generation:020d}-{value.get('checkpoint_sha256')}.json":
                raise Proto17HLT14Error("checkpoint filename/generation differs")
            try:
                if validate_checkpoint(value) != value:
                    raise Proto17HLT14Error("checkpoint differs from canonical construction")
            except Proto17ConstructionError as error:
                raise Proto17HLT14Error("checkpoint construction differs") from error
            checkpoints.append(value)
        return checkpoints

    def _verify_states(self, root: Path, checkpoint_value: Mapping[str, object]) -> None:
        for key in MEMBER_KEYS:
            address = _sha(checkpoint_value["states"][key], "checkpoint state")
            path = self._state_path(root, address)
            if not path.is_file() or sha256(path.read_bytes()).hexdigest() != address:
                raise Proto17HLT14Error("checkpoint state object is absent or drifted")
            _decode_canonical(path.read_bytes(), "state object")

    @staticmethod
    def _accepted_ledger(old: Mapping[str, object], target: Mapping[str, object]) -> dict[str, object]:
        result = deepcopy(old)
        result["last_accepted_time_hex"] = time_identity(target)["binary64_hex"]
        result["accepted_macro_step_count"] += 1
        result["last_accepted_macro_step_temporal_retry_count"] = result["current_macro_step_temporal_retry_count"]
        result["current_macro_step_temporal_retry_count"] = 0
        return validate_ledger(result, accepted_time=target)

    def _semantic_accepted_record(self, current: Mapping[str, object], key: str, sequence: int, parent_tip: str) -> tuple[dict[str, object], dict[str, object], dict[str, object]]:
        old_cursor = current["cursors"][key]
        old_ledger = current["ledgers"][key]
        target = current["active_event_target_time"]
        new_ledger = self._accepted_ledger(old_ledger, target)
        next_serial = old_cursor["attempt_serial"] + 1
        update = {**old_cursor,
            "accepted_boundary_time": target,
            "TDG6_ledger_sha256": digest(new_ledger),
            "cursor_generation": old_cursor["cursor_generation"] + 1,
            "attempt_serial": next_serial,
            "attempt_id": f"{current['campaign_id']}:{key}:{current['committed_common_event_index']}:{next_serial}",
            "journal_tip_sha256": parent_tip,
            "cursor_chain_parent_sha256": old_cursor["cursor_chain_sha256"],
        }
        update.pop("cursor_chain_sha256", None)
        from .proto17_pure_construction import cursor
        new_cursor = cursor(update)
        payload = {
            "member_key": key,
            "predecessor_cursor_chain_sha256": old_cursor["cursor_chain_sha256"],
            "successor_cursor_chain_sha256": new_cursor["cursor_chain_sha256"],
            "accepted_state_sha256": current["states"][key],
            "TDG6_ledger_sha256": new_cursor["TDG6_ledger_sha256"],
            "executed_plan": {"synthetic": True},
            "cursor": new_cursor,
            "sequence": sequence,
        }
        bare = {"record_kind": "ACCEPTED_FINE_COMMIT", "journal_parent_sha256": parent_tip, "payload": payload}
        return {**bare, "record_sha256": digest(bare)}, new_cursor, new_ledger

    def replay_six_semantic_accepted_records(self, root: Path, *, fault_hook: FaultHook | None = None) -> dict[str, object]:
        root = _safe_root(root)
        current = self.open(root)
        if current["checkpoint_sha256"] != self.genesis["checkpoint_sha256"] or self._records(root):
            raise Proto17HLT14Error("semantic accepted prefix requires untouched generation zero")
        base = deepcopy(current)
        working = deepcopy(current)
        for key in MEMBER_KEYS:
            sequence = len(self._records(root)) + 1
            parent = working["journal_tip_sha256"]
            record, new_cursor, new_ledger = self._semantic_accepted_record(working, key, sequence, parent)
            _atomic(self._record_path(root, record), _frame(record), fault_hook, "accepted_journal")
            working["cursors"][key] = new_cursor
            working["ledgers"][key] = new_ledger
            working["attempt_high_water"][key] = new_cursor["attempt_serial"]
            working["cursor_generation_high_water"][key] = new_cursor["cursor_generation"]
            working["journal_tip_sha256"] = record["record_sha256"]
        current = checkpoint(
            protocol_artifact_id=base["protocol_artifact_id"], campaign_id=base["campaign_id"],
            campaign_generation=base["campaign_generation"] + 1,
            committed_common_event_index=base["committed_common_event_index"],
            active_event_target_time=base["active_event_target_time"], cursors=working["cursors"],
            ledgers=working["ledgers"], states=base["states"], journal_tip_sha256=working["journal_tip_sha256"],
            attempt_high_water=working["attempt_high_water"], cursor_generation_high_water=working["cursor_generation_high_water"],
            parent_checkpoint_sha256=base["checkpoint_sha256"],
        )
        self._write_checkpoint(root, current, fault_hook)
        return current

    def commit_common_event(self, root: Path, *, fault_hook: FaultHook | None = None) -> dict[str, object]:
        root = _safe_root(root)
        current = self.open(root)
        records = self._records(root)
        try:
            edge = construct_common_event(current, prior_journal_records=records)
        except Proto17ConstructionError as error:
            raise Proto17HLT14Error("common event premises differ") from error
        receipt = edge["receipt"]
        _atomic(self._record_path(root, receipt), _frame(receipt), fault_hook, "common_event_receipt")
        self._write_checkpoint(root, edge["checkpoint"], fault_hook)
        return edge["checkpoint"]

    def _replay_committed(self, root: Path) -> tuple[dict[str, object], list[dict[str, Any]], int]:
        root = _safe_root(root)
        checkpoints = self._load_checkpoints(root)
        records = self._records(root)
        if checkpoints[0] != self.genesis:
            raise Proto17HLT14Error("store genesis differs from injected authority")
        for value in checkpoints:
            self._verify_states(root, value)
        # Semantic replay is intentionally restricted to the HLT14 synthetic
        # accepted prefix plus the PROTO17 receipt.  No production history is
        # admitted by this runtime layer.
        expected = deepcopy(self.genesis)
        position = 0
        for next_checkpoint in checkpoints[1:]:
            if position >= len(records):
                raise Proto17HLT14Error("checkpoint lacks its semantic journal record")
            record = records[position]
            if record["record_kind"] == "ACCEPTED_FINE_COMMIT":
                base = deepcopy(expected)
                working = deepcopy(expected)
                for key in MEMBER_KEYS:
                    if position >= len(records):
                        raise Proto17HLT14Error("accepted prefix is incomplete")
                    candidate = records[position]
                    predicted_record, predicted_cursor, predicted_ledger = self._semantic_accepted_record(
                        working, key, position + 1, working["journal_tip_sha256"]
                    )
                    if candidate != predicted_record:
                        raise Proto17HLT14Error("accepted record semantic payload differs")
                    working["cursors"][key] = predicted_cursor
                    working["ledgers"][key] = predicted_ledger
                    working["attempt_high_water"][key] = predicted_cursor["attempt_serial"]
                    working["cursor_generation_high_water"][key] = predicted_cursor["cursor_generation"]
                    working["journal_tip_sha256"] = candidate["record_sha256"]
                    position += 1
                expected = checkpoint(protocol_artifact_id=base["protocol_artifact_id"], campaign_id=base["campaign_id"], campaign_generation=base["campaign_generation"] + 1, committed_common_event_index=base["committed_common_event_index"], active_event_target_time=base["active_event_target_time"], cursors=working["cursors"], ledgers=working["ledgers"], states=base["states"], journal_tip_sha256=working["journal_tip_sha256"], attempt_high_water=working["attempt_high_water"], cursor_generation_high_water=working["cursor_generation_high_water"], parent_checkpoint_sha256=base["checkpoint_sha256"])
            elif record["record_kind"] == "COMMON_EVENT_COMMIT":
                edge = construct_common_event(expected, prior_journal_records=records[:position])
                if record != edge["receipt"]:
                    raise Proto17HLT14Error("common-event receipt differs")
                expected = edge["checkpoint"]
            else:
                raise Proto17HLT14Error("HLT14 synthetic replay record kind differs")
            if next_checkpoint != expected:
                raise Proto17HLT14Error("checkpoint differs from semantic replay")
            if record["record_kind"] == "COMMON_EVENT_COMMIT":
                position += 1
        return expected, records, position

    def open(self, root: Path) -> dict[str, object]:
        expected, records, position = self._replay_committed(root)
        if records[position:]:
            raise Proto17HLT14RecoveryStop("open refuses every uncheckpointed journal suffix")
        return expected

    def recover(self, root: Path, *, fault_hook: FaultHook | None = None) -> dict[str, object]:
        root = _safe_root(root)
        checkpoints = self._load_checkpoints(root)
        expected, records, position = self._replay_committed(root)
        suffix = records[position:]
        if not suffix:
            return expected
        if len(suffix) != 1 or suffix[0]["record_kind"] != "COMMON_EVENT_COMMIT":
            raise Proto17HLT14RecoveryStop("recovery suffix is not one common receipt")
        edge = construct_common_event(expected, prior_journal_records=records[:position])
        if suffix[0] != edge["receipt"]:
            raise Proto17HLT14RecoveryStop("recovery common receipt differs")
        if checkpoints[-1] != expected:
            raise Proto17HLT14RecoveryStop("latest checkpoint differs before suffix recovery")
        self._write_checkpoint(root, edge["checkpoint"], fault_hook)
        return edge["checkpoint"]


__all__ = ["FaultHook", "Proto17HLT14Error", "Proto17HLT14RecoveryStop", "Proto17HLT14Runtime"]
