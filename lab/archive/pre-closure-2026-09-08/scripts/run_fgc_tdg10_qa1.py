#!/usr/bin/env python3
"""Run one retry-3 SSPRK3 exact complete-C diagnostic qualification."""

from __future__ import annotations

import argparse
import ctypes
from dataclasses import asdict, is_dataclass
import errno
from fractions import Fraction
from hashlib import sha256
import importlib
import json
from math import isfinite
import os
from pathlib import Path
import secrets
import stat
import sys
from typing import Mapping, NoReturn, Sequence

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

from scripts import run_fgc_tdg9_ti2 as ti2_runner  # noqa: E402
from recursive_horizons.fgc.evolution import tdg6_temporal_admission_runtime as tdg6  # noqa: E402
from recursive_horizons.fgc.evolution import tdg9_ar1_authority as ar1  # noqa: E402
from recursive_horizons.fgc.evolution import tdg10_qa1_authority as authority  # noqa: E402
from recursive_horizons.fgc.evolution.hlt16_campaign_store import (  # noqa: E402
    HLT16CampaignStore,
)
from recursive_horizons.fgc.evolution.hlt16_member_codec import ARRAY_NAMES  # noqa: E402
from recursive_horizons.fgc.evolution.hlt16_state_store import (  # noqa: E402
    HLT16StateStoreError,
    decode_payload,
    encode_payload,
)
from recursive_horizons.fgc.evolution.numerical_engine import COMPARATOR_METHOD  # noqa: E402
from recursive_horizons.fgc.evolution.proto19_gr0_static_factory import (  # noqa: E402
    build_static_gr0_shells,
)
from recursive_horizons.fgc.evolution.tdg10_exact_complete_c_admission import (  # noqa: E402
    ExactCompleteCResourceExhausted,
    ExactCompleteCRouteDisagreement,
)
from recursive_horizons.fgc.evolution.tdg10_exact_complete_c_runtime import (  # noqa: E402
    TDG10ExactCompleteCRuntimeClosed,
)


RUNNER_ID = "FGC-1-TDG10-QA1-RUN1"
RAW_SCHEMA = "FGC-1-TDG10-QA1-raw-v1"
DIAGNOSTIC_ENDPOINT_SCHEMA = "FGC-1-TDG10-QA1-diagnostic-fine-endpoint-v1"
EXACT_COMPLETE_C_RUNTIME_MODULE = (
    "recursive_horizons.fgc.evolution.tdg10_exact_complete_c_runtime"
)
EXACT_COMPLETE_C_ASSESS_NAME = "assess_tdg10_exact_complete_c"
INCONCLUSIVE_EXCEPTIONS = (
    ExactCompleteCResourceExhausted,
    ExactCompleteCRouteDisagreement,
    TDG10ExactCompleteCRuntimeClosed,
)
PASS_CLASS = "completed_exact_all_channel_pass_no_state_advance"
NONPASS_CLASS = "completed_exact_all_channel_nonpass_no_state_advance"
INCONCLUSIVE_CLASS = "exact_route_or_resource_inconclusive_no_state_advance"
PREMISE_STOP_CLASS = "shadow_or_scientific_premise_stop_no_state_advance"
INVALID_CLASS = "invalid_provenance_or_implementation"
TERMINAL_CLASSES = frozenset(
    {
        PASS_CLASS,
        NONPASS_CLASS,
        INCONCLUSIVE_CLASS,
        PREMISE_STOP_CLASS,
        INVALID_CLASS,
    }
)
PUBLISHABLE_TERMINAL_CLASSES = TERMINAL_CLASSES - {INVALID_CLASS}
_OUTPUT_LEAF_LIMIT = 512 * 1024 * 1024
_HEX = frozenset("0123456789abcdef")
_RENAME_EXCL = 0x00000004
_FORBIDDEN_OUTPUT_NAMES = frozenset(
    {
        "checkpoints",
        "journal",
        "members",
        "cursor",
        "ledger",
        "campaign",
        "common-event",
        "common_event",
    }
)


class QA1RunnerError(RuntimeError):
    """Typed fail-closed QA1 runner error."""

    def __init__(self, owner: str, code: str, detail: object) -> None:
        self.owner = str(owner)
        self.code = str(code)
        self.detail = " ".join(str(detail).replace(str(ROOT), "<repo>").split())[:640]
        super().__init__(f"{self.owner}/{self.code}: {self.detail}")


def _fail(owner: str, code: str, detail: object) -> NoReturn:
    raise QA1RunnerError(owner, code, detail)


def _canonical(value: object) -> bytes:
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode()


def _pretty(value: object) -> bytes:
    return (
        json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n"
    ).encode()


def _unique(items: list[tuple[str, object]]) -> dict[str, object]:
    answer: dict[str, object] = {}
    for key, value in items:
        if key in answer:
            _fail("serialization", "duplicate_JSON_key", key)
        answer[key] = value
    return answer


def _fingerprint_exact(value: object) -> object:
    if isinstance(value, Fraction):
        return {
            "numerator": str(value.numerator),
            "denominator": str(value.denominator),
        }
    if type(value) is float:
        if not isfinite(value):
            _fail("serialization", "nonfinite_binary64", value)
        return {"binary64_hex": value.hex()}
    if is_dataclass(value) and not isinstance(value, type):
        return _fingerprint_exact(asdict(value))
    if isinstance(value, Mapping):
        return {str(key): _fingerprint_exact(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_fingerprint_exact(item) for item in value]
    if value is None or type(value) in {str, int, bool}:
        return value
    _fail("serialization", "unsupported_value", type(value).__name__)


def _fraction(value: object, *, label: str) -> Fraction:
    if isinstance(value, Fraction) and not isinstance(value, bool):
        return value
    if isinstance(value, Mapping) and set(value) == {"numerator", "denominator"}:
        numerator = value["numerator"]
        denominator = value["denominator"]
        if not isinstance(numerator, str) or not isinstance(denominator, str):
            _fail("serialization", "fraction_types", label)
        try:
            answer = Fraction(int(numerator), int(denominator))
        except (TypeError, ValueError, ZeroDivisionError) as exc:
            raise QA1RunnerError("serialization", "fraction_parse", label) from exc
        if (
            str(answer.numerator) != numerator
            or str(answer.denominator) != denominator
        ):
            _fail("serialization", "fraction_roundtrip", label)
        return answer
    _fail("serialization", "nonexact_fraction", label)


def _hex_digest(value: object, *, label: str) -> str:
    if (
        not isinstance(value, str)
        or len(value) != 64
        or value.lower() != value
        or any(character not in _HEX for character in value)
    ):
        _fail("serialization", "digest", label)
    return value


def _require_keys(
    value: object, expected: set[str], *, label: str
) -> Mapping[str, object]:
    if not isinstance(value, Mapping) or set(value) != expected:
        _fail("status", "terminal_schema", (label, sorted(expected)))
    return value


def _nonclaims() -> dict[str, object]:
    return {
        "diagnostic_qualification_only": True,
        "diagnostic_fine_endpoint_is_accepted_state": False,
        "new_protocol_id_authorized": False,
        "target_protocol": authority.TARGET_PROTOCOL,
        "production_SSPRK3_comparator_earned": False,
        "independent_method_agreement_earned": False,
        "state_advance_authorized": False,
        "campaign_state_write_authorized": False,
        "PDE_state_commit_authorized": False,
        "PDE_state_committed": False,
        "fine_path_commit_authorized": False,
        "fine_path_committed": False,
        "temporal_retry_admission_called": False,
        "candidate_execution_authorized": False,
        "candidate_branch_opened": False,
        "mechanism_result_authorized": False,
        "physical_result_authorized": False,
        "common_event_authorized": False,
        "GR0_calibration_authorized": False,
        "SGBL_authorized": False,
        "FGCQR_authorized": False,
        "retained_EFT_evolution_authorized": False,
        "physical_transition_claim_authorized": False,
        "external_publication_authorized": False,
        "push_authorized": False,
        "retry_4_authorized": False,
        "retry_5_authorized": False,
        "successor_remedy_selected": False,
        "common_event_completed": False,
        "GR0_calibration_completed": False,
        "mechanism_result_earned": False,
        "physical_result_earned": False,
        "production_method_earned": False,
    }


def _execution_authority(root: Path, commit: str) -> authority.QA1Authority:
    try:
        return authority.authorize(
            root,
            authority.read_leaf(root, authority.CONFIG_PATH),
            authority.read_leaf(root, authority.RESULT_PATH),
            commit,
        )
    except Exception as exc:
        raise QA1RunnerError("authority", "rejected", exc) from exc


def _snapshot_store(root: Path) -> tuple[int, str]:
    try:
        return ti2_runner._snapshot_store(root)
    except Exception as exc:
        raise QA1RunnerError("provenance", "store_snapshot", exc) from exc


def _restore_replay(root: Path) -> ti2_runner.RestoredReplay:
    replay = dict(authority.REPLAY)
    if int(replay["predecessor_generation"]) == authority.FORBIDDEN_RETRY3_GENERATION:
        _fail("replay", "retry3_must_not_restore_generation_10", 10)
    if int(replay["predecessor_generation"]) != authority.PREDECESSOR_GENERATION:
        _fail("replay", "generation_not_9", replay["predecessor_generation"])
    if int(replay["journal_sequence"]) != authority.HISTORICAL_RETRY3_REJECTION_SEQUENCE:
        _fail("replay", "historical_rejection_not_sequence_11", replay["journal_sequence"])
    store = HLT16CampaignStore(root / ar1.PREF2_STORE_PATH)
    shells = build_static_gr0_shells(root)
    try:
        return ti2_runner._restore_replay(root, store, shells, replay)
    except ti2_runner.TI1RunnerError as exc:
        raise QA1RunnerError(exc.owner, exc.code, exc.detail) from exc


def _prepare_shadow(
    restored: ti2_runner.RestoredReplay,
) -> tdg6.TDG6PreparedGR0Compositor:
    try:
        return ti2_runner._prepare_shadow(restored)
    except tdg6.TDG6RefinementPathStop:
        raise
    except ti2_runner.TI1RunnerError as exc:
        raise QA1RunnerError(exc.owner, exc.code, exc.detail) from exc


def _assess_prepared_exact_complete_c(prepared: object) -> object:
    """Isolated adapter import; adjust EXACT_COMPLETE_C_* if the public name differs."""

    try:
        module = importlib.import_module(EXACT_COMPLETE_C_RUNTIME_MODULE)
        assess = getattr(module, EXACT_COMPLETE_C_ASSESS_NAME)
    except (ModuleNotFoundError, AttributeError) as exc:
        raise QA1RunnerError("exact", "adapter_unavailable", exc) from exc
    return assess(
        prepared,
        maximum_candidates_D01=authority.PER_CHANNEL_MAXIMUM_CANDIDATES_D01,
        maximum_candidates_D12=authority.PER_CHANNEL_MAXIMUM_CANDIDATES_D12,
        refinement_depth=authority.PRIMARY_REFINEMENT_DEPTH,
        expected_owned_row_count=authority.OWNED_ROW_COUNT,
    )


def expected_replay_receipt() -> dict[str, object]:
    return authority.expected_raw_replay_receipt()


def _replay_receipt(restored: ti2_runner.RestoredReplay) -> dict[str, object]:
    fingerprint = restored.fingerprint
    if not isinstance(fingerprint, Mapping):
        _fail("replay", "fingerprint_type", type(fingerprint).__name__)
    receipt = {
        "member_key": fingerprint.get("member_key"),
        "retry": authority.RETRY,
        "predecessor_generation": authority.PREDECESSOR_GENERATION,
        "checkpoint_journal_sequence": authority.CHECKPOINT_JOURNAL_SEQUENCE,
        "historical_retry3_rejection_sequence": (
            authority.HISTORICAL_RETRY3_REJECTION_SEQUENCE
        ),
        "attempted_width_hex": authority.WIDTH_HEX,
        "accepted_time_hex": fingerprint.get("accepted_time_hex"),
        "state_sha256": fingerprint.get("state_sha256"),
        "descriptor_sha256": fingerprint.get("descriptor_sha256"),
        "transaction_sha256": fingerprint.get("transaction_sha256"),
        "historical_journal_sha256": restored.historical_journal_sha256,
    }
    expected = expected_replay_receipt()
    if receipt != expected:
        _fail("replay", "receipt_identity", receipt)
    return receipt


def _channel_records(assessment: object) -> list[tuple[str, object]]:
    channels = getattr(assessment, "channels", None)
    if channels is None:
        channels = getattr(assessment, "channel_evidence", None)
    if isinstance(assessment, Mapping) and channels is None:
        channels = assessment.get("channels", assessment.get("channel_evidence"))
    if isinstance(channels, Mapping):
        return [(name, channels[name]) for name in authority.CHANNEL_ORDER]
    if isinstance(channels, (list, tuple)):
        records: list[tuple[str, object]] = []
        for index, item in enumerate(channels):
            if isinstance(item, tuple) and len(item) == 2:
                records.append((str(item[0]), item[1]))
            elif isinstance(item, Mapping) and "channel" in item:
                records.append((str(item["channel"]), item))
            elif hasattr(item, "channel"):
                records.append((str(item.channel), item))
            else:
                records.append((authority.CHANNEL_ORDER[index], item))
        return records
    _fail("exact", "channel_records", type(assessment).__name__)


def _nested_evidence(item: object) -> object:
    if isinstance(item, Mapping):
        return item.get("evidence", item)
    return getattr(item, "evidence", item)


def _d12_upper(item: object) -> Fraction:
    if isinstance(item, Mapping) and "d12_upper" in item:
        return _fraction(item["d12_upper"], label="d12_upper")
    evidence = _nested_evidence(item)
    d12 = (
        evidence.get("d12")
        if isinstance(evidence, Mapping)
        else getattr(evidence, "d12", None)
    )
    if d12 is None:
        d12 = item.get("d12") if isinstance(item, Mapping) else getattr(item, "d12", None)
    if d12 is not None:
        upper = d12.get("upper") if isinstance(d12, Mapping) else getattr(d12, "upper", None)
        if upper is not None:
            return _fraction(upper, label="d12.upper")
    _fail("exact", "d12_upper_absent", type(item).__name__)


def _admission_passed(item: object) -> bool:
    if isinstance(item, Mapping) and "admission_passed" in item:
        return item["admission_passed"] is True
    if hasattr(item, "admission_passed"):
        return item.admission_passed is True
    evidence = _nested_evidence(item)
    decision = (
        evidence.get("decision")
        if isinstance(evidence, Mapping)
        else getattr(evidence, "decision", None)
    )
    if decision is None:
        decision = getattr(item, "decision", None)
    if decision is not None:
        flag = (
            decision.get("admission_passed")
            if isinstance(decision, Mapping)
            else getattr(decision, "admission_passed", False)
        )
        return flag is True
    _fail("exact", "admission_flag_absent", type(item).__name__)


def _classification_of(item: object) -> str:
    if isinstance(item, Mapping) and "classification" in item:
        return str(item["classification"])
    if hasattr(item, "classification"):
        return str(item.classification)
    evidence = _nested_evidence(item)
    decision = (
        evidence.get("decision")
        if isinstance(evidence, Mapping)
        else getattr(evidence, "decision", None)
    )
    if decision is None:
        decision = getattr(item, "decision", None)
    if decision is not None:
        value = (
            decision.get("classification")
            if isinstance(decision, Mapping)
            else getattr(decision, "classification", None)
        )
        if value is not None:
            return str(value)
    _fail("exact", "classification_absent", type(item).__name__)


def serialize_exact_assessment(assessment: object) -> dict[str, object]:
    records = _channel_records(assessment)
    if tuple(name for name, _item in records) != authority.CHANNEL_ORDER:
        _fail("exact", "channel_order", [name for name, _item in records])
    channels: list[dict[str, object]] = []
    failed: list[str] = []
    for name, item in records:
        passed = _admission_passed(item)
        if not passed:
            failed.append(name)
        channels.append(
            {
                "channel": name,
                "admission_passed": passed,
                "classification": _classification_of(item),
                "d12_upper": _fingerprint_exact(_d12_upper(item)),
            }
        )
    declared_complete = getattr(assessment, "complete_admission_passed", None)
    if declared_complete is None and isinstance(assessment, Mapping):
        declared_complete = assessment.get("complete_admission_passed")
    complete = not failed
    if declared_complete is not None and declared_complete is not complete:
        _fail("exact", "complete_pass_identity", declared_complete)
    declared_failed = getattr(assessment, "failed_channels", None)
    if declared_failed is None and isinstance(assessment, Mapping):
        declared_failed = assessment.get("failed_channels")
    if declared_failed is not None and tuple(declared_failed) != tuple(failed):
        _fail("exact", "failed_channels_identity", declared_failed)
    return {
        "channel_order": list(authority.CHANNEL_ORDER),
        "channel_count": authority.CHANNEL_COUNT,
        "complete_admission_passed": complete,
        "failed_channels": failed,
        "admission_is_all_of": True,
        "interval_owner": authority.INTERVAL_OWNER,
        "independent_route_required": True,
        "channels": channels,
    }


def reduce_qa1_terminal(
    *,
    complete_admission_passed: bool,
    failed_channels: Sequence[str],
) -> str:
    if complete_admission_passed is True and tuple(failed_channels) == ():
        return PASS_CLASS
    return NONPASS_CLASS


def _terminal_base(
    authority_commit: str,
    classification: str,
    before: tuple[int, str],
    after: tuple[int, str],
) -> dict[str, object]:
    if classification not in TERMINAL_CLASSES:
        _fail("reduction", "terminal_class", classification)
    return {
        "schema": RAW_SCHEMA,
        "artifact_id": authority.ARTIFACT_ID,
        "runner_id": RUNNER_ID,
        "classification": classification,
        "authority_commit": authority_commit,
        "tableau_selector": authority.TABLEAU_SELECTOR,
        "tableau_runtime_selector": COMPARATOR_METHOD,
        "actual_spatial_operator": authority.ACTUAL_SPATIAL_OPERATOR,
        "retry": authority.RETRY,
        "predecessor_generation": authority.PREDECESSOR_GENERATION,
        "checkpoint_journal_sequence": authority.CHECKPOINT_JOURNAL_SEQUENCE,
        "historical_retry3_rejection_sequence": (
            authority.HISTORICAL_RETRY3_REJECTION_SEQUENCE
        ),
        "attempted_width_hex": authority.WIDTH_HEX,
        "channel_order": list(authority.CHANNEL_ORDER),
        "store_snapshot_before": {"leaf_count": before[0], "sha256": before[1]},
        "store_snapshot_after": {"leaf_count": after[0], "sha256": after[1]},
        "store_unchanged": before == after,
        **_nonclaims(),
    }


def _manifest(
    authority_commit: str,
    *,
    output_leaves: Sequence[str] | None = None,
) -> dict[str, object]:
    leaves = list(output_leaves or ("manifest.json", "terminal.json"))
    return {
        "schema": RAW_SCHEMA,
        "artifact_id": authority.ARTIFACT_ID,
        "runner_id": RUNNER_ID,
        "authority_commit": authority_commit,
        "tableau_selector": authority.TABLEAU_SELECTOR,
        "actual_spatial_operator": authority.ACTUAL_SPATIAL_OPERATOR,
        "retry": authority.RETRY,
        "predecessor_generation": authority.PREDECESSOR_GENERATION,
        "target_protocol": authority.TARGET_PROTOCOL,
        "diagnostic_qualification_only": True,
        "diagnostic_fine_endpoint_is_accepted_state": False,
        "state_advance_authorized": False,
        "output_leaves": leaves,
    }


def _array(value: object, *, label: str) -> np.ndarray:
    array = np.asarray(value, dtype=np.float64)
    if array.dtype != np.dtype("<f8") or not array.flags.c_contiguous:
        array = np.ascontiguousarray(array, dtype=np.float64)
    if not array.shape or not np.isfinite(array).all():
        _fail("diagnostic", "array_invalid", label)
    return array


def _reconstruct_fine_tracers(
    member: object, prepared: object
) -> tuple[object, object]:
    fine = prepared.fine
    clone = tdg6._clone_tracers(member.tracers)
    attempts = tuple(fine.attempts)
    if len(attempts) != 4:
        _fail("diagnostic", "fine_attempt_count", len(attempts))
    for attempt in attempts:
        payload = getattr(attempt, "preaccept_payload", None)
        if payload is None:
            _fail("diagnostic", "preaccept_payload_absent", type(attempt).__name__)
        clone.commit_advance(*payload)
    snapshot = tdg6._tracer_snapshot(clone)
    expected = fine.final_tracer_snapshot
    if snapshot != expected:
        _fail("diagnostic", "tracer_snapshot_mismatch", snapshot)
    if not np.array_equal(
        np.asarray(clone.positions, dtype=np.float64),
        np.asarray(fine.final_tracer_positions, dtype=np.float64),
    ) or not np.array_equal(
        np.asarray(clone.proper_times, dtype=np.float64),
        np.asarray(fine.final_tracer_proper_times, dtype=np.float64),
    ):
        _fail("diagnostic", "tracer_array_mismatch", "positions_or_proper_times")
    return clone, snapshot


def _executed_four_quarter_plan(prepared: object) -> dict[str, object]:
    fine = prepared.fine
    attempts = tuple(fine.attempts)
    boundaries = [float(fine.initial_time)]
    for attempt in attempts:
        proposal = attempt.proposal
        boundaries.append(float(proposal.final_time))
    if len(boundaries) != 5:
        _fail("diagnostic", "four_quarter_boundaries", len(boundaries))
    return {
        "level": "fine",
        "step_count": 4,
        "initial_time": _fingerprint_exact(float(fine.initial_time)),
        "final_time": _fingerprint_exact(float(fine.final_time)),
        "boundaries": [_fingerprint_exact(value) for value in boundaries],
        "quarter_widths": [
            _fingerprint_exact(right - left)
            for left, right in zip(boundaries, boundaries[1:])
        ],
    }


def _payload_arrays(member: object, accepted: object, tracers: object) -> dict[str, np.ndarray]:
    arrays = {
        "u": _array(accepted.state.u, label="u"),
        "p": _array(accepted.state.p, label="p"),
        "q": _array(accepted.state.q, label="q"),
        "grid_coordinates": _array(
            member.initial.grid.coordinates, label="grid_coordinates"
        ),
        "tracer_labels": _array(tracers.labels, label="tracer_labels"),
        "tracer_positions": _array(tracers.positions, label="tracer_positions"),
        "tracer_proper_times": _array(
            tracers.proper_times, label="tracer_proper_times"
        ),
        "event_proper_times": _array(
            np.stack(tracers.event_proper_times), label="event_proper_times"
        ),
        "event_fields": _array(
            np.stack(tracers.event_fields), label="event_fields"
        ),
    }
    if tuple(arrays) != ARRAY_NAMES:
        _fail("diagnostic", "array_order", tuple(arrays))
    return arrays


def serialize_diagnostic_endpoint(
    restored: ti2_runner.RestoredReplay,
    prepared: object,
    exact: Mapping[str, object],
) -> tuple[dict[str, bytes], dict[str, object]]:
    member = restored.member
    after = ti2_runner._member_fingerprint(
        member, str(restored.fingerprint["descriptor_sha256"])
    )
    if after != restored.fingerprint:
        _fail("diagnostic", "member_mutated", after)
    clone, snapshot = _reconstruct_fine_tracers(member, prepared)
    accepted = prepared.fine.final_accepted
    arrays = _payload_arrays(member, accepted, clone)
    try:
        raw, semantic, manifest = encode_payload(arrays)
    except Exception as exc:
        raise QA1RunnerError("diagnostic", "payload_encode", exc) from exc
    semantic = _hex_digest(semantic, label="payload_semantic")
    payload_raw_sha256 = sha256(raw).hexdigest()
    d12_vector = [
        {
            "channel": item["channel"],
            "upper": item["d12_upper"],
        }
        for item in exact["channels"]
        if isinstance(item, Mapping)
    ]
    body = {
        "schema": DIAGNOSTIC_ENDPOINT_SCHEMA,
        "artifact_id": authority.ARTIFACT_ID,
        "classification": "diagnostic_nonaccepted_complete_fine_endpoint",
        "diagnostic_qualification_only": True,
        "accepted_state": False,
        "proto15_cursor": False,
        "campaign_descriptor": False,
        "journal": False,
        "checkpoint": False,
        "ledger": False,
        "payload_semantic_sha256": semantic,
        "payload_raw_sha256": payload_raw_sha256,
        "payload_manifest": [dict(item) for item in manifest],
        "predecessor": {
            "member_key": authority.MEMBER_KEY,
            "generation": authority.PREDECESSOR_GENERATION,
            "checkpoint_sha256": authority.CHECKPOINT_SHA256,
            "descriptor_sha256": authority.MEMBER_DESCRIPTOR_SHA256,
            "physical_state_sha256": authority.PHYSICAL_STATE_SHA256,
            "accepted_time": _fingerprint_exact(
                float.fromhex(authority.ACCEPTED_TIME_HEX)
            ),
            "transaction_sha256": authority.TRANSACTION_SHA256,
        },
        "accepted_step_tableau": authority.TABLEAU_SELECTOR,
        "actual_spatial_operator": authority.ACTUAL_SPATIAL_OPERATOR,
        "endpoint_time": _fingerprint_exact(float(accepted.time)),
        "endpoint_step_index": int(accepted.step_index),
        "endpoint_transaction_serial": int(accepted.transaction_serial),
        "executed_four_quarter_plan": _executed_four_quarter_plan(prepared),
        "monitor_state": _fingerprint_exact(prepared.fine.final_monitor_state),
        "causal_state": _fingerprint_exact(prepared.fine.final_causal_state),
        "tracer_snapshot": _fingerprint_exact(snapshot),
        "exact_D12_upper_vector": d12_vector,
        **{
            key: False
            for key, value in _nonclaims().items()
            if value is False
        },
        "diagnostic_qualification_only": True,
        "target_protocol": authority.TARGET_PROTOCOL,
    }
    descriptor_sha256 = sha256(_canonical(body)).hexdigest()
    body["descriptor_sha256"] = descriptor_sha256
    payload_relative = f"payloads/{semantic}.npz"
    descriptor_relative = f"fine-endpoint/{descriptor_sha256}.json"
    extras = {
        payload_relative: raw,
        descriptor_relative: _pretty(body),
    }
    record = {
        "present": True,
        "accepted_state": False,
        "payload_relative": payload_relative,
        "payload_semantic_sha256": semantic,
        "payload_raw_sha256": payload_raw_sha256,
        "payload_manifest": [dict(item) for item in manifest],
        "descriptor_relative": descriptor_relative,
        "descriptor_sha256": descriptor_sha256,
    }
    return extras, record


def _open_output_parent(root: Path) -> tuple[int, str]:
    relative = Path(authority.OUTPUT_NAMESPACE)
    descriptor = os.open(
        root,
        os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0),
    )
    try:
        for part in relative.parent.parts:
            try:
                child = os.open(
                    part,
                    os.O_RDONLY
                    | getattr(os, "O_DIRECTORY", 0)
                    | getattr(os, "O_NOFOLLOW", 0),
                    dir_fd=descriptor,
                )
            except FileNotFoundError:
                os.mkdir(part, mode=0o755, dir_fd=descriptor)
                child = os.open(
                    part,
                    os.O_RDONLY
                    | getattr(os, "O_DIRECTORY", 0)
                    | getattr(os, "O_NOFOLLOW", 0),
                    dir_fd=descriptor,
                )
            os.close(descriptor)
            descriptor = child
        return descriptor, relative.name
    except Exception:
        os.close(descriptor)
        raise


def _write_leaf(directory_fd: int, name: str, raw: bytes) -> None:
    leaf = os.open(
        name,
        os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0),
        0o644,
        dir_fd=directory_fd,
    )
    try:
        pending = memoryview(raw)
        while pending:
            written = os.write(leaf, pending)
            if written <= 0:
                _fail("output", "short_write", name)
            pending = pending[written:]
        os.fsync(leaf)
    finally:
        os.close(leaf)


def _cleanup_stage(
    stage_fd: int | None,
    stage_name: str,
    parent_fd: int,
    extras: Mapping[str, bytes],
) -> None:
    if stage_fd is None:
        return
    for name in ("manifest.json", "terminal.json"):
        try:
            os.unlink(name, dir_fd=stage_fd)
        except FileNotFoundError:
            pass
    for relative in extras:
        path = Path(relative)
        if len(path.parts) == 2:
            try:
                child = os.open(
                    path.parts[0],
                    os.O_RDONLY
                    | getattr(os, "O_DIRECTORY", 0)
                    | getattr(os, "O_NOFOLLOW", 0),
                    dir_fd=stage_fd,
                )
            except FileNotFoundError:
                continue
            try:
                try:
                    os.unlink(path.parts[1], dir_fd=child)
                except FileNotFoundError:
                    pass
            finally:
                os.close(child)
            try:
                os.rmdir(path.parts[0], dir_fd=stage_fd)
            except OSError:
                pass
    os.close(stage_fd)
    try:
        os.rmdir(stage_name, dir_fd=parent_fd)
    except FileNotFoundError:
        pass


def _publish(
    root: Path,
    manifest: Mapping[str, object],
    terminal: Mapping[str, object],
    *,
    extra_files: Mapping[str, bytes] | None = None,
) -> None:
    extras = dict(extra_files or {})
    parent_fd, target_name = _open_output_parent(root)
    stage_name = f"{authority.STAGING_PREFIX}{os.getpid()}-{secrets.token_hex(8)}"
    stage_fd: int | None = None
    try:
        try:
            os.stat(target_name, dir_fd=parent_fd, follow_symlinks=False)
        except FileNotFoundError:
            pass
        else:
            _fail("output", "namespace_exists", authority.OUTPUT_NAMESPACE)
        os.mkdir(stage_name, mode=0o700, dir_fd=parent_fd)
        stage_fd = os.open(
            stage_name,
            os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0),
            dir_fd=parent_fd,
        )
        _write_leaf(stage_fd, "manifest.json", _pretty(manifest))
        _write_leaf(stage_fd, "terminal.json", _pretty(terminal))
        created_dirs: set[str] = set()
        for relative, raw in extras.items():
            path = Path(relative)
            if (
                path.is_absolute()
                or len(path.parts) != 2
                or any(part in {"", ".", ".."} for part in path.parts)
                or path.parts[0] in _FORBIDDEN_OUTPUT_NAMES
            ):
                _fail("output", "unsafe_extra_leaf", relative)
            directory = path.parts[0]
            if directory not in created_dirs:
                os.mkdir(directory, mode=0o700, dir_fd=stage_fd)
                created_dirs.add(directory)
            child = os.open(
                directory,
                os.O_RDONLY
                | getattr(os, "O_DIRECTORY", 0)
                | getattr(os, "O_NOFOLLOW", 0),
                dir_fd=stage_fd,
            )
            try:
                _write_leaf(child, path.parts[1], raw)
                os.fsync(child)
            finally:
                os.close(child)
        os.fsync(stage_fd)
        if os.uname().sysname != "Darwin":
            _fail("output", "exclusive_adoption_unavailable", os.uname().sysname)
        rename = getattr(ctypes.CDLL(None, use_errno=True), "renameatx_np", None)
        if rename is None:
            _fail("output", "exclusive_adoption_unavailable", "renameatx_np")
        rename.argtypes = (
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.c_uint,
        )
        rename.restype = ctypes.c_int
        if rename(
            parent_fd,
            os.fsencode(stage_name),
            parent_fd,
            os.fsencode(target_name),
            _RENAME_EXCL,
        ):
            error = ctypes.get_errno()
            if error in {errno.EEXIST, errno.ENOTEMPTY}:
                _fail("output", "namespace_arrived", authority.OUTPUT_NAMESPACE)
            _fail("output", "exclusive_adoption_failed", error)
        os.close(stage_fd)
        stage_fd = None
        os.fsync(parent_fd)
    except Exception:
        _cleanup_stage(stage_fd, stage_name, parent_fd, extras)
        stage_fd = None
        raise
    finally:
        os.close(parent_fd)


def _publish_terminal(
    root: Path,
    manifest: Mapping[str, object],
    terminal: Mapping[str, object],
    *,
    store_before: tuple[int, str],
    store_after: tuple[int, str],
    extra_files: Mapping[str, bytes] | None = None,
) -> None:
    if store_after != store_before:
        _fail("provenance", "campaign_store_mutated", (store_before, store_after))
    classification = terminal.get("classification")
    if classification not in PUBLISHABLE_TERMINAL_CLASSES:
        _fail("output", "invalid_not_published", classification)
    if extra_files and classification != PASS_CLASS:
        _fail("output", "diagnostic_without_all_pass", classification)
    if classification == PASS_CLASS and not extra_files:
        _fail("output", "all_pass_missing_diagnostic", classification)
    _publish(root, manifest, terminal, extra_files=extra_files)


def _identity(value: os.stat_result) -> tuple[int, int, int, int, int, int]:
    return (
        value.st_dev,
        value.st_ino,
        value.st_size,
        value.st_mtime_ns,
        value.st_ctime_ns,
        value.st_nlink,
    )


def _read_output_leaf(directory_fd: int, name: str, before: os.stat_result) -> bytes:
    if (
        stat.S_ISLNK(before.st_mode)
        or not stat.S_ISREG(before.st_mode)
        or before.st_nlink != 1
        or before.st_size < 0
        or before.st_size > _OUTPUT_LEAF_LIMIT
    ):
        _fail("status", "unsafe_output_leaf", name)
    descriptor = -1
    try:
        descriptor = os.open(
            name,
            os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0),
            dir_fd=directory_fd,
        )
        active = os.fstat(descriptor)
        expected = _identity(before)
        if _identity(active) != expected:
            _fail("status", "output_leaf_substituted", name)
        chunks: list[bytes] = []
        remaining = before.st_size
        while remaining:
            block = os.read(descriptor, min(1 << 20, remaining))
            if not block:
                _fail("status", "output_leaf_short_read", name)
            chunks.append(block)
            remaining -= len(block)
        if os.read(descriptor, 1):
            _fail("status", "output_leaf_grew", name)
        if _identity(os.fstat(descriptor)) != expected:
            _fail("status", "output_leaf_changed", name)
        after = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
        if _identity(after) != expected:
            _fail("status", "output_leaf_substituted", name)
        return b"".join(chunks)
    except OSError as exc:
        raise QA1RunnerError("status", "unsafe_output_leaf", name) from exc
    finally:
        if descriptor != -1:
            os.close(descriptor)


def _parse_output_json(raw: bytes, *, name: str) -> dict[str, object]:
    def reject_constant(token: str) -> NoReturn:
        _fail("status", "nonfinite_json_constant", (name, token))

    try:
        value = json.loads(
            raw,
            object_pairs_hook=_unique,
            parse_constant=reject_constant,
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise QA1RunnerError("status", "noncanonical_output", name) from exc
    if not isinstance(value, dict) or raw != _pretty(value):
        _fail("status", "noncanonical_output", name)
    return value


def _snapshot_output(root: Path) -> tuple[dict[str, object], dict[str, bytes], dict[str, str]] | None:
    flags = (
        os.O_RDONLY
        | getattr(os, "O_DIRECTORY", 0)
        | getattr(os, "O_NOFOLLOW", 0)
    )
    root_before = root.lstat()
    if stat.S_ISLNK(root_before.st_mode) or not stat.S_ISDIR(root_before.st_mode):
        _fail("status", "unsafe_repository_root", root)
    descriptors: list[int] = []
    identities: list[tuple[int, int, int, int, int, int]] = []
    edge_names: list[str] = []
    try:
        current = os.open(root, flags)
        descriptors.append(current)
        identities.append(_identity(os.fstat(current)))
        if identities[0] != _identity(root_before):
            _fail("status", "repository_root_substituted", root)
        for part in Path(authority.OUTPUT_NAMESPACE).parts:
            try:
                before = os.stat(part, dir_fd=current, follow_symlinks=False)
            except FileNotFoundError:
                return None
            if stat.S_ISLNK(before.st_mode) or not stat.S_ISDIR(before.st_mode):
                _fail("status", "unsafe_namespace", part)
            child = os.open(part, flags, dir_fd=current)
            if _identity(os.fstat(child)) != _identity(before):
                os.close(child)
                _fail("status", "output_directory_substituted", part)
            descriptors.append(child)
            identities.append(_identity(before))
            edge_names.append(part)
            current = child
        output_fd = descriptors[-1]
        names = tuple(sorted(os.listdir(output_fd)))
        if any(name in _FORBIDDEN_OUTPUT_NAMES for name in names):
            _fail("status", "campaign_artifact_present", names)
        values: dict[str, object] = {}
        blobs: dict[str, bytes] = {}
        hashes: dict[str, str] = {}

        def record_file(directory_fd: int, relative: str, name: str) -> None:
            before = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
            raw = _read_output_leaf(directory_fd, name, before)
            blobs[relative] = raw
            hashes[relative] = sha256(raw).hexdigest()
            if relative.endswith(".json"):
                values[relative] = _parse_output_json(raw, name=relative)

        json_names = ("manifest.json", "terminal.json")
        if not set(json_names).issubset(names):
            _fail("status", "partial_or_foreign_namespace", names)
        for name in names:
            before = os.stat(name, dir_fd=output_fd, follow_symlinks=False)
            if stat.S_ISLNK(before.st_mode):
                _fail("status", "unsafe_namespace", name)
            if stat.S_ISREG(before.st_mode):
                if name not in json_names:
                    _fail("status", "partial_or_foreign_namespace", names)
                record_file(output_fd, name, name)
            elif stat.S_ISDIR(before.st_mode):
                if name not in {"payloads", "fine-endpoint"}:
                    _fail("status", "partial_or_foreign_namespace", names)
                child = os.open(name, flags, dir_fd=output_fd)
                try:
                    children = tuple(sorted(os.listdir(child)))
                    if len(children) != 1:
                        _fail("status", "partial_or_foreign_namespace", children)
                    record_file(child, f"{name}/{children[0]}", children[0])
                finally:
                    os.close(child)
            else:
                _fail("status", "unsafe_namespace", name)
        for index, (descriptor, expected) in enumerate(
            zip(descriptors, identities, strict=True)
        ):
            if _identity(os.fstat(descriptor)) != expected:
                _fail("status", "output_directory_changed", index)
        if _identity(root.lstat()) != identities[0]:
            _fail("status", "repository_root_substituted", root)
        return values, blobs, hashes
    except OSError as exc:
        raise QA1RunnerError(
            "status", "output_snapshot_unreadable", authority.OUTPUT_NAMESPACE
        ) from exc
    finally:
        for descriptor in reversed(descriptors):
            os.close(descriptor)


def _validate_replay_receipt(value: object) -> dict[str, object]:
    expected = expected_replay_receipt()
    if not isinstance(value, Mapping) or set(value) != set(expected):
        _fail("status", "replay_receipt_schema", (sorted(expected), value))
    receipt = {key: value[key] for key in expected}
    if receipt != expected:
        _fail("status", "replay_receipt_identity", receipt)
    return receipt


def _validate_exact_summary(value: object) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        _fail("status", "exact_schema", type(value).__name__)
    channels = value.get("channels")
    if not isinstance(channels, list) or len(channels) != authority.CHANNEL_COUNT:
        _fail("status", "exact_channel_count", type(channels).__name__)
    names = tuple(
        item.get("channel") for item in channels if isinstance(item, Mapping)
    )
    if names != authority.CHANNEL_ORDER:
        _fail("status", "exact_channel_order", names)
    failed = [
        item["channel"]
        for item in channels
        if isinstance(item, Mapping) and item.get("admission_passed") is not True
    ]
    if (
        value.get("channel_order") != list(authority.CHANNEL_ORDER)
        or value.get("channel_count") != authority.CHANNEL_COUNT
        or value.get("admission_is_all_of") is not True
        or value.get("interval_owner") != authority.INTERVAL_OWNER
        or value.get("independent_route_required") is not True
        or list(value.get("failed_channels", [])) != failed
        or value.get("complete_admission_passed") is not (not failed)
    ):
        _fail("status", "exact_all_of_identity", value.get("failed_channels"))
    for item in channels:
        if not isinstance(item, Mapping):
            _fail("status", "exact_channel_type", type(item).__name__)
        _fraction(item.get("d12_upper"), label=f"{item.get('channel')}.d12_upper")
    return value


def _validate_terminal(
    value: object, *, authority_commit: str
) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        _fail("status", "terminal_schema", "terminal is not a mapping")
    classification = value.get("classification")
    if classification not in PUBLISHABLE_TERMINAL_CLASSES:
        _fail("status", "terminal_class_not_publishable", classification)
    sealed = (
        authority.SEALED_STORE_LEAF_COUNT,
        authority.SEALED_STORE_SNAPSHOT_SHA256,
    )
    expected_base = _terminal_base(
        authority_commit, str(classification), sealed, sealed
    )
    if any(value.get(name) != expected for name, expected in expected_base.items()):
        _fail("status", "terminal_identity", classification)
    base_keys = set(expected_base)
    if classification == PREMISE_STOP_CLASS:
        _require_keys(
            value,
            base_keys | {"typed_stop", "replay_receipt"},
            label="premise stop terminal",
        )
        stop = _require_keys(
            value["typed_stop"], {"type", "detail"}, label="premise stop"
        )
        if (
            not isinstance(stop["type"], str)
            or not stop["type"]
            or not isinstance(stop["detail"], str)
            or not stop["detail"]
            or len(str(stop["detail"])) > 640
        ):
            _fail("status", "premise_stop_terminal", stop)
        _validate_replay_receipt(value["replay_receipt"])
        return value
    if classification == INCONCLUSIVE_CLASS:
        _require_keys(
            value,
            base_keys
            | {
                "typed_inconclusive",
                "replay_receipt",
                "shadow_path_count",
                "shadow_proposal_count",
                "SSPRK3_stage_and_endpoint_record_count",
            },
            label="inconclusive terminal",
        )
        closed = _require_keys(
            value["typed_inconclusive"],
            {"type", "reason", "detail"},
            label="typed inconclusive",
        )
        if (
            closed["type"]
            not in {
                "ExactCompleteCResourceExhausted",
                "ExactCompleteCRouteDisagreement",
                "TDG10ExactCompleteCRuntimeClosed",
            }
            or not isinstance(closed["reason"], str)
            or not closed["reason"]
        ):
            _fail("status", "inconclusive_terminal", closed)
        _validate_replay_receipt(value["replay_receipt"])
        if (
            value.get("shadow_path_count") != 7
            or value.get("shadow_proposal_count") != 7
            or value.get("SSPRK3_stage_and_endpoint_record_count") != 28
        ):
            _fail("status", "inconclusive_shadow_budget", classification)
        return value
    extra = {
        "replay_receipt",
        "shadow_path_count",
        "shadow_proposal_count",
        "SSPRK3_stage_and_endpoint_record_count",
        "exact_complete_C",
    }
    if classification == PASS_CLASS:
        extra.add("diagnostic_fine_endpoint")
    _require_keys(value, base_keys | extra, label="completed terminal")
    _validate_replay_receipt(value["replay_receipt"])
    exact = _validate_exact_summary(value["exact_complete_C"])
    expected = reduce_qa1_terminal(
        complete_admission_passed=bool(exact["complete_admission_passed"]),
        failed_channels=tuple(exact["failed_channels"]),
    )
    if (
        classification != expected
        or value.get("shadow_path_count") != 7
        or value.get("shadow_proposal_count") != 7
        or value.get("SSPRK3_stage_and_endpoint_record_count") != 28
    ):
        _fail("status", "completed_terminal_counts", classification)
    if classification == PASS_CLASS:
        diagnostic = _require_keys(
            value["diagnostic_fine_endpoint"],
            {
                "present",
                "accepted_state",
                "payload_relative",
                "payload_semantic_sha256",
                "payload_raw_sha256",
                "payload_manifest",
                "descriptor_relative",
                "descriptor_sha256",
            },
            label="diagnostic fine endpoint",
        )
        if diagnostic.get("present") is not True:
            _fail("status", "diagnostic_absent", classification)
        if diagnostic.get("accepted_state") is not False:
            _fail("status", "diagnostic_accepted_state", classification)
        _hex_digest(
            diagnostic.get("payload_semantic_sha256"),
            label="payload_semantic_sha256",
        )
        _hex_digest(
            diagnostic.get("payload_raw_sha256"), label="payload_raw_sha256"
        )
        _hex_digest(
            diagnostic.get("descriptor_sha256"), label="descriptor_sha256"
        )
    return value


def _descriptor_address(descriptor: Mapping[str, object]) -> str:
    body = {
        key: value
        for key, value in descriptor.items()
        if key != "descriptor_sha256"
    }
    return sha256(_canonical(body)).hexdigest()


def _inspect_output(root: Path, *, authority_commit: str) -> dict[str, object]:
    snapshot = _snapshot_output(root)
    if snapshot is None:
        return {"state": "absent"}
    values, blobs, hashes = snapshot
    terminal = _validate_terminal(
        values["terminal.json"], authority_commit=authority_commit
    )
    classification = str(terminal["classification"])
    expected_leaves = ["manifest.json", "terminal.json"]
    if classification == PASS_CLASS:
        diagnostic = terminal["diagnostic_fine_endpoint"]
        expected_leaves.extend(
            sorted(
                (
                    str(diagnostic["payload_relative"]),
                    str(diagnostic["descriptor_relative"]),
                )
            )
        )
        payload_relative = str(diagnostic["payload_relative"])
        descriptor_relative = str(diagnostic["descriptor_relative"])
        if payload_relative not in blobs or descriptor_relative not in blobs:
            _fail("status", "diagnostic_files_absent", classification)
        if hashes[payload_relative] != diagnostic["payload_raw_sha256"]:
            _fail("status", "payload_raw_hash", payload_relative)
        if hashes[descriptor_relative] != sha256(blobs[descriptor_relative]).hexdigest():
            _fail("status", "descriptor_file_hash", descriptor_relative)
        descriptor = values[descriptor_relative]
        if not isinstance(descriptor, Mapping):
            _fail("status", "descriptor_type", type(descriptor).__name__)
        recomputed = _descriptor_address(descriptor)
        stored_address = _hex_digest(
            descriptor.get("descriptor_sha256"), label="descriptor.descriptor_sha256"
        )
        terminal_address = _hex_digest(
            diagnostic.get("descriptor_sha256"), label="terminal.descriptor_sha256"
        )
        if (
            recomputed != stored_address
            or recomputed != terminal_address
            or Path(descriptor_relative).name != f"{recomputed}.json"
        ):
            _fail("status", "descriptor_address", descriptor_relative)
        try:
            _arrays, semantic, manifest = decode_payload(
                blobs[payload_relative],
                expected_semantic_sha256=str(diagnostic["payload_semantic_sha256"]),
                expected_manifest=diagnostic.get("payload_manifest"),
            )
        except (HLT16StateStoreError, TypeError, ValueError) as exc:
            raise QA1RunnerError("status", "payload_decode", exc) from exc
        if (
            semantic != diagnostic["payload_semantic_sha256"]
            or semantic != descriptor.get("payload_semantic_sha256")
            or hashes[payload_relative] != diagnostic["payload_raw_sha256"]
            or hashes[payload_relative] != descriptor.get("payload_raw_sha256")
            or list(manifest) != list(diagnostic.get("payload_manifest") or ())
            or list(manifest) != list(descriptor.get("payload_manifest") or ())
            or Path(payload_relative).name != f"{semantic}.npz"
        ):
            _fail("status", "payload_identity", payload_relative)
        if (
            descriptor.get("accepted_state") is not False
            or descriptor.get("campaign_descriptor") is not False
            or descriptor.get("proto15_cursor") is not False
        ):
            _fail("status", "descriptor_scope", descriptor_relative)
    elif any(name.startswith("payloads/") or name.startswith("fine-endpoint/") for name in blobs):
        _fail("status", "diagnostic_without_all_pass", classification)
    if values["manifest.json"] != _manifest(
        authority_commit, output_leaves=expected_leaves
    ):
        _fail("status", "manifest_identity", values["manifest.json"])
    return {
        "state": "terminal",
        "hashes": hashes,
        "classification": classification,
    }


def status(root: Path, *, authority_commit: str) -> dict[str, object]:
    _execution_authority(root.resolve(), authority_commit)
    before = _snapshot_store(root.resolve())
    if before != (
        authority.SEALED_STORE_LEAF_COUNT,
        authority.SEALED_STORE_SNAPSHOT_SHA256,
    ):
        _fail("status", "sealed_store_identity", before)
    output = _inspect_output(root.resolve(), authority_commit=authority_commit)
    return {
        "artifact_id": authority.ARTIFACT_ID,
        "authority_commit": authority_commit,
        "output": output,
        "store": {"leaf_count": before[0], "sha256": before[1]},
        "safe_to_run": output["state"] == "absent",
        "state_advance_authorized": False,
        "status_read_only": True,
    }


def _shadow_budget() -> dict[str, int]:
    return {
        "shadow_path_count": 7,
        "shadow_proposal_count": 7,
        "SSPRK3_stage_and_endpoint_record_count": 28,
    }


def run(root: Path, *, authority_commit: str) -> dict[str, object]:
    repository = root.resolve()
    receipt = _execution_authority(repository, authority_commit)
    if (
        receipt.tableau_selector != authority.TABLEAU_SELECTOR
        or receipt.retry != authority.RETRY
        or receipt.predecessor_generation != authority.PREDECESSOR_GENERATION
        or receipt.state_advance_authorized
        or receipt.new_protocol_id_authorized
        or not receipt.diagnostic_qualification_only
    ):
        _fail("authority", "receipt_semantics", receipt)
    authority.require_output_absent(repository)
    before = _snapshot_store(repository)
    sealed = (
        authority.SEALED_STORE_LEAF_COUNT,
        authority.SEALED_STORE_SNAPSHOT_SHA256,
    )
    if before != sealed:
        _fail("provenance", "sealed_store_identity", (before, sealed))
    replay_receipt: dict[str, object] = {}
    try:
        restored = _restore_replay(repository)
        replay_receipt = _replay_receipt(restored)
        prepared = _prepare_shadow(restored)
    except tdg6.TDG6RefinementPathStop as exc:
        after = _snapshot_store(repository)
        if after != before:
            _fail("provenance", "campaign_store_mutated", (before, after))
        terminal = {
            **_terminal_base(authority_commit, PREMISE_STOP_CLASS, before, after),
            "typed_stop": {
                "type": type(exc).__name__,
                "detail": " ".join(str(exc).split())[:640],
            },
            "replay_receipt": replay_receipt,
        }
        _publish_terminal(
            repository,
            _manifest(authority_commit),
            terminal,
            store_before=before,
            store_after=after,
        )
        return terminal

    try:
        assessment = _assess_prepared_exact_complete_c(prepared)
    except INCONCLUSIVE_EXCEPTIONS as exc:
        after = _snapshot_store(repository)
        if after != before:
            _fail("provenance", "campaign_store_mutated", (before, after))
        evidence = getattr(exc, "evidence", None)
        reason = getattr(evidence, "reason", type(exc).__name__)
        terminal = {
            **_terminal_base(authority_commit, INCONCLUSIVE_CLASS, before, after),
            "typed_inconclusive": {
                "type": type(exc).__name__,
                "reason": str(reason),
                "detail": " ".join(str(exc).split())[:640],
            },
            "replay_receipt": replay_receipt,
            **_shadow_budget(),
        }
        _publish_terminal(
            repository,
            _manifest(authority_commit),
            terminal,
            store_before=before,
            store_after=after,
        )
        return terminal

    serialized = serialize_exact_assessment(assessment)
    classification = reduce_qa1_terminal(
        complete_admission_passed=bool(serialized["complete_admission_passed"]),
        failed_channels=tuple(serialized["failed_channels"]),
    )
    extras: dict[str, bytes] | None = None
    diagnostic = None
    if classification == PASS_CLASS:
        extras, diagnostic = serialize_diagnostic_endpoint(
            restored, prepared, serialized
        )
    after = _snapshot_store(repository)
    if after != before:
        _fail("provenance", "campaign_store_mutated", (before, after))
    after_fingerprint = ti2_runner._member_fingerprint(
        restored.member, str(restored.fingerprint["descriptor_sha256"])
    )
    if after_fingerprint != restored.fingerprint:
        _fail("provenance", "member_mutated", after_fingerprint)
    terminal = {
        **_terminal_base(authority_commit, classification, before, after),
        "replay_receipt": replay_receipt,
        **_shadow_budget(),
        "exact_complete_C": serialized,
    }
    if diagnostic is not None:
        terminal["diagnostic_fine_endpoint"] = diagnostic
    leaves = ["manifest.json", "terminal.json"]
    if extras:
        leaves.extend(sorted(extras))
    _publish_terminal(
        repository,
        _manifest(authority_commit, output_leaves=leaves),
        terminal,
        store_before=before,
        store_after=after,
        extra_files=extras,
    )
    return terminal


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--authority-commit", required=True)
    parser.add_argument("--status", action="store_true")
    arguments = parser.parse_args(argv)
    try:
        result = (
            status(ROOT, authority_commit=arguments.authority_commit)
            if arguments.status
            else run(ROOT, authority_commit=arguments.authority_commit)
        )
    except QA1RunnerError as exc:
        print(
            json.dumps(
                {
                    "artifact_id": authority.ARTIFACT_ID,
                    "classification": INVALID_CLASS,
                    "owner": exc.owner,
                    "code": exc.code,
                    "detail": exc.detail,
                },
                sort_keys=True,
                separators=(",", ":"),
            )
        )
        return 2
    print(json.dumps(result, sort_keys=True, separators=(",", ":"), default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
