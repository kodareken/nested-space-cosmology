"""Focused adversarial controls for the immutable TDG11-IMP1 ledger codec."""

from __future__ import annotations

import ast
from dataclasses import FrozenInstanceError
from fractions import Fraction
from hashlib import sha256
import json
import math
from pathlib import Path
import sys
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src")]

from recursive_horizons.evidence_io import (  # noqa: E402
    canonical_json_bytes,
    load_canonical_json,
)
from recursive_horizons.fgc.evolution.numerical_engine import (  # noqa: E402
    COMPARATOR_METHOD,
    PRIMARY_METHOD,
)
from recursive_horizons.fgc.evolution.tdg6_temporal_admission_design import (  # noqa: E402
    TDG6_COMPLETE_STATE_CHANNELS,
)
from recursive_horizons.fgc.evolution.tdg6_temporal_admission_runtime import (  # noqa: E402
    TDG6TemporalLedger,
)
from recursive_horizons.fgc.evolution.tdg11_imp1_ledger import (  # noqa: E402
    IMP1_ACCEPTED_SUBSTEP_COUNT,
    IMP1_CHANNEL_COUNT,
    IMP1_CHECKPOINT_KEYS,
    IMP1_COMPLETE_STATE_CHANNELS,
    IMP1_LEDGER_SCHEMA,
    IMP1_METHOD_STAGE_RECORD_COUNT,
    IMP1_REJECTION_EVENT_TYPE,
    IMP1_REJECTION_KEYS,
    IMP1_ZERO_DEBIT,
    MAX_RATIONAL_BITS,
    TDG11IMP1Error,
    TDG11IMP1Ledger,
    TDG11IMP1ResourceError,
    TDG11IMP1TemporalRejection,
    accept_imp1_step,
    append_imp1_rejection,
    combined_diagnostic_debit_vector,
    encode_inherited_snapshot,
    imp1_checkpoint_extension,
    inherited_tdg6_ledger,
    integer_text,
    parse_integer_text,
    restore_imp1_checkpoint_extension,
    seed_imp1_ledger,
)
from recursive_horizons.fgc.evolution import tdg11_imp1_ledger as ledger_mod  # noqa: E402


MODULE_PATH = ROOT / "src/recursive_horizons/fgc/evolution/tdg11_imp1_ledger.py"
Q = Fraction
ORIGIN_TIME = 23 / 16
START_INDEX = 10
START_SERIAL = 100
FORBIDDEN_IMPORT_MARKERS = (
    "binder",
    "scripts",
    "tdg11_msel1_pref1",
    "publish_exclusive",
    "git_read",
    "hlt16",
    "campaign",
)


def _sha(tag: str) -> str:
    return sha256(tag.encode("ascii")).hexdigest()


def _tdg6_rejection() -> str:
    return json.dumps(
        {"event_type": "rejected_TDG6_temporal_admission", "synthetic": True},
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )


def _inherited(*, active_retry: bool = False) -> TDG6TemporalLedger:
    current = 1 if active_retry else 0
    last_accepted = 0 if active_retry else 1
    return TDG6TemporalLedger(
        last_accepted_time=ORIGIN_TIME,
        accepted_macro_step_count=2,
        cumulative_temporal_retry_count=1,
        current_macro_step_temporal_retry_count=current,
        last_accepted_macro_step_temporal_retry_count=last_accepted,
        accumulated_debit_vector=(0.5, 0.25) + (0.0,) * 16,
        serialized_temporal_rejections=(_tdg6_rejection(),),
    )


def _seed(*, method: str = PRIMARY_METHOD, inherited: TDG6TemporalLedger | None = None):
    return seed_imp1_ledger(
        _inherited() if inherited is None else inherited,
        method=method,
        state_sha256=_sha("state"),
        step_index=START_INDEX,
        transaction_serial=START_SERIAL,
        origin_receipt_sha256=_sha("receipt"),
    )


def _zeros() -> tuple[Fraction, ...]:
    return IMP1_ZERO_DEBIT


def _debit(*pairs: tuple[int, Fraction]) -> tuple[Fraction, ...]:
    values = list(_zeros())
    for index, value in pairs:
        values[index] = value
    return tuple(values)


def _half(width: float) -> float:
    return math.ldexp(width, -1)


def _rejection(
    ledger: TDG11IMP1Ledger,
    *,
    current: int | None = None,
    cumulative: int | None = None,
    method: str | None = None,
    event: str | None = None,
    serial: int | None = None,
    channels: list[str] | None = None,
    extra: dict[str, object] | None = None,
    omit: str | None = None,
    width: float | None = None,
    next_cap: float | None = None,
    time_hex: str | None = None,
    state: str | None = None,
    step_index: int | None = None,
    event_target: float = 32.0,
    preparation: str | None = None,
) -> dict[str, object]:
    if width is None:
        if ledger.current_macro_step_temporal_retry_count > 0:
            previous = json.loads(ledger.serialized_temporal_rejections[-1])
            width = float.fromhex(str(previous["next_cap_hex"]))
        else:
            width = 0.125
    if next_cap is None:
        next_cap = _half(width)
    mapping: dict[str, object] = {
        "event_type": IMP1_REJECTION_EVENT_TYPE if event is None else event,
        "method": ledger.method if method is None else method,
        "time_hex": ledger.last_accepted_time.hex() if time_hex is None else time_hex,
        "event_target_hex": event_target.hex(),
        "attempted_width_hex": width.hex(),
        "next_cap_hex": next_cap.hex(),
        "accepted_state_sha256": ledger.current_state_sha256
        if state is None
        else state,
        "accepted_step_index": (
            ledger.current_step_index if step_index is None else step_index
        ),
        "accepted_transaction_serial": (
            ledger.current_transaction_serial if serial is None else serial
        ),
        "assessment_sha256": _sha("assessment"),
        "preparation_sha256": _sha("preparation")
        if preparation is None
        else preparation,
        "channel_order": list(IMP1_COMPLETE_STATE_CHANNELS),
        "failed_channels": ["u:alpha"] if channels is None else channels,
        "retry_count_for_current_macro_step": (
            ledger.current_macro_step_temporal_retry_count + 1
            if current is None
            else current
        ),
        "cumulative_temporal_retry_count": (
            ledger.cumulative_temporal_retry_count + 1
            if cumulative is None
            else cumulative
        ),
    }
    if extra:
        mapping.update(extra)
    if omit is not None:
        mapping.pop(omit)
    return mapping


def _recompute_checkpoint_digests(mapping: dict[str, object]) -> dict[str, object]:
    updated = dict(mapping)
    snapshot = updated["inherited_snapshot"]
    updated["inherited_snapshot_sha256"] = sha256(
        str(snapshot).encode("ascii")
    ).hexdigest()
    updated["accumulated_debit_sha256"] = sha256(
        canonical_json_bytes(updated["accumulated_debit_vector"])
    ).hexdigest()
    updated["serialized_temporal_rejections_sha256"] = sha256(
        canonical_json_bytes(updated["serialized_temporal_rejections"])
    ).hexdigest()
    return updated


def _rewrite_rejection(text: str, **changes: object) -> str:
    payload = json.loads(text)
    payload.update(changes)
    return TDG11IMP1TemporalRejection.from_mapping(payload).canonical_text()


def _history_ledger() -> TDG11IMP1Ledger:
    seeded = _seed()
    first = append_imp1_rejection(seeded, _rejection(seeded))
    second = append_imp1_rejection(first, _rejection(first))
    accepted = _accept(second, debit=_debit((0, Q(1, 6))))
    gapped = _accept(accepted, final_time=1.625, state=_sha("gap-state"))
    return append_imp1_rejection(gapped, _rejection(gapped, preparation=_sha("prep-2")))


def _accept(
    ledger: TDG11IMP1Ledger,
    *,
    final_time: float = 1.5,
    debit: tuple[Fraction, ...] | None = None,
    step_index: int | None = None,
    transaction_serial: int | None = None,
    state: str | None = None,
) -> TDG11IMP1Ledger:
    return accept_imp1_step(
        ledger,
        final_time=final_time,
        state_sha256=_sha("next-state") if state is None else state,
        step_index=(
            ledger.current_step_index + IMP1_ACCEPTED_SUBSTEP_COUNT
            if step_index is None
            else step_index
        ),
        transaction_serial=(
            ledger.current_transaction_serial
            + IMP1_METHOD_STAGE_RECORD_COUNT[ledger.method]
            if transaction_serial is None
            else transaction_serial
        ),
        debit_vector=_zeros() if debit is None else debit,
        assessment_sha256=_sha("pass"),
    )


def _roundtrip(ledger: TDG11IMP1Ledger) -> TDG11IMP1Ledger:
    mapping = imp1_checkpoint_extension(ledger)
    restored = restore_imp1_checkpoint_extension(mapping)
    raw = canonical_json_bytes(mapping)
    from_bytes = restore_imp1_checkpoint_extension(load_canonical_json(raw))
    if restored != ledger or from_bytes != ledger:
        raise AssertionError("IMP1 checkpoint round-trip lost ledger identity")
    return restored


class TDG11IMP1LedgerTests(unittest.TestCase):
    def test_rounded_subnormal_half_is_not_an_exact_retry_cap(self):
        # 3 * the smallest subnormal has no representable exact half. IEEE
        # rounding produces a positive number, which is not the promised half.
        width = float.fromhex("0x0.0000000000003p-1022")
        with self.assertRaisesRegex(TDG11IMP1Error, "exact binary64 half"):
            append_imp1_rejection(_seed(), _rejection(_seed(), width=width))

    def test_exported_api_and_closed_field_contracts(self) -> None:
        public = set(ledger_mod.__all__)
        for name in (
            "TDG11IMP1Ledger",
            "TDG11IMP1TemporalRejection",
            "seed_imp1_ledger",
            "append_imp1_rejection",
            "accept_imp1_step",
            "imp1_checkpoint_extension",
            "restore_imp1_checkpoint_extension",
            "combined_diagnostic_debit_vector",
            "MAX_RATIONAL_BITS",
            "MAX_REJECTION_RECORDS",
            "MAX_CHECKPOINT_BYTES",
        ):
            self.assertIn(name, public)
            self.assertTrue(hasattr(ledger_mod, name))
        record = TDG11IMP1TemporalRejection.from_mapping(_rejection(_seed()))
        self.assertEqual(set(record.as_mapping()), set(IMP1_REJECTION_KEYS))
        self.assertEqual(
            set(imp1_checkpoint_extension(_seed())), set(IMP1_CHECKPOINT_KEYS)
        )
        self.assertEqual(IMP1_CHANNEL_COUNT, 18)
        self.assertEqual(IMP1_COMPLETE_STATE_CHANNELS, TDG6_COMPLETE_STATE_CHANNELS)
        self.assertEqual(IMP1_METHOD_STAGE_RECORD_COUNT[PRIMARY_METHOD], 20)
        self.assertEqual(IMP1_METHOD_STAGE_RECORD_COUNT[COMPARATOR_METHOD], 16)
        self.assertEqual(IMP1_ACCEPTED_SUBSTEP_COUNT, 4)
        self.assertEqual(IMP1_LEDGER_SCHEMA, "FGC-1-TDG11-IMP1-ledger-v1")
        self.assertEqual(
            IMP1_REJECTION_EVENT_TYPE, "rejected_TDG11_IMP1_temporal_admission"
        )
        self.assertIn("preparation_sha256", IMP1_REJECTION_KEYS)
        self.assertIn("event_target_hex", IMP1_REJECTION_KEYS)

    def test_seed_preserves_inherited_snapshot_and_refuses_active_retry(self) -> None:
        inherited = _inherited()
        ledger = _seed(inherited=inherited)
        frozen = inherited_tdg6_ledger(ledger)
        self.assertEqual(frozen, inherited)
        self.assertEqual(
            frozen.last_accepted_time.hex(), inherited.last_accepted_time.hex()
        )
        self.assertEqual(
            tuple(value.hex() for value in frozen.accumulated_debit_vector),
            tuple(value.hex() for value in inherited.accumulated_debit_vector),
        )
        self.assertEqual(
            frozen.serialized_temporal_rejections,
            inherited.serialized_temporal_rejections,
        )
        self.assertEqual(
            ledger.last_accepted_time.hex(), inherited.last_accepted_time.hex()
        )
        self.assertEqual(ledger.current_time, ledger.last_accepted_time)
        self.assertEqual(ledger.accepted_macro_step_count, 0)
        self.assertEqual(ledger.accumulated_debit_vector, _zeros())
        self.assertEqual(ledger.origin_state_sha256, ledger.current_state_sha256)
        self.assertEqual(ledger.origin_step_index, START_INDEX)
        self.assertEqual(ledger.origin_transaction_serial, START_SERIAL)
        self.assertIs(ledger.historical_origin_authenticated, False)
        self.assertEqual(
            encode_inherited_snapshot(inherited), ledger.inherited_snapshot
        )
        with self.assertRaises(TDG11IMP1Error):
            _seed(inherited=_inherited(active_retry=True))
        with self.assertRaises(TypeError):
            seed_imp1_ledger(
                object(),  # type: ignore[arg-type]
                method=PRIMARY_METHOD,
                state_sha256=_sha("state"),
                step_index=0,
                transaction_serial=0,
                origin_receipt_sha256=_sha("receipt"),
            )

    def test_roundtrip_after_reject_and_accept_retains_history(self) -> None:
        seeded = _seed()
        original = seeded
        rejected = append_imp1_rejection(seeded, _rejection(seeded))
        self.assertIsNot(rejected, seeded)
        self.assertEqual(seeded, original)
        self.assertEqual(rejected.current_macro_step_temporal_retry_count, 1)
        self.assertEqual(rejected.cumulative_temporal_retry_count, 1)
        self.assertEqual(rejected.accepted_macro_step_count, 0)
        self.assertEqual(
            rejected.last_accepted_time.hex(), seeded.last_accepted_time.hex()
        )
        self.assertEqual(rejected.accumulated_debit_vector, _zeros())
        self.assertEqual(len(rejected.serialized_temporal_rejections), 1)
        restored_reject = _roundtrip(rejected)
        self.assertEqual(restored_reject, rejected)
        accepted = _accept(
            restored_reject,
            debit=_debit((0, Q(1, 6)), (1, Q(1, 3))),
        )
        self.assertEqual(accepted.accepted_macro_step_count, 1)
        self.assertEqual(accepted.current_macro_step_temporal_retry_count, 0)
        self.assertEqual(accepted.last_accepted_macro_step_temporal_retry_count, 1)
        self.assertEqual(accepted.cumulative_temporal_retry_count, 1)
        self.assertEqual(len(accepted.serialized_temporal_rejections), 1)
        self.assertEqual(accepted.current_step_index, START_INDEX + 4)
        self.assertEqual(accepted.current_transaction_serial, START_SERIAL + 20)
        self.assertEqual(accepted.accumulated_debit_vector[0], Q(1, 6))
        self.assertEqual(accepted.accumulated_debit_vector[1], Q(1, 3))
        restored_accept = _roundtrip(accepted)
        self.assertEqual(restored_accept, accepted)
        self.assertEqual(restored_accept.inherited_snapshot, seeded.inherited_snapshot)

    def test_nonzero_inherited_debit_history_and_combined_diagnostic(self) -> None:
        inherited = _inherited()
        ledger = _accept(
            append_imp1_rejection(
                _seed(inherited=inherited), _rejection(_seed(inherited=inherited))
            ),
            debit=_debit((0, Q(1, 3))),
        )
        frozen = inherited_tdg6_ledger(ledger)
        self.assertEqual(frozen.accumulated_debit_vector[0].hex(), (0.5).hex())
        self.assertEqual(frozen.accumulated_debit_vector[1].hex(), (0.25).hex())
        self.assertEqual(
            frozen.serialized_temporal_rejections,
            inherited.serialized_temporal_rejections,
        )
        self.assertEqual(frozen.last_accepted_time, ORIGIN_TIME)
        self.assertNotEqual(ledger.last_accepted_time, frozen.last_accepted_time)
        combined = combined_diagnostic_debit_vector(ledger)
        self.assertEqual(combined[0], Q(1, 2) + Q(1, 3))
        self.assertEqual(combined[1], Q(1, 4))
        self.assertEqual(combined[2], Q(0))

    def test_ssp_method_requires_sixteen_stage_records(self) -> None:
        ledger = _seed(method=COMPARATOR_METHOD)
        accepted = _accept(ledger)
        self.assertEqual(accepted.current_transaction_serial, START_SERIAL + 16)
        self.assertEqual(accepted.current_step_index, START_INDEX + 4)
        with self.assertRaises(TDG11IMP1Error):
            _accept(ledger, transaction_serial=START_SERIAL + 20)

    def test_bool_and_int_aliases_rejected(self) -> None:
        ledger = _seed()
        with self.assertRaises(TypeError):
            _seed().__class__(
                **{
                    **{field: getattr(ledger, field) for field in ledger.__slots__},
                    "current_step_index": True,
                }
            )
        with self.assertRaises(TypeError):
            _accept(ledger, step_index=True)  # type: ignore[arg-type]
        with self.assertRaises(TypeError):
            accept_imp1_step(
                ledger,
                final_time=1.5,
                state_sha256=_sha("next-state"),
                step_index=START_INDEX + 4,
                transaction_serial=START_SERIAL + 20,
                debit_vector=(1,) + (Q(0),) * 17,  # type: ignore[arg-type]
                assessment_sha256=_sha("pass"),
            )
        mapping = imp1_checkpoint_extension(ledger)
        mapping["historical_origin_authenticated"] = 0
        with self.assertRaises(TDG11IMP1Error):
            restore_imp1_checkpoint_extension(mapping)
        mapping = imp1_checkpoint_extension(ledger)
        mapping["complete_rejection_evidence_present"] = 1
        with self.assertRaises(TDG11IMP1Error):
            restore_imp1_checkpoint_extension(mapping)
        with self.assertRaises(TDG11IMP1Error):
            parse_integer_text(True)  # type: ignore[arg-type]
        with self.assertRaises(TDG11IMP1Error):
            parse_integer_text(1)  # type: ignore[arg-type]

    def test_negative_noncanonical_and_oversized_rationals(self) -> None:
        ledger = _seed()
        with self.assertRaises(TDG11IMP1Error):
            _accept(ledger, debit=_debit((0, Q(-1, 3))))
        with self.assertRaises(TDG11IMP1Error):
            _accept(ledger, debit=_debit((0, Q(1, 9))))
        with self.assertRaises(TDG11IMP1Error):
            _accept(ledger, debit=_debit((0, Q(1, 5))))
        with self.assertRaises(TDG11IMP1ResourceError):
            _accept(ledger, debit=_debit((0, Q(2**32768, 1))))
        with self.assertRaises(TDG11IMP1ResourceError):
            _accept(ledger, debit=_debit((0, Q(1, 2**32768))))
        self.assertEqual(ledger.accumulated_debit_vector, _zeros())
        mapping = imp1_checkpoint_extension(_accept(ledger, debit=_debit((0, Q(1, 6)))))
        mapping["accumulated_debit_vector"][0] = {"numerator": "2", "denominator": "4"}
        mapping["accumulated_debit_sha256"] = sha256(
            canonical_json_bytes(mapping["accumulated_debit_vector"])
        ).hexdigest()
        with self.assertRaises(TDG11IMP1Error):
            restore_imp1_checkpoint_extension(mapping)
        mapping = imp1_checkpoint_extension(ledger)
        mapping["accumulated_debit_vector"][0] = {"numerator": "-1", "denominator": "2"}
        mapping["accumulated_debit_sha256"] = sha256(
            canonical_json_bytes(mapping["accumulated_debit_vector"])
        ).hexdigest()
        with self.assertRaises(TDG11IMP1Error):
            restore_imp1_checkpoint_extension(mapping)
        huge = _accept(ledger, debit=_debit((0, Q(2 ** (MAX_RATIONAL_BITS - 1), 1))))
        restored = _roundtrip(huge)
        self.assertEqual(
            restored.accumulated_debit_vector[0], Q(2 ** (MAX_RATIONAL_BITS - 1), 1)
        )

    def test_wrong_channel_order_method_event_and_serial(self) -> None:
        ledger = _seed()
        with self.assertRaises(TDG11IMP1Error):
            append_imp1_rejection(
                ledger,
                _rejection(ledger, event="rejected_TDG6_temporal_admission"),
            )
        with self.assertRaises(TDG11IMP1Error):
            append_imp1_rejection(ledger, _rejection(ledger, method=COMPARATOR_METHOD))
        with self.assertRaises(TDG11IMP1Error):
            append_imp1_rejection(ledger, _rejection(ledger, serial=START_SERIAL + 1))
        with self.assertRaises(TDG11IMP1Error):
            append_imp1_rejection(
                ledger,
                _rejection(ledger, channels=["q:chi", "u:alpha"]),
            )
        mapping = imp1_checkpoint_extension(ledger)
        mapping["channel_order"] = list(reversed(IMP1_COMPLETE_STATE_CHANNELS))
        with self.assertRaises(TDG11IMP1Error):
            restore_imp1_checkpoint_extension(mapping)
        with self.assertRaises(TDG11IMP1Error):
            _accept(ledger, step_index=START_INDEX + 3)
        with self.assertRaises(TDG11IMP1Error):
            _accept(ledger, transaction_serial=START_SERIAL + 19)
        with self.assertRaises(TDG11IMP1Error):
            _accept(ledger, final_time=ORIGIN_TIME)
        with self.assertRaises(TDG11IMP1Error):
            _accept(ledger, final_time=ORIGIN_TIME - 0.125)

    def test_missing_and_extra_keys(self) -> None:
        ledger = _seed()
        with self.assertRaises(TDG11IMP1Error):
            append_imp1_rejection(ledger, _rejection(ledger, omit="failed_channels"))
        with self.assertRaises(TDG11IMP1Error):
            append_imp1_rejection(ledger, _rejection(ledger, omit="preparation_sha256"))
        with self.assertRaises(TDG11IMP1Error):
            append_imp1_rejection(ledger, _rejection(ledger, omit="event_target_hex"))
        with self.assertRaises(TDG11IMP1Error):
            append_imp1_rejection(
                ledger, _rejection(ledger, extra={"physical_classification": False})
            )
        mapping = imp1_checkpoint_extension(ledger)
        missing = dict(mapping)
        missing.pop("channel_order")
        with self.assertRaises(TDG11IMP1Error):
            restore_imp1_checkpoint_extension(missing)
        extra = dict(mapping)
        extra["physical_result"] = True
        with self.assertRaises(TDG11IMP1Error):
            restore_imp1_checkpoint_extension(extra)
        with self.assertRaises(TDG11IMP1Error):
            TDG11IMP1TemporalRejection.from_mapping(_rejection(ledger, omit="time_hex"))

    def test_changed_hashes_fail_restore(self) -> None:
        ledger = _accept(append_imp1_rejection(_seed(), _rejection(_seed())))
        mapping = imp1_checkpoint_extension(ledger)
        mutated = dict(mapping)
        mutated["accumulated_debit_sha256"] = _sha("tampered-debit")
        with self.assertRaises(TDG11IMP1Error):
            restore_imp1_checkpoint_extension(mutated)
        mutated = dict(mapping)
        mutated["serialized_temporal_rejections_sha256"] = _sha("tampered-history")
        with self.assertRaises(TDG11IMP1Error):
            restore_imp1_checkpoint_extension(mutated)
        mutated = dict(mapping)
        mutated["inherited_snapshot_sha256"] = _sha("tampered-snapshot")
        with self.assertRaises(TDG11IMP1Error):
            restore_imp1_checkpoint_extension(mutated)
        mutated = dict(mapping)
        debit = list(mutated["accumulated_debit_vector"])
        debit[0] = {"numerator": "1", "denominator": "2"}
        mutated["accumulated_debit_vector"] = debit
        with self.assertRaises(TDG11IMP1Error):
            restore_imp1_checkpoint_extension(mutated)

    def test_no_implicit_tdg6_time_advancement(self) -> None:
        inherited = _inherited()
        seeded = _seed(inherited=inherited)
        rejected = append_imp1_rejection(seeded, _rejection(seeded))
        accepted = _accept(rejected, final_time=1.5)
        for ledger in (seeded, rejected, accepted):
            frozen = inherited_tdg6_ledger(ledger)
            self.assertEqual(frozen.last_accepted_time.hex(), ORIGIN_TIME.hex())
            self.assertEqual(frozen.accepted_macro_step_count, 2)
            self.assertEqual(ledger.inherited_snapshot, seeded.inherited_snapshot)
        self.assertEqual(accepted.last_accepted_time.hex(), (1.5).hex())
        mapping = imp1_checkpoint_extension(accepted)
        mapping["inherited_last_accepted_time_hex"] = (1.5).hex()
        with self.assertRaises(TDG11IMP1Error):
            restore_imp1_checkpoint_extension(mapping)

    def test_immutable_failures(self) -> None:
        ledger = _seed()
        record = TDG11IMP1TemporalRejection.from_mapping(_rejection(ledger))
        with self.assertRaises(FrozenInstanceError):
            ledger.accepted_macro_step_count = 1  # type: ignore[misc]
        with self.assertRaises(FrozenInstanceError):
            record.event_type = "rejected_TDG6_temporal_admission"  # type: ignore[misc]
        rejected = append_imp1_rejection(ledger, record)
        self.assertEqual(ledger.cumulative_temporal_retry_count, 0)
        self.assertEqual(rejected.cumulative_temporal_retry_count, 1)
        accepted = _accept(rejected)
        self.assertEqual(rejected.last_accepted_time.hex(), ORIGIN_TIME.hex())
        self.assertEqual(accepted.last_accepted_time.hex(), (1.5).hex())

    def test_unsafe_state_promotion_data_fail_closed(self) -> None:
        ledger = _seed()
        mapping = imp1_checkpoint_extension(ledger)
        mapping["historical_origin_authenticated"] = True
        with self.assertRaises(TDG11IMP1Error):
            restore_imp1_checkpoint_extension(mapping)
        mapping = imp1_checkpoint_extension(ledger)
        mapping["physical_classification"] = True
        with self.assertRaises(TDG11IMP1Error):
            restore_imp1_checkpoint_extension(mapping)
        mapping = imp1_checkpoint_extension(ledger)
        mapping["complete_rejection_evidence_present"] = False
        with self.assertRaises(TDG11IMP1Error):
            restore_imp1_checkpoint_extension(mapping)
        mapping = imp1_checkpoint_extension(ledger)
        mapping["source_and_CFL_retry_counters_owned_elsewhere"] = False
        with self.assertRaises(TDG11IMP1Error):
            restore_imp1_checkpoint_extension(mapping)
        snapshot = json.loads(ledger.inherited_snapshot)
        snapshot["ledger_class"] = "TDG11IMP1Ledger"
        mapping = imp1_checkpoint_extension(ledger)
        mapping["inherited_snapshot"] = canonical_json_bytes(snapshot).decode("ascii")
        mapping["inherited_snapshot_sha256"] = sha256(
            mapping["inherited_snapshot"].encode("ascii")
        ).hexdigest()
        with self.assertRaises(TDG11IMP1Error):
            restore_imp1_checkpoint_extension(mapping)

    def test_resource_ceilings_fail_before_successor(self) -> None:
        ledger = _seed()
        with patch.object(ledger_mod, "MAX_REJECTION_RECORDS", 1):
            once = append_imp1_rejection(ledger, _rejection(ledger))
            with self.assertRaises(TDG11IMP1ResourceError):
                append_imp1_rejection(once, _rejection(once))
            self.assertEqual(once.cumulative_temporal_retry_count, 1)
        with patch.object(ledger_mod, "MAX_CHECKPOINT_BYTES", 64):
            with self.assertRaises(TDG11IMP1ResourceError):
                imp1_checkpoint_extension(ledger)
        self.assertEqual(ledger.serialized_temporal_rejections, ())

    def test_integer_codec_keeps_global_limit(self) -> None:
        before = sys.get_int_max_str_digits()
        self.assertEqual(parse_integer_text(integer_text(0)), 0)
        self.assertEqual(parse_integer_text(integer_text(1)), 1)
        self.assertEqual(parse_integer_text(integer_text(10**100)), 10**100)
        self.assertEqual(
            parse_integer_text(integer_text(2**32767)),
            2**32767,
        )
        self.assertEqual(
            parse_integer_text(integer_text(2**32768 - 1)),
            2**32768 - 1,
        )
        with self.assertRaises(TDG11IMP1ResourceError):
            integer_text(2**32768)
        for bad in ("-0", "+1", "01", "1.0", "", "١"):
            with self.subTest(bad=bad), self.assertRaises(TDG11IMP1Error):
                parse_integer_text(bad)
        self.assertEqual(sys.get_int_max_str_digits(), before)
        _roundtrip(_seed())
        self.assertEqual(sys.get_int_max_str_digits(), before)

    def test_preparation_and_event_target_bind_and_half_cap(self) -> None:
        ledger = _seed()
        record = TDG11IMP1TemporalRejection.from_mapping(_rejection(ledger))
        self.assertEqual(record.preparation_sha256, _sha("preparation"))
        self.assertEqual(record.event_target_hex, (32.0).hex())
        self.assertEqual(record.next_cap_hex, _half(0.125).hex())
        with self.assertRaises(TDG11IMP1Error):
            TDG11IMP1TemporalRejection.from_mapping(_rejection(ledger, next_cap=0.125))
        with self.assertRaises(TDG11IMP1Error):
            TDG11IMP1TemporalRejection.from_mapping(
                _rejection(ledger, event_target=ORIGIN_TIME)
            )
        with self.assertRaises(TDG11IMP1Error):
            TDG11IMP1TemporalRejection.from_mapping(
                _rejection(ledger, event_target=ORIGIN_TIME + 0.0625)
            )
        with self.assertRaises(TDG11IMP1Error):
            TDG11IMP1TemporalRejection.from_mapping(
                _rejection(ledger, preparation="not-a-digest")
            )
        with self.assertRaises(TDG11IMP1Error):
            TDG11IMP1TemporalRejection.from_mapping(
                _rejection(ledger, preparation=_sha("preparation").upper())
            )
        first = append_imp1_rejection(ledger, record)
        with self.assertRaises(TDG11IMP1Error):
            append_imp1_rejection(
                first, _rejection(first, event_target=16.0, width=0.0625)
            )
        with self.assertRaises(TDG11IMP1Error):
            append_imp1_rejection(first, _rejection(first, width=0.125))
        second = append_imp1_rejection(first, _rejection(first))
        self.assertEqual(second.current_macro_step_temporal_retry_count, 2)
        self.assertEqual(
            json.loads(second.serialized_temporal_rejections[1])["attempted_width_hex"],
            _half(0.125).hex(),
        )

    def test_semantic_history_rejects_earlier_mutations_with_recomputed_hashes(
        self,
    ) -> None:
        ledger = _history_ledger()
        self.assertEqual(ledger.accepted_macro_step_count, 2)
        self.assertEqual(ledger.current_macro_step_temporal_retry_count, 1)
        self.assertEqual(ledger.last_accepted_macro_step_temporal_retry_count, 0)
        self.assertEqual(ledger.cumulative_temporal_retry_count, 3)
        self.assertEqual(ledger.current_step_index, START_INDEX + 8)
        mapping = imp1_checkpoint_extension(ledger)
        original = list(mapping["serialized_temporal_rejections"])
        first = original[0]
        second = original[1]

        mutated_retry = dict(mapping)
        mutated_retry["serialized_temporal_rejections"] = [
            _rewrite_rejection(first, cumulative_temporal_retry_count=2),
            second,
            original[2],
        ]
        with self.assertRaises(TDG11IMP1Error):
            restore_imp1_checkpoint_extension(
                _recompute_checkpoint_digests(mutated_retry)
            )
        mutated_local = dict(mapping)
        mutated_local["serialized_temporal_rejections"] = [
            first,
            _rewrite_rejection(second, retry_count_for_current_macro_step=1),
            original[2],
        ]
        with self.assertRaises(TDG11IMP1Error):
            restore_imp1_checkpoint_extension(
                _recompute_checkpoint_digests(mutated_local)
            )

        mutated_time = dict(mapping)
        mutated_time["serialized_temporal_rejections"] = [
            _rewrite_rejection(first, time_hex=(1.5).hex()),
            second,
            original[2],
        ]
        with self.assertRaises(TDG11IMP1Error):
            restore_imp1_checkpoint_extension(
                _recompute_checkpoint_digests(mutated_time)
            )

        mutated_state = dict(mapping)
        mutated_state["serialized_temporal_rejections"] = [
            _rewrite_rejection(first, accepted_state_sha256=_sha("forged-early-state")),
            second,
            original[2],
        ]
        with self.assertRaises(TDG11IMP1Error):
            restore_imp1_checkpoint_extension(
                _recompute_checkpoint_digests(mutated_state)
            )

        mutated_width = dict(mapping)
        mutated_width["serialized_temporal_rejections"] = [
            first,
            _rewrite_rejection(
                second,
                attempted_width_hex=(0.125).hex(),
                next_cap_hex=_half(0.125).hex(),
            ),
            original[2],
        ]
        with self.assertRaises(TDG11IMP1Error):
            restore_imp1_checkpoint_extension(
                _recompute_checkpoint_digests(mutated_width)
            )

        mutated_counter = dict(mapping)
        mutated_counter["last_accepted_macro_step_temporal_retry_count"] = 2
        with self.assertRaises(TDG11IMP1Error):
            restore_imp1_checkpoint_extension(
                _recompute_checkpoint_digests(mutated_counter)
            )

        mutated_index = dict(mapping)
        mutated_index["current_step_index"] = START_INDEX + 4
        mutated_index["accepted_macro_step_count"] = 1
        with self.assertRaises(TDG11IMP1Error):
            restore_imp1_checkpoint_extension(
                _recompute_checkpoint_digests(mutated_index)
            )

    def test_unaccepted_identity_and_non_ascii_snapshot_fail_closed(self) -> None:
        ledger = _seed()
        mapping = imp1_checkpoint_extension(ledger)
        mapping["last_accepted_time_hex"] = (1.5).hex()
        with self.assertRaises(TDG11IMP1Error):
            restore_imp1_checkpoint_extension(_recompute_checkpoint_digests(mapping))
        mapping = imp1_checkpoint_extension(ledger)
        mapping["current_state_sha256"] = _sha("other-state")
        with self.assertRaises(TDG11IMP1Error):
            restore_imp1_checkpoint_extension(_recompute_checkpoint_digests(mapping))
        mapping = imp1_checkpoint_extension(ledger)
        mapping["inherited_snapshot"] = "not-ascii-\u00e9"
        mapping["inherited_snapshot_sha256"] = _sha("ignore")
        with self.assertRaises(TDG11IMP1Error):
            restore_imp1_checkpoint_extension(mapping)

    def test_module_imports_are_closed(self) -> None:
        names: set[str] = set()
        for node in ast.walk(ast.parse(MODULE_PATH.read_text())):
            if isinstance(node, ast.Import):
                names.update(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                names.add(node.module)
        joined = " ".join(sorted(names))
        for marker in FORBIDDEN_IMPORT_MARKERS:
            self.assertNotIn(marker, joined)
        self.assertIn("recursive_horizons.evidence_io", names)


if __name__ == "__main__":
    unittest.main()
