"""Immutable TDG11-IMP1 temporal ledger and closed checkpoint codec.

The module freezes a complete TDG6 sibling snapshot and then owns a distinct
IMP1 time, identity, retry, exact-rational debit, and rejection history.  It
does not prepare, admit, or commit a scientific step, write a store, or claim
historical origin authentication.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
import math
from numbers import Real
import re
from typing import Final, Mapping, Sequence

from recursive_horizons.evidence_io import (
    CanonicalJSONError,
    canonical_json_bytes,
    load_canonical_json,
)

from .numerical_engine import COMPARATOR_METHOD, METHODS, PRIMARY_METHOD
from .tdg6_temporal_admission_design import (
    TDG6_COMPLETE_STATE_CHANNELS,
    TDG6_METHOD_STAGE_RECORDS,
)
from .tdg6_temporal_admission_runtime import TDG6TemporalLedger


MAX_RATIONAL_BITS: Final[int] = 32768
MAX_REJECTION_RECORDS: Final[int] = 65536
MAX_CHECKPOINT_BYTES: Final[int] = 16 * 1024 * 1024
MAX_INTEGER_TEXT: Final[int] = 10000
IMP1_ACCEPTED_SUBSTEP_COUNT: Final[int] = 4
IMP1_CHANNEL_COUNT: Final[int] = len(TDG6_COMPLETE_STATE_CHANNELS)
IMP1_COMPLETE_STATE_CHANNELS: Final[tuple[str, ...]] = TDG6_COMPLETE_STATE_CHANNELS
IMP1_LEDGER_SCHEMA: Final[str] = "FGC-1-TDG11-IMP1-ledger-v1"
IMP1_LEDGER_SCHEMA_VERSION: Final[int] = 1
IMP1_REJECTION_EVENT_TYPE: Final[str] = "rejected_TDG11_IMP1_temporal_admission"
IMP1_RUNTIME_METHODS: Final[tuple[str, ...]] = METHODS
IMP1_METHOD_LABEL: Final[dict[str, str]] = {
    PRIMARY_METHOD: "RK4",
    COMPARATOR_METHOD: "SSPRK3",
}
IMP1_STAGE_RECORDS_PER_SUBSTEP: Final[dict[str, int]] = {
    method: TDG6_METHOD_STAGE_RECORDS[label]
    for method, label in IMP1_METHOD_LABEL.items()
}
IMP1_METHOD_STAGE_RECORD_COUNT: Final[dict[str, int]] = {
    method: IMP1_ACCEPTED_SUBSTEP_COUNT * count
    for method, count in IMP1_STAGE_RECORDS_PER_SUBSTEP.items()
}
IMP1_ZERO_DEBIT: Final[tuple[Fraction, ...]] = (Fraction(0),) * IMP1_CHANNEL_COUNT

IMP1_REJECTION_KEYS: Final[tuple[str, ...]] = (
    "accepted_state_sha256",
    "accepted_step_index",
    "accepted_transaction_serial",
    "assessment_sha256",
    "attempted_width_hex",
    "channel_order",
    "cumulative_temporal_retry_count",
    "event_target_hex",
    "event_type",
    "failed_channels",
    "method",
    "next_cap_hex",
    "preparation_sha256",
    "retry_count_for_current_macro_step",
    "time_hex",
)
IMP1_INHERITED_SNAPSHOT_KEYS: Final[tuple[str, ...]] = (
    "accepted_macro_step_count",
    "accumulated_debit_hex",
    "channel_order",
    "cumulative_temporal_retry_count",
    "current_macro_step_temporal_retry_count",
    "last_accepted_macro_step_temporal_retry_count",
    "last_accepted_time_hex",
    "ledger_class",
    "physical_classification",
    "serialized_temporal_rejections",
)
IMP1_CHECKPOINT_KEYS: Final[tuple[str, ...]] = (
    "accepted_macro_step_count",
    "accepted_substep_count",
    "accumulated_debit_sha256",
    "accumulated_debit_vector",
    "channel_order",
    "complete_rejection_evidence_present",
    "cumulative_temporal_retry_count",
    "current_macro_step_temporal_retry_count",
    "current_state_sha256",
    "current_step_index",
    "current_transaction_serial",
    "historical_origin_authenticated",
    "inherited_last_accepted_time_hex",
    "inherited_snapshot",
    "inherited_snapshot_sha256",
    "last_accepted_macro_step_temporal_retry_count",
    "last_accepted_time_hex",
    "method",
    "method_stage_record_count",
    "origin_receipt_sha256",
    "origin_state_sha256",
    "origin_step_index",
    "origin_transaction_serial",
    "physical_classification",
    "schema",
    "schema_version",
    "serialized_temporal_rejections",
    "serialized_temporal_rejections_sha256",
    "source_and_CFL_retry_counters_owned_elsewhere",
)

_SHA256 = re.compile(r"\A[0-9a-f]{64}\Z")
_INTEGER_TEXT = re.compile(r"\A-?(?:0|[1-9][0-9]*)\Z")
_REJECTION_KEY_SET = frozenset(IMP1_REJECTION_KEYS)
_SNAPSHOT_KEY_SET = frozenset(IMP1_INHERITED_SNAPSHOT_KEYS)
_CHECKPOINT_KEY_SET = frozenset(IMP1_CHECKPOINT_KEYS)


class TDG11IMP1Error(ValueError):
    """Malformed IMP1 ledger, rejection, or checkpoint data."""


class TDG11IMP1ResourceError(TDG11IMP1Error):
    """A prospective resource ceiling was reached before any successor existed."""


def _fail(message: str) -> None:
    raise TDG11IMP1Error(message)


def _digest(payload: bytes) -> str:
    return sha256(payload).hexdigest()


def _exact_mapping(
    value: object, names: frozenset[str], label: str
) -> dict[str, object]:
    if type(value) is not dict or set(value) != names:
        _fail(f"{label} is incomplete or broadened")
    return value


def _text(value: object, label: str) -> str:
    if type(value) is not str:
        raise TypeError(f"{label} must be text")
    return value


def _bool(value: object, expected: bool, label: str) -> bool:
    if type(value) is not bool or value is not expected:
        _fail(f"{label} crossed its numerical scope")
    return value


def _digest_text(value: object, label: str) -> str:
    text = _text(value, label)
    if _SHA256.fullmatch(text) is None:
        _fail(f"{label} must be a SHA-256 digest")
    return text


def _nonnegative_int(value: object, label: str) -> int:
    if type(value) is not int or value < 0:
        raise TypeError(f"{label} must be a nonnegative integer")
    return value


def _positive_int(value: object, label: str) -> int:
    number = _nonnegative_int(value, label)
    if number < 1:
        _fail(f"{label} must be a positive integer")
    return number


def _method(value: object) -> str:
    if type(value) is not str or value not in IMP1_RUNTIME_METHODS:
        _fail("unknown IMP1 method")
    return value


def _channel_order(value: object, label: str) -> tuple[str, ...]:
    if type(value) is not list or any(type(item) is not str for item in value):
        raise TypeError(f"{label} must be a list of channel names")
    order = tuple(value)
    if order != IMP1_COMPLETE_STATE_CHANNELS:
        _fail("IMP1 channel order differs")
    return order


def _canonical_bytes(value: object, *, label: str) -> bytes:
    try:
        raw = canonical_json_bytes(value)
    except CanonicalJSONError as error:
        if "byte bound" in str(error):
            raise TDG11IMP1ResourceError(
                f"{label} exceeds MAX_CHECKPOINT_BYTES"
            ) from error
        raise TDG11IMP1Error(f"{label} is not canonical JSON") from error
    if len(raw) > MAX_CHECKPOINT_BYTES:
        raise TDG11IMP1ResourceError(f"{label} exceeds MAX_CHECKPOINT_BYTES")
    return raw


def _load_canonical_object(text: object, *, label: str) -> object:
    payload = _text(text, label)
    try:
        raw = payload.encode("ascii")
    except UnicodeEncodeError as error:
        raise TDG11IMP1Error(f"{label} is not canonical JSON") from error
    if len(raw) > MAX_CHECKPOINT_BYTES:
        raise TDG11IMP1ResourceError(f"{label} exceeds MAX_CHECKPOINT_BYTES")
    try:
        return load_canonical_json(raw)
    except CanonicalJSONError as error:
        raise TDG11IMP1Error(f"{label} is not canonical JSON") from error


def _canonical_hex(value: object, label: str) -> str:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{label} must be a finite real scalar")
    try:
        number = float(value)
    except (OverflowError, ValueError, TypeError) as error:
        raise TDG11IMP1Error(f"{label} must be finite") from error
    if not math.isfinite(number):
        _fail(f"{label} must be finite")
    text = number.hex()
    if float.fromhex(text).hex() != text:
        _fail(f"{label} is not canonical binary64 hex")
    return text


def _parse_hex(value: object, label: str, *, positive: bool = False) -> float:
    text = _text(value, label)
    try:
        number = float.fromhex(text)
    except (OverflowError, ValueError) as error:
        raise TDG11IMP1Error(f"{label} is not canonical binary64 hex") from error
    if not math.isfinite(number) or number.hex() != text:
        _fail(f"{label} is not canonical binary64 hex")
    if positive and not (number > 0.0):
        _fail(f"{label} must be positive")
    return number


def _same_binary64(left: float, right: float) -> bool:
    return left.hex() == right.hex()


def _binary64_half(width: float) -> float:
    half = math.ldexp(width, -1)
    if (
        not math.isfinite(half)
        or not (half > 0.0)
        or 2 * Fraction.from_float(half) != Fraction.from_float(width)
    ):
        _fail("IMP1 next cap must be the exact binary64 half attempted width")
    return half


def integer_text(value: int) -> str:
    """Encode a bounded integer without changing Python's global digit limit."""

    if type(value) is not int:
        raise TypeError("integer text requires an int")
    if abs(value).bit_length() > MAX_RATIONAL_BITS:
        raise TDG11IMP1ResourceError("IMP1 integer exceeds MAX_RATIONAL_BITS")
    if value == 0:
        return "0"
    remaining = abs(value)
    groups: list[int] = []
    while remaining:
        remaining, group = divmod(remaining, 1_000_000_000)
        groups.append(group)
    digits = str(groups.pop())
    digits += "".join(format(group, "09d") for group in reversed(groups))
    return ("-" if value < 0 else "") + digits


def parse_integer_text(value: object) -> int:
    """Parse canonical decimal text without an unbounded ``int(text)`` conversion."""

    if (
        type(value) is not str
        or len(value) > MAX_INTEGER_TEXT
        or value == "-0"
        or _INTEGER_TEXT.fullmatch(value) is None
    ):
        _fail("noncanonical exact integer text")
    negative = value.startswith("-")
    digits = value[1:] if negative else value
    result = 0
    for start in range(0, len(digits), 8):
        part = digits[start : start + 8]
        result = result * (10 ** len(part)) + int(part)
        if result.bit_length() > MAX_RATIONAL_BITS:
            raise TDG11IMP1ResourceError("IMP1 integer exceeds MAX_RATIONAL_BITS")
    return -result if negative else result


def _rational_bits(value: Fraction) -> int:
    return max(abs(value.numerator).bit_length(), value.denominator.bit_length())


def _require_rational_bits(value: Fraction, label: str) -> Fraction:
    if _rational_bits(value) > MAX_RATIONAL_BITS:
        raise TDG11IMP1ResourceError(f"{label} exceeds MAX_RATIONAL_BITS")
    return value


def _production_denominator(denominator: int, label: str) -> None:
    if type(denominator) is not int or denominator <= 0:
        _fail(f"{label} denominator must be a positive integer")
    rest = denominator >> ((denominator & -denominator).bit_length() - 1)
    if rest % 3 == 0:
        rest //= 3
    if rest != 1:
        _fail(f"{label} denominator must have only primes 2 and 3 with v3<=1")


def _exact_fraction(value: object, label: str) -> Fraction:
    if type(value) is not Fraction:
        raise TypeError(f"{label} must be an exact Fraction")
    if value < 0:
        _fail(f"{label} must be nonnegative")
    _require_rational_bits(value, label)
    _production_denominator(value.denominator, label)
    return value


def _wire_fraction(value: Fraction) -> dict[str, str]:
    _require_rational_bits(value, "IMP1 rational")
    return {
        "numerator": integer_text(value.numerator),
        "denominator": integer_text(value.denominator),
    }


def _read_fraction(value: object, label: str) -> Fraction:
    item = _exact_mapping(value, frozenset({"numerator", "denominator"}), label)
    numerator = parse_integer_text(item["numerator"])
    denominator = parse_integer_text(item["denominator"])
    if denominator <= 0:
        _fail(f"{label} denominator must be positive")
    result = Fraction(numerator, denominator)
    if (result.numerator, result.denominator) != (numerator, denominator):
        _fail(f"{label} is not a reduced canonical rational")
    if result < 0:
        _fail(f"{label} must be nonnegative")
    _require_rational_bits(result, label)
    _production_denominator(result.denominator, label)
    return result


def _debit_vector(value: object, label: str) -> tuple[Fraction, ...]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise TypeError(f"{label} must be an 18-element debit vector")
    if len(value) != IMP1_CHANNEL_COUNT:
        _fail("IMP1 debit vector must have 18 channels")
    return tuple(
        _exact_fraction(item, f"{label}[{index}]") for index, item in enumerate(value)
    )


def _wired_debit(values: Sequence[Fraction]) -> list[dict[str, str]]:
    return [_wire_fraction(item) for item in values]


def _failed_channels(value: object, order: tuple[str, ...]) -> tuple[str, ...]:
    if type(value) is not list or any(type(item) is not str for item in value):
        raise TypeError("failed_channels must be a list of channel names")
    if not value:
        _fail("IMP1 rejection omitted failed channels")
    index = {name: position for position, name in enumerate(order)}
    previous = -1
    seen: set[str] = set()
    for name in value:
        position = index.get(name)
        if position is None or name in seen or position <= previous:
            _fail("IMP1 rejection channels are reordered or repeated")
        seen.add(name)
        previous = position
    return tuple(value)


def _inherited_snapshot_mapping(inherited: TDG6TemporalLedger) -> dict[str, object]:
    return {
        "ledger_class": "TDG6TemporalLedger",
        "last_accepted_time_hex": _canonical_hex(
            inherited.last_accepted_time, "inherited last_accepted_time"
        ),
        "accepted_macro_step_count": inherited.accepted_macro_step_count,
        "cumulative_temporal_retry_count": inherited.cumulative_temporal_retry_count,
        "current_macro_step_temporal_retry_count": (
            inherited.current_macro_step_temporal_retry_count
        ),
        "last_accepted_macro_step_temporal_retry_count": (
            inherited.last_accepted_macro_step_temporal_retry_count
        ),
        "channel_order": list(IMP1_COMPLETE_STATE_CHANNELS),
        "accumulated_debit_hex": [
            _canonical_hex(item, f"inherited debit[{index}]")
            for index, item in enumerate(inherited.accumulated_debit_vector)
        ],
        "serialized_temporal_rejections": list(
            inherited.serialized_temporal_rejections
        ),
        "physical_classification": False,
    }


def encode_inherited_snapshot(inherited: TDG6TemporalLedger) -> str:
    """Freeze every TDG6 time, debit-bit, count, and rejection string."""

    if not isinstance(inherited, TDG6TemporalLedger):
        raise TypeError("inherited must be TDG6TemporalLedger")
    return _canonical_bytes(
        _inherited_snapshot_mapping(inherited),
        label="IMP1 inherited snapshot",
    ).decode("ascii")


def decode_inherited_snapshot(snapshot: str) -> TDG6TemporalLedger:
    """Restore the frozen TDG6 sibling without promoting it to IMP1 state."""

    mapping = _exact_mapping(
        _load_canonical_object(snapshot, label="IMP1 inherited snapshot"),
        _SNAPSHOT_KEY_SET,
        "IMP1 inherited snapshot",
    )
    if mapping["ledger_class"] != "TDG6TemporalLedger":
        _fail("IMP1 inherited snapshot is not a TDG6TemporalLedger")
    _bool(
        mapping["physical_classification"], False, "inherited physical_classification"
    )
    _channel_order(mapping["channel_order"], "inherited channel_order")
    debit_hex = mapping["accumulated_debit_hex"]
    if type(debit_hex) is not list or len(debit_hex) != IMP1_CHANNEL_COUNT:
        _fail("inherited debit vector must have 18 channels")
    rejections = mapping["serialized_temporal_rejections"]
    if type(rejections) is not list or any(
        type(item) is not str for item in rejections
    ):
        raise TypeError("inherited serialized rejections must be text")
    restored = TDG6TemporalLedger(
        last_accepted_time=_parse_hex(
            mapping["last_accepted_time_hex"], "inherited last_accepted_time_hex"
        ),
        accepted_macro_step_count=_nonnegative_int(
            mapping["accepted_macro_step_count"], "inherited accepted_macro_step_count"
        ),
        cumulative_temporal_retry_count=_nonnegative_int(
            mapping["cumulative_temporal_retry_count"],
            "inherited cumulative_temporal_retry_count",
        ),
        current_macro_step_temporal_retry_count=_nonnegative_int(
            mapping["current_macro_step_temporal_retry_count"],
            "inherited current_macro_step_temporal_retry_count",
        ),
        last_accepted_macro_step_temporal_retry_count=_nonnegative_int(
            mapping["last_accepted_macro_step_temporal_retry_count"],
            "inherited last_accepted_macro_step_temporal_retry_count",
        ),
        accumulated_debit_vector=tuple(
            _parse_hex(item, f"inherited debit[{index}]")
            for index, item in enumerate(debit_hex)
        ),
        serialized_temporal_rejections=tuple(rejections),
    )
    if encode_inherited_snapshot(restored) != snapshot:
        _fail("IMP1 inherited snapshot is not a closed TDG6 encoding")
    return restored


@dataclass(frozen=True, slots=True)
class TDG11IMP1TemporalRejection:
    """Explicit field contract for one IMP1 temporal-admission rejection."""

    event_type: str
    method: str
    time_hex: str
    event_target_hex: str
    attempted_width_hex: str
    next_cap_hex: str
    accepted_state_sha256: str
    accepted_step_index: int
    accepted_transaction_serial: int
    assessment_sha256: str
    preparation_sha256: str
    channel_order: tuple[str, ...]
    failed_channels: tuple[str, ...]
    retry_count_for_current_macro_step: int
    cumulative_temporal_retry_count: int

    def __post_init__(self) -> None:
        if self.event_type != IMP1_REJECTION_EVENT_TYPE:
            _fail("IMP1 rejection evidence has the wrong event type")
        object.__setattr__(self, "method", _method(self.method))
        time = _parse_hex(self.time_hex, "time_hex")
        target = _parse_hex(self.event_target_hex, "event_target_hex")
        if not (target > time):
            _fail("IMP1 event target must be a strictly later finite instant")
        attempted = _parse_hex(
            self.attempted_width_hex, "attempted_width_hex", positive=True
        )
        next_cap = _parse_hex(self.next_cap_hex, "next_cap_hex", positive=True)
        half = _binary64_half(attempted)
        if not _same_binary64(next_cap, half):
            _fail("IMP1 next cap must be the exact binary64 half attempted width")
        span = Fraction.from_float(target) - Fraction.from_float(time)
        if Fraction.from_float(attempted) > span:
            _fail("IMP1 attempted width exceeds the exact event-target interval")
        object.__setattr__(
            self,
            "accepted_state_sha256",
            _digest_text(self.accepted_state_sha256, "accepted_state_sha256"),
        )
        object.__setattr__(
            self,
            "accepted_step_index",
            _nonnegative_int(self.accepted_step_index, "accepted_step_index"),
        )
        object.__setattr__(
            self,
            "accepted_transaction_serial",
            _nonnegative_int(
                self.accepted_transaction_serial, "accepted_transaction_serial"
            ),
        )
        object.__setattr__(
            self,
            "assessment_sha256",
            _digest_text(self.assessment_sha256, "assessment_sha256"),
        )
        object.__setattr__(
            self,
            "preparation_sha256",
            _digest_text(self.preparation_sha256, "preparation_sha256"),
        )
        if type(self.channel_order) not in (tuple, list):
            raise TypeError("channel_order must be a sequence of channel names")
        order = _channel_order(list(self.channel_order), "channel_order")
        object.__setattr__(self, "channel_order", order)
        if type(self.failed_channels) not in (tuple, list):
            raise TypeError("failed_channels must be a sequence of channel names")
        object.__setattr__(
            self,
            "failed_channels",
            _failed_channels(list(self.failed_channels), order),
        )
        object.__setattr__(
            self,
            "retry_count_for_current_macro_step",
            _positive_int(
                self.retry_count_for_current_macro_step,
                "retry_count_for_current_macro_step",
            ),
        )
        object.__setattr__(
            self,
            "cumulative_temporal_retry_count",
            _positive_int(
                self.cumulative_temporal_retry_count,
                "cumulative_temporal_retry_count",
            ),
        )
        if (
            self.cumulative_temporal_retry_count
            < self.retry_count_for_current_macro_step
        ):
            _fail("IMP1 cumulative retry count is inconsistent")

    def as_mapping(self) -> dict[str, object]:
        return {
            "event_type": self.event_type,
            "method": self.method,
            "time_hex": self.time_hex,
            "event_target_hex": self.event_target_hex,
            "attempted_width_hex": self.attempted_width_hex,
            "next_cap_hex": self.next_cap_hex,
            "accepted_state_sha256": self.accepted_state_sha256,
            "accepted_step_index": self.accepted_step_index,
            "accepted_transaction_serial": self.accepted_transaction_serial,
            "assessment_sha256": self.assessment_sha256,
            "preparation_sha256": self.preparation_sha256,
            "channel_order": list(self.channel_order),
            "failed_channels": list(self.failed_channels),
            "retry_count_for_current_macro_step": self.retry_count_for_current_macro_step,
            "cumulative_temporal_retry_count": self.cumulative_temporal_retry_count,
        }

    @classmethod
    def from_mapping(
        cls, mapping: Mapping[str, object]
    ) -> "TDG11IMP1TemporalRejection":
        item = _exact_mapping(mapping, _REJECTION_KEY_SET, "IMP1 rejection")
        return cls(
            event_type=item["event_type"],
            method=item["method"],
            time_hex=item["time_hex"],
            event_target_hex=item["event_target_hex"],
            attempted_width_hex=item["attempted_width_hex"],
            next_cap_hex=item["next_cap_hex"],
            accepted_state_sha256=item["accepted_state_sha256"],
            accepted_step_index=item["accepted_step_index"],
            accepted_transaction_serial=item["accepted_transaction_serial"],
            assessment_sha256=item["assessment_sha256"],
            preparation_sha256=item["preparation_sha256"],
            channel_order=item["channel_order"],
            failed_channels=item["failed_channels"],
            retry_count_for_current_macro_step=item[
                "retry_count_for_current_macro_step"
            ],
            cumulative_temporal_retry_count=item["cumulative_temporal_retry_count"],
        )

    def canonical_text(self) -> str:
        return _canonical_bytes(self.as_mapping(), label="IMP1 rejection").decode(
            "ascii"
        )


def _parse_rejection_text(record: str) -> TDG11IMP1TemporalRejection:
    parsed = TDG11IMP1TemporalRejection.from_mapping(
        _load_canonical_object(record, label="IMP1 rejection")
    )
    if parsed.canonical_text() != record:
        _fail("IMP1 rejection evidence is not canonical JSON")
    return parsed


def _boundary_macro(step_index: int, origin_step_index: int) -> int:
    delta = step_index - origin_step_index
    if delta < 0 or delta % IMP1_ACCEPTED_SUBSTEP_COUNT != 0:
        _fail("IMP1 rejection step index is not aligned to origin")
    return delta // IMP1_ACCEPTED_SUBSTEP_COUNT


def _validate_rejection_history(
    *,
    method: str,
    origin_time: float,
    origin_state: str,
    origin_step_index: int,
    origin_serial: int,
    current_time: float,
    current_state: str,
    current_step_index: int,
    current_serial: int,
    accepted_macro_step_count: int,
    current_retry: int,
    last_accepted_retry: int,
    stage_records: int,
    records: tuple[TDG11IMP1TemporalRejection, ...],
) -> None:
    groups: list[tuple[int, list[TDG11IMP1TemporalRejection]]] = []
    for position, record in enumerate(records):
        if record.method != method:
            _fail("IMP1 rejection evidence has the wrong method")
        if record.cumulative_temporal_retry_count != position + 1:
            _fail("IMP1 rejection cumulative ordinal does not match history")
        macro = _boundary_macro(record.accepted_step_index, origin_step_index)
        if macro > accepted_macro_step_count:
            _fail("IMP1 rejection boundary macro index exceeds accepted count")
        expected_serial = origin_serial + stage_records * macro
        if record.accepted_transaction_serial != expected_serial:
            _fail("IMP1 rejection serial is not aligned to origin")
        rec_time = _parse_hex(record.time_hex, "time_hex")
        if rec_time < origin_time or rec_time > current_time:
            _fail("IMP1 rejection time is outside the origin/current interval")
        if macro == 0:
            if not _same_binary64(rec_time, origin_time):
                _fail("IMP1 origin-boundary rejection time differs")
            if record.accepted_state_sha256 != origin_state:
                _fail("IMP1 origin-boundary rejection state differs")
            if record.accepted_step_index != origin_step_index:
                _fail("IMP1 origin-boundary rejection step index differs")
            if record.accepted_transaction_serial != origin_serial:
                _fail("IMP1 origin-boundary rejection serial differs")
        elif macro == accepted_macro_step_count:
            if not _same_binary64(rec_time, current_time):
                _fail("IMP1 active-retry boundary time differs")
            if record.accepted_state_sha256 != current_state:
                _fail("IMP1 active-retry boundary state differs")
            if record.accepted_step_index != current_step_index:
                _fail("IMP1 active-retry boundary step index differs")
            if record.accepted_transaction_serial != current_serial:
                _fail("IMP1 active-retry boundary serial differs")
        elif not (origin_time < rec_time < current_time):
            _fail("IMP1 historical rejection time is not between origin and current")
        if groups and groups[-1][0] == macro:
            groups[-1][1].append(record)
            continue
        if groups and macro <= groups[-1][0]:
            _fail("IMP1 rejection macro groups do not advance")
        groups.append((macro, [record]))

    for macro, items in groups:
        head = items[0]
        for index, record in enumerate(items):
            if record.retry_count_for_current_macro_step != index + 1:
                _fail("IMP1 per-macro retry numbering is not contiguous from 1")
            if (
                record.time_hex != head.time_hex
                or record.accepted_state_sha256 != head.accepted_state_sha256
                or record.accepted_step_index != head.accepted_step_index
                or record.accepted_transaction_serial
                != head.accepted_transaction_serial
                or record.event_target_hex != head.event_target_hex
            ):
                _fail("IMP1 retry group boundary identity or event target differs")
            if index == 0:
                continue
            previous_cap = _parse_hex(items[index - 1].next_cap_hex, "next_cap_hex")
            attempted = _parse_hex(
                record.attempted_width_hex, "attempted_width_hex", positive=True
            )
            if attempted > previous_cap:
                _fail("IMP1 attempted width exceeds the previous next cap")
    for left, right in zip(groups, groups[1:]):
        earlier = left[1][0]
        later = right[1][0]
        earlier_time = _parse_hex(earlier.time_hex, "time_hex")
        later_time = _parse_hex(later.time_hex, "time_hex")
        if not (later_time > earlier_time):
            _fail("IMP1 later rejection groups must advance time")
        if later.accepted_step_index <= earlier.accepted_step_index:
            _fail("IMP1 later rejection groups must advance step index")

    completed = {
        macro: len(items)
        for macro, items in groups
        if macro < accepted_macro_step_count
    }
    active_items = [
        items for macro, items in groups if macro == accepted_macro_step_count
    ]
    active_size = len(active_items[0]) if active_items else 0
    if len(active_items) > 1:
        _fail("IMP1 active retry group is repeated")
    if current_retry != active_size:
        _fail("IMP1 current retry count does not match the final history group")
    expected_last_accepted = (
        0
        if accepted_macro_step_count == 0
        else completed.get(accepted_macro_step_count - 1, 0)
    )
    if last_accepted_retry != expected_last_accepted:
        _fail("IMP1 last-accepted retry count does not match the final history group")
    if current_retry and (
        not active_items
        or not _same_binary64(
            _parse_hex(active_items[0][-1].time_hex, "time_hex"), current_time
        )
        or active_items[0][-1].accepted_state_sha256 != current_state
        or active_items[0][-1].accepted_step_index != current_step_index
        or active_items[0][-1].accepted_transaction_serial != current_serial
    ):
        _fail("IMP1 active-retry latest boundary does not match the live identity")


@dataclass(frozen=True, slots=True)
class TDG11IMP1Ledger:
    """Immutable IMP1 ledger with a frozen TDG6 sibling snapshot."""

    method: str
    inherited_snapshot: str
    origin_receipt_sha256: str
    origin_state_sha256: str
    origin_step_index: int
    origin_transaction_serial: int
    last_accepted_time: float
    current_state_sha256: str
    current_step_index: int
    current_transaction_serial: int
    accepted_macro_step_count: int
    cumulative_temporal_retry_count: int
    current_macro_step_temporal_retry_count: int
    last_accepted_macro_step_temporal_retry_count: int
    accumulated_debit_vector: tuple[Fraction, ...]
    serialized_temporal_rejections: tuple[str, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "method", _method(self.method))
        snapshot = _text(self.inherited_snapshot, "inherited_snapshot")
        inherited = decode_inherited_snapshot(snapshot)
        if inherited.current_macro_step_temporal_retry_count != 0:
            _fail("IMP1 inherited snapshot has an active temporal retry")
        object.__setattr__(
            self,
            "origin_receipt_sha256",
            _digest_text(self.origin_receipt_sha256, "origin_receipt_sha256"),
        )
        object.__setattr__(
            self,
            "origin_state_sha256",
            _digest_text(self.origin_state_sha256, "origin_state_sha256"),
        )
        object.__setattr__(
            self,
            "origin_step_index",
            _nonnegative_int(self.origin_step_index, "origin_step_index"),
        )
        object.__setattr__(
            self,
            "origin_transaction_serial",
            _nonnegative_int(
                self.origin_transaction_serial, "origin_transaction_serial"
            ),
        )
        object.__setattr__(
            self,
            "last_accepted_time",
            _parse_hex(
                _canonical_hex(self.last_accepted_time, "last_accepted_time"),
                "last_accepted_time",
            ),
        )
        object.__setattr__(
            self,
            "current_state_sha256",
            _digest_text(self.current_state_sha256, "current_state_sha256"),
        )
        object.__setattr__(
            self,
            "current_step_index",
            _nonnegative_int(self.current_step_index, "current_step_index"),
        )
        object.__setattr__(
            self,
            "current_transaction_serial",
            _nonnegative_int(
                self.current_transaction_serial, "current_transaction_serial"
            ),
        )
        for name in (
            "accepted_macro_step_count",
            "cumulative_temporal_retry_count",
            "current_macro_step_temporal_retry_count",
            "last_accepted_macro_step_temporal_retry_count",
        ):
            object.__setattr__(self, name, _nonnegative_int(getattr(self, name), name))
        stage_records = IMP1_METHOD_STAGE_RECORD_COUNT[self.method]
        expected_index = (
            self.origin_step_index
            + IMP1_ACCEPTED_SUBSTEP_COUNT * self.accepted_macro_step_count
        )
        expected_serial = (
            self.origin_transaction_serial
            + stage_records * self.accepted_macro_step_count
        )
        if self.current_step_index != expected_index:
            _fail("IMP1 step index is not origin plus four substeps per accepted macro")
        if self.current_transaction_serial != expected_serial:
            _fail(
                "IMP1 serial is not origin plus method-owned stage records per accepted macro"
            )
        origin_time = inherited.last_accepted_time
        if self.accepted_macro_step_count == 0:
            if not _same_binary64(self.last_accepted_time, origin_time):
                _fail("IMP1 unaccepted time must equal the inherited origin")
            if self.current_state_sha256 != self.origin_state_sha256:
                _fail("IMP1 unaccepted state must equal the origin state")
            if self.last_accepted_macro_step_temporal_retry_count != 0:
                _fail("IMP1 unaccepted last-accepted retry count must be zero")
        elif not (self.last_accepted_time > origin_time):
            _fail("IMP1 accepted time must be strictly after the inherited origin")
        if (
            self.current_macro_step_temporal_retry_count
            > self.cumulative_temporal_retry_count
            or self.last_accepted_macro_step_temporal_retry_count
            > self.cumulative_temporal_retry_count
        ):
            _fail("IMP1 per-step retry count exceeds its cumulative ledger")
        object.__setattr__(
            self,
            "accumulated_debit_vector",
            _debit_vector(self.accumulated_debit_vector, "accumulated_debit_vector"),
        )
        if (
            self.accepted_macro_step_count == 0
            and self.accumulated_debit_vector != IMP1_ZERO_DEBIT
        ):
            _fail("IMP1 unaccepted debit vector must be all zero")
        if type(self.serialized_temporal_rejections) not in (tuple, list):
            raise TypeError("serialized_temporal_rejections must be a sequence of text")
        if any(type(item) is not str for item in self.serialized_temporal_rejections):
            raise TypeError("IMP1 serialized rejection must be text")
        if len(self.serialized_temporal_rejections) > MAX_REJECTION_RECORDS:
            raise TDG11IMP1ResourceError(
                "IMP1 rejection history exceeds MAX_REJECTION_RECORDS"
            )
        if (
            len(self.serialized_temporal_rejections)
            != self.cumulative_temporal_retry_count
        ):
            _fail("IMP1 retry count differs from serialized evidence")
        parsed = tuple(
            _parse_rejection_text(item) for item in self.serialized_temporal_rejections
        )
        object.__setattr__(
            self,
            "serialized_temporal_rejections",
            tuple(item.canonical_text() for item in parsed),
        )
        _validate_rejection_history(
            method=self.method,
            origin_time=origin_time,
            origin_state=self.origin_state_sha256,
            origin_step_index=self.origin_step_index,
            origin_serial=self.origin_transaction_serial,
            current_time=self.last_accepted_time,
            current_state=self.current_state_sha256,
            current_step_index=self.current_step_index,
            current_serial=self.current_transaction_serial,
            accepted_macro_step_count=self.accepted_macro_step_count,
            current_retry=self.current_macro_step_temporal_retry_count,
            last_accepted_retry=self.last_accepted_macro_step_temporal_retry_count,
            stage_records=stage_records,
            records=parsed,
        )

    @property
    def current_time(self) -> float:
        return self.last_accepted_time

    @property
    def historical_origin_authenticated(self) -> bool:
        return False

    @property
    def accepted_substep_count(self) -> int:
        return IMP1_ACCEPTED_SUBSTEP_COUNT

    @property
    def method_stage_record_count(self) -> int:
        return IMP1_METHOD_STAGE_RECORD_COUNT[self.method]


def inherited_tdg6_ledger(ledger: TDG11IMP1Ledger) -> TDG6TemporalLedger:
    """Return the frozen TDG6 sibling; it is never the live IMP1 identity."""

    if not isinstance(ledger, TDG11IMP1Ledger):
        raise TypeError("ledger must be TDG11IMP1Ledger")
    return decode_inherited_snapshot(ledger.inherited_snapshot)


def _checked_ledger(ledger: object) -> TDG11IMP1Ledger:
    if not isinstance(ledger, TDG11IMP1Ledger):
        raise TypeError("ledger must be TDG11IMP1Ledger")
    return ledger


def _successor(ledger: TDG11IMP1Ledger, **changes: object) -> TDG11IMP1Ledger:
    payload = {
        "method": ledger.method,
        "inherited_snapshot": ledger.inherited_snapshot,
        "origin_receipt_sha256": ledger.origin_receipt_sha256,
        "origin_state_sha256": ledger.origin_state_sha256,
        "origin_step_index": ledger.origin_step_index,
        "origin_transaction_serial": ledger.origin_transaction_serial,
        "last_accepted_time": ledger.last_accepted_time,
        "current_state_sha256": ledger.current_state_sha256,
        "current_step_index": ledger.current_step_index,
        "current_transaction_serial": ledger.current_transaction_serial,
        "accepted_macro_step_count": ledger.accepted_macro_step_count,
        "cumulative_temporal_retry_count": ledger.cumulative_temporal_retry_count,
        "current_macro_step_temporal_retry_count": (
            ledger.current_macro_step_temporal_retry_count
        ),
        "last_accepted_macro_step_temporal_retry_count": (
            ledger.last_accepted_macro_step_temporal_retry_count
        ),
        "accumulated_debit_vector": ledger.accumulated_debit_vector,
        "serialized_temporal_rejections": ledger.serialized_temporal_rejections,
    }
    payload.update(changes)
    successor = TDG11IMP1Ledger(**payload)
    imp1_checkpoint_extension(successor)
    return successor


def seed_imp1_ledger(
    inherited: TDG6TemporalLedger,
    *,
    method: str,
    state_sha256: str,
    step_index: int,
    transaction_serial: int,
    origin_receipt_sha256: str,
) -> TDG11IMP1Ledger:
    """Create an IMP1 ledger from a frozen TDG6 sibling.

    The snapshot is bitwise-preserving.  An active inherited temporal retry is
    refused.  Origin authentication is not claimed; PROTO19 owns that later.
    """

    if not isinstance(inherited, TDG6TemporalLedger):
        raise TypeError("inherited must be TDG6TemporalLedger")
    if inherited.current_macro_step_temporal_retry_count != 0:
        _fail("IMP1 seed requires no active inherited temporal retry")
    snapshot = encode_inherited_snapshot(inherited)
    ledger = TDG11IMP1Ledger(
        method=_method(method),
        inherited_snapshot=snapshot,
        origin_receipt_sha256=origin_receipt_sha256,
        origin_state_sha256=state_sha256,
        origin_step_index=step_index,
        origin_transaction_serial=transaction_serial,
        last_accepted_time=_parse_hex(
            _canonical_hex(
                inherited.last_accepted_time, "inherited last_accepted_time"
            ),
            "inherited last_accepted_time",
        ),
        current_state_sha256=state_sha256,
        current_step_index=step_index,
        current_transaction_serial=transaction_serial,
        accepted_macro_step_count=0,
        cumulative_temporal_retry_count=0,
        current_macro_step_temporal_retry_count=0,
        last_accepted_macro_step_temporal_retry_count=0,
        accumulated_debit_vector=IMP1_ZERO_DEBIT,
        serialized_temporal_rejections=(),
    )
    imp1_checkpoint_extension(ledger)
    return ledger


def _rejection_record(record: object) -> TDG11IMP1TemporalRejection:
    if isinstance(record, TDG11IMP1TemporalRejection):
        return TDG11IMP1TemporalRejection.from_mapping(record.as_mapping())
    if isinstance(record, Mapping):
        return TDG11IMP1TemporalRejection.from_mapping(record)
    raise TypeError("record must be TDG11IMP1TemporalRejection or a closed mapping")


def append_imp1_rejection(
    ledger: TDG11IMP1Ledger,
    record: TDG11IMP1TemporalRejection | Mapping[str, object],
) -> TDG11IMP1Ledger:
    """Return a new ledger with one canonical rejection appended.

    The runtime owns durable-sink ordering and retry exhaustion.  History is
    never truncated.
    """

    current = _checked_ledger(ledger)
    if current.cumulative_temporal_retry_count >= MAX_REJECTION_RECORDS:
        raise TDG11IMP1ResourceError(
            "IMP1 rejection history exceeds MAX_REJECTION_RECORDS"
        )
    item = _rejection_record(record)
    expected_current = current.current_macro_step_temporal_retry_count + 1
    expected_cumulative = current.cumulative_temporal_retry_count + 1
    if item.method != current.method:
        _fail("IMP1 rejection evidence has the wrong method")
    if item.time_hex != current.last_accepted_time.hex():
        _fail("IMP1 rejection time differs from the accepted boundary")
    if item.accepted_state_sha256 != current.current_state_sha256:
        _fail("IMP1 rejection accepted state differs")
    if item.accepted_step_index != current.current_step_index:
        _fail("IMP1 rejection accepted step index differs")
    if item.accepted_transaction_serial != current.current_transaction_serial:
        _fail("IMP1 rejection accepted transaction serial differs")
    if item.retry_count_for_current_macro_step != expected_current:
        _fail("IMP1 rejection current retry count is not incremented")
    if item.cumulative_temporal_retry_count != expected_cumulative:
        _fail("IMP1 rejection cumulative retry count is not incremented")
    if current.current_macro_step_temporal_retry_count > 0:
        previous = _parse_rejection_text(current.serialized_temporal_rejections[-1])
        if item.event_target_hex != previous.event_target_hex:
            _fail("IMP1 retry group event target differs")
        attempted = _parse_hex(
            item.attempted_width_hex, "attempted_width_hex", positive=True
        )
        previous_cap = _parse_hex(previous.next_cap_hex, "next_cap_hex", positive=True)
        if attempted > previous_cap:
            _fail("IMP1 attempted width exceeds the previous next cap")
    serialized = item.canonical_text()
    return _successor(
        current,
        cumulative_temporal_retry_count=expected_cumulative,
        current_macro_step_temporal_retry_count=expected_current,
        serialized_temporal_rejections=(
            *current.serialized_temporal_rejections,
            serialized,
        ),
    )


def accept_imp1_step(
    ledger: TDG11IMP1Ledger,
    *,
    final_time: object,
    state_sha256: str,
    step_index: int,
    transaction_serial: int,
    debit_vector: Sequence[Fraction],
    assessment_sha256: str,
) -> TDG11IMP1Ledger:
    """Accept one fine four-substep IMP1 macro step into a new ledger.

    ``step_index`` must advance by four accepted substeps and
    ``transaction_serial`` by the method-owned 20 RK4 or 16 SSPRK3 stage
    records.  Only the current IMP1 retry count is reset.
    """

    current = _checked_ledger(ledger)
    _digest_text(assessment_sha256, "assessment_sha256")
    final = _parse_hex(_canonical_hex(final_time, "final_time"), "final_time")
    if not (final > current.last_accepted_time):
        _fail("IMP1 accepted time must be a strictly later finite instant")
    next_index = _nonnegative_int(step_index, "step_index")
    next_serial = _nonnegative_int(transaction_serial, "transaction_serial")
    if next_index != current.current_step_index + IMP1_ACCEPTED_SUBSTEP_COUNT:
        _fail("IMP1 accepted step must cover exactly four substeps")
    expected_serial = (
        current.current_transaction_serial
        + IMP1_METHOD_STAGE_RECORD_COUNT[current.method]
    )
    if next_serial != expected_serial:
        _fail("IMP1 accepted serial must cover the method-owned stage records")
    added = _debit_vector(debit_vector, "debit_vector")
    accumulated: list[Fraction] = []
    for index, (old, extra) in enumerate(
        zip(current.accumulated_debit_vector, added, strict=True)
    ):
        total = old + extra
        _require_rational_bits(total, f"accumulated debit[{index}]")
        _production_denominator(total.denominator, f"accumulated debit[{index}]")
        accumulated.append(total)
    return _successor(
        current,
        last_accepted_time=final,
        current_state_sha256=_digest_text(state_sha256, "state_sha256"),
        current_step_index=next_index,
        current_transaction_serial=next_serial,
        accepted_macro_step_count=current.accepted_macro_step_count + 1,
        current_macro_step_temporal_retry_count=0,
        last_accepted_macro_step_temporal_retry_count=(
            current.current_macro_step_temporal_retry_count
        ),
        accumulated_debit_vector=tuple(accumulated),
        serialized_temporal_rejections=current.serialized_temporal_rejections,
    )


def combined_diagnostic_debit_vector(
    ledger: TDG11IMP1Ledger,
) -> tuple[Fraction, ...]:
    """Exact diagnostic upper bound: inherited binary64 bits plus IMP1 debit.

    The inherited float vector is not rewritten.  The sum is comparison
    evidence only.
    """

    current = _checked_ledger(ledger)
    inherited = inherited_tdg6_ledger(current)
    combined: list[Fraction] = []
    for index, (old, extra) in enumerate(
        zip(
            inherited.accumulated_debit_vector,
            current.accumulated_debit_vector,
            strict=True,
        )
    ):
        total = Fraction.from_float(old) + extra
        combined.append(_require_rational_bits(total, f"combined debit[{index}]"))
    return tuple(combined)


def imp1_checkpoint_extension(ledger: TDG11IMP1Ledger) -> dict[str, object]:
    """Encode the closed IMP1 checkpoint mapping with redundant identities."""

    current = _checked_ledger(ledger)
    inherited = inherited_tdg6_ledger(current)
    debit = _wired_debit(current.accumulated_debit_vector)
    rejections = list(current.serialized_temporal_rejections)
    mapping: dict[str, object] = {
        "schema": IMP1_LEDGER_SCHEMA,
        "schema_version": IMP1_LEDGER_SCHEMA_VERSION,
        "method": current.method,
        "inherited_snapshot": current.inherited_snapshot,
        "inherited_snapshot_sha256": _digest(
            current.inherited_snapshot.encode("ascii")
        ),
        "inherited_last_accepted_time_hex": inherited.last_accepted_time.hex(),
        "origin_receipt_sha256": current.origin_receipt_sha256,
        "origin_state_sha256": current.origin_state_sha256,
        "origin_step_index": current.origin_step_index,
        "origin_transaction_serial": current.origin_transaction_serial,
        "last_accepted_time_hex": current.last_accepted_time.hex(),
        "current_state_sha256": current.current_state_sha256,
        "current_step_index": current.current_step_index,
        "current_transaction_serial": current.current_transaction_serial,
        "accepted_macro_step_count": current.accepted_macro_step_count,
        "accepted_substep_count": IMP1_ACCEPTED_SUBSTEP_COUNT,
        "method_stage_record_count": IMP1_METHOD_STAGE_RECORD_COUNT[current.method],
        "cumulative_temporal_retry_count": current.cumulative_temporal_retry_count,
        "current_macro_step_temporal_retry_count": (
            current.current_macro_step_temporal_retry_count
        ),
        "last_accepted_macro_step_temporal_retry_count": (
            current.last_accepted_macro_step_temporal_retry_count
        ),
        "channel_order": list(IMP1_COMPLETE_STATE_CHANNELS),
        "accumulated_debit_vector": debit,
        "accumulated_debit_sha256": _digest(
            _canonical_bytes(debit, label="IMP1 debit")
        ),
        "serialized_temporal_rejections": rejections,
        "serialized_temporal_rejections_sha256": _digest(
            _canonical_bytes(rejections, label="IMP1 rejections")
        ),
        "historical_origin_authenticated": False,
        "physical_classification": False,
        "complete_rejection_evidence_present": True,
        "source_and_CFL_retry_counters_owned_elsewhere": True,
    }
    if set(mapping) != _CHECKPOINT_KEY_SET:
        _fail("IMP1 checkpoint metadata is incomplete or broadened")
    _canonical_bytes(mapping, label="IMP1 checkpoint")
    return mapping


def restore_imp1_checkpoint_extension(mapping: Mapping[str, object]) -> TDG11IMP1Ledger:
    """Restore and revalidate a closed IMP1 checkpoint mapping."""

    item = _exact_mapping(mapping, _CHECKPOINT_KEY_SET, "IMP1 checkpoint")
    if item["schema"] != IMP1_LEDGER_SCHEMA:
        _fail("IMP1 checkpoint schema differs")
    if item["schema_version"] != IMP1_LEDGER_SCHEMA_VERSION:
        _fail("IMP1 checkpoint schema differs")
    method = _method(item["method"])
    snapshot = _text(item["inherited_snapshot"], "inherited_snapshot")
    try:
        snapshot_bytes = snapshot.encode("ascii")
    except UnicodeEncodeError as error:
        raise TDG11IMP1Error("IMP1 inherited snapshot is not canonical JSON") from error
    if _digest(snapshot_bytes) != _digest_text(
        item["inherited_snapshot_sha256"], "inherited_snapshot_sha256"
    ):
        _fail("IMP1 inherited snapshot hash differs")
    inherited = decode_inherited_snapshot(snapshot)
    if inherited.last_accepted_time.hex() != _text(
        item["inherited_last_accepted_time_hex"], "inherited_last_accepted_time_hex"
    ):
        _fail("IMP1 inherited time identity differs")
    _bool(
        item["historical_origin_authenticated"],
        False,
        "historical_origin_authenticated",
    )
    _bool(item["physical_classification"], False, "physical_classification")
    _bool(
        item["complete_rejection_evidence_present"],
        True,
        "complete_rejection_evidence_present",
    )
    _bool(
        item["source_and_CFL_retry_counters_owned_elsewhere"],
        True,
        "source_and_CFL_retry_counters_owned_elsewhere",
    )
    if item["accepted_substep_count"] != IMP1_ACCEPTED_SUBSTEP_COUNT:
        _fail("IMP1 accepted substep count differs")
    if item["method_stage_record_count"] != IMP1_METHOD_STAGE_RECORD_COUNT[method]:
        _fail("IMP1 method-owned stage record count differs")
    _channel_order(item["channel_order"], "channel_order")
    debit_wire = item["accumulated_debit_vector"]
    if type(debit_wire) is not list or len(debit_wire) != IMP1_CHANNEL_COUNT:
        _fail("IMP1 debit vector must have 18 channels")
    debit = tuple(
        _read_fraction(entry, f"accumulated_debit_vector[{index}]")
        for index, entry in enumerate(debit_wire)
    )
    if _digest(
        _canonical_bytes(_wired_debit(debit), label="IMP1 debit")
    ) != _digest_text(item["accumulated_debit_sha256"], "accumulated_debit_sha256"):
        _fail("IMP1 debit hash differs")
    rejections = item["serialized_temporal_rejections"]
    if type(rejections) is not list or any(
        type(entry) is not str for entry in rejections
    ):
        raise TypeError("IMP1 serialized rejections must be text")
    if _digest(_canonical_bytes(rejections, label="IMP1 rejections")) != _digest_text(
        item["serialized_temporal_rejections_sha256"],
        "serialized_temporal_rejections_sha256",
    ):
        _fail("IMP1 rejection hash differs")
    restored = TDG11IMP1Ledger(
        method=method,
        inherited_snapshot=snapshot,
        origin_receipt_sha256=item["origin_receipt_sha256"],
        origin_state_sha256=item["origin_state_sha256"],
        origin_step_index=item["origin_step_index"],
        origin_transaction_serial=item["origin_transaction_serial"],
        last_accepted_time=_parse_hex(
            item["last_accepted_time_hex"], "last_accepted_time_hex"
        ),
        current_state_sha256=item["current_state_sha256"],
        current_step_index=item["current_step_index"],
        current_transaction_serial=item["current_transaction_serial"],
        accepted_macro_step_count=item["accepted_macro_step_count"],
        cumulative_temporal_retry_count=item["cumulative_temporal_retry_count"],
        current_macro_step_temporal_retry_count=item[
            "current_macro_step_temporal_retry_count"
        ],
        last_accepted_macro_step_temporal_retry_count=item[
            "last_accepted_macro_step_temporal_retry_count"
        ],
        accumulated_debit_vector=debit,
        serialized_temporal_rejections=tuple(rejections),
    )
    encoded = imp1_checkpoint_extension(restored)
    if _canonical_bytes(encoded, label="IMP1 checkpoint") != _canonical_bytes(
        dict(item), label="IMP1 checkpoint"
    ):
        _fail("IMP1 checkpoint identities are not redundant")
    return restored


__all__ = [
    "IMP1_ACCEPTED_SUBSTEP_COUNT",
    "IMP1_CHANNEL_COUNT",
    "IMP1_CHECKPOINT_KEYS",
    "IMP1_COMPLETE_STATE_CHANNELS",
    "IMP1_INHERITED_SNAPSHOT_KEYS",
    "IMP1_LEDGER_SCHEMA",
    "IMP1_LEDGER_SCHEMA_VERSION",
    "IMP1_METHOD_LABEL",
    "IMP1_METHOD_STAGE_RECORD_COUNT",
    "IMP1_REJECTION_EVENT_TYPE",
    "IMP1_REJECTION_KEYS",
    "IMP1_RUNTIME_METHODS",
    "IMP1_STAGE_RECORDS_PER_SUBSTEP",
    "IMP1_ZERO_DEBIT",
    "MAX_CHECKPOINT_BYTES",
    "MAX_INTEGER_TEXT",
    "MAX_RATIONAL_BITS",
    "MAX_REJECTION_RECORDS",
    "TDG11IMP1Error",
    "TDG11IMP1Ledger",
    "TDG11IMP1ResourceError",
    "TDG11IMP1TemporalRejection",
    "accept_imp1_step",
    "append_imp1_rejection",
    "combined_diagnostic_debit_vector",
    "decode_inherited_snapshot",
    "encode_inherited_snapshot",
    "imp1_checkpoint_extension",
    "inherited_tdg6_ledger",
    "integer_text",
    "parse_integer_text",
    "restore_imp1_checkpoint_extension",
    "seed_imp1_ledger",
]
