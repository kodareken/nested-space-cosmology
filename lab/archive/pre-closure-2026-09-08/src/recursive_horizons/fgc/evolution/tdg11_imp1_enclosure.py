"""Exact Bernstein sufficient-proof enclosure for prospective TDG11-IMP1.

This module is a new evaluator of the same C1 corrected cubics selected by
MSEL1-PREF1.  It is not a historical MSEL1 wire-equivalent label, not a
scientific runner, and not a global PDE-error or raw-medium difference bound.

Hermite segments are converted directly to Bernstein controls

``(y0, y0 + h f0 / 3, y1 - h f1 / 3, y1)``.

Exact dyadic de Casteljau restriction and subtraction then form the complete
two D01 and four D12 difference cubics per owned row.  Those difference
cubics are bounded: the global lower bound is the exact sample maximum at
``0, 1/4, 1/2, 3/4, 1`` and the upper bound is the convex hull of the
Bernstein controls.  A nonzero cubic cannot vanish at all five samples, so an
interval that merely contains zero cannot be advertised as exact-zero or
enclosure-debit-only.

Adaptive de Casteljau bisection may tighten only certified bounds.  Pieces
are pruned only when their certified upper is at most the current global
lower.  Ties are deterministic.  Every original row and piece remains in the
hashes and polynomial counts.

Sufficient pass is ``8 U12^2 <= L01^2``.  Sufficient failure is
``8 L12^2 > U01^2``.  Zero cases and the TDG6 five-case classifier are
retained.  If the comparison is still unresolved after the bounded
refinement budget, the sealed dual-localizer instrument
``assess_rational_complete_c_rows`` is invoked.  Each agreed exact interval
is intersected with the retained Bernstein interval; an empty intersection
is a typed disagreement.  Classification is derived afresh from the
intersection.  Fallback and exact errors never become an admission.

The public ledger debit is ``Bernstein_U12 + R_fine`` with
``R_fine = max_rows(sum_steps(|delta|))``.  The gate debit is
``gate_U12 + R_fine``.  The nonnegative extra enclosure debit
``Bernstein_U12 - gate_U12`` is recorded separately and is not summed into
the ledger.  Direct Bernstein bounds and the public debit stay in
``(1/3) Z[1/2]`` for binary64-derived records.  The cubic ``6 t^2 - 5 t^3``
has exact maximum ``32/25`` and Bernstein controls ``(0, 0, 2, 1)``: the
odd-prime exact location must not pollute the accumulated public debit.

Raw visibility uses the initial exact Bernstein bounds only and an explicit
unresolved label when needed.  The expensive fallback is never run solely
for the raw path.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
import heapq
import json
from typing import Final, Iterable, Sequence

from .tdg6_temporal_admission_design import (
    CertifiedMagnitudeInterval,
    TDG6_COMPLETE_STATE_CHANNELS,
    TDG6ChannelAdmission,
    classify_tdg6_channel,
)
from .tdg11_msel1_reconstruction import (
    RecordedFamilyError,
    ValidatedRecordedFamily,
    reconstruct_channel,
    validate_recorded_family,
)
from .tdg11_rational_complete_c import (
    RationalCompleteCResourceExhausted,
    RationalCompleteCRouteDisagreement,
    assess_rational_complete_c_rows,
)


Q = Fraction
IMP1_ENCLOSURE_SCHEMA_VERSION: Final[int] = 1
IMP1_ENCLOSURE_EVALUATOR_ID: Final[str] = (
    "tdg11_imp1_exact_bernstein_with_dual_fallback_v1"
)
IMP1_ORDER_SQUARED_MULTIPLIER: Final[int] = 8
_ROW_HASH_DOMAIN: Final[bytes] = b"TDG11-IMP1-BERNSTEIN-ROW-STREAM-v1\n"
_COEFFICIENT_HASH_DOMAIN: Final[bytes] = b"TDG11-IMP1-BERNSTEIN-CUBIC-STREAM-v1\n"
_COMBINED_HASH_DOMAIN: Final[bytes] = b"TDG11-IMP1-BERNSTEIN-COMBINED-v1\n"
_CHANNEL_HASH_DOMAIN: Final[bytes] = b"TDG11-IMP1-CHANNEL-ASSESSMENT-v1\n"
_FAMILY_HASH_DOMAIN: Final[bytes] = b"TDG11-IMP1-FAMILY-ASSESSMENT-v1\n"
_SAMPLE_POINTS: Final[tuple[Fraction, ...]] = (
    Q(0),
    Q(1, 4),
    Q(1, 2),
    Q(3, 4),
    Q(1),
)
_SCOPE_FALSE_FIELDS: Final[tuple[str, ...]] = (
    "absolute_tolerance_used",
    "physical_signal_used_for_normalization",
    "global_PDE_error_certified",
    "raw_medium_difference_bound_certified",
    "production_authority",
)

Cubic = tuple[Fraction, Fraction, Fraction, Fraction]
Segment = tuple[Fraction, Fraction, Fraction, Fraction, Fraction]


def _sha256_digest(value: object, *, name: str) -> str:
    if (
        not isinstance(value, str)
        or len(value) != 64
        or any(character not in "0123456789abcdef" for character in value)
    ):
        raise ValueError(f"{name} must be a lowercase SHA-256 digest")
    return value


def _strict_bool(value: object, *, name: str) -> bool:
    if type(value) is not bool:
        raise TypeError(f"{name} must be a built-in bool")
    return value


def _positive_integer(value: object, *, name: str, allow_zero: bool = False) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{name} must be an integer")
    if value < (0 if allow_zero else 1):
        qualifier = "nonnegative" if allow_zero else "positive"
        raise ValueError(f"{name} must be {qualifier}")
    return value


def _rational(value: object, *, name: str) -> Fraction:
    if type(value) is int:
        return Q(value)
    if type(value) is Fraction:
        if value.denominator <= 0:
            raise ValueError(f"{name} denominator must be positive")
        return value
    raise TypeError(f"{name} must be an exact Fraction or built-in int")


def _exact_nonnegative(value: object, *, name: str) -> Fraction:
    answer = _rational(value, name=name)
    if answer < 0:
        raise ValueError(f"{name} must be nonnegative")
    return answer


def _fixed_sequence(value: object, *, length: int, name: str) -> tuple[object, ...]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise TypeError(f"{name} must be a sequence")
    answer = tuple(value)
    if len(answer) != length:
        raise ValueError(f"{name} must contain exactly {length} entries")
    return answer


_DECIMAL_CHUNK_BASE: Final[int] = 1_000_000_000
_INT_ENCODING_MARKERS: Final[tuple[str, ...]] = (
    "integer string conversion",
    "int_max_str_digits",
)


def _odd_part(denominator: int) -> int:
    if denominator <= 0:
        return denominator
    return denominator >> ((denominator & -denominator).bit_length() - 1)


def _require_direct_ring(
    value: Fraction, *, name: str, allow_three: bool
) -> Fraction:
    odd = _odd_part(value.denominator)
    if allow_three:
        if odd not in (1, 3):
            raise ValueError(
                f"{name} is not in the direct Bernstein ring (1/3)Z[1/2]"
            )
    elif odd != 1:
        raise ValueError(f"{name} is not an exact dyadic rational")
    return value


def _bit_size(value: Fraction) -> int:
    return max(abs(value.numerator).bit_length(), value.denominator.bit_length())


def _int_to_decimal_text(value: int) -> str:
    """Format an int without Python's global decimal-digit conversion limit."""

    if value == 0:
        return "0"
    sign = "-" if value < 0 else ""
    remaining = abs(value)
    chunks: list[int] = []
    base = _DECIMAL_CHUNK_BASE
    while remaining:
        remaining, group = divmod(remaining, base)
        chunks.append(group)
    digits = str(chunks.pop())
    digits += "".join(format(group, "09d") for group in reversed(chunks))
    return sign + digits


def _canonical_text(value: Fraction) -> str:
    return (
        _int_to_decimal_text(value.numerator)
        + "/"
        + _int_to_decimal_text(value.denominator)
    )


def _is_int_encoding_error(exc: BaseException) -> bool:
    if not isinstance(exc, ValueError):
        return False
    text = str(exc)
    return any(marker in text for marker in _INT_ENCODING_MARKERS)


def _canonical_json(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def _digest(domain: bytes, payload: object) -> str:
    return sha256(domain + _canonical_json(payload).encode("ascii")).hexdigest()


class IMP1EnclosureResourceStop(RuntimeError):
    """Hard prospective resource cap; never a scientific nonpass."""

    def __init__(self, reason: str, *, detail: str) -> None:
        if reason not in {"maximum_owned_rows", "maximum_rational_bits"}:
            raise ValueError("unknown IMP1 enclosure resource reason")
        if not isinstance(detail, str) or not detail:
            raise ValueError("resource stop detail must be a nonempty string")
        self.reason = reason
        self.detail = detail
        super().__init__(f"tdg11_imp1_enclosure_resource_stop:{reason}")


class IMP1EnclosureFallbackStop(RuntimeError):
    """Typed fail-closed record for dual-localizer fallback errors."""

    def __init__(self, reason: str, *, detail: str) -> None:
        if reason not in {
            "candidate_ceiling_exhausted",
            "root_isolation_inconclusive",
            "route_disagreement",
            "integer_encoding_exhausted",
        }:
            raise ValueError("unknown IMP1 enclosure fallback reason")
        if not isinstance(detail, str) or not detail:
            raise ValueError("fallback stop detail must be a nonempty string")
        self.reason = reason
        self.detail = detail
        super().__init__(f"tdg11_imp1_enclosure_fallback_stop:{reason}")


class IMP1EnclosureDisagreement(RuntimeError):
    """Bernstein and dual-localizer intervals have empty intersection."""

    def __init__(self, *, level: str, detail: str) -> None:
        if level not in {"D01", "D12"}:
            raise ValueError("disagreement level must be D01 or D12")
        if not isinstance(detail, str) or not detail:
            raise ValueError("disagreement detail must be a nonempty string")
        self.level = level
        self.detail = detail
        super().__init__(f"tdg11_imp1_enclosure_disagreement:{level}")


@dataclass(frozen=True, slots=True)
class IMP1EnclosureLimits:
    """Prospective implementation caps; not scientific thresholds."""

    maximum_owned_rows: int = 2044
    maximum_rational_bits: int = 32768
    maximum_subdivision_depth: int = 12
    maximum_subdivisions_per_pair: int = 32768
    fallback_maximum_candidates_D01: int = 16352
    fallback_maximum_candidates_D12: int = 32704
    fallback_refinement_depth: int = 160

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "maximum_owned_rows",
            _positive_integer(self.maximum_owned_rows, name="maximum_owned_rows"),
        )
        object.__setattr__(
            self,
            "maximum_rational_bits",
            _positive_integer(
                self.maximum_rational_bits, name="maximum_rational_bits"
            ),
        )
        object.__setattr__(
            self,
            "maximum_subdivision_depth",
            _positive_integer(
                self.maximum_subdivision_depth,
                name="maximum_subdivision_depth",
                allow_zero=True,
            ),
        )
        object.__setattr__(
            self,
            "maximum_subdivisions_per_pair",
            _positive_integer(
                self.maximum_subdivisions_per_pair,
                name="maximum_subdivisions_per_pair",
                allow_zero=True,
            ),
        )
        object.__setattr__(
            self,
            "fallback_maximum_candidates_D01",
            _positive_integer(
                self.fallback_maximum_candidates_D01,
                name="fallback_maximum_candidates_D01",
            ),
        )
        object.__setattr__(
            self,
            "fallback_maximum_candidates_D12",
            _positive_integer(
                self.fallback_maximum_candidates_D12,
                name="fallback_maximum_candidates_D12",
            ),
        )
        object.__setattr__(
            self,
            "fallback_refinement_depth",
            _positive_integer(
                self.fallback_refinement_depth,
                name="fallback_refinement_depth",
                allow_zero=True,
            ),
        )

    def as_mapping(self) -> dict[str, int]:
        return {
            "maximum_owned_rows": self.maximum_owned_rows,
            "maximum_rational_bits": self.maximum_rational_bits,
            "maximum_subdivision_depth": self.maximum_subdivision_depth,
            "maximum_subdivisions_per_pair": self.maximum_subdivisions_per_pair,
            "fallback_maximum_candidates_D01": self.fallback_maximum_candidates_D01,
            "fallback_maximum_candidates_D12": self.fallback_maximum_candidates_D12,
            "fallback_refinement_depth": self.fallback_refinement_depth,
        }


DEFAULT_IMP1_ENCLOSURE_LIMITS: Final[IMP1EnclosureLimits] = IMP1EnclosureLimits()


def _require_limits(value: object) -> IMP1EnclosureLimits:
    if not isinstance(value, IMP1EnclosureLimits):
        raise TypeError("limits must be IMP1EnclosureLimits")
    return IMP1EnclosureLimits(
        maximum_owned_rows=value.maximum_owned_rows,
        maximum_rational_bits=value.maximum_rational_bits,
        maximum_subdivision_depth=value.maximum_subdivision_depth,
        maximum_subdivisions_per_pair=value.maximum_subdivisions_per_pair,
        fallback_maximum_candidates_D01=value.fallback_maximum_candidates_D01,
        fallback_maximum_candidates_D12=value.fallback_maximum_candidates_D12,
        fallback_refinement_depth=value.fallback_refinement_depth,
    )


def _check_bits(
    value: Fraction, *, limits: IMP1EnclosureLimits, name: str
) -> Fraction:
    if _bit_size(value) > limits.maximum_rational_bits:
        raise IMP1EnclosureResourceStop("maximum_rational_bits", detail=name)
    return value


def _admit_rational(
    value: object,
    *,
    name: str,
    limits: IMP1EnclosureLimits,
    allow_three: bool,
    nonnegative: bool = False,
) -> Fraction:
    exact = _rational(value, name=name)
    if nonnegative and exact < 0:
        raise ValueError(f"{name} must be nonnegative")
    _check_bits(exact, limits=limits, name=name)
    return _require_direct_ring(exact, name=name, allow_three=allow_three)


def _cubic(value: object, *, name: str) -> Cubic:
    raw = _fixed_sequence(value, length=4, name=name)
    return tuple(_rational(item, name=f"{name}[{index}]") for index, item in enumerate(raw))  # type: ignore[return-value]


def _parse_segment(
    value: object, *, name: str, limits: IMP1EnclosureLimits
) -> Segment:
    raw = _fixed_sequence(value, length=5, name=name)
    exact = tuple(
        _admit_rational(
            item,
            name=f"{name}[{index}]",
            limits=limits,
            allow_three=index in (0, 2),
        )
        for index, item in enumerate(raw)
    )
    if exact[4] <= 0:
        raise ValueError(f"{name} width must be positive")
    return exact  # type: ignore[return-value]


def _hermite_bernstein(segment: Segment) -> Cubic:
    y0, f0, y1, f1, width = segment
    return (y0, y0 + width * f0 / 3, y1 - width * f1 / 3, y1)


def _split_bernstein(controls: Cubic) -> tuple[Cubic, Cubic]:
    b0, b1, b2, b3 = controls
    c01 = (b0 + b1) / 2
    c12 = (b1 + b2) / 2
    c23 = (b2 + b3) / 2
    c012 = (c01 + c12) / 2
    c123 = (c12 + c23) / 2
    midpoint = (c012 + c123) / 2
    return (b0, c01, c012, midpoint), (midpoint, c123, c23, b3)


def _subtract_bernstein(left: Cubic, right: Cubic) -> Cubic:
    return (left[0] - right[0], left[1] - right[1], left[2] - right[2], left[3] - right[3])


def _bernstein_value(controls: Cubic, point: Fraction) -> Fraction:
    b0, b1, b2, b3 = controls
    one_minus = 1 - point
    return (
        one_minus**3 * b0
        + 3 * one_minus**2 * point * b1
        + 3 * one_minus * point**2 * b2
        + point**3 * b3
    )


def _identically_zero(controls: Cubic) -> bool:
    return controls[0] == 0 and controls[1] == 0 and controls[2] == 0 and controls[3] == 0


def _hull_upper(controls: Cubic) -> Fraction:
    return max(abs(item) for item in controls)


def _sample_lower(controls: Cubic) -> Fraction:
    samples = tuple(_bernstein_value(controls, point) for point in _SAMPLE_POINTS)
    if _identically_zero(controls):
        if any(item != 0 for item in samples):
            raise ValueError("zero Bernstein cubic produced a nonzero sample")
        return Q(0)
    magnitude = max(abs(item) for item in samples)
    if magnitude == 0:
        raise ValueError("nonzero cubic vanished at all five dyadic samples")
    return magnitude


def _certify_controls(
    controls: Cubic, *, limits: IMP1EnclosureLimits, name: str
) -> Cubic:
    certified: list[Fraction] = []
    for index, item in enumerate(controls):
        labeled = f"{name}[{index}]"
        certified.append(
            _admit_rational(
                item,
                name=labeled,
                limits=limits,
                allow_three=True,
            )
        )
    return tuple(certified)  # type: ignore[return-value]


def _certify_bound(
    value: Fraction, *, limits: IMP1EnclosureLimits, name: str
) -> Fraction:
    return _admit_rational(
        value,
        name=name,
        limits=limits,
        allow_three=True,
        nonnegative=True,
    )


def _feed_coefficient_hash(hasher: object, ordinal: int, controls: Cubic) -> None:
    line = (
        f"{ordinal}|"
        + "|".join(_canonical_text(item) for item in controls)
        + "\n"
    )
    hasher.update(line.encode("ascii"))  # type: ignore[attr-defined]


def _combined_coefficient_hash(outer_digest: str, finest_digest: str) -> str:
    return sha256(
        _COMBINED_HASH_DOMAIN + f"{outer_digest}\n{finest_digest}\n".encode("ascii")
    ).hexdigest()


def _row_line(
    ordinal: int,
    outer: Segment,
    medium: Sequence[Segment],
    fine: Sequence[Segment],
) -> str:
    values: list[Fraction] = []
    for segment in (outer, *medium, *fine):
        values.extend(segment)
    encoded = "|".join(_canonical_text(item) for item in values)
    return f"{ordinal}|{encoded}\n"


def _contraction_flags(
    *,
    outer_lower: Fraction,
    outer_upper: Fraction,
    finest_lower: Fraction,
    finest_upper: Fraction,
    limits: IMP1EnclosureLimits | None = None,
) -> tuple[Fraction, Fraction, bool, bool, bool, bool]:
    left = IMP1_ORDER_SQUARED_MULTIPLIER * finest_upper * finest_upper
    right = outer_lower * outer_lower
    if limits is not None:
        _check_bits(left, limits=limits, name="sufficient_pass_left")
        _check_bits(right, limits=limits, name="sufficient_pass_right")
    exact_zero = outer_upper == 0 and finest_upper == 0
    if exact_zero:
        return left, right, False, False, False, True
    if left <= right:
        return left, right, True, False, False, False
    failure_left = IMP1_ORDER_SQUARED_MULTIPLIER * finest_lower * finest_lower
    failure_right = outer_upper * outer_upper
    if limits is not None:
        _check_bits(failure_left, limits=limits, name="sufficient_failure_left")
        _check_bits(failure_right, limits=limits, name="sufficient_failure_right")
    if failure_left > failure_right:
        return left, right, False, True, False, False
    return left, right, False, False, True, False


def _interval_mapping(lower: Fraction, upper: Fraction) -> dict[str, str]:
    return {"lower": _canonical_text(lower), "upper": _canonical_text(upper)}


def _decision_mapping(decision: TDG6ChannelAdmission) -> dict[str, object]:
    order_passed = decision.order_threshold_passed
    return {
        "classification": decision.classification,
        "admission_passed": decision.admission_passed,
        "temporal_retry_permitted": decision.temporal_retry_permitted,
        "minimum_observed_order": _canonical_text(decision.minimum_observed_order),
        "outer_difference": _interval_mapping(
            decision.outer_difference.lower, decision.outer_difference.upper
        ),
        "finest_difference": _interval_mapping(
            decision.finest_difference.lower, decision.finest_difference.upper
        ),
        "finest_pair_debit": _canonical_text(decision.finest_pair_debit),
        "order_threshold_resolved": decision.order_threshold_resolved,
        "order_threshold_passed": order_passed,
    }


@dataclass(frozen=True, slots=True)
class IMP1CubicFamilyBounds:
    """Exact sample/hull enclosure of one Bernstein cubic family."""

    lower: Fraction
    upper: Fraction
    polynomial_count: int
    coefficient_stream_sha256: str
    subdivisions_used: int
    maximum_depth_reached: int
    identically_zero: bool

    def __post_init__(self) -> None:
        lower = _exact_nonnegative(self.lower, name="lower")
        upper = _exact_nonnegative(self.upper, name="upper")
        if lower > upper:
            raise ValueError("cubic-family lower bound exceeds upper bound")
        object.__setattr__(self, "lower", lower)
        object.__setattr__(self, "upper", upper)
        _positive_integer(self.polynomial_count, name="polynomial_count")
        _positive_integer(
            self.subdivisions_used, name="subdivisions_used", allow_zero=True
        )
        _positive_integer(
            self.maximum_depth_reached,
            name="maximum_depth_reached",
            allow_zero=True,
        )
        _sha256_digest(
            self.coefficient_stream_sha256, name="coefficient stream hash"
        )
        identically_zero = _strict_bool(
            self.identically_zero, name="identically_zero"
        )
        if identically_zero is not (lower == 0 and upper == 0):
            raise ValueError("identically-zero flag differs from the certified bounds")


@dataclass(frozen=True, slots=True)
class IMP1PairBounds:
    """Retained Bernstein D01/D12 enclosure, flags, and stream provenance."""

    d01_lower: Fraction
    d01_upper: Fraction
    d12_lower: Fraction
    d12_upper: Fraction
    row_count: int
    polynomial_count_d01: int
    polynomial_count_d12: int
    row_stream_sha256: str
    d01_coefficient_stream_sha256: str
    d12_coefficient_stream_sha256: str
    combined_coefficient_stream_sha256: str
    subdivisions_used: int
    maximum_depth_reached: int
    sufficient_pass: bool
    sufficient_failure: bool
    threshold_inconclusive: bool
    exact_zero: bool

    def __post_init__(self) -> None:
        d01_lower = _exact_nonnegative(self.d01_lower, name="d01_lower")
        d01_upper = _exact_nonnegative(self.d01_upper, name="d01_upper")
        d12_lower = _exact_nonnegative(self.d12_lower, name="d12_lower")
        d12_upper = _exact_nonnegative(self.d12_upper, name="d12_upper")
        if d01_lower > d01_upper or d12_lower > d12_upper:
            raise ValueError("difference lower bound exceeds upper bound")
        object.__setattr__(self, "d01_lower", d01_lower)
        object.__setattr__(self, "d01_upper", d01_upper)
        object.__setattr__(self, "d12_lower", d12_lower)
        object.__setattr__(self, "d12_upper", d12_upper)
        row_count = _positive_integer(self.row_count, name="row_count")
        if self.polynomial_count_d01 != 2 * row_count:
            raise ValueError("D01 polynomial count must equal twice row_count")
        if self.polynomial_count_d12 != 4 * row_count:
            raise ValueError("D12 polynomial count must equal four times row_count")
        _positive_integer(
            self.subdivisions_used, name="subdivisions_used", allow_zero=True
        )
        _positive_integer(
            self.maximum_depth_reached,
            name="maximum_depth_reached",
            allow_zero=True,
        )
        _sha256_digest(self.row_stream_sha256, name="row stream hash")
        d01_hash = _sha256_digest(
            self.d01_coefficient_stream_sha256, name="D01 coefficient hash"
        )
        d12_hash = _sha256_digest(
            self.d12_coefficient_stream_sha256, name="D12 coefficient hash"
        )
        combined = _sha256_digest(
            self.combined_coefficient_stream_sha256, name="combined coefficient hash"
        )
        if combined != _combined_coefficient_hash(d01_hash, d12_hash):
            raise ValueError("combined coefficient hash differs from component hashes")
        (
            _left,
            _right,
            expected_pass,
            expected_failure,
            expected_inconclusive,
            expected_zero,
        ) = _contraction_flags(
            outer_lower=d01_lower,
            outer_upper=d01_upper,
            finest_lower=d12_lower,
            finest_upper=d12_upper,
        )
        if _strict_bool(self.sufficient_pass, name="sufficient_pass") is not expected_pass:
            raise ValueError("sufficient-pass flag differs from the certified bounds")
        if (
            _strict_bool(self.sufficient_failure, name="sufficient_failure")
            is not expected_failure
        ):
            raise ValueError("sufficient-failure flag differs from the certified bounds")
        if (
            _strict_bool(self.threshold_inconclusive, name="threshold_inconclusive")
            is not expected_inconclusive
        ):
            raise ValueError(
                "threshold-inconclusive flag differs from the certified bounds"
            )
        if _strict_bool(self.exact_zero, name="exact_zero") is not expected_zero:
            raise ValueError("exact-zero flag differs from the certified bounds")

    def as_mapping(self) -> dict[str, object]:
        return {
            "d01_lower": _canonical_text(self.d01_lower),
            "d01_upper": _canonical_text(self.d01_upper),
            "d12_lower": _canonical_text(self.d12_lower),
            "d12_upper": _canonical_text(self.d12_upper),
            "row_count": self.row_count,
            "polynomial_count_d01": self.polynomial_count_d01,
            "polynomial_count_d12": self.polynomial_count_d12,
            "row_stream_sha256": self.row_stream_sha256,
            "d01_coefficient_stream_sha256": self.d01_coefficient_stream_sha256,
            "d12_coefficient_stream_sha256": self.d12_coefficient_stream_sha256,
            "combined_coefficient_stream_sha256": self.combined_coefficient_stream_sha256,
            "subdivisions_used": self.subdivisions_used,
            "maximum_depth_reached": self.maximum_depth_reached,
            "sufficient_pass": self.sufficient_pass,
            "sufficient_failure": self.sufficient_failure,
            "threshold_inconclusive": self.threshold_inconclusive,
            "exact_zero": self.exact_zero,
        }


@dataclass(frozen=True, slots=True)
class IMP1PathAssessment:
    """One raw or corrected path: retained Bernstein bounds, gate, and debit."""

    bernstein: IMP1PairBounds
    gate_d01_lower: Fraction
    gate_d01_upper: Fraction
    gate_d12_lower: Fraction
    gate_d12_upper: Fraction
    decision: TDG6ChannelAdmission
    sufficient_pass_left: Fraction
    sufficient_pass_right: Fraction
    sufficient_contraction_pass: bool
    sufficient_contraction_failure: bool
    threshold_inconclusive: bool
    exact_zero: bool
    accumulation_bounds: tuple[Fraction, Fraction, Fraction]
    public_fine_debit: Fraction
    gate_debit: Fraction
    extra_enclosure_debit: Fraction
    limits: IMP1EnclosureLimits
    fallback_called: bool
    fallback_calls: int
    fallback_d01_lower: Fraction | None = None
    fallback_d01_upper: Fraction | None = None
    fallback_d12_lower: Fraction | None = None
    fallback_d12_upper: Fraction | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.bernstein, IMP1PairBounds):
            raise TypeError("bernstein must be IMP1PairBounds")
        if not isinstance(self.decision, TDG6ChannelAdmission):
            raise TypeError("decision must be TDG6ChannelAdmission")
        gate_d01_lower = _exact_nonnegative(self.gate_d01_lower, name="gate_d01_lower")
        gate_d01_upper = _exact_nonnegative(self.gate_d01_upper, name="gate_d01_upper")
        gate_d12_lower = _exact_nonnegative(self.gate_d12_lower, name="gate_d12_lower")
        gate_d12_upper = _exact_nonnegative(self.gate_d12_upper, name="gate_d12_upper")
        if gate_d01_lower > gate_d01_upper or gate_d12_lower > gate_d12_upper:
            raise ValueError("gate lower bound exceeds upper bound")
        if (
            gate_d01_lower < self.bernstein.d01_lower
            or gate_d01_upper > self.bernstein.d01_upper
            or gate_d12_lower < self.bernstein.d12_lower
            or gate_d12_upper > self.bernstein.d12_upper
        ):
            raise ValueError("gate interval is not inside the retained Bernstein interval")
        object.__setattr__(self, "gate_d01_lower", gate_d01_lower)
        object.__setattr__(self, "gate_d01_upper", gate_d01_upper)
        object.__setattr__(self, "gate_d12_lower", gate_d12_lower)
        object.__setattr__(self, "gate_d12_upper", gate_d12_upper)
        accumulation = _fixed_sequence(
            self.accumulation_bounds, length=3, name="accumulation_bounds"
        )
        exact_accumulation = tuple(
            _exact_nonnegative(item, name=f"accumulation_bounds[{index}]")
            for index, item in enumerate(accumulation)
        )
        object.__setattr__(self, "accumulation_bounds", exact_accumulation)
        public = _exact_nonnegative(self.public_fine_debit, name="public_fine_debit")
        gate_debit = _exact_nonnegative(self.gate_debit, name="gate_debit")
        extra = _exact_nonnegative(
            self.extra_enclosure_debit, name="extra_enclosure_debit"
        )
        object.__setattr__(self, "public_fine_debit", public)
        object.__setattr__(self, "gate_debit", gate_debit)
        object.__setattr__(self, "extra_enclosure_debit", extra)
        r_fine = exact_accumulation[2]
        if public != self.bernstein.d12_upper + r_fine:
            raise ValueError("public fine debit differs from Bernstein_U12 + R_fine")
        if gate_debit != gate_d12_upper + r_fine:
            raise ValueError("gate debit differs from gate_U12 + R_fine")
        if extra != self.bernstein.d12_upper - gate_d12_upper:
            raise ValueError("extra enclosure debit differs from Bernstein_U12 - gate_U12")
        limits = _require_limits(self.limits)
        object.__setattr__(self, "limits", limits)
        for labeled, item in (
            ("gate_d01_lower", gate_d01_lower),
            ("gate_d01_upper", gate_d01_upper),
            ("gate_d12_lower", gate_d12_lower),
            ("gate_d12_upper", gate_d12_upper),
            ("public_fine_debit", public),
            ("gate_debit", gate_debit),
            ("extra_enclosure_debit", extra),
        ):
            _check_bits(item, limits=limits, name=labeled)
        for index, item in enumerate(exact_accumulation):
            _check_bits(item, limits=limits, name=f"accumulation_bounds[{index}]")
        _require_direct_ring(public, name="public_fine_debit", allow_three=True)
        for index, item in enumerate(exact_accumulation):
            _require_direct_ring(
                item, name=f"accumulation_bounds[{index}]", allow_three=True
            )
        left = _exact_nonnegative(self.sufficient_pass_left, name="sufficient_pass_left")
        right = _exact_nonnegative(
            self.sufficient_pass_right, name="sufficient_pass_right"
        )
        object.__setattr__(self, "sufficient_pass_left", left)
        object.__setattr__(self, "sufficient_pass_right", right)
        expected_left = IMP1_ORDER_SQUARED_MULTIPLIER * gate_d12_upper * gate_d12_upper
        expected_right = gate_d01_lower * gate_d01_lower
        _check_bits(expected_left, limits=limits, name="sufficient_pass_left")
        _check_bits(expected_right, limits=limits, name="sufficient_pass_right")
        _check_bits(left, limits=limits, name="sufficient_pass_left")
        _check_bits(right, limits=limits, name="sufficient_pass_right")
        if left != expected_left:
            raise ValueError("sufficient-pass left side differs from 8 U12^2")
        if right != expected_right:
            raise ValueError("sufficient-pass right side differs from L01^2")
        (
            _left,
            _right,
            expected_pass,
            expected_failure,
            expected_inconclusive,
            expected_zero,
        ) = _contraction_flags(
            outer_lower=gate_d01_lower,
            outer_upper=gate_d01_upper,
            finest_lower=gate_d12_lower,
            finest_upper=gate_d12_upper,
        )
        if (
            _strict_bool(
                self.sufficient_contraction_pass, name="sufficient_contraction_pass"
            )
            is not expected_pass
        ):
            raise ValueError("sufficient-pass flag differs from the gate bounds")
        if (
            _strict_bool(
                self.sufficient_contraction_failure,
                name="sufficient_contraction_failure",
            )
            is not expected_failure
        ):
            raise ValueError("sufficient-failure flag differs from the gate bounds")
        if (
            _strict_bool(self.threshold_inconclusive, name="threshold_inconclusive")
            is not expected_inconclusive
        ):
            raise ValueError(
                "threshold-inconclusive flag differs from the gate bounds"
            )
        if _strict_bool(self.exact_zero, name="exact_zero") is not expected_zero:
            raise ValueError("exact-zero flag differs from the gate bounds")
        expected_decision = classify_tdg6_channel(
            CertifiedMagnitudeInterval(gate_d01_lower, gate_d01_upper),
            CertifiedMagnitudeInterval(gate_d12_lower, gate_d12_upper),
        )
        if self.decision != expected_decision:
            raise ValueError(
                "channel decision differs from classify_tdg6_channel on the gate"
            )
        fallback_called = _strict_bool(self.fallback_called, name="fallback_called")
        fallback_calls = _positive_integer(
            self.fallback_calls, name="fallback_calls", allow_zero=True
        )
        fallback_fields = (
            self.fallback_d01_lower,
            self.fallback_d01_upper,
            self.fallback_d12_lower,
            self.fallback_d12_upper,
        )
        if fallback_called:
            if fallback_calls != 1:
                raise ValueError("fallback call count must be one when fallback ran")
            if any(item is None for item in fallback_fields):
                raise ValueError("fallback bounds are required when fallback ran")
            fallback_d01_lower = _exact_nonnegative(
                self.fallback_d01_lower, name="fallback_d01_lower"
            )
            fallback_d01_upper = _exact_nonnegative(
                self.fallback_d01_upper, name="fallback_d01_upper"
            )
            fallback_d12_lower = _exact_nonnegative(
                self.fallback_d12_lower, name="fallback_d12_lower"
            )
            fallback_d12_upper = _exact_nonnegative(
                self.fallback_d12_upper, name="fallback_d12_upper"
            )
            if (
                fallback_d01_lower > fallback_d01_upper
                or fallback_d12_lower > fallback_d12_upper
            ):
                raise ValueError("fallback lower bound exceeds upper bound")
            for labeled, item in (
                ("fallback_d01_lower", fallback_d01_lower),
                ("fallback_d01_upper", fallback_d01_upper),
                ("fallback_d12_lower", fallback_d12_lower),
                ("fallback_d12_upper", fallback_d12_upper),
            ):
                _check_bits(item, limits=limits, name=labeled)
            object.__setattr__(self, "fallback_d01_lower", fallback_d01_lower)
            object.__setattr__(self, "fallback_d01_upper", fallback_d01_upper)
            object.__setattr__(self, "fallback_d12_lower", fallback_d12_lower)
            object.__setattr__(self, "fallback_d12_upper", fallback_d12_upper)
            _require_gate_equals_intersection(
                bernstein=self.bernstein,
                gate_d01=(gate_d01_lower, gate_d01_upper),
                gate_d12=(gate_d12_lower, gate_d12_upper),
                fallback_d01=(fallback_d01_lower, fallback_d01_upper),
                fallback_d12=(fallback_d12_lower, fallback_d12_upper),
            )
        else:
            if fallback_calls != 0:
                raise ValueError("fast-only decisions must record zero fallback calls")
            if any(item is not None for item in fallback_fields):
                raise ValueError("fast-only decisions must not invent fallback bounds")
            if extra != 0:
                raise ValueError("fast-only extra enclosure debit must be zero")
            if (
                gate_d01_lower != self.bernstein.d01_lower
                or gate_d01_upper != self.bernstein.d01_upper
                or gate_d12_lower != self.bernstein.d12_lower
                or gate_d12_upper != self.bernstein.d12_upper
            ):
                raise ValueError("fast-only gate must equal the Bernstein enclosure")

    def as_mapping(self) -> dict[str, object]:
        if self.fallback_called:
            _require_gate_equals_intersection(
                bernstein=self.bernstein,
                gate_d01=(self.gate_d01_lower, self.gate_d01_upper),
                gate_d12=(self.gate_d12_lower, self.gate_d12_upper),
                fallback_d01=(self.fallback_d01_lower, self.fallback_d01_upper),
                fallback_d12=(self.fallback_d12_lower, self.fallback_d12_upper),
            )
        mapping: dict[str, object] = {
            "bernstein": self.bernstein.as_mapping(),
            "gate": {
                "d01": _interval_mapping(self.gate_d01_lower, self.gate_d01_upper),
                "d12": _interval_mapping(self.gate_d12_lower, self.gate_d12_upper),
            },
            "decision": _decision_mapping(self.decision),
            "sufficient_pass_left": _canonical_text(self.sufficient_pass_left),
            "sufficient_pass_right": _canonical_text(self.sufficient_pass_right),
            "sufficient_contraction_pass": self.sufficient_contraction_pass,
            "sufficient_contraction_failure": self.sufficient_contraction_failure,
            "threshold_inconclusive": self.threshold_inconclusive,
            "exact_zero": self.exact_zero,
            "accumulation_bounds": [
                _canonical_text(item) for item in self.accumulation_bounds
            ],
            "public_fine_debit": _canonical_text(self.public_fine_debit),
            "gate_debit": _canonical_text(self.gate_debit),
            "extra_enclosure_debit": _canonical_text(self.extra_enclosure_debit),
            "limits": self.limits.as_mapping(),
            "fallback_called": self.fallback_called,
            "fallback_calls": self.fallback_calls,
        }
        if self.fallback_called:
            mapping["fallback"] = {
                "d01": _interval_mapping(
                    self.fallback_d01_lower, self.fallback_d01_upper  # type: ignore[arg-type]
                ),
                "d12": _interval_mapping(
                    self.fallback_d12_lower, self.fallback_d12_upper  # type: ignore[arg-type]
                ),
            }
        return mapping


def _require_scope_false(**flags: object) -> None:
    for name, value in flags.items():
        if name not in _SCOPE_FALSE_FIELDS or value is not False:
            raise ValueError("IMP1 enclosure evidence crossed its scope")


@dataclass(frozen=True, slots=True)
class IMP1ChannelAssessment:
    """Raw visibility plus corrected Bernstein/gate admission for one channel."""

    channel: str
    method: str
    family_sha256: str
    row_count: int
    limits: IMP1EnclosureLimits
    raw: IMP1PairBounds
    raw_unresolved: bool
    corrected: IMP1PathAssessment
    evaluator_id: str = IMP1_ENCLOSURE_EVALUATOR_ID
    schema_version: int = IMP1_ENCLOSURE_SCHEMA_VERSION
    absolute_tolerance_used: bool = False
    physical_signal_used_for_normalization: bool = False
    global_PDE_error_certified: bool = False
    raw_medium_difference_bound_certified: bool = False
    production_authority: bool = False
    channel_sha256: str = ""

    def __post_init__(self) -> None:
        if self.channel not in TDG6_COMPLETE_STATE_CHANNELS:
            raise ValueError("channel is not a TDG6 complete-state channel")
        if not isinstance(self.method, str) or not self.method:
            raise ValueError("method must be a nonempty string")
        _sha256_digest(self.family_sha256, name="family hash")
        _positive_integer(self.row_count, name="row_count")
        limits = _require_limits(self.limits)
        object.__setattr__(self, "limits", limits)
        if not isinstance(self.raw, IMP1PairBounds):
            raise TypeError("raw must be IMP1PairBounds")
        if not isinstance(self.corrected, IMP1PathAssessment):
            raise TypeError("corrected must be IMP1PathAssessment")
        if self.raw.row_count != self.row_count or self.corrected.bernstein.row_count != self.row_count:
            raise ValueError("raw and corrected row counts must match the channel")
        if self.corrected.limits != limits:
            raise ValueError("channel limits differ from the path limits")
        if self.raw.subdivisions_used != 0 or self.raw.maximum_depth_reached != 0:
            raise ValueError("raw visibility must retain initial Bernstein bounds only")
        for labeled, item in (
            ("raw.d01_lower", self.raw.d01_lower),
            ("raw.d01_upper", self.raw.d01_upper),
            ("raw.d12_lower", self.raw.d12_lower),
            ("raw.d12_upper", self.raw.d12_upper),
        ):
            _check_bits(item, limits=limits, name=labeled)
        raw_unresolved = _strict_bool(self.raw_unresolved, name="raw_unresolved")
        if raw_unresolved is not self.raw.threshold_inconclusive:
            raise ValueError("raw unresolved label differs from the Bernstein flags")
        if (
            type(self.schema_version) is not int
            or self.schema_version != IMP1_ENCLOSURE_SCHEMA_VERSION
            or self.evaluator_id != IMP1_ENCLOSURE_EVALUATOR_ID
        ):
            raise ValueError("IMP1 enclosure evidence crossed its scope")
        _require_scope_false(
            absolute_tolerance_used=self.absolute_tolerance_used,
            physical_signal_used_for_normalization=self.physical_signal_used_for_normalization,
            global_PDE_error_certified=self.global_PDE_error_certified,
            raw_medium_difference_bound_certified=self.raw_medium_difference_bound_certified,
            production_authority=self.production_authority,
        )
        expected = _digest(_CHANNEL_HASH_DOMAIN, self._payload())
        digest = self.channel_sha256
        if digest == "":
            object.__setattr__(self, "channel_sha256", expected)
        else:
            if digest != expected:
                raise ValueError("channel assessment hash does not match the payload")

    def _payload(self) -> dict[str, object]:
        return {
            "channel": self.channel,
            "method": self.method,
            "family_sha256": self.family_sha256,
            "row_count": self.row_count,
            "limits": self.limits.as_mapping(),
            "raw": self.raw.as_mapping(),
            "raw_unresolved": self.raw_unresolved,
            "corrected": self.corrected.as_mapping(),
            "evaluator_id": self.evaluator_id,
            "schema_version": self.schema_version,
            "absolute_tolerance_used": self.absolute_tolerance_used,
            "physical_signal_used_for_normalization": (
                self.physical_signal_used_for_normalization
            ),
            "global_PDE_error_certified": self.global_PDE_error_certified,
            "raw_medium_difference_bound_certified": (
                self.raw_medium_difference_bound_certified
            ),
            "production_authority": self.production_authority,
        }

    def as_mapping(self) -> dict[str, object]:
        payload = self._payload()
        payload["channel_sha256"] = self.channel_sha256
        payload["admission_passed"] = self.corrected.decision.admission_passed
        payload["public_fine_debit"] = _canonical_text(self.corrected.public_fine_debit)
        payload["gate_debit"] = _canonical_text(self.corrected.gate_debit)
        payload["extra_enclosure_debit"] = _canonical_text(
            self.corrected.extra_enclosure_debit
        )
        payload["accumulation_bounds"] = [
            _canonical_text(item) for item in self.corrected.accumulation_bounds
        ]
        return payload


@dataclass(frozen=True, slots=True)
class IMP1FamilyAssessment:
    """All eighteen complete-state channel records for one recorded family."""

    method: str
    family_sha256: str
    row_count: int
    limits: IMP1EnclosureLimits
    channels: tuple[IMP1ChannelAssessment, ...]
    admission_passed: bool
    public_fine_debits: tuple[Fraction, ...]
    evaluator_id: str = IMP1_ENCLOSURE_EVALUATOR_ID
    schema_version: int = IMP1_ENCLOSURE_SCHEMA_VERSION
    absolute_tolerance_used: bool = False
    physical_signal_used_for_normalization: bool = False
    global_PDE_error_certified: bool = False
    raw_medium_difference_bound_certified: bool = False
    production_authority: bool = False
    assessment_sha256: str = ""

    def __post_init__(self) -> None:
        if not isinstance(self.method, str) or not self.method:
            raise ValueError("method must be a nonempty string")
        _sha256_digest(self.family_sha256, name="family hash")
        _positive_integer(self.row_count, name="row_count")
        limits = _require_limits(self.limits)
        object.__setattr__(self, "limits", limits)
        if not isinstance(self.channels, tuple):
            raise TypeError("channels must be a tuple")
        if len(self.channels) != len(TDG6_COMPLETE_STATE_CHANNELS):
            raise ValueError("family assessment must contain exactly 18 channel records")
        names: list[str] = []
        expected_debits: list[Fraction] = []
        for index, record in enumerate(self.channels):
            if not isinstance(record, IMP1ChannelAssessment):
                raise TypeError("channel records must be IMP1ChannelAssessment")
            names.append(record.channel)
            if record.method != self.method:
                raise ValueError("channel method differs from the family")
            if record.family_sha256 != self.family_sha256:
                raise ValueError("channel family hash differs from the family")
            if record.row_count != self.row_count:
                raise ValueError("channel row count differs from the family")
            if record.limits != limits:
                raise ValueError("channel limits differ from the family")
            if record.evaluator_id != IMP1_ENCLOSURE_EVALUATOR_ID:
                raise ValueError("channel evaluator differs from IMP1")
            expected_debits.append(record.corrected.public_fine_debit)
        if tuple(names) != TDG6_COMPLETE_STATE_CHANNELS:
            raise ValueError("family channels are not the TDG6 complete-state order")
        if not isinstance(self.public_fine_debits, tuple):
            raise TypeError("public_fine_debits must be a tuple")
        if len(self.public_fine_debits) != len(TDG6_COMPLETE_STATE_CHANNELS):
            raise ValueError("public_fine_debits must contain exactly 18 entries")
        exact_debits = tuple(
            _exact_nonnegative(item, name=f"public_fine_debits[{index}]")
            for index, item in enumerate(self.public_fine_debits)
        )
        object.__setattr__(self, "public_fine_debits", exact_debits)
        if exact_debits != tuple(expected_debits):
            raise ValueError("public_fine_debits differ from the channel records")
        expected_passed = all(
            record.corrected.decision.admission_passed for record in self.channels
        )
        if _strict_bool(self.admission_passed, name="admission_passed") is not expected_passed:
            raise ValueError("family admission flag differs from the channel records")
        if (
            type(self.schema_version) is not int
            or self.schema_version != IMP1_ENCLOSURE_SCHEMA_VERSION
            or self.evaluator_id != IMP1_ENCLOSURE_EVALUATOR_ID
        ):
            raise ValueError("IMP1 enclosure evidence crossed its scope")
        _require_scope_false(
            absolute_tolerance_used=self.absolute_tolerance_used,
            physical_signal_used_for_normalization=self.physical_signal_used_for_normalization,
            global_PDE_error_certified=self.global_PDE_error_certified,
            raw_medium_difference_bound_certified=self.raw_medium_difference_bound_certified,
            production_authority=self.production_authority,
        )
        expected = _digest(_FAMILY_HASH_DOMAIN, self._payload())
        digest = self.assessment_sha256
        if digest == "":
            object.__setattr__(self, "assessment_sha256", expected)
        elif digest != expected:
            raise ValueError("family assessment hash does not match the payload")

    def _payload(self) -> dict[str, object]:
        return {
            "method": self.method,
            "family_sha256": self.family_sha256,
            "row_count": self.row_count,
            "limits": self.limits.as_mapping(),
            "channels": [record.as_mapping() for record in self.channels],
            "admission_passed": self.admission_passed,
            "public_fine_debits": [
                _canonical_text(item) for item in self.public_fine_debits
            ],
            "evaluator_id": self.evaluator_id,
            "schema_version": self.schema_version,
            "absolute_tolerance_used": self.absolute_tolerance_used,
            "physical_signal_used_for_normalization": (
                self.physical_signal_used_for_normalization
            ),
            "global_PDE_error_certified": self.global_PDE_error_certified,
            "raw_medium_difference_bound_certified": (
                self.raw_medium_difference_bound_certified
            ),
            "production_authority": self.production_authority,
        }

    def as_mapping(self) -> dict[str, object]:
        payload = self._payload()
        payload["assessment_sha256"] = self.assessment_sha256
        return payload


def _piece_bounds(
    controls: Cubic, *, limits: IMP1EnclosureLimits, name: str
) -> tuple[Fraction, Fraction]:
    certified = _certify_controls(controls, limits=limits, name=name)
    lower = _certify_bound(_sample_lower(certified), limits=limits, name=f"{name}.lower")
    upper = _certify_bound(_hull_upper(certified), limits=limits, name=f"{name}.upper")
    if lower > upper:
        raise ValueError(f"{name} sample lower exceeds convex-hull upper")
    return lower, upper


def _family_initial_bounds(
    cubics: Sequence[Cubic], *, limits: IMP1EnclosureLimits, name: str
) -> tuple[Fraction, Fraction, str]:
    if not cubics:
        raise ValueError(f"{name} must contain at least one cubic")
    hasher = sha256(_COEFFICIENT_HASH_DOMAIN)
    lower: Fraction | None = None
    upper: Fraction | None = None
    for ordinal, controls in enumerate(cubics):
        _feed_coefficient_hash(hasher, ordinal, controls)
        piece_lower, piece_upper = _piece_bounds(
            controls, limits=limits, name=f"{name}[{ordinal}]"
        )
        lower = piece_lower if lower is None else max(lower, piece_lower)
        upper = piece_upper if upper is None else max(upper, piece_upper)
    assert lower is not None and upper is not None
    return lower, upper, hasher.hexdigest()


def _heap_item(
    *,
    upper: Fraction,
    depth: int,
    family: int,
    ordinal: int,
    serial: int,
    controls: Cubic,
    lower: Fraction,
) -> tuple[object, ...]:
    return (-upper, depth, family, ordinal, serial, controls, lower)


def _current_upper(
    heap: list[tuple[object, ...]], *, lower: Fraction
) -> Fraction:
    while heap and -heap[0][0] <= lower:  # type: ignore[operator]
        heapq.heappop(heap)
    if not heap:
        return lower
    return -heap[0][0]  # type: ignore[operator]


def _pop_splittable(
    heap: list[tuple[object, ...]],
    *,
    lower: Fraction,
    max_depth: int,
) -> tuple[object, ...] | None:
    deferred: list[tuple[object, ...]] = []
    chosen: tuple[object, ...] | None = None
    while heap:
        item = heapq.heappop(heap)
        upper = -item[0]  # type: ignore[operator]
        if upper <= lower:
            continue
        if item[1] >= max_depth:  # type: ignore[operator]
            deferred.append(item)
            continue
        chosen = item
        break
    for item in deferred:
        heapq.heappush(heap, item)
    return chosen


def _seed_heap(
    cubics: Sequence[Cubic],
    *,
    family: int,
    limits: IMP1EnclosureLimits,
    name: str,
    serial: int,
) -> tuple[list[tuple[object, ...]], int]:
    heap: list[tuple[object, ...]] = []
    for ordinal, controls in enumerate(cubics):
        piece_lower, piece_upper = _piece_bounds(
            controls, limits=limits, name=f"{name}[{ordinal}]"
        )
        heapq.heappush(
            heap,
            _heap_item(
                upper=piece_upper,
                depth=0,
                family=family,
                ordinal=ordinal,
                serial=serial,
                controls=controls,
                lower=piece_lower,
            ),
        )
        serial += 1
    return heap, serial


def _split_piece(
    item: tuple[object, ...],
    heap: list[tuple[object, ...]],
    *,
    lower: Fraction,
    limits: IMP1EnclosureLimits,
    serial: int,
) -> tuple[Fraction, int, int]:
    family = item[2]  # type: ignore[index]
    depth = item[1]  # type: ignore[index]
    ordinal = item[3]  # type: ignore[index]
    controls = item[5]  # type: ignore[index]
    left, right = _split_bernstein(controls)  # type: ignore[arg-type]
    child_depth = depth + 1  # type: ignore[operator]
    for child in (left, right):
        child_lower, child_upper = _piece_bounds(
            child,
            limits=limits,
            name=f"refine[{family}][{ordinal}].child",
        )
        lower = max(lower, child_lower)
        heapq.heappush(
            heap,
            _heap_item(
                upper=child_upper,
                depth=child_depth,
                family=family,  # type: ignore[arg-type]
                ordinal=ordinal,  # type: ignore[arg-type]
                serial=serial,
                controls=child,
                lower=child_lower,
            ),
        )
        serial += 1
    return lower, child_depth, serial


def _refine_family(
    cubics: Sequence[Cubic],
    *,
    limits: IMP1EnclosureLimits,
    lower: Fraction,
    upper: Fraction,
) -> tuple[Fraction, Fraction, int, int]:
    heap, serial = _seed_heap(
        cubics, family=0, limits=limits, name="cubics", serial=0
    )
    subdivisions = 0
    maximum_depth = 0
    while lower < upper:
        if subdivisions >= limits.maximum_subdivisions_per_pair:
            break
        chosen = _pop_splittable(
            heap, lower=lower, max_depth=limits.maximum_subdivision_depth
        )
        if chosen is None:
            break
        lower, child_depth, serial = _split_piece(
            chosen, heap, lower=lower, limits=limits, serial=serial
        )
        maximum_depth = max(maximum_depth, child_depth)
        subdivisions += 1
        upper = max(lower, _current_upper(heap, lower=lower))
    return lower, max(lower, _current_upper(heap, lower=lower)), subdivisions, maximum_depth


def _refine_pair(
    d01: Sequence[Cubic],
    d12: Sequence[Cubic],
    *,
    limits: IMP1EnclosureLimits,
    d01_lower: Fraction,
    d01_upper: Fraction,
    d12_lower: Fraction,
    d12_upper: Fraction,
) -> tuple[Fraction, Fraction, Fraction, Fraction, int, int]:
    d01_heap, serial = _seed_heap(
        d01, family=0, limits=limits, name="D01", serial=0
    )
    d12_heap, serial = _seed_heap(
        d12, family=1, limits=limits, name="D12", serial=serial
    )
    lowers = [d01_lower, d12_lower]
    uppers = [d01_upper, d12_upper]
    heaps = [d01_heap, d12_heap]
    subdivisions = 0
    maximum_depth = 0

    while True:
        (
            _left,
            _right,
            _passed,
            _failed,
            inconclusive,
            _zero,
        ) = _contraction_flags(
            outer_lower=lowers[0],
            outer_upper=uppers[0],
            finest_lower=lowers[1],
            finest_upper=uppers[1],
        )
        if not inconclusive:
            break
        if subdivisions >= limits.maximum_subdivisions_per_pair:
            break
        first = _pop_splittable(
            heaps[0],
            lower=lowers[0],
            max_depth=limits.maximum_subdivision_depth,
        )
        second = _pop_splittable(
            heaps[1],
            lower=lowers[1],
            max_depth=limits.maximum_subdivision_depth,
        )
        if first is None and second is None:
            break
        if first is not None and second is not None:
            chosen = first if first[:5] <= second[:5] else second
            leftover = second if chosen is first else first
            heapq.heappush(heaps[leftover[2]], leftover)  # type: ignore[index]
        else:
            chosen = first if first is not None else second
        assert chosen is not None
        family = int(chosen[2])  # type: ignore[arg-type]
        lowers[family], child_depth, serial = _split_piece(
            chosen,
            heaps[family],
            lower=lowers[family],
            limits=limits,
            serial=serial,
        )
        maximum_depth = max(maximum_depth, child_depth)
        subdivisions += 1
        uppers[0] = max(lowers[0], _current_upper(heaps[0], lower=lowers[0]))
        uppers[1] = max(lowers[1], _current_upper(heaps[1], lower=lowers[1]))
    uppers[0] = max(lowers[0], _current_upper(heaps[0], lower=lowers[0]))
    uppers[1] = max(lowers[1], _current_upper(heaps[1], lower=lowers[1]))
    return lowers[0], uppers[0], lowers[1], uppers[1], subdivisions, maximum_depth


def _parse_rows(
    rows: Iterable[object],
    *,
    expected_row_count: int,
    limits: IMP1EnclosureLimits,
) -> tuple[
    tuple[Cubic, ...],
    tuple[Cubic, ...],
    str,
    str,
    str,
    tuple[tuple[Segment, tuple[Segment, ...], tuple[Segment, ...]], ...],
]:
    if expected_row_count > limits.maximum_owned_rows:
        raise IMP1EnclosureResourceStop(
            "maximum_owned_rows",
            detail=f"expected_row_count={expected_row_count}",
        )
    row_hasher = sha256(_ROW_HASH_DOMAIN)
    d01_hasher = sha256(_COEFFICIENT_HASH_DOMAIN)
    d12_hasher = sha256(_COEFFICIENT_HASH_DOMAIN)
    d01_cubics: list[Cubic] = []
    d12_cubics: list[Cubic] = []
    parsed_rows: list[tuple[Segment, tuple[Segment, ...], tuple[Segment, ...]]] = []
    row_count = 0
    for row_index, row in enumerate(rows):
        if row_index >= expected_row_count:
            raise ValueError("rows exceed expected_row_count")
        raw_outer, raw_medium, raw_fine = _fixed_sequence(
            row, length=3, name=f"rows[{row_index}]"
        )
        medium = _fixed_sequence(raw_medium, length=2, name=f"rows[{row_index}].medium")
        fine = _fixed_sequence(raw_fine, length=4, name=f"rows[{row_index}].fine")
        outer = _parse_segment(
            raw_outer, name=f"rows[{row_index}].outer", limits=limits
        )
        medium_segments = tuple(
            _parse_segment(
                item, name=f"rows[{row_index}].medium[{index}]", limits=limits
            )
            for index, item in enumerate(medium)
        )
        fine_segments = tuple(
            _parse_segment(
                item, name=f"rows[{row_index}].fine[{index}]", limits=limits
            )
            for index, item in enumerate(fine)
        )
        width = outer[4]
        if any(item[4] != width / 2 for item in medium_segments):
            raise ValueError("medium widths must equal half the outer width")
        if any(item[4] != width / 4 for item in fine_segments):
            raise ValueError("fine widths must equal one quarter of the outer width")
        if outer[0] != medium_segments[0][0] or outer[0] != fine_segments[0][0]:
            raise ValueError("initial values must agree across outer, medium, and fine")
        if outer[1] != medium_segments[0][1] or outer[1] != fine_segments[0][1]:
            raise ValueError("initial RHS must agree across outer, medium, and fine")
        if medium_segments[0][2] != medium_segments[1][0]:
            raise ValueError("medium endpoint values must be continuous")
        if medium_segments[0][3] != medium_segments[1][1]:
            raise ValueError("medium endpoint RHS must be continuous")
        for index in range(3):
            if fine_segments[index][2] != fine_segments[index + 1][0]:
                raise ValueError("fine endpoint values must be continuous")
            if fine_segments[index][3] != fine_segments[index + 1][1]:
                raise ValueError("fine endpoint RHS must be continuous")
        row_hasher.update(
            _row_line(row_index, outer, medium_segments, fine_segments).encode("ascii")
        )
        outer_controls = _certify_controls(
            _hermite_bernstein(outer),
            limits=limits,
            name=f"rows[{row_index}].outer.bernstein",
        )
        medium_controls = tuple(
            _certify_controls(
                _hermite_bernstein(item),
                limits=limits,
                name=f"rows[{row_index}].medium[{index}].bernstein",
            )
            for index, item in enumerate(medium_segments)
        )
        fine_controls = tuple(
            _certify_controls(
                _hermite_bernstein(item),
                limits=limits,
                name=f"rows[{row_index}].fine[{index}].bernstein",
            )
            for index, item in enumerate(fine_segments)
        )
        outer_halves = _split_bernstein(outer_controls)
        for half, child in enumerate(medium_controls):
            difference = _certify_controls(
                _subtract_bernstein(outer_halves[half], child),
                limits=limits,
                name=f"rows[{row_index}].D01[{half}]",
            )
            _feed_coefficient_hash(d01_hasher, len(d01_cubics), difference)
            d01_cubics.append(difference)
        for medium_index, parent in enumerate(medium_controls):
            parent_halves = _split_bernstein(parent)
            for half in (0, 1):
                child_index = 2 * medium_index + half
                difference = _certify_controls(
                    _subtract_bernstein(parent_halves[half], fine_controls[child_index]),
                    limits=limits,
                    name=f"rows[{row_index}].D12[{child_index}]",
                )
                _feed_coefficient_hash(d12_hasher, len(d12_cubics), difference)
                d12_cubics.append(difference)
        parsed_rows.append((outer, medium_segments, fine_segments))
        row_count += 1
    if row_count != expected_row_count:
        raise ValueError("rows must contain exactly expected_row_count refinement rows")
    if row_count > limits.maximum_owned_rows:
        raise IMP1EnclosureResourceStop(
            "maximum_owned_rows",
            detail=f"row_count={row_count}",
        )
    return (
        tuple(d01_cubics),
        tuple(d12_cubics),
        row_hasher.hexdigest(),
        d01_hasher.hexdigest(),
        d12_hasher.hexdigest(),
        tuple(parsed_rows),
    )


def _make_pair_bounds(
    *,
    d01_lower: Fraction,
    d01_upper: Fraction,
    d12_lower: Fraction,
    d12_upper: Fraction,
    row_count: int,
    row_digest: str,
    d01_digest: str,
    d12_digest: str,
    subdivisions_used: int,
    maximum_depth_reached: int,
) -> IMP1PairBounds:
    (
        _left,
        _right,
        sufficient_pass,
        sufficient_failure,
        inconclusive,
        exact_zero,
    ) = _contraction_flags(
        outer_lower=d01_lower,
        outer_upper=d01_upper,
        finest_lower=d12_lower,
        finest_upper=d12_upper,
    )
    return IMP1PairBounds(
        d01_lower=d01_lower,
        d01_upper=d01_upper,
        d12_lower=d12_lower,
        d12_upper=d12_upper,
        row_count=row_count,
        polynomial_count_d01=2 * row_count,
        polynomial_count_d12=4 * row_count,
        row_stream_sha256=row_digest,
        d01_coefficient_stream_sha256=d01_digest,
        d12_coefficient_stream_sha256=d12_digest,
        combined_coefficient_stream_sha256=_combined_coefficient_hash(
            d01_digest, d12_digest
        ),
        subdivisions_used=subdivisions_used,
        maximum_depth_reached=maximum_depth_reached,
        sufficient_pass=sufficient_pass,
        sufficient_failure=sufficient_failure,
        threshold_inconclusive=inconclusive,
        exact_zero=exact_zero,
    )


def _interval_intersection(
    left_lower: Fraction,
    left_upper: Fraction,
    right_lower: Fraction,
    right_upper: Fraction,
) -> tuple[Fraction, Fraction] | None:
    lower = max(left_lower, right_lower)
    upper = min(left_upper, right_upper)
    if lower > upper:
        return None
    return (lower, upper)


def _require_gate_equals_intersection(
    *,
    bernstein: IMP1PairBounds,
    gate_d01: tuple[Fraction, Fraction],
    gate_d12: tuple[Fraction, Fraction],
    fallback_d01: tuple[object, object],
    fallback_d12: tuple[object, object],
) -> None:
    d01 = _interval_intersection(
        bernstein.d01_lower,
        bernstein.d01_upper,
        _exact_nonnegative(fallback_d01[0], name="fallback_d01_lower"),
        _exact_nonnegative(fallback_d01[1], name="fallback_d01_upper"),
    )
    d12 = _interval_intersection(
        bernstein.d12_lower,
        bernstein.d12_upper,
        _exact_nonnegative(fallback_d12[0], name="fallback_d12_lower"),
        _exact_nonnegative(fallback_d12[1], name="fallback_d12_upper"),
    )
    if d01 is None:
        raise ValueError("empty Bernstein/fallback intersection")
    if d12 is None:
        raise ValueError("empty Bernstein/fallback intersection")
    if gate_d01 != d01 or gate_d12 != d12:
        raise ValueError("gate does not equal Bernstein/fallback intersection")


def _intersect(
    *,
    bernstein_lower: Fraction,
    bernstein_upper: Fraction,
    exact_lower: Fraction,
    exact_upper: Fraction,
    level: str,
) -> tuple[Fraction, Fraction]:
    answer = _interval_intersection(
        bernstein_lower, bernstein_upper, exact_lower, exact_upper
    )
    if answer is None:
        raise IMP1EnclosureDisagreement(
            level=level,
            detail=(
                f"bernstein=[{_canonical_text(bernstein_lower)},"
                f"{_canonical_text(bernstein_upper)}] "
                f"exact=[{_canonical_text(exact_lower)},{_canonical_text(exact_upper)}]"
            ),
        )
    return answer


def _run_fallback(
    parsed_rows: Sequence[tuple[Segment, tuple[Segment, ...], tuple[Segment, ...]]],
    *,
    limits: IMP1EnclosureLimits,
) -> tuple[Fraction, Fraction, Fraction, Fraction]:
    try:
        evidence = assess_rational_complete_c_rows(
            parsed_rows,
            expected_row_count=len(parsed_rows),
            maximum_candidates_D01=limits.fallback_maximum_candidates_D01,
            maximum_candidates_D12=limits.fallback_maximum_candidates_D12,
            refinement_depth=limits.fallback_refinement_depth,
        )
    except RationalCompleteCResourceExhausted as exc:
        raise IMP1EnclosureFallbackStop(
            exc.evidence.reason, detail=exc.evidence.detail
        ) from exc
    except RationalCompleteCRouteDisagreement as exc:
        raise IMP1EnclosureFallbackStop(
            exc.evidence.reason, detail=exc.evidence.detail
        ) from exc
    except ValueError as exc:
        if _is_int_encoding_error(exc):
            raise IMP1EnclosureFallbackStop(
                "integer_encoding_exhausted",
                detail="sealed fallback integer encoding",
            ) from exc
        raise
    bounds = (
        evidence.d01.lower,
        evidence.d01.upper,
        evidence.d12.lower,
        evidence.d12.upper,
    )
    for name, item in zip(
        (
            "fallback_d01_lower",
            "fallback_d01_upper",
            "fallback_d12_lower",
            "fallback_d12_upper",
        ),
        bounds,
        strict=True,
    ):
        _check_bits(item, limits=limits, name=name)
    return bounds


def _path_from_bounds(
    bernstein: IMP1PairBounds,
    *,
    limits: IMP1EnclosureLimits,
    gate_d01_lower: Fraction,
    gate_d01_upper: Fraction,
    gate_d12_lower: Fraction,
    gate_d12_upper: Fraction,
    accumulation_bounds: tuple[Fraction, Fraction, Fraction],
    fallback_called: bool,
    fallback_bounds: tuple[Fraction, Fraction, Fraction, Fraction] | None,
) -> IMP1PathAssessment:
    (
        left,
        right,
        sufficient_pass,
        sufficient_failure,
        inconclusive,
        exact_zero,
    ) = _contraction_flags(
        outer_lower=gate_d01_lower,
        outer_upper=gate_d01_upper,
        finest_lower=gate_d12_lower,
        finest_upper=gate_d12_upper,
        limits=limits,
    )
    decision = classify_tdg6_channel(
        CertifiedMagnitudeInterval(gate_d01_lower, gate_d01_upper),
        CertifiedMagnitudeInterval(gate_d12_lower, gate_d12_upper),
    )
    r_fine = accumulation_bounds[2]
    public = bernstein.d12_upper + r_fine
    gate_debit = gate_d12_upper + r_fine
    extra = bernstein.d12_upper - gate_d12_upper
    for labeled, item in (
        ("public_fine_debit", public),
        ("gate_debit", gate_debit),
        ("extra_enclosure_debit", extra),
        ("gate_d01_lower", gate_d01_lower),
        ("gate_d01_upper", gate_d01_upper),
        ("gate_d12_lower", gate_d12_lower),
        ("gate_d12_upper", gate_d12_upper),
    ):
        _check_bits(item, limits=limits, name=labeled)
    fallback_kwargs: dict[str, object] = {
        "fallback_called": fallback_called,
        "fallback_calls": 1 if fallback_called else 0,
    }
    if fallback_called:
        assert fallback_bounds is not None
        fallback_kwargs.update(
            {
                "fallback_d01_lower": fallback_bounds[0],
                "fallback_d01_upper": fallback_bounds[1],
                "fallback_d12_lower": fallback_bounds[2],
                "fallback_d12_upper": fallback_bounds[3],
            }
        )
    return IMP1PathAssessment(
        bernstein=bernstein,
        gate_d01_lower=gate_d01_lower,
        gate_d01_upper=gate_d01_upper,
        gate_d12_lower=gate_d12_lower,
        gate_d12_upper=gate_d12_upper,
        decision=decision,
        sufficient_pass_left=left,
        sufficient_pass_right=right,
        sufficient_contraction_pass=sufficient_pass,
        sufficient_contraction_failure=sufficient_failure,
        threshold_inconclusive=inconclusive,
        exact_zero=exact_zero,
        accumulation_bounds=accumulation_bounds,
        public_fine_debit=public,
        gate_debit=gate_debit,
        extra_enclosure_debit=extra,
        limits=limits,
        **fallback_kwargs,  # type: ignore[arg-type]
    )


def bound_exact_cubics(
    cubics: Sequence[object],
    *,
    limits: IMP1EnclosureLimits = DEFAULT_IMP1_ENCLOSURE_LIMITS,
    refine: bool = False,
) -> IMP1CubicFamilyBounds:
    """Enclose the global absolute maximum of Bernstein cubics.

    The default path is the exact five-sample lower bound and convex-hull
    upper bound.  Optional bounded de Casteljau refinement stays inside
    ``(1/3) Z[1/2]`` and therefore cannot replace an odd-prime exact maximum
    such as ``32/25`` for Bernstein ``(0, 0, 2, 1)``.
    """

    limits = _require_limits(limits)
    if isinstance(cubics, (str, bytes)) or not isinstance(cubics, Sequence):
        raise TypeError("cubics must be a sequence")
    certified = tuple(
        _certify_controls(
            _cubic(item, name=f"cubics[{index}]"),
            limits=limits,
            name=f"cubics[{index}]",
        )
        for index, item in enumerate(cubics)
    )
    if not certified:
        raise ValueError("cubics must contain at least one cubic")
    lower, upper, digest = _family_initial_bounds(
        certified, limits=limits, name="cubics"
    )
    subdivisions = 0
    maximum_depth = 0
    if refine and lower < upper:
        lower, upper, subdivisions, maximum_depth = _refine_family(
            certified, limits=limits, lower=lower, upper=upper
        )
    identically_zero = lower == 0 and upper == 0
    return IMP1CubicFamilyBounds(
        lower=lower,
        upper=upper,
        polynomial_count=len(certified),
        coefficient_stream_sha256=digest,
        subdivisions_used=subdivisions,
        maximum_depth_reached=maximum_depth,
        identically_zero=identically_zero,
    )


def assess_imp1_rows(
    rows: Iterable[object],
    *,
    expected_row_count: int,
    limits: IMP1EnclosureLimits = DEFAULT_IMP1_ENCLOSURE_LIMITS,
    accumulation_bounds: tuple[object, object, object] = (0, 0, 0),
    refine: bool = True,
    allow_fallback: bool = True,
) -> IMP1PathAssessment:
    """Assess one exact Hermite 1/2/4 row stream with the IMP1 Bernstein law."""

    limits = _require_limits(limits)
    expected = _positive_integer(expected_row_count, name="expected_row_count")
    accumulation = tuple(
        _admit_rational(
            item,
            name=f"accumulation_bounds[{index}]",
            limits=limits,
            allow_three=True,
            nonnegative=True,
        )
        for index, item in enumerate(
            _fixed_sequence(
                accumulation_bounds, length=3, name="accumulation_bounds"
            )
        )
    )
    d01, d12, row_digest, d01_digest, d12_digest, parsed = _parse_rows(
        rows, expected_row_count=expected, limits=limits
    )
    d01_lower, d01_upper, _hashed_d01 = _family_initial_bounds(
        d01, limits=limits, name="D01"
    )
    d12_lower, d12_upper, _hashed_d12 = _family_initial_bounds(
        d12, limits=limits, name="D12"
    )
    if d01_digest != _hashed_d01 or d12_digest != _hashed_d12:
        raise ValueError("Bernstein coefficient hashes diverged during construction")
    subdivisions = 0
    maximum_depth = 0
    (
        _left,
        _right,
        _passed,
        _failed,
        inconclusive,
        _zero,
    ) = _contraction_flags(
        outer_lower=d01_lower,
        outer_upper=d01_upper,
        finest_lower=d12_lower,
        finest_upper=d12_upper,
        limits=limits,
    )
    if refine and inconclusive:
        d01_lower, d01_upper, d12_lower, d12_upper, subdivisions, maximum_depth = (
            _refine_pair(
                d01,
                d12,
                limits=limits,
                d01_lower=d01_lower,
                d01_upper=d01_upper,
                d12_lower=d12_lower,
                d12_upper=d12_upper,
            )
        )
        (
            _left,
            _right,
            _passed,
            _failed,
            inconclusive,
            _zero,
        ) = _contraction_flags(
            outer_lower=d01_lower,
            outer_upper=d01_upper,
            finest_lower=d12_lower,
            finest_upper=d12_upper,
            limits=limits,
        )
    bernstein = _make_pair_bounds(
        d01_lower=d01_lower,
        d01_upper=d01_upper,
        d12_lower=d12_lower,
        d12_upper=d12_upper,
        row_count=expected,
        row_digest=row_digest,
        d01_digest=d01_digest,
        d12_digest=d12_digest,
        subdivisions_used=subdivisions,
        maximum_depth_reached=maximum_depth,
    )
    fallback_called = False
    fallback_bounds: tuple[Fraction, Fraction, Fraction, Fraction] | None = None
    gate = (d01_lower, d01_upper, d12_lower, d12_upper)
    if inconclusive and allow_fallback:
        fallback_bounds = _run_fallback(parsed, limits=limits)
        gate_d01 = _intersect(
            bernstein_lower=d01_lower,
            bernstein_upper=d01_upper,
            exact_lower=fallback_bounds[0],
            exact_upper=fallback_bounds[1],
            level="D01",
        )
        gate_d12 = _intersect(
            bernstein_lower=d12_lower,
            bernstein_upper=d12_upper,
            exact_lower=fallback_bounds[2],
            exact_upper=fallback_bounds[3],
            level="D12",
        )
        gate = (gate_d01[0], gate_d01[1], gate_d12[0], gate_d12[1])
        fallback_called = True
    return _path_from_bounds(
        bernstein,
        limits=limits,
        gate_d01_lower=gate[0],
        gate_d01_upper=gate[1],
        gate_d12_lower=gate[2],
        gate_d12_upper=gate[3],
        accumulation_bounds=accumulation,  # type: ignore[arg-type]
        fallback_called=fallback_called,
        fallback_bounds=fallback_bounds,
    )


def _require_family(family: object) -> ValidatedRecordedFamily:
    if not isinstance(family, ValidatedRecordedFamily):
        raise TypeError("family must be ValidatedRecordedFamily")
    replayed = validate_recorded_family(
        family.paths,
        method=family.method,
        arithmetic_id=family.arithmetic_id,
        expected_owned_row_count=family.owned_row_count,
    )
    if replayed.family_sha256 != family.family_sha256:
        raise RecordedFamilyError("family hash does not match recorded content")
    return replayed


def assess_imp1_channel(
    family: ValidatedRecordedFamily,
    channel: str,
    *,
    limits: IMP1EnclosureLimits = DEFAULT_IMP1_ENCLOSURE_LIMITS,
) -> IMP1ChannelAssessment:
    """Replay one recorded family channel and return the IMP1 assessment."""

    limits = _require_limits(limits)
    replayed = _require_family(family)
    if replayed.owned_row_count > limits.maximum_owned_rows:
        raise IMP1EnclosureResourceStop(
            "maximum_owned_rows",
            detail=f"owned_row_count={replayed.owned_row_count}",
        )
    rebuilt = reconstruct_channel(replayed, channel)
    raw = assess_imp1_rows(
        rebuilt.raw_rows,
        expected_row_count=rebuilt.row_count,
        limits=limits,
        accumulation_bounds=rebuilt.accumulation_bounds,
        refine=False,
        allow_fallback=False,
    )
    corrected = assess_imp1_rows(
        rebuilt.corrected_rows,
        expected_row_count=rebuilt.row_count,
        limits=limits,
        accumulation_bounds=rebuilt.accumulation_bounds,
        refine=True,
        allow_fallback=True,
    )
    return IMP1ChannelAssessment(
        channel=channel,
        method=replayed.method,
        family_sha256=replayed.family_sha256,
        row_count=rebuilt.row_count,
        limits=limits,
        raw=raw.bernstein,
        raw_unresolved=raw.bernstein.threshold_inconclusive,
        corrected=corrected,
    )


def assess_imp1_family(
    family: ValidatedRecordedFamily,
    *,
    limits: IMP1EnclosureLimits = DEFAULT_IMP1_ENCLOSURE_LIMITS,
) -> IMP1FamilyAssessment:
    """Assess all eighteen complete-state channels of one recorded family."""

    limits = _require_limits(limits)
    replayed = _require_family(family)
    records = tuple(
        assess_imp1_channel(replayed, channel, limits=limits)
        for channel in TDG6_COMPLETE_STATE_CHANNELS
    )
    debits = tuple(record.corrected.public_fine_debit for record in records)
    return IMP1FamilyAssessment(
        method=replayed.method,
        family_sha256=replayed.family_sha256,
        row_count=replayed.owned_row_count,
        limits=limits,
        channels=records,
        admission_passed=all(
            record.corrected.decision.admission_passed for record in records
        ),
        public_fine_debits=debits,
    )


__all__ = [
    "DEFAULT_IMP1_ENCLOSURE_LIMITS",
    "IMP1_ENCLOSURE_EVALUATOR_ID",
    "IMP1_ENCLOSURE_SCHEMA_VERSION",
    "IMP1ChannelAssessment",
    "IMP1CubicFamilyBounds",
    "IMP1EnclosureDisagreement",
    "IMP1EnclosureFallbackStop",
    "IMP1EnclosureLimits",
    "IMP1EnclosureResourceStop",
    "IMP1FamilyAssessment",
    "IMP1PairBounds",
    "IMP1PathAssessment",
    "assess_imp1_channel",
    "assess_imp1_family",
    "assess_imp1_rows",
    "bound_exact_cubics",
]
