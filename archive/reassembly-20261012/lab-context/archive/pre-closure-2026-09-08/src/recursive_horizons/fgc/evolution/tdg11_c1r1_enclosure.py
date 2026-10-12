"""C1R1 exact Bernstein enclosure of the selected C1 numerical object.

This module is a prospective representation of the same corrected cubics
evaluated by TDG11-IMP1.  It does not change the mathematical object, the
eighteen-channel order, ``p >= 3/2``, public debit, or resource law.

Returned ``IMP1*`` dataclasses, hash domains and evaluator identifiers are
the IMP1 reference-wire format.  Using them is a compatibility claim that
the C1R1 representation produced the same certified numbers; it is not a
claim that the sealed IMP1 implementation ran.  Implementation provenance
is ``C1R1_IMPLEMENTATION_ID``.

Inputs outside the dyadic module ``(1/3)Z[1/2]`` (or the dyadic-only
slots) take the unchanged IMP1 reference path.  They are not made
admissible, omitted, or given a relaxed limit.  A one-pass iterator is
replayed from its consumed prefix so the reference path sees the same
logical stream; at most ``expected_row_count + 1`` rows are retained,
matching IMP1's extra-row detection budget.  Dual-localizer fallback and
TDG6 classification remain the sealed independent instruments.

Adaptive refinement, heap ties and squared bit-size witnesses follow the
IMP1 decision order.  Coefficient streams are hashed once from reduced
integer/exponent text; they are not re-certified through Fraction gcd.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
import heapq
from typing import Final, Iterable, Iterator, Sequence

from .tdg6_temporal_admission_design import (
    CertifiedMagnitudeInterval,
    TDG6_COMPLETE_STATE_CHANNELS,
    classify_tdg6_channel,
)
from .tdg11_c1r1_ring import (
    RING_ZERO,
    Cubic,
    DirectRingConversionError,
    RingElement,
    cmp,
    cubic_from_inputs,
    div_pow2,
    from_input,
    hermite_bernstein,
    hull_upper,
    int_to_decimal_text,
    mul_int,
    sample_lower,
    split_bernstein,
    subtract_bernstein,
)
from .tdg11_imp1_enclosure import (
    DEFAULT_IMP1_ENCLOSURE_LIMITS,
    IMP1_ENCLOSURE_EVALUATOR_ID,
    IMP1_ENCLOSURE_SCHEMA_VERSION,
    IMP1_ORDER_SQUARED_MULTIPLIER,
    IMP1ChannelAssessment,
    IMP1CubicFamilyBounds,
    IMP1EnclosureDisagreement,
    IMP1EnclosureFallbackStop,
    IMP1EnclosureLimits,
    IMP1EnclosureResourceStop,
    IMP1FamilyAssessment,
    IMP1PairBounds,
    IMP1PathAssessment,
    assess_imp1_channel as _imp1_reference_assess_channel,
    assess_imp1_family as _imp1_reference_assess_family,
    assess_imp1_rows as _imp1_reference_assess_rows,
    bound_exact_cubics as _imp1_reference_bound_exact_cubics,
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
C1R1_IMPLEMENTATION_ID: Final[str] = "tdg11_c1r1_integer_exponent_direct_ring_v1"
C1R1_MATHEMATICAL_OBJECT: Final[str] = (
    "exact_accumulation_reconstruction_with_debit"
)
C1R1_REFERENCE_WIRE_EVALUATOR_ID: Final[str] = IMP1_ENCLOSURE_EVALUATOR_ID
C1R1_REFERENCE_WIRE_SCHEMA_VERSION: Final[int] = IMP1_ENCLOSURE_SCHEMA_VERSION

# IMP1-owned hash domains.  Emitting them is reference-wire compatibility,
# not IMP1 execution provenance.
_IMP1_REFERENCE_WIRE_ROW_HASH_DOMAIN: Final[bytes] = (
    b"TDG11-IMP1-BERNSTEIN-ROW-STREAM-v1\n"
)
_IMP1_REFERENCE_WIRE_COEFFICIENT_HASH_DOMAIN: Final[bytes] = (
    b"TDG11-IMP1-BERNSTEIN-CUBIC-STREAM-v1\n"
)
_IMP1_REFERENCE_WIRE_COMBINED_HASH_DOMAIN: Final[bytes] = (
    b"TDG11-IMP1-BERNSTEIN-COMBINED-v1\n"
)
_INT_ENCODING_MARKERS: Final[tuple[str, ...]] = (
    "integer string conversion",
    "int_max_str_digits",
)

Segment = tuple[Fraction, Fraction, Fraction, Fraction, Fraction]
RingSegment = tuple[RingElement, RingElement, RingElement, RingElement, RingElement]


class _BoundedRowReplay:
    """Replay a one-pass row stream from its consumed prefix.

    IMP1 pulls at most ``expected_row_count + 1`` rows (the extra row only
    to detect overflow).  This buffer uses that same finite budget so a
    mid-stream ``DirectRingConversionError`` can hand the original logical
    sequence to the sealed reference path without unbounded materialization
    and without changing ``maximum_owned_rows``.
    """

    def __init__(self, source: Iterable[object], *, bound: int) -> None:
        if type(bound) is not int or bound < 1:
            raise ValueError("row replay bound must be a positive integer")
        self._source: Iterator[object] = iter(source)
        self._buffer: list[object] = []
        self._bound = bound

    def __iter__(self) -> Iterator[object]:
        yield from self._buffer
        for item in self._source:
            if len(self._buffer) < self._bound:
                self._buffer.append(item)
            yield item


def _replayable_rows(rows: object, *, bound: int) -> Iterable[object]:
    if isinstance(rows, Sequence) and not isinstance(rows, (str, bytes)):
        return rows
    return _BoundedRowReplay(rows, bound=bound)  # type: ignore[arg-type]


def _positive_integer(value: object, *, name: str, allow_zero: bool = False) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{name} must be an integer")
    if value < (0 if allow_zero else 1):
        qualifier = "nonnegative" if allow_zero else "positive"
        raise ValueError(f"{name} must be {qualifier}")
    return value


def _fixed_sequence(value: object, *, length: int, name: str) -> tuple[object, ...]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise TypeError(f"{name} must be a sequence")
    answer = tuple(value)
    if len(answer) != length:
        raise ValueError(f"{name} must contain exactly {length} entries")
    return answer


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
    value: RingElement, *, limits: IMP1EnclosureLimits, name: str
) -> RingElement:
    if value.bit_size() > limits.maximum_rational_bits:
        raise IMP1EnclosureResourceStop("maximum_rational_bits", detail=name)
    return value


def _certify(
    value: RingElement,
    *,
    limits: IMP1EnclosureLimits,
    name: str,
    allow_three: bool = True,
    nonnegative: bool = False,
) -> RingElement:
    if nonnegative and not value.is_nonnegative():
        raise ValueError(f"{name} must be nonnegative")
    if allow_three:
        if not value.in_direct_ring():
            raise ValueError(
                f"{name} is not in the direct Bernstein ring (1/3)Z[1/2]"
            )
    elif not value.is_dyadic():
        raise ValueError(f"{name} is not an exact dyadic rational")
    return _check_bits(value, limits=limits, name=name)


def _certify_cubic(
    controls: Cubic, *, limits: IMP1EnclosureLimits, name: str
) -> Cubic:
    certified: list[RingElement] = []
    for index, item in enumerate(controls):
        certified.append(
            _certify(item, limits=limits, name=f"{name}[{index}]")
        )
    return tuple(certified)  # type: ignore[return-value]


def _is_int_encoding_error(exc: BaseException) -> bool:
    if not isinstance(exc, ValueError):
        return False
    text = str(exc)
    return any(marker in text for marker in _INT_ENCODING_MARKERS)


def _feed_coefficient_hash(hasher: object, ordinal: int, controls: Cubic) -> None:
    line = (
        f"{ordinal}|"
        + "|".join(item.canonical_text() for item in controls)
        + "\n"
    )
    hasher.update(line.encode("ascii"))  # type: ignore[attr-defined]


def _combined_coefficient_hash(outer_digest: str, finest_digest: str) -> str:
    return sha256(
        _IMP1_REFERENCE_WIRE_COMBINED_HASH_DOMAIN
        + f"{outer_digest}\n{finest_digest}\n".encode("ascii")
    ).hexdigest()


def _row_line(
    ordinal: int,
    outer: RingSegment,
    medium: Sequence[RingSegment],
    fine: Sequence[RingSegment],
) -> str:
    values: list[RingElement] = []
    for segment in (outer, *medium, *fine):
        values.extend(segment)
    encoded = "|".join(item.canonical_text() for item in values)
    return f"{ordinal}|{encoded}\n"


def _max_ring(left: RingElement, right: RingElement) -> RingElement:
    return left if cmp(left, right) >= 0 else right


def _contraction_flags(
    *,
    outer_lower: RingElement,
    outer_upper: RingElement,
    finest_lower: RingElement,
    finest_upper: RingElement,
    limits: IMP1EnclosureLimits | None = None,
) -> tuple[Fraction, Fraction, bool, bool, bool, bool]:
    left = mul_int(finest_upper.square(), IMP1_ORDER_SQUARED_MULTIPLIER)
    right = outer_lower.square()
    if limits is not None:
        _check_bits(left, limits=limits, name="sufficient_pass_left")
        _check_bits(right, limits=limits, name="sufficient_pass_right")
    exact_zero = outer_upper.is_zero() and finest_upper.is_zero()
    if exact_zero:
        return left.as_fraction(), right.as_fraction(), False, False, False, True
    if cmp(left, right) <= 0:
        return left.as_fraction(), right.as_fraction(), True, False, False, False
    failure_left = mul_int(finest_lower.square(), IMP1_ORDER_SQUARED_MULTIPLIER)
    failure_right = outer_upper.square()
    if limits is not None:
        _check_bits(failure_left, limits=limits, name="sufficient_failure_left")
        _check_bits(failure_right, limits=limits, name="sufficient_failure_right")
    if cmp(failure_left, failure_right) > 0:
        return left.as_fraction(), right.as_fraction(), False, True, False, False
    return left.as_fraction(), right.as_fraction(), False, False, True, False


def _check_bits_fraction(
    value: Fraction, *, limits: IMP1EnclosureLimits, name: str
) -> Fraction:
    size = max(abs(value.numerator).bit_length(), value.denominator.bit_length())
    if size > limits.maximum_rational_bits:
        raise IMP1EnclosureResourceStop("maximum_rational_bits", detail=name)
    return value


def _contraction_flags_fraction(
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
        _check_bits_fraction(left, limits=limits, name="sufficient_pass_left")
        _check_bits_fraction(right, limits=limits, name="sufficient_pass_right")
    exact_zero = outer_upper == 0 and finest_upper == 0
    if exact_zero:
        return left, right, False, False, False, True
    if left <= right:
        return left, right, True, False, False, False
    failure_left = IMP1_ORDER_SQUARED_MULTIPLIER * finest_lower * finest_lower
    failure_right = outer_upper * outer_upper
    if limits is not None:
        _check_bits_fraction(
            failure_left, limits=limits, name="sufficient_failure_left"
        )
        _check_bits_fraction(
            failure_right, limits=limits, name="sufficient_failure_right"
        )
    if failure_left > failure_right:
        return left, right, False, True, False, False
    return left, right, False, False, True, False


def _piece_bounds(
    controls: Cubic, *, limits: IMP1EnclosureLimits, name: str
) -> tuple[RingElement, RingElement]:
    lower = _certify(
        sample_lower(controls),
        limits=limits,
        name=f"{name}.lower",
        nonnegative=True,
    )
    upper = _certify(
        hull_upper(controls),
        limits=limits,
        name=f"{name}.upper",
        nonnegative=True,
    )
    if cmp(lower, upper) > 0:
        raise ValueError(f"{name} sample lower exceeds convex-hull upper")
    return lower, upper


def _family_initial_bounds(
    cubics: Sequence[Cubic],
    *,
    limits: IMP1EnclosureLimits,
    name: str,
    hasher: object | None = None,
) -> tuple[RingElement, RingElement, str | None]:
    if not cubics:
        raise ValueError(f"{name} must contain at least one cubic")
    lower: RingElement | None = None
    upper: RingElement | None = None
    for ordinal, controls in enumerate(cubics):
        if hasher is not None:
            _feed_coefficient_hash(hasher, ordinal, controls)
        piece_lower, piece_upper = _piece_bounds(
            controls, limits=limits, name=f"{name}[{ordinal}]"
        )
        lower = piece_lower if lower is None else _max_ring(lower, piece_lower)
        upper = piece_upper if upper is None else _max_ring(upper, piece_upper)
    assert lower is not None and upper is not None
    digest = hasher.hexdigest() if hasher is not None else None  # type: ignore[union-attr]
    return lower, upper, digest


def _heap_item(
    *,
    upper: RingElement,
    depth: int,
    family: int,
    ordinal: int,
    serial: int,
    controls: Cubic,
    lower: RingElement,
) -> tuple[object, ...]:
    return (
        -upper.as_fraction(),
        depth,
        family,
        ordinal,
        serial,
        controls,
        lower.as_fraction(),
    )


def _current_upper(
    heap: list[tuple[object, ...]], *, lower: RingElement
) -> RingElement:
    lower_fraction = lower.as_fraction()
    while heap and -heap[0][0] <= lower_fraction:  # type: ignore[operator]
        heapq.heappop(heap)
    if not heap:
        return lower
    return from_input(-heap[0][0], name="heap_upper", allow_three=True)  # type: ignore[arg-type]


def _pop_splittable(
    heap: list[tuple[object, ...]],
    *,
    lower: RingElement,
    max_depth: int,
) -> tuple[object, ...] | None:
    lower_fraction = lower.as_fraction()
    deferred: list[tuple[object, ...]] = []
    chosen: tuple[object, ...] | None = None
    while heap:
        item = heapq.heappop(heap)
        upper = -item[0]  # type: ignore[operator]
        if upper <= lower_fraction:
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
    lower: RingElement,
    limits: IMP1EnclosureLimits,
    serial: int,
) -> tuple[RingElement, int, int]:
    family = item[2]
    depth = item[1]
    ordinal = item[3]
    controls = item[5]
    left, right = split_bernstein(controls)  # type: ignore[arg-type]
    child_depth = depth + 1  # type: ignore[operator]
    for child in (left, right):
        _certify_cubic(
            child,
            limits=limits,
            name=f"refine[{family}][{ordinal}].child",
        )
        child_lower, child_upper = _piece_bounds(
            child,
            limits=limits,
            name=f"refine[{family}][{ordinal}].child",
        )
        lower = _max_ring(lower, child_lower)
        heapq.heappush(
            heap,
            _heap_item(
                upper=child_upper,
                depth=child_depth,  # type: ignore[arg-type]
                family=family,  # type: ignore[arg-type]
                ordinal=ordinal,  # type: ignore[arg-type]
                serial=serial,
                controls=child,
                lower=child_lower,
            ),
        )
        serial += 1
    return lower, child_depth, serial  # type: ignore[return-value]


def _refine_family(
    cubics: Sequence[Cubic],
    *,
    limits: IMP1EnclosureLimits,
    lower: RingElement,
    upper: RingElement,
) -> tuple[RingElement, RingElement, int, int]:
    heap, serial = _seed_heap(
        cubics, family=0, limits=limits, name="cubics", serial=0
    )
    subdivisions = 0
    maximum_depth = 0
    while cmp(lower, upper) < 0:
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
        upper = _max_ring(lower, _current_upper(heap, lower=lower))
    return lower, _max_ring(lower, _current_upper(heap, lower=lower)), subdivisions, maximum_depth


def _refine_pair(
    d01: Sequence[Cubic],
    d12: Sequence[Cubic],
    *,
    limits: IMP1EnclosureLimits,
    d01_lower: RingElement,
    d01_upper: RingElement,
    d12_lower: RingElement,
    d12_upper: RingElement,
) -> tuple[RingElement, RingElement, RingElement, RingElement, int, int]:
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
        uppers[0] = _max_ring(lowers[0], _current_upper(heaps[0], lower=lowers[0]))
        uppers[1] = _max_ring(lowers[1], _current_upper(heaps[1], lower=lowers[1]))
    uppers[0] = _max_ring(lowers[0], _current_upper(heaps[0], lower=lowers[0]))
    uppers[1] = _max_ring(lowers[1], _current_upper(heaps[1], lower=lowers[1]))
    return lowers[0], uppers[0], lowers[1], uppers[1], subdivisions, maximum_depth


def _parse_segment(
    value: object, *, name: str, limits: IMP1EnclosureLimits
) -> tuple[RingSegment, Segment]:
    raw = _fixed_sequence(value, length=5, name=name)
    ring_values: list[RingElement] = []
    for index, item in enumerate(raw):
        labeled = f"{name}[{index}]"
        converted = from_input(
            item,
            name=labeled,
            allow_three=index in (0, 2),
        )
        _certify(
            converted,
            limits=limits,
            name=labeled,
            allow_three=index in (0, 2),
        )
        ring_values.append(converted)
    if cmp(ring_values[4], RING_ZERO) <= 0:
        raise ValueError(f"{name} width must be positive")
    exact = tuple(item.as_fraction() for item in ring_values)
    return tuple(ring_values), exact  # type: ignore[return-value]


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
    row_hasher = sha256(_IMP1_REFERENCE_WIRE_ROW_HASH_DOMAIN)
    d01_hasher = sha256(_IMP1_REFERENCE_WIRE_COEFFICIENT_HASH_DOMAIN)
    d12_hasher = sha256(_IMP1_REFERENCE_WIRE_COEFFICIENT_HASH_DOMAIN)
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
        outer_ring, outer = _parse_segment(
            raw_outer, name=f"rows[{row_index}].outer", limits=limits
        )
        medium_ring: list[RingSegment] = []
        medium_segments: list[Segment] = []
        for index, item in enumerate(medium):
            ring_segment, fraction_segment = _parse_segment(
                item, name=f"rows[{row_index}].medium[{index}]", limits=limits
            )
            medium_ring.append(ring_segment)
            medium_segments.append(fraction_segment)
        fine_ring: list[RingSegment] = []
        fine_segments: list[Segment] = []
        for index, item in enumerate(fine):
            ring_segment, fraction_segment = _parse_segment(
                item, name=f"rows[{row_index}].fine[{index}]", limits=limits
            )
            fine_ring.append(ring_segment)
            fine_segments.append(fraction_segment)
        width = outer_ring[4]
        medium_width = div_pow2(width, 1)
        fine_width = div_pow2(width, 2)
        if any(item[4] != medium_width for item in medium_ring):
            raise ValueError("medium widths must equal half the outer width")
        if any(item[4] != fine_width for item in fine_ring):
            raise ValueError("fine widths must equal one quarter of the outer width")
        if outer_ring[0] != medium_ring[0][0] or outer_ring[0] != fine_ring[0][0]:
            raise ValueError("initial values must agree across outer, medium, and fine")
        if outer_ring[1] != medium_ring[0][1] or outer_ring[1] != fine_ring[0][1]:
            raise ValueError("initial RHS must agree across outer, medium, and fine")
        if medium_ring[0][2] != medium_ring[1][0]:
            raise ValueError("medium endpoint values must be continuous")
        if medium_ring[0][3] != medium_ring[1][1]:
            raise ValueError("medium endpoint RHS must be continuous")
        for index in range(3):
            if fine_ring[index][2] != fine_ring[index + 1][0]:
                raise ValueError("fine endpoint values must be continuous")
            if fine_ring[index][3] != fine_ring[index + 1][1]:
                raise ValueError("fine endpoint RHS must be continuous")
        row_hasher.update(
            _row_line(row_index, outer_ring, medium_ring, fine_ring).encode("ascii")
        )
        outer_controls = _certify_cubic(
            hermite_bernstein(*outer_ring),
            limits=limits,
            name=f"rows[{row_index}].outer.bernstein",
        )
        medium_controls = tuple(
            _certify_cubic(
                hermite_bernstein(*item),
                limits=limits,
                name=f"rows[{row_index}].medium[{index}].bernstein",
            )
            for index, item in enumerate(medium_ring)
        )
        fine_controls = tuple(
            _certify_cubic(
                hermite_bernstein(*item),
                limits=limits,
                name=f"rows[{row_index}].fine[{index}].bernstein",
            )
            for index, item in enumerate(fine_ring)
        )
        outer_halves = split_bernstein(outer_controls)
        for half_index, child in enumerate(medium_controls):
            difference = _certify_cubic(
                subtract_bernstein(outer_halves[half_index], child),
                limits=limits,
                name=f"rows[{row_index}].D01[{half_index}]",
            )
            _feed_coefficient_hash(d01_hasher, len(d01_cubics), difference)
            d01_cubics.append(difference)
        for medium_index, parent in enumerate(medium_controls):
            parent_halves = split_bernstein(parent)
            for half_index in (0, 1):
                child_index = 2 * medium_index + half_index
                difference = _certify_cubic(
                    subtract_bernstein(
                        parent_halves[half_index], fine_controls[child_index]
                    ),
                    limits=limits,
                    name=f"rows[{row_index}].D12[{child_index}]",
                )
                _feed_coefficient_hash(d12_hasher, len(d12_cubics), difference)
                d12_cubics.append(difference)
        parsed_rows.append((outer, tuple(medium_segments), tuple(fine_segments)))
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
    d01_lower: RingElement,
    d01_upper: RingElement,
    d12_lower: RingElement,
    d12_upper: RingElement,
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
        d01_lower=d01_lower.as_fraction(),
        d01_upper=d01_upper.as_fraction(),
        d12_lower=d12_lower.as_fraction(),
        d12_upper=d12_upper.as_fraction(),
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


def _canonical_text(value: Fraction) -> str:
    return int_to_decimal_text(value.numerator) + "/" + int_to_decimal_text(
        value.denominator
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
        _check_bits_fraction(item, limits=limits, name=name)
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
    ) = _contraction_flags_fraction(
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
        _check_bits_fraction(item, limits=limits, name=labeled)
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


def _assess_rows_ring(
    rows: Iterable[object],
    *,
    expected_row_count: int,
    limits: IMP1EnclosureLimits,
    accumulation_bounds: tuple[object, object, object],
    refine: bool,
    allow_fallback: bool,
) -> IMP1PathAssessment:
    expected = _positive_integer(expected_row_count, name="expected_row_count")
    accumulation = tuple(
        _certify(
            from_input(
                item,
                name=f"accumulation_bounds[{index}]",
                allow_three=True,
            ),
            limits=limits,
            name=f"accumulation_bounds[{index}]",
            nonnegative=True,
        )
        for index, item in enumerate(
            _fixed_sequence(accumulation_bounds, length=3, name="accumulation_bounds")
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
    gate = (
        d01_lower.as_fraction(),
        d01_upper.as_fraction(),
        d12_lower.as_fraction(),
        d12_upper.as_fraction(),
    )
    if inconclusive and allow_fallback:
        fallback_bounds = _run_fallback(parsed, limits=limits)
        gate_d01 = _intersect(
            bernstein_lower=gate[0],
            bernstein_upper=gate[1],
            exact_lower=fallback_bounds[0],
            exact_upper=fallback_bounds[1],
            level="D01",
        )
        gate_d12 = _intersect(
            bernstein_lower=gate[2],
            bernstein_upper=gate[3],
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
        accumulation_bounds=tuple(item.as_fraction() for item in accumulation),
        fallback_called=fallback_called,
        fallback_bounds=fallback_bounds,
    )


def bound_c1r1_cubics(
    cubics: Sequence[object],
    *,
    limits: IMP1EnclosureLimits = DEFAULT_IMP1_ENCLOSURE_LIMITS,
    refine: bool = False,
) -> IMP1CubicFamilyBounds:
    """Enclose the global absolute maximum of Bernstein cubics via C1R1."""

    limits = _require_limits(limits)
    if isinstance(cubics, (str, bytes)) or not isinstance(cubics, Sequence):
        raise TypeError("cubics must be a sequence")
    try:
        certified = tuple(
            _certify_cubic(
                cubic_from_inputs(item, name=f"cubics[{index}]"),
                limits=limits,
                name=f"cubics[{index}]",
            )
            for index, item in enumerate(cubics)
        )
    except DirectRingConversionError:
        return _imp1_reference_bound_exact_cubics(
            cubics, limits=limits, refine=refine
        )
    if not certified:
        raise ValueError("cubics must contain at least one cubic")
    hasher = sha256(_IMP1_REFERENCE_WIRE_COEFFICIENT_HASH_DOMAIN)
    lower, upper, digest = _family_initial_bounds(
        certified, limits=limits, name="cubics", hasher=hasher
    )
    assert digest is not None
    subdivisions = 0
    maximum_depth = 0
    if refine and cmp(lower, upper) < 0:
        lower, upper, subdivisions, maximum_depth = _refine_family(
            certified, limits=limits, lower=lower, upper=upper
        )
    identically_zero_flag = lower.is_zero() and upper.is_zero()
    return IMP1CubicFamilyBounds(
        lower=lower.as_fraction(),
        upper=upper.as_fraction(),
        polynomial_count=len(certified),
        coefficient_stream_sha256=digest,
        subdivisions_used=subdivisions,
        maximum_depth_reached=maximum_depth,
        identically_zero=identically_zero_flag,
    )


def assess_c1r1_rows(
    rows: Iterable[object],
    *,
    expected_row_count: int,
    limits: IMP1EnclosureLimits = DEFAULT_IMP1_ENCLOSURE_LIMITS,
    accumulation_bounds: tuple[object, object, object] = (0, 0, 0),
    refine: bool = True,
    allow_fallback: bool = True,
) -> IMP1PathAssessment:
    """Assess one exact Hermite 1/2/4 row stream with the C1R1 ring."""

    limits = _require_limits(limits)
    expected = _positive_integer(expected_row_count, name="expected_row_count")
    stream = _replayable_rows(rows, bound=expected + 1)
    try:
        return _assess_rows_ring(
            stream,
            expected_row_count=expected,
            limits=limits,
            accumulation_bounds=accumulation_bounds,
            refine=refine,
            allow_fallback=allow_fallback,
        )
    except DirectRingConversionError:
        return _imp1_reference_assess_rows(
            stream,
            expected_row_count=expected,
            limits=limits,
            accumulation_bounds=accumulation_bounds,
            refine=refine,
            allow_fallback=allow_fallback,
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


def assess_c1r1_channel(
    family: ValidatedRecordedFamily,
    channel: str,
    *,
    limits: IMP1EnclosureLimits = DEFAULT_IMP1_ENCLOSURE_LIMITS,
) -> IMP1ChannelAssessment:
    """Replay one recorded family channel and return the IMP1-wire assessment."""

    limits = _require_limits(limits)
    replayed = _require_family(family)
    if replayed.owned_row_count > limits.maximum_owned_rows:
        raise IMP1EnclosureResourceStop(
            "maximum_owned_rows",
            detail=f"owned_row_count={replayed.owned_row_count}",
        )
    rebuilt = reconstruct_channel(replayed, channel)
    try:
        raw = assess_c1r1_rows(
            rebuilt.raw_rows,
            expected_row_count=rebuilt.row_count,
            limits=limits,
            accumulation_bounds=rebuilt.accumulation_bounds,
            refine=False,
            allow_fallback=False,
        )
        corrected = assess_c1r1_rows(
            rebuilt.corrected_rows,
            expected_row_count=rebuilt.row_count,
            limits=limits,
            accumulation_bounds=rebuilt.accumulation_bounds,
            refine=True,
            allow_fallback=True,
        )
    except DirectRingConversionError:
        return _imp1_reference_assess_channel(replayed, channel, limits=limits)
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


def assess_c1r1_family(
    family: ValidatedRecordedFamily,
    *,
    limits: IMP1EnclosureLimits = DEFAULT_IMP1_ENCLOSURE_LIMITS,
) -> IMP1FamilyAssessment:
    """Assess all eighteen complete-state channels with the C1R1 ring."""

    limits = _require_limits(limits)
    replayed = _require_family(family)
    try:
        records = tuple(
            assess_c1r1_channel(replayed, channel, limits=limits)
            for channel in TDG6_COMPLETE_STATE_CHANNELS
        )
    except DirectRingConversionError:
        return _imp1_reference_assess_family(replayed, limits=limits)
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
    "C1R1_IMPLEMENTATION_ID",
    "C1R1_MATHEMATICAL_OBJECT",
    "C1R1_REFERENCE_WIRE_EVALUATOR_ID",
    "C1R1_REFERENCE_WIRE_SCHEMA_VERSION",
    "assess_c1r1_channel",
    "assess_c1r1_family",
    "assess_c1r1_rows",
    "bound_c1r1_cubics",
]
