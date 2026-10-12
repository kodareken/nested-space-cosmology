"""Rational dual-route complete-C adapter for prospective TDG11.

This module is a standard-library instrument, not an authorization, binder,
runtime, or campaign calculator.  It consumes an iterable of exact rational
cubic-Hermite refinement rows and never restores state.  Input segments are
``Fraction`` or built-in ``int`` only; binary64 Hermite reconstruction is
not used.

One input row is ``(outer, medium, fine)``.  ``outer`` is one segment,
``medium`` two segments, and ``fine`` four segments.  A segment is the plain
five-tuple ``(left_value, left_rhs, right_value, right_rhs, width)``.  Exact
cubics are formed from these rationals, then the public TDG9 operations
``restrict_exact_cubic_to_half`` and ``subtract_exact_cubics``.  Each row
contributes two ``D01`` cubics and four ``D12`` cubics.  Widths must be the
exact positive uniform triple ``h``, ``h/2``, ``h/4``.  The three levels must
share the same initial value and initial RHS, and values/RHS must be
continuous at every intra-level endpoint.

The global absolute maximum of each complete-C family is enclosed by both
existing exact localizers ``localize_absolute_maximum`` and
``localize_absolute_maximum_independently_v2``.  The routes must agree on
polynomial counts, candidate counts, the independent stationary-count digest,
survivor keys, overlapping survivor intervals, and global interval overlap.
Disagreement or resource exhaustion fails closed as a typed exception, never
as a scientific pass or failure.  Compact evidence keeps exact rational
bounds, counts, and digests, and does not retain cubics or candidates.

Channel admission is exactly ``classify_tdg6_channel`` applied to the
dual-route intersection through ``CertifiedMagnitudeInterval``.  Separate
sufficient-pass, sufficient-failure, threshold-inconclusive, and exact-zero
flags are retained so a merely failed ``8 U12^2 <= L01^2`` test is not
advertised as a proven failure.  Coefficient and row hashes use explicit
TDG11 domains, distinct from the TDG10 replay domains.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
from typing import Final, Iterable, Sequence

from .tdg6_temporal_admission_design import (
    CertifiedMagnitudeInterval,
    TDG6ChannelAdmission,
    classify_tdg6_channel,
)
from .tdg9_exact_temporal_arithmetic import (
    restrict_exact_cubic_to_half,
    subtract_exact_cubics,
)
from .tdg9_local_extrema import (
    EVALUATOR_ID as PRIMARY_LOCALIZER_ID,
    LocalCubic,
    localize_absolute_maximum,
)
from .tdg9_local_extrema import _stationary_intervals as _primary_stationary_intervals
from .tdg9_local_extrema_independent_v2 import (
    EVALUATOR_ID as INDEPENDENT_LOCALIZER_ID,
    IndependentLocalCubicV2,
    RootIsolationInconclusive,
    localize_absolute_maximum_independently_v2,
)


Q = Fraction
TDG11_RATIONAL_COMPLETE_C_SCHEMA_VERSION: Final[int] = 1
TDG11_RATIONAL_COMPLETE_C_EVALUATOR_ID: Final[str] = (
    "tdg11_rational_complete_c_dual_route_v1"
)
TDG11_ORDER_SQUARED_MULTIPLIER: Final[int] = 8
_ROW_HASH_DOMAIN: Final[bytes] = b"TDG11-RATIONAL-COMPLETE-C-ROW-STREAM-v1\n"
_COEFFICIENT_HASH_DOMAIN: Final[bytes] = b"TDG11-RATIONAL-COMPLETE-C-CUBIC-STREAM-v1\n"
_COMBINED_HASH_DOMAIN: Final[bytes] = b"TDG11-RATIONAL-COMPLETE-C-COMBINED-v1\n"
_SURVIVOR_KEY_HASH_DOMAIN: Final[bytes] = (
    b"TDG11-RATIONAL-COMPLETE-C-SURVIVOR-KEY-STREAM-v1\n"
)
_INDEPENDENT_STATIONARY_COUNT_DOMAIN: Final[bytes] = (
    b"TDG9-LOC2-STATIONARY-COUNT-STREAM-v1\n"
)
_CLOSED_REASONS: Final[frozenset[str]] = frozenset(
    {
        "candidate_ceiling_exhausted",
        "root_isolation_inconclusive",
        "route_disagreement",
    }
)
_LOCALIZATION_CLASSES: Final[frozenset[str]] = frozenset(
    {
        "unique_maximum",
        "nonunique_or_interval_inconclusive",
    }
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


def _exact_nonnegative(value: object, *, name: str) -> Fraction:
    if type(value) is int:
        answer = Q(value)
    elif type(value) is Fraction:
        answer = value
    else:
        raise TypeError(f"{name} must be an exact Fraction or built-in int")
    if answer < 0:
        raise ValueError(f"{name} must be nonnegative")
    return answer


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


def _rational(value: object, *, name: str) -> Fraction:
    if type(value) is int:
        return Q(value)
    if type(value) is Fraction:
        return value
    raise TypeError(f"{name} must be an exact Fraction or built-in int")


def _parse_segment(value: object, *, name: str) -> Segment:
    raw = _fixed_sequence(value, length=5, name=name)
    exact = tuple(
        _rational(item, name=f"{name}[{index}]") for index, item in enumerate(raw)
    )
    if exact[4] <= 0:
        raise ValueError(f"{name} width must be positive")
    return exact  # type: ignore[return-value]


def _rational_hermite_coefficients(segment: Segment) -> Cubic:
    left, left_rhs, right, right_rhs, width = segment
    left_slope = width * left_rhs
    right_slope = width * right_rhs
    return (
        left,
        left_slope,
        -3 * left - 2 * left_slope + 3 * right - right_slope,
        2 * left + left_slope - 2 * right + right_slope,
    )


def _feed_coefficient_hash(hasher: object, ordinal: int, coefficients: Cubic) -> None:
    line = (
        f"{ordinal}|"
        + "|".join(f"{item.numerator}/{item.denominator}" for item in coefficients)
        + "\n"
    )
    hasher.update(line.encode("ascii"))  # type: ignore[attr-defined]


def _combined_coefficient_hash(outer_digest: str, finest_digest: str) -> str:
    return sha256(
        _COMBINED_HASH_DOMAIN + f"{outer_digest}\n{finest_digest}\n".encode("ascii")
    ).hexdigest()


def _metadata_text(metadata: tuple[tuple[str, object], ...]) -> str:
    parts: list[str] = []
    for key, value in metadata:
        if not isinstance(key, str):
            raise TypeError("metadata keys must be strings")
        if isinstance(value, bool) or not isinstance(value, (str, int)):
            raise TypeError("complete-C metadata values must be str or int")
        parts.append(f"{key}={value}")
    return ",".join(parts)


def _survivor_key(item: object) -> tuple[object, ...]:
    return (
        item.polynomial_ordinal,  # type: ignore[attr-defined]
        item.location,  # type: ignore[attr-defined]
        item.location_ordinal,  # type: ignore[attr-defined]
        tuple(item.metadata),  # type: ignore[attr-defined]
    )


def _survivor_key_digest(keys: Sequence[tuple[object, ...]]) -> str:
    hasher = sha256(_SURVIVOR_KEY_HASH_DOMAIN)
    for ordinal, location, location_ordinal, metadata in keys:
        hasher.update(
            (
                f"{ordinal}|{location}|{location_ordinal}|{_metadata_text(metadata)}\n"  # type: ignore[arg-type]
            ).encode("ascii")
        )
    return hasher.hexdigest()


def _intervals_overlap(
    left_lower: Fraction,
    left_upper: Fraction,
    right_lower: Fraction,
    right_upper: Fraction,
) -> bool:
    return max(left_lower, right_lower) <= min(left_upper, right_upper)


def _candidate_intervals_overlap(left: object, right: object) -> bool:
    return _intervals_overlap(
        left.parameter_lower,  # type: ignore[attr-defined]
        left.parameter_upper,  # type: ignore[attr-defined]
        right.parameter_lower,  # type: ignore[attr-defined]
        right.parameter_upper,  # type: ignore[attr-defined]
    ) and _intervals_overlap(
        left.absolute_lower,  # type: ignore[attr-defined]
        left.absolute_upper,  # type: ignore[attr-defined]
        right.absolute_lower,  # type: ignore[attr-defined]
        right.absolute_upper,  # type: ignore[attr-defined]
    )


def _primary_stationary_count_digest(
    cubics: Sequence[LocalCubic],
    *,
    refinement_depth: int,
) -> str:
    digest = sha256(_INDEPENDENT_STATIONARY_COUNT_DOMAIN)
    for ordinal, cubic in enumerate(cubics):
        roots = _primary_stationary_intervals(
            cubic.coefficients,
            depth=refinement_depth,
        )
        digest.update(f"{ordinal}|{len(roots)}\n".encode("ascii"))
    return digest.hexdigest()


def _row_line(
    ordinal: int,
    outer: Segment,
    medium: Sequence[Segment],
    fine: Sequence[Segment],
) -> str:
    values: list[Fraction] = []
    for segment in (outer, *medium, *fine):
        values.extend(segment)
    encoded = "|".join(f"{item.numerator}/{item.denominator}" for item in values)
    return f"{ordinal}|{encoded}\n"


@dataclass(frozen=True, slots=True)
class RationalCompleteCDifferenceEvidence:
    """Compact dual-route enclosure of one rational complete-C family."""

    lower: Fraction
    upper: Fraction
    polynomial_count: int
    candidate_count: int
    survivor_count: int
    localization_classification: str
    coefficient_stream_sha256: str
    survivor_key_stream_sha256: str
    primary_stationary_count_stream_sha256: str
    independent_stationary_count_stream_sha256: str
    primary_evaluator_id: str
    independent_evaluator_id: str
    maximum_candidates: int
    routes_agree: bool = True

    def __post_init__(self) -> None:
        lower = _exact_nonnegative(self.lower, name="lower")
        upper = _exact_nonnegative(self.upper, name="upper")
        if lower > upper:
            raise ValueError("difference lower bound exceeds upper bound")
        object.__setattr__(self, "lower", lower)
        object.__setattr__(self, "upper", upper)
        _positive_integer(self.polynomial_count, name="polynomial_count")
        candidate_count = _positive_integer(
            self.candidate_count, name="candidate_count"
        )
        survivor_count = _positive_integer(self.survivor_count, name="survivor_count")
        _positive_integer(self.maximum_candidates, name="maximum_candidates")
        if candidate_count > self.maximum_candidates:
            raise ValueError("candidate count exceeds the declared ceiling")
        if (
            not 2 * self.polynomial_count
            <= candidate_count
            <= 4 * self.polynomial_count
        ):
            raise ValueError("candidate count is incompatible with cubic extrema")
        if survivor_count > candidate_count:
            raise ValueError("survivor count exceeds candidate count")
        if self.localization_classification not in _LOCALIZATION_CLASSES:
            raise ValueError("unknown localization classification")
        if self.localization_classification == "unique_maximum":
            if survivor_count != 1:
                raise ValueError("unique maximum must retain one survivor")
        elif survivor_count < 2:
            raise ValueError("nonunique localization must retain co-survivors")
        _sha256_digest(self.coefficient_stream_sha256, name="coefficient stream hash")
        _sha256_digest(self.survivor_key_stream_sha256, name="survivor key stream hash")
        _sha256_digest(
            self.primary_stationary_count_stream_sha256,
            name="primary stationary-count hash",
        )
        _sha256_digest(
            self.independent_stationary_count_stream_sha256,
            name="independent stationary-count hash",
        )
        if (
            self.primary_stationary_count_stream_sha256
            != self.independent_stationary_count_stream_sha256
        ):
            raise ValueError("stationary-count digests disagree")
        if self.primary_evaluator_id != PRIMARY_LOCALIZER_ID:
            raise ValueError("unexpected primary localizer identifier")
        if self.independent_evaluator_id != INDEPENDENT_LOCALIZER_ID:
            raise ValueError("unexpected independent localizer identifier")
        if self.routes_agree is not True:
            raise ValueError("compact difference evidence cannot record disagreement")


@dataclass(frozen=True, slots=True)
class RationalCompleteCAdmissionEvidence:
    """Completed dual-route rational complete-C assessment for one channel."""

    d01: RationalCompleteCDifferenceEvidence
    d12: RationalCompleteCDifferenceEvidence
    decision: TDG6ChannelAdmission
    row_count: int
    expected_row_count: int
    row_stream_sha256: str
    combined_coefficient_stream_sha256: str
    sufficient_pass_left: Fraction
    sufficient_pass_right: Fraction
    sufficient_contraction_pass: bool
    sufficient_contraction_failure: bool
    threshold_inconclusive: bool
    maximum_candidates_D01: int
    maximum_candidates_D12: int
    refinement_depth: int
    schema_version: int = TDG11_RATIONAL_COMPLETE_C_SCHEMA_VERSION
    evaluator_id: str = TDG11_RATIONAL_COMPLETE_C_EVALUATOR_ID
    absolute_tolerance_used: bool = False
    physical_signal_used_for_normalization: bool = False
    declared_cubic_is_exact_PDE_history: bool = False
    production_authority: bool = False

    def __post_init__(self) -> None:
        if not isinstance(self.d01, RationalCompleteCDifferenceEvidence):
            raise TypeError("d01 must be RationalCompleteCDifferenceEvidence")
        if not isinstance(self.d12, RationalCompleteCDifferenceEvidence):
            raise TypeError("d12 must be RationalCompleteCDifferenceEvidence")
        if not isinstance(self.decision, TDG6ChannelAdmission):
            raise TypeError("decision must be TDG6ChannelAdmission")
        row_count = _positive_integer(self.row_count, name="row_count")
        expected = _positive_integer(self.expected_row_count, name="expected_row_count")
        if row_count != expected:
            raise ValueError("row_count must equal expected_row_count")
        _positive_integer(
            self.refinement_depth, name="refinement_depth", allow_zero=True
        )
        if self.d01.polynomial_count != 2 * row_count:
            raise ValueError("D01 polynomial count must equal twice row_count")
        if self.d12.polynomial_count != 4 * row_count:
            raise ValueError("D12 polynomial count must equal four times row_count")
        if self.d01.maximum_candidates != self.maximum_candidates_D01:
            raise ValueError("D01 candidate ceiling differs from the request")
        if self.d12.maximum_candidates != self.maximum_candidates_D12:
            raise ValueError("D12 candidate ceiling differs from the request")
        _sha256_digest(self.row_stream_sha256, name="row stream hash")
        combined = _sha256_digest(
            self.combined_coefficient_stream_sha256,
            name="combined coefficient hash",
        )
        expected_combined = _combined_coefficient_hash(
            self.d01.coefficient_stream_sha256,
            self.d12.coefficient_stream_sha256,
        )
        if combined != expected_combined:
            raise ValueError("combined coefficient hash differs from component hashes")
        left = _exact_nonnegative(
            self.sufficient_pass_left, name="sufficient_pass_left"
        )
        right = _exact_nonnegative(
            self.sufficient_pass_right, name="sufficient_pass_right"
        )
        object.__setattr__(self, "sufficient_pass_left", left)
        object.__setattr__(self, "sufficient_pass_right", right)
        if left != TDG11_ORDER_SQUARED_MULTIPLIER * self.d12.upper**2:
            raise ValueError("sufficient-pass left side differs from 8 U12^2")
        if right != self.d01.lower**2:
            raise ValueError("sufficient-pass right side differs from L01^2")
        exact_zero = self.d01.upper == 0 and self.d12.upper == 0
        expected_pass = (not exact_zero) and left <= right
        expected_failure = (
            (not exact_zero)
            and (not expected_pass)
            and (TDG11_ORDER_SQUARED_MULTIPLIER * self.d12.lower**2 > self.d01.upper**2)
        )
        expected_inconclusive = (
            (not exact_zero) and (not expected_pass) and (not expected_failure)
        )
        if self.sufficient_contraction_pass is not expected_pass:
            raise ValueError("sufficient-pass flag differs from the certified bounds")
        if self.sufficient_contraction_failure is not expected_failure:
            raise ValueError(
                "sufficient-failure flag differs from the certified bounds"
            )
        if self.threshold_inconclusive is not expected_inconclusive:
            raise ValueError(
                "threshold-inconclusive flag differs from the certified bounds"
            )
        expected_decision = classify_tdg6_channel(
            CertifiedMagnitudeInterval(self.d01.lower, self.d01.upper),
            CertifiedMagnitudeInterval(self.d12.lower, self.d12.upper),
        )
        if self.decision != expected_decision:
            raise ValueError(
                "channel decision differs from classify_tdg6_channel on the "
                "certified complete-C intervals"
            )
        if (
            type(self.schema_version) is not int
            or self.schema_version != TDG11_RATIONAL_COMPLETE_C_SCHEMA_VERSION
            or self.evaluator_id != TDG11_RATIONAL_COMPLETE_C_EVALUATOR_ID
            or self.absolute_tolerance_used is not False
            or self.physical_signal_used_for_normalization is not False
            or self.declared_cubic_is_exact_PDE_history is not False
            or self.production_authority is not False
        ):
            raise ValueError("rational complete-C evidence crossed its scope")


@dataclass(frozen=True, slots=True)
class RationalCompleteCClosedEvidence:
    """Typed fail-closed record for disagreement or exhausted exact work."""

    reason: str
    level: str
    row_count: int
    maximum_candidates_D01: int
    maximum_candidates_D12: int
    refinement_depth: int
    detail: str

    def __post_init__(self) -> None:
        if self.reason not in _CLOSED_REASONS:
            raise ValueError("unknown rational complete-C closed reason")
        if self.level not in {"D01", "D12"}:
            raise ValueError("closed evidence level must be D01 or D12")
        _positive_integer(self.row_count, name="row_count", allow_zero=True)
        _positive_integer(self.maximum_candidates_D01, name="maximum_candidates_D01")
        _positive_integer(self.maximum_candidates_D12, name="maximum_candidates_D12")
        _positive_integer(
            self.refinement_depth, name="refinement_depth", allow_zero=True
        )
        if not isinstance(self.detail, str) or not self.detail:
            raise ValueError("closed evidence detail must be a nonempty string")


class RationalCompleteCResourceExhausted(RuntimeError):
    """Raised when a localizer ceiling or root-isolation budget is exhausted."""

    def __init__(self, evidence: RationalCompleteCClosedEvidence) -> None:
        if evidence.reason not in {
            "candidate_ceiling_exhausted",
            "root_isolation_inconclusive",
        }:
            raise ValueError("resource evidence must use a resource reason")
        self.evidence = evidence
        super().__init__(f"tdg11_rational_complete_c_exhausted:{evidence.reason}")


class RationalCompleteCRouteDisagreement(RuntimeError):
    """Raised when the two exact localizers fail a required agreement check."""

    def __init__(self, evidence: RationalCompleteCClosedEvidence) -> None:
        if evidence.reason != "route_disagreement":
            raise ValueError("disagreement evidence must use route_disagreement")
        self.evidence = evidence
        super().__init__("tdg11_rational_complete_c_closed:route_disagreement")


def _complete_c_metadata(
    *, level: str, row_index: int, subinterval: int
) -> dict[str, object]:
    return {
        "level": level,
        "row_index": row_index,
        "subinterval": subinterval,
    }


def _route_agreement_facts(
    first: object,
    second: object,
    *,
    primary_stationary_digest: str,
    independent_stationary_digest: str,
) -> dict[str, bool]:
    first_items = {
        _survivor_key(item): item
        for item in first.candidates  # type: ignore[attr-defined]
    }
    second_items = {
        _survivor_key(item): item
        for item in second.candidates  # type: ignore[attr-defined]
    }
    same_keys = (
        len(first_items) == len(first.candidates)  # type: ignore[attr-defined]
        and len(second_items) == len(second.candidates)  # type: ignore[attr-defined]
        and tuple(sorted(first_items)) == tuple(sorted(second_items))
    )
    survivor_overlaps = same_keys and all(
        _candidate_intervals_overlap(first_items[key], second_items[key])
        for key in first_items
    )
    global_overlap = _intervals_overlap(
        first.global_absolute_lower,  # type: ignore[attr-defined]
        first.global_absolute_upper,  # type: ignore[attr-defined]
        second.global_absolute_lower,  # type: ignore[attr-defined]
        second.global_absolute_upper,  # type: ignore[attr-defined]
    )
    return {
        "localization_classification_agrees": (
            first.classification == second.classification  # type: ignore[attr-defined]
            and first.classification in _LOCALIZATION_CLASSES  # type: ignore[attr-defined]
        ),
        "polynomial_count_agrees": (
            first.polynomial_count == second.polynomial_count  # type: ignore[attr-defined]
        ),
        "candidate_count_agrees": (
            first.candidate_count == second.candidate_count  # type: ignore[attr-defined]
        ),
        "stationary_count_digest_agrees": (
            primary_stationary_digest == independent_stationary_digest
        ),
        "survivor_keys_agree": same_keys,
        "survivor_intervals_overlap": survivor_overlaps,
        "global_intervals_overlap": global_overlap,
    }


def _facts_text(facts: dict[str, bool]) -> str:
    return ",".join(
        f"{name}={'true' if value else 'false'}" for name, value in facts.items()
    )


def _closed_evidence(
    *,
    reason: str,
    level: str,
    row_count: int,
    maximum_candidates_D01: int,
    maximum_candidates_D12: int,
    refinement_depth: int,
    detail: str,
) -> RationalCompleteCClosedEvidence:
    return RationalCompleteCClosedEvidence(
        reason=reason,
        level=level,
        row_count=row_count,
        maximum_candidates_D01=maximum_candidates_D01,
        maximum_candidates_D12=maximum_candidates_D12,
        refinement_depth=refinement_depth,
        detail=detail,
    )


def _localize_family(
    coefficients: Sequence[tuple[Cubic, dict[str, object]]],
    *,
    level: str,
    row_count: int,
    coefficient_digest: str,
    maximum_candidates: int,
    maximum_candidates_D01: int,
    maximum_candidates_D12: int,
    refinement_depth: int,
) -> RationalCompleteCDifferenceEvidence:
    primary_cubics = tuple(
        LocalCubic(item, metadata) for item, metadata in coefficients
    )
    independent_cubics = tuple(
        IndependentLocalCubicV2(item, metadata) for item, metadata in coefficients
    )
    try:
        primary = localize_absolute_maximum(
            primary_cubics,
            maximum_candidates=maximum_candidates,
            refinement_depth=refinement_depth,
        )
    except RuntimeError as exc:
        if str(exc) != "candidate_ceiling_exhausted":
            raise
        raise RationalCompleteCResourceExhausted(
            _closed_evidence(
                reason="candidate_ceiling_exhausted",
                level=level,
                row_count=row_count,
                maximum_candidates_D01=maximum_candidates_D01,
                maximum_candidates_D12=maximum_candidates_D12,
                refinement_depth=refinement_depth,
                detail="primary_localizer",
            )
        ) from exc
    try:
        independent = localize_absolute_maximum_independently_v2(
            independent_cubics,
            maximum_candidates=maximum_candidates,
        )
    except RootIsolationInconclusive as exc:
        raise RationalCompleteCResourceExhausted(
            _closed_evidence(
                reason="root_isolation_inconclusive",
                level=level,
                row_count=row_count,
                maximum_candidates_D01=maximum_candidates_D01,
                maximum_candidates_D12=maximum_candidates_D12,
                refinement_depth=refinement_depth,
                detail=f"{exc.reason}:{exc.detail}",
            )
        ) from exc
    except RuntimeError as exc:
        if str(exc) != "candidate_ceiling_exhausted":
            raise
        raise RationalCompleteCResourceExhausted(
            _closed_evidence(
                reason="candidate_ceiling_exhausted",
                level=level,
                row_count=row_count,
                maximum_candidates_D01=maximum_candidates_D01,
                maximum_candidates_D12=maximum_candidates_D12,
                refinement_depth=refinement_depth,
                detail="independent_localizer_v2",
            )
        ) from exc

    if primary.tolerance_used is not False:
        raise ValueError("primary exact localizer used an absolute tolerance")
    if (
        independent.tolerance_used is not False
        or independent.boundary_clipping_used is not False
    ):
        raise ValueError(
            "independent exact localizer used a tolerance or boundary clipping"
        )

    primary_stationary = _primary_stationary_count_digest(
        primary_cubics, refinement_depth=refinement_depth
    )
    independent_stationary = independent.stationary_count_stream_sha256
    facts = _route_agreement_facts(
        primary,
        independent,
        primary_stationary_digest=primary_stationary,
        independent_stationary_digest=independent_stationary,
    )
    if not all(facts.values()):
        raise RationalCompleteCRouteDisagreement(
            _closed_evidence(
                reason="route_disagreement",
                level=level,
                row_count=row_count,
                maximum_candidates_D01=maximum_candidates_D01,
                maximum_candidates_D12=maximum_candidates_D12,
                refinement_depth=refinement_depth,
                detail=_facts_text(facts),
            )
        )

    certified_lower = max(
        primary.global_absolute_lower, independent.global_absolute_lower
    )
    certified_upper = min(
        primary.global_absolute_upper, independent.global_absolute_upper
    )
    survivor_keys = tuple(sorted(_survivor_key(item) for item in primary.candidates))
    return RationalCompleteCDifferenceEvidence(
        lower=certified_lower,
        upper=certified_upper,
        polynomial_count=primary.polynomial_count,
        candidate_count=primary.candidate_count,
        survivor_count=len(primary.candidates),
        localization_classification=primary.classification,
        coefficient_stream_sha256=coefficient_digest,
        survivor_key_stream_sha256=_survivor_key_digest(survivor_keys),
        primary_stationary_count_stream_sha256=primary_stationary,
        independent_stationary_count_stream_sha256=independent_stationary,
        primary_evaluator_id=primary.evaluator_id,
        independent_evaluator_id=independent.evaluator_id,
        maximum_candidates=maximum_candidates,
    )


def _contraction_flags(
    *,
    outer_lower: Fraction,
    outer_upper: Fraction,
    finest_lower: Fraction,
    finest_upper: Fraction,
) -> tuple[Fraction, Fraction, bool, bool, bool]:
    left = TDG11_ORDER_SQUARED_MULTIPLIER * finest_upper**2
    right = outer_lower**2
    exact_zero = outer_upper == 0 and finest_upper == 0
    if exact_zero:
        return left, right, False, False, False
    if left <= right:
        return left, right, True, False, False
    if TDG11_ORDER_SQUARED_MULTIPLIER * finest_lower**2 > outer_upper**2:
        return left, right, False, True, False
    return left, right, False, False, True


def assess_rational_complete_c_rows(
    rows: Iterable[object],
    *,
    expected_row_count: int,
    maximum_candidates_D01: int,
    maximum_candidates_D12: int,
    refinement_depth: int = 160,
) -> RationalCompleteCAdmissionEvidence:
    """Enclose streamed rational complete-C D01/D12 maxima and classify.

    The input iterable is consumed once and must yield exactly
    ``expected_row_count`` rows.  Compact evidence never retains the
    per-candidate localizer output.  Work is bounded by the explicit D01/D12
    candidate ceilings and by the primary localizer refinement depth.
    """

    expected_rows = _positive_integer(expected_row_count, name="expected_row_count")
    depth_limit = _positive_integer(
        refinement_depth, name="refinement_depth", allow_zero=True
    )
    d01_ceiling = _positive_integer(
        maximum_candidates_D01, name="maximum_candidates_D01"
    )
    d12_ceiling = _positive_integer(
        maximum_candidates_D12, name="maximum_candidates_D12"
    )
    row_hasher = sha256(_ROW_HASH_DOMAIN)
    d01_hasher = sha256(_COEFFICIENT_HASH_DOMAIN)
    d12_hasher = sha256(_COEFFICIENT_HASH_DOMAIN)
    d01_coefficients: list[tuple[Cubic, dict[str, object]]] = []
    d12_coefficients: list[tuple[Cubic, dict[str, object]]] = []
    row_count = 0

    for row_index, row in enumerate(rows):
        if row_index >= expected_rows:
            raise ValueError("rows exceed expected_row_count")
        raw_outer, raw_medium, raw_fine = _fixed_sequence(
            row, length=3, name=f"rows[{row_index}]"
        )
        medium = _fixed_sequence(raw_medium, length=2, name=f"rows[{row_index}].medium")
        fine = _fixed_sequence(raw_fine, length=4, name=f"rows[{row_index}].fine")
        outer = _parse_segment(raw_outer, name=f"rows[{row_index}].outer")
        medium_segments = tuple(
            _parse_segment(item, name=f"rows[{row_index}].medium[{index}]")
            for index, item in enumerate(medium)
        )
        fine_segments = tuple(
            _parse_segment(item, name=f"rows[{row_index}].fine[{index}]")
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

        outer_coefficients = _rational_hermite_coefficients(outer)
        medium_coefficients = tuple(
            _rational_hermite_coefficients(item) for item in medium_segments
        )
        fine_coefficients = tuple(
            _rational_hermite_coefficients(item) for item in fine_segments
        )
        for half, child in enumerate(medium_coefficients):
            difference = subtract_exact_cubics(
                restrict_exact_cubic_to_half(outer_coefficients, half=half),
                child,
            )
            _feed_coefficient_hash(d01_hasher, len(d01_coefficients), difference)
            d01_coefficients.append(
                (
                    difference,
                    _complete_c_metadata(
                        level="D01", row_index=row_index, subinterval=half
                    ),
                )
            )
        for medium_index, parent in enumerate(medium_coefficients):
            for half in (0, 1):
                child_index = 2 * medium_index + half
                difference = subtract_exact_cubics(
                    restrict_exact_cubic_to_half(parent, half=half),
                    fine_coefficients[child_index],
                )
                _feed_coefficient_hash(d12_hasher, len(d12_coefficients), difference)
                d12_coefficients.append(
                    (
                        difference,
                        _complete_c_metadata(
                            level="D12",
                            row_index=row_index,
                            subinterval=child_index,
                        ),
                    )
                )
        row_count += 1

    if row_count != expected_rows:
        raise ValueError("rows must contain exactly expected_row_count refinement rows")

    d01 = _localize_family(
        d01_coefficients,
        level="D01",
        row_count=row_count,
        coefficient_digest=d01_hasher.hexdigest(),
        maximum_candidates=d01_ceiling,
        maximum_candidates_D01=d01_ceiling,
        maximum_candidates_D12=d12_ceiling,
        refinement_depth=depth_limit,
    )
    d12 = _localize_family(
        d12_coefficients,
        level="D12",
        row_count=row_count,
        coefficient_digest=d12_hasher.hexdigest(),
        maximum_candidates=d12_ceiling,
        maximum_candidates_D01=d01_ceiling,
        maximum_candidates_D12=d12_ceiling,
        refinement_depth=depth_limit,
    )
    (
        pass_left,
        pass_right,
        sufficient_pass,
        sufficient_failure,
        inconclusive,
    ) = _contraction_flags(
        outer_lower=d01.lower,
        outer_upper=d01.upper,
        finest_lower=d12.lower,
        finest_upper=d12.upper,
    )
    decision = classify_tdg6_channel(
        CertifiedMagnitudeInterval(d01.lower, d01.upper),
        CertifiedMagnitudeInterval(d12.lower, d12.upper),
    )
    return RationalCompleteCAdmissionEvidence(
        d01=d01,
        d12=d12,
        decision=decision,
        row_count=row_count,
        expected_row_count=expected_rows,
        row_stream_sha256=row_hasher.hexdigest(),
        combined_coefficient_stream_sha256=_combined_coefficient_hash(
            d01.coefficient_stream_sha256,
            d12.coefficient_stream_sha256,
        ),
        sufficient_pass_left=pass_left,
        sufficient_pass_right=pass_right,
        sufficient_contraction_pass=sufficient_pass,
        sufficient_contraction_failure=sufficient_failure,
        threshold_inconclusive=inconclusive,
        maximum_candidates_D01=d01_ceiling,
        maximum_candidates_D12=d12_ceiling,
        refinement_depth=depth_limit,
    )


__all__ = [
    "RationalCompleteCAdmissionEvidence",
    "RationalCompleteCClosedEvidence",
    "RationalCompleteCDifferenceEvidence",
    "RationalCompleteCResourceExhausted",
    "RationalCompleteCRouteDisagreement",
    "TDG11_RATIONAL_COMPLETE_C_EVALUATOR_ID",
    "TDG11_RATIONAL_COMPLETE_C_SCHEMA_VERSION",
    "assess_rational_complete_c_rows",
]
