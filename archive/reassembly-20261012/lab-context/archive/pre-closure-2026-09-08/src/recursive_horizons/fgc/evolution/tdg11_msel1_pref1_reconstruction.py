"""Independent recorded-RK family validator for TDG11-MSEL1-PREF1.

This module inspects already-collected ``StepProposal`` records.  It never
calls a proposer, RHS, projector, or guard, never accepts a state, and never
emits a scientific pass label.  Byte-for-byte replay of the declared RK
schedule is a coherence check on the recorded object.  Exact accumulation and
embedded-defect reductions are formed here from the declared tableaux and the
recorded binary64 values; they are not copied from a runner or from the
production reconstruction owner.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
from math import isfinite
from typing import Final, Sequence

import numpy as np

from .numerical_engine import (
    COMPARATOR_METHOD,
    PRIMARY_METHOD,
    EvolutionRHS,
    EvolutionState,
    StepProposal,
)
from .tdg11_compensated_rk import (
    ARITHMETIC_ID as COMPENSATED_ARITHMETIC_ID,
    CompensatedArithmeticStop,
    compensated_update,
)


ORIGINAL_ARITHMETIC_ID: Final[str] = "tdg11_recorded_legacy_binary64_rk_v1"
_ARITHMETIC_IDS: Final[frozenset[str]] = frozenset(
    {ORIGINAL_ARITHMETIC_ID, COMPENSATED_ARITHMETIC_ID}
)
_OWNED: Final[slice] = slice(1, -4)
_BLOCKS: Final[tuple[str, str, str]] = ("u", "p", "q")
_DERIVATIVES: Final[tuple[str, str, str]] = ("du", "dp", "dq")
_FIELDS: Final[tuple[str, ...]] = ("alpha", "v", "lambda", "R", "phi", "chi")
_CHANNELS: Final[tuple[str, ...]] = tuple(
    f"{block}:{field}" for block in _BLOCKS for field in _FIELDS
)
_PATH_LENGTHS: Final[tuple[int, int, int]] = (1, 2, 4)
_RK4_NAMES: Final[tuple[str, ...]] = (
    "rk4_k1",
    "rk4_k2",
    "rk4_k3",
    "rk4_k4",
    "candidate_endpoint",
)
_SSP_NAMES: Final[tuple[str, ...]] = (
    "ssprk3_s0",
    "ssprk3_s1",
    "ssprk3_s2",
    "candidate_endpoint",
)
_RK4_WEIGHTS: Final[tuple[Fraction, ...]] = (
    Fraction(1, 6),
    Fraction(1, 3),
    Fraction(1, 3),
    Fraction(1, 6),
)
_SSP_WEIGHTS: Final[tuple[Fraction, ...]] = (
    Fraction(1, 6),
    Fraction(1, 6),
    Fraction(2, 3),
)
_HASH_DOMAIN: Final[bytes] = b"TDG11-MSEL1-RECORDED-FAMILY-v1\n"

Segment = tuple[Fraction, Fraction, Fraction, Fraction, Fraction]
PathSegments = tuple[Segment, ...]
RefinementRow = tuple[Segment, tuple[Segment, ...], tuple[Segment, ...]]


class IndependentRecordedFamilyError(ValueError):
    """The recorded 1/2/4 RK family is not a coherent reconstruction object."""


def _fail(message: str) -> None:
    raise IndependentRecordedFamilyError(message)


def _exact(value: object) -> Fraction:
    return Fraction.from_float(float(value))


def _positive_int(value: object, *, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{name} must be an integer")
    if value < 1:
        _fail(f"{name} must be positive")
    return value


def _bits(value: object) -> bytes:
    if type(value) not in (float, np.float64):
        _fail("times and widths must be binary64 values, not coercible labels")
    number = float(value)
    if not isfinite(number):
        _fail("times and widths must be finite binary64 values")
    return np.float64(number).tobytes()


def _same_bits(left: object, right: object) -> bool:
    return _bits(left) == _bits(right)


def _require_native(array: object, *, name: str) -> np.ndarray:
    if type(array) is not np.ndarray or array.dtype != np.dtype(np.float64):
        _fail(f"{name} must be a finite native binary64 array")
    if array.ndim != 2 or not array.flags.c_contiguous:
        _fail(f"{name} must be a finite native binary64 array")
    if not np.all(np.isfinite(array)):
        _fail(f"{name} must be a finite native binary64 array")
    return array


def _owned_bytes(array: np.ndarray) -> bytes:
    return np.ascontiguousarray(array[_OWNED]).tobytes()


def _pack_state(state: EvolutionState, *, owned: bool) -> bytes:
    pack = _owned_bytes if owned else (lambda item: item.tobytes())
    return b"".join(pack(getattr(state, name)) for name in _BLOCKS)


def _pack_rhs(rhs: EvolutionRHS, *, owned: bool) -> bytes:
    pack = _owned_bytes if owned else (lambda item: item.tobytes())
    return b"".join(pack(getattr(rhs, name)) for name in _DERIVATIVES)


def _check_state(state: object, *, owned_row_count: int, name: str) -> EvolutionState:
    if not isinstance(state, EvolutionState):
        _fail(f"{name} must be EvolutionState")
    if state.shape[1] != len(_FIELDS):
        _fail(f"{name} must contain six fields")
    if state.shape[0] - 5 != owned_row_count:
        _fail("owned row count must equal point_count - 5")
    for block in _BLOCKS:
        _require_native(getattr(state, block), name=f"{name}.{block}")
    return state


def _check_rhs(rhs: object, *, owned_row_count: int, name: str) -> EvolutionRHS:
    if not isinstance(rhs, EvolutionRHS):
        _fail(f"{name} must be EvolutionRHS")
    if rhs.shape[1] != len(_FIELDS):
        _fail(f"{name} must contain six fields")
    if rhs.shape[0] - 5 != owned_row_count:
        _fail("owned row count must equal point_count - 5")
    for derivative in _DERIVATIVES:
        _require_native(getattr(rhs, derivative), name=f"{name}.{derivative}")
    return rhs


def _freeze_paths(paths: object) -> tuple[tuple[StepProposal, ...], ...]:
    if (
        isinstance(paths, (str, bytes))
        or not isinstance(paths, Sequence)
        or len(paths) != 3
    ):
        _fail("paths must be a tuple of three StepProposal tuples")
    frozen: list[tuple[StepProposal, ...]] = []
    for path, length in zip(paths, _PATH_LENGTHS, strict=True):
        if isinstance(path, (str, bytes)) or not isinstance(path, Sequence):
            _fail("paths must be a tuple of three StepProposal tuples")
        items = tuple(path)
        if len(items) != length or any(
            not isinstance(item, StepProposal) for item in items
        ):
            _fail("paths must contain StepProposal tuples of lengths 1, 2, and 4")
        frozen.append(items)
    return tuple(frozen)


def _feed_array(hasher: object, array: np.ndarray) -> None:
    hasher.update(np.asarray(array.shape, dtype=np.int64).tobytes())  # type: ignore[attr-defined]
    hasher.update(np.ascontiguousarray(array, dtype=np.float64).tobytes())  # type: ignore[attr-defined]


def _family_digest(
    paths: tuple[tuple[StepProposal, ...], ...],
    *,
    method: str,
    arithmetic_id: str,
    owned_row_count: int,
) -> str:
    hasher = sha256(_HASH_DOMAIN)
    hasher.update(f"{method}\n{arithmetic_id}\n{owned_row_count}\n".encode("ascii"))
    for path in paths:
        for proposal in path:
            hasher.update(proposal.method.encode("ascii") + b"\n")
            hasher.update(_bits(proposal.initial_time))
            hasher.update(_bits(proposal.final_time))
            for state in (proposal.initial_state, proposal.candidate_state):
                for block in _BLOCKS:
                    _feed_array(hasher, getattr(state, block))
            for stage in proposal.stages:
                hasher.update(stage.stage_name.encode("ascii") + b"\n")
                hasher.update(_bits(stage.time))
                for block in _BLOCKS:
                    _feed_array(hasher, getattr(stage.state, block))
                for derivative in _DERIVATIVES:
                    _feed_array(hasher, getattr(stage.rhs, derivative))
    return hasher.hexdigest()


def _compensated_state(
    base: EvolutionState,
    terms: tuple[tuple[Fraction, EvolutionRHS], ...],
) -> EvolutionState:
    try:
        return EvolutionState(
            *(
                compensated_update(
                    getattr(base, block),
                    tuple(
                        (coefficient, getattr(rhs, derivative))
                        for coefficient, rhs in terms
                    ),
                )
                for block, derivative in zip(_BLOCKS, _DERIVATIVES, strict=True)
            )
        )
    except CompensatedArithmeticStop as error:
        raise IndependentRecordedFamilyError("compensated replay failed") from error


def _legacy_rk4(
    base: EvolutionState, h: float, stages: tuple[EvolutionRHS, ...]
) -> tuple[EvolutionState, ...]:
    k1, k2, k3, k4 = stages
    return (
        base.linear_combination(((h / 2.0, k1),)),
        base.linear_combination(((h / 2.0, k2),)),
        base.linear_combination(((h, k3),)),
        base.linear_combination(
            ((h / 6.0, k1), (h / 3.0, k2), (h / 3.0, k3), (h / 6.0, k4))
        ),
    )


def _legacy_ssp(
    base: EvolutionState, h: float, stages: tuple[EvolutionRHS, ...]
) -> tuple[EvolutionState, ...]:
    k0, k1, k2 = stages
    y1 = base.linear_combination(((h, k0),))
    euler_y1 = y1.linear_combination(((h, k1),))
    y2 = EvolutionState(
        0.75 * base.u + 0.25 * euler_y1.u,
        0.75 * base.p + 0.25 * euler_y1.p,
        0.75 * base.q + 0.25 * euler_y1.q,
    )
    euler_y2 = y2.linear_combination(((h, k2),))
    candidate = EvolutionState(
        base.u / 3 + 2 * euler_y2.u / 3,
        base.p / 3 + 2 * euler_y2.p / 3,
        base.q / 3 + 2 * euler_y2.q / 3,
    )
    return y1, y2, candidate


def _compensated_rk4(
    base: EvolutionState, h: Fraction, stages: tuple[EvolutionRHS, ...]
) -> tuple[EvolutionState, ...]:
    k1, k2, k3, k4 = stages
    return (
        _compensated_state(base, ((h * Fraction(1, 2), k1),)),
        _compensated_state(base, ((h * Fraction(1, 2), k2),)),
        _compensated_state(base, ((h * Fraction(1), k3),)),
        _compensated_state(
            base,
            (
                (h * Fraction(1, 6), k1),
                (h * Fraction(1, 3), k2),
                (h * Fraction(1, 3), k3),
                (h * Fraction(1, 6), k4),
            ),
        ),
    )


def _compensated_ssp(
    base: EvolutionState, h: Fraction, stages: tuple[EvolutionRHS, ...]
) -> tuple[EvolutionState, ...]:
    k0, k1, k2 = stages
    return (
        _compensated_state(base, ((h * Fraction(1), k0),)),
        _compensated_state(base, ((h * Fraction(1, 4), k0), (h * Fraction(1, 4), k1))),
        _compensated_state(
            base,
            (
                (h * Fraction(1, 6), k0),
                (h * Fraction(1, 6), k1),
                (h * Fraction(2, 3), k2),
            ),
        ),
    )


def _replay_states(
    proposal: StepProposal, *, method: str, arithmetic_id: str, h: float
) -> tuple[EvolutionState, ...]:
    base = proposal.initial_state
    if method == PRIMARY_METHOD:
        rhs_records = tuple(stage.rhs for stage in proposal.stages[:4])
        if arithmetic_id == ORIGINAL_ARITHMETIC_ID:
            return _legacy_rk4(base, h, rhs_records)
        return _compensated_rk4(base, Fraction.from_float(h), rhs_records)
    rhs_records = tuple(stage.rhs for stage in proposal.stages[:3])
    if arithmetic_id == ORIGINAL_ARITHMETIC_ID:
        return _legacy_ssp(base, h, rhs_records)
    return _compensated_ssp(base, Fraction.from_float(h), rhs_records)


def _validate_proposal(
    proposal: StepProposal,
    *,
    method: str,
    arithmetic_id: str,
    owned_row_count: int,
) -> float:
    if proposal.method != method:
        _fail("proposal method differs from the recorded family method")
    names = _RK4_NAMES if method == PRIMARY_METHOD else _SSP_NAMES
    if len(proposal.stages) != len(names):
        _fail("proposal must retain 5 RK4 or 4 SSPRK3 stage records")
    if tuple(stage.stage_name for stage in proposal.stages) != names:
        _fail("stage names do not match the recorded method")
    start = proposal.initial_time
    final = proposal.final_time
    width = final - start
    if not isfinite(width) or width <= 0.0:
        _fail("step width must be positive")
    middle = start + width / 2.0
    if not (start < middle < final) or not isfinite(middle):
        _fail("stage midpoint is not stage-safe")
    expected_times = (
        (start, middle, middle, final, final)
        if method == PRIMARY_METHOD
        else (start, final, middle, final)
    )
    for stage, expected in zip(proposal.stages, expected_times, strict=True):
        if not _same_bits(stage.time, expected):
            _fail("stage time does not match the recorded lattice")
        _check_state(stage.state, owned_row_count=owned_row_count, name="stage.state")
        _check_rhs(stage.rhs, owned_row_count=owned_row_count, name="stage.rhs")
    _check_state(
        proposal.initial_state, owned_row_count=owned_row_count, name="initial_state"
    )
    _check_state(
        proposal.candidate_state,
        owned_row_count=owned_row_count,
        name="candidate_state",
    )
    if _pack_state(proposal.stages[0].state, owned=False) != _pack_state(
        proposal.initial_state, owned=False
    ):
        _fail("first-stage state differs from the initial state")
    if _pack_state(proposal.stages[-1].state, owned=False) != _pack_state(
        proposal.candidate_state, owned=False
    ):
        _fail("endpoint record does not match the candidate")
    replayed = _replay_states(
        proposal, method=method, arithmetic_id=arithmetic_id, h=width
    )
    recorded = tuple(stage.state for stage in proposal.stages[1:])
    for expected, actual, label in zip(replayed, recorded, names[1:], strict=True):
        if _pack_state(expected, owned=True) != _pack_state(actual, owned=True):
            _fail(f"owned {label} bytes differ from recorded arithmetic replay")
    return width


def _validate_path(
    path: tuple[StepProposal, ...],
    *,
    method: str,
    arithmetic_id: str,
    owned_row_count: int,
    start_bits: bytes,
    final_bits: bytes,
    expected_width: float,
) -> None:
    previous_state: bytes | None = None
    previous_rhs: bytes | None = None
    previous_final: bytes | None = None
    for index, proposal in enumerate(path):
        if index == 0:
            if _bits(proposal.initial_time) != start_bits:
                _fail("path does not share the family initial time")
        elif previous_final != _bits(proposal.initial_time):
            _fail("path times are not contiguous")
        if previous_state is not None and previous_state != _pack_state(
            proposal.initial_state, owned=False
        ):
            _fail("path endpoint states are not contiguous")
        if previous_rhs is not None and previous_rhs != _pack_rhs(
            proposal.stages[0].rhs, owned=False
        ):
            _fail("path endpoint RHS values are not contiguous")
        width = _validate_proposal(
            proposal,
            method=method,
            arithmetic_id=arithmetic_id,
            owned_row_count=owned_row_count,
        )
        if not _same_bits(width, expected_width):
            _fail("path widths are not bitwise uniform")
        previous_state = _pack_state(proposal.candidate_state, owned=False)
        previous_rhs = _pack_rhs(proposal.stages[-1].rhs, owned=False)
        previous_final = _bits(proposal.final_time)
    if previous_final != final_bits:
        _fail("path does not share the family final time")


def _validate_family_contents(
    frozen: tuple[tuple[StepProposal, ...], ...],
    *,
    method: str,
    arithmetic_id: str,
    owned: int,
) -> None:
    sample = frozen[0][0]
    start_bits = _bits(sample.initial_time)
    final_bits = _bits(sample.final_time)
    outer_width = sample.final_time - sample.initial_time
    expected_widths = (outer_width, outer_width / 2.0, outer_width / 4.0)
    if any(width <= 0.0 or not isfinite(width) for width in expected_widths):
        _fail("h, h/2, and h/4 must be positive")
    initial_state = None
    initial_rhs = None
    for path, width in zip(frozen, expected_widths, strict=True):
        _validate_path(
            path,
            method=method,
            arithmetic_id=arithmetic_id,
            owned_row_count=owned,
            start_bits=start_bits,
            final_bits=final_bits,
            expected_width=width,
        )
        if initial_state is None:
            initial_state = _pack_state(path[0].initial_state, owned=False)
            initial_rhs = _pack_rhs(path[0].stages[0].rhs, owned=False)
        elif _pack_state(path[0].initial_state, owned=False) != initial_state:
            _fail("paths do not share a bitwise initial full state")
        elif initial_rhs != _pack_rhs(path[0].stages[0].rhs, owned=False):
            _fail("paths do not share a bitwise initial RHS")


def _require_segment(value: object, *, name: str) -> None:
    if not isinstance(value, tuple) or len(value) != 5:
        _fail(f"{name} must be a five-tuple segment")
    for item in value:
        if type(item) not in (int, Fraction):
            raise TypeError(f"{name} must contain exact rationals")


def _require_bounds(value: object, *, name: str) -> None:
    if not isinstance(value, tuple) or len(value) != 3:
        _fail(f"{name} must be the (outer, medium, fine) triple")
    for item in value:
        if type(item) not in (int, Fraction) or item < 0:
            _fail(f"{name} must contain nonnegative exact rationals")


@dataclass(frozen=True, slots=True)
class IndependentRecordedFamily:
    """Immutable 1/2/4 recorded-RK family after schedule-coherent replay."""

    paths: tuple[tuple[StepProposal, ...], ...]
    method: str
    arithmetic_id: str
    owned_row_count: int
    family_sha256: str

    def __post_init__(self) -> None:
        paths = _freeze_paths(self.paths)
        object.__setattr__(self, "paths", paths)
        if self.method not in (PRIMARY_METHOD, COMPARATOR_METHOD):
            _fail("unrecognized method")
        if self.arithmetic_id not in _ARITHMETIC_IDS:
            _fail("unrecognized arithmetic")
        owned = _positive_int(self.owned_row_count, name="owned_row_count")
        digest = self.family_sha256
        if (
            not isinstance(digest, str)
            or len(digest) != 64
            or any(character not in "0123456789abcdef" for character in digest)
        ):
            _fail("family hash must be a lowercase SHA-256 digest")
        _validate_family_contents(
            paths, method=self.method, arithmetic_id=self.arithmetic_id, owned=owned
        )
        expected = _family_digest(
            paths,
            method=self.method,
            arithmetic_id=self.arithmetic_id,
            owned_row_count=owned,
        )
        if digest != expected:
            _fail("family hash does not match recorded content")


@dataclass(frozen=True, slots=True)
class IndependentChannelReconstruction:
    """Exact rational raw/corrected Hermite rows and per-level defect bounds."""

    channel: str
    row_count: int
    raw_rows: tuple[RefinementRow, ...]
    corrected_rows: tuple[RefinementRow, ...]
    accumulation_bounds: tuple[Fraction, Fraction, Fraction]
    embedded_defect_bounds: tuple[Fraction, Fraction, Fraction]

    def __post_init__(self) -> None:
        if self.channel not in _CHANNELS:
            _fail("channel is not a TDG6 complete-state channel")
        row_count = _positive_int(self.row_count, name="row_count")
        if len(self.raw_rows) != row_count or len(self.corrected_rows) != row_count:
            _fail("row count differs from the reconstructed row stream")
        for name, rows in (
            ("raw_rows", self.raw_rows),
            ("corrected_rows", self.corrected_rows),
        ):
            for index, row in enumerate(rows):
                if not isinstance(row, tuple) or len(row) != 3:
                    _fail(f"{name}[{index}] must be (outer, medium, fine)")
                outer, medium, fine = row
                _require_segment(outer, name=f"{name}[{index}].outer")
                if not isinstance(medium, tuple) or len(medium) != 2:
                    _fail(f"{name}[{index}].medium must contain two segments")
                if not isinstance(fine, tuple) or len(fine) != 4:
                    _fail(f"{name}[{index}].fine must contain four segments")
                for inner_index, segment in enumerate(medium):
                    _require_segment(
                        segment, name=f"{name}[{index}].medium[{inner_index}]"
                    )
                for inner_index, segment in enumerate(fine):
                    _require_segment(
                        segment, name=f"{name}[{index}].fine[{inner_index}]"
                    )
        _require_bounds(self.accumulation_bounds, name="accumulation_bounds")
        _require_bounds(self.embedded_defect_bounds, name="embedded_defect_bounds")


def validate_recorded_family_independently(
    paths: object,
    *,
    method: str,
    arithmetic_id: str,
    owned_row_count: int,
) -> IndependentRecordedFamily:
    """Replay recorded RK arithmetic and freeze a 1/2/4 proposal family."""

    if method not in (PRIMARY_METHOD, COMPARATOR_METHOD):
        _fail("unrecognized method")
    if arithmetic_id not in _ARITHMETIC_IDS:
        _fail("unrecognized arithmetic")
    owned = _positive_int(owned_row_count, name="owned_row_count")
    frozen = _freeze_paths(paths)
    _validate_family_contents(
        frozen, method=method, arithmetic_id=arithmetic_id, owned=owned
    )
    digest = _family_digest(
        frozen, method=method, arithmetic_id=arithmetic_id, owned_row_count=owned
    )
    return IndependentRecordedFamily(
        paths=frozen,
        method=method,
        arithmetic_id=arithmetic_id,
        owned_row_count=owned,
        family_sha256=digest,
    )


def _require_unchanged_family(family: IndependentRecordedFamily) -> None:
    # NumPy's read-only flag can be reversed by a caller. The frozen dataclass
    # therefore does not authenticate later reads of the recorded arrays.
    if (
        _family_digest(
            family.paths,
            method=family.method,
            arithmetic_id=family.arithmetic_id,
            owned_row_count=family.owned_row_count,
        )
        != family.family_sha256
    ):
        _fail("recorded family changed after validation")


def _channel_parts(channel: str) -> tuple[str, str]:
    if channel not in _CHANNELS:
        _fail("channel is not a TDG6 complete-state channel")
    block, field = channel.split(":")
    return block, field


def _sample(array: np.ndarray, *, field: str, row: int) -> Fraction:
    return _exact(array[_OWNED][row, _FIELDS.index(field)])


def _embedded_defect(
    *,
    method: str,
    h: Fraction,
    ks: tuple[Fraction, ...],
    endpoint: Fraction,
) -> Fraction:
    if method == PRIMARY_METHOD:
        return abs(h * (ks[3] - endpoint) / 6)
    return abs(h * (-ks[0] - ks[1] + 2 * ks[2]) / 3)


def _level_segments(row: RefinementRow, level: int) -> PathSegments:
    if level == 0:
        return (row[0],)
    return row[level]


def _absolute_delta_sum(
    raw_segments: PathSegments, corrected_segments: PathSegments
) -> Fraction:
    total = Fraction(0)
    gap = Fraction(0)
    if not raw_segments or len(raw_segments) != len(corrected_segments):
        _fail("raw and corrected segment streams differ")
    if raw_segments[0][0] != corrected_segments[0][0]:
        _fail("corrected reconstruction must start at the recorded initial value")
    for raw, corrected in zip(raw_segments, corrected_segments, strict=True):
        if raw[4] != corrected[4]:
            _fail("raw and corrected segments must share the recorded width")
        if raw[1] != corrected[1] or raw[3] != corrected[3]:
            _fail("corrected segments must retain the recorded endpoint RHS slopes")
        delta = (raw[2] - raw[0]) - (corrected[2] - corrected[0])
        gap += delta
        if raw[2] - corrected[2] != gap:
            _fail("recorded-minus-corrected identity y-z=sum(delta) failed")
        total += abs(delta)
    return total


def _reconstruct_path(
    path: tuple[StepProposal, ...],
    *,
    block: str,
    derivative: str,
    field: str,
    method: str,
    row_count: int,
) -> tuple[tuple[PathSegments, ...], tuple[PathSegments, ...], Fraction, Fraction]:
    weights = _RK4_WEIGHTS if method == PRIMARY_METHOD else _SSP_WEIGHTS
    raw_pending: list[list[Segment]] = [[] for _ in range(row_count)]
    corr_pending: list[list[Segment]] = [[] for _ in range(row_count)]
    abs_delta = [Fraction(0)] * row_count
    abs_embedded = [Fraction(0)] * row_count
    z = [
        _sample(getattr(path[0].initial_state, block), field=field, row=row)
        for row in range(row_count)
    ]
    gap = [Fraction(0)] * row_count
    for proposal in path:
        h = _exact(proposal.final_time - proposal.initial_time)
        for row in range(row_count):
            y_old = _sample(
                getattr(proposal.initial_state, block), field=field, row=row
            )
            y_new = _sample(
                getattr(proposal.candidate_state, block), field=field, row=row
            )
            ks = tuple(
                _sample(getattr(stage.rhs, derivative), field=field, row=row)
                for stage in proposal.stages[: len(weights)]
            )
            left_rhs = _sample(
                getattr(proposal.stages[0].rhs, derivative), field=field, row=row
            )
            right_rhs = _sample(
                getattr(proposal.stages[-1].rhs, derivative), field=field, row=row
            )
            increment = h * sum(
                weight * stage for weight, stage in zip(weights, ks, strict=True)
            )
            if y_old - z[row] != gap[row]:
                _fail("recorded-minus-corrected identity y-z=sum(delta) failed")
            delta = y_new - y_old - increment
            z_next = z[row] + increment
            gap[row] += delta
            if y_new - z_next != gap[row]:
                _fail("recorded-minus-corrected identity y-z=sum(delta) failed")
            raw_pending[row].append((y_old, left_rhs, y_new, right_rhs, h))
            corr_pending[row].append((z[row], left_rhs, z_next, right_rhs, h))
            abs_delta[row] += abs(delta)
            abs_embedded[row] += _embedded_defect(
                method=method, h=h, ks=ks, endpoint=right_rhs
            )
            z[row] = z_next
    return (
        tuple(tuple(item) for item in raw_pending),
        tuple(tuple(item) for item in corr_pending),
        max(abs_delta),
        max(abs_embedded),
    )


def reconstruct_channel_independently(
    family: IndependentRecordedFamily, channel: str
) -> IndependentChannelReconstruction:
    """Return exact raw/corrected Hermite rows and per-level defect bounds."""

    if not isinstance(family, IndependentRecordedFamily):
        raise TypeError("family must be IndependentRecordedFamily")
    _require_unchanged_family(family)
    block, field = _channel_parts(channel)
    derivative = f"d{block}"
    row_count = family.owned_row_count
    raw_levels: list[tuple[PathSegments, ...]] = []
    corr_levels: list[tuple[PathSegments, ...]] = []
    accumulation: list[Fraction] = []
    embedded: list[Fraction] = []
    for path in family.paths:
        raw_path, corr_path, acc, emb = _reconstruct_path(
            path,
            block=block,
            derivative=derivative,
            field=field,
            method=family.method,
            row_count=row_count,
        )
        raw_levels.append(raw_path)
        corr_levels.append(corr_path)
        accumulation.append(acc)
        embedded.append(emb)
    raw_rows = tuple(
        (raw_levels[0][row][0], raw_levels[1][row], raw_levels[2][row])
        for row in range(row_count)
    )
    corrected_rows = tuple(
        (corr_levels[0][row][0], corr_levels[1][row], corr_levels[2][row])
        for row in range(row_count)
    )
    for level in range(3):
        recomputed = max(
            _absolute_delta_sum(
                _level_segments(raw, level), _level_segments(corrected, level)
            )
            for raw, corrected in zip(raw_rows, corrected_rows, strict=True)
        )
        if recomputed != accumulation[level]:
            _fail("fine debit inputs were not retained in the accumulation bound")
        if len(_level_segments(raw_rows[0], level)) != _PATH_LENGTHS[level]:
            _fail("every recorded step must contribute to the level debit")
    result = IndependentChannelReconstruction(
        channel=channel,
        row_count=row_count,
        raw_rows=raw_rows,
        corrected_rows=corrected_rows,
        accumulation_bounds=(accumulation[0], accumulation[1], accumulation[2]),
        embedded_defect_bounds=(embedded[0], embedded[1], embedded[2]),
    )
    _require_unchanged_family(family)
    return result


__all__ = [
    "ORIGINAL_ARITHMETIC_ID",
    "IndependentChannelReconstruction",
    "IndependentRecordedFamily",
    "IndependentRecordedFamilyError",
    "reconstruct_channel_independently",
    "validate_recorded_family_independently",
]
