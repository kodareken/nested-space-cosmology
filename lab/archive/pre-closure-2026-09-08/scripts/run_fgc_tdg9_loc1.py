#!/usr/bin/env python3
"""Run the bounded read-only TDG9 LOC1 localization diagnostic.

The runner reconstructs only the ten sealed PREF1 failure occurrences at RCV3
retry depths 3, 4, and 5.  It imports neither the AR1 runner nor the persisted
retry wrapper.  It never opens the campaign store for writing and publishes
only canonical ``manifest.json`` and ``terminal.json`` in a new namespace.
"""

from __future__ import annotations

import argparse
from copy import deepcopy
import ctypes
from dataclasses import asdict, dataclass, is_dataclass
from fractions import Fraction
from hashlib import sha256
import errno
import json
import math
import os
from pathlib import Path
import secrets
import stat
import sys
from typing import Any, Iterable, Iterator, Mapping, NoReturn, Sequence

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.fgc.evolution import proto15_runtime as p15  # noqa: E402
from recursive_horizons.fgc.evolution import tdg5_stage_complete_refinement_runtime as tdg5  # noqa: E402
from recursive_horizons.fgc.evolution import tdg6_temporal_admission_runtime as tdg6  # noqa: E402
from recursive_horizons.fgc.evolution import tdg7_binary64_subdivision_lattice as lattice  # noqa: E402
from recursive_horizons.fgc.evolution import tdg9_ar1_authority as ar1  # noqa: E402
from recursive_horizons.fgc.evolution import tdg9_exact_temporal_arithmetic as pref1_exact  # noqa: E402
from recursive_horizons.fgc.evolution import tdg9_exact_temporal_arithmetic_independent as pref1_exact_independent  # noqa: E402
from recursive_horizons.fgc.evolution import tdg9_loc1_authority as authority  # noqa: E402
from recursive_horizons.fgc.evolution import tdg9_local_extrema as primary  # noqa: E402
from recursive_horizons.fgc.evolution import tdg9_local_extrema_independent as independent  # noqa: E402
from recursive_horizons.fgc.evolution.hlt16_campaign_store import HLT16CampaignStore  # noqa: E402
from recursive_horizons.fgc.evolution.hlt16_member_codec import HLT16MemberSnapshot, restore_member  # noqa: E402
from recursive_horizons.fgc.evolution.proto19_gr0_static_factory import build_static_gr0_shells  # noqa: E402


RUNNER_ID = "FGC-1-TDG9-LOC1-RUN1"
RAW_SCHEMA = "FGC-1-TDG9-LOC1-raw-v1"
_BLOCKS = ("u", "p", "q")
_FIELDS = ("alpha", "v", "lambda", "R", "phi", "chi")
_MAX_DETAIL = 640


class LOC1RunnerError(RuntimeError):
    def __init__(self, owner: str, code: str, detail: object) -> None:
        self.owner, self.code = str(owner), str(code)
        self.detail = " ".join(str(detail).replace(str(ROOT), "<repo>").split())[:_MAX_DETAIL]
        super().__init__(f"{self.owner}/{self.code}: {self.detail}")


def _fail(owner: str, code: str, detail: object) -> NoReturn:
    raise LOC1RunnerError(owner, code, detail)


def _canonical(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def _pretty(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + "\n").encode()


def _pairs(items: list[tuple[str, Any]]) -> dict[str, Any]:
    answer: dict[str, Any] = {}
    for key, value in items:
        if key in answer:
            raise ValueError("duplicate JSON key")
        answer[key] = value
    return answer


def _read(root: Path, relative: str, maximum: int = 16 * 1024 * 1024) -> bytes:
    candidate = Path(relative)
    if candidate.is_absolute() or any(part in {"", ".", ".."} for part in candidate.parts):
        _fail("provenance", "unsafe_path", relative)
    path = root / candidate
    before = path.lstat()
    if stat.S_ISLNK(before.st_mode) or not stat.S_ISREG(before.st_mode) or before.st_size > maximum:
        _fail("provenance", "unsafe_leaf", relative)
    descriptor = os.open(path, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
    try:
        active = os.fstat(descriptor)
        identity = (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns)
        if identity != (active.st_dev, active.st_ino, active.st_size, active.st_mtime_ns):
            _fail("provenance", "raced_leaf", relative)
        raw = b""
        while len(raw) <= maximum:
            block = os.read(descriptor, min(1 << 20, maximum + 1 - len(raw)))
            if not block:
                break
            raw += block
        after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    if identity != (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns) or len(raw) != before.st_size:
        _fail("provenance", "changed_leaf", relative)
    return raw


def _authority(root: Path, commit: str) -> authority.LOC1Authority:
    try:
        return authority.authorize(root, _read(root, authority.CONFIG_PATH), _read(root, authority.RESULT_PATH), commit)
    except Exception as exc:
        raise LOC1RunnerError("authority", "rejected", exc) from exc


def _pref1_failure_map(root: Path) -> dict[tuple[int, str], Mapping[str, Any]]:
    raw = _read(root, authority.PREF1_RESULT_PATH)
    if sha256(raw).hexdigest() != authority.PREF1_RESULT_SHA256:
        _fail("provenance", "PREF1_result_hash", "tracked bytes differ")
    try:
        value = json.loads(raw, object_pairs_hook=_pairs)
        retries = value["artifact_payload"]["independent_replay"]["retries"]
    except (KeyError, TypeError, json.JSONDecodeError, ValueError) as exc:
        raise LOC1RunnerError("provenance", "PREF1_shape", exc) from exc
    answer: dict[tuple[int, str], Mapping[str, Any]] = {}
    for retry in retries:
        for channel in retry["channels"]:
            if channel["classification"] == "sufficient_contraction_failure":
                answer[(int(retry["retry"]), str(channel["channel"]))] = channel
    if tuple(answer) != authority.FAILED_OCCURRENCES:
        _fail("selection", "PREF1_failures", tuple(answer))
    return answer


def _journal(root: Path, replay: Mapping[str, object]) -> Mapping[str, object]:
    relative = f"{ar1.PREF2_STORE_PATH}/journal/{int(replay['journal_sequence']):020d}-{replay['journal_sha256']}.journal"
    raw = _read(root, relative)
    if sha256(raw).hexdigest() != replay["journal_raw_sha256"]:
        _fail("replay", "journal_hash", replay["retry"])
    value = json.loads(raw, object_pairs_hook=_pairs)
    evidence = value.get("payload", {}).get("evidence")
    if not isinstance(evidence, Mapping):
        _fail("replay", "journal_shape", replay["retry"])
    return evidence


def _plan(cursor: p15.Proto15Cursor, replay: Mapping[str, object]) -> Any:
    pending = cursor.payload.get("retry_successor_payload_or_none")
    if not isinstance(pending, Mapping):
        _fail("replay", "plan_absent", replay["retry"])
    persisted = p15._plan_mapping(pending["successor_plan"])
    plan = lattice.plan_forward_proto14_subdivision(float.fromhex(str(persisted["current_hex"])), float.fromhex(str(persisted["event_target_hex"])), float.fromhex(str(persisted["requested_cap_hex"])), minimum_width=float.fromhex(str(persisted["minimum_width_hex_or_none"])))
    if p15._plan_mapping(plan) != persisted or plan.macro_width.hex() != replay["attempted_width_hex"]:
        _fail("replay", "plan_drift", replay["retry"])
    return plan


def _prepare(root: Path, store: HLT16CampaignStore, shells: Mapping[str, object], replay: Mapping[str, object]) -> tuple[Any, Mapping[str, object]]:
    generation = int(replay["predecessor_generation"])
    checkpoint = store.authenticated_checkpoint_at_generation(generation)
    if checkpoint.sha256 != replay["checkpoint_sha256"]:
        _fail("replay", "checkpoint_hash", generation)
    relative = f"{ar1.PREF2_STORE_PATH}/checkpoints/{generation:020d}-{checkpoint.sha256}.json"
    if sha256(_read(root, relative)).hexdigest() != replay["checkpoint_raw_sha256"]:
        _fail("replay", "checkpoint_raw_hash", generation)
    member_state = checkpoint.members[ar1.MEMBER_KEY]
    member = deepcopy(shells[ar1.MEMBER_KEY])
    persisted = store.load_state(member_state.descriptor_sha256)
    restore_member(member, HLT16MemberSnapshot(persisted.arrays, persisted.descriptor["metadata"]))
    identity = persisted.descriptor["metadata"]["runtime_identity"]
    if identity["generation"] == 0:
        member.source_retry_count += member_state.source_total
        member.CFL_retry_count += member_state.cfl_total
    member.temporal_ledger = p15._ledger_from_mapping(member_state.ledger)
    cursor = p15.Proto15Cursor(dict(member_state.cursor))
    cursor.validate()
    ledger = member.temporal_ledger
    if ledger is None or cursor.mode != "RETRY_PENDING" or member_state.pending_owner != "temporal" or ledger.current_macro_step_temporal_retry_count != replay["prior_retry_count"]:
        _fail("replay", "predecessor_selection", replay["retry"])
    plan = _plan(cursor, replay)
    prepared = tdg6.prepare_tdg6_gr0_compositor(method=member.integrator_id, time=plan.current, step_size=plan.macro_width, state=member.state, rhs=member.operator, projector=member.projector, transaction=member.transaction, tracers=member.tracers, coordinates=member.initial.grid.coordinates, temporal_ledger=ledger, previous_step_index=member.step_index, previous_transaction_serial=member.transaction_serial)
    historical = _journal(root, replay)
    seen: list[Mapping[str, object]] = []
    try:
        tdg6.require_tdg6_temporal_admission(prepared, transaction=member.transaction, tracers=member.tracers, temporal_ledger=ledger, current_time=member.time, current_state=member.state, current_step_index=member.step_index, current_transaction_serial=member.transaction_serial, durable_rejection_sink=lambda item: seen.append(dict(item)))
    except tdg6.TDG6TemporalRetryRequired as rejected:
        if rejected.evidence.retry_count_for_current_macro_step != replay["retry"]:
            _fail("replay", "retry_count", replay["retry"])
    except tdg6.TDG6TemporalRetryExhausted as exc:
        raise LOC1RunnerError("replay", "unexpected_exhaustion", exc) from exc
    else:
        _fail("replay", "historical_rejection_admitted", replay["retry"])
    if len(seen) != 1 or _canonical(p15._json_safe(seen[0])) != _canonical(historical):
        _fail("replay", "historical_identity", replay["retry"])
    return prepared, historical


@dataclass(frozen=True, slots=True)
class _Segment:
    left: np.ndarray
    left_rhs: np.ndarray
    right: np.ndarray
    right_rhs: np.ndarray
    width: float


def _owned(value: object) -> np.ndarray:
    if type(value) is not np.ndarray or value.dtype.str != "<f8" or value.shape != (3, authority.OWNED_ROW_COUNT, 6) or not value.flags.c_contiguous or not bool(np.isfinite(value).all()):
        _fail("replay", "owned_array", type(value).__name__)
    return value


def _segment(proposal: Any) -> _Segment:
    start_rhs, end_rhs = tdg5._proposal_endpoint_records(proposal)
    width = proposal.final_time - proposal.initial_time
    if type(width) is not float or not math.isfinite(width) or width <= 0:
        _fail("replay", "segment_width", width)
    return _Segment(_owned(tdg5._owned_state(proposal.initial_state)), _owned(tdg5._owned_rhs(start_rhs)), _owned(tdg5._owned_state(proposal.candidate_state)), _owned(tdg5._owned_rhs(end_rhs)), width)


def _surface(prepared: Any) -> tuple[tuple[_Segment, ...], ...]:
    paths = tuple(tuple(_segment(attempt.proposal) for attempt in path.attempts) for path in (prepared.outer, prepared.medium, prepared.fine))
    if tuple(map(len, paths)) != (1, 2, 4):
        _fail("replay", "path_shape", tuple(map(len, paths)))
    return paths


def _binary_segment(item: _Segment, block: int, row: int, field: int) -> tuple[float, float, float, float, float]:
    values = (float(item.left[block, row, field]), float(item.left_rhs[block, row, field]), float(item.right[block, row, field]), float(item.right_rhs[block, row, field]), item.width)
    if any(type(value) is not float or not math.isfinite(value) for value in values):
        _fail("replay", "binary_segment", row)
    return values


def _rows(surface: tuple[tuple[_Segment, ...], ...], channel: str) -> Iterator[object]:
    block_name, field_name = channel.split(":")
    block, field = _BLOCKS.index(block_name), _FIELDS.index(field_name)
    outer, medium, fine = surface
    for row in range(authority.OWNED_ROW_COUNT):
        yield (_binary_segment(outer[0], block, row, field), tuple(_binary_segment(item, block, row, field) for item in medium), tuple(_binary_segment(item, block, row, field) for item in fine))


def _row_hash(rows: Iterable[object]) -> str:
    digest = sha256(b"TDG9-AR1-BINARY64-HERMITE-ROWS-v1\n")
    for ordinal, row in enumerate(rows):
        outer, medium, fine = row  # type: ignore[misc]
        values = [value for segment in (outer, *medium, *fine) for value in segment]
        digest.update((f"{ordinal}|" + "|".join(value.hex() for value in values) + "\n").encode())
    return digest.hexdigest()


def _exact_hermite(segment: Sequence[float]) -> tuple[Fraction, Fraction, Fraction, Fraction]:
    left, left_rhs, right, right_rhs, width = (Fraction(*value.as_integer_ratio()) for value in segment)
    slope_left, slope_right = width * left_rhs, width * right_rhs
    return left, slope_left, -3 * left - 2 * slope_left + 3 * right - slope_right, 2 * left + slope_left - 2 * right + slope_right


def _restrict(coefficients: tuple[Fraction, ...], half: int) -> tuple[Fraction, Fraction, Fraction, Fraction]:
    a0, a1, a2, a3 = coefficients
    if half == 0:
        return a0, a1 / 2, a2 / 4, a3 / 8
    return a0 + a1 / 2 + a2 / 4 + a3 / 8, a1 / 2 + a2 / 2 + 3 * a3 / 8, a2 / 4 + 3 * a3 / 8, a3 / 8


def _difference(parent: Sequence[float], child: Sequence[float], half: int) -> tuple[Fraction, Fraction, Fraction, Fraction]:
    return _component_coefficients(parent, child, half)["complete_C"]


def _component_coefficients(
    parent: Sequence[float],
    child: Sequence[float],
    half: int,
) -> dict[str, tuple[Fraction, Fraction, Fraction, Fraction]]:
    """Split one restricted-parent minus child cubic into value and slope parts."""

    restricted = _restrict(_exact_hermite(parent), half)
    child_coefficients = _exact_hermite(child)
    complete = tuple(
        left - right
        for left, right in zip(restricted, child_coefficients, strict=True)
    )
    parent_left = restricted[0]
    parent_right = sum(restricted, Fraction(0))
    parent_left_slope = restricted[1]
    parent_right_slope = restricted[1] + 2 * restricted[2] + 3 * restricted[3]
    child_left = Fraction(*child[0].as_integer_ratio())
    child_right = Fraction(*child[2].as_integer_ratio())
    child_width = Fraction(*child[4].as_integer_ratio())
    state_left = parent_left - child_left
    state_right = parent_right - child_right
    rhs_left = parent_left_slope - child_width * Fraction(*child[1].as_integer_ratio())
    rhs_right = parent_right_slope - child_width * Fraction(*child[3].as_integer_ratio())
    value = (
        state_left,
        Fraction(0),
        -3 * state_left + 3 * state_right,
        2 * state_left - 2 * state_right,
    )
    slope = (
        Fraction(0),
        rhs_left,
        -2 * rhs_left - rhs_right,
        rhs_left + rhs_right,
    )
    if tuple(left + right for left, right in zip(value, slope, strict=True)) != complete:
        _fail("localization", "component_identity", half)
    return {"value_V": value, "slope_S": slope, "complete_C": complete}  # type: ignore[return-value]


def _region(row: int) -> str:
    for start, end, name in authority.ROW_REGIONS:
        if start <= row <= end:
            return name
    _fail("localization", "row_region", row)


def _cubics(
    rows: Sequence[object],
    level: str,
    component: str = "complete_C",
) -> tuple[list[primary.LocalCubic], list[independent.IndependentLocalCubic]]:
    if component not in {"value_V", "slope_S", "complete_C"}:
        _fail("localization", "component", component)
    first: list[primary.LocalCubic] = []
    second: list[independent.IndependentLocalCubic] = []
    ordinal = 0
    for row, raw in enumerate(rows):
        outer, medium, fine = raw  # type: ignore[misc]
        pairs = ((outer, medium[index], index, 0, index) for index in range(2)) if level == "D01" else ((medium[index // 2], fine[index], index, index // 2, index) for index in range(4))
        for parent, child, subinterval, parent_index, child_index in pairs:
            half = subinterval if level == "D01" else subinterval % 2
            coefficients = _component_coefficients(parent, child, half)[component]
            metadata = {"component": component, "level": level, "subinterval": subinterval, "parent_index": parent_index, "child_index": child_index, "owned_row": row, "global_index": row + 1, "radius_hex": float((row + 1) / 16).hex(), "row_region": _region(row), "ordinal": ordinal}
            first.append(primary.LocalCubic(coefficients, metadata))
            second.append(independent.IndependentLocalCubic(coefficients, metadata))
            ordinal += 1
    return first, second


def _fraction(value: Fraction) -> dict[str, str]:
    return {"numerator": str(value.numerator), "denominator": str(value.denominator)}


def _json_exact(value: object) -> object:
    if isinstance(value, Fraction):
        return _fraction(value)
    if is_dataclass(value) and not isinstance(value, type):
        return _json_exact(asdict(value))
    if isinstance(value, Mapping):
        return {str(key): _json_exact(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_exact(item) for item in value]
    if value is None or isinstance(value, (str, int, bool)):
        return value
    _fail("serialization", "nonexact_value", type(value).__name__)


def _candidate_key(item: object) -> tuple[object, ...]:
    return (item.polynomial_ordinal, item.location, item.location_ordinal, tuple(item.metadata))  # type: ignore[attr-defined]


def _overlap(left: object, right: object) -> bool:
    return max(left.global_absolute_lower, right.global_absolute_lower) <= min(left.global_absolute_upper, right.global_absolute_upper)  # type: ignore[attr-defined]


def _candidate_intervals_overlap(left: object, right: object) -> bool:
    return (
        max(left.parameter_lower, right.parameter_lower)  # type: ignore[attr-defined]
        <= min(left.parameter_upper, right.parameter_upper)  # type: ignore[attr-defined]
        and max(left.absolute_lower, right.absolute_lower)  # type: ignore[attr-defined]
        <= min(left.absolute_upper, right.absolute_upper)  # type: ignore[attr-defined]
    )


def _require_route_agreement(first: object, second: object, detail: object) -> None:
    first_items = {_candidate_key(item): item for item in first.candidates}  # type: ignore[attr-defined]
    second_items = {_candidate_key(item): item for item in second.candidates}  # type: ignore[attr-defined]
    if (
        first.polynomial_count != second.polynomial_count  # type: ignore[attr-defined]
        or first.candidate_count != second.candidate_count  # type: ignore[attr-defined]
        or first.classification != second.classification  # type: ignore[attr-defined]
        or len(first.candidates) != len(second.candidates)  # type: ignore[attr-defined]
        or len(first_items) != len(first.candidates)  # type: ignore[attr-defined]
        or len(second_items) != len(second.candidates)  # type: ignore[attr-defined]
        or tuple(sorted(first_items)) != tuple(sorted(second_items))
        or not _overlap(first, second)
        or any(
            not _candidate_intervals_overlap(first_items[key], second_items[key])
            for key in first_items
        )
    ):
        _fail("localization", "independent_disagreement", detail)


def _ulp_record(value: float) -> dict[str, object]:
    exact = Fraction(*value.as_integer_ratio())
    lower = Fraction(*math.nextafter(value, -math.inf).as_integer_ratio())
    upper = Fraction(*math.nextafter(value, math.inf).as_integer_ratio())
    return {"hex": value.hex(), "exact": _fraction(exact), "downward": _fraction(exact - lower), "upward": _fraction(upper - exact)}


def _hermite_inputs(rows: Sequence[object], level: str, candidate: object) -> dict[str, object]:
    metadata = dict(candidate.metadata)  # type: ignore[attr-defined]
    row, subinterval = int(metadata["owned_row"]), int(metadata["subinterval"])
    outer, medium, fine = rows[row]  # type: ignore[misc]
    if level == "D01":
        parent, child, half = outer, medium[subinterval], subinterval
    else:
        parent, child, half = medium[subinterval // 2], fine[subinterval], subinterval % 2
    restricted = _restrict(_exact_hermite(parent), half)
    components = _component_coefficients(parent, child, half)
    parent_left = restricted[0]
    parent_right = sum(restricted, Fraction(0))
    parent_left_slope = restricted[1]
    parent_right_slope = restricted[1] + 2 * restricted[2] + 3 * restricted[3]
    child_left = Fraction(*child[0].as_integer_ratio())
    child_right = Fraction(*child[2].as_integer_ratio())
    child_width = Fraction(*child[4].as_integer_ratio())
    child_left_slope = child_width * Fraction(*child[1].as_integer_ratio())
    child_right_slope = child_width * Fraction(*child[3].as_integer_ratio())
    state_left = parent_left - child_left
    state_right = parent_right - child_right
    rhs_left = parent_left_slope - child_left_slope
    rhs_right = parent_right_slope - child_right_slope
    parent_update = parent_right - parent_left
    child_update = Fraction(*child[2].as_integer_ratio()) - Fraction(*child[0].as_integer_ratio())
    return {"restricted_parent_endpoint_data": {"left_value": _fraction(parent_left), "right_value": _fraction(parent_right), "left_normalized_slope": _fraction(parent_left_slope), "right_normalized_slope": _fraction(parent_right_slope)}, "child_endpoint_data": {"left_value": _fraction(child_left), "right_value": _fraction(child_right), "left_normalized_slope": _fraction(child_left_slope), "right_normalized_slope": _fraction(child_right_slope)}, "endpoint_state_difference": {"left": _fraction(state_left), "right": _fraction(state_right)}, "width_scaled_endpoint_RHS_difference": {"left": _fraction(rhs_left), "right": _fraction(rhs_right)}, "component_polynomials": {name: [_fraction(item) for item in coefficients] for name, coefficients in components.items()}, "complete_difference_polynomial": [_fraction(item) for item in components["complete_C"]], "complete_equals_value_plus_slope": True, "updates": {"restricted_parent": _fraction(parent_update), "child": _fraction(child_update)}, "binary64_contributors": {"parent_left": _ulp_record(parent[0]), "parent_left_RHS": _ulp_record(parent[1]), "parent_right": _ulp_record(parent[2]), "parent_right_RHS": _ulp_record(parent[3]), "parent_width": _ulp_record(parent[4]), "child_left": _ulp_record(child[0]), "child_left_RHS": _ulp_record(child[1]), "child_right": _ulp_record(child[2]), "child_right_RHS": _ulp_record(child[3]), "child_width": _ulp_record(child[4])}, "ULP_used_as_tolerance": False}


def _contraction_assessment(outer: object, finest: object) -> dict[str, object]:
    left_lower = outer.global_absolute_lower**2  # type: ignore[attr-defined]
    left_upper = outer.global_absolute_upper**2  # type: ignore[attr-defined]
    right_lower = 8 * finest.global_absolute_lower**2  # type: ignore[attr-defined]
    right_upper = 8 * finest.global_absolute_upper**2  # type: ignore[attr-defined]
    exact_zero = outer.global_absolute_upper == 0 and finest.global_absolute_upper == 0  # type: ignore[attr-defined]
    classification = _classify_contraction_bounds(
        left_lower,
        left_upper,
        right_lower,
        right_upper,
        exact_zero=exact_zero,
    )
    sufficient_pass = classification == "sufficient_contraction_pass"
    sufficient_failure = classification == "sufficient_contraction_failure"
    return {
        "classification": classification,
        "sufficient_contraction_pass_certified": sufficient_pass,
        "sufficient_contraction_failure_certified": sufficient_failure,
        "threshold_inconclusive": classification == "threshold_inconclusive",
        "left_D01_squared": {"lower": _fraction(left_lower), "upper": _fraction(left_upper)},
        "right_8_D12_squared": {"lower": _fraction(right_lower), "upper": _fraction(right_upper)},
    }


def _classify_contraction_bounds(
    left_lower: Fraction,
    left_upper: Fraction,
    right_lower: Fraction,
    right_upper: Fraction,
    *,
    exact_zero: bool,
) -> str:
    if exact_zero:
        return "exact_zero"
    if right_upper <= left_lower:
        return "sufficient_contraction_pass"
    if right_lower > left_upper:
        return "sufficient_contraction_failure"
    return "threshold_inconclusive"


def _component_ownership(components: Mapping[str, Mapping[str, object]]) -> dict[str, object]:
    classes = {
        name: str(value["contraction"]["classification"])  # type: ignore[index]
        for name, value in components.items()
    }
    value_failure = classes["value_V"] == "sufficient_contraction_failure"
    slope_failure = classes["slope_S"] == "sufficient_contraction_failure"
    complete_failure = classes["complete_C"] == "sufficient_contraction_failure"
    value_pass = classes["value_V"] == "sufficient_contraction_pass"
    slope_pass = classes["slope_S"] == "sufficient_contraction_pass"
    value_clear = value_pass or classes["value_V"] == "exact_zero"
    slope_clear = slope_pass or classes["slope_S"] == "exact_zero"
    if complete_failure and value_failure and slope_clear:
        classification = "endpoint_state_owned_failure"
    elif complete_failure and slope_failure and value_clear:
        classification = "width_scaled_RHS_owned_failure"
    elif complete_failure and value_failure and slope_failure:
        classification = "both_components_independently_fail"
    elif complete_failure and value_pass and slope_pass:
        classification = "mixed_only_failure"
    elif complete_failure:
        classification = "component_ownership_inconclusive"
    else:
        classification = "complete_not_failure"
    cancellation_by_level: dict[str, bool] = {}
    reinforcement_by_level: dict[str, bool] = {}
    for level in ("D01", "D12"):
        value = components["value_V"][level]  # type: ignore[index]
        slope = components["slope_S"][level]  # type: ignore[index]
        complete = components["complete_C"][level]  # type: ignore[index]
        cancellation_by_level[level] = complete.global_absolute_upper < max(  # type: ignore[attr-defined]
            value.global_absolute_lower, slope.global_absolute_lower  # type: ignore[attr-defined]
        )
        reinforcement_by_level[level] = complete.global_absolute_lower > max(  # type: ignore[attr-defined]
            value.global_absolute_upper, slope.global_absolute_upper  # type: ignore[attr-defined]
        )
    return {
        "classification": classification,
        "value_only_failure": complete_failure and value_failure and slope_clear,
        "slope_only_failure": complete_failure and slope_failure and value_clear,
        "both_components_independently_fail": value_failure and slope_failure,
        "mixed_only_failure": complete_failure and value_pass and slope_pass,
        "cancellation_certified_by_level": cancellation_by_level,
        "reinforcement_certified_by_level": reinforcement_by_level,
    }


def _localize_occurrence(rows: Sequence[object], retry: int, channel: str, pref1: Mapping[str, Any]) -> dict[str, object]:
    row_hash = _row_hash(rows)
    if row_hash != pref1["row_stream_sha256"]:
        _fail("replay", "PREF1_row_hash", (retry, channel))
    # Reproduce the sealed primary evidence identity before localization.
    sealed = pref1_exact.assess_exact_temporal_refinement(rows, max_depth=ar1.MAX_DEPTH, max_nodes=ar1.MAX_NODES_PER_CHANNEL_PER_EVALUATOR)
    sealed_hash = sha256(_canonical(_json_exact(sealed))).hexdigest()
    sealed_independent = pref1_exact_independent.assess_exact_temporal_refinement_independently(rows, max_depth=ar1.MAX_DEPTH, max_nodes=ar1.MAX_NODES_PER_CHANNEL_PER_EVALUATOR)
    sealed_independent_hash = sha256(_canonical(_json_exact(sealed_independent))).hexdigest()
    if (
        sealed_hash != pref1["primary_evidence_sha256"]
        or sealed_independent_hash != pref1["independent_evidence_sha256"]
        or sealed.classification != "sufficient_contraction_failure"
        or sealed_independent.classification != sealed.classification
        or sealed.combined_coefficient_stream_sha256 != sealed_independent.combined_coefficient_stream_sha256
    ):
        _fail("replay", "PREF1_evidence_hash", (retry, channel))
    component_evidence: dict[str, dict[str, object]] = {}
    levels: dict[str, dict[str, object]] = {"D01": {}, "D12": {}}
    for component in ("value_V", "slope_S", "complete_C"):
        component_evidence[component] = {}
        for level, expected in (("D01", 2 * authority.OWNED_ROW_COUNT), ("D12", 4 * authority.OWNED_ROW_COUNT)):
            first_cubics, second_cubics = _cubics(rows, level, component)
            if len(first_cubics) != expected or len(second_cubics) != expected:
                _fail("localization", "polynomial_count", (retry, channel, component, level))
            first = primary.localize_absolute_maximum(
                first_cubics,
                maximum_candidates=4 * expected,
                refinement_depth=authority.PRIMARY_REFINEMENT_DEPTH,
            )
            second = independent.localize_absolute_maximum_independently(
                second_cubics,
                maximum_candidates=4 * expected,
                refinement_bits=authority.INDEPENDENT_REFINEMENT_BITS,
            )
            _require_route_agreement(first, second, (retry, channel, component, level))
            component_evidence[component][level] = first
            levels[level][component] = {
                "polynomial_count": expected,
                "primary": _json_exact(first),
                "independent": _json_exact(second),
                "co_maximizer_count": len(first.candidates),
                "maximizer_input_decomposition": [
                    _hermite_inputs(rows, level, item) for item in first.candidates
                ],
            }
        contraction = _contraction_assessment(
            component_evidence[component]["D01"],
            component_evidence[component]["D12"],
        )
        levels["D01"][component]["contraction"] = contraction
        levels["D12"][component]["contraction"] = contraction
        component_evidence[component]["contraction"] = contraction
    if component_evidence["complete_C"]["contraction"]["classification"] != sealed.classification:  # type: ignore[index]
        _fail("localization", "complete_classification_drift", (retry, channel))
    ownership = _component_ownership(component_evidence)
    return {
        "retry": retry,
        "channel": channel,
        "PREF1_row_stream_sha256": row_hash,
        "PREF1_primary_evidence_sha256": sealed_hash,
        "PREF1_independent_evidence_sha256": sealed_independent_hash,
        "levels": [
            {"level": level, "base_cubic_count": 2 * authority.OWNED_ROW_COUNT if level == "D01" else 4 * authority.OWNED_ROW_COUNT, "components": levels[level]}
            for level in ("D01", "D12")
        ],
        "component_ownership": ownership,
    }


def _snapshot_store(root: Path) -> tuple[int, str]:
    base = root / ar1.PREF2_STORE_PATH
    metadata = base.lstat()
    if stat.S_ISLNK(metadata.st_mode) or not stat.S_ISDIR(metadata.st_mode):
        _fail("provenance", "store_root_unsafe", ar1.PREF2_STORE_PATH)
    digest = sha256(b"TDG9-LOC1-STORE-SNAPSHOT-v1\n")
    count = 0
    stack = [base]
    while stack:
        directory = stack.pop()
        with os.scandir(directory) as entries:
            items = sorted(entries, key=lambda item: item.name)
        for item in items:
            metadata = item.stat(follow_symlinks=False)
            if item.is_symlink():
                _fail("provenance", "store_symlink", item.path)
            if stat.S_ISDIR(metadata.st_mode):
                stack.append(Path(item.path))
            elif stat.S_ISREG(metadata.st_mode):
                relative = Path(item.path).relative_to(base).as_posix()
                raw = _read(base, relative, 128 * 1024 * 1024)
                digest.update((relative + "\0" + sha256(raw).hexdigest() + "\n").encode())
                count += 1
            else:
                _fail("provenance", "store_foreign_type", item.path)
    return count, digest.hexdigest()


def _open_output_parent(root: Path) -> tuple[int, str]:
    relative = Path(authority.OUTPUT_NAMESPACE)
    if relative.is_absolute() or any(part in {"", ".", ".."} for part in relative.parts):
        _fail("output", "unsafe_namespace", authority.OUTPUT_NAMESPACE)
    try:
        descriptor = os.open(
            root,
            os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0),
        )
    except OSError as exc:
        raise LOC1RunnerError("output", "unsafe_root", exc.errno) from exc
    try:
        for part in relative.parent.parts:
            try:
                child = os.open(
                    part,
                    os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0),
                    dir_fd=descriptor,
                )
            except FileNotFoundError:
                try:
                    os.mkdir(part, mode=0o755, dir_fd=descriptor)
                except FileExistsError:
                    pass
                try:
                    child = os.open(
                        part,
                        os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_NOFOLLOW", 0),
                        dir_fd=descriptor,
                    )
                except OSError as exc:
                    raise LOC1RunnerError("output", "unsafe_parent_component", (part, exc.errno)) from exc
            except OSError as exc:
                raise LOC1RunnerError("output", "unsafe_parent_component", (part, exc.errno)) from exc
            os.close(descriptor)
            descriptor = child
        return descriptor, relative.name
    except Exception:
        os.close(descriptor)
        raise


def _publish(root: Path, manifest: Mapping[str, object], terminal: Mapping[str, object]) -> None:
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
        for name, value in (("manifest.json", manifest), ("terminal.json", terminal)):
            descriptor = os.open(
                name,
                os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0),
                0o644,
                dir_fd=stage_fd,
            )
            try:
                pending = memoryview(_pretty(value))
                while pending:
                    written = os.write(descriptor, pending)
                    if written <= 0:
                        _fail("output", "short_write", name)
                    pending = pending[written:]
                os.fsync(descriptor)
            finally:
                os.close(descriptor)
        os.fsync(stage_fd)
        if os.uname().sysname != "Darwin":
            _fail("output", "exclusive_adoption_unavailable", os.uname().sysname)
        libc = ctypes.CDLL(None, use_errno=True)
        renamex = getattr(libc, "renameatx_np", None)
        if renamex is None:
            _fail("output", "exclusive_adoption_unavailable", "renameatx_np")
        renamex.argtypes = (ctypes.c_int, ctypes.c_char_p, ctypes.c_int, ctypes.c_char_p, ctypes.c_uint)
        renamex.restype = ctypes.c_int
        if renamex(parent_fd, os.fsencode(stage_name), parent_fd, os.fsencode(target_name), 0x00000004) != 0:
            error = ctypes.get_errno()
            if error in {errno.EEXIST, errno.ENOTEMPTY}:
                _fail("output", "namespace_arrived", authority.OUTPUT_NAMESPACE)
            _fail("output", "exclusive_adoption_failed", error)
        os.close(stage_fd)
        stage_fd = None
        os.fsync(parent_fd)
    except Exception:
        if stage_fd is not None:
            for name in ("manifest.json", "terminal.json"):
                try:
                    os.unlink(name, dir_fd=stage_fd)
                except FileNotFoundError:
                    pass
            os.close(stage_fd)
        try:
            os.rmdir(stage_name, dir_fd=parent_fd)
        except FileNotFoundError:
            pass
        raise
    finally:
        os.close(parent_fd)


def run(root: Path, *, authority_commit: str) -> dict[str, object]:
    repository = root.resolve()
    receipt = _authority(repository, authority_commit)
    failures = _pref1_failure_map(repository)
    before = _snapshot_store(repository)
    store = HLT16CampaignStore(repository / ar1.PREF2_STORE_PATH)
    shells = build_static_gr0_shells(repository)
    coordinates = shells[authority.MEMBER_KEY].initial.grid.coordinates
    expected_coordinates = np.arange(authority.POINT_COUNT, dtype=np.float64) * float.fromhex(authority.GRID_SPACING_HEX)
    if (
        type(coordinates) is not np.ndarray
        or coordinates.dtype.str != "<f8"
        or coordinates.shape != (authority.POINT_COUNT,)
        or not np.array_equal(coordinates, expected_coordinates)
        or float(coordinates[-1]).hex() != authority.OUTER_RADIUS_HEX
    ):
        _fail("selection", "grid_identity", getattr(coordinates, "shape", None))
    prepared_by_retry: dict[int, Any] = {}
    historical_hashes: dict[int, str] = {}
    for replay in ar1.REPLAYS:
        prepared, historical = _prepare(repository, store, shells, replay)
        prepared_by_retry[int(replay["retry"])] = prepared
        historical_hashes[int(replay["retry"])] = sha256(_canonical(historical)).hexdigest()
    occurrences: list[dict[str, object]] = []
    base_cubic_count = 0
    candidate_counts = {"value_V": 0, "slope_S": 0, "complete_C": 0}
    for retry, channel in authority.FAILED_OCCURRENCES:
        rows = tuple(_rows(_surface(prepared_by_retry[retry]), channel))
        result = _localize_occurrence(rows, retry, channel, failures[(retry, channel)])
        occurrences.append(result)
        for level in result["levels"]:  # type: ignore[assignment]
            base_cubic_count += int(level["base_cubic_count"])  # type: ignore[index]
            for component in candidate_counts:
                candidate_counts[component] += int(level["components"][component]["primary"]["candidate_count"])  # type: ignore[index]
    total_candidate_count = sum(candidate_counts.values())
    if (
        base_cubic_count != receipt.total_polynomials
        or candidate_counts["complete_C"] > authority.COMPLETE_C_MAXIMUM_CANDIDATES
        or candidate_counts["value_V"] > authority.VALUE_V_MAXIMUM_CANDIDATES
        or candidate_counts["slope_S"] > authority.SLOPE_S_MAXIMUM_CANDIDATES
        or total_candidate_count > receipt.total_VSC_maximum_candidates
    ):
        _fail("localization", "global_budget", (base_cubic_count, candidate_counts))
    after = _snapshot_store(repository)
    if before != after:
        _fail("provenance", "store_mutated", (before, after))
    manifest = {"schema": RAW_SCHEMA, "artifact_id": authority.ARTIFACT_ID, "runner_id": RUNNER_ID, "authority_commit": authority_commit, "PREF1_commit": authority.PREF1_COMMIT, "failed_occurrences": [{"retry": retry, "channel": channel} for retry, channel in authority.FAILED_OCCURRENCES], "base_cubic_count": authority.TOTAL_POLYNOMIALS, "component_cubic_counts": {"complete_C": authority.COMPONENT_CUBIC_COUNT, "value_V": authority.COMPONENT_CUBIC_COUNT, "slope_S": authority.COMPONENT_CUBIC_COUNT, "total_VSC": authority.EXECUTED_COMPONENT_CUBIC_COUNT}, "candidate_ceilings": {"complete_C": authority.COMPLETE_C_MAXIMUM_CANDIDATES, "value_V": authority.VALUE_V_MAXIMUM_CANDIDATES, "slope_S": authority.SLOPE_S_MAXIMUM_CANDIDATES, "total_VSC": authority.TOTAL_VSC_MAXIMUM_CANDIDATES}, "output_leaves": ["manifest.json", "terminal.json"]}
    terminal = {"schema": RAW_SCHEMA, "artifact_id": authority.ARTIFACT_ID, "runner_id": RUNNER_ID, "classification": "bounded_exact_localization_completed", "authority_commit": authority_commit, "historical_TDG6_evidence_sha256": historical_hashes, "base_cubic_count": base_cubic_count, "executed_component_cubic_count": authority.EXECUTED_COMPONENT_CUBIC_COUNT, "component_candidate_counts": candidate_counts, "total_VSC_candidate_count": total_candidate_count, "occurrences": occurrences, "store_snapshot_before": {"leaf_count": before[0], "sha256": before[1]}, "store_snapshot_after": {"leaf_count": after[0], "sha256": after[1]}, "store_unchanged": True, "ULP_used_as_tolerance": False, "fourth_width_executed": False, "resource_escalation_used": False, "stage2_source_decomposition_executed": False, "SSPRK3_comparator_executed": False, "campaign_store_mutated": False, "PDE_state_committed": False, "continuation_authorized": False, "candidate_branch_opened": False, "physical_result_earned": False}
    _publish(repository, manifest, terminal)
    return terminal


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--authority-commit", required=True)
    arguments = parser.parse_args(argv)
    result = run(ROOT, authority_commit=arguments.authority_commit)
    print(json.dumps({"artifact_id": authority.ARTIFACT_ID, "classification": result["classification"], "base_cubic_count": result["base_cubic_count"], "campaign_store_mutated": False, "candidate_branch_opened": False}, sort_keys=True, separators=(",", ":")))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
