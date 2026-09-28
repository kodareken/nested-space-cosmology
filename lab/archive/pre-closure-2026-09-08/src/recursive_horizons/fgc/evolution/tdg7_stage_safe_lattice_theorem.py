"""Independent exact theorem controls for TDG7's stage-safe binary64 lattice.

This module deliberately does *not* import the prospective TDG7 design helper.
It reconstructs IEEE-754 binary64 values from their bits, derives the ``8Q``
selection with :class:`fractions.Fraction`, and checks the actual RK4/SSPRK3
time expressions by parsing the numerical-engine source syntax tree.

The theorem is restricted to the forward PROTO14 coordinate envelope.  It is
not a runtime adapter, a campaign runner, or a physics result.
"""

from __future__ import annotations

import ast
from dataclasses import dataclass
from fractions import Fraction
import math
import struct
from typing import Final, Iterable


Q = Fraction
TIME_LOWER: Final[Fraction] = Q(23, 16)
TIME_UPPER: Final[Fraction] = Q(32)


class IndependentTDG7TheoremError(ValueError):
    """One exact stage-safe lattice assumption is not satisfied."""

    def __init__(self, reason: str) -> None:
        self.reason = reason
        super().__init__(f"tdg7_independent_theorem_error:{reason}")


def _bits(value: float) -> int:
    return struct.unpack(">Q", struct.pack(">d", value))[0]


def _power_of_two(exponent: int) -> Fraction:
    return Q(1 << exponent) if exponent >= 0 else Q(1, 1 << (-exponent))


def _binary64(value: object, *, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (float, int)):
        raise IndependentTDG7TheoremError(f"nonfinite_coordinate:{name}")
    try:
        answer = float(value)
    except (OverflowError, ValueError) as exc:
        raise IndependentTDG7TheoremError(
            f"nonfinite_coordinate:{name}"
        ) from exc
    if not math.isfinite(answer):
        raise IndependentTDG7TheoremError(f"nonfinite_coordinate:{name}")
    return answer


def binary64_fraction(value: object, *, name: str = "value") -> Fraction:
    """Return the exact IEEE-754 binary64 rational encoded by *value*.

    This is intentionally derived from sign/exponent/significand bits rather
    than delegated to the TDG7 freeze implementation or a decimal conversion.
    """

    number = _binary64(value, name=name)
    word = _bits(number)
    sign = -1 if word >> 63 else 1
    exponent = (word >> 52) & 0x7FF
    significand = word & ((1 << 52) - 1)
    if exponent == 0x7FF:  # Defensive: _binary64 already rejects this.
        raise IndependentTDG7TheoremError(f"nonfinite_coordinate:{name}")
    if exponent == 0:
        return sign * Q(significand, 1 << 1074)
    return sign * Q((1 << 52) | significand) * _power_of_two(exponent - 1023 - 52)


def binary64_ulp_fraction(value: object, *, name: str = "value") -> Fraction:
    """Derive the upward binary64 ULP from the encoded exponent bits."""

    number = _binary64(value, name=name)
    word = _bits(number)
    exponent = (word >> 52) & 0x7FF
    if exponent == 0x7FF:
        raise IndependentTDG7TheoremError(f"nonfinite_coordinate:{name}")
    if exponent == 0:
        return Q(1, 1 << 1074)
    return _power_of_two(exponent - 1023 - 52)


def _exact_float(value: Fraction, *, reason: str = "internal_exactness_failure") -> float:
    answer = float(value)
    if not math.isfinite(answer) or binary64_fraction(answer) != value:
        raise IndependentTDG7TheoremError(reason)
    return answer


@dataclass(frozen=True, slots=True)
class HistoricalTDG6GuardCertificate:
    """Exact reconstruction of TDG6's former separately rounded guards."""

    current: float
    requested_width: float
    expected_final: float
    nominal_steps: tuple[float, float, float]
    boundaries: tuple[
        tuple[float, float],
        tuple[float, float, float],
        tuple[float, float, float, float, float],
    ]
    adjacent_differences: tuple[tuple[float, ...], tuple[float, ...], tuple[float, ...]]
    failed_counts: tuple[int, ...]


def reconstruct_historical_tdg6_guard(
    current: object, requested_width: object
) -> HistoricalTDG6GuardCertificate:
    """Recreate the historical independent ``start + i*(width/count)`` guards."""

    start = _binary64(current, name="current")
    width = _binary64(requested_width, name="requested_width")
    if width <= 0.0:
        raise IndependentTDG7TheoremError("nonpositive_requested_cap")
    expected = start + width
    nominal_steps: list[float] = []
    all_boundaries: list[tuple[float, ...]] = []
    all_differences: list[tuple[float, ...]] = []
    failed: list[int] = []
    for count in (1, 2, 4):
        step = width / count
        boundaries = tuple(start + index * step for index in range(count + 1))
        differences = tuple(right - left for left, right in zip(boundaries, boundaries[1:]))
        nominal_steps.append(step)
        all_boundaries.append(boundaries)
        all_differences.append(differences)
        if _bits(boundaries[-1]) != _bits(expected) or any(
            _bits(delta) != _bits(step) for delta in differences
        ):
            failed.append(count)
    return HistoricalTDG6GuardCertificate(
        current=start,
        requested_width=width,
        expected_final=expected,
        nominal_steps=(nominal_steps[0], nominal_steps[1], nominal_steps[2]),
        boundaries=(all_boundaries[0], all_boundaries[1], all_boundaries[2]),  # type: ignore[arg-type]
        adjacent_differences=(all_differences[0], all_differences[1], all_differences[2]),
        failed_counts=tuple(failed),
    )


def _stage_path(boundaries: Iterable[float]) -> tuple[tuple[float, float, float], ...]:
    points = tuple(boundaries)
    return tuple(
        (left, left + (right - left) / 2.0, right)
        for left, right in zip(points, points[1:])
    )


@dataclass(frozen=True, slots=True)
class IndependentBinary64LatticePlan:
    """The independently derived stage-safe 1/2/4 shared partition."""

    current: float
    target: float
    requested_cap: float
    quantum: Fraction
    macro_quantum: Fraction
    selection_budget: Fraction
    macro_width: Fraction
    fine_width: Fraction
    conservative_reduction: Fraction
    boundaries: tuple[float, float, float, float, float]
    target_limited: bool

    @property
    def outer_boundaries(self) -> tuple[float, float]:
        return (self.boundaries[0], self.boundaries[4])

    @property
    def medium_boundaries(self) -> tuple[float, float, float]:
        return (self.boundaries[0], self.boundaries[2], self.boundaries[4])

    @property
    def fine_boundaries(self) -> tuple[float, float, float, float, float]:
        return self.boundaries

    @property
    def outer_stage_times(self) -> tuple[tuple[float, float, float], ...]:
        return _stage_path(self.outer_boundaries)

    @property
    def medium_stage_times(self) -> tuple[tuple[float, float, float], ...]:
        return _stage_path(self.medium_boundaries)

    @property
    def fine_stage_times(self) -> tuple[tuple[float, float, float], ...]:
        return _stage_path(self.fine_boundaries)


@dataclass(frozen=True, slots=True)
class StageAbscissaCertificate:
    """Exactness evidence for the engine's actual ``start + dt/2`` operations."""

    stage_triplet_count: int
    all_stage_times_exact: bool
    all_midpoints_strictly_interior: bool
    all_endpoints_exact: bool


def derive_stage_safe_lattice(
    current: object,
    target: object,
    requested_cap: object,
    *,
    minimum_width: object | None = None,
) -> IndependentBinary64LatticePlan:
    """Independently derive the largest stage-safe ``8Q``-aligned macro step."""

    x_float = _binary64(current, name="current")
    target_float = _binary64(target, name="target")
    cap_float = _binary64(requested_cap, name="requested_cap")
    x = binary64_fraction(x_float)
    t = binary64_fraction(target_float)
    cap = binary64_fraction(cap_float)
    if not (TIME_LOWER <= x < t <= TIME_UPPER):
        raise IndependentTDG7TheoremError(
            "target_not_ahead" if t <= x else "outside_frozen_time_envelope"
        )
    if cap <= 0:
        raise IndependentTDG7TheoremError("nonpositive_requested_cap")
    minimum: Fraction | None = None
    if minimum_width is not None:
        minimum_float = _binary64(minimum_width, name="minimum_width")
        minimum = binary64_fraction(minimum_float)
        if minimum <= 0:
            raise IndependentTDG7TheoremError("below_minimum_aligned_width")

    quantum = max(binary64_ulp_fraction(x_float), binary64_ulp_fraction(target_float))
    macro_quantum = 8 * quantum
    if x % quantum != 0 or t % quantum != 0:
        raise IndependentTDG7TheoremError("endpoint_not_q_aligned")
    remaining = t - x
    if remaining % macro_quantum != 0:
        raise IndependentTDG7TheoremError("target_not_stage_lattice_aligned")
    budget = min(cap, remaining)
    ticks = budget // macro_quantum
    if ticks < 1:
        raise IndependentTDG7TheoremError("no_positive_aligned_width")
    macro = ticks * macro_quantum
    fine = macro / 4
    if minimum is not None and macro < minimum:
        raise IndependentTDG7TheoremError("below_minimum_aligned_width")
    if fine % (2 * quantum) != 0:
        raise IndependentTDG7TheoremError("internal_exactness_failure")

    exact_boundaries = tuple(x + index * fine for index in range(5))
    boundaries = tuple(_exact_float(value) for value in exact_boundaries)
    plan = IndependentBinary64LatticePlan(
        current=x_float,
        target=target_float,
        requested_cap=cap_float,
        quantum=quantum,
        macro_quantum=macro_quantum,
        selection_budget=budget,
        macro_width=macro,
        fine_width=fine,
        conservative_reduction=budget - macro,
        boundaries=boundaries,  # type: ignore[arg-type]
        target_limited=(macro == remaining),
    )
    verify_stage_safe_plan(plan)
    return plan


def verify_stage_safe_plan(plan: IndependentBinary64LatticePlan) -> StageAbscissaCertificate:
    """Verify every shared boundary and actual stage expression exactly."""

    if not isinstance(plan, IndependentBinary64LatticePlan):
        raise TypeError("plan must be IndependentBinary64LatticePlan")
    x = binary64_fraction(plan.current, name="current")
    t = binary64_fraction(plan.target, name="target")
    cap = binary64_fraction(plan.requested_cap, name="requested_cap")
    q = plan.quantum
    if not (TIME_LOWER <= x < t <= TIME_UPPER):
        raise IndependentTDG7TheoremError("internal_exactness_failure")
    if cap <= 0:
        raise IndependentTDG7TheoremError("internal_exactness_failure")
    if q <= 0 or q != max(binary64_ulp_fraction(plan.current), binary64_ulp_fraction(plan.target)):
        raise IndependentTDG7TheoremError("internal_exactness_failure")
    if x % q != 0 or t % q != 0 or (t - x) % (8 * q) != 0:
        raise IndependentTDG7TheoremError("internal_exactness_failure")
    if plan.macro_quantum != 8 * q or plan.macro_width <= 0 or plan.fine_width <= 0:
        raise IndependentTDG7TheoremError("internal_exactness_failure")
    if plan.macro_width % (8 * q) != 0 or plan.fine_width != plan.macro_width / 4:
        raise IndependentTDG7TheoremError("internal_exactness_failure")
    if (
        plan.fine_width % (2 * q) != 0
        or plan.selection_budget != min(cap, t - x)
        or not (0 < plan.macro_width <= plan.selection_budget <= cap)
        or plan.macro_width > t - x
    ):
        raise IndependentTDG7TheoremError("internal_exactness_failure")
    if plan.macro_width != (plan.selection_budget // (8 * q)) * (8 * q):
        raise IndependentTDG7TheoremError("internal_exactness_failure")
    if plan.conservative_reduction != plan.selection_budget - plan.macro_width:
        raise IndependentTDG7TheoremError("internal_exactness_failure")
    exact_boundaries = tuple(x + index * plan.fine_width for index in range(5))
    if (
        len(plan.boundaries) != 5
        or any(binary64_fraction(value) != expected for value, expected in zip(plan.boundaries, exact_boundaries, strict=True))
        or binary64_fraction(plan.boundaries[-1]) > t
        or plan.target_limited != (plan.macro_width == t - x)
    ):
        raise IndependentTDG7TheoremError("internal_exactness_failure")

    count = 0
    for path in (plan.outer_stage_times, plan.medium_stage_times, plan.fine_stage_times):
        for left, midpoint, right in path:
            left_exact = binary64_fraction(left)
            right_exact = binary64_fraction(right)
            midpoint_exact = binary64_fraction(midpoint)
            expected_midpoint = left_exact + (right_exact - left_exact) / 2
            if (
                left_exact >= right_exact
                or binary64_fraction(right - left) != right_exact - left_exact
                or midpoint_exact != expected_midpoint
                or not left < midpoint < right
                or binary64_fraction(left + (right - left)) != right_exact
            ):
                raise IndependentTDG7TheoremError("internal_exactness_failure")
            count += 1
    return StageAbscissaCertificate(
        stage_triplet_count=count,
        all_stage_times_exact=True,
        all_midpoints_strictly_interior=True,
        all_endpoints_exact=True,
    )


@dataclass(frozen=True, slots=True)
class EngineStageContract:
    """AST-level evidence of the immutable engine's actual time abscissae."""

    stage_time_kinds: tuple[tuple[str, str], ...]
    candidate_endpoint_uses_final_time: bool


@dataclass(frozen=True, slots=True)
class HistoricalTDG6SourceContract:
    """AST-level evidence for the exact historical subdivision guard shape."""

    step_is_width_over_count: bool
    boundaries_are_start_plus_index_times_step: bool
    expected_final_is_start_plus_width: bool
    macro_endpoint_guard_present: bool
    adjacent_uniformity_guard_present: bool


_CONTROL_FLOW_NODES = (
    ast.If,
    ast.For,
    ast.AsyncFor,
    ast.While,
    ast.Try,
    ast.With,
    ast.AsyncWith,
    ast.Match,
)


def _time_kind(node: ast.AST) -> str | None:
    if isinstance(node, ast.Name):
        return node.id if node.id in {"start", "final_time"} else None
    if (
        isinstance(node, ast.BinOp)
        and isinstance(node.op, ast.Add)
        and isinstance(node.left, ast.Name)
        and node.left.id == "start"
        and isinstance(node.right, ast.Name)
        and node.right.id == "dt"
    ):
        return "end"
    if (
        isinstance(node, ast.BinOp)
        and isinstance(node.op, ast.Add)
        and isinstance(node.left, ast.Name)
        and node.left.id == "start"
        and isinstance(node.right, ast.BinOp)
        and isinstance(node.right.op, ast.Div)
        and isinstance(node.right.left, ast.Name)
        and node.right.left.id == "dt"
        and isinstance(node.right.right, ast.Constant)
        and node.right.right.value == 2
    ):
        return "midpoint"
    return None


def parse_immutable_engine_stage_contract(engine_source: str) -> EngineStageContract:
    """Parse, rather than grep, the source-time expressions in ``propose_step``."""

    if not isinstance(engine_source, str):
        raise TypeError("engine_source must be text")
    tree = ast.parse(engine_source)
    proposal = next(
        (node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == "propose_step"),
        None,
    )
    if proposal is None:
        raise IndependentTDG7TheoremError("engine_propose_step_missing")
    method_branch = next(
        (
            node
            for node in proposal.body
            if isinstance(node, ast.If)
            and isinstance(node.test, ast.Compare)
            and isinstance(node.test.left, ast.Name)
            and node.test.left.id == "method"
            and len(node.test.ops) == 1
            and isinstance(node.test.ops[0], ast.Eq)
            and len(node.test.comparators) == 1
            and isinstance(node.test.comparators[0], ast.Name)
            and node.test.comparators[0].id == "PRIMARY_METHOD"
        ),
        None,
    )
    if method_branch is None or not method_branch.orelse:
        raise IndependentTDG7TheoremError("engine_method_branch_missing")

    def live_calls(
        statements: list[ast.stmt], *, reject_nested_control: bool
    ) -> list[tuple[str, str, str]]:
        calls: list[tuple[str, str, str]] = []
        for statement in statements:
            if reject_nested_control and any(
                isinstance(node, _CONTROL_FLOW_NODES)
                for node in ast.walk(statement)
            ):
                raise IndependentTDG7TheoremError(
                    "engine_stage_control_flow_changed"
                )
            for node in ast.walk(statement):
                if (
                    not isinstance(node, ast.Call)
                    or not isinstance(node.func, ast.Name)
                    or node.func.id != "evaluate"
                ):
                    continue
                if (
                    len(node.args) != 3
                    or not isinstance(node.args[0], ast.Constant)
                    or not isinstance(node.args[0].value, str)
                    or not isinstance(node.args[2], ast.Name)
                ):
                    raise IndependentTDG7TheoremError(
                        "engine_stage_abscissa_contract_changed"
                    )
                kind = _time_kind(node.args[1])
                if kind is None:
                    raise IndependentTDG7TheoremError(
                        "engine_stage_abscissa_contract_changed"
                    )
                calls.append((node.args[0].value, kind, node.args[2].id))
        return calls

    primary_expected = [
        ("rk4_k1", "start", "state"),
        ("rk4_k2", "midpoint", "y2"),
        ("rk4_k3", "midpoint", "y3"),
        ("rk4_k4", "end", "y4"),
    ]
    comparator_expected = [
        ("ssprk3_s0", "start", "state"),
        ("ssprk3_s1", "end", "y1"),
        ("ssprk3_s2", "midpoint", "y2"),
    ]
    branch_index = proposal.body.index(method_branch)
    prebranch_calls = live_calls(
        proposal.body[:branch_index], reject_nested_control=False
    )
    postbranch_calls = live_calls(
        proposal.body[branch_index + 1 :], reject_nested_control=True
    )
    candidate_expected = [("candidate_endpoint", "final_time", "candidate")]
    if (
        prebranch_calls
        or live_calls(method_branch.body, reject_nested_control=True)
        != primary_expected
        or live_calls(method_branch.orelse, reject_nested_control=True)
        != comparator_expected
        or postbranch_calls != candidate_expected
    ):
        raise IndependentTDG7TheoremError("engine_stage_abscissa_contract_changed")
    expected = primary_expected + comparator_expected + candidate_expected
    return EngineStageContract(
        stage_time_kinds=tuple((name, kind) for name, kind, _ in expected),
        candidate_endpoint_uses_final_time=True,
    )


def _assigned_value(statement: ast.stmt, name: str) -> ast.AST | None:
    if (
        isinstance(statement, ast.Assign)
        and len(statement.targets) == 1
        and isinstance(statement.targets[0], ast.Name)
        and statement.targets[0].id == name
    ):
        return statement.value
    return None


def _same_expression(left: ast.AST, source: str) -> bool:
    expected = ast.parse(source).body[0]
    if not isinstance(expected, ast.Expr):
        raise RuntimeError("internal AST template must be an expression")
    return ast.dump(left, include_attributes=False) == ast.dump(
        expected.value, include_attributes=False
    )


def _is_value_error_message(statements: list[ast.stmt], message: str) -> bool:
    return (
        len(statements) == 1
        and isinstance(statements[0], ast.Raise)
        and isinstance(statements[0].exc, ast.Call)
        and isinstance(statements[0].exc.func, ast.Name)
        and statements[0].exc.func.id == "ValueError"
        and len(statements[0].exc.args) == 1
        and isinstance(statements[0].exc.args[0], ast.Constant)
        and statements[0].exc.args[0].value == message
    )


def parse_historical_tdg6_subdivision_contract(
    runtime_source: str,
) -> HistoricalTDG6SourceContract:
    """Require TDG6's historical step, boundaries, and both guard shapes.

    This is deliberately stricter than a source hash.  The independent
    arithmetic reconstruction is only evidence for the historical runtime if
    the immutable source actually contains the reconstructed expressions.
    """

    if not isinstance(runtime_source, str):
        raise TypeError("runtime_source must be text")
    tree = ast.parse(runtime_source)
    function = next(
        (
            node
            for node in tree.body
            if isinstance(node, ast.FunctionDef)
            and node.name == "_subdivision_boundaries"
        ),
        None,
    )
    if function is None or len(function.body) != 6:
        raise IndependentTDG7TheoremError("historical_tdg6_subdivision_contract_changed")
    step, boundaries, expected_final, endpoint_guard, adjacent_guard, returned = function.body
    step_ok = (_assigned_value(step, "step") is not None and _same_expression(
        _assigned_value(step, "step"), "width / count"  # type: ignore[arg-type]
    ))
    boundary_ok = (_assigned_value(boundaries, "boundaries") is not None and _same_expression(
        _assigned_value(boundaries, "boundaries"),
        "tuple(start + index * step for index in range(count + 1))",  # type: ignore[arg-type]
    ))
    final_ok = (_assigned_value(expected_final, "expected_final") is not None and _same_expression(
        _assigned_value(expected_final, "expected_final"), "start + width"  # type: ignore[arg-type]
    ))
    endpoint_ok = (
        isinstance(endpoint_guard, ast.If)
        and _same_expression(
            endpoint_guard.test,
            "not _same_binary64(boundaries[-1], expected_final)",
        )
        and _is_value_error_message(
            endpoint_guard.body,
            "TDG6 subdivision does not share the macro endpoint",
        )
    )
    adjacent_ok = (
        isinstance(adjacent_guard, ast.For)
        and isinstance(adjacent_guard.target, ast.Tuple)
        and [
            item.id for item in adjacent_guard.target.elts if isinstance(item, ast.Name)
        ]
        == ["left", "right"]
        and _same_expression(adjacent_guard.iter, "zip(boundaries, boundaries[1:])")
        and len(adjacent_guard.body) == 1
        and isinstance(adjacent_guard.body[0], ast.If)
        and _same_expression(
            adjacent_guard.body[0].test,
            "not _same_binary64(right - left, step)",
        )
        and _is_value_error_message(
            adjacent_guard.body[0].body,
            "TDG6 subdivision is not bitwise uniform",
        )
    )
    returned_ok = isinstance(returned, ast.Return) and isinstance(returned.value, ast.Name) and returned.value.id == "boundaries"
    if not (step_ok and boundary_ok and final_ok and endpoint_ok and adjacent_ok and returned_ok):
        raise IndependentTDG7TheoremError("historical_tdg6_subdivision_contract_changed")
    return HistoricalTDG6SourceContract(
        step_is_width_over_count=True,
        boundaries_are_start_plus_index_times_step=True,
        expected_final_is_start_plus_width=True,
        macro_endpoint_guard_present=True,
        adjacent_uniformity_guard_present=True,
    )


def cal11_stage_safe_witness_certificate() -> dict[str, object]:
    """Recompute CAL11's historical cap, old failure, and independent 8Q plan."""

    current = 23.0 / 16.0
    target = 24.0 / 16.0
    cfl = float.fromhex("0x1.0000000000000p-3")
    spacing = float.fromhex("0x1.0000000000000p-4")
    speed = float.fromhex("0x1.33345e70b5188p+0")
    requested = min(target - current, cfl * spacing / speed)
    historical = reconstruct_historical_tdg6_guard(current, requested)
    plan = derive_stage_safe_lattice(current, target, requested)
    stage = verify_stage_safe_plan(plan)
    old_macro = binary64_fraction(historical.adjacent_differences[0][0])
    old_fine = old_macro / 4
    old_midpoint = binary64_fraction(current) + old_fine / 2
    old_midpoint_float = float(old_midpoint)
    if (
        requested.hex() != "0x1.aaa90b0fb5c26p-8"
        or historical.failed_counts != (1, 2, 4)
        or plan.quantum != Q(1, 1 << 52)
        or plan.macro_quantum != Q(1, 1 << 49)
        or plan.macro_width != Q(3664984285035, 1 << 49)
        or plan.fine_width != Q(3664984285035, 1 << 51)
        or plan.conservative_reduction != Q(531, 1 << 59)
        or plan.fine_width / plan.quantum != 7329968570070
        or binary64_fraction(old_midpoint_float) == old_midpoint
        or not stage.all_stage_times_exact
    ):
        raise IndependentTDG7TheoremError("cal11_witness_mismatch")
    return {
        "requested_cap_hex": requested.hex(),
        "historical_failed_counts": list(historical.failed_counts),
        "historical_expected_endpoint_hex": historical.expected_final.hex(),
        "event_quantum_hex": _exact_float(plan.quantum).hex(),
        "macro_quantum_hex": _exact_float(plan.macro_quantum).hex(),
        "macro_width_hex": _exact_float(plan.macro_width).hex(),
        "fine_width_hex": _exact_float(plan.fine_width).hex(),
        "fine_width_over_Q": int(plan.fine_width / plan.quantum),
        "conservative_reduction": str(plan.conservative_reduction),
        "selected_boundaries_hex": [value.hex() for value in plan.boundaries],
        "selected_fine_midpoints_hex": [value[1].hex() for value in plan.fine_stage_times],
        "endpoint_only_four_q_counterexample_derived_from_historical_one_step": True,
        "endpoint_only_four_q_midpoint_exact": False,
        "endpoint_only_four_q_midpoint_binary64_hex": old_midpoint_float.hex(),
        "stage_triplet_count": stage.stage_triplet_count,
    }


def deterministic_stage_safe_property_certificate() -> dict[str, int | bool]:
    """Small deterministic cross-binade and event-lattice theorem controls."""

    plans = 0
    triplets = 0
    for numerator in range(23, 512):
        current = numerator / 16.0
        target = (numerator + 1) / 16.0
        q = max(math.ulp(current), math.ulp(target))
        for factor in (8, 9, 31, 32, 255, 1_000_003):
            plan = derive_stage_safe_lattice(current, target, factor * 8.0 * q)
            triplets += verify_stage_safe_plan(plan).stage_triplet_count
            plans += 1
    for pivot in (2.0, 4.0, 8.0, 16.0):
        q = math.ulp(pivot)
        for low in range(-32, 1, 8):
            for high in range(8, 40, 8):
                current, target = pivot + low * q, pivot + high * q
                if not (float(TIME_LOWER) <= current < target <= float(TIME_UPPER)):
                    continue
                plan = derive_stage_safe_lattice(current, target, 127.0 * 8.0 * q)
                triplets += verify_stage_safe_plan(plan).stage_triplet_count
                plans += 1
    return {
        "successful_plan_count": plans,
        "stage_triplet_count": triplets,
        "all_exact_and_strictly_interior": True,
    }
