"""Prospective orbit-local Picard-Lindelof continuation of the affine SGB-L IVP.

This owner does not edit the Picard, Taylor, family, principal, continuity,
cell-admission, or trap instruments.  It starts from the authenticated
Misner-Sharp ``(C,k)`` prefix ending at ``r0=10749/1024`` and asks whether a
frozen tube policy proves a unique local affine ``(lambda,k)`` orbit on the
remaining declared support.

The policy is declared before evaluation: 16 equal continuation steps, tube
radius ``rho=1/8``, rational-bit cap 16384, and no adaptive bisection.
Existence is a strict Picard self-map of the orbit tube.  Uniqueness is the
interval Jacobian Lipschitz bound on that tube.  Compactness is evaluated on
the certified orbit image, not the whole product domain.

Nominal ``A_chi=3`` is immutable.  ``A_chi=1/8`` is a named control only.
Aggregate health, ``FRZ1``, ``PREF1``, execution, holdout, trapping, and
physics flags remain false.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
from json import dumps
from numbers import Integral
from types import MappingProxyType
from typing import Any, Mapping

from .exact_interval import Interval, interval
from .sgb1_ctl1_initial_health import (
    CERTIFICATE_CHART_KIND,
    DECLARED_C_DOMAIN,
    DECLARED_K_DOMAIN,
    DECLARED_LAMBDA_DOMAIN,
    DECLARED_MAX_RATIONAL_BITS,
    SGBLExactInitialSlice,
    SGBLInitialHealthStop,
    SGBLODECellEnclosure,
    SGBLValidatedODERecord,
    sgbl_centered_mean_value_compactness_rhs,
    sgbl_centered_mean_value_rhs,
    sgbl_constraint_rhs_state_jacobian,
    sgbl_exact_minkowski_rhs_reduction,
    sgbl_intersect_algebraic_compactness_enclosures,
    sgbl_interval_affine_coefficients,
    sgbl_misner_sharp_compactness_interval,
    sgbl_misner_sharp_compactness_invariant_residual,
    sgbl_validated_constraint_ode,
)
from .sgb1_ctl1_trap_refinement2 import (
    PREDECESSOR_NOMINAL_CONTRACT_SHA256 as IMPORTED_TRAP_REFINEMENT_HASH,
)
from .sgb1_ctl1_trap_taylor import (
    PREDECESSOR_TRAP_REFINEMENT2_CONTRACT_SHA256 as IMPORTED_TRAP_REFINEMENT2_HASH,
    PREDECESSOR_TRAP_REFINEMENT_CONTRACT_SHA256 as IMPORTED_TRAP_REFINEMENT_HASH_TAYLOR,
)


Q = Fraction
INSTRUMENT_ID = "FGC-1-SGB1-CTL1-SOL1"
NOMINAL_CHI_AMPLITUDE = Q(3)
SMALL_AMPLITUDE_CONTROL = Q(1, 8)
PREFIX_RADIUS = Q(10749, 1024)
PREFIX_LAST_CELL_LEFT = Q(2687, 256)
PREFIX_CELL_COUNT = 9
PREFIX_OBSTRUCTION = "picard_strict_self_map_failed"
PREFIX_C_UPPER = Q(17053332284082708107, 2**64)
PREFIX_C_MARGIN = Q(1393411789626843509, 2**64)
PREFIX_D_MARGIN = Q(326806051576253549, 2**62)
PREFIX_RHS_EVALUATIONS = 492
PREFIX_COMPACTNESS_RHS_EVALUATIONS = 123
DECLARED_SUPPORT_MAXIMUM = Q(14)
DECLARED_CONTINUATION_STEPS = 16
DECLARED_TUBE_RADIUS = Q(1, 8)
DECLARED_MAX_BISECTION_DEPTH = 0
DECLARED_MAX_RATIONAL_BIT_LENGTH = DECLARED_MAX_RATIONAL_BITS
DECLARED_DYADIC_DENOMINATOR = 2**64
DECLARED_COMPACTNESS_STRICT_UPPER = 1
DECLARED_STEP_WIDTH = (
    DECLARED_SUPPORT_MAXIMUM - PREFIX_RADIUS
) / DECLARED_CONTINUATION_STEPS
PREDECESSOR_TRAP_REFINEMENT_CONTRACT_SHA256 = (
    "a82363a29b641c538a0e15f14ec2f4c99a49182a9a01c84ae7b1d14a953834c8"
)
PREDECESSOR_TRAP_REFINEMENT2_CONTRACT_SHA256 = (
    "2cc8ba39a756afc433f1b04ce6fa3e26f45d795f2ba4c617ec080e7342cc396d"
)
PREDECESSOR_TRAP_TAYLOR_CONTRACT_SHA256 = (
    "7de1421d200f9d276b95a7e541c0dd8b4b75eddddbe966db6611316b969f820b"
)
PREDECESSOR_TRAP_BARRIER_CONTRACT_SHA256 = (
    "a75d564c0ea0658802c138e0575edf386375ebba35728661bd8a11c579e7d689"
)
AFFINE_STATE = "(lambda,k)"
PICARD_OPERATOR = "P(y)=y0+[0,h]*f(T)"
LIPSCHITZ_NORM = "||Df/d(lambda,k)||_inf"
ZERO = Interval.singleton(0)
ONE = Interval.singleton(1)
FORBIDDEN_HEALTH_IMPORTS = (
    "WeakCouplingThresholds",
    "weak_coupling_health_certificate",
    "esf_reference_cluster_certificate",
    "canonical_health_monitor_values",
    "background_from_spherical_state",
    "FGCQRActionParameters",
    "solve_initial_data",
    "run_fgc_gr0_calibration_v1",
    "pro20_ev1_runtime",
    "build_artifact_catalog",
)
TUBE_CLASSIFICATIONS = frozenset(
    {
        "unique_local_affine_orbit",
        "on_orbit_compactness_not_below_one",
        "jacobian_diagonal_contains_zero",
        "interval_inconclusive",
        "domain_error",
        "resource_limit",
    }
)
INCONCLUSIVE_REASONS = frozenset(
    {
        "prefix_endpoint_not_inside_declared_tube",
        "picard_strict_self_map_failed",
        "lipschitz_jacobian_not_bounded",
        "residual_enclosure_misses_origin",
        "compactness_invariant_misses_origin",
        "mean_value_inconsistent",
        "zero_in_interval_reciprocal",
        "incomplete_remaining_support",
        "wrapping",
    }
)
RESOURCE_REASONS = frozenset({"max_rational_bit_length"})
CONTRACT_STOP_REASONS = frozenset(
    {
        "changed_policy",
        "wrong_predecessor",
        "tampered_contract",
        "gapped_inventory",
        "claim_mutation",
        "family_mutation",
        "witness_mutation",
    }
)
MISSING_THEOREMS = MappingProxyType(
    {
        "prefix_endpoint_not_inside_declared_tube": (
            "authenticated_prefix_endpoint_contained_in_the_declared_orbit_tube"
        ),
        "picard_strict_self_map_failed": (
            "strict_Picard_self_map_of_the_affine_constraint_IVP_on_the_declared_orbit_tube"
        ),
        "jacobian_diagonal_contains_zero": (
            "affine_H_L_and_M_k_diagonals_excluding_zero_on_every_declared_orbit_tube"
        ),
        "on_orbit_compactness_not_below_one": (
            "on_orbit_compactness_C_strictly_below_one_along_the_certified_affine_orbit"
        ),
        "lipschitz_jacobian_not_bounded": (
            "interval_Jacobian_Lipschitz_bound_of_the_affine_constraint_vector_field"
        ),
        "residual_enclosure_misses_origin": (
            "affine_H_and_M_residual_enclosures_containing_the_origin_on_every_orbit_tube"
        ),
        "compactness_invariant_misses_origin": (
            "algebraic_and_propagated_orbit_compactness_overlapping_with_invariant_residual_zero"
        ),
        "mean_value_inconsistent": (
            "centered_mean_value_RHS_intersecting_the_natural_interval_extension"
        ),
        "zero_in_interval_reciprocal": (
            "interval_reciprocals_of_r_and_lambda_remaining_defined_on_every_orbit_tube"
        ),
        "incomplete_remaining_support": (
            "unique_affine_orbit_tubes_tiling_the_remaining_declared_support"
        ),
        "wrapping": (
            "orbit_local_Picard_Lindelof_tube_containing_the_authenticated_prefix_endpoint"
        ),
        "physical_domain_escape": (
            "Picard_image_remaining_inside_the_declared_physical_state_domain"
        ),
        "max_rational_bit_length": "declared_orbit_tube_rational_bit_budget",
    }
)
FALSE_CLAIMS = (
    "SGBL_branch_owned_and_healthy",
    "FRZ1",
    "PREF1",
    "holdout_authorized",
    "execution_authorized",
    "whole_domain_nagumo_invariant",
    "trajectory_trapping_claimed",
    "actual_orbit_traps",
    "sgbl_model_rejected",
    "copied_gr0_or_fgcqr_health_evidence",
    "adapted_tube_radius",
    "adapted_step_count",
    "replaced_nominal_family_by_small_amplitude",
    "spacetime_ibvp_theorem",
    "quasilinear_existence",
)

if IMPORTED_TRAP_REFINEMENT_HASH != PREDECESSOR_TRAP_REFINEMENT_CONTRACT_SHA256:
    raise RuntimeError("SOL1 requires the preserved TRAP-REFINEMENT contract hash")
if IMPORTED_TRAP_REFINEMENT_HASH_TAYLOR != PREDECESSOR_TRAP_REFINEMENT_CONTRACT_SHA256:
    raise RuntimeError("SOL1 requires Taylor to bind the same TRAP-REFINEMENT hash")
if IMPORTED_TRAP_REFINEMENT2_HASH != PREDECESSOR_TRAP_REFINEMENT2_CONTRACT_SHA256:
    raise RuntimeError("SOL1 requires the preserved TRAP-REFINEMENT2 contract hash")


def _fraction(name: str, value: object, *, nonnegative: bool = False) -> Fraction:
    if type(value) not in (int, Fraction):
        raise TypeError(f"{name} must be a Fraction or built-in int")
    result = Fraction(value)
    if nonnegative and result < 0:
        raise ValueError(f"{name} must be nonnegative")
    return result


def _positive_int(name: str, value: object, *, allow_zero: bool = False) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise ValueError(f"{name} must be an integer")
    result = int(value)
    if result < 0 or (result == 0 and not allow_zero):
        raise ValueError(f"{name} must be a positive integer")
    return result


def _as_interval(value: object) -> Interval:
    if type(value) is Interval:
        return value
    if type(value) in (int, Fraction):
        return Interval.singleton(Fraction(value))
    raise TypeError("value must be an Interval or exact rational")


def _fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def _interval_text(value: Interval) -> tuple[str, str]:
    return (_fraction_text(value.lower), _fraction_text(value.upper))


def _interval_bits(value: Interval) -> int:
    return max(
        abs(value.lower.numerator).bit_length(),
        value.lower.denominator.bit_length(),
        abs(value.upper.numerator).bit_length(),
        value.upper.denominator.bit_length(),
    )


def _floor_dyadic(value: Fraction, denominator: int) -> Fraction:
    return Fraction((value.numerator * denominator) // value.denominator, denominator)


def _ceil_dyadic(value: Fraction, denominator: int) -> Fraction:
    return -_floor_dyadic(-value, denominator)


def _outward(value: Interval, denominator: int = DECLARED_DYADIC_DENOMINATOR) -> Interval:
    """Valid outward rounding onto the declared dyadic grid.  This is a resource."""

    return Interval(_floor_dyadic(value.lower, denominator), _ceil_dyadic(value.upper, denominator))


def _intersect(left: Interval, right: Interval) -> Interval | None:
    lower = max(left.lower, right.lower)
    upper = min(left.upper, right.upper)
    if lower > upper:
        return None
    return Interval(lower, upper)


class SGBLOrbitLocalStop(ValueError):
    """Typed contract stop of the prospective orbit-local Picard-Lindelof owner."""

    def __init__(
        self,
        reason: str,
        message: str,
        payload: Mapping[str, Any] | None = None,
    ) -> None:
        if reason not in CONTRACT_STOP_REASONS:
            raise ValueError("unknown orbit-local contract stop reason")
        self.reason = reason
        self.payload = dict(payload or {})
        super().__init__(message)


@dataclass
class _TubeBudget:
    rhs_evaluations: int = 0
    compactness_rhs_evaluations: int = 0
    max_bits: int = 0


def _observe(budget: _TubeBudget, limit: int, *values: Interval) -> None:
    for value in values:
        budget.max_bits = max(budget.max_bits, _interval_bits(value))
        if budget.max_bits > limit:
            raise SGBLInitialHealthStop(
                "resource_limit",
                "orbit-tube endpoints exceeded the declared rational bit cap",
                {
                    "observed": budget.max_bits,
                    "limit": limit,
                    "resource_reason": "max_rational_bit_length",
                },
            )


def _outward_all(budget: _TubeBudget, limit: int, *values: Interval) -> tuple[Interval, ...]:
    rounded = tuple(_outward(value) for value in values)
    _observe(budget, limit, *rounded)
    return rounded


@dataclass(frozen=True, slots=True, kw_only=True)
class SGBLOrbitTubePolicy:
    """Declared orbit-tube continuation policy.  Bind requires the frozen tuple."""

    continuation_steps: int = DECLARED_CONTINUATION_STEPS
    tube_radius: Fraction = DECLARED_TUBE_RADIUS
    max_bisection_depth: int = DECLARED_MAX_BISECTION_DEPTH
    max_rational_bit_length: int = DECLARED_MAX_RATIONAL_BIT_LENGTH
    prefix_radius: Fraction = PREFIX_RADIUS
    support_maximum: Fraction = DECLARED_SUPPORT_MAXIMUM
    lambda_domain: Interval = DECLARED_LAMBDA_DOMAIN
    k_domain: Interval = DECLARED_K_DOMAIN
    compactness_domain: Interval = DECLARED_C_DOMAIN

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "continuation_steps",
            _positive_int("continuation_steps", self.continuation_steps),
        )
        object.__setattr__(self, "tube_radius", _fraction("tube_radius", self.tube_radius))
        object.__setattr__(
            self,
            "max_bisection_depth",
            _positive_int(
                "max_bisection_depth", self.max_bisection_depth, allow_zero=True
            ),
        )
        object.__setattr__(
            self,
            "max_rational_bit_length",
            _positive_int("max_rational_bit_length", self.max_rational_bit_length),
        )
        object.__setattr__(
            self, "prefix_radius", _fraction("prefix_radius", self.prefix_radius)
        )
        object.__setattr__(
            self, "support_maximum", _fraction("support_maximum", self.support_maximum)
        )
        if type(self.lambda_domain) is not Interval or type(self.k_domain) is not Interval:
            raise TypeError("orbit-tube domains must be exact intervals")
        if type(self.compactness_domain) is not Interval:
            raise TypeError("orbit-tube compactness domain must be an exact interval")
        if self.tube_radius <= 0:
            raise ValueError("tube radius must be positive")
        if self.support_maximum <= self.prefix_radius:
            raise ValueError("remaining support must have positive length")

    @property
    def step_width(self) -> Fraction:
        return (self.support_maximum - self.prefix_radius) / self.continuation_steps


DECLARED_ORBIT_TUBE_POLICY = SGBLOrbitTubePolicy()
if DECLARED_ORBIT_TUBE_POLICY.step_width != DECLARED_STEP_WIDTH:
    raise RuntimeError("declared orbit-tube step width must be (14-10749/1024)/16")
if DECLARED_STEP_WIDTH != Q(3587, 16384):
    raise RuntimeError("declared remaining-support step must be the exact 3587/16384")
if 1 - PREFIX_C_UPPER != PREFIX_C_MARGIN:
    raise RuntimeError("frozen prefix C_margin must be the exact 1-C_upper remainder")


def sgbl_declared_orbit_step_grid(
    *,
    origin: Fraction | int = PREFIX_RADIUS,
    terminus: Fraction | int = DECLARED_SUPPORT_MAXIMUM,
    steps: int = DECLARED_CONTINUATION_STEPS,
) -> tuple[Fraction, ...]:
    """Exact equal partition of the remaining declared support."""

    origin = _fraction("origin", origin)
    terminus = _fraction("terminus", terminus)
    steps = _positive_int("steps", steps)
    width = (terminus - origin) / steps
    nodes = tuple(origin + index * width for index in range(steps + 1))
    if nodes[-1] != terminus:
        raise SGBLOrbitLocalStop(
            "changed_policy",
            "orbit-tube grid must terminate at the declared support maximum",
            {"terminus": _fraction_text(nodes[-1])},
        )
    return nodes


def sgbl_orbit_tube_around(state: Interval, radius: Fraction | int) -> Interval:
    """Closed infinity-norm tube of the declared radius about the midpoint."""

    state = _as_interval(state)
    radius = _fraction("radius", radius)
    if radius <= 0:
        raise ValueError("tube radius must be positive")
    midpoint = state.midpoint()
    return Interval(midpoint - radius, midpoint + radius)


def sgbl_interval_matrix_infinity_norm_2x2(
    j11: Interval,
    j12: Interval,
    j21: Interval,
    j22: Interval,
) -> Fraction:
    """Infinity-norm Lipschitz bound of a 2x2 interval Jacobian."""

    return max(
        j11.abs_upper() + j12.abs_upper(),
        j21.abs_upper() + j22.abs_upper(),
    )


def sgbl_require_affine_jacobian_diagonals(
    hamiltonian_lambda_r: Interval,
    momentum_k_r: Interval,
) -> None:
    """Refuse a tube whose affine H/M Jacobian diagonal contains zero."""

    hamiltonian_lambda_r = _as_interval(hamiltonian_lambda_r)
    momentum_k_r = _as_interval(momentum_k_r)
    if hamiltonian_lambda_r.contains_zero() or momentum_k_r.contains_zero():
        raise SGBLInitialHealthStop(
            "interval_inconclusive",
            "affine constraint Jacobian diagonal contains zero",
            {
                "obstruction": "jacobian_diagonal_contains_zero",
                "hamiltonian_lambda_r": (
                    hamiltonian_lambda_r.lower,
                    hamiltonian_lambda_r.upper,
                ),
                "momentum_k_r": (momentum_k_r.lower, momentum_k_r.upper),
            },
        )


def sgbl_classify_on_orbit_compactness(compactness: Interval) -> str:
    """Classify C along a certified orbit image.  The product box is not used."""

    compactness = _as_interval(compactness)
    if compactness.upper < DECLARED_COMPACTNESS_STRICT_UPPER:
        return "unique_local_affine_orbit"
    return "on_orbit_compactness_not_below_one"


def _require_declared_policy(policy: object) -> SGBLOrbitTubePolicy:
    if type(policy) is not SGBLOrbitTubePolicy:
        raise SGBLOrbitLocalStop(
            "changed_policy",
            "orbit-tube policy must be SGBLOrbitTubePolicy",
        )
    if policy != DECLARED_ORBIT_TUBE_POLICY:
        raise SGBLOrbitLocalStop(
            "changed_policy",
            "orbit-tube policy is not the frozen declaration",
            {
                "continuation_steps": policy.continuation_steps,
                "tube_radius": _fraction_text(policy.tube_radius),
                "max_bisection_depth": policy.max_bisection_depth,
                "max_rational_bit_length": policy.max_rational_bit_length,
            },
        )
    if policy.max_bisection_depth != 0:
        raise SGBLOrbitLocalStop(
            "changed_policy",
            "orbit-local continuation does not adaptively bisect",
        )
    return policy


def _require_predecessor_hashes(
    trap_refinement: object,
    trap_refinement2: object,
    trap_taylor: object,
    trap_barrier: object,
) -> tuple[str, str, str, str]:
    required = (
        (
            trap_refinement,
            PREDECESSOR_TRAP_REFINEMENT_CONTRACT_SHA256,
            "TRAP-REFINEMENT",
        ),
        (
            trap_refinement2,
            PREDECESSOR_TRAP_REFINEMENT2_CONTRACT_SHA256,
            "TRAP-REFINEMENT2",
        ),
        (
            trap_taylor,
            PREDECESSOR_TRAP_TAYLOR_CONTRACT_SHA256,
            "TRAP-TAYLOR",
        ),
        (
            trap_barrier,
            PREDECESSOR_TRAP_BARRIER_CONTRACT_SHA256,
            "TRAP-BARRIER",
        ),
    )
    values: list[str] = []
    for supplied, expected, name in required:
        if type(supplied) is not str or supplied != expected:
            raise SGBLOrbitLocalStop(
                "wrong_predecessor",
                f"SOL1 does not bind the preserved {name} contract hash",
                {"supplied": supplied, "required": expected},
            )
        values.append(supplied)
    return values[0], values[1], values[2], values[3]


def _require_nominal_family(spec: object) -> SGBLExactInitialSlice:
    if type(spec) is not SGBLExactInitialSlice:
        raise TypeError("spec must be SGBLExactInitialSlice")
    if spec.chi_amplitude != NOMINAL_CHI_AMPLITUDE:
        raise SGBLOrbitLocalStop(
            "family_mutation",
            "nominal bind requires the frozen A_chi=3 family",
            {
                "chi_amplitude": _fraction_text(spec.chi_amplitude),
                "expected": _fraction_text(NOMINAL_CHI_AMPLITUDE),
            },
        )
    nominal = SGBLExactInitialSlice()
    if spec != nominal:
        raise SGBLOrbitLocalStop(
            "family_mutation",
            "nominal bind cannot retune the matched family fixtures",
        )
    return spec


def _require_false_claims(claims: Mapping[str, Any] | None) -> dict[str, bool]:
    flags = {name: False for name in FALSE_CLAIMS}
    if claims is None:
        return flags
    if not isinstance(claims, Mapping):
        raise TypeError("claims must be a mapping")
    for name, value in claims.items():
        if name not in flags:
            raise SGBLOrbitLocalStop(
                "claim_mutation",
                "unknown promoting claim is not part of the SOL1 contract",
                {"claim": name},
            )
        if value is not False:
            raise SGBLOrbitLocalStop(
                "claim_mutation",
                "SOL1 contract cannot promote aggregate, freeze, or physics claims",
                {"claim": name, "value": value},
            )
        flags[name] = False
    return flags


def sgbl_authenticate_ck_prefix(
    ode: SGBLValidatedODERecord,
    *,
    spec: SGBLExactInitialSlice | None = None,
) -> SGBLODECellEnclosure:
    """Bind the live ``(C,k)`` prefix to the frozen initial-health identities."""

    if type(ode) is not SGBLValidatedODERecord:
        raise TypeError("ode must be SGBLValidatedODERecord")
    spec = _require_nominal_family(spec if spec is not None else ode.spec)
    if ode.spec != spec:
        raise SGBLOrbitLocalStop(
            "wrong_predecessor",
            "prefix ODE family does not match the frozen A_chi=3 slice",
        )
    if ode.chart_kind != CERTIFICATE_CHART_KIND or not ode.chart_coordinates:
        raise SGBLOrbitLocalStop(
            "wrong_predecessor",
            "SOL1 starts from the authenticated (C,k) certificate prefix",
            {"chart_kind": ode.chart_kind},
        )
    if ode.obstruction != PREFIX_OBSTRUCTION:
        raise SGBLOrbitLocalStop(
            "wrong_predecessor",
            "authenticated prefix obstruction is not picard_strict_self_map_failed",
            {"obstruction": ode.obstruction},
        )
    if len(ode.cells) != PREFIX_CELL_COUNT:
        raise SGBLOrbitLocalStop(
            "wrong_predecessor",
            "authenticated prefix cell count is not the frozen nine-cell inventory",
            {"cells": len(ode.cells)},
        )
    if ode.coverage_left != spec.support_minimum or ode.coverage_right != PREFIX_RADIUS:
        raise SGBLOrbitLocalStop(
            "wrong_predecessor",
            "authenticated prefix coverage is not [10, 10749/1024]",
            {
                "coverage_left": _fraction_text(ode.coverage_left),
                "coverage_right": _fraction_text(ode.coverage_right),
            },
        )
    if ode.support_compactness.upper != PREFIX_C_UPPER:
        raise SGBLOrbitLocalStop(
            "wrong_predecessor",
            "authenticated prefix C_upper is not the frozen (C,k) identity",
            {"c_upper": _fraction_text(ode.support_compactness.upper)},
        )
    if ode.rhs_evaluations != PREFIX_RHS_EVALUATIONS:
        raise SGBLOrbitLocalStop(
            "wrong_predecessor",
            "authenticated prefix affine-RHS count is not the frozen identity",
            {"rhs_evaluations": ode.rhs_evaluations},
        )
    if ode.compactness_rhs_evaluations != PREFIX_COMPACTNESS_RHS_EVALUATIONS:
        raise SGBLOrbitLocalStop(
            "wrong_predecessor",
            "authenticated prefix C_r count is not the frozen identity",
            {"compactness_rhs_evaluations": ode.compactness_rhs_evaluations},
        )
    last = ode.cells[-1]
    if last.radius.lower != PREFIX_LAST_CELL_LEFT or last.radius.upper != PREFIX_RADIUS:
        raise SGBLOrbitLocalStop(
            "wrong_predecessor",
            "last prefix cell is not [2687/256, 10749/1024]",
            {
                "left": _fraction_text(last.radius.lower),
                "right": _fraction_text(last.radius.upper),
            },
        )
    if last.denominator_margin != PREFIX_D_MARGIN:
        raise SGBLOrbitLocalStop(
            "wrong_predecessor",
            "authenticated prefix D_margin is not the frozen identity",
            {"d_margin": _fraction_text(last.denominator_margin)},
        )
    return last


def _state_at_prefix_radius(ode: SGBLValidatedODERecord) -> tuple[Interval, Interval]:
    """Enclose ``(lambda,k)`` at ``r0`` by the unique covering cell of a graph."""

    if type(ode) is not SGBLValidatedODERecord:
        raise TypeError("ode must be SGBLValidatedODERecord")
    covering = [
        cell
        for cell in ode.cells
        if cell.radius.lower <= PREFIX_RADIUS <= cell.radius.upper
    ]
    if not covering:
        raise SGBLOrbitLocalStop(
            "wrong_predecessor",
            "ODE inventory does not cover the frozen prefix radius",
        )
    cell = covering[0]
    if PREFIX_RADIUS == cell.radius.upper:
        return cell.lambda_right, cell.k_right
    if PREFIX_RADIUS == cell.radius.lower:
        return cell.lambda_left, cell.k_left
    return cell.lambda_box, cell.k_box


def _missing_theorem(reason: str | None) -> str | None:
    if reason is None:
        return None
    return MISSING_THEOREMS.get(reason, reason)


@dataclass(frozen=True, slots=True, kw_only=True)
class SGBLOrbitTubeEnclosure:
    """One frozen-radius Picard-Lindelof tube of the affine ``(lambda,k)`` IVP."""

    radius: Interval
    lambda_left: Interval
    k_left: Interval
    lambda_tube: Interval
    k_tube: Interval
    lambda_image: Interval
    k_image: Interval
    lambda_right: Interval
    k_right: Interval
    compactness: Interval
    direct_compactness: Interval
    invariant_residual: Interval
    hamiltonian_lambda_r: Interval
    momentum_k_r: Interval
    lipschitz: Fraction
    max_speed: Fraction
    contraction_constant: Fraction
    initial_inside_tube: bool
    picard_self_map: bool
    jacobian_diagonals_exclude_zero: bool
    residual_contains_origin: bool
    compactness_invariant_contains_zero: bool
    unique_local_affine_orbit: bool
    classification: str
    obstruction: str | None
    missing_theorem: str | None
    rhs_evaluations: int
    compactness_rhs_evaluations: int
    max_rational_bit_length_observed: int
    depth: int
    chart_kind: str = AFFINE_STATE

    def __post_init__(self) -> None:
        if self.classification not in TUBE_CLASSIFICATIONS:
            raise ValueError("unknown orbit-tube classification")
        if self.depth != 0:
            raise SGBLOrbitLocalStop(
                "changed_policy",
                "orbit-local tubes do not adaptively bisect",
                {"depth": self.depth},
            )
        unique = (
            self.initial_inside_tube
            and self.picard_self_map
            and self.jacobian_diagonals_exclude_zero
            and self.residual_contains_origin
            and self.compactness_invariant_contains_zero
            and self.classification
            in {"unique_local_affine_orbit", "on_orbit_compactness_not_below_one"}
        )
        if self.unique_local_affine_orbit != unique:
            raise SGBLOrbitLocalStop(
                "tampered_contract",
                "unique_local_affine_orbit does not match the Picard-Lindelof predicates",
            )
        if (
            self.classification == "unique_local_affine_orbit"
            and self.compactness.upper >= DECLARED_COMPACTNESS_STRICT_UPPER
        ):
            raise SGBLOrbitLocalStop(
                "tampered_contract",
                "a unique-orbit pass cannot retain C>=1 along the certified image",
            )
        if self.classification == "on_orbit_compactness_not_below_one" and not unique:
            raise SGBLOrbitLocalStop(
                "tampered_contract",
                "on-orbit compactness nonpass requires a unique certified orbit",
            )


def sgbl_validate_orbit_tube_inventory(
    tubes: tuple[SGBLOrbitTubeEnclosure, ...],
    *,
    origin: Fraction,
    terminus: Fraction | None = None,
) -> tuple[Fraction, Fraction]:
    """Require a strictly ordered, nonoverlapping, gap-free tube inventory."""

    origin = _fraction("origin", origin)
    if terminus is not None:
        terminus = _fraction("terminus", terminus)
    if not tubes:
        if terminus is not None and terminus != origin:
            raise ValueError("empty orbit inventory cannot cover a nonempty interval")
        return origin, origin
    if any(type(tube) is not SGBLOrbitTubeEnclosure for tube in tubes):
        raise TypeError("tubes must be SGBLOrbitTubeEnclosure")
    if tubes[0].radius.lower != origin:
        raise ValueError("orbit inventory does not start at the requested left boundary")
    previous_right = origin
    previous: SGBLOrbitTubeEnclosure | None = None
    for tube in tubes:
        if tube.radius.upper <= tube.radius.lower:
            raise ValueError("orbit tube radius is not a positive-length interval")
        if tube.radius.lower < previous_right:
            raise ValueError("orbit inventory contains overlapping tubes")
        if tube.radius.lower > previous_right:
            raise ValueError("orbit inventory contains a radial gap")
        if previous is not None and (
            tube.lambda_left != previous.lambda_right or tube.k_left != previous.k_right
        ):
            raise ValueError("orbit inventory does not concatenate exact neighboring endpoints")
        previous_right = tube.radius.upper
        previous = tube
    if terminus is not None and previous_right != terminus:
        raise ValueError("orbit inventory does not reach the required right boundary")
    return origin, previous_right


def sgbl_orbit_inventory_span(
    tubes: tuple[SGBLOrbitTubeEnclosure, ...],
    *,
    origin: Fraction = PREFIX_RADIUS,
) -> tuple[Fraction, Fraction]:
    try:
        return sgbl_validate_orbit_tube_inventory(tubes, origin=origin)
    except ValueError as exc:
        raise SGBLOrbitLocalStop(
            "gapped_inventory",
            "SOL1 refuses a dropped or gapped orbit-tube inventory",
            {"detail": str(exc)},
        ) from exc


def _classify_tube(
    *,
    unique: bool,
    compactness: Interval,
    obstruction: str | None,
    resource_reason: str | None,
    domain: bool,
) -> tuple[str, str | None]:
    if resource_reason is not None:
        return "resource_limit", resource_reason
    if domain:
        return "domain_error", obstruction or "physical_domain_escape"
    if obstruction == "jacobian_diagonal_contains_zero":
        return "jacobian_diagonal_contains_zero", obstruction
    if unique:
        if compactness.upper < DECLARED_COMPACTNESS_STRICT_UPPER:
            return "unique_local_affine_orbit", None
        return "on_orbit_compactness_not_below_one", "on_orbit_compactness_not_below_one"
    return "interval_inconclusive", obstruction or "picard_strict_self_map_failed"


def sgbl_orbit_tube_picard_lindelof(
    r0: Fraction | int,
    r1: Fraction | int,
    lambda_left: Interval | Fraction | int,
    k_left: Interval | Fraction | int,
    spec: SGBLExactInitialSlice,
    *,
    policy: SGBLOrbitTubePolicy | None = None,
    budget: _TubeBudget | None = None,
) -> SGBLOrbitTubeEnclosure:
    """Prove or refuse a unique affine orbit on one frozen-radius tube."""

    if type(spec) is not SGBLExactInitialSlice:
        raise TypeError("spec must be SGBLExactInitialSlice")
    policy = policy or DECLARED_ORBIT_TUBE_POLICY
    if type(policy) is not SGBLOrbitTubePolicy:
        raise TypeError("policy must be SGBLOrbitTubePolicy")
    if policy.max_bisection_depth != 0:
        raise SGBLOrbitLocalStop(
            "changed_policy",
            "orbit-local continuation does not adaptively bisect",
        )
    r0 = _fraction("r0", r0)
    r1 = _fraction("r1", r1)
    lambda_left = _as_interval(lambda_left)
    k_left = _as_interval(k_left)
    budget = budget or _TubeBudget()
    width = r1 - r0
    radius = interval(r0, r1)
    if width <= 0 or not radius.strictly_positive():
        return _failed_tube(
            radius=radius,
            lambda_left=lambda_left,
            k_left=k_left,
            policy=policy,
            budget=budget,
            classification="domain_error",
            obstruction="physical_domain_escape",
            domain=True,
        )
    lambda_tube = _outward(sgbl_orbit_tube_around(lambda_left, policy.tube_radius))
    k_tube = _outward(sgbl_orbit_tube_around(k_left, policy.tube_radius))
    try:
        _observe(budget, policy.max_rational_bit_length, lambda_tube, k_tube, radius)
    except SGBLInitialHealthStop:
        return _failed_tube(
            radius=radius,
            lambda_left=lambda_left,
            k_left=k_left,
            policy=policy,
            budget=budget,
            classification="resource_limit",
            obstruction="max_rational_bit_length",
            resource_reason="max_rational_bit_length",
            lambda_tube=lambda_tube,
            k_tube=k_tube,
        )
    initial_inside = lambda_left.subset_of(lambda_tube) and k_left.subset_of(k_tube)
    if not initial_inside:
        return _failed_tube(
            radius=radius,
            lambda_left=lambda_left,
            k_left=k_left,
            policy=policy,
            budget=budget,
            classification="interval_inconclusive",
            obstruction="prefix_endpoint_not_inside_declared_tube",
            lambda_tube=lambda_tube,
            k_tube=k_tube,
            initial_inside_tube=False,
        )
    if not lambda_tube.subset_of(policy.lambda_domain) or not k_tube.subset_of(
        policy.k_domain
    ):
        return _failed_tube(
            radius=radius,
            lambda_left=lambda_left,
            k_left=k_left,
            policy=policy,
            budget=budget,
            classification="domain_error",
            obstruction="physical_domain_escape",
            domain=True,
            lambda_tube=lambda_tube,
            k_tube=k_tube,
            initial_inside_tube=True,
        )
    if not lambda_tube.strictly_positive():
        return _failed_tube(
            radius=radius,
            lambda_left=lambda_left,
            k_left=k_left,
            policy=policy,
            budget=budget,
            classification="domain_error",
            obstruction="physical_domain_escape",
            domain=True,
            lambda_tube=lambda_tube,
            k_tube=k_tube,
            initial_inside_tube=True,
        )
    try:
        coefficients = sgbl_interval_affine_coefficients(
            radius=radius,
            radial_metric=lambda_tube,
            angular_extrinsic_curvature=k_tube,
            spec=spec,
        )
        budget.rhs_evaluations += 1
        hamiltonian_constant = coefficients["hamiltonian_constant"]
        hamiltonian_lambda_r = coefficients["hamiltonian_lambda_r"]
        momentum_constant = coefficients["momentum_constant"]
        momentum_k_r = coefficients["momentum_k_r"]
        (
            hamiltonian_constant,
            hamiltonian_lambda_r,
            momentum_constant,
            momentum_k_r,
        ) = _outward_all(
            budget,
            policy.max_rational_bit_length,
            hamiltonian_constant,
            hamiltonian_lambda_r,
            momentum_constant,
            momentum_k_r,
        )
        sgbl_require_affine_jacobian_diagonals(hamiltonian_lambda_r, momentum_k_r)
        natural_lambda_r = _outward(-hamiltonian_constant / hamiltonian_lambda_r)
        natural_k_r = _outward(-momentum_constant / momentum_k_r)
        centered_lambda_r, centered_k_r = sgbl_centered_mean_value_rhs(
            radius, lambda_tube, k_tube, spec
        )
        budget.rhs_evaluations += 3
        jacobian = sgbl_constraint_rhs_state_jacobian(
            radius, lambda_tube, k_tube, spec
        )
        budget.rhs_evaluations += 2
        centered_lambda_r, centered_k_r = _outward_all(
            budget,
            policy.max_rational_bit_length,
            centered_lambda_r,
            centered_k_r,
        )
        lambda_r = _intersect(natural_lambda_r, centered_lambda_r)
        k_r = _intersect(natural_k_r, centered_k_r)
        if lambda_r is None or k_r is None:
            return _failed_tube(
                radius=radius,
                lambda_left=lambda_left,
                k_left=k_left,
                policy=policy,
                budget=budget,
                classification="interval_inconclusive",
                obstruction="mean_value_inconsistent",
                lambda_tube=lambda_tube,
                k_tube=k_tube,
                initial_inside_tube=True,
                hamiltonian_lambda_r=hamiltonian_lambda_r,
                momentum_k_r=momentum_k_r,
                jacobian_diagonals_exclude_zero=True,
            )
        lambda_r, k_r = _outward_all(budget, policy.max_rational_bit_length, lambda_r, k_r)
        lipschitz = sgbl_interval_matrix_infinity_norm_2x2(
            jacobian["d_lambda_r_d_lambda"],
            jacobian["d_lambda_r_d_k"],
            jacobian["d_k_r_d_lambda"],
            jacobian["d_k_r_d_k"],
        )
        max_speed = max(lambda_r.abs_upper(), k_r.abs_upper())
        step = interval(0, width)
        lambda_image = _outward(lambda_left + step * lambda_r)
        k_image = _outward(k_left + step * k_r)
        lambda_right = _outward(lambda_left + Interval.singleton(width) * lambda_r)
        k_right = _outward(k_left + Interval.singleton(width) * k_r)
        _observe(
            budget,
            policy.max_rational_bit_length,
            lambda_image,
            k_image,
            lambda_right,
            k_right,
        )
        hamiltonian_residual = _outward(hamiltonian_constant + hamiltonian_lambda_r * lambda_r)
        momentum_residual = _outward(momentum_constant + momentum_k_r * k_r)
        residual_ok = hamiltonian_residual.contains_zero() and momentum_residual.contains_zero()
        self_map = lambda_image.strictly_inside(lambda_tube) and k_image.strictly_inside(
            k_tube
        )
        compactness_r = sgbl_centered_mean_value_compactness_rhs(
            radius, lambda_tube, k_tube, spec
        )
        budget.compactness_rhs_evaluations += 1
        compactness_r = _outward(compactness_r)
        if not lambda_image.strictly_positive():
            return _failed_tube(
                radius=radius,
                lambda_left=lambda_left,
                k_left=k_left,
                policy=policy,
                budget=budget,
                classification="domain_error",
                obstruction="physical_domain_escape",
                domain=True,
                lambda_tube=lambda_tube,
                k_tube=k_tube,
                lambda_image=lambda_image,
                k_image=k_image,
                lambda_right=lambda_right,
                k_right=k_right,
                initial_inside_tube=True,
                picard_self_map=self_map,
                hamiltonian_lambda_r=hamiltonian_lambda_r,
                momentum_k_r=momentum_k_r,
                jacobian_diagonals_exclude_zero=True,
                lipschitz=lipschitz,
                max_speed=max_speed,
            )
        compactness_left = sgbl_misner_sharp_compactness_interval(
            Interval.singleton(r0), lambda_left, k_left
        )
        propagated = _outward(compactness_left + step * compactness_r)
        enclosures = sgbl_intersect_algebraic_compactness_enclosures(
            radius,
            lambda_image,
            k_image,
            propagated=propagated,
        )
        direct = _outward(enclosures["natural"])
        compactness = _outward(enclosures["correlated"])
        invariant = _outward(
            sgbl_misner_sharp_compactness_invariant_residual(
                radius, lambda_image, k_image, compactness
            )
        )
        _observe(budget, policy.max_rational_bit_length, compactness, direct, invariant)
        invariant_ok = invariant.contains_zero()
        unique = initial_inside and self_map and residual_ok and invariant_ok
        classification, obstruction = _classify_tube(
            unique=unique,
            compactness=compactness,
            obstruction=None
            if unique
            else (
                "residual_enclosure_misses_origin"
                if not residual_ok
                else (
                    "compactness_invariant_misses_origin"
                    if not invariant_ok
                    else "picard_strict_self_map_failed"
                )
            ),
            resource_reason=None,
            domain=False,
        )
        contraction = width * lipschitz
        return SGBLOrbitTubeEnclosure(
            radius=radius,
            lambda_left=lambda_left,
            k_left=k_left,
            lambda_tube=lambda_tube,
            k_tube=k_tube,
            lambda_image=lambda_image,
            k_image=k_image,
            lambda_right=lambda_right,
            k_right=k_right,
            compactness=compactness,
            direct_compactness=direct,
            invariant_residual=invariant,
            hamiltonian_lambda_r=hamiltonian_lambda_r,
            momentum_k_r=momentum_k_r,
            lipschitz=lipschitz,
            max_speed=max_speed,
            contraction_constant=contraction,
            initial_inside_tube=True,
            picard_self_map=self_map,
            jacobian_diagonals_exclude_zero=True,
            residual_contains_origin=residual_ok,
            compactness_invariant_contains_zero=invariant_ok,
            unique_local_affine_orbit=unique,
            classification=classification,
            obstruction=obstruction,
            missing_theorem=_missing_theorem(obstruction),
            rhs_evaluations=budget.rhs_evaluations,
            compactness_rhs_evaluations=budget.compactness_rhs_evaluations,
            max_rational_bit_length_observed=budget.max_bits,
            depth=0,
        )
    except SGBLInitialHealthStop as exc:
        if exc.reason == "resource_limit":
            return _failed_tube(
                radius=radius,
                lambda_left=lambda_left,
                k_left=k_left,
                policy=policy,
                budget=budget,
                classification="resource_limit",
                obstruction=str(
                    exc.payload.get("resource_reason") or "max_rational_bit_length"
                ),
                resource_reason="max_rational_bit_length",
                lambda_tube=lambda_tube,
                k_tube=k_tube,
                initial_inside_tube=True,
            )
        if exc.reason == "domain_error":
            return _failed_tube(
                radius=radius,
                lambda_left=lambda_left,
                k_left=k_left,
                policy=policy,
                budget=budget,
                classification="domain_error",
                obstruction=str(exc.payload.get("obstruction") or "physical_domain_escape"),
                domain=True,
                lambda_tube=lambda_tube,
                k_tube=k_tube,
                initial_inside_tube=True,
            )
        reason = str(exc.payload.get("obstruction") or "zero_in_interval_reciprocal")
        classification = (
            "jacobian_diagonal_contains_zero"
            if reason == "jacobian_diagonal_contains_zero"
            else "interval_inconclusive"
        )
        return _failed_tube(
            radius=radius,
            lambda_left=lambda_left,
            k_left=k_left,
            policy=policy,
            budget=budget,
            classification=classification,
            obstruction=reason,
            lambda_tube=lambda_tube,
            k_tube=k_tube,
            initial_inside_tube=True,
            jacobian_diagonals_exclude_zero=reason != "jacobian_diagonal_contains_zero",
        )
    except ZeroDivisionError:
        return _failed_tube(
            radius=radius,
            lambda_left=lambda_left,
            k_left=k_left,
            policy=policy,
            budget=budget,
            classification="interval_inconclusive",
            obstruction="zero_in_interval_reciprocal",
            lambda_tube=lambda_tube,
            k_tube=k_tube,
            initial_inside_tube=True,
        )


def _failed_tube(
    *,
    radius: Interval,
    lambda_left: Interval,
    k_left: Interval,
    policy: SGBLOrbitTubePolicy,
    budget: _TubeBudget,
    classification: str,
    obstruction: str | None,
    resource_reason: str | None = None,
    domain: bool = False,  # classification already encodes domain vs wrapping
    lambda_tube: Interval | None = None,
    k_tube: Interval | None = None,
    lambda_image: Interval | None = None,
    k_image: Interval | None = None,
    lambda_right: Interval | None = None,
    k_right: Interval | None = None,
    initial_inside_tube: bool = False,
    picard_self_map: bool = False,
    hamiltonian_lambda_r: Interval | None = None,
    momentum_k_r: Interval | None = None,
    jacobian_diagonals_exclude_zero: bool = False,
    lipschitz: Fraction = Q(0),
    max_speed: Fraction = Q(0),
) -> SGBLOrbitTubeEnclosure:
    if lambda_tube is None:
        lambda_tube = sgbl_orbit_tube_around(lambda_left, policy.tube_radius)
    if k_tube is None:
        k_tube = sgbl_orbit_tube_around(k_left, policy.tube_radius)
    if lambda_image is None:
        lambda_image = lambda_left
    if k_image is None:
        k_image = k_left
    if lambda_right is None:
        lambda_right = lambda_left
    if k_right is None:
        k_right = k_left
    if hamiltonian_lambda_r is None:
        hamiltonian_lambda_r = ZERO
    if momentum_k_r is None:
        momentum_k_r = ZERO
    compactness = Interval.singleton(2)
    if radius.upper <= radius.lower:
        radius = Interval(radius.lower, radius.lower + policy.step_width)
    return SGBLOrbitTubeEnclosure(
        radius=radius,
        lambda_left=lambda_left,
        k_left=k_left,
        lambda_tube=lambda_tube,
        k_tube=k_tube,
        lambda_image=lambda_image,
        k_image=k_image,
        lambda_right=lambda_right,
        k_right=k_right,
        compactness=compactness,
        direct_compactness=compactness,
        invariant_residual=ONE,
        hamiltonian_lambda_r=hamiltonian_lambda_r,
        momentum_k_r=momentum_k_r,
        lipschitz=lipschitz,
        max_speed=max_speed,
        contraction_constant=lipschitz * max(radius.upper - radius.lower, Q(0)),
        initial_inside_tube=initial_inside_tube,
        picard_self_map=picard_self_map,
        jacobian_diagonals_exclude_zero=jacobian_diagonals_exclude_zero,
        residual_contains_origin=False,
        compactness_invariant_contains_zero=False,
        unique_local_affine_orbit=False,
        classification=classification,
        obstruction=obstruction if classification != "resource_limit" else resource_reason,
        missing_theorem=_missing_theorem(
            resource_reason if classification == "resource_limit" else obstruction
        ),
        rhs_evaluations=budget.rhs_evaluations,
        compactness_rhs_evaluations=budget.compactness_rhs_evaluations,
        max_rational_bit_length_observed=budget.max_bits,
        depth=0,
    )


def sgbl_picard_lindelof_constant_ode_control(
    *,
    policy: SGBLOrbitTubePolicy | None = None,
) -> SGBLOrbitTubeEnclosure:
    """Exact control: ``y'=0`` has a unique constant orbit in the frozen tube."""

    policy = policy or DECLARED_ORBIT_TUBE_POLICY
    if type(policy) is not SGBLOrbitTubePolicy:
        raise TypeError("policy must be SGBLOrbitTubePolicy")
    r0 = PREFIX_RADIUS
    r1 = r0 + policy.step_width
    lambda_left = ONE
    k_left = ZERO
    lambda_tube = sgbl_orbit_tube_around(lambda_left, policy.tube_radius)
    k_tube = sgbl_orbit_tube_around(k_left, policy.tube_radius)
    compactness = sgbl_misner_sharp_compactness_interval(
        interval(r0, r1), lambda_left, k_left
    )
    invariant = sgbl_misner_sharp_compactness_invariant_residual(
        interval(r0, r1), lambda_left, k_left, compactness
    )
    if not lambda_left.strictly_inside(lambda_tube) or not compactness.upper < 1:
        raise SGBLOrbitLocalStop(
            "tampered_contract",
            "constant-ODE control must remain a unique untrapped orbit",
        )
    return SGBLOrbitTubeEnclosure(
        radius=interval(r0, r1),
        lambda_left=lambda_left,
        k_left=k_left,
        lambda_tube=lambda_tube,
        k_tube=k_tube,
        lambda_image=lambda_left,
        k_image=k_left,
        lambda_right=lambda_left,
        k_right=k_left,
        compactness=compactness,
        direct_compactness=compactness,
        invariant_residual=invariant,
        hamiltonian_lambda_r=ONE,
        momentum_k_r=ONE,
        lipschitz=Q(0),
        max_speed=Q(0),
        contraction_constant=Q(0),
        initial_inside_tube=True,
        picard_self_map=True,
        jacobian_diagonals_exclude_zero=True,
        residual_contains_origin=True,
        compactness_invariant_contains_zero=invariant.contains_zero(),
        unique_local_affine_orbit=True,
        classification="unique_local_affine_orbit",
        obstruction=None,
        missing_theorem=None,
        rhs_evaluations=0,
        compactness_rhs_evaluations=0,
        max_rational_bit_length_observed=_interval_bits(lambda_tube),
        depth=0,
    )


def sgbl_orbit_local_minkowski_buffer_control(
    spec: SGBLExactInitialSlice | None = None,
    *,
    policy: SGBLOrbitTubePolicy | None = None,
) -> SGBLOrbitTubeEnclosure:
    """Inherited Minkowski-buffer control: exact ``(lambda_r,k_r)=(0,0)``."""

    spec = spec or SGBLExactInitialSlice()
    if type(spec) is not SGBLExactInitialSlice:
        raise TypeError("spec must be SGBLExactInitialSlice")
    policy = policy or DECLARED_ORBIT_TUBE_POLICY
    if type(policy) is not SGBLOrbitTubePolicy:
        raise TypeError("policy must be SGBLOrbitTubePolicy")
    reduction = sgbl_exact_minkowski_rhs_reduction(spec)
    if reduction["lambda_r"] != 0 or reduction["k_r"] != 0 or reduction["compactness"] != 0:
        raise SGBLOrbitLocalStop(
            "tampered_contract",
            "Minkowski-buffer control requires vanishing affine RHS and C=0",
        )
    r0 = spec.support_minimum
    r1 = r0 + policy.step_width
    control = sgbl_picard_lindelof_constant_ode_control(policy=policy)
    return SGBLOrbitTubeEnclosure(
        radius=interval(r0, r1),
        lambda_left=control.lambda_left,
        k_left=control.k_left,
        lambda_tube=control.lambda_tube,
        k_tube=control.k_tube,
        lambda_image=control.lambda_image,
        k_image=control.k_image,
        lambda_right=control.lambda_right,
        k_right=control.k_right,
        compactness=control.compactness,
        direct_compactness=control.direct_compactness,
        invariant_residual=control.invariant_residual,
        hamiltonian_lambda_r=control.hamiltonian_lambda_r,
        momentum_k_r=control.momentum_k_r,
        lipschitz=control.lipschitz,
        max_speed=control.max_speed,
        contraction_constant=control.contraction_constant,
        initial_inside_tube=control.initial_inside_tube,
        picard_self_map=control.picard_self_map,
        jacobian_diagonals_exclude_zero=control.jacobian_diagonals_exclude_zero,
        residual_contains_origin=control.residual_contains_origin,
        compactness_invariant_contains_zero=control.compactness_invariant_contains_zero,
        unique_local_affine_orbit=control.unique_local_affine_orbit,
        classification=control.classification,
        obstruction=control.obstruction,
        missing_theorem=control.missing_theorem,
        rhs_evaluations=control.rhs_evaluations,
        compactness_rhs_evaluations=control.compactness_rhs_evaluations,
        max_rational_bit_length_observed=control.max_rational_bit_length_observed,
        depth=0,
    )


def _continue_tubes(
    spec: SGBLExactInitialSlice,
    *,
    lambda_left: Interval,
    k_left: Interval,
    policy: SGBLOrbitTubePolicy,
) -> tuple[SGBLOrbitTubeEnclosure, ...]:
    nodes = sgbl_declared_orbit_step_grid(
        origin=policy.prefix_radius,
        terminus=policy.support_maximum,
        steps=policy.continuation_steps,
    )
    tubes: list[SGBLOrbitTubeEnclosure] = []
    current_lambda = lambda_left
    current_k = k_left
    for index in range(policy.continuation_steps):
        budget = _TubeBudget()
        tube = sgbl_orbit_tube_picard_lindelof(
            nodes[index],
            nodes[index + 1],
            current_lambda,
            current_k,
            spec,
            policy=policy,
            budget=budget,
        )
        tubes.append(tube)
        if not tube.unique_local_affine_orbit:
            break
        if tube.classification != "unique_local_affine_orbit":
            break
        current_lambda = tube.lambda_right
        current_k = tube.k_right
    inventory = tuple(tubes)
    sgbl_orbit_inventory_span(inventory, origin=policy.prefix_radius)
    return inventory


def _aggregate_classification(
    tubes: tuple[SGBLOrbitTubeEnclosure, ...],
    *,
    tiles_remaining_support: bool,
) -> tuple[str, str | None]:
    if not tubes:
        return "interval_inconclusive", "prefix_endpoint_not_inside_declared_tube"
    last = tubes[-1]
    if last.classification == "resource_limit":
        return "resource_limit", last.obstruction
    if last.classification == "domain_error":
        return "domain_error", last.obstruction
    if last.classification == "jacobian_diagonal_contains_zero":
        return "jacobian_diagonal_contains_zero", last.obstruction
    if last.classification == "on_orbit_compactness_not_below_one":
        return "on_orbit_compactness_not_below_one", last.obstruction
    if (
        tiles_remaining_support
        and all(tube.classification == "unique_local_affine_orbit" for tube in tubes)
        and all(tube.compactness.upper < 1 for tube in tubes)
    ):
        return "unique_local_affine_orbit", None
    if last.obstruction == "prefix_endpoint_not_inside_declared_tube":
        return "interval_inconclusive", "prefix_endpoint_not_inside_declared_tube"
    if not tiles_remaining_support:
        return "interval_inconclusive", last.obstruction or "incomplete_remaining_support"
    return "interval_inconclusive", last.obstruction or "picard_strict_self_map_failed"


def _contract_payload(
    spec: SGBLExactInitialSlice,
    policy: SGBLOrbitTubePolicy,
    tubes: tuple[SGBLOrbitTubeEnclosure, ...],
    *,
    classification: str,
    obstruction: str | None,
    predecessor_hashes: Mapping[str, str],
    flags: Mapping[str, bool],
    prefix_identities: Mapping[str, Any],
) -> dict[str, Any]:
    coverage_left, coverage_right = sgbl_orbit_inventory_span(
        tubes, origin=policy.prefix_radius
    )
    certified = tuple(tube for tube in tubes if tube.unique_local_affine_orbit)
    compactness_upper = (
        max(tube.compactness.upper for tube in certified) if certified else Q(2)
    )
    return {
        "INSTRUMENT_ID": INSTRUMENT_ID,
        "chi_amplitude": _fraction_text(spec.chi_amplitude),
        "affine_state": AFFINE_STATE,
        "picard_operator": PICARD_OPERATOR,
        "lipschitz_norm": LIPSCHITZ_NORM,
        "prefix_radius": _fraction_text(policy.prefix_radius),
        "support_maximum": _fraction_text(policy.support_maximum),
        "continuation_steps": policy.continuation_steps,
        "step_width": _fraction_text(policy.step_width),
        "tube_radius": _fraction_text(policy.tube_radius),
        "max_bisection_depth": policy.max_bisection_depth,
        "max_rational_bit_length": policy.max_rational_bit_length,
        "classification": classification,
        "obstruction": obstruction,
        "missing_theorem": _missing_theorem(obstruction),
        "tube_count": len(tubes),
        "unique_tubes": sum(1 for tube in tubes if tube.unique_local_affine_orbit),
        "coverage_left": _fraction_text(coverage_left),
        "coverage_right": _fraction_text(coverage_right),
        "tiles_remaining_support": coverage_right == policy.support_maximum
        and bool(tubes)
        and all(tube.classification == "unique_local_affine_orbit" for tube in tubes),
        "compactness_upper": _fraction_text(compactness_upper),
        "sampled_nodes_are_not_the_certificate": True,
        "prefix": dict(prefix_identities),
        "tubes": [
            {
                "radius": list(_interval_text(tube.radius)),
                "classification": tube.classification,
                "obstruction": tube.obstruction,
                "unique_local_affine_orbit": tube.unique_local_affine_orbit,
                "initial_inside_tube": tube.initial_inside_tube,
                "picard_self_map": tube.picard_self_map,
                "jacobian_diagonals_exclude_zero": tube.jacobian_diagonals_exclude_zero,
                "lipschitz": _fraction_text(tube.lipschitz),
                "contraction_constant": _fraction_text(tube.contraction_constant),
                "compactness": list(_interval_text(tube.compactness)),
            }
            for tube in tubes
        ],
        **dict(predecessor_hashes),
        **{name: bool(flags[name]) for name in FALSE_CLAIMS},
    }


def sgbl_orbit_local_contract_sha256(payload: Mapping[str, Any]) -> str:
    """SHA-256 of the compact nonpromoting SOL1 payload."""

    raw = dumps(dict(payload), sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return sha256(raw.encode("ascii")).hexdigest()


@dataclass(frozen=True, slots=True, kw_only=True)
class SGBLOrbitLocalRecord:
    """Frozen-policy orbit-local Picard-Lindelof continuation record."""

    slice: SGBLExactInitialSlice
    policy: SGBLOrbitTubePolicy
    tubes: tuple[SGBLOrbitTubeEnclosure, ...]
    classification: str
    obstruction: str | None
    missing_theorem: str | None
    coverage_left: Fraction
    coverage_right: Fraction
    compactness_upper: Fraction
    compactness_margin: Fraction
    predecessor_trap_refinement_contract_sha256: str
    predecessor_trap_refinement2_contract_sha256: str
    predecessor_trap_taylor_contract_sha256: str
    predecessor_trap_barrier_contract_sha256: str
    contract_payload: Mapping[str, Any]
    contract_sha256: str

    def __post_init__(self) -> None:
        if type(self.slice) is not SGBLExactInitialSlice:
            raise TypeError("slice must be SGBLExactInitialSlice")
        if type(self.policy) is not SGBLOrbitTubePolicy:
            raise TypeError("policy must be SGBLOrbitTubePolicy")
        if self.classification not in TUBE_CLASSIFICATIONS:
            raise ValueError("unknown orbit-local classification")
        left, right = sgbl_validate_orbit_tube_inventory(
            self.tubes, origin=self.policy.prefix_radius
        )
        if self.coverage_left != left or self.coverage_right != right:
            raise SGBLOrbitLocalStop(
                "tampered_contract",
                "coverage does not match the orbit-tube inventory",
            )
        if self.compactness_margin != 1 - self.compactness_upper:
            raise SGBLOrbitLocalStop(
                "tampered_contract",
                "compactness_margin must be the exact C<1 remainder",
            )
        if any(dict(self.contract_payload)[name] is not False for name in FALSE_CLAIMS):
            raise SGBLOrbitLocalStop(
                "claim_mutation",
                "nonpromoting contract payload cannot set aggregate or physics claims true",
            )
        digest = sgbl_orbit_local_contract_sha256(dict(self.contract_payload))
        if self.contract_sha256 != digest:
            raise SGBLOrbitLocalStop(
                "tampered_contract",
                "contract hash does not match the nonpromoting payload",
            )

    @property
    def tiles_remaining_support(self) -> bool:
        return (
            bool(self.tubes)
            and self.coverage_left == self.policy.prefix_radius
            and self.coverage_right == self.policy.support_maximum
            and all(tube.classification == "unique_local_affine_orbit" for tube in self.tubes)
        )

    @property
    def unique_affine_orbit_on_remaining_support(self) -> bool:
        return (
            self.classification == "unique_local_affine_orbit"
            and self.tiles_remaining_support
            and self.compactness_upper < 1
            and self.compactness_margin > 0
        )

    @property
    def sampled_nodes_are_not_the_certificate(self) -> bool:
        return True

    @property
    def SGBL_branch_owned_and_healthy(self) -> bool:
        return False

    @property
    def FRZ1(self) -> bool:
        return False

    @property
    def PREF1(self) -> bool:
        return False

    @property
    def holdout_authorized(self) -> bool:
        return False

    @property
    def execution_authorized(self) -> bool:
        return False

    @property
    def actual_orbit_traps(self) -> bool:
        return False

    @property
    def sgbl_model_rejected(self) -> bool:
        return False

    @property
    def whole_domain_nagumo_invariant(self) -> bool:
        return False


def _build_record(
    spec: SGBLExactInitialSlice,
    policy: SGBLOrbitTubePolicy,
    tubes: tuple[SGBLOrbitTubeEnclosure, ...],
    *,
    predecessor_hashes: Mapping[str, str],
    flags: Mapping[str, bool],
    prefix_identities: Mapping[str, Any],
) -> SGBLOrbitLocalRecord:
    coverage_left, coverage_right = sgbl_orbit_inventory_span(
        tubes, origin=policy.prefix_radius
    )
    tiles = (
        bool(tubes)
        and coverage_right == policy.support_maximum
        and all(tube.classification == "unique_local_affine_orbit" for tube in tubes)
    )
    classification, obstruction = _aggregate_classification(
        tubes, tiles_remaining_support=tiles
    )
    certified = tuple(tube for tube in tubes if tube.unique_local_affine_orbit)
    compactness_upper = max((tube.compactness.upper for tube in certified), default=Q(2))
    payload = _contract_payload(
        spec,
        policy,
        tubes,
        classification=classification,
        obstruction=obstruction,
        predecessor_hashes=predecessor_hashes,
        flags=flags,
        prefix_identities=prefix_identities,
    )
    return SGBLOrbitLocalRecord(
        slice=spec,
        policy=policy,
        tubes=tubes,
        classification=classification,
        obstruction=obstruction,
        missing_theorem=_missing_theorem(obstruction),
        coverage_left=coverage_left,
        coverage_right=coverage_right,
        compactness_upper=compactness_upper,
        compactness_margin=1 - compactness_upper,
        predecessor_trap_refinement_contract_sha256=predecessor_hashes[
            "predecessor_trap_refinement_contract_sha256"
        ],
        predecessor_trap_refinement2_contract_sha256=predecessor_hashes[
            "predecessor_trap_refinement2_contract_sha256"
        ],
        predecessor_trap_taylor_contract_sha256=predecessor_hashes[
            "predecessor_trap_taylor_contract_sha256"
        ],
        predecessor_trap_barrier_contract_sha256=predecessor_hashes[
            "predecessor_trap_barrier_contract_sha256"
        ],
        contract_payload=MappingProxyType(payload),
        contract_sha256=sgbl_orbit_local_contract_sha256(payload),
    )


def sgbl_bind_orbit_local_picard_lindelof(
    spec: SGBLExactInitialSlice | None = None,
    *,
    policy: SGBLOrbitTubePolicy | None = None,
    prefix_ode: SGBLValidatedODERecord | None = None,
    predecessor_trap_refinement_contract_sha256: str | None = None,
    predecessor_trap_refinement2_contract_sha256: str | None = None,
    predecessor_trap_taylor_contract_sha256: str | None = None,
    predecessor_trap_barrier_contract_sha256: str | None = None,
    claims: Mapping[str, Any] | None = None,
) -> SGBLOrbitLocalRecord:
    """Authenticate the ``(C,k)`` prefix and evaluate the frozen orbit-tube policy."""

    spec = _require_nominal_family(spec or SGBLExactInitialSlice())
    frozen = _require_declared_policy(
        DECLARED_ORBIT_TUBE_POLICY if policy is None else policy
    )
    hashes = _require_predecessor_hashes(
        PREDECESSOR_TRAP_REFINEMENT_CONTRACT_SHA256
        if predecessor_trap_refinement_contract_sha256 is None
        else predecessor_trap_refinement_contract_sha256,
        PREDECESSOR_TRAP_REFINEMENT2_CONTRACT_SHA256
        if predecessor_trap_refinement2_contract_sha256 is None
        else predecessor_trap_refinement2_contract_sha256,
        PREDECESSOR_TRAP_TAYLOR_CONTRACT_SHA256
        if predecessor_trap_taylor_contract_sha256 is None
        else predecessor_trap_taylor_contract_sha256,
        PREDECESSOR_TRAP_BARRIER_CONTRACT_SHA256
        if predecessor_trap_barrier_contract_sha256 is None
        else predecessor_trap_barrier_contract_sha256,
    )
    flags = _require_false_claims(claims)
    ode = prefix_ode if prefix_ode is not None else sgbl_validated_constraint_ode(spec)
    last = sgbl_authenticate_ck_prefix(ode, spec=spec)
    tubes = _continue_tubes(
        spec,
        lambda_left=last.lambda_right,
        k_left=last.k_right,
        policy=frozen,
    )
    prefix_identities = {
        "chart_kind": CERTIFICATE_CHART_KIND,
        "coverage_right": _fraction_text(PREFIX_RADIUS),
        "last_cell_left": _fraction_text(PREFIX_LAST_CELL_LEFT),
        "cell_count": PREFIX_CELL_COUNT,
        "obstruction": PREFIX_OBSTRUCTION,
        "c_upper": _fraction_text(PREFIX_C_UPPER),
        "c_margin": _fraction_text(PREFIX_C_MARGIN),
        "d_margin": _fraction_text(PREFIX_D_MARGIN),
        "rhs_evaluations": PREFIX_RHS_EVALUATIONS,
        "compactness_rhs_evaluations": PREFIX_COMPACTNESS_RHS_EVALUATIONS,
    }
    return _build_record(
        spec,
        frozen,
        tubes,
        predecessor_hashes={
            "predecessor_trap_refinement_contract_sha256": hashes[0],
            "predecessor_trap_refinement2_contract_sha256": hashes[1],
            "predecessor_trap_taylor_contract_sha256": hashes[2],
            "predecessor_trap_barrier_contract_sha256": hashes[3],
        },
        flags=flags,
        prefix_identities=prefix_identities,
    )


def sgbl_orbit_local_small_amplitude_control(
    spec: SGBLExactInitialSlice | None = None,
    *,
    policy: SGBLOrbitTubePolicy | None = None,
    prefix_ode: SGBLValidatedODERecord | None = None,
) -> SGBLOrbitLocalRecord:
    """Named ``A_chi=1/8`` control on the same frozen remaining-support tubes.

    This control never replaces the nominal family and never opens aggregate
    health.  It uses the same ``r0``, step count, ``rho``, and bit cap.
    """

    if spec is None:
        spec = SGBLExactInitialSlice(chi_amplitude=SMALL_AMPLITUDE_CONTROL)
    if type(spec) is not SGBLExactInitialSlice:
        raise TypeError("spec must be SGBLExactInitialSlice")
    if spec.chi_amplitude != SMALL_AMPLITUDE_CONTROL:
        raise SGBLOrbitLocalStop(
            "family_mutation",
            "small-amplitude control must remain A_chi=1/8",
            {"chi_amplitude": _fraction_text(spec.chi_amplitude)},
        )
    frozen = _require_declared_policy(
        DECLARED_ORBIT_TUBE_POLICY if policy is None else policy
    )
    hashes = _require_predecessor_hashes(
        PREDECESSOR_TRAP_REFINEMENT_CONTRACT_SHA256,
        PREDECESSOR_TRAP_REFINEMENT2_CONTRACT_SHA256,
        PREDECESSOR_TRAP_TAYLOR_CONTRACT_SHA256,
        PREDECESSOR_TRAP_BARRIER_CONTRACT_SHA256,
    )
    ode = prefix_ode if prefix_ode is not None else sgbl_validated_constraint_ode(spec)
    if ode.spec != spec:
        raise SGBLOrbitLocalStop(
            "family_mutation",
            "small-amplitude prefix ODE does not match A_chi=1/8",
        )
    lambda_left, k_left = _state_at_prefix_radius(ode)
    tubes = _continue_tubes(
        spec,
        lambda_left=lambda_left,
        k_left=k_left,
        policy=frozen,
    )
    return _build_record(
        spec,
        frozen,
        tubes,
        predecessor_hashes={
            "predecessor_trap_refinement_contract_sha256": hashes[0],
            "predecessor_trap_refinement2_contract_sha256": hashes[1],
            "predecessor_trap_taylor_contract_sha256": hashes[2],
            "predecessor_trap_barrier_contract_sha256": hashes[3],
        },
        flags={name: False for name in FALSE_CLAIMS},
        prefix_identities={
            "control": "A_chi=1/8",
            "prefix_radius": _fraction_text(PREFIX_RADIUS),
            "nominal_family_replaced": False,
        },
    )


def sgbl_orbit_local_health_gate(record: SGBLOrbitLocalRecord) -> dict[str, Any]:
    """Aggregate flags stay false even after a local unique-orbit tube."""

    if type(record) is not SGBLOrbitLocalRecord:
        raise TypeError("record must be SGBLOrbitLocalRecord")
    return {
        "classification": record.classification,
        "unique_affine_orbit_on_remaining_support": (
            record.unique_affine_orbit_on_remaining_support
        ),
        "tiles_remaining_support": record.tiles_remaining_support,
        "sampled_nodes_are_not_the_certificate": True,
        "contract_sha256": record.contract_sha256,
        **{name: False for name in FALSE_CLAIMS},
    }


__all__ = [
    "AFFINE_STATE",
    "DECLARED_CONTINUATION_STEPS",
    "DECLARED_ORBIT_TUBE_POLICY",
    "DECLARED_STEP_WIDTH",
    "DECLARED_TUBE_RADIUS",
    "FALSE_CLAIMS",
    "FORBIDDEN_HEALTH_IMPORTS",
    "INSTRUMENT_ID",
    "LIPSCHITZ_NORM",
    "NOMINAL_CHI_AMPLITUDE",
    "PICARD_OPERATOR",
    "PREDECESSOR_TRAP_BARRIER_CONTRACT_SHA256",
    "PREDECESSOR_TRAP_REFINEMENT2_CONTRACT_SHA256",
    "PREDECESSOR_TRAP_REFINEMENT_CONTRACT_SHA256",
    "PREDECESSOR_TRAP_TAYLOR_CONTRACT_SHA256",
    "PREFIX_C_MARGIN",
    "PREFIX_C_UPPER",
    "PREFIX_CELL_COUNT",
    "PREFIX_D_MARGIN",
    "PREFIX_LAST_CELL_LEFT",
    "PREFIX_OBSTRUCTION",
    "PREFIX_RADIUS",
    "SGBLOrbitLocalRecord",
    "SGBLOrbitLocalStop",
    "SGBLOrbitTubeEnclosure",
    "SGBLOrbitTubePolicy",
    "SMALL_AMPLITUDE_CONTROL",
    "sgbl_authenticate_ck_prefix",
    "sgbl_bind_orbit_local_picard_lindelof",
    "sgbl_classify_on_orbit_compactness",
    "sgbl_declared_orbit_step_grid",
    "sgbl_interval_matrix_infinity_norm_2x2",
    "sgbl_orbit_inventory_span",
    "sgbl_orbit_local_contract_sha256",
    "sgbl_orbit_local_health_gate",
    "sgbl_orbit_local_minkowski_buffer_control",
    "sgbl_orbit_local_small_amplitude_control",
    "sgbl_orbit_tube_around",
    "sgbl_orbit_tube_picard_lindelof",
    "sgbl_picard_lindelof_constant_ode_control",
    "sgbl_require_affine_jacobian_diagonals",
    "sgbl_validate_orbit_tube_inventory",
]
