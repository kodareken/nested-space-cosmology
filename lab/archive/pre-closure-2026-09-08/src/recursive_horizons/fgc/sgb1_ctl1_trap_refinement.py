"""Prospective product-box ``(lambda,k,C)`` continuous-no-trap resource ladder.

This owner does not change the ``(C,k)`` certificate chart, the family, the
bump, the physical domains, or the strict ``C<1`` gate.  It freezes a finite
interval-Picard resource ladder and then evaluates
``sgbl_validated_lambda_k_constraint_ode`` at every declared level.

A local pass may say only that one prospectively declared level proves
continuous no-initial-trap on the product-box graph.  Aggregate health,
``FRZ1``, ``PREF1``, execution, and holdout remain false.  Sampled nodes are
not the certificate.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
from json import dumps
from numbers import Integral
from time import perf_counter
from types import MappingProxyType
from typing import Any, Mapping

from .exact_interval import Interval
from .sgb1_ctl1_initial_health import (
    CERTIFICATE_CHART_KIND,
    DECLARED_BASE_CELLS,
    DECLARED_C_DOMAIN,
    DECLARED_J_DOMAIN,
    DECLARED_K_DOMAIN,
    DECLARED_LAMBDA_DOMAIN,
    DECLARED_MAX_PICARD_ITERATIONS,
    DECLARED_MAX_RATIONAL_BITS,
    ODE_INCONCLUSIVE_REASONS,
    SGBLExactInitialSlice,
    SGBLInitialHealthStop,
    SGBLODECellEnclosure,
    SGBLValidatedODEPolicy,
    SGBLValidatedODERecord,
    sgbl_intersect_algebraic_compactness_enclosures,
    sgbl_validate_ode_cell_inventory,
    sgbl_validated_lambda_k_constraint_ode,
)


Q = Fraction
INSTRUMENT_ID = "FGC-1-SGB1-CTL1-TRAP-REFINEMENT"
PRODUCT_BOX_CHART_KIND = "lambda_k"
DECLARED_LADDER_DEPTHS = (8, 10, 12)
DECLARED_COMPACTNESS_STRICT_UPPER = 1
DECLARED_AFFINE_EVALUATIONS_PER_RHS_CALL = 4
DECLARED_RHS_CALLS_PER_PARTITION = 1 + DECLARED_MAX_PICARD_ITERATIONS
NOMINAL_CHI_AMPLITUDE = Q(3)
SMALL_AMPLITUDE_CONTROL = Q(1, 8)
LEVEL_CLASSIFICATIONS = frozenset(
    {
        "continuous_no_initial_trapped_sphere",
        "interval_inconclusive",
        "resource_limit",
        "domain_error",
    }
)
RESOURCE_REASONS = frozenset(
    {
        "max_cells",
        "max_rhs_evaluations",
        "max_rational_bit_length",
    }
)
REFINEMENT_STOP_REASONS = frozenset(
    {
        "changed_ladder",
        "reordered_levels",
        "skipped_level",
        "tampered_contract",
        "gapped_inventory",
    }
)
MISSING_THEOREMS = MappingProxyType(
    {
        "picard_strict_self_map_failed": (
            "strict_interval_Picard_self_map_of_the_affine_constraint_ODE_on_the_declared_refinement_tree"
        ),
        "jacobian_diagonal_contains_zero": (
            "affine_H_L_and_M_k_diagonals_excluding_zero_on_every_declared_cell_box"
        ),
        "physical_domain_escape": (
            "Picard_image_remaining_inside_the_declared_physical_state_domain"
        ),
        "compactness_enclosure_not_below_one": (
            "correlated_interval_compactness_C_strictly_below_one_on_the_validated_augmented_graph"
        ),
        "residual_enclosure_misses_origin": (
            "affine_H_and_M_residual_enclosures_containing_the_origin_on_every_cell"
        ),
        "exterior_compactness_not_below_one": (
            "analytic_exterior_C_equals_two_M_over_r_with_maximum_strictly_below_one"
        ),
        "zero_in_interval_reciprocal": (
            "interval_reciprocals_of_r_and_lambda_remaining_defined_on_every_cell_box"
        ),
        "mean_value_inconsistent": (
            "centered_mean_value_RHS_intersecting_the_natural_interval_extension"
        ),
        "compactness_mean_value_inconsistent": (
            "centered_mean_value_compactness_RHS_intersecting_the_natural_chain_rule_extension"
        ),
        "algebraic_compactness_mean_value_inconsistent": (
            "centered_mean_value_algebraic_C_intersecting_the_natural_interval_extension"
        ),
        "compactness_invariant_misses_origin": (
            "algebraic_and_propagated_compactness_enclosures_overlapping_with_invariant_residual_zero"
        ),
        "chart_denominator_not_strictly_positive": (
            "chart_denominator_D_equals_one_plus_J_squared_over_r_to_the_fourth_minus_C_strictly_positive"
        ),
        "sqrt_branch_not_positive": (
            "positive_square_root_branch_of_D_producing_strictly_positive_lambda"
        ),
        "chart_roundtrip_misses_origin": (
            "Misner_Sharp_chart_roundtrip_residuals_containing_the_origin"
        ),
        "incomplete_support_coverage": (
            "complete_compact_support_covering_by_product_box_Picard_cells"
        ),
    }
)


def _nonnegative_int(name: str, value: object, *, allow_zero: bool = False) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral):
        raise ValueError(f"{name} must be an integer")
    result = int(value)
    if result < 0 or (result == 0 and not allow_zero):
        raise ValueError(f"{name} must be a positive integer")
    return result


def _fraction_text(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def _interval_text(value: Interval) -> tuple[str, str]:
    return (_fraction_text(value.lower), _fraction_text(value.upper))


def _intersect(left: Interval, right: Interval) -> Interval | None:
    lower = max(left.lower, right.lower)
    upper = min(left.upper, right.upper)
    if lower > upper:
        return None
    return Interval(lower, upper)


def sgbl_declared_product_box_cell_cap(max_bisection_depth: int) -> int:
    """Worst-case partition calls on the 16-base-cell binary tree."""

    depth = _nonnegative_int(
        "max_bisection_depth", max_bisection_depth, allow_zero=True
    )
    return DECLARED_BASE_CELLS * (2 ** (depth + 1) - 1)


def sgbl_declared_product_box_rhs_cap(max_bisection_depth: int) -> int:
    """Worst-case affine evaluations: seed plus Picard iterations, four evals each."""

    return (
        sgbl_declared_product_box_cell_cap(max_bisection_depth)
        * DECLARED_RHS_CALLS_PER_PARTITION
        * DECLARED_AFFINE_EVALUATIONS_PER_RHS_CALL
    )


def sgbl_declared_product_box_bit_cap(max_bisection_depth: int) -> int:
    """A priori bit budget.  Depth 8 matches the existing owner; later levels double."""

    depth = _nonnegative_int(
        "max_bisection_depth", max_bisection_depth, allow_zero=True
    )
    extra = max(0, (depth - DECLARED_LADDER_DEPTHS[0]) // 2)
    return DECLARED_MAX_RATIONAL_BITS * (2 ** extra)


class SGBLTrapRefinementStop(ValueError):
    """Typed contract stop of the prospective product-box ladder."""

    def __init__(
        self,
        reason: str,
        message: str,
        payload: Mapping[str, Any] | None = None,
    ) -> None:
        if reason not in REFINEMENT_STOP_REASONS:
            raise ValueError("unknown trap-refinement stop reason")
        self.reason = reason
        self.payload = dict(payload or {})
        super().__init__(message)


@dataclass(frozen=True, slots=True, kw_only=True)
class SGBLProductBoxLadderLevel:
    """One prospectively declared product-box Picard resource level."""

    name: str
    max_bisection_depth: int
    max_cells: int
    max_rhs_evaluations: int
    max_rational_bit_length: int
    base_cells: int = DECLARED_BASE_CELLS
    max_picard_iterations: int = DECLARED_MAX_PICARD_ITERATIONS
    lambda_domain: Interval = DECLARED_LAMBDA_DOMAIN
    k_domain: Interval = DECLARED_K_DOMAIN
    compactness_domain: Interval = DECLARED_C_DOMAIN
    angular_momentum_domain: Interval = DECLARED_J_DOMAIN

    def __post_init__(self) -> None:
        object.__setattr__(self, "name", str(self.name))
        object.__setattr__(
            self,
            "max_bisection_depth",
            _nonnegative_int(
                "max_bisection_depth", self.max_bisection_depth, allow_zero=True
            ),
        )
        object.__setattr__(self, "max_cells", _nonnegative_int("max_cells", self.max_cells))
        object.__setattr__(
            self,
            "max_rhs_evaluations",
            _nonnegative_int("max_rhs_evaluations", self.max_rhs_evaluations),
        )
        object.__setattr__(
            self,
            "max_rational_bit_length",
            _nonnegative_int("max_rational_bit_length", self.max_rational_bit_length),
        )
        object.__setattr__(
            self, "base_cells", _nonnegative_int("base_cells", self.base_cells)
        )
        object.__setattr__(
            self,
            "max_picard_iterations",
            _nonnegative_int("max_picard_iterations", self.max_picard_iterations),
        )
        if self.base_cells != DECLARED_BASE_CELLS:
            raise ValueError("product-box ladder must keep the declared 16 base cells")
        if self.max_picard_iterations != DECLARED_MAX_PICARD_ITERATIONS:
            raise ValueError("product-box ladder must keep the declared Picard iteration cap")
        if type(self.lambda_domain) is not Interval or type(self.k_domain) is not Interval:
            raise TypeError("ladder domains must be exact intervals")
        if type(self.compactness_domain) is not Interval:
            raise TypeError("ladder compactness domain must be an exact interval")
        if (
            self.lambda_domain != DECLARED_LAMBDA_DOMAIN
            or self.k_domain != DECLARED_K_DOMAIN
            or self.compactness_domain != DECLARED_C_DOMAIN
        ):
            raise ValueError("product-box ladder must keep the declared lambda/k/C domains")
        if self.compactness_domain.upper <= DECLARED_COMPACTNESS_STRICT_UPPER:
            raise ValueError("compactness resource domain must not assume C<=1")
        if not (self.compactness_domain.lower < 0 < self.compactness_domain.upper):
            raise ValueError("compactness resource domain must contain a neighborhood of C=0")

    def policy(self) -> SGBLValidatedODEPolicy:
        return SGBLValidatedODEPolicy(
            lambda_domain=self.lambda_domain,
            k_domain=self.k_domain,
            compactness_domain=self.compactness_domain,
            angular_momentum_domain=self.angular_momentum_domain,
            base_cells=self.base_cells,
            max_bisection_depth=self.max_bisection_depth,
            max_cells=self.max_cells,
            max_rhs_evaluations=self.max_rhs_evaluations,
            max_picard_iterations=self.max_picard_iterations,
            max_rational_bit_length=self.max_rational_bit_length,
        )

    def as_contract_mapping(self) -> dict[str, int | str | tuple[str, str]]:
        return {
            "name": self.name,
            "max_bisection_depth": self.max_bisection_depth,
            "max_cells": self.max_cells,
            "max_rhs_evaluations": self.max_rhs_evaluations,
            "max_rational_bit_length": self.max_rational_bit_length,
            "base_cells": self.base_cells,
            "max_picard_iterations": self.max_picard_iterations,
            "lambda_domain": _interval_text(self.lambda_domain),
            "k_domain": _interval_text(self.k_domain),
            "compactness_domain": _interval_text(self.compactness_domain),
        }


def _declare_product_box_resource_ladder() -> tuple[SGBLProductBoxLadderLevel, ...]:
    levels = []
    for depth in DECLARED_LADDER_DEPTHS:
        levels.append(
            SGBLProductBoxLadderLevel(
                name=f"depth_{depth}",
                max_bisection_depth=depth,
                max_cells=sgbl_declared_product_box_cell_cap(depth),
                max_rhs_evaluations=sgbl_declared_product_box_rhs_cap(depth),
                max_rational_bit_length=sgbl_declared_product_box_bit_cap(depth),
            )
        )
    return tuple(levels)


DECLARED_PRODUCT_BOX_RESOURCE_LADDER = _declare_product_box_resource_ladder()


def _require_declared_ladder(
    ladder: object,
) -> tuple[SGBLProductBoxLadderLevel, ...]:
    if type(ladder) is not tuple:
        raise SGBLTrapRefinementStop(
            "changed_ladder",
            "resource ladder must be the frozen declared tuple",
        )
    if not ladder:
        raise SGBLTrapRefinementStop(
            "skipped_level",
            "declared ladder levels were skipped",
        )
    if not all(type(level) is SGBLProductBoxLadderLevel for level in ladder):
        raise SGBLTrapRefinementStop(
            "changed_ladder",
            "resource ladder levels must be SGBLProductBoxLadderLevel",
        )
    depths = tuple(level.max_bisection_depth for level in ladder)
    declared_set = set(DECLARED_LADDER_DEPTHS)
    if len(ladder) < len(DECLARED_LADDER_DEPTHS) and set(depths) <= declared_set:
        raise SGBLTrapRefinementStop(
            "skipped_level",
            "a declared ladder level was skipped",
            {"depths": depths},
        )
    if (
        len(ladder) == len(DECLARED_LADDER_DEPTHS)
        and set(depths) == declared_set
        and depths != DECLARED_LADDER_DEPTHS
    ):
        raise SGBLTrapRefinementStop(
            "reordered_levels",
            "declared ladder levels were reordered",
            {"depths": depths},
        )
    if ladder != DECLARED_PRODUCT_BOX_RESOURCE_LADDER:
        raise SGBLTrapRefinementStop(
            "changed_ladder",
            "resource ladder is not the frozen declaration",
            {"depths": depths},
        )
    return ladder


def sgbl_product_box_inventory_span(
    cells: tuple[SGBLODECellEnclosure, ...],
    spec: SGBLExactInitialSlice,
) -> tuple[Fraction, Fraction]:
    """Return the gap-free inventory span, refusing dropped or gapped cells."""

    if type(spec) is not SGBLExactInitialSlice:
        raise TypeError("spec must be SGBLExactInitialSlice")
    if any(type(cell) is not SGBLODECellEnclosure for cell in cells):
        raise TypeError("cells must be SGBLODECellEnclosure")
    try:
        return sgbl_validate_ode_cell_inventory(cells, origin=spec.support_minimum)
    except ValueError as exc:
        raise SGBLTrapRefinementStop(
            "gapped_inventory",
            "product-box ladder refuses a dropped or gapped cell inventory",
            {"detail": str(exc)},
        ) from exc


def sgbl_product_box_cells_cover_support(
    cells: tuple[SGBLODECellEnclosure, ...],
    spec: SGBLExactInitialSlice,
) -> bool:
    """True only when the validated inventory tiles the whole compact support."""

    if not cells:
        return False
    _left, right = sgbl_product_box_inventory_span(cells, spec)
    return right == spec.support_maximum


def sgbl_product_box_shared_coverage_monotone(
    coarser: tuple[SGBLODECellEnclosure, ...],
    finer: tuple[SGBLODECellEnclosure, ...],
) -> bool:
    """True when overlapping finer product-box enclosures contract inside coarser ones."""

    for left in coarser:
        for right in finer:
            if _intersect(left.radius, right.radius) is None:
                continue
            if (
                _intersect(left.lambda_box, right.lambda_box) is None
                or _intersect(left.k_box, right.k_box) is None
                or _intersect(left.compactness, right.compactness) is None
            ):
                return False
            if right.radius.subset_of(left.radius) and not (
                right.lambda_box.subset_of(left.lambda_box)
                and right.k_box.subset_of(left.k_box)
                and right.compactness.subset_of(left.compactness)
            ):
                return False
    return True


def _resource_reason_from_payload(payload: Mapping[str, Any]) -> str:
    if "evaluations" in payload:
        return "max_rhs_evaluations"
    if "cells" in payload:
        return "max_cells"
    if "observed" in payload:
        return "max_rational_bit_length"
    return "max_rhs_evaluations"


def _support_compactness_upper(ode: SGBLValidatedODERecord | None) -> Fraction:
    if ode is None or not ode.cells:
        return Q(2)
    return ode.support_compactness.upper


def _level_proved(
    *,
    classification: str,
    obstruction: str | None,
    complete_support_coverage: bool,
    compactness_upper: Fraction,
    exterior_compactness_upper: Fraction,
    compactness_margin: Fraction,
    ode: SGBLValidatedODERecord | None,
    resource_reason: str | None,
) -> bool:
    if classification != "continuous_no_initial_trapped_sphere":
        return False
    if (
        obstruction is not None
        or resource_reason is not None
        or not complete_support_coverage
        or compactness_upper >= DECLARED_COMPACTNESS_STRICT_UPPER
        or exterior_compactness_upper >= DECLARED_COMPACTNESS_STRICT_UPPER
        or compactness_margin <= 0
        or ode is None
        or not ode.cells
        or ode.chart_kind != PRODUCT_BOX_CHART_KIND
        or ode.chart_coordinates
        or ode.obstruction is not None
        or not ode.tiles_compact_support
        or not ode.covers_complete_support
    ):
        return False
    return all(
        cell.picard_strict_self_map
        and cell.residual_contains_origin
        and cell.compactness_invariant_contains_zero
        and cell.chart_kind == PRODUCT_BOX_CHART_KIND
        and not cell.chart_coordinates
        and cell.compactness.upper < DECLARED_COMPACTNESS_STRICT_UPPER
        for cell in ode.cells
    )


def _exterior_upper(
    spec: SGBLExactInitialSlice,
    ode: SGBLValidatedODERecord,
) -> tuple[Fraction, str | None]:
    end_radius = Interval.singleton(spec.support_maximum)
    try:
        end_enclosures = sgbl_intersect_algebraic_compactness_enclosures(
            end_radius,
            ode.lambda_end,
            ode.k_end,
            propagated=ode.compactness_end,
        )
    except SGBLInitialHealthStop as exc:
        if exc.reason in {"resource_limit", "domain_error"}:
            raise
        return max(ode.compactness_end.upper, Q(0)), str(
            exc.payload.get("obstruction") or "compactness_invariant_misses_origin"
        )
    return max(end_enclosures["correlated"].upper, Q(0)), None


@dataclass(frozen=True, slots=True, kw_only=True)
class SGBLProductBoxLadderLevelRecord:
    """Outcome of one prospectively declared product-box resource level."""

    level: SGBLProductBoxLadderLevel
    classification: str
    obstruction: str | None
    complete_support_coverage: bool
    compactness_upper: Fraction
    exterior_compactness_upper: Fraction
    compactness_margin: Fraction
    coverage_left: Fraction
    coverage_right: Fraction
    ode_cell_count: int
    rhs_evaluations: int
    compactness_rhs_evaluations: int
    max_rational_bit_length_observed: int
    picard_strict_inclusions: int
    continuous_no_initial_trapped_sphere: bool
    resource_reason: str | None
    missing_theorem: str | None
    diagnostic_wall_seconds: float
    ode: SGBLValidatedODERecord | None = None

    def __post_init__(self) -> None:
        if type(self.level) is not SGBLProductBoxLadderLevel:
            raise TypeError("level must be SGBLProductBoxLadderLevel")
        if self.classification not in LEVEL_CLASSIFICATIONS:
            raise ValueError("unknown product-box ladder classification")
        if type(self.compactness_upper) is not Fraction:
            raise TypeError("compactness_upper must be a Fraction")
        if type(self.exterior_compactness_upper) is not Fraction:
            raise TypeError("exterior_compactness_upper must be a Fraction")
        if type(self.compactness_margin) is not Fraction:
            raise TypeError("compactness_margin must be a Fraction")
        if type(self.coverage_left) is not Fraction or type(self.coverage_right) is not Fraction:
            raise TypeError("coverage endpoints must be Fractions")
        if self.compactness_margin != (
            1 - max(self.compactness_upper, self.exterior_compactness_upper)
        ):
            raise ValueError("compactness_margin must be the exact C<1 remainder")
        if self.resource_reason is not None and self.resource_reason not in RESOURCE_REASONS:
            raise ValueError("unknown resource reason")
        if self.classification == "resource_limit" and self.resource_reason is None:
            raise ValueError("resource_limit requires a typed resource reason")
        if self.classification != "resource_limit" and self.resource_reason is not None:
            raise ValueError("resource reasons cannot label a non-resource classification")
        if self.obstruction is not None and self.obstruction not in ODE_INCONCLUSIVE_REASONS:
            raise ValueError("unknown product-box obstruction")
        proved = _level_proved(
            classification=self.classification,
            obstruction=self.obstruction,
            complete_support_coverage=self.complete_support_coverage,
            compactness_upper=self.compactness_upper,
            exterior_compactness_upper=self.exterior_compactness_upper,
            compactness_margin=self.compactness_margin,
            ode=self.ode,
            resource_reason=self.resource_reason,
        )
        if self.continuous_no_initial_trapped_sphere != proved:
            raise ValueError("continuous-no-trap bit does not match the product-box evidence")
        if proved and self.missing_theorem is not None:
            raise ValueError("a proved product-box level cannot retain a missing theorem")
        if not proved and self.classification == "continuous_no_initial_trapped_sphere":
            raise ValueError("unproved product-box level cannot use the pass classification")
        if self.ode is not None and self.ode.chart_kind != PRODUCT_BOX_CHART_KIND:
            raise ValueError("product-box ladder must keep chart_kind lambda_k")
        if self.ode is not None and self.ode.chart_coordinates:
            raise ValueError("product-box ladder is not the (C,k) certificate chart")
        if self.ode is None:
            if self.complete_support_coverage:
                raise ValueError("a resource or domain stop cannot claim complete support tiling")
            if self.ode_cell_count != 0:
                raise ValueError("missing ODE record cannot report produced cells")
        else:
            left, right = sgbl_validate_ode_cell_inventory(
                self.ode.cells, origin=self.ode.spec.support_minimum
            )
            if self.coverage_left != left or self.coverage_right != right:
                raise ValueError("coverage span does not match the gap-free inventory")
            if self.ode_cell_count != len(self.ode.cells):
                raise ValueError("ode_cell_count does not match the returned inventory")
            tiles = self.ode.tiles_compact_support
            if self.complete_support_coverage != tiles:
                raise ValueError("complete_support_coverage does not match compact-support tiling")
            if self.continuous_no_initial_trapped_sphere and not tiles:
                raise ValueError("a local pass requires complete support tiling")

    @property
    def sampled_nodes_are_not_the_certificate(self) -> bool:
        return True

    @property
    def chart_kind(self) -> str:
        return PRODUCT_BOX_CHART_KIND

    @property
    def product_box_is_not_the_certificate_chart(self) -> bool:
        return True

    def as_contract_mapping(self) -> dict[str, Any]:
        return {
            "name": self.level.name,
            "max_bisection_depth": self.level.max_bisection_depth,
            "classification": self.classification,
            "obstruction": self.obstruction,
            "complete_support_coverage": self.complete_support_coverage,
            "coverage_left": _fraction_text(self.coverage_left),
            "coverage_right": _fraction_text(self.coverage_right),
            "compactness_upper": _fraction_text(self.compactness_upper),
            "exterior_compactness_upper": _fraction_text(self.exterior_compactness_upper),
            "compactness_margin": _fraction_text(self.compactness_margin),
            "ode_cell_count": self.ode_cell_count,
            "rhs_evaluations": self.rhs_evaluations,
            "compactness_rhs_evaluations": self.compactness_rhs_evaluations,
            "max_rational_bit_length_observed": self.max_rational_bit_length_observed,
            "picard_strict_inclusions": self.picard_strict_inclusions,
            "continuous_no_initial_trapped_sphere": self.continuous_no_initial_trapped_sphere,
            "resource_reason": self.resource_reason,
            "missing_theorem": self.missing_theorem,
            "sampled_nodes_are_not_the_certificate": True,
            "chart_kind": PRODUCT_BOX_CHART_KIND,
        }


def sgbl_evaluate_product_box_ladder_level(
    spec: SGBLExactInitialSlice,
    level: SGBLProductBoxLadderLevel,
) -> SGBLProductBoxLadderLevelRecord:
    """Evaluate one declared or attack resource level.  Caps are not raised on failure."""

    if type(spec) is not SGBLExactInitialSlice:
        raise TypeError("spec must be SGBLExactInitialSlice")
    if type(level) is not SGBLProductBoxLadderLevel:
        raise TypeError("level must be SGBLProductBoxLadderLevel")
    started = perf_counter()
    ode: SGBLValidatedODERecord | None = None
    obstruction: str | None = None
    resource_reason: str | None = None
    classification = "interval_inconclusive"
    rhs_evaluations = 0
    compactness_rhs_evaluations = 0
    observed_bits = 0
    try:
        ode = sgbl_validated_lambda_k_constraint_ode(spec, policy=level.policy())
    except SGBLInitialHealthStop as exc:
        if exc.reason == "resource_limit":
            classification = "resource_limit"
            resource_reason = _resource_reason_from_payload(exc.payload)
            rhs_evaluations = int(exc.payload.get("evaluations") or 0)
            observed_bits = int(exc.payload.get("observed") or 0)
        elif exc.reason == "domain_error":
            classification = "domain_error"
        else:
            obstruction = str(
                exc.payload.get("obstruction") or "zero_in_interval_reciprocal"
            )
    wall = perf_counter() - started
    if ode is not None:
        obstruction = ode.obstruction
        rhs_evaluations = ode.rhs_evaluations
        compactness_rhs_evaluations = ode.compactness_rhs_evaluations
        observed_bits = ode.max_rational_bit_length_observed
        if ode.chart_kind != PRODUCT_BOX_CHART_KIND or ode.chart_coordinates:
            raise SGBLTrapRefinementStop(
                "tampered_contract",
                "product-box ladder received a non-lambda_k ODE record",
            )
    cells = ode.cells if ode is not None else ()
    if ode is not None:
        coverage_left, coverage_right = sgbl_product_box_inventory_span(cells, spec)
        complete = sgbl_product_box_cells_cover_support(cells, spec)
        if complete != ode.tiles_compact_support:
            raise SGBLTrapRefinementStop(
                "tampered_contract",
                "product-box tiling bit does not match the validated inventory",
            )
    else:
        coverage_left = spec.support_minimum
        coverage_right = spec.support_minimum
        complete = False
    compactness_upper = _support_compactness_upper(ode)
    exterior_upper = compactness_upper
    if (
        complete
        and ode is not None
        and obstruction is None
        and resource_reason is None
        and classification not in {"resource_limit", "domain_error"}
    ):
        exterior_upper, exterior_reason = _exterior_upper(spec, ode)
        if exterior_reason is not None:
            obstruction = exterior_reason
        elif exterior_upper >= DECLARED_COMPACTNESS_STRICT_UPPER:
            obstruction = "exterior_compactness_not_below_one"
    margin = 1 - max(compactness_upper, exterior_upper)
    if classification not in {"resource_limit", "domain_error"}:
        classification = "interval_inconclusive"
    proved = _level_proved(
        classification="continuous_no_initial_trapped_sphere",
        obstruction=obstruction,
        complete_support_coverage=complete,
        compactness_upper=compactness_upper,
        exterior_compactness_upper=exterior_upper,
        compactness_margin=margin,
        ode=ode,
        resource_reason=resource_reason,
    )
    if proved:
        classification = "continuous_no_initial_trapped_sphere"
    if classification == "resource_limit":
        missing = "declared_product_box_resource_budget"
    elif classification == "domain_error":
        missing = "declared_product_box_physical_domain"
    elif obstruction is not None:
        missing = MISSING_THEOREMS.get(obstruction, obstruction)
    elif not complete:
        missing = MISSING_THEOREMS["incomplete_support_coverage"]
    else:
        missing = None
    return SGBLProductBoxLadderLevelRecord(
        level=level,
        classification=classification,
        obstruction=obstruction,
        complete_support_coverage=complete,
        compactness_upper=compactness_upper,
        exterior_compactness_upper=exterior_upper,
        compactness_margin=margin,
        coverage_left=coverage_left,
        coverage_right=coverage_right,
        ode_cell_count=len(cells),
        rhs_evaluations=rhs_evaluations,
        compactness_rhs_evaluations=compactness_rhs_evaluations,
        max_rational_bit_length_observed=observed_bits,
        picard_strict_inclusions=sum(1 for cell in cells if cell.picard_strict_self_map),
        continuous_no_initial_trapped_sphere=proved,
        resource_reason=resource_reason,
        missing_theorem=missing,
        diagnostic_wall_seconds=wall,
        ode=ode,
    )


def _shared_coverage_monotone(
    levels: tuple[SGBLProductBoxLadderLevelRecord, ...],
) -> bool:
    for earlier, later in zip(levels, levels[1:]):
        earlier_cells = earlier.ode.cells if earlier.ode is not None else ()
        later_cells = later.ode.cells if later.ode is not None else ()
        if not sgbl_product_box_shared_coverage_monotone(earlier_cells, later_cells):
            return False
    return True


def _compactness_upper_nonincreasing(
    levels: tuple[SGBLProductBoxLadderLevelRecord, ...],
) -> bool:
    previous: Fraction | None = None
    for record in levels:
        if record.ode is None or not record.ode.cells:
            continue
        if previous is not None and record.compactness_upper > previous:
            return False
        previous = record.compactness_upper
    return True


def _contract_payload(
    spec: SGBLExactInitialSlice,
    ladder: tuple[SGBLProductBoxLadderLevel, ...],
    levels: tuple[SGBLProductBoxLadderLevelRecord, ...],
    *,
    shared_coverage_monotone: bool,
    compactness_upper_nonincreasing: bool,
) -> dict[str, Any]:
    proved_names = tuple(
        record.level.name
        for record in levels
        if record.continuous_no_initial_trapped_sphere
    )
    return {
        "INSTRUMENT_ID": INSTRUMENT_ID,
        "chart_kind": PRODUCT_BOX_CHART_KIND,
        "certificate_chart_kind": CERTIFICATE_CHART_KIND,
        "product_box_is_not_the_certificate_chart": True,
        "chi_amplitude": _fraction_text(spec.chi_amplitude),
        "phi_amplitude": _fraction_text(spec.phi_amplitude),
        "center": _fraction_text(spec.center),
        "half_width": _fraction_text(spec.half_width),
        "compactness_gate": "strict_C_lt_1",
        "compactness_strict_upper": _fraction_text(Q(DECLARED_COMPACTNESS_STRICT_UPPER)),
        "base_cells": DECLARED_BASE_CELLS,
        "declared_depths": list(DECLARED_LADDER_DEPTHS),
        "ladder": [level.as_contract_mapping() for level in ladder],
        "levels": [record.as_contract_mapping() for record in levels],
        "evaluated_every_declared_level": len(levels) == len(DECLARED_LADDER_DEPTHS),
        "shared_coverage_monotone": shared_coverage_monotone,
        "compactness_upper_nonincreasing": compactness_upper_nonincreasing,
        "proved_level_names": list(proved_names),
        "any_declared_level_proves_continuous_no_initial_trap": bool(proved_names),
        "SGBL_branch_owned_and_healthy": False,
        "FRZ1": False,
        "PREF1": False,
        "holdout_authorized": False,
        "execution_authorized": False,
        "sampled_nodes_are_not_the_certificate": True,
    }


def sgbl_trap_refinement_contract_sha256(payload: Mapping[str, Any]) -> str:
    raw = dumps(dict(payload), sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return sha256(raw.encode("ascii")).hexdigest()


@dataclass(frozen=True, slots=True, kw_only=True)
class SGBLProductBoxLadderRecord:
    """Complete evaluation of the frozen product-box resource ladder.

    A true ``any_declared_level_proves_continuous_no_initial_trap`` bit is only
    a local product-box fact.  It is not ``SGBL_branch_owned_and_healthy``,
    ``FRZ1``, ``PREF1``, or holdout.
    """

    slice: SGBLExactInitialSlice
    ladder: tuple[SGBLProductBoxLadderLevel, ...]
    levels: tuple[SGBLProductBoxLadderLevelRecord, ...]
    shared_coverage_monotone: bool
    compactness_upper_nonincreasing: bool
    proved_level_names: tuple[str, ...]
    any_declared_level_proves_continuous_no_initial_trap: bool
    contract_payload: Mapping[str, Any]
    contract_sha256: str

    def __post_init__(self) -> None:
        if type(self.slice) is not SGBLExactInitialSlice:
            raise TypeError("slice must be SGBLExactInitialSlice")
        if self.ladder != DECLARED_PRODUCT_BOX_RESOURCE_LADDER:
            raise SGBLTrapRefinementStop(
                "changed_ladder",
                "aggregate record must keep the frozen resource ladder",
            )
        if len(self.levels) != len(DECLARED_LADDER_DEPTHS):
            raise SGBLTrapRefinementStop(
                "skipped_level",
                "aggregate record must evaluate every declared level",
            )
        if tuple(record.level for record in self.levels) != DECLARED_PRODUCT_BOX_RESOURCE_LADDER:
            raise SGBLTrapRefinementStop(
                "reordered_levels",
                "aggregate record levels must follow the frozen ladder order",
            )
        expected_monotone = _shared_coverage_monotone(self.levels)
        expected_nonincreasing = _compactness_upper_nonincreasing(self.levels)
        if self.shared_coverage_monotone != expected_monotone:
            raise SGBLTrapRefinementStop(
                "tampered_contract",
                "shared-coverage monotone bit does not match the returned cells",
            )
        if self.compactness_upper_nonincreasing != expected_nonincreasing:
            raise SGBLTrapRefinementStop(
                "tampered_contract",
                "compactness-upper monotone bit does not match the returned levels",
            )
        proved = tuple(
            record.level.name
            for record in self.levels
            if record.continuous_no_initial_trapped_sphere
        )
        if self.proved_level_names != proved:
            raise SGBLTrapRefinementStop(
                "tampered_contract",
                "proved level names do not match the level records",
            )
        if self.any_declared_level_proves_continuous_no_initial_trap != bool(proved):
            raise SGBLTrapRefinementStop(
                "tampered_contract",
                "aggregate pass bit does not match proved levels",
            )
        if proved and not self.shared_coverage_monotone:
            raise SGBLTrapRefinementStop(
                "tampered_contract",
                "a product-box pass cannot be claimed on nonmonotone shared coverage",
            )
        expected_payload = _contract_payload(
            self.slice,
            self.ladder,
            self.levels,
            shared_coverage_monotone=self.shared_coverage_monotone,
            compactness_upper_nonincreasing=self.compactness_upper_nonincreasing,
        )
        if dict(self.contract_payload) != expected_payload:
            raise SGBLTrapRefinementStop(
                "tampered_contract",
                "contract payload does not match the frozen ladder record",
            )
        expected_hash = sgbl_trap_refinement_contract_sha256(expected_payload)
        if self.contract_sha256 != expected_hash:
            raise SGBLTrapRefinementStop(
                "tampered_contract",
                "contract hash does not match the nonpromoting payload",
            )
        if any(
            dict(self.contract_payload)[name] is not False
            for name in (
                "SGBL_branch_owned_and_healthy",
                "FRZ1",
                "PREF1",
                "holdout_authorized",
                "execution_authorized",
            )
        ):
            raise SGBLTrapRefinementStop(
                "tampered_contract",
                "nonpromoting contract payload cannot set aggregate health true",
            )

    @property
    def sampled_nodes_are_not_the_certificate(self) -> bool:
        return True

    @property
    def product_box_is_not_the_certificate_chart(self) -> bool:
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
    def evaluated_every_declared_level(self) -> bool:
        return len(self.levels) == len(DECLARED_LADDER_DEPTHS)

    @property
    def diagnostic_wall_seconds(self) -> tuple[float, ...]:
        return tuple(record.diagnostic_wall_seconds for record in self.levels)


def sgbl_evaluate_product_box_resource_ladder(
    spec: SGBLExactInitialSlice | None = None,
    *,
    ladder: tuple[SGBLProductBoxLadderLevel, ...] | None = None,
) -> SGBLProductBoxLadderRecord:
    """Evaluate every frozen level.  A pass does not stop the remaining levels."""

    if spec is None:
        spec = SGBLExactInitialSlice()
    elif type(spec) is not SGBLExactInitialSlice:
        raise TypeError("spec must be SGBLExactInitialSlice")
    if spec.chi_amplitude != NOMINAL_CHI_AMPLITUDE and spec.chi_amplitude != SMALL_AMPLITUDE_CONTROL:
        raise ValueError(
            "product-box ladder slice must be nominal A_chi=3 or the named small-amplitude control"
        )
    frozen = _require_declared_ladder(
        DECLARED_PRODUCT_BOX_RESOURCE_LADDER if ladder is None else ladder
    )
    records = []
    for level in frozen:
        records.append(sgbl_evaluate_product_box_ladder_level(spec, level))
    levels = tuple(records)
    monotone = _shared_coverage_monotone(levels)
    nonincreasing = _compactness_upper_nonincreasing(levels)
    proved = tuple(
        record.level.name for record in levels if record.continuous_no_initial_trapped_sphere
    )
    payload = _contract_payload(
        spec,
        frozen,
        levels,
        shared_coverage_monotone=monotone,
        compactness_upper_nonincreasing=nonincreasing,
    )
    return SGBLProductBoxLadderRecord(
        slice=spec,
        ladder=frozen,
        levels=levels,
        shared_coverage_monotone=monotone,
        compactness_upper_nonincreasing=nonincreasing,
        proved_level_names=proved,
        any_declared_level_proves_continuous_no_initial_trap=bool(proved),
        contract_payload=MappingProxyType(payload),
        contract_sha256=sgbl_trap_refinement_contract_sha256(payload),
    )


def sgbl_trap_refinement_health_gate(
    record: SGBLProductBoxLadderRecord,
) -> dict[str, Any]:
    """Aggregate flags stay false even if a declared product-box level locally passes."""

    if type(record) is not SGBLProductBoxLadderRecord:
        raise TypeError("record must be SGBLProductBoxLadderRecord")
    return {
        "any_declared_level_proves_continuous_no_initial_trap": (
            record.any_declared_level_proves_continuous_no_initial_trap
        ),
        "proved_level_names": record.proved_level_names,
        "evaluated_every_declared_level": record.evaluated_every_declared_level,
        "shared_coverage_monotone": record.shared_coverage_monotone,
        "sampled_nodes_are_not_the_certificate": True,
        "product_box_is_not_the_certificate_chart": True,
        "SGBL_branch_owned_and_healthy": False,
        "FRZ1": False,
        "PREF1": False,
        "holdout_authorized": False,
        "execution_authorized": False,
        "copied_gr0_or_fgcqr_health_evidence": False,
    }


__all__ = [
    "DECLARED_AFFINE_EVALUATIONS_PER_RHS_CALL",
    "DECLARED_COMPACTNESS_STRICT_UPPER",
    "DECLARED_LADDER_DEPTHS",
    "DECLARED_PRODUCT_BOX_RESOURCE_LADDER",
    "DECLARED_RHS_CALLS_PER_PARTITION",
    "INSTRUMENT_ID",
    "NOMINAL_CHI_AMPLITUDE",
    "PRODUCT_BOX_CHART_KIND",
    "SMALL_AMPLITUDE_CONTROL",
    "SGBLProductBoxLadderLevel",
    "SGBLProductBoxLadderLevelRecord",
    "SGBLProductBoxLadderRecord",
    "SGBLTrapRefinementStop",
    "sgbl_declared_product_box_bit_cap",
    "sgbl_declared_product_box_cell_cap",
    "sgbl_declared_product_box_rhs_cap",
    "sgbl_evaluate_product_box_ladder_level",
    "sgbl_evaluate_product_box_resource_ladder",
    "sgbl_product_box_cells_cover_support",
    "sgbl_product_box_inventory_span",
    "sgbl_product_box_shared_coverage_monotone",
    "sgbl_trap_refinement_contract_sha256",
    "sgbl_trap_refinement_health_gate",
]
