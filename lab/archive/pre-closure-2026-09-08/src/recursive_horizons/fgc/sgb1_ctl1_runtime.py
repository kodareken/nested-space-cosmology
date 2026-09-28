"""Branch-owned SGB-L source and synthetic stage runtime.

The linear branch is ``f(phi)=alpha_gb*phi`` with ``beta=eta=0`` and
independent chi.  The physical residual is the same complete six-row MHG
map already owned by :mod:`sgb1_ctl1_source`.  This module adds a
binary64/grid adapter, method-owned spatial reconstruction of the mixed
and second-radial jets, and a fail-closed member transaction.

``EvolutionState`` stores only ``(u,p,q)=(value,dt,dr)``.  The source
algebra needs the complete lower two-jet
``(u,p,q,p_r,q_r)=(value,dt,dr,dtr,drr)``.  Those last two slots are not
filled with zeros and are not borrowed from an FGC-QR reference
projection.  They are either supplied as declared arrays or reconstructed
by the bound method-owned SBP operator acting on ``p`` and ``q``.  An
``EvolutionState`` without that completion is a typed missing-interface
stop, not a runnable smaller source.

RK4/SSPRK3 here are the numerical-engine *time* methods.  The family's
identically named radial constraint integrators are a different owner.
This slice is synthetic: small annular grids and point batches only.  It
does not open ``SGBL_branch_owned_and_healthy`` or authorize physical
execution.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from fractions import Fraction
from math import isfinite
from numbers import Integral, Real
from types import MappingProxyType
from typing import Any, Mapping, Sequence

import numpy as np

from .evolution.floating_jet import FloatJet2, float_jet, scalar_primal
from .exact_linear_algebra import inverse
from .evolution.numerical_engine import (
    COMPARATOR_METHOD,
    METHODS,
    PRIMARY_METHOD,
    EvolutionRHS,
    EvolutionState,
    SBPFirstDerivative,
    UniformRadialGrid,
    accept_step,
    array_content_sha256,
    propose_step,
)
from .modified_harmonic_reference import modified_harmonic_full_residuals
from .reference_connection import flat_spherical_annulus_reference
from .sgb1_ctl1_family import (
    CONSTRAINT_INTEGRATOR_NAMES,
    sgbl_branch_health_gate,
    sgbl_tensor_constraint_pair,
)
from .sgb1_ctl1_principal import sgbl_principal_point_facts
from .sgb1_ctl1_source import (
    HAT_NORMAL_FACTOR,
    REFERENCE_RADIAL_MINIMUM,
    SGBLSourceInputs,
    SGBLSourceJacobianSolveStop,
    SOURCE_ACCELERATION_ORDER,
    SOURCE_EQUATION_ORDER,
    SOURCE_FIELD_ORDER,
    TILDE_NORMAL_FACTOR,
    sgbl_source_coefficients,
    sgbl_source_residual,
    sgbl_source_solve,
    sgbl_source_state,
)
from .spherical_reduction import ADM_FIELD_ORDER, SphericalState, _state_from_adm_pg_jets


Q = Fraction
ADM_SOURCE_CHART = "adm_accelerations"
METRIC_DTT_CHART = "metric_dtt"
LOWER_JET_GROUPS = ("u", "p", "q", "p_r", "q_r")
LOWER_JET_SLOTS = ("value", "dt", "dr", "dtr", "drr")
FIELD_COUNT = len(SOURCE_FIELD_ORDER)
SOURCE_INPUT_ATTR = {
    "alpha": "alpha",
    "shift": "shift",
    "lambda": "radial_metric",
    "R": "areal_radius",
    "areal_radius": "areal_radius",
    "phi": "phi",
    "chi": "chi",
}
TIME_METHOD_SPATIAL_ORDER = {
    PRIMARY_METHOD: 4,
    COMPARATOR_METHOD: 2,
}
MAX_SYNTHETIC_GRID_POINTS = 17
_ZERO6 = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
_SOURCE_REFERENCE = flat_spherical_annulus_reference(
    radial_domain_minimum=REFERENCE_RADIAL_MINIMUM,
)
MISSING_LOWER_TWO_JET_INTERFACE = (
    "Complete SGB-L lower two-jets are "
    "(u,p,q,p_r,q_r)=(value,dt,dr,dtr,drr) in ADM field order "
    f"{SOURCE_FIELD_ORDER}. EvolutionState declares only (u,p,q). "
    "p_r and q_r are mixed and second-radial derivatives; they are not "
    "recoverable from (u,p,q) without a declared completion. Supply the "
    "arrays explicitly, or bind a method-owned SBP operator that reconstructs "
    "p_r=D_r p and q_r=D_r q. Zeros are not a completion, and FGC-QR "
    "reference-map projection is not a completion."
)


def _fraction(name: str, value: object) -> Fraction:
    if type(value) not in (int, Fraction):
        raise TypeError(f"{name} must be a Fraction or built-in int")
    return Fraction(value)


def _finite(name: str, value: object, *, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise TypeError(f"{name} must be a finite real scalar")
    result = float(value)
    if not isfinite(result):
        raise ValueError(f"{name} must be finite")
    if positive and result <= 0.0:
        raise ValueError(f"{name} must be positive")
    return result


def _radius_fraction(value: object) -> Fraction:
    if type(value) in (int, Fraction):
        result = Fraction(value)
    else:
        result = Fraction.from_float(_finite("coordinate_radius", value))
    if result <= 0:
        raise ValueError("coordinate_radius must be positive")
    return result


def _six_float(name: str, value: object) -> np.ndarray:
    array = np.asarray(value, dtype=np.float64)
    if array.shape != (FIELD_COUNT,) or not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be a finite length-{FIELD_COUNT} vector")
    frozen = np.ascontiguousarray(array).copy()
    frozen.setflags(write=False)
    return frozen


def _frozen_array(value: object) -> np.ndarray:
    frozen = np.ascontiguousarray(np.asarray(value, dtype=np.float64)).copy()
    frozen.setflags(write=False)
    return frozen


def _grid_array(name: str, value: object, shape: tuple[int, int]) -> np.ndarray:
    array = np.asarray(value, dtype=np.float64)
    if array.shape != shape or not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be a finite array of shape {shape}")
    frozen = np.ascontiguousarray(array).copy()
    frozen.setflags(write=False)
    return frozen


def _residual_vector(objects: Sequence[object]) -> np.ndarray:
    return np.asarray([scalar_primal(item) for item in objects], dtype=np.float64)


def _action_mapping(point: SGBLSourceInputs) -> dict[str, Fraction]:
    return {
        "planck_mass": point.planck_mass,
        "scalar_mass": point.scalar_mass,
        "quartic_coupling": point.quartic_coupling,
        "alpha_gb": point.alpha_gb,
        "beta": Q(0),
        "eta": Q(0),
    }


def _require_action_identity(
    *,
    alpha_gb: Fraction,
    beta: Fraction,
    eta: Fraction,
    chart: str,
    allow_zero_coupling_control: bool = False,
) -> None:
    if chart != ADM_SOURCE_CHART:
        raise SGBLRuntimeStop(
            "source_chart_error",
            "SGB-L runtime source chart is ADM accelerations, not metric dtt",
            {"chart": chart},
        )
    if beta != 0 or eta != 0:
        raise SGBLRuntimeStop(
            "action_identity_error",
            "SGB-L action identity is f=alpha_gb*phi with beta=eta=0",
            {"beta": beta, "eta": eta},
        )
    if alpha_gb == 0 and not allow_zero_coupling_control:
        raise SGBLRuntimeStop(
            "action_identity_error",
            "alpha_gb=0 is only an algebraic coupling-limit control",
        )


def _fixture(point: SGBLSourceInputs) -> dict[str, object]:
    return {
        "model_id": "SGB-L" if point.alpha_gb else "GR-0",
        "action_parameters": _action_mapping(point),
    }


class SGBLRuntimeStop(RuntimeError):
    """Typed source, CFL, health, or missing-interface stop.

    The transaction restores the pre-stage member.  A stop is not branch
    health, defocusing evidence, or a Jacobian-floor certificate.
    """

    def __init__(
        self,
        reason: str,
        message: str,
        payload: Mapping[str, Any] | None = None,
    ) -> None:
        self.reason = reason
        self.payload = dict(payload or {})
        super().__init__(message)


@dataclass(frozen=True, slots=True, kw_only=True)
class SGBLRuntimeAction:
    """Declared linear-branch couplings. Health is not a field."""

    planck_mass: Fraction | int
    scalar_mass: Fraction | int
    quartic_coupling: Fraction | int
    alpha_gb: Fraction | int
    beta: Fraction | int = 0
    eta: Fraction | int = 0
    chart: str = ADM_SOURCE_CHART
    allow_zero_coupling_control: bool = False

    def __post_init__(self) -> None:
        for name in (
            "planck_mass",
            "scalar_mass",
            "quartic_coupling",
            "alpha_gb",
            "beta",
            "eta",
        ):
            object.__setattr__(self, name, _fraction(name, getattr(self, name)))
        for name in ("planck_mass", "scalar_mass", "quartic_coupling"):
            if getattr(self, name) <= 0:
                raise ValueError(f"{name} must be positive")
        _require_action_identity(
            alpha_gb=self.alpha_gb,
            beta=self.beta,
            eta=self.eta,
            chart=self.chart,
            allow_zero_coupling_control=self.allow_zero_coupling_control,
        )

    @property
    def independent_chi(self) -> bool:
        return True

    @property
    def linear_gb_branch(self) -> bool:
        return self.beta == 0 and self.eta == 0 and (
            self.alpha_gb != 0 or self.allow_zero_coupling_control
        )

    def as_mapping(self) -> Mapping[str, Fraction]:
        return MappingProxyType(
            {
                "planck_mass": self.planck_mass,
                "scalar_mass": self.scalar_mass,
                "quartic_coupling": self.quartic_coupling,
                "alpha_gb": self.alpha_gb,
                "beta": self.beta,
                "eta": self.eta,
            }
        )


def sgbl_runtime_action_from_source(point: SGBLSourceInputs) -> SGBLRuntimeAction:
    if type(point) is not SGBLSourceInputs:
        raise TypeError("point must be SGBLSourceInputs")
    return SGBLRuntimeAction(
        planck_mass=point.planck_mass,
        scalar_mass=point.scalar_mass,
        quartic_coupling=point.quartic_coupling,
        alpha_gb=point.alpha_gb,
        allow_zero_coupling_control=point.alpha_gb == 0,
    )


def sgbl_floating_source_state(
    point: SGBLSourceInputs,
    accelerations: Sequence[Real] = _ZERO6,
    *,
    chart: str = ADM_SOURCE_CHART,
) -> SphericalState:
    """Map an exact lower jet and floating ADM accelerations to FloatJet2."""

    if type(point) is not SGBLSourceInputs:
        raise TypeError("point must be SGBLSourceInputs")
    _require_action_identity(
        alpha_gb=point.alpha_gb,
        beta=Q(0),
        eta=Q(0),
        chart=chart,
        allow_zero_coupling_control=point.alpha_gb == 0,
    )
    acceleration = _six_float("accelerations", accelerations)
    jets = {
        "alpha": float_jet(point.alpha, dtt=acceleration[0]),
        "shift": float_jet(point.shift, dtt=acceleration[1]),
        "lambda": float_jet(point.radial_metric, dtt=acceleration[2]),
        "areal_radius": float_jet(point.areal_radius, dtt=acceleration[3]),
        "phi": float_jet(point.phi, dtt=acceleration[4]),
        "chi": float_jet(point.chi, dtt=acceleration[5]),
    }
    return _state_from_adm_pg_jets(_fixture(point), jets)


def sgbl_floating_source_residual(
    point: SGBLSourceInputs,
    accelerations: Sequence[Real] = _ZERO6,
    *,
    chart: str = ADM_SOURCE_CHART,
) -> np.ndarray:
    """Complete six-row MHG residual in binary64 at the given accelerations."""

    state = sgbl_floating_source_state(point, accelerations, chart=chart)
    try:
        residual = modified_harmonic_full_residuals(
            state,
            reference=_SOURCE_REFERENCE,
            coordinate_radius=point.coordinate_radius,
            tilde_normal_factor=TILDE_NORMAL_FACTOR,
            hat_normal_factor=HAT_NORMAL_FACTOR,
        )["full_residual_vector"]
    except (ArithmeticError, TypeError, ValueError, ZeroDivisionError) as exc:
        raise SGBLRuntimeStop(
            "source_chart_domain",
            "floating MHG residual evaluation failed",
            {"error": str(exc)},
        ) from exc
    vector = _residual_vector(residual)
    if vector.shape != (FIELD_COUNT,) or not np.all(np.isfinite(vector)):
        raise SGBLRuntimeStop(
            "nonfinite_residual_or_jacobian",
            "floating MHG residual is nonfinite or has the wrong shape",
        )
    return vector


@dataclass(frozen=True, slots=True, kw_only=True)
class SGBLFloatingSourceCoefficients:
    """Floating affine map ``R(a)=R0+J a`` plus the exact ADM inverse when it exists.

    Invertibility of the stored matrix is not kinetic health.
    """

    constant: np.ndarray
    jacobian: np.ndarray
    inverse: np.ndarray | None
    classification: str
    exact_constant: tuple[Fraction, ...] | None = None
    exact_jacobian: tuple[tuple[Fraction, ...], ...] | None = None
    exact_inverse: tuple[tuple[Fraction, ...], ...] | None = None
    exact_determinant: Fraction | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "constant", _six_float("constant", self.constant))
        jacobian = np.asarray(self.jacobian, dtype=np.float64)
        if jacobian.shape != (FIELD_COUNT, FIELD_COUNT) or not np.all(
            np.isfinite(jacobian)
        ):
            raise ValueError("jacobian must be a finite 6-by-6 array")
        frozen = np.ascontiguousarray(jacobian).copy()
        frozen.setflags(write=False)
        object.__setattr__(self, "jacobian", frozen)
        if self.inverse is not None:
            inverse = np.asarray(self.inverse, dtype=np.float64)
            if inverse.shape != (FIELD_COUNT, FIELD_COUNT) or not np.all(
                np.isfinite(inverse)
            ):
                raise ValueError("inverse must be a finite 6-by-6 array")
            frozen_inverse = np.ascontiguousarray(inverse).copy()
            frozen_inverse.setflags(write=False)
            object.__setattr__(self, "inverse", frozen_inverse)
        if self.classification not in {
            "exact_inverse",
            "floating_affine_solve",
            "singular_source_jacobian",
            "interval_inconclusive",
        }:
            raise ValueError("unknown floating-source classification")

    @property
    def equation_order(self) -> tuple[str, ...]:
        return SOURCE_EQUATION_ORDER

    @property
    def acceleration_order(self) -> tuple[str, ...]:
        return SOURCE_ACCELERATION_ORDER

    @property
    def SGBL_branch_owned_and_healthy(self) -> bool:
        return False

    def evaluate(self, accelerations: Sequence[Real]) -> np.ndarray:
        acceleration = _six_float("accelerations", accelerations)
        return self.constant + self.jacobian @ acceleration


def sgbl_floating_source_coefficients(
    point: SGBLSourceInputs,
    *,
    chart: str = ADM_SOURCE_CHART,
) -> SGBLFloatingSourceCoefficients:
    """Seven-evaluation floating ADM Jacobian of the complete MHG residual."""

    exact = sgbl_source_coefficients(point)
    constant = sgbl_floating_source_residual(point, _ZERO6, chart=chart)
    columns = []
    for index in range(FIELD_COUNT):
        seed = [0.0] * FIELD_COUNT
        seed[index] = 1.0
        columns.append(
            sgbl_floating_source_residual(point, seed, chart=chart) - constant
        )
    jacobian = np.column_stack(columns)
    if exact.jacobian_determinant == 0:
        return SGBLFloatingSourceCoefficients(
            constant=constant,
            jacobian=jacobian,
            inverse=None,
            classification="singular_source_jacobian",
            exact_constant=exact.constant,
            exact_jacobian=exact.jacobian,
            exact_inverse=None,
            exact_determinant=exact.jacobian_determinant,
        )
    try:
        floating_inverse = np.linalg.inv(jacobian)
    except np.linalg.LinAlgError as exc:
        raise SGBLRuntimeStop(
            "interval_inconclusive",
            "floating ADM Jacobian inverse is inconclusive",
            {"error": str(exc)},
        ) from exc
    return SGBLFloatingSourceCoefficients(
        constant=constant,
        jacobian=jacobian,
        inverse=floating_inverse,
        classification="exact_inverse",
        exact_constant=exact.constant,
        exact_jacobian=exact.jacobian,
        exact_inverse=tuple(tuple(row) for row in inverse(exact.jacobian)),
        exact_determinant=exact.jacobian_determinant,
    )


def sgbl_floating_source_solve(
    point: SGBLSourceInputs,
    *,
    chart: str = ADM_SOURCE_CHART,
) -> tuple[np.ndarray, np.ndarray, SGBLFloatingSourceCoefficients]:
    """Solve the floating affine map and re-evaluate the actual residual.

    Exact ``det J=0`` remains ``singular_source_jacobian``.  A finite
    residual after the floating solve is a witness, not a fitted floor and
    not branch health.
    """

    coefficients = sgbl_floating_source_coefficients(point, chart=chart)
    if coefficients.classification == "singular_source_jacobian":
        raise SGBLSourceJacobianSolveStop(
            sgbl_source_coefficients(point)
        )
    if coefficients.inverse is None:
        raise SGBLRuntimeStop(
            "interval_inconclusive",
            "floating ADM source inverse is not available",
        )
    accelerations = -coefficients.inverse @ coefficients.constant
    if not np.all(np.isfinite(accelerations)):
        raise SGBLRuntimeStop(
            "nonfinite_residual_or_jacobian",
            "floating ADM accelerations are nonfinite",
        )
    residual = sgbl_floating_source_residual(point, accelerations, chart=chart)
    return accelerations, residual, coefficients


def sgbl_exact_floating_source_comparison(
    point: SGBLSourceInputs,
) -> dict[str, Any]:
    """Independent exact-versus-floating residual, Jacobian and solve witness."""

    exact_residual = np.asarray(
        [float(entry) for entry in sgbl_source_residual(point)], dtype=np.float64
    )
    exact_solve = np.asarray(
        [float(entry) for entry in sgbl_source_solve(point)], dtype=np.float64
    )
    exact_coefficients = sgbl_source_coefficients(point)
    floating_residual = sgbl_floating_source_residual(point)
    accelerations, solved_residual, coefficients = sgbl_floating_source_solve(point)
    exact_jacobian = np.asarray(
        [[float(entry) for entry in row] for row in exact_coefficients.jacobian],
        dtype=np.float64,
    )
    return {
        "residual_infinity_difference": float(
            np.max(np.abs(floating_residual - exact_residual))
        ),
        "jacobian_infinity_difference": float(
            np.max(np.abs(coefficients.jacobian - exact_jacobian))
        ),
        "solve_infinity_difference": float(
            np.max(np.abs(accelerations - exact_solve))
        ),
        "solved_residual_infinity": float(np.max(np.abs(solved_residual))),
        "exact_determinant": exact_coefficients.jacobian_determinant,
        "classification": coefficients.classification,
        "chart": ADM_SOURCE_CHART,
        "independent_chi": True,
        "SGBL_branch_owned_and_healthy": False,
    }


def sgbl_declared_polynomial_grid(
    point: SGBLSourceInputs,
    grid: UniformRadialGrid,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Extend one exact two-jet as a radial polynomial on a small annulus.

    Field values are quadratic in ``r-r0``; first time and radial derivatives
    are linear.  Mixed and second-radial slots stay the exact point values.
    This is declared data, not a zero-jet fill.
    """

    if type(point) is not SGBLSourceInputs:
        raise TypeError("point must be SGBLSourceInputs")
    if not isinstance(grid, UniformRadialGrid):
        raise TypeError("grid must be UniformRadialGrid")
    _check_synthetic_annulus(grid)
    radii = grid.coordinates
    origin = float(point.coordinate_radius)
    u = np.empty((grid.point_count, FIELD_COUNT), dtype=np.float64)
    p = np.empty_like(u)
    q = np.empty_like(u)
    p_r = np.empty_like(u)
    q_r = np.empty_like(u)
    for column, name in enumerate(SOURCE_FIELD_ORDER):
        jet = getattr(point, SOURCE_INPUT_ATTR[name])
        value = float(jet.value)
        dt = float(jet.dt)
        dr = float(jet.dr)
        dtr = float(jet.dtr)
        drr = float(jet.drr)
        displacement = radii - origin
        u[:, column] = value + dr * displacement + 0.5 * drr * displacement**2
        p[:, column] = dt + dtr * displacement
        q[:, column] = dr + drr * displacement
        p_r[:, column] = dtr
        q_r[:, column] = drr
    return u, p, q, p_r, q_r


def _check_synthetic_annulus(grid: UniformRadialGrid) -> None:
    if grid.point_count > MAX_SYNTHETIC_GRID_POINTS:
        raise SGBLRuntimeStop(
            "resource_limit",
            "synthetic SGB-L grids stay at or below the small-mesh cap",
            {
                "point_count": grid.point_count,
                "maximum": MAX_SYNTHETIC_GRID_POINTS,
            },
        )
    if grid.minimum <= float(REFERENCE_RADIAL_MINIMUM):
        raise SGBLRuntimeStop(
            "source_chart_domain",
            "annular SGB-L source requires r greater than the frozen 1/2 minimum",
            {"minimum": grid.minimum},
        )


def complete_lower_jets_from_evolution_state(
    state: EvolutionState,
    *,
    p_r: object | None = None,
    q_r: object | None = None,
    spatial_operator: SBPFirstDerivative | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Return declared or method-owned ``(p_r,q_r)``; never invent zeros."""

    if not isinstance(state, EvolutionState):
        raise TypeError("state must be EvolutionState")
    if state.shape[1] != FIELD_COUNT:
        raise ValueError("evolution state must use the six ADM fields")
    declared = p_r is not None or q_r is not None
    if declared:
        if p_r is None or q_r is None:
            raise SGBLRuntimeStop(
                "missing_lower_two_jets",
                MISSING_LOWER_TWO_JET_INTERFACE,
                {"supplied": {"p_r": p_r is not None, "q_r": q_r is not None}},
            )
        mixed = _grid_array("p_r", p_r, state.shape)
        radial_second = _grid_array("q_r", q_r, state.shape)
        return mixed, radial_second
    if spatial_operator is None:
        raise SGBLRuntimeStop(
            "missing_lower_two_jets",
            MISSING_LOWER_TWO_JET_INTERFACE,
            {"groups": LOWER_JET_GROUPS, "evolution_state_groups": ("u", "p", "q")},
        )
    if not isinstance(spatial_operator, SBPFirstDerivative):
        raise TypeError("spatial_operator must be SBPFirstDerivative")
    if spatial_operator.grid.point_count != state.shape[0]:
        raise ValueError("spatial operator grid does not match the evolution state")
    mixed = spatial_operator.differentiate(state.p)
    radial_second = spatial_operator.differentiate(state.q)
    return _grid_array("p_r", mixed, state.shape), _grid_array(
        "q_r", radial_second, state.shape
    )


def _floating_state_from_adm_slots(
    *,
    values: Sequence[Real],
    dt: Sequence[Real],
    dr: Sequence[Real],
    dtr: Sequence[Real],
    drr: Sequence[Real],
    dtt: Sequence[Real],
    action: SGBLRuntimeAction,
) -> SphericalState:
    jets = {
        adm_name: FloatJet2(
            _finite(f"{adm_name}.value", values[index]),
            dt=_finite(f"{adm_name}.dt", dt[index]),
            dr=_finite(f"{adm_name}.dr", dr[index]),
            dtt=_finite(f"{adm_name}.dtt", dtt[index]),
            dtr=_finite(f"{adm_name}.dtr", dtr[index]),
            drr=_finite(f"{adm_name}.drr", drr[index]),
        )
        for index, adm_name in enumerate(ADM_FIELD_ORDER)
    }
    return _state_from_adm_pg_jets(
        {
            "model_id": "SGB-L",
            "action_parameters": dict(action.as_mapping()),
        },
        jets,
    )


def sgbl_grid_source_residual(
    *,
    u: object,
    p: object,
    q: object,
    p_r: object,
    q_r: object,
    radii: object,
    action: SGBLRuntimeAction,
    accelerations: object | None = None,
    chart: str = ADM_SOURCE_CHART,
) -> np.ndarray:
    """Pointwise complete-MHG residual on a declared two-jet batch."""

    _require_action_identity(
        alpha_gb=action.alpha_gb,
        beta=action.beta,
        eta=action.eta,
        chart=chart,
        allow_zero_coupling_control=action.allow_zero_coupling_control,
    )
    radius = np.asarray(radii, dtype=np.float64)
    if radius.ndim != 1 or radius.size == 0 or not np.all(np.isfinite(radius)):
        raise ValueError("radii must be a finite one-dimensional array")
    if np.any(radius <= float(REFERENCE_RADIAL_MINIMUM)):
        raise SGBLRuntimeStop(
            "source_chart_domain",
            "annular SGB-L source requires r greater than the frozen 1/2 minimum",
        )
    shape = (radius.size, FIELD_COUNT)
    fields = {
        "u": _grid_array("u", u, shape),
        "p": _grid_array("p", p, shape),
        "q": _grid_array("q", q, shape),
        "p_r": _grid_array("p_r", p_r, shape),
        "q_r": _grid_array("q_r", q_r, shape),
    }
    if accelerations is None:
        dtt = np.zeros(shape, dtype=np.float64)
    else:
        dtt = _grid_array("accelerations", accelerations, shape)
    residual = np.empty(shape, dtype=np.float64)
    for index in range(radius.size):
        try:
            state = _floating_state_from_adm_slots(
                values=fields["u"][index],
                dt=fields["p"][index],
                dr=fields["q"][index],
                dtr=fields["p_r"][index],
                drr=fields["q_r"][index],
                dtt=dtt[index],
                action=action,
            )
            vector = modified_harmonic_full_residuals(
                state,
                reference=_SOURCE_REFERENCE,
                coordinate_radius=_radius_fraction(radius[index]),
                tilde_normal_factor=TILDE_NORMAL_FACTOR,
                hat_normal_factor=HAT_NORMAL_FACTOR,
            )["full_residual_vector"]
        except SGBLRuntimeStop:
            raise
        except (ArithmeticError, TypeError, ValueError, ZeroDivisionError) as exc:
            raise SGBLRuntimeStop(
                "source_chart_domain",
                "grid MHG residual evaluation failed",
                {"index": index, "error": str(exc)},
            ) from exc
        residual[index] = _residual_vector(vector)
    if not np.all(np.isfinite(residual)):
        raise SGBLRuntimeStop(
            "nonfinite_residual_or_jacobian",
            "grid MHG residual is nonfinite",
        )
    return residual


def sgbl_grid_source_solve(
    *,
    u: object,
    p: object,
    q: object,
    p_r: object,
    q_r: object,
    radii: object,
    action: SGBLRuntimeAction,
    chart: str = ADM_SOURCE_CHART,
) -> dict[str, Any]:
    """Affine 6x6 ADM solve at every declared two-jet sample."""

    radius = np.asarray(radii, dtype=np.float64)
    shape = (int(radius.size), FIELD_COUNT)
    constant = sgbl_grid_source_residual(
        u=u,
        p=p,
        q=q,
        p_r=p_r,
        q_r=q_r,
        radii=radii,
        action=action,
        accelerations=None,
        chart=chart,
    )
    jacobian = np.empty((shape[0], FIELD_COUNT, FIELD_COUNT), dtype=np.float64)
    for column in range(FIELD_COUNT):
        seed = np.zeros(shape, dtype=np.float64)
        seed[:, column] = 1.0
        evaluated = sgbl_grid_source_residual(
            u=u,
            p=p,
            q=q,
            p_r=p_r,
            q_r=q_r,
            radii=radii,
            action=action,
            accelerations=seed,
            chart=chart,
        )
        jacobian[:, :, column] = evaluated - constant
    accelerations = np.empty(shape, dtype=np.float64)
    inverse = np.empty_like(jacobian)
    classifications = []
    for index in range(shape[0]):
        matrix = jacobian[index]
        try:
            inverse[index] = np.linalg.inv(matrix)
        except np.linalg.LinAlgError as exc:
            raise SGBLRuntimeStop(
                "source_inversion_inconclusive",
                "floating grid ADM source inversion is inconclusive",
                {"index": index},
            ) from exc
        accelerations[index] = -inverse[index] @ constant[index]
        if not np.all(np.isfinite(accelerations[index])):
            raise SGBLRuntimeStop(
                "nonfinite_residual_or_jacobian",
                "grid ADM accelerations are nonfinite",
                {"index": index},
            )
        classifications.append("floating_affine_candidate_with_residual_witness")
    solved_residual = sgbl_grid_source_residual(
        u=u,
        p=p,
        q=q,
        p_r=p_r,
        q_r=q_r,
        radii=radii,
        action=action,
        accelerations=accelerations,
        chart=chart,
    )
    residual_infinity = np.max(np.abs(solved_residual), axis=1)
    coefficient_scale = np.maximum(
        1.0,
        np.max(np.abs(constant), axis=1)
        + np.max(np.sum(np.abs(jacobian), axis=2), axis=1)
        * np.max(np.abs(accelerations), axis=1),
    )
    scaled_residual = residual_infinity / coefficient_scale
    return {
        "accelerations": _frozen_array(accelerations),
        "constant": _frozen_array(constant),
        "jacobian": _frozen_array(jacobian),
        "inverse": _frozen_array(inverse),
        "solved_residual": _frozen_array(solved_residual),
        "solved_residual_infinity": float(np.max(residual_infinity)),
        "scaled_residual_witness": float(np.max(scaled_residual)),
        "source_solve_admission_qualified": False,
        "classifications": tuple(classifications),
        "chart": ADM_SOURCE_CHART,
        "SGBL_branch_owned_and_healthy": False,
    }


def sgbl_coordinate_speed_upper(state: EvolutionState) -> float:
    """Null-coordinate speed ``|shift|+lapse/lambda``; not a borrowed cone certificate."""

    if not isinstance(state, EvolutionState) or state.shape[1] != FIELD_COUNT:
        raise ValueError("state must be a six-field EvolutionState")
    lapse = state.u[:, 0]
    shift = state.u[:, 1]
    radial_metric = state.u[:, 2]
    if np.any(lapse <= 0.0) or np.any(radial_metric <= 0.0):
        raise SGBLRuntimeStop(
            "health_premise_failed",
            "lapse and radial metric must stay positive",
            {
                "minimum_lapse": float(np.min(lapse)),
                "minimum_radial_metric": float(np.min(radial_metric)),
            },
        )
    speeds = np.abs(shift) + lapse / radial_metric
    return float(np.max(speeds))


@dataclass(frozen=True, slots=True)
class SGBLRuntimeMonitor:
    accepted_stage_count: int = 0
    last_accepted_time: float = 0.0
    last_accepted_stage_index: int = -1
    first_failed_premise: str | None = None
    first_failed_stage_index: int | None = None
    first_failed_time: float | None = None

    def __post_init__(self) -> None:
        for name in ("accepted_stage_count", "last_accepted_stage_index"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, Integral):
                raise TypeError(f"{name} must be an integer")
        _finite("last_accepted_time", self.last_accepted_time)


@dataclass(frozen=True, slots=True)
class SGBLCausalState:
    accepted_time: float = 0.0
    accumulated_distance: float = 0.0
    previous_speed_upper: float = 0.0

    def __post_init__(self) -> None:
        for name in ("accepted_time", "accumulated_distance", "previous_speed_upper"):
            value = _finite(name, getattr(self, name))
            if value < 0.0:
                raise ValueError(f"{name} must be nonnegative")
            object.__setattr__(self, name, value)


@dataclass(frozen=True, slots=True)
class SGBLRuntimeCounters:
    rhs_evaluations: int = 0
    source_solves: int = 0
    accepted_steps: int = 0
    accepted_stages: int = 0
    rollbacks: int = 0

    def __post_init__(self) -> None:
        for name in (
            "rhs_evaluations",
            "source_solves",
            "accepted_steps",
            "accepted_stages",
            "rollbacks",
        ):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, Integral) or int(value) < 0:
                raise ValueError(f"{name} must be a nonnegative integer")
            object.__setattr__(self, name, int(value))


@dataclass(frozen=True, slots=True)
class SGBLTracerState:
    labels: tuple[str, ...] = ()
    positions: np.ndarray | None = None

    def __post_init__(self) -> None:
        labels = tuple(str(label) for label in self.labels)
        object.__setattr__(self, "labels", labels)
        if self.positions is None:
            empty = np.zeros((0,), dtype=np.float64)
            empty.setflags(write=False)
            object.__setattr__(self, "positions", empty)
            return
        array = np.asarray(self.positions, dtype=np.float64)
        if not np.all(np.isfinite(array)):
            raise ValueError("tracer positions must be finite")
        frozen = np.ascontiguousarray(array).copy()
        frozen.setflags(write=False)
        object.__setattr__(self, "positions", frozen)


def sgbl_runtime_health_facts(
    *,
    origin: SGBLSourceInputs | None = None,
    solved_state: SphericalState | None = None,
) -> dict[str, Any]:
    """Constraint, center, principal and cone facts. Health remains closed."""

    gate = dict(sgbl_branch_health_gate())
    gate["source_runtime"] = "implemented_synthetic_not_a_health_certificate"
    payload: dict[str, Any] = {
        "health_gate": MappingProxyType(gate),
        "center": MappingProxyType(
            {
                "evolving_center_after_matter_arrives": "unqualified",
                "initial_empty_minkowski_buffer_is_not_this_runtime": True,
                "source_annulus_minimum_is_not_a_centre_chart": True,
            }
        ),
        "cone_certificate": "unqualified",
        "constraint": None,
        "principal": None,
        "SGBL_branch_owned_and_healthy": False,
        "physical_execution_authorized": False,
    }
    state = solved_state
    if state is None and origin is not None:
        state = sgbl_source_state(origin, sgbl_source_solve(origin))
    if state is not None:
        hamiltonian, momentum = sgbl_tensor_constraint_pair(state)
        payload["constraint"] = MappingProxyType(
            {
                "H": hamiltonian,
                "M": momentum,
                "source": "unredefined RED1 physical projections",
                "not_defining_rhs_zeros": True,
            }
        )
        facts = sgbl_principal_point_facts(state)
        payload["principal"] = facts
        payload["cone_certificate"] = facts.cone_certificate
        payload["SGBL_branch_owned_and_healthy"] = facts.SGBL_branch_owned_and_healthy
    return payload


@dataclass(frozen=True, slots=True, kw_only=True)
class SGBLRuntimeMember:
    """Finite synthetic member. Forbidden health bits are derived and false."""

    time: float
    method: str
    action: SGBLRuntimeAction
    grid: UniformRadialGrid
    state: EvolutionState
    monitor: SGBLRuntimeMonitor = SGBLRuntimeMonitor()
    causal: SGBLCausalState = SGBLCausalState()
    tracers: SGBLTracerState = SGBLTracerState()
    counters: SGBLRuntimeCounters = SGBLRuntimeCounters()
    origin: SGBLSourceInputs | None = None
    declared_p_r: np.ndarray | None = None
    declared_q_r: np.ndarray | None = None

    def __post_init__(self) -> None:
        if self.method not in METHODS:
            raise ValueError("method must be a numerical-engine time method")
        if self.method in CONSTRAINT_INTEGRATOR_NAMES:
            raise ValueError("radial constraint integrator names are not time methods")
        if not isinstance(self.action, SGBLRuntimeAction):
            raise TypeError("action must be SGBLRuntimeAction")
        if not isinstance(self.grid, UniformRadialGrid):
            raise TypeError("grid must be UniformRadialGrid")
        if not isinstance(self.state, EvolutionState):
            raise TypeError("state must be EvolutionState")
        _check_synthetic_annulus(self.grid)
        if self.state.shape != (self.grid.point_count, FIELD_COUNT):
            raise ValueError("state shape must match the synthetic ADM grid")
        object.__setattr__(self, "time", _finite("time", self.time))
        if self.declared_p_r is not None:
            object.__setattr__(
                self,
                "declared_p_r",
                _grid_array("declared_p_r", self.declared_p_r, self.state.shape),
            )
        if self.declared_q_r is not None:
            object.__setattr__(
                self,
                "declared_q_r",
                _grid_array("declared_q_r", self.declared_q_r, self.state.shape),
            )
        if (self.declared_p_r is None) != (self.declared_q_r is None):
            raise SGBLRuntimeStop(
                "missing_lower_two_jets",
                MISSING_LOWER_TWO_JET_INTERFACE,
            )

    @property
    def spatial_order(self) -> int:
        return TIME_METHOD_SPATIAL_ORDER[self.method]

    @property
    def not_a_radial_constraint_integrator(self) -> bool:
        return True

    @property
    def SGBL_branch_owned_and_healthy(self) -> bool:
        return False

    @property
    def physical_execution_authorized(self) -> bool:
        return False

    @property
    def physical_fingerprint(self) -> str:
        positions = (
            self.tracers.positions
            if self.tracers.positions is not None
            else np.zeros((0,), dtype=np.float64)
        )
        extras = []
        if self.declared_p_r is not None:
            extras.append(self.declared_p_r)
        if self.declared_q_r is not None:
            extras.append(self.declared_q_r)
        return array_content_sha256(
            np.asarray([self.time], dtype=np.float64),
            self.state.u,
            self.state.p,
            self.state.q,
            *extras,
            np.asarray(
                [
                    self.causal.accepted_time,
                    self.causal.accumulated_distance,
                    self.causal.previous_speed_upper,
                ],
                dtype=np.float64,
            ),
            positions,
            np.asarray(
                [
                    self.monitor.accepted_stage_count,
                    self.monitor.last_accepted_time,
                    self.monitor.last_accepted_stage_index,
                    self.counters.rhs_evaluations,
                    self.counters.source_solves,
                    self.counters.accepted_steps,
                    self.counters.accepted_stages,
                    self.counters.rollbacks,
                ],
                dtype=np.float64,
            ),
        )


class SGBLEvolutionOperator:
    """Method-owned spatial operator plus complete-MHG ADM source."""

    def __init__(self, member: SGBLRuntimeMember) -> None:
        if not isinstance(member, SGBLRuntimeMember):
            raise TypeError("member must be SGBLRuntimeMember")
        self.grid = member.grid
        self.action = member.action
        self.method = member.method
        self.derivative = SBPFirstDerivative(member.grid, member.spatial_order)

    def __call__(self, _time: float, state: EvolutionState) -> EvolutionRHS:
        if not isinstance(state, EvolutionState) or state.shape != (
            self.grid.point_count,
            FIELD_COUNT,
        ):
            raise ValueError("SGB-L RHS state shape differs")
        mixed, radial_second = complete_lower_jets_from_evolution_state(
            state,
            spatial_operator=self.derivative,
        )
        source = sgbl_grid_source_solve(
            u=state.u,
            p=state.p,
            q=state.q,
            p_r=mixed,
            q_r=radial_second,
            radii=self.grid.coordinates,
            action=self.action,
        )
        speed = sgbl_coordinate_speed_upper(state)
        reduction = state.q - self.derivative.differentiate(state.u)
        return EvolutionRHS(
            state.p,
            source["accelerations"],
            mixed,
            {
                "source_residual_infinity": source["solved_residual_infinity"],
                "scaled_residual_witness": source["scaled_residual_witness"],
                "source_solve_admission_qualified": False,
                "source_classification": source["classifications"],
                "chart": source["chart"],
                "coordinate_speed_upper": speed,
                "reduction_constraint_infinity": float(
                    np.max(np.abs(reduction), initial=0.0)
                ),
                "minimum_lapse": float(np.min(state.u[:, 0])),
                "minimum_radial_metric": float(np.min(state.u[:, 2])),
                "spatial_order": self.derivative.order,
                "time_method": self.method,
                "not_a_radial_constraint_integrator": True,
                "SGBL_branch_owned_and_healthy": False,
            },
        )


@dataclass(frozen=True, slots=True)
class SGBLStageReceipt:
    time: float
    step_index: int
    accepted_stages: int
    source_residual_infinity: float
    scaled_residual_witness: float
    coordinate_speed_upper: float
    physical_fingerprint: str
    method: str
    spatial_order: int

    @property
    def SGBL_branch_owned_and_healthy(self) -> bool:
        return False

    @property
    def source_solve_admission_qualified(self) -> bool:
        return False

    @property
    def synthetic_transition_only(self) -> bool:
        return True


class SGBLRuntime:
    """Capture/restore transaction around one synthetic time step."""

    def __init__(self, member: SGBLRuntimeMember) -> None:
        if not isinstance(member, SGBLRuntimeMember):
            raise TypeError("member must be SGBLRuntimeMember")
        self._member = member
        self.last_stop: SGBLRuntimeStop | None = None

    @property
    def member(self) -> SGBLRuntimeMember:
        return self._member

    def capture(self) -> SGBLRuntimeMember:
        return self._member

    def restore(self, member: SGBLRuntimeMember) -> None:
        if not isinstance(member, SGBLRuntimeMember):
            raise TypeError("member must be SGBLRuntimeMember")
        self._member = member

    def _operator(self) -> SGBLEvolutionOperator:
        return SGBLEvolutionOperator(self._member)

    def attempt_step(self, step_size: Real) -> SGBLStageReceipt:
        original = self._member
        dt = _finite("step_size", step_size, positive=True)
        try:
            speed = sgbl_coordinate_speed_upper(original.state)
            if speed * dt > original.grid.spacing:
                raise SGBLRuntimeStop(
                    "cfl_violation",
                    "synthetic CFL bound would be exceeded",
                    {
                        "coordinate_speed_upper": speed,
                        "step_size": dt,
                        "spacing": original.grid.spacing,
                    },
                )
            operator = self._operator()
            proposal = propose_step(
                method=original.method,
                time=original.time,
                step_size=dt,
                state=original.state,
                rhs=operator,
            )
            endpoint = proposal.stages[-1]
            if endpoint.rhs.diagnostics["minimum_lapse"] <= 0.0:
                raise SGBLRuntimeStop(
                    "health_premise_failed",
                    "candidate lapse is not positive",
                )
            if endpoint.rhs.diagnostics["minimum_radial_metric"] <= 0.0:
                raise SGBLRuntimeStop(
                    "health_premise_failed",
                    "candidate radial metric is not positive",
                )
            residual = float(endpoint.rhs.diagnostics["source_residual_infinity"])
            scaled_residual = float(
                endpoint.rhs.diagnostics["scaled_residual_witness"]
            )
            if not isfinite(residual) or not isfinite(scaled_residual):
                raise SGBLRuntimeStop(
                    "health_premise_failed",
                    "source residual witness is nonfinite",
                )
            accepted = accept_step(
                proposal,
                previous_step_index=original.counters.accepted_steps,
                previous_transaction_serial=original.monitor.last_accepted_stage_index
                + 1
                if original.monitor.last_accepted_stage_index >= 0
                else 0,
            )
            stage_count = len(proposal.stages)
            updated = replace(
                original,
                time=proposal.final_time,
                state=proposal.candidate_state,
                monitor=SGBLRuntimeMonitor(
                    accepted_stage_count=original.monitor.accepted_stage_count
                    + stage_count,
                    last_accepted_time=proposal.final_time,
                    last_accepted_stage_index=(
                        original.monitor.last_accepted_stage_index + stage_count
                    ),
                ),
                causal=SGBLCausalState(
                    accepted_time=proposal.final_time,
                    accumulated_distance=(
                        original.causal.accumulated_distance + speed * dt
                    ),
                    previous_speed_upper=speed,
                ),
                counters=SGBLRuntimeCounters(
                    rhs_evaluations=original.counters.rhs_evaluations + stage_count,
                    source_solves=original.counters.source_solves + stage_count,
                    accepted_steps=original.counters.accepted_steps + 1,
                    accepted_stages=original.counters.accepted_stages + stage_count,
                    rollbacks=original.counters.rollbacks,
                ),
                declared_p_r=None,
                declared_q_r=None,
            )
            self._member = updated
            self.last_stop = None
            return SGBLStageReceipt(
                time=updated.time,
                step_index=accepted.step_index,
                accepted_stages=stage_count,
                source_residual_infinity=residual,
                scaled_residual_witness=scaled_residual,
                coordinate_speed_upper=speed,
                physical_fingerprint=updated.physical_fingerprint,
                method=updated.method,
                spatial_order=updated.spatial_order,
            )
        except SGBLRuntimeStop as stop:
            self._member = original
            self.last_stop = stop
            raise
        except SGBLSourceJacobianSolveStop as stop:
            wrapped = SGBLRuntimeStop(
                "singular_source_jacobian",
                "singular source Jacobian",
                {"reason": stop.reason},
            )
            self._member = original
            self.last_stop = wrapped
            raise wrapped from stop


def sgbl_synthetic_member_from_source_point(
    point: SGBLSourceInputs,
    *,
    method: str,
    grid: UniformRadialGrid | None = None,
) -> SGBLRuntimeMember:
    """Declared polynomial embedding of one exact source point on a small annulus."""

    if method not in METHODS:
        raise ValueError("method must be a numerical-engine time method")
    mesh = grid or UniformRadialGrid(2.0, 3.0, 9)
    u, p, q, p_r, q_r = sgbl_declared_polynomial_grid(point, mesh)
    return SGBLRuntimeMember(
        time=0.0,
        method=method,
        action=sgbl_runtime_action_from_source(point),
        grid=mesh,
        state=EvolutionState(u, p, q),
        origin=point,
        declared_p_r=p_r,
        declared_q_r=q_r,
    )


__all__ = [
    "ADM_SOURCE_CHART",
    "FIELD_COUNT",
    "LOWER_JET_GROUPS",
    "LOWER_JET_SLOTS",
    "MAX_SYNTHETIC_GRID_POINTS",
    "METRIC_DTT_CHART",
    "MISSING_LOWER_TWO_JET_INTERFACE",
    "SGBLCausalState",
    "SGBLEvolutionOperator",
    "SGBLFloatingSourceCoefficients",
    "SGBLRuntime",
    "SGBLRuntimeAction",
    "SGBLRuntimeCounters",
    "SGBLRuntimeMember",
    "SGBLRuntimeMonitor",
    "SGBLRuntimeStop",
    "SGBLStageReceipt",
    "SGBLTracerState",
    "TIME_METHOD_SPATIAL_ORDER",
    "complete_lower_jets_from_evolution_state",
    "sgbl_coordinate_speed_upper",
    "sgbl_declared_polynomial_grid",
    "sgbl_exact_floating_source_comparison",
    "sgbl_floating_source_coefficients",
    "sgbl_floating_source_residual",
    "sgbl_floating_source_solve",
    "sgbl_floating_source_state",
    "sgbl_grid_source_residual",
    "sgbl_grid_source_solve",
    "sgbl_runtime_action_from_source",
    "sgbl_runtime_health_facts",
    "sgbl_synthetic_member_from_source_point",
]
