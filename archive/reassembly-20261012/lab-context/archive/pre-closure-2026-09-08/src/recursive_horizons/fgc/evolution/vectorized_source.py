"""Vectorized grid evaluation of the unchanged ACT1/VAR1 REF1 equations.

The exact tensor implementation is written for scalar two-jets.  Repeating it
point by point would make the declared three-resolution experiment
unnecessarily opaque and slow.  This adapter supplies an array scalar and
array two-jet with the same arithmetic interface, then sends the complete
grid through the same curvature, stress, Gauss--Bonnet, and modified-harmonic
tensor routines.

All six candidate ADM accelerations are seeded in one leading batch axis.  As
the Horndeski-class equations are quasilinear, residual differences recover
the six-by-six kinetic block at every radius.  A batched linear solve obtains
the acceleration update, and an independent second residual evaluation checks
the accepted root.  NUM1 compares this route with the scalar binary64 and
exact-rational backends before it may be used by COL1.

The centre is intentionally excluded from this annular evaluator.  Its
accelerations are supplied by the parity-regular limit adapter whose
convergence is tested separately against CTR1.  No small nonzero radius is
silently substituted for ``r=0``.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import isfinite
from numbers import Real
from typing import Any, Sequence

import numpy as np

from ..action import ActionParameters, ModelID
from ..modified_harmonic_reference import (
    ModifiedHarmonicGaugeData,
    _independent_metric_vector,
    physical_connection_data,
)
from ..spherical_reduction import Jet2, SphericalState, residuals
from .floating_jet import scalar_primal
from .health_monitor import FIELD_ORDER


FIELD_COUNT = len(FIELD_ORDER)
FIELD_INDEX = {name: index for index, name in enumerate(FIELD_ORDER)}
ADM_CENTER_PARITIES = np.asarray((1, -1, 1, -1, 1, 1), dtype=np.int64)


def _coerce_action(value: object) -> ActionParameters:
    if isinstance(value, ActionParameters):
        return value
    required = ("planck_mass", "scalar_mass", "quartic_coupling", "beta", "eta")
    if all(hasattr(value, name) for name in required):
        return ActionParameters(
            model_id=ModelID.FGC_QR,
            planck_mass=float(getattr(value, "planck_mass")),
            scalar_mass=float(getattr(value, "scalar_mass")),
            quartic_coupling=float(getattr(value, "quartic_coupling")),
            pulse_width=1.0,
            ricci_coupling=float(getattr(value, "beta")),
            quadratic_gb_coupling=float(getattr(value, "eta")),
        )
    raise TypeError("action must be ActionParameters or the frozen FGCQR action record")


def _finite_array(name: str, value: object, *, ndim: int | None = None) -> np.ndarray:
    answer = np.asarray(value, dtype=np.float64)
    if ndim is not None and answer.ndim != ndim:
        raise ValueError(f"{name} must have {ndim} dimensions")
    if not np.all(np.isfinite(answer)):
        raise ValueError(f"{name} must contain only finite values")
    return np.ascontiguousarray(answer)


class BatchScalar:
    """A finite ndarray with tensor-code-compatible scalar arithmetic."""

    __array_priority__ = 1000

    def __init__(self, value: object) -> None:
        if isinstance(value, BatchScalar):
            self.data = value.data
            return
        if isinstance(value, bool):
            raise TypeError("batch scalar cannot be boolean")
        data = np.asarray(value, dtype=np.float64)
        if not np.all(np.isfinite(data)):
            raise ValueError("batch scalar must be finite")
        self.data = data

    @staticmethod
    def coerce(value: object) -> "BatchScalar":
        return value if isinstance(value, BatchScalar) else BatchScalar(value)

    def _binary(self, other: object, operation) -> "BatchScalar":
        return BatchScalar(operation(self.data, self.coerce(other).data))

    def __add__(self, other: object) -> "BatchScalar":
        return self._binary(other, np.add)

    __radd__ = __add__

    def __neg__(self) -> "BatchScalar":
        return BatchScalar(-self.data)

    def __sub__(self, other: object) -> "BatchScalar":
        return self._binary(other, np.subtract)

    def __rsub__(self, other: object) -> "BatchScalar":
        return self.coerce(other) - self

    def __mul__(self, other: object) -> "BatchScalar":
        return self._binary(other, np.multiply)

    __rmul__ = __mul__

    def __truediv__(self, other: object) -> "BatchScalar":
        denominator = self.coerce(other).data
        if np.any(denominator == 0.0):
            raise ZeroDivisionError("batch scalar division by zero")
        return BatchScalar(self.data / denominator)

    def __rtruediv__(self, other: object) -> "BatchScalar":
        return self.coerce(other) / self

    def __pow__(self, exponent: object) -> "BatchScalar":
        if isinstance(exponent, bool) or not isinstance(exponent, Real):
            return NotImplemented
        return BatchScalar(np.power(self.data, float(exponent)))

    def __abs__(self) -> "BatchScalar":
        return BatchScalar(np.abs(self.data))

    def __eq__(self, other: object) -> bool:
        try:
            left, right = np.broadcast_arrays(self.data, self.coerce(other).data)
        except (TypeError, ValueError):
            return False
        return bool(np.array_equal(left, right))

    def __ne__(self, other: object) -> bool:
        return not self == other

    # Comparisons are domain predicates, not elementwise numerical outputs.
    # A nonpositive member must fail the corresponding strict-positive test.
    def __le__(self, other: object) -> bool:
        return bool(np.any(self.data <= self.coerce(other).data))

    def __lt__(self, other: object) -> bool:
        return bool(np.all(self.data < self.coerce(other).data))

    def __ge__(self, other: object) -> bool:
        return bool(np.any(self.data >= self.coerce(other).data))

    def __gt__(self, other: object) -> bool:
        return bool(np.all(self.data > self.coerce(other).data))

    def __repr__(self) -> str:
        return f"BatchScalar(shape={self.data.shape})"


@dataclass(frozen=True)
class BatchJet2(Jet2):
    """Array-valued implementation of RED1's public ``Jet2`` protocol."""

    def __post_init__(self) -> None:
        for name in ("value", "dt", "dr", "dtt", "dtr", "drr"):
            object.__setattr__(self, name, BatchScalar(getattr(self, name)))

    @classmethod
    def constant(cls, value: object) -> "BatchJet2":
        return cls(BatchScalar(value))

    @staticmethod
    def coerce(value: object) -> "BatchJet2":
        if isinstance(value, BatchJet2):
            return value
        if isinstance(value, Jet2):
            return BatchJet2(
                value.value, value.dt, value.dr, value.dtt, value.dtr, value.drr
            )
        return BatchJet2.constant(value)

    def __add__(self, other: object) -> "BatchJet2":
        value = self.coerce(other)
        return BatchJet2(
            self.value + value.value,
            self.dt + value.dt,
            self.dr + value.dr,
            self.dtt + value.dtt,
            self.dtr + value.dtr,
            self.drr + value.drr,
        )

    __radd__ = __add__

    def __neg__(self) -> "BatchJet2":
        return BatchJet2(
            -self.value, -self.dt, -self.dr, -self.dtt, -self.dtr, -self.drr
        )

    def __sub__(self, other: object) -> "BatchJet2":
        return self + (-self.coerce(other))

    def __rsub__(self, other: object) -> "BatchJet2":
        return self.coerce(other) - self

    def __mul__(self, other: object) -> "BatchJet2":
        value = self.coerce(other)
        return BatchJet2(
            self.value * value.value,
            self.dt * value.value + self.value * value.dt,
            self.dr * value.value + self.value * value.dr,
            self.dtt * value.value
            + 2 * self.dt * value.dt
            + self.value * value.dtt,
            self.dtr * value.value
            + self.dt * value.dr
            + self.dr * value.dt
            + self.value * value.dtr,
            self.drr * value.value
            + 2 * self.dr * value.dr
            + self.value * value.drr,
        )

    __rmul__ = __mul__

    def compose(self, value: object, first: object, second: object) -> "BatchJet2":
        q0 = BatchScalar(value)
        q1 = BatchScalar(first)
        q2 = BatchScalar(second)
        return BatchJet2(
            q0,
            q1 * self.dt,
            q1 * self.dr,
            q1 * self.dtt + q2 * self.dt * self.dt,
            q1 * self.dtr + q2 * self.dt * self.dr,
            q1 * self.drr + q2 * self.dr * self.dr,
        )

    def reciprocal(self) -> "BatchJet2":
        if np.any(self.value.data == 0.0):
            raise ZeroDivisionError("batch jet reciprocal crossed zero")
        return self.compose(
            1 / self.value,
            -1 / (self.value * self.value),
            2 / (self.value**3),
        )

    def __truediv__(self, other: object) -> "BatchJet2":
        return self * self.coerce(other).reciprocal()

    def __rtruediv__(self, other: object) -> "BatchJet2":
        return self.coerce(other) / self


def _unchecked_state(
    *,
    h_tt: BatchJet2,
    h_tr: BatchJet2,
    h_rr: BatchJet2,
    areal_radius: BatchJet2,
    phi: BatchJet2,
    chi: BatchJet2,
    action: ActionParameters,
) -> SphericalState:
    """Construct after vectorized domain checks without scalar bool coercion."""

    state = object.__new__(SphericalState)
    values = {
        "h_tt": h_tt,
        "h_tr": h_tr,
        "h_rr": h_rr,
        "areal_radius": areal_radius,
        "phi": phi,
        "chi": chi,
        "planck_mass": Fraction.from_float(float(action.planck_mass)),
        "beta": Fraction.from_float(float(action.ricci_coupling)),
        "mu": Fraction.from_float(float(action.scalar_mass)),
        "g4": Fraction.from_float(float(action.quartic_coupling)),
        "alpha": Fraction.from_float(float(action.linear_gb_coupling)),
        "eta": Fraction.from_float(float(action.quadratic_gb_coupling)),
        "branch": action.model_id.value,
    }
    for name, value in values.items():
        object.__setattr__(state, name, value)
    return state


def _batch_state_from_adm(
    u: np.ndarray,
    p: np.ndarray,
    q: np.ndarray,
    dtt: np.ndarray,
    p_r: np.ndarray,
    q_r: np.ndarray,
    *,
    action: ActionParameters,
) -> SphericalState:
    # u,p,q,p_r,q_r have (point,field); dtt has (batch,point,field).
    batch = dtt.shape[0]
    points = u.shape[0]

    def jet(field: str) -> BatchJet2:
        index = FIELD_INDEX[field]
        return BatchJet2(
            np.broadcast_to(u[None, :, index], (batch, points)),
            dt=np.broadcast_to(p[None, :, index], (batch, points)),
            dr=np.broadcast_to(q[None, :, index], (batch, points)),
            dtt=dtt[:, :, index],
            dtr=np.broadcast_to(p_r[None, :, index], (batch, points)),
            drr=np.broadcast_to(q_r[None, :, index], (batch, points)),
        )

    lapse = jet("alpha")
    shift = jet("shift")
    radial = jet("lambda")
    lam2 = radial * radial
    return _unchecked_state(
        h_tt=-(lapse * lapse) + lam2 * shift * shift,
        h_tr=lam2 * shift,
        h_rr=lam2,
        areal_radius=jet("R"),
        phi=jet("phi"),
        chi=jet("chi"),
        action=action,
    )


def _zero2() -> list[list[object]]:
    return [[Fraction(0) for _ in range(4)] for _ in range(4)]


def _zero3() -> list[list[list[object]]]:
    return [[[Fraction(0) for _ in range(4)] for _ in range(4)] for _ in range(4)]


def _zero4() -> list[list[list[list[object]]]]:
    return [
        [[[Fraction(0) for _ in range(4)] for _ in range(4)] for _ in range(4)]
        for _ in range(4)
    ]


def _flat_reference_arrays(radius: BatchScalar) -> tuple[Any, Any]:
    gamma = _zero3()
    derivative = _zero4()
    _, radial, theta, phi = range(4)
    gamma[radial][theta][theta] = -radius
    gamma[radial][phi][phi] = -radius
    gamma[theta][radial][theta] = gamma[theta][theta][radial] = 1 / radius
    gamma[phi][radial][phi] = gamma[phi][phi][radial] = 1 / radius
    derivative[radial][radial][theta][theta] = Fraction(-1)
    derivative[radial][radial][phi][phi] = Fraction(-1)
    derivative[radial][theta][radial][theta] = derivative[radial][theta][theta][radial] = -1 / radius**2
    derivative[radial][phi][radial][phi] = derivative[radial][phi][phi][radial] = -1 / radius**2
    derivative[theta][theta][phi][phi] = Fraction(1)
    derivative[theta][phi][theta][phi] = Fraction(-1)
    derivative[theta][phi][phi][theta] = Fraction(-1)
    freeze = lambda value: tuple(
        tuple(tuple(item for item in lower) for lower in upper) for upper in value
    )
    frozen_gamma = freeze(gamma)
    frozen_derivative = tuple(
        tuple(
            tuple(tuple(item for item in lower) for lower in upper)
            for upper in direction
        )
        for direction in derivative
    )
    return frozen_gamma, frozen_derivative


def _batch_auxiliary_with_derivative(
    physical: Any,
    normal_factor: int,
) -> tuple[Any, Any]:
    if normal_factor <= 1:
        raise ValueError("auxiliary normal factor must exceed one")
    inverse = physical.inverse_metric
    derivative = physical.inverse_metric_derivative
    inverse_tt = BatchScalar.coerce(inverse[0][0])
    if np.any(inverse_tt.data >= 0.0):
        raise ValueError("coordinate-time covector must remain timelike")
    normal_outer = tuple(
        tuple(-inverse[a][0] * inverse[b][0] / inverse[0][0] for b in range(4))
        for a in range(4)
    )
    auxiliary = tuple(
        tuple(
            inverse[a][b] - (normal_factor - 1) * normal_outer[a][b]
            for b in range(4)
        )
        for a in range(4)
    )
    normal_derivative = _zero3()
    for direction in range(4):
        denominator = inverse[0][0]
        denominator_derivative = derivative[direction][0][0]
        for a in range(4):
            for b in range(4):
                numerator = inverse[a][0] * inverse[b][0]
                numerator_derivative = (
                    derivative[direction][a][0] * inverse[b][0]
                    + inverse[a][0] * derivative[direction][b][0]
                )
                normal_derivative[direction][a][b] = (
                    -numerator_derivative / denominator
                    + numerator * denominator_derivative / denominator**2
                )
    auxiliary_derivative = tuple(
        tuple(
            tuple(
                derivative[direction][a][b]
                - (normal_factor - 1) * normal_derivative[direction][a][b]
                for b in range(4)
            )
            for a in range(4)
        )
        for direction in range(4)
    )
    return auxiliary, auxiliary_derivative


def _batch_gauge(
    state: SphericalState,
    radius: BatchScalar,
    *,
    tilde_normal_factor: int,
) -> ModifiedHarmonicGaugeData:
    physical = physical_connection_data(state)
    reference, reference_derivative = _flat_reference_arrays(radius)
    gamma_difference = _zero3()
    gamma_difference_derivative = _zero4()
    for a in range(4):
        for b in range(4):
            for c in range(4):
                gamma_difference[a][b][c] = physical.christoffel[a][b][c] - reference[a][b][c]
                for direction in range(4):
                    gamma_difference_derivative[direction][a][b][c] = (
                        physical.christoffel_derivative[direction][a][b][c]
                        - reference_derivative[direction][a][b][c]
                    )
    tilde, tilde_derivative = _batch_auxiliary_with_derivative(
        physical, tilde_normal_factor
    )
    constraint = tuple(
        -sum(
            (tilde[rho][sigma] * gamma_difference[a][rho][sigma] for rho in range(4) for sigma in range(4)),
            BatchScalar(0.0),
        )
        for a in range(4)
    )
    partial = tuple(
        tuple(
            -sum(
                (
                    tilde_derivative[direction][rho][sigma] * gamma_difference[a][rho][sigma]
                    + tilde[rho][sigma] * gamma_difference_derivative[direction][a][rho][sigma]
                    for rho in range(4)
                    for sigma in range(4)
                ),
                BatchScalar(0.0),
            )
            for a in range(4)
        )
        for direction in range(4)
    )
    covariant = tuple(
        tuple(
            partial[direction][a]
            + sum(
                (physical.christoffel[a][direction][source] * constraint[source] for source in range(4)),
                BatchScalar(0.0),
            )
            for a in range(4)
        )
        for direction in range(4)
    )
    constraint_down = tuple(
        sum((physical.metric[a][b] * constraint[b] for b in range(4)), BatchScalar(0.0))
        for a in range(4)
    )
    freeze3 = tuple(
        tuple(tuple(item for item in lower) for lower in upper)
        for upper in gamma_difference
    )
    freeze4 = tuple(
        tuple(
            tuple(tuple(item for item in lower) for lower in upper)
            for upper in direction
        )
        for direction in gamma_difference_derivative
    )
    return ModifiedHarmonicGaugeData(
        reference_id="flat_spherical_annulus",
        coordinate_radius=radius,  # type: ignore[arg-type]
        physical=physical,
        gamma_difference=freeze3,  # type: ignore[arg-type]
        gamma_difference_derivative=freeze4,  # type: ignore[arg-type]
        tilde_inverse_metric=tilde,
        tilde_inverse_metric_derivative=tilde_derivative,
        constraint_up=constraint,  # type: ignore[arg-type]
        constraint_down=constraint_down,  # type: ignore[arg-type]
        partial_constraint_up=partial,  # type: ignore[arg-type]
        covariant_constraint_derivative=covariant,  # type: ignore[arg-type]
    )


def _batch_extension(
    state: SphericalState,
    gauge: ModifiedHarmonicGaugeData,
    *,
    hat_normal_factor: int,
) -> tuple[tuple[object, ...], ...]:
    hat, _ = _batch_auxiliary_with_derivative(
        gauge.physical, hat_normal_factor
    )
    effective_planck = state.planck_mass**2 + state.beta * state.phi.value**2
    if np.any(BatchScalar.coerce(effective_planck).data <= 0.0):
        raise ValueError("modified-harmonic grid extension requires positive F")

    def projector(alpha: int, beta: int, mu: int, nu: int) -> object:
        return (
            int(alpha == mu) * hat[nu][beta]
            + int(alpha == nu) * hat[mu][beta]
            - int(alpha == beta) * hat[mu][nu]
        ) / 2

    extension_up = _zero2()
    for mu in range(4):
        for nu in range(4):
            extension_up[mu][nu] = effective_planck * sum(
                (
                    projector(alpha, beta, mu, nu)
                    * gauge.covariant_constraint_derivative[beta][alpha]
                    for alpha in range(4)
                    for beta in range(4)
                ),
                BatchScalar(0.0),
            )
    metric = gauge.physical.metric
    extension_down = _zero2()
    for a in range(4):
        for b in range(4):
            extension_down[a][b] = sum(
                (
                    metric[a][mu] * metric[b][nu] * extension_up[mu][nu]
                    for mu in range(4)
                    for nu in range(4)
                ),
                BatchScalar(0.0),
            )
    return tuple(tuple(row) for row in extension_down)


def batch_ref1_residual(
    u: object,
    p: object,
    q: object,
    dtt: object,
    p_r: object,
    q_r: object,
    radii: object,
    *,
    action: ActionParameters,
    tilde_normal_factor: int = 4,
    hat_normal_factor: int = 9,
) -> np.ndarray:
    """Return ``(batch, point, equation)`` complete REF1 residuals."""

    arrays = {
        name: _finite_array(name, value, ndim=2)
        for name, value in (("u", u), ("p", p), ("q", q), ("p_r", p_r), ("q_r", q_r))
    }
    if len({value.shape for value in arrays.values()}) != 1:
        raise ValueError("grid state arrays must share one shape")
    points, fields = arrays["u"].shape
    if fields != FIELD_COUNT:
        raise ValueError("grid state must use the canonical six ADM fields")
    acceleration = _finite_array("dtt", dtt, ndim=3)
    if acceleration.shape[1:] != (points, fields):
        raise ValueError("dtt must have (batch,point,field) shape")
    coordinate = _finite_array("radii", radii, ndim=1)
    if coordinate.shape != (points,) or np.any(coordinate <= 0.0):
        raise ValueError("annular radii must be strictly positive and match the grid")
    action = _coerce_action(action)
    if tilde_normal_factor <= 1 or hat_normal_factor <= tilde_normal_factor:
        raise ValueError("auxiliary normal factors must satisfy 1<tilde<hat")
    if np.any(arrays["u"][:, FIELD_INDEX["alpha"]] <= 0.0):
        raise ValueError("grid lapse must be positive")
    if np.any(arrays["u"][:, FIELD_INDEX["lambda"]] <= 0.0):
        raise ValueError("grid radial metric must be positive")
    if np.any(arrays["u"][:, FIELD_INDEX["R"]] <= 0.0):
        raise ValueError("annular areal radius must be positive")
    effective_planck = (
        action.planck_mass**2
        + action.ricci_coupling * arrays["u"][:, FIELD_INDEX["phi"]] ** 2
    )
    if np.any(effective_planck <= 0.0):
        raise ValueError("grid effective Planck coefficient must be positive")

    state = _batch_state_from_adm(
        arrays["u"], arrays["p"], arrays["q"], acceleration,
        arrays["p_r"], arrays["q_r"], action=action
    )
    radius = BatchScalar(coordinate[None, :])
    gauge = _batch_gauge(state, radius, tilde_normal_factor=tilde_normal_factor)
    extension = _batch_extension(
        state, gauge, hat_normal_factor=hat_normal_factor
    )
    original = residuals(state)
    full_metric = tuple(
        tuple(original["metric"][a][b] + extension[a][b] for b in range(4))
        for a in range(4)
    )
    vector = _independent_metric_vector(full_metric) + (original["phi"], original["chi"])
    rows = []
    for equation, value in enumerate(vector):
        scalar = BatchScalar.coerce(value).data
        try:
            broadcast = np.broadcast_to(scalar, (acceleration.shape[0], points))
        except ValueError as exc:
            raise ValueError(f"residual equation {equation} has incompatible batch shape") from exc
        rows.append(broadcast)
    answer = np.moveaxis(np.asarray(rows, dtype=np.float64), 0, -1)
    if answer.shape != (acceleration.shape[0], points, FIELD_COUNT) or not np.all(np.isfinite(answer)):
        raise ValueError("batch REF1 residual has invalid shape or values")
    return answer


@dataclass(frozen=True, slots=True)
class GridAccelerationResult:
    accelerations: np.ndarray
    residuals: np.ndarray
    residual_infinity: float
    initial_residual_infinity: float
    condition_infinity_maximum: float
    branch_displacement_infinity: float
    iterations: int
    residual_decreased_monotonically: bool
    affine_in_accelerations_verified: bool

    def __post_init__(self) -> None:
        acceleration = _finite_array("accelerations", self.accelerations, ndim=2)
        residual = _finite_array("residuals", self.residuals, ndim=2)
        if acceleration.shape != residual.shape or acceleration.shape[1] != FIELD_COUNT:
            raise ValueError("grid acceleration result shape differs")
        for name, value in (("accelerations", acceleration), ("residuals", residual)):
            immutable = value.copy()
            immutable.setflags(write=False)
            object.__setattr__(self, name, immutable)
        for name in (
            "residual_infinity",
            "initial_residual_infinity",
            "condition_infinity_maximum",
            "branch_displacement_infinity",
        ):
            value = float(getattr(self, name))
            if not isfinite(value) or value < 0.0:
                raise ValueError(f"{name} must be finite and nonnegative")
            object.__setattr__(self, name, value)
        if self.iterations not in {0, 1}:
            raise ValueError("quasilinear grid solve must take zero or one Newton update")


@dataclass(frozen=True, slots=True)
class CenterAccelerationLimit:
    acceleration: np.ndarray
    three_point_acceleration: np.ndarray
    estimator_infinity: float

    def __post_init__(self) -> None:
        values = _finite_array("center acceleration", self.acceleration, ndim=1)
        comparison = _finite_array(
            "three-point center acceleration", self.three_point_acceleration, ndim=1
        )
        if values.shape != (FIELD_COUNT,) or comparison.shape != (FIELD_COUNT,):
            raise ValueError("center acceleration vectors must contain six fields")
        for name, value in (("acceleration", values), ("three_point_acceleration", comparison)):
            immutable = value.copy()
            immutable.setflags(write=False)
            object.__setattr__(self, name, immutable)
        estimator = float(self.estimator_infinity)
        if not isfinite(estimator) or estimator < 0.0:
            raise ValueError("center acceleration estimator must be finite and nonnegative")
        object.__setattr__(self, "estimator_infinity", estimator)


def regular_center_acceleration_limit(
    positive_radius_accelerations: object,
) -> CenterAccelerationLimit:
    """Take the parity-regular annular acceleration limit at ``r=0``.

    Four positive, equally spaced grid values determine the even centre value
    through the degree-six polynomial weights ``(8/5,-4/5,8/35,-1/35)``.
    The independent three-point degree-four estimate uses
    ``(3/2,-3/5,1/10)``.  Odd primitive accelerations vanish exactly.  This
    never evaluates the singular spherical reference connection at zero.
    """

    values = _finite_array(
        "positive_radius_accelerations", positive_radius_accelerations, ndim=2
    )
    if values.shape[0] < 4 or values.shape[1] != FIELD_COUNT:
        raise ValueError("center limit requires four positive rows in six-field order")
    four = (
        8.0 / 5.0 * values[0]
        - 4.0 / 5.0 * values[1]
        + 8.0 / 35.0 * values[2]
        - 1.0 / 35.0 * values[3]
    )
    three = 1.5 * values[0] - 0.6 * values[1] + 0.1 * values[2]
    four[ADM_CENTER_PARITIES == -1] = 0.0
    three[ADM_CENTER_PARITIES == -1] = 0.0
    return CenterAccelerationLimit(
        acceleration=four,
        three_point_acceleration=three,
        estimator_infinity=float(np.max(np.abs(four - three), initial=0.0)),
    )


def adm_accelerations_from_base_metric(
    adm_u: Sequence[Real],
    adm_p: Sequence[Real],
    base_metric_accelerations: Sequence[Real],
) -> np.ndarray:
    """Invert the exact ADM-to-``h_AB`` second-time-derivative chain rule."""

    u = _finite_array("adm_u", adm_u, ndim=1)
    p = _finite_array("adm_p", adm_p, ndim=1)
    h = _finite_array("base_metric_accelerations", base_metric_accelerations, ndim=1)
    if u.shape != (FIELD_COUNT,) or p.shape != (FIELD_COUNT,) or h.shape != (FIELD_COUNT,):
        raise ValueError("ADM/base acceleration vectors must contain six fields")
    lapse = u[FIELD_INDEX["alpha"]]
    shift = u[FIELD_INDEX["shift"]]
    radial = u[FIELD_INDEX["lambda"]]
    lapse_t = p[FIELD_INDEX["alpha"]]
    shift_t = p[FIELD_INDEX["shift"]]
    radial_t = p[FIELD_INDEX["lambda"]]
    if lapse <= 0.0 or radial <= 0.0:
        raise ValueError("ADM acceleration inversion requires positive lapse and lambda")
    h_tt_tt, h_tr_tt, h_rr_tt, radius_tt, phi_tt, chi_tt = h
    radial_tt = (h_rr_tt - 2.0 * radial_t**2) / (2.0 * radial)
    shift_tt = (
        h_tr_tt
        - shift * h_rr_tt
        - 4.0 * radial * radial_t * shift_t
    ) / radial**2
    lapse_tt = (
        -h_tt_tt
        - 2.0 * lapse_t**2
        + h_rr_tt * shift**2
        + 8.0 * radial * radial_t * shift * shift_t
        + 2.0 * radial**2 * (shift_t**2 + shift * shift_tt)
    ) / (2.0 * lapse)
    answer = np.asarray(
        (lapse_tt, shift_tt, radial_tt, radius_tt, phi_tt, chi_tt),
        dtype=np.float64,
    )
    if not np.all(np.isfinite(answer)):
        raise ValueError("ADM acceleration inversion produced nonfinite values")
    return answer


def adm_lower_jets_from_base_state(
    state: SphericalState,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Invert a scalar ``h_AB`` two-jet into ADM ``u,p,q,p_r,q_r`` vectors."""

    if not isinstance(state, SphericalState):
        raise TypeError("state must be SphericalState")

    def scalar(value: object) -> float:
        if isinstance(value, BatchScalar):
            if value.data.shape != ():
                raise ValueError("ADM scalar inversion received a batch value")
            answer = float(value.data)
        else:
            answer = float(scalar_primal(value))
        if not isfinite(answer):
            raise ValueError("ADM scalar inversion received a nonfinite value")
        return answer

    htt = state.h_tt
    htr = state.h_tr
    hrr = state.h_rr
    hrr0 = scalar(hrr.value)
    if hrr0 <= 0.0:
        raise ValueError("ADM inversion requires h_rr>0")
    radial = np.sqrt(hrr0)
    shift = scalar(htr.value) / hrr0
    lapse_squared = -scalar(htt.value) + hrr0 * shift**2
    if lapse_squared <= 0.0:
        raise ValueError("ADM inversion requires positive lapse squared")
    lapse = np.sqrt(lapse_squared)

    hrr_t, hrr_r = scalar(hrr.dt), scalar(hrr.dr)
    htr_t, htr_r = scalar(htr.dt), scalar(htr.dr)
    htt_t, htt_r = scalar(htt.dt), scalar(htt.dr)
    radial_t = hrr_t / (2.0 * radial)
    radial_r = hrr_r / (2.0 * radial)
    shift_t = (htr_t - hrr_t * shift) / hrr0
    shift_r = (htr_r - hrr_r * shift) / hrr0
    lapse_t = (
        -htt_t + hrr_t * shift**2 + 2.0 * hrr0 * shift * shift_t
    ) / (2.0 * lapse)
    lapse_r = (
        -htt_r + hrr_r * shift**2 + 2.0 * hrr0 * shift * shift_r
    ) / (2.0 * lapse)

    hrr_tr, hrr_rr = scalar(hrr.dtr), scalar(hrr.drr)
    htr_tr, htr_rr = scalar(htr.dtr), scalar(htr.drr)
    htt_tr, htt_rr = scalar(htt.dtr), scalar(htt.drr)
    radial_tr = (hrr_tr - 2.0 * radial_t * radial_r) / (2.0 * radial)
    radial_rr = (hrr_rr - 2.0 * radial_r**2) / (2.0 * radial)
    shift_tr = (
        htr_tr
        - hrr_tr * shift
        - hrr_t * shift_r
        - hrr_r * shift_t
    ) / hrr0
    shift_rr = (
        htr_rr - hrr_rr * shift - 2.0 * hrr_r * shift_r
    ) / hrr0
    lapse_tr = (
        -htt_tr
        - 2.0 * lapse_t * lapse_r
        + hrr_tr * shift**2
        + 2.0 * hrr_t * shift * shift_r
        + 2.0 * hrr_r * shift * shift_t
        + 2.0 * hrr0 * (shift_t * shift_r + shift * shift_tr)
    ) / (2.0 * lapse)
    lapse_rr = (
        -htt_rr
        - 2.0 * lapse_r**2
        + hrr_rr * shift**2
        + 4.0 * hrr_r * shift * shift_r
        + 2.0 * hrr0 * (shift_r**2 + shift * shift_rr)
    ) / (2.0 * lapse)

    u = np.asarray(
        (
            lapse,
            shift,
            radial,
            scalar(state.areal_radius.value),
            scalar(state.phi.value),
            scalar(state.chi.value),
        )
    )
    p = np.asarray(
        (
            lapse_t,
            shift_t,
            radial_t,
            scalar(state.areal_radius.dt),
            scalar(state.phi.dt),
            scalar(state.chi.dt),
        )
    )
    q = np.asarray(
        (
            lapse_r,
            shift_r,
            radial_r,
            scalar(state.areal_radius.dr),
            scalar(state.phi.dr),
            scalar(state.chi.dr),
        )
    )
    p_r = np.asarray(
        (
            lapse_tr,
            shift_tr,
            radial_tr,
            scalar(state.areal_radius.dtr),
            scalar(state.phi.dtr),
            scalar(state.chi.dtr),
        )
    )
    q_r = np.asarray(
        (
            lapse_rr,
            shift_rr,
            radial_rr,
            scalar(state.areal_radius.drr),
            scalar(state.phi.drr),
            scalar(state.chi.drr),
        )
    )
    if not all(np.all(np.isfinite(value)) for value in (u, p, q, p_r, q_r)):
        raise ValueError("ADM lower-jet inversion produced nonfinite values")
    return u, p, q, p_r, q_r


def solve_grid_accelerations(
    u: object,
    p: object,
    q: object,
    p_r: object,
    q_r: object,
    radii: object,
    *,
    action: ActionParameters,
    warm_start: object | None = None,
    residual_tolerance: Real = 1.0e-12,
    condition_number_maximum: Real = 1.0e10,
    branch_displacement_maximum: Real = 1.0 / 16.0,
) -> GridAccelerationResult:
    """Solve all positive-radius ADM accelerations in one batched transaction."""

    base = _finite_array("u", u, ndim=2)
    points, fields = base.shape
    if fields != FIELD_COUNT:
        raise ValueError("u must contain six canonical ADM fields")
    initial = (
        np.zeros_like(base)
        if warm_start is None
        else _finite_array("warm_start", warm_start, ndim=2)
    )
    if initial.shape != base.shape:
        raise ValueError("warm_start shape differs from u")
    tolerance = float(residual_tolerance)
    condition_limit = float(condition_number_maximum)
    displacement_limit = float(branch_displacement_maximum)
    if not all(isfinite(value) and value > 0.0 for value in (tolerance, condition_limit, displacement_limit)):
        raise ValueError("grid solver thresholds must be finite and positive")

    seeded = np.broadcast_to(initial[None, :, :], (FIELD_COUNT + 1, points, fields)).copy()
    for field in range(FIELD_COUNT):
        seeded[field + 1, :, field] += 1.0
    residual_seed = batch_ref1_residual(
        base, p, q, seeded, p_r, q_r, radii, action=action
    )
    initial_residual = residual_seed[0]
    jacobian = np.empty((points, FIELD_COUNT, FIELD_COUNT), dtype=np.float64)
    for field in range(FIELD_COUNT):
        jacobian[:, :, field] = residual_seed[field + 1] - initial_residual
    try:
        inverse = np.linalg.inv(jacobian)
    except np.linalg.LinAlgError as exc:
        raise ValueError("grid kinetic block is singular") from exc
    condition = np.max(np.sum(np.abs(jacobian), axis=2), axis=1) * np.max(
        np.sum(np.abs(inverse), axis=2), axis=1
    )
    if not np.all(np.isfinite(condition)) or np.any(condition >= condition_limit):
        raise ValueError("grid kinetic condition limit reached")
    initial_norm = float(np.max(np.abs(initial_residual), initial=0.0))
    if initial_norm <= tolerance:
        candidate = initial
        iterations = 0
        displacement = 0.0
    else:
        correction = np.linalg.solve(
            jacobian, -initial_residual[..., None]
        )[..., 0]
        candidate = initial + correction
        displacement = float(np.max(np.abs(correction), initial=0.0))
        if displacement >= displacement_limit:
            raise ValueError("grid branch displacement limit reached")
        iterations = 1
    verified = batch_ref1_residual(
        base, p, q, candidate[None, :, :], p_r, q_r, radii, action=action
    )[0]
    final_norm = float(np.max(np.abs(verified), initial=0.0))
    monotonic = final_norm < initial_norm if iterations else final_norm <= initial_norm
    if final_norm > tolerance or not monotonic:
        raise ValueError("grid acceleration residual did not reach the strict tolerance")

    # A second seed scale catches accidental nonlinear dependence or a broken
    # batch axis without rerunning a pointwise finite-difference loop.
    half_seeded = np.broadcast_to(initial[None, :, :], (FIELD_COUNT + 1, points, fields)).copy()
    for field in range(FIELD_COUNT):
        half_seeded[field + 1, :, field] += 0.5
    half_residual = batch_ref1_residual(
        base, p, q, half_seeded, p_r, q_r, radii, action=action
    )
    half_jacobian = np.empty_like(jacobian)
    for field in range(FIELD_COUNT):
        half_jacobian[:, :, field] = 2.0 * (half_residual[field + 1] - half_residual[0])
    affine_error = float(np.max(np.abs(half_jacobian - jacobian), initial=0.0))
    affine_scale = max(1.0, float(np.max(np.abs(jacobian), initial=0.0)))
    affine_verified = affine_error <= 4096.0 * np.finfo(np.float64).eps * affine_scale
    if not affine_verified:
        raise ValueError("REF1 grid residual is not affine in ADM accelerations")
    return GridAccelerationResult(
        accelerations=candidate,
        residuals=verified,
        residual_infinity=final_norm,
        initial_residual_infinity=initial_norm,
        condition_infinity_maximum=float(np.max(condition, initial=0.0)),
        branch_displacement_infinity=displacement,
        iterations=iterations,
        residual_decreased_monotonically=monotonic,
        affine_in_accelerations_verified=True,
    )


__all__ = [
    "BatchJet2",
    "BatchScalar",
    "CenterAccelerationLimit",
    "ADM_CENTER_PARITIES",
    "FIELD_COUNT",
    "FIELD_INDEX",
    "GridAccelerationResult",
    "adm_accelerations_from_base_metric",
    "adm_lower_jets_from_base_state",
    "batch_ref1_residual",
    "regular_center_acceleration_limit",
    "solve_grid_accelerations",
]
