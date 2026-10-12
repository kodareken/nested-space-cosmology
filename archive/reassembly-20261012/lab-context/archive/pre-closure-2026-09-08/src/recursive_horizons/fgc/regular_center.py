"""Exact regular-centre series for the spherical REF1 formulation.

The inherited REF1 evaluator is deliberately defined on a positive-radius
annulus.  Spherical coordinates are singular at ``r=0`` even when the
geometry is smooth, so extending that evaluator by substituting a small
floating radius would not define a centre formulation.  This module instead
uses exact truncated Laurent series and the regular variables

``shift = r*V`` and ``areal_radius = r*A``.

The primitive fields ``alpha, lambda, A, V, phi, chi`` and all of their time
derivatives are even in ``r``.  Elementary flatness is imposed coefficient by
coefficient at the centre through ``A=lambda`` for the value, first time
derivative, and second time derivative.  The complete unredefined ACT1/VAR1
residual is evaluated through the sealed tensor evaluator, while the affected
flat-spherical reference-connection descendants are rederived directly in
Laurent arithmetic.  Only coefficients within the declared exact series
window are retained; the public certificate inspects a much smaller safe
window whose coefficients cannot depend on discarded terms.

This is a local centre formulation and exact regression surface.  It is not a
constraint-compatible initial-data family, an existence theorem, a boundary
condition, a numerical evolution, or a retained-EFT authorization.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from numbers import Integral
from typing import Any, Mapping, Sequence

from .modified_harmonic import _trace_reversal_projector
from .modified_harmonic_reference import physical_connection_data
from .spherical_reduction import residuals


Q = Fraction
N = 4
SERIES_MIN_POWER = -12
SERIES_MAX_POWER = 20
CERTIFIED_MIN_POWER = -4
CERTIFIED_MAX_POWER = 4
REGULAR_PRIMITIVE_ORDER = ("alpha", "lambda", "A", "V", "phi", "chi")
TIME_LAYER_ORDER = ("value", "dt", "dtt")
REGULAR_EQUATION_ORDER = (
    "metric_tt",
    "metric_tr_over_r",
    "metric_rr",
    "metric_angular_over_r_squared",
    "scalar_phi",
    "scalar_chi",
)


def _fraction(name: str, value: object) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (Fraction, Integral)):
        raise TypeError(f"{name} must be an exact rational")
    return value if isinstance(value, Fraction) else Q(int(value))


@dataclass(frozen=True, slots=True, eq=False)
class LaurentSeries:
    """One exact Laurent series on a fixed finite coefficient window."""

    coefficients: tuple[Fraction, ...]

    def __post_init__(self) -> None:
        width = SERIES_MAX_POWER - SERIES_MIN_POWER + 1
        if len(self.coefficients) != width:
            raise ValueError(f"Laurent series requires exactly {width} coefficients")
        object.__setattr__(
            self,
            "coefficients",
            tuple(_fraction(f"coefficients[{index}]", value) for index, value in enumerate(self.coefficients)),
        )

    @classmethod
    def zero(cls) -> "LaurentSeries":
        return cls((Q(0),) * (SERIES_MAX_POWER - SERIES_MIN_POWER + 1))

    @classmethod
    def constant(cls, value: int | Fraction) -> "LaurentSeries":
        return cls.monomial(value, 0)

    @classmethod
    def monomial(cls, value: int | Fraction, power: int) -> "LaurentSeries":
        coefficient = _fraction("monomial coefficient", value)
        if power < SERIES_MIN_POWER or power > SERIES_MAX_POWER:
            if coefficient == 0:
                return cls.zero()
            raise ValueError("monomial power lies outside the exact series window")
        values = [Q(0)] * (SERIES_MAX_POWER - SERIES_MIN_POWER + 1)
        values[power - SERIES_MIN_POWER] = coefficient
        return cls(tuple(values))

    @classmethod
    def from_mapping(cls, values: Mapping[int, int | Fraction]) -> "LaurentSeries":
        output = cls.zero()
        coefficients = list(output.coefficients)
        for power, value in values.items():
            if isinstance(power, bool) or not isinstance(power, int):
                raise TypeError("Laurent powers must be integers")
            coefficient = _fraction(f"coefficient[{power}]", value)
            if power < SERIES_MIN_POWER or power > SERIES_MAX_POWER:
                if coefficient:
                    raise ValueError("nonzero Laurent coefficient lies outside the exact window")
                continue
            coefficients[power - SERIES_MIN_POWER] += coefficient
        return cls(tuple(coefficients))

    @staticmethod
    def _coerce(value: object) -> "LaurentSeries" | Any:
        if isinstance(value, LaurentSeries):
            return value
        if isinstance(value, bool) or not isinstance(value, (Fraction, Integral)):
            return NotImplemented
        return LaurentSeries.constant(value if isinstance(value, Fraction) else int(value))

    def coefficient(self, power: int) -> Fraction:
        if power < SERIES_MIN_POWER or power > SERIES_MAX_POWER:
            return Q(0)
        return self.coefficients[power - SERIES_MIN_POWER]

    def coefficient_map(
        self, *, minimum: int = CERTIFIED_MIN_POWER, maximum: int = CERTIFIED_MAX_POWER
    ) -> dict[int, Fraction]:
        if minimum > maximum:
            raise ValueError("coefficient-map bounds are reversed")
        return {
            power: self.coefficient(power)
            for power in range(minimum, maximum + 1)
            if self.coefficient(power) != 0
        }

    def is_zero(self) -> bool:
        return all(value == 0 for value in self.coefficients)

    def valuation(self) -> int | None:
        for offset, value in enumerate(self.coefficients):
            if value:
                return SERIES_MIN_POWER + offset
        return None

    def leading_coefficient(self) -> Fraction:
        power = self.valuation()
        return Q(0) if power is None else self.coefficient(power)

    def has_negative_power(self, *, minimum: int = SERIES_MIN_POWER) -> bool:
        return any(self.coefficient(power) != 0 for power in range(minimum, 0))

    def parity(self, expected: str, *, minimum: int = 0, maximum: int = CERTIFIED_MAX_POWER) -> bool:
        if expected not in {"even", "odd"}:
            raise ValueError("series parity must be even or odd")
        rejected = 1 if expected == "even" else 0
        return all(self.coefficient(power) == 0 for power in range(minimum, maximum + 1) if power % 2 == rejected)

    def __neg__(self) -> "LaurentSeries":
        return LaurentSeries(tuple(-value for value in self.coefficients))

    def __add__(self, other: object) -> "LaurentSeries" | Any:
        value = self._coerce(other)
        if value is NotImplemented:
            return NotImplemented
        return LaurentSeries(tuple(a + b for a, b in zip(self.coefficients, value.coefficients)))

    __radd__ = __add__

    def __sub__(self, other: object) -> "LaurentSeries" | Any:
        value = self._coerce(other)
        if value is NotImplemented:
            return NotImplemented
        return self + (-value)

    def __rsub__(self, other: object) -> "LaurentSeries" | Any:
        value = self._coerce(other)
        return NotImplemented if value is NotImplemented else value - self

    def __mul__(self, other: object) -> "LaurentSeries" | Any:
        value = self._coerce(other)
        if value is NotImplemented:
            return NotImplemented
        if self.is_zero() or value.is_zero():
            return LaurentSeries.zero()
        output = [Q(0)] * len(self.coefficients)
        left = tuple(
            (SERIES_MIN_POWER + index, coefficient)
            for index, coefficient in enumerate(self.coefficients)
            if coefficient
        )
        right = tuple(
            (SERIES_MIN_POWER + index, coefficient)
            for index, coefficient in enumerate(value.coefficients)
            if coefficient
        )
        for left_power, left_value in left:
            for right_power, right_value in right:
                power = left_power + right_power
                if power < SERIES_MIN_POWER:
                    raise ArithmeticError(
                        "Laurent product crossed the exact negative-power guard edge"
                    )
                if SERIES_MIN_POWER <= power <= SERIES_MAX_POWER:
                    output[power - SERIES_MIN_POWER] += left_value * right_value
        return LaurentSeries(tuple(output))

    __rmul__ = __mul__

    def reciprocal(self) -> "LaurentSeries":
        leading_power = self.valuation()
        if leading_power is None:
            raise ZeroDivisionError("zero Laurent series has no reciprocal")
        leading = self.coefficient(leading_power)
        first_output_power = -leading_power
        if not SERIES_MIN_POWER <= first_output_power <= SERIES_MAX_POWER:
            raise ValueError("reciprocal leading power lies outside the exact window")
        maximum_index = SERIES_MAX_POWER - first_output_power
        inverse_coefficients = [Q(0)] * (maximum_index + 1)
        inverse_coefficients[0] = 1 / leading
        for order in range(1, maximum_index + 1):
            convolution = sum(
                self.coefficient(leading_power + source_order)
                * inverse_coefficients[order - source_order]
                for source_order in range(1, order + 1)
            )
            inverse_coefficients[order] = -convolution / leading
        return LaurentSeries.from_mapping(
            {
                first_output_power + order: coefficient
                for order, coefficient in enumerate(inverse_coefficients)
            }
        )

    def __truediv__(self, other: object) -> "LaurentSeries" | Any:
        value = self._coerce(other)
        if value is NotImplemented:
            return NotImplemented
        return self * value.reciprocal()

    def __rtruediv__(self, other: object) -> "LaurentSeries" | Any:
        value = self._coerce(other)
        return NotImplemented if value is NotImplemented else value / self

    def __pow__(self, exponent: object) -> "LaurentSeries" | Any:
        if isinstance(exponent, bool) or not isinstance(exponent, Integral):
            return NotImplemented
        power = int(exponent)
        if power == 0:
            return LaurentSeries.constant(1)
        if power < 0:
            return self.reciprocal() ** (-power)
        result = LaurentSeries.constant(1)
        base = self
        while power:
            if power & 1:
                result = result * base
            power //= 2
            if power:
                base = base * base
        return result

    def derivative(self) -> "LaurentSeries":
        if self.coefficient(SERIES_MIN_POWER):
            raise ArithmeticError(
                "Laurent derivative crossed the exact negative-power guard edge"
            )
        return LaurentSeries.from_mapping(
            {
                power - 1: power * self.coefficient(power)
                for power in range(SERIES_MIN_POWER, SERIES_MAX_POWER + 1)
                if power and self.coefficient(power) and SERIES_MIN_POWER <= power - 1
            }
        )

    def evaluate(self, radius: int | Fraction) -> Fraction:
        value = _fraction("radius", radius)
        if value == 0 and any(self.coefficient(power) for power in range(SERIES_MIN_POWER, 0)):
            raise ZeroDivisionError("cannot evaluate a singular Laurent series at zero")
        return sum(
            (self.coefficient(power) * value**power for power in range(SERIES_MIN_POWER, SERIES_MAX_POWER + 1)),
            Q(0),
        )

    def _comparison(self, other: object, operation: str) -> bool | Any:
        value = self._coerce(other)
        if value is NotImplemented:
            return NotImplemented
        difference = self - value
        leading = difference.leading_coefficient()
        if operation == "lt":
            return leading < 0
        if operation == "le":
            return leading <= 0
        if operation == "gt":
            return leading > 0
        if operation == "ge":
            return leading >= 0
        raise AssertionError("unknown comparison")

    def __eq__(self, other: object) -> bool:
        value = self._coerce(other)
        return False if value is NotImplemented else self.coefficients == value.coefficients

    def __lt__(self, other: object) -> bool | Any:
        return self._comparison(other, "lt")

    def __le__(self, other: object) -> bool | Any:
        return self._comparison(other, "le")

    def __gt__(self, other: object) -> bool | Any:
        return self._comparison(other, "gt")

    def __ge__(self, other: object) -> bool | Any:
        return self._comparison(other, "ge")


@dataclass(frozen=True, slots=True)
class SeriesJet2:
    """A two-jet whose six entries are exact Laurent series."""

    value: LaurentSeries
    dt: LaurentSeries
    dr: LaurentSeries
    dtt: LaurentSeries
    dtr: LaurentSeries
    drr: LaurentSeries

    @classmethod
    def from_time_layers(
        cls,
        value: LaurentSeries,
        dt: LaurentSeries,
        dtt: LaurentSeries,
    ) -> "SeriesJet2":
        return cls(value, dt, value.derivative(), dtt, dt.derivative(), value.derivative().derivative())

    @classmethod
    def constant(cls, value: int | Fraction) -> "SeriesJet2":
        series = LaurentSeries.constant(value)
        zero = LaurentSeries.zero()
        return cls(series, zero, zero, zero, zero, zero)

    @staticmethod
    def _coerce(value: object) -> "SeriesJet2" | Any:
        if isinstance(value, SeriesJet2):
            return value
        if isinstance(value, bool) or not isinstance(value, (Fraction, Integral)):
            return NotImplemented
        return SeriesJet2.constant(value if isinstance(value, Fraction) else int(value))

    def __neg__(self) -> "SeriesJet2":
        return SeriesJet2(*(-getattr(self, name) for name in ("value", "dt", "dr", "dtt", "dtr", "drr")))

    def __add__(self, other: object) -> "SeriesJet2" | Any:
        value = self._coerce(other)
        if value is NotImplemented:
            return NotImplemented
        return SeriesJet2(*(getattr(self, name) + getattr(value, name) for name in ("value", "dt", "dr", "dtt", "dtr", "drr")))

    __radd__ = __add__

    def __sub__(self, other: object) -> "SeriesJet2" | Any:
        value = self._coerce(other)
        return NotImplemented if value is NotImplemented else self + (-value)

    def __rsub__(self, other: object) -> "SeriesJet2" | Any:
        value = self._coerce(other)
        return NotImplemented if value is NotImplemented else value - self

    def __mul__(self, other: object) -> "SeriesJet2" | Any:
        value = self._coerce(other)
        if value is NotImplemented:
            return NotImplemented
        a, b = self, value
        return SeriesJet2(
            a.value * b.value,
            a.dt * b.value + a.value * b.dt,
            a.dr * b.value + a.value * b.dr,
            a.dtt * b.value + 2 * a.dt * b.dt + a.value * b.dtt,
            a.dtr * b.value + a.dt * b.dr + a.dr * b.dt + a.value * b.dtr,
            a.drr * b.value + 2 * a.dr * b.dr + a.value * b.drr,
        )

    __rmul__ = __mul__


@dataclass(frozen=True, slots=True)
class RegularSeriesState:
    h_tt: SeriesJet2
    h_tr: SeriesJet2
    h_rr: SeriesJet2
    areal_radius: SeriesJet2
    phi: SeriesJet2
    chi: SeriesJet2
    planck_mass: Fraction
    beta: Fraction
    mu: Fraction
    g4: Fraction
    alpha: Fraction
    eta: Fraction
    branch: str


def _coefficient_layer(name: str, value: object) -> dict[int, Fraction]:
    if not isinstance(value, Mapping):
        raise TypeError(f"{name} must be a power-to-rational mapping")
    output: dict[int, Fraction] = {}
    for power, coefficient in value.items():
        if isinstance(power, bool) or not isinstance(power, int) or power < 0:
            raise ValueError(f"{name} powers must be nonnegative integers")
        if power > CERTIFIED_MAX_POWER + 4:
            raise ValueError(f"{name} power exceeds the declared profile window")
        exact = _fraction(f"{name}[{power}]", coefficient)
        if exact:
            output[power] = exact
    return output


def validate_regular_profile(profile: Mapping[str, Any]) -> dict[str, Any]:
    """Validate parity, elementary flatness, action, and exact profile scope."""

    if not isinstance(profile, Mapping):
        raise TypeError("regular-centre profile must be a mapping")
    expected = {"profile_id", "fields", "action"}
    if set(profile) != expected:
        raise ValueError("regular-centre profile keys differ")
    profile_id = profile["profile_id"]
    if not isinstance(profile_id, str) or not profile_id:
        raise ValueError("regular-centre profile requires a nonempty identifier")
    fields = profile["fields"]
    if not isinstance(fields, Mapping) or set(fields) != set(REGULAR_PRIMITIVE_ORDER):
        raise ValueError("regular-centre primitive field set differs")
    parsed: dict[str, dict[str, dict[int, Fraction]]] = {}
    for field in REGULAR_PRIMITIVE_ORDER:
        layers = fields[field]
        if not isinstance(layers, Mapping) or set(layers) != set(TIME_LAYER_ORDER):
            raise ValueError(f"{field} must define value, dt, and dtt layers")
        parsed[field] = {}
        for layer in TIME_LAYER_ORDER:
            coefficients = _coefficient_layer(f"fields.{field}.{layer}", layers[layer])
            if any(power % 2 for power in coefficients):
                raise ValueError(f"{field}.{layer} violates even centre parity")
            parsed[field][layer] = coefficients

    for field in ("alpha", "lambda", "A"):
        if parsed[field]["value"].get(0, Q(0)) <= 0:
            raise ValueError(f"{field}(0) must be strictly positive")
    for layer in TIME_LAYER_ORDER:
        if parsed["A"][layer].get(0, Q(0)) != parsed["lambda"][layer].get(0, Q(0)):
            raise ValueError(f"elementary flatness A=lambda fails in {layer} layer")

    action = profile["action"]
    required_action = {"branch", "planck_mass", "beta", "mu", "g4", "alpha", "eta"}
    if not isinstance(action, Mapping) or set(action) != required_action:
        raise ValueError("regular-centre action keys differ")
    exact_action = {
        key: _fraction(f"action.{key}", action[key])
        for key in ("planck_mass", "beta", "mu", "g4", "alpha", "eta")
    }
    branch = action["branch"]
    if branch != "FGC-QR" or exact_action != {
        "planck_mass": Q(2),
        "beta": -Q(1, 4),
        "mu": Q(3),
        "g4": Q(1, 2),
        "alpha": Q(0),
        "eta": Q(1, 2),
    }:
        raise ValueError("regular-centre control must use the frozen ACT1 FGC-QR action")
    phi0 = parsed["phi"]["value"].get(0, Q(0))
    effective_planck = exact_action["planck_mass"] ** 2 + exact_action["beta"] * phi0**2
    if effective_planck <= 0:
        raise ValueError("regular-centre profile violates F>0")
    return {
        "profile_id": profile_id,
        "fields": parsed,
        "action": {"branch": branch, **exact_action},
        "parity_valid": True,
        "elementary_flatness_value_dt_dtt_valid": True,
        "effective_planck_at_center": effective_planck,
    }


def _series_layer(parsed: Mapping[str, Any], field: str, layer: str) -> LaurentSeries:
    return LaurentSeries.from_mapping(parsed["fields"][field][layer])


def regular_series_state(profile: Mapping[str, Any]) -> tuple[RegularSeriesState, dict[str, SeriesJet2]]:
    parsed = validate_regular_profile(profile)
    primitive: dict[str, SeriesJet2] = {}
    for field in REGULAR_PRIMITIVE_ORDER:
        primitive[field] = SeriesJet2.from_time_layers(
            _series_layer(parsed, field, "value"),
            _series_layer(parsed, field, "dt"),
            _series_layer(parsed, field, "dtt"),
        )
    radius = SeriesJet2.from_time_layers(
        LaurentSeries.monomial(1, 1), LaurentSeries.zero(), LaurentSeries.zero()
    )
    shift = radius * primitive["V"]
    areal_radius = radius * primitive["A"]
    lam2 = primitive["lambda"] * primitive["lambda"]
    state = RegularSeriesState(
        h_tt=-(primitive["alpha"] * primitive["alpha"]) + lam2 * shift * shift,
        h_tr=lam2 * shift,
        h_rr=lam2,
        areal_radius=areal_radius,
        phi=primitive["phi"],
        chi=primitive["chi"],
        **parsed["action"],
    )
    return state, primitive


def _zero2() -> list[list[LaurentSeries]]:
    return [[LaurentSeries.zero() for _ in range(N)] for _ in range(N)]


def _zero3() -> list[list[list[LaurentSeries]]]:
    return [[[LaurentSeries.zero() for _ in range(N)] for _ in range(N)] for _ in range(N)]


def _zero4() -> list[list[list[list[LaurentSeries]]]]:
    return [[[[LaurentSeries.zero() for _ in range(N)] for _ in range(N)] for _ in range(N)] for _ in range(N)]


def _freeze(value: Any) -> Any:
    if isinstance(value, list):
        return tuple(_freeze(item) for item in value)
    return value


def _flat_reference_series() -> tuple[Any, Any]:
    """Return flat-spherical ``bar Gamma`` and its coordinate derivative."""

    radius = LaurentSeries.monomial(1, 1)
    inverse_radius = radius.reciprocal()
    gamma = _zero3()
    derivative = _zero4()
    _, radial, theta, phi = range(N)
    gamma[radial][theta][theta] = -radius
    gamma[radial][phi][phi] = -radius
    gamma[theta][radial][theta] = gamma[theta][theta][radial] = inverse_radius
    gamma[phi][radial][phi] = gamma[phi][phi][radial] = inverse_radius
    derivative[radial][radial][theta][theta] = LaurentSeries.constant(-1)
    derivative[radial][radial][phi][phi] = LaurentSeries.constant(-1)
    derivative[radial][theta][radial][theta] = derivative[radial][theta][theta][radial] = -inverse_radius**2
    derivative[radial][phi][radial][phi] = derivative[radial][phi][phi][radial] = -inverse_radius**2
    derivative[theta][theta][phi][phi] = LaurentSeries.constant(1)
    derivative[theta][phi][theta][phi] = LaurentSeries.constant(-1)
    derivative[theta][phi][phi][theta] = LaurentSeries.constant(-1)
    return _freeze(gamma), _freeze(derivative)


def _generic_auxiliary_inverse_with_derivative(physical: Any, factor: Fraction) -> tuple[Any, Any]:
    if factor <= 1:
        raise ValueError("auxiliary normal factor must exceed one")
    inverse = physical.inverse_metric
    inverse_derivative = physical.inverse_metric_derivative
    normal_outer = tuple(
        tuple(-inverse[a][0] * inverse[b][0] / inverse[0][0] for b in range(N))
        for a in range(N)
    )
    auxiliary = tuple(
        tuple(inverse[a][b] - (factor - 1) * normal_outer[a][b] for b in range(N))
        for a in range(N)
    )
    normal_outer_derivative = _zero3()
    for direction in range(N):
        denominator = inverse[0][0]
        denominator_derivative = inverse_derivative[direction][0][0]
        for a in range(N):
            for b in range(N):
                numerator = inverse[a][0] * inverse[b][0]
                numerator_derivative = (
                    inverse_derivative[direction][a][0] * inverse[b][0]
                    + inverse[a][0] * inverse_derivative[direction][b][0]
                )
                normal_outer_derivative[direction][a][b] = (
                    -numerator_derivative / denominator
                    + numerator * denominator_derivative / denominator**2
                )
    auxiliary_derivative = tuple(
        tuple(
            tuple(
                inverse_derivative[direction][a][b]
                - (factor - 1) * normal_outer_derivative[direction][a][b]
                for b in range(N)
            )
            for a in range(N)
        )
        for direction in range(N)
    )
    return auxiliary, auxiliary_derivative


def _series_gauge_data(state: RegularSeriesState, physical: Any, *, tilde_normal_factor: Fraction) -> dict[str, Any]:
    reference, reference_derivative = _flat_reference_series()
    gamma_difference = _zero3()
    gamma_difference_derivative = _zero4()
    for a in range(N):
        for b in range(N):
            for c in range(N):
                gamma_difference[a][b][c] = physical.christoffel[a][b][c] - reference[a][b][c]
                for direction in range(N):
                    gamma_difference_derivative[direction][a][b][c] = (
                        physical.christoffel_derivative[direction][a][b][c]
                        - reference_derivative[direction][a][b][c]
                    )
    tilde, tilde_derivative = _generic_auxiliary_inverse_with_derivative(physical, tilde_normal_factor)
    constraint = tuple(
        -sum(
            (tilde[rho][sigma] * gamma_difference[a][rho][sigma] for rho in range(N) for sigma in range(N)),
            LaurentSeries.zero(),
        )
        for a in range(N)
    )
    partial_constraint = tuple(
        tuple(
            -sum(
                (
                    tilde_derivative[direction][rho][sigma] * gamma_difference[a][rho][sigma]
                    + tilde[rho][sigma] * gamma_difference_derivative[direction][a][rho][sigma]
                    for rho in range(N)
                    for sigma in range(N)
                ),
                LaurentSeries.zero(),
            )
            for a in range(N)
        )
        for direction in range(N)
    )
    covariant = tuple(
        tuple(
            partial_constraint[direction][a]
            + sum(
                (physical.christoffel[a][direction][source] * constraint[source] for source in range(N)),
                LaurentSeries.zero(),
            )
            for a in range(N)
        )
        for direction in range(N)
    )
    return {
        "reference_christoffel": reference,
        "reference_christoffel_derivative": reference_derivative,
        "gamma_difference": _freeze(gamma_difference),
        "gamma_difference_derivative": _freeze(gamma_difference_derivative),
        "tilde_inverse_metric": tilde,
        "tilde_inverse_metric_derivative": tilde_derivative,
        "constraint_up": constraint,
        "partial_constraint_up": partial_constraint,
        "covariant_constraint_derivative": covariant,
    }


def _series_extension(state: RegularSeriesState, physical: Any, gauge: Mapping[str, Any], *, hat_normal_factor: Fraction) -> tuple[Any, Any]:
    hat, _ = _generic_auxiliary_inverse_with_derivative(physical, hat_normal_factor)
    effective_planck = state.planck_mass**2 + state.beta * state.phi.value**2
    if effective_planck <= 0:
        raise ValueError("regular-centre extension requires F>0 near the centre")
    extension_up = _zero2()
    for mu in range(N):
        for nu in range(N):
            extension_up[mu][nu] = effective_planck * sum(
                (
                    _trace_reversal_projector(hat, alpha, beta, mu, nu)
                    * gauge["covariant_constraint_derivative"][beta][alpha]
                    for alpha in range(N)
                    for beta in range(N)
                ),
                LaurentSeries.zero(),
            )
    extension_down = _zero2()
    metric = physical.metric
    for a in range(N):
        for b in range(N):
            extension_down[a][b] = sum(
                (
                    metric[a][mu] * metric[b][nu] * extension_up[mu][nu]
                    for mu in range(N)
                    for nu in range(N)
                ),
                LaurentSeries.zero(),
            )
    return _freeze(extension_up), _freeze(extension_down)


def _all_series(value: Any) -> tuple[LaurentSeries, ...]:
    if isinstance(value, LaurentSeries):
        return (value,)
    if isinstance(value, (tuple, list)):
        output: list[LaurentSeries] = []
        for item in value:
            output.extend(_all_series(item))
        return tuple(output)
    raise TypeError(f"expected Laurent tensor, received {type(value).__name__}")


def _finite_tensor(value: Any) -> bool:
    return all(not series.has_negative_power() for series in _all_series(value))


def _parity_tensor(value: Sequence[LaurentSeries], parities: Sequence[str]) -> bool:
    return all(series.parity(parity) for series, parity in zip(value, parities))


def _ricci_squared(ricci: Any, inverse: Any) -> LaurentSeries:
    return sum(
        (
            inverse[a][c] * inverse[b][d] * ricci[a][b] * ricci[c][d]
            for a in range(N)
            for b in range(N)
            for c in range(N)
            for d in range(N)
        ),
        LaurentSeries.zero(),
    )


def regular_center_series_certificate(profile: Mapping[str, Any]) -> dict[str, Any]:
    """Evaluate and classify the complete regularized REF1 centre series."""

    parsed = validate_regular_profile(profile)
    state, primitive = regular_series_state(profile)
    unredefined = residuals(state)  # type: ignore[arg-type]
    physical = physical_connection_data(state)  # type: ignore[arg-type]
    gauge = _series_gauge_data(state, physical, tilde_normal_factor=Q(4))
    _, extension = _series_extension(state, physical, gauge, hat_normal_factor=Q(9))
    full_metric = tuple(
        tuple(unredefined["metric"][a][b] + extension[a][b] for b in range(N))
        for a in range(N)
    )
    radius = LaurentSeries.monomial(1, 1)
    raw_vector = (
        full_metric[0][0],
        full_metric[0][1],
        full_metric[1][1],
        full_metric[2][2],
        unredefined["phi"],
        unredefined["chi"],
    )
    regular_vector = (
        raw_vector[0],
        raw_vector[1] / radius,
        raw_vector[2],
        raw_vector[3] / radius**2,
        raw_vector[4],
        raw_vector[5],
    )
    regular_gauge = (
        gauge["constraint_up"][0],
        gauge["constraint_up"][1] / radius,
    )
    ricci_squared = _ricci_squared(unredefined["ricci"], physical.inverse_metric)
    kretschmann = unredefined["GB"] + 4 * ricci_squared - unredefined["R"] ** 2

    # Exact identities that remove the explicit reference 1/r and 1/r^2
    # components before a centre value is taken.
    A = primitive["A"].value
    R = state.areal_radius.value
    identity_r = R.derivative() / R - 1 / radius
    identity_r_expected = A.derivative() / A
    identity_rr = identity_r.derivative()
    identity_rr_expected = (A.derivative() / A).derivative()
    identity_t = state.areal_radius.dt / R
    identity_t_expected = primitive["A"].dt / A
    if identity_r != identity_r_expected or identity_rr != identity_rr_expected or identity_t != identity_t_expected:
        raise ValueError("regular reference quotient identities failed exactly")

    all_gamma_difference_finite = _finite_tensor(gauge["gamma_difference"])
    all_gamma_derivative_difference_finite = _finite_tensor(gauge["gamma_difference_derivative"])
    all_regular_equations_finite = _finite_tensor(regular_vector)
    all_regular_gauge_finite = _finite_tensor(regular_gauge)
    curvature_finite = _finite_tensor((unredefined["R"], ricci_squared, kretschmann, unredefined["GB"]))
    expected_equation_parity = ("even", "even", "even", "even", "even", "even")
    if not all(
        (
            all_gamma_difference_finite,
            all_gamma_derivative_difference_finite,
            all_regular_equations_finite,
            all_regular_gauge_finite,
            curvature_finite,
            _parity_tensor(regular_vector, expected_equation_parity),
            regular_gauge[0].parity("even"),
            regular_gauge[1].parity("even"),
        )
    ):
        raise ValueError("regular-centre Laurent cancellation, parity, or curvature gate failed")

    center_system = tuple(
        {
            "equation": name,
            "center_coefficient": series.coefficient(0),
            "next_even_coefficient": series.coefficient(2),
            "negative_coefficients": series.coefficient_map(minimum=SERIES_MIN_POWER, maximum=-1),
        }
        for name, series in zip(REGULAR_EQUATION_ORDER, regular_vector)
    )
    return {
        "profile_id": parsed["profile_id"],
        "series_window": {
            "computed_minimum_power": SERIES_MIN_POWER,
            "computed_maximum_power": SERIES_MAX_POWER,
            "certified_minimum_power": CERTIFIED_MIN_POWER,
            "certified_maximum_power": CERTIFIED_MAX_POWER,
            "all_computed_negative_powers_checked_for_zero": True,
            "structural_lowest_possible_final_pole": -4,
        },
        "regular_variables": {
            "areal_radius_equals_r_times_A": True,
            "shift_equals_r_times_V": True,
            "primitive_fields_even": tuple(REGULAR_PRIMITIVE_ORDER),
            "radial_derivatives_odd": True,
            "A0_equals_lambda0_in_value_dt_dtt_layers": parsed["elementary_flatness_value_dt_dtt_valid"],
        },
        "reference_regularization": {
            "reference_choice": "flat_spherical_connection_rederived_by_exact_punctured_series_limits",
            "R_r_over_R_minus_one_over_r_equals_A_r_over_A": identity_r == identity_r_expected,
            "radial_derivative_identity_exact": identity_rr == identity_rr_expected,
            "R_t_over_R_equals_A_t_over_A": identity_t == identity_t_expected,
            "all_connection_differences_have_no_certified_negative_powers": all_gamma_difference_finite,
            "all_connection_derivative_differences_have_no_certified_negative_powers": all_gamma_derivative_difference_finite,
            "affected_REF1_descendants_rederived": True,
        },
        "analytic_regularity_basis": {
            "base_two_geometry": "finite_even_or_odd_ADM_coefficients_and_nonzero_alpha0_lambda0",
            "mixed_warp_curvature": "nabla_A_nabla_B_R_over_R_has_a_finite_even_center_limit",
            "sphere_warp_curvature": "one_minus_grad_R_squared_over_R_squared_is_finite_iff_elementary_flatness_holds",
            "scalar_angular_hessian": "grad_R_dot_grad_scalar_over_R_is_finite_for_even_scalars",
            "reference_mixed_connection": "R_r_over_R_minus_one_over_r_equals_A_r_over_A",
            "reference_radial_derivative": "d_r_of_reference_mixed_connection_equals_d_r_A_r_over_A",
            "ACT1_VAR1_closure": "Einstein_nonminimal_scalar_stress_and_double_dual_GB_terms_are_algebraic_contractions_of_the_finite_basis",
            "REF1_closure": "auxiliary_metrics_projectors_and_nabla_C_are_algebraic_or_first_derivative_constructions_from_the_finite_difference_basis",
            "quadratic_curvature_sets_worst_final_pole_to_r_minus_4_before_cancellation": True,
        },
        "centre_taylor_system": center_system,
        "regular_equation_series": {
            name: series.coefficient_map()
            for name, series in zip(REGULAR_EQUATION_ORDER, regular_vector)
        },
        "regular_gauge_series": {
            "C_t": regular_gauge[0].coefficient_map(),
            "C_r_over_r": regular_gauge[1].coefficient_map(),
        },
        "curvature_series": {
            "Ricci_scalar": unredefined["R"].coefficient_map(),
            "Ricci_squared": ricci_squared.coefficient_map(),
            "Kretschmann": kretschmann.coefficient_map(),
            "Gauss_Bonnet": unredefined["GB"].coefficient_map(),
        },
        "checks": {
            "all_removable_reference_singularities_cancelled": all_gamma_difference_finite and all_gamma_derivative_difference_finite,
            "all_regularized_REF1_equations_finite": all_regular_equations_finite,
            "regularized_REF1_equations_even": _parity_tensor(regular_vector, expected_equation_parity),
            "regularized_gauge_components_finite_and_even": all_regular_gauge_finite and regular_gauge[0].parity("even") and regular_gauge[1].parity("even"),
            "Ricci_Ricci2_Kretschmann_GB_finite": curvature_finite,
        },
        "center_limits": tuple(series.coefficient(0) for series in regular_vector),
        "gauge_center_limits": tuple(series.coefficient(0) for series in regular_gauge),
        "nonclaims": {
            "compatible_initial_data_family_constructed": False,
            "smooth_solution_exists": False,
            "constraint_preserving_boundary_supplied": False,
            "evolution_authorized": False,
            "collapse_solution_derived": False,
            "affine_null_defocusing_derived": False,
            "retained_EFT_evolution_authorized": False,
        },
    }


def analytic_defect_ledger(profile: Mapping[str, Any]) -> dict[str, Any]:
    """Return leading analytic singular coefficients before validation.

    This helper intentionally accepts only the field portion of a profile and
    does not call :func:`validate_regular_profile`; mutation tests use it to
    show why conical or scalar-parity defects must fail closed.
    """

    fields = profile.get("fields") if isinstance(profile, Mapping) else None
    if not isinstance(fields, Mapping):
        raise TypeError("defect ledger requires profile fields")

    def coefficient(field: str, power: int) -> Fraction:
        layers = fields.get(field)
        if not isinstance(layers, Mapping):
            raise ValueError(f"defect ledger lacks {field}")
        values = layers.get("value")
        if not isinstance(values, Mapping):
            raise ValueError(f"defect ledger lacks {field}.value")
        return _fraction(f"{field}.value[{power}]", values.get(power, Q(0)))

    A0 = coefficient("A", 0)
    lambda0 = coefficient("lambda", 0)
    if A0 == 0 or lambda0 == 0:
        raise ValueError("defect ledger requires nonzero A0 and lambda0")
    phi1 = coefficient("phi", 1)
    chi1 = coefficient("chi", 1)
    return {
        "Ricci_scalar_r_minus_2_from_conical_defect": 2 * (1 - A0**2 / lambda0**2) / A0**2,
        "box_phi_r_minus_1_from_odd_scalar_term": 2 * phi1 / lambda0**2,
        "box_chi_r_minus_1_from_odd_scalar_term": 2 * chi1 / lambda0**2,
        "conical_defect_absent_iff_A0_squared_equals_lambda0_squared": A0**2 == lambda0**2,
        "positive_orientation_selects_A0_equals_lambda0": A0 > 0 and lambda0 > 0,
        "scalar_odd_defects_absent": phi1 == 0 and chi1 == 0,
    }


def _polynomial_derivative_value(
    coefficients: Mapping[int, Fraction], radius: Fraction, derivative_order: int
) -> Fraction:
    if derivative_order not in {0, 1, 2}:
        raise ValueError("only radial derivatives through second order are supported")
    total = Q(0)
    for power, coefficient in coefficients.items():
        if power < derivative_order:
            continue
        factor = Q(1)
        for offset in range(derivative_order):
            factor *= power - offset
        total += coefficient * factor * radius ** (power - derivative_order)
    return total


def annular_regularized_evaluation(
    profile: Mapping[str, Any],
    *,
    coordinate_radius: int | Fraction,
    reference_radial_minimum: int | Fraction,
) -> dict[str, Any]:
    """Evaluate the unchanged annular REF1 code on one regular polynomial jet."""

    from .modified_harmonic_reference import modified_harmonic_full_residuals
    from .reference_connection import flat_spherical_annulus_reference
    from .spherical_reduction import Jet2, state_from_generalized_adm_pg_fixture

    parsed = validate_regular_profile(profile)
    radius = _fraction("coordinate_radius", coordinate_radius)
    minimum = _fraction("reference_radial_minimum", reference_radial_minimum)
    if not 0 < minimum < radius:
        raise ValueError("annular control requires 0 < reference minimum < radius")

    def jet(field: str, *, multiply_by_r: bool = False) -> Jet2:
        layers = parsed["fields"][field]
        shifted = {
            layer: {power + 1: value for power, value in layers[layer].items()}
            if multiply_by_r
            else dict(layers[layer])
            for layer in TIME_LAYER_ORDER
        }
        return Jet2(
            _polynomial_derivative_value(shifted["value"], radius, 0),
            _polynomial_derivative_value(shifted["dt"], radius, 0),
            _polynomial_derivative_value(shifted["value"], radius, 1),
            _polynomial_derivative_value(shifted["dtt"], radius, 0),
            _polynomial_derivative_value(shifted["dt"], radius, 1),
            _polynomial_derivative_value(shifted["value"], radius, 2),
        )

    jets = {
        "alpha": jet("alpha"),
        "shift": jet("V", multiply_by_r=True),
        "lambda": jet("lambda"),
        "areal_radius": jet("A", multiply_by_r=True),
        "phi": jet("phi"),
        "chi": jet("chi"),
    }
    state_mapping = {
        name: {
            component: getattr(value, component)
            for component in ("value", "dt", "dr", "dtt", "dtr", "drr")
        }
        for name, value in jets.items()
    }
    action = parsed["action"]
    fixture = {
        "model_id": action["branch"],
        "action_parameters": {
            "planck_mass": action["planck_mass"],
            "ricci_coupling": action["beta"],
            "scalar_mass": action["mu"],
            "quartic_coupling": action["g4"],
            "linear_gb_coupling": action["alpha"],
            "quadratic_gb_coupling": action["eta"],
        },
        "state": state_mapping,
    }
    state = state_from_generalized_adm_pg_fixture(fixture)
    reference = flat_spherical_annulus_reference(radial_domain_minimum=minimum)
    full = modified_harmonic_full_residuals(
        state,
        reference=reference,
        coordinate_radius=radius,
        tilde_normal_factor=Q(4),
        hat_normal_factor=Q(9),
    )
    vector = full["full_residual_vector"]
    regular = (
        vector[0],
        vector[1] / radius,
        vector[2],
        vector[3] / radius**2,
        vector[4],
        vector[5],
    )
    gauge = (
        full["gauge"].constraint_up[0],
        full["gauge"].constraint_up[1] / radius,
    )
    return {
        "coordinate_radius": radius,
        "regular_equation_order": REGULAR_EQUATION_ORDER,
        "regular_equations": regular,
        "regular_gauge_order": ("C_t", "C_r_over_r"),
        "regular_gauge": gauge,
    }


def first_grid_point_convergence_certificate(
    profile: Mapping[str, Any],
    *,
    radii: Sequence[int | Fraction] = (Q(1, 16), Q(1, 32), Q(1, 64)),
    formal_certificate: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Compare first positive grid points with the exact centre Taylor limit.

    Every regularized field is even, so its point value differs from its
    centre value by ``O(h^2)``.  The frozen exact controls require the error to
    shrink by at least ``15/4`` when ``h`` is halved.  This is a strict finite
    check of the expected second-order centre approach, not a numerical PDE
    convergence claim.
    """

    exact_radii = tuple(_fraction(f"radii[{index}]", value) for index, value in enumerate(radii))
    if len(exact_radii) < 3 or any(value <= 0 for value in exact_radii):
        raise ValueError("centre convergence requires at least three positive radii")
    if any(exact_radii[index + 1] * 2 != exact_radii[index] for index in range(len(exact_radii) - 1)):
        raise ValueError("centre convergence radii must halve exactly")
    formal = (
        regular_center_series_certificate(profile)
        if formal_certificate is None
        else formal_certificate
    )
    if formal.get("profile_id") != profile.get("profile_id"):
        raise ValueError("formal centre certificate belongs to a different profile")
    names = REGULAR_EQUATION_ORDER + ("C_t", "C_r_over_r")
    limits = tuple(formal["center_limits"]) + tuple(formal["gauge_center_limits"])
    minimum = exact_radii[-1] / 4
    evaluations = tuple(
        annular_regularized_evaluation(
            profile,
            coordinate_radius=radius,
            reference_radial_minimum=minimum,
        )
        for radius in exact_radii
    )
    values = tuple(
        tuple(evaluation["regular_equations"]) + tuple(evaluation["regular_gauge"])
        for evaluation in evaluations
    )
    records: list[dict[str, Any]] = []
    for component, (name, limit) in enumerate(zip(names, limits)):
        errors = tuple(abs(row[component] - limit) for row in values)
        ratios: list[Fraction | None] = []
        for coarse, fine in zip(errors, errors[1:]):
            if fine == 0:
                ratios.append(None)
            else:
                ratios.append(coarse / fine)
        passed = all(
            (fine == 0 and coarse == 0)
            or (fine > 0 and coarse / fine >= Q(15, 4))
            for coarse, fine in zip(errors, errors[1:])
        )
        records.append(
            {
                "component": name,
                "exact_center_limit": limit,
                "values": tuple(row[component] for row in values),
                "absolute_errors": errors,
                "coarse_to_fine_error_ratios": tuple(ratios),
                "at_least_second_order_control_passed": passed,
            }
        )
    if not all(record["at_least_second_order_control_passed"] for record in records):
        raise ValueError("first-grid regular-centre convergence control failed")
    return {
        "radii": exact_radii,
        "halving_sequence": True,
        "minimum_required_error_ratio": Q(15, 4),
        "records": tuple(records),
        "all_regular_equation_and_gauge_controls_passed": True,
        "interpretation": "exact_annular_point_values_approach_formal_center_limits_at_least_quadratically",
        "not_a_PDE_discretization_convergence_claim": True,
    }
