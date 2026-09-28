"""Outcome-blind local sensitivity of the complete spherical Q expression.

The input is a box of four ADM metric two-jets and a radial tangent. All
curvature contractions use the full four-dimensional Ricci tensor, including
the term that the exact-null Hess(R) shortcut would omit off the null cone.
Interval first tangents bound every partial derivative throughout the box;
the mean-value bound does not fit an error to the measured Q.

This pure arithmetic instrument establishes neither the supplied input-error
radii nor null/affine transport, a trajectory, or a physical DEF1 result.
Returned numbers are not self-authenticating evidence records.
"""

from __future__ import annotations

from fractions import Fraction
from typing import Mapping

from .exact_interval import Interval
from .interval_tangent import IntervalFirstTangent, primal_and_tangent
from .spherical_reduction import (
    Jet2,
    SphericalState,
    _ricci,
    _warped_2plus2_metric_connection,
    warped_2plus2_curvature,
)


Q = Fraction
METRIC_FIELDS = ("alpha", "shift", "lambda", "R")
JET_COMPONENTS = ("value", "dt", "dr", "dtt", "dtr", "drr")
INPUT_NAMES = tuple(f"{field}.{part}" for field in METRIC_FIELDS
                    for part in JET_COMPONENTS) + ("k.t", "k.r")
_INPUT_KEYS = frozenset(INPUT_NAMES)
MAXIMUM_INPUT_RATIONAL_BITS = 4096
MAXIMUM_OUTPUT_RATIONAL_BITS = 262144


class QGeometryPremiseError(ValueError):
    """The declared annular metric/input box is incomplete or unsafe."""


class QGeometryResourceError(QGeometryPremiseError):
    """A prospective exact-arithmetic size limit was exceeded."""


def _exact(value: object, name: str) -> Fraction:
    if type(value) not in (int, Fraction):
        raise TypeError(f"{name} must be a Fraction or built-in int")
    result = Q(value)
    if max(result.numerator.bit_length(), result.denominator.bit_length()) > MAXIMUM_INPUT_RATIONAL_BITS:
        raise QGeometryResourceError(f"{name} exceeds the input rational limit")
    return result


def _inventory(value: object, label: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping) or set(value) != _INPUT_KEYS:
        raise QGeometryPremiseError(f"{label} requires exactly the 26 named inputs")
    return value


def _box(value: object) -> dict[str, Interval]:
    given = _inventory(value, "geometry box")
    result = {}
    for name in INPUT_NAMES:
        item = given[name]
        if type(item) is Interval:
            result[name] = Interval(_exact(item.lower, name), _exact(item.upper, name))
        else:
            result[name] = Interval.singleton(_exact(item, name))
    for name in ("alpha.value", "lambda.value", "R.value", "k.t"):
        if not result[name].strictly_positive():
            raise QGeometryPremiseError(f"{name} must be positive on the whole box")
    return result


def _bounded(interval: Interval, name: str) -> Interval:
    if max(abs(item.numerator).bit_length() for item in (interval.lower, interval.upper)) > MAXIMUM_OUTPUT_RATIONAL_BITS or max(item.denominator.bit_length() for item in (interval.lower, interval.upper)) > MAXIMUM_OUTPUT_RATIONAL_BITS:
        raise QGeometryResourceError(f"{name} exceeds the output rational limit")
    return interval


def _evaluate(values: Mapping[str, object]) -> dict[str, object]:
    jets = {
        field: Jet2(**{part: values[f"{field}.{part}"] for part in JET_COMPONENTS})
        for field in METRIC_FIELDS
    }
    alpha, shift, radial, radius = (jets[field] for field in METRIC_FIELDS)
    radial_squared = radial * radial
    # The scalar slots are unused by curvature. This is a metric-only use of
    # the geometric carrier, not an assertion that this metric solves GR-0.
    state = SphericalState(
        h_tt=-(alpha * alpha) + radial_squared * shift * shift,
        h_tr=radial_squared * shift,
        h_rr=radial_squared,
        areal_radius=radius,
        phi=Jet2.constant(0), chi=Jet2.constant(0),
        planck_mass=Q(1), mu=Q(1), g4=Q(1),
        alpha=Q(0), beta=Q(0), eta=Q(0), branch="GR-0",
    )
    metric, inverse, _connection = _warped_2plus2_metric_connection(state)
    ricci = _ricci(warped_2plus2_curvature(state), inverse)
    tangent = (values["k.t"], values["k.r"])
    ricci_kk = sum(ricci[a][b] * tangent[a] * tangent[b]
                   for a in range(2) for b in range(2))
    null_residual = sum(metric[a][b] * tangent[a] * tangent[b]
                        for a in range(2) for b in range(2))
    theta = 2 * (radius.dt * tangent[0] + radius.dr * tangent[1]) / radius.value
    return {
        "theta": theta,
        "ricci_kk": ricci_kk,
        "q": -(theta * theta) / 2 - ricci_kk,
        "metric_null_residual": null_residual,
    }


def evaluate_q_box(box: Mapping[str, object]) -> dict[str, Interval]:
    """Enclose theta, full R_kk, Q and metric nullness on a supplied box.

    Q here is its smooth geometric expression extended off the null cone.
    The actual Raychaudhuri interpretation separately requires a metric-null,
    affinely parametrized, radial twist-free congruence. An interval merely
    containing zero is not a proof of nullness.
    """

    given = _box(box)
    try:
        values = _evaluate(given)
    except (ValueError, ZeroDivisionError) as error:
        raise QGeometryPremiseError("could not certify the metric box's annular Lorentzian domain") from error
    return {name: _bounded(primal_and_tangent(value)[0], name)
            for name, value in values.items()}


def q_sensitivity_enclosures(box: Mapping[str, object]) -> dict[str, Interval]:
    """Enclose dQ/dx_i everywhere in the box, not just at its center."""

    given = _box(box)
    result = {}
    for selected in INPUT_NAMES:
        seeded = {
            name: IntervalFirstTangent(value, Interval.singleton(int(name == selected)))
            for name, value in given.items()
        }
        try:
            derivative = primal_and_tangent(_evaluate(seeded)["q"])[1]
        except (ValueError, ZeroDivisionError) as error:
            raise QGeometryPremiseError("could not certify the sensitivity box's Lorentzian domain") from error
        result[selected] = _bounded(derivative, selected)
    return result


def q_input_error_bound(
    nominal: Mapping[str, object], input_radii: Mapping[str, object]
) -> Fraction:
    """Conditional mean-value bound sum_i sup_box|dQ/dx_i| epsilon_i.

    Both mappings must specify all 26 exact inputs, including explicit zero
    radii. The box is nominal +/- epsilon, established before Q is evaluated.
    For any true input inside that box, this bounds |Q(true)-Q(nominal)|.
    Correlated errors are allowed; no independence or favorable cancellation
    is assumed. The caller still owns proof that its true input is enclosed.
    """

    given = _inventory(nominal, "nominal geometry")
    radii = _inventory(input_radii, "input radii")
    center = {name: _exact(given[name], name) for name in INPUT_NAMES}
    error = {name: _exact(radii[name], name) for name in INPUT_NAMES}
    if any(value < 0 for value in error.values()):
        raise QGeometryPremiseError("input radii must be nonnegative")
    box = {name: Interval(center[name] - error[name], center[name] + error[name])
           for name in INPUT_NAMES}
    derivatives = q_sensitivity_enclosures(box)
    bound = sum((derivatives[name].abs_upper() * error[name] for name in INPUT_NAMES), Q(0))
    return _bounded(Interval.singleton(bound), "Q input-error bound").upper
