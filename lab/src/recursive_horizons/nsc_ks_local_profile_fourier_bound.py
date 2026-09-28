"""Local Fourier enclosures for compact Chebyshev/cutoff profile products.

The earlier profile owner obtains retained coefficients from a sampled DFT
and bounds both aliasing and the omitted band with one global derivative-L1
estimate.  That estimate is valid but unusably large on the flat cutoff
ramps of the current Chebyshev-32 history.

This module instead performs two directed operations:

* each retained coefficient is integrated directly with ``acb.integral``;
  the exponentially flat endpoint strips are replaced by their limiting
  values and their omitted contribution is enclosed explicitly;
* the difference between the analytic profile and the retained Fourier
  polynomial is bounded on local panels by midpoint Taylor expansion and a
  directed derivative remainder.

No fitted Fourier decay or sampled maximum is used.  The returned tail
bounds are uniform on the complete numerical period and can be consumed by
``polynomial_residual_bounds``.
"""
from dataclasses import dataclass
from hashlib import sha256
import json
from math import factorial

from flint import acb, arb, ctx

from . import nsc_ks_ball_geometry as G
from .nsc_ks_chebyshev_jet_bound import (
    actual_profile_functions,
    chebyshev_endpoint_derivative_bounds,
    monomial_interval_series,
)
from .nsc_local_incoming_family import LocalIncomingFamily


@dataclass(frozen=True)
class LocalProfileFourierBound:
    key: tuple
    coefficients: tuple
    uniform_tail_bounds: tuple
    coefficient_flat_strip_error: object
    retained_index: int
    taylor_order: int
    local_cells: int
    maximum_subdivision_depth: int
    period_origin: object
    period_length: object


def _arb_fraction(value):
    return arb(value.numerator) / arb(value.denominator)


def _pack_arb(value, digits=60):
    value = arb(value)
    return {
        "mid": value.mid().str(digits),
        "radius_upper": value.rad().upper().str(digits),
    }


def _unpack_arb(value, name="Arb ball"):
    if (not isinstance(value, dict) or set(value) != {"mid", "radius_upper"}
            or not isinstance(value["mid"], str)
            or not isinstance(value["radius_upper"], str)):
        raise ValueError("malformed serialized " + name)
    midpoint, radius = arb(value["mid"]), arb(value["radius_upper"])
    if not midpoint.is_finite() or not radius.is_finite() or not radius >= 0:
        raise ValueError("finite nonnegative serialized " + name + " required")
    return midpoint + radius * arb(0, 1)


def serialize_local_profile_fourier(rows):
    """Lossless-for-proof outward serialization of local Fourier balls."""
    result = {}
    for key, row in rows.items():
        result[f"{key[0]}_{key[1]}"] = {
            "key": list(key),
            "coefficients": [
                {"real": _pack_arb(value.real), "imag": _pack_arb(value.imag)}
                for value in row.coefficients
            ],
            "uniform_tail_bounds": [_pack_arb(value)
                                    for value in row.uniform_tail_bounds],
            "coefficient_flat_strip_error": _pack_arb(
                row.coefficient_flat_strip_error),
            "retained_index": row.retained_index,
            "taylor_order": row.taylor_order,
            "local_cells": row.local_cells,
            "maximum_subdivision_depth": row.maximum_subdivision_depth,
            "period_origin": _pack_arb(row.period_origin),
            "period_length": _pack_arb(row.period_length),
        }
    return result


def profile_payload_digest(value):
    """SHA-256 of the canonical JSON of a serialized local Fourier payload."""
    if not isinstance(value, dict) or not value:
        raise ValueError("serialized local Fourier payload required")
    return sha256(json.dumps(
        value, sort_keys=True, separators=(",", ":"),
        allow_nan=False).encode("utf-8")).hexdigest()


def restore_local_profile_fourier(value, *, expected_sha256=None):
    """Restore a conservative superset of serialized coefficient balls.

    A supplied ``expected_sha256`` must match ``profile_payload_digest``.
    A mismatched hash is a restore failure, not a scientific verdict.
    """
    if not isinstance(value, dict) or not value:
        raise ValueError("serialized local Fourier payload required")
    digest = profile_payload_digest(value)
    if expected_sha256 is not None:
        if (not isinstance(expected_sha256, str)
                or len(expected_sha256) != 64 or digest != expected_sha256):
            raise ValueError("mismatched profile payload")
    result = {}
    required = {
        "key", "coefficients", "uniform_tail_bounds",
        "coefficient_flat_strip_error", "retained_index", "taylor_order",
        "local_cells", "maximum_subdivision_depth", "period_origin",
        "period_length",
    }
    for label, packed in value.items():
        if not isinstance(packed, dict) or set(packed) != required:
            raise ValueError("malformed serialized local Fourier row")
        raw_key = packed["key"]
        if (not isinstance(raw_key, list) or len(raw_key) != 2
                or any(isinstance(item, bool) or not isinstance(item, int)
                       or item < 0 for item in raw_key)
                or sum(raw_key) < 1):
            raise ValueError("malformed serialized profile key")
        key = tuple(raw_key)
        if label != f"{key[0]}_{key[1]}" or key in result:
            raise ValueError("serialized profile key mismatch")
        integers = {}
        for name in ("retained_index", "taylor_order", "local_cells",
                     "maximum_subdivision_depth"):
            item = packed[name]
            if (isinstance(item, bool) or not isinstance(item, int)
                    or item < (0 if name == "maximum_subdivision_depth" else 1)):
                raise ValueError("malformed serialized " + name)
            integers[name] = item
        raw_coefficients = packed["coefficients"]
        if (not isinstance(raw_coefficients, list)
                or len(raw_coefficients) != 2 * integers["retained_index"] + 1):
            raise ValueError("serialized coefficient band does not match retained index")
        coefficients = []
        for item in raw_coefficients:
            if not isinstance(item, dict) or set(item) != {"real", "imag"}:
                raise ValueError("malformed serialized complex coefficient")
            coefficients.append(acb(
                _unpack_arb(item["real"], "coefficient real ball"),
                _unpack_arb(item["imag"], "coefficient imaginary ball")))
        raw_tails = packed["uniform_tail_bounds"]
        if not isinstance(raw_tails, list) or not raw_tails:
            raise ValueError("serialized derivative tail bounds required")
        tails = tuple(_unpack_arb(item, "tail bound") for item in raw_tails)
        if any(not item >= 0 for item in tails):
            raise ValueError("nonnegative serialized tail bounds required")
        flat = _unpack_arb(
            packed["coefficient_flat_strip_error"], "flat-strip error")
        if not flat >= 0:
            raise ValueError("nonnegative serialized flat-strip error required")
        period_origin = _unpack_arb(packed["period_origin"], "period origin")
        period_length = _unpack_arb(packed["period_length"], "period length")
        if not period_length.lower() > 0:
            raise ValueError("positive serialized period length required")
        result[key] = LocalProfileFourierBound(
            key, tuple(coefficients), tails, flat,
            integers["retained_index"], integers["taylor_order"],
            integers["local_cells"], integers["maximum_subdivision_depth"],
            period_origin, period_length,
        )
    return result


def _validate_keys(keys):
    keys = tuple(keys)
    if (not keys or len(set(keys)) != len(keys)
            or any(not isinstance(key, tuple) or len(key) != 2
                   or min(key) < 0 or sum(key) < 1 for key in keys)):
        raise ValueError("distinct nonnegative positive-degree (p,q) keys required")
    return keys


def _chebyshev_acb(profile, z):
    """Evaluate the owned finite Chebyshev polynomial on one complex ball."""
    offset, scale = profile._derivatives[0].mapparms()
    x = acb(float(offset)) + acb(float(scale)) * z
    if len(profile.coefficients) == 1:
        return acb(profile.coefficients[0])
    b1 = b2 = acb(0)
    for coefficient in reversed(profile.coefficients[1:]):
        b1, b2 = acb(coefficient) + 2 * x * b1 - b2, b1
    return acb(profile.coefficients[0]) + x * b1 - b2


def _transition_eta(z, center, inner, width, *, left):
    distance = acb(center) - z if left else z - acb(center)
    u = (distance - acb(inner)) / acb(width)
    q = 1 / (1 - u) - 1 / u
    return 1 / (1 + q.exp())


def _flat_endpoint_error(family, keys, epsilon):
    """Integral error for replacing the inner/outer flat strips by 1/0."""
    w_profile, u_profile = actual_profile_functions(family)
    w_bound = chebyshev_endpoint_derivative_bounds(w_profile, 1)["polynomial_bounds"][0]
    u_bound = chebyshev_endpoint_derivative_bounds(u_profile, 1)["polynomial_bounds"][0]
    q = 1 / (1 - epsilon) - 1 / epsilon
    tiny = q.exp() / (1 + q.exp())
    width = arb(w_profile.outer) - arb(w_profile.inner)
    result = {}
    for key in keys:
        degree = sum(key)
        polynomial = (_arb_fraction(w_bound) ** key[0]
                      * _arb_fraction(u_bound) ** key[1])
        # Two inner strips: 1-eta^degree <= degree*(1-eta).
        # Two outer strips: eta^degree <= eta <= tiny.
        result[key] = (2 * epsilon * width * polynomial * (degree + 1) * tiny).upper()
    return result


def _positive_coefficients(family, keys, origin, length, retained_index,
                           flat_denominator, *, rel_bits, abs_bits):
    """Direct nonnegative Fourier coefficients with endpoint-strip error."""
    w_profile, u_profile = actual_profile_functions(family)
    center = arb(w_profile.center)
    inner = arb(w_profile.inner)
    outer = arb(w_profile.outer)
    width = outer - inner
    epsilon = arb(1) / flat_denominator
    flat = _flat_endpoint_error(family, keys, epsilon)
    result = {key: [] for key in keys}

    def polynomial(z, key):
        return _chebyshev_acb(w_profile, z) ** key[0] * _chebyshev_acb(
            u_profile, z) ** key[1]

    for mode in range(retained_index + 1):
        wave = 2 * arb.pi() * mode / length

        def phase(z):
            return (-acb(0, wave) * (z - acb(origin))).exp()

        for key in keys:
            degree = sum(key)
            total = acb.integral(
                lambda z, _analytic: polynomial(z, key) * phase(z),
                acb(center - inner - epsilon * width),
                acb(center + inner + epsilon * width),
                rel_tol=arb(2) ** -rel_bits,
                abs_tol=arb(2) ** -abs_bits,
            )
            for left, right, is_left in (
                (center - outer + epsilon * width,
                 center - inner - epsilon * width, True),
                (center + inner + epsilon * width,
                 center + outer - epsilon * width, False),
            ):
                total += acb.integral(
                    lambda z, _analytic, side=is_left: (
                        _transition_eta(z, center, inner, width, left=side) ** degree
                        * polynomial(z, key) * phase(z)
                    ),
                    acb(left), acb(right),
                    rel_tol=arb(2) ** -rel_bits,
                    abs_tol=arb(2) ** -abs_bits,
                )
            error = flat[key] / length
            result[key].append(total / length + acb(arb(0, error), arb(0, error)))
    return result, {key: (flat[key] / length).upper() for key in keys}


def _panel_layout(origin, length, center, inner, outer, *, exterior_panels,
                  transition_panels, interior_panels):
    endpoints = (origin, center - outer, center - inner,
                 center + inner, center + outer, origin + length)
    counts = (exterior_panels, transition_panels, interior_panels,
              transition_panels, exterior_panels)
    intervals = []
    for left, right, count in zip(endpoints[:-1], endpoints[1:], counts):
        if not right >= left:
            raise ValueError("profile support must lie inside the numerical period")
        for index in range(count):
            intervals.append((left + (right - left) * index / count,
                              left + (right - left) * (index + 1) / count, 0))
    return intervals


def _partial_derivative(positive, waves, relative_z, order):
    value = positive[0] if order == 0 else acb(0)
    for mode in range(1, len(positive)):
        phase = acb(0, waves[mode] * relative_z).exp()
        value += ((acb(0, waves[mode]) ** order) * positive[mode] * phase
                  + (-acb(0, waves[mode])) ** order
                  * positive[mode].conjugate() / phase)
    return value


def _local_tail_bounds(family, keys, positive, origin, length, *, taylor_order,
                       exterior_panels, transition_panels, interior_panels,
                       max_derivative, max_depth):
    w_profile, u_profile = actual_profile_functions(family)
    highest = max_derivative + taylor_order + 1
    endpoint_cache = (
        chebyshev_endpoint_derivative_bounds(w_profile, max(highest + 1, 9)),
        chebyshev_endpoint_derivative_bounds(u_profile, max(highest + 1, 9)),
    )
    retained_index = len(next(iter(positive.values()))) - 1
    waves = tuple(2 * arb.pi() * mode / length for mode in range(retained_index + 1))
    intervals = _panel_layout(
        origin, length, arb(w_profile.center), arb(w_profile.inner), arb(w_profile.outer),
        exterior_panels=exterior_panels, transition_panels=transition_panels,
        interior_panels=interior_panels,
    )
    maxima = {key: [arb(0) for _ in range(max_derivative + 1)] for key in keys}
    cells = depth_seen = 0
    while intervals:
        left, right, depth = intervals.pop()
        middle = (left + right) / 2
        halfwidth = (right - left) / 2
        try:
            center_series = monomial_interval_series(
                family, middle, keys, highest - 1, endpoint_cache)
            full_series = monomial_interval_series(
                family, left.union(right), keys, highest, endpoint_cache)
        except G.SubdivisionNeeded:
            if depth >= max_depth:
                raise
            split = (left + right) / 2
            intervals.extend(((left, split, depth + 1),
                              (split, right, depth + 1)))
            depth_seen = max(depth_seen, depth + 1)
            continue
        relative = middle - origin
        for key in keys:
            center_coefficients = G._coeffs(center_series[key], highest - 1)
            full_coefficients = G._coeffs(full_series[key], highest)
            for target in range(max_derivative + 1):
                bound = arb(0)
                for offset in range(taylor_order + 1):
                    order = target + offset
                    analytic = center_coefficients[order] * factorial(order)
                    difference = acb(analytic) - _partial_derivative(
                        positive[key], waves, relative, order)
                    bound += (difference.abs_upper() * halfwidth ** offset
                              / factorial(offset))
                remainder_order = target + taylor_order + 1
                analytic_remainder = abs(
                    full_coefficients[remainder_order]
                    * factorial(remainder_order)).upper()
                fourier_remainder = sum((
                    2 * abs(waves[mode]) ** remainder_order
                    * positive[key][mode].abs_upper()
                    for mode in range(1, retained_index + 1)), arb(0))
                bound += ((analytic_remainder + fourier_remainder)
                          * halfwidth ** (taylor_order + 1)
                          / factorial(taylor_order + 1))
                maxima[key][target] = arb.max(maxima[key][target], bound.upper())
        cells += 1
    return {key: tuple(value.upper() for value in row)
            for key, row in maxima.items()}, cells, depth_seen


def enclose_profile_fourier_local(
        family, keys, period_origin, period_length, *, retained_index=256,
        taylor_order=7, exterior_panels=512, transition_panels=256,
        interior_panels=512, max_derivative=2, flat_denominator=256,
        relative_tolerance_bits=45, absolute_tolerance_bits=65,
        max_depth=8, bits=90):
    """Return direct retained coefficients and local uniform tail bounds."""
    if not isinstance(family, LocalIncomingFamily):
        raise TypeError("LocalIncomingFamily required")
    keys = _validate_keys(keys)
    integer_values = (retained_index, taylor_order, exterior_panels,
                      transition_panels, interior_panels, max_derivative,
                      flat_denominator, relative_tolerance_bits,
                      absolute_tolerance_bits, max_depth, bits)
    if any(isinstance(value, bool) or not isinstance(value, int) or value < 1
           for value in integer_values):
        raise ValueError("positive integer enclosure settings required")
    if taylor_order < 1 or max_derivative > 2:
        raise ValueError("positive Taylor order and derivative orders through two required")
    with ctx.workprec(bits):
        origin, length = arb(period_origin), arb(period_length)
        if not length > 0:
            raise ValueError("positive numerical period required")
        w_profile, _ = actual_profile_functions(family)
        if not (origin < arb(w_profile.center) - arb(w_profile.outer)
                and origin + length > arb(w_profile.center) + arb(w_profile.outer)):
            raise ValueError("compact support must lie strictly inside the period")
        positive, flat_error = _positive_coefficients(
            family, keys, origin, length, retained_index, flat_denominator,
            rel_bits=relative_tolerance_bits, abs_bits=absolute_tolerance_bits)
        tails, cells, depth = _local_tail_bounds(
            family, keys, positive, origin, length, taylor_order=taylor_order,
            exterior_panels=exterior_panels, transition_panels=transition_panels,
            interior_panels=interior_panels, max_derivative=max_derivative,
            max_depth=max_depth)
        result = {}
        for key in keys:
            coefficients = tuple(
                positive[key][-mode].conjugate() if mode < 0 else positive[key][mode]
                for mode in range(-retained_index, retained_index + 1))
            result[key] = LocalProfileFourierBound(
                key, coefficients, tails[key], flat_error[key], retained_index,
                taylor_order, cells, depth, origin, length)
        return result
