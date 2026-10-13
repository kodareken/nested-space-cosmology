"""Binary64 replay envelope for the discrete Euler diagnostic, not physics.

The finite-term action, Fourier operator and physical acceptance criteria are
unchanged. This bounds differences in an uncancelled arithmetic expression
whose small remainder depends on the FFT and elementwise evaluation order.
"""
from functools import lru_cache
from math import ceil, log2

import numpy as np
import sympy as sp

from .nsc_covariant_operator import smooth_metric
from .nsc_finite_terms import _density_functions, _symbolic_basis, finite_response, periodic_derivative


def _gamma(operations):
    # Twice binary64 unit roundoff is conservative for scalar elementary
    # operations and the faithfully rounded transcendental metric inputs.
    accumulated = int(operations) * np.finfo(np.float64).eps
    if not 0 <= accumulated < 1:
        raise ValueError("invalid binary64 operation budget")
    return float(np.nextafter(accumulated / (1 - accumulated), np.inf))


def _outward(lower, upper):
    return np.nextafter(lower, -np.inf), np.nextafter(upper, np.inf)


def _interval(expression, boxes):
    """Outward interval evaluation of the owned rational Euler jet formula."""
    if expression in boxes:
        return boxes[expression]
    if expression.is_Number:
        value = float(expression)
        return _outward(np.asarray(value), np.asarray(value))
    children = [_interval(arg, boxes) for arg in expression.args]
    if expression.is_Add:
        lower, upper = children[0]
        for lo, hi in children[1:]:
            lower, upper = _outward(lower + lo, upper + hi)
        return lower, upper
    if expression.is_Mul:
        lower, upper = children[0]
        for lo, hi in children[1:]:
            terms = np.broadcast_arrays(lower * lo, lower * hi, upper * lo, upper * hi)
            lower, upper = _outward(np.minimum.reduce(terms), np.maximum.reduce(terms))
        return lower, upper
    if expression.is_Pow and expression.exp.is_Integer:
        lower, upper = children[0]
        exponent = int(expression.exp)
        if exponent < 0:
            if np.any((lower <= 0) & (upper >= 0)):
                raise ValueError("Euler jet enclosure crosses a denominator zero")
            lower, upper = _outward(1 / upper, 1 / lower)
            exponent = -exponent
        lo, hi = np.asarray(1.), np.asarray(1.)
        for _ in range(exponent):
            terms = np.broadcast_arrays(lo * lower, lo * upper, hi * lower, hi * upper)
            lo, hi = _outward(np.minimum.reduce(terms), np.maximum.reduce(terms))
        return lo, hi
    raise ValueError("unexpected nonrational Euler jet operation")


def _majorant(expression, boxes):
    """Absolute term sum before cancellation, with denominator lower bounds."""
    if expression in boxes:
        lower, upper = boxes[expression]
        return np.maximum(abs(lower), abs(upper))
    if expression.is_Number:
        return abs(float(expression))
    if expression.is_Add:
        result = 0.
        for arg in expression.args:
            result = np.nextafter(result + _majorant(arg, boxes), np.inf)
        return result
    if expression.is_Mul:
        result = 1.
        for arg in expression.args:
            result = np.nextafter(result * _majorant(arg, boxes), np.inf)
        return result
    if expression.is_Pow and expression.exp.is_Integer:
        exponent = int(expression.exp)
        if exponent >= 0:
            return _majorant(expression.base, boxes) ** exponent
        lower, upper = _interval(expression.base, boxes)
        if np.any((lower <= 0) & (upper >= 0)):
            raise ValueError("Euler majorant crosses a denominator zero")
        return np.minimum(abs(lower), abs(upper)) ** exponent
    raise ValueError("unexpected Euler majorant operation")


def _fft_error(values, length, order):
    """Two radix FFTs, derivative multiplier and normalization, in max norm.

    Each radix-2 dependency path costs at most eight real operations per
    complex butterfly; 32*log2(n)+16 charges both transforms, twiddle error,
    multiplier evaluation and normalization. Fourier triangle inequality
    gives sum(|k|**order)*max|values|, independently of a cancelled result.
    """
    n = len(values)
    frequencies = 2 * np.pi * np.fft.fftfreq(n, length / n)
    frequencies[n // 2] = 0
    mass = float(np.sum(abs(frequencies) ** order)) * float(np.max(abs(values)))
    return float(np.nextafter(_gamma(32 * ceil(log2(n)) + 16) * mass, np.inf))


@lru_cache(maxsize=3)
def euler_replay_roundoff_envelope(points):
    """Pairwise arithmetic envelope for the original 32/64/128 smooth cells.

    Input construction has fewer than 64 scalar operations per node,
    including trig/exp/sqrt error. Derivative input errors propagate through
    the rational jet coefficients using intervals. Both coefficient
    evaluations and both subsequent FFT differentiations are charged. The
    final absolute maximum is 1-Lipschitz in the componentwise max norm.
    No stored residual, observed drift or physical tolerance sets this bound.
    """
    if isinstance(points, bool) or points not in (32, 64, 128):
        raise ValueError("Euler portability is restricted to the declared three grids")
    metric = smooth_metric(points, general=True)
    # Exact first-derivative Fourier collocation matrix has off-diagonal
    # (pi/L)*(-1)**j*cot(pi*j/n); the retained second derivative is its square.
    angle = np.pi * np.arange(1, points) / points
    first_norm = (np.pi / metric.length) * float(np.sum(abs(1 / np.tan(angle))))
    first_norm *= 1 + _gamma(8 * points + 32)
    second_norm = first_norm * first_norm
    jets, mass, _, _, _, _, densities = _symbolic_basis()
    values, boxes = [], {}
    for symbols, field in zip(jets, (metric.lapse, metric.radial_scale, metric.sphere_radius)):
        raw_pair_error = 2 * _gamma(64) * float(np.max(abs(field)))
        for order, symbol in enumerate(symbols):
            value = field if order == 0 else periodic_derivative(field, metric.length, order)
            operator_norm = (1., first_norm, second_norm)[order]
            pair_error = operator_norm * raw_pair_error
            if order:
                pair_error += 2 * _fft_error(field, metric.length, order)
            values.append(value)
            boxes[symbol] = _outward(value - pair_error, value + pair_error)
    boxes[mass] = (np.asarray(1.), np.asarray(1.))
    evaluated = _density_functions()[0][4](*values, 1.)
    envelopes = []
    for field_index, symbols in enumerate(jets):
        coefficients, errors = [], []
        for order, symbol in enumerate(symbols):
            expression = sp.diff(densities[4], symbol)
            coefficient = np.broadcast_to(evaluated[1 + 3 * field_index + order], (points,))
            lower, upper = _interval(expression, boxes)
            operation_budget = 4 * int(sp.count_ops(expression)) + 32
            error = np.maximum(abs(lower - coefficient), abs(upper - coefficient))
            error += 2 * _gamma(operation_budget) * _majorant(expression, boxes)
            coefficients.append(coefficient)
            errors.append(float(np.max(error)))
        p0, p1, p2 = coefficients
        uncancelled = (float(np.max(abs(p0))) + first_norm * float(np.max(abs(p1)))
                       + second_norm * float(np.max(abs(p2))))
        error = errors[0] + first_norm * errors[1] + second_norm * errors[2]
        error += 2 * (_fft_error(p1, metric.length, 1) + _fft_error(p2, metric.length, 2))
        error += 2 * _gamma(8) * uncancelled
        envelopes.append(4 * np.pi * error)
    bound = float(np.nextafter(max(envelopes), np.inf))
    if not np.isfinite(bound) or bound <= 0:
        raise ValueError("nonfinite Euler arithmetic envelope")
    return bound


@lru_cache(maxsize=3)
def euler_summary_consistency_envelope(points):
    """Return Euler scale and rounding in the independent residual identity.

    The full and C2-only vectors share the same determinant source. Their
    exact difference is dx*delta_E4*G_E4 because the boxR vector is zero.
    The reverse triangle inequality then connects their reported max norms.
    gamma_32 charges both six-term contractions, spacing and source additions
    on uncancelled absolute operands. This is arithmetic, not a physics bound.
    """
    if isinstance(points, bool) or points not in (32, 64, 128):
        raise ValueError("Euler consistency is restricted to the declared three grids")
    from .nsc_normalization_flow import (
        axial_metric, c2_only_coefficient_shift, determinant_source_change,
        finite_coefficient_shift,
    )
    metric = smooth_metric(points, general=True)
    response = finite_response(axial_metric(metric))
    source = determinant_source_change(metric, 2., 1., 2.)["metric_gradients"]
    full, c2 = finite_coefficient_shift(2.), c2_only_coefficient_shift(2.)
    operands = abs(response["gradients"])
    full_mass = metric.spacing * np.einsum("k,kij->ij", abs(full), operands)
    c2_mass = metric.spacing * np.einsum("k,kij->ij", abs(c2), operands)
    uncancelled = float(np.max(2 * abs(source) + full_mass + c2_mass))
    arithmetic = float(np.nextafter(2 * _gamma(32) * uncancelled, np.inf))
    return abs(metric.spacing * full[4]), arithmetic


def validate_euler_summary(row):
    """Reject a diagnostic inconsistent with its unchanged physical summaries."""
    points = row["points"]
    scale, arithmetic = euler_summary_consistency_envelope(points)
    values = [row[key] for key in ("discrete_E4_el_maximum",
              "c2_only_gradient_residual", "full_gradient_residual")]
    if any(type(value) is not float or not np.isfinite(value) or value < 0 for value in values):
        raise ValueError("invalid Euler diagnostic or residual summary")
    diagnostic, c2_residual, full_residual = values
    if abs(scale * diagnostic - c2_residual) > full_residual + arithmetic:
        raise ValueError("Euler diagnostic is inconsistent with its residual summaries")
