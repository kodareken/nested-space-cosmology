"""Endpoint comparison for a quiet monotone source and certified radius image.

Source positivity/monotonicity and the exact reference image are caller-proved
assumptions. Derivatives use the original affine source label. This component
does not certify a radius or clock coefficient or evaluate a sampled field.
"""
from fractions import Fraction as Q

from .nsc_atlas_transmission_bound import _rational

A = Q(4, 5)


def prefix_kernel_integral(x, radius):
    """Exact transverse kernel integral at the fixed reference scale 4/5."""
    x = _rational(x, 'x', positive=True)
    radius = _rational(radius, 'radius', positive=True)
    if not x < radius < A:
        raise ValueError('require 0 < x < radius < 4/5')
    return A*x/(2*(A-x))*(1/x**2-1/radius**2)


def _ordered_pair(value, name):
    if not isinstance(value, tuple) or len(value) != 2:
        raise ValueError(name+' requires an ordered pair of nonnegative rationals')
    lower = _rational(value[0], name+' lower')
    upper = _rational(value[1], name+' upper')
    if lower > upper:
        raise ValueError(name+' bounds are not ordered')
    return lower, upper


def monotone_wall_enclosure(*, left, right, radius_lower, radius_upper, g, g_x, budget):
    """Propagate proved endpoint source bounds through the positive Volterra bound.

    The source must be nonnegative and nondecreasing over the entire interval,
    quiet at its left edge, with endpoint value/derivative in g/g_x. The radius
    image must enclose the exact reference wall throughout the interval.
    Checking these scalar inputs alone does not establish those hypotheses.
    """
    if not callable(getattr(budget, 'check', None)):
        raise ValueError('require a shared cumulative budget')
    budget.check(512)
    left = _rational(left, 'left', positive=True)
    right = _rational(right, 'right', positive=True)
    radius_lower = _rational(radius_lower, 'radius lower', positive=True)
    radius_upper = _rational(radius_upper, 'radius upper', positive=True)
    if not left < right < radius_lower <= radius_upper < A:
        raise ValueError('require 0 < left < right < radius lower <= radius upper < 4/5')
    g_lo, g_hi = _ordered_pair(g, 'g')
    gx_lo, gx_hi = _ordered_pair(g_x, 'g_x')
    b_upper = prefix_kernel_integral(left, radius_upper)
    b_lo = prefix_kernel_integral(right, radius_lower)
    b_hi = prefix_kernel_integral(right, radius_upper)
    duration = right-left
    q = duration*b_upper
    if q >= 1:
        raise ValueError('clipped Volterra norm is not strictly below one')
    margin = 1-q
    budget.check()
    return dict(
        duration=duration, B_upper=b_upper, B_endpoint=dict(lower=b_lo, upper=b_hi), q=q,
        wall_value=dict(lower=margin*g_lo, upper=g_hi),
        wall_derivative=dict(lower=gx_lo-b_hi*g_hi, upper=gx_hi-margin*b_lo*g_lo),
        source_minus_wall=dict(lower=Q(0), upper=q*g_hi),
        source_derivative_minus_wall_derivative=dict(lower=margin*b_lo*g_lo, upper=b_hi*g_hi),
        scope='leading scalar and affine-label derivative from proved monotone source and reference image; '
              'no finite-amplitude radius or clock coefficient')
