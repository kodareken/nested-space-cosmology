"""Certified Fourier coefficients of the compact axial potential profiles.

Arb DFT arithmetic is combined with continuous derivative L1 bounds. The
classical alias identity and integration by parts bound coefficients outside
the sampled/retained sets. No empirical Fourier decay fit is used.
"""
from dataclasses import dataclass
from math import factorial

from flint import acb, arb, ctx

from . import nsc_ks_ball_geometry as G


@dataclass(frozen=True)
class ProfileFourierBound:
    key: tuple
    coefficients: object
    derivative_l1_upper: object
    alias_error_upper: object
    uniform_tail_bounds: object
    period_origin: object
    period_length: object
    derivative_order: int
    quadrature_points: int
    retained_index: int
    derivative_cells: int


def alias_and_tail_bounds(derivative_l1, length, derivative_order, quadrature_points,
                          retained_index, max_derivative=2):
    """Bounds with Fourier convention c_k=L^-1 integral f exp(-2 pi i k z/L)."""
    p, M, K = derivative_order, quadrature_points, retained_index
    if any(isinstance(v, bool) or not isinstance(v, int) for v in (p, M, K, max_derivative)):
        raise ValueError('integer derivative and Fourier counts required')
    if not 0 < K < M/2 or p <= max_derivative+1:
        raise ValueError('retained band below Nyquist and integrable derivative tails required')
    L, B = arb(length), arb(derivative_l1)
    if not L > 0 or not B >= 0:
        raise ValueError('positive period and nonnegative derivative bound required')
    scale = B/L*(L/(2*arb.pi()))**p
    alias = (2*scale*p/((p-1)*arb(M-K)**p)).upper()
    tails = tuple((2*scale*(2*arb.pi()/L)**j*arb(K)**(j-p+1)/(p-j-1)).upper()
                  for j in range(max_derivative+1))
    return alias, tails


def _common_window(model):
    profiles = [profile for d in model.metric.directions for profile in (d.w, d.U)]
    windows = {(float(p.center), float(p.inner), float(p.outer)) for p in profiles}
    if len(windows) != 1:
        raise ValueError('the approved common axial cutoff is required')
    return next(iter(windows))


def _derivative_l1(model, keys, order, transition_panels):
    center, inner, outer = map(arb, _common_window(model))
    sums = {key: arb(0) for key in keys}
    intervals = []
    for left, right, count in ((-outer,-inner,transition_panels), (-inner,inner,4),
                                (inner,outer,transition_panels)):
        for i in range(count):
            intervals.append((center+left+(right-left)*i/count,
                              center+left+(right-left)*(i+1)/count, 0))
    cells = 0
    while intervals:
        left, right, depth = intervals.pop()
        try:
            with G._JetWork(order):
                w, U = model.profile_series(left.union(right),order)
                values = {key: G._coeffs(w**key[0]*U**key[1],order)[order]*factorial(order)
                          for key in keys}
        except G.SubdivisionNeeded:
            if depth >= 12:
                raise
            middle = (left+right)/2
            intervals.extend(((left,middle,depth+1),(middle,right,depth+1)))
            continue
        for key, value in values.items():
            sums[key] += (right-left)*abs(value).upper()
        cells += 1
    return {key: value.upper() for key,value in sums.items()}, cells


def enclose_profile_fourier(model, keys, period_origin, period_length, *,
                            derivative_order=8, transition_panels=256,
                            quadrature_points=16384, retained_index=2048,
                            max_derivative=2, bits=120):
    """Return retained coefficient balls and continuous tails for w^p U^q."""
    if (not keys or len(set(keys)) != len(keys) or
            any(len(k)!=2 or min(k)<0 or sum(k)<1 for k in keys)):
        raise ValueError('distinct positive-total-degree profile powers required')
    if transition_panels < 1 or not isinstance(transition_panels,int):
        raise ValueError('positive transition panel count required')
    with ctx.workprec(bits):
        origin, length = arb(period_origin), arb(period_length)
        center, inner, outer = map(arb,_common_window(model))
        if not origin < center-outer or not origin+length > center+outer:
            raise ValueError('compact support must lie strictly inside the numerical period')
        l1, cells = _derivative_l1(model,keys,derivative_order,transition_panels)
        errors = {key:alias_and_tail_bounds(l1[key],length,derivative_order,
                    quadrature_points,retained_index,max_derivative) for key in keys}
        samples = {key:[] for key in keys}
        for node in range(quadrature_points):
            z = origin+length*node/quadrature_points
            if z <= center-outer or z >= center+outer:
                values = {key:arb(0) for key in keys}
            else:
                w, U = model.profile_series(z,0)
                values = {key:w[0]**key[0]*U[0]**key[1] for key in keys}
            for key in keys:
                samples[key].append(acb(values[key]))
        results = {}
        for key in keys:
            transformed = acb.dft(samples[key])
            alias, tails = errors[key]
            error = arb((0,0),alias.man_exp())
            coefficients = tuple(transformed[k % quadrature_points]/quadrature_points+acb(error,error)
                                 for k in range(-retained_index,retained_index+1))
            results[key] = ProfileFourierBound(key,coefficients,l1[key],alias,tails,
                origin,length,derivative_order,quadrature_points,retained_index,cells)
        return results
