"""Directed defect of the homogeneous reference in an anchored joint segment.

Reuses the existing Bernstein and background Taylor owners. Anchored endpoint
differences and polynomial conversion are formed in ball arithmetic, not by
converting the DOP853 coefficients to binary64 powers first.
"""
from math import comb

import numpy as np
from flint import acb, arb, ctx

from .nsc_ks_ball_trajectory import complex_ball, exact_upper
from .nsc_ks_trajectory import TrajectorySegment
from .nsc_ks_upstream_error import (
    _background_power_polynomials, _poly_add, _poly_convolution,
    _vector_bernstein_norm,
)


def reference_segment_defect(segment, energies, weights, mass, angular, *, bits=120):
    """Bound integral of weighted reference defect over one decreasing cell.

    Every source column has its own declared energy and original weight.
    The returned weighted Frobenius bound also bounds the max-spin-row norm
    used by propagate_difference_error. This is numerical field error only.
    """
    if not isinstance(segment, TrajectorySegment):
        raise TypeError('anchored joint trajectory segment required')
    energies, weights = np.asarray(energies, float), np.asarray(weights, float)
    if (energies.ndim != 1 or not len(energies) or weights.shape != energies.shape
            or not np.isfinite(energies).all() or not np.isfinite(weights).all()
            or np.any(weights <= 0) or not np.isfinite([mass, angular]).all()):
        raise ValueError('finite labels, channel and original positive weights required')
    n = len(energies)
    if len(segment.start) < 2*n:
        raise ValueError('reference columns missing from joint segment')
    with ctx.workprec(bits):
        # Background owner uses increasing u; joint trajectory uses x=1-u.
        polys, errors, width = _background_power_polynomials(
            segment.rho_end, segment.rho_start, 8, bits=bits)
        geometry = {name: [sum((c[j]*comb(j,k)*(-1)**k
                                     for j in range(k,len(c))), arb(0))
                                for k in range(len(c))]
                    for name,c in polys.items()}
        h = arb(segment.rho_end)-arb(segment.rho_start)
        columns = []
        for source in range(n):
            ys = []
            for spin in range(2):
                index = spin*n+source
                c = [acb(0) for _ in range(8)]
                c[0] = complex_ball(segment.start[index])
                c[1] = complex_ball(segment.end[index])-c[0]
                for j in range(1,7):
                    p,q = (j+2)//2,(j+1)//2
                    v = complex_ball(segment.coefficients[j-1,index])
                    for k in range(q+1):c[p+k] += (-1)**k*comb(q,k)*v
                ys.append([v*arb(float(weights[source])) for v in c])
            y0,y1 = ys
            r0,r1 = [[(j+1)*c[j+1]/h for j in range(7)] for c in ys]
            energy = float(energies[source])
            _poly_add(r0,_poly_convolution(geometry['inv_a2'],y0),acb(0,energy))
            _poly_add(r1,_poly_convolution(geometry['inv_a2'],y1),-acb(0,energy))
            _poly_add(r0,_poly_convolution(geometry['inv_a'],y1),acb(0,float(mass)))
            _poly_add(r1,_poly_convolution(geometry['inv_a'],y0),acb(0,float(mass)))
            _poly_add(r0,_poly_convolution(geometry['inv_ar'],y1),-float(angular))
            _poly_add(r1,_poly_convolution(geometry['inv_ar'],y0),float(angular))
            tail = (abs(arb(energy))*errors['inv_a2']
                    +abs(arb(float(mass)))*errors['inv_a']
                    +abs(arb(float(angular)))*errors['inv_ar'])
            columns.append((_vector_bernstein_norm(r0,r1)
                            +tail*_vector_bernstein_norm(y0,y1)).upper())
        supremum = sum((v*v for v in columns),arb(0)).sqrt().upper()
        return {'weighted_frobenius_defect_sup': exact_upper(supremum),
                'weighted_row_residual_integral': exact_upper(width*supremum),
                'numerical_reference_only': True,
                'physical_source_error_included': False}
