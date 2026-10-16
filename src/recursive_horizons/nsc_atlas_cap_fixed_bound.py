"""Local metric certificate for exact fixed shapes and binary coefficients.

Bernstein hulls enclose each existing physical polynomial piece. The radial
source envelope bounds the exact constrained solution, not its numerical
evaluation or a time evolution. No subdivision or profile query is used.
"""
from fractions import Fraction as Q
from math import comb

import numpy as np
from flint import arb,ctx,fmpq

from .nsc_atlas_background_cap import _budget
from .nsc_atlas_cap_background_bound import _envelope_step
from .nsc_atlas_cap_fixed_basis import FixedCapBasis,exact_number
from .nsc_atlas_polar_cap import _scalar,_vector


def _ball(value):
    value=Q(value)
    return arb(fmpq(value.numerator,value.denominator))


def _bernstein_controls(power,cut=Q(1)):
    """Exact controls on xi in [0,cut], at the declared polynomial degree."""
    n=len(power)-1
    scaled=tuple(Q(c)*cut**i for i,c in enumerate(power))
    return tuple(sum(scaled[k]*Q(comb(j,k),comb(n,k)) for k in range(j+1))
                 for j in range(n+1))


def _hull(power,cut):
    return max(abs(_ball(v)).upper() for v in _bernstein_controls(power,cut))


def fixed_metric_bound(basis,chi,canonical_density,Z,*,alpha,budget):
    """Sufficient all-radius positivity for one exact fixed-space field pair."""
    budget=_budget(budget);before=budget.evaluations
    if not isinstance(basis,FixedCapBasis):
        raise ValueError('fixed metric certificate requires a fixed local basis')
    chi=_vector(chi,'chi');density=_vector(canonical_density,'canonical density')
    if chi.shape!=(basis.dimension,) or density.shape!=chi.shape:
        raise ValueError('fixed metric coefficient dimensions differ')
    Z=exact_number(Z,'wall coordinate')
    cut=(Z-basis.exact_knots[-2])/(basis.exact_knots[-1]-basis.exact_knots[-2])
    if not Q(1,4)<=cut<=Q(3,4):
        raise ValueError('fixed metric wall cut must lie in [1/4,3/4]')
    alpha=_scalar(alpha,'alpha',nonnegative=True)
    exact=(tuple(Q(float(v)) for v in chi),tuple(Q(float(v)) for v in density))
    budget.check(64*basis.dimension)
    with ctx.workprec(160):
        zb,al=_ball(Z),_ball(alpha)
        if not 0<zb.tanh()<1:
            raise ValueError('metric interval is not resolved inside (0,1)')
        U=arb(0);floors=[]
        for cell in range(basis.cells):
            lo,hi=basis.exact_knots[cell:cell+2];h=hi-lo
            local_cut=cut if cell==basis.cells-1 else Q(1)
            hi=min(hi,Z)
            idx=basis.indices(cell);shapes=basis.exact_shapes(cell)
            pieces=[tuple(sum(row[j]*coeff[int(i)] for j,i in enumerate(idx)) for row in shapes)
                    for coeff in exact]
            first=[tuple(p*poly[p]/h for p in range(1,len(poly))) for poly in pieces]
            lc=_hull(first[0],local_cut)
            lob,hib=_ball(lo),_ball(hi)
            rl,rh=lob.tanh(),hib.tanh()
            Ahi,Alo=1/hib.cosh()**2,1/lob.cosh()**2
            if cell==0:
                lk=_hull(first[1],local_cut)
                fmax=(lc*hib/rh.lower()).upper()
                pmax=(lk*hib/Ahi.lower()).upper()
            else:
                fmax=(_hull(pieces[0],local_cut)/rl.lower()).upper()
                pmax=(_hull(pieces[1],local_cut)/Ahi.lower()).upper()
            hmax=(lc/Ahi.lower()+fmax).upper()
            argument=(al/6).sqrt()*fmax
            if not argument<=100:
                raise ValueError('fixed metric potential bound leaves its finite domain')
            potential=(fmax*fmax if al==0 or fmax==0 else
                       fmax*fmax*(argument.sinh()/argument)**2*(argument.cosh()**2+1)/2)
            source=(rh*rh*potential/2+Alo*(pmax*pmax+hmax*hmax)/4).upper()
            U,ratio=_envelope_step(U,source,rl,rh)
            floor=(Ahi-2*al*ratio).lower()
            budget.check(1024)
            if not floor.is_finite() or not floor>0:
                raise ValueError('sufficient fixed metric bound is not positive; physical sign is undetermined')
            floors.append(floor)
        lower=min(floors)
        mantissa,exponent=lower.man_exp()
        downward=float(np.nextafter(float(lower),-np.inf))
        mass_upper=0. if U==0 else float(np.nextafter(float(U.upper()),np.inf))
        if not np.isfinite(downward) or downward<=0 or not np.isfinite(mass_upper):
            raise ValueError('fixed metric certificate is not representable')
        result=dict(A_lower=dict(mantissa=str(mantissa),exponent=int(exponent)),
            A_lower_float_down=downward,mass_upper=mass_upper,precision_bits=ctx.prec,
            interval_count=basis.cells,wall_cut_fraction=float(cut),
            coefficient_evaluations=budget.evaluations-before,
            scope='exact rational fixed shapes with binary coefficient inputs; floating field evaluation, radial ODE and evolution errors excluded')
    budget.check()
    return result
