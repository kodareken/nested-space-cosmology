"""Local exact-C2 metric bound from original Hermite and departure data.

Arb bounds the mathematical Hermite interpolant of the original binary
inputs. Numerical PPoly conversion, field evaluation and ODE errors are
excluded. Existing profile knots define the partition; no adaptive scan.
"""
import numpy as np
from flint import arb, ctx

from .nsc_atlas_background_cap import _budget
from .nsc_atlas_cap_background_jets import RetainedCapBackground
from .nsc_atlas_cap_polynomial_basis import _count
from .nsc_atlas_cap_variational_reference import _exact_ball
from .nsc_atlas_polar_cap import _scalar, _vector
from .nsc_atlas_characteristic_background import _onset_bernstein


def _hermite_controls(v0,d0,s0,v1,d1,s1,h):
    """Exact quintic Bernstein controls, then derivative controls."""
    values=[v0,v0+h*d0/5,v0+2*h*d0/5+h*h*s0/20,
            v1-2*h*d1/5+h*h*s1/20,v1-h*d1/5,v1]
    first=[5*(b-a)/h for a,b in zip(values[:-1],values[1:])]
    second=[4*(b-a)/h for a,b in zip(first[:-1],first[1:])]
    return values,first,second


def _background_controls(background,knots,values,first,second):
    """Original-input controls, including the one-sided active onset piece."""
    controls=[_hermite_controls(values[i],first[i],second[i],values[i+1],first[i+1],second[i+1],
                                knots[i+1]-knots[i]) for i in range(len(knots)-1)]
    onset=background.onset
    if onset is not None:
        i=onset['index'];h=knots[i+1]-knots[i]
        v=list(map(arb,_onset_bernstein(_exact_ball(onset['quartic_coefficient']),
                                      values[i+1],first[i+1],second[i+1],h)))
        d=[7*(b-a)/h for a,b in zip(v[:-1],v[1:])]
        s=[6*(b-a)/h for a,b in zip(d[:-1],d[1:])]
        controls[i]=(v,d,s)
    return controls


def _upper_max(values):
    return max(value.upper() for value in values)


def _envelope_step(U,b,left,right):
    """Outward upper envelope and its largest n/r ratio on one interval."""
    following=(U+b*(right-left).upper()).upper()
    ratio=b.upper() if left==0 else _upper_max([U/left.lower(),following/right.lower()])
    return following,ratio


def _piece_ranges(knots,controls,lower,upper,budget):
    # All intersected pieces contribute to a maximum, not a sum. The whole
    # control hull is a conservative enclosure even for a point overlap.
    result=[arb(0),arb(0),arb(0)]
    for i,piece in enumerate(controls):
        budget.check(4)
        if upper < knots[i] or lower > knots[i+1]:
            continue
        for order,coefficients in enumerate(piece):
            result[order]=_upper_max([result[order],*(abs(c) for c in coefficients)])
            budget.check(len(coefficients))
    return result


def background_metric_bound(background,T,Z,chi_coefficients,density_coefficients,*,alpha,budget):
    """Sufficient all-radius A>0 for one exact-C2 background plus departures."""
    budget=_budget(budget)
    before=budget.evaluations
    if not isinstance(background,RetainedCapBackground):
        raise ValueError('require original retained Hermite data')
    T,Z,alpha=_scalar(T,'T'),_scalar(Z,'Z',positive=True),_scalar(alpha,'alpha',nonnegative=True)
    chi=_vector(chi_coefficients,'chi_coefficients')
    k=_vector(density_coefficients,'density_coefficients')
    count=_count(len(chi))
    if k.shape != chi.shape:
        raise ValueError('matching departure coefficient vectors required')
    data=background.inputs
    with ctx.workprec(160):
        tb,zb,al=map(_exact_ball,(T,Z,alpha))
        knots=[_exact_ball(x) for x in data['labels']]
        if not tb+zb<=knots[-1] or (not background.profile.quiet_extension and not tb-zb>=knots[0]):
            raise ValueError('exact slice exits the supplied Hermite domain; an inward slice must be explicit')
        R=zb.tanh()
        if not R>0 or not R<1:
            raise ValueError('metric interval is not resolved inside (0,1)')
        values,first,second=([_exact_ball(x) for x in data[key]] for key in ('values','first','second'))
        controls=_background_controls(background,knots,values,first,second)
        budget.check(96*len(controls)+64*count)
        cb,kb=([abs(_exact_ball(x)) for x in a] for a in (chi,k))
        slopes=[arb(p*(p+1)//2) for p in range(1,2*count,2)]
        lc=sum(c*s for c,s in zip(cb,slopes))/zb
        lk=sum(c*s for c,s in zip(kb,slopes))/zb
        mapped=np.r_[data['labels']-T,T-data['labels']]
        cuts=np.unique(np.r_[0.,mapped[(mapped>0)&(mapped<Z)],Z])
        U=arb(0)
        floors=[]
        for left,right in zip(cuts[:-1],cuts[1:]):
            lo,hi=_exact_ball(left),_exact_ball(right)
            rl,rh=lo.tanh(),hi.tanh()
            Ahi,Alo=1/hi.cosh()**2,1/lo.cosh()**2
            plus=_piece_ranges(knots,controls,tb+lo,tb+hi,budget)
            minus=_piece_ranges(knots,controls,tb-hi,tb-lo,budget)
            if lo==0:
                Lchi=(plus[1]+minus[1]+lc).upper()
                Lk=(plus[2]+minus[2]+lk).upper()
                fmax=(Lchi*hi/rh.lower()).upper()
                pmax=(Lk*hi/Ahi.lower()).upper()
                hmax=(Lchi/Ahi.lower()+fmax).upper()
            else:
                factors=[min(arb(1),(s*hi/zb).upper()) for s in slopes]
                chib=(plus[0]+minus[0]+sum(c*b for c,b in zip(cb,factors))).upper()
                kbnd=(plus[1]+minus[1]+sum(c*b for c,b in zip(kb,factors))).upper()
                qnum=(plus[1]+minus[1]+lc).upper()
                fmax=(chib/rl.lower()).upper()
                pmax=(kbnd/Ahi.lower()).upper()
                hmax=(qnum/Ahi.lower()+fmax).upper()
            q=(al/6).sqrt()*fmax
            if not q<=100:
                raise ValueError('potential bound leaves the finite reference domain')
            potential=fmax*fmax if al==0 or fmax==0 else fmax*fmax*(q.sinh()/q)**2*(q.cosh()**2+1)/2
            source=(rh*rh*potential/2+Alo*(pmax*pmax+hmax*hmax)/4).upper()
            U,ratio=_envelope_step(U,source,rl,rh)
            floor=(Ahi-2*al*ratio).lower()
            budget.check(160+24*count)
            if not floor.is_finite() or not floor>0:
                raise ValueError('sufficient local metric bound is not positive; physical metric sign is undetermined')
            floors.append(floor)
        lower=min(floors)
        mantissa,exponent=lower.man_exp()
        down=float(np.nextafter(float(lower),-np.inf))
        mass_upper=0. if U==0 else float(np.nextafter(float(U.upper()),np.inf))
        if not np.isfinite(down) or down<=0 or not np.isfinite(mass_upper):
            raise ValueError('metric certificate is not representable')
        result=dict(A_lower=dict(mantissa=str(mantissa),exponent=int(exponent)),
            A_lower_float_down=down,precision_bits=ctx.prec,interval_count=len(cuts)-1,
            mass_upper=mass_upper,evaluated_points=budget.evaluations-before,
            scope='exact C2 Hermite interpolant of original binary inputs plus odd Legendre departures; numerical interpolation, evaluation and ODE errors excluded')
    budget.check()
    return result
