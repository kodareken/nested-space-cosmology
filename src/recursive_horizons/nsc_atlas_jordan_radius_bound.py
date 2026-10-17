"""Finite-onset material-radius conversion from the canonical scalar."""
from fractions import Fraction as Q

from flint import ctx

from .nsc_atlas_geometric_departure_bound import _round_upper
from .nsc_atlas_geometric_remainder_bound import endpoint_comparison
from .nsc_atlas_onset_geometry_bound import dyadic
from .nsc_atlas_reference_gauge_bound import _ball


def interval_product(a,b):
    products=[x*y for x in a for y in b]
    return min(products),max(products)


def first_jordan_interval(first):
    def pair(row):
        lo,hi=(Q(row[k]) for k in ('lower','upper'))
        if lo>hi:raise ValueError('ordered first-response interval required')
        return lo,hi
    factor=(Q(2),Q(2))
    for name in ('a','v','d'):
        value=pair(first['constants'][name])
        if value[0]<=0:raise ValueError('positive reference event-shift factor required')
        factor=interval_product(factor,value)
    shift=interval_product(factor,pair(first['outputs']['eta']))
    rho=pair(first['outputs']['rho_proper_time'])
    return dict(lower=rho[0]+shift[0],upper=rho[1]+shift[1])


def frame_half_hessian(*,radius_lower,radius_upper,X_upper,radius_difference,
                       X_difference,alpha_max,budget):
    """Positive analytic bound for (R*cosh(sqrt(alpha/6)*X/R))''/2.

    R and X are affine along the comparison segment, not along the actual
    solution family. Quotient derivatives of f=X/R are retained.
    """
    budget.check(4096)
    values=(radius_lower,radius_upper,X_upper,radius_difference,X_difference,alpha_max)
    if any(isinstance(v,(float,bool)) for v in values):raise ValueError('exact frame inputs required')
    r,h,x,cr,cx,al=map(Q,values)
    if not 0<r<=h or any(v<0 for v in (x,cr,cx,al)):
        raise ValueError('positive radius and nonnegative frame bounds required')
    f=x/r
    d1=cx/r+x*cr/r**2
    d2=2*x*cr**2/r**3+2*cx*cr/r**2
    q1=f*f/6+al*f*d1/3
    q2=2*f*d1/3+al*(d1*d1+f*d2)/3
    with ctx.workprec(160):
        ch=_round_upper((_ball(al/6).sqrt()*_ball(f)).cosh())
        raw=cr*ch*q1/2+h*ch*q1*q1/24+h*ch*q2/4
        bound=_round_upper(_ball(raw))
    budget.check()
    return dict(f_upper=f,f_first=d1,f_second=d2,q_first=q1,q_second=q2,
                cosh_upper=ch,half_hessian=bound,unrounded_positive_bound=raw,
                conformal_factor_lower=1/(ch*ch))


def jordan_radius_bound(tube,geometry,scalar,remainder,first,*,budget):
    budget.check(2048)
    al=Q(tube['alpha_max'])
    if al<=0:raise ValueError('positive amplitude interval required')
    for record in (geometry,scalar,remainder):
        if Q(record['alpha_max'])!=al or Q(record['duration'])!=Q(tube['duration']):
            raise ValueError('common original amplitude and cell required')
    if Q(scalar['z_star'])!=Q(tube['z_star']):raise ValueError('common collar required')
    def ceiling(name):
        row=tube['constants'][name]['upper']
        return Q(row['mantissa'])*Q(2)**row['exponent']
    # These stored floor constants are point dyadics from the strict tube.
    for name in ('wall_radius_floor','radius_upper','X_upper'):
        row=tube['constants'][name]
        if row['lower']!=row['upper']:raise ValueError('point-valued original tube constant required')
    inputs=dict(radius_lower=ceiling('wall_radius_floor'),radius_upper=ceiling('radius_upper'),
        X_upper=ceiling('X_upper'),radius_difference=Q(geometry['coefficients']['radius']),
        X_difference=Q(scalar['coefficients']['X']),alpha_max=al)
    frame=frame_half_hessian(**inputs,budget=budget)
    first_jordan=first_jordan_interval(first)
    ce=Q(remainder['coefficients']['radius'])
    if ce<0:raise ValueError('nonnegative Einstein remainder required')
    coefficients=dict(radius=ce+frame['half_hessian'],delta_T=Q(remainder['coefficients']['delta_T']))
    first_values=dict(radius_fixed_source=first_jordan,clock_fixed_source=first['outputs']['clock_fixed_source'])
    endpoint=endpoint_comparison(first_values,coefficients,al,budget=budget)
    budget.check(1024)
    return dict(alpha_max=al,duration=Q(tube['duration']),frame_inputs=inputs,frame=frame,
        first_jordan=first_jordan,einstein_radius_remainder=ce,coefficients=coefficients,
        maximum_amplitude_remainders={k:al*al*v for k,v in coefficients.items()},
        endpoint=endpoint,common_strict_sign_alpha_max=min(e['strict_sign_alpha_max'] for e in endpoint.values()),
        strict_endpoint_sign_on_full_interval=all(e['strict_on_full_interval'] for e in endpoint.values()),
        scope='Jordan areal radius at fixed source label on original onset; canonical f=X/R_E; '
              'same Einstein-polar slicing clock anchored to physical central proper time in quiet core; '
              'not wall proper time, later arrival, source selection or observational identification')
