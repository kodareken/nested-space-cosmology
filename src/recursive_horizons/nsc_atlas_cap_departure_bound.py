"""Radial cap differences on the already certified common z collar.

The quiet-core mass and lapse integrals are compared algebraically.
N=(m-r^3/2)/alpha and D=delta/alpha have nonzero leading profiles.
No radial solution, field sample, or full geometric remainder is computed.
"""
from fractions import Fraction as Q

from flint import ctx

from .nsc_atlas_coupled_tube_bound import potential_majorants
from .nsc_atlas_geometric_departure_bound import _stored, _u, _l, _round_upper
from .nsc_atlas_onset_geometry_bound import dyadic
from .nsc_atlas_reference_gauge_bound import _ball


def potential_remainder_coefficient(field, alpha):
    """Upper coefficient of alpha in Vcap(alpha,f)-f^2, including alpha=0."""
    if not field >= 0 or not alpha >= 0:
        raise ValueError('require nonnegative field and coupling envelopes')
    q = (alpha/6).sqrt()*field
    return _u(field**4*(4*(4*q).cosh()+(2*q).cosh())/36)


def radial_comparison(values, *, budget):
    """Exact positive product comparison, with caller-certified envelopes."""
    budget.check(512)
    if any(isinstance(v,(bool,float)) for v in values.values()):
        raise ValueError('require exact rational cap inputs')
    v = {k:Q(x) for k,x in values.items()}
    if any(x<0 for x in v.values()):
        raise ValueError('require nonnegative cap envelopes')
    r, R, A, Amax, L = (v[k] for k in ('r_min','r_max','A0_min','A0_max','length_max'))
    if not 0<r<=R<1 or not 0<A<=Amax<=1 or L<=0:
        raise ValueError('require positive common cap floors and finite collar')
    f = v['chi_max']/r
    p = (v['E_plus_max']+v['E_minus_max'])/(2*A)
    h = p+f
    if v['potential_max']<f*f or v['potential_remainder']<Q(5,36)*f**4:
        raise ValueError('potential envelopes omit the zero-alpha limit')
    cr = Amax*v['C_Z']
    ca = 2*R*cr
    cf = v['C_chi']/r+v['chi_max']*cr/r**2
    cp = (v['C_E_plus']+v['C_E_minus'])/(2*A)+p*ca/A
    ch = cp+cf
    qc = (p*p+h*h)/r**2
    cqc = 2*(p*cp+h*ch)/r**2+2*(p*p+h*h)*cr/r**3
    cv = v['potential_remainder']+2*f*cf
    jac = L*Amax
    cj = Amax*v['C_Z']+L*ca
    b = R**2*v['potential_max']/2+R**2*Amax*qc/4
    cb = (2*R*cr*v['potential_max']+R**2*cv)/2
    cb += ((2*R*cr*Amax+R**2*ca)*qc+R**2*Amax*cqc)/4
    attenuation = jac*R*qc*v['N_max']/2
    cn = cj*b+jac*cb+attenuation
    cd = (cj*R*qc+jac*(cr*qc+R*cqc))/2
    # These are second-order residuals after removing the vacuum part at
    # the actual radius, not the complete geometric Taylor remainder.
    mu2 = cd+v['mu_max']*v['D_max']**2/2
    a2 = 2*(cn/r+v['N_max']*cr/r**2)
    budget.check()
    return dict(C_r=cr,C_A0=ca,C_f=cf,C_p=cp,C_h=ch,Qc_max=qc,C_Qc=cqc,
        C_potential=cv,jacobian_max=jac,C_jacobian=cj,density_max=b,
        C_density=cb,attenuation_contribution=attenuation,N=cn,D=cd,
        mu_second_order=mu2,A_vacuum_subtracted_second_order=a2)


def cap_departure_bound(tube, scalar, *, budget):
    if not callable(getattr(budget,'check',None)):
        raise ValueError('require a shared cumulative budget')
    budget.check(4096)
    if not tube['strict'] or not all(tube['component_checks'].values()) or not all(tube['margin_checks'].values()):
        raise ValueError('require original strict cap tube')
    for key in ('alpha_max','duration','z_star'):
        if Q(tube[key])!=Q(scalar[key]):
            raise ValueError('scalar/cap '+key+' mismatch')
    al = Q(tube['alpha_max'])
    if al<=0:raise ValueError('require positive maximum coupling')
    with ctx.workprec(160):
        t = {k:_stored(v) for k,v in tube['constants'].items()}
        r, R = t['cap_inner_radius'],t['radius_upper']
        A, Amax = _l(1-R*R),_u(1-r*r)
        length = _u(R.atanh()-_ball(Q(tube['z_star'])))
        budget.check(2048)
        f = t['cap_field_upper']
        _, potential, _ = potential_majorants(f,_ball(al))
        remainder = potential_remainder_coefficient(f,_ball(al))
        values=dict(r_min=dyadic(r),r_max=_round_upper(R),
            A0_min=dyadic(A),A0_max=_round_upper(Amax),length_max=_round_upper(length),
            chi_max=Q(tube['box']['chi']),E_plus_max=Q(tube['box']['E_plus']),
            E_minus_max=Q(tube['box']['E_minus']),N_max=_round_upper(t['cap_mass_upper']),
            D_max=_round_upper(t['cap_lapse_upper']/_ball(al)),mu_max=_round_upper(t['cap_mu_upper']),
            potential_max=_round_upper(potential),potential_remainder=_round_upper(remainder),
            C_Z=Q(scalar['comparison']['wall_z']),C_chi=Q(scalar['coefficients']['chi']),
            C_E_plus=Q(scalar['coefficients']['E_plus']),C_E_minus=Q(scalar['coefficients']['E_minus']))
        result=radial_comparison(values,budget=budget)
        rounded={k:_round_upper(_ball(v)) for k,v in result.items()}
        coefficients={k:rounded[k] for k in ('N','D','mu_second_order','A_vacuum_subtracted_second_order')}
        budget.check(1024)
        return dict(alpha_max=al,duration=Q(tube['duration']),z_star=Q(tube['z_star']),
            inputs=values,comparison=rounded,coefficients=coefficients,
            maximum_amplitude_differences={k:al*coefficients[k] for k in ('N','D')},
            maximum_amplitude_second_order={
                'mass_excess':al*al*coefficients['N'],'lapse':al*al*coefficients['D'],
                'mu':al*al*coefficients['mu_second_order'],
                'A_vacuum_subtracted':al*al*coefficients['A_vacuum_subtracted_second_order']},
            precision_bits=160,coefficient_rounding_bits=128,
            scope='uniform normalized cap N,D differences at common eta on original onset; '
                  'second-order lapse and mass/A after subtracting vacuum at actual radius; '
                  'not the full radius/clock/total-mass Taylor remainder or a later arrival')
