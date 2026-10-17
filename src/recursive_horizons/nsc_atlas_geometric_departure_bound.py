"""Amplitude-proportional geometric bounds from the existing positive tube.

All coefficient and scalar bounds stay frozen at the certified tube. Only
three geometric errors participate in feedback. No scalar amplitude
derivative, new field solution, or quadratic remainder is computed here.
"""
from fractions import Fraction as Q

from flint import arb, ctx

from .nsc_atlas_coupled_tube_bound import wall_difference_bounds, LONGITUDINAL_BOUND
from .nsc_atlas_reference_gauge_bound import _ball
from .nsc_atlas_onset_geometry_bound import dyadic

FEEDBACK = ('J_r', 'J_L', 'delta_r_U')
OUTPUTS = (*FEEDBACK, 'delta_L_U', 'delta_T', 'radius', 'log_C', 'r_x', 'L_x')


def positive_resolvent(matrix, forcing, weights, *, budget):
    """Exact rational positive inverse, certified by a weighted row norm."""
    budget.check(2048)
    n = len(weights)
    if n != 3 or len(matrix) != n or len(forcing) != n or any(len(r)!=n for r in matrix):
        raise ValueError('require a three-component comparison')
    def exact(x):
        if isinstance(x,(bool,float)):
            raise ValueError('require explicit rational comparison inputs')
        return Q(x)
    M = [[exact(x) for x in r] for r in matrix]
    f, b = [exact(x) for x in forcing], [exact(x) for x in weights]
    if any(x<0 for r in M for x in r) or any(x<0 for x in f) or any(x<=0 for x in b):
        raise ValueError('require nonnegative matrix/forcing and positive weights')
    q = max(sum(M[i][j]*b[j] for j in range(n))/b[i] for i in range(n))
    if q >= 1:
        raise ValueError('weighted comparison norm is not below one')
    rows = [[Q(i==j)-M[i][j] for j in range(n)]+[f[i]] for i in range(n)]
    for j in range(n):
        pivot = rows[j][j]
        if pivot <= 0:
            raise ValueError('positive resolvent pivot failed')
        rows[j] = [x/pivot for x in rows[j]]
        for i in range(n):
            if i != j:
                scale = rows[i][j]
                rows[i] = [x-scale*y for x,y in zip(rows[i],rows[j])]
    result = [r[-1] for r in rows]
    if any(x<0 for x in result) or any(result[i]!=f[i]+sum(M[i][j]*result[j] for j in range(n)) for i in range(n)):
        raise ValueError('positive resolvent exact residual failed')
    budget.check()
    return dict(solution=result, weighted_norm=q)


def _stored(row, end='upper'):
    v=row[end]
    return _ball(Q(v['mantissa'])*Q(2)**v['exponent'])


def _u(v):
    if not v.is_finite():raise ValueError('nonfinite frozen coefficient')
    return v.upper()


def _l(v):
    if not v.is_finite():raise ValueError('nonfinite frozen coefficient')
    return v.lower()


def _round_upper(v):
    # Directed dyadic rounding limits rational record size; it is not a
    # floating estimate. All matrix entries and forcing coefficients are >=0.
    value=dyadic(_u(v)); scale=2**128
    if value<0:raise ValueError('negative affine coefficient')
    return Q(-((-value.numerator*scale)//value.denominator),scale)


def frozen_context(source, gauge, tube):
    """Reconstruct only fixed coefficients, retaining the receipt's floors."""
    c, margins = tube['constants'], tube['margins']
    if not tube['strict'] or not all(tube['component_checks'].values()) or not all(tube['margin_checks'].values()):
        raise ValueError('require the original strict tube')
    al=Q(source['alpha_max'])
    if not al>0 or al!=Q(tube['alpha_max']) or Q(source['duration'])!=Q(tube['duration']):
        raise ValueError('source/tube amplitude or duration mismatch')
    left,right=(Q(gauge['inputs'][k]) for k in ('left','right'))
    if right-left!=Q(tube['duration']):raise ValueError('gauge/tube duration mismatch')
    Amax=_ball(al)
    s={k:_u(_ball(Q(v))/Amax) for k,v in source['forcing'].items()
       if k in ('r','log_C','log_C_x','r_times_r_x_minus_x')}
    fixed={k:_stored(v) for k,v in c.items()}
    rlo,rhi,Rlo=(fixed[k] for k in ('radius_floor','radius_upper','wall_radius_floor'))
    R0lo,R0hi=(_ball(Q(gauge['inputs'][k])) for k in ('radius_lower','radius_upper'))
    xl,xr,a=_ball(left),_ball(right),_ball(Q(4,5))
    W,Wxhi=fixed['W'],fixed['W_x_upper']
    Wxlo=_stored(gauge['wall']['W_x'],'lower')
    eL=fixed['lapse_value_error']
    Cwlo=_l(_stored(gauge['wall']['C_wall'],'lower')*(-eL).exp())
    Cwhi=_u(_stored(gauge['wall']['C_wall'])*eL.exp())
    Upmin,Upmax=_l((Wxlo/Cwhi).sqrt()),_u((Wxhi/Cwlo).sqrt())
    Vpmin,Vpmax=_l(1/(Wxhi*Cwhi).sqrt()),_u(1/(Wxlo*Cwlo).sqrt())
    Amin,Hmin=_stored(margins['cap_A'],'lower'),_stored(margins['wall_H'],'lower')
    for name,val in dict(rlo=rlo,Rlo=Rlo,Upmin=Upmin,Vpmin=Vpmin,Wxlo=Wxlo,Amin=Amin,Hmin=Hmin).items():
        if not val>0:raise ValueError('positive frozen '+name+' required')
    rx0hi=_u(xr*(a-R0lo)/(R0lo*(a-xr)))
    rxhi=_u(rx0hi+fixed['r_x_error'])
    fw=_u(fixed['X_upper']/Rlo)
    root=_u(((Amax/6).sqrt()*fw).cosh())
    Sw=_u(2*root**3)
    drt=_u(fw*fw*root/12)
    dSw=_u(6*root*root*drt)
    gapmax=_u((1-Rlo*Rlo)/Hmin)
    bbar=_u((2*Vpmax*rxhi+gapmax)/2)
    framegap=_u(root*Amax*fw**3/18)
    cx=_u(root*fixed['Q_upper']+framegap*rxhi)
    eminus=_ball(Q(tube['box']['E_minus']))
    cm=_u(gapmax*(root*eminus/(1-rhi*rhi)+framegap))
    nc=_u(Vpmax*cx+cm/2)
    u_wall=_u(fw*root)
    mu_coefficient=_u(fixed['cap_lapse_upper']/Amax*fixed['cap_mu_upper'])
    pabs,k,Cmax=(fixed[key] for key in ('p_r_upper','K_upper','C_upper'))
    vmax=_ball(LONGITUDINAL_BOUND)
    fext=_u(fixed['X_upper']/rlo)
    FU=_u((_ball(Q(tube['box']['P']))+fext*pabs)/rlo)
    Fx=_u((fixed['Q_upper']+fext*vmax)/rlo)
    return dict(source=s,al=al,dt=_ball(Q(tube['duration'])),rlo=rlo,rhi=rhi,Rlo=Rlo,
        R0lo=R0lo,R0hi=R0hi,W=W,Wxhi=Wxhi,Wxlo=Wxlo,
        Upmin=Upmin,Upmax=Upmax,Vpmin=Vpmin,Vpmax=Vpmax,Amin=Amin,Hmin=Hmin,
        rx0hi=rx0hi,Sw=Sw,drt=drt,dSw=dSw,u_wall=u_wall,nc=nc,bbar=bbar,
        mu_coefficient=mu_coefficient,N=fixed['cap_mass_upper'],vmax=vmax,pabs=pabs,k=k,
        Cmax=Cmax,kinetic=_u(FU*Fx))


def affine_image(c, vector, scale):
    """Positive affine image with frozen coefficients; scale is forcing only.

    scale=0 removes explicit source and alpha terms, not the tube values.
    The unit-vector evaluations then enclose homogeneous coefficients.
    """
    jr,jl,pr=map(_ball,vector);t=_ball(scale)
    s=c['source'];rlo,rhi,W,k=c['rlo'],c['rhi'],c['W'],c['k']
    er=_u(t*s['r']+jr);eL=_u(t*s['log_C']+jl)
    V=_u((t*s['r_times_r_x_minus_x']+W*c['Cmax']*eL/4+c['vmax']*er)/rlo)
    def kernel(jet,span):
        return _u(c['vmax']*jet/rlo**2+c['pabs']*span*V/rlo**2
                  +c['Cmax']*span*eL/(4*rlo**2)+2*k*span*er/rlo)
    KI,KS=kernel(jr,W),kernel(pr,arb(1))
    Lx=_u(t*s['log_C_x']+2*KI+W*t*c['kinetic'])
    wall=wall_difference_bounds(eps_r=er,eps_L=eL,Verror=V,Lxerror=Lx,
        alpha=t,N=c['N'],mu_error=t*c['mu_coefficient'],Rlo=c['Rlo'],Rhi=rhi,
        R0lo=c['R0lo'],R0hi=c['R0hi'],Amin=c['Amin'],Hmin=c['Hmin'],
        Upmin=c['Upmin'],Upmax=c['Upmax'],Vpmin=c['Vpmin'],Vpmax=c['Vpmax'],
        Wxmin=c['Wxlo'],rx0max=c['rx0hi'],Sw=c['Sw'],root_change=t*c['drt'],
        Sw_change=t*c['dSw'],u_wall=c['u_wall'],nc=c['nc'],beta_bar_max=c['bbar'])
    dt=c['dt']
    return dict(J_r=_u(dt*(c['Wxhi']*wall['r_U']+rhi*KI+k*W*er)),
        J_L=_u(dt*(c['Wxhi']*wall['L_U']+2*KI+W*t*c['kinetic'])),
        delta_r_U=_u(wall['r_U']+dt*(rhi*KS+k*er)),
        delta_L_U=_u(wall['L_U']+dt*(2*KS+t*c['kinetic'])),
        delta_T=_u(dt*wall['T_x']),radius=er,log_C=eL,r_x=V,L_x=Lx)


def geometric_departure_bound(source,gauge,tube,*,budget):
    if not callable(getattr(budget,'check',None)):
        raise ValueError('require a shared cumulative budget')
    budget.check(8192)
    with ctx.workprec(160):
        c=frozen_context(source,gauge,tube)
        budget.check(4096)
        forcing={k:_round_upper(v) for k,v in affine_image(c,(0,0,0),1).items()}
        columns=[]
        for j in range(3):
            budget.check(4096)
            columns.append({k:_round_upper(v) for k,v in affine_image(c,tuple(int(i==j) for i in range(3)),0).items()})
        rows={k:[column[k] for column in columns] for k in OUTPUTS}
        weights=[Q(tube['box'][k]) for k in FEEDBACK]
        answer=positive_resolvent([rows[k] for k in FEEDBACK],
                                  [forcing[k] for k in FEEDBACK],weights,budget=budget)
        u=answer['solution']
        coefficients={k:forcing[k]+sum(rows[k][j]*u[j] for j in range(3)) for k in OUTPUTS}
        assert [coefficients[k] for k in FEEDBACK]==u
        budget.check(2048)
        return dict(feedback=list(FEEDBACK),alpha_max=Q(tube['alpha_max']),
            duration=Q(tube['duration']),matrix_rows=rows,forcing=forcing,weights=weights,
            weighted_norm=answer['weighted_norm'],coefficients=coefficients,
            endpoint_bounds={k:Q(tube['alpha_max'])*v for k,v in coefficients.items()},
            precision_bits=160,coefficient_rounding_bits=128,
            scope='uniform |geometry(alpha)-vacuum| <= alpha times the reported coefficient; '
                  'original onset cell and coupling interval; no scalar first difference, '
                  'second derivative, quadratic remainder or stored-trajectory error')
