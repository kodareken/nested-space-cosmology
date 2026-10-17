"""Uniform quadratic geometric remainder on the original onset cell.

Second jets differentiate algebraic maps over certified boxes. No second
scalar evolution is used. The existing positive geometric matrix closes
the comparison after source, scalar and normalized cap bounds are supplied.
"""
from dataclasses import dataclass, field
from fractions import Fraction as Q

from flint import arb, ctx

from .nsc_atlas_coupled_tube_bound import wall_difference_bounds
from .nsc_atlas_geometric_departure_bound import (
    FEEDBACK, OUTPUTS, frozen_context, positive_resolvent, _stored, _u, _l, _round_upper,
)
from .nsc_atlas_onset_geometry_bound import dyadic, _stored as interval
from .nsc_atlas_reference_gauge_bound import _ball


def ball(value):
    return value if isinstance(value, arb) else _ball(Q(value))


@dataclass(frozen=True)
class Jet:
    value: arb
    first: arb = field(default_factory=lambda: arb(0))
    second: arb = field(default_factory=lambda: arb(0))

    def __post_init__(self):
        for name in ('value','first','second'):
            value=ball(getattr(self,name))
            if not value.is_finite():raise ValueError('nonfinite second jet')
            object.__setattr__(self,name,value)

    @staticmethod
    def constant(value):
        return value if isinstance(value,Jet) else Jet(ball(value))

    def __add__(self,other):
        b=Jet.constant(other)
        return Jet(self.value+b.value,self.first+b.first,self.second+b.second)
    __radd__=__add__

    def __neg__(self):return Jet(-self.value,-self.first,-self.second)
    def __sub__(self,other):return self+-Jet.constant(other)
    def __rsub__(self,other):return Jet.constant(other)+-self

    def __mul__(self,other):
        b=Jet.constant(other);a=self
        return Jet(a.value*b.value,a.first*b.value+a.value*b.first,
                   a.second*b.value+2*a.first*b.first+a.value*b.second)
    __rmul__=__mul__

    def inverse(self):
        if not self.value>0:raise ValueError('positive reciprocal jet required')
        return Jet(1/self.value,-self.first/self.value**2,
                   2*(self.first*self.first)/self.value**3-self.second/self.value**2)

    def __truediv__(self,other):return self*Jet.constant(other).inverse()
    def __rtruediv__(self,other):return Jet.constant(other)*self.inverse()

    def __pow__(self,n):
        if type(n) is not int or n<0:raise ValueError('nonnegative integer jet power required')
        result=Jet(1)
        for _ in range(n):result=result*self
        return result

    def exp(self):
        e=self.value.exp()
        # Arb's generic power can be indeterminate on a zero-crossing
        # base. Interval multiplication encloses its integer square.
        return Jet(e,e*self.first,e*(self.second+self.first*self.first))


def jet_box(lo,hi,direction):
    lo,hi,direction=ball(lo),ball(hi),ball(direction)
    if not hi>=lo or not direction>=0:raise ValueError('ordered box and nonnegative direction required')
    return Jet(lo.union(hi),(-direction).union(direction),0)


def even_product(alpha_max,f):
    """Even hyperbolic jets for the declared physical product alpha*f^2/6.

    Actual alpha lies in [0,alpha_max]; an outward ball may extend slightly
    below zero. Nonnegativity follows from this constructor's exact domain,
    not from clamping an arbitrary interval argument to cosh(sqrt(q)).
    """
    if isinstance(alpha_max,(float,bool)):raise ValueError('require exact coupling endpoint')
    alpha_max=Q(alpha_max)
    if alpha_max<0:raise ValueError('nonnegative coupling endpoint required')
    alpha=Jet(arb(0).union(_ball(alpha_max)),1,0)
    q=alpha*f*f/6
    qmax=_u(_ball(alpha_max)*abs(f.value).upper()**2/6)
    chmax=_u(qmax.sqrt().cosh())
    def composition(den1,den2):
        value=arb(1).union(chmax)
        first=_ball(Q(1,den1)).union(_u(chmax/den1))
        second=_ball(Q(1,den2)).union(_u(chmax/den2))
        return Jet(value,first*q.first,second*(q.first*q.first)+first*q.second)
    return alpha,composition(2,12),composition(6,60),q


def wall_ports(*,R,L,qr,X,Qx,Em,N,D,Wx,alpha_max,qL=0,Wxx=0):
    """Signed canonical wall maps; affine qL/boost terms remain explicit."""
    Wx=ball(Wx)
    if not Wx>0:raise ValueError('positive fixed wall gauge required')
    j=Wx.sqrt();boost=(-L/2).exp()
    Up,Vp=boost*j,boost/j
    f=X/R
    alpha,ch,sh,_=even_product(alpha_max,f)
    Sw=2*ch**3;A0=1-R*R;A=A0-2*alpha*N/R;mu=(alpha*D).exp()
    H=2*Vp*qr+Sw*R;gap=A/H;beta=(2*Vp*qr+gap)/2
    K=-3*ch/4-alpha*f*sh/(2*R*ch)*(Vp*Qx-gap*Em/(2*A0)-f*beta)-2*beta/R-ch**3
    ports=dict(r_U=(Sw*R-gap)/(2*Up),L_U=(Jet.constant(qL)-2*K/Vp-Jet.constant(Wxx)/Wx)/Wx,
               T_x=(H+gap)/(2*A*mu*Vp))
    guards={k:dyadic(_l(v.value)) for k,v in dict(A0=A0,A=A,H=H,Udot=Up,Vdot=Vp).items()}
    if any(v<=0 for v in guards.values()):raise ValueError('wall Hessian box loses a positive floor')
    return ports,guards


def bulk_kernel(r,pr,qr,L):return (pr*qr+L.exp()/4)/(r*r)


def nonlinear_bounds(c,gauge,tube,geometry,scalar,cap,*,budget):
    budget.check(16384)
    t={k:_stored(v) for k,v in tube['constants'].items()}
    g={k:Q(v) for k,v in geometry['coefficients'].items()}
    f={k:Q(v) for k,v in scalar['coefficients'].items()}
    n={k:Q(v) for k,v in cap['coefficients'].items()}
    eL=t['lapse_value_error']
    Llo=_l(_stored(gauge['vacuum_wedge']['C'],'lower').log()-eL)
    Lhi=_u(_stored(gauge['vacuum_wedge']['C']).log()+eL)
    K=bulk_kernel(jet_box(c['rlo'],c['rhi'],g['radius']),
        jet_box(-c['pabs'],c['pabs'],g['delta_r_U']),
        jet_box(-c['vmax'],c['vmax'],g['r_x']),jet_box(Llo,Lhi,g['log_C']))
    cw=gauge['wall']['C_wall']
    lwlo=_l(_stored(cw,'lower').log()-eL);lwhi=_u(_stored(cw).log()+eL)
    xl=_ball(Q(gauge['inputs']['left']));a=_ball(Q(4,5))
    qr0lo=_l(xl*(a-c['R0hi'])/(c['R0hi']*(a-xl)))
    qrlo=_l(qr0lo-t['r_x_error']);qrhi=_u(c['rx0hi']+t['r_x_error'])
    if not qrlo>0:raise ValueError('positive wall longitudinal radius required')
    ports,guards=wall_ports(R=jet_box(c['Rlo'],c['rhi'],g['radius']),
        L=jet_box(lwlo,lwhi,g['log_C']),qr=jet_box(qrlo,qrhi,g['r_x']),
        X=jet_box(-t['X_upper'],t['X_upper'],f['X']),
        Qx=jet_box(-t['Q_upper'],t['Q_upper'],f['Q']),
        Em=jet_box(-_ball(Q(tube['box']['E_minus'])),_ball(Q(tube['box']['E_minus'])),f['E_minus']),
        N=jet_box(0,t['cap_mass_upper'],n['N']),
        D=jet_box(0,_u(t['cap_lapse_upper']/_ball(c['al'])),n['D']),
        Wx=interval(gauge['wall']['W_x']),alpha_max=c['al'])
    hessians={k:_round_upper(abs(v.second).upper()/2) for k,v in ports.items()}
    hessians['K']=_round_upper(abs(K.second).upper()/2)
    budget.check(4096)
    rlo,W=c['rlo'],c['W'];cr=_ball(g['radius'])
    fmax=_u(t['X_upper']/rlo)
    fu=_u((_ball(Q(tube['box']['P']))+fmax*c['pabs'])/rlo)
    fx=_u((t['Q_upper']+fmax*c['vmax'])/rlo)
    cf=_u(_ball(f['X'])/rlo+t['X_upper']*cr/rlo**2)
    cfu=_u((_ball(f['P'])+fmax*_ball(g['delta_r_U'])+c['pabs']*cf)/rlo+fu*cr/rlo)
    cfx=_u((_ball(f['Q'])+fmax*_ball(g['r_x'])+c['vmax']*cf)/rlo+fx*cr/rlo)
    kin=_u(fu*cfx+fx*cfu)
    fq=_u((W*c['Cmax']*_ball(g['log_C'])**2/8+cr*_ball(g['r_x']))/rlo)
    cross=_u(cr*_ball(Q(scalar['exterior_inputs']['kernel_pointwise'])))
    terms={k:_round_upper(v) for k,v in dict(radius_longitudinal=fq,kinetic=kin,radius_kernel_cross=cross).items()}
    details={k:_round_upper(v) for k,v in dict(f_upper=fmax,f_U_upper=fu,f_x_upper=fx,
                f_difference=cf,f_U_difference=cfu,f_x_difference=cfx).items()}
    budget.check()
    return dict(hessians=hessians,terms=terms,field_comparison=details,wall_floors=guards)


def remainder_image(c,source2,nonlinear,vector):
    """Positive image after dividing all remainder errors by alpha^2."""
    jr,jl,pr=map(_ball,vector)
    s={k:_ball(Q(v)) for k,v in source2.items()}
    h={k:_ball(Q(v)) for k,v in nonlinear['hessians'].items()}
    f={k:_ball(Q(v)) for k,v in nonlinear['terms'].items()}
    rlo,rhi,W,k=c['rlo'],c['rhi'],c['W'],c['k']
    er=_u(s['r']+jr);eL=_u(s['log_C']+jl)
    V=_u((s['r_times_r_x_minus_x']+W*c['Cmax']*eL/4+c['vmax']*er)/rlo+f['radius_longitudinal'])
    def kernel(jet,span):
        return _u(c['vmax']*jet/rlo**2+c['pabs']*span*V/rlo**2
            +c['Cmax']*span*eL/(4*rlo**2)+2*k*span*er/rlo+span*h['K'])
    KI,KS=kernel(jr,W),kernel(pr,arb(1))
    Lx=_u(s['log_C_x']+2*KI+W*f['kinetic'])
    wall=wall_difference_bounds(eps_r=er,eps_L=eL,Verror=V,Lxerror=Lx,
        alpha=arb(0),N=c['N'],mu_error=arb(0),Rlo=c['Rlo'],Rhi=rhi,
        R0lo=c['R0lo'],R0hi=c['R0hi'],Amin=c['Amin'],Hmin=c['Hmin'],
        Upmin=c['Upmin'],Upmax=c['Upmax'],Vpmin=c['Vpmin'],Vpmax=c['Vpmax'],
        Wxmin=c['Wxlo'],rx0max=c['rx0hi'],Sw=c['Sw'],root_change=arb(0),
        Sw_change=arb(0),u_wall=c['u_wall'],nc=c['nc'],beta_bar_max=c['bbar'])
    wr,wl,wt=(_u(wall[key]+h[key]) for key in ('r_U','L_U','T_x'))
    dt=c['dt']
    return dict(J_r=_u(dt*(c['Wxhi']*wr+rhi*KI+k*W*er+W*f['radius_kernel_cross'])),
        J_L=_u(dt*(c['Wxhi']*wl+2*KI+W*f['kinetic'])),
        delta_r_U=_u(wr+dt*(rhi*KS+k*er+f['radius_kernel_cross'])),
        delta_L_U=_u(wl+dt*(2*KS+f['kinetic'])),delta_T=_u(dt*wt),
        radius=er,log_C=eL,r_x=V,L_x=Lx)


def endpoint_comparison(first,coefficients,alpha_max,*,budget):
    budget.check(512)
    al=Q(alpha_max);result={}
    if al<=0:raise ValueError('positive endpoint coupling required')
    for key,name in (('radius','radius_fixed_source'),('delta_T','clock_fixed_source')):
        lo,hi=(Q(first[name][end]) for end in ('lower','upper'))
        c2=Q(coefficients[key])
        if not lo<=hi<0 or c2<0:raise ValueError('require certified negative first interval and positive remainder')
        ratio=al*c2/(-hi)
        # Equality only proves nonpositivity. A reduced range uses half
        # the threshold so the sign remains strict at its closed endpoint.
        sign_alpha=al if ratio<1 else min(al,(-hi)/(2*c2))
        result[key]=dict(first_lower=lo,first_upper=hi,remainder=c2,
            relative_error_at_maximum=ratio,strict_on_full_interval=ratio<1,
            strict_sign_alpha_max=sign_alpha,
            normalized_lower_at_maximum=lo-al*c2,normalized_upper_at_maximum=hi+al*c2,
            shift_lower_at_maximum=al*lo-al*al*c2,shift_upper_at_maximum=al*hi+al*al*c2)
    budget.check()
    return result


def geometric_remainder_bound(source,gauge,tube,geometry,scalar,cap,source2,first,*,budget):
    if not callable(getattr(budget,'check',None)):raise ValueError('require a shared cumulative budget')
    budget.check(8192)
    for record in (geometry,scalar,cap,source2):
        if Q(record['alpha_max'])!=Q(tube['alpha_max']) or Q(record['duration'])!=Q(tube['duration']):
            raise ValueError('original amplitude/domain mismatch')
    for record in (scalar,cap):
        if Q(record['z_star'])!=Q(tube['z_star']):raise ValueError('common collar mismatch')
    for key in ('left','right'):
        if Q(source2[key])!=Q(gauge['inputs'][key]):raise ValueError('source remainder cell mismatch')
    if geometry['feedback']!=list(FEEDBACK):raise ValueError('original feedback order changed')
    with ctx.workprec(160):
        c=frozen_context(source,gauge,tube)
        nonlinear=nonlinear_bounds(c,gauge,tube,geometry,scalar,cap,budget=budget)
        budget.check(4096)
        forcing={k:_round_upper(v) for k,v in remainder_image(c,source2['coefficients'],nonlinear,(0,0,0)).items()}
        rows={k:[Q(x) for x in geometry['matrix_rows'][k]] for k in OUTPUTS}
        weights=[Q(x) for x in geometry['weights']]
        answer=positive_resolvent([rows[k] for k in FEEDBACK],[forcing[k] for k in FEEDBACK],weights,budget=budget)
        if answer['weighted_norm']!=Q(geometry['weighted_norm']):raise ValueError('frozen matrix norm changed')
        u=answer['solution']
        coefficients={k:forcing[k]+sum(rows[k][j]*u[j] for j in range(3)) for k in OUTPUTS}
        assert [coefficients[k] for k in FEEDBACK]==u
        endpoint=endpoint_comparison(first['outputs'],coefficients,c['al'],budget=budget)
        budget.check(2048)
        return dict(alpha_max=c['al'],duration=Q(tube['duration']),feedback=list(FEEDBACK),
            nonlinear=nonlinear,forcing=forcing,matrix_rows=rows,weights=weights,
            weighted_norm=answer['weighted_norm'],coefficients=coefficients,
            maximum_amplitude_remainders={k:c['al']**2*v for k,v in coefficients.items()},
            endpoint=endpoint,precision_bits=160,coefficient_rounding_bits=128,
            common_strict_sign_alpha_max=min(v['strict_sign_alpha_max'] for v in endpoint.values()),
            strict_endpoint_sign_on_full_interval=all(v['strict_on_full_interval'] for v in endpoint.values()),
            scope='uniform |G(alpha)-G0-alpha*G1|<=alpha^2*C on original onset; '
                  'endpoint Einstein radius and centre-anchored clock at fixed source label; '
                  'strict signs concern positive alpha only; '
                  'not absolute contraction, evaporation, later arrival or a stored-trajectory error')
