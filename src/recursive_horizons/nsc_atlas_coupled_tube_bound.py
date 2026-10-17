"""One fixed vacuum-relative bootstrap image for the sourced continuum system.

All inputs are uniform certificates. The image is conditional on a smooth
solution in the declared box, with compatible quiet initial data. Strict
improvement is an a priori estimate, not a construction/existence theorem.
"""
from fractions import Fraction as Q

from flint import arb,ctx

from .nsc_atlas_background_cap import _budget
from .nsc_atlas_reference_gauge_bound import _ball,_record
from .nsc_atlas_reference_initial_bound import _from_record

DECLARED_BOX={
    'J_r':Q('1e-6'),'J_L':Q('1e-6'),'J_X':Q('1e-7'),'delta_T':Q('1e-6'),
    'E_plus':Q('.004'),'E_minus':Q('1e-8'),'chi':Q('1e-5'),
    'delta_r_U':Q('1e-5'),'delta_L_U':Q('.001'),'P':Q('1e-7')}
LONGITUDINAL_BOUND=Q(101,100)


def _u(value):
    if not value.is_finite():raise ValueError('nonfinite upper bound')
    return value.upper()


def _l(value):
    if not value.is_finite():raise ValueError('nonfinite lower bound')
    return value.lower()


def _positive(**values):
    for name,value in values.items():
        if not value>0:raise ValueError('positive coefficient not certified: '+name)


def potential_majorants(field,alpha):
    """Positive hyperbolic majorants including the regular alpha-zero limit."""
    ch=_u(((alpha/6).sqrt()*field).cosh())
    potential=_u(field**2*ch**2*(ch**2+1)/2)
    reaction=_u(field**3/3*(ch**4/6+ch**2*(ch**2+ch+1)/2))
    return ch,potential,reaction


def wall_difference_bounds(*,eps_r,eps_L,Verror,Lxerror,alpha,N,mu_error,
                           Rlo,Rhi,R0lo,R0hi,Amin,Hmin,Upmin,Upmax,
                           Vpmin,Vpmax,Wxmin,rx0max,Sw,root_change,
                           Sw_change,u_wall,nc,beta_bar_max):
    """Finite differences from the zero-field wall on the same exact gauge."""
    dUp=_u(Upmax*eps_L/2);dVp=_u(Vpmax*eps_L/2)
    dSR=_u(Sw*eps_r+R0hi*Sw_change)
    dRp=_u(2*Vpmax*Verror+2*rx0max*dVp)
    dH=_u(dRp+dSR)
    Aref=_u(1-R0lo**2)
    dA=_u(2*Rhi*eps_r+2*alpha*N/Rlo)
    dgap=_u(dA/Hmin+Aref*dH/Hmin**2)
    reference_numerator=_u(2*R0hi+Aref/Hmin)
    dpr=_u((dSR+dgap)/(2*Upmin)+reference_numerator*dUp/(2*Upmin**2))
    dbbar=_u((dRp+dgap)/2)
    bbar0=_u(1/(5*R0lo**2)-R0lo/4)
    # The coefficient 2 is required by the actual -2*beta_bar/R curvature term.
    dk=_u(3*root_change/4+alpha*u_wall*(nc+u_wall*beta_bar_max)/(2*Rlo)
          +2*dbbar/Rlo+2*bbar0*eps_r/Rlo**2+Sw_change/2)
    k0=_u(2/(5*R0lo**3)+_ball(Q(5,4)))
    dpL=_u((2*(dk/Vpmin+k0*dVp/Vpmin**2)+Lxerror)/Wxmin)
    dbminus=_u((dH+dgap)/2)
    bminus0=_u(1/(5*R0lo**2)+3*R0hi/4)
    # Both lapses have mu>=1. Split the reciprocal product before bounding it.
    dTx=_u(dbminus/(Amin*Vpmin)+bminus0*(dA/(Amin**2*Vpmin)
             +mu_error/(Amin*Vpmin)+dVp/(Amin*Vpmin**2)))
    return dict(r_U=dpr,L_U=dpL,T_x=dTx,k_plus=dk,A=dA,H=dH,A_over_H=dgap,
                beta_bar=dbbar,Udot=dUp,Vdot=dVp)


def coupled_vacuum_tube_bound(source,gauge,*,z_star,box,budget):
    budget=_budget(budget);budget.check(8192)
    if set(box)!=set(DECLARED_BOX):raise ValueError('require all ten declared norm bounds')
    exact={k:Q(v) for k,v in box.items()}
    if any(v<=0 for v in exact.values()):raise ValueError('box bounds must be positive')
    xlq,xrq=(Q(gauge['inputs'][k]) for k in ('left','right'))
    if Q(source['duration'])!=xrq-xlq:raise ValueError('source and gauge cells differ')
    with ctx.workprec(160):
        b={k:_ball(v) for k,v in exact.items()}
        s={k:_u(_ball(Q(v))) for k,v in source['forcing'].items()}
        al=_u(_ball(Q(source['alpha_max'])));dt=_u(_ball(xrq-xlq))
        xl,xr,a=_ball(xlq),_ball(xrq),_ball(Q(4,5))
        R0lo=_l(_ball(Q(gauge['inputs']['radius_lower'])))
        R0hi=_u(_ball(Q(gauge['inputs']['radius_upper'])))
        W=_u(_from_record(gauge['wall']['W']))
        Wx=_from_record(gauge['wall']['W_x']);Wxlo,Wxhi=_l(Wx),_u(Wx)
        er=_u(s['r']+b['J_r']);eL=_u(s['log_C']+b['J_L'])
        X=_u(s['X']+b['J_X'])
        rlo,rhi=_l(xl-er),_u(R0hi+er)
        Rlo,Rhi=_l(R0lo-er),rhi
        zs=_ball(Q(z_star));rstar=_l(zs.tanh())
        _positive(rlo=rlo,rstar=rstar,collar=Rlo-rstar,horizon=1-Rhi,
                  source_wall_gap=Rlo-xr,W=W,Wx=Wxlo)
        Cmax=_u(eL.exp())
        Cw0=_from_record(gauge['wall']['C_wall'])
        Cwlo=_l(Cw0.lower()*(-eL).exp());Cwhi=_u(Cw0.upper()*eL.exp())
        Upmin,Upmax=_l((Wxlo/Cwhi).sqrt()),_u((Wxhi/Cwlo).sqrt())
        Vpmin,Vpmax=_l(1/(Wxhi*Cwhi).sqrt()),_u(1/(Wxlo*Cwlo).sqrt())
        p0=_u((a-xl)/(4*xl));pabs=_u(p0+b['delta_r_U'])
        vmax=_ball(LONGITUDINAL_BOUND)
        k=_u((pabs*vmax+Cmax/4)/rlo**2)
        V=_u((s['r_times_r_x_minus_x']+W*Cmax*eL/4+vmax*er)/rlo)
        def kernel_difference(jet,span):
            return _u(vmax*jet/rlo**2+pabs*span*V/rlo**2+Cmax*span*eL/(4*rlo**2)
                      +2*k*span*er/rlo)
        KI,KS=kernel_difference(b['J_r'],W),kernel_difference(b['delta_r_U'],arb(1))
        Qmax=_u(s['Q']+k*W*X)
        fext=_u(X/rlo)
        FU=_u((b['P']+fext*pabs)/rlo);Fx=_u((Qmax+fext*vmax)/rlo)
        kinetic=_u(al*FU*Fx)
        Lxerror=_u(s['log_C_x']+2*KI+W*kinetic)
        # Radial cap constraints: N'+a_c N=b_c with a_c,b_c>=0 and zero core.
        A0min=_l(1-Rhi**2);A0max=_u(1-rstar**2);A0wallmax=_u(1-Rlo**2)
        fcap=_u(b['chi']/rstar)
        pcap=_u((b['E_plus']+b['E_minus'])/(2*A0min));hcap=_u(pcap+fcap)
        _,potential,reaction=potential_majorants(fcap,al)
        width=_u(Rhi-rstar);energy=_u(pcap**2+hcap**2)
        N=_u(width*(Rhi**2*potential/2+A0max*energy/4))
        delta=_u(al*width*energy/(2*rstar));mu=_u(delta.exp())
        Amin=_l(A0min-2*al*N/rstar);Amax=A0wallmax
        gamma=_u(al*mu*(2*N/rstar**2+Rhi*potential+4*N/A0min))
        reaction_source=_u(al*mu*(Rhi*reaction+(2*N/rstar**2+Rhi*potential)*fcap))
        # The wall uses its actual sourced radius, not a reference A0 maximum.
        rx0lo=_l(xl*(a-R0hi)/(R0hi*(a-xl)))
        rx0hi=_u(xr*(a-R0lo)/(R0lo*(a-xr)))
        rxlo,rxhi=_l(rx0lo-V),_u(rx0hi+V)
        fw=_u(X/Rlo);root=_u(((al/6).sqrt()*fw).cosh());Sw=_u(2*root**3)
        drt=_u(al*fw**2*root/12);dSw=_u(6*root**2*drt)
        Hlo=_l(2*Vpmin*rxlo+2*Rlo);Hhi=_u(2*Vpmax*rxhi+Sw*Rhi)
        _positive(Amin=Amin,wall_r_x=rxlo,Hmin=Hlo)
        gapmax=_u(Amax/Hlo);gapmin=_l(Amin/Hhi)
        beta_minus=_u((Hhi+gapmax)/2);beta_bar=_u((2*Vpmax*rxhi+gapmax)/2)
        Tx=_u(beta_minus/(Amin*Vpmin));beta=mu
        speed_min=_l((Hlo-gapmax)/2)
        inward_margin=_l(gapmin/(A0wallmax*Vpmax))
        outward_margin=_l(Hlo/(A0wallmax*Vpmax))
        _positive(expanding_wall=speed_min,inward_characteristic_gap=inward_margin,
                  outward_characteristic_gap=outward_margin)
        J=_u(Sw*Rhi*al*fw**3/18)
        frame_gap=_u(root*al*fw**3/18);u_wall=_u(fw*root)
        cx=_u(root*Qmax+frame_gap*rxhi)
        cm=_u(gapmax*(root*b['E_minus']/A0min+frame_gap))
        nc=_u(Vpmax*cx+cm/2)
        differences=wall_difference_bounds(eps_r=er,eps_L=eL,Verror=V,Lxerror=Lxerror,
            alpha=al,N=N,mu_error=_u(delta*mu),Rlo=Rlo,Rhi=Rhi,R0lo=R0lo,R0hi=R0hi,
            Amin=Amin,Hmin=Hlo,Upmin=Upmin,Upmax=Upmax,
            Vpmin=Vpmin,Vpmax=Vpmax,Wxmin=Wxlo,rx0max=rx0hi,
            Sw=Sw,root_change=drt,Sw_change=dSw,u_wall=u_wall,nc=nc,beta_bar_max=beta_bar)
        Ewall=_u(A0wallmax/Hlo*(2*Vpmax*Qmax+J))
        Pwall=_u((Amax*b['E_minus']/(Hlo*A0min)+J)/(2*Upmin))
        image=dict(
            J_r=_u(dt*(Wxhi*differences['r_U']+rhi*KI+k*W*er)),
            J_L=_u(dt*(Wxhi*differences['L_U']+2*KI+W*kinetic)),
            J_X=_u(dt*(Wxhi*Pwall+k*W*X)),
            delta_T=_u(dt*differences['T_x']),
            E_plus=_u(Ewall+dt*Tx*(gamma*b['E_plus']+A0max*reaction_source)),
            E_minus=_u(dt*Tx*(gamma*b['E_minus']+A0max*reaction_source)),
            chi=_u(X+dt*Tx*beta*(b['E_plus']+b['E_minus'])/2),
            delta_r_U=_u(differences['r_U']+dt*(rhi*KS+k*er)),
            delta_L_U=_u(differences['L_U']+dt*(2*KS+kinetic)),
            P=_u(Pwall+dt*k*X))
        margins=dict(longitudinal=_l(vmax-1-V),core_travel=_l(Rlo.atanh()-zs-dt*Tx*beta),
                     cap_A=Amin,wall_H=Hlo,inward_characteristic=inward_margin,
                     outward_characteristic=outward_margin,expanding_wall=speed_min)
        checks={k:bool(value<b[k]) for k,value in image.items()}
        guards={k:bool(value>0) for k,value in margins.items()}
        constants=dict(radius_floor=rlo,radius_upper=rhi,wall_radius_floor=Rlo,cap_inner_radius=rstar,
            W=W,W_x_upper=Wxhi,C_upper=Cmax,p_r_upper=pabs,K_upper=k,
            radius_value_error=er,lapse_value_error=eL,X_upper=X,Q_upper=Qmax,
            r_x_error=V,K_integrated_difference=KI,K_pointwise_difference=KS,L_x_error=Lxerror,
            cap_field_upper=fcap,cap_mass_upper=N,cap_lapse_upper=delta,cap_mu_upper=mu,
            cap_beta_upper=beta,cap_beta_z_upper=gamma,cap_reaction_upper=reaction_source,
            wall_clock_rate_upper=Tx,wall_incoming_upper=Ewall,wall_transverse_scalar_upper=Pwall)
        budget.check()
        return dict(precision_bits=160,box={k:str(v) for k,v in exact.items()},
                    alpha_max=source['alpha_max'],duration=source['duration'],z_star=str(Q(z_star)),
                    image={k:_record(v) for k,v in image.items()},
                    ratios={k:_record(v/b[k]) for k,v in image.items()},
                    margins={k:_record(v) for k,v in margins.items()},
                    constants={k:_record(v) for k,v in constants.items()},
                    wall_differences={k:_record(v) for k,v in differences.items()},
                    component_checks=checks,margin_checks=guards,strict=all(checks.values()) and all(guards.values()),
                    scope='conditional uniform bootstrap improvement; smooth solution construction and derivative continuation not supplied')
