"""Regular two-radius Euclidean saddle and its open Lorentzian continuation.

This is the fixed-q, source-free subsector u=(a,b,b,b) of the combined
parent action. The massive scalar is exactly at phi_star(q). No membrane,
charge jump, incoming pulse or independently chosen exit state is inserted.
Existence of a profile does not prove its full physical fluctuation index.
"""
from dataclasses import dataclass
from fractions import Fraction
import math
import time

import numpy as np
from scipy.integrate import solve_bvp,solve_ivp,simpson,quad
from scipy.linalg import eigh

from .nsc_atlas_sphere_flux_completion import (
    log_radius_potential,log_radius_gradient,log_radius_hessian,
    log_radius_kinetic_metric,symmetric_radius_branches,action_reduction_record,
    rolling_parent_state,
)

EMBED=np.array([[1.,0.],[0.,1.],[0.,1.],[0.,1.]])
METRIC=EMBED.T@log_radius_kinetic_metric()@EMBED
INVERSE=np.linalg.inv(METRIC)


def expand(fields):return np.asarray(fields,dtype=float)@EMBED.T
def potential(fields):return log_radius_potential(expand(fields))
def gradient(fields):return log_radius_gradient(expand(fields))@EMBED
def hessian(fields):return np.einsum('ai,...ab,bj->...ij',EMBED,log_radius_hessian(expand(fields)),EMBED)


def reference_points():
    rows=symmetric_radius_branches()['branches']
    saddle=next(row for row in rows if row['branch']=='1+3' and row['homogeneous_negative_count']==1)
    false=next(row for row in rows if row['branch']=='isotropic' and row['V']>0 and row['homogeneous_negative_count']==0)
    ads=next(row for row in rows if row['branch']=='1+3' and row['V']<0 and row['homogeneous_negative_count']==0)
    return false,saddle,ads


def euclidean_rhs(x,y,length):
    fields=y[:2].T;vel=y[2:4].T;R=y[4];P=y[5]
    value=potential(fields);force=gradient(fields)@INVERSE
    kin=np.einsum('ni,ij,nj->n',vel,METRIC,vel)
    return length*np.vstack([vel.T,(force-3*P[:,None]/R[:,None]*vel).T,P,-R*(kin+value)/3])


def pole_boundaries(left,right,length,cutoff):
    d=length*cutoff;vl=potential(left[:2]);vr=potential(right[:2])
    fl=INVERSE@gradient(left[:2]);fr=INVERSE@gradient(right[:2])
    return np.r_[left[4]-d+vl*d**3/18,left[5]-1+vl*d*d/6,
        left[2:4]-d*fl/4,right[4]-d+vr*d**3/18,right[2:4]+d*fr/4]


@dataclass
class BounceSolution:
    solution: object
    cutoff: float
    attempts: list

    @property
    def length(self):return float(self.solution.p[0])

    def pole_fields(self):
        delta=self.length*self.cutoff
        return tuple(point-INVERSE@gradient(point)*delta**2/8
            for point in (self.solution.y[:2,0],self.solution.y[:2,-1]))


def solve_bounce(*,tolerance=2e-6,cutoff=1e-4,mesh=161,previous=None,budget_seconds=60.):
    if not (0<tolerance<.01 and 0<cutoff<.001 and mesh>=65 and budget_seconds>0):
        raise ValueError('positive bounded tolerance/cutoff/budget and mesh>=65 required')
    started=time.perf_counter();false,saddle,_=reference_points()
    center=np.asarray(saddle['log_radii'])[[0,1]]
    masses,vectors=eigh(hessian(center),METRIC);direction=vectors[:,0]
    if direction[0]>0:direction=-direction
    length=math.pi*math.sqrt(3/saddle['V'])
    x=np.linspace(cutoff,1-cutoff,mesh);attempts=[]
    for seed in ((None,) if previous is not None else (.4,.8)):
        if previous is not None:
            length=previous.length;guess=previous.solution.sol(x)
        else:
            fields=center[:,None]+seed*direction[:,None]*np.cos(math.pi*x)
            vel=-seed*direction[:,None]*math.pi*np.sin(math.pi*x)/length
            guess=np.vstack([fields,vel,np.sin(math.pi*x)*length/math.pi,np.cos(math.pi*x)])
        def rhs(at,values,p):
            if time.perf_counter()-started>budget_seconds:raise TimeoutError('bounded coupled bounce allowance exhausted')
            if np.max(abs(values[:2]))>5 or p[0]<=0:raise FloatingPointError('uncontrolled Newton trial; no profile admitted')
            return euclidean_rhs(at,values,p[0])
        try:
            with np.errstate(over='raise',invalid='raise',divide='raise'):
                sol=solve_bvp(rhs,lambda a,b,p:pole_boundaries(a,b,p[0],cutoff),
                    x,guess,p=[length],tol=tolerance,max_nodes=5000)
            field_span=np.ptp(sol.y[:2],axis=1)
            admitted=bool(sol.success and max(field_span)>.03 and min(sol.y[4])>0 and sol.p[0]>0)
            attempts.append(dict(seed_amplitude=seed,status=int(sol.status),message=sol.message,
                maximum_field_span=float(max(field_span)),nonconstant_profile=admitted,nodes=len(sol.x)))
            if admitted:return BounceSolution(sol,cutoff,attempts)
        except (ValueError,FloatingPointError) as exc:attempts.append(dict(seed_amplitude=seed,error=str(exc)))
    raise RuntimeError('No nonconstant profile admitted in bounded seeds: '+str(attempts))


def bounce_diagnostics(bounce):
    sol=bounce.solution;L=bounce.length;cut=bounce.cutoff
    x=np.linspace(cut,1-cut,4001);y=sol.sol(x);dy=sol.sol(x,1)
    fields=y[:2].T;vel=y[2:4].T;R,P=y[4:];V=potential(fields)
    kin=np.einsum('ni,ij,nj->n',vel,METRIC,vel)
    C=P*P-1-R*R*(kin/2-V)/3
    residual=dy-euclidean_rhs(x,y,L)
    actual_Rpp=dy[5]/L
    R3=-4*actual_Rpp/R+2*(1-P*P)/R**2
    ham=R3-kin-2*V
    pole_n,pole_s=bounce.pole_fields();delta=L*cut;vn=potential(pole_n);vs=potential(pole_s)
    trace=-2*math.pi**2*(L*simpson(R**3*V,x=x)+(vn+vs)*delta**4/4)
    full=2*math.pi**2*(L*simpson(R**3*(kin/2+V)-3*R*(P*P+1),x=x)
        -6*delta**2+7*(vn+vs)*delta**4/12)
    false,saddle,_=reference_points();false_action=-24*math.pi**2/false['V']
    common_scale=false['V']/3+abs(R3)+abs(kin+2*V)
    Q=1-R*R*kin/6
    transverse=np.array([0.,1.,-1.,0.])/math.sqrt(2)
    H=log_radius_hessian(expand(fields));m_trans=np.einsum('i,nij,j->n',transverse,H,transverse)/3
    direction=np.asarray(saddle['negative_mode'])[[0,1]]
    projection=(fields-np.asarray(saddle['log_radii'])[[0,1]])@METRIC@direction
    crossings=int(np.sum(projection[:-1]*projection[1:]<0))
    equator=int(np.argmin(abs(P)))
    return dict(length=L,cutoff=cut,mesh_nodes=len(sol.x),attempts=bounce.attempts,
        pole_fields=[pole_n.tolist(),pole_s.tolist()],pole_potentials=[float(vn),float(vs)],
        max_Euclidean_constraint=float(max(abs(C))),max_common_Cauchy_Hamiltonian=float(max(abs(ham))),
        max_common_Cauchy_relative_Hamiltonian=float(max(abs(ham)/common_scale)),
        common_Cauchy_Hamiltonian_RMS=float(np.sqrt(simpson(R*R*ham*ham,x=x)/simpson(R*R,x=x))),
        max_scaled_ODE_residual=float(np.max(abs(residual)/(1+abs(euclidean_rhs(x,y,L))))),
        right_pole_slope_defect=float(y[5,-1]+1-vs*delta**2/6),
        action_trace=trace,action_full=full,action_identity_difference=full-trace,
        action_difference_times_g=trace-false_action,minimum_Q=float(min(Q)),
        sampled_transverse_shape_mass_squared_minimum=float(min(m_trans)),
        projected_barrier_crossings=crossings,maximum_radius=float(max(R)),
        equator_x=float(x[equator]),
        full_physical_negative_mode_count_proved=False,
        scope='Regular nonconstant profile in the1+3 invariant scalar sector. Scalar-metric fluctuation index, internal KK modes and decay-rate interpretation are not certified.')


def lorentzian_rhs(t,y):
    fields=y[:2];vel=y[2:4];a,P=y[4:]
    V=potential(fields);kin=vel@METRIC@vel
    return np.r_[vel,-3*P/a*vel-INVERSE@gradient(fields),P,-a*(kin-V)/3]


def continue_pole(bounce,*,pole_index=0,maximum_time=30.,tolerance=2e-11,budget_seconds=30.):
    if pole_index not in (0,1) or maximum_time<=0:raise ValueError('one regular pole and positive time allowance required')
    started=time.perf_counter();pole=bounce.pole_fields()[pole_index];force=INVERSE@gradient(pole)
    V0=potential(pole);t0=1e-5
    initial=np.r_[pole-force*t0*t0/8,-force*t0/4,t0+V0*t0**3/18,1+V0*t0*t0/6]
    def rhs(t,y):
        if time.perf_counter()-started>budget_seconds:raise TimeoutError('bounded Lorentzian allowance exhausted')
        return lorentzian_rhs(t,y)
    def turning(t,y):return y[5]
    turning.terminal=True;turning.direction=-1
    def acceleration_zero(t,y):return potential(y[:2])-y[2:4]@METRIC@y[2:4]
    acceleration_zero.direction=0
    sol=solve_ivp(rhs,(t0,maximum_time),initial,method='DOP853',rtol=tolerance,atol=tolerance/100,
        max_step=.01,events=(turning,acceleration_zero),dense_output=True)
    if not sol.success:raise RuntimeError(sol.message)
    times=np.linspace(t0,sol.t[-1],501);y=sol.sol(times);fields=y[:2].T;vel=y[2:4].T
    kin=np.einsum('ni,ij,nj->n',vel,METRIC,vel);V=potential(fields);rho=kin/2+V;p=kin/2-V;a,P=y[4:]
    constraint=P*P-1-a*a*rho/3
    loss=quad(lambda t:3*sol.sol(t)[5]/sol.sol(t)[4]*(sol.sol(t)[2:4]@METRIC@sol.sol(t)[2:4]),
        t0,sol.t[-1],epsabs=1e-10)[0]
    candidates=np.r_[t0,sol.t[-1],sol.t_events[1]]
    Pmax=float(max(sol.sol(candidates)[5]));mean_log=(fields[:,0]+3*fields[:,1])/4
    bulk_scale=a*np.exp(-6*mean_log)
    false,_,ads=reference_points()
    false_delta=fields[-1]-np.asarray(false['log_radii'])[[0,1]]
    ads_delta=fields[-1]-np.asarray(ads['log_radii'])[[0,1]]
    return sol,dict(pole_index=pole_index,exit_fields=pole.tolist(),exit_density=float(V0),
        exit_excitation_above_known_AdS_minimum=float(V0-ads['V']),
        initial_Ricci_scalar_limit=float(4*V0),initial_Kretschmann_limit=float(8*V0*V0/3),
        proper_endpoint=float(sol.t[-1]),endpoint_reason='first_scale_factor_turn' if len(sol.t_events[0]) else 'declared_time_allowance',
        scale_at_endpoint=float(a[-1]),density_at_endpoint=float(rho[-1]),pressure_at_endpoint=float(p[-1]),
        maximum_Friedmann_constraint=float(max(abs(constraint))),density_work_balance_residual=float(rho[-1]-rho[0]+loss),
        minimum_curvature_fraction=1/Pmax**2,maximum_bulk_ruler_scale_sampled=float(max(bulk_scale)),
        log_radii_at_endpoint=fields[-1].tolist(),incoming_material_flux=0.,
        endpoint_field_distance_to_false_minimum=float(np.sqrt(false_delta@METRIC@false_delta)),
        endpoint_field_distance_to_known_AdS_minimum=float(np.sqrt(ads_delta@METRIC@ads_delta)),
        scope='Open-FRW region from the regular pole, no independent exit/matter state. Finite curvature at a=0 is checked by the regular series. Initial geometric excitation/work is not measured incoming parent material energy.')


def connecting_geometry(bounce):
    """The same profile in the connecting chart and its common Cauchy cut.

    ds4²=dξ²+R_E(ξ)²[-dT²+cosh²(T)dOmega2²]. At T=0, K=Pi=0.
    Its contracting/expanding cuts are not a separately prepared false-vacuum
    past or a proof of black-hole nucleation. Material normal flux is zero.
    """
    from scipy.optimize import brentq
    sol=bounce.solution;cut=bounce.cutoff;L=bounce.length
    equator=brentq(lambda x:float(sol.sol(x)[5]),cut,1-cut,xtol=2e-13)
    y=sol.sol(equator);R=float(y[4]);u=expand(y[:2]);s=float(np.mean(u))
    Rpp=float(sol.sol(equator,1)[5]/L)
    V0=action_reduction_record().reference_volume
    area4=4*math.pi*R**2
    area16=(4*math.pi*R**2*math.exp(-12*s))*(V0*math.exp(3*sum(u)))
    x=np.linspace(cut,1-cut,4001);values=sol.sol(x);vel=values[2:4].T
    kin=np.einsum('ni,ij,nj->n',vel,METRIC,vel);rho=kin/2+potential(values[:2].T)
    energy=4*math.pi*L*simpson(values[4]**2*rho,x=x)
    volume=4*math.pi*L*simpson(values[4]**2,x=x)
    return dict(metric='ds4²=dξ²+R_E(ξ)²[-dT²+cosh²T dOmega2²]',
        equator_x=equator,equator_radius=R,equator_spatial_second_derivative=Rpp,
        T0_common_Cauchy_volume=volume,T0_common_Cauchy_normal_energy=energy,
        T0_Kij=0.,T0_normal_field_momenta=0.,normal_material_flux=0.,
        equator_null_expansions={str(T):[2*math.tanh(T)/R]*2 for T in (-1.,0.,1.)},
        equator_geodesic_radius_second_proper_derivative=1/R,
        dimensional_screen_area_ratio=area16/(V0*area4),
        area_identity='A14=4piR_E² exp(-12s)*V0 exp(12s)=V0 A2; kappa16=V0 kappa4.',
        null_expansion_identity='theta16=exp(6s) theta4 for the screen S2 x M12; signs agree.',
        scope='One inherited source profile and regular overlapping charts. The turning equator is a spatial maximum, not a persisting neck minimum. Area/entropy covariance is not a measured incoming energy supply.')


def evolve_rolling_parent(*,initial_s,initial_v,efolds=40.,q=2.,tolerance=2e-10,budget_seconds=20.):
    """A declared member of the parent family; no finite state is selected.

    y=(s,ds/dN,bulk proper duration). Einstein a_E/a_E0=exp(N).
    Integrates exact full-potential equations, not the asymptotic power law.
    """
    if not (0<efolds<=60 and 0<tolerance<.001 and budget_seconds>0):
        raise ValueError('bounded e-fold/tolerance/compute domain required')
    rolling_parent_state(initial_s,initial_v,q);started=time.perf_counter()
    def rhs(N,y):
        if time.perf_counter()-started>budget_seconds:raise TimeoutError('rolling-parent allowance exhausted')
        row=rolling_parent_state(y[0],y[1],q)
        return np.array([y[1],row['v_N'],1/row['H_scaled']])
    sol=solve_ivp(rhs,(0,efolds),[initial_s,initial_v,0.],method='DOP853',
        rtol=tolerance,atol=tolerance/100,max_step=.1,dense_output=True)
    if not sol.success:raise RuntimeError(sol.message)
    return sol


def flat_horizon_pp_record():
    """Test the flat Nariai idea in the ACTUAL source, in an affine frame.

    For ds16=-dT²+A(T)²dx³+b(T)²gamma12, an external affine null ray
    has dT/dlambda=P/A. The internal unit screen is parallel transported,
    and R(k,e_a,k,e_a)=-b''(lambda)/b. Finite scalar curvature does not
    establish a C² extension. No numerical horizon fitting is performed.
    """
    from .nsc_atlas_sphere_flux_completion import all_volume_extrema,volume_thresholds
    top=next(row for row in all_volume_extrema(1.5) if row['kind']=='maximum')
    alpha=-3*top['canonical_mass_squared']/top['vacuum']
    p=2*alpha/(math.sqrt(9+4*alpha)+3)
    p1_threshold=Fraction(175,22)**5*Fraction(28,55)
    reference_A5=Fraction(7)**5;merger_A5=volume_thresholds()['critical_A_fifth']
    assert p1_threshold<reference_A5<merger_A5
    product=p*(p-1);external_H2=math.exp(12*top['s'])*top['vacuum']/3
    internal_section=(5/8)*math.exp(-2*top['s'])
    return dict(q=1.5,barrier_s=top['s'],barrier_V=top['vacuum'],
        minus_mass_squared_over_H_squared=alpha,departure_exponent=p,
        exact_reference_A_fifth=reference_A5,exact_p1_threshold_A_fifth=p1_threshold,
        exact_merger_A_fifth=merger_A5,certified_range='0<p<1',
        affine_relation='d(lambda)=exp(-12s)*a_E*dt_E/P; s-stop=C_lambda*lambda^p+... ',
        internal_pp_coefficient_per_C=-product,external_pp_coefficient_per_C=6*product,
        pp_power=p-2,leading_null_Ricci_trace=12*(-product)+2*(6*product),
        full_Ricci_limit=16+(25/56)*math.exp(-12*top['s']),
        full_Kretschmann_limit=24*external_H2**2+48*internal_section**2,
        departing_horizon_is_C2=False,finite_double_tidal_integral=True,
        divergent_single_tidal_integral=True,
        scope='Nonconstant, exactly homogeneous common-radius flat-patch departure only. '
            'A finite-distortion pp null singularity is not a smooth classical horizon. '
            'The exact static top is regular but has no transition; closed Cauchy '
            'preparations are different solutions, not coordinate repairs of this one.')


def critical_parent_response(carrier_scale=.25,budget_seconds=30.):
    """Recompute the conditioned linear critical history and its weak CK root.

    This uses the fixed global-de Sitter barrier background, not a nonlinear
    critical solution or a natural/quantum state-selection rule. It inserts
    no balanced perturbation, solves no Hamiltonian constraint and evolves
    no GR Cauchy data. All returned coefficients multiply unnormalized U_l.
    """
    from scipy.special import eval_gegenbauer,roots_jacobi
    from .nsc_atlas_parent_circulation import Parent
    from .nsc_atlas_sphere_flux_completion import all_volume_extrema,combined_record

    if (isinstance(carrier_scale,(bool,np.bool_)) or not math.isfinite(float(carrier_scale)) or
        carrier_scale<0 or not math.isfinite(float(budget_seconds)) or budget_seconds<=0):
        raise ValueError('finite nonnegative carrier scale and positive allowance required')
    started=time.perf_counter();scale=float(carrier_scale)
    model=combined_record();K=float(model.K);extrema=all_volume_extrema(1.5)
    top=next(r for r in extrema if r['kind']=='maximum');low=next(r for r in extrema if r['kind']=='minimum')
    H2=top['vacuum']/3;H=math.sqrt(H2);mu2=400*math.exp(-12*top['s'])
    alpha=-top['canonical_mass_squared']/H2;p=(-3+math.sqrt(9+4*alpha))/2
    retained=Parent.retained()
    A=scale*math.sqrt(24/(retained.mass_squared/retained.H**2+3))
    B=scale*retained.H**1.5*math.sqrt(12)/math.pi
    initial=np.array([0.,A/2,0.,-B*math.sqrt(.6)/(12*H),0.,-B*math.sqrt(.6)/(4*H)])
    def guard():
        if time.perf_counter()-started>budget_seconds:raise TimeoutError('critical-parent linear-response allowance exhausted')
    def carrier_rhs(T,y):
        guard();return np.r_[y[3:],-3*np.tanh(T)*y[3:]-(mu2/H2+np.array([0,3,8])/np.cosh(T)**2)*y[:3]]
    carrier=solve_ivp(carrier_rhs,(0.,12.),initial,method='DOP853',rtol=2e-11,atol=2e-13,max_step=.02,dense_output=True)
    def unit_rhs(T,y):
        guard();return [y[1],-3*math.tanh(T)*y[1]-mu2/H2*y[0]]
    unit=solve_ivp(unit_rhs,(0.,12.),[0.,1/H],method='DOP853',rtol=2e-11,atol=2e-13,max_step=.02,dense_output=True)
    if not carrier.success or not unit.success:raise RuntimeError('massive linear source history failed')
    prefactor=6*mu2/(math.sqrt(K)*H2)
    def source(T,index):
        if index==2:return prefactor*unit.sol(T)[0:1]**2
        q=carrier.sol(T)[:3]
        if index==1:return 2*prefactor*unit.sol(T)[0]*q
        a,b,c=q
        return prefactor*np.array([a*a+b*b+c*c,2*b*(a+c),b*b+2*a*c+c*c,2*b*c,c*c])
    audit_T=np.array([0.,.031,.113,.37,1.,2.,5.,10.,12.])
    def response(channels,source_index,cutoff=12.,audit=False):
        labels=np.array(channels);n=len(labels);eigen=labels*(labels+2)
        def rhs(T,y):
            guard();return np.vstack([y[n:],-3*np.tanh(T)*y[n:]-(eigen[:,None]/np.cosh(T)**2-alpha)*y[:n]+source(T,source_index)])
        def jac(T,y):
            out=np.zeros((2*n,2*n,len(T)))
            for i in range(n):out[i,n+i]=1.;out[n+i,i]=alpha-eigen[i]/np.cosh(T)**2;out[n+i,n+i]=-3*np.tanh(T)
            return out
        def bc(a,b):return np.r_[[a[n+i] if l%2==0 else a[i] for i,l in enumerate(labels)],b[n:]+(3+p)*b[:n]]
        def bcjac(a,b):
            da=np.zeros((2*n,2*n));db=da.copy()
            for i,l in enumerate(labels):da[i,n+i if l%2==0 else i]=1.;db[n+i,i]=3+p;db[n+i,n+i]=1.
            return da,db
        mesh=np.linspace(0.,cutoff,int(cutoff*60)+1)
        sol=solve_bvp(rhs,bc,mesh,np.zeros((2*n,len(mesh))),fun_jac=jac,bc_jac=bcjac,tol=2e-9,max_nodes=30000)
        if not sol.success:raise RuntimeError('critical response: '+sol.message)
        check=np.linspace(0.,cutoff,4001);value=sol.sol(check);residual=sol.sol(check,1)-rhs(check,value)
        info=dict(channels=labels.tolist(),cutoff=cutoff,nodes=len(sol.x),
            max_scaled_ODE_residual=float(np.max(abs(residual)/(1+abs(rhs(check,value))))),
            boundary_max=float(np.max(abs(bc(sol.sol(0.),sol.sol(cutoff))))),
            maximum_tail_source=float(np.max(abs(source(np.array([cutoff]),source_index)))))
        if audit:
            value=sol.sol(audit_T);derivative=sol.sol(audit_T,1)
            info['audit']=dict(z=value[:n].T.tolist(),z_T=value[n:].T.tolist(),
                z_T_from_polynomial=derivative[:n].T.tolist(),z_TT=derivative[n:].T.tolist())
        return sol.sol(0.),info
    base,base_info=response(range(5),0,audit=True)
    shorter,short_info=response(range(5),0,cutoff=10.)
    first,first_info=response((0,1,2),1,audit=True)
    second,second_info=response((0,),2,audit=True)
    Sbar=math.sqrt(K)*(top['s']-low['s'])
    S0=base[:5].copy();S0[0]+=Sbar;Pi0=H*base[5:]
    S1=np.zeros(5);Pi1=np.zeros(5);S2=np.zeros(5)
    S1[:3]=first[:3];Pi1[:3]=H*first[3:];S2[0]=second[0]
    S0[[1,3]]=0.;S1[[1,3,4]]=0.;Pi0[[0,2,4]]=0.;Pi1[[0,2,3,4]]=0.
    u,w=roots_jacobi(24,1.5,1.5);w/=sum(w)
    Y=np.array([eval_gegenbauer(l,1,u) for l in range(5)])
    DY=np.array([np.zeros_like(u) if l==0 else 2*eval_gegenbauer(l-1,2,u) for l in range(5)])
    mean=lambda x:float(w@x)
    C0=mean((Pi0@Y)*(S0@DY))
    C1=mean((Pi1@Y)*(S0@DY)+(Pi0@Y)*(S1@DY));C2=mean((Pi1@Y)*(S1@DY))
    linear=A+C1;discriminant=linear**2-4*C2*C0;d2=C1*C1-4*C2*C0;d1=2*A*C1
    candidates=[0.,1.]
    if d2>0 and 0<-d1/(2*d2)<1:candidates.append(-d1/(2*d2))
    minimum_discriminant=min(A*A+d1*x+d2*x*x for x in candidates)
    result=dict(model=dict(q=1.5,carrier_scale=scale,K=K,s_top=top['s'],sstar=low['s'],
            V_top=top['vacuum'],mu2=mu2,alpha=alpha,p=p,canonical_S_bar=Sbar,source_prefactor=prefactor),
        H=H,H2=H2,carrier_A=A,carrier_B=B,
        normalization='Unnormalized U_l=C_l^1; canonical S=sqrt84(s-sstar), properPiS=H*z_T. Delta is constant proper/conformal initialPiChi; h0_T(0)=1/H.',
        response_basis=dict(S0_U=S0.tolist(),Pi0_U=Pi0.tolist(),S1_U=S1.tolist(),Pi1_U=Pi1.tolist(),S2_U=S2.tolist()),
        CK_coefficients=dict(C0=C0,C1=C1,C2=C2,A=A,linear=linear,discriminant=discriminant,
            minimum_reduced_discriminant_along_carrier_scaling=minimum_discriminant),
        controls=dict(full_source_interval=[-12.,12.],carrier_IVP_rtol=2e-11,carrier_IVP_atol=2e-13,
            carrier_max_step=.02,BVP_tolerance=2e-9,BVP_max_nodes=30000,base=base_info,
            shorter_base=short_info,first_response=first_info,second_response=second_info,
            cutoff_initial_S_max_spread=float(np.max(abs(base[:5]-shorter[:5]))),
            cutoff_initial_PiS_max_spread=float(H*np.max(abs(base[5:]-shorter[5:])))),
        audit=dict(T=audit_T.tolist(),carrier_q=carrier.sol(audit_T)[:3].T.tolist(),
            unit_q0=unit.sol(audit_T)[0].tolist(),carrier_initial=initial.tolist()),
        root_selection='Weak branchdelta(lambda)=-2lambda^3C0/[A+lambda²C1+sqrt((A+lambda²C1)²-4lambda^4C2C0)], delta~-lambda³C0/A.',
        scope='Bounded two-sided linear reference conditioning and source/CK closure, not natural/quantum state selection, exact nonlinear criticality, a joint Hamiltonian solve or a formed child.')
    if scale and (discriminant<0 or minimum_discriminant<=0):
        result.update(admitted=False,reason='No real continuous weak source/CK branch.')
    else:
        delta=-2*C0/(linear+math.sqrt(discriminant)) if scale and C0 else 0.
        S=S0+delta*S1+delta**2*S2;Pi=Pi0+delta*Pi1
        PiChi=np.array([-B*math.sqrt(.6)/12,0.,-B*math.sqrt(.6)/4]);PiChi[0]+=delta
        result.update(admitted=True,delta_PiChi=delta,
            other_root=(-linear-math.sqrt(discriminant))/(2*C2) if C2 else None,
            canonical_S_critical_U=S.tolist(),canonical_PiS_critical_U=Pi.tolist(),updated_PiChi_conformal_U=PiChi.tolist(),
            CK_polynomial_residual=C2*delta**2+linear*delta+C0,
            CK_direct_weighted_residual=mean((Pi@Y)*(S@DY))+A*mean(PiChi@Y[:3]))
    guard();result['elapsed_seconds']=time.perf_counter()-started;result['budget_seconds']=float(budget_seconds)
    return result
