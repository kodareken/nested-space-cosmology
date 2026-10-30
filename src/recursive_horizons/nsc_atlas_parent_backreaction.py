"""Finite quadratic metric, volume and cubic scalar on the inherited S3 parent.

Gaussian-normal observers are anchored to the retained initial slice. The
complete order-epsilon^4 matter-excitation work identity subtracts the
constant vacuum stress; it is not the fourth-order Einstein solution.
"""
from dataclasses import dataclass
from concurrent.futures import ProcessPoolExecutor
import math
import multiprocessing
import time

import numpy as np
from scipy.integrate import quad, solve_ivp

from .nsc_atlas_cubic_operators import cubic_basis,numeric_operator_matrices
from .nsc_atlas_horizon_response import solve_source
from .nsc_atlas_parent_circulation import VOLUME,initial_constraint,evolve,differential_coefficient
from .nsc_atlas_s3_polynomials import inner_product,spherical_components
from .nsc_atlas_volume_feedback import EIGENVALUES,gram_matrix,solve_feedback

C1=math.sqrt(2)/math.pi
C2=math.sqrt(12)/math.pi


def matter_source(z):
    q1,p1,q2,p2=z[2:6]
    A=C1*C1*q1*q1;B=C2*C2*q2*q2;E=C1*C2*q1*q2
    Ap=C1*C1*p1*p1;Bp=C2*C2*p2*p2;Ep=C1*C2*p1*p2
    F=np.array([A/4+B/24,A,B/8,2*E,B])
    P=np.array([Ap/4+Bp/24,Ap,Bp/8,2*Ep,Bp])
    grad=np.array([3*A/4+B/3,-A,B/2,-4*E,-4*B])
    gamma=np.array([A/4+B/9,A,B/2,2*E,4*B/3])
    Hess=np.array([A/2,B/8,2*E/3,B/3])
    Q=np.r_[gamma,Hess,-q1*q2/2,q2*q2]
    momentum=np.array([C1*C1*p1*q1/2,C2*C2*p2*q2/16,
                      C1*C2*(2*p1*q2+p2*q1)/3,C2*C2*p2*q2/2])
    return F,P,grad,Q,momentum


def tensor_gram(G):
    T=np.zeros((11,11));T[:5,:5]=3*G
    T[:5,5:9]=-G[:,1:]*EIGENVALUES[1:]
    T[5:9,:5]=T[:5,5:9].T
    T[5:9,5:9]=G[1:,1:]*(EIGENVALUES[1:]*(EIGENVALUES[1:]-2))
    T[9,9]=80/(3*math.pi**4)
    T[10,10]=4/(3*math.pi**4)
    return T


def initial_metric(parent):
    c=np.array(initial_constraint(parent)['coefficients'])
    U=np.array([c[0]+c[1]/4+c[2]/2+c[3]/24,c[1],c[2]+c[3]/8,c[4],c[3]])
    return np.r_[4*U,np.zeros(6)]


@dataclass
class Backreaction:
    volume: object
    solution: object
    operators: dict
    G: np.ndarray
    GP: np.ndarray
    GT: np.ndarray

    @property
    def source(self):return self.volume.source

    def fields(self,u):
        y=self.solution.sol(u);scale=math.cosh(u)**-1.5;H=self.source.parent.H
        cubic=[]
        for i,j in ((22,30),(38,46)):
            cubic.append((scale*y[i:i+8],H*scale*(y[j:j+8]-1.5*math.tanh(u)*y[i:i+8])))
        return y[:11],H*y[11:22],cubic


def energy_and_work(volume,ops,G,GP,GT,u,y):
    p=volume.source.parent;H=p.H;z=volume.source.state(u)
    a=z[0];h=H*math.tanh(u);k=y[:11];kd=H*y[11:22]
    scale=math.cosh(u)**-1.5
    c=scale*(y[22:30]+y[38:46])
    cd=H*scale*(y[30:38]+y[46:54]-1.5*math.tanh(u)*(y[22:30]+y[38:46]))
    f=np.r_[C1*z[2],C2*z[4],np.zeros(6)]
    fd=np.r_[C1*z[3],C2*z[5],np.zeros(6)]
    s,su=volume.coefficients(u);sd=H*su
    F,P,D,Q,_=matter_source(z)
    tr=ops['trace']@k;trd=ops['trace']@kd
    rho2=(P+p.mass_squared*F+D/a**2)/2
    pressure2=(P-p.mass_squared*F-D/(3*a**2))/2
    vel_cross=float(fd@GP@cd)
    mass_cross=p.mass_squared*float(f@GP@c)
    grad_cross=-float(f@GP@(ops['laplacian']@c))/a**2
    metric_grad=float(k@GT@Q)/a**2
    exchange=volume.d*p.mass_squared*float(s@G@F)/2
    srho=volume.kinetic/2*(float(sd@G@sd)+volume.mass_squared*float(s@G@s)
                          +float(s@G@(EIGENVALUES*s))/a**2)
    spress=volume.kinetic/2*(float(sd@G@sd)-volume.mass_squared*float(s@G@s)
                          -float(s@G@(EIGENVALUES*s))/(3*a**2))
    rho4=vel_cross+mass_cross+grad_cross-metric_grad/2-exchange+srho
    pressure4=vel_cross-mass_cross-grad_cross/3+metric_grad/6+exchange+spress
    density_volume=VOLUME*a**3
    E4=density_volume*(rho4+float(tr@G@rho2)/2)
    W4=-density_volume*(3*h*(pressure4+float(tr@G@pressure2)/2)
        +float(kd@GT@Q)/(2*a**2)
        +float(trd@G@(P-D/a**2-p.mass_squared*F))/4)
    return E4,W4,dict(volume_mode_energy=srho*density_volume,
        scalar_cross_and_interaction_energy=(rho4-srho)*density_volume,
        perturbed_spatial_measure_energy=float(tr@G@rho2)*density_volume/2)


def solve_backreaction(*,volume=None,method='DOP853',rtol=2e-11,atol=2e-13):
    volume=volume or solve_feedback();source=volume.source;p=source.parent
    ops=numeric_operator_matrices();G=gram_matrix();basis=cubic_basis()
    GP=np.array([[float(inner_product(a,b)) for b in basis] for a in basis]);GT=tensor_gram(G)
    H=p.H;nu2=source.nu_squared
    lap=ops['laplacian'];mult=ops['multiplication'];metric=ops['metric_delta_laplacian']
    z0=np.zeros(55);z0[:11]=initial_metric(p)
    def rhs(u,y):
        cosh=math.cosh(u);sech2=cosh**-2;tanh=math.tanh(u)
        z=source.state(u);q=z[[2,4]];qu=z[[3,5]]/H
        F,_,_,Q,_=matter_source(z)
        k=y[:11];ku=y[11:22]
        Bp=np.r_[0.,ku[5:9]]
        kacc=np.r_[p.kappa*p.mass_squared/H**2*F+2*p.kappa*sech2*Q[:5]
            -6*tanh*ku[:5]+tanh*EIGENVALUES*Bp-(EIGENVALUES-4)*sech2*k[:5],
            -3*tanh*ku[5:9]+sech2*k[1:5]+2*p.kappa*sech2*Q[5:9],
            -3*tanh*ku[9]-p.kappa*q[0]*q[1]*sech2,
            -3*tanh*ku[10]-8*sech2*k[10]+2*p.kappa*q[1]**2*sech2]
        gravity=sech2*np.einsum('aij,i,j->a',metric,k,q)
        gravity-=.5*np.einsum('aij,i,j->a',mult,ops['trace']@ku,qu)
        s,_=volume.coefficients(u)
        breathing=volume.d*p.mass_squared/H**2*np.einsum('aij,i,j->a',mult,s,q)
        bg,bgu,bs,bsu=y[22:30],y[30:38],y[38:46],y[46:54]
        def acc(v,force):return -nu2*v-sech2*(.75*v-lap@v)+cosh**1.5*force
        _,work,_=energy_and_work(volume,ops,G,GP,GT,u,y)
        return np.r_[ku,kacc,bgu,acc(bg,gravity),bsu,acc(bs,breathing),work/H]
    sol=solve_ivp(rhs,(0.,source.endpoint_u),z0,method=method,rtol=rtol,atol=atol,
                  dense_output=True,max_step=.008)
    if not sol.success or sol.nfev>350000:
        raise RuntimeError('bounded metric/cubic response failed')
    return Backreaction(volume,sol,ops,G,GP,GT)


def constraints(response,u):
    p=response.source.parent;z=response.source.state(u);a=z[0];H=p.H
    h=H*math.tanh(u);y=response.solution.sol(u);k=y[:11];kd=H*y[11:22]
    F,P,D,Q,m= matter_source(z)
    Bdot=np.r_[0.,kd[5:9]]
    Hamiltonian=2*(EIGENVALUES-3)*k[:5]/a**2+2*h*(3*kd[:5]-EIGENVALUES*Bdot)
    Hamiltonian-=p.kappa*(P+D/a**2+p.mass_squared*F)
    momentum=kd[5:9]-kd[1:5]-p.kappa*m
    J=z[3]*z[4]-z[2]*z[5]
    vector=-5*kd[9]-p.kappa*J
    return Hamiltonian,momentum,vector


def anchor_point(parent):
    q=parent.initial();ratio=C1*q[2]/(C2*q[4])
    if not 0<ratio<1:raise ValueError('retained initial maximum changed branch')
    return np.array([ratio,math.sqrt((1-ratio*ratio)/2),0.,math.sqrt((1-ratio*ratio)/2)])


def observer_record(response,u):
    p=response.source.parent;X0=anchor_point(p)
    from .nsc_atlas_s3_polynomials import X
    vals=np.array([float(a.subs(dict(zip(X,X0)))) for a in cubic_basis()])
    z=response.source.state(u);_,_,pairs=response.fields(u)
    f=C1*X0[0]*z[2]+C2*X0[1]*X0[3]*z[4]
    fd=C1*X0[0]*z[3]+C2*X0[1]*X0[3]*z[5]
    fdd=-3*p.H*math.tanh(u)*fd-p.mass_squared*f
    fdd-=(3*C1*X0[0]*z[2]+8*C2*X0[1]*X0[3]*z[4])/z[0]**2
    w=p.reference_frequency;complex_base=complex(f,-fd/w)
    phase_speed=(complex(fd,-fdd/w)/complex_base).imag
    terms=[]
    for c,cd in pairs:
        correction=complex(float(vals@c),-float(vals@cd)/w)
        phase=(correction/complex_base).imag
        terms.append(dict(field_over_epsilon_cubed=float(vals@c),
            proper_derivative_over_epsilon_cubed=float(vals@cd),
            relative_phase_over_epsilon_squared=phase,
            H_phase_event_proper_shift_over_epsilon_squared=-p.H*phase/phase_speed))
    return dict(anchor_initial_embedding=X0.tolist(),Ht=u,proper_readout=u/p.H,
        phase_speed_over_H=phase_speed/p.H,
        metric_feedback=terms[0],volume_feedback=terms[1],
        total_H_phase_event_proper_shift_over_epsilon_squared=sum(t['H_phase_event_proper_shift_over_epsilon_squared'] for t in terms),
        interpretation='Normal geodesic through the initial scalar maximum, proper duration, and supplied scalar quadrature clock; conditional comparison to the unchanged leading wave.')


def harmonic_grams():
    basis=cubic_basis();sectors={}
    pieces=[spherical_components(p) for p in basis]
    for l in (1,2,3,4,5,6):
        vectors=[p.get(l,0) for p in pieces]
        sectors[l]=np.array([[float(inner_product(a,b)) for b in vectors] for a in vectors])
    return sectors


def cross_observable_record(response,observer):
    from types import SimpleNamespace
    p=response.source.parent
    def previous_guide_state(t):
        u=p.H*t;z=response.source.state(u)
        # Cdot=-beta and C(0)=0. The original held guides are transported
        # into normal coordinates by x_old=x_normal+epsilon²*C*W.
        z[7]=-response.solution.sol(u)[9]
        return z
    light=differential_coefficient(p,SimpleNamespace(sol=previous_guide_state))
    clock=observer['total_H_phase_event_proper_shift_over_epsilon_squared']
    return dict(H_guided_differential_light_delay_over_epsilon_squared=light['H_delay_over_epsilon_squared'],
        H_scalar_phase_event_shift_over_epsilon_squared=clock,
        scalar_clock_shift_over_guided_light_delay=clock/light['H_delay_over_epsilon_squared'],
        amplitude_and_SI_ruler_cancel=True,
        same_event_required=False,
        guide_protocol='Existing radius1/8, separationpi/6, phase-averaged opposite directions; transport the held-guide worldlines under x_old=x_normal+epsilon²*C*W',
        scope='Conditional two-protocol relation for the same fixed two-mode preparation; not a universal source-independent constant or observed sky comparison.')


def independent_solution(volume):
    # Return only the picklable SciPy trajectory, not immutable operator maps.
    return solve_backreaction(volume=volume,method='RK45',rtol=8e-11,atol=3e-13).solution


def calculation():
    begin=time.perf_counter();volume=solve_feedback()
    with ProcessPoolExecutor(max_workers=1,mp_context=multiprocessing.get_context('spawn')) as pool:
        pending=pool.submit(independent_solution,volume)
        main=solve_backreaction(volume=volume)
        other=Backreaction(volume,pending.result(),main.operators,main.G,main.GP,main.GT)
    event=main.source.parent.H*evolve().y[8,-1]
    points=np.linspace(0,main.source.endpoint_u,257)
    residuals=[constraints(main,float(u)) for u in points]
    hmax=max(float(np.max(np.abs(r[0]))) for r in residuals)
    mmax=max(float(np.max(np.abs(r[1]))) for r in residuals)
    vmax=max(abs(r[2]) for r in residuals)
    k_gap=0.;c_gap=0.
    for u in points:
        k,kd,pairs=main.fields(float(u));k2,kd2,pairs2=other.fields(float(u))
        k_gap=max(k_gap,float(np.max(np.abs(k-k2))))
        c_gap=max(c_gap,max(float(np.max(np.abs(a[0]-b[0]))) for a,b in zip(pairs,pairs2)))
    E0=energy_and_work(volume,main.operators,main.G,main.GP,main.GT,0.,main.solution.y[:,0])[0]
    balances=[]
    for u in (event,main.source.endpoint_u):
        E,W,parts=energy_and_work(volume,main.operators,main.G,main.GP,main.GT,u,main.solution.sol(u))
        integrated,error=quad(lambda v:energy_and_work(volume,main.operators,main.G,main.GP,main.GT,
            v,main.solution.sol(v))[1]/main.source.parent.H,0.,u,
            epsabs=5e-9,epsrel=5e-11,points=np.linspace(0,u,65),limit=2000)
        balances.append(dict(Ht=u,initial_excitation_energy_coefficient=E0,
            final_excitation_energy_coefficient=E,parts=parts,
            integrated_work_from_evolution=float(main.solution.sol(u)[54]),
            independent_integrated_work=integrated,quadrature_error_estimate=error,
            conservation_absolute_residual=abs(E-E0-integrated)))
    grams=harmonic_grams();harmonics=[]
    for u in (event,2.,4.,8.):
        _,_,pairs=main.fields(u);c=pairs[0][0]+pairs[1][0]
        row=dict(Ht=u,cubic_field_L2_squared_by_degree={str(l):max(0.,float(c@G@c)*VOLUME) for l,G in grams.items()})
        row['scope']='Order-epsilon^6 squared-field coefficient; not a separately conserved mode energy.'
        harmonics.append(row)
    elapsed=time.perf_counter()-begin
    if elapsed>180:raise RuntimeError('bounded coupled perturbative calculation exceeded 180 seconds')
    observer=observer_record(main,event);observer_control=observer_record(other,event)
    observer['two_integrator_absolute_H_phase_shift_gap']=abs(observer['total_H_phase_event_proper_shift_over_epsilon_squared']
        -observer_control['total_H_phase_event_proper_shift_over_epsilon_squared'])
    return dict(schema='nsc-atlas-parent-backreaction-v1',
        gauge='Gaussian normal to the retained initial slice; three spatial directions with offset evolution',
        inherited_data=dict(metric_initial='k=4u gamma from the existing nonlinear constraint-data branch; normal derivative zero',
            volume_initial='s2=sdot2=0',scalar_initial='phi3=normal_derivative(phi3)=0',new_source_or_fit=False),
        equations=dict(metric_orders='All scalar, compensated vector and TT metric sectors at epsilon^2',
            volume_order='Existing sourced s2 at epsilon^2',scalar_order='Complete phi3 at epsilon^3, separated by metric/volume drive',
            matter_energy_order='All excitation-stress and spatial-measure terms at epsilon^4, subtracting constant vacuum V_star',
            fourth_order_Einstein_solution=False,finite_amplitude_remainder_bound=False),
        numerical_controls=dict(Hamiltonian_max_absolute_residual=hmax,scalar_momentum_max_absolute_residual=mmax,
            vector_momentum_max_absolute_residual=vmax,metric_two_integrator_max_gap=k_gap,
            cubic_two_integrator_max_gap=c_gap,methods=['DOP853','RK45'],nfev=[main.solution.nfev,other.solution.nfev],
            independent_integrations_run_concurrently=True),
        energy_balances=balances,observer=observer,harmonic_generation=harmonics,
        cross_observable=cross_observable_record(main,observer),
        complete_finite_amplitude_concentration_or_child_formation=False,
        SI_scale_selected=False,charge_asymmetry_selected=False,elapsed_seconds=elapsed)
