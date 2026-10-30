"""Causal internal-volume response to the retained, inherited scalar state.

The perturbation s=s_star+epsilon^2*c uses the same action and two-mode
source as the smooth-parent handoff.  No charged preparation, new frequency,
SI ruler or fitted amplitude is selected.  This closes the leading volume
equation; the complete second-order metric and cubic scalar remain separate.
"""
from dataclasses import dataclass
import math

import numpy as np
from scipy.integrate import quad_vec, solve_ivp

from .nsc_atlas_horizon_response import solve_source
from .nsc_atlas_parent_circulation import VOLUME, evolve, energy
from .nsc_atlas_parent_scale_selector import VolumeParent
from .nsc_atlas_inherited_charge import harmonic_sector
from .nsc_atlas_s3_polynomials import X, inner_product

EIGENVALUES=np.array([0.,8.,8.,15.,24.])


def spatial_basis():
    import sympy as s
    x,y,_,w=X
    return [s.S.One,x*x-s.Rational(1,4),y*y+w*w-s.Rational(1,2),x*y*w,
            y*y*w*w-(y*y+w*w)/8+s.Rational(1,48)]


def gram_matrix():
    # Exact algebra is used by the independent tests; rational entries here
    # are converted only after the normalized sphere pairings are evaluated.
    S=spatial_basis()
    return np.array([[float(inner_product(a,b)) for b in S] for a in S])


def source_coefficients(source,u):
    z=source.state(float(u));a=math.sqrt(2)/math.pi*z[2]
    b=math.sqrt(12)/math.pi*z[4]
    return np.array([a*a/4+b*b/24,a*a,b*b/8,2*a*b,b*b])


def horizon_means(u):
    c2=1/math.cosh(u)**2;t2=math.tanh(u)**2
    return np.array([1.,c2/3-1/4,c2/3+t2-1/2,0.,
                     c2*t2/3-(c2/3+t2)/8+1/48])


@dataclass
class VolumeResponse:
    source: object
    solution: object
    stiffness: float
    kinetic: float
    d: int

    @property
    def mass_squared(self):return self.stiffness/self.kinetic

    @property
    def nu_squared(self):return self.mass_squared/self.source.parent.H**2-2.25

    def coefficients(self,u):
        y=self.solution.sol(u);scale=math.cosh(u)**-1.5
        c=scale*y[:5]
        cu=scale*(y[5:10]-1.5*math.tanh(u)*y[:5])
        return c,cu


def solve_feedback(*,source=None,method='DOP853',rtol=2e-11,atol=2e-13):
    source=source or solve_source();model=VolumeParent()
    minimum=next(r for r in model.vacuum_roots() if r['kind']=='minimum')
    stiffness=minimum['stiffness'];K=model.K;H=source.parent.H
    nu2=stiffness/(K*H*H)-2.25
    if nu2<=0:raise ValueError('the retained response must have an oscillatory stable floor')
    factor=model.d*source.parent.mass_squared/(2*K*H*H)
    def rhs(u,z):
        cosh=math.cosh(u)
        force=factor*cosh**1.5*source_coefficients(source,u)
        return np.r_[z[5:10],force-(nu2+(EIGENVALUES+.75)/cosh**2)*z[:5]]
    sol=solve_ivp(rhs,(0.,source.endpoint_u),np.zeros(10),method=method,
        rtol=rtol,atol=atol,dense_output=True,max_step=.008)
    if not sol.success or sol.nfev>250000:
        raise RuntimeError('bounded five-harmonic volume response failed')
    return VolumeResponse(source,sol,stiffness,K,model.d)


def sample(response,u,G):
    p=response.source.parent;H=p.H;a=math.cosh(u)/H;h=H*math.tanh(u)
    c,cu=response.coefficients(u);v=H*cu
    F=source_coefficients(response.source,u);J=response.d*p.mass_squared*F/2
    mass2=response.mass_squared
    kinetic=float(v@G@v);potential=mass2*float(c@G@c)
    gradient=float(c@G@(EIGENVALUES*c))/a**2
    rho=response.kinetic*(kinetic+potential+gradient)/2
    pressure=response.kinetic*(kinetic-potential-gradient/3)/2
    proper_volume=VOLUME*a**3
    density_energy=proper_volume*rho
    exchange=proper_volume*float(J@G@v)
    work=-3*h*proper_volume*pressure
    acceleration=J/response.kinetic-3*h*v-(mass2+EIGENVALUES/a**2)*c
    direct=response.kinetic*proper_volume*(
        float(v@G@acceleration)+mass2*float(c@G@v)
        +float(v@G@(EIGENVALUES*c))/a**2-h*gradient)
    direct+=3*h*density_energy
    inverse=1/(response.stiffness+response.kinetic*EIGENVALUES/a**2)
    static=-proper_volume*float(J@G@(inverse*J))/2
    local=-proper_volume*float(J@G@J)/(2*response.stiffness)
    center=np.array([1.,-1/4,1/2,0.,-5/48])
    central=float(center@c);horizon=float(horizon_means(u)@c)
    return dict(Ht=u,proper_readout=u/H,
        volume_coefficients_over_epsilon_squared=c.tolist(),
        proper_derivative_coefficients_over_epsilon_squared=v.tolist(),
        spatial_mean_volume_shift_over_epsilon_squared=float(c[0]),
        spatial_RMS_volume_shift_over_epsilon_squared=math.sqrt(max(0.,float(c@G@c))),
        central_volume_shift_over_epsilon_squared=central,
        horizon_mean_volume_shift_over_epsilon_squared=horizon,
        central_delta_log_gauge_coupling_over_epsilon_squared=-7*central,
        horizon_mean_delta_log_gauge_coupling_over_epsilon_squared=-7*horizon,
        central_delta_log_KK_increment_over_base_mass_squared_over_epsilon_squared=-2*central,
        horizon_mean_delta_log_KK_increment_over_base_mass_squared_over_epsilon_squared=-2*horizon,
        central_delta_log_sqrt_KK_increment_over_base_mass_over_epsilon_squared=-central,
        horizon_mean_delta_log_sqrt_KK_increment_over_base_mass_over_epsilon_squared=-horizon,
        induced_volume_energy_over_epsilon_fourth=density_energy,
        induced_volume_density_over_epsilon_fourth=rho,
        induced_volume_pressure_over_epsilon_fourth=pressure,
        source_work_rate_over_epsilon_fourth=exchange,
        pressure_work_rate_over_epsilon_fourth=work,
        energy_identity_absolute_residual=abs(direct-work-exchange),
        static_relaxation_energy_contribution_over_epsilon_fourth=static,
        zero_gradient_local_relaxation_energy_over_epsilon_fourth=local,
        dynamic_response_is_static_minimizer=False)


def memory_bound(response,G):
    """Analytic future bound on the exact rescaled oscillator's quadratures.

    Source amplitudes are bounded by the original modes' nonincreasing
    positive energies. Numerical integration error is controlled separately.
    """
    source=response.source;p=source.parent;U=source.endpoint_u
    nu=math.sqrt(response.nu_squared);b=response.solution.sol(U)
    A=b[:5]*math.cos(nu*U)-b[5:10]/nu*math.sin(nu*U)
    B=b[:5]*math.sin(nu*U)+b[5:10]/nu*math.cos(nu*U)
    initial=p.initial()[[2,4]]
    qbound=np.abs(initial)*np.sqrt((source.nu_squared+np.array([3.75,8.75]))/source.nu_squared)
    a=math.sqrt(2)/math.pi*qbound[0];v=math.sqrt(12)/math.pi*qbound[1]
    Fbound=np.array([a*a/4+v*v/24,a*a,v*v/8,2*a*v,v*v])
    factor=response.d*p.mass_squared/(2*response.kinetic*p.H*p.H)
    forcing=factor*Fbound*2**1.5*(2/3)*math.exp(-1.5*U)/nu
    integral=(EIGENVALUES+.75)/nu*(1-math.tanh(U))
    tail=np.expm1(integral)*np.sqrt(A*A+B*B)+np.exp(integral)*forcing
    norm=math.sqrt(float(A@G@A+B@G@B))
    delta=math.sqrt(float(np.linalg.eigvalsh(G)[-1])*2*float(tail@tail))
    C=response.kinetic*response.mass_squared*VOLUME/(2*p.H**3)
    return dict(A=A.tolist(),B=B.tolist(),each_quadrature_future_error_upper_bound=tail.tolist(),
        phase_averaged_late_energy_estimate=C*norm*norm,
        late_energy_future_error_upper_bound=C*(2*norm*delta+delta*delta),
        bound_scope='Analytic future bound at exact endpoint data; numerical solve error is separate.')


def integrated_balance(response,G,end):
    def rates(u):
        r=sample(response,float(u),G)
        return np.array([r['pressure_work_rate_over_epsilon_fourth'],
                         r['source_work_rate_over_epsilon_fourth']])/response.source.parent.H
    work,error=quad_vec(rates,0.,end,epsabs=2e-10,epsrel=2e-11,
        points=np.linspace(0,end,33),limit=1500)
    final=sample(response,end,G)['induced_volume_energy_over_epsilon_fourth']
    return dict(endpoint_Ht=end,final_energy_over_epsilon_fourth=final,
        integrated_pressure_work_over_epsilon_fourth=float(work[0]),
        integrated_source_work_over_epsilon_fourth=float(work[1]),
        quadrature_error_estimate=float(error),absolute_residual=abs(final-float(sum(work))))


def calculation():
    import time
    begin=time.perf_counter();source=solve_source()
    main=solve_feedback(source=source)
    other=solve_feedback(source=source,method='RK45',rtol=1e-10,atol=2e-13)
    G=gram_matrix();phase=evolve();event=source.parent.H*phase.y[8,-1]
    events=[0.,.25,.5,1.,event,2.,4.,8.]
    rows=[sample(main,u,G) for u in events]
    check=np.linspace(0,8,257)
    gap=max(float(np.max(np.abs(main.coefficients(u)[0]-other.coefficients(u)[0]))) for u in check)
    nu=math.sqrt(main.nu_squared);memory=memory_bound(main,G)
    late_energy=memory['phase_averaged_late_energy_estimate']
    balances=[integrated_balance(main,G,u) for u in (event,source.endpoint_u)]
    initial_scalar=energy(source.parent.initial(),source.parent)
    sector=harmonic_sector();dM=-12*sector['zero_mode_mass_squared']-14*sector['KK_mass_squared_increment']
    c_real=main.d**2*source.parent.mass_squared**2/(8*main.stiffness)
    c_charged=dM*dM/(2*main.stiffness)
    duration=time.perf_counter()-begin
    if duration>90:raise RuntimeError('volume-feedback calculation exceeded its 90-second budget')
    return dict(schema='nsc-atlas-volume-feedback-v1',
        action=dict(kinetic=main.kinetic,volume_stiffness=main.stiffness,
            volume_mass_squared=main.mass_squared,volume_mass_over_H=math.sqrt(main.mass_squared)/source.parent.H,
            volume_late_angular_rate_per_Ht=nu,
            equation='s2ddot+3h*s2dot-a^-2 Delta(s2)+Vpp*s2/K=d*mu^2*phi1^2/(2K)',
            initial_data='s2=normal_derivative(s2)=0 on the existing retained slice',
            gauge_invariance='delta_xi s2=-xi^mu partial_mu s_star=0; metric corrections first affect this volume equation at epsilon^4',
            source='unchanged q1 Y1+q2 Y2, inherited on both sides of the same smooth horizon'),
        spectral=dict(basis=['1','X1^2-1/4','X2^2+X4^2-1/2','X1 X2 X4',
            'X2^2 X4^2-(X2^2+X4^2)/8+1/48'],eigenvalues=EIGENVALUES.tolist(),Gram=G.tolist()),
        samples=rows,shared_300_phase_pass_endpoint_Ht=event,
        integrated_energy_balances=balances,
        numerical_control=dict(methods=['DOP853','RK45'],sample_count=len(check),
            coefficient_maximum_absolute_gap=gap,nfev=[main.solution.nfev,other.solution.nfev],
            formal_perturbation_error='s=s_star+epsilon^2*s2+O(epsilon^4); no finite-amplitude remainder bound'),
        causal_memory=dict(rescaled_quadrature_and_tail=memory,
            homogeneous_volume_tail_power_in_a=-1.5,instantaneous_source_tail_power_in_a=-3.,
            phase_averaged_late_comoving_volume_energy_over_epsilon_fourth_estimate=late_energy,
            phase_averaged_volume_energy_over_initial_scalar_energy_over_epsilon_squared_estimate=late_energy/initial_scalar,
            coefficient_future_error_upper_bound=memory['late_energy_future_error_upper_bound']/initial_scalar,
            late_phase_average_is_not_an_ordinary_finite_time_window_average=True,
            selected_present_matter_fraction=False),
        Schur_interaction=dict(static_operator='L_s=-K Delta_physical+Vpp >0',
            source='J=d*mu^2*phi1^2/2 for the real scalar',
            energy='Delta E=-1/2 <J,L_s^-1 J> <0 for nonzero source',
            real_zero_gradient_quartic_hat=c_real,complex_zero_gradient_quartic_hat=c_charged,
            physical_canonical_quartic_is_ggrav_times_hat=True,
            real_and_complex_normalizations_differ=True,
            dynamic_retarded_kernel_is_instantaneous_negative_quartic=False,
            full_gravitational_or_gauged_bound_state_proved=False,
            metastable_vacuum_does_not_prove_global_stability=True),
        scale_and_charge=dict(SI_scale_selected=False,nonzero_mean_charge_generated=False,
            mass_gap_ratio_log_response='delta ln[(M1^2-mu^2)/mu^2]=-2 epsilon^2 s2',
            gauge_log_response='delta ln g_KK=-7 epsilon^2 s2',
            coupling_to_sqrt_mass_increment_ratio_log_slope=7.,
            gauge_gap_identity_unchanged='g_KK^2/(kappa_nat DeltaM_KK^2)=13/12'),
        elapsed_seconds=duration)
