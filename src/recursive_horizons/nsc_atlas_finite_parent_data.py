"""Finite-amplitude conformal initial constraints for the unchanged S3 source.

The amplitude is a declared benchmark, not selected by the action. Positive
conformal geometry is solved nonlinearly; no spherical averaging or scaled
perturbative metric is substituted for the inherited nonspherical field.
"""
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
import math
import multiprocessing
import time

import numpy as np
from scipy.sparse.linalg import LinearOperator,gmres
from threadpoolctl import threadpool_limits

from .nsc_atlas_parent_circulation import Parent,VOLUME
from .nsc_atlas_s3_spectral import build_surface_basis,harmonic_values,sphere_quadrature


def source_profile(points,parent):
    q=parent.initial();a=math.sqrt(2)/math.pi*q[2];b=math.sqrt(12)/math.pi*q[4]
    x,y,_,w=points.T
    f=a*x+b*y*w
    D=a*a+b*b*(y*y+w*w)-(a*x+2*b*y*w)**2
    if float(np.min(D)) < -1e-12:raise ArithmeticError('intrinsic gradient norm became negative')
    return f,np.maximum(0.,D)


def existence_certificate(parent=None):
    p=parent or Parent.retained();z=p.mass_squared/p.H**2
    # These identities follow directly from the retained equal-energy q_i.
    aa=12/(z+3);bb=72/(z+8)
    mean_D=9/(z+3)+24/(z+8)
    mean_mass=z*(aa/4+bb/24)
    maximum_mass=z*(aa+bb)**2/(4*bb)
    mean_A=6-mean_D;mean_B=6+mean_mass;max_B=6+maximum_mass
    ratio_cubed=(mean_A/6)**3*max_B/(4*mean_B)
    if not z>33:raise ValueError('the uniform finite-amplitude sufficient condition m^2>33 is not met')
    return dict(mass_squared_over_H_squared=z,
        symmetry_group_order=8,minimum_orbit_size=2,
        normalized_a_squared=aa,normalized_b_squared=bb,
        mean_A_at_epsilon1=mean_A,mean_B_at_epsilon1=mean_B,maximum_B_at_epsilon1=max_B,
        constant_trial_over_bubbling_threshold=ratio_cubed**(1/3),
        uniform_amplitude_interval=[0.,1.],
        universal_coercivity_lower_bound='38/11',
        uniform_trial_threshold_ratio_cubed_upper_bound='61/88',
        strict_numerator_polynomial='49(z+3)^2(z+8)-z(7z+26)^2=322z^2+2117z+3528>0',
        theorem='For m^2>33 and every 0<=epsilon<=1 the symmetry-invariant weighted Yamabe quotient attains below the minimum two-point-orbit concentration threshold. Rescaling its positive minimizer solves the full initial constraint.',
        unproved=['Uniqueness','Identification of every exact minimizer with the numerical continuation branch',
                  'Nonlinear future concentration or daughter formation'])


def gh_initial_geometry(psi,reference_grad_log_psi):
    """Compatible initial H_a=0 data in a static-reference orthonormal frame.

    Pi_ab=-Dbar_t g_ab and Phi_iab=Dbar_i g_ab. A future gauge driver,
    regular spatial differentiation and nonlinear GH evolution are separate.
    """
    p=np.asarray(reference_grad_log_psi,dtype=float)
    if not math.isfinite(psi) or psi<=0 or p.shape!=(3,) or not np.all(np.isfinite(p)):
        raise ValueError('positive Psi and a finite three-component reference gradient required')
    g=np.diag([-1.,psi**4,psi**4,psi**4]);Pi=np.zeros((4,4));Phi=np.zeros((3,4,4))
    Pi[0,1:]=Pi[1:,0]=2*p
    for i in range(3):Phi[i,1:,1:]=4*p[i]*psi**4*np.eye(3)
    return g,Pi,Phi


@dataclass
class FiniteParentData:
    parent: object
    basis: object
    epsilon: float
    coefficients: np.ndarray
    continuation: list

    def values(self,points):
        return harmonic_values(points,self.basis.labels)@self.coefficients


def solve_initial_data(degree,*,epsilon=1.,tolerance=2e-11):
    if not 0<=epsilon<=1 or degree<4 or degree>20:
        raise ValueError('this finite campaign admits 0 <= epsilon <= 1 and 4 <= degree <= 20')
    p=Parent.retained();basis=build_surface_basis(degree)
    f,D=source_profile(basis.points,p);m2=p.mass_squared/p.H**2
    c=np.zeros(len(basis.labels));c[0]=math.sqrt(VOLUME)
    linear=8*basis.eigenvalues+6
    preconditioner=1/np.abs(8*basis.eigenvalues-24)
    if not np.all(np.isfinite(preconditioner)):raise ArithmeticError('ell1 must be excluded by source symmetry')
    stages=[]
    for eps in (np.linspace(0,epsilon,5)[1:] if epsilon else [0.]):
        A=p.kappa*eps*eps*D;B=6+p.kappa*eps*eps*m2*f*f
        def residual(x):
            psi=basis.evaluate(x)
            return linear*x-basis.project(A*psi+B*psi**5)
        inner=0
        for iteration in range(16):
            psi=basis.evaluate(c);R=residual(c);norm=float(np.linalg.norm(R))
            if norm<tolerance:break
            multiplier=A+5*B*psi**4
            J=LinearOperator((len(c),len(c)),matvec=lambda v:linear*v-basis.project(multiplier*basis.evaluate(v)),dtype=float)
            M=LinearOperator(J.shape,matvec=lambda v:preconditioner*v,dtype=float)
            def count(_):
                nonlocal inner
                inner+=1
            step,info=gmres(J,-R,M=M,rtol=5e-11,atol=2e-14,restart=min(60,len(c)),maxiter=10,
                            callback=count,callback_type='pr_norm')
            if info:raise RuntimeError(f'finite constraint Jacobian solve failed: {info}')
            for halve in range(14):
                trial=c+step*2.**(-halve)
                if min(basis.evaluate(trial))>0 and np.linalg.norm(residual(trial))<norm*(1-1e-4*2.**(-halve)):
                    c=trial;break
            else:raise RuntimeError('positive-branch Newton step failed')
        else:raise RuntimeError('finite constraint Newton budget exceeded')
        psi=basis.evaluate(c);V6=float(basis.weights@psi**6)
        excitation=float(basis.weights@(p.H**2*psi**2*D+p.mass_squared*psi**6*f*f))*eps**2/2
        numerator=8*float(c@(basis.eigenvalues*c))+float(basis.weights@((6-A)*psi**2))
        denominator=float(basis.weights@(B*psi**6))**(1/3)
        max_B=6+eps**2*(existence_certificate(p)['maximum_B_at_epsilon1']-6)
        bubble=2**(2/3)*6*VOLUME**(2/3)/max_B**(1/3)
        stages.append(dict(epsilon=float(eps),Newton_steps=iteration,
            Krylov_iterations=inner,projected_residual_norm=float(np.linalg.norm(residual(c))),
            proper_volume_over_reference_round_volume=V6/VOLUME,
            volume_equivalent_length_over_reference_length=(V6/VOLUME)**(1/3),
            mean_excitation_density_over_vacuum=excitation/(p.vacuum*V6),
            computed_quotient_over_concentration_threshold=numerator/(denominator*bubble)))
    return FiniteParentData(p,basis,epsilon,c,stages)


def positivity_lower_bound(data):
    c=data.coefficients;labels=np.array(data.basis.labels)
    value=c[0]/math.sqrt(VOLUME)
    for l in range(1,data.basis.degree+1):
        value-=(l+1)*np.linalg.norm(c[labels[:,0]==l])/math.sqrt(VOLUME)
    return float(value)


def validate(data,*,chunk=8192):
    L=data.basis.degree;p=data.parent;eps=data.epsilon;c=data.coefficients
    # The conformal residual is a polynomial of degree <=5L+4. This rule
    # integrates its square, not only its retained Galerkin projection.
    points,weights=sphere_quadrature(5*L+5,5*L+5,10*L+9)
    R2=physicalR2=V6=E0=0.;minimum=math.inf;maximum=0.;maxR=0.
    for lo in range(0,len(points),chunk):
        x=points[lo:lo+chunk];w=weights[lo:lo+chunk]
        Y=harmonic_values(x,data.basis.labels);psi=Y@c
        lap=Y@(-data.basis.eigenvalues*c);f,D=source_profile(x,p)
        if np.min(psi)<=0:raise ArithmeticError('finite geometry lost positivity')
        R=-8*lap+6*psi-6*psi**5-p.kappa*eps**2*(D*psi+p.mass_squared/p.H**2*f*f*psi**5)
        R2+=float(w@(R*R));physicalR2+=float(w@(R*R*psi**-4))
        V6+=float(w@psi**6)
        E0+=float(w@(p.H**2*psi**2*D+p.mass_squared*psi**6*f*f))
        minimum=min(minimum,float(min(psi)));maximum=max(maximum,float(max(psi)))
        maxR=max(maxR,float(max(abs(R))))
    volume=V6/p.H**3;energy=eps**2*E0/(2*p.H**3)
    q=p.initial();a=math.sqrt(2)/math.pi*q[2];b=math.sqrt(12)/math.pi*q[4]
    fmax=b/2+a*a/(2*b)
    rho_anchor=eps**2*p.mass_squared*fmax*fmax/2
    return dict(degree=L,basis_modes=len(c),epsilon=eps,
        solve_quadrature_points=len(data.basis.points),independent_validation_points=len(points),
        full_conformal_PDE_RMS_residual=math.sqrt(R2/VOLUME),
        full_conformal_PDE_sample_max_residual=maxR,
        physical_Hamiltonian_RMS_over_H_squared=math.sqrt(physicalR2/V6),
        analytic_positivity_lower_bound=positivity_lower_bound(data),
        sampled_Psi_min=minimum,sampled_Psi_max=maximum,
        proper_volume_over_reference_round_volume=V6/VOLUME,
        volume_equivalent_length_over_reference_length=(V6/VOLUME)**(1/3),
        mean_excitation_density_over_vacuum=energy/(p.vacuum*volume),
        excitation_energy_over_reference_round_vacuum_energy=energy/(p.vacuum*VOLUME/p.H**3),
        scalar_maximum_anchor_density_over_vacuum=rho_anchor/p.vacuum,
        initial_normal_expansion_derivative_over_H_squared_at_anchor=3*(1+rho_anchor/p.vacuum),
        continuation=data.continuation,
        quadrature_scope='Conformal residual-square polynomial exactly integrated by the rule before floating-point error; physical rational-weight residual is a numerical estimate.',
        geometry_scope='Numerically solved time-symmetric finite initial slice; no evolved concentration, horizon radius, or daughter claimed.')


def level_result(degree):
    begin=time.perf_counter()
    with threadpool_limits(limits=2):
        data=solve_initial_data(degree);result=validate(data)
    result['elapsed_seconds']=time.perf_counter()-begin
    return result


def calculation():
    begin=time.perf_counter()
    with ProcessPoolExecutor(max_workers=3,mp_context=multiprocessing.get_context('spawn')) as pool:
        levels=list(pool.map(level_result,(8,12,16)))
    elapsed=time.perf_counter()-begin
    if elapsed>240:raise RuntimeError('finite initial-data campaign exhausted its 240-second budget')
    p=Parent.retained();q=p.initial();a=math.sqrt(2)/math.pi*q[2];b=math.sqrt(12)/math.pi*q[4]
    return dict(schema='nsc-atlas-finite-parent-data-v1',
        finite_amplitude_existence=existence_certificate(p),
        model=dict(source='Unchanged q1 Y1+q2 Y2, fixed action, s=s_star, all normal velocities and extrinsic curvature zero',
            amplitude=1.,amplitude_status='Declared equal-reference-energy benchmark, not an inherited amplitude-selection theorem or an observational fit',
            conformal_constraint='-8 DeltaPsi+6Psi-6Psi^5-kappa eps^2[|Df|^2 Psi+(mu^2/H^2)f^2 Psi^5]=0',
            principal_linear_operator_lower_bound=6-p.kappa*(a*a+b*b),
            principal_bound_does_not_prove_nonlinear_uniqueness=True,
            momentum_constraint_exactly_zero=True,SI_scale_selected=False),
        levels=levels,
        convergence=dict(full_residual_ratios=[levels[i]['full_conformal_PDE_RMS_residual']/levels[i+1]['full_conformal_PDE_RMS_residual'] for i in range(2)],
            fine_middle_volume_spread=abs(levels[-1]['proper_volume_over_reference_round_volume']-levels[-2]['proper_volume_over_reference_round_volume']),
            fine_middle_mean_density_spread=abs(levels[-1]['mean_excitation_density_over_vacuum']-levels[-2]['mean_excitation_density_over_vacuum'])),
        formation=dict(initial_trace_expansion=0.,
            initial_trace_acceleration='theta_dot=kappa[V_star+mu^2 eps^2 f^2/2]>0; gradients do not enter R_nn for canonical scalar at zero normal velocities',
            comparison_between_preparations_is_not_time_evolution=True,
            nonlinear_evolution_solved=False),
        covariant_GH_handoff=dict(reference='Static product -dt^2+a0^2 gamma, expressed in its orthonormal frame',
            initial_gauge_source='H_a=0',
            p='p_i=Dbar_i lnPsi',metric='g00=-1,g0i=0,gij=Psi^4 deltaij',
            Pi='Pi00=Piij=0, Pi0i=2 p_i',Phi='Phi_ijk=4 p_i Psi^4 delta_jk; other initial components zero',
            constraint='g^bc DeltaGamma^a_bc=0; allPi=0 would fail when p is nonzero',
            numerical_gradient_evaluator_and_nonlinear_GH_evolution_implemented=False),
        compute=dict(independent_resolutions=3,processes=3,BLAS_threads_per_process=2,elapsed_seconds=elapsed,
            saved_arrays=False))
