"""Finite conformal data for a scalar kick and its York gravitational seed.

The background action and scalar profile are unchanged.  The freely
specified conformal TT datum is zero.  A positive numerical continuation
is not a global existence, uniqueness or state-selection theorem.
"""
from dataclasses import dataclass
import math
import time

import numpy as np
from scipy.sparse.linalg import LinearOperator,gmres

from .nsc_atlas_finite_parent_data import FiniteParentData,solve_initial_data,source_profile,positivity_lower_bound
from .nsc_atlas_kick_york import solve_york_kick
from .nsc_atlas_s3_spectral import harmonic_values,sphere_quadrature


@dataclass
class KickedParentData(FiniteParentData):
    kick_amplitude: float
    kick: object
    metadata: dict


def _inverse_source(kick,points,*,chunk=8192):
    values=np.empty(len(points));r=kick.parent.H
    for start in range(0,len(points),chunk):
        x=points[start:start+chunk];f=kick.momentum_profile(x)
        values[start:start+chunk]=(kick.squared_norm(x)+kick.parent.kappa*f*f)/r**2
    return values


def solve_kicked_initial_data(degree,*,epsilon=.25,kick_amplitude=0.,source_harmonic=1,
                             tolerance=2e-11,budget_seconds=120):
    amplitude=float(kick_amplitude);started=time.perf_counter()
    if not math.isfinite(amplitude) or abs(amplitude)>1:
        raise ValueError('this finite preparation admits |kick_amplitude|<=1')
    base=solve_initial_data(degree,epsilon=epsilon,tolerance=tolerance)
    kick=solve_york_kick(epsilon=epsilon,source_harmonic=source_harmonic,parent=base.parent)
    basis=base.basis;p=base.parent;c=base.coefficients.copy()
    profile,D=source_profile(basis.points,p)
    A=p.kappa*epsilon**2*D;B=6+p.kappa*epsilon**2*p.mass_squared/p.H**2*profile*profile
    Q=_inverse_source(kick,basis.points)
    linear=8*basis.eigenvalues+6;preconditioner=1/np.abs(8*basis.eigenvalues-24)
    stages=[]
    for fraction in ([1.] if abs(amplitude)<=.125 else [.25,.5,.75,1.]):
        square=amplitude**2*fraction
        def residual(x):
            psi=basis.evaluate(x)
            return linear*x-basis.project(A*psi+B*psi**5+square*Q*psi**-7)
        calls=0
        for iteration in range(16):
            if time.perf_counter()-started>budget_seconds:raise TimeoutError('kicked initial constraint budget exceeded')
            psi=basis.evaluate(c);R=residual(c);norm=float(np.linalg.norm(R))
            if norm<tolerance:break
            multiplier=A+5*B*psi**4-7*square*Q*psi**-8
            J=LinearOperator((len(c),len(c)),matvec=lambda v:linear*v-basis.project(multiplier*basis.evaluate(v)),dtype=float)
            M=LinearOperator(J.shape,matvec=lambda v:preconditioner*v,dtype=float)
            def counted(_):
                nonlocal calls
                calls+=1
            step,info=gmres(J,-R,M=M,rtol=5e-11,atol=2e-14,restart=min(60,len(c)),
                maxiter=10,callback=counted,callback_type='pr_norm')
            if info:raise RuntimeError('kicked initial Jacobian failed: '+str(info))
            for halve in range(14):
                trial=c+step*2.**(-halve)
                if np.min(basis.evaluate(trial))>0 and np.linalg.norm(residual(trial))<norm*(1-1e-4*2.**(-halve)):
                    c=trial;break
            else:raise RuntimeError('kicked positive-branch line search failed')
        else:raise RuntimeError('kicked initial Newton budget exceeded')
        stages.append(dict(kick_squared=float(square),Newton_steps=iteration,Krylov_iterations=calls,
            projected_residual_norm=float(np.linalg.norm(residual(c)))))
    psi=basis.evaluate(c)
    metadata=dict(kick_amplitude=amplitude,source_harmonic=source_harmonic,stages=stages,
        sampled_psi_min=float(np.min(psi)),sampled_psi_max=float(np.max(psi)),
        psi_change_coefficient_norm=float(np.linalg.norm(c-base.coefficients)),
        elapsed_seconds=time.perf_counter()-started,york=kick.diagnostics,
        branch_scope='Positive numerical continuation of the existing parent; no global uniqueness or exact finite-branch existence theorem.',
        TT_choice='Zero freely specified conformal TT datum; reference-longitudinal norm minimum, not a physical radiation selection.',
        scalar_normal_rate='Pi_chi=kick_amplitude*f_phys/psi^6; Pi_S=0',
        extrinsic_curvature='K_ij=kick_amplitude*psi^-2*Abar_ij; traceK=0')
    volume=float(basis.weights@psi**6)/(2*math.pi**2)
    reference_volume=base.continuation[-1]['proper_volume_over_reference_round_volume']
    metadata.update(solved_volume_over_reference_round=volume,
        zero_kick_volume_over_reference_round=reference_volume,
        volume_change_from_zero_kick=volume-reference_volume)
    answer=KickedParentData(p,basis,float(epsilon),c,base.continuation,amplitude,kick,metadata)
    metadata['analytic_finite_psi_lower_bound']=positivity_lower_bound(answer)
    return answer


def validate_kicked_data(data,*,chunk=8192):
    """Independent dense quadrature of the complete rational constraint."""
    degree=data.basis.degree;p=data.parent;amplitude=data.kick_amplitude
    points,weights=sphere_quadrature(5*degree+5,5*degree+5,10*degree+9)
    R2=physical2=V6=kinetic=shear=0.;minimum=math.inf;maximum=0.
    for start in range(0,len(points),chunk):
        x=points[start:start+chunk];w=weights[start:start+chunk]
        Y=harmonic_values(x,data.basis.labels);psi=Y@data.coefficients
        if np.min(psi)<=0:raise ValueError('kicked data lost positivity')
        lap=Y@(-data.basis.eigenvalues*data.coefficients)
        profile,D=source_profile(x,p);f=data.kick.momentum_profile(x);tensor2=data.kick.squared_norm(x)
        R=-8*lap+6*psi-p.kappa*data.epsilon**2*D*psi
        R-=(6+p.kappa*data.epsilon**2*p.mass_squared/p.H**2*profile**2)*psi**5
        R-=amplitude**2*(tensor2+p.kappa*f*f)/p.H**2*psi**-7
        R2+=float(w@(R*R));physical2+=float(w@(R*R*psi**-4));V6+=float(w@psi**6)
        kinetic+=float(w@(amplitude**2*f*f*psi**-6))/2
        shear+=float(w@(amplitude**2*tensor2*psi**-6))
        minimum=min(minimum,float(np.min(psi)));maximum=max(maximum,float(np.max(psi)))
    V3=2*math.pi**2
    return dict(degree=degree,kick_amplitude=amplitude,source_harmonic=data.kick.source_harmonic,
        full_conformal_RMS_estimate=math.sqrt(R2/V3),
        physical_Hamiltonian_RMS_over_H_squared_estimate=math.sqrt(physical2/V6),
        sampled_psi_min=minimum,sampled_psi_max=maximum,
        proper_volume_over_initial_round_reference=V6/V3,
        initial_normal_kinetic_over_round_vacuum_energy=kinetic/(p.vacuum*V3),
        integrated_shear_constraint_over_round_vacuum= shear/(2*p.kappa*p.vacuum*V3),
        independent_quadrature_points=len(points),
        quadrature_scope='Independent full rational-constraint estimate; inverse powers prevent polynomial-exactness claims.',
        york=data.kick.diagnostics)


def build_kicked_gh_state(system,*,data=None):
    """Represent the conformal kick in the existing GH packed state.

    The represented F=psi^4 is used in K and Pi before their own projection.
    Finite projections still need independent ADM/GH diagnostics.
    """
    if system.backend!='hopf':raise ValueError('gravitating kicks currently require the regular Hopf implementation')
    from .nsc_atlas_parent_gh import PAIRS
    from .nsc_atlas_s3_hopf import build_hopf_frame_basis
    degree=max(8,min(system.basis.degree,20))
    if data is None:
        data=solve_kicked_initial_data(degree,epsilon=system.epsilon,
            kick_amplitude=system.kick_amplitude,source_harmonic=system.kick_source)
    elif (data.basis.degree!=degree or data.epsilon!=system.epsilon or
            data.kick_amplitude!=system.kick_amplitude or data.kick.source_harmonic!=system.kick_source):
        raise ValueError('prepared kick data do not match the represented GH preparation')
    order=4*degree+system.basis.degree
    basis=build_hopf_frame_basis(system.basis.degree,(math.ceil(order/4)+2,2*math.ceil((order+1)/2)))
    psi=np.empty(len(basis.points))
    for start in range(0,len(psi),8192):psi[start:start+8192]=data.values(basis.points[start:start+8192])
    Fcoef=basis.project(psi**4);Fcoef[~np.all(system.basis.characters==[1,1],axis=1)]=0.
    F=basis.evaluate(Fcoef)
    if np.min(F)<=0:raise ValueError('represented kicked geometry lost positivity')
    C=np.zeros((24,len(system.basis.labels)));C[0,0]=-math.sqrt(2*math.pi**2)
    for k,(i,j) in enumerate(PAIRS):
        if i==j and i:C[k]=Fcoef
    profile,_=source_profile(basis.points,system.parent);C[21]=system.epsilon*basis.project(profile)
    C[21,~np.all(system.basis.characters==[-1,-1],axis=1)]=0.
    tensor=np.empty((len(F),3,3));f=data.kick.momentum_profile(basis.points)
    for start in range(0,len(F),8192):tensor[start:start+8192]=data.kick.values(basis.points[start:start+8192])
    for k,(i,j) in enumerate(PAIRS):
        if i and j:C[10+k]=basis.project(-2*system.kick_amplitude*F**-.5*tensor[:,i-1,j-1])
    C[23]=basis.project(system.kick_amplitude*f*F**-1.5)
    for i,Di in enumerate(system.basis.generators):
        k=10+PAIRS.index((0,i+1))
        C[k]=basis.project(-.5*basis.evaluate(system.r*(Di@Fcoef))/F)
    for k,char in enumerate(system.chars):C[k,~np.all(system.basis.characters==char,axis=1)]=0.
    system.initial_kick_metadata=data.metadata
    lower=Fcoef[0]/math.sqrt(2*math.pi**2)
    labels=np.array([label[0] for label in system.basis.labels])
    for level in range(1,system.basis.degree+1):
        lower-=(level+1)*np.linalg.norm(Fcoef[labels==level])/math.sqrt(2*math.pi**2)
    system.initial_kick_metadata['represented_F_lower_bound']=float(lower)
    return np.r_[system.pack(C),system.X0,[1.,0.,0.,0.],0.,0.]
