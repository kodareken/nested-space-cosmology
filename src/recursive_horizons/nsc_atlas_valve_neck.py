"""Time-symmetric one-mode reflection neck in the retained scalar action.

This changes the declared occupation, not the action: chi=epsilon*A*X1,
kappa*A**2=24/(z+3), z=mu**2/H**2. Its ROUND-REFERENCE excitation energy
equals that of the earlier two-mode preparation. Physical energy after
solving the conformal constraint is calculated separately.

The SO(3)-invariant, reflection-even conformal factor is expanded in even
Gegenbauer U_l. No spacetime evolution or received child energy is computed.
Numerical residuals and refinement are not a certified continuum enclosure.
"""
from dataclasses import dataclass
import math
import time

import numpy as np
from scipy.optimize import brentq
from scipy.special import eval_gegenbauer

from .nsc_atlas_parent_circulation import Parent, VOLUME


def _degree(value):
    if isinstance(value,(bool,np.bool_)) or not isinstance(value,(int,np.integer)) or value<4 or value>64 or value%2:
        raise ValueError('even degree in [4,64] required')
    return int(value)


def _epsilon(value):
    if isinstance(value,(bool,np.bool_)):
        raise ValueError('finite epsilon in [0,1] required')
    value=float(value)
    if not math.isfinite(value) or not 0<=value<=1:
        raise ValueError('finite epsilon in [0,1] required')
    return value


def zonal_rule(degree,order):
    """S3 integral rule and normalized even U_l; weights sum to 2*pi**2.

    For u=cos(chi), dV=4*pi*sqrt(1-u**2)du. The Gauss-U rule is exact
    for polynomials of degree <=2*order-1 before floating-point error.
    """
    degree=_degree(degree)
    if isinstance(order,(bool,np.bool_)) or not isinstance(order,(int,np.integer)) or order<degree+1:
        raise ValueError('integer rule order >=degree+1 required')
    theta=np.arange(1,order+1)*math.pi/(order+1);u=np.cos(theta)
    weights=4*math.pi**2/(order+1)*np.sin(theta)**2
    labels=np.arange(0,degree+1,2)
    Y=np.column_stack([eval_gegenbauer(l,1,u)/math.sqrt(VOLUME) for l in labels])
    return u,weights,labels,Y


@dataclass(frozen=True)
class NeckInitialData:
    parent: object
    degree: int
    epsilon: float
    coefficients: np.ndarray
    continuation: tuple

    def values(self,u,derivative=0):
        u=np.asarray(u,dtype=float)
        if not np.all(np.isfinite(u)) or np.any(abs(u)>1):
            raise ValueError('finite u in [-1,1] required')
        if derivative not in (0,1,2):raise ValueError('derivative must be 0,1,2')
        labels=np.arange(0,self.degree+1,2)
        values=np.zeros(u.shape,dtype=float)
        factors=(1.,2.,8.)
        for coefficient,l in zip(self.coefficients,labels):
            if l>=derivative:
                values+=coefficient*factors[derivative]*eval_gegenbauer(l-derivative,1+derivative,u)/math.sqrt(VOLUME)
        return values


def solve_necks(degree,epsilons=(0.,.25,.5,.75,1.),*,max_epsilon_step=.25,tolerance=2e-12,parent=None):
    """Continue the positive Galerkin branch through requested amplitudes.

    Requested amplitudes may be unsorted; each transition is subdivided
    into steps no larger than max_epsilon_step, starting from exact vacuum.
    Returned arrays are copies, so later continuation cannot alter them.
    """
    degree=_degree(degree);requested=tuple(_epsilon(e) for e in epsilons)
    if not requested:raise ValueError('at least one requested epsilon required')
    max_epsilon_step=float(max_epsilon_step);tolerance=float(tolerance)
    if not math.isfinite(max_epsilon_step) or not 0<max_epsilon_step<=.25:
        raise ValueError('max_epsilon_step must be in (0,.25]')
    if not math.isfinite(tolerance) or tolerance<=0:raise ValueError('positive finite tolerance required')
    p=parent or Parent.retained();z=p.mass_squared/p.H**2;aa=24/(z+3)
    u,w,labels,Y=zonal_rule(degree,5*degree+5)
    linear=8*labels*(labels+2)+6
    c=np.zeros(len(labels));c[0]=math.sqrt(VOLUME);previous=0.;stages=[];out=[]
    for target in requested:
        count=max(1,math.ceil(abs(target-previous)/max_epsilon_step))
        for eps in np.linspace(previous,target,count+1)[1:]:
            a=eps*eps*aa;A=a*(1-u*u);B=6+a*z*u*u
            def residual(coef):
                psi=Y@coef
                return linear*coef-Y.T@(w*(A*psi+B*psi**5))
            for iteration in range(30):
                psi=Y@c;R=residual(c);norm=float(np.linalg.norm(R))
                if norm<tolerance:break
                J=np.diag(linear)-Y.T@((w*(A+5*B*psi**4))[:,None]*Y)
                step=np.linalg.solve(J,-R)
                for half in range(18):
                    trial=c+step*2.**-half
                    if np.min(Y@trial)>0 and np.linalg.norm(residual(trial))<norm*(1-1e-4*2.**-half):
                        c=trial;break
                else:raise RuntimeError('positive-branch neck Newton line search failed')
            else:raise RuntimeError('neck Newton iteration budget exhausted')
            stages.append(dict(epsilon=float(eps),Newton_steps=iteration,projected_residual_norm=float(np.linalg.norm(residual(c)))))
        frozen=c.copy();frozen.setflags(write=False)
        out.append(NeckInitialData(p,degree,target,frozen,tuple(dict(s) for s in stages)))
        previous=target
    return out


def neck_diagnostics(data,*,validation_order=None):
    """Constraint, stress and hemisphere accounting on one initial slice.

    The independent quadrature order integrates the polynomial conformal
    residual squared exactly before roundoff. Rational physical residuals
    are quadrature diagnostics. Positive sign margins are finite numerical
    evidence; no interval enclosure is asserted.
    """
    p=data.parent;L=data.degree;eps=data.epsilon;z=p.mass_squared/p.H**2
    aa=24/(z+3);a=eps*eps*aa
    order=8*L+9 if validation_order is None else validation_order
    if order<5*L+3:raise ValueError('validation order must be >=5*degree+3')
    u,w,labels,Y=zonal_rule(L,order);psi=Y@data.coefficients
    if np.min(psi)<=0:raise ValueError('sampled conformal factor lost positivity')
    lap=Y@(-labels*(labels+2)*data.coefficients)
    R=-8*lap+6*psi-a*(1-u*u)*psi-(6+a*z*u*u)*psi**5
    psi0=float(data.values(0.));psi2=float(data.values(0.,2))
    direct=psi0**-4*(2-4*psi2/psi0-a)
    reduced=3-(1+a/2)*psi0**-4
    grad2=p.H**2*eps**2*aa/p.kappa*psi0**-4
    V6=float(w@psi**6);volume=V6/p.H**3
    Echi=eps**2*aa/(2*p.kappa*p.H)*float(w@((1-u*u)*psi**2+z*u*u*psi**6))
    Evac=p.vacuum*volume;reference_volume=VOLUME/p.H**3;reference_energy=p.vacuum*reference_volume
    lower=(data.coefficients[0]-sum((int(l)+1)*abs(c) for l,c in zip(labels[1:],data.coefficients[1:])))/math.sqrt(VOLUME)
    return dict(degree=L,epsilon=eps,modes=len(labels),solve_quadrature_order=5*L+5,
        independent_validation_order=int(order),projected_residual_norm=data.continuation[-1]['projected_residual_norm'],
        full_conformal_PDE_RMS=math.sqrt(float(w@(R*R))/VOLUME),
        full_conformal_PDE_sample_max=float(np.max(abs(R))),
        physical_Hamiltonian_RMS_over_H_squared=math.sqrt(float(w@(R*R*psi**-4))/V6),
        physical_Hamiltonian_sample_max_over_H_squared=float(np.max(abs(R*psi**-5))),
        momentum_constraint='exactly zero: K_ij=0 and all normal scalar momenta=0',
        analytic_positivity_lower_bound=float(lower),sampled_psi_min=float(np.min(psi)),sampled_psi_max=float(np.max(psi)),
        psi_neck=psi0,psi_uu_neck=psi2,psi_u_neck=float(data.values(0.,1)),
        initial_theta_plus=0.,initial_theta_minus=0.,
        dot_theta_pm_over_H_squared_direct=direct,dot_theta_pm_over_H_squared_constraint=reduced,
        slope_identity_residual=abs(direct-reduced),
        dot_trace_expansion_over_H_squared=3.,
        dot_radial_expansion_over_H_squared=3-reduced,
        dot_each_tangential_expansion_over_H_squared=reduced/2,
        psi_neck_strict_trapping_threshold=((1+a/2)/3)**.25,
        neck_areal_radius=psi0**2/p.H,H_times_neck_areal_radius=psi0**2,
        neck_gradient_energy_density=grad2/2,neck_vacuum_energy_density=p.vacuum,
        neck_total_energy_density=p.vacuum+grad2/2,
        neck_radial_pressure=grad2/2-p.vacuum,
        neck_tangential_pressure=-grad2/2-p.vacuum,
        neck_isotropic_pressure=-grad2/6-p.vacuum,
        proper_total_volume=volume,proper_half_volume=volume/2,
        proper_volume_over_reference_round_volume=V6/VOLUME,
        physical_excitation_energy=Echi,physical_half_excitation_energy=Echi/2,
        physical_vacuum_energy=Evac,physical_half_vacuum_energy=Evac/2,
        physical_total_energy=Evac+Echi,physical_half_total_energy=(Evac+Echi)/2,
        physical_excitation_over_reference_round_vacuum_energy=Echi/reference_energy,
        supplied_round_reference_excitation_over_reference_vacuum=eps**2,
        energy_normalization='E=integral T(n,n)dV_h on the initial slice; positive vacuum included separately. Each X1 hemisphere has exactly half by reflection. These are initial regional energies, not received child energies.',
        units='Retained action units with kappa=1; no SI scale selected.',
        numerical_scope='Galerkin approximation with independent full residual and resolution comparison; not a certified continuum enclosure.')


def critical_amplitude(degree,*,parent=None,xtol=2e-12):
    """Numerical zero of the constraint-reduced neck slope in [.75,1]."""
    p=parent or Parent.retained()
    def objective(eps):
        data=solve_necks(degree,(eps,),parent=p)[0]
        aa=24/(p.mass_squared/p.H**2+3)
        return 3-(1+eps*eps*aa/2)*float(data.values(0.))**-4
    root=brentq(objective,.75,1.,xtol=xtol,rtol=1e-14)
    result=neck_diagnostics(solve_necks(degree,(root,),parent=p)[0])
    result['epsilon_crit']=root
    result['root_definition']='Zero of 3-(1+kappa epsilon^2 A^2/2)psi_neck^(-4) on the continued Galerkin branch.'
    return result


def _zonal_values(coefficients,labels,u,derivative=0):
    u=np.asarray(u,float)
    if np.any(~np.isfinite(u)) or np.any(abs(u)>1) or derivative not in (0,1,2):
        raise ValueError('finite sphere coordinates and derivative0,1,2 required')
    factor=(1.,2.,8.)[derivative]
    return sum((c*factor*eval_gegenbauer(int(l)-derivative,1+derivative,u)/math.sqrt(VOLUME)
        for c,l in zip(coefficients,labels) if l>=derivative),np.zeros_like(u))


@dataclass(frozen=True)
class PocketInitialData:
    """Joint constrained state; source amplitudes/compression remain declared.

    S is canonical and measured from sstar. No membrane or charge jump is
    inserted. The current phase is a benchmark, not a selected quantum state.
    """
    parent: object
    degree: int
    epsilon: float
    carrier_scale: float
    compression: float
    coefficients: np.ndarray
    labels: np.ndarray
    S_coefficients: np.ndarray
    S_labels: np.ndarray
    sstar: float
    carrier_A: float
    carrier_B: float
    PiS_conformal: float
    continuation: tuple

    def values(self,u,derivative=0):
        return _zonal_values(self.coefficients,self.labels,u,derivative)

    def S_values(self,u,derivative=0):
        return _zonal_values(self.S_coefficients,self.S_labels,u,derivative)

    def radial_shear_over_v(self,u):
        # div L gradF=(4/3)grad(Delta+3/R0²)F. The scalar source is
        # -PiS_conformal*gradS, with no l1 component in the even S profile.
        denominator=self.S_labels*(self.S_labels+2)-3
        scalar=self.parent.kappa*self.PiS_conformal*_zonal_values(
            self.S_coefficients/denominator,self.S_labels,u,2)
        C=-self.parent.kappa*self.carrier_A*self.carrier_B*math.sqrt(3/5)
        return scalar+C*np.asarray(u)/6

    def state_values(self,u):
        u=np.atleast_1d(np.asarray(u,float))
        if u.ndim!=1:raise ValueError('one-dimensional sphere sample points required')
        psi=self.values(u)
        if min(psi)<=0:raise ValueError('positive conformal metric required')
        v=1-u*u;cover=self.radial_shear_over_v(u);scale=psi**-6
        state=np.zeros((8,len(u)));state[0]=psi*psi/self.parent.H
        state[2]=-.5*v*cover*scale;state[3]=1.5*cover*scale
        state[4]=self.S_values(u);state[5]=self.carrier_A*u
        state[6]=self.PiS_conformal*scale
        state[7]=-self.carrier_B*math.sqrt(3/5)*(u*u-1/6)*scale
        return state

    @property
    def positivity_lower_bound(self):
        return float((self.coefficients[0]-sum((int(l)+1)*abs(c)
            for c,l in zip(self.coefficients[1:],self.labels[1:])))/math.sqrt(VOLUME))


def solve_pocket_data(degree=40,*,epsilon=.02,carrier_scale=0.,compression=0.,
                      tolerance=2e-12,budget_seconds=15.):
    """New barrier/incident-carrier preparation using the same fixed action.

    The radion baseline is the instantaneous force top at the actual chi,
    plus a declared even l2 perturbation. Source position/momentum use the
    retained carrier normalization. A declared inward radion momentum may
    break time-reflection symmetry; it is not inferred by this routine.
    """
    from .nsc_atlas_sphere_flux_completion import (
        all_volume_extrema,conditional_volume_barrier,charged_potential,combined_record)
    degree=_degree(degree)
    if (not all(math.isfinite(float(x)) for x in (epsilon,carrier_scale,compression,tolerance,budget_seconds))
        or abs(epsilon)>.05 or carrier_scale<0 or compression<0 or tolerance<=0 or budget_seconds<=0):
        raise ValueError('bounded finite preparation, positive tolerance and allowance required')
    started=time.perf_counter();reference=Parent.retained();rows=all_volume_extrema(1.5)
    top=next(x for x in rows if x['kind']=='maximum');low=next(x for x in rows if x['kind']=='minimum')
    parent=Parent(top['vacuum'],400*math.exp(-12*top['s']),reference.offset)
    A=carrier_scale*math.sqrt(24/(reference.mass_squared/reference.H**2+3))
    B=carrier_scale*reference.H**1.5*math.sqrt(12)/math.pi;PiS=-math.sqrt(84)*compression
    # Source field is even; its derivative/source is represented once and
    # used identically in the Hamiltonian, York tensor and physical data.
    u,w,Slabels,Yeven=zonal_rule(degree,6*degree+1)
    stop,_=conditional_volume_barrier(A*u)
    S=math.sqrt(84)*(stop+epsilon*(u*u-.25)-low['s'])
    Scoef=Yeven.T@(w*S);S=Yeven@Scoef
    S_u=_zonal_values(Scoef,Slabels,u,1);chi=A*u;s=low['s']+S/math.sqrt(84)
    V=charged_potential(s,combined_record().phi_star_reference+chi,1.5)
    # The scalar and carrier York tensors have opposite reflection parity.
    # Their cross term requires all conformal harmonics when both are present.
    labels=np.arange(degree+1) if compression and carrier_scale else Slabels
    Y=np.column_stack([eval_gegenbauer(int(l),1,u)/math.sqrt(VOLUME) for l in labels])
    C=-A*B*math.sqrt(3/5);v=1-u*u
    Ar=PiS*v*_zonal_values(Scoef/(Slabels*(Slabels+2)-3),Slabels,u,2)+C*u*v/6
    At=-Ar/2;PiC=-B*math.sqrt(3/5)*(u*u-1/6)
    alpha=(S_u*S_u+A*A)*v;beta=2*V/parent.H**2
    Q=(Ar*Ar+2*At*At+PiS*PiS+PiC*PiC)/parent.H**2
    linear=8*labels*(labels+2)+6
    coef=np.zeros(len(labels));coef[0]=math.sqrt(VOLUME)
    def residual(c):
        psi=Y@c
        return linear*c-Y.T@(w*(alpha*psi+beta*psi**5+Q*psi**-7))
    for iteration in range(30):
        if time.perf_counter()-started>budget_seconds:raise TimeoutError('joint pocket constraint allowance exhausted')
        R=residual(coef);norm=float(np.linalg.norm(R))
        if norm<tolerance:break
        psi=Y@coef
        jac=np.diag(linear)-Y.T@((w*(alpha+5*beta*psi**4-7*Q*psi**-8))[:,None]*Y)
        step=np.linalg.solve(jac,-R)
        for half in range(18):
            trial=coef+step*2.**-half
            if min(Y@trial)>0 and np.linalg.norm(residual(trial))<norm:
                coef=trial;break
        else:raise RuntimeError('joint pocket positive Newton branch was not admitted')
    else:raise RuntimeError('joint pocket Newton iteration allowance exhausted')
    for array in (coef,labels,Scoef,Slabels):array.setflags(write=False)
    data=PocketInitialData(parent,degree,float(epsilon),float(carrier_scale),float(compression),
        coef,labels,Scoef,Slabels,low['s'],A,B,PiS,
        (dict(Newton_steps=iteration,projected_residual_norm=norm),))
    if data.positivity_lower_bound<=0:raise RuntimeError('whole-sphere coefficient positivity control failed')
    return data
