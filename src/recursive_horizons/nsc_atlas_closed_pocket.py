"""Joint closed Cauchy preparation over the unchanged spherical field engine.

Only initial data and their adapter live here. Historical neck and physical
PDE symbols keep their original sources and authentication contracts.
"""
from dataclasses import dataclass
import math
import time

import numpy as np
from scipy.optimize import brentq
from scipy.special import eval_gegenbauer

from .nsc_atlas_parent_circulation import Parent,VOLUME
from .nsc_atlas_valve_neck import _degree,zonal_rule
from .nsc_atlas_valve_spherical import SphericalValve


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


class PocketValve(SphericalValve):
    """Use a new constrained state with the original physical RHS unchanged."""
    def __init__(self,modes=64,*,initial_data,reflection_symmetry=False,budget_seconds=60.):
        if not isinstance(initial_data,PocketInitialData):
            raise ValueError('jointly constrained pocket preparation required')
        if initial_data.carrier_scale and reflection_symmetry:
            raise ValueError('directed carrier momentum requires both reflection parities')
        super().__init__(modes,epsilon=0.,initial_degree=initial_data.degree,
            reflection_symmetry=reflection_symmetry,budget_seconds=budget_seconds)
        if not math.isclose(self.sstar,initial_data.sstar,rel_tol=0,abs_tol=2e-13):
            raise ValueError('initial data and physical engine use different radion origins')
        self.initial_data=initial_data;self.parent=initial_data.parent
        self.H=self.parent.H;self.kappa=self.parent.kappa;self.epsilon=initial_data.epsilon
        self.initial=initial_data.state_values(self.u)
