"""Core/reservoir scale separation at fixed inherited parent vacuum.

New initial preparations, not scaled copies of the evolved diffuse lens.
The existing scalar action, core shape and spherical constraint are reused.
Scaled coordinates keep tiny physical supports numerically resolved.
"""
from dataclasses import dataclass
from functools import cached_property
import math

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from .nsc_atlas_scalar_core import CoreAtlasEmbedding
from .nsc_atlas_core_reservoir import compact_bump


@dataclass
class ScaleHierarchy:
    frequency: float = 2000.
    weak_amplitude: float = .08266217132541569
    depletion: float = .5

    def __post_init__(self):
        if not (math.isfinite(self.frequency) and 2000<=self.frequency<=1e10
                and .001<=self.weak_amplitude<=.1 and 0<self.depletion<=1):
            raise ValueError('outside the declared hierarchy family')

    @cached_property
    def core(self):
        return CoreAtlasEmbedding(frequency=self.frequency,weak_amplitude=self.weak_amplitude)

    @property
    def reference(self): return self.core.reference

    @property
    def radius(self): return self.core.cutoff_radius

    @property
    def curvature(self): return self.reference.cap_geometry['parent_curvature']

    @property
    def ambient_coefficient(self):
        s=self.reference
        return math.sqrt(1-s.scalar.epsilon)*(s.parent_flux-s.fluxes[0])

    @property
    def excitation(self): return self.ambient_coefficient**2/2

    @property
    def vacuum(self):
        s=self.reference
        return s.scalar.kappa*s.scalar.vacuum(s.parent_flux)

    def shape(self,x):
        z=self.weak_amplitude;w=self.frequency;k=self.reference.scalar.kappa
        f,F=self.core.field(np.asarray(x)/(w*z));normal=math.sqrt(2/k)*z*z
        return f/normal,F/(normal*w*z)

    def core_coefficients(self,*,rtol=2e-11):
        """Exact N=N0+N1/omega+N2/omega^2 for the specified shape.

        M=zeta*N/omega, x=omega*zeta*R. The coefficient of N
        in the radial constraint is independent of omega.
        """
        z=self.weak_amplitude;k=self.reference.scalar.kappa;c=self.curvature
        d=self.ambient_coefficient
        def rhs(x,y):
            A,Ap=self.shape(x);coef=-z**4*x*Ap*Ap
            source=np.array([x*x*(A*A+z*z*Ap*Ap)/2,
                x*x*math.sqrt(k/2)*d*A/(z*z),-c*x**4*Ap*Ap/2])
            return coef*y+source
        sol=solve_ivp(rhs,(0,20),np.zeros(3),method='DOP853',rtol=rtol,
            atol=(rtol*1e-3,rtol*1e-3,rtol*1e-3),max_step=.08,dense_output=True)
        if not sol.success: raise RuntimeError(sol.message)
        return sol

    def core_data(self,*,rtol=2e-11):
        sol=self.core_coefficients(rtol=rtol);w=self.frequency;z=self.weak_amplitude
        weights=np.array([1.,1/w,1/w**2]);N=float(sol.y[:,-1]@weights)
        M=z*N/w
        half=brentq(lambda x:float(sol.sol(x)@weights)-N/2,0,20,xtol=1e-12)
        samples=np.linspace(1e-4,20,1201);n=weights@sol.sol(samples)
        beta=1-self.curvature*(samples/(w*z))**2-2*z*z*n/samples
        return dict(mass=M,half_mass_radius=half/(w*z),cutoff_radius=self.radius,
            coefficients=sol.y[:,-1].tolist(),scaled_mass=N,
            minimum_core_beta_squared=float(np.min(beta)),
            scaled_core_compactness=float(np.max(2*z*z*n/samples))),sol

    def reservoir(self,outer,core_mass,*,depletion=None,rtol=2e-11,proper=False):
        """Solve the exact constraint on x=(R-Rc)/(Rd-Rc), y=M/Mc.

        The independent proper-distance check evolves R itself, instead
        of assuming R(x) is the independent radial coordinate.
        """
        a=self.depletion if depletion is None else depletion
        Rc=self.radius;W=outer-Rc;w=self.frequency;d=self.ambient_coefficient
        k=self.reference.scalar.kappa;c=self.curvature
        if W<=0 or outer>=self.reference.outer_radius:
            raise ValueError('reservoir must fit between the core and the material wall')
        def values(x,y):
            R=Rc+W*x;b,bp=compact_bump(x,0.,1.)
            F=-a*d*bp/(w*W)
            B=1-c*R*R-2*core_mass*y/R
            if B<=0:raise RuntimeError('reservoir radial chart is not positive')
            # Cancellation-safe U-U0; do not subtract nearly equal vacua.
            du=d*d*(a*a*b*b-2*a*b)/2
            mass_R=k*R*R*(B*F*F/2+du)/2
            return B,mass_R
        if not proper:
            def rhs(x,y):return (W*values(x,y[0])[1]/core_mass,)
            sol=solve_ivp(rhs,(0,1),(1.,),method='DOP853',rtol=rtol,atol=rtol*1e-3,
                max_step=.01,dense_output=True)
            if not sol.success:raise RuntimeError(sol.message)
            x=np.linspace(0,1,801);y=sol.sol(x)[0]
            B=np.array([values(xx,yy)[0] for xx,yy in zip(x,y)])
            return dict(endpoint_mass_ratio=float(sol.y[0,-1]),minimum_beta_squared=float(np.min(B)),
                maximum_mass_ratio=float(np.max(y))),sol
        # Independent Hamilton geometry: beta^2=1-c*R^2-A*h,
        # A=2*Mc/Rc. Evolve h from R^(3)+6H^2=2*kappa*rho,
        # not from the preceding mass_R function. This scaled metric
        # deficit avoids subtracting beta^2 from one for a small core.
        A=2*core_mass/Rc
        def rhs(u,y):
            x,h=y;R=Rc+W*x;b,bp=compact_bump(x,0.,1.)
            F=-a*d*bp/(w*W);B=1-c*R*R-A*h;beta=math.sqrt(B)
            drho=B*F*F/2+d*d*(a*a*b*b-2*a*b)/2
            return beta,W*beta*(-h/R+k*R*drho/A)
        def stop(u,y):return y[0]-1
        stop.terminal=True;stop.direction=1
        sol=solve_ivp(rhs,(0,2),(0.,1.),method='DOP853',rtol=rtol,atol=rtol*1e-3,
            max_step=.01,events=stop)
        if not sol.success or len(sol.t_events[0])!=1:raise RuntimeError('proper-distance endpoint missing')
        return dict(endpoint_mass_ratio=float(sol.y_events[0][0,1]*outer/Rc),
            proper_width=float(W*sol.t_events[0][0])),sol

    def solve(self,*,rtol=2e-11,independent=True):
        core,sol=self.core_data(rtol=rtol);M=core['mass'];Rc=self.radius;k=self.reference.scalar.kappa
        ideal=(Rc**3+6*M/(k*self.excitation))**(1/3)
        # The capacity bound lies below any admissible solution. A bracket
        # expansion is a root search with a fixed depletion, not a grid scan.
        lo=ideal;hi=min(2*ideal,self.reference.outer_radius*.95)
        def residual(R):return self.reservoir(R,M,rtol=rtol)[0]['endpoint_mass_ratio']
        flo=residual(lo);fhi=residual(hi)
        while fhi>0 and hi<self.reference.outer_radius*.95:
            hi=min(2*hi,self.reference.outer_radius*.95);fhi=residual(hi)
        if flo<0 or fhi>0:raise RuntimeError('conservation root not bracketed inside the parent')
        # Solve in dimensionless radius so a tiny support does not meet
        # an unrelated absolute SI-scale root tolerance.
        ratio=brentq(lambda z:residual(z*ideal),1,hi/ideal,xtol=1e-12)
        Rd=ratio*ideal;reservoir,_=self.reservoir(Rd,M,rtol=rtol)
        proper=self.reservoir(Rd,M,rtol=rtol,proper=True)[0] if independent else None
        J=depletion_integral(self.depletion)
        asymptotic=4*self.vacuum/(k*self.ambient_coefficient**2*J)
        causal=self.reference.parent_causal_domain();a0=causal['initial_scale_factor']
        lower=math.asin(Rd/a0)+causal['future_conformal_distance_upper_bound']
        upper=causal['initial_parent_angle']-causal['future_conformal_distance_upper_bound']
        if lower>=upper:raise RuntimeError('no protected parent annulus')
        return dict(frequency=self.frequency,weak_amplitude=self.weak_amplitude,depletion=self.depletion,
            core=core,reservoir_outer_radius=Rd,reservoir_to_core_cutoff=Rd/Rc,
            ideal_capacity_lower_radius=ideal,lower_bound_margin=Rd-ideal,
            reservoir_constraint=reservoir,independent_proper_distance=proper,
            Lambda_reservoir_radius_cubed_over_core_mass=self.vacuum*Rd**3/M,
            asymptotic_Lambda_reservoir_radius_cubed_over_core_mass=asymptotic,
            core_radius_over_parent_radius=Rc/self.reference.cap_geometry['parent_maximum_radius'],
            reservoir_radius_over_parent_radius=Rd/self.reference.cap_geometry['parent_maximum_radius'],
            parent_vacuum=self.vacuum,parent_initial_density=self.reference.densities[2],
            parent_initial_curvature=self.curvature,
            parent_Lambda_reference_radius_squared=self.vacuum/self.curvature,
            protected_parent_angles=[lower,upper],
            future_conformal_distance_bound=causal['future_conformal_distance_upper_bound'],
            source='Existing coherent C2-cutoff core plus a C4 depleted ambient scalar annulus; Pi=0 and no added pulse. Outer radius is solved, depletion is supplied.')


def depletion_integral(a):
    # Exact beta integrals of x^2(2*a*b-a^2*b^2), b=1024*x^5(1-x)^5.
    first=1024*math.factorial(7)*math.factorial(5)/math.factorial(13)
    second=1024**2*math.factorial(12)*math.factorial(10)/math.factorial(23)
    return 2*a*first-a*a*second
