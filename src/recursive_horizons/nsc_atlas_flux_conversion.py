"""Selection limits and an explicit scalar/four-form conversion candidate.

The homogeneous evolution starts AFTER specified preparation. It does not
replace P002's membrane solution or calculate a nucleation probability.
Classical numerical variables are normalized; no SI ruler is selected.
"""
from dataclasses import dataclass
from fractions import Fraction
import math

import numpy as np
from scipy.integrate import solve_ivp

from .nsc_atlas_flux_transition import FluxLaw


def sector_minimum(parent_flux, charge, *, kappa=1, bare=0, inward_only=True):
    """Exact rational arithmetic floor on f0-n*e; accessibility is not assumed.

    Inward-only means n>=0. With both orientations n is any integer.
    A transition rate can restrict this set but cannot add points to it.
    """
    f,e,k,b=map(lambda x:Fraction(str(x)),(parent_flux,charge,kappa,bare))
    if e<=0 or k<=0:
        raise ValueError("positive charge and gravity required")
    ratio=f/e
    lo=ratio.numerator//ratio.denominator
    candidates={lo,lo+1}
    if inward_only:
        candidates={max(0,n) for n in candidates}
    energies={n:b+k*(f-n*e)**2/2 for n in candidates}
    minimum=min(energies.values())
    winners=sorted(n for n,value in energies.items() if value==minimum)
    return dict(flux_residue_mod_charge=str(f-e*lo),
                nearest_generations=winners,
                residual_fluxes=[str(f-n*e) for n in winners],
                minimum_vacuum=str(minimum),
                minimum_vacuum_float=float(minimum),
                domain="n>=0" if inward_only else "integer n",
                interpretation="Kinematic minimum, not a selected terminal state or probability")


def isolated_bounce(law: FluxLaw, parent_flux):
    """Ordinary positive-cap, zero-Kottler-mass BT reference, not a P002 rate.

    Returns the leading Euclidean ACTION only, without its determinant
    prefactor. For normalized P002 variables multiply by Lref^2/kappa_phys.
    """
    lo,li=law.vacuum(parent_flux),law.vacuum(parent_flux-law.charge)
    k=law.k
    if not (lo>li>=0 and lo-li>3*k*k):
        raise ValueError("requires downward nonnegative vacua and positive-cap orientation")
    B=(lo-li)/(6*k)+k/2
    R=1/math.sqrt(li/3+B*B)
    bi,bo=B*R,(B-k)*R
    action=4*math.pi**2*k*R**3/(law.kappa*(1+bi)*(1+bo))
    return dict(parent_Lambda=lo,child_Lambda=li,radius=R,beta_inside=bi,beta_outside=bo,
                action=action,angular_residual=bi-bo-k*R,
                inside_constraint=bi*bi+li*R*R/3-1,
                outside_constraint=bo*bo+lo*R*R/3-1,
                scope="Isolated zero-mass vacuum bounce; no determinant or finite-parent nucleation rate")


@dataclass(frozen=True)
class ScalarFlux:
    """U(phi,q)=rho_b + m^2 phi^2/2 +(q+mu phi)^2/2.

    epsilon=m^2/(m^2+mu^2), omega^2=m^2+mu^2. This explicitly adds a
    propagating scalar and its positive mass; no axion protection is asserted.
    """
    epsilon: float=.25
    omega: float=20.
    bare_density: float=0.
    kappa: float=1.

    def __post_init__(self):
        if (not all(math.isfinite(x) for x in (self.epsilon,self.omega,self.bare_density,self.kappa))
                or not 0<=self.epsilon<=1 or min(self.omega,self.kappa)<=0):
            raise ValueError("finite data, 0<=epsilon<=1 and positive frequency/gravity required")

    @property
    def mu(self):
        return self.omega*math.sqrt(1-self.epsilon)

    def minimum(self,q):
        if not math.isfinite(q):
            raise ValueError("finite displacement flux required")
        return -self.mu*q/self.omega**2

    def vacuum(self,q):
        return self.bare_density+self.epsilon*q*q/2

    def potential(self,phi,q):
        return self.bare_density+self.epsilon*self.omega**2*phi*phi/2+(q+self.mu*phi)**2/2

    def wall_pressure(self,phi,outside_flux,charge):
        """Normal pressure jump for continuous phi and normal derivative.

        With no direct scalar membrane coupling, kinetic projections cancel.
        sigma*K_mean equals this, not the separately relaxed vacuum drop.
        """
        if not all(math.isfinite(x) for x in (phi,outside_flux,charge)):
            raise ValueError("finite wall field and flux data required")
        return charge*(outside_flux+self.mu*phi-charge/2)

    def prepare(self,parent_flux=1.5,charge=.5,*,parent_displacement=0.,parent_velocity=0.):
        if not all(math.isfinite(x) for x in (parent_flux,charge,parent_displacement,parent_velocity)) or charge<0:
            raise ValueError("finite parent flux and nonnegative charge required")
        child=parent_flux-charge
        old,new=self.minimum(parent_flux)+parent_displacement,self.minimum(child)
        oscillation=(parent_velocity**2+self.omega**2*(old-new)**2)/2
        released=charge*(self.epsilon*parent_flux+self.mu*parent_displacement)-charge*charge/2
        parent_energy=self.vacuum(parent_flux)+(parent_velocity**2+self.omega**2*parent_displacement**2)/2
        return dict(parent_flux=parent_flux,charge=charge,child_flux=child,
                    inherited_phi=old,child_minimum=new,displacement=old-new,
                    inherited_velocity=parent_velocity,parent_vacuum=self.vacuum(parent_flux),
                    parent_displacement=parent_displacement,parent_total_energy=parent_energy,
                    child_vacuum=self.vacuum(child),oscillation_energy=oscillation,
                    released_to_transition=released,
                    actual_wall_pressure=self.wall_pressure(old,parent_flux,charge),
                    required_external_work=max(0.,-released),
                    energy_balance_residual=parent_energy-self.vacuum(child)-oscillation-released,
                    interpretation="Conditional rapid transfer with continuous phi and velocity; boundary work and parent excitation are explicit. No nucleation or homogeneous wall passage is derived.")


def evolve_prepared(*,model=None,efolds=2.,rtol=2e-10,atol=2e-12,clock="proper"):
    """Exact homogeneous scalar-gravity evolution; dust averaging is a check only."""
    model=ScalarFlux() if model is None else model
    if not 0<efolds<=2 or clock not in ("proper","efolds"):
        raise ValueError("declared comparison is 0<N<=2 with a named clock")
    p=model.prepare()
    if p["child_vacuum"]<=0 or p["required_external_work"]>0:
        raise ValueError("reference requires positive vacuum and no unaccounted energy injection")
    w,k,rv=model.omega,model.kappa,p["child_vacuum"]
    H0=math.sqrt(k*(rv+p["oscillation_energy"])/3)
    calls=0
    if clock=="proper":
        # chi,velocity,N,H,integrated expansion work
        def rhs(t,y):
            nonlocal calls
            calls+=1
            if calls>30000:
                raise RuntimeError("finite scalar evolution budget exceeded")
            x,v,N,H,W=y
            return v,-3*H*v-w*w*x,H,-k*v*v/2,3*H*v*v
        def stop(t,y):
            return y[2]-efolds
        stop.terminal=True
        stop.direction=1
        end=efolds/math.sqrt(k*rv/3)*1.01
        result=solve_ivp(rhs,(0,end),(p["displacement"],0.,0.,H0,0.),method="DOP853",
                         events=stop,rtol=rtol,atol=atol,dense_output=True,max_step=math.pi/(8*w))
        if not result.success or len(result.t_events[0])!=1:
            raise RuntimeError("scalar evolution did not reach its finite endpoint")
    else:
        # Independently enforce Friedmann and use N as clock: chi,velocity,proper time.
        def rhs(N,y):
            nonlocal calls
            calls+=1
            if calls>30000:
                raise RuntimeError("finite scalar evolution budget exceeded")
            x,v,t=y
            H=math.sqrt(k*(rv+(v*v+w*w*x*x)/2)/3)
            return v/H,-3*v-w*w*x/H,1/H
        result=solve_ivp(rhs,(0,efolds),(p["displacement"],0.,0.),method="DOP853",
                         rtol=rtol,atol=atol,dense_output=True,max_step=.002)
        if not result.success:
            raise RuntimeError(result.message)
    return result,dict(rhs_evaluations=calls,preparation=p,initial_H=H0,
                       initial_omega_over_H=w/H0,clock=clock)
