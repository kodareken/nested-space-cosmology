"""Compensated scalar structure on the existing compact atlas.

Exact spherical constraints and instantaneous mass transport, not an
evolution of the isolated envelope cooling run. All masses are geometric
(G times physical mass); kappa=8*pi*G and c=1.
"""
from dataclasses import dataclass
from functools import cached_property
import math

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq

from .nsc_atlas_scalar_core import CoreAtlasEmbedding


def compact_bump(R, left, right):
    """C4 bump and analytic radial derivative, normalized to unit peak."""
    x=np.clip((np.asarray(R)-left)/(right-left),0,1)
    return 1024*x**5*(1-x)**5,5120*x**4*(1-x)**4*(1-2*x)/(right-left)


def mass_transport(R, kappa, v, beta, Pi, F, U):
    """Einstein/Hayward identities in the orthonormal (n,s) frame.

    F=dphi/dR on the initial increasing-radius chart; s(phi)=beta*F.
    fixed_radius_rate is (n-v*s/beta)(m), not unit observer proper time.
    """
    rho=(Pi*Pi+beta*beta*F*F)/2+U
    pressure=(Pi*Pi+beta*beta*F*F)/2-U
    radial=kappa*R*R*(rho-v*Pi*F)/2
    normal=kappa*R*R*(Pi*beta*beta*F-v*pressure)/2
    fixed=normal-v*radial
    return dict(density=rho,radial_pressure=pressure,
                outward_energy_current=-Pi*beta*F,
                radial_mass_derivative=radial,normal_mass_derivative=normal,
                fixed_radius_rate=fixed)


@dataclass
class CoreReservoir:
    frequency: float = 2000.
    weak_amplitude: float = .08266217132541569
    pulse_fraction: float = .02
    core_enabled: bool = True
    # The core representative is inherited from the previous finite cooling
    # record. These annuli are a declared preparation, not predicted radii.
    reservoir_outer_factor: float = 2.
    pulse_inner_factor: float = 3.
    outer_factor: float = 4.

    def __post_init__(self):
        if not (math.isfinite(self.frequency) and 1000<=self.frequency<=4000
                and math.isfinite(self.weak_amplitude) and .05<=self.weak_amplitude<=.1
                and 0<=self.pulse_fraction<=.05
                and 1<self.reservoir_outer_factor<self.pulse_inner_factor<self.outer_factor<=5):
            raise ValueError('outside the admitted finite preparation family')

    @cached_property
    def core(self):
        return CoreAtlasEmbedding(frequency=self.frequency,weak_amplitude=self.weak_amplitude,
                                  with_core=self.core_enabled)

    @property
    def reference(self):return self.core.reference

    @property
    def radius(self):return self.core.cutoff_radius

    @property
    def outer(self):return self.outer_factor*self.radius

    @property
    def H(self):return self.reference.expansions[1]

    @property
    def ambient(self):
        s=self.reference
        return s.field-s.scalar.minimum(s.parent_flux)

    @property
    def breaks(self):
        return self.radius*np.array([0.,1.,self.reservoir_outer_factor,
                                    self.pulse_inner_factor,self.outer_factor])

    def fields(self,R,depletion):
        """Field displacement, radial derivative and pulse-support flag."""
        r=self.radius
        f,F=self.core.field(R)
        b,bp=compact_bump(R,r,self.reservoir_outer_factor*r)
        f=f-depletion*self.ambient*b;F=F-depletion*self.ambient*bp
        b,bp=compact_bump(R,self.pulse_inner_factor*r,self.outer)
        phase=self.frequency*(np.asarray(R)-self.pulse_inner_factor*r)
        amp=self.pulse_fraction*self.ambient
        f=f+amp*b*np.sin(phase)
        F=F+amp*(bp*np.sin(phase)+self.frequency*b*np.cos(phase))
        outgoing=(np.asarray(R)>self.pulse_inner_factor*r)&(np.asarray(R)<self.outer)
        return f,F,outgoing

    def state(self,R,M,depletion):
        s=self.reference;kap=s.scalar.kappa;U0=s.densities[2];v=self.H*R
        beta2=1+v*v-kap*U0*R*R/3-2*M/np.maximum(R,1e-100)
        if np.any(beta2<=0):raise RuntimeError('initial radial chart lost positivity')
        beta=np.sqrt(beta2);f,F,outgoing=self.fields(R,depletion)
        # One outward characteristic initially vanishes. The potential is
        # retained, so this does not assert a purely outgoing future massive wave.
        Pi=np.where(outgoing,-beta*F,0.)
        U=s.scalar.potential(s.field+f,s.parent_flux)
        values=mass_transport(R,kap,v,beta,Pi,F,U)
        values.update(beta=beta,field_displacement=f,field_derivative=F,Pi=Pi,
                      radial_K=-self.H-kap*R*Pi*F/2,angular_K=-self.H,
                      excess_mass_derivative=values['radial_mass_derivative']-kap*U0*R*R/2)
        return values

    def integrate(self,depletion,*,rtol=2e-10,substeps=150):
        def rhs(R,y):return (float(self.state(R,y[0],depletion)['excess_mass_derivative']),)
        result=solve_ivp(rhs,(0,self.outer),(0.,),method='DOP853',rtol=rtol,atol=2e-15,
                         dense_output=True,max_step=self.radius/substeps)
        if not result.success:raise RuntimeError(result.message)
        return result

    def solve(self,*,rtol=2e-10,substeps=150):
        if self.outer>=self.reference.outer_radius:
            raise ValueError('preparation reaches the material-wall radius')
        unbalanced=self.integrate(0.,rtol=rtol,substeps=substeps).y[0,-1]
        depleted=self.integrate(1.,rtol=rtol,substeps=substeps).y[0,-1]
        if unbalanced< -1e-14 or depleted>1e-14:
            raise ValueError('declared reservoir cannot compensate this supplied core/pulse')
        if abs(unbalanced)<1e-16:amplitude=0.
        else:amplitude=brentq(lambda a:self.integrate(a,rtol=rtol,substeps=substeps).y[0,-1],
                             0,1,xtol=2e-12)
        solution=self.integrate(amplitude,rtol=rtol,substeps=substeps)
        return ReservoirSolution(self,amplitude,solution,unbalanced,depleted)


@dataclass
class ReservoirSolution:
    model: CoreReservoir
    depletion: float
    solution: object
    unbalanced_mass: float
    full_depletion_mass: float

    def state(self,R):
        return self.model.state(np.asarray(R),self.solution.sol(R)[0],self.depletion)

    def summary(self):
        c=self.model;s=c.reference;kap=s.scalar.kappa;U0=s.densities[2]
        breaks=c.breaks;mass=self.solution.sol(breaks)[0]
        def integral(func):
            return sum(quad(func,a,b,epsabs=2e-12,limit=300)[0] for a,b in zip(breaks[:-1],breaks[1:]))
        pieces={
            'gradient':integral(lambda R:kap*R*R*self.state(R)['beta']**2*self.state(R)['field_derivative']**2/4),
            'kinetic':integral(lambda R:kap*R*R*self.state(R)['Pi']**2/4),
            'potential_excess':integral(lambda R:kap*R*R*(s.scalar.potential(s.field+self.state(R)['field_displacement'],s.parent_flux)-U0)/2),
            'momentum':integral(lambda R:-kap*R*R*c.H*R*self.state(R)['Pi']*self.state(R)['field_derivative']/2),
        }
        length=integral(lambda R:1/float(self.state(R)['beta']))
        reference_length=math.asin(math.sqrt(s.cap_geometry['parent_curvature'])*c.outer)/math.sqrt(s.cap_geometry['parent_curvature'])
        proper_excess=integral(lambda R:4*math.pi*R*R*(self.state(R)['density']-U0)/self.state(R)['beta'])
        curvature=s.cap_geometry['parent_curvature']
        volume_change=integral(lambda R:4*math.pi*R*R*(1/self.state(R)['beta']-1/math.sqrt(1-curvature*R*R)))
        causal=s.parent_causal_domain();a0=causal['initial_scale_factor']
        theta_b=math.asin(c.outer/a0);eta=causal['future_conformal_distance_upper_bound']
        theta_low=theta_b+eta;theta_high=causal['initial_parent_angle']-eta
        if theta_low>=theta_high:raise ValueError('no reference future annulus remains protected')
        angular_volume=lambda x:x-math.sin(2*x)/2
        protected_volume=2*math.pi*a0**3*(angular_volume(theta_high)-angular_volume(theta_low))
        actual_parent_volume=2*math.pi*a0**3*angular_volume(causal['initial_parent_angle'])+volume_change
        # Independent proper-distance formulation evolves R and beta, using
        # Hamilton + momentum directly rather than the M_R equation.
        def rhs(ell,y):
            R,beta=y
            if R<1e-13:return (1.,0.)
            f,F,out=c.fields(R,self.depletion);Pi=-beta*F if out else 0.
            rho=(Pi*Pi+beta*beta*F*F)/2+s.scalar.potential(s.field+f,s.parent_flux)
            Lt=-c.H;Lr=-c.H-kap*R*Pi*F/2
            R3=2*kap*rho-4*Lr*Lt-2*Lt*Lt
            return (beta,(1-beta*beta)/(2*R)-R*R3/4)
        independent=solve_ivp(rhs,(0,length),(0.,1.),method='DOP853',rtol=2e-11,
                              atol=2e-13,max_step=c.radius/200)
        if not independent.success:raise RuntimeError(independent.message)
        Rb,bb=independent.y[:,-1]
        independent_M=Rb*(1+c.H*c.H*Rb*Rb-bb*bb)/2-kap*U0*Rb**3/6
        R=np.linspace(c.radius/1000,c.outer,3001);row=self.state(R);M=self.solution.sol(R)[0]
        v=c.H*R;Mp=row['excess_mass_derivative'];beta=row['beta']
        Bp=2*v*c.H-2*kap*U0*R/3-2*Mp/R+2*M/R**2
        R3=-2*Bp/R+2*(1-beta**2)/R**2
        Lr=row['radial_K'];Lt=row['angular_K']
        ham=R3+4*Lr*Lt+2*Lt*Lt-2*kap*row['density']
        mom=2*beta/R*(Lr-Lt)+kap*row['Pi']*beta*row['field_derivative']
        pulse=(R>breaks[3])&(R<breaks[4]);j=int(np.argmax(row['outward_energy_current']))
        # This compares identical initial exterior data, not global proper
        # volumes or an assertion of complete wall evolution.
        exterior_mass=float(mass[-1]);wall_trial=c.core.wall(exterior_mass)
        capacity=kap*(U0-s.scalar.vacuum(s.parent_flux))*(breaks[2]**3-breaks[1]**3)/6
        return dict(depletion_fraction=self.depletion,ambient_displacement=c.ambient,
            frequency=c.frequency,weak_amplitude=c.weak_amplitude,pulse_fraction=c.pulse_fraction,
            annulus_boundaries=breaks.tolist(),mass_at_boundaries=mass.tolist(),
            core_geometric_excess=float(mass[1]),reservoir_geometric_excess=float(mass[2]-mass[1]),
            pulse_geometric_excess=float(mass[4]-mass[3]),total_geometric_excess=exterior_mass,
            unbalanced_geometric_excess=self.unbalanced_mass,
            full_depletion_geometric_excess=self.full_depletion_mass,
            mass_integral_parts=pieces,mass_integral_gap=sum(pieces.values())-exterior_mass,
            proper_energy_excess=proper_excess,local_proper_radius=length,
            parent_proper_length_change=length-reference_length,
            parent_proper_volume_change=volume_change,
            protected_reference_annulus=dict(initial_support_angle=theta_b,
                future_conformal_distance_bound=eta,permanent_inner_angle=theta_low,
                permanent_outer_angle=theta_high,initial_protected_proper_volume=protected_volume,
                actual_initial_parent_volume=actual_parent_volume,
                fraction_of_actual_initial_parent_volume=protected_volume/actual_parent_volume,
                scope='Domain of dependence of identical smooth exterior initial data. The central perturbed ball and complete material-wall future are excluded. No bound on null propagation inside the altered core is assumed.'),
            minimum_absolute_density=float(np.min(row['density'])),vacuum_floor=s.scalar.vacuum(s.parent_flux),
            minimum_beta_squared=float(np.min(beta*beta)),
            maximum_local_compactness_excess=float(np.max(2*M/R)),
            reservoir_absolute_capacity_bound=capacity,
            minimum_reservoir_outer_radius_bound=float((breaks[1]**3+6*mass[1]/(kap*(U0-s.scalar.vacuum(s.parent_flux))))**(1/3)),
            maximum_sampled_Hamilton_residual=float(np.max(abs(ham))),
            maximum_sampled_momentum_residual=float(np.max(abs(mom))),
            independent_radius_gap=float(Rb-c.outer),independent_mass_gap=float(independent_M-exterior_mass),
            outer_wall_speed_gap=wall_trial['velocity']-s.walls[1]['velocity'],
            outer_wall_acceleration_gap=wall_trial['acceleration']-s.walls[1]['acceleration'],
            peak_outward_current=dict(radius=float(R[j]),current=float(row['outward_energy_current'][j]),
                normal_mass_rate=float(row['normal_mass_derivative'][j]),
                fixed_areal_radius_rate=float(row['fixed_radius_rate'][j]),
                fixed_areal_radius_rate_per_proper_time=float(row['fixed_radius_rate'][j]/math.sqrt(1-(v[j]/beta[j])**2))),
            maximum_fixed_radius_rate_in_pulse=float(np.max(row['fixed_radius_rate'][pulse])),
            exterior_field_jets=list(map(float,c.fields(c.outer,self.depletion)[:2])),
            vacuum_Lambda_reference_radius_squared=float(kap*s.scalar.vacuum(s.parent_flux)/s.cap_geometry['parent_curvature']),
            scope='Exact nonlinear initial constraints and instantaneous Einstein mass transport. Supplied core/pulse and solved depletion; no claim that the isolated cooling run evolved into this state, no full membrane history or observed abundance.')
