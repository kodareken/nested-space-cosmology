"""A compact, three-region Einstein/scalar/membrane initial atlas.

Kottler formulas are reused ONLY to construct instantaneous constraint data
and wall jets: the uniform scalar is excited in the two outer regions and
cannot be continued as vacuum matter. The two walls comove with the initial
slice normal and share one middle integration constant.
"""
from dataclasses import dataclass,field
from functools import cached_property
import math

import numpy as np
from scipy.integrate import quad,solve_ivp

from .nsc_atlas_flux_conversion import ScalarFlux,sector_minimum
from .nsc_atlas_flux_transition import FluxLaw,KottlerRegion,ChargedWall
from .nsc_atlas_geometric_observer import WallGeometry,ScalarWallObserver
from .nsc_atlas_scalar_preparation import ScalarPreparation


@dataclass(frozen=True)
class SharedScalarSlice:
    scalar: ScalarFlux = field(default_factory=ScalarFlux)
    parent_flux: float = 1.5
    charge: float = .5
    tension: float = .2
    middle_mass: float = .12
    inner_radius: float = 1.
    outer_radius: float = 1.9

    def __post_init__(self):
        values=(self.parent_flux,self.charge,self.tension,self.middle_mass,self.inner_radius,self.outer_radius)
        if (not all(map(math.isfinite,values)) or min(self.charge,self.tension,self.middle_mass,self.inner_radius)<=0
                or self.outer_radius<=self.inner_radius):
            raise ValueError('finite data, positive charge/tension/middle mass and ordered positive radii required')

    @property
    def fluxes(self):
        return (self.parent_flux-2*self.charge,self.parent_flux-self.charge,self.parent_flux)

    @property
    def field(self):
        return self.scalar.minimum(self.fluxes[0])

    @property
    def densities(self):
        return tuple(self.scalar.potential(self.field,q) for q in self.fluxes)

    @property
    def k(self):
        return self.scalar.kappa*self.tension/2

    @cached_property
    def regions(self):
        return tuple(KottlerRegion(M,self.scalar.kappa*U) for M,U in zip((0.,self.middle_mass,0.),self.densities))

    @cached_property
    def reference_law(self):
        # This is an algebraic representation of U(phi0,q), not a new bare term.
        frozen_scalar_term=self.scalar.bare_density+self.scalar.epsilon*self.scalar.omega**2*self.field**2/2
        return FluxLaw(self.scalar.kappa,self.scalar.kappa*frozen_scalar_term,self.charge,self.tension)

    @cached_property
    def reference_walls(self):
        return tuple(ChargedWall(self.reference_law,self.fluxes[j+1]+self.scalar.mu*self.field,
                                self.regions[j],self.regions[j+1]) for j in (0,1))

    def wall(self,index):
        if index not in (0,1):raise ValueError('wall index must be inner0 or outer1')
        R=(self.inner_radius,self.outer_radius)[index]
        ref=self.reference_walls[index]
        beta_signed=(ref.outside.mass-ref.inside.mass)/(self.k*R**2)
        beta_signed+=(ref.outside.Lambda-ref.inside.Lambda+3*self.k**2)*R/(6*self.k)
        v2=beta_signed**2-ref.inside.f(R)
        if beta_signed-self.k*R<=0 or v2<=0:
            raise ValueError('initial wall left the expanding ordinary-normal branch')
        v=math.sqrt(v2)
        ref_geometry=ref.geometry(R,v)
        acceleration=ref_geometry['acceleration']
        potential_A,potential_D,potential_C=ref.coefficients()
        acceleration_R=10*potential_A**2/R**6+2*potential_D/R**3+potential_C
        sides=[]
        for j,beta in ((index,ref_geometry['beta_inside']),(index+1,ref_geometry['beta_outside'])):
            region=self.regions[j]
            A=(acceleration+region.f_prime(R)/2)/beta
            B=beta/R
            Bd=v/R*(A-B)
            Ad_raw=v/beta*(acceleration_R-2*region.mass/R**3-region.Lambda/3-A*A)
            # Pdot=0 is exact for the uniform Pi=0 initial state. Use its
            # cancellation-safe form, checking the independent metric jet.
            Ad=-2*Bd
            mass=R*(1+v*v-beta*beta)/2
            pressure=-self.densities[j]
            independent_acceleration=-mass/R**2-self.scalar.kappa*R*pressure/2+A*beta
            # Initial curvature from the Kottler metric derivatives; the scalar
            # stress has exactly this vacuum form only at the initial instant.
            Ruu=(-4*region.mass/R**3-2*region.Lambda/3)/2+region.f_prime(R)/R
            sides.append(dict(region=j,beta=beta,A=A,B=B,A_dot=Ad,B_dot=Bd,
                temporal_derivative_identity_residual=Ad_raw+2*Bd,
                mass=mass,mass_offset=mass-region.Lambda*R**3/6,
                independent_areal_acceleration=independent_acceleration,
                Ruu=Ruu,timelike_residual=1-(beta*beta-v*v)/region.f(R)))
        P=self.scalar.wall_pressure(self.field,self.fluxes[index+1],self.charge)
        inner,outer=sides
        meanK=(inner['A']+outer['A']+2*inner['B']+2*outer['B'])/2
        return dict(index=index,radius=R,velocity=v,expansion=v/R,acceleration=acceleration,
            pressure=P,inside=inner,outside=outer,
            angular_residual=inner['beta']-outer['beta']-self.k*R,
            temporal_residual=outer['A']-inner['A']+self.k,
            force_residual=self.tension*meanK-P,
            pressure_time_derivative=self.tension*(inner['A_dot']+2*inner['B_dot']),
            reference_energy_residual=v*v+ref.potential(R))

    @cached_property
    def walls(self):
        return (self.wall(0),self.wall(1))

    @property
    def expansions(self):
        return tuple(w['expansion'] for w in self.walls)

    @cached_property
    def cap_geometry(self):
        H1,H2=self.expansions
        ci=self.regions[0].Lambda/3-H1*H1
        co=self.regions[2].Lambda/3-H2*H2
        if ci<0:
            a=math.sqrt(-ci);Li=math.asinh(a*self.inner_radius)/a
        elif ci>0:
            a=math.sqrt(ci);Li=math.asin(a*self.inner_radius)/a
        else:Li=self.inner_radius
        if co<=0 or co*self.outer_radius**2>=1:
            raise ValueError('ordinary compact outer cap requires positive spatial curvature')
        a=math.sqrt(co)
        Lo=(math.pi-math.asin(a*self.outer_radius))/a
        return dict(inner_curvature=ci,parent_curvature=co,inner_length=Li,parent_length=Lo,
                    parent_maximum_radius=1/a)

    def middle_velocity(self,R):
        H1,H2=self.expansions
        z=(np.asarray(R)-self.inner_radius)/(self.outer_radius-self.inner_radius)
        step=6*z**5-15*z**4+10*z**3
        step_prime=30*z*z*(1-z)**2/(self.outer_radius-self.inner_radius)
        H=H1+(H2-H1)*step
        return R*H,H+R*(H2-H1)*step_prime

    def middle(self,R):
        R=np.asarray(R)
        if np.any((R<self.inner_radius)|(R>self.outer_radius)):
            raise ValueError('middle radius outside the constructed annulus')
        v,vp=self.middle_velocity(R)
        lam=self.regions[1].Lambda
        f=1-2*self.middle_mass/R-lam*R*R/3
        fp=2*self.middle_mass/R**2-2*lam*R/3
        beta=np.sqrt(f+v*v)
        Rll=fp/2+v*vp
        Lr,Lt=-vp,-v/R
        R3=-4*Rll/R+2*(1-beta*beta)/(R*R)
        Hconstraint=R3+(Lr+2*Lt)**2-Lr*Lr-2*Lt*Lt-2*self.scalar.kappa*self.densities[1]
        Lt_r=-vp/R+v/R**2
        momentum=-2*beta*Lt_r+2*beta/R*(Lr-Lt)
        return dict(radius=R,velocity=v,beta=beta,slice_radial_curvature=Lr,
                    slice_angular_curvature=Lt,spatial_scalar_curvature=R3,
                    Hamilton_residual=Hconstraint,momentum_residual=momentum)

    @cached_property
    def lengths_and_volumes(self):
        cap=self.cap_geometry
        Lm=quad(lambda r:1/float(self.middle(r)['beta']),self.inner_radius,self.outer_radius,epsabs=1e-12)[0]
        Vm=quad(lambda r:4*math.pi*r*r/float(self.middle(r)['beta']),self.inner_radius,self.outer_radius,epsabs=1e-11)[0]
        def cap_radius(l,c):
            if c>0:return math.sin(math.sqrt(c)*l)/math.sqrt(c)
            if c<0:return math.sinh(math.sqrt(-c)*l)/math.sqrt(-c)
            return l
        Vi=quad(lambda l:4*math.pi*cap_radius(l,cap['inner_curvature'])**2,0,cap['inner_length'],epsabs=1e-11)[0]
        Vo=quad(lambda l:4*math.pi*cap_radius(l,cap['parent_curvature'])**2,0,cap['parent_length'],epsabs=1e-9)[0]
        return dict(middle_length=Lm,total_length=cap['inner_length']+Lm+cap['parent_length'],
            volumes=[Vi,Vm,Vo],total_volume=Vi+Vm+Vo,
            wall_positions=[cap['inner_length'],cap['inner_length']+Lm])

    @property
    def lapse_collar_width(self):
        c=self.cap_geometry;l=self.lengths_and_volumes
        return min(.1,c['inner_length']/4,l['middle_length']/4,c['parent_length']/4)

    def lapse(self,ell):
        """Positive initial lapse realizing full temporal Israel at fixed wall labels.

        The lapse equals one on both walls and near both regular poles. This
        is a gauge extension of the initial metric jets, not a new source.
        """
        ell=np.asarray(ell)
        logN=np.zeros_like(ell,dtype=float)
        d=self.lapse_collar_width
        for w,position in zip(self.walls,self.lengths_and_volumes['wall_positions']):
            z=(ell-position)/d
            shape=np.where(abs(z)<1,z*np.maximum(0,1-z*z)**4,0.)
            A=np.where(z<=0,w['inside']['A'],w['outside']['A'])
            logN+=A*d*shape
        return np.exp(logN)

    def observer(self,index):
        w=self.walls[index];s=w['inside']
        geo=WallGeometry(w['radius'],w['velocity'],w['acceleration'],s['A'],s['B'],s['A_dot'],s['B_dot'])
        model=ScalarPreparation(self.scalar,self.fluxes[index+1],self.charge,self.tension)
        return ScalarWallObserver(model).reconstruct(geo,measured_Ruu=s['Ruu'])

    def tangent_coefficients(self,index):
        """Conditional weak-development coefficients at a separate local wall clock."""
        w=self.walls[index]
        S=self.scalar.mu*self.charge
        forces=[self.scalar.omega**2*self.field+self.scalar.mu*q for q in self.fluxes]
        mean=(forces[index]+forces[index+1])/2
        flux=S*mean/2
        side_rows=[]
        for j,key in ((index,'inside'),(index+1,'outside')):
            pressure=(mean*mean+S*S/4+forces[j]*mean)/2
            mass=self.scalar.kappa*w['radius']**2/6*(w[key]['beta']*flux-w['velocity']*pressure)
            side_rows.append(dict(region=j,normal_pressure_coefficient=pressure,
                initial_mass_offset=self.regions[j].mass,
                initial_density_subtracted_mass_coefficient=mass))
        return dict(inward_flux_coefficient=flux,pressure_drop_coefficient=flux,
                    wall_field_displacement_coefficient=-mean/2,
                    normal_field_derivative_coefficient=-S/2,sides=side_rows)

    def checks(self):
        c=self.cap_geometry;l=self.lengths_and_volumes
        R=np.linspace(self.inner_radius,self.outer_radius,501)
        mid=self.middle(R)
        endpoints=[self.middle(self.inner_radius),self.middle(self.outer_radius)]
        seam=[]
        for j,key in ((0,'outside'),(1,'inside')):
            edge=endpoints[j];w=self.walls[j]
            seam.append(dict(beta_gap=float(edge['beta']-w[key]['beta']),
                velocity_gap=float(edge['velocity']-w['velocity']),
                radial_extrinsic_gap=float(edge['slice_radial_curvature']+w['expansion']),
                angular_extrinsic_gap=float(edge['slice_angular_curvature']+w['expansion'])))
        def ode(_,y):
            v=float(self.middle_velocity(y[0])[0])
            return math.sqrt(self.regions[1].f(y[0])+v*v),4*math.pi*y[0]**2
        # Independent proper-distance integration, using beta rather than the
        # radius-coordinate quadrature which defines the target length.
        sol=solve_ivp(ode,(0,l['middle_length']),(self.inner_radius,0.),method='DOP853',rtol=2e-12,atol=2e-13)
        if not sol.success:raise RuntimeError(sol.message)
        fmin=min(self.regions[1].f(self.inner_radius),self.regions[1].f(self.outer_radius))
        if fmin<=0:raise ValueError('middle static seed collar is not positive')
        # f''<0 on this positive-mass, positive-density interval, so its
        # endpoint minimum is a lower bound for f throughout the annulus.
        bound=max(abs(w[key]['A']) for w in self.walls for key in ('inside','outside'))*self.lapse_collar_width*4096/19683
        return dict(middle_Hamilton_max=float(np.max(abs(mid['Hamilton_residual']))),
            middle_momentum_max=float(np.max(abs(mid['momentum_residual']))),
            cap_Hamilton_residuals=[6*(c['inner_curvature']+self.expansions[0]**2)-2*self.scalar.kappa*self.densities[0],
                6*(c['parent_curvature']+self.expansions[1]**2)-2*self.scalar.kappa*self.densities[2]],
            middle_f_lower_bound=fmin,seams=seam,
            independent_radius_endpoint_gap=float(sol.y[0,-1]-self.outer_radius),
            independent_volume_gap=float(sol.y[1,-1]-l['volumes'][1]),
            lapse_global_bounds=[math.exp(-bound),math.exp(bound)],
            Hamilton_surface_coefficient=4*self.k,
            required_surface_coefficient=2*self.scalar.kappa*self.tension,
            scalar_value_jump=0.,scalar_normal_jump=0.,scalar_momentum=0.)

    def parent_causal_domain(self):
        """A smooth closed-FLRW domain fixed by the actual parent initial cap.

        Its positive vacuum prevents recollapse and bounds the whole future
        conformal distance. No continuation of the material walls is assumed
        when constructing this smooth domain of dependence.
        """
        c=self.cap_geometry
        a0=1/math.sqrt(c['parent_curvature'])
        ground=self.scalar.vacuum(self.parent_flux)
        hmin=math.sqrt(self.scalar.kappa*ground/3)
        if hmin*a0<=1:
            raise ValueError('this future shielding bound requires a0*h_min>1')
        boundary=c['parent_length']/a0
        future=math.asin(1/(hmin*a0))
        protected=boundary-future
        if protected<=0:raise ValueError('parent pole is not protected by this bound')
        volume=lambda x:x-math.sin(2*x)/2
        return dict(initial_scale_factor=a0,initial_H=self.expansions[1],
            initial_parent_angle=boundary,ground_density=ground,
            inherited_vacuum_Lambda_a0_squared=self.scalar.kappa*ground*a0*a0,
            H_positive_lower_bound=math.sqrt(hmin*hmin-1/(a0*a0)),
            H_upper_bound=math.sqrt(self.scalar.kappa*self.densities[2]/3),
            future_conformal_distance_upper_bound=future,
            forever_protected_angle_lower_bound=protected,
            protected_initial_parent_volume_fraction_lower_bound=volume(protected)/volume(boundary),
            scope='Smooth Einstein-scalar evolution on the parent cap domain of dependence. A causal global completion would agree there, but the walls and their interacting future are not constructed by this bound.')

    def evolve_parent(self,*,efolds=1.,clock='proper',rtol=2e-11,atol=2e-13):
        if not 0<efolds<=1 or clock not in ('proper','efolds'):
            raise ValueError('declared parent evolution is 0<N<=1 in a named clock')
        d=self.parent_causal_domain()
        a0,H0,rv=d['initial_scale_factor'],d['initial_H'],d['ground_density']
        chi0=self.field-self.scalar.minimum(self.parent_flux)
        E0=self.scalar.omega**2*chi0*chi0/2
        omega,kap=self.scalar.omega,self.scalar.kappa
        calls=0
        if clock=='proper':
            # chi, Pi, N, H, dissipated scalar work, conformal distance
            def rhs(t,y):
                nonlocal calls
                calls+=1
                if calls>30000:raise RuntimeError('parent ODE budget exceeded')
                chi,Pi,N,H,_,_=y
                a=a0*math.exp(N)
                return (Pi,-3*H*Pi-omega**2*chi,H,1/(a*a)-kap*Pi*Pi/2,
                        3*H*Pi*Pi,1/a)
            def endpoint(t,y):return y[2]-efolds
            endpoint.terminal=True
            endpoint.direction=1
            sol=solve_ivp(rhs,(0,1.05*efolds/d['H_positive_lower_bound']),
                (chi0,0.,0.,H0,0.,0.),method='DOP853',events=endpoint,dense_output=True,
                rtol=rtol,atol=atol,max_step=.02)
            if not sol.success or not len(sol.t_events[0]):raise RuntimeError('parent endpoint not reached')
            t=float(sol.t_events[0][0]);y=sol.y_events[0][0]
            chi,Pi,N,H,work,conformal=map(float,y)
            sample=sol.sol(np.linspace(0,t,257))
            cc,pp,nn,hh,ww=sample[:5]
        else:
            # chi, Pi, H, work, conformal distance, proper time
            def rhs(N,y):
                nonlocal calls
                calls+=1
                if calls>30000:raise RuntimeError('parent ODE budget exceeded')
                chi,Pi,H,_,_,_=y
                a=a0*math.exp(N)
                if H<=0:raise RuntimeError('expanding parent branch lost')
                return (Pi/H,-3*Pi-omega**2*chi/H,(1/(a*a)-kap*Pi*Pi/2)/H,
                        3*Pi*Pi,1/(a*H),1/H)
            sol=solve_ivp(rhs,(0,efolds),(chi0,0.,H0,0.,0.,0.),method='DOP853',
                dense_output=True,rtol=rtol,atol=atol,max_step=.008)
            if not sol.success:raise RuntimeError(sol.message)
            chi,Pi,H,work,conformal,t=map(float,sol.y[:,-1]);N=efolds
            nn=np.linspace(0,efolds,257)
            sample=sol.sol(nn);cc,pp,hh,ww=sample[:4]
        a=a0*math.exp(N)
        E=(Pi*Pi+omega**2*chi*chi)/2
        energies=(pp*pp+omega**2*cc*cc)/2
        constraint=3*(hh*hh+np.exp(-2*nn)/a0**2)-kap*(rv+energies)
        balance=energies+ww-E0
        return dict(clock=clock,efolds=N,proper_time=t,scale_factor=a,H=H,
            scalar_displacement=chi,scalar_velocity=Pi,initial_excitation=E0,
            excitation_density=E,integrated_scalar_expansion_work=work,
            comoving_scalar_energy_ratio=E*math.exp(3*N)/E0,
            Omega_vacuum=kap*rv/(3*H*H),Omega_scalar_excitation=kap*E/(3*H*H),
            Omega_spatial_curvature=-1/(a*a*H*H),
            conformal_distance=conformal,
            minimum_untouched_parent_angle=d['initial_parent_angle']-conformal,
            max_Friedmann_residual=float(np.max(abs(constraint))),
            max_scalar_energy_balance_residual=float(np.max(abs(balance))),
            minimum_sampled_H=float(np.min(hh)),rhs_calls=calls,
            scope='Actual homogeneous scalar plus closed-FLRW evolution on the protected parent domain; no matter-clustering or wall-transmission solution is inferred.')


def alternative_floor_test():
    """Reject the locally equivalent B law on the new low-density region.

    Its normalized flux lattice is 3/5-n and its ground density is
    (25/32)*(3/5-n)^2. The canonical kinetic/gradient terms are nonnegative.
    """
    # The helper's coefficient is the arithmetic weight after q/e rescaling;
    # the physical Einstein coupling of candidate B remains kappa=1.
    floor=sector_minimum('3/5',1,kappa='25/16',bare=0,inward_only=False)
    required=1/32
    return dict(exact_lattice_floor=floor,normalized_density_polynomial='(25/32)*(3/5-n)^2; physical kappa=1',
        required_grandchild_density=required,
        density_gap=floor['minimum_vacuum_float']-required,
        largest_possible_Hamilton_residual=2*(required-floor['minimum_vacuum_float']),
        fixed_next_sector_density=49/32,fixed_next_sector_Hamilton_residual=2*(required-49/32),
        scope='Fixed candidate B, common normalization, unchanged charge lattice and canonical real scalar. This is a source-independent rejection for this constructed third-region geometry, not an empirical rejection of B.')
