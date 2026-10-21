"""Self-gravitating real-scalar envelopes and a sourced atlas initial slice.

The Schrödinger--Newton equation is a leading classical envelope reduction,
not an extra complex field or a claim of quantum state reduction. Local
cooling in a de Sitter limit and exact spherical initial constraints have
separate domains; their combination is not a full membrane formation history.
"""
from dataclasses import dataclass
from functools import cached_property
import math
import time

import numpy as np
from scipy.fft import dst,idst
from scipy.integrate import quad,solve_bvp,solve_ivp
from scipy.optimize import brentq

from .nsc_atlas_flux_conversion import ScalarFlux
from .nsc_atlas_flux_transition import KottlerRegion,ChargedWall
from .nsc_atlas_shared_scalar_slice import SharedScalarSlice


@dataclass
class CoreProfile:
    solution: object
    scale: float
    radius: float

    def fields(self,r):
        r=np.asarray(r,dtype=float)
        if np.any(r<0):raise ValueError('negative radius')
        bounded=np.minimum(r,self.radius)
        y=self.solution.sol(self.scale*bounded)
        S=self.scale**2*y[0];Sd=self.scale**3*y[1]
        s=self.scale**2*y[2];sd=self.scale**3*y[3]
        outside=r>self.radius
        if np.any(outside):
            R=self.radius;tail_power=self.mass_norm/2-1
            factor=np.exp(-(r-R))*(np.maximum(r,R)/R)**tail_power
            S=np.where(outside,S*factor,S)
            Sd=np.where(outside,S*(-1+tail_power/np.maximum(r,R)),Sd)
            s=np.where(outside,-1+self.mass_norm/np.maximum(r,R),s)
            sd=np.where(outside,-self.mass_norm/np.maximum(r,R)**2,sd)
        return S,Sd,s,sd

    @cached_property
    def mass_norm(self):
        return -self.radius**2*self.scale**3*float(self.solution.sol(self.scale*self.radius)[3])

    def summary(self):
        R=self.radius;N=quad(lambda r:r*r*self.fields(r)[0]**2,0,R,epsabs=1e-11)[0]
        T=quad(lambda r:.5*r*r*self.fields(r)[1]**2,0,R,epsabs=1e-11)[0]
        W=quad(lambda r:-.25*r*r*(self.fields(r)[2]+1)*self.fields(r)[0]**2,0,R,epsabs=1e-11)[0]
        radii={str(f):brentq(lambda r:quad(lambda x:x*x*self.fields(x)[0]**2,0,r,epsabs=1e-11)[0]-f*N,0,R)
               for f in (.5,.95,.99)}
        return dict(central_amplitude=float(self.fields(0)[0]),radial_mass_norm=self.mass_norm,
            integrated_mass=N,mass_integral_gap=N-self.mass_norm,kinetic=T,gravitational_energy=W,
            virial_residual=2*T+W,eigenvalue_identity=T+2*W+N/2,
            mass_fraction_radii=radii,
            isolated_mass_radius_frequency_squared={f:self.mass_norm*value/2 for f,value in radii.items()},
            outer_radius=R,
            maximum_collocation_residual=float(np.max(self.solution.rms_residuals)))


def solve_core(*,tolerance=1e-8,outer=25.):
    if not 15<=outer<=30:raise ValueError('admitted stationary-core interval is 15..30')
    r=np.linspace(0,outer,401)
    S=np.exp(-r*r/6);s=-1+3.5/np.sqrt(r*r+4)
    guess=np.array([S,-r*S/3,s,-3.5*r/(r*r+4)**1.5])
    def rhs(r,y,p):return np.array([y[1],-y[2]*y[0],y[3],-y[0]**2])
    def bc(a,b,p):
        k=np.exp(p[0]/2)
        return np.array([a[0]-1,a[1],a[3],b[2]+outer*b[3]+k*k,
            b[1]+(k+1/outer+outer*b[3]/(2*k))*b[0]])
    result=solve_bvp(rhs,bc,r,guess,p=[0.],S=np.diag([0.,-2.,0.,-2.]),
        tol=tolerance,max_nodes=6000)
    if not result.success:raise RuntimeError('stationary scalar core: '+result.message)
    scale=math.exp(-result.p[0]/2)
    profile=CoreProfile(result,scale,outer/scale)
    if np.min(profile.fields(np.linspace(0,profile.radius,1001))[0])<=0:
        raise RuntimeError('stationary solve left the nodeless branch')
    return profile


@dataclass(frozen=True)
class CoolingRun:
    shape: str = 'gaussian'
    points: int = 2047
    outer: float = 64.
    step: float = .02
    end: float = 600.
    h: float = .015309310892394859
    width: float = 6.
    gravity: bool = True
    absorber: bool = True


def evolve_core(config=CoolingRun(),*,retain=False):
    """Unitary split evolution with an explicitly accounted outer absorber."""
    start=time.process_time();wall=time.perf_counter();p=solve_core()
    c=config;n=c.points;L=c.outer;dt=c.step
    if (c.shape not in ('gaussian','supergaussian','stationary') or n<255 or n>4095
            or not 24<=L<=128 or not 0<dt<=.05 or not 0<c.end<=600 or not 0<=c.h<=.02):
        raise ValueError('cooling calculation outside admitted finite domain')
    steps=round(c.end/dt)
    if not math.isclose(steps*dt,c.end,abs_tol=1e-10):raise ValueError('endpoint must be a whole step')
    dr=L/(n+1);r=np.arange(1,n+1)*dr;k=np.pi*np.arange(1,n+1)/L
    if c.shape=='stationary':psi=p.fields(r)[0]
    else:psi=np.exp(-.5*(r/c.width)**(2 if c.shape=='gaussian' else 4))
    u=np.asarray(r*psi,dtype=complex);u*=math.sqrt(p.mass_norm/(dr*np.sum(abs(u)**2)))
    initial_u=u.copy();kinetic=np.exp(-.5j*k*k*dt)
    gamma=.2*np.maximum(0,(r-.7*L)/(.3*L))**4 if c.absorber else np.zeros_like(r)
    absorb=np.exp(-gamma*dt);tide=-.5*c.h*c.h*r*r
    def potential(u):
        density=abs(u)**2
        if not c.gravity:return np.zeros_like(r)
        mass=dr*np.cumsum(density);tail=dr*np.cumsum(density/r)
        # Symmetric discrete shell kernel -dr/(2 max(r_i,r_j)).
        return -.5*(mass/r+tail[-1]-tail)
    def energy(u,V=None):
        density=abs(u)**2
        if V is None:V=potential(u)
        kinetic_energy=.5*dr*np.sum(k*k*abs(dst(u,type=1,norm='ortho'))**2)
        return float(kinetic_energy+.5*dr*np.sum(V*density)+dr*np.sum(tide*density))
    cut=20.;j=int(round(cut/dr))-1
    if j<2 or j+2>=n or not math.isclose(r[j],cut,abs_tol=1e-12):
        raise ValueError('admitted mesh must place the independent flux cut at radius20')
    inside=lambda u:float(dr*(np.sum(abs(u[:j])**2)+.5*abs(u[j])**2))
    def flux(u):
        derivative=(u[j-2]-8*u[j-1]+8*u[j+1]-u[j+2])/(12*dr)
        return float(np.imag(np.conj(u[j])*derivative))
    initial_mass=float(dr*np.sum(abs(u)**2));initial_energy=energy(u);initial_inside=inside(u)
    removed_mass=0.;removed_energy=0.;outward=0.;oldflux=flux(u)
    snapshots=[];average=np.zeros_like(r);average_count=0
    stride=max(1,steps//120);sample_stride=max(1,round(.2/dt))
    for iteration in range(steps):
        u*=np.exp(-.5j*dt*(potential(u)+tide))
        u=idst(dst(u,type=1,norm='ortho')*kinetic,type=1,norm='ortho')
        V=potential(u);u*=np.exp(-.5j*dt*(V+tide))
        before=float(dr*np.sum(abs(u)**2))
        if c.absorber:
            Ebefore=energy(u,V);u*=absorb
            removed_energy+=Ebefore-energy(u)
        removed_mass+=before-float(dr*np.sum(abs(u)**2))
        newflux=flux(u);outward+=dt*(oldflux+newflux)/2;oldflux=newflux
        tau=(iteration+1)*dt
        if tau>=c.end-100 and iteration%sample_stride==0:
            average+=abs(u/r)**2;average_count+=1
        if (iteration+1)%stride==0 or iteration+1==steps:
            amp=abs(u/r)
            snapshots.append(dict(tau=tau,mass=float(dr*np.sum(abs(u)**2)),
                inner_mass=inside(u),central_amplitude=float((4*amp[0]-amp[1])/3),
                mass_balance=float(dr*np.sum(abs(u)**2))+removed_mass-initial_mass,
                energy_balance=energy(u)+removed_energy-initial_energy))
    amplitude=abs(u/r);mean_density=average/max(1,average_count)
    core=r<=10
    core_mass=float(dr*np.sum((r*r*mean_density)[core]))
    # Use mass, not a fit to the profile, to fix the reference scale.
    def retained_norm(scale):
        return float(dr*np.sum((r[core]*scale**2*p.fields(scale*r[core])[0])**2))
    scale=brentq(lambda z:retained_norm(z)-core_mass,.1,2)
    target=scale**2*p.fields(scale*r)[0]
    shape_error=float(np.sqrt(np.sum(r[core]**2*(np.sqrt(mean_density[core])-target[core])**2)
                              /np.sum(r[core]**2*target[core]**2)))
    initial_core_norm=float(dr*np.sum(abs(initial_u[core])**2))
    initial_scale=brentq(lambda z:retained_norm(z)-initial_core_norm,.1,2)
    initial_target=initial_scale**2*p.fields(initial_scale*r)[0]
    initial_gap=float(np.sqrt(np.sum(r[core]**2*(abs(initial_u[core]/r[core])-initial_target[core])**2)
        /np.sum(r[core]**2*initial_target[core]**2)))
    result=dict(config=c.__dict__,initial_mass=initial_mass,final_mass=float(dr*np.sum(abs(u)**2)),
        removed_mass=removed_mass,initial_energy=initial_energy,final_energy=energy(u),removed_energy=removed_energy,
        mass_balance=float(dr*np.sum(abs(u)**2))+removed_mass-initial_mass,
        energy_balance=energy(u)+removed_energy-initial_energy,
        maximum_sampled_energy_balance=max(abs(row['energy_balance']) for row in snapshots),
        outward_flux_at20=outward,inside20_mass_change=inside(u)-initial_inside,
        independent_flux_balance=inside(u)-initial_inside+outward,
        averaging_interval=[max(0,c.end-100),c.end],mean_core_mass_inside10=core_mass,
        mass_selected_scale=scale,initial_profile_gap=initial_gap,late_profile_gap=shape_error,
        final_central_amplitude=float((4*amplitude[0]-amplitude[1])/3),
        mean_central_density=float((4*mean_density[0]-mean_density[1])/3),
        cpu_seconds=time.process_time()-start,wall_seconds=time.perf_counter()-wall)
    if c.shape=='gaussian' and not c.gravity and c.h==0:
        A=math.sqrt(initial_mass/(math.sqrt(math.pi)*c.width**3/4))
        analytic=A*c.width**3/(c.width**2+1j*c.end)**1.5*np.exp(-r*r/(2*(c.width**2+1j*c.end)))
        result['analytic_free_field_error_inside20']=float(np.max(abs(u[r<=20]/r[r<=20]-analytic[r<=20])))
    if retain:return result,dict(radius=r,initial_amplitude=abs(initial_u/r),final_amplitude=amplitude,
        mean_amplitude=np.sqrt(mean_density),mass_selected_reference=target,snapshots=snapshots,complex_radial_field=u)
    return result


@dataclass
class CoreAtlasEmbedding:
    frequency: float = 2000.
    weak_amplitude: float = .1
    with_core: bool = True
    phase_sign: int = 1

    @cached_property
    def reference(self):return SharedScalarSlice(scalar=ScalarFlux(omega=self.frequency))

    @cached_property
    def profile(self):return solve_core()

    @property
    def cutoff_radius(self):return 20/(self.frequency*self.weak_amplitude)

    def field(self,R):
        if self.phase_sign not in (-1,1):raise ValueError('the admitted coherent turning-point phases are +1 and -1')
        rho=np.asarray(R)*self.frequency*self.weak_amplitude
        S,Sd,_,_=self.profile.fields(rho)
        z=np.clip((rho-15)/5,0,1)
        window=1-(6*z**5-15*z**4+10*z**3)
        window_prime=-6*z*z*(1-z)**2
        amplitude=math.sqrt(2/self.reference.scalar.kappa)*self.weak_amplitude**2*int(self.with_core)*self.phase_sign
        return (amplitude*S*window,
            amplitude*(Sd*window+S*window_prime)*self.frequency*self.weak_amplitude)

    def wall(self,mass):
        s=self.reference;R=s.outer_radius;k=s.k
        ref=ChargedWall(s.reference_law,s.parent_flux+s.scalar.mu*s.field,
            s.regions[1],KottlerRegion(mass,s.regions[2].Lambda))
        beta=(mass-s.middle_mass)/(k*R*R)+(ref.outside.Lambda-ref.inside.Lambda+3*k*k)*R/(6*k)
        v2=beta*beta-ref.inside.f(R)
        if v2<=0 or beta-k*R<=0:raise ValueError('core left ordinary expanding wall branch')
        row=ref.geometry(R,math.sqrt(v2));row['expansion']=math.sqrt(v2)/R
        return row

    def slicing(self,R,H):
        # Retain a homogeneous parent slice expansion. In the zero-core
        # control this recovers the original parent data, including K.
        return H*R,H+np.zeros_like(R)

    def mass_solution(self,H):
        s=self.reference;k=s.scalar.kappa;U0=s.densities[2];w=self.frequency
        slope=s.scalar.omega**2*s.field+s.scalar.mu*s.parent_flux
        def rhs(R,y):
            if R==0:return (0.,)
            f,fp=self.field(R);v,_=self.slicing(R,H)
            beta2=1+v*v-k*U0*R*R/3-2*y[0]/R
            if beta2<=0:raise RuntimeError('core radial chart lost positivity')
            return (k*R*R*(beta2*fp*fp/4+slope*f/2+w*w*f*f/4),)
        result=solve_ivp(rhs,(0,self.cutoff_radius),(0.,),method='DOP853',rtol=2e-11,atol=2e-14,
            dense_output=True,max_step=self.cutoff_radius/120)
        if not result.success:raise RuntimeError(result.message)
        return result,rhs

    def solve(self):
        s=self.reference
        if not self.cutoff_radius<s.outer_radius:raise ValueError('core support reaches the material wall')
        estimate=self.weak_amplitude*self.profile.mass_norm/(2*self.frequency)
        def residual(M):
            sol,_=self.mass_solution(self.wall(M)['expansion'])
            return sol.y[0,-1]-M
        M=brentq(residual,0,3*estimate,xtol=2e-14)
        wall=self.wall(M);H=wall['expansion'];sol,rhs=self.mass_solution(H)
        c=s.densities[2]*s.scalar.kappa/3-H*H
        B=lambda R:1-2*M/R-c*R*R
        Rmax=brentq(B,s.outer_radius,2/math.sqrt(c),xtol=1e-12)
        Rp=self.cutoff_radius
        def bound_integral(bottom,power):
            limit=math.sqrt(Rmax-bottom);bp=2*M/Rmax**2-2*c*Rmax
            def integrand(u):
                if u<1e-7:return 2*Rmax**power/math.sqrt(-bp)
                R=Rmax-u*u
                return 2*u*R**power/math.sqrt(max(B(R),1e-30))
            return quad(integrand,0,limit,epsabs=2e-10)[0]
        def core_values(R):
            if R==0:return 1.,0.,0.
            mass=float(sol.sol(R)[0]);mp=float(rhs(R,[mass])[0]);v,vp=self.slicing(R,H)
            beta2=1+v*v-s.scalar.kappa*s.densities[2]*R*R/3-2*mass/R
            bprime=2*v*vp-2*s.scalar.kappa*s.densities[2]*R/3-2*mp/R+2*mass/R**2
            R3=-2*bprime/R+2*(1-beta2)/R**2
            f,fp=self.field(R)
            density=.5*beta2*fp*fp+s.scalar.potential(s.field+f,s.parent_flux)
            constraint=R3+4*v*vp/R+2*v*v/R**2-2*s.scalar.kappa*density
            return beta2,constraint,density
        core_length=quad(lambda R:1/math.sqrt(core_values(R)[0]),0,Rp,epsabs=1e-11)[0]
        proper_excess=quad(lambda R:4*math.pi*R*R*(core_values(R)[2]-s.densities[2])
            /math.sqrt(core_values(R)[0]),0,Rp,epsabs=2e-12)[0]
        force=s.scalar.omega**2*s.field+s.scalar.mu*s.parent_flux
        pieces=dict(
            gradient=quad(lambda R:s.scalar.kappa*R*R*core_values(R)[0]*self.field(R)[1]**2/4,
                0,Rp,epsabs=2e-13)[0],
            inherited_cross_term=quad(lambda R:s.scalar.kappa*R*R*force*self.field(R)[0]/2,
                0,Rp,epsabs=2e-13)[0],
            quadratic=quad(lambda R:s.scalar.kappa*R*R*self.frequency**2*self.field(R)[0]**2/4,
                0,Rp,epsabs=2e-13)[0])
        length=core_length+bound_integral(Rp,0)+bound_integral(s.outer_radius,0)
        samples=np.linspace(Rp/1000,Rp,401)
        residuals=[abs(core_values(float(R))[1]) for R in samples]
        # Independent proper-radius system, including its geometric mass.
        def ell_rhs(ell,y):
            R,mass=y
            if R<=0:return (1.,0.)
            v,_=self.slicing(R,H)
            beta=math.sqrt(1+v*v-s.scalar.kappa*s.densities[2]*R*R/3-2*mass/R)
            return (beta,float(rhs(R,[mass])[0])*beta)
        check=solve_ivp(ell_rhs,(0,core_length),(0.,0.),method='DOP853',rtol=2e-11,atol=2e-13)
        if not check.success:raise RuntimeError(check.message)
        def middle(R):
            z=(np.asarray(R)-s.inner_radius)/(s.outer_radius-s.inner_radius)
            H1=s.expansions[0]
            step=6*z**5-15*z**4+10*z**3
            step_r=30*z*z*(1-z)**2/(s.outer_radius-s.inner_radius)
            expansion=H1+(H-H1)*step
            v=R*expansion;vp=expansion+R*(H-H1)*step_r
            B=s.regions[1].f(R)+v*v
            Br=s.regions[1].f_prime(R)+2*v*vp
            R3=-2*Br/R+2*(1-B)/R**2
            constraint=R3+4*v*vp/R+2*v*v/R**2-2*s.scalar.kappa*s.densities[1]
            return np.sqrt(B),v,vp,constraint
        left=middle(s.inner_radius);right=middle(s.outer_radius)
        middle_length=quad(lambda R:1/float(middle(R)[0]),s.inner_radius,s.outer_radius,epsabs=1e-12)[0]
        middle_res=max(abs(float(middle(R)[3])) for R in np.linspace(s.inner_radius,s.outer_radius,101))
        seams=dict(inner_beta=float(left[0]-s.walls[0]['outside']['beta']),
            outer_beta=float(right[0]-wall['beta_inside']),
            inner_expansion=float(left[2]-s.expansions[0]),outer_expansion=float(right[2]-H))
        outer_A=[wall['K_inside']-2*wall['beta_inside']/s.outer_radius,
                 wall['K_outside']-2*wall['beta_outside']/s.outer_radius]
        bound=.1*4096/19683*max(abs(x) for x in outer_A+
            [s.walls[0]['inside']['A'],s.walls[0]['outside']['A']])
        return dict(frequency=self.frequency,weak_amplitude=self.weak_amplitude,coherent_phase_sign=self.phase_sign,
            core_cutoff_radius=Rp,geometric_core_mass=M,leading_isolated_mass=estimate,
            mass_difference_from_isolated_estimate=M-estimate,
            geometric_mass_in_energy_units=8*math.pi*M/s.scalar.kappa,
            proper_slice_energy_excess=proper_excess,geometric_mass_parts=pieces,
            independent_field_mass_gap=sum(pieces.values())-M,
            outer_wall=wall,outer_wall_speed_change=wall['velocity']-s.walls[1]['velocity'],
            parent_curvature_parameter=c,parent_maximum_radius=Rmax,parent_pole_to_wall_length=length,
            core_proper_radius=core_length,minimum_core_beta_squared=min(core_values(float(R))[0] for R in samples),
            maximum_sampled_Hamilton_residual=max(residuals),
            independent_radius_gap=float(check.y[0,-1]-Rp),independent_mass_gap=float(check.y[1,-1]-M),
            middle_constraint_residual=middle_res,middle_seams=seams,middle_proper_length=middle_length,
            total_pole_to_pole_length=s.cap_geometry['inner_length']+middle_length+length,
            initial_lapse=dict(collar_width=.1,outer_one_sided_log_derivatives=outer_A,
                global_bounds=[math.exp(-bound),math.exp(bound)],
                construction='Same positive compact log-lapse collars as the shared-slice owner, with the updated outer A jets; unit lapse at the core and both poles.'),
            fixed_middle_mass=s.middle_mass,inner_wall_unchanged=True,
            scalar_value_and_first_derivative_at_core_cutoff=list(map(float,self.field(Rp))),
            scope='Exact spherical initial-constraint equations with a supplied coherent core profile from the envelope eigenproblem. This is not its full nonlinear formation history or a retroactive preparation of the P003 packet.')
