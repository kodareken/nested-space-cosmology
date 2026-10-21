"""Finite parent scalar Cauchy data with gravity and material-wall matching.

This constructs time-symmetric initial data, not their finite-time evolution.
The smooth branch solves a scalar second-normal jet with the Hamiltonian
constraint. A distinct uniform seed supplies a conditional weak transfer limit;
second-order smooth corner compatibility is not imposed on that second state.
"""
from dataclasses import dataclass, field
import math

import numpy as np
from scipy.integrate import quad, solve_bvp
from scipy.special import j0

from .nsc_atlas_scalar_preparation import ScalarPreparation


@dataclass(frozen=True)
class ParentScalarSource:
    prepared: ScalarPreparation = field(default_factory=ScalarPreparation)
    amplitude: float = .02
    center: float = .3
    half_width: float = .12
    collar_width: float = .08

    def __post_init__(self):
        if (not all(map(math.isfinite,(self.amplitude,self.center,self.half_width,self.collar_width)))
                or not 0 <= self.amplitude <= .02
                or not 0 < self.collar_width < self.center-self.half_width
                or self.half_width <= 0):
            raise ValueError("requires 0<=amplitude<=.02 and disjoint positive packet/collar widths")

    def packet(self, r, parent_length):
        """Specified C5 compact displacement and seed-distance derivative."""
        z = (parent_length-np.asarray(r)-self.center)/self.half_width
        b = np.maximum(0.,1-z*z)
        return self.amplitude*b**6, self.amplitude*12*z*b**5/self.half_width

    def collar(self, r, parent_length, second_jet):
        d = parent_length-np.asarray(r)
        z = d/self.collar_width
        b = np.maximum(0.,1-z)
        on = (z>=0)&(z<1)
        return (np.where(on,second_jet*d*d*b**6/2,0.),
                np.where(on,-second_jet*(d*b**6-3*d*d*b**5/self.collar_width),0.))

    def required_second_jet(self, wall_factor, seed_normal_field):
        m = self.prepared
        return (m.scalar.mu*m.charge*(wall_factor**4-1)
                +2*m.k*seed_normal_field*(wall_factor**2-1))

    def solve(self, background, *, tolerance=1e-8, match_scalar_jet=True):
        if not 1e-10 <= tolerance <= 1e-6:
            raise ValueError("declared source tolerance range is 1e-10..1e-6")
        m = self.prepared
        Li,Lo = map(float,background.p)
        if self.center+self.half_width >= Lo:
            raise ValueError("parent source must avoid the pole")
        kappa,w = m.scalar.kappa,m.scalar.omega
        calls = 0

        def rhs(x,y,params):
            nonlocal calls
            calls += 1
            if calls > 3000:
                raise RuntimeError("source collocation RHS budget exceeded")
            old = background.sol(x)
            out = np.empty_like(y)
            for side,L,rv in ((0,Li,m.child_vacuum),(1,Lo,m.parent_vacuum)):
                a,ap,z,zp = old[4*side:4*side+4]
                if side:
                    f,fp = self.packet(L*x,Lo)
                    g,gp = self.collar(L*x,Lo,params[0])
                    f,fp = f+g,fp+gp
                else:
                    f,fp = 0.,0.
                psi,p = y[2*side:2*side+2]
                U0,U = rv+w*w*z*z/2,rv+w*w*(z+f)**2/2
                F = kappa*(zp*zp-(zp+fp)**2)*psi/8+kappa*(U0*psi-U*psi**5)/4
                regular = np.zeros_like(x)
                mask = x>0
                regular[mask] = L*ap[mask]/a[mask]-1/x[mask]
                out[2*side] = L*p
                out[2*side+1] = L*F-2*regular*p
            return out

        def boundary(left,right,params):
            W = right[0]
            required = self.required_second_jet(W,background.y[3,-1]) if match_scalar_jet else 0.
            return np.array([left[1],left[3],W-right[2],
                             right[1]+right[3]-m.k/2*(W**3-W),params[0]-required])

        grid = np.unique(np.r_[np.linspace(0,1,401),
            np.linspace(1-(self.center+self.half_width)/Lo,1-(self.center-self.half_width)/Lo,121),
            np.linspace(1-self.collar_width/Lo,1,61)])
        initial = np.zeros((4,len(grid)))
        initial[[0,2]] = 1
        singular = np.zeros((4,4))
        singular[1,1] = singular[3,3] = -2
        solution = solve_bvp(rhs,boundary,grid,initial,p=[0.],S=singular,
                             tol=tolerance,bc_tol=min(tolerance,1e-11),max_nodes=10000)
        if not solution.success or np.min(solution.y[[0,2]]) <= .5:
            raise RuntimeError("positive near-reference source branch failed: "+solution.message)
        return solution,dict(nodes=len(solution.x),rhs_calls=calls,iterations=solution.niter,
            maximum_collocation_residual=float(np.max(solution.rms_residuals)),
            maximum_boundary_residual=float(np.max(np.abs(boundary(solution.y[:,0],solution.y[:,-1],solution.p)))))

    def wall(self, background, source):
        m = self.prepared
        v = background.y[:,-1]
        w = source.y[:,-1]
        R = w[0]**2*v[0]
        beta_i = v[1]+2*v[0]*w[1]/w[0]
        beta_o = -(v[5]+2*v[4]*w[3]/w[2])
        phi = m.scalar.minimum(m.child_flux)+v[2]
        pressure = m.scalar.wall_pressure(phi,m.parent_flux,m.charge)
        mean_A = pressure/m.tension-(beta_i+beta_o)/R
        A_i,A_o = mean_A+m.k/2,mean_A-m.k/2
        rows = []
        for beta,A,j,W,p,rv,sign,jet in (
            (beta_i,A_i,0,w[0],w[1],m.child_vacuum,1,0.),
            (beta_o,A_o,4,w[2],w[3],m.parent_vacuum,-1,float(source.p[0]))):
            normal = sign*v[j+3]/W**2
            second = (m.scalar.omega**2*v[j+2]-3*v[j+1]*v[j+3]/v[j]
                      +jet-2*p*v[j+3]/W)/W**4
            U = rv+m.scalar.omega**2*v[j+2]**2/2
            radial_pressure = normal*normal/2-U
            mass = R*(1-beta*beta)/2
            acceleration = -mass/R**2-m.scalar.kappa*R*radial_pressure/2+A*beta
            field_acceleration = second+2*beta*normal/R-m.scalar.omega**2*v[j+2]+A*normal
            rows.append(dict(beta=float(beta),angular_curvature=float(beta/R),
                normal_acceleration=float(A),normal_field=float(normal),second_normal_field=float(second),
                radial_pressure=float(radial_pressure),Misner_Sharp_mass=float(mass),
                wall_areal_acceleration=float(acceleration),wall_field_acceleration=float(field_acceleration),
                outward_flux_slope=float(-normal*field_acceleration)))
        inside,outside = rows
        scalar_compatibility = (outside["second_normal_field"]-inside["second_normal_field"]
            -3*m.k*inside["normal_field"]-m.scalar.mu*m.charge)
        return dict(radius=float(R),wall_factor=float(w[0]),second_jet=float(source.p[0]),
            pressure=float(pressure),inside=inside,outside=outside,
            angular_Israel_gap=float((beta_o-beta_i)/R+m.k),
            temporal_Israel_gap=float(A_o-A_i+m.k),
            mean_force_gap=float(m.tension*(mean_A+(beta_i+beta_o)/R)-pressure),
            scalar_second_jet_gap=float(scalar_compatibility),
            areal_acceleration_gap=outside["wall_areal_acceleration"]-inside["wall_areal_acceleration"],
            field_acceleration_gap=outside["wall_field_acceleration"]-inside["wall_field_acceleration"],
            outward_flux_slope_gap=outside["outward_flux_slope"]-inside["outward_flux_slope"])

    def diagnostics(self, background, source):
        """Recover the spatial scalar curvature from the physical areal metric."""
        m = self.prepared
        x = np.unique(np.r_[np.linspace(1e-3,1,1001),
            np.linspace(1-(self.center+self.half_width)/background.p[1],1,501)])
        old = background.sol(x)
        y,dy = source.sol(x),source.sol(x,1)
        weighted,unweighted = [],[]
        for side,L,rv in ((0,background.p[0],m.child_vacuum),(1,background.p[1],m.parent_vacuum)):
            a,ap,z,zp = old[4*side:4*side+4]
            psi,p = y[2*side:2*side+2]
            pp = dy[2*side+1]/L
            f,fp = (0.,0.)
            if side:
                f,fp = self.packet(L*x,L)
                g,gp = self.collar(L*x,L,source.p[0])
                f,fp = f+g,fp+gp
            U0,U = rv+m.scalar.omega**2*z*z/2,rv+m.scalar.omega**2*(z+f)**2/2
            app = -m.scalar.kappa*a*(zp*zp+U0)/3
            R = psi*psi*a
            beta = ap+2*a*p/psi
            Rll = (app+2*ap*p/psi+2*a*pp/psi-2*a*p*p/psi**2)/psi**2
            rho = (zp+fp)**2/(2*psi**4)+U
            residual_times_R2 = -4*Rll*R+2*(1-beta*beta)-2*m.scalar.kappa*rho*R*R
            weighted.append(float(np.max(np.abs(residual_times_R2))))
            unweighted.append(float(np.max(np.abs(residual_times_R2/(R*R)))))
        return dict(weighted_Hamilton_residual=weighted,Hamilton_residual=unweighted,
                    minimum_factor=float(np.min(y[[0,2]])),maximum_factor=float(np.max(y[[0,2]])),
                    momentum_constraint="Exactly zero: K_ij=0 and scalar Pi=0.")

    def energy(self, background, source):
        """Proper scalar excess on the sourced metric, including interference.

        The baseline field in this subtraction uses the SAME sourced metric.
        It is not an ADM energy or a difference between entire universe energies.
        """
        m = self.prepared
        L = float(background.p[1])
        def integrand(r):
            psi = source.sol(r/L)[2]
            a,_,z,zp = background.sol(r/L)[4:8]
            f,fp = self.packet(r,L)
            g,gp = self.collar(r,L,source.p[0])
            delta,dp = f+g,fp+gp
            weight = 4*math.pi*a*a*psi**6
            return np.array([weight*(dp*dp/psi**4+m.scalar.omega**2*delta*delta)/2,
                             weight*(zp*dp/psi**4+m.scalar.omega**2*z*delta)])
        intervals = ((L-self.center-self.half_width,L-self.center+self.half_width),
                     (L-self.collar_width,L))
        quad_values = [sum(quad(lambda r: integrand(r)[i],a,b,epsabs=2e-10,epsrel=2e-10)[0]
                           for a,b in intervals) for i in (0,1)]
        return dict(positive_quadratic=float(quad_values[0]),interference=float(quad_values[1]),
                    total_scalar_excess=float(sum(quad_values)))

    def mass_checks(self, background, source):
        """Independent integrated Hamiltonian constraint from each regular pole."""
        m = self.prepared
        wall = self.wall(background,source)
        checks = []
        for side,L,rv,name in ((0,float(background.p[0]),m.child_vacuum,"inside"),
                               (1,float(background.p[1]),m.parent_vacuum,"outside")):
            def integrand(r):
                a,ap,z,zp = background.sol(r/L)[4*side:4*side+4]
                psi,p = source.sol(r/L)[2*side:2*side+2]
                f,fp = (0.,0.)
                if side:
                    f,fp = self.packet(r,L)
                    g,gp = self.collar(r,L,source.p[0])
                    f,fp = f+g,fp+gp
                R = a*psi*psi
                beta = ap+2*a*p/psi
                rho = (zp+fp)**2/(2*psi**4)+rv+m.scalar.omega**2*(z+f)**2/2
                return m.scalar.kappa*R*R*rho*beta*psi*psi/2
            cuts = [0.,L]
            if side:
                cuts.extend([L-self.center-self.half_width,L-self.center+self.half_width,L-self.collar_width])
            cuts = sorted(cuts)
            value = sum(quad(integrand,a,b,epsabs=2e-10,epsrel=2e-10,limit=200)[0]
                        for a,b in zip(cuts[:-1],cuts[1:]))
            checks.append(dict(side=name,integrated_mass=float(value),
                               geometric_mass=wall[name]["Misner_Sharp_mass"],
                               gap=float(value-wall[name]["Misner_Sharp_mass"])))
        return checks


def uniform_excited_parent(prepared=None):
    """Exact round Cauchy caps with a uniform field at the CHILD minimum.

    The parent is displaced from its minimum. Its subsequent spacetime is
    not de Sitter. The initial scalar corner fails C2 compatibility; the
    incoming limit applies to bounded-geometry causal weak developments.
    """
    m = ScalarPreparation() if prepared is None else prepared
    phi = m.scalar.minimum(m.child_flux)
    rhoi = m.child_vacuum
    rhop = m.scalar.potential(phi,m.parent_flux)
    hi,hp = math.sqrt(m.scalar.kappa*rhoi/3),math.sqrt(m.scalar.kappa*rhop/3)
    Bi = (hp*hp-hi*hi+m.k*m.k)/(2*m.k)
    Bo = Bi-m.k
    if min(hi,hp,Bo)<=0:
        raise ValueError('uniform preparation requires the ordinary positive two-cap branch')
    R = 1/math.hypot(hi,Bi)
    Li = math.asin(hi*R)/hi
    theta_o = math.pi-math.asin(hp*R)
    Lo = theta_o/hp
    volume = 2*math.pi/hp**3*(theta_o-math.sin(2*theta_o)/2)
    S = m.scalar.mu*m.charge
    pressure = m.scalar.wall_pressure(phi,m.parent_flux,m.charge)
    coeff = S*S/4
    return dict(field_value=phi,child_density=rhoi,parent_initial_density=rhop,
        parent_ground_density=m.parent_vacuum,parent_excitation_density=rhop-m.parent_vacuum,
        parent_proper_excitation=(rhop-m.parent_vacuum)*volume,
        radius=R,child_length=Li,parent_length=Lo,child_h=hi,parent_initial_h=hp,
        beta_inside=Bi*R,beta_outside=Bo*R,initial_wall_areal_acceleration=1/R,
        pressure=pressure,field_force_jump=S,
        angular_Israel_gap=(Bi-Bo-m.k)*R,
        mean_force_gap=m.tension*3*(Bi+Bo)/2-pressure,
        child_constraint_gap=(Bi*R)**2+hi*hi*R*R-1,
        parent_constraint_gap=(Bo*R)**2+hp*hp*R*R-1,
        scalar_initial_C2_gap=-S,
        conditional_inward_flux_over_time_squared=coeff,
        conditional_child_energy_per_area_over_time_cubed=S*S/12,
        conditional_child_energy_over_time_cubed=4*math.pi*R*R*S*S/12,
        conditional_pressure_drop_over_time_squared=coeff,
        earliest_central_front_time=math.asinh(hi/Bi)/hi,
        scope='Exact initial constraints. Flux/energy coefficients are conditional short-time limits of a causal weak development with bounded one-sided geometry; no finite-time coupled evolution or error bar is supplied.')


def planar_interface_control(time, *, force_jump, omega):
    """Exact massive 1+1 retarded step-source control, not the curved spacetime.

    Child x<0 starts at its minimum, and parent x>0 has force -S.
    Reflection plus homogeneous forcing fixes the wall value. The spatial
    derivative follows independently from the retarded Bessel kernel.
    """
    if not all(map(math.isfinite,(time,force_jump,omega))) or min(time,omega)<=0:
        raise ValueError('positive finite time and frequency required')
    def inward(t):
        integral = quad(lambda u:j0(omega*u),0,t,epsabs=1e-14)[0]
        return force_jump**2*math.sin(omega*t)*integral/(4*omega)
    integral = quad(lambda u:j0(omega*u),0,time,epsabs=1e-14)[0]
    field = -force_jump*math.sin(omega*time/2)**2/omega**2
    dt = -force_jump*math.sin(omega*time)/(2*omega)
    dx = -force_jump*integral/2
    return dict(time=time,wall_displacement=field,wall_time_derivative=dt,
                wall_normal_derivative=dx,inward_flux=dt*dx,
                child_energy_per_area=quad(inward,0,time,epsabs=1e-14)[0],
                pressure_drop=-force_jump*field,
                scope='Exact planar fixed-interface massive-field control only.')


def planar_volume_energies(time, *, force_jump, omega, order=24):
    """Independent retarded-volume integral of scalar energy inside x<0.

    The two nested Bessel integrals evaluate the actual field, and separate
    single integrals its derivatives. This does not integrate boundary flux.
    """
    if not 0<omega*time<=1 or not 8<=order<=64:
        raise ValueError('control requires 0<omega*time<=1 and 8<=order<=64')
    nodes,weights=np.polynomial.legendre.leggauss(order)
    nodes,weights=(nodes+1)/2,weights/2
    values=[]
    parent_values=[]
    homogeneous=-2*force_jump*math.sin(omega*time/2)**2/omega**2
    homogeneous_dt=-force_jump*math.sin(omega*time)/omega
    for r in time*nodes:
        u=r+(time-r)*nodes
        z=r+(u[:,None]-r)*nodes[None,:]
        kernel=j0(omega*np.sqrt(np.maximum(0.,u[:,None]**2-z*z)))
        phi=-force_jump/2*(time-r)*np.dot(weights,(u-r)*(kernel@weights))
        dt=-force_jump/2*(time-r)*np.dot(weights,j0(omega*np.sqrt(np.maximum(0.,time*time-u*u))))
        dx=-force_jump/2*(time-r)*np.dot(weights,j0(omega*np.sqrt(np.maximum(0.,u*u-r*r))))
        values.append((dt*dt+dx*dx+omega*omega*phi*phi)/2)
        parent_phi=homogeneous-phi
        parent_dt=homogeneous_dt-dt
        parent_values.append((parent_dt**2+dx*dx+omega*omega*parent_phi**2)/2+force_jump*parent_phi)
    child=float(time*np.dot(weights,values))
    parent=float(time*np.dot(weights,parent_values))
    return dict(child=child,parent_excess=parent,net=child+parent)


def planar_child_energy_integral(time, *, force_jump, omega, order=24):
    """Child component of the independent two-sided retarded volume integral."""
    return planar_volume_energies(time,force_jump=force_jump,omega=omega,order=order)['child']
