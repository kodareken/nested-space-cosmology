"""A charged scalar layer solved with gravity and its Lorentzian preparation.

Regular O(4) boundary data define a candidate saddle; no fluctuation index,
quantum probability or observational abundance is inferred from convergence.
The two caps use separate pole-to-wall coordinates with opposite normals.
"""
from dataclasses import dataclass, field
import math

import numpy as np
from scipy.integrate import solve_bvp, solve_ivp
from scipy.integrate import quad

from .nsc_atlas_flux_conversion import ScalarFlux


@dataclass(frozen=True)
class ScalarPreparation:
    scalar: ScalarFlux = field(default_factory=ScalarFlux)
    parent_flux: float = 1.5
    charge: float = .5
    tension: float = .2

    def __post_init__(self):
        if (not all(math.isfinite(x) for x in (self.parent_flux,self.charge,self.tension))
                or min(self.charge,self.tension)<=0):
            raise ValueError("positive finite charge and tension required")
        if not self.parent_vacuum>self.child_vacuum>0:
            raise ValueError("this compact branch requires two positive descending vacua")

    @property
    def child_flux(self):
        return self.parent_flux-self.charge

    @property
    def child_vacuum(self):
        return self.scalar.vacuum(self.child_flux)

    @property
    def parent_vacuum(self):
        return self.scalar.vacuum(self.parent_flux)

    @property
    def shift(self):
        return self.scalar.minimum(self.child_flux)-self.scalar.minimum(self.parent_flux)

    @property
    def k(self):
        return self.scalar.kappa*self.tension/2

    def vacuum_guess(self):
        """Pure-vacuum caps initialize the solve; the answer retains scalar stress."""
        ki,ko=self.scalar.kappa*self.child_vacuum,self.scalar.kappa*self.parent_vacuum
        B=(ko-ki)/(6*self.k)+self.k/2
        if B-self.k<=0:
            raise ValueError("ordinary outer-normal branch required")
        R=1/math.sqrt(ki/3+B*B)
        hi,ho=math.sqrt(ki/3),math.sqrt(ko/3)
        return R,math.asin(hi*R)/hi,(math.pi-math.asin(ho*R))/ho

    def boundary(self,left,right,lengths):
        return np.array([left[0],left[1]-1,left[3],left[4],left[5]-1,left[7],
                         right[0]-right[4],right[2]-right[6]+self.shift,
                         right[3]+right[7],right[1]+right[5]-self.k*right[0]])

    def solve(self,*,tolerance=1e-8):
        if not 1e-10<=tolerance<=1e-5:
            raise ValueError("declared collocation tolerance range is 1e-10..1e-5")
        _,Li,Lo=self.vacuum_guess()
        w,kappa=self.scalar.omega,self.scalar.kappa
        x=np.linspace(0,1,321)
        y=np.empty((8,len(x)))
        for j,L,rv,sign in ((0,Li,self.child_vacuum,-1),(4,Lo,self.parent_vacuum,1)):
            h=math.sqrt(kappa*rv/3)
            y[j]=np.sin(h*L*x)/h
            y[j+1]=np.cos(h*L*x)
            y[j+2]=sign*self.shift/2*np.exp(-w*L*(1-x))
            y[j+3]=w*y[j+2]

        calls=0
        def rhs(x,y,p):
            nonlocal calls
            calls+=1
            if calls>3000:
                raise RuntimeError("finite boundary-value budget exceeded")
            out=np.empty_like(y)
            for j,L,rv in ((0,p[0],self.child_vacuum),(4,p[1],self.parent_vacuum)):
                a,ap,z,zp=y[j:j+4]
                regular=np.zeros_like(x)
                mask=x>0
                regular[mask]=L*ap[mask]/a[mask]-1/x[mask]
                out[j]=L*ap
                out[j+1]=-kappa*L*a*(zp*zp+rv+w*w*z*z/2)/3
                out[j+2]=L*zp
                out[j+3]=L*w*w*z-3*regular*zp
            return out

        S=np.zeros((8,8))
        S[3,3]=S[7,7]=-3
        solution=solve_bvp(rhs,self.boundary,x,y,p=[Li,Lo],S=S,
                           tol=tolerance,bc_tol=min(tolerance,1e-11),max_nodes=12000)
        if not solution.success or min(solution.p)<=0:
            raise RuntimeError(solution.message)
        wall=solution.y[:,-1]
        if not (wall[1]>0>wall[5] and wall[0]>0):
            raise RuntimeError("solution left the declared ordinary-cap branch")
        return solution,dict(nodes=len(solution.x),rhs_calls=calls,iterations=solution.niter,
                             max_collocation_residual=float(np.max(solution.rms_residuals)))

    def diagnostics(self,solution):
        # Independent Einstein constraint, not imposed as an evolution equation.
        x=np.linspace(1e-4,1,1025)
        y=solution.sol(x)
        dy=solution.sol(x,1)
        constraints=[]
        scalar_equation=[]
        tangent_equation=[]
        for j,L,rv in ((0,solution.p[0],self.child_vacuum),(4,solution.p[1],self.parent_vacuum)):
            a,ap,z,zp=y[j:j+4]
            U=rv+self.scalar.omega**2*z*z/2
            C=ap*ap-1-self.scalar.kappa*a*a/3*(zp*zp/2-U)
            constraints.append(float(np.max(np.abs(C))))
            scalar_equation.append(float(np.max(np.abs(dy[j+3]/L+3*ap/a*zp-self.scalar.omega**2*z))))
            # Multiply by a^2 to avoid a removable coordinate cancellation at a=0.
            E=2*a*dy[j+1]/L+ap*ap-1+self.scalar.kappa*a*a*(zp*zp/2+U)
            tangent_equation.append(float(np.max(np.abs(E))))
        v=solution.y[:,-1]
        R=v[0]
        phi=self.scalar.minimum(self.child_flux)+v[2]
        Ki,Ko=3*v[1]/R,-3*v[5]/R
        pressure=self.scalar.wall_pressure(phi,self.parent_flux,self.charge)
        return dict(radius=float(R),lengths=[float(t) for t in solution.p],
                    phi_wall=float(phi),child_wall_displacement=float(v[2]),
                    parent_wall_displacement=float(v[6]),child_normal_derivative=float(v[3]),
                    parent_normal_derivative=float(-v[7]),K_inside=float(Ki),K_outside=float(Ko),
                    wall_pressure=float(pressure),
                    force_residual=float(self.tension*(Ki+Ko)/2-pressure),
                    full_Israel_component_residual=float((Ko-Ki)/3+self.k),
                    scalar_value_jump=float(v[2]-v[6]+self.shift),
                    scalar_normal_jump=float(v[3]+v[7]),
                    Hamilton_max=constraints,scalar_equation_max=scalar_equation,
                    weighted_tangential_Einstein_max=tangent_equation,
                    max_boundary_residual=float(np.max(np.abs(self.boundary(solution.y[:,0],v,solution.p)))))

    def pole_transfer(self,solution,*,side=0,rtol=2e-11):
        """Resolve the regular mode's logarithmic gain independently of tiny tails.

        chi(0)=chi(wall)*exp(-loggain); ratio w=chi'/chi satisfies a Riccati ODE.
        The decaying sign of chi is in the supplied nonzero wall value.
        """
        if side not in (0,1):
            raise ValueError("side must be child0 or parent1")
        j,L=4*side,float(solution.p[side])
        amplitude=float(solution.y[j+2,-1])
        if amplitude==0:
            return dict(pole_displacement=0.,wall_displacement=0.,log_gain=0.,
                        logarithmic_derivative_gap=0.,scalar_energy_at_pole=0.)
        w=self.scalar.omega
        h=1e-6/w
        def rhs(rho,y):
            a,ap=solution.sol(rho/L)[j:j+2]
            return w*w-3*ap/a*y[0]-y[0]**2,y[0]
        result=solve_ivp(rhs,(h,L),(w*w*h/4,w*w*h*h/8),method="DOP853",
                         rtol=rtol,atol=rtol/50,max_step=.1/w)
        if not result.success:
            raise RuntimeError(result.message)
        slope,gain=result.y[:,-1]
        pole=amplitude*math.exp(-gain)
        return dict(pole_displacement=pole,wall_displacement=amplitude,log_gain=float(gain),
                    logarithmic_derivative_gap=float(slope-solution.y[j+3,-1]/amplitude),
                    scalar_energy_at_pole=.5*w*w*pole*pole,
                    ode_evaluations=result.nfev)

    def lorentzian_wall(self,radius,eta):
        if not (radius>0 and math.isfinite(eta) and abs(eta)<=1):
            raise ValueError("declared wall continuation is |eta|<=1 at positive radius")
        return dict(eta=eta,proper_time=radius*eta,areal_radius=radius*math.cosh(eta),
                    radius_proper_speed=math.sinh(eta),
                    radius_proper_acceleration=math.cosh(eta)/radius,
                    worldvolume_scalar_curvature=6/radius**2,
                    scalar_rest_frame_flux=0.)

    def scalar_reduction(self,solution):
        """Actual regular-cap Dirichlet-to-Neumann maps, at the solved geometry.

        This is an elliptic saddle-sector Schur reduction, not a Lorentzian
        retarded kernel and not an additional bulk cosmological constant.
        """
        v=solution.y[:,-1]
        if self.shift==0:
            return dict(zero_mixing=True,field_layer_energy_per_wall_area=0.)
        si,so=v[3]/v[2],v[7]/v[6]
        predicted_i=-self.shift*so/(si+so)
        predicted_o=self.shift*si/(si+so)
        R=v[0]
        action=math.pi**2*R**3*self.shift**2*si*so/(si+so)
        integrated_action=0.
        energy=0.
        for j,L in ((0,solution.p[0]),(4,solution.p[1])):
            def density(x,power):
                a,ap,z,zp=solution.sol(x)[j:j+4]
                return a**power*(zp*zp+self.scalar.omega**2*z*z)/2*L
            integrated_action+=2*math.pi**2*quad(lambda x:density(x,3),0,1,epsabs=2e-12)[0]
            energy+=quad(lambda x:density(x,2),0,1,epsabs=2e-12)[0]/R**2
        return dict(child_DtN=float(si),parent_DtN=float(so),
                    Schur_coefficient=float(si*so/(si+so)),
                    wall_displacement_reconstruction_gap=float(max(abs(predicted_i-v[2]),abs(predicted_o-v[6]))),
                    boundary_excitation_action=float(action),
                    direct_bulk_excitation_action=float(integrated_action),
                    bulk_boundary_action_gap=float(integrated_action-action),
                    field_layer_energy_per_wall_area=float(energy))

    def child_response(self,pole_amplitude,*,end=1.):
        """Open de Sitter response with a separate rigorous energy comparison.

        The normalized scalar is computed on the vacuum background. Actual
        homogeneous backreaction is bounded by its nonincreasing initial
        excitation; the small numerical amplitude is not rounded to zero.
        """
        if not 0<end<=1:
            raise ValueError("declared child proper interval is 0<tau<=1")
        w,kappa=self.scalar.omega,self.scalar.kappa
        H=math.sqrt(kappa*self.child_vacuum/3)
        t0=1e-6/w
        def rhs(t,y):
            return y[1],-3*H/math.tanh(H*t)*y[1]-w*w*y[0]
        solution=solve_ivp(rhs,(t0,end),(1-w*w*t0*t0/8,-w*w*t0/4),method="DOP853",
                           rtol=2e-12,atol=2e-14,max_step=.05/w,dense_output=True)
        if not solution.success:
            raise RuntimeError(solution.message)
        u,v=solution.y[:,-1]
        E0=.5*w*w*pole_amplitude**2
        E=.5*pole_amplitude**2*(v*v+w*w*u*u)
        return dict(proper_time=end,normalized_field=float(u),normalized_field_velocity=float(v),
                    excitation_density=float(E),initial_excitation_density=E0,
                    tidal_trace_increment=float(kappa*pole_amplitude**2*(w*w*u*u/2-v*v)),
                    vacuum_tidal_trace=kappa*self.child_vacuum,
                    vacuum_scale_factor=math.sinh(H*end)/H,
                    max_extra_H_squared_at_fixed_a=kappa*E0/3,
                    spatial_curvature=-1,
                    approximation="Scalar on open de Sitter, with actual Friedmann source correction bounded by kappa*E0/3 at fixed scale factor. Not an exactly flat FRW quench.")
