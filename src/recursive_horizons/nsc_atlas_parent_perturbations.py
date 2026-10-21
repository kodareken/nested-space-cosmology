"""Unaveraged canonical-scalar/metric transfer on the protected closed parent.

All results are linear in a declared perturbation of the same action. The
Newtonian-gauge system avoids division by the oscillating scalar velocity.
Its S3 harmonics are computational extensions, not membrane eigenmodes.
"""
from dataclasses import dataclass,field
from functools import cached_property
import math

import numpy as np
from scipy.integrate import solve_ivp,quad
from scipy.special import eval_chebyu

from .nsc_atlas_shared_scalar_slice import SharedScalarSlice


@dataclass(frozen=True)
class ParentModes:
    atlas: SharedScalarSlice = field(default_factory=SharedScalarSlice)

    @property
    def initial(self):
        s=self.atlas
        return dict(a=s.cap_geometry['parent_maximum_radius'],H=s.expansions[1],
            chi=s.field-s.scalar.minimum(s.parent_flux),Pi=0.,
            vacuum=s.scalar.vacuum(s.parent_flux),kappa=s.scalar.kappa,omega=s.scalar.omega)

    def solve(self,modes=(2,8,32,128),*,clock='proper',efolds=1.,rtol=2e-10,atol=2e-12):
        n=np.asarray(modes,dtype=int)
        if (n.ndim!=1 or len(n)==0 or np.any(n<2) or np.any(n>256)
                or np.any(n!=np.asarray(modes)) or len(set(n))!=len(n)):
            raise ValueError('distinct integer harmonics 2<=n<=256 are required')
        if clock not in ('proper','efolds') or not 0<efolds<=1:
            raise ValueError('declared domain is one e-fold, in proper or expansion clock')
        p=self.initial;kap,w=p['kappa'],p['omega'];lam=n*(n+2)
        E0=w*w*p['chi']**2/2
        z=np.zeros((3,len(n),2))
        # Column 0: unit initial comoving excitation-density contrast, Y=0.
        # Column 1: zero density but unit canonical phase displacement Y/(w chi).
        # Physical perturbations multiply these linear coefficients by a small
        # declared amplitude; unit transfer columns are not finite solutions.
        z[0,:,0]=p['chi']/2
        z[1,:,1]=w*p['chi']
        z[2,:,0]=-kap*p['a']**2*E0/(2*(lam-3))
        state0=np.r_[p['chi'],0.,0.,p['H'],0.,0.,0.,z.ravel()]
        calls=0
        def rhs(clock_value,y):
            nonlocal calls
            calls+=1
            if calls>40000:raise RuntimeError('parent perturbation RHS budget exceeded')
            chi,Pi,N,H,_,_,_=y[:7]
            a=p['a']*math.exp(N)
            X,Y,Phi=y[7:].reshape(3,len(n),2)
            Phid=-H*Phi+kap*Pi*X/2
            Yd=-3*H*Y-(w*w+lam[:,None]/a**2)*X+4*Pi*Phid-2*w*w*chi*Phi
            rates=self.atlas.parent_rates(chi,Pi,N,H)
            full=np.r_[rates,1.,np.array([Y,Yd,Phid]).ravel()]
            if clock=='efolds':
                if H<=0:raise RuntimeError('expanding branch lost')
                full/=H
            return full
        def endpoint(t,y):return y[2]-efolds
        endpoint.terminal=True;endpoint.direction=1
        end=efolds if clock=='efolds' else 1.05*efolds/self.atlas.parent_causal_domain()['H_positive_lower_bound']
        sol=solve_ivp(rhs,(0,end),state0,method='DOP853',dense_output=True,
            events=endpoint if clock=='proper' else None,rtol=rtol,atol=atol,
            max_step=.02 if clock=='proper' else .006)
        if not sol.success or (clock=='proper' and not len(sol.t_events[0])):
            raise RuntimeError('linear parent endpoint failed: '+sol.message)
        return ModeEvolution(self,n,sol,clock,calls)


@dataclass
class ModeEvolution:
    model: ParentModes
    modes: np.ndarray
    solution: object
    clock: str
    rhs_calls: int

    def sample(self,values):
        values=np.atleast_1d(values)
        if np.any(values<0) or np.any(values>self.solution.t[-1]+1e-12):
            raise ValueError('sample outside solved interval')
        y=self.solution.sol(values)
        chi,Pi,N,H,work,eta,t=y[:7]
        p=self.model.initial;a=p['a']*np.exp(N);kap,w=p['kappa'],p['omega']
        E=(Pi*Pi+w*w*chi*chi)/2
        X,Y,Phi=y[7:].reshape(3,len(self.modes),2,-1)
        lam=(self.modes*(self.modes+2))[:,None,None];D=lam-3
        drho=Pi*Y-Pi*Pi*Phi+w*w*chi*X
        drho_com=drho+3*H*Pi*X
        poisson=-2*D*Phi/(kap*a*a)
        constraint=D*Phi/(a*a)+kap*drho_com/2
        # Normalize by the fixed initial excitation, never an instantaneous
        # scalar velocity or a coefficient which can cross zero.
        E0=w*w*p['chi']**2/2
        return dict(proper_time=t,N=N,a=a,H=H,chi=chi,Pi=Pi,E=E,work=work,conformal=eta,
            X=X,Y=Y,Phi=Phi,Delta=poisson/E,Delta_from_matter=drho_com/E,
            Hamilton_residual=constraint,Hamilton_normalized=constraint/(kap*E0),
            Friedmann=3*(H*H+1/(a*a))-kap*(p['vacuum']+E),
            energy_balance=E+work-E0)

    def summary(self):
        history=self.sample(np.linspace(0,self.solution.t[-1],513))
        end=self.sample([self.solution.t[-1]])
        p=self.model.initial
        rows=[]
        for j,n in enumerate(self.modes):
            initial_Phi=-p['kappa']*p['a']**2*(p['omega']**2*p['chi']**2/2)/(2*(n*(n+2)-3))
            rows.append(dict(n=int(n),eigenvalue=int(n*(n+2)),
                initial_physical_wavenumber=math.sqrt(n*(n+2))/p['a'],
                density_seed_final_Delta=float(end['Delta'][j,0,0]),
                phase_seed_final_Delta=float(end['Delta'][j,1,0]),
                density_seed_potential_ratio=float(end['Phi'][j,0,0]/initial_Phi)))
        return dict(clock=self.clock,proper_time=float(end['proper_time'][0]),efolds=float(end['N'][0]),
            rows=rows,rhs_calls=self.rhs_calls,
            max_normalized_Hamilton=float(np.max(abs(history['Hamilton_normalized']))),
            max_density_constraint_gap=float(np.max(abs(history['Delta']-history['Delta_from_matter']))),
            max_Friedmann=float(np.max(abs(history['Friedmann']))),
            max_energy_balance=float(np.max(abs(history['energy_balance']))),
            scalar_velocity_zero_crossings=int(np.sum(np.diff(np.sign(history['Pi'][1:]))!=0)))


@dataclass(frozen=True)
class ParentPacket:
    """A C7 compensated radial potential supported away from both walls.

    Initial n=0 and n=1 projections vanish by two weighted moments. The
    exact profile has compact support; finite spectral approximations have
    small tails which must be measured, not declared identically zero.
    """
    model: ParentModes = field(default_factory=ParentModes)
    support_angle: float = 1.4
    power: int = 8

    def __post_init__(self):
        if not 0<self.support_angle<math.pi or not isinstance(self.power,int) or self.power<6:
            raise ValueError('a regular finite cap and at least C5 support are required')
        domain=self.model.atlas.parent_causal_domain()
        if self.support_angle+2*domain['future_conformal_distance_upper_bound']>=domain['initial_parent_angle']:
            raise ValueError('the whole packet causal cone must remain inside the protected parent')

    @cached_property
    def polynomial(self):
        c=math.cos(self.support_angle)
        weight=lambda x:(x-c)**self.power*math.sqrt(max(0,1-x*x))
        norm=quad(weight,c,1,epsabs=1e-14)[0]
        mean=quad(lambda x:x*weight(x),c,1,epsabs=1e-14)[0]/norm
        var=quad(lambda x:(x-mean)**2*weight(x),c,1,epsabs=1e-14)[0]/norm
        third=quad(lambda x:(x-mean)**3*weight(x),c,1,epsabs=1e-14)[0]/norm
        # Polynomial in x=cos(theta); all derivatives are regular at the pole.
        P=np.polynomial.Polynomial
        return P([-c,1])**self.power*(P([-mean,1])**2-(third/var)*P([-mean,1])-var)

    @cached_property
    def normalization(self):
        p=self.model.initial;f=self.polynomial
        required=(2/(p['kappa']*p['a']**2))*(-3*f.deriv()(1)+3*f(1))
        E0=p['omega']**2*p['chi']**2/2
        return E0/required

    def profiles(self,theta):
        theta=np.asarray(theta);x=np.cos(theta);f=self.polynomial
        mask=theta<self.support_angle
        Phi=np.where(mask,f(x),0.)*self.normalization
        lap=np.where(mask,(1-x*x)*f.deriv(2)(x)-3*x*f.deriv()(x),0.)*self.normalization
        p=self.model.initial;drho=2*(lap+3*Phi)/(p['kappa']*p['a']**2)
        return dict(Phi=Phi,laplacian=lap,delta_rho=drho,X=drho/(p['omega']**2*p['chi']),
            Delta=drho/(p['omega']**2*p['chi']**2/2))

    def coefficients(self,max_mode=128):
        if max_mode<2 or max_mode>256:raise ValueError('packet cutoff must be 2..256')
        n=np.arange(max_mode+1)
        # A fixed high-order quadrature resolves all admitted Fourier factors;
        # independent lower/higher cutoffs test the retained spatial response.
        x,w=np.polynomial.legendre.leggauss(max(400,3*max_mode))
        theta=(x+1)*self.support_angle/2;weight=w*self.support_angle/2
        F=self.profiles(theta)['Phi']
        coeff=2*(n+1)/math.pi*(np.sin(np.outer(n+1,theta))@(weight*F*np.sin(theta)))
        return coeff

    @staticmethod
    def radial_harmonics(modes,theta):
        theta=np.asarray(theta);modes=np.asarray(modes)
        return eval_chebyu(modes[:,None],np.cos(theta)[None,:])/(modes[:,None]+1)

    def evolve(self,max_mode=128,*,clock='proper',rtol=2e-10,atol=2e-12):
        coefficients=self.coefficients(max_mode);modes=np.arange(2,max_mode+1)
        evolution=self.model.solve(modes,clock=clock,rtol=rtol,atol=atol)
        p=self.model.initial;E0=p['omega']**2*p['chi']**2/2
        initial_delta=-2*(modes*(modes+2)-3)*coefficients[2:]/(p['kappa']*p['a']**2*E0)
        return evolution,coefficients,initial_delta

    def reconstruct(self,evolution,initial_delta,clock_value,theta):
        s=evolution.sample([clock_value]);Q=self.radial_harmonics(evolution.modes,theta)
        return dict(Delta=(initial_delta*s['Delta'][:,0,0])@Q,
                    Phi=(initial_delta*s['Phi'][:,0,0])@Q)


def jeans_comparison(model=ParentModes()):
    """Attributed averaged short-wavelength comparison, not evolution input."""
    p=model.initial;E=p['omega']**2*p['chi']**2/2
    kJ=(2*p['kappa']*E*p['omega']**2)**.25
    lv=math.sqrt(3/(p['kappa']*p['vacuum']))
    return dict(initial_excitation_to_vacuum_ratio=E/p['vacuum'],
        initial_H_over_frequency=p['H']/p['omega'],
        initial_physical_Jeans_wavenumber=kJ,
        initial_Jeans_wavelength_over_vacuum_radius=2*math.pi/(kJ*lv),
        initial_Jeans_wavelength_over_a0=2*math.pi/(kJ*p['a']),
        initial_curvature_over_Jeans_wavenumber_squared=1/(p['a']**2*kJ**2),
        formula='k_J,physical^4=2*kappa*E*omega^2; leading averaged, short-curvature-scale estimate only')


def fixed_metric_control(evolution):
    """Deliberately remove gravitational response, retaining the same background.

    This Klein--Gordon test-field control violates the perturbed Einstein
    constraints and is not another physical solution of the coupled atlas.
    It diagnoses whether the chosen initial field alone produces the growth.
    """
    n=evolution.modes;p=evolution.model.initial
    lam=n*(n+2);w=p['omega']
    initial=np.r_[np.full(len(n),p['chi']/2),np.zeros(len(n))]
    def rhs(t,z):
        _,_,N,H=evolution.solution.sol(t)[:4]
        a=p['a']*math.exp(N);X,Y=z.reshape(2,len(n))
        rates=np.r_[Y,-3*H*Y-(w*w+lam/a**2)*X]
        return rates/H if evolution.clock=='efolds' else rates
    sol=solve_ivp(rhs,(0,evolution.solution.t[-1]),initial,method='DOP853',
        rtol=2e-11,atol=2e-13,max_step=.01 if evolution.clock=='proper' else .003)
    if not sol.success:raise RuntimeError('fixed-metric control failed')
    chi,Pi,_,H=evolution.solution.y[:4,-1]
    E=(Pi*Pi+w*w*chi*chi)/2
    X,Y=sol.y[:,-1].reshape(2,len(n))
    return dict(Delta=(Pi*Y+(w*w*chi+3*H*Pi)*X)/E,
        rhs_calls=sol.nfev,scope='Fixed-background test field only; the gravitational response is deliberately omitted.')
