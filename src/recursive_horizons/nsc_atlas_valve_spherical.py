"""Smooth spherical continuation of the declared one-mode valve preparation.

The action is unchanged. Spherical symmetry is exact for this NEW occupation,
not an average of the old two-mode state. Synchronous proper time and regular
pole variables are used. A future-trapped neck is not identified with a
regular child, an event horizon, or a received-energy fraction.
"""
import math
import time

import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq
from scipy.special import eval_gegenbauer

from .nsc_atlas_parent_scale_selector import VolumeParent
from .nsc_atlas_valve_neck import solve_necks


def gegenbauer_jets(u, modes):
    u=np.atleast_1d(np.asarray(u,dtype=float))
    if np.any(~np.isfinite(u)) or np.any(abs(u)>1):
        raise ValueError('finite u in [-1,1] required')
    values=[]
    for derivative,factor in enumerate((1.,2.,8.)):
        values.append(np.column_stack([
            factor*eval_gegenbauer(l-derivative,1+derivative,u)
            if l>=derivative else np.zeros_like(u) for l in range(modes)]))
    return values


class SphericalValve:
    """Exact SO(3) reduction of Einstein gravity and the two canonical fields.

    h=q² dchi²+b² sin²chi dOmega², u=cos(chi), q=b+(1-u²)c.
    K_rad=K_tan+(1-u²)d enforces regular metric/curvature poles.
    State rows are b,c,K_tan,d,S,chi,Pi_S,Pi_chi. No shell is present.
    """
    def __init__(self,modes=32,*,epsilon=1.,initial_degree=40,budget_seconds=60.):
        if isinstance(modes,bool) or not isinstance(modes,int) or not 12<=modes<=128:
            raise ValueError('integer modes in [12,128] required')
        self.modes=modes
        self.initial_data=solve_necks(initial_degree,(epsilon,))[0]
        self.parent=self.initial_data.parent;self.H=self.parent.H
        self.kappa=self.parent.kappa;self.epsilon=float(epsilon)
        self.sstar=next(x['s'] for x in VolumeParent().vacuum_roots() if x['kind']=='minimum')
        self.budget_seconds=float(budget_seconds)
        if not math.isfinite(self.budget_seconds) or self.budget_seconds<=0:
            raise ValueError('positive finite compute budget required')
        angles=np.pi*np.arange(1,modes+1)/(modes+1)
        self.u=np.cos(angles)
        Y,Yu,Yuu=gegenbauer_jets(self.u,modes)
        self.inverse=np.linalg.inv(Y);self.D=Yu@self.inverse;self.D2=Yuu@self.inverse
        self.weights=4*np.pi**2/(modes+1)*np.sin(angles)**2
        psi=self.initial_data.values(self.u)
        self.initial=np.zeros((8,modes));self.initial[0]=psi**2/self.H
        z=self.parent.mass_squared/self.H**2
        self.initial[5]=self.epsilon*math.sqrt(24/(self.kappa*(z+3)))*self.u
        self.calls=0;self.started=None

    def potential(self,S,chi):
        s=self.sstar+S/math.sqrt(84)
        a=7*np.exp(-12*s);b=-7.5*np.exp(-14*s);c=25/32*np.exp(-24*s)
        mu2=400*np.exp(-12*s)
        return (a+b+c+.5*mu2*chi**2,
                (-12*a-14*b-24*c-6*mu2*chi**2)/math.sqrt(84),mu2*chi)

    def coefficients(self,state):
        state=np.asarray(state,dtype=float)
        if state.shape!=(8,self.modes) or np.any(~np.isfinite(state)):
            raise ValueError('finite 8 by modes state required')
        return state@self.inverse.T

    def jets(self,state,u=None):
        if u is None:
            values=np.asarray(state);first=(self.D@values.T).T;second=(self.D2@values.T).T
        else:
            Y,Yu,Yuu=gegenbauer_jets(u,self.modes);coef=self.coefficients(state)
            values=coef@Y.T;first=coef@Yu.T;second=coef@Yuu.T
        return values,first,second

    def pieces(self,state,u=None):
        u=self.u if u is None else np.atleast_1d(u)
        y,yu,yuu=self.jets(state,None if u is self.u else u)
        b,c,kb,d,S,chi,PS,PC=y
        bu,cu,Su,chi_u=yu[[0,1,4,5]];buu,Suu,chi_uu=yuu[[0,4,5]]
        v=1-u*u;q=b+v*c;qu=bu-2*u*c+v*cu;ka=kb+v*d;K=ka+2*kb
        if np.any(q<=0) or np.any(b<=0):
            raise ValueError('positive spatial metric ended; change foliation before continuing')
        rr=u*b-v*bu
        radial=(-b-3*u*bu+v*buu+(qu/q)*rr)/(q*q*b)
        angular=(b*b+2*b*c+2*u*b*bu+v*(c*c-bu*bu))/(q*q*b*b)
        Rr=-2*radial;Rt=angular-radial
        Rdiff=(b/q*(-2*b*c-u*c*bu-2*c*c-u*b*cu+bu*bu+v*cu*bu)-b*buu-c*c+bu*bu)/(q*q*b*b)
        gradient_square=v*(Su*Su+chi_u*chi_u)/(q*q)
        W,WS,WC=self.potential(S,chi)
        drift=(v*(2*bu/b-qu/q)-3*u)/(q*q)
        lapS=v*Suu/(q*q)+drift*Su;lapC=v*chi_uu/(q*q)+drift*chi_u
        pi_square=PS*PS+PC*PC
        rho=.5*(pi_square+gradient_square)+W
        pr=.5*(pi_square+gradient_square)-W
        pt=.5*(pi_square-gradient_square)-W
        # j is physical energy flux in the increasing-chi unit direction.
        j=np.sqrt(v)/q*(PS*Su+PC*chi_u)
        return dict(y=y,yu=yu,u=u,q=q,b=b,ka=ka,kb=kb,K=K,Rr=Rr,Rt=Rt,Rdiff=Rdiff,
                    grad2=gradient_square,W=W,WS=WS,WC=WC,lapS=lapS,lapC=lapC,
                    rho=rho,pr=pr,pt=pt,j=j,rr=rr,radius=np.sqrt(v)*b)

    def rhs(self,t,flat):
        self.calls+=1
        if self.started is not None and self.calls%50==0 and time.perf_counter()-self.started>self.budget_seconds:
            raise TimeoutError('bounded spherical evolution compute budget exceeded')
        y=np.asarray(flat).reshape(8,self.modes);z=self.pieces(y)
        b,c,kb,d,S,chi,PS,PC=y
        v=1-self.u*self.u;Su,chi_u=z['yu'][[4,5]]
        out=np.array([-kb*b,-kb*c-d*z['q'],z['Rt']+z['K']*kb-self.kappa*z['W'],
                      z['Rdiff']+z['K']*d-self.kappa*(Su*Su+chi_u*chi_u)/z['q']**2,
                      PS,PC,z['K']*PS+z['lapS']-z['WS'],z['K']*PC+z['lapC']-z['WC']])
        # Exact reflection parity; this projects roundoff only, not a physical mode.
        parity=np.array([1,1,1,1,1,-1,1,-1])[:,None]
        return (.5*(out+parity*out[:,::-1])).ravel()

    def evolve(self,endpoint=.3,*,rtol=3e-11,atol=1e-13,max_step=.002):
        if not math.isfinite(endpoint) or not 0<endpoint<=2:
            raise ValueError('finite proper endpoint in (0,2] required')
        self.started=time.perf_counter();self.calls=0
        sol=solve_ivp(self.rhs,(0.,endpoint),self.initial.ravel(),method='DOP853',rtol=rtol,
                      atol=atol,max_step=max_step,dense_output=True)
        if not sol.success:raise RuntimeError(sol.message)
        return sol

    def expansions(self,state,u):
        z=self.pieces(state,u)
        if np.any(abs(z['u'])>=1):raise ValueError('closed spheres require |u|<1')
        curvature=2*z['rr']/(z['q']*z['radius'])
        return -2*z['kb']+curvature,-2*z['kb']-curvature

    def diagnostic(self,t,state,*,validation_order=None):
        order=validation_order or 3*self.modes+7
        angles=np.pi*np.arange(1,order+1)/(order+1)
        u=np.cos(angles);w=4*np.pi**2/(order+1)*np.sin(angles)**2
        z=self.pieces(state,u);measure=z['q']*z['b']**2
        volume=float(w@measure)
        average=lambda f:float(w@(measure*f))/volume
        Ham=z['Rr']+2*z['Rt']+z['K']**2-z['ka']**2-2*z['kb']**2-2*self.kappa*z['rho']
        kbu=z['yu'][2];bu=z['yu'][0]
        geom=2*((u/np.sqrt(1-u*u)-np.sqrt(1-u*u)*bu/z['b'])*(z['ka']-z['kb'])+np.sqrt(1-u*u)*kbu)/z['q']
        Mom=geom-self.kappa*z['j']
        neck=self.pieces(state,[0.]);tp,tm=self.expansions(state,[0.])
        r0=float(neck['radius'][0]);b0=float(neck['b'][0]);kb0=float(neck['kb'][0])
        norm_dr=-(kb0*b0)**2 # r_chi=0 exactly on the reflection neck.
        pressure=average((z['pr']+2*z['pt'])/3);density=average(z['rho'])
        work=float(w@(measure*(z['ka']*z['pr']+2*z['kb']*z['pt'])))
        return dict(proper_time=float(t),modes=self.modes,epsilon=self.epsilon,
            theta_plus_over_H=float(tp[0]/self.H),theta_minus_over_H=float(tm[0]/self.H),
            closed_neck_future_trapped=bool(tp[0]<0 and tm[0]<0),neck_radius=r0,
            neck_radial_scale=float(neck['q'][0]),neck_rho=float(neck['rho'][0]),
            neck_radial_pressure=float(neck['pr'][0]),neck_tangential_pressure=float(neck['pt'][0]),
            neck_timelike_energy_flux=float(neck['j'][0]),
            neck_radius_gradient_norm_squared=norm_dr,
            neck_geometric_mass=0.5*r0*(1-norm_dr),
            volume=volume,half_volume=volume/2,energy=volume*density,half_energy=volume*density/2,
            mean_density=density,mean_isotropic_pressure=pressure,mean_w=pressure/density,
            pressure_work_rate=work,half_pressure_work_rate=work/2,
            Hamiltonian_RMS_over_H_squared=math.sqrt(average(Ham*Ham))/self.H**2,
            momentum_RMS_over_H_squared=math.sqrt(average(Mom*Mom))/self.H**2,
            min_radial_scale=float(min(z['q'])),min_angular_scale=float(min(z['b'])),
            initial_constraint_degree=self.initial_data.degree,validation_order=order,
            scope='Numerical finite continuation of explicitly supplied one-mode occupation; no regular child or inherited preparation selection.')

    def trapped_band(self,state):
        """Locate the first marginal sphere on each side of the neck.

        Returns the connected future-trapped band containing u=0, not an
        event horizon or a new universe. A root is accepted only after both
        expansions are checked in the band and marginality is re-evaluated.
        """
        tp,tm=self.expansions(state,[0.])
        if tp[0]>=0 or tm[0]>=0:return None
        def edge(u):
            a,b=self.expansions(state,[u]);return float(max(a[0],b[0]))
        samples=np.linspace(0.,.999,257)
        values=np.maximum(*self.expansions(state,samples))
        idx=np.flatnonzero(values>=0)
        if not len(idx):return None
        i=int(idx[0]);root=brentq(edge,float(samples[i-1]),float(samples[i]),xtol=2e-13)
        def integrands(u):
            z=self.pieces(state,[u]);measure=4*np.pi*z['q'][0]*z['b'][0]**2*math.sqrt(1-u*u)
            return measure,measure*z['rho'][0],measure*(z['ka'][0]*z['pr'][0]+2*z['kb'][0]*z['pt'][0])
        values=[2*quad(lambda u:integrands(u)[i],0,root,epsabs=2e-11,epsrel=2e-11)[0] for i in range(3)]
        z=self.pieces(state,[root]);a,b=self.expansions(state,[root]);radius=float(z['radius'][0])
        return dict(u_edge=root,volume=values[0],energy=values[1],pressure_work_rate=values[2],
                    marginal_expansion=float(max(a[0],b[0])),other_expansion=float(min(a[0],b[0])),
                    boundary_radius=radius,boundary_rho=float(z['rho'][0]),
                    # At +u the outward direction opposes increasing chi.
                    outward_energy_flux=float(-z['j'][0]),boundary_radial_scale=float(z['q'][0]),
                    scope='Connected trapped band; moving-boundary sweeping, flux and pressure work must be distinguished.')
