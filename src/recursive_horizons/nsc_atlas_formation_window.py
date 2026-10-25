"""Ambient formation budget in the same averaged inherited parent.

u=1/A makes the whole expanding future a regular finite interval. This
linear calculation does not install a core or assert nonlinear collapse.
"""
from dataclasses import dataclass
import math

import numpy as np
from scipy.integrate import solve_ivp,quad
from scipy.special import roots_legendre,spherical_jn,beta

from .nsc_atlas_scale_hierarchy import ScaleHierarchy


@dataclass(frozen=True)
class FormationBackground:
    vacuum: float=9/32
    excitation: float=3/8
    curvature: float=.0227765023676067
    frequency: float=1e9

    @classmethod
    def inherited(cls):
        m=ScaleHierarchy(1e9,.001,.5)
        return cls(m.vacuum,m.reference.scalar.kappa*m.excitation,m.curvature,m.frequency)

    def H(self,u):
        F=self.vacuum/3+self.excitation*np.asarray(u)**3/3-self.curvature*np.asarray(u)**2
        if np.any(F<=0):raise ValueError('background leaves the expanding positive-H branch')
        return np.sqrt(F)

    def minimum_H_squared(self,u_initial):
        critical=2*self.curvature/self.excitation if self.excitation else -1.
        points=[0.,u_initial]+([critical] if 0<critical<u_initial else [])
        return min(float(self.H(u)**2) for u in points)


def transfer(k=0.,*,initial_scale=1.,rtol=2e-11,background=None,samples=401):
    """Columns: delta_initial and (delta_dot/H)_initial.

    The auxiliary momentum is p=A^2 delta_dot. No k=0 physical global
    mode is inferred: k=0 supplies the positive comparison solution.
    """
    b=background or FormationBackground.inherited()
    if not (math.isfinite(k) and k>=0 and .001<=initial_scale<=1.):
        raise ValueError('outside the declared local-mode/initial-epoch domain')
    ui=1/initial_scale;Hi=float(b.H(ui));nu=k*k/(2*b.frequency)
    if b.minimum_H_squared(ui)<=0:raise ValueError('future contains a turnaround')
    # Rescale momentum by a constant to resolve oscillatory high-k cases.
    pscale=max(1.,nu,initial_scale**2*Hi)
    def rhs(u,y):
        H=float(b.H(u));Y=y.reshape(2,2)
        mat=np.array([[0.,-u*pscale/H],[(nu*nu*u-b.excitation/2)/(H*pscale),0.]])
        return (mat@Y).ravel()
    start=np.diag([1.,initial_scale**2*Hi/pscale])
    sol=solve_ivp(rhs,(ui,0),start.ravel(),method='DOP853',rtol=rtol,atol=rtol*.01,dense_output=True)
    if not sol.success:raise RuntimeError(sol.message)
    u=np.linspace(ui,0,samples);Y=sol.sol(u).T.reshape(-1,2,2)
    determinant=np.linalg.det(Y);target=np.linalg.det(start)
    # Restore true p in the second row, while keeping initial column units.
    Y[:,1,:]*=pscale
    return dict(k=float(k),initial_scale=initial_scale,
        future_density_response=Y[-1,0,:].tolist(),
        future_momentum_response=Y[-1,1,:].tolist(),
        sampled_max_density_response=float(np.max(abs(Y[:,0,0]))),
        sampled_max_velocity_response=float(np.max(abs(Y[:,0,1]))),
        determinant_relative_error=float(np.max(abs(determinant/target-1))),
        minimum_H_squared=b.minimum_H_squared(ui),
        initial_H_over_frequency=Hi/b.frequency,
        initial_k_over_frequency_A=k/(b.frequency*initial_scale),
        Jeans_crossing_scale=k**4/(2*b.frequency**2*b.excitation) if b.excitation>0 else None,
        function_evaluations=sol.nfev),dict(u=u,response=Y,solution=sol,pscale=pscale)


def dust_by_quadrature(initial_scale=1.,background=None):
    """Independent H and H*integral growing/decaying basis, no ODE solve."""
    b=background or FormationBackground.inherited();u=1/initial_scale
    H=float(b.H(u));Hu=(b.excitation*u*u-2*b.curvature*u)/(2*H)
    # Growing solution H(u) int_u^infinity v/H(v)^3 dv. Overall
    # normalization cancels in the initial/final fundamental matrix.
    I=quad(lambda v:v/float(b.H(v))**3,u,np.inf,epsabs=1e-11,epsrel=2e-12)[0]
    I0=quad(lambda v:v/float(b.H(v))**3,0,np.inf,epsabs=1e-11,epsrel=2e-12)[0]
    G=H*I;Gu=Hu*I-u/H**2
    # p=-H delta_u/u, and initial v=delta_dot/H=-u delta_u.
    initial=np.array([[H,G],[-u*Hu,-u*Gu]])
    final=np.array([math.sqrt(b.vacuum/3),math.sqrt(b.vacuum/3)*I0])
    return final@np.linalg.inv(initial)


def density_phase(k,*,initial_scale=1.,rtol=2e-11,background=None):
    """Independent continuity/Euler phase system; initial current is zero.

    delta_dot=k^2 theta/(omega A^2),
    theta_dot=-k^2 delta/(4 omega A^2)+omega D delta/(2 A k^2).
    """
    b=background or FormationBackground.inherited()
    if k<=0:raise ValueError('phase check needs nonzero k')
    def rhs(u,y):
        H=float(b.H(u));delta,theta=y
        return (-k*k*u*theta/(b.frequency*H),
            (k*k*u/(4*b.frequency)-b.frequency*b.excitation/(2*k*k))*delta/H)
    sol=solve_ivp(rhs,(1/initial_scale,0),(1.,0.),method='DOP853',rtol=rtol,atol=rtol*.01)
    if not sol.success:raise RuntimeError(sol.message)
    return float(sol.y[0,-1])


def elementary_gain_bound(initial_scale=1.,background=None):
    """Volterra/Gronwall upper bound, no density-transfer integration."""
    b=background or FormationBackground.inherited();ui=1/initial_scale
    minimum=b.minimum_H_squared(ui)
    exponent=b.excitation*ui**3/(12*minimum)
    return dict(log_density_gain_bound=exponent,
        density_gain_bound=math.exp(exponent) if exponent<700 else None,
        explanation='|delta(t)| <= exp[D ui^3/(12 Hmin^2)] |delta_i| for zero initial current, in the linear envelope at every k. A loose early-time bound may overflow; use the exact positive dust comparison instead.')


def packet_tail_fraction(Q):
    """Analytic bound beyond Q for the p=4 compact autocorrelation seed.

    G/G0=10395*j5(q)/q^5; bound its explicit sin/cos polynomial.
    The normalization is Parseval's exact integral of |grad g|^2.
    """
    if Q<1:raise ValueError('tail estimate requires Q>=1')
    coeff=np.array([1.,15.,105.,420.,945.,945.])
    square=np.convolve(coeff,coeff)
    integral=4*np.pi*16*beta(2.5,7)/beta(1.5,5)**2
    tail=10395.**2*sum(c*Q**(-7-j)/(7+j) for j,c in enumerate(square))
    return float(tail/integral)


def compact_packet(length,*,initial_scale=.01,order=128,qmax=40.,rtol=2e-10,background=None):
    """delta_i=-Laplacian(g*g) normalized to central value one.

    g=(1-r^2/L^2)^4_+; support(delta_i)=2L and its Fourier weight is
    nonnegative. Finite q quadrature is NOT renormalized to hide its tail.
    """
    b=background or FormationBackground.inherited();z,w=roots_legendre(order)
    q=(z+1)*qmax/2;w=w*qmax/2
    norm=4*np.pi*16*beta(2.5,7)/beta(1.5,5)**2
    # At small q the power series avoids division of tiny spherical j5.
    G=10395*spherical_jn(5,q)/q**5
    weights=w*q**4*G*G/norm;k=q/length;nu=k*k/(2*b.frequency);scale=np.maximum(1.,nu)
    def rhs(u,y):
        d,p=y.reshape(2,-1);H=float(b.H(u))
        return np.array([-u*scale*p/H,(nu*nu*u-b.excitation/2)*d/(H*scale)]).ravel()
    start=np.r_[np.ones(order),np.zeros(order)]
    sol=solve_ivp(rhs,(1/initial_scale,0),start,method='DOP853',rtol=rtol,atol=rtol*.01,dense_output=True)
    if not sol.success:raise RuntimeError(sol.message)
    us=np.r_[np.geomspace(1/initial_scale,.05,181),0.]
    density=sol.sol(us)[:order];central=weights@density
    dust=dust_by_quadrature(initial_scale,b)[0];tail=packet_tail_fraction(qmax)
    return dict(length=length,initial_scale=initial_scale,quadrature_order=order,qmax=qmax,
        initial_normalization=float(np.sum(weights)),exact_normalization_integral=float(norm),
        central_future_gain=float(central[-1]),initial_omitted_spectral_weight_bound=tail,
        all_future_pointwise_tail_error_bound=float(dust*tail),
        maximum_retained_k_over_frequency_A=float(qmax/(length*b.frequency*initial_scale)),
        initial_H_over_frequency=float(b.H(1/initial_scale)/b.frequency),
        function_evaluations=sol.nfev),dict(u=us,central=central,solution=sol,weights=weights)
