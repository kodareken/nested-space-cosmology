"""Uniform inherited-scalar background bounds and a rational light-cone map.

The exact homogeneous source belongs to the protected shared parent. The
comparison bounds its density, Hubble rate and specified FLRW distances;
it does not average perturbations or select an ancestral quantum state.
"""
from __future__ import annotations

from fractions import Fraction as F
import math
import time

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq

from .nsc_atlas_flux_conversion import ScalarFlux
from .nsc_atlas_shared_scalar_slice import SharedScalarSlice


def exact_initial_geometry():
    """Independent unsquared outer junction, using exact inherited densities."""
    R, k, mass = F(19, 10), F(1, 10), F(3, 25)
    Ui, Uo, vacuum = F(7, 32), F(21, 32), F(9, 32)
    beta = -mass/(k*R**2)+(Uo-Ui+3*k*k)*R/(6*k)
    f = 1-2*mass/R-Ui*R*R/3
    H02 = (beta*beta-f)/(R*R)
    curvature = Uo/3-H02
    assert 0 < curvature < vacuum/3
    assert beta-k*R > 0
    return dict(H02=H02, curvature=curvature, vacuum=vacuum,
                excitation=Uo-vacuum, initial_density=Uo)


def uniform_bound(omega=20):
    w = F(str(omega))
    if w <= 0:
        raise ValueError("positive frequency required")
    initial = exact_initial_geometry()
    # Conservative rational constants, verified against the exact geometry.
    hmin, hmax, H0max, Kmax = F(133, 500), F(117, 250), F(443, 1000), F(57, 2500)
    assert hmin*hmin < initial['vacuum']/3-initial['curvature']
    assert hmax*hmax > initial['initial_density']/3
    assert H0max*H0max > initial['H02'] and Kmax > initial['curvature']
    delta = 3*hmax/(2*w)
    variation = H0max-hmin+Kmax/hmin
    if delta >= 1:
        raise ValueError("frequency does not give positive corrected energy")
    gamma = 3*variation/(2*w*(1-delta))
    if gamma >= 1:
        raise ValueError("rational exponential comparison requires Gamma<1")
    lower = (1-gamma)/(1+delta)
    upper = 1/((1-gamma)*(1-delta))
    error = max(1-lower, upper-1)
    hubble_squared_fraction = (initial['excitation']/3)/initial['H02']*error
    if hubble_squared_fraction >= 1:
        raise ValueError("positive relative Hubble comparison lost")
    return dict(omega=w,delta=delta,total_H_variation=variation,Gamma=gamma,
                lower=lower,upper=upper,density_relative_error=error,
                hubble_squared_relative_error=hubble_squared_fraction,**initial)


def _mul(a, b):
    c = [F(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i+j] += x*y
    return c


def _integral(poly, lo, hi):
    return sum(c*(hi**(i+1)-lo**(i+1))/(i+1) for i,c in enumerate(poly))


def _sinc_bounds(x2, terms=6):
    if not 0 <= x2 <= 1 or terms < 2 or terms % 2:
        raise ValueError("squared angle in [0,1] and an even truncation required")
    upper = sum(F((-1)**n, math.factorial(2*n+1))*x2**n for n in range(terms+1))
    lower = upper-F(1,math.factorial(2*terms+3))*x2**(terms+1)
    return lower, upper


def dust_lightcone(*, observer_scale=F(2), redshift=F(1,4), terms=12):
    """Exact rational enclosure of Lambda*D_L^2 in the dust comparison.

    u=1/A makes H_d^2=P(u) cubic. Integrate its binomial polynomial exactly
    and bound the entire remaining series, then enclose sinc with its
    alternating series. No floating quadrature defines these endpoints.
    """
    Ao, z = F(observer_scale), F(redshift)
    if Ao <= 1 or z < 0 or Ao/(1+z) < 1 or not isinstance(terms,int) or terms < 4:
        raise ValueError("future comoving endpoints and nonnegative redshift required")
    d = exact_initial_geometry()
    K, vacuum, matter = d['curvature'], d['vacuum'], d['excitation']
    lo, hi = 1/Ao, (1+z)/Ao
    P = [vacuum/3, F(0), -K, matter/3]
    evaluate = lambda u: sum(c*u**i for i,c in enumerate(P))
    # The declared observer protocol lies on this monotone cubic branch.
    if lo*(matter*lo-2*K) <= 0:
        raise ValueError("positive polynomial derivative not established")
    Pmin, Pmax = evaluate(lo), evaluate(hi)
    B = (Pmin+Pmax)/2
    q = (Pmax-Pmin)/(Pmax+Pmin)
    assert Pmin > 0 and 0 <= q < 1
    x = [1-P[0]/B, F(0), -P[2]/B, -P[3]/B]
    power = [F(1)]
    value = F(0)
    for n in range(terms+1):
        value += F(math.comb(2*n,n),4**n)*_integral(power,lo,hi)
        power = _mul(power,x)
    tail = (hi-lo)*q**(terms+1)/(1-q)
    I_lo, I_hi = value-tail, value+tail
    if I_lo < 0:
        raise ValueError("integral lower bound is not positive")
    eta2lo, eta2hi = K*I_lo**2/B, K*I_hi**2/B
    if eta2hi > 1:
        raise ValueError("this sine enclosure requires squared angle below one")
    sinc_lo = _sinc_bounds(eta2hi)[0]
    sinc_hi = _sinc_bounds(eta2lo)[1]
    factor = vacuum*Ao**2*(1+z)**2/K
    Q_lo = factor*eta2lo*sinc_lo*sinc_lo
    Q_hi = factor*eta2hi*sinc_hi*sinc_hi
    return dict(observer_scale=Ao,redshift=z,terms=terms,binomial_q=q,
                integral_tail=tail,eta_squared_upper=eta2hi,
                lower=Q_lo,upper=Q_hi)


def scalar_lightcone_bound(omega=20):
    b, light = uniform_bound(omega), dust_lightcone()
    r = b['hubble_squared_relative_error']
    # eta_d/sqrt(1-r)<1<pi/2, so sin is increasing and concave on both rays.
    assert light['eta_squared_upper']/(1-r) < 1
    return dict(lower=light['lower']/(1+r),upper=light['upper']/(1-r),
                relative_squared_distance_error=r/(1-r))


def outward_decimal_interval(lower, upper, places=15):
    """Authoritative decimal strings rounded outwards from exact rationals."""
    scale = 10**places
    lo = lower.numerator*scale//lower.denominator
    hi = -((-upper.numerator*scale)//upper.denominator)
    def render(n):
        whole, part = divmod(n,scale)
        return str(whole)+'.'+str(part).zfill(places)
    return [render(lo), render(hi)]


def solve_parent(*, clock='proper', omega=20, rtol=2e-11, short_cycles=None):
    """Independent proper-time/constraint-eliminated expansion clocks.

    The old parent rates remain unchanged. Only compact readouts are returned.
    """
    started = time.process_time()
    atlas = SharedScalarSlice(scalar=ScalarFlux(omega=omega))
    domain = atlas.parent_causal_domain()
    a0, H0, vacuum = (domain[k] for k in ('initial_scale_factor','initial_H','ground_density'))
    chi0 = atlas.field-atlas.scalar.minimum(atlas.parent_flux)
    E0 = omega*omega*chi0*chi0/2
    target = math.log(2)
    calls = 0
    def charge():
        nonlocal calls
        calls += 1
        if calls > 60000 or time.process_time()-started > 90:
            raise RuntimeError('finite parent comparison allowance exhausted')
    if clock == 'proper':
        def rhs(t,y):
            charge()
            return atlas.parent_rates(*y[:4])
        def endpoint(t,y):return y[2]-target
        endpoint.terminal=True;endpoint.direction=1
        end = 2*math.pi*short_cycles/omega if short_cycles else 1.1*target/domain['H_positive_lower_bound']
        sol = solve_ivp(rhs,(0,end),(chi0,0.,0.,H0,0.,0.),method='DOP853',
            events=None if short_cycles else endpoint,dense_output=True,
            rtol=rtol,atol=rtol*.01,max_step=math.pi/(8*omega))
        if not sol.success or (not short_cycles and not len(sol.t_events[0])):
            raise RuntimeError(sol.message)
        points=np.linspace(0,sol.t[-1],1025)
        chi,Pi,N,H,work,eta=sol.sol(points)
        def sample_by_N(n):
            if abs(n-sol.y[2,-1])<1e-12:
                return sol.t[-1],sol.y[:,-1]
            t=brentq(lambda t:sol.sol(t)[2]-n,0,sol.t[-1],xtol=1e-13)
            return t,sol.sol(t)
    elif clock == 'efolds' and short_cycles is None:
        def rhs(N,y):
            charge()
            chi,Pi,t,eta,work=y
            E=(Pi*Pi+omega*omega*chi*chi)/2
            H=math.sqrt((vacuum+E)/3-math.exp(-2*N)/(a0*a0))
            return Pi/H,-3*Pi-omega*omega*chi/H,1/H,math.exp(-N)/(a0*H),3*Pi*Pi
        sol=solve_ivp(rhs,(0,target),(chi0,0.,0.,0.,0.),method='DOP853',dense_output=True,
            rtol=rtol,atol=rtol*.01,max_step=.003)
        if not sol.success:raise RuntimeError(sol.message)
        N=np.linspace(0,target,1025)
        chi,Pi,points,eta,work=sol.sol(N)
        E=(Pi*Pi+omega*omega*chi*chi)/2
        H=np.sqrt((vacuum+E)/3-np.exp(-2*N)/(a0*a0))
        def sample_by_N(n):
            x,p,t,eta,work=sol.sol(n)
            h=math.sqrt((vacuum+(p*p+omega*omega*x*x)/2)/3-math.exp(-2*n)/(a0*a0))
            return t,np.array([x,p,n,h,work,eta])
    else:
        raise ValueError('proper or expansion clock required; short control uses proper time')
    E=(Pi*Pi+omega*omega*chi*chi)/2
    scaled=np.exp(3*N)*E/E0
    S=np.exp(3*N)*chi*Pi
    corrected=scaled+1.5*H*S/E0
    b=uniform_bound(omega)
    if np.min(scaled)<float(b['lower']) or np.max(scaled)>float(b['upper']):
        raise RuntimeError('unaveraged source escaped the proved energy interval')
    result=dict(clock=clock,omega=omega,short_cycles=short_cycles,rhs_calls=calls,
        process_cpu_seconds=time.process_time()-started,final_scale=float(math.exp(N[-1])),
        final_proper_time=float(points[-1]),final_comoving_energy_ratio=float(scaled[-1]),
        sampled_comoving_energy_extrema=[float(np.min(scaled)),float(np.max(scaled))],
        sampled_corrected_energy_extrema=[float(np.min(corrected)),float(np.max(corrected))],
        maximum_Friedmann_residual=float(np.max(abs(3*(H*H+np.exp(-2*N)/(a0*a0))-vacuum-E))),
        maximum_energy_work_residual=float(np.max(abs(E+work-E0))))
    if short_cycles is None:
        te,ye=sample_by_N(math.log(8/5))
        to,yo=sample_by_N(target)
        optical_eta=yo[5]-ye[5]
        Q=vacuum*(a0*2*1.25*math.sin(optical_eta))**2
        if clock=='proper':
            alternative=quad(lambda t:1/(a0*math.exp(sol.sol(t)[2])),te,to,
                epsabs=1e-13,epsrel=1e-12)[0]
        else:
            alternative=quad(lambda n:1/(a0*math.exp(n)*sample_by_N(n)[1][3]),
                math.log(8/5),target,epsabs=1e-13,epsrel=1e-12)[0]
        def energy_excess(n):
            _,y=sample_by_N(n)
            return (y[1]**2+omega*omega*y[0]**2)/2-vacuum
        Neq=brentq(energy_excess,0,target,xtol=1e-13)
        result.update(Lambda_DL_squared=Q,photon_conformal_distance=optical_eta,
            independent_photon_integral_gap=abs(optical_eta-alternative),
            energy_equality_scale=math.exp(Neq))
        prediction=scalar_lightcone_bound(omega)
        if not float(prediction['lower'])<Q<float(prediction['upper']):
            raise RuntimeError('unaveraged photon result escaped the analytic interval')
    return result
