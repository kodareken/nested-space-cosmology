"""Short-shell Fourier response in the owned unitary KS radial frame.

Real energies label original input fibers and output Fourier arguments.
Neither source occupations nor energy quadrature weights enter this ODE.
No physical preparation accuracy is inferred from its numerical tolerances.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp

from .nsc_compatible_history_geometry import CompatibleRadiusDirection
from .nsc_ks_spacetime_variation import chart_coordinates
from .nsc_lorentzian import geometry
from .nsc_common_time_bulk_split import CommonTimeBulkSplit
from .nsc_transmitting_dirac_domain import TransmittingDiracSeamDomain,MODE_TO_CURRENT,I2,S1,S2,S3
from .nsc_transmitting_resolvent import profile


def canonical_map(rho):
    beta=float(geometry(rho)[0]);r=float(np.sqrt(1+rho*rho))
    domain=TransmittingDiracSeamDomain(lapse=1.,radial_scale=1.,shift=beta,radius=r)
    _,inverse,_=domain.trace_map(np.ones(1))
    return inverse[0]/r


def archived_amplitudes(x,phi,energies,index):
    """Spin-major characteristic archive -> unweighted KS source amplitudes."""
    x=np.asarray(x);E=np.asarray(energies);phi=np.asarray(phi)
    if E.shape!=(4,) or phi.shape!=(2*len(x),12) or not 0<=index<len(x):
        raise ValueError('four energies and three archived source columns required')
    current=MODE_TO_CURRENT@phi[[index,len(x)+index]]
    mapped=(canonical_map(float(x[index]))@current).reshape(2,4,3).transpose(1,0,2)
    return mapped*np.exp(1j*E*chart_coordinates(float(x[index]))[1])[:,None,None]


def pulse_fourier(omega,center,halfwidth,order):
    """Integral exp(+i omega z) w(z) dz; inverse measure is d omega/(2 pi)."""
    omega=np.asarray(omega)
    if np.iscomplexobj(omega) or not np.isfinite(omega).all() or halfwidth<=0 or order not in (96,192):
        raise ValueError('real finite Fourier arguments and fixed96/192 rules required')
    nodes,weights=leggauss(order)
    bump=np.array([profile(float(u))[0] for u in nodes])
    return halfwidth*np.exp(1j*omega*center)*np.sum(
        np.exp(1j*omega[...,None]*halfwidth*nodes)*weights*bump,axis=-1)


def ks_generator(rho,energies,mass,angular):
    beta=float(geometry(rho)[0]);a=np.sqrt(beta*beta-1);r=np.sqrt(1+rho*rho)
    E=np.asarray(energies)
    if np.iscomplexobj(E) or not np.isfinite(E).all():raise ValueError('finite real Fourier energies required')
    H=-mass*S1+angular/r*S2-E[...,None,None]/a*S3
    return 1j/a*H


def local_frame_checks(radii,energies,mass,angular):
    """Replay owned frame connection and pure-radius forcing identities."""
    generator=forcing=unitary=0.
    for rho in radii:
        beta,bp,_=map(float,geometry(rho));a=np.sqrt(beta*beta-1);r=np.sqrt(1+rho*rho)
        R=canonical_map(rho);F=np.linalg.inv(R);v=S2-beta*I2
        connection=(-beta*bp*I2-bp*S2)/(2*a*a)
        M=CommonTimeBulkSplit.local_potential(1.,r,mass,angular)
        forcing=max(forcing,float(np.max(abs(R@np.linalg.solve(v,S3)@F-S2/a))))
        for E in energies:
            G=ks_generator(rho,E,mass,angular)
            pg=np.linalg.solve(v,1j*(E*I2-M)+.5*bp*I2)
            mapped=connection-1j*E*beta/(a*a)*I2+F@G@R
            generator=max(generator,float(np.max(abs(pg-mapped))))
            unitary=max(unitary,float(np.max(abs(G+G.conj().T))))
    return {'PG_KS_generator':generator,'radius_forcing_map':forcing,'KS_antiHermitian_generator':unitary}


def short_response(upstream,energies,mass,angular,initial,what,*,rtol,atol):
    """One coupled solve for four input columns and sixteen Fourier responses."""
    E=np.asarray(energies);A0=np.asarray(initial,complex);W=np.asarray(what,complex)
    if E.shape!=(4,) or A0.shape!=(4,2,3) or W.shape!=(4,4) or upstream<=1:
        raise ValueError('four source fibers, sixteen responses and upstream rho>1 required')
    if not all(np.isfinite(v).all() for v in (E,A0,W)):raise ValueError('finite response data required')
    # This constant callback produces exactly s*chi(s), reusing the frozen
    # supplied radius basis rather than implementing another cutoff.
    normal=CompatibleRadiusDirection(lambda z,n:1. if n==0 else 0.,lambda z,n:0.,.007,.03)
    T1=chart_coordinates(1.)[0]
    def rhs(rho,state):
        A=state[:24].reshape(4,2,3);B=state[24:].reshape(4,4,2,3)
        G=ks_generator(rho,E,mass,angular)
        beta=float(geometry(rho)[0]);a=np.sqrt(beta*beta-1);r=np.sqrt(1+rho*rho)
        s=chart_coordinates(float(rho))[0]-T1
        coefficient=-1j*angular*normal.value(s,0.)/(a*r*r)
        dA=np.einsum('iab,ibc->iac',G,A)
        dB=np.einsum('oab,oibc->oiac',G,B)+coefficient*W[:,:,None,None]*np.einsum('ab,ibc->iac',S2,A)[None,:,:,:]
        return np.concatenate((dA.ravel(),dB.ravel()))
    initial_state=np.concatenate((A0.ravel(),np.zeros(96,complex)))
    solution=solve_ivp(rhs,(float(upstream),1.),initial_state,method='DOP853',rtol=rtol,atol=atol)
    if not solution.success:raise ArithmeticError(solution.message)
    end=solution.y[:,-1]
    return end[:24].reshape(4,2,3),end[24:].reshape(4,4,2,3),{'function_evaluations':int(solution.nfev),'accepted_steps':len(solution.t)-1,'rtol':rtol,'atol':atol}
