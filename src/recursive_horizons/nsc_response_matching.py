"""Ultrastatic Euclidean/retarded Dirac response in an explicit radius channel.

The proper-time and Abel regulators are kept distinct at finite cutoff.
Their periodic-minus-antiperiodic continuum limit is a matching control,
not a completed finite-cutoff gravitational influence functional.
"""
from dataclasses import replace
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import eigh, solve
from scipy.special import erfc
from .nsc_covariant_operator import SIGMA1, SIGMA2, euclidean_operator
from .nsc_influence import canonical_hamiltonian, geometric_vertex


def _ultrastatic(metric):
    if not np.all(metric.lapse == 1):
        raise ValueError('matching control requires unit lapse')


def radius_profile(metric):
    return np.exp(2*(np.cos(2*np.pi*(metric.x+.6)/metric.length)-1))


def radius_path(metric, amplitude, coordinate='inverse_radius'):
    s=radius_profile(metric)
    if coordinate=='inverse_radius':
        radius=metric.sphere_radius/(1+amplitude*s)
    elif coordinate=='log_radius':
        radius=metric.sphere_radius*np.exp(amplitude*s)
    else: raise ValueError('unknown metric coordinate')
    return replace(metric,sphere_radius=radius)


def canonical_data(metric,kappa=1):
    _ultrastatic(metric)
    h=canonical_hamiltonian(metric,kappa)
    e,u=eigh(h,driver='evd')
    if min(abs(e))<1e-9:raise ValueError('zero-mode state requires separate treatment')
    v=geometric_vertex(metric,kappa)
    w=np.kron(SIGMA1,np.diag(kappa*radius_profile(metric)**2/metric.sphere_radius))
    return {'energies':e,'vertex_squared':abs(u.conj().T@v@u)**2,
            'contact_diagonal':np.diag(u.conj().T@w@u).real}


def static_proper_time(data,cutoff):
    """One representative angular block, with analytic frequency integration."""
    if not np.isfinite(cutoff) or cutoff<=0:raise ValueError('positive cutoff required')
    e=data['energies'];v2=data['vertex_squared']
    energy=.5*np.sum(cutoff/np.sqrt(np.pi)*np.exp(-(e/cutoff)**2)
                       -abs(e)*erfc(abs(e)/cutoff))
    fp=-.5*np.sign(e)*erfc(abs(e)/cutoff)
    delta=e[:,None]-e[None,:];mid=(e[:,None]+e[None,:])/2
    near=(abs(delta)<1e-9)&(e[:,None]*e[None,:]>0)
    divided=np.divide(fp[:,None]-fp[None,:],delta,out=np.zeros_like(delta),where=~near)
    divided[near]=np.exp(-(mid[near]/cutoff)**2)/(np.sqrt(np.pi)*cutoff)
    bubble=np.sum(divided*v2)
    contact=np.dot(fp,data['contact_diagonal'])
    return {'energy':float(energy),'inverse_radius_hessian':float(bubble),
            'log_radius_contact':float(contact),'log_radius_hessian':float(bubble+contact)}


def retarded_susceptibility(data,z,abel_cutoff=None):
    """Upper-half-plane vacuum Kubo response, H_int=+J V, no contact term.

    The optional positive transition weight is an Abel summation regulator.
    It is not identified with the finite proper-time action regulator.
    """
    if not np.isfinite(z) or np.imag(z)<0:raise ValueError('upper-half-plane energy required')
    e=data['energies'];lo=e<0;hi=e>0
    gap=e[hi][None,:]-e[lo][:,None]
    weight=data['vertex_squared'][np.ix_(lo,hi)].copy()
    if abel_cutoff is not None:
        if not np.isfinite(abel_cutoff) or abel_cutoff<=0:raise ValueError('positive Abel cutoff required')
        weight*=np.exp(-(e[lo][:,None]**2+e[hi][None,:]**2)/abel_cutoff**2)
    if np.any(abs(z-gap)<1e-12) or np.any(abs(z+gap)<1e-12):
        raise ValueError('real-axis pole requires a distributional prescription')
    return complex(np.sum(weight*(1/(z-gap)-1/(z+gap))))


def euclidean_kernel(metric,nu,cutoff,kappa=1,frequency_points=64):
    """One angular block's proper-time quadratic kernel at external |nu|.

    Differentiate Tr g(D), g(d)=E1(d²/Lambda²)/2, between the two
    internal frequencies omega±nu/2. The inverse-radius coordinate is affine
    in D, so its explicit second operator variation vanishes.
    """
    _ultrastatic(metric)
    if not np.isfinite(nu) or nu<0 or not np.isfinite(cutoff) or cutoff<=0:
        raise ValueError('nonnegative frequency and positive cutoff required')
    if frequency_points<16:raise ValueError('frequency quadrature needs at least16 nodes')
    v=-np.kron(SIGMA2,np.diag(kappa*radius_profile(metric)/metric.sphere_radius))
    nodes,weights=leggauss(frequency_points);extent=8*cutoff+nu
    value=0.
    for omega,weight in zip((nodes+1)*extent/2,weights*extent/2):
        e,u=eigh(euclidean_operator(metric,omega-nu/2,kappa),driver='evd')
        f,t=eigh(euclidean_operator(metric,omega+nu/2,kappa),driver='evd')
        if min(min(abs(e)),min(abs(f)))<1e-9:raise ValueError('unresolved zero mode')
        ep=-np.exp(-(e/cutoff)**2)/e;fp=-np.exp(-(f/cutoff)**2)/f
        delta=e[:,None]-f[None,:];near=(abs(delta)<1e-8)&(e[:,None]*f[None,:]>0)
        divided=np.divide(ep[:,None]-fp[None,:],delta,out=np.zeros_like(delta),where=~near)
        mid=(e[:,None]+f[None,:])/2
        divided[near]=np.exp(-(mid[near]/cutoff)**2)*(1/mid[near]**2+2/cutoff**2)
        value+=weight*np.sum(divided*abs(u.conj().T@v@t)**2)/np.pi
    return float(value)


def rational_frequency_bubble(metric,nu,kappa=1,frequency_points=96):
    """Independent unregulated finite-space frequency loop from matrix solves.

    omega=t/(1-t) integrates the infinite frequency interval. This is the
    ordinary finite Dirac determinant, not the finite proper-time action.
    """
    _ultrastatic(metric)
    if not np.isfinite(nu) or nu<0:raise ValueError('nonnegative Euclidean frequency required')
    v=-np.kron(SIGMA2,np.diag(kappa*radius_profile(metric)/metric.sphere_radius))
    nodes,weights=leggauss(frequency_points)
    value=0j
    for t,weight in zip((nodes+1)/2,weights/2):
        omega=t/(1-t)
        left=solve(euclidean_operator(metric,omega-nu/2,kappa),v,assume_a='her')
        right=solve(euclidean_operator(metric,omega+nu/2,kappa),v,assume_a='her')
        value+=weight*np.trace(left@right)/((1-t)**2*np.pi)
    if abs(value.imag)>1e-9:raise RuntimeError('frequency loop lost its real symmetry')
    return float(value.real)


def abel_extrapolation(periodic,antiperiodic,nu,cutoff):
    cuts=[cutoff,cutoff/np.sqrt(2),cutoff/2]
    values=[(retarded_susceptibility(periodic,1j*nu,c)
             -retarded_susceptibility(antiperiodic,1j*nu,c)).real for c in cuts]
    return {'cutoffs':cuts,'differences':values,
            'extrapolated':float((8*values[0]-6*values[1]+values[2])/3),
            'remainder_bound_proved':False}
