"""Normalized finite Gaussian fermionic influence functional for geometry.

H is derived from the covariant static Dirac operator. Bare finite-space
sources and time-smeared matter kernels are distinct from a completed
renormalized gravitational action. No extra gravitational weight is added.
"""
from __future__ import annotations

import numpy as np
import warnings
from scipy.linalg import eigh, expm, lu_factor, LinAlgWarning

from .nsc_covariant_operator import SIGMA1, SIGMA3, euclidean_operator


def canonical_hamiltonian(metric,kappa=1):
    """D_omega=N^-1/2 beta(omega+iH)N^-1/2 in the finite representation."""
    beta=np.kron(SIGMA3,np.eye(metric.points))
    h0=-1j*beta@euclidean_operator(metric,0.,kappa)
    root=np.sqrt(np.r_[metric.lapse,metric.lapse])
    return root[:,None]*h0*root[None,:]


def geometric_vertex(metric,kappa=1):
    """Exact H tangent for r(t,x)=r0/(1+J(t)s(x)), with N,q fixed."""
    profile=np.exp(2*(np.cos(2*np.pi*(metric.x+.6)/metric.length)-1))
    return np.kron(SIGMA1,np.diag(metric.lapse*kappa*profile/metric.sphere_radius))


def ground_covariance(h):
    e,u=eigh(h,driver='evd')
    if min(abs(e))<1e-10:raise ValueError('zero-mode occupation must be specified')
    return e,u,(u[:,e<0]@u[:,e<0].conj().T)


def thermal_covariance(h,inverse_temperature,chemical_potential=0.):
    from scipy.special import expit
    if not np.isfinite(inverse_temperature) or inverse_temperature<0:
        raise ValueError('finite nonnegative inverse temperature required')
    e,u=eigh(h,driver='evd')
    c=(u*expit(-inverse_temperature*(e-chemical_potential))[None,:])@u.conj().T
    return c


def _covariance(c):
    c=np.asarray(c,dtype=complex)
    if c.ndim!=2 or c.shape[0]!=c.shape[1] or not np.allclose(c,c.conj().T,rtol=0,atol=1e-11):
        raise ValueError('Hermitian one-particle covariance required')
    values=np.linalg.eigvalsh(c)
    if min(values)<-1e-11 or max(values)>1+1e-11:
        raise ValueError('fermionic covariance requires 0<=C<=I')
    return c


def influence(c,u_plus,u_minus):
    """Tr(U_plus rho U_minus†)=det(I-C+C U_minus† U_plus)."""
    c=_covariance(c);plus=np.asarray(u_plus);minus=np.asarray(u_minus)
    if plus.shape!=c.shape or minus.shape!=c.shape:raise ValueError('history dimensions differ')
    identity=np.eye(len(c))
    for unitary in (plus,minus):
        if not np.allclose(unitary.conj().T@unitary,identity,atol=1e-10,rtol=0):
            raise ValueError('unitary branch evolution required')
    kernel=identity-c+c@minus.conj().T@plus
    # Sum LU log-diagonals so tiny nonzero overlaps retain their action even
    # if the displayed amplitude underflows. Exact singularity has no log.
    with warnings.catch_warnings(record=True) as notices:
        warnings.simplefilter('always',LinAlgWarning)
        lu,pivots=lu_factor(kernel,check_finite=True)
    diagonal=np.diag(lu)
    singular=bool(np.any(diagonal==0))
    if notices and not singular:raise ArithmeticError('unexpected influence LU warning')
    if singular:
        z=0j;log_abs=-np.inf;phase=None
    else:
        log_abs=float(np.sum(np.log(abs(diagonal))))
        angle=float(np.sum(np.angle(diagonal))+np.pi*np.count_nonzero(pivots!=np.arange(len(c))))
        phase=float(np.angle(np.exp(1j*angle)))
        z=complex(np.exp(log_abs)*np.exp(1j*phase))
    return {'amplitude':z,'log_modulus':float(log_abs),
            'principal_phase':phase,
            'principal_action':complex(phase,-log_abs) if not singular else None,
            'amplitude_underflow':bool(not singular and z==0)}


def heisenberg(h,a,time):
    u=expm(-1j*time*h)
    return u.conj().T@a@u


def mean(c,a):return complex(np.trace(c@a))


def connected(c,a,b):
    """Wick covariance <A B>-<A><B> for number-conserving bilinears."""
    return complex(np.trace(c@a@(np.eye(len(c))-c)@b))


def retarded_response(h,c,a,b,time,source_time=0.):
    """For H_int=+J B: chi_AB=-i theta(t-s) <[A(t),B(s)]>."""
    if time<=source_time:return 0.
    at=heisenberg(h,a,time);bs=heisenberg(h,b,source_time)
    return float((-1j*np.trace(c@(at@bs-bs@at))).real)


def noise_kernel(h,c,a,times):
    c=_covariance(c)
    operators=[heisenberg(h,a,float(t)) for t in times]
    return np.array([[connected(c,at,bs).real for bs in operators] for at in operators])


def gaussian_smeared_vertex(h,vertex,width):
    """Integral dt exp(-t²/(2width²)) V_H(t), not a time-ordered propagator."""
    if not np.isfinite(width) or width<=0:raise ValueError('positive finite width required')
    e,u=eigh(h,driver='evd');w=u.conj().T@vertex@u
    frequency=e[:,None]-e[None,:]
    smeared=np.sqrt(2*np.pi)*width*np.exp(-.5*(width*frequency)**2)*w
    return u@smeared@u.conj().T


def response_spectrum(h,c,vertex):
    """Positive-frequency weights for a stationary Gaussian covariance.

    noise_plus is the positive-frequency symmetric-noise delta weight/pi;
    minus_im_chi_plus is the matching -Im chi_R weight/pi. The distributions
    themselves need time/frequency smearing before pointwise comparison.
    """
    c=_covariance(c);e,u=eigh(h,driver='evd');rotated=u.conj().T@c@u
    if np.max(abs(rotated-np.diag(np.diag(rotated))))>1e-9:
        raise ValueError('spectral stationary response requires [H,C]=0')
    n=np.diag(rotated).real;v=u.conj().T@vertex@u
    rows=[]
    for i in range(len(e)):
        for j in range(i+1,len(e)):
            frequency=e[j]-e[i]
            if frequency<1e-10:continue
            weight=float(abs(v[i,j])**2)
            rows.append({'frequency':float(frequency),
                         'noise_plus':float((n[i]*(1-n[j])+n[j]*(1-n[i]))*weight),
                         'minus_im_chi_plus':float((n[i]-n[j])*weight)})
    return rows


def zero_frequency_noise(h,c,vertex):
    """Stationary elastic noise, separate from positive-frequency FDT weights."""
    c=_covariance(c);e,u=eigh(h,driver='evd');rotated=u.conj().T@c@u
    if np.max(abs(rotated-np.diag(np.diag(rotated))))>1e-9:raise ValueError('stationary covariance required')
    n=np.diag(rotated).real;v=u.conj().T@vertex@u
    same=abs(e[:,None]-e[None,:])<1e-10
    return float(np.sum(n[:,None]*(1-n[None,:])*abs(v)**2*same))


def with_real_branch_action(amplitude,action_plus,action_minus):
    """Real local action differences preserve normalization and overlap modulus.

    This demonstrates the phase/matching freedom; the coefficients of such
    terms must come from the common action, not from a desired source value.
    """
    if not np.isfinite(action_plus) or not np.isfinite(action_minus):raise ValueError('finite real actions required')
    return complex(amplitude*np.exp(1j*(action_plus-action_minus)))


def pulse_unitary(h,vertex,epsilon=.02,width=.5,steps=256):
    """Exact-unitary midpoint products in the free interaction frame.

    Independent of the previous occupied-mode DOP853 solver. The finite
    Gaussian window is +/-7width; quadrature error is controlled by steps.
    """
    if (isinstance(steps,bool) or not isinstance(steps,int) or steps<8
            or not np.isfinite(width) or width<=0 or not np.isfinite(epsilon)):
        raise ValueError('finite amplitude, positive width and at least eight steps required')
    e,basis=eigh(h,driver='evd');v=basis.conj().T@vertex@basis
    dt=14*width/steps;unitary=np.eye(len(h),dtype=complex)
    for index in range(steps):
        time=-7*width+(index+.5)*dt
        phase=np.exp(1j*e*time)
        vi=phase[:,None]*v*phase.conj()[None,:]
        unitary=expm(-1j*epsilon*np.exp(-.5*(time/width)**2)*dt*vi)@unitary
    return basis@unitary@basis.conj().T


def branch_unitary(h,kicks,final_time):
    """Ordered source probes, each (time, Hermitian vertex, amplitude).

    Coincident probes retain input order. They test functional derivatives,
    not an assertion that a delta-function metric is a physical geometry.
    """
    if not np.isfinite(final_time) or final_time<0:raise ValueError('finite final time required')
    if not np.allclose(h,h.conj().T,atol=1e-12,rtol=0):raise ValueError('Hermitian H required')
    unitary=np.eye(len(h),dtype=complex);previous=0.
    for time,vertex,amplitude in sorted(kicks,key=lambda item:item[0]):
        if (not np.isfinite(time) or not 0<=time<=final_time or not np.isfinite(amplitude)
                or vertex.shape!=h.shape or not np.allclose(vertex,vertex.conj().T,atol=1e-12,rtol=0)):
            raise ValueError('valid time, amplitude and Hermitian vertex required')
        unitary=expm(-1j*amplitude*vertex)@expm(-1j*h*(time-previous))@unitary
        previous=time
    return expm(-1j*h*(final_time-previous))@unitary


def mixed_response_from_influence(h,c,a,b,observation_time,source_time,step=2e-4):
    """Return d²Re Gamma/(dJ_delta,A dJ_common,B), expected -chi_AB.

    A future common source cancels between branches by unitarity. Causality
    is tested by the ordered branch calculation, not a hard-coded theta.
    """
    final_time=max(observation_time,source_time)+.3
    values={}
    for sa in (-1,1):
        for sb in (-1,1):
            plus=branch_unitary(h,[(observation_time,a,sa*step/2),(source_time,b,sb*step)],final_time)
            minus=branch_unitary(h,[(observation_time,a,-sa*step/2),(source_time,b,sb*step)],final_time)
            action=influence(c,plus,minus)['principal_action']
            if action is None:raise ValueError('response crossed a zero-overlap history')
            values[sa,sb]=action.real
    return float((values[1,1]-values[1,-1]-values[-1,1]+values[-1,-1])/(4*step**2))


def static_susceptibility(h,vertex):
    """Ground-state static response for affine H+J V; contacts are separate."""
    e,u,_=ground_covariance(h);v=u.conj().T@vertex@u
    occupied=e<0;empty=e>0
    gaps=e[empty][None,:]-e[occupied][:,None]
    return float(-2*np.sum(abs(v[np.ix_(occupied,empty)])**2/gaps))


def fock_annihilators(modes):
    """Independent small-Fock-space oracle; never used for the large model."""
    if not isinstance(modes,int) or not 1<=modes<=5:raise ValueError('oracle limited to 1--5 modes')
    matrices=[]
    for mode in range(modes):
        a=np.zeros((2**modes,2**modes),complex)
        for state in range(2**modes):
            if state&(1<<mode):
                parity=(state&((1<<mode)-1)).bit_count()
                a[state^(1<<mode),state]=(-1)**parity
        matrices.append(a)
    return matrices


def second_quantize(matrix,annihilators):
    return sum((matrix[i,j]*annihilators[i].conj().T@annihilators[j]
                for i in range(len(matrix)) for j in range(len(matrix))),
               np.zeros_like(annihilators[0]))


def gaussian_fock_state(c,annihilators):
    n,u=eigh(c)
    if min(n)<=0 or max(n)>=1:raise ValueError('mixed-state oracle requires strict occupations')
    k=(u*np.log((1-n)/n)[None,:])@u.conj().T
    density=expm(-second_quantize(k,annihilators))
    return density/np.trace(density)
