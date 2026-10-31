"""Phase/amplitude gradient split of an auxiliary real-scalar quadrature.

All derivatives are raw E_i derivatives in the already normalized spatial
reference frame.  p is a coordinate-time derivative, not normal momentum;
d2q[i,j]=E_i E_j q is ordered.  No evolution equation is evaluated here.
The quadrature Z=chi-i*pi/omega_ref is not an internal or baryon charge.
"""
import math

import numpy as np


def _array(value,name,shape):
    raw=np.asarray(value)
    if raw.dtype.kind not in 'iuf' or raw.shape!=shape or not np.all(np.isfinite(raw)):
        raise ValueError(name+' needs finite real values with shape '+str(shape))
    return np.asarray(raw,dtype=float)


def wave_phase_diagnostics(g,dg,q,p,dq,d2q,dp,*,omega_ref):
    """Return nodewise reference and physical scalar gradient diagnostics.

    Trailing shapes: g=(4,4), dg=(3,4,4), q,p=(2,), dq,dp=(3,2),
    d2q=(3,3,2), with common optional leading batch dimensions.  Field
    ordering is (canonical volume S, real parent scalar chi).  The fixed
    positive omega_ref defines the auxiliary quadrature and must be stated
    when interpreting its phase.  No factor of the reference radius is
    inserted: callers supply E_i jets in physical reference-frame units.

    phase_gradient_contribution means A^2 |Dtheta|^2, not |Dtheta|^2.
    phase_current is Im(conj(Z)*E_i Z)=A^2 E_i theta.  At exact Z=0 both
    split contributions and the MASKED identity residual are set to zero;
    the full gradient remains available, and defined/undefined_count report
    where the amplitude/phase decomposition has no pointwise meaning.
    """
    if isinstance(omega_ref,(bool,np.bool_)):
        raise ValueError('omega_ref must be finite and positive')
    omega=float(omega_ref)
    if not math.isfinite(omega) or omega<=0:
        raise ValueError('omega_ref must be finite and positive')
    raw=np.asarray(g)
    if raw.ndim<2 or raw.shape[-2:]!=(4,4):raise ValueError('g needs trailing shape (4,4)')
    batch=raw.shape[:-2]
    if raw.size==0:raise ValueError('at least one point is required')
    g=_array(g,'g',batch+(4,4));dg=_array(dg,'dg',batch+(3,4,4))
    q=_array(q,'q',batch+(2,));p=_array(p,'coordinate rates p',batch+(2,))
    dq=_array(dq,'dq',batch+(3,2));dp=_array(dp,'dp',batch+(3,2))
    d2q=_array(d2q,'ordered d2q',batch+(3,3,2))
    for value,name in ((g,'g'),(dg,'dg')):
        if np.max(abs(value-np.swapaxes(value,-1,-2)))>2e-12*max(1.,float(np.max(abs(value)))):
            raise ValueError(name+' must be symmetric in metric indices')
    h=g[...,1:,1:]
    if np.min(np.linalg.eigvalsh(h))<=0:raise ValueError('spatial metric must be positive definite')
    try:gi=np.linalg.inv(g);hi=np.linalg.inv(h)
    except np.linalg.LinAlgError as error:raise ValueError('metric inverse failed') from error
    if not np.all(np.isfinite(gi)) or np.any(gi[...,0,0]>=0):
        raise ValueError('time slices must have a positive real lapse')
    N=1/np.sqrt(-gi[...,0,0])
    beta=np.einsum('...ij,...j->...i',hi,g[...,0,1:])
    dgi00=-np.einsum('...a,...iab,...b->...i',gi[...,0,:],dg,gi[..., :,0])
    dN=N[...,None]**3*dgi00/2
    shift_covariant_rate=dg[..., :,0,1:]-np.einsum('...ijk,...k->...ij',dg[..., :,1:,1:],beta)
    dbeta=np.einsum('...jl,...il->...ij',hi,shift_covariant_rate)

    chi=q[...,1];dchi=dq[..., :,1]
    pi=(p[...,1]-np.einsum('...j,...j->...',beta,dchi))/N
    dpi=(dp[..., :,1]-np.einsum('...ij,...j->...i',dbeta,dchi)
         -np.einsum('...j,...ij->...i',beta,d2q[..., :, :,1]))/N[...,None]
    dpi-=pi[...,None]*dN/N[...,None]

    imaginary=-pi/omega;dimaginary=-dpi/omega
    amplitude=np.hypot(chi,imaginary)
    defined=amplitude!=0
    # Normalize first to avoid division by A^2 or A^4 near a small nonzero node.
    real_hat=np.divide(chi,amplitude,out=np.zeros_like(amplitude),where=defined)
    imaginary_hat=np.divide(imaginary,amplitude,out=np.zeros_like(amplitude),where=defined)
    amp_projection=real_hat[...,None]*dchi+imaginary_hat[...,None]*dimaginary
    phase_projection=real_hat[...,None]*dimaginary-imaginary_hat[...,None]*dchi
    full=np.sum(dchi*dchi+dimaginary*dimaginary,axis=-1)
    amp=np.sum(amp_projection**2,axis=-1)
    phase=np.sum(phase_projection**2,axis=-1)
    physical_full=(np.einsum('...ij,...i,...j->...',hi,dchi,dchi)
                   +np.einsum('...ij,...i,...j->...',hi,dimaginary,dimaginary))
    physical_amp=np.einsum('...ij,...i,...j->...',hi,amp_projection,amp_projection)
    physical_phase=np.einsum('...ij,...i,...j->...',hi,phase_projection,phase_projection)
    result=dict(amplitude_squared=amplitude*amplitude,
                reference_gradient_squared=full,
                amplitude_gradient_squared=amp,
                phase_gradient_contribution=phase,
                phase_current=amplitude[...,None]*phase_projection,
                identity_residual=np.where(defined,full-amp-phase,0.),
                defined=defined,undefined_count=int(np.count_nonzero(~defined)),
                normal_momentum=pi,normal_momentum_gradient=dpi,
                physical_gradient_squared=physical_full,
                physical_amplitude_gradient_squared=physical_amp,
                physical_phase_gradient_contribution=physical_phase,
                physical_identity_residual=np.where(defined,physical_full-physical_amp-physical_phase,0.))
    if any(not np.all(np.isfinite(value)) for value in result.values()):
        raise ValueError('derived wave diagnostics are not finite')
    return result
