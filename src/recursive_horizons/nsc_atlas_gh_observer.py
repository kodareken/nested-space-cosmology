"""Invariant point readouts for a specified timelike scalar-field observer.

The fields are canonical (S, chi).  Derivatives are components in the
coordinate-time / physical-reference-spatial frame used by the caller;
no metric, evolution equation or observer trajectory is supplied here.
"""
import math

import numpy as np


def _array(value,name,shape):
    raw=np.asarray(value)
    if raw.dtype.kind not in 'iuf' or raw.shape!=shape or not np.all(np.isfinite(raw)):
        raise ValueError(name+' needs finite real values with shape '+str(shape))
    return np.asarray(raw,dtype=float)


def observer_readout(g,q,p,dq,u,*,potential,omega_ref):
    """Measure total canonical-scalar density and chi along a timelike u.

    Trailing shapes are g=(4,4), q,p=(2,), dq=(3,2), u=(4,), with common
    optional leading batch dimensions.  p is partial_t q and dq_i=E_i q;
    the E_i are already in physical reference-frame derivative units.
    potential is the full W(S,chi), including its vacuum contribution once,
    and may be scalar or broadcast to the batch shape.

    The supplied future-timelike u is normalized for every readout.  In
    particular density_total=T_ab u_unit^a u_unit^b is not the ADM normal
    density unless these two observers coincide.  input_norm_error reports
    abs(g(u,u)+1) before normalization; it is not an evolution error bound.

    The auxiliary quadrature is Z=chi-i*u_unit(chi)/omega_ref.  Its wrapped
    phase lies in [-pi,pi].  At exactly Z=0 phase_wrapped=0 by convention
    and defined=False: this is an undefined phase, not a metric failure.
    No winding, internal charge or matter--antimatter count is inferred.
    """
    if isinstance(omega_ref,(bool,np.bool_)):
        raise ValueError('omega_ref must be finite and positive')
    omega=float(omega_ref)
    if not math.isfinite(omega) or omega<=0:
        raise ValueError('omega_ref must be finite and positive')
    raw=np.asarray(g)
    if raw.ndim<2 or raw.shape[-2:]!=(4,4):
        raise ValueError('g needs trailing shape (4,4)')
    batch=raw.shape[:-2]
    if raw.size==0:raise ValueError('at least one point is required')
    g=_array(g,'g',batch+(4,4));q=_array(q,'q',batch+(2,))
    p=_array(p,'coordinate rates p',batch+(2,));dq=_array(dq,'dq',batch+(3,2))
    u=_array(u,'observer velocity u',batch+(4,))
    if np.max(abs(g-np.swapaxes(g,-1,-2)))>2e-12*max(1.,float(np.max(abs(g)))):
        raise ValueError('g must be symmetric in metric indices')
    if np.min(np.linalg.eigvalsh(g[...,1:,1:]))<=0:
        raise ValueError('spatial metric must be positive definite')
    try:gi=np.linalg.inv(g)
    except np.linalg.LinAlgError as error:raise ValueError('metric inverse failed') from error
    if not np.all(np.isfinite(gi)) or np.any(gi[...,0,0]>=0):
        raise ValueError('time slices must have a positive real lapse')
    W=np.asarray(potential)
    if W.dtype.kind not in 'iuf' or not np.all(np.isfinite(W)):
        raise ValueError('potential needs finite real values')
    try:W=np.broadcast_to(np.asarray(W,dtype=float),batch)
    except ValueError as error:raise ValueError('potential must broadcast to the batch shape') from error

    norm=np.einsum('...a,...ab,...b->...',u,g,u)
    if not np.all(np.isfinite(norm)) or np.any(norm>=0) or np.any(u[...,0]<=0):
        raise ValueError('observer velocity must be future timelike')
    unit=u/np.sqrt(-norm)[...,None]
    gradients=np.concatenate((p[...,None,:],dq),axis=-2)
    proper_rates=np.einsum('...a,...aA->...A',unit,gradients)
    kinetic_contraction=np.einsum('...ab,...aA,...bA->...',gi,gradients,gradients)
    density=np.sum(proper_rates**2,axis=-1)+.5*kinetic_contraction+W
    chi_rate=proper_rates[...,1]
    imaginary=-chi_rate/omega
    amplitude=np.hypot(q[...,1],imaginary)
    defined=amplitude!=0
    phase=np.where(defined,np.arctan2(imaginary,q[...,1]),0.)
    result=dict(unit_velocity=unit,input_norm_error=abs(norm+1),density_total=density,
                chi_proper_derivative=chi_rate,quadrature_amplitude=amplitude,
                phase_wrapped=phase,defined=defined)
    if any(not np.all(np.isfinite(value)) for value in result.values()):
        raise ValueError('derived observer readouts are not finite')
    return result
