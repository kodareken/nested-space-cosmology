"""Canonical finite-band matter fluctuations on a supplied evolving metric.

The canonical fields are (S, chi), with S=sqrt(84)*(s-s_star).  This module
projects the quadratic *Hamiltonian*, including its shift term, in a real
reference-orthonormal spatial basis.  The canonical momentum coefficients
are those of the momentum density relative to the reference measure.
They are not coefficients of a normal time derivative, and the projected
kinetic matrix is not the inverse of a truncated Lagrangian mass matrix.

These are Gaussian matter-probe equations on a supplied classical metric.
No metric fluctuations, quantum state, entropy law or Hadamard continuum
limit are selected by this finite projection.
"""
import math

import numpy as np


KINETIC=84.0


def _real_array(value,name,shape):
    raw=np.asarray(value)
    if raw.dtype.kind not in 'iuf' or raw.shape!=shape or not np.all(np.isfinite(raw)):
        raise ValueError(name+' needs finite real values with shape '+str(shape))
    return np.ascontiguousarray(raw,dtype=float)


def potential_hessian(background,*,sstar):
    """Return W_AB for the existing fixed two-field action.

    ``background[...,0]`` is canonical S and ``background[...,1]`` is chi.
    The returned trailing shape is (2,2), with no vacuum subtraction:
    differentiating the full W removes a constant automatically.
    """
    raw=np.asarray(background)
    if raw.ndim<1 or raw.shape[-1]!=2:
        raise ValueError('background needs trailing shape (2,)')
    background=_real_array(background,'canonical background',raw.shape)
    if isinstance(sstar,(bool,np.bool_)):
        raise ValueError('sstar must be a finite real number')
    try:sstar=float(sstar)
    except (TypeError,ValueError) as error:
        raise ValueError('sstar must be a finite real number') from error
    if not math.isfinite(sstar):raise ValueError('sstar must be a finite real number')
    s=sstar+background[...,0]/math.sqrt(KINETIC)
    chi=background[...,1]
    with np.errstate(over='ignore',invalid='ignore'):
        e12=np.exp(-12*s);e14=np.exp(-14*s);e24=np.exp(-24*s)
        mass2=400*e12
        Vss=1008*e12-1470*e14+450*e24
        result=np.empty(background.shape[:-1]+(2,2))
        result[...,0,0]=(Vss+72*mass2*chi*chi)/KINETIC
        result[...,0,1]=-12*mass2*chi/math.sqrt(KINETIC)
        result[...,1,0]=result[...,0,1]
        result[...,1,1]=mass2
    if not np.all(np.isfinite(result)):
        raise ValueError('potential Hessian is not finite')
    return result


def _adm_nodes(g):
    """Extract ADM data; reject invalid slices instead of repairing them."""
    scale=max(1.,float(np.max(abs(g))))
    if np.max(abs(g-np.swapaxes(g,-1,-2)))>2e-12*scale:
        raise ValueError('g must be symmetric')
    h=g[:,1:,1:]
    a,b,c=h[:,0,0],h[:,0,1],h[:,0,2]
    d,e,f=h[:,1,1],h[:,1,2],h[:,2,2]
    with np.errstate(over='ignore',invalid='ignore',divide='ignore'):
        cof=np.empty_like(h)
        cof[:,0,0]=d*f-e*e;cof[:,0,1]=c*e-b*f;cof[:,0,2]=b*e-c*d
        cof[:,1,0]=cof[:,0,1];cof[:,1,1]=a*f-c*c;cof[:,1,2]=b*c-a*e
        cof[:,2,0]=cof[:,0,2];cof[:,2,1]=cof[:,1,2];cof[:,2,2]=a*d-b*b
        det=a*cof[:,0,0]+b*cof[:,0,1]+c*cof[:,0,2]
        valid=(a>0)&(cof[:,2,2]>0)&(det>0)&np.isfinite(det)
        if not np.all(valid):raise ValueError('spatial metric must be positive definite')
        inverse=cof/det[:,None,None]
        shift=np.einsum('nij,nj->ni',inverse,g[:,0,1:])
        lapse2=-g[:,0,0]+np.einsum('ni,ni->n',g[:,0,1:],shift)
        if not np.all(np.isfinite(lapse2)) or np.any(lapse2<=0):
            raise ValueError('metric requires a positive finite lapse')
        volume=np.sqrt(det);lapse=np.sqrt(lapse2)
    if not np.all(np.isfinite(inverse)) or not np.all(np.isfinite(shift)):
        raise ValueError('derived metric data are not finite')
    return lapse,shift,volume,inverse


def quadratic_hamiltonian(g,background,Y,DY,weights,fields,*,sstar):
    """Project H2=.5*q^T*K*q+p^T*B*q+.5*p^T*A*p.

    Shapes are g=(nodes,4,4), background=(nodes,2), Y=(nodes,modes),
    DY=(nodes,3,modes), weights=(nodes,), fields=(modes,).  Each retained
    mode has field 0 (S) or 1 (chi); a spatial basis column may therefore
    appear twice, once for each field.  DY already contains the physical
    inverse-reference-radius factor and the supplied weights are unit-S3
    weights.  Physical harmonics r^(3/2)*Y cancel the common r^-3 measure.

    For z=(q,p), G=[[K,B^T],[B,A]] and zdot=J*G*z.  In particular
    qdot=A*p+B*q: the sign and transpose of the shift are fixed here.
    A is the direct projection of N/sqrt(det h), not inv(project(sqrt(h)/N)).
    The quadrature and finite spatial band are supplied by the caller.
    """
    raw=np.asarray(g)
    if raw.ndim!=3 or raw.shape[1:]!=(4,4) or raw.shape[0]==0:
        raise ValueError('g needs shape (nodes,4,4) with nonzero nodes')
    nodes=len(raw);g=_real_array(g,'g',(nodes,4,4))
    background=_real_array(background,'canonical background',(nodes,2))
    raw=np.asarray(Y)
    if raw.ndim!=2 or raw.shape[0]!=nodes or raw.shape[1]==0:
        raise ValueError('Y needs shape (nodes,modes) with nonzero modes')
    modes=raw.shape[1];Y=_real_array(Y,'Y',(nodes,modes))
    DY=_real_array(DY,'physical basis derivatives',(nodes,3,modes))
    weights=_real_array(weights,'quadrature weights',(nodes,))
    if np.any(weights<0) or not np.sum(weights)>0:
        raise ValueError('quadrature weights must be nonnegative with positive sum')
    raw=np.asarray(fields)
    if raw.shape!=(modes,) or raw.dtype.kind not in 'iu' or np.any((raw!=0)&(raw!=1)):
        raise ValueError('fields must assign integer 0 or 1 to each mode')
    fields=np.asarray(raw,dtype=int)
    lapse,shift,volume,inverse=_adm_nodes(g)
    mass=potential_hessian(background,sstar=sstar)
    A=np.zeros((modes,modes));B=np.zeros_like(A);K=np.zeros_like(A)
    kinetic_weight=weights*lapse/volume
    stiffness_weight=weights*lapse*volume
    indices=[np.flatnonzero(fields==field) for field in (0,1)]
    for field,idx in enumerate(indices):
        if not len(idx):continue
        y=Y[:,idx];dy=DY[:,:,idx]
        A[np.ix_(idx,idx)]=y.T@(kinetic_weight[:,None]*y)
        shift_derivative=np.einsum('ni,nij->nj',shift,dy)
        B[np.ix_(idx,idx)]=y.T@(weights[:,None]*shift_derivative)
        # Contract the three spatial directions pointwise, then use one
        # stacked Gram product.  Nine separate mode-space products would
        # repeat the dominant nodes*modes^2 work for every (i,j) pair.
        contracted=np.einsum('nij,njm->nim',inverse,dy)
        weighted=stiffness_weight[:,None,None]*contracted
        block=dy.reshape(3*nodes,len(idx)).T@weighted.reshape(3*nodes,len(idx))
        block+=y.T@((stiffness_weight*mass[:,field,field])[:,None]*y)
        K[np.ix_(idx,idx)]=block
    if len(indices[0]) and len(indices[1]):
        cross=Y[:,indices[0]].T@((stiffness_weight*mass[:,0,1])[:,None]*Y[:,indices[1]])
        K[np.ix_(indices[0],indices[1])]=cross
        K[np.ix_(indices[1],indices[0])]=cross.T
    # Only remove roundoff antisymmetry of the integrated quadratic form.
    # No metric, mass eigenvalue or Hamiltonian mode is clipped or repaired.
    asymmetry=max(float(np.max(abs(A-A.T))),float(np.max(abs(K-K.T))))
    A=(A+A.T)/2;K=(K+K.T)/2
    G=np.block([[K,B.T],[B,A]])
    if not np.all(np.isfinite(G)):
        raise ValueError('projected quadratic Hamiltonian is not finite')
    mass_min=.5*(mass[:,0,0]+mass[:,1,1]
        -np.hypot(mass[:,0,0]-mass[:,1,1],2*mass[:,0,1]))
    return dict(G=G,A=A,K=K,B=B,lapse=lapse,shift=shift,
        volume_weight=volume,mass_hessian=mass,
        minimum_nodal_mass_eigenvalue=float(np.min(mass_min)),
        integration_antisymmetry_before_roundoff_average=asymmetry,
        phase_dimension=2*modes,quadrature_nodes=nodes)
