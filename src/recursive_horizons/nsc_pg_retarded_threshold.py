"""Subgap covariance contour in terms of the bounded whole-line resolvent.

This is algebraically the same horizon state as the direct R_h split. It
avoids cancellation between O(10^9) horizon-basis Gram terms at large angular
label. No state component or local counterterm is added or removed.
"""
from dataclasses import dataclass
import numpy as np
from numpy.polynomial.legendre import leggauss

from .nsc_pg_threshold_projection import horizon_vectors


def coherence_coefficients(z,kappa,basis,dual_basis):
    u,v,w=horizon_vectors(basis);ud,vd,wd=horizon_vectors(dual_basis)
    exp_minus=np.exp(-2*np.pi*z/kappa)
    f=exp_minus/(1+exp_minus);s=np.exp(-np.pi*z/kappa)/(1+exp_minus);d=f-.5
    outer=lambda a,b:np.outer(a,b.conj())
    smooth=-2*d*outer(w,wd)+1j*s*(outer(w,ud)-outer(u,wd))
    reflected=-1j*s*outer(v,wd)
    return d,smooth,reflected


@dataclass
class RetardedThresholdIntegral:
    gram_positive_energy: np.ndarray
    centered_positive_energy: np.ndarray
    analytic_integrals: np.ndarray
    smooth_integral: np.ndarray


def retarded_threshold_integral(left,right,kappa,node,*,points=24,height=None):
    """Use G_E=(R^+-R^{+dagger})/i and the unchanged horizon coherence.

    node(z) returns projection, projection_conjugate, and, off the real axis,
    reflection and packet_resolvent. The Green part is bounded by the owned
    self-adjoint resolvent; its state coefficient remains the same f-1/2.
    """
    height=kappa/4 if height is None else height
    if not 0<left<right or not 0<height<kappa/2 or points<4:
        raise ValueError('closed-channel interval and pole-free contour required')
    x,weights=leggauss(points);smooth=np.zeros((8,8),complex)
    for t,wgt in zip(x,weights):
        z=complex(left+(t+1)*(right-left)/2);data=node(z)
        _,part,_=coherence_coefficients(z,kappa,data['projection'],data['projection_conjugate'])
        smooth+=wgt*(right-left)/2*part
    integrals=np.zeros((3,8,8),complex)
    for a,b in ((complex(left),left+1j*height),(left+1j*height,right+1j*height),(right+1j*height,complex(right))):
        for t,wgt in zip(x,weights):
            z=a+(t+1)*(b-a)/2;data=node(z)
            d,_,reflected=coherence_coefficients(z,kappa,data['projection'],data['projection_conjugate'])
            G=data['packet_resolvent'];R=data['reflection'].item()
            integrals+=wgt*(b-a)/2*np.array([G,d*G,R*reflected])
    gram=-1j*(integrals[0]-integrals[0].conj().T)
    centered=-1j*(integrals[1]-integrals[1].conj().T)+smooth+integrals[2]+integrals[2].conj().T
    return RetardedThresholdIntegral(gram,centered,integrals,smooth)
