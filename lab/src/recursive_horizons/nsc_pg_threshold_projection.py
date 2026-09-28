"""Closed-infinity-channel spectral integration using the owned reflection.

The physical subgap integrand is split into a smooth finite-horizon part and
a part linear in R(E). Cauchy's theorem moves only the analytic latter term
into the upper frequency half-plane. No average reflection, covariance fill,
new physical regulator, or metric history is introduced.
"""
from dataclasses import dataclass

import numpy as np
from numpy.polynomial.legendre import leggauss


def horizon_vectors(projected_basis):
    if np.shape(projected_basis)!=(8,2):raise ValueError('four full-spin packet basis projections required')
    u=np.zeros(8,complex);v=np.array(projected_basis[:,1],complex);w=np.zeros(8,complex)
    u[6:]=projected_basis[6:,0]
    w[:6]=projected_basis[:6,0]
    return u,v,w


def closed_fiber_split(z,kappa,projected_basis,conjugate_projected_basis):
    """Analytic coefficients G0,G1,K0,K1 with real-E value G0+R G1+h.c.

    Products use the conjugate-energy solution's adjoint. Conjugating the
    solution at z itself would make the contour integrand non-holomorphic.
    K is centered covariance: F(C_H-P/2)F^dagger. Below the parent mass gap
    P has rank two and the incoming-infinity column is absent.
    """
    z=complex(z)
    if kappa<=0 or z.real<=0 or abs(z.imag)>=kappa/2:
        raise ValueError('positive real frequency in the pole-free horizon strip required')
    u,v,w=horizon_vectors(projected_basis)
    ud,vd,wd=horizon_vectors(conjugate_projected_basis)
    outer=lambda a,b:np.outer(a,b.conj())
    exp_minus=np.exp(-2*np.pi*z/kappa)
    f=exp_minus/(1+exp_minus);s=np.exp(-np.pi*z/kappa)/(1+exp_minus);d=f-.5
    G0=outer(u,ud)+outer(v,vd)+outer(w,wd)
    G1=outer(v,ud)
    K0=d*(outer(u,ud)+outer(v,vd)-outer(w,wd))+1j*s*(outer(w,ud)-outer(u,wd))
    K1=d*outer(v,ud)-1j*s*outer(v,wd)
    return np.array([G0,G1,K0,K1])


@dataclass
class ThresholdIntegral:
    left: float
    right: float
    contour_height: float
    points_per_segment: int
    gram_positive_energy: np.ndarray
    centered_positive_energy: np.ndarray
    analytic_integral: np.ndarray
    smooth_integral: np.ndarray


def contour_integral(left,right,kappa,basis_projection,reflection,*,points=16,height=None):
    """Integrate the subgap Gram/centered-C kernels, before dE/(2pi).

    basis_projection(z) and reflection(z) must belong to the same transmitted
    Dirac domain. right is the parent mass threshold; callers must not apply
    this rank-two identity to an open infinity channel. Height is a numerical
    contour parameter inside the existing thermal analytic strip.
    """
    height=kappa/4 if height is None else float(height)
    if not 0<left<right or not 0<height<kappa/2 or points<4:
        raise ValueError('subgap interval and pole-free numerical contour required')
    nodes,weights=leggauss(points)
    def coefficients(z):
        plus=basis_projection(complex(z))
        minus=plus if np.imag(z)==0 else basis_projection(np.conj(z))
        return closed_fiber_split(z,kappa,plus,minus)
    smooth=np.zeros((2,8,8),complex)
    for x,w in zip(nodes,weights):
        z=left+(x+1)*(right-left)/2
        c=coefficients(z);smooth+=w*(right-left)/2*c[[0,2]]
    analytic=np.zeros_like(smooth)
    path=((complex(left),left+1j*height),(left+1j*height,right+1j*height),(right+1j*height,complex(right)))
    for a,b in path:
        for x,w in zip(nodes,weights):
            z=a+(x+1)*(b-a)/2
            c=coefficients(z);analytic+=w*(b-a)/2*reflection(z)*c[[1,3]]
    result=smooth+analytic+analytic.swapaxes(-1,-2).conj()
    return ThresholdIntegral(left,right,height,points,result[0],result[1],analytic,smooth)
