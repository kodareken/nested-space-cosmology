"""Controlled repeated-throat spatial chain, with links fixed by its stencil.

This finite/periodic intervention repeats rho in [-R,R] with w=kappa/r(rho).
It is a spatial motif, not a solved nesting of Lorentzian cosmological rooms.
R and Omega are declared inputs. No independent constant mass Phi is inserted.
"""
from __future__ import annotations
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import eigh
from scipy.optimize import brentq


def geometric_motif(radius: float, intervals: int, kappa: float = 1.):
    if not np.isfinite(radius) or radius <= 0 or intervals < 4 or not np.isfinite(kappa) or kappa < 0:
        raise ValueError("positive radius and resolution, nonnegative finite kappa required")
    h = 2*radius/intervals
    edges = -radius + (np.arange(intervals)+.5)*h
    w = kappa/np.sqrt(1+edges**2)
    if np.max(h*w) >= 2:
        raise ValueError("refine grid so h*w<2; retain the first-order transfer branch")
    a, b = -1/h + w/2, 1/h + w/2
    forward = np.diag(a) + np.diag(b[:-1],1)
    local = np.block([[np.zeros_like(forward),forward.T], [forward,np.zeros_like(forward)]])
    link = np.zeros_like(local)
    link[-1,0] = b[-1]
    return {"H":local, "B":link, "h":h,
            "zero_energy_transfer": float(np.prod(-a/b)),
            "continuum_zero_transfer": float(np.exp(-2*kappa*np.arcsinh(radius)))}


def bloch_gap(radius, intervals, phases=33):
    motif = geometric_motif(radius, intervals)
    minimum = np.inf
    location = None
    for phase in np.linspace(0,np.pi,phases):
        matrix = motif['H'] + np.exp(1j*phase)*motif['B'] + np.exp(-1j*phase)*motif['B'].T
        gap = float(np.min(np.abs(eigh(matrix,eigvals_only=True))))
        if gap < minimum:
            minimum, location = gap, float(phase)
    return {"intervals":intervals,"spacing":motif['h'],"phase_samples":phases,
            "minimum_sampled_bloch_gap":minimum,"minimum_phase":location,
            "zero_transfer":motif['zero_energy_transfer']}


def continuum_discriminant(energy, radius, kappa=1.):
    def rhs(rho, state):
        w=kappa/np.sqrt(1+rho*rho)
        generator=np.array([[-w,energy],[-energy,w]])
        return (generator@state.reshape(2,2)).ravel()
    solution=solve_ivp(rhs,(-radius,radius),np.eye(2).ravel(),method='DOP853',
                       rtol=2e-11,atol=2e-12,max_step=.15)
    if not solution.success: raise RuntimeError(solution.message)
    transfer=solution.y[:,-1].reshape(2,2)
    return float(np.trace(transfer)),float(np.linalg.det(transfer))


def continuum_first_band_edge(radius):
    # Zero energy has discriminant 2*cosh(2*asinh(R))>2. Find the first
    # sampled crossing of +2 and refine it; lattice convergence is independent.
    previous=0.
    for energy in np.linspace(0,1,129)[1:]:
        trace,_=continuum_discriminant(energy,radius)
        if trace <= 2:
            root=brentq(lambda e: continuum_discriminant(e,radius)[0]-2,previous,energy,xtol=1e-12)
            final,det=continuum_discriminant(root,radius)
            return {"band_edge":root,"discriminant_residual":final-2,"monodromy_determinant":det,
                    "bracket":[float(previous),float(energy)]}
        previous=float(energy)
    raise RuntimeError("first band edge was not bracketed; no gap value inferred")
