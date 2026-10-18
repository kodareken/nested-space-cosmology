"""Characteristic Dirac transfer and metric insertions on a supplied null chart.

The physical coefficient is W=k*sqrt(C)/r from nsc_atlas_dirac_radar.
2*dV eta_plus=-i*W*eta_minus and 2*dU eta_minus=-i*W*eta_plus
give the midpoint Cayley gate below. Positive characteristic flux is the
Euclidean norm of sqrt(du)*eta_plus, sqrt(dv)*eta_minus wires.

Admission: reuse the radial Dirac action and Cayley transform; supply the
missing causal transfer and its metric derivative on the finite wall diamond.
W and profiles are supplied, real, finite and resolved; no state is selected.
Controls cover plane waves, a C1 wall-crossing exact solution, conservation,
causality, source finite differences, sign/units and inherited CAR differences.
Control allowance: 60 process-CPU seconds, <=32 cells per direction.
First application: only N=64/128, k=1..8, nested and smooth-cap coefficients;
180 process-CPU seconds, 8 workers, BLAS=1, <=256 MiB working matrices.
Stop on a failed control, tolerance or budget. No P001 evolution is invoked.
The finite-rank preparation difference is conditional on a compatible parent;
the vacuum polarization and the gravitational/junction solve remain absent.
"""
from __future__ import annotations

import math

import numpy as np
from numba import njit, prange

from .nsc_boundary_state import GaussianBoundaryState
from .nsc_influence import _covariance


@njit(cache=True)
def _transfer(w, du, dv, profiles):
    nu, nv = w.shape
    size = nu+nv
    history = np.eye(size, dtype=np.complex128)
    vertices = np.zeros((len(profiles), size, size), dtype=np.complex128)
    weight = math.sqrt(du*dv)/4
    for i in range(nu):
        for j in range(nv):
            p, m = i, nu+j
            gamma = w[i,j]*weight
            denominator = 1+gamma*gamma
            cosine = (1-gamma*gamma)/denominator
            sine = -2j*gamma/denominator
            # Every insertion is a SINGLE angular copy. Magnetic
            # degeneracy multiplies a final expectation exactly once.
            for q in range(len(profiles)):
                factor = 2*gamma*profiles[q,i,j]/denominator
                if factor != 0:
                    for a in range(size):
                        up = history[p,a].conjugate()
                        um = history[m,a].conjugate()
                        for b in range(size):
                            vertices[q,a,b] += factor*(
                                up*history[m,b]+um*history[p,b])
            for a in range(size):
                up, um = history[p,a], history[m,a]
                history[p,a] = cosine*up+sine*um
                history[m,a] = sine*up+cosine*um
    return history, vertices


@njit(cache=True, parallel=True)
def _batch(w, du, dv, profiles):
    count, nu, nv = w.shape
    size = nu+nv
    histories = np.empty((count,size,size), dtype=np.complex128)
    vertices = np.empty((count,len(profiles),size,size), dtype=np.complex128)
    for sector in prange(count):
        h, v = _transfer(w[sector],du,dv,profiles)
        histories[sector] = h
        vertices[sector] = v
    return histories, vertices


def _real_array(value, name, ndim):
    raw = np.asarray(value)
    if raw.dtype.kind not in 'iuf' or raw.ndim != ndim:
        raise ValueError(name+' must be real with the declared dimension')
    result = np.asarray(raw, dtype=np.float64)
    if not np.all(np.isfinite(result)):
        raise ValueError(name+' must be finite')
    return np.ascontiguousarray(result)


def _inputs(potential, du, dv, profiles, ndim):
    w = _real_array(potential,'potential',ndim)
    if not all(w.shape):
        raise ValueError('a nonempty characteristic grid is required')
    for width in (du,dv):
        if isinstance(width,bool) or not np.isscalar(width) or not np.isrealobj(width):
            raise ValueError('positive finite real characteristic widths required')
        if not np.isfinite(width) or width <= 0:
            raise ValueError('positive finite real characteristic widths required')
    p = (np.empty((0,)+w.shape[-2:]) if profiles is None
         else _real_array(profiles,'metric profiles',3))
    if p.shape[1:] != w.shape[-2:]:
        raise ValueError('metric profiles must match the potential cells')
    return w,float(du),float(dv),p


def null_transfer(potential, du, dv, *, profiles=None):
    """Incoming-to-outgoing unitary and i*U^dagger*dU/dlambda.

    The metric deformation is delta ln(N)=delta ln(a)=lambda*f,
    delta r=0, with support off the null boundary: W -> W*exp(lambda*f).
    The returned vertices are integrated generators in units hbar, for
    ONE radial/angular copy, with no magnetic degeneracy or chosen state.
    Wires are ordered (eta_plus on V=0, eta_minus on U=0) on entry,
    (eta_plus on V=end, eta_minus on U=end) on exit.
    """
    w,du,dv,p = _inputs(potential,du,dv,profiles,2)
    return _transfer(w,du,dv,p)


def batch_null_transfer(potentials, du, dv, *, profiles=None):
    """Numba parallel sectors; caller selects the physical-core thread count."""
    w,du,dv,p = _inputs(potentials,du,dv,profiles,3)
    return _batch(w,du,dv,p)


def ordered_metric_response(vertices, profiles):
    """-i[J_later,J_earlier] for two strictly timelike ordered supports.

    Compact disjoint supports remove the cross contact for these two
    conformal-radial insertions. Overlap and general gauge response need
    additional contacts and the metric/observer solve.
    """
    j = np.asarray(vertices,dtype=complex)
    p = _real_array(profiles,'metric profiles',3)
    if (len(p) != 2 or j.ndim != 3 or j.shape[0] != 2 or
            j.shape[1] != j.shape[2] or not np.isfinite(j).all()):
        raise ValueError('two square finite metric vertices and matching profiles required')
    if j.shape[1] != sum(p.shape[1:]):
        raise ValueError('metric vertices and null grids must match')
    if not np.allclose(j,j.conj().transpose(0,2,1),rtol=0,atol=2e-11):
        raise ValueError('Hermitian metric insertions required')
    later,earlier = [np.argwhere(f != 0) for f in p]
    if not len(later) or not len(earlier):
        raise ValueError('both metric supports must be nonempty')
    if not np.all(np.max(earlier,axis=0) < np.min(later,axis=0)):
        raise ValueError('the entire earlier support must precede the later in U and V')
    return -1j*(j[0]@j[1]-j[1]@j[0])


def reset_inherited_mode(parent_covariance, incoming_mode, occupation):
    """Retain a supplied parent's complement, reset ONE smooth mode's occupation.

    P C_parent P + nu |f><f| is a CAR covariance for 0<=nu<=1.
    This is a declared preparation operation, not a choice of root vacuum.
    A continuum Hadamard statement additionally requires a compatible parent
    and smoothness of C_parent*f; those properties are not inferred here.
    """
    c = _covariance(parent_covariance)
    f = _normalized_mode(incoming_mode,len(c))
    nu = _occupation(occupation)
    projector = np.eye(len(c))-np.outer(f,f.conj())
    return projector@c@projector+nu*np.outer(f,f.conj())


def inherit_boundary_state(incoming, transfer, *, outgoing_parent_mask):
    """Propagate the full supplied covariance, including cross-leg correlations.

    The actual characteristic transfer supplies the map missing from a formal
    inheritance rule. The outgoing partition must be identified geometrically;
    on the radar square, future V=end is exterior and future U=end is cap.
    This determines conditional inheritance, not the initial root state.
    """
    if not isinstance(incoming,GaussianBoundaryState):
        raise ValueError('a supplied correlated Gaussian boundary state is required')
    u = np.asarray(transfer,dtype=complex)
    if u.shape != incoming.covariance.shape or not np.isfinite(u).all():
        raise ValueError('transfer must match every incoming characteristic mode')
    if not np.allclose(u.conj().T@u,np.eye(len(u)),rtol=0,atol=2e-11):
        raise ValueError('a flux-unitary characteristic transfer is required')
    return GaussianBoundaryState(u@incoming.covariance@u.conj().T,
                                 outgoing_parent_mask)


def _normalized_mode(value, size):
    f = np.asarray(value,dtype=complex)
    if f.shape != (size,) or not np.isfinite(f).all():
        raise ValueError('a finite normalized mode on all incoming wires is required')
    if abs(float(np.vdot(f,f).real)-1) > 2e-11:
        raise ValueError('normalize the incoming mode in characteristic flux')
    return f


def _occupation(value):
    if isinstance(value,bool) or not np.isscalar(value) or not np.isrealobj(value):
        raise ValueError('a real occupation in [0,1] is required')
    if not np.isfinite(value) or not 0 <= value <= 1:
        raise ValueError('a real occupation in [0,1] is required')
    return float(value)


def occupation_response(response, incoming_mode, *, angular_k, occupation=1.):
    """Change in <J_later>/d(lambda_earlier)/hbar for one full multiplet.

    The expectation difference is independent of an unchosen baseline
    parent covariance. It does not calculate that baseline's response.
    k and -k are separate sectors, each of degeneracy 2*abs(k).
    """
    r = np.asarray(response,dtype=complex)
    if r.ndim != 2 or r.shape[0] != r.shape[1] or not np.isfinite(r).all():
        raise ValueError('a square finite response is required')
    if not np.allclose(r,r.conj().T,rtol=0,atol=2e-11):
        raise ValueError('a Hermitian retarded response is required')
    if isinstance(angular_k,bool) or not isinstance(angular_k,int) or not angular_k:
        raise ValueError('a nonzero integer angular Dirac eigenvalue is required')
    f = _normalized_mode(incoming_mode,len(r))
    return float(2*abs(angular_k)*_occupation(occupation)*np.vdot(f,r@f).real)


def smooth_bump(value, center, half_width):
    """C-infinity compact bump, peak one; no bandwidth cutoff is imposed."""
    if not np.isfinite(center) or not np.isfinite(half_width) or half_width <= 0:
        raise ValueError('finite center and positive finite half width required')
    z = (np.asarray(value,dtype=float)-center)/half_width
    if not np.isfinite(z).all():
        raise ValueError('finite bump arguments required')
    out = np.zeros_like(z)
    inside = np.abs(z) < 1
    out[inside] = np.exp(1-1/(1-z[inside]**2))
    return out


def reference_profiles(duration, cells):
    """Declared exterior metric source B and later cap metric detector A."""
    if not np.isfinite(duration) or duration <= 0 or cells < 32:
        raise ValueError('positive duration and at least 32 cells required')
    labels = (np.arange(cells)+.5)*duration/cells
    u,v = np.meshgrid(labels,labels,indexing='ij')
    profiles = []
    for cu,cv in ((.80,.70),(.35,.55)):
        profiles.append(smooth_bump(u,cu*duration,.04*duration)*
                        smooth_bump(v,cv*duration,.04*duration))
    return labels,np.asarray(profiles)


def exterior_prepared_mode(duration, cells):
    """One explicitly modulated incoming exterior mode, not a vacuum spectrum."""
    labels = (np.arange(cells)+.5)*duration/cells
    incoming = smooth_bump(labels,.6*duration,.3*duration)*np.exp(-12j*labels/duration)
    incoming *= math.sqrt(duration/cells)
    incoming /= np.linalg.norm(incoming)
    return np.concatenate((np.zeros(cells,dtype=complex),incoming))
