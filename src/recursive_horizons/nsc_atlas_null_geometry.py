"""Direct spacetime curvature from ordered regular-frame metric jets.

The frame is E_0=partial_t, [E_i,E_j]=2*r*epsilon_ijk E_k, and
[E_0,E_i]=0.  ``r=0`` also supports a commuting-frame control.  Derivatives
are component derivatives in this frame, not covariant derivatives.

Conventions: nabla_{E_c} E_b = Gamma[a,c,b] E_a and
Riemann[a,b,c,d] = g_{ae} (R(E_c,E_d) E_b)^e.  Thus de Sitter has
R_abcd=H**2*(g_ac*g_bd-g_ad*g_bc), and the optical tidal matrix is
R(e_A,k,e_B,k).  No Einstein equation or stress tensor is used: metric
second time derivatives must come from the actual metric history.
"""
from __future__ import annotations

import math

import numpy as np


def _array(value, name, shape):
    raw = np.asarray(value)
    if raw.dtype.kind not in 'iuf' or raw.shape != shape or not np.all(np.isfinite(raw)):
        raise ValueError(f'{name} must be a finite real array of shape {shape}')
    out = np.ascontiguousarray(raw, dtype=float)
    tolerance = 2e-12 * max(1., float(np.max(np.abs(out), initial=0.)))
    if np.max(np.abs(out - out.swapaxes(-1, -2)), initial=0.) > tolerance:
        raise ValueError(f'{name} must be symmetric in its metric indices')
    return out


def _frame_brackets(frame_scale):
    if isinstance(frame_scale, (bool, np.bool_)):
        raise ValueError('frame_scale must be finite and nonnegative')
    r = float(frame_scale)
    if not math.isfinite(r) or r < 0:
        raise ValueError('frame_scale must be finite and nonnegative')
    C = np.zeros((4, 4, 4))  # C[c,d,e] = C^e_cd
    for i, j, k in ((1, 2, 3), (2, 3, 1), (3, 1, 2)):
        C[i, j, k] = 2*r
        C[j, i, k] = -2*r
    return C


def _first_jets(g, v, dg, frame_scale):
    raw = np.asarray(g)
    if raw.ndim < 2 or raw.shape[-2:] != (4, 4):
        raise ValueError('g needs trailing shape (4,4)')
    batch = raw.shape[:-2]
    g = _array(g, 'g', batch + (4, 4))
    v = _array(v, 'v', batch + (4, 4))
    dg = _array(dg, 'dg', batch + (3, 4, 4))
    count = int(np.prod(batch)) if batch else 1
    if not count:
        raise ValueError('the metric batch must not be empty')
    g = g.reshape(count, 4, 4)
    if np.min(np.linalg.eigvalsh(g[:, 1:, 1:])) <= 0:
        raise ValueError('the spatial metric must be positive definite')
    try:
        inverse = np.linalg.inv(g)
    except np.linalg.LinAlgError as exc:
        raise ValueError('the spacetime metric must be invertible') from exc
    if np.any(inverse[:, 0, 0] >= 0):
        raise ValueError('the metric must have Lorentzian signature and spacelike time slices')
    D = np.empty((count, 4, 4, 4))
    D[:, 0] = v.reshape(count, 4, 4)
    D[:, 1:] = dg.reshape(count, 3, 4, 4)
    return batch, g, inverse, D, _frame_brackets(frame_scale)


def _koszul(g, D, C):
    count = len(g)
    lower = np.empty((count, 4, 4, 4))
    for a in range(4):
        for c in range(4):
            for b in range(4):
                value = D[:, c, b, a] + D[:, b, c, a] - D[:, a, c, b]
                for e in range(4):
                    value = (value - g[:, c, e]*C[b, a, e]
                             + g[:, b, e]*C[a, c, e] + g[:, a, e]*C[c, b, e])
                lower[:, a, c, b] = value/2
    return lower


def first_connection(g, v, dg, *, frame_scale):
    """Return inverse metric and full connection from the first metric jets.

    Inputs have optional common batch axes followed by g,v=(4,4),
    dg=(3,4,4).  Gamma's final indices are (raised, derivative, vector).
    """
    batch, g, inverse, D, C = _first_jets(g, v, dg, frame_scale)
    Gamma = np.einsum('nae,necb->nacb', inverse, _koszul(g, D, C))
    return {'inverse': inverse.reshape(batch + (4, 4)),
            'Gamma': Gamma.reshape(batch + (4, 4, 4))}


def geometry_from_jets(g, v, dg, d2g, dv, accg, *, frame_scale):
    """Return direct metric connection, covariant Riemann and Ricci tensors.

    Trailing shapes: g,v,accg=(4,4), dg,dv=(3,4,4),
    d2g=(3,3,4,4).  Here v=partial_t g, dg_i=E_i g,
    d2g_ij=E_i(E_j g), dv_i=E_i(partial_t g), accg=partial_t**2 g.
    All inputs share g's optional leading batch axes.  The spatial second
    jets are ordered and MUST NOT be symmetrized.  ``commutator_error``
    reports the maximum absolute jet identity residual per point; no
    repair or curvature symmetrization masks an inconsistent supplied jet.
    """
    batch, g, inverse, D, C = _first_jets(g, v, dg, frame_scale)
    count = len(g)
    d2g = _array(d2g, 'ordered d2g', batch + (3, 3, 4, 4)).reshape(count, 3, 3, 4, 4)
    dv = _array(dv, 'dv', batch + (3, 4, 4)).reshape(count, 3, 4, 4)
    accg = _array(accg, 'accg', batch + (4, 4)).reshape(count, 4, 4)
    DD = np.empty((count, 4, 4, 4, 4))
    DD[:, 0, 0] = accg
    DD[:, 0, 1:] = dv
    DD[:, 1:, 0] = dv
    DD[:, 1:, 1:] = d2g
    commutator = DD - DD.swapaxes(1, 2) - np.einsum('cde,neab->ncdab', C, D)
    commutator_error = np.max(np.abs(commutator), axis=(1, 2, 3, 4))

    lower = _koszul(g, D, C)
    dlower = np.empty((count, 4, 4, 4, 4))
    for a in range(4):
        for c in range(4):
            for b in range(4):
                value = DD[:, :, c, b, a] + DD[:, :, b, c, a] - DD[:, :, a, c, b]
                for e in range(4):
                    value = (value - D[:, :, c, e]*C[b, a, e]
                             + D[:, :, b, e]*C[a, c, e] + D[:, :, a, e]*C[c, b, e])
                dlower[:, :, a, c, b] = value/2
    dinverse = -np.einsum('nap,nmpq,nqe->nmae', inverse, D, inverse)
    Gamma = np.einsum('nae,necb->nacb', inverse, lower)
    dGamma = (np.einsum('nmae,necb->nmacb', dinverse, lower)
              + np.einsum('nae,nmecb->nmacb', inverse, dlower))
    raised = np.empty((count, 4, 4, 4, 4))
    for c in range(4):
        for d in range(4):
            raised[:, :, :, c, d] = (
                dGamma[:, c, :, d, :] - dGamma[:, d, :, c, :]
                + Gamma[:, :, c, :] @ Gamma[:, :, d, :]
                - Gamma[:, :, d, :] @ Gamma[:, :, c, :]
                - np.einsum('e,naeb->nab', C[c, d], Gamma))
    Riemann = np.einsum('nae,nebcd->nabcd', g, raised)
    Ricci = np.einsum('nabad->nbd', raised)
    return {'inverse': inverse.reshape(batch + (4, 4)),
            'Gamma': Gamma.reshape(batch + (4, 4, 4)),
            'Riemann': Riemann.reshape(batch + (4, 4, 4, 4)),
            'Ricci': Ricci.reshape(batch + (4, 4)),
            'commutator_error': commutator_error.reshape(batch)}


def optical_tidal(Riemann, screen, k):
    """Contract R(e_A,k,e_B,k), with broadcastable optional batch axes.

    ``screen`` ends in (number_of_screen_vectors,4) and k ends in (4,).
    These are frame components; callers select and normalize the null
    tangent and its spacelike orthonormal screen using the actual metric.
    No normalization or trace subtraction is applied here.
    """
    Riemann, screen, k = map(np.asarray, (Riemann, screen, k))
    if Riemann.ndim < 4 or Riemann.shape[-4:] != (4, 4, 4, 4):
        raise ValueError('Riemann needs trailing shape (4,4,4,4)')
    if screen.ndim < 2 or screen.shape[-1] != 4 or k.ndim < 1 or k.shape[-1] != 4:
        raise ValueError('screen and k need trailing shapes (screen_count,4) and (4,)')
    for array in (Riemann, screen, k):
        if array.dtype.kind not in 'iuf' or not np.all(np.isfinite(array)):
            raise ValueError('curvature, screen and tangent must be finite real arrays')
    return np.einsum('...abcd,...Aa,...b,...Bc,...d->...AB', Riemann, screen, k, screen, k)
