"""Retained PG covariance from physical source-mode integrals.

The packet matrices are compressions of the whole transmitted field.  They
are never used as a closed finite-system initial state for the CTP determinant.
"""
import numpy as np

from .nsc_pg_packet_modes import SIGNED_PACKET_MAP
from .nsc_pg_state_covariance import equal_time_ctp_blocks


def positive_panel(projection, weights, source=None, projector=None):
    """Integrate positive-energy Gram and centered covariance with dE/(2pi)."""
    F = np.asarray(projection, complex)
    w = np.asarray(weights, float).ravel() / (2 * np.pi)
    if F.shape != (len(w), 8, 3) or not np.isfinite(F).all():
        raise ValueError('evaluated eight-packet, three-source mode columns required')
    if not np.isfinite(w).all() or np.any(w <= 0):
        raise ValueError('positive finite real-energy quadrature weights required')
    if (source is None) != (projector is None):
        raise ValueError('source and open-fiber projector must accompany one another')
    if source is None:
        P = np.broadcast_to(np.eye(3), (len(w), 3, 3))
        C = np.broadcast_to(np.diag([0., 1., 0.]), P.shape)
    else:
        C, P = np.asarray(source, complex), np.asarray(projector, complex)
        if C.shape != (len(w), 3, 3) or P.shape != C.shape:
            raise ValueError('physical source-fiber covariance required; seed slots are excluded')
        if not np.isfinite(C).all() or not np.isfinite(P).all():
            raise ValueError('finite source covariance and projector required')
        if max(np.max(abs(C-C.swapaxes(-1,-2).conj())),
               np.max(abs(P-P.swapaxes(-1,-2).conj())),
               np.max(abs(P@P-P))) > 3e-11:
            raise ValueError('source must be Hermitian on its canonical open fiber')
        if min(np.linalg.eigvalsh(C).min(), np.linalg.eigvalsh(P-C).min()) < -3e-11:
            raise ValueError('source covariance violates its CAR interval')
    gram = np.einsum('n,nai,nij,nbj->ab', w, F, P, F.conj())
    centered = np.einsum('n,nai,nij,nbj->ab', w, F, C-.5*P, F.conj())
    return gram, centered


def signed_covariance(positive, opposite):
    """Combine actual opposite-angular positive-E maps using the locked identity.

    Angular-zero groups pass the same physical family twice; their two signed
    energy halves remain distinct.  No covariance completion is performed.
    """
    G, K = (np.asarray(x, complex) for x in positive)
    Go, Ko = (np.asarray(x, complex) for x in opposite)
    if any(x.shape != (8, 8) or not np.isfinite(x).all() for x in (G,K,Go,Ko)):
        raise ValueError('both computed positive-energy families required')
    S = SIGNED_PACKET_MAP
    gram = G + S @ Go.conj() @ S
    centered = K - S @ Ko.conj() @ S
    covariance = .5*gram + centered
    return {'covariance': covariance, 'car_gram': gram,
            **equal_time_ctp_blocks(covariance, gram)}


def covariance_residuals(result):
    C, G = result['covariance'], result['car_gram']
    eigenvalues = np.linalg.eigvalsh(C)
    return {
        'CAR': float(np.linalg.norm(G-np.eye(len(G)), 2)),
        'Hermiticity': float(np.linalg.norm(C-C.conj().T, 2)),
        'eigenvalue_min': float(eigenvalues.min()),
        'eigenvalue_max': float(eigenvalues.max()),
        'probe_bulk_correlation': float(np.linalg.norm(C[:4,4:], 2)),
        'CTP_CAR': float(np.linalg.norm(1j*(result['greater']-result['lesser'])-G, 2)),
        'CTP_Keldysh': float(np.linalg.norm(result['keldysh']-result['greater']-result['lesser'], 2)),
    }

