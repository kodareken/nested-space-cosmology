"""Gaussian covariance blocks from the physical PG source-mode map.

Probe coordinates are observables of the full field. The Q complement is
retained as a spatial kernel and is never replaced by a closed probe system.
"""
import numpy as np


def spectral_probe_bulk_blocks(source_covariance,source_projector,probe_modes,field_modes,probe_embedding,*,other_field_modes=None,other_probe_embedding=None):
    """Return C_JJ(E), C_JQ(E;rho), C_QQ(E;rho,rho) before dE/(2pi).

    field_modes is Phi_E(rho), probe_modes=J^dagger Phi_E and
    probe_embedding=J(rho). Each input comes from the same normalized mode
    owner; this is not a map from a stored seed covariance to spatial slots.
    """
    C=np.asarray(source_covariance,complex);P=np.asarray(source_projector,complex)
    F=np.asarray(probe_modes,complex);phi=np.asarray(field_modes,complex);J=np.asarray(probe_embedding,complex)
    if C.shape!=P.shape or C.ndim!=2 or C.shape[0]!=C.shape[1] or F.shape[1]!=C.shape[0] or phi.shape[1]!=C.shape[0] or J.shape!=(phi.shape[0],F.shape[0]):
        raise ValueError('source, field and probe maps have incompatible domains')
    if not all(np.isfinite(a).all() for a in (C,P,F,phi,J)):raise ValueError('finite evaluated mode maps required')
    if np.linalg.norm(P-P.conj().T)>3e-11 or np.linalg.norm(P@P-P)>3e-11:
        raise ValueError('source support must be its canonical open-fiber projector')
    if np.linalg.norm(C-C.conj().T)>3e-11:raise ValueError('source covariance is not Hermitian')
    if min(np.linalg.eigvalsh(C).min(),np.linalg.eigvalsh(P-C).min()) < -3e-11:
        raise ValueError('source covariance violates CAR in its declared fiber')
    bulk=phi-J@F
    result={'probe':F@C@F.conj().T,'probe_bulk':F@C@bulk.conj().T,
            'bulk_diagonal':bulk@C@bulk.conj().T,'bulk_mode_field':bulk,
            'field_reconstruction_residual':float(np.linalg.norm(phi-(J@F+bulk)))}
    if other_field_modes is not None or other_probe_embedding is not None:
        if other_field_modes is None or other_probe_embedding is None:raise ValueError('both maps at the other spatial point are required')
        other=np.asarray(other_field_modes)-np.asarray(other_probe_embedding)@F
        result['bulk_pair']=bulk@C@other.conj().T
    return result


def equal_time_ctp_blocks(covariance,car_gram):
    """Keep the computed CAR Gram in the numerical lesser/greater relation."""
    C=np.asarray(covariance,complex);G=np.asarray(car_gram,complex)
    if C.shape!=G.shape or C.ndim!=2 or C.shape[0]!=C.shape[1]:raise ValueError('one common probe algebra required')
    return {'lesser':1j*C,'greater':-1j*(G-C),'keldysh':-1j*(G-2*C)}
