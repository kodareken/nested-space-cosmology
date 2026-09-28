"""Common interior KS Cauchy restriction and PG/KS Gaussian trace pairing.

Inputs are already resolved GLOBAL source-mode fields and their horizon/
infinity covariance. A same-event spin frame is used only after restriction
of those solutions; it is not declared to be a Cauchy propagator or U_L0.
"""
import numpy as np

from .nsc_ks_spacetime_variation import chart_coordinates, CompactKSHarmonic
from .nsc_transmitting_dirac_domain import TransmittingDiracSeamDomain, S1, S2, S3
from .nsc_lorentzian import geometry


def restrict_resolved_modes(rho, energies, pg_fields):
    """Global fields -> canonical KS plane-wave amplitudes, k=-E.

    Phi_PG(E,rho) exp(-iE tau) = F_frame(rho) Psi_KS(E,rho) exp(-iE z).
    The flat half-density and the actual KS normal supply the frame factor.
    No stored seed covariance or freely chosen transport interval is accepted.
    """
    if np.iscomplexobj(energies) and np.any(np.imag(energies) != 0):
        raise ValueError('real physical frequencies required; analytic contour columns are not covariance states')
    E, phi = np.asarray(energies, float), np.asarray(pg_fields, complex)
    if E.ndim != 1 or phi.shape != (len(E), 2, 3) or not np.isfinite(phi).all() or not np.isfinite(E).all():
        raise ValueError('resolved three-source mode columns, not a two-moment/seed covariance, required')
    beta = float(geometry(rho)[0]); radius = np.sqrt(1+rho*rho)
    domain = TransmittingDiracSeamDomain(1., 1., beta, radius,
                                       surface_id=f'reference-KS-Cauchy-section-rho={rho:g}')
    _, inverse, _ = domain.trace_map(np.ones(1))
    # trace_inverse acts on psi; supplied PG fields are r*psi at q_PG=1.
    canonical = np.einsum('ab,ebs->eas', inverse[0]/radius, phi)
    T, S = chart_coordinates(float(rho))
    canonical *= np.exp(1j*E*S)[:, None, None]
    return canonical, {'KS_time': T, 'clock_shift': S, 'axial_scale': domain.induced_axial_scale}


def incoming_KS_state(energies, pg_fields, source_covariance, source_projectors, *, rho=1.):
    """Physical Gaussian restriction on the upstream interior slice.

    The supplied response metric is identically unperturbed for rho>=1;
    both characteristics point downstream. This slice is a Cauchy surface
    for the future interior collar, not for the entire exterior spacetime.
    The output is a Fourier-symbol sample; it does not close the continuum
    into a few-mode dynamics.
    """
    if rho < 1 or float(geometry(rho)[0]) <= 1:
        raise ValueError('unperturbed upstream slice inside the trapped collar required')
    F, chart = restrict_resolved_modes(rho, energies, pg_fields)
    C, P = np.asarray(source_covariance, complex), np.asarray(source_projectors, complex)
    if C.shape != (len(F), 3, 3) or P.shape != C.shape or not np.isfinite(C).all() or not np.isfinite(P).all():
        raise ValueError('the same resolved horizon/infinity source fibers are required')
    adj = lambda a: a.swapaxes(-1, -2).conj()
    if max(np.max(abs(C-adj(C))), np.max(abs(P-adj(P))), np.max(abs(P@P-P))) > 3e-11:
        raise ValueError('Hermitian source covariance and orthogonal source projector required')
    source_eig = np.linalg.eigvalsh(C)
    if source_eig.min() < -3e-11 or source_eig.max() > 1+3e-11:
        raise ValueError('source covariance must satisfy the inherited CAR interval')
    initial = F@C@adj(F)
    gram = F@adj(F)
    complement = P-adj(F)@F
    correlations = F@C@complement
    eig = np.linalg.eigvalsh(initial)
    residuals = {
        'source_coisometry': float(np.max(abs(gram-np.eye(2)))),
        'closed_source_projection': float(np.max(abs(F@P-F))),
        'source_complement_projector': float(np.max(abs(complement@complement-complement))),
        'Hermiticity': float(np.max(abs(initial-adj(initial)))),
        'CAR_lower': max(0., -float(eig.min())),
        'CAR_upper': max(0., float(eig.max())-1.),
    }
    return {'canonical_momenta': -np.asarray(energies), 'mode_map': F, 'covariance': initial,
            'source_complement_correlations': correlations, 'source_complement': complement,
            'eigenvalues': eig, 'chart': chart, 'residuals': residuals,
            'measure': 'dk/(2*pi)=dE/(2*pi), with reversed integration orientation',
            'physical_scope': 'restriction to the upstream KS Cauchy slice and its future interior domain'}


def reference_KS_vertices(data, *, omega=.4):
    """Whole Gaussian weak kernel in the common reference KS coordinates.

    This independently reconstructs the archived PG metric-variation kernel.
    It is a reference-geometry trace identity, not a new stationary history.
    The midpoint canonical momentum is required by the symmetric derivative.
    """
    E = np.asarray(data['energies'], float)
    repeated = np.repeat(E, 3); kmid = -.5*(repeated[:, None]+repeated[None, :])
    difference = repeated[:, None]-repeated[None, :]
    mass, angular = float(data['mass'].item()), float(data['angular'].item())
    result = np.zeros((4, len(repeated), len(repeated)), complex)
    direction = CompactKSHarmonic(omega)
    for rho, weight, field in zip(data['fields/rho'], data['fields/weights'], data['fields/fields']):
        F, chart = restrict_resolved_modes(float(rho), E, field)
        F = F.transpose(1, 0, 2).reshape(2, -1)
        radius, axial = np.sqrt(1+rho*rho), chart['axial_scale']
        density = F.conj().T@F
        s1, s2, s3 = (F.conj().T@s@F for s in (S1, S2, S3))
        vertex = np.array([-mass*s1+angular/radius*s2+kmid/axial*s3,
                           -kmid*density, -kmid/axial**2*s3, -angular/radius**2*s2])
        phase = direction.data(float(rho))[0]*np.exp(1j*difference*chart['clock_shift'])
        result += weight/axial*phase*vertex
    return result
