"""Unintegrated state-minus-reference moments on the SAME rho=1 KS slice.

Global PG source columns are restricted by the existing common-Cauchy owner.
The homogeneous fourth-order Bloch subtraction is reused at q=3*pi/4;
its conformal-clock recursion agrees with the proper-KS reference symbol.
It is a formal subtraction, never a Gaussian state or occupation fit.

The stored degeneracy already includes both angular signs when lambda!=0.
For the actual sign inventory, negative frequencies are imported by the owned
opposite-angular antiunitary identity.  Summing positive-E kernels over that
inventory therefore has prefactor copies*deg/(4*pi^2*r^2*a*n_signs).
No spectral integration, local restoration, full stress or constraints are
supplied here. The separately owned massless LLL allocation is excluded.
"""
import numpy as np

from .nsc_common_ks_trace import incoming_KS_state
from .nsc_compact_ctp_neck import _adiabatic_bloch_function
from .nsc_transmitting_dirac_domain import S1, S2, S3


RHO = 1.
RADIUS = np.sqrt(2.)
AXIAL = np.sqrt(3*np.pi/2-4.)
CONFORMAL_Q = 3*np.pi/4
SIGMA = np.array([S1, S2, S3])


def incoming_adiabatic_reference(energies, mass, signed_angular):
    """Existing superadiabatic P through order four, vectorized over E."""
    raw = np.asarray(energies)
    if raw.ndim != 1 or not len(raw) or np.iscomplexobj(raw):
        raise ValueError('positive real source frequencies required')
    E = np.asarray(raw, float)
    if (not np.isfinite(E).all() or np.any(E <= 0) or mass < 0
            or not np.isfinite([mass, signed_angular]).all()):
        raise ValueError('finite positive source frequencies and existing channel labels required')
    if mass == 0 and signed_angular == 0:
        raise ValueError('massless LLL remains in its existing analytic owner')
    # Same generated expressions as adiabatic_bloch; no new recursion or
    # clock prescription. Frequency dependence broadcasts in its five terms.
    values = _adiabatic_bloch_function()(CONFORMAL_Q, mass, signed_angular, E)
    terms = np.asarray(values, float).reshape(5, 3, len(E))
    bloch = terms.sum(axis=0).T
    reference = (np.eye(2)+np.einsum('ni,ijk->njk', bloch, SIGMA))/2
    return reference, terms.transpose(0, 2, 1)


def incoming_group_factor(channel, angular_signs):
    """Factor multiplying sum_sign integral_0^infty dE of these kernels.

    Before signed folding: copies*(deg/n_signs)/(4*pi*r^2)
    times dk/(2*pi*a).  The imported negative-frequency equality supplies
    the factor two. Group degeneracy is not assigned independently per sign.
    """
    signs = tuple(angular_signs)
    angular = float(channel['angular_eigenvalue'])
    expected = {1} if angular == 0 else {-1, 1}
    if len(signs) != len(expected) or set(signs) != expected:
        raise ValueError('the complete actual angular-sign inventory is required once')
    copies, degeneracy = channel['copy_count'], channel['degeneracy']
    if (isinstance(copies, bool) or isinstance(degeneracy, bool)
            or int(copies) != copies or int(degeneracy) != degeneracy
            or min(copies, degeneracy) <= 0):
        raise ValueError('positive integer multiplicities from the retained channel ledger required')
    return float(copies*degeneracy/(4*np.pi**2*RADIUS**2*AXIAL*len(signs)))


def incoming_state_moments(energies, pg_fields, source_covariance, source_projectors,
                           mass, signed_angular, *, covariance_tolerance=3e-8):
    """Per-E rho, parallel pressure, covariant T01 and sphere pressure kernels.

    Inputs must be resolved global source columns at rho=1 and the unchanged
    C_H/incoming source fibers. No seed covariance, probe-covariance inverse,
    new normal jets, spectral cutoff or small-value filter is accepted.
    """
    if covariance_tolerance <= 0 or not np.isfinite(covariance_tolerance):
        raise ValueError('positive finite covariance tolerance required')
    reference, reference_bloch_orders = incoming_adiabatic_reference(energies, mass, signed_angular)
    E = np.asarray(energies, float)
    state = incoming_KS_state(E, pg_fields, source_covariance, source_projectors, rho=RHO)
    Csource, Psource = np.asarray(source_covariance, complex), np.asarray(source_projectors, complex)
    support = float(np.max(abs(Csource@Psource-Csource)))
    complement_lower = max(0., -float(np.linalg.eigvalsh(Psource-Csource).min()))
    diagnostics = dict(state['residuals'], source_support=support,
                       source_complement_CAR_lower=complement_lower)
    if max(diagnostics.values()) > covariance_tolerance:
        raise ValueError('resolved canonical source fails the declared covariance/current tolerance')
    C = state['covariance']
    difference = C-reference
    momentum = -E/AXIAL
    H = -mass*S1+signed_angular/RADIUS*S2+momentum[:, None, None]*S3
    vertices = np.stack((H, momentum[:, None, None]*S3,
                         -momentum[:, None, None]*np.eye(2),
                         np.broadcast_to(signed_angular/(2*RADIUS)*S2, H.shape)), axis=1)
    complex_kernels = np.einsum('nij,nvji->nv', difference, vertices)
    kernels = complex_kernels.real
    diagnostics.update(
        reference_trace=float(np.max(abs(np.trace(reference, axis1=1, axis2=2)-1))),
        reference_Hermiticity=float(np.max(abs(reference-reference.swapaxes(-1, -2).conj()))),
        kernel_imaginary=float(np.max(abs(complex_kernels.imag))))
    # Arithmetic cancellation indicator only. It neither clips the kernels
    # nor bounds upstream modal error or a missing differentiated UV tail.
    scale = np.einsum('nij,nvji->nv', abs(C)+abs(reference), abs(vertices)).real
    roundoff = 16*np.finfo(float).eps*scale
    mode_map = state['mode_map']
    gram = mode_map@mode_map.swapaxes(-1, -2).conj()
    per_energy_diagnostics = {
        'source_coisometry': np.max(abs(gram-np.eye(2)), axis=(1, 2)),
        'covariance_eigenvalues': state['eigenvalues'],
        'reference_trace': abs(np.trace(reference, axis1=1, axis2=2)-1),
        'reference_eigenvalues_not_a_state': np.linalg.eigvalsh(reference),
        'kernel_imaginary': np.max(abs(complex_kernels.imag), axis=1),
    }
    return {
        'energies': E.copy(), 'canonical_momenta': -E,
        'covariance': C, 'reference': reference,
        'reference_bloch_orders': reference_bloch_orders,
        'difference': difference, 'hamiltonian_normal': H,
        'kernels': kernels, 'kernel_order': ('rho', 'p_parallel', 'T01', 'p_perp'),
        'roundoff_indicator': roundoff,
        'roundoff_dominates_kernel': abs(kernels) <= roundoff,
        'source_complement_correlations': state['source_complement_correlations'],
        'diagnostics': diagnostics,
        'per_energy_diagnostics': per_energy_diagnostics,
        'intrinsic_geometry': {'rho': RHO, 'radius': RADIUS, 'axial': AXIAL,
                               'lapse': 1., 'shift': 0., 'frame': 'existing KS normal'},
        'reference_scope': 'formal fourth-order subtraction, not a Gaussian covariance',
        'scope': 'unintegrated incoming state moments; local allocation and full spectral source not evaluated',
        'roundoff_scope': 'arithmetic indicator, not a bound on modal, spectral, or renormalized-stress error',
    }
