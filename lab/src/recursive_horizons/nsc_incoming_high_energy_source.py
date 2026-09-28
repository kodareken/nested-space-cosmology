"""Stable incoming subtraction for the archived massive middle-mode recipe.

Only arithmetic changes: the order-16 Riccati columns, their source law and
the existing fourth-order subtraction remain unchanged.  The vacuum Bloch
z difference is rationalized before subtracting the higher reference orders.
No clipping, fitted occupation, extra term, or new spectral cutoff is used.
This is an unintegrated source observable, not a certified infinite tail.
"""
import numpy as np

from .nsc_common_ks_trace import restrict_resolved_modes
from .nsc_incoming_state_moments import incoming_state_moments, RADIUS, AXIAL, SIGMA
from .nsc_paired_horizon_preparation import source_covariance as horizon_source_covariance
from .nsc_pg_archived_high_energy_modes import archived_middle_boundary_modes
from .nsc_pg_high_energy import riccati_coefficients


def incoming_high_energy_source(energies, metadata, source_covariance, *, repo_root=None):
    """Return stable rho/p_parallel/T01/p_perp kernels in the same convention.

    ``metadata`` authenticates the original middle-panel producer; frequencies
    must remain inside that panel. ``source_covariance`` is required to be
    exactly its unchanged horizon/incoming law. The first archived source
    column is zero, so horizon coherence is retained in that law but projects
    to zero in this PARTICULAR vacuum-mode approximation, not in the theory.
    """
    modes = archived_middle_boundary_modes(energies, metadata, repo_root=repo_root)
    E = modes.energies
    mass = float(metadata['channel']['compact_mass'])
    angular = metadata['angular_sign']*float(metadata['channel']['angular_eigenvalue'])
    config = metadata['config']
    expected = np.array([horizon_source_covariance(e, config['surface_gravity'], config['omega'], mass) for e in E])
    source = np.asarray(source_covariance, complex)
    if source.shape != expected.shape or not np.array_equal(source, expected):
        raise ValueError('the exact unchanged horizon/incoming source covariance is required')
    projector = np.array([np.diag([1., 1., float(e > mass)]) for e in E])
    direct = incoming_state_moments(E, modes.mode_at_one, source, projector, mass, angular)
    coefficients, _ = riccati_coefficients([1.], mass, angular, 16)
    powers = (1/(2*E[:, None]))**np.arange(1, 17)
    series = powers@coefficients[0]
    u = AXIAL**2*abs(series)**2
    denominator = 1+u
    bloch = np.stack((2*AXIAL*series.imag/denominator,
                      -2*AXIAL*series.real/denominator,
                      1-2*u/denominator), axis=1)
    x = AXIAL**2*(mass**2+angular**2/RADIUS**2)/E**2
    root = np.sqrt(1+x)
    relativistic_deficit = x/(root*(1+root))
    momentum = -E/AXIAL
    energy = E/AXIAL*root
    zeroth = np.stack((np.full(len(E), mass)/energy,
                       np.full(len(E), -angular/RADIUS)/energy,
                       1/root), axis=1)
    difference0 = np.stack((bloch[:, 0]-zeroth[:, 0], bloch[:, 1]-zeroth[:, 1],
                            relativistic_deficit-2*u/denominator), axis=1)
    orders = direct['reference_bloch_orders']
    higher = orders[1:].sum(axis=0)
    vacuum_difference = difference0-higher
    # Read f directly, rather than subtracting the rounded number (1-f).
    outgoing = source[:, 0, 0].real
    incoming = source[:, 2, 2].real
    thermal_bloch = -(outgoing+incoming)[:, None]*bloch
    difference = vacuum_difference+thermal_bloch
    trace_difference = incoming-outgoing
    h = np.stack((-np.full(len(E), mass), np.full(len(E), angular/RADIUS), momentum), axis=1)

    def moments(delta_bloch, delta_trace):
        return np.stack((np.sum(h*delta_bloch, axis=1), momentum*delta_bloch[:, 2],
                         -momentum*delta_trace, angular/(2*RADIUS)*delta_bloch[:, 1]), axis=1)

    kernels = moments(difference, trace_difference)
    vacuum = moments(vacuum_difference, np.zeros(len(E)))
    thermal = moments(thermal_bloch, trace_difference)
    canonical, _ = restrict_resolved_modes(1., E, modes.mode_at_one)
    partner = canonical[:, :, 1]
    outer = partner[:, :, None]*partner[:, None, :].conj()
    field_bloch = np.einsum('nij,vji->nv', outer, SIGMA).real
    # Account for subtraction operand sizes, without pretending that this
    # covers the mode approximation, omitted reflection, or UV integration.
    operands = abs(bloch)+abs(zeroth)+np.sum(abs(orders[1:]), axis=0)+abs(thermal_bloch)
    operands[:, 2] = (abs(relativistic_deficit)+abs(2*u/denominator)
                      +np.sum(abs(orders[1:, :, 2]), axis=0)+abs(thermal_bloch[:, 2]))
    delta_roundoff = 16*np.finfo(float).eps*operands
    roundoff = np.stack((np.sum(abs(h)*delta_roundoff, axis=1),
                         abs(momentum)*delta_roundoff[:, 2],
                         16*np.finfo(float).eps*abs(momentum)*(outgoing+incoming),
                         abs(angular)/(2*RADIUS)*delta_roundoff[:, 1]), axis=1)
    return {
        'energies': E.copy(), 'kernels': kernels,
        'kernel_order': direct['kernel_order'],
        'vacuum_kernels': vacuum, 'thermal_kernels': thermal,
        'vacuum_bloch': bloch, 'vacuum_reference_bloch_difference': vacuum_difference,
        'state_reference_trace_difference': trace_difference,
        'source_covariance': source.copy(),
        'source_horizon_coherence': source[:, 0, 1].copy(),
        'projected_horizon_coherence': np.zeros((len(E), 2, 2), complex),
        'roundoff_indicator': roundoff,
        'roundoff_dominates_kernel': abs(kernels) <= roundoff,
        'direct_kernels': direct['kernels'],
        'diagnostics': {
            'source_covariance_change': 0.,
            'canonical_partner_bloch_residual': np.max(abs(field_bloch-bloch), axis=1),
            'zeroth_reference_residual': np.max(abs(orders[0]-zeroth), axis=1),
            'direct_kernel_difference': np.max(abs(kernels-direct['kernels']), axis=1),
            'decomposition_residual': np.max(abs(kernels-vacuum-thermal), axis=1),
            'producer_current_residual': modes.current_coisometry_residual,
            'producer_riccati_residual': np.max(modes.riccati_equation_residual, axis=1),
            'producer_order12_16_field_difference': modes.order12_16_field_difference,
        },
        'provenance': modes.provenance,
        'scope': 'same archived middle-mode approximation; unintegrated source with stable vacuum subtraction',
        'coherence_scope': 'source C_H coherence retained; zero first producer column makes its projection zero',
        'error_scope': 'roundoff indicator only; no thermal/reflection approximation or infinite stress-tail bound',
        'reference_is_gaussian_state': False, 'new_cutoff_or_filter': False,
    }
