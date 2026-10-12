"""Coherent matter differences without subtracting two background contractions."""
import numpy as np

from .nsc_evolved_incoming_constraints import (
    source_column_matter, raw_ks_vertex_coefficients,
)
from .nsc_ks_difference_envelope import KSDifferenceIncoming


def source_fixed_matter_difference(prepared, *, axial_scale, radius, multiplicity,
                                  reference='joint'):
    """Contract the exact cross and quadratic terms of F=F_ref+delta F.

    F_ref is the homogeneous solution from this object's A_up, evaluated
    jointly or by the assembly's preparation-keyed independent cache. Its
    variation is zero. No free incoming covariance or stress is accepted.
    """
    if not isinstance(prepared, KSDifferenceIncoming):
        raise TypeError('joint reference/difference prepared state required')
    prepared.validate()
    # The frozen object constructor binds its read-only difference arrays.
    if reference == 'joint':
        Aref = prepared.reference_amplitudes
    elif reference == 'cached':
        from .nsc_ks_local_constraints import homogeneous_reference_amplitudes
        Aref = homogeneous_reference_amplitudes(prepared)
    else:
        raise ValueError('reference must be joint or cached')
    # Retain the difference between both numerical integrations. Dropping it
    # would silently remove drift and break equality with the raw assembly.
    difference = prepared.envelope_difference + (prepared.reference_amplitudes - Aref)[None]
    parameters = dict(mass=prepared.mass, angular=prepared.angular,
                      axial_scale=axial_scale, radius=radius, multiplicity=multiplicity)
    current = source_column_matter(
        prepared.weighted_columns, prepared.weighted_axial_columns,
        prepared.source_covariance, prepared.weighted_column_tangents,
        prepared.weighted_axial_tangents, **parameters)
    phase_weight = np.exp(-1j * prepared.source_energies * prepared.z[:, None])[:, None, :] * prepared.column_weights
    ref = phase_weight * Aref[None]
    refz = phase_weight * (-1j * prepared.source_energies * Aref)[None]
    delta = phase_weight * difference
    deltaz = phase_weight * (prepared.envelope_difference_z - 1j * prepared.source_energies * difference)
    C = prepared.source_covariance

    def adjoint(value):
        return value.swapaxes(-1, -2).conj()

    cross = delta @ C @ adjoint(ref)
    density = cross + adjoint(cross) + delta @ C @ adjoint(delta)
    bilinear = deltaz @ C @ adjoint(ref) + refz @ C @ adjoint(delta) + deltaz @ C @ adjoint(delta)
    momentum = (bilinear - adjoint(bilinear)) / (2j)
    vertices = raw_ks_vertex_coefficients(
        np.array([1., 0., axial_scale, radius]), prepared.mass, prepared.angular, envelopes=np.ones(4))
    gradient = -multiplicity * (
        np.einsum('bij,...ji->...b', vertices.multiplication[:2], density)
        + np.einsum('bij,...ji->...b', vertices.momentum[:2], momentum))
    imag = float(np.max(abs(gradient.imag)))
    if imag > 3e-11:
        raise ArithmeticError('non-real matter difference exceeds algebra tolerance')
    return {
        'action_gradient_change': gradient.real,
        'action_gradient_tangent': current['action_gradient_tangent'],
        'current_action_gradient': current['action_gradient'],
        'density_change': density, 'momentum_density_change': momentum,
        'imaginary_residual': imag,
        'reference_source': prepared.fixed_preparation_digest,
        'reference_evaluation': reference,
        'reference_integration_difference_retained': float(np.max(abs(prepared.reference_amplitudes - Aref))),
        'reference_tangent': 'zero by fixed upstream preparation',
        'cross_and_quadratic_terms_included': True,
        'stress_drift_subtracted': False,
        'physical_constraint_status': 'OPEN', 'spectral_and_field_error_bound': None,
    }
