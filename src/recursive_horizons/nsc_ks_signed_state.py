"""Opposite-angular negative-energy map of a source-fixed KS difference state.

The inhomogeneous envelope law on a real radius history is

    X_rho = (S3 / a^2) X_z + (i / a) (-m S1 + ell / r S2 - E / a S3) X.

It is mapped by (E, ell, X) -> (-E, -ell, S3 conjugate(X)). The partner is
the actual opposite-sign family, not a factor-two folding of one evolution.
After inhomogeneous evolution, F(E, z) F(E, z)^dagger is not assumed to be
I; the negative covariance uses the transformed field and the supplied
C_negative, never I - C_positive(z).
"""
from __future__ import annotations

from types import MappingProxyType

import numpy as np

from .nsc_evolved_incoming_state import FixedSourcePreparation
from .nsc_ks_difference_envelope import KSDifferenceIncoming
from .nsc_ks_source_envelope import _fixed_preparation_digest
from .nsc_lorentzian import geometry
from .nsc_transmitting_dirac_domain import S1, S2, S3


SOURCE_COMPLEMENT_TOLERANCE = 3e-13


def s3_conjugate(field):
    """Apply S3 K on the two-component spinor in the second-to-last axis."""
    value = np.asarray(field, complex)
    if value.ndim < 2 or value.shape[-2] != 2:
        raise ValueError(
            'spinor axis of length 2 required in the second-to-last position')
    if not np.isfinite(value).all():
        raise ValueError('finite spinor field required')
    return S3 @ np.conjugate(value)


def source_complement_residual(positive_covariance, negative_covariance):
    """max |C_neg - (I - conjugate(C_pos))| of the matching source law."""
    positive = np.asarray(positive_covariance, complex)
    negative = np.asarray(negative_covariance, complex)
    if (positive.ndim != 2 or positive.shape != negative.shape
            or positive.shape[0] != positive.shape[1] or positive.shape[0] < 1):
        raise ValueError('matching square source covariances required')
    if not (np.isfinite(positive).all() and np.isfinite(negative).all()):
        raise ValueError('finite source covariances required')
    complement = np.eye(positive.shape[0], dtype=complex) - positive.conj()
    return float(np.max(np.abs(negative - complement)))


def actual_radius_envelope_rhs(X, Xz, rho, energies, mass, angular, radius):
    """Full actual-radius KS envelope PDE, including the z derivative.

    X and Xz have shape (nz, 2, nsrc). radius is a positive scalar or (nz,).
    """
    field = np.asarray(X, complex)
    axial = np.asarray(Xz, complex)
    if field.ndim != 3 or field.shape[1] != 2 or axial.shape != field.shape:
        raise ValueError('matching (nz, 2, nsrc) envelope and z derivative required')
    if not (np.isfinite(field).all() and np.isfinite(axial).all()):
        raise ValueError('finite envelope and z derivative required')
    rho = float(rho)
    mass = float(mass)
    angular = float(angular)
    raw_energies = np.asarray(energies)
    raw_radius = np.asarray(radius)
    if np.iscomplexobj(raw_energies) and np.any(raw_energies.imag != 0):
        raise ValueError('finite real energies required')
    if np.iscomplexobj(raw_radius) and np.any(raw_radius.imag != 0):
        raise ValueError('real radius required')
    labels = np.asarray(raw_energies.real, float)
    if (labels.ndim != 1 or labels.shape != (field.shape[2],)
            or not np.isfinite(labels).all() or not np.isfinite([rho, mass, angular]).all()):
        raise ValueError('finite real energies, mass, angular and rho required')
    radius = np.asarray(raw_radius.real, float)
    if radius.ndim == 0:
        radius = np.full(field.shape[0], float(radius))
    elif radius.shape != (field.shape[0],):
        raise ValueError('radius must be a positive scalar or shape (nz,)')
    if not np.isfinite(radius).all() or np.any(radius <= 0):
        raise ValueError('positive finite actual radius required')
    beta = float(geometry(rho)[0])
    axial_scale = np.sqrt(beta * beta - 1.)
    if not np.isfinite(axial_scale) or axial_scale <= 0:
        raise ArithmeticError('strictly trapped a>0 required')
    term_z = np.einsum('ab,zbs->zas', S3, axial) / (axial_scale * axial_scale)
    term_m = -mass * np.einsum('ab,zbs->zas', S1, field)
    term_ell = (angular / radius)[:, None, None] * np.einsum('ab,zbs->zas', S2, field)
    term_E = -(labels[None, None, :] / axial_scale) * np.einsum('ab,zbs->zas', S3, field)
    return term_z + (1j / axial_scale) * (term_m + term_ell + term_E)


def _carrier_phase(energies, z):
    return np.exp(-1j * np.asarray(energies) * np.asarray(z)[:, None])[:, None, :]


def _reconstruct_carrier(reference, difference, difference_z, energies, z):
    phase = _carrier_phase(energies, z)
    total = np.asarray(reference)[None] + np.asarray(difference)
    columns = phase * total
    axial = phase * (np.asarray(difference_z) - 1j * np.asarray(energies) * total)
    return columns, axial, phase


def negative_angular_partner(prepared, negative_source):
    """Map a positive-energy KSDifferenceIncoming to the opposite-angular family.

    ``negative_source`` is the caller's original FixedSourcePreparation. Its
    covariance is stored as supplied. The I-conjugate(C_src) relation is a
    diagnostic of the matching source law; it does not replace C_negative or
    set envelope differences to zero.
    """
    if not isinstance(prepared, KSDifferenceIncoming):
        raise TypeError(
            'KSDifferenceIncoming required; the signed map owns the joint '
            'reference/difference incoming state')
    if not isinstance(negative_source, FixedSourcePreparation):
        raise TypeError(
            'original FixedSourcePreparation required for the negative source; '
            'a covariance matrix or frozen incoming C0 is not accepted')
    prepared.validate()
    energies = np.asarray(prepared.source_energies)
    if np.any(energies <= 0):
        raise ValueError('positive input energies required for the negative partner')
    if not np.array_equal(negative_source.energies, -energies):
        raise ValueError('negative-source energies must equal -positive energies')
    if not np.array_equal(negative_source.column_weights, prepared.column_weights):
        raise ValueError('negative-source column weights must equal the positive weights')
    residual = source_complement_residual(
        prepared.source_covariance, negative_source.covariance)
    if residual > SOURCE_COMPLEMENT_TOLERANCE:
        raise ValueError('negative source differs from the fixed horizon source law')
    angular = -float(prepared.angular)
    if angular == 0.0:
        angular = 0.0
    initial = s3_conjugate(prepared.initial_columns)
    reference = s3_conjugate(prepared.reference_amplitudes)
    difference = s3_conjugate(prepared.envelope_difference)
    difference_z = s3_conjugate(prepared.envelope_difference_z)
    columns, axial, _ = _reconstruct_carrier(
        reference, difference, difference_z, negative_source.energies, prepared.z)
    # S3 conjugation commutes with real z/history derivatives. Mapping the
    # complete tangents directly avoids cancelling and restoring large E*Y
    # carrier terms just to recover the same axial derivative.
    column_tangents = s3_conjugate(prepared.column_tangents)
    axial_tangents = s3_conjugate(prepared.axial_tangents)
    digest = _fixed_preparation_digest(
        negative_source, initial, prepared.mass, angular, prepared.rho_up,
        prepared.binding.rho_sigma)
    diagnostics = dict(prepared.diagnostics)
    diagnostics.update({
        'fixed_preparation_digest': digest,
        'positive_preparation_digest': prepared.fixed_preparation_digest,
        'signed_map': 'S3 conjugation with opposite angular sign',
        'energy_folding_factor': 1,
        'supplied_negative_covariance_preserved': True,
        'covariance_replaced_by_rounded_complement': False,
        'field_isometry_assumed': False,
        'source_complement_residual': residual,
        'boundary_state_signed_relation': residual,
        'same_real_history': True,
        'physical_source_error': 'OPEN',
        'physical_gate': 'OPEN',
    })
    coverage = dict(prepared.coverage)
    coverage['angular_multiplicity_added'] = False
    coverage['quadrature_applied_once'] = True
    coverage['physical_source_error'] = 'OPEN'
    return KSDifferenceIncoming(
        prepared.z, columns, column_tangents, axial, axial_tangents,
        negative_source.covariance, negative_source.column_weights,
        negative_source.energies, prepared.mass, angular, initial, prepared.rho_up,
        prepared.binding, digest, MappingProxyType(coverage), diagnostics,
        reference, difference, difference_z)
