"""Value-only incoming residual: the history is evolved, its tangents are not.

Merit is the joint sup-norm of the two seam residuals on the supplied nodes:

    m(g) = max_nodes max(|E_N|, |E_beta|).

It is the acceptance norm for later search steps. It is not an error bound,
and a nodal merit is not the continuous EXISTENCE test.
"""
from dataclasses import replace
from hashlib import sha256

import numpy as np

from .nsc_compatible_history_geometry import (
    CompatibleIncomingMetric, CompatibleRadiusDirection)
from .nsc_dirac_source_phase import formal_source_phase_coefficient
from .nsc_evolved_incoming_constraints import (
    compatible_history_slots, surface_geometry_response)
from .nsc_ks_cutoff_bridge_evaluator import source_edge_from_phase
from .nsc_ks_difference_envelope import KSDifferenceIncoming
from .nsc_ks_source_envelope import _digest_arrays
from .nsc_ks_history_evaluator import fixed_target_union
from .nsc_ks_local_constraints import _family_matter
from .nsc_local_incoming_family import LocalIncomingFamily

MERIT_NAME = 'nodal joint sup-norm'
MERIT_DEFINITION = 'm(g)=max over nodes of max(|E_N|, |E_beta|)'
TARGET_COUNTS = (47, 129, 257)


def target_nodes(family, count):
    """47 is the owned solve/verification union. 129 and 257 are Chebyshev on I."""
    if count not in TARGET_COUNTS:
        raise ValueError('target count must be 47, 129 or 257')
    if not isinstance(family, LocalIncomingFamily):
        raise TypeError('owned local incoming family required')
    if count == 47:
        nodes = fixed_target_union(family.collocation_nodes(16), family.collocation_nodes(33))
    else:
        nodes = np.asarray(family.collocation_nodes(count), float)
    if len(nodes) != count:
        raise ValueError('target node count changed')
    return nodes


def merit(gradient):
    """Declared joint norm. Both components are inside the same maximum."""
    values = np.asarray(gradient, float)
    if values.ndim != 2 or values.shape[1] != 2 or not np.isfinite(values).all():
        raise ValueError('nodal (N, beta) residual required')
    return float(np.max(np.abs(values)))


def one_direction_metric(family):
    """The physical (w, U) profile as one direction, with no basis tangents."""
    if not isinstance(family, LocalIncomingFamily):
        raise TypeError('owned local incoming family required')
    w_function, u_function = family.functions
    direction = CompatibleRadiusDirection(
        w_function, u_function, family.normal_inner, family.normal_outer)
    return CompatibleIncomingMetric((1.0,), (direction,))


def assert_zero_tangents(prepared):
    if not isinstance(prepared, KSDifferenceIncoming):
        raise TypeError('value evaluation requires the difference-envelope state')
    prepared.validate()
    if int(prepared.column_tangents.shape[0]) != 0 or int(prepared.axial_tangents.shape[0]) != 0:
        raise ValueError('value evaluation carries no tangent directions')
    substitution = prepared.coverage.get('outgoing_momentum_substitution')
    if substitution not in (None, False, 'not used'):
        raise ValueError('source energies are labels; k=-E is not a column momentum')
    return True


def matter_change(prepared, channel, coefficients):
    """Contract one signed family. The tangent axis stays empty."""
    assert_zero_tangents(prepared)
    z = np.asarray(prepared.z, float)
    matter = _family_matter(prepared, channel, coefficients, 0, z)
    tangent = np.asarray(matter['tangent'])
    if tangent.shape[0] != 0:
        raise ValueError('zero-direction matter tangent was not empty')
    change = np.asarray(matter['correction'], float)
    if change.shape != (len(z), 2) or not np.isfinite(change).all():
        raise ValueError('matter change must be a finite nodal (N, beta) value')
    return change


def geometry_change(family, z, coefficients):
    slots, slot_tangents = compatible_history_slots(family.metric(), z, 0)
    if slot_tangents.shape[0] != 0:
        raise ValueError('value geometry does not carry slot tangents')
    value = surface_geometry_response(slots, slot_tangents, coefficients)['action_gradient_change']
    return np.asarray(value, float)


def edge_change(family, z, signed_angular_square_weights, axial_scale, rho_up, gauss_nodes=192):
    """Unit-angular edge of the one-direction profile. The tangent row is dropped."""
    phase = formal_source_phase_coefficient(
        one_direction_metric(family), z, angular=1.0, rho_up=float(rho_up),
        gauss_nodes=gauss_nodes)
    gradient, tangent = source_edge_from_phase(
        phase, signed_angular_square_weights, axial_scale)
    if np.asarray(tangent).shape[0] != 1:
        raise ValueError('one-direction edge tangent has the wrong count')
    return np.asarray(gradient, float)


def assemble_residual(baseline, geometry, matter, edge):
    """Baseline, geometry, matter and edge, each once."""
    base = np.asarray(baseline, float)
    parts = [np.asarray(geometry, float), np.asarray(matter, float), np.asarray(edge, float)]
    if base.shape == (2,):
        base = np.broadcast_to(base, parts[0].shape).copy()
    shapes = {base.shape, *(part.shape for part in parts)}
    if len(shapes) != 1 or base.ndim != 2 or base.shape[1] != 2:
        raise ValueError('baseline, geometry, matter and edge must share nodal (N, beta) shape')
    if any(not np.isfinite(part).all() for part in (base, *parts)):
        raise ValueError('residual parts must be finite')
    return base + parts[0] + parts[1] + parts[2]


def _with_preparation(prepared, covariance, weights):
    """Rebind covariance or weights and the digest that commits them."""
    source_digest = _digest_arrays(covariance, weights, prepared.source_energies)
    preparation = sha256(
        source_digest.encode()
        + np.ascontiguousarray(prepared.initial_columns).tobytes()
        + repr((float(prepared.mass), float(prepared.angular), float(prepared.rho_up),
                float(prepared.binding.rho_sigma))).encode()).hexdigest()
    diagnostics = dict(prepared.diagnostics)
    diagnostics['fixed_preparation_digest'] = preparation
    return replace(
        prepared, source_covariance=np.asarray(covariance),
        column_weights=np.asarray(weights, float),
        fixed_preparation_digest=preparation, diagnostics=diagnostics)


def mutated_matter(prepared, channel, coefficients, mutation):
    """Controls that must not reproduce the honest matter change."""
    if mutation == 'drop-delta-C':
        return np.zeros((len(prepared.z), 2), float)
    if mutation == 'drop-coherence':
        covariance = np.diag(np.diag(np.asarray(prepared.source_covariance, complex)))
        prepared = _with_preparation(prepared, covariance, prepared.column_weights)
    elif mutation == 'double-weights':
        prepared = _with_preparation(
            prepared, prepared.source_covariance, np.asarray(prepared.column_weights, float) * 2)
    elif mutation == 'k=-E':
        coverage = dict(prepared.coverage)
        coverage['outgoing_momentum_substitution'] = 'k=-E'
        prepared = replace(prepared, coverage=coverage)
        assert_zero_tangents(prepared)
        raise AssertionError('k=-E substitution must be rejected before contraction')
    else:
        raise ValueError('unknown value-evaluator mutation')
    return matter_change(prepared, channel, coefficients)
