"""Same-action source-cutoff edge on the evolved local constraint evaluator.

The raw Dirac source sum is retained. This adds the derived ordered-kernel
edge and its complete history tangent once per represented signed family.
The uncomputed finite-source tail and quadrature errors remain OPEN.
"""
from collections.abc import Mapping

import numpy as np

from .nsc_dirac_source_phase import formal_source_phase_coefficient
from .nsc_ks_history_evaluator import KSHistoryEvaluator, subset_node_indices, _snapshot
from .nsc_ks_profile_identity import profile_identity


def source_edge_from_phase(phase, signed_angular_square_weights, axial_scale):
    """Contract a unit-angular geometric coefficient with the original ledger.

    This mathematical factorization uses f_s proportional to ell^2; it does
    not substitute ell=1 into the physical Dirac source. There is no extra
    factor for energy signs or spinor dimension.
    """
    weights = np.asarray(signed_angular_square_weights, float)
    a = float(axial_scale)
    if (weights.shape != (2,) or not np.isfinite(weights).all() or np.any(weights < 0)
            or not np.isfinite(a) or a <= 0):
        raise ValueError('original nonnegative (+E,-E) angular-square weights and positive a required')
    if phase.signs != (1, -1) or phase.angular != 1.:
        raise ValueError('unit-angular phase coefficient in (+E,-E) order required')
    fz = np.asarray(phase.f_z)
    dfz = np.asarray(phase.delta_f_z)
    if fz.shape != (2, len(phase.z)) or dfz.ndim != 3 or dfz.shape[::2] != fz.shape:
        raise ValueError('phase derivatives must retain sign, direction and target axes')
    plus, minus = weights[:, None]*fz
    dplus, dminus = weights[:, None, None]*dfz
    gradient = np.column_stack(((plus-minus)/a, -(plus+minus)))/(2*np.pi)
    tangent = np.stack(((dplus-dminus)/a, -(dplus+dminus)), axis=-1)/(2*np.pi)
    return gradient, tangent


def signed_family_angular_square_weights(family_records):
    """Use each signed family once, regardless of its source-panel batch count."""
    if not isinstance(family_records, Mapping) or not family_records:
        raise ValueError('bound signed family records required')
    weights = np.zeros(2)
    seen = set()
    for record in family_records.values():
        key = (int(record['group']), int(record['angular_sign']), int(record['energy_sign']))
        if key in seen or key[2] not in (1, -1):
            raise ValueError('unique explicitly signed source family required')
        seen.add(key)
        ell = float(record['signed_angular'])
        multiplicity = float(record['multiplicity'])
        if not np.isfinite([ell, multiplicity]).all() or multiplicity <= 0:
            raise ValueError('finite original angular eigenvalue and positive multiplicity required')
        weights[0 if key[2] == 1 else 1] += multiplicity*ell*ell
    return weights


class KSCutoffBridgeEvaluator:
    """Wrap the actual source-fixed evaluator; never accept a matter dictionary."""

    def __init__(self, raw_evaluator, *, phase_gauss_nodes=32):
        if not isinstance(raw_evaluator, KSHistoryEvaluator):
            raise TypeError('owned KSHistoryEvaluator required')
        if (isinstance(phase_gauss_nodes, (bool, np.bool_))
                or not isinstance(phase_gauss_nodes, (int, np.integer)) or phase_gauss_nodes < 1):
            raise ValueError('positive integer phase quadrature count required')
        self._raw = raw_evaluator
        self._phase_nodes = int(phase_gauss_nodes)
        self._phase_cache = {}

    @property
    def source_identity(self):
        return self._raw.source_identity

    @property
    def target_nodes(self):
        return self._raw.target_nodes

    def evaluate(self, family, z=None):
        requested = self.target_nodes if z is None else np.asarray(z, float)
        raw = self._raw.evaluate(family, requested)
        identity = profile_identity(family)
        if raw['history_identity'] != identity:
            raise ValueError('raw state and phase must belong to the same analytic history')
        phase = self._phase_cache.get(identity)
        if phase is None:
            phase = formal_source_phase_coefficient(family, self.target_nodes,
                angular=1., rho_up=self._raw._rho_up, gauss_nodes=self._phase_nodes)
            self._phase_cache[identity] = phase
        if (phase.profile_identity != identity or phase.rho_up != self._raw._rho_up
                or not np.array_equal(phase.z, self.target_nodes)):
            raise ValueError('phase/preparation restriction binding mismatch')
        weights = signed_family_angular_square_weights(raw['family_records'])
        edge, tangent = source_edge_from_phase(phase, weights, self._raw._coefficients['a'])
        indices = subset_node_indices(self.target_nodes, requested)
        edge, tangent = edge[indices], tangent[:, indices]
        if tangent.shape != raw['history_jacobian'].shape:
            raise ValueError('full retarded phase derivative must match the state direction order')
        result = dict(raw)
        budget = dict(raw['error_budget'])
        # A bound for the old raw summation is not a bound for this new
        # ordered-limit assembly, including its phase quadrature and tail.
        budget['changed_history_tail'] = None
        result.update({
            'raw_action_gradient': _snapshot(raw['action_gradient']),
            'raw_history_jacobian': _snapshot(raw['history_jacobian']),
            'source_cutoff_edge_gradient': _snapshot(edge),
            'source_cutoff_edge_history_tangent': _snapshot(tangent),
            'action_gradient': _snapshot(raw['action_gradient']+edge),
            'history_jacobian': _snapshot(raw['history_jacobian']+tangent),
            'source_cutoff_phase_binding': phase.binding,
            'source_cutoff_angular_square_weights': _snapshot(weights),
            'source_cutoff_finite_tail_error_bound': None,
            'source_cutoff_phase_quadrature_error_bound': None,
            'error_budget': _snapshot(budget),
        })
        scope = dict(raw['scope'])
        scope.update({
            'adapter': 'KSCutoffBridgeEvaluator',
            'source_cutoff_coincidence_bridge': 'derived analytic identity; numerical phase and finite tail OPEN',
            'source_cutoff_edge_included_once_per_signed_family': True,
            'source_cutoff_edge_full_retarded_derivative': True,
            'source_cutoff_phase_quadrature_error_bound': None,
            'source_cutoff_finite_tail_error_bound': None,
            'physical_constraint_status': 'OPEN',
            'physical_EXISTENCE_certificate': False,
            'new_action_term': False,
            'Gamma_rest_assigned': False,
            'local_reference_coefficients_changed': False,
        })
        scope['unresolved'] = tuple(raw['scope']['unresolved'])+(
            'numerical source-cutoff phase quadrature', 'finite-source-cutoff raw tail')
        result['scope'] = _snapshot(scope)
        return result
