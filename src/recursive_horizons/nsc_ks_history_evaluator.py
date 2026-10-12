"""Source-fixed production evaluator for local incoming Newton on I=S(1)+[.12,.18].

Each candidate LocalIncomingFamily is evolved from frozen KSUpstreamBatch
objects through the owned energy propagator and batched KS constraints.
C_up, C_src, weights, mass, ell and rho_up stay fixed. The physical gate
remains OPEN: interpolation, field, source-cutoff and coincidence errors
are not supplied here.
"""
from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import replace
from types import MappingProxyType

import numpy as np

from .nsc_evolved_incoming_state import FixedSourcePreparation
from .nsc_ks_batched_constraints import KSConstraintAccumulator, KSUpstreamBatch
from .nsc_ks_difference_envelope import KSDifferenceIncoming
from .nsc_ks_energy_propagator import evolve_energy_propagator
from .nsc_ks_profile_identity import profile_identity
from .nsc_ks_signed_state import (
    SOURCE_COMPLEMENT_TOLERANCE, negative_angular_partner, s3_conjugate,
    source_complement_residual,
)
from .nsc_ks_source_envelope import (
    RHO_SIGMA, _finite_real, physical_incoming_interval, usual_axial_support,
)
from .nsc_local_history_newton import (
    FIXED_AXIAL_INNER, FIXED_AXIAL_OUTER, FIXED_NORMAL_INNER, FIXED_NORMAL_OUTER,
    PRODUCTION_LOCAL_GATE,
)
from .nsc_local_incoming_constraints import (
    _channel_record, local_error_budget,
)
from .nsc_local_incoming_family import LocalIncomingFamily


_FOREIGN = frozenset({
    'history', 'HistoryBinding', 'node_fields', 'incoming_C0', 'C0',
    'covariance', 'matter', 'optimizer', 'Gamma_rest', 'metric_step',
    'tangents',
})
_NODE_VECTORS = (
    'action_gradient', 'baseline_action_gradient',
    'local_reference_gradient_change', 'incoming_slots',
)
_NODE_TANGENTS = (
    'history_jacobian', 'local_reference_gradient_tangent',
    'incoming_slot_tangents',
)
_NODE_MAPS = (
    'batch_corrections', 'batch_current_matter', 'batch_reference_matter',
    'family_corrections',
)
_TANGENT_MAPS = (
    'batch_matter_tangents', 'family_matter_tangents',
)


def _reject_foreign(kwargs):
    for name in kwargs:
        if name in _FOREIGN or name.lower() in {key.lower() for key in _FOREIGN}:
            raise ValueError(
                f'{name} is not an evaluator argument; frozen C0, Gamma_rest '
                'and metric steps are rejected')
    if kwargs:
        raise TypeError(f'unsupported evaluator argument {next(iter(kwargs))}')


def _copy_array(value):
    array = np.array(value, copy=True)
    array.setflags(write=False)
    return array


def _snapshot(value):
    """Own nested numerical inputs and receipts independently of caller buffers."""
    if isinstance(value, np.ndarray):
        return _copy_array(value)
    if isinstance(value, Mapping):
        return MappingProxyType({key: _snapshot(item) for key, item in value.items()})
    if isinstance(value, (list, tuple)):
        return tuple(_snapshot(item) for item in value)
    return value


def _require_local_family(family):
    if not isinstance(family, LocalIncomingFamily):
        raise TypeError(
            'LocalIncomingFamily required; frozen C0 and matter dictionaries '
            'are not this evaluator')
    expected = physical_incoming_interval()
    if not np.allclose(family.interval, expected, atol=3e-12, rtol=0):
        raise ValueError('incoming interval is the fixed I=S(1)+[.12,.18]')
    if not np.allclose(
            [family.axial_inner, family.axial_outer, family.normal_inner, family.normal_outer],
            [FIXED_AXIAL_INNER, FIXED_AXIAL_OUTER, FIXED_NORMAL_INNER, FIXED_NORMAL_OUTER],
            atol=3e-15, rtol=0):
        raise ValueError('axial cutoff and normal support are fixed solver-domain data')
    if family.radius_lower_bound() <= 0:
        raise ValueError('family radius_lower_bound must be positive before operator evolution')
    return family, len(family.coefficients[0])


def fixed_target_union(solve_nodes, verification_nodes):
    """Strictly increasing union of solver and held-out verification nodes on I."""
    solve = np.array(_finite_real(solve_nodes, 'solve nodes'), dtype=float, copy=True)
    verify = np.array(_finite_real(verification_nodes, 'verification nodes'), dtype=float, copy=True)
    if solve.ndim != 1 or verify.ndim != 1 or len(solve) < 1 or len(verify) < 1:
        raise ValueError('nonempty one-dimensional solve and verification nodes required')
    if np.any(np.diff(solve) <= 0) or np.any(np.diff(verify) <= 0):
        raise ValueError('increasing solve and verification nodes required')
    union = np.unique(np.concatenate([solve, verify]))
    if union.ndim != 1 or len(union) < 1 or np.any(np.diff(union) <= 0):
        raise ValueError('strictly increasing target union required')
    interval = physical_incoming_interval()
    if np.any(union < interval[0] - 1e-12) or np.any(union > interval[1] + 1e-12):
        raise ValueError('fixed target nodes must lie in I=S(1)+[.12,.18]')
    union.setflags(write=False)
    return union


def subset_node_indices(owned, requested):
    """Exact membership of requested nodes in the frozen union; offsets are absent."""
    owned = np.asarray(owned, float)
    requested = np.asarray(requested, float)
    if owned.ndim != 1 or requested.ndim != 1 or len(requested) < 1:
        raise ValueError('nonempty one-dimensional requested nodes required')
    if len(requested) != len(np.unique(requested)):
        raise ValueError('duplicate requested nodes')
    indices = []
    for value in requested:
        hits = np.flatnonzero(owned == value)
        if hits.size == 0:
            raise ValueError(
                'requested incoming node is absent from the fixed solver/verification union')
        if hits.size > 1:
            raise ValueError('duplicate nodes in the fixed target grid')
        indices.append(int(hits[0]))
    return np.asarray(indices, dtype=int)


def history_operator_cache_key(
        family, *, mass, angular, rho_up, z_grid, target_z, energy_interval,
        degree, axial_support, rtol, atol, max_step, integrator='dop853',
        step_control='primal', tangent_rtol=None, tangent_atol=None):
    """Bind analytic profile/directions and the numerical channel, not amplitudes alone."""
    return (
        profile_identity(family, include_normal_window=True),
        tuple(np.asarray(z_grid, float).tolist()),
        tuple(np.asarray(target_z, float).tolist()),
        tuple(float(v) for v in energy_interval),
        int(degree),
        float(mass),
        float(angular),
        float(rho_up),
        tuple(float(v) for v in axial_support),
        float(rtol),
        float(atol),
        float(max_step),
        str(integrator),
        str(step_control),
        None if tangent_rtol is None else float(tangent_rtol),
        None if tangent_atol is None else float(tangent_atol),
        'tangents=all',
        float(RHO_SIGMA),
    )


def signed_symmetry_applies(positive, negative):
    """True only for the owned S3/opposite-angular negative of one frozen batch."""
    if not isinstance(positive, KSUpstreamBatch) or not isinstance(negative, KSUpstreamBatch):
        return False
    if positive.energy_sign <= 0 or negative.energy_sign >= 0:
        return False
    if (positive.group != negative.group
            or positive.original_panel != negative.original_panel
            or tuple(positive.rows) != tuple(negative.rows)
            or float(positive.mass) != float(negative.mass)
            or float(positive.rho_up) != float(negative.rho_up)):
        return False
    expected_angular = 0.0 if positive.angular == 0.0 else -float(positive.angular)
    if float(negative.angular) != expected_angular:
        return False
    if not np.array_equal(negative.source.energies, -np.asarray(positive.source.energies)):
        return False
    if not np.array_equal(negative.source.column_weights, positive.source.column_weights):
        return False
    residual = source_complement_residual(
        positive.source.covariance, negative.source.covariance)
    if residual > SOURCE_COMPLEMENT_TOLERANCE:
        return False
    return np.array_equal(negative.initial_columns, s3_conjugate(positive.initial_columns))


def group_signed_operator_families(entries):
    """One operator evolution per (mass, angular); map partners instead of a second solve.

    A negative batch is reused from a positive solve only when signed_symmetry_applies.
    Opposite-angular positive-energy families are different operator channels.
    """
    if not isinstance(entries, Sequence) or isinstance(entries, (str, bytes)):
        raise TypeError('signed-family entries required')
    remaining_neg = {
        index for index, item in enumerate(entries) if item[1].energy_sign < 0
    }
    pairs = []
    used_pos = set()
    for index, (_, batch, _) in enumerate(entries):
        if batch.energy_sign <= 0:
            continue
        for other in list(remaining_neg):
            if signed_symmetry_applies(batch, entries[other][1]):
                pairs.append((index, other))
                used_pos.add(index)
                remaining_neg.remove(other)
                break
    mapped_neg = {other for _, other in pairs}
    groups = {}

    def bucket(mass, angular):
        key = (float(mass), float(angular))
        if key not in groups:
            groups[key] = {
                'mass': key[0], 'angular': key[1],
                'applies': [], 'maps': [],
            }
        return groups[key]

    for pos_index, neg_index in pairs:
        pos_id, pos_batch, pos_channel = entries[pos_index]
        neg_id, neg_batch, neg_channel = entries[neg_index]
        group = bucket(pos_batch.mass, pos_batch.angular)
        group['applies'].append((pos_id, pos_batch, pos_channel))
        group['maps'].append((pos_id, neg_id, neg_batch, neg_channel))
    for index, item in enumerate(entries):
        if index in used_pos or index in mapped_neg:
            continue
        batch_id, batch, channel = item
        group = bucket(batch.mass, batch.angular)
        group['applies'].append((batch_id, batch, channel))
    return tuple(groups[key] for key in sorted(groups))


def _normalize_entries(batches):
    if not isinstance(batches, Sequence) or isinstance(batches, (str, bytes)):
        raise TypeError(
            'batches must be a sequence of (KSUpstreamBatch, channel) or '
            '(id, KSUpstreamBatch, channel)')
    entries = []
    for item in batches:
        if not isinstance(item, Sequence) or isinstance(item, (str, bytes)):
            raise TypeError(
                'each entry is (batch, channel) or (id, batch, channel)')
        if len(item) == 2:
            batch, channel = item
            batch_id = None
        elif len(item) == 3:
            batch_id, batch, channel = item
        else:
            raise TypeError(
                'each entry is (batch, channel) or (id, batch, channel)')
        if not isinstance(batch, KSUpstreamBatch):
            raise TypeError(
                'bound KSUpstreamBatch required; prepare A_up/C_src once and freeze them')
        # Copy the already prepared inputs, without repeating preparation.
        # A caller retaining the original NumPy buffers cannot alter this owner.
        batch = replace(batch, source=FixedSourcePreparation(
            batch.source.covariance, batch.source.column_weights, batch.source.energies),
            continuation=_snapshot(batch.continuation), provenance=_snapshot(batch.provenance))
        if isinstance(channel, Mapping) and {'incoming_C0', 'C0'} & set(channel):
            raise ValueError('a frozen incoming C0 is not a production KS batch')
        channel = _channel_record('_', {'_': channel})
        if channel['ell0_pure_radius_identity'] or batch.angular == 0.0:
            raise ValueError(
                'ell=0 radius-response identities are declared analytically; '
                'they are not evolved as numerical batches')
        if (channel['group'] != batch.group
                or channel['angular_sign'] != batch.angular_sign):
            raise ValueError('channel group/sign does not match the bound upstream batch')
        if float(channel['compact_mass']) != float(batch.mass):
            raise ValueError('channel mass does not match the bound upstream batch')
        if float(channel['angular_eigenvalue']) != abs(float(batch.angular)):
            raise ValueError('channel angular eigenvalue does not match the bound batch')
        if batch_id is None:
            batch_id = f'{batch.panel_name}:{batch.energy_sign:+d}:{batch.rows[0]}:{batch.rows[1]}'
        entries.append((str(batch_id), batch, _snapshot(channel)))
    ids = [item[0] for item in entries]
    if not ids:
        raise ValueError('at least one bound KSUpstreamBatch is required')
    if len(ids) != len(set(ids)):
        raise ValueError('duplicate batch IDs')
    return tuple(entries)


def _source_identity(entries):
    items = []
    for batch_id, batch, channel in entries:
        items.append((
            batch_id,
            batch.source.digest,
            batch.preparation_digest,
            str(batch.panel_name),
            tuple(int(v) for v in batch.rows),
            float(batch.mass),
            float(batch.angular),
            float(batch.rho_up),
            int(batch.energy_sign),
            int(batch.group),
            int(batch.angular_sign),
            int(channel['copy_count']),
            int(channel['degeneracy']),
        ))
    return tuple(items)


def _applied_energies(groups):
    labels = []
    for group in groups:
        for _, batch, _ in group['applies']:
            labels.append(np.asarray(batch.source.energies, float))
    if not labels:
        raise ValueError('no operator-applied source energies; mapped partners still need a positive parent')
    return np.concatenate(labels)


def _subset_mapping(raw, requested, indices):
    out = dict(raw)
    out['z'] = _copy_array(requested)
    for key in _NODE_VECTORS:
        if key in out and out[key] is not None:
            out[key] = _copy_array(np.asarray(out[key])[indices])
    for key in _NODE_TANGENTS:
        if key in out and out[key] is not None:
            out[key] = _copy_array(np.asarray(out[key])[:, indices])
    for key in _NODE_MAPS:
        if key in out and isinstance(out[key], Mapping):
            out[key] = {
                name: _copy_array(np.asarray(value)[indices])
                for name, value in out[key].items()
            }
    for key in _TANGENT_MAPS:
        if key in out and isinstance(out[key], Mapping):
            out[key] = {
                name: _copy_array(np.asarray(value)[:, indices])
                for name, value in out[key].items()
            }
    return out


def _freeze_mapping(raw):
    return {key: _snapshot(value) for key, value in raw.items()}


class KSHistoryEvaluator:
    """Production evaluate(family, z) adapter. Finite retained-source CONTROL only."""

    def __init__(
            self, batches, *, baseline_gradient, coefficients, energy_interval,
            interpolant_degree, z_grid, solve_nodes, verification_nodes,
            rtol, atol, max_step, integrator='dop853', step_control='primal',
            tangent_rtol=None, tangent_atol=None, axial_support=None, ell0_channels=(),
            error_budget=None, coverage=None, interval=None, **kwargs):
        _reject_foreign(kwargs)
        entries = _normalize_entries(batches)
        rho = {float(batch.rho_up) for _, batch, _ in entries}
        if len(rho) != 1:
            raise ValueError('all bound batches must share one frozen rho_up')
        groups = group_signed_operator_families(entries)
        interval_e = tuple(float(v) for v in energy_interval)
        if len(interval_e) != 2 or not np.isfinite(interval_e).all() or interval_e[0] >= interval_e[1]:
            raise ValueError('finite positive-length energy interval required')
        if isinstance(interpolant_degree, (bool, np.bool_)) or int(interpolant_degree) != interpolant_degree:
            raise ValueError('integer interpolation degree >= 1 required')
        degree = int(interpolant_degree)
        if degree < 1:
            raise ValueError('integer interpolation degree >= 1 required')
        applied = _applied_energies(groups)
        if np.any((applied < interval_e[0]) | (applied > interval_e[1])):
            raise ValueError(
                'operator-applied source energies must lie in the interpolant interval; '
                'mapped negative partners are not interpolated')
        target = fixed_target_union(solve_nodes, verification_nodes)
        support = usual_axial_support() if axial_support is None else tuple(
            float(v) for v in axial_support)
        incoming = physical_incoming_interval() if interval is None else tuple(
            float(v) for v in interval)
        if not np.allclose(incoming, physical_incoming_interval(), atol=3e-12, rtol=0):
            raise ValueError('incoming interval is the fixed I=S(1)+[.12,.18]')
        ell0 = []
        declared = set()
        for channel in ell0_channels:
            record = _channel_record('_', {'_': channel})
            if record['group'] not in (0, 13, 23) or not record['ell0_pure_radius_identity']:
                raise ValueError('only the owned ell=0 groups 0,13,23 have this radius-response identity')
            if record['group'] in declared:
                raise ValueError('channel response already included')
            if any(batch.group == record['group'] for _, batch, _ in entries):
                raise ValueError('channel response already included')
            declared.add(record['group'])
            ell0.append(record)
        self._entries = entries
        self._groups = groups
        self._baseline = _copy_array(baseline_gradient)
        self._coefficients = _snapshot(coefficients)
        self._energy_interval = interval_e
        self._degree = degree
        self._z_grid = _copy_array(_finite_real(z_grid, 'computational z-grid'))
        self._solve_nodes = _copy_array(solve_nodes)
        self._verification_nodes = _copy_array(verification_nodes)
        self.target_nodes = target
        self._rtol = float(rtol)
        self._atol = float(atol)
        self._max_step = float(max_step)
        if min(self._rtol, self._atol, self._max_step) <= 0 or not np.isfinite(
                [self._rtol, self._atol, self._max_step]).all():
            raise ValueError('explicit finite positive rtol, atol and max_step required')
        if integrator not in ('dop853', 'cf4'):
            raise ValueError("integrator must be 'dop853' or 'cf4'")
        if step_control not in ('joint', 'primal'):
            raise ValueError("step_control must be 'joint' or 'primal'")
        if step_control == 'joint' and (tangent_rtol is not None or tangent_atol is not None):
            raise ValueError('separate tangent tolerances require primal step control')
        if step_control == 'primal':
            tangent_rtol = self._rtol if tangent_rtol is None else float(tangent_rtol)
            tangent_atol = self._atol if tangent_atol is None else float(tangent_atol)
            if min(tangent_rtol, tangent_atol) <= 0 or not np.isfinite(
                    [tangent_rtol, tangent_atol]).all():
                raise ValueError('finite positive tangent tolerances required')
        self._integrator = str(integrator)
        self._step_control = str(step_control)
        self._tangent_rtol = tangent_rtol
        self._tangent_atol = tangent_atol
        self._axial_support = support
        self._ell0 = _snapshot(ell0)
        self._error_budget = _snapshot(local_error_budget() if error_budget is None else error_budget)
        self._coverage = _snapshot(coverage)
        self._interval = incoming
        self._rho_up = next(iter(rho))
        self._source_identity = _source_identity(entries)
        self._operator_cache = {}
        self._assembly_cache = {}
        self.operator_evolution_count = 0
        self.signed_operator_reuse = tuple(
            {
                'mass': group['mass'],
                'angular': group['angular'],
                'applied_batch_ids': tuple(item[0] for item in group['applies']),
                'mapped_partner_ids': tuple(item[1] for item in group['maps']),
                'signed_symmetry_maps': len(group['maps']),
            }
            for group in groups
        )

    @property
    def source_identity(self):
        return self._source_identity

    def evaluate(self, family, z=None):
        family, n = _require_local_family(family)
        requested = self.target_nodes if z is None else np.asarray(z, float)
        if requested.ndim != 1 or len(requested) < 1:
            raise ValueError('nonempty requested incoming z required')
        indices = subset_node_indices(self.target_nodes, requested)
        raw = self._evaluate_union(family, n)
        result = _subset_mapping(raw, requested, indices)
        if not np.array_equal(result['z'], requested):
            raise ValueError('evaluator z must equal the requested nodes exactly')
        ndir = 2 * n
        if result['history_jacobian'].shape != (ndir, len(requested), 2):
            raise ValueError('full retarded history_jacobian of shape (2n, nz, 2) required')
        if result['local_reference_gradient_tangent'].shape != (ndir, len(requested), 2):
            raise ValueError('geometric tangent must have shape (2n, nz, 2)')
        return result

    def _evaluate_union(self, family, n):
        cache_id = profile_identity(family, include_normal_window=True)
        cached = self._assembly_cache.get(cache_id)
        if cached is not None:
            return cached
        acc = KSConstraintAccumulator(
            family, self._baseline, self._coefficients, self._interval,
            self._error_budget, coverage=self._coverage)
        for group in self._groups:
            operator = self._evolve_operator(family, group, n)
            states = {}
            for batch_id, batch, channel in group['applies']:
                prepared = operator.apply(batch.source, batch.initial_columns)
                self._require_reconstructed_source(prepared, batch)
                if not np.array_equal(prepared.z, self.target_nodes):
                    raise ValueError('operator targets must be the frozen solver/verification union')
                if prepared.column_tangents.shape[0] != 2 * n:
                    raise ValueError('every new g must carry all 2n retarded history directions')
                states[batch_id] = prepared
                acc.add(prepared, batch, channel, batch_id=batch_id)
            for pos_id, neg_id, neg_batch, neg_channel in group['maps']:
                partner = negative_angular_partner(states[pos_id], neg_batch.source)
                self._require_frozen_negative(partner, neg_batch)
                acc.add(partner, neg_batch, neg_channel, batch_id=neg_id)
        for channel in self._ell0:
            acc.declare_radius_invariant_channel(channel)
        raw = self._annotate(acc.finalize(), family)
        self._assembly_cache[cache_id] = raw
        return raw

    def _evolve_operator(self, family, group, n):
        key = history_operator_cache_key(
            family, mass=group['mass'], angular=group['angular'],
            rho_up=self._rho_up, z_grid=self._z_grid, target_z=self.target_nodes,
            energy_interval=self._energy_interval, degree=self._degree,
            axial_support=self._axial_support, rtol=self._rtol, atol=self._atol,
            max_step=self._max_step, integrator=self._integrator,
            step_control=self._step_control, tangent_rtol=self._tangent_rtol,
            tangent_atol=self._tangent_atol)
        cached = self._operator_cache.get(key)
        if cached is not None:
            return cached
        self.operator_evolution_count += 1
        operator = evolve_energy_propagator(
            self._energy_interval, self._degree, family, self._z_grid,
            self.target_nodes, group['mass'], group['angular'], self._rho_up,
            axial_support=self._axial_support, rtol=self._rtol, atol=self._atol,
            max_step=self._max_step, tangents='all', integrator=self._integrator,
            step_control=self._step_control, tangent_rtol=self._tangent_rtol,
            tangent_atol=self._tangent_atol)
        if operator.tangent.shape[1] != 2 * n:
            raise ValueError('every new g must carry all 2n retarded history directions')
        self._operator_cache[key] = operator
        return operator

    def _require_reconstructed_source(self, prepared, batch):
        if not isinstance(prepared, KSDifferenceIncoming):
            raise TypeError('KSDifferenceIncoming required from the energy propagator')
        reconstructed = FixedSourcePreparation(
            prepared.source_covariance, prepared.column_weights, prepared.source_energies)
        if reconstructed.digest != batch.source.digest:
            raise ValueError('propagator replaced the frozen source')
        if not np.array_equal(prepared.source_covariance, batch.source.covariance):
            raise ValueError('original C_src columns were not reconstructed')
        if not np.array_equal(prepared.column_weights, batch.source.column_weights):
            raise ValueError('original source weights were not reconstructed')
        if not np.array_equal(prepared.source_energies, batch.source.energies):
            raise ValueError('original source energy labels were not reconstructed')
        if not np.array_equal(prepared.initial_columns, batch.initial_columns):
            raise ValueError('original A_up was not reconstructed')
        if prepared.fixed_preparation_digest != batch.preparation_digest:
            raise ValueError('reconstructed preparation digest differs from the frozen batch')
        if prepared.diagnostics.get('source_interpolated') is True:
            raise ValueError('source values must not be interpolated')

    def _require_frozen_negative(self, partner, batch):
        if not np.array_equal(partner.source_covariance, batch.source.covariance):
            raise ValueError('explicit negative C_src was replaced')
        if not np.array_equal(partner.column_weights, batch.source.column_weights):
            raise ValueError('negative-source column weights must equal the positive weights')
        if not np.array_equal(partner.initial_columns, batch.initial_columns):
            raise ValueError('signed map A_up differs from the bound negative batch')
        if partner.fixed_preparation_digest != batch.preparation_digest:
            raise ValueError('signed map changed the frozen negative preparation')
        if partner.diagnostics.get('energy_folding_factor', 1) != 1:
            raise ValueError('energy folding is not performed; both signs remain explicit')
        if partner.diagnostics.get('field_isometry_assumed') is True:
            raise ValueError('field isometry is not assumed')

    def _annotate(self, raw, family):
        out = _freeze_mapping(raw)
        identity = profile_identity(family, include_normal_window=True)
        scope = dict(out['scope'])
        scope.update({
            'full_retarded_state_derivative': True,
            'finite_retained_source_control': True,
            'source_cutoff_coincidence_bridge': 'OPEN',
            'operator_interpolation_error_bound': None,
            'node_evolution_error_bound': None,
            'fresh_source_evolution_from_fixed_upstream': True,
            'positive_result_doubled': False,
            'energy_folding_factor': 1,
            'production_local_gate': PRODUCTION_LOCAL_GATE,
            'physical_NONEXISTENCE_certificate': False,
            'adapter': 'KSHistoryEvaluator',
            'integrator': self._integrator,
            'step_control': self._step_control,
            'tangent_rtol': self._tangent_rtol,
            'tangent_atol': self._tangent_atol,
        })
        unresolved = list(scope.get('unresolved', ()))
        if 'source-cutoff/coincidence bridge' not in unresolved:
            unresolved.append('source-cutoff/coincidence bridge')
        scope['unresolved'] = unresolved
        out['scope'] = _snapshot(scope)
        out['constraint_order'] = ('N', 'beta')
        out['profile_identity'] = identity
        out['history_identity'] = identity
        out['source_identity'] = self._source_identity
        out['full_retarded_state_derivative'] = True
        out['physical_constraint_status'] = 'OPEN'
        out['physical_EXISTENCE_certificate'] = False
        out['physical_NONEXISTENCE_certificate'] = False
        out['certified_constraint_residual'] = None
        out['full_source_error_bound'] = None
        out['missing_error_budget'] = None
        out['constant_frozen_matter_used_as_physical_state'] = False
        out['reference_state_tangent_subtracted'] = False
        out['incoming_state_derivative_removed'] = False
        out['frozen_incoming_state_derivative'] = False
        return out
