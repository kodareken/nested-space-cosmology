"""Bounded-memory KS N,beta assembly from disjoint three-column source batches.

E_B[g] = E_B[g_ref] + Delta E_B^{local+reference}[j_Sigma g]
         + sum_batches (G_B[C_Sigma[g]] - G_B[C_ref]),  B = N, beta.

Upstream batches come from ReferenceSourcePanel.fixed_upstream and ignore
trial history. Baseline and geometry are added once. Both energy signs stay
explicit. Physical status is OPEN.
"""
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from types import MappingProxyType

import numpy as np

from .nsc_evolved_incoming_constraints import (
    compatible_history_slots, surface_geometry_response,
)
from .nsc_evolved_incoming_state import FixedSourcePreparation
from .nsc_ks_difference_envelope import KSDifferenceIncoming
from .nsc_ks_local_constraints import (
    _MISSING_BOUNDS, _family_matter, _history_family, _incoming_intrinsic,
)
from .nsc_ks_signed_state import (
    SOURCE_COMPLEMENT_TOLERANCE, negative_angular_partner, s3_conjugate,
    source_complement_residual,
)
from .nsc_ks_source_envelope import RHO_SIGMA, _fixed_preparation_digest
from .nsc_ks_source_inventory import ReferenceSourcePanel
from .nsc_local_incoming_constraints import (
    _channel_record, _coverage_record, _finite, _positive_interval,
    _require_error_budget,
)


def _mapping(value):
    return value if isinstance(value, Mapping) else {}


def _panel_identity(panel):
    prov, nested = _mapping(panel.provenance), _mapping(_mapping(panel.provenance).get('source'))
    spec = prov.get('rows', nested.get('rows'))
    count = int(len(panel.energies))
    if spec is not None and (len(spec)!=2 or any(isinstance(v,(bool,np.bool_)) or int(v)!=v for v in spec)):
        raise ValueError('integer original panel row bounds required')
    rows = (0, count) if spec is None else (int(spec[0]), int(spec[1]))
    if rows[1] - rows[0] != count or rows[0] < 0 or rows[1] <= rows[0]:
        raise ValueError('original panel row range must match the bound energy fibers')
    for mapping, key in ((prov, 'parent_panel'), (nested, 'parent_panel'), (prov, 'positive_panel')):
        if key in mapping:
            name = str(mapping[key])
            break
    else:
        name = str(panel.name)
    if name.endswith('/negative-partner'):
        name = name[:-17]
    cut = name.find('/rows')
    return name[:cut] if cut >= 0 else name, rows


def _energy_sign(energies):
    labels = np.asarray(energies, float)
    if np.any(labels == 0) or (np.any(labels > 0) and np.any(labels < 0)):
        raise ValueError('one energy sign per upstream batch; mixed-sign panels are not folded')
    return 1 if np.all(labels > 0) else -1


def _ell0_note(flag):
    return ('included; pure-radius response identically zero by ell=0 Dirac '
            'r-independence; baseline retained; the displayed discrete '
            'current-minus-analytic-reference also includes numerical drift' if flag else None)


@dataclass(frozen=True)
class KSUpstreamBatch:
    """History-independent fixed-upstream batch: panel rows, C_src, A_up, m, ell, rho_up."""

    panel_name: str
    original_panel: str
    rows: object
    group: int
    angular_sign: int
    energy_sign: int
    mass: float
    angular: float
    rho_up: float
    source: object
    initial_columns: object
    preparation_digest: str
    continuation: object
    provenance: object

    def __post_init__(self):
        if not isinstance(self.source, FixedSourcePreparation):
            raise TypeError('FixedSourcePreparation required; a covariance or frozen C0 is not a batch')
        if self.group not in range(1, 33):
            raise ValueError('retained real-field groups are 1 through 32; LLL columns are not invented here')
        if self.angular_sign not in (-1, 1) or self.energy_sign not in (-1, 1):
            raise ValueError('actual angular and energy signs required')
        if self.angular == 0.0 and self.angular_sign != 1:
            raise ValueError('ell=0 keeps angular_sign +1; LLL columns are not invented')
        if len(self.rows)!=2 or any(isinstance(v,(bool,np.bool_)) or int(v)!=v for v in self.rows):
            raise ValueError('integer original panel row bounds required')
        rows = tuple(int(v) for v in self.rows)
        if len(rows) != 2 or rows[1] <= rows[0] or rows[0] < 0:
            raise ValueError('original panel row range [start, end) required')
        if self.source.covariance.shape[0]!=3*(rows[1]-rows[0]):
            raise ValueError('three complete source columns per original energy row required')
        energies=self.source.energies.reshape(-1,3)
        weights=self.source.column_weights.reshape(-1,3)
        if not np.all(energies==energies[:,:1]) or not np.all(weights==weights[:,:1]):
            raise ValueError('source energy labels and weights must repeat across each coherent three-column fiber')
        fibers=np.arange(len(self.source.energies))//3
        if np.any(self.source.covariance[fibers[:,None]!=fibers[None,:]]!=0):
            raise ValueError('cross-energy coherences cannot be split as independent panel fibers')
        if not np.isfinite([self.mass,self.angular,self.rho_up]).all() or self.mass<0 or self.rho_up<1.03:
            raise ValueError('finite owned mass/angular/rho-up parameters required')
        if self.angular!=0 and np.sign(self.angular)!=self.angular_sign:
            raise ValueError('signed angular label does not match angular_sign')
        initial = np.array(self.initial_columns, complex, copy=True)
        if (not np.isfinite(initial).all() or initial.shape != (2, self.source.covariance.shape[0])):
            raise ValueError('finite A_up of shape (2, nsrc) required')
        initial.setflags(write=False)
        if _energy_sign(self.source.energies) != int(self.energy_sign):
            raise ValueError('bound energy sign must match the source labels')
        digest = _fixed_preparation_digest(
            self.source, initial, self.mass, self.angular, self.rho_up, RHO_SIGMA)
        if digest != self.preparation_digest:
            raise ValueError('upstream preparation digest does not match source/A_up/m/ell/rho')
        if self.continuation.get('energy_folding_factor', 1) != 1:
            raise ValueError('energy folding is not performed; both signs remain explicit')
        object.__setattr__(self, 'rows', rows)
        object.__setattr__(self, 'group', int(self.group))
        object.__setattr__(self, 'angular_sign', int(self.angular_sign))
        object.__setattr__(self, 'energy_sign', int(self.energy_sign))
        object.__setattr__(self, 'mass', float(self.mass))
        object.__setattr__(self, 'angular', 0.0 if float(self.angular) == 0.0 else float(self.angular))
        object.__setattr__(self, 'rho_up', float(self.rho_up))
        object.__setattr__(self, 'initial_columns', initial)
        object.__setattr__(self, 'continuation', MappingProxyType(dict(self.continuation or {})))
        object.__setattr__(self, 'provenance', MappingProxyType(dict(self.provenance or {})))


def bind_upstream_batch(panel, *, rho_up, rtol, atol, max_step):
    """Wrap panel.fixed_upstream as an immutable history-independent batch."""
    if not isinstance(panel, ReferenceSourcePanel):
        raise TypeError('ReferenceSourcePanel required')
    if panel.group not in range(1, 33):
        raise ValueError('retained real-field groups are 1 through 32; LLL columns are not invented here')
    original, rows = _panel_identity(panel)
    source, initial, report = panel.fixed_upstream(
        rho_up=rho_up, rtol=rtol, atol=atol, max_step=max_step)
    digest = _fixed_preparation_digest(source, initial, panel.mass, panel.angular, rho_up, RHO_SIGMA)
    return KSUpstreamBatch(
        panel.name, original, rows, panel.group, panel.angular_sign, _energy_sign(panel.energies),
        panel.mass, panel.angular, float(rho_up), source, initial, digest, {
            'panel': report['panel'], 'rho_up': report['rho_up'], 'nfev': report['nfev'],
            'source_digest': report['source_digest'],
            'continuation_error_bound': report['continuation_error_bound'],
            'physical_source_error_bound': report['physical_source_error_bound'],
            'physical_gate': report['physical_gate'],
            'incoming_state_frozen': report['incoming_state_frozen'],
            'energy_folding_factor': 1, 'history_independent': True,
        }, {'parent_panel': original, 'rows': list(rows), 'panel': panel.name,
            'source': _mapping(panel.provenance), 'energy_folding_factor': 1})


def map_negative_batch(batch, negative_source):
    """Opposite-angular negative batch from the owned signed map; no ODE and no factor two."""
    if not isinstance(batch, KSUpstreamBatch):
        raise TypeError('KSUpstreamBatch required')
    if batch.energy_sign <= 0:
        raise ValueError('positive-energy upstream batch required for the negative map')
    if not isinstance(negative_source, FixedSourcePreparation):
        raise TypeError('original FixedSourcePreparation required for the negative source')
    if not np.array_equal(negative_source.energies, -np.asarray(batch.source.energies)):
        raise ValueError('negative-source energies must equal -positive energies')
    if not np.array_equal(negative_source.column_weights, batch.source.column_weights):
        raise ValueError('negative-source column weights must equal the positive weights')
    residual = source_complement_residual(batch.source.covariance, negative_source.covariance)
    if residual > SOURCE_COMPLEMENT_TOLERANCE:
        raise ValueError('negative source differs from the fixed horizon source law')
    angular = 0.0 if batch.angular == 0.0 else -float(batch.angular)
    initial = s3_conjugate(batch.initial_columns)
    digest = _fixed_preparation_digest(
        negative_source, initial, batch.mass, angular, batch.rho_up, RHO_SIGMA)
    return KSUpstreamBatch(
        batch.panel_name + '/negative-partner', batch.original_panel, batch.rows, batch.group,
        1 if angular == 0.0 else -batch.angular_sign, -1, batch.mass, angular, batch.rho_up,
        negative_source, initial, digest, {
            'map': 'S3 conjugation with opposite angular sign', 'energy_folding_factor': 1,
            'positive_preparation_digest': batch.preparation_digest,
            'source_complement_residual': residual, 'history_independent': True,
            'field_isometry_assumed': False, 'physical_gate': 'OPEN',
        }, {'positive_panel': batch.panel_name, 'parent_panel': batch.original_panel,
            'rows': list(batch.rows), 'energy_folding_factor': 1})


class KSConstraintAccumulator:
    """Sum disjoint KSDifferenceIncoming matter corrections; geometry once at finalize."""

    def __init__(self, provider, baseline_gradient, coefficients, interval, error_budget,
                 coverage=None):
        self._history = _history_family(provider)
        self._baseline, self._coefficients = baseline_gradient, coefficients
        self._interval = _positive_interval(interval)
        self._bounds, extra = _require_error_budget(error_budget)
        declared = {k: extra[k] for k in (
            'sampled_family_ids', 'missing_family_ids', 'sampled_energy_ranges',
            'missing_energy_ranges', 'energy_coverage') if k in extra}
        self._coverage = coverage if coverage is not None else declared or None
        self._certificate = extra.get('independent_source_error_certificate')
        axial, radius = _incoming_intrinsic()
        if not np.allclose([coefficients['a'], coefficients['r']], [axial, radius], atol=3e-12, rtol=0):
            raise ValueError(
                'constraint coefficients and prepared states must use the same intrinsic identification')
        self._z = self._ndir = None
        self._corrections, self._currents, self._references, self._tangents = {}, {}, {}, {}
        self._records, self._occupied, self._families = {}, {}, {}
        self._analytic_zero_groups = {}
        self._finalized = False

    def declare_radius_invariant_channel(self, channel):
        """Reuse the owned ell=0 identity; its baseline is still counted once."""
        if self._finalized:
            raise ValueError('cannot declare channels after finalization')
        channel=_channel_record('_',{'_':channel})
        group=channel['group']
        if group not in (0,13,23) or channel['angular_eigenvalue']!=0 or channel['angular_sign']!=1:
            raise ValueError('only the owned ell=0 groups 0,13,23 have this radius-response identity')
        if group in self._analytic_zero_groups or any(r['group']==group for r in self._records.values()):
            raise ValueError('channel response already included')
        self._analytic_zero_groups[group]={
            'group':group,'compact_mass':channel['compact_mass'],'angular_eigenvalue':0.,
            'copy_count':channel['copy_count'],'degeneracy':channel['degeneracy'],
            'matter_change_exact':(0.,0.),'matter_history_tangent_exact':0.,
            'all_source_energies':True,'baseline_retained':True,
            'reason':'pure-radius KS generator has no r dependence at ell=0; intrinsic Sigma vertex unchanged'}

    def add(self, prepared, batch, channel, batch_id=None):
        if self._finalized:
            raise ValueError('cannot add batches after finalization')
        if isinstance(prepared, dict):
            raise ValueError(
                'a frozen incoming C0, PG HistoryBinding, or node_fields dictionary '
                'is not a KS difference batch')
        if not isinstance(prepared, KSDifferenceIncoming):
            raise TypeError(
                'KSDifferenceIncoming required; the batched assembly owns the joint '
                'reference/difference incoming state')
        if not isinstance(batch, KSUpstreamBatch):
            raise TypeError('bound KSUpstreamBatch required')
        prepared.validate()
        prepared.require_history(self._history)
        source = FixedSourcePreparation(
            prepared.source_covariance, prepared.column_weights, prepared.source_energies)
        if source.digest != batch.source.digest:
            raise ValueError('prepared source differs from the bound upstream batch')
        if (prepared.fixed_preparation_digest != batch.preparation_digest
                or not np.array_equal(prepared.initial_columns, batch.initial_columns)):
            raise ValueError('prepared A_up/preparation digest differs from the bound upstream batch')
        channel = _channel_record('_', {'_': channel})
        if channel['group'] in self._analytic_zero_groups:
            raise ValueError('channel response already included by the exact ell=0 identity')
        if channel['group'] != batch.group or channel['angular_sign'] != batch.angular_sign:
            raise ValueError('channel group/sign does not match the bound upstream batch')
        batch_id = batch_id or f'{batch.panel_name}:{batch.energy_sign:+d}:{batch.rows[0]}:{batch.rows[1]}'
        if batch_id in self._records:
            raise ValueError('duplicate batch IDs')
        occupied = self._occupied.setdefault((batch.group, batch.angular_sign, batch.original_panel, batch.energy_sign), [])
        for other_start, other_end, other_id in occupied:
            if batch.rows[0] < other_end and other_start < batch.rows[1]:
                raise ValueError(
                    'overlapping original panel rows for the same energy sign: '
                    f'{batch.original_panel} {batch.rows} overlaps {(other_start, other_end)} ({other_id})')
        z, ndir = np.asarray(prepared.z), int(prepared.column_tangents.shape[0])
        if self._z is None:
            self._z, self._ndir = z, ndir
        elif not np.array_equal(z, self._z) or ndir != self._ndir:
            raise ValueError('mismatched incoming z or direction order across prepared batches')
        matter = _family_matter(prepared, channel, self._coefficients, self._ndir, self._z)
        family = (channel['group'], channel['angular_sign'], batch.energy_sign)
        previous = self._families.get(family)
        if previous is not None and previous['multiplicity'] != matter['multiplicity']:
            raise ValueError('original signed-family multiplicity must be applied once')
        if previous is None:
            self._families[family] = {
                'group': channel['group'], 'angular_sign': channel['angular_sign'],
                'energy_sign': batch.energy_sign, 'compact_mass': channel['compact_mass'],
                'angular_eigenvalue': channel['angular_eigenvalue'],
                'signed_angular': matter['signed_angular'], 'copy_count': channel['copy_count'],
                'degeneracy': channel['degeneracy'],
                'n_actual_angular_signs': channel['n_actual_angular_signs'],
                'multiplicity': matter['multiplicity'],
                'ell0_pure_radius_identity': matter['ell0_pure_radius_identity'],
                'batch_ids': [batch_id],
            }
        else:
            previous['batch_ids'].append(batch_id)
        self._corrections[batch_id] = matter['correction']
        self._currents[batch_id] = matter['current']
        self._references[batch_id] = matter['reference']
        self._tangents[batch_id] = matter['tangent']
        self._records[batch_id] = {
            'group': channel['group'], 'angular_sign': channel['angular_sign'],
            'energy_sign': batch.energy_sign, 'original_panel': batch.original_panel,
            'panel_name': batch.panel_name, 'rows': batch.rows,
            'compact_mass': channel['compact_mass'],
            'angular_eigenvalue': channel['angular_eigenvalue'],
            'signed_angular': matter['signed_angular'], 'copy_count': channel['copy_count'],
            'degeneracy': channel['degeneracy'],
            'n_actual_angular_signs': channel['n_actual_angular_signs'],
            'multiplicity': matter['multiplicity'],
            'preparation_digest': matter['preparation_digest'],
            'source_digest': matter['source_digest'],
            'computational_count': matter['computational_count'],
            'stable_difference': True,
            'ell0_pure_radius_identity': matter['ell0_pure_radius_identity'],
            'ell0_documentation': _ell0_note(matter['ell0_pure_radius_identity']),
            'energy_folding_factor': 1, 'field_isometry_assumed': False,
            'history_independent_upstream': True,
        }
        occupied.append((batch.rows[0], batch.rows[1], batch_id))
        return batch_id

    def finalize(self):
        if not self._records:
            raise ValueError(
                'at least one prepared batch required; omitted families remain explicit OPEN')
        self._finalized = True
        z, ndir = self._z, self._ndir
        slots, slot_tangents = compatible_history_slots(self._history, z, ndir)
        geometric = surface_geometry_response(slots, slot_tangents, self._coefficients)
        baseline = _finite(self._baseline, 'baseline complete-assembly approximant', real=True)
        if baseline.shape == (2,):
            nodes = np.broadcast_to(baseline, (len(z), 2)).copy()
        elif baseline.shape == (len(z), 2):
            nodes = np.array(baseline, dtype=float, copy=True)
        else:
            raise ValueError('baseline N,beta assembly approximant required with shape (2,) or (nz, 2)')
        family_corrections, family_tangents, family_records = {}, {}, {}
        for key, record in self._families.items():
            ids = record['batch_ids']
            family_id = f'{key[0]}_{key[1]}_E{key[2]:+d}'
            family_corrections[family_id] = sum(self._corrections[i] for i in ids)
            family_tangents[family_id] = sum(self._tangents[i] for i in ids)
            family_records[family_id] = {
                **record, 'batch_ids': tuple(ids),
                'ell0_documentation': _ell0_note(record['ell0_pure_radius_identity']),
            }
        unresolved = [label for name, label in _MISSING_BOUNDS if self._bounds[name] is None]
        unresolved.append('complete source/error certificate')
        scope = MappingProxyType({
            'interval': list(self._interval), 'production_local_gate': 'I=S(1)+[.12,.18]',
            'global_matching': False, 'extended_stationarity': 'OPEN',
            'metric_timestep': False, 'observations': False,
            'physical_constraint_status': 'OPEN', 'physical_EXISTENCE_certificate': False,
            'certified_constraint_residual': None, 'full_source_error_bound': None,
            'baseline_and_geometry_counted_once': True,
            'baseline_and_geometry_added_per_batch': False, 'reference_tangent': 0,
            'reference_state_tangent_subtracted': False,
            'constant_frozen_matter_used_as_physical_state': False,
            'pde_axial_derivatives': True, 'interpolated_axial_derivatives': False,
            'homogeneous_ks_generator_reference': True,
            'computational_meshes_need_not_match': True,
            'finite_energy_samples_are_complete_spectrum': False,
            'dense_all_energy_covariance_allocated': False, 'energy_folding_factor': 1,
            'positive_result_doubled': False, 'field_isometry_assumed': False,
            'history_independent_upstream_batches': True,
            'independent_source_error_certificate': self._certificate,
            **_coverage_record(family_records, self._coverage),
            'sampled_batch_ids': tuple(self._records),
            'ell0_pure_radius_identities': {
                i: rec['ell0_documentation'] for i, rec in self._records.items()
                if rec['ell0_pure_radius_identity']
            },
            'unresolved': unresolved,
            'analytic_zero_response_groups': dict(self._analytic_zero_groups),
        })
        return {
            'constraint_order': ('N', 'beta'), 'z': np.array(z, copy=True),
            'action_gradient': nodes + geometric['action_gradient_change'] + sum(self._corrections.values()),
            'history_jacobian': geometric['action_gradient_tangent'] + sum(self._tangents.values()),
            'baseline_action_gradient': nodes,
            'local_reference_gradient_change': geometric['action_gradient_change'],
            'local_reference_gradient_tangent': geometric['action_gradient_tangent'],
            'incoming_slots': slots, 'incoming_slot_tangents': slot_tangents,
            'batch_corrections': self._corrections, 'batch_current_matter': self._currents,
            'batch_reference_matter': self._references, 'batch_matter_tangents': self._tangents,
            'batch_records': self._records, 'family_corrections': family_corrections,
            'family_matter_tangents': family_tangents, 'family_records': family_records,
            'error_budget': self._bounds, 'interval': self._interval, 'scope': scope,
            'physical_constraint_status': 'OPEN', 'physical_EXISTENCE_certificate': False,
        }


def assemble_batched_incoming(
        entries, provider, baseline_gradient, coefficients, interval, error_budget,
        coverage=None):
    """Assemble disjoint KSDifferenceIncoming batches with original ledger multiplicity."""
    acc = KSConstraintAccumulator(
        provider, baseline_gradient, coefficients, interval, error_budget, coverage=coverage)
    if not isinstance(entries, Sequence) or isinstance(entries, (str, bytes)):
        raise TypeError('each entry is (prepared, batch, channel) or (id, prepared, batch, channel)')
    for entry in entries:
        if not isinstance(entry, Sequence) or isinstance(entry, (str, bytes)):
            raise TypeError('each entry is (prepared, batch, channel) or (id, prepared, batch, channel)')
        if len(entry) == 3:
            acc.add(*entry)
        elif len(entry) == 4:
            acc.add(entry[1], entry[2], entry[3], batch_id=entry[0])
        else:
            raise TypeError('each entry is (prepared, batch, channel) or (id, prepared, batch, channel)')
    return acc.finalize()


def signed_incoming_pair(prepared, batch, negative_source):
    """Map one evolved positive batch to its explicit negative partner; weights once."""
    partner = negative_angular_partner(prepared, negative_source)
    return prepared, batch, partner, map_negative_batch(batch, negative_source)
