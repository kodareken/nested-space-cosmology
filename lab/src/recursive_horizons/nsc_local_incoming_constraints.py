"""Bounded local incoming N,beta assembly from exact-phase families.

The local residual on a declared interval I is the owned source-fixed split

    E_B[g] = E_B[g_ref] + Delta E_B^{local+reference}[j_Sigma g]
             + sum_families (G_B[C_Sigma[g]] - G_B[C_ref]),  B = N, beta.

Baseline and geometry enter once. Current matter uses weighted exact-phase
PDE columns F, F_z and tangents dF, dF_z. The computational reference is the
stationary phase of the same initial fields through the fixed restriction R;
its tangent is identically zero. No free C0 or reference-stress dictionary
is accepted. Finite sampled families and energies are not a complete spectrum.
Physical constraint status remains OPEN until an independent source/error
certificate is supplied outside this residual.
"""
from __future__ import annotations

from collections.abc import Mapping, Sequence
from types import MappingProxyType

import numpy as np

from .nsc_evolved_incoming_constraints import (
    compatible_history_slots, source_column_matter, surface_geometry_response,
)
from .nsc_evolved_incoming_state import FixedSourcePreparation
from .nsc_exact_phase_prepared_state import ExactPhaseIncoming, _fixed_preparation_digest
from .nsc_local_incoming_family import LocalIncomingFamily
from .nsc_transmitting_dirac_domain import TransmittingDiracSeamDomain


STATIONARITY_TOLERANCE = 3e-11
ERROR_BUDGET_COMPONENTS = (
    'baseline_covered_regions',
    'baseline_low_subgap',
    'changed_history_finite_energy',
    'changed_history_tail',
    'preparation_field_axial_derivative',
    'coefficient_arithmetic',
    'between_node_remainder',
)
_ALLOWED_COVERAGE_ENERGY = {None, 'OPEN', 'partial', False, 'sampled'}


def _finite(value, name, *, real=False):
    if value is None:
        raise ValueError('explicit ' + name + ' required; missing data are not zero')
    raw = np.asarray(value)
    if real and np.iscomplexobj(raw) and np.any(raw.imag != 0):
        raise ValueError('real ' + name + ' required')
    result = np.asarray(raw.real if real else raw, dtype=float if real else complex)
    if not np.isfinite(result).all():
        raise ValueError('finite ' + name + ' required')
    return result


def _positive_interval(interval, name='interval'):
    if interval is None:
        raise ValueError(name + ' required; a missing interval cannot pass a local residual')
    bounds = _finite(interval, name, real=True)
    if bounds.shape != (2,) or not np.isfinite(bounds).all() or float(bounds[1]) <= float(bounds[0]):
        raise ValueError('connected positive-length ' + name + ' required')
    return (float(bounds[0]), float(bounds[1]))


def _component_bound(value, name):
    """Length-2 N,beta bound, or None. Missing/None is not converted to zero."""
    if value is None:
        return None
    bound = _finite(value, name, real=True)
    if bound.shape != (2,) or np.any(bound < 0):
        raise ValueError(name + ' must be a nonnegative componentwise N,beta bound or None')
    return np.array(bound, dtype=float, copy=True)


def local_error_budget(
        baseline_covered_regions=None,
        baseline_low_subgap=None,
        changed_history_finite_energy=None,
        changed_history_tail=None,
        preparation_field_axial_derivative=None,
        coefficient_arithmetic=None,
        between_node_remainder=None):
    """Minimal error-budget mapping. Every component is an N,beta pair or None."""
    return {
        'baseline_covered_regions': _component_bound(
            baseline_covered_regions, 'baseline_covered_regions'),
        'baseline_low_subgap': _component_bound(
            baseline_low_subgap, 'baseline_low_subgap'),
        'changed_history_finite_energy': _component_bound(
            changed_history_finite_energy, 'changed_history_finite_energy'),
        'changed_history_tail': _component_bound(
            changed_history_tail, 'changed_history_tail'),
        'preparation_field_axial_derivative': _component_bound(
            preparation_field_axial_derivative, 'preparation_field_axial_derivative'),
        'coefficient_arithmetic': _component_bound(
            coefficient_arithmetic, 'coefficient_arithmetic'),
        'between_node_remainder': _component_bound(
            between_node_remainder, 'between_node_remainder'),
    }


def _require_error_budget(error_budget):
    if error_budget is None:
        raise ValueError('error_budget required; missing bounds cannot become zero')
    if not isinstance(error_budget, Mapping):
        raise TypeError('error_budget must be a mapping of componentwise N,beta bounds or None')
    unknown = sorted(set(error_budget) - set(ERROR_BUDGET_COMPONENTS) - {
        'independent_source_error_certificate', 'sampled_family_ids',
        'missing_family_ids', 'sampled_energy_ranges', 'missing_energy_ranges',
        'energy_coverage',
    })
    if unknown:
        raise ValueError('unsupported error-budget keys: ' + ', '.join(unknown))
    missing = [name for name in ERROR_BUDGET_COMPONENTS if name not in error_budget]
    if missing:
        raise ValueError(
            'error_budget missing mandatory components ' + ', '.join(missing)
            + '; omitted bounds are not zero')
    normalized = local_error_budget(**{
        name: error_budget[name] for name in ERROR_BUDGET_COMPONENTS
    })
    extra = {
        key: error_budget[key]
        for key in error_budget
        if key not in ERROR_BUDGET_COMPONENTS
    }
    return normalized, extra


def _metric_provider(provider):
    if isinstance(provider, LocalIncomingFamily):
        return provider.metric()
    return provider


def _as_family_map(prepared_families):
    if prepared_families is None:
        raise ValueError('prepared_families required')
    if isinstance(prepared_families, Mapping):
        items = list(prepared_families.items())
    elif isinstance(prepared_families, Sequence) and not isinstance(prepared_families, (str, bytes)):
        items = []
        for entry in prepared_families:
            if isinstance(entry, Mapping) and 'id' in entry and 'prepared' in entry:
                items.append((entry['id'], entry['prepared']))
            elif isinstance(entry, Sequence) and len(entry) == 2:
                items.append((entry[0], entry[1]))
            else:
                raise TypeError(
                    'prepared families must be a mapping or a sequence of (id, ExactPhaseIncoming)')
    else:
        raise TypeError('prepared families must be a mapping of family IDs to ExactPhaseIncoming')
    ids = [entry[0] for entry in items]
    if len(ids) != len(set(ids)):
        raise ValueError('duplicate family IDs')
    if not ids:
        raise ValueError(
            'at least one prepared family required; omitted families remain explicit OPEN')
    families = {}
    for family_id, prepared in items:
        if isinstance(prepared, dict) and (
                'stress_approximant' in prepared or 'parts' in prepared
                or 'incoming_C0' in prepared or 'C0' in prepared):
            raise ValueError(
                'a frozen incoming C0 or legacy matter dictionary is not an exact-phase family')
        if not isinstance(prepared, ExactPhaseIncoming):
            raise TypeError(
                'ExactPhaseIncoming required for each family; frozen C0, interpolated '
                'axial columns, or a bare evolved state without PDE F_z are not this assembly')
        prepared.validate()
        families[family_id] = prepared
    return families


def _channel_record(family_id, channels):
    if not isinstance(channels, Mapping):
        raise TypeError('channels must map each family ID to its retained-ledger record')
    if family_id not in channels:
        raise ValueError('channel record required for family ' + repr(family_id))
    channel = dict(channels[family_id])
    group = channel.get('index', channel.get('group'))
    if group is None:
        raise ValueError('each channel must bind an actual group/index')
    if 'angular_sign' not in channel:
        raise ValueError('each prepared family must be bound to its actual angular_sign')
    sign = channel['angular_sign']
    if isinstance(sign, (bool, np.bool_)) or int(sign) != sign or int(sign) not in (-1, 1):
        raise ValueError('angular_sign must be +1 or -1')
    copies, degeneracy = channel['copy_count'], channel['degeneracy']
    if any(isinstance(value, (bool, np.bool_)) or int(value) != value or value <= 0
           for value in (copies, degeneracy)):
        raise ValueError('positive integer copies/degeneracy from the retained ledger required')
    angular = float(channel['angular_eigenvalue'])
    signs = (1,) if angular == 0 else (-1, 1)
    if int(sign) not in signs:
        raise ValueError('owned angular sign required; no energy-sign folding is performed')
    return {
        **channel,
        'index': int(group),
        'group': int(group),
        'angular_sign': int(sign),
        'compact_mass': float(channel['compact_mass']),
        'angular_eigenvalue': angular,
        'copy_count': int(copies),
        'degeneracy': int(degeneracy),
        'n_actual_angular_signs': len(signs),
        'multiplicity': int(copies) * int(degeneracy) / len(signs),
        'ell0_pure_radius_identity': angular == 0.,
    }


def _fixed_source(prepared):
    source = FixedSourcePreparation(
        prepared.state.source_covariance,
        prepared.state.column_weights,
        prepared.state.source_energies,
    )
    if source.digest != prepared.state.history.source_digest:
        raise ValueError('prepared source arrays do not match the bound source digest')
    return source


def _initial_fields(prepared):
    fields = prepared.diagnostics.get('node_fields')
    if fields is None:
        raise ValueError(
            'ExactPhaseIncoming diagnostics[node_fields][0] required to rebuild the '
            'preparation digest and the stationary computational reference')
    phi = np.asarray(fields)
    if phi.ndim != 3 or phi.shape[0] < 1:
        raise ValueError('node_fields must cover every time node')
    return np.asarray(phi[0], complex)


def _verify_preparation_digest(prepared, channel, source, phi, x):
    """Rebuild the private m/ell digest from initial fields and the actual source."""
    mass = channel['compact_mass']
    angular = channel['angular_sign'] * channel['angular_eigenvalue']
    digest = _fixed_preparation_digest(
        source, phi, x, mass, angular, float(prepared.state.history.times[0]))
    if digest != prepared.fixed_preparation_digest:
        raise ValueError(
            'channel mass/sign does not match the bound preparation digest rebuilt '
            'from the initial fields and actual source')
    if digest != prepared.diagnostics.get('fixed_preparation_digest'):
        raise ValueError('fixed-preparation digest is not bound to the stored diagnostics')
    expected_source = channel.get('expected_source_digest')
    if expected_source is not None and source.digest != expected_source:
        raise ValueError('prepared source differs from its declared family source')
    return digest


def stationary_computational_reference(prepared):
    """Stationary Fref, Fref_z from the same initial fields, fixed R, times and labels.

    This is a computational subtraction reference. It is not a free C0, a
    physical incoming covariance, or a reference-stress dictionary.
    """
    if not isinstance(prepared, ExactPhaseIncoming):
        raise TypeError('ExactPhaseIncoming required')
    prepared.validate()
    phi = _initial_fields(prepared)
    history = prepared.state.history
    x = np.asarray(history.spatial_grid)
    rho1 = int(history.rho1_index)
    take = np.array([rho1, len(x) + rho1], dtype=np.intp)
    if phi.ndim != 2 or phi.shape[0] != 2 * len(x) or phi.shape[1] != prepared.state.columns.shape[-1]:
        raise ValueError('initial node fields must have shape (2*len(x), nsrc)')
    trace = np.asarray(history.restriction) @ phi[take]
    times = np.asarray(history.times, float)
    energies = np.asarray(prepared.state.source_energies, float)
    phase = np.exp(-1j * (times - times[0])[:, None, None] * energies)
    reference = trace[None] * phase
    reference_z = -1j * energies * reference
    return reference, reference_z


def _family_matter(prepared, channel, coefficients, ndir):
    source = _fixed_source(prepared)
    phi = _initial_fields(prepared)
    x = prepared.state.history.spatial_grid
    digest = _verify_preparation_digest(prepared, channel, source, phi, x)
    F = np.asarray(prepared.weighted_columns)
    Fz = np.asarray(prepared.weighted_axial_columns)
    dF = np.asarray(prepared.weighted_column_tangents)
    dFz = np.asarray(prepared.weighted_axial_tangents)
    if dF.shape[0] != ndir or dFz.shape[0] != ndir:
        raise ValueError('prepared direction count must match the assembled history')
    reference, reference_z = stationary_computational_reference(prepared)
    weights = np.asarray(prepared.state.column_weights)
    zeros = np.zeros((ndir, *reference.shape), complex)
    parameters = dict(
        mass=channel['compact_mass'],
        angular=channel['angular_sign'] * channel['angular_eigenvalue'],
        axial_scale=float(coefficients['a']),
        radius=float(coefficients['r']),
        multiplicity=channel['multiplicity'],
    )
    current = source_column_matter(F, Fz, prepared.state.source_covariance, dF, dFz, **parameters)
    reference_matter = source_column_matter(
        reference * weights, reference_z * weights,
        prepared.state.source_covariance, zeros, zeros, **parameters)
    correction = current['action_gradient'] - reference_matter['action_gradient']
    # Reference tangent is identically zero; do not subtract a reference-history derivative.
    tangent = current['action_gradient_tangent']
    return {
        'current': current['action_gradient'],
        'reference': reference_matter['action_gradient'],
        'correction': correction,
        'tangent': tangent,
        'multiplicity': parameters['multiplicity'],
        'signed_angular': parameters['angular'],
        'preparation_digest': digest,
        'source_digest': source.digest,
        'coherent_source_retained': True,
        'reference_tangent_zero': True,
        'pde_axial_derivatives': True,
        'outgoing_momentum_replaced_by_source_energy': False,
        'ell0_pure_radius_identity': channel['ell0_pure_radius_identity'],
        'weights_applied_once': True,
    }


def _coverage_record(family_ids, coverage):
    sampled = tuple(family_ids)
    record = {
        'sampled_family_ids': sampled,
        'missing_family_ids': 'OPEN',
        'sampled_energy_ranges': 'OPEN',
        'missing_energy_ranges': 'OPEN',
        'energy_coverage': 'OPEN',
        'finite_energy_samples_are_complete_spectrum': False,
        'all_retained_ids_cannot_certify_energy_coverage': True,
    }
    if coverage is None:
        return record
    if not isinstance(coverage, Mapping):
        raise TypeError('coverage must be a mapping declared by the caller')
    declared_ids = coverage.get('sampled_family_ids', sampled)
    declared = tuple(declared_ids)
    if set(declared) != set(sampled):
        raise ValueError(
            'declared sampled family coverage must match the assembled family IDs')
    record['sampled_family_ids'] = sampled
    if 'missing_family_ids' in coverage:
        missing = coverage['missing_family_ids']
        if missing is None:
            record['missing_family_ids'] = 'OPEN'
        elif isinstance(missing, str):
            if missing != 'OPEN':
                raise ValueError('missing family IDs are an explicit list or OPEN')
            record['missing_family_ids'] = 'OPEN'
        else:
            record['missing_family_ids'] = tuple(missing)
    energy = coverage.get('energy_coverage', 'OPEN')
    if energy not in _ALLOWED_COVERAGE_ENERGY:
        raise ValueError(
            'finite sampled families cannot certify energy coverage or a complete spectrum')
    if coverage.get('finite_energy_samples_are_complete_spectrum') or coverage.get(
            'complete_spectrum') is True:
        raise ValueError(
            'finite energy samples are not a complete spectrum; even all 33 family IDs '
            'cannot certify energy coverage')
    if 'sampled_energy_ranges' in coverage:
        ranges = coverage['sampled_energy_ranges']
        record['sampled_energy_ranges'] = 'OPEN' if ranges in (None, False, 'OPEN') else ranges
    if 'missing_energy_ranges' in coverage:
        missing_e = coverage['missing_energy_ranges']
        if missing_e is None or missing_e is False:
            record['missing_energy_ranges'] = 'OPEN'
        elif isinstance(missing_e, str):
            if missing_e != 'OPEN':
                raise ValueError('missing energy ranges are an explicit list or OPEN')
            record['missing_energy_ranges'] = 'OPEN'
        else:
            record['missing_energy_ranges'] = tuple(missing_e)
    record['energy_coverage'] = 'OPEN'
    record['finite_energy_samples_are_complete_spectrum'] = False
    return record


def assemble_local_incoming(
        prepared_families, provider, baseline_gradient, coefficients,
        channels, interval, error_budget, coverage=None):
    """Assemble raw local N,beta residuals from exact-phase families.

    prepared_families: mapping (or (id, state) sequence) of ExactPhaseIncoming.
    provider: LocalIncomingFamily or CompatibleIncomingMetric history.
    baseline_gradient: complete-assembly N,beta approximant, shape (2,) or (nz, 2).
    coefficients: owned surface-coefficient mapping with D[0]=0.
    channels: mapping family ID -> retained group, mass, ell, sign, copies, degeneracy.
    interval: connected positive I, production gate I=S(1)+[.12,.18].
    error_budget: mandatory componentwise N,beta bounds or None; missing is not 0.
    coverage: optional caller declaration of sampled/missing family IDs and energies.
    """
    families = _as_family_map(prepared_families)
    provider = _metric_provider(provider)
    interval = _positive_interval(interval)
    bounds, budget_extra = _require_error_budget(error_budget)
    coverage_record = _coverage_record(families, coverage if coverage is not None else {
        key: budget_extra[key]
        for key in (
            'sampled_family_ids', 'missing_family_ids', 'sampled_energy_ranges',
            'missing_energy_ranges', 'energy_coverage',
        )
        if key in budget_extra
    } or None)

    first = next(iter(families.values()))
    times = np.asarray(first.state.history.times)
    x = np.asarray(first.state.history.spatial_grid)
    z = np.asarray(first.state.z)
    rho1_metric = np.asarray(first.state.history.rho1_metric)
    restriction = np.asarray(first.state.history.restriction)
    ndir = int(first.state.column_tangents.shape[0])
    first.require_history(provider, times, x)

    N, beta, q, r = rho1_metric
    intrinsic = TransmittingDiracSeamDomain(lapse=N, radial_scale=q, shift=beta, radius=r)
    if not np.allclose(
            [coefficients['a'], coefficients['r']],
            [intrinsic.induced_axial_scale, r], atol=3e-12, rtol=0):
        raise ValueError(
            'constraint coefficients and prepared states must use the same intrinsic identification')

    seen_groups = {}
    records = {}
    corrections = {}
    currents = {}
    references = {}
    tangents = {}
    for family_id, prepared in families.items():
        if not np.array_equal(prepared.state.history.times, times):
            raise ValueError('mismatched times across prepared families')
        if not np.array_equal(prepared.state.history.spatial_grid, x):
            raise ValueError('mismatched grids across prepared families')
        if float(np.max(np.abs(prepared.state.history.rho1_metric - rho1_metric))) > 3e-12:
            raise ValueError('mismatched intrinsic rho1 data across prepared families')
        if float(np.max(np.abs(prepared.state.history.restriction - restriction))) > 3e-12:
            raise ValueError('mismatched intrinsic restriction across prepared families')
        if not np.array_equal(prepared.state.z, z):
            raise ValueError('mismatched incoming z across prepared families')
        prepared.require_history(provider, times, x)
        channel = _channel_record(family_id, channels)
        key = (channel['group'], channel['angular_sign'])
        if key in seen_groups:
            raise ValueError('duplicate family IDs for the same group and angular_sign')
        seen_groups[key] = family_id
        matter = _family_matter(prepared, channel, coefficients, ndir)
        # Each signed family owns its source and numerical energy sampling.
        # Equality of source/weight digests across different angular signs
        # is not implied by sharing a mass/group label.
        records[family_id] = {
            'group': channel['group'],
            'angular_sign': channel['angular_sign'],
            'compact_mass': channel['compact_mass'],
            'angular_eigenvalue': channel['angular_eigenvalue'],
            'signed_angular': matter['signed_angular'],
            'copy_count': channel['copy_count'],
            'degeneracy': channel['degeneracy'],
            'n_actual_angular_signs': channel['n_actual_angular_signs'],
            'multiplicity': matter['multiplicity'],
            'preparation_digest': matter['preparation_digest'],
            'source_digest': matter['source_digest'],
            'ell0_pure_radius_identity': matter['ell0_pure_radius_identity'],
            'ell0_documentation': (
                'included; pure-radius response identically zero by ell=0 Dirac '
                'r-independence; baseline retained; the displayed discrete '
                'current-minus-analytic-reference also includes numerical drift'
                if matter['ell0_pure_radius_identity'] else None
            ),
        }
        corrections[family_id] = matter['correction']
        currents[family_id] = matter['current']
        references[family_id] = matter['reference']
        tangents[family_id] = matter['tangent']

    slots, slot_tangents = compatible_history_slots(provider, z, ndir)
    geometric = surface_geometry_response(slots, slot_tangents, coefficients)
    baseline = _finite(baseline_gradient, 'baseline complete-assembly approximant', real=True)
    if baseline.shape == (2,):
        baseline_nodes = np.broadcast_to(baseline, (len(z), 2)).copy()
    elif baseline.shape == (len(z), 2):
        baseline_nodes = np.array(baseline, dtype=float, copy=True)
    else:
        raise ValueError('baseline N,beta assembly approximant required with shape (2,) or (nz, 2)')

    matter_sum = sum(corrections.values())
    residual = baseline_nodes + geometric['action_gradient_change'] + matter_sum
    jacobian = geometric['action_gradient_tangent'] + sum(tangents.values())
    ell0 = {
        family_id: record['ell0_documentation']
        for family_id, record in records.items()
        if record['ell0_pure_radius_identity']
    }
    scope = {
        'interval': list(interval),
        'production_local_gate': 'I=S(1)+[.12,.18]',
        'global_matching': False,
        'extended_stationarity': 'OPEN',
        'metric_timestep': False,
        'observations': False,
        'physical_constraint_status': 'OPEN',
        'physical_EXISTENCE_certificate': False,
        'certified_constraint_residual': None,
        'full_source_error_bound': None,
        'baseline_and_geometry_counted_once': True,
        'reference_tangent': 0,
        'reference_state_tangent_subtracted': False,
        'constant_frozen_matter_used_as_physical_state': False,
        'pde_axial_derivatives': True,
        'interpolated_axial_derivatives': False,
        'finite_energy_samples_are_complete_spectrum': False,
        'independent_source_error_certificate': budget_extra.get(
            'independent_source_error_certificate'),
        **coverage_record,
        'ell0_pure_radius_identities': ell0,
        'unresolved': [
            'baseline low-energy/subgap source accuracy' if bounds['baseline_low_subgap'] is None
            else None,
            'changed-history finite-energy coverage' if bounds['changed_history_finite_energy'] is None
            else None,
            'changed-history tail' if bounds['changed_history_tail'] is None else None,
            'preparation/field/axial-derivative error'
            if bounds['preparation_field_axial_derivative'] is None else None,
            'coefficient interval and arithmetic' if bounds['coefficient_arithmetic'] is None
            else None,
            'between-node remainder on I' if bounds['between_node_remainder'] is None else None,
            'complete source/error certificate',
        ],
    }
    scope['unresolved'] = [item for item in scope['unresolved'] if item is not None]
    return {
        'constraint_order': ('N', 'beta'),
        'z': np.array(z, copy=True),
        'action_gradient': residual,
        'history_jacobian': jacobian,
        'baseline_action_gradient': baseline_nodes,
        'local_reference_gradient_change': geometric['action_gradient_change'],
        'local_reference_gradient_tangent': geometric['action_gradient_tangent'],
        'incoming_slots': slots,
        'incoming_slot_tangents': slot_tangents,
        'family_corrections': corrections,
        'family_current_matter': currents,
        'family_reference_matter': references,
        'family_matter_tangents': tangents,
        'family_records': records,
        'error_budget': bounds,
        'interval': interval,
        'scope': MappingProxyType(scope),
        'stationarity_tolerance': STATIONARITY_TOLERANCE,
        'physical_constraint_status': 'OPEN',
        'physical_EXISTENCE_certificate': False,
    }


def assess_local_residual(assembly, interval=None, error_budget=None):
    """Max sampled |N,beta| plus declared bounds on a connected positive I.

    Endpoints of I must be sampled. None bounds cannot become zero. A
    combined value <= 3e-11 is only a conditional numeric tolerance, never a
    physical EXISTENCE certificate. Zero residual without complete bounds is
    not a PASS.
    """
    if not isinstance(assembly, Mapping) or 'action_gradient' not in assembly:
        raise TypeError('assemble_local_incoming result required')
    z = _finite(assembly['z'], 'incoming z', real=True)
    residual = _finite(assembly['action_gradient'], 'raw N,beta residuals', real=True)
    if residual.shape != (len(z), 2):
        raise ValueError('raw residuals must have shape (nz, 2) in N,beta order')
    interval = _positive_interval(
        interval if interval is not None else assembly.get('interval'))
    if error_budget is None:
        stored = assembly.get('error_budget')
        if stored is None:
            raise ValueError('error_budget required; missing bounds cannot become zero')
        bounds = local_error_budget(**{
            name: stored[name] for name in ERROR_BUDGET_COMPONENTS
        }) if isinstance(stored, Mapping) else stored
        if not isinstance(bounds, Mapping):
            raise TypeError('error_budget mapping required')
        bounds = local_error_budget(**{name: bounds[name] for name in ERROR_BUDGET_COMPONENTS})
    else:
        bounds, _ = _require_error_budget(error_budget)

    lo, hi = interval
    left = np.abs(z - lo) <= 3e-12
    right = np.abs(z - hi) <= 3e-12
    if not np.any(left) or not np.any(right):
        raise ValueError(
            'sample coverage of both interval endpoints required; missing coverage '
            'cannot pass a local residual')
    on_interval = (z >= lo - 3e-12) & (z <= hi + 3e-12)
    sampled = residual[on_interval]
    sampled_max = np.max(np.abs(sampled), axis=0)
    complete = all(bounds[name] is not None for name in ERROR_BUDGET_COMPONENTS)
    combined = None
    within = False
    if complete:
        combined = np.array(sampled_max, dtype=float, copy=True)
        for name in ERROR_BUDGET_COMPONENTS:
            combined = combined + bounds[name]
        within = bool(np.all(combined <= STATIONARITY_TOLERANCE))
    return {
        'constraint_order': ('N', 'beta'),
        'interval': list(interval),
        'sample_count_on_interval': int(np.count_nonzero(on_interval)),
        'endpoint_coverage': True,
        'sampled_maximum_absolute_residual': sampled_max,
        'declared_component_bounds': bounds,
        'declared_bounds_complete': complete,
        'combined_residual_bound': combined,
        'conditional_numeric_tolerance': within,
        'stationarity_tolerance': STATIONARITY_TOLERANCE,
        'physical_EXISTENCE_certificate': False,
        'physical_constraint_status': 'OPEN',
        'zero_residual_without_bounds_is_not_pass': True,
        'missing_bounds_are_not_zero': not complete,
    }
