"""KS-envelope local incoming N,beta assembly.

E_B[g] = E_B[g_ref] + Delta E_B^{local+reference}[j_Sigma g]
         + sum_families (G_B[C_Sigma[g]] - G_B[C_ref]),  B = N, beta.

Families share target z/history; each may have its own numerical mesh and
source digest. Reference: homogeneous ks_generator ODE of A_up, E, m, ell,
rho_up to rho1, then Fref=Aref exp(-iEz), Fref_z=-iE Fref, zero tangent,
cached by full prep and solver settings. No C0/HistoryBinding/node_fields.
Physical status stays OPEN. Baseline and geometry are counted once.
"""
from collections.abc import Mapping, Sequence
from types import MappingProxyType

import numpy as np
from scipy.integrate import solve_ivp

from .nsc_compatible_history_geometry import CompatibleIncomingMetric
from .nsc_evolved_incoming_constraints import (
    compatible_history_slots, source_column_matter, surface_geometry_response,
)
from .nsc_evolved_incoming_state import (
    AmplitudeOnlyMetric, CachedCompatibleIncomingMetric, FixedSourcePreparation,
    HistoryBinding,
)
from .nsc_ks_source_envelope import KSEnvelopeIncoming, RHO_SIGMA, _fixed_preparation_digest
from .nsc_ks_difference_envelope import KSDifferenceIncoming
from .nsc_ks_matter_difference import source_fixed_matter_difference
from .nsc_local_incoming_constraints import (
    _channel_record, _coverage_record, _finite, _positive_interval,
    _require_error_budget,
)
from .nsc_local_incoming_family import LocalIncomingFamily
from .nsc_pg_ks_metric_pullback import ks_to_pg, reference_chart
from .nsc_retarded_radial_response import ks_generator
from .nsc_transmitting_dirac_domain import TransmittingDiracSeamDomain


_REFERENCE_CACHE = {}
_LEGACY = frozenset({
    'stress_approximant', 'parts', 'incoming_C0', 'C0', 'node_fields',
    'history', 'HistoryBinding', 'covariance', 'matter',
})
_MISSING_BOUNDS = (
    ('baseline_low_subgap', 'baseline low-energy/subgap source accuracy'),
    ('changed_history_finite_energy', 'changed-history finite-energy coverage'),
    ('changed_history_tail', 'changed-history tail'),
    ('preparation_field_axial_derivative', 'preparation/field/axial-derivative error'),
    ('coefficient_arithmetic', 'coefficient interval and arithmetic'),
    ('between_node_remainder', 'between-node remainder on I'),
)


def _history_family(provider):
    family = provider.inner if isinstance(provider, AmplitudeOnlyMetric) else provider
    if isinstance(family, CachedCompatibleIncomingMetric):
        family = family.original
    if isinstance(family, LocalIncomingFamily):
        family = family.metric()
    if isinstance(family, HistoryBinding) or not isinstance(family, CompatibleIncomingMetric):
        raise TypeError(
            'LocalIncomingFamily or CompatibleIncomingMetric required; HistoryBinding and PG node_fields are not this assembly')
    return family


def _as_family_map(prepared_families):
    if prepared_families is None:
        raise ValueError('prepared_families required')
    if isinstance(prepared_families, Mapping):
        items = list(prepared_families.items())
    elif isinstance(prepared_families, Sequence) and not isinstance(prepared_families, (str, bytes)):
        items = [(entry[0], entry[1]) for entry in prepared_families]
    else:
        raise TypeError('prepared families must be a mapping or a sequence of (id, KSEnvelopeIncoming)')
    ids = [entry[0] for entry in items]
    if len(ids) != len(set(ids)) or not ids:
        raise ValueError(
            'duplicate family IDs' if ids else
            'at least one prepared family required; omitted families remain explicit OPEN')
    families = {}
    for family_id, prepared in items:
        if isinstance(prepared, HistoryBinding) or (
                isinstance(prepared, dict) and _LEGACY.intersection(prepared)):
            raise ValueError(
                'a frozen incoming C0, PG HistoryBinding, or node_fields dictionary '
                'is not a KS envelope family')
        if not isinstance(prepared, KSEnvelopeIncoming):
            raise TypeError(
                'KSEnvelopeIncoming required for each family; frozen C0, PG '
                'HistoryBinding, or node_fields are not this assembly')
        prepared.validate()
        if float(prepared.binding.rho_sigma) != RHO_SIGMA:
            raise ValueError('incoming surface is the owned rho_sigma=1 slice')
        families[family_id] = prepared
    return families


def _family_source(prepared, channel):
    source = FixedSourcePreparation(
        prepared.source_covariance, prepared.column_weights, prepared.source_energies)
    digest = _fixed_preparation_digest(
        source, prepared.initial_columns, channel['compact_mass'],
        channel['angular_sign'] * channel['angular_eigenvalue'],
        prepared.rho_up, prepared.binding.rho_sigma)
    if digest != prepared.fixed_preparation_digest:
        raise ValueError(
            'channel mass/sign does not match the bound preparation digest rebuilt from A_up and the actual source')
    expected = channel.get('expected_source_digest')
    if expected is not None and source.digest != expected:
        raise ValueError('prepared source differs from its declared family source')
    return source, digest


def _incoming_intrinsic():
    _, _, axial, _, _ = reference_chart(RHO_SIGMA)
    N, beta, q, r = ks_to_pg(RHO_SIGMA, 1., 0., axial, np.sqrt(2.))
    domain = TransmittingDiracSeamDomain(
        lapse=float(N), radial_scale=float(q), shift=float(beta), radius=float(r))
    return float(domain.induced_axial_scale), float(r)


def homogeneous_reference_amplitudes(prepared):
    """Endpoint Aref of the homogeneous ks_generator radial ODE, cached by prep."""
    if not isinstance(prepared, KSEnvelopeIncoming):
        raise TypeError('KSEnvelopeIncoming required')
    prepared.validate()
    key = (prepared.fixed_preparation_digest, float(prepared.binding.rtol),
           float(prepared.binding.atol), float(prepared.binding.max_step))
    cached = _REFERENCE_CACHE.get(key)
    if cached is not None:
        return cached
    initial = np.asarray(prepared.initial_columns, complex)
    energies = np.asarray(prepared.source_energies, float)

    def rhs(rho, state):
        return np.einsum(
            'sab,bs->as', ks_generator(rho, energies, prepared.mass, prepared.angular),
            state.reshape(initial.shape)).ravel()

    solution = solve_ivp(
        rhs, (float(prepared.rho_up), float(prepared.binding.rho_sigma)),
        initial.ravel(), method='DOP853', rtol=float(prepared.binding.rtol),
        atol=float(prepared.binding.atol), max_step=float(prepared.binding.max_step))
    if not solution.success:
        raise ArithmeticError(solution.message)
    Aref = np.array(solution.y[:, -1].reshape(initial.shape), complex, copy=True)
    if not np.isfinite(Aref).all():
        raise ArithmeticError('finite homogeneous KS reference amplitudes required')
    Aref.setflags(write=False)
    _REFERENCE_CACHE[key] = Aref
    return Aref


def stationary_computational_reference(prepared, z=None):
    """Fref = Aref exp(-i E z), Fref_z = -i E Fref. Tangent is identically zero."""
    Aref = homogeneous_reference_amplitudes(prepared)
    nodes = np.asarray(prepared.z if z is None else z, float)
    energies = np.asarray(prepared.source_energies, float)
    phase = np.exp(-1j * energies * nodes[:, None])
    reference = Aref[None, :, :] * phase[:, None, :]
    return reference, -1j * energies * reference


def _family_matter(prepared, channel, coefficients, ndir, z):
    source, digest = _family_source(prepared, channel)
    F, Fz = np.asarray(prepared.weighted_columns), np.asarray(prepared.weighted_axial_columns)
    dF, dFz = np.asarray(prepared.weighted_column_tangents), np.asarray(prepared.weighted_axial_tangents)
    if dF.shape[0] != ndir or dFz.shape[0] != ndir:
        raise ValueError('prepared direction count must match the assembled history')
    if not np.array_equal(np.asarray(prepared.z), z):
        raise ValueError('mismatched incoming z across prepared families')

    reference, reference_z = stationary_computational_reference(prepared, z)
    zeros = np.zeros((ndir, *reference.shape), complex)
    parameters = dict(
        mass=channel['compact_mass'],
        angular=channel['angular_sign'] * channel['angular_eigenvalue'],
        axial_scale=float(coefficients['a']), radius=float(coefficients['r']),
        multiplicity=channel['multiplicity'])
    reference_matter = source_column_matter(
        reference * prepared.column_weights, reference_z * prepared.column_weights,
        prepared.source_covariance, zeros, zeros, **parameters)
    if isinstance(prepared, KSDifferenceIncoming):
        stable = source_fixed_matter_difference(prepared, axial_scale=parameters['axial_scale'],
            radius=parameters['radius'], multiplicity=parameters['multiplicity'], reference='cached')
        current_value, tangent = stable['current_action_gradient'], stable['action_gradient_tangent']
        correction = stable['action_gradient_change']
    else:
        current = source_column_matter(F, Fz, prepared.source_covariance, dF, dFz, **parameters)
        current_value, tangent = current['action_gradient'], current['action_gradient_tangent']
        correction = current_value - reference_matter['action_gradient']
    return {
        'current': current_value,
        'reference': reference_matter['action_gradient'],
        'correction': correction,
        'tangent': tangent,
        'multiplicity': parameters['multiplicity'],
        'signed_angular': parameters['angular'],
        'preparation_digest': digest,
        'source_digest': source.digest,
        'ell0_pure_radius_identity': channel['ell0_pure_radius_identity'],
        'computational_count': int(len(prepared.binding.computational_z)),
        'stable_difference': isinstance(prepared, KSDifferenceIncoming),
    }


def assemble_local_incoming(
        prepared_families, provider, baseline_gradient, coefficients,
        channels, interval, error_budget, coverage=None):
    """Raw N,beta residuals and full retarded derivative from KSEnvelopeIncoming.

    Subclasses are accepted by isinstance. Baseline and geometry enter once.
    Physical constraint status is always OPEN.
    """
    families = _as_family_map(prepared_families)
    history = _history_family(provider)
    interval = _positive_interval(interval)
    bounds, budget_extra = _require_error_budget(error_budget)
    extra = {k: budget_extra[k] for k in (
        'sampled_family_ids', 'missing_family_ids', 'sampled_energy_ranges',
        'missing_energy_ranges', 'energy_coverage') if k in budget_extra}
    coverage_record = _coverage_record(families, coverage if coverage is not None else extra or None)
    first = next(iter(families.values()))
    z, ndir = np.asarray(first.z), int(first.column_tangents.shape[0])
    first.require_history(history)
    axial, radius = _incoming_intrinsic()
    if not np.allclose([coefficients['a'], coefficients['r']], [axial, radius], atol=3e-12, rtol=0):
        raise ValueError(
            'constraint coefficients and prepared states must use the same intrinsic identification')

    seen, records, corrections, currents, references, tangents = {}, {}, {}, {}, {}, {}
    for family_id, prepared in families.items():
        prepared.require_history(history)
        channel = _channel_record(family_id, channels)
        key = (channel['group'], channel['angular_sign'])
        if key in seen:
            raise ValueError('duplicate family IDs for the same group and angular_sign')
        seen[key] = family_id
        matter = _family_matter(prepared, channel, coefficients, ndir, z)
        ell0 = matter['ell0_pure_radius_identity']
        records[family_id] = {
            'group': channel['group'], 'angular_sign': channel['angular_sign'],
            'compact_mass': channel['compact_mass'],
            'angular_eigenvalue': channel['angular_eigenvalue'],
            'signed_angular': matter['signed_angular'],
            'copy_count': channel['copy_count'], 'degeneracy': channel['degeneracy'],
            'n_actual_angular_signs': channel['n_actual_angular_signs'],
            'multiplicity': matter['multiplicity'],
            'preparation_digest': matter['preparation_digest'],
            'source_digest': matter['source_digest'],
            'computational_count': matter['computational_count'],
            'stable_difference': matter['stable_difference'],
            'ell0_pure_radius_identity': ell0,
            'ell0_documentation': (
                'included; pure-radius response identically zero by ell=0 Dirac '
                'r-independence; baseline retained; the displayed discrete '
                'current-minus-analytic-reference also includes numerical drift' if ell0 else None),
        }
        corrections[family_id] = matter['correction']
        currents[family_id] = matter['current']
        references[family_id] = matter['reference']
        tangents[family_id] = matter['tangent']

    slots, slot_tangents = compatible_history_slots(history, z, ndir)
    geometric = surface_geometry_response(slots, slot_tangents, coefficients)
    baseline = _finite(baseline_gradient, 'baseline complete-assembly approximant', real=True)
    if baseline.shape == (2,):
        baseline_nodes = np.broadcast_to(baseline, (len(z), 2)).copy()
    elif baseline.shape == (len(z), 2):
        baseline_nodes = np.array(baseline, dtype=float, copy=True)
    else:
        raise ValueError('baseline N,beta assembly approximant required with shape (2,) or (nz, 2)')
    residual = baseline_nodes + geometric['action_gradient_change'] + sum(corrections.values())
    jacobian = geometric['action_gradient_tangent'] + sum(tangents.values())
    unresolved = [label for name, label in _MISSING_BOUNDS if bounds[name] is None]
    unresolved.append('complete source/error certificate')
    scope = MappingProxyType({
        'interval': list(interval), 'production_local_gate': 'I=S(1)+[.12,.18]',
        'global_matching': False, 'extended_stationarity': 'OPEN',
        'metric_timestep': False, 'observations': False,
        'physical_constraint_status': 'OPEN', 'physical_EXISTENCE_certificate': False,
        'certified_constraint_residual': None, 'full_source_error_bound': None,
        'baseline_and_geometry_counted_once': True, 'reference_tangent': 0,
        'reference_state_tangent_subtracted': False,
        'constant_frozen_matter_used_as_physical_state': False,
        'pde_axial_derivatives': True, 'interpolated_axial_derivatives': False,
        'homogeneous_ks_generator_reference': True,
        'reference_cached_by_preparation_and_solver_settings': True,
        'computational_meshes_need_not_match': True,
        'finite_energy_samples_are_complete_spectrum': False,
        'independent_source_error_certificate': budget_extra.get(
            'independent_source_error_certificate'),
        **coverage_record,
        'ell0_pure_radius_identities': {
            family_id: record['ell0_documentation']
            for family_id, record in records.items() if record['ell0_pure_radius_identity']
        },
        'unresolved': unresolved,
    })
    return {
        'constraint_order': ('N', 'beta'), 'z': np.array(z, copy=True),
        'action_gradient': residual, 'history_jacobian': jacobian,
        'baseline_action_gradient': baseline_nodes,
        'local_reference_gradient_change': geometric['action_gradient_change'],
        'local_reference_gradient_tangent': geometric['action_gradient_tangent'],
        'incoming_slots': slots, 'incoming_slot_tangents': slot_tangents,
        'family_corrections': corrections, 'family_current_matter': currents,
        'family_reference_matter': references, 'family_matter_tangents': tangents,
        'family_records': records, 'error_budget': bounds, 'interval': interval,
        'scope': scope, 'physical_constraint_status': 'OPEN',
        'physical_EXISTENCE_certificate': False,
    }
