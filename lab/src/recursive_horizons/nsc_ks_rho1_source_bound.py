"""Physical rho=1 source-column coverage and bound skeleton.

This owner traces the authenticated V5 inventory, the retained source
panels, and the existing upstream/low-subgap classifiers. It preserves
coherent three-column fibers, original weights, signed multiplicities,
degeneracies, E/-E pairing and the analytic ell=0 zero radius-history
response while retaining those groups' baseline. The 102 panels / 29356
rows are a panel-level low/subgap campaign subset, not a complete physical
column-accuracy universe: 59 panels straddle their V5 joined threshold and
V5 action bounds do not certify source-column reconstruction. This v1 core
never fills a finite physical N,beta budget.
"""
from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

import numpy as np

from .nsc_ks_residual_error import nonnegative
from .nsc_ks_source_inventory import ReferenceSourcePanel
from .nsc_ks_upstream_error import (
    ELL0_GROUPS,
    IDENTITIES_NOT_ERROR_BOUNDS,
    classify_low_subgap_windows,
    combine_upstream_row_error,
    family14_signed_multiplicity,
    nine_component_upstream_feed,
    pack_upper,
    restore_packed_upper,
)
from .nsc_transmitting_dirac_domain import S3


SCHEMA = 'NSC-KS-RHO1-SOURCE-BOUND-v1'
INVENTORY_POSITIVE_PANELS = 164
INVENTORY_POSITIVE_ROWS = 49372
V5_ACTION_COVERED_PANELS = 62
V5_ACTION_COVERED_ROWS = 20016
REMAINING_RECONSTRUCTION_PANELS = 102
REMAINING_RECONSTRUCTION_ROWS = 29356
BELOW_JOINED_ROWS_IN_REMAINING_PANELS = 27828
AT_OR_ABOVE_JOINED_ROWS_IN_REMAINING_PANELS = 1528
STRADDLING_REMAINING_PANELS = 59
HISTORICAL_COARSE_AGGREGATE = 4.998e-9
HERMITIAN_TOLERANCE = 3e-13
SAMPLE_ROW_LIMIT = 3
FORBIDDEN_BOUND_METHODS = (
    'independent_column_max',
    'drop_coherence',
    'drop_minus_energy',
    'drop_weights',
    'drop_degeneracy',
    'historical_coarse_4.998e-9',
    'pauli_identity_as_bound',
    'nested_node_as_bound',
)
FINE_PREFIXES = (('low_ref/', 'low/'), ('mid_ref/', 'mid/'))


def coverage_identity(group, sign, panel, row_index=None):
    """Exact family/channel/panel/row key; row_index None is the whole panel."""
    if isinstance(group, bool) or int(group) != group or int(group) not in range(0, 33):
        raise ValueError('retained group 0 through 32 required')
    if sign not in (-1, 1):
        raise ValueError('actual angular sign required')
    if panel is None:
        if int(group) != 0:
            raise ValueError('only the analytic LLL group may omit a mode panel')
        key = (0, 1, None)
    else:
        if not isinstance(panel, str) or not panel:
            raise ValueError('named source panel required')
        key = (int(group), int(sign), str(panel))
    if row_index is None:
        return key
    if isinstance(row_index, bool) or int(row_index) != row_index or int(row_index) < 0:
        raise ValueError('nonnegative row index required')
    return key + (int(row_index),)


def row_identities(group, sign, panel, n_rows):
    if isinstance(n_rows, bool) or int(n_rows) != n_rows or int(n_rows) < 1:
        raise ValueError('positive selected row count required')
    return tuple(coverage_identity(group, sign, panel, index) for index in range(int(n_rows)))


def signed_channel_multiplicity(channel):
    """copies*degeneracy split across actual angular signs, never across E/-E."""
    return family14_signed_multiplicity(channel)


def radius_history_response(group):
    return 'exact_zero' if int(group) in ELL0_GROUPS else 'nontrivial'


def serialize_rational(value):
    number = nonnegative(value, 'rational')
    return {'numerator': str(number.numerator), 'denominator': str(number.denominator)}


def restore_rational(item):
    if not isinstance(item, dict) or set(item) != {'numerator', 'denominator'}:
        raise ValueError('rational numerator/denominator required')
    return nonnegative(Fraction(int(item['numerator']), int(item['denominator'])),
                       'rational')


def serialize_directed_upper(value):
    """Directed dyadic upper of a nonnegative bound; None stays missing."""
    if value is None:
        return None
    from flint import arb, ctx
    number = nonnegative(value, 'directed upper')
    with ctx.workprec(80):
        return pack_upper(arb(number.numerator) / arb(number.denominator))


def restore_directed_upper(item):
    if item is None:
        return None
    return restore_packed_upper(item)


def coverage_digest(identities):
    payload = json.dumps(list(identities), sort_keys=True, separators=(',', ':')).encode()
    return sha256(payload).hexdigest()


def _coarse_counterpart(panel):
    for fine, coarse in FINE_PREFIXES:
        if panel.startswith(fine):
            return coarse + panel[len(fine):]
    return None


def reject_double_counting(panel_names):
    names = tuple(panel_names)
    if any(name is None or not isinstance(name, str) or not name for name in names):
        raise ValueError('named source panels required for reconstruction coverage')
    if len(names) != len(set(names)):
        raise ValueError('duplicate rho=1 source panel')
    present = set(names)
    for name in names:
        coarse = _coarse_counterpart(name)
        if coarse is not None and coarse in present:
            raise ValueError('fine panel must replace, never add, its coarse counterpart')
    return names


def reject_forbidden_bound_method(method):
    if method in FORBIDDEN_BOUND_METHODS or method == 'historical_coarse_4.998e-9':
        raise ValueError('forbidden bound method: ' + str(method))
    if method is None or not isinstance(method, str) or not method:
        raise ValueError('explicit reconstruction method required')
    return method


def reject_historical_coarse(bound):
    if bound is None:
        return None
    number = nonnegative(bound, 'physical source bound')
    if abs(float(number) - HISTORICAL_COARSE_AGGREGATE) <= 5e-13:
        raise ValueError('historical 4.998e-9 coarse aggregate is not the current remaining bound')
    return number


def coherent_covariance_required(covariance):
    blocks = np.asarray(covariance, complex)
    if blocks.ndim == 2:
        blocks = blocks[None]
    if blocks.ndim != 3 or blocks.shape[-2:] != (3, 3) or len(blocks) == 0:
        raise ValueError('three-column coherent source covariance required')
    if not np.isfinite(blocks).all():
        raise ValueError('finite source covariance required')
    hermitian = float(np.max(np.abs(blocks - np.swapaxes(blocks, -1, -2).conj())))
    if hermitian > HERMITIAN_TOLERANCE:
        raise ValueError('Hermitian source covariance required')
    eig = np.linalg.eigvalsh(blocks)
    if eig.min() < -HERMITIAN_TOLERANCE or eig.max() > 1 + HERMITIAN_TOLERANCE:
        raise ValueError('source CAR interval violated')
    return blocks


def covariance_preservation_report(covariance):
    blocks = coherent_covariance_required(covariance)
    diagonal = np.zeros_like(blocks)
    for index in range(3):
        diagonal[..., index, index] = blocks[..., index, index]
    return {
        'hermitian_residual': float(np.max(np.abs(blocks - np.swapaxes(blocks, -1, -2).conj()))),
        'offdiag_max': float(np.max(np.abs(blocks - diagonal))),
        'coherence_retained': True,
        'independent_column_max_used': False,
    }


def assert_coherence_retained(covariance, *, require_offdiag=False):
    report = covariance_preservation_report(covariance)
    if require_offdiag and report['offdiag_max'] == 0:
        raise ValueError('coherent off-diagonal source blocks required')
    return report


def pair_negative_energy(panel, config):
    """Owned antiunitary E/-E map; same weights, no signed-energy folding."""
    if not isinstance(panel, ReferenceSourcePanel):
        raise TypeError('ReferenceSourcePanel required')
    partner = panel.negative_partner(config)
    if partner.provenance.get('energy_folding_factor') != 1:
        raise ValueError('E/-E pairing must not fold signed energy')
    if not np.array_equal(partner.weights, panel.weights):
        raise ValueError('negative partner must retain original weights')
    if np.max(np.abs(partner.energies + panel.energies)) != 0:
        raise ValueError('negative partner must negate the original energies')
    residual = float(np.max(np.abs(partner.covariance - (np.eye(3) - panel.covariance.conj()))))
    if residual > HERMITIAN_TOLERANCE:
        raise ValueError('negative partner must use the owned source complement')
    if panel.angular_magnitude and partner.angular_sign != -panel.angular_sign:
        raise ValueError('nonzero angular channels pair with opposite sign')
    if not panel.angular_magnitude and partner.angular_sign != 1:
        raise ValueError('ell=0 negative partner keeps the unique plus sign')
    mapped = S3 @ panel.amplitudes_at_one.conj()
    if float(np.max(np.abs(partner.amplitudes_at_one - mapped))) > 0:
        raise ValueError('negative partner amplitudes must be S3 conjugate')
    return partner


def assert_ell0_baseline_retained(group, *, baseline_retained, panel=None,
                                  radius_response=None):
    group = int(group)
    response = radius_history_response(group) if radius_response is None else radius_response
    if group in ELL0_GROUPS:
        if response != 'exact_zero':
            raise ValueError('ell=0 radius-history response is exactly zero')
        if not baseline_retained:
            raise ValueError('ell=0 exact zero response retains the nonzero baseline')
        if group == 0 and panel is not None:
            raise ValueError('group 0 has no synthetic mode panel')
    elif response == 'exact_zero':
        raise ValueError('non-ell=0 groups have a nontrivial radius-history response')
    return response


def attach_reconstruction_bound(entry, bound, method):
    """Attach an explicit directed bound, or keep None. Never invent a budget."""
    reject_forbidden_bound_method(method)
    if bound is None:
        updated = dict(entry)
        updated['reconstruction_error_upper'] = None
        updated['method'] = method
        return updated
    reject_historical_coarse(bound)
    updated = dict(entry)
    updated['reconstruction_error_upper'] = nonnegative(bound, 'reconstruction error')
    updated['method'] = method
    return updated


def combine_directed_n_beta(contributions, *, declared_universe_size):
    """Validate provisional rows while keeping v1 incapable of closure.

    A count alone cannot authenticate the 164-panel / 49372-row physical
    column universe, and zero-filled rows cannot prove reconstruction accuracy.
    A successor must bind every contribution to a row manifest and enclosure
    payload before summation. This v1 API therefore always returns OPEN.
    """
    if (isinstance(declared_universe_size, bool)
            or int(declared_universe_size) != declared_universe_size
            or int(declared_universe_size) < 1):
        raise ValueError('positive declared universe size required')
    universe = int(declared_universe_size)
    rows = list(contributions)
    if len(rows) != universe:
        raise ValueError('contribution count does not match the declared universe')
    identities = []
    for row in rows:
        reject_forbidden_bound_method(row.get('method', ''))
        if not row.get('coherence_retained') or not row.get('weights_retained'):
            raise ValueError('covariance, coherence and weights must be retained')
        if not row.get('minus_energy_retained'):
            raise ValueError('minus-energy sector must be retained')
        if not row.get('degeneracy_retained', True):
            raise ValueError('channel degeneracy must be retained')
        identities.append(tuple(row['identity']))
    if len(identities) != len(set(identities)):
        raise ValueError('double-counted rho=1 source contribution')
    finite_inputs = all(
        row.get('N') is not None and row.get('beta') is not None for row in rows)
    for row in rows:
        if row.get('N') is not None:
            nonnegative(row['N'], 'N contribution')
        if row.get('beta') is not None:
            nonnegative(row['beta'], 'beta contribution')
    return {
        'N': None,
        'beta': None,
        'status': (
            'OPEN: v1 has no authenticated row-manifest and enclosure-payload '
            'closure for the complete physical source universe'),
        'declared_universe_size': universe,
        'independent_column_max_used': False,
        'finite_inputs_do_not_close_v1': finite_inputs,
    }


def physical_source_error_for_upstream(physical_n_beta, continuation_n_beta=None):
    """Feed the existing upstream owner. Missing physical input stays OPEN."""
    if physical_n_beta is None:
        combined = combine_upstream_row_error(
            None if continuation_n_beta is None else continuation_n_beta[0], None)
        feed = nine_component_upstream_feed(continuation_n_beta, None)
        return {'combined_row_error': combined, 'nine_component_feed': feed}
    raise ValueError(
        'v1 has no authenticated physical rho=1 column-reconstruction bound')


def classify_rho1_windows(inventory_record, panels, v5_coverage):
    """Reuse the existing low/subgap classifier; reconstruction stays unbounded."""
    return classify_low_subgap_windows(inventory_record, panels, v5_coverage)


def _joined(coverage, group):
    row = coverage[str(int(group))]
    return list(row['joined'])


def build_rho1_coverage(windows, *, inventory_record, v5_coverage, channels=None):
    """Exact panel/row coverage from the existing classified windows."""
    ell0 = tuple(inventory_record['radius_response_exactly_zero_groups'])
    if ell0 != ELL0_GROUPS and ell0 != list(ELL0_GROUPS):
        raise ValueError('ell=0 radius-response groups changed')
    reconstruction = []
    remaining = []
    covered_action = []
    lll = 0
    for window in windows:
        group = int(window['group'])
        sign = window.get('sign')
        panel = window.get('panel')
        response = window.get('radius_history_response')
        baseline = window.get('baseline_action_region')
        assert_ell0_baseline_retained(
            group, baseline_retained=True, panel=panel, radius_response=response)
        if group == 0:
            if panel is not None or window.get('kind') != 'lll_analytic':
                raise ValueError('group 0 must remain the analytic LLL owner')
            lll += 1
            continue
        if panel is None:
            if window.get('kind') != 'covered_baseline_action':
                raise ValueError('mode-column windows require a panel name')
            continue
        if sign not in (-1, 1):
            raise ValueError('actual signed reconstruction panel required')
        if group in ELL0_GROUPS and sign != 1:
            raise ValueError('ell=0 channels have a unique plus sign')
        if channels is not None:
            channel = channels[group]
            expected = (1,) if float(channel['angular_eigenvalue']) == 0 else (-1, 1)
            if sign not in expected:
                raise ValueError('panel sign is outside the retained channel inventory')
        if window.get('v5_joined') != _joined(v5_coverage, group):
            raise ValueError('panel V5 joined threshold does not match the V5 owner')
        n_rows = int(window['n_rows'])
        if n_rows < 1:
            raise ValueError('positive selected row count required')
        item = {
            'group': group,
            'sign': int(sign),
            'panel': panel,
            'n_rows': n_rows,
            'kind': window['kind'],
            'energy_interval': list(window['energy_interval']),
            'v5_joined': list(window['v5_joined']),
            'radius_history_response': response,
            'baseline_action_region': baseline,
            'source_reconstruction': window.get('source_reconstruction'),
            'error_bound': window.get('error_bound'),
            'baseline_retained': True,
        }
        reconstruction.append(item)
        if baseline == 'not_covered':
            remaining.append(item)
        elif baseline == 'covered_v5_joined':
            covered_action.append(item)
        else:
            raise ValueError('reconstruction panel must be V5-covered or remaining')
        if item['source_reconstruction'] != 'unbounded' or item['error_bound'] is not None:
            raise ValueError('V5 action coverage is not a rho=1 column reconstruction bound')
    if lll != 1:
        raise ValueError('analytic LLL owner must appear once')
    reject_double_counting(item['panel'] for item in reconstruction)
    reject_double_counting(item['panel'] for item in remaining)
    remaining_rows = sum(item['n_rows'] for item in remaining)
    covered_rows = sum(item['n_rows'] for item in covered_action)
    recon_rows = sum(item['n_rows'] for item in reconstruction)
    if (len(reconstruction) != INVENTORY_POSITIVE_PANELS
            or recon_rows != INVENTORY_POSITIVE_ROWS):
        raise ValueError('source inventory panel/row coverage changed')
    if (len(remaining) != REMAINING_RECONSTRUCTION_PANELS
            or remaining_rows != REMAINING_RECONSTRUCTION_ROWS):
        raise ValueError('remaining rho=1 reconstruction coverage changed')
    if (len(covered_action) != V5_ACTION_COVERED_PANELS
            or covered_rows != V5_ACTION_COVERED_ROWS):
        raise ValueError('V5 action-covered panel/row coverage changed')
    if remaining_rows + covered_rows != recon_rows:
        raise ValueError('covered and remaining rows do not partition the inventory')
    remaining_identities = [
        list(coverage_identity(item['group'], item['sign'], item['panel']))
        + [item['n_rows']] for item in remaining]
    covered_identities = [
        list(coverage_identity(item['group'], item['sign'], item['panel']))
        + [item['n_rows']] for item in covered_action]
    by_kind = {}
    for item in remaining:
        by_kind[item['kind']] = by_kind.get(item['kind'], 0) + 1
    ell0_remaining = [item for item in remaining if item['group'] in ELL0_GROUPS]
    family_channels = sorted({(item['group'], item['sign']) for item in reconstruction})
    return {
        'inventory_positive_panels': len(reconstruction),
        'inventory_positive_rows': recon_rows,
        'v5_action_covered_panels': len(covered_action),
        'v5_action_covered_rows': covered_rows,
        'remaining_reconstruction_panels': len(remaining),
        'remaining_reconstruction_rows': remaining_rows,
        'remaining_by_kind': by_kind,
        'angular_families_with_columns': len(family_channels),
        'ell0_zero_response_groups': list(ELL0_GROUPS),
        'ell0_remaining_reconstruction_panels': len(ell0_remaining),
        'ell0_remaining_reconstruction_rows': sum(item['n_rows'] for item in ell0_remaining),
        'remaining_panel_identities': remaining_identities,
        'v5_action_panel_identities': covered_identities,
        'remaining_identity_digest': coverage_digest(remaining_identities),
        'v5_action_identity_digest': coverage_digest(covered_identities),
        'double_counted': False,
        'v5_action_is_not_column_reconstruction': True,
        'baseline_ell0_allocations_retained': True,
        'identities_that_are_not_error_bounds': list(IDENTITIES_NOT_ERROR_BOUNDS),
        'remaining_physical_rho1_source_error_bound': None,
        'remaining_low_subgap_source_error_bound': None,
        'directed_N_beta_contribution_bound': {
            'N': None, 'beta': None,
            'status': (
                'OPEN: complete 164-panel physical column reconstruction '
                'unevaluated; 102-panel low/subgap subset identified'),
        },
        'remaining_windows': remaining,
        'v5_action_windows': covered_action,
    }


def panel_from_arrays(name, meta_panel, arrays, *, rows=None):
    energies = np.asarray(arrays[name + '/energies'])
    weights = np.asarray(arrays[name + '/weights'])
    covariance = np.asarray(arrays[name + '/covariance'])
    amplitudes = np.asarray(arrays[name + '/amplitudes_at_one'])
    if rows is not None:
        rows = np.asarray(rows, int)
        energies, weights = energies[rows], weights[rows]
        covariance, amplitudes = covariance[rows], amplitudes[rows]
    return ReferenceSourcePanel(
        name, int(meta_panel['group']), int(meta_panel['angular_sign']),
        float(meta_panel['mass']), float(meta_panel['angular_magnitude']),
        energies, weights, covariance, amplitudes, meta_panel.get('provenance', {}))


def load_authenticated_inventory_payload(root, inventory_record):
    spec = inventory_record['payload']
    path = Path(root) / spec['path']
    digest = sha256(path.read_bytes()).hexdigest()
    if digest != spec['sha256']:
        raise ValueError('source inventory payload changed: ' + spec['path'])
    with np.load(path, allow_pickle=False) as loaded:
        arrays = {key: np.array(loaded[key], copy=True) for key in loaded.files}
    meta = json.loads(arrays['metadata_json'].tobytes())
    return arrays, meta, digest


def audit_row_threshold_partition(arrays, coverage):
    """Count authenticated rows relative to each family V5 joined threshold."""
    if not isinstance(coverage, dict):
        raise TypeError('rho=1 coverage mapping required')
    rows = []
    totals = {'below_joined': 0, 'at_or_above_joined': 0, 'straddling_panels': 0}
    for kind in ('remaining_windows', 'v5_action_windows'):
        for window in coverage[kind]:
            energies = np.asarray(arrays[window['panel'] + '/energies'], float)
            if (energies.ndim != 1 or len(energies) != int(window['n_rows'])
                    or not np.isfinite(energies).all()):
                raise ValueError('authenticated source-panel energy rows changed')
            joined = float(window['v5_joined'][0])
            below = int(np.count_nonzero(energies < joined))
            at_or_above = int(len(energies) - below)
            if kind == 'v5_action_windows' and below:
                raise ValueError('V5 action panel contains rows below its joined threshold')
            if kind == 'remaining_windows':
                totals['below_joined'] += below
                totals['at_or_above_joined'] += at_or_above
                totals['straddling_panels'] += int(below > 0 and at_or_above > 0)
            rows.append([
                int(window['group']), int(window['sign']), window['panel'],
                int(len(energies)), below, at_or_above, float(joined).hex(), kind,
            ])
    if (totals['below_joined'] != BELOW_JOINED_ROWS_IN_REMAINING_PANELS
            or totals['at_or_above_joined'] !=
            AT_OR_ABOVE_JOINED_ROWS_IN_REMAINING_PANELS
            or totals['straddling_panels'] != STRADDLING_REMAINING_PANELS):
        raise ValueError('row-level V5 threshold partition changed')
    return {
        'remaining_panel_campaign_panels': REMAINING_RECONSTRUCTION_PANELS,
        'remaining_panel_campaign_rows': REMAINING_RECONSTRUCTION_ROWS,
        'rows_below_joined_in_remaining_panels': totals['below_joined'],
        'rows_at_or_above_joined_in_remaining_panels': totals['at_or_above_joined'],
        'straddling_remaining_panels': totals['straddling_panels'],
        'physical_column_accuracy_universe_panels': INVENTORY_POSITIVE_PANELS,
        'physical_column_accuracy_universe_rows': INVENTORY_POSITIVE_ROWS,
        'v5_action_is_not_column_reconstruction': True,
        'remaining_panel_campaign_is_not_complete_physical_source_universe': True,
        'row_partition_digest': coverage_digest(rows),
    }


def inspect_authenticated_sample(arrays, meta, remaining_windows, config, *,
                                 panel_names=None, rows_per_panel=2):
    """Load a few remaining fibers. This is not a 102-panel reconstruction."""
    if rows_per_panel < 1 or rows_per_panel > SAMPLE_ROW_LIMIT:
        raise ValueError('sample must stay inside the declared small-row limit')
    remaining = {item['panel']: item for item in remaining_windows}
    names = list(panel_names) if panel_names is not None else [
        'group14/low32_1', 'group13/low']
    reports = []
    for name in names:
        if name not in remaining:
            raise ValueError('sample panel is not a remaining reconstruction window: ' + name)
        window = remaining[name]
        n_take = min(int(rows_per_panel), int(window['n_rows']))
        if name == 'group13/low':
            n_take = min(1, n_take)
        indices = tuple(range(n_take))
        panel = panel_from_arrays(name, meta['panels'][name], arrays, rows=indices)
        if len(panel.energies) != n_take or n_take < 1:
            raise ValueError('sample rows missing from authenticated payload')
        if int(panel.group) != int(window['group']) or int(panel.angular_sign) != int(window['sign']):
            raise ValueError('sample panel family/channel mismatch')
        report = covariance_preservation_report(panel.covariance)
        partner = pair_negative_energy(panel, config)
        response = assert_ell0_baseline_retained(
            panel.group, baseline_retained=True, panel=name,
            radius_response=window['radius_history_response'])
        reports.append({
            'panel': name,
            'group': int(panel.group),
            'sign': int(panel.angular_sign),
            'row_indices': list(indices),
            'n_rows': int(n_take),
            'panel_rows': int(window['n_rows']),
            'energies': [float(value) for value in panel.energies],
            'weights': [float(value) for value in panel.weights],
            'hermitian_residual': report['hermitian_residual'],
            'offdiag_max': report['offdiag_max'],
            'coherence_retained': True,
            'minus_energy_residual': float(np.max(np.abs(
                partner.covariance - (np.eye(3) - panel.covariance.conj())))),
            'energy_folding_factor': 1,
            'radius_history_response': response,
            'baseline_retained': True,
            'reconstruction_error_upper': None,
            'remaining_panel': True,
        })
    if sum(item['n_rows'] for item in reports) > SAMPLE_ROW_LIMIT:
        raise ValueError('sample exceeded the declared small-row limit')
    return {
        'panels': reports,
        'inspected_sample_panels': len(reports),
        'inspected_sample_rows': sum(item['n_rows'] for item in reports),
        'reconstruction_bounds_evaluated_panels': 0,
        'reconstruction_bounds_evaluated_rows': 0,
        'full_remaining_campaign_run': False,
    }


def validate_rho1_source_bound_record(record):
    if record.get('schema') != SCHEMA:
        raise ValueError('unexpected rho=1 source-bound schema')
    status = record.get('status')
    if not isinstance(status, str) or not status.startswith('OPEN'):
        raise ValueError('rho=1 source-bound core must remain OPEN')
    coverage = record['coverage']
    if coverage['remaining_reconstruction_panels'] != REMAINING_RECONSTRUCTION_PANELS:
        raise ValueError('remaining reconstruction panel count changed')
    if coverage['remaining_reconstruction_rows'] != REMAINING_RECONSTRUCTION_ROWS:
        raise ValueError('remaining reconstruction row count changed')
    if coverage['inventory_positive_panels'] != INVENTORY_POSITIVE_PANELS:
        raise ValueError('inventory panel count changed')
    if coverage['inventory_positive_rows'] != INVENTORY_POSITIVE_ROWS:
        raise ValueError('inventory row count changed')
    if coverage['v5_action_covered_panels'] != V5_ACTION_COVERED_PANELS:
        raise ValueError('V5 action-covered panel count changed')
    if coverage['v5_action_covered_rows'] != V5_ACTION_COVERED_ROWS:
        raise ValueError('V5 action-covered row count changed')
    partition = record.get('row_threshold_partition')
    if not isinstance(partition, dict):
        raise ValueError('authenticated row-level threshold partition required')
    expected_partition = {
        'remaining_panel_campaign_panels': REMAINING_RECONSTRUCTION_PANELS,
        'remaining_panel_campaign_rows': REMAINING_RECONSTRUCTION_ROWS,
        'rows_below_joined_in_remaining_panels':
            BELOW_JOINED_ROWS_IN_REMAINING_PANELS,
        'rows_at_or_above_joined_in_remaining_panels':
            AT_OR_ABOVE_JOINED_ROWS_IN_REMAINING_PANELS,
        'straddling_remaining_panels': STRADDLING_REMAINING_PANELS,
        'physical_column_accuracy_universe_panels': INVENTORY_POSITIVE_PANELS,
        'physical_column_accuracy_universe_rows': INVENTORY_POSITIVE_ROWS,
        'v5_action_is_not_column_reconstruction': True,
        'remaining_panel_campaign_is_not_complete_physical_source_universe': True,
    }
    for key, expected in expected_partition.items():
        if partition.get(key) != expected:
            raise ValueError('row-level V5 threshold partition changed: ' + key)
    for key in ('remaining_physical_rho1_source_error_bound',
                'remaining_low_subgap_source_error_bound'):
        if record.get(key) is not None or coverage.get(key) is not None:
            raise ValueError('v1 may not invent a missing physical rho=1 bound')
    record_directed = record.get('directed_N_beta_contribution_bound')
    coverage_directed = coverage.get('directed_N_beta_contribution_bound')
    for directed in (record_directed, coverage_directed):
        if (not isinstance(directed, dict) or directed.get('N') is not None
                or directed.get('beta') is not None):
            raise ValueError('v1 slice has no finite directed N,beta reconstruction budget')
    if record_directed != coverage_directed:
        raise ValueError('directed N,beta OPEN ledgers disagree')
    if record.get('historical_coarse_4p998e_minus9_is_current_remaining_bound'):
        raise ValueError('historical 4.998e-9 coarse aggregate is not the current remaining bound')
    sample = record['sample']
    if sample['full_remaining_campaign_run']:
        raise ValueError('v1 slice must not run all 102 remaining reconstruction panels')
    if sample['reconstruction_bounds_evaluated_panels'] != 0:
        raise ValueError('v1 slice must not claim evaluated reconstruction bounds')
    if sample.get('reconstruction_bounds_evaluated_rows') != 0:
        raise ValueError('v1 slice must not claim evaluated reconstruction bounds')
    if any(item.get('reconstruction_error_upper') is not None
           for item in sample.get('panels', ())):
        raise ValueError('sample inspection may not supply reconstruction bounds')
    if sample['inspected_sample_rows'] > SAMPLE_ROW_LIMIT:
        raise ValueError('sample exceeded the declared small-row limit')
    if record.get('redo_group12_or_group32_order24_middle_bounds'):
        raise ValueError('already-certified middle windows must be reused')
    if (record.get('physical_EXISTENCE_certificate')
            or record.get('physical_NONEXISTENCE_certificate')):
        raise ValueError('OPEN source-bound core cannot issue a physical certificate')
    upstream = record.get('upstream_physical_slot')
    if (not isinstance(upstream, dict)
            or upstream.get('physical_rho1_source_error_upper') is not None
            or upstream.get('combined_upstream_bound_N_beta') is not None
            or not isinstance(upstream.get('combined_status'), str)
            or not upstream['combined_status'].startswith('OPEN')):
        raise ValueError('v1 upstream physical source slot must remain OPEN')
    if not coverage.get('v5_action_is_not_column_reconstruction'):
        raise ValueError('V5 action coverage is not column reconstruction')
    if not coverage.get('baseline_ell0_allocations_retained'):
        raise ValueError('ell=0 baseline allocations must be retained')
    if coverage.get('double_counted'):
        raise ValueError('coverage claims a double count')
    if coverage_digest(coverage['remaining_panel_identities']) != coverage['remaining_identity_digest']:
        raise ValueError('remaining identity digest changed')
    if coverage_digest(coverage['v5_action_panel_identities']) != coverage['v5_action_identity_digest']:
        raise ValueError('V5 action identity digest changed')
    return True
