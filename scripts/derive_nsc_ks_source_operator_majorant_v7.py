#!/usr/bin/env python3
"""Join accepted v6 with the 108 disjoint replayed subgap rows.

v6 is the accepted v4 aggregate plus 18 threshold rows. v5 is that same v4
aggregate plus these 108 subgap rows. The two totals are not added. Original
weights are applied once. Saved subgap witnesses are checked by tree hash,
payload hash and archive identity; the 54-row proof replay is not repeated.
NumPy 2.5.1 drops the negative zero at each negative fiber's C[0, 1].real.
Restoring that one sign bit recovers the digest stored by v6. Values stay identical.
"""
import argparse, hashlib, json, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
sys.path.insert(0, str(ROOT / 'scripts'))
from derive_nsc_ks_source_operator_majorant import digest, read_bound, source_error_moments
from derive_nsc_ks_source_operator_majorant_v4 import require_same_slice
from derive_nsc_ks_source_operator_majorant_v5 import (
    bind_subgap_inputs, require_saved_witness, subgap_rows)
from flint import arb, ctx
from recursive_horizons.evidence_io import publish_exclusive_file
from recursive_horizons.nsc_evolved_incoming_state import _digest_arrays
from recursive_horizons.nsc_ks_ball_trajectory import exact_upper, restored_upper
from recursive_horizons.nsc_ks_evaluation_binding import implementation_hashes
from recursive_horizons.nsc_ks_massive_cubic_uv import require_original_history
from recursive_horizons.nsc_ks_retained_upstream_archive import RetainedUpstreamArchive
from recursive_horizons.nsc_ks_source_envelope import RHO_SIGMA
from recursive_horizons.nsc_ks_source_operator_majorant import history_source_insertion
from recursive_horizons.nsc_subgap_source_window import (
    EXPECTED_POSITIVE_ROWS, EXPECTED_SIGNED_ROWS, SCHEMA, checkpoint_relative, load_checkpoints,
    load_failures)
from recursive_horizons.nsc_threshold18_replay import (
    ADJUDICATION, CHECKPOINT_RELATIVE, PANEL, ROWS, _batch, require_signed_zero_only, verify_closure)
PRIOR = 'results/development/nsc-ks-source-operator-majorant-v6.json'
V5 = 'results/development/nsc-ks-source-operator-majorant-v5.json'
OUTPUT = 'results/development/nsc-ks-source-operator-majorant-v7.json'
WINDOW = 'results/development/nsc-subgap-source-window-v1'
V6_SHA256 = '15615d1ee0206ed0fb6b8c3a10df8d860c581b2ebc206bef2242492f0154d517'
V5_SHA256 = '21f9eb1214126c513e80470fa67b9e2075252beeb93a76491769607d067122ca'
SUBGAP_TREE_SHA256 = 'c310e9023e945fad5985651cc0811e0f737bdff6d8eab352e091b52ac05ee161'
PRIOR_COVERED, PRIOR_REMAINING = 1060, 108
ADDED_SIGNED = 108
THRESHOLD_SIGNED = 18
COVERED, REMAINING = 1168, 0
FAMILY_SIGNED = 1168
OWNERS = (
    'scripts/derive_nsc_ks_source_operator_majorant_v7.py',
    'scripts/derive_nsc_ks_source_operator_majorant_v6.py',
    'scripts/derive_nsc_ks_source_operator_majorant_v5.py',
    'scripts/derive_nsc_ks_source_operator_majorant_v4.py',
    'scripts/derive_nsc_ks_source_operator_majorant.py',
    'src/recursive_horizons/nsc_threshold18_replay.py',
    'src/recursive_horizons/nsc_subgap_source_window.py',
    'src/recursive_horizons/nsc_ks_source_operator_majorant.py',
    'tests/test_nsc_ks_source_operator_majorant_v7.py')


def subgap_tree_sha256(root=None):
    """Directory entries and file bytes. Directory nodes contribute no marker byte."""
    root = Path(root) if root is not None else ROOT / WINDOW
    digest = hashlib.sha256()
    for path in sorted(root.rglob('*'), key=lambda item: item.relative_to(root).as_posix()):
        digest.update(path.relative_to(root).as_posix().encode())
        digest.update(b'\0')
        if path.is_file():
            digest.update(path.read_bytes())
        digest.update(b'\0')
    if digest.hexdigest() != SUBGAP_TREE_SHA256:
        raise ValueError('subgap witness tree changed')
    return digest.hexdigest()


def accepted_v6():
    raw = (ROOT / PRIOR).read_bytes()
    if hashlib.sha256(raw).hexdigest() != V6_SHA256:
        raise ValueError('accepted v6 source aggregate changed')
    record = read_bound(PRIOR)
    coverage = record.get('source_coverage') if isinstance(record, dict) else None
    equivalence = record.get('negative_source_equivalence') if isinstance(record, dict) else None
    closure = record.get('threshold18_closure') if isinstance(record, dict) else None
    historical = (closure or {}).get('historical_objects') if isinstance(closure, dict) else None
    if (record.get('schema') != 'NSC-KS-SOURCE-OPERATOR-MAJORANT-v6'
            or record.get('family') != [14, 1] or record.get('physical_local_gate') != 'OPEN'
            or record.get('physical_upstream_budget_component') is not None
            or record.get('source_weights_applied_once') is not True
            or record.get('source_quadrature_error_included') is not False
            or record.get('numerical_field_error_included') is not False
            or record.get('value_assembly_arithmetic_included') is not False
            or record.get('negative_bounds_supplied_separately') is not True
            or not isinstance(coverage, dict)
            or coverage.get('covered_signed_rows') != PRIOR_COVERED
            or coverage.get('remaining_signed_rows_in_family') != PRIOR_REMAINING
            or coverage.get('all_source_families') is not False
            or coverage.get('complete_middle_window_replayed') is not True
            or coverage.get('complete_direct_source_window_replayed') is not True
            or coverage.get('v4_signed_rows_reused') != 1042
            or coverage.get('threshold18_signed_rows_added') != THRESHOLD_SIGNED
            or coverage.get('ode_rerun') is not False
            or coverage.get('failed_rows_entered_as_zero') is not False
            or not isinstance(equivalence, dict)
            or equivalence.get('relation') != 'signed_zero_only'
            or equivalence.get('nonzero_numeric_differences') != 0
            or equivalence.get('positive_digests_match') is not True
            or equivalence.get('batch_signed_zeros') != 16
            or equivalence.get('original_raw_source_digest') == equivalence.get('local_raw_source_digest')
            or equivalence.get('original_raw_preparation_digest') == equivalence.get('local_raw_preparation_digest')
            or not isinstance(historical, list) or len(historical) != 1
            or historical[0].get('blob') != '31818fdd6b8cc9c4d8b9d33f41e0f241b2dd1227'
            or historical[0].get('sha256') == historical[0].get('head_sha256')):
        raise ValueError('accepted v6 family14_1 coverage required')
    return record


def prior_rows(record):
    coverage = record.get('source_coverage') if isinstance(record, dict) else None
    rows = []
    for row in (coverage or {}).get('rows') or []:
        if not isinstance(row, dict) or row.get('covariance_error_upper') is None:
            raise ValueError('missing v6 source bound')
        if row.get('energy_sign') not in (-1, 1) or row.get('bound_owner') is None:
            raise ValueError('missing v6 source bound')
        rows.append({**row, 'angular_sign': row['energy_sign'], 'owner': row['bound_owner'],
                     'epsilon': row['covariance_error_upper']})
    if len(rows) != PRIOR_COVERED:
        raise ValueError('accepted v6 family14_1 coverage required')
    return rows


def _summary(record):
    payload = record.get('payload') if isinstance(record, dict) else None
    if not isinstance(payload, dict) or not isinstance(payload.get('path'), str):
        raise ValueError('missing subgap window file')
    return {
        'panel': record['panel'], 'row': record['row'],
        'checkpoint': payload['path'].rsplit('/', 1)[0],
        'energy_hex': record['energy_hex'], 'negative_energy_hex': record['negative_energy_hex'],
        'transport_profile': record['transport_profile'], 'tube': record['tube'],
        'defect_subdivisions': record['defect_subdivisions'],
        'positive_covariance_error_upper': record['positive_covariance_error_upper'],
        'negative_covariance_error_upper': record['negative_covariance_error_upper'],
        'signed_comparisons': record['signed_comparisons'],
        'quadrature_weight_applied': False, 'negative_error_copied_from_positive': False,
        'physical_upstream_budget_component': None, 'physical_local_gate': 'OPEN'}


def ready_subgap_window():
    """Read the finished witness tree. This does not call the proof replay."""
    if EXPECTED_POSITIVE_ROWS != 54 or EXPECTED_SIGNED_ROWS != ADDED_SIGNED or SCHEMA != 'NSC-SUBGAP-SOURCE-WINDOW-v1':
        raise ValueError('subgap window census changed')
    root = ROOT / WINDOW
    subgap_tree_sha256(root)
    freeze_path = root / 'proof-freeze.json'
    try:
        freeze = json.loads(freeze_path.read_bytes())
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError('missing subgap window file') from error
    proof = freeze.get('sha256') if isinstance(freeze, dict) else None
    proof_path = freeze.get('proof') if isinstance(freeze, dict) else None
    if (not isinstance(proof, str) or not isinstance(proof_path, str)
            or hashlib.sha256((ROOT / proof_path).read_bytes()).hexdigest() != proof
            or freeze.get('positive_rows') != 54 or freeze.get('signed_rows') != ADDED_SIGNED):
        raise ValueError('frozen proof code changed')
    if load_failures(root):
        raise ValueError('a failed row must not be entered as zero')
    loaded = load_checkpoints(root)
    if len(loaded) != 54:
        raise ValueError('subgap census does not add up')
    rows = []
    for _dirname, record, payload in loaded:
        require_saved_witness(record)
        described = record.get('payload')
        if (record.get('proof_sha256') != proof
                or record.get('replay_uses_saved_phase_and_bloch_witness') is not True
                or record.get('preparation_jost_rerun_on_replay') is not False
                or record.get('bloch_ode_rerun_on_replay') is not False
                or not isinstance(described, dict)
                or described.get('sha256') != hashlib.sha256(payload).hexdigest()
                or described.get('bytes') != len(payload)):
            raise ValueError('replayed subgap window required')
        summary = _summary(record)
        if summary['checkpoint'] != checkpoint_relative(record['panel'], record['row']):
            raise ValueError('subgap witness does not match the established control')
        rows.append(summary)
    return {
        'schema': SCHEMA, 'mode': 'check', 'new_rows_captured': 0,
        'replay_uses_saved_phase_and_bloch_witness': True,
        'preparation_jost_rerun_on_replay': False, 'quadrature_weight_applied': False,
        'covariance_errors_unweighted': True, 'group': 14, 'angular_sign': 1,
        'failed_rows_entered_as_zero': False, 'missing_rows_entered_as_zero': False,
        'finite_window_closes_gate': False, 'physical_upstream_budget_component': None,
        'physical_local_gate': 'OPEN', 'negative_bounds_supplied_separately': True,
        'source_columns_replaced': False, 'all_source_families': False,
        'rows': rows, 'failed': [], 'missing': [],
        'completed_positive_rows': 54, 'completed_signed_rows': ADDED_SIGNED,
        'failed_positive_rows': 0, 'missing_positive_rows': 0, 'coverage_complete': True}


def add_disjoint_rows(rows, added):
    keys = {(row['panel'], row['row'], row['energy_sign']) for row in rows}
    if len(keys) != len(rows) or len(keys) != PRIOR_COVERED:
        raise ValueError('duplicate bounded source row')
    for item in added:
        key = (item['panel'], item['row'], item['energy_sign'])
        if key in keys:
            raise ValueError('duplicate subgap source row')
        if item.get('epsilon') is None:
            raise ValueError('missing subgap source bound')
        keys.add(key)
    if len(added) != ADDED_SIGNED or len(keys) != COVERED:
        raise ValueError('source census changed')
    return list(rows) + list(added)


def _row_key(row):
    return (row['panel'], row['row'], row['energy_sign'])


def require_disjoint_complements(base, added, threshold):
    """v5's extra rows and v6's extra rows meet only on the shared v4 keys."""
    raw = (ROOT / V5).read_bytes()
    if hashlib.sha256(raw).hexdigest() != V5_SHA256:
        raise ValueError('accepted v5 subgap aggregate changed')
    try:
        record = json.loads(raw)
    except json.JSONDecodeError as error:
        raise ValueError('accepted v5 subgap aggregate changed') from error
    coverage = record.get('source_coverage') if isinstance(record, dict) else None
    if (not isinstance(coverage, dict) or coverage.get('covered_signed_rows') != 1150
            or coverage.get('remaining_signed_rows_in_family') != THRESHOLD_SIGNED
            or coverage.get('subgap_signed_rows_added') != ADDED_SIGNED
            or coverage.get('subgap_positive_rows_witnessed') != 54
            or coverage.get('subgap_positive_rows_failed') != 0
            or coverage.get('subgap_positive_rows_missing') != 0
            or coverage.get('v4_signed_rows_reused') != 1042
            or coverage.get('failed_rows_entered_as_zero') is not False
            or record.get('source_weights_applied_once') is not True):
        raise ValueError('accepted v5 subgap aggregate changed')
    v5_rows = {_row_key(row): row for row in coverage.get('rows') or []}
    v6_keys = {_row_key(row) for row in base}
    added_map = {_row_key(row): row for row in added}
    threshold_keys = {_row_key(row) for row in threshold}
    if (len(v5_rows) != 1150 or len(v6_keys) != PRIOR_COVERED or len(added_map) != ADDED_SIGNED
            or len(threshold_keys) != THRESHOLD_SIGNED):
        raise ValueError('source census changed')
    if not threshold_keys <= v6_keys or not threshold_keys.isdisjoint(added_map):
        raise ValueError('threshold rows overlap the subgap complement')
    if set(v5_rows) - v6_keys != set(added_map):
        raise ValueError('subgap rows are not the disjoint v5 complement')
    if v6_keys - set(v5_rows) != threshold_keys:
        raise ValueError('threshold rows are not the disjoint v6 complement')
    if len(v6_keys | set(added_map)) != COVERED:
        raise ValueError('overlapping v5 and v6 totals were added')
    for key, item in added_map.items():
        stored = v5_rows[key]
        if (stored.get('covariance_error_upper') != item.get('epsilon')
                or stored.get('source_digest') != item.get('source_digest')
                or stored.get('preparation_digest') != item.get('preparation_digest')
                or stored.get('energy_hex') != item.get('energy_hex')):
            raise ValueError('subgap identity does not match the replayed witness row')
    return v5_rows


def signed_zero_spell(batch):
    """v6's negative digest is the live batch plus one negative zero at each C[0, 1].real."""
    if batch.energy_sign >= 0:
        raise ValueError('negative batch required')
    covariance = np.array(batch.source.covariance, copy=True)
    fibers = covariance.shape[0] // 3
    if fibers < 1 or covariance.shape != (3 * fibers, 3 * fibers):
        raise ValueError('negative batch fiber census changed')
    view = covariance.view(np.float64)
    for index in range(fibers):
        flat = ((3 * index) * covariance.shape[1] + (3 * index + 1)) * 2
        slot = view.ravel()[flat]
        if slot != 0 or np.signbit(slot):
            raise ValueError('negative covariance corner is not a positive zero')
        view.ravel()[flat] = np.copysign(slot, -1.0)
    compared = require_signed_zero_only(batch.source.covariance, covariance, 'restored signed zeros')
    if compared['signed_zero_differences'] != fibers or compared['raw_equal'] or compared['max_abs'] != 0.0:
        raise ValueError('signed-zero restoration changed a numeric value')
    signed = _digest_arrays(covariance, batch.source.column_weights, batch.source.energies)
    if signed == batch.source.digest:
        raise ValueError('signed zero was already the live digest')
    initial = np.ascontiguousarray(np.asarray(batch.initial_columns, complex))
    preparation = hashlib.sha256(
        signed.encode() + initial.tobytes() + repr((
            float(batch.mass), float(batch.angular), float(batch.rho_up), float(RHO_SIGMA))).encode()
    ).hexdigest()
    return signed, preparation


def verify_threshold_identity(archive, prior):
    """Exact values match. NumPy 2.5.1 no longer keeps the negative-zero bytes."""
    equivalence = prior['negative_source_equivalence']
    positive, negative = _batch(archive)
    if (negative.source.digest != equivalence['original_raw_source_digest']
            or negative.preparation_digest != equivalence['original_raw_preparation_digest']
            or positive.source.digest != equivalence['positive_source_digest']
            or positive.preparation_digest != equivalence['positive_preparation_digest']):
        raise ValueError('live source digest does not match the recorded exact-value identity')
    signed_source, signed_preparation = signed_zero_spell(negative)
    if (signed_source != equivalence['local_raw_source_digest']
            or signed_preparation != equivalence['local_raw_preparation_digest']):
        raise ValueError('restoring signed zeros did not recover the local low32 digest')
    rows_root = ROOT / CHECKPOINT_RELATIVE / 'rows'
    closure = None
    for index in ROWS:
        directory = rows_root / f'row-{PANEL.replace("/", "__")}-{index:05d}'
        record = json.loads((directory / 'record.json').read_bytes())
        payload = (directory / 'witness.npz').read_bytes()
        described = record.get('payload') if isinstance(record, dict) else None
        if (not isinstance(described, dict) or described.get('sha256') != hashlib.sha256(payload).hexdigest()
                or described.get('bytes') != len(payload)):
            raise ValueError('missing payload')
        if closure is None:
            closure = verify_closure(ROOT, record)
        if (record.get('negative_endpoint_adjudication') != ADJUDICATION
                or record.get('negative_density_endpoint') != '(nx,-ny,-nz)'
                or record.get('weight_applied_to_error') is not False):
            raise ValueError('imported dyadic replay required')
        local_index = index - 32
        sl = slice(3 * local_index, 3 * local_index + 3)
        with np.load(directory / 'witness.npz', allow_pickle=False) as saved:
            witness = {name: np.array(saved[name]) for name in saved.files}
        for name, local, imported in (
                ('positive_columns', positive.initial_columns[:, sl], witness['positive_columns']),
                ('negative_columns', negative.initial_columns[:, sl], witness['negative_columns']),
                ('positive_covariance', positive.source.covariance[sl, sl], witness['positive_covariance']),
                ('negative_covariance', negative.source.covariance[sl, sl], witness['negative_covariance'])):
            compared = require_signed_zero_only(local, imported, f'{index}:{name}')
            if compared['signed_zero_differences'] or not compared['raw_equal'] or compared['max_abs'] != 0.0:
                raise ValueError('numerical source arrays differ: ' + name)
        stored = {_row_key(row): row for row in prior_rows(prior)
                  if row['panel'] == PANEL and row['row'] == index}
        comparisons = record.get('signed_comparisons')
        if not isinstance(comparisons, list) or [item.get('energy_sign') for item in comparisons] != [1, -1]:
            raise ValueError('separate signed source bounds required')
        for item in comparisons:
            row = stored.get((PANEL, index, item.get('energy_sign')))
            if row is None or row['epsilon'] != item.get('signed_total_upper'):
                raise ValueError('threshold dyadic does not match the accepted v6 row')
            if item.get('weight_applied_to_error') is not False:
                raise ValueError('source weight already applied')
    resolution = equivalence.get('propagator_resolution')
    historical = closure.get('historical_objects') if isinstance(closure, dict) else None
    if (not isinstance(resolution, dict) or historical != [resolution]
            or resolution.get('blob') != '31818fdd6b8cc9c4d8b9d33f41e0f241b2dd1227'
            or resolution.get('sha256') == resolution.get('head_sha256')
            or closure.get('source_hashes_checked') != 67
            or closure.get('archive_digests_checked') != 64
            or closure.get('base_commit') != 'eebbe0c3b0bfe76be0bae4a64d3b52ba8da7cfb9'):
        raise ValueError('historical propagator object changed')
    return equivalence, closure, negative.source.digest, negative.preparation_digest


def rebind_signed_zero_spelling(rows, entries):
    """Replace a proved signed-zero digest with the live exact-value digest."""
    batches = {}
    spells = {}
    for batch, _channel in entries:
        if batch.energy_sign < 0:
            spells[id(batch)] = signed_zero_spell(batch)
        for index, _energy in enumerate(batch.source.energies[::3]):
            key = (batch.original_panel, batch.rows[0] + index, batch.energy_sign)
            if key in batches:
                raise ValueError('duplicate original source row')
            batches[key] = batch
    rebound = []
    count = 0
    for row in rows:
        key = (row['panel'], row['row'], row['energy_sign'])
        if key not in batches:
            raise ValueError('bounded row is not in the original source')
        batch = batches[key]
        live = (batch.source.digest, batch.preparation_digest)
        current = (row['source_digest'], row['preparation_digest'])
        if current == live:
            rebound.append(row)
            continue
        if batch.energy_sign > 0 or current != spells[id(batch)]:
            raise ValueError('source bound does not match its original preparation')
        count += 1
        rebound.append({**row, 'source_digest': live[0], 'preparation_digest': live[1]})
    if count < 1:
        raise ValueError('signed-zero row census changed')
    return rebound, count


def indexed_keys(entries):
    indexed = set()
    for batch, _channel in entries:
        for index, _energy in enumerate(batch.source.energies[::3]):
            key = (batch.original_panel, batch.rows[0] + index, batch.energy_sign)
            if key in indexed:
                raise ValueError('duplicate original source row')
            indexed.add(key)
    return indexed


def refuse_mutated_identity(archive, entries, rows):
    altered = [dict(item) for item in rows]
    altered[-1] = {**altered[-1], 'preparation_digest': '0' * 64}
    try:
        source_error_moments(archive, entries, altered)
    except ValueError as error:
        if 'original preparation' not in str(error):
            raise
    else:
        raise ValueError('mutated preparation was accepted')
    negative = [dict(item) for item in rows]
    negative[-1] = {**negative[-1], 'epsilon': {'exponent': 0, 'mantissa': '-1'}}
    try:
        source_error_moments(archive, entries, negative)
    except ValueError as error:
        if 'nonnegative' not in str(error):
            raise
    else:
        raise ValueError('negative covariance bound was accepted')
    try:
        source_error_moments(archive, entries, list(rows) + [rows[-1]])
    except ValueError as error:
        if 'duplicate' not in str(error):
            raise
    else:
        raise ValueError('duplicate source row was accepted')


def exercise_pure_refusals(base):
    """Reject overlap, a resumed window, a pre-weighted error and a copied sign."""
    try:
        add_disjoint_rows(base, base[:1])
    except ValueError as error:
        if 'duplicate' not in str(error):
            raise
    else:
        raise ValueError('duplicate subgap row was accepted')
    window = {
        'schema': SCHEMA, 'mode': 'resume', 'new_rows_captured': 1,
        'replay_uses_saved_phase_and_bloch_witness': True,
        'preparation_jost_rerun_on_replay': False, 'quadrature_weight_applied': False,
        'covariance_errors_unweighted': True, 'group': 14, 'angular_sign': 1,
        'failed_rows_entered_as_zero': False, 'missing_rows_entered_as_zero': False,
        'finite_window_closes_gate': False, 'physical_upstream_budget_component': None,
        'physical_local_gate': 'OPEN', 'negative_bounds_supplied_separately': True,
        'source_columns_replaced': False, 'all_source_families': False,
        'rows': [], 'failed': [], 'missing': [],
        'completed_positive_rows': 0, 'completed_signed_rows': 0,
        'failed_positive_rows': 0, 'missing_positive_rows': 0}
    try:
        subgap_rows(window)
    except ValueError as error:
        if 'replayed' not in str(error):
            raise
    else:
        raise ValueError('resumed subgap window was accepted')
    weighted = dict(window)
    weighted['mode'] = 'check'
    weighted['new_rows_captured'] = 0
    weighted['quadrature_weight_applied'] = True
    try:
        subgap_rows(weighted)
    except ValueError as error:
        if 'weight' not in str(error):
            raise
    else:
        raise ValueError('weighted subgap error was accepted')
    copied = {'negative_complement_applied': True, 'negative_error_copied_from_positive': True,
              'quadrature_weight_applied': False, 'method_too_loose': False,
              'covariance_error_unweighted': True, 'physical_local_gate': 'OPEN'}
    try:
        require_saved_witness(copied)
    except ValueError as error:
        if 'covariance complement' not in str(error):
            raise
    else:
        raise ValueError('copied negative error was accepted')


def calculate():
    prior = accepted_v6()
    base = prior_rows(prior)
    exercise_pure_refusals(base)
    archive = RetainedUpstreamArchive(ROOT)
    rho, entries = require_same_slice(prior, archive)
    window = ready_subgap_window()
    added = subgap_rows(window)
    equivalence, closure, _live_source, _live_preparation = verify_threshold_identity(archive, prior)
    threshold = [{'panel': PANEL, 'row': index, 'energy_sign': sign}
                 for index in ROWS for sign in (1, -1)]
    require_disjoint_complements(base, added, threshold)
    rebound_base, _base_rebound = rebind_signed_zero_spelling(base, entries)
    rows, rebound_count = rebind_signed_zero_spelling(add_disjoint_rows(base, added), entries)
    family, identity = require_original_history(ROOT)
    if identity != prior.get('profile_identity'):
        raise ValueError('profile identity changed')
    channel = archive.meta['channels'][14]
    if channel['copy_count'] * channel['degeneracy'] / 2 != 12:
        raise ValueError('multiplicity changed')
    with ctx.workprec(192):
        _base_moment, _base_energy, base_covered, base_missing = source_error_moments(
            archive, entries, rebound_base)
        if len(base_covered) != PRIOR_COVERED or base_missing != PRIOR_REMAINING:
            raise ValueError('source census changed')
        refuse_mutated_identity(archive, entries, rows)
        moment, energy, covered, missing = source_error_moments(archive, entries, rows)
        universe = indexed_keys(entries)
        covered_keys = {_row_key(row) for row in covered}
        missing_keys = sorted(universe - covered_keys)
        if (len(universe) != FAMILY_SIGNED or len(covered) != COVERED or missing != REMAINING
                or missing_keys or len(covered_keys) != COVERED):
            raise ValueError('source census changed')
        bound = history_source_insertion(
            family, rho, moment, energy,
            mass=arb(channel['compact_mass']).union(arb.pi() / 2),
            absolute_angular=arb(channel['angular_eigenvalue']).union(arb(5).sqrt()),
            multiplicity=12)
        previous = prior['continuous_action_error_upper']
        if not (restored_upper(bound['N']) > restored_upper(previous['N'])
                and restored_upper(bound['beta']) > restored_upper(previous['beta'])):
            raise ValueError('joined source bound did not increase')
    return {
        'schema': 'NSC-KS-SOURCE-OPERATOR-MAJORANT-v7',
        'status': 'continuous source-error bound for 1168 signed rows of family14_1; full upstream component OPEN',
        'profile_identity': identity, 'rho_up_hex': rho.hex(), 'family': [14, 1],
        'multiplicity_per_signed_family': 12,
        'weighted_error_moment': exact_upper(moment),
        'weighted_energy_error_moment': exact_upper(energy),
        'continuous_action_error_upper': bound,
        'negative_source_equivalence': equivalence,
        'threshold18_closure': closure,
        'source_coverage': {
            'covered_signed_rows': COVERED, 'remaining_signed_rows_in_family': REMAINING,
            'missing_keys': [list(key) for key in missing_keys],
            'rows': covered, 'all_source_families': False,
            'complete_middle_window_replayed': True,
            'complete_direct_source_window_replayed': True,
            'family14_1_signed_rows_complete': True,
            'v6_signed_rows_reused': PRIOR_COVERED,
            'v4_signed_rows_inside_v6': 1042,
            'threshold18_signed_rows_retained': THRESHOLD_SIGNED,
            'subgap_signed_rows_added': ADDED_SIGNED,
            'subgap_positive_rows_witnessed': 54,
            'subgap_positive_rows_failed': 0,
            'subgap_positive_rows_missing': 0,
            'v5_covered_total_not_added': True,
            'subgap_tree_sha256': SUBGAP_TREE_SHA256,
            'subgap_proof_replay_repeated': False,
            'negative_rows_rebound_to_cleared_digest': rebound_count,
            'negative_signed_zero_corner': 'C[0, 1].real',
            'ode_rerun': False, 'failed_rows_entered_as_zero': False},
        'sign_audit': ADJUDICATION,
        'negative_bloch_convention': '(nx,-ny,-nz)',
        'covariance_complement_included': True, 'negative_bounds_supplied_separately': True,
        'source_weights_applied_once': True, 'source_quadrature_error_included': False,
        'numerical_field_error_included': False, 'value_assembly_arithmetic_included': False,
        'physical_upstream_budget_component': None, 'physical_local_gate': 'OPEN',
        'source_hashes': implementation_hashes(ROOT, owners=OWNERS),
        'input_hashes': bind_subgap_inputs(dict(prior['input_hashes']), window)}


def encoded_result(result):
    return (json.dumps(result, sort_keys=True, indent=2, allow_nan=False) + '\n').encode()


def _snapshot(path):
    path = Path(path)
    digest = hashlib.sha256()
    if path.is_file():
        digest.update(path.read_bytes())
        return digest.hexdigest()
    for item in sorted(path.rglob('*'), key=lambda candidate: candidate.relative_to(path).as_posix()):
        digest.update(item.relative_to(path).as_posix().encode())
        digest.update(b'\0')
        if item.is_symlink():
            digest.update(b'link')
        elif item.is_file():
            digest.update(item.read_bytes())
        else:
            digest.update(b'dir')
        digest.update(b'\0')
    return digest.hexdigest()


HISTORIC = (
    PRIOR, V5, 'results/development/nsc-ks-source-operator-majorant-v4.json',
    WINDOW, 'results/development/nsc-threshold18-import-v1',
    'results/development/nsc-direct-source-window-v1', '.source-history',
    'src/recursive_horizons/nsc_ks_energy_propagator.py')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--record', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    before = {path: _snapshot(ROOT / path) for path in HISTORIC}
    result = calculate()
    encoded = encoded_result(result)
    if {path: _snapshot(ROOT / path) for path in HISTORIC} != before:
        raise ValueError('historic bytes changed')
    if args.record:
        if (ROOT / OUTPUT).exists():
            if (ROOT / OUTPUT).read_bytes() != encoded:
                raise ValueError('existing v7 record does not match recomputation')
        else:
            publish_exclusive_file(ROOT, OUTPUT, encoded)
    elif (ROOT / OUTPUT).read_bytes() != encoded:
        raise ValueError('source insertion or dependencies changed')
    action = result['continuous_action_error_upper']
    coverage = result['source_coverage']
    equivalence = result['negative_source_equivalence']
    print(json.dumps({
        'covered_signed_rows': coverage['covered_signed_rows'],
        'remaining_signed_rows_in_family': coverage['remaining_signed_rows_in_family'],
        'missing_keys': coverage['missing_keys'],
        'N': float(restored_upper(action['N'])),
        'beta': float(restored_upper(action['beta'])),
        'N_mantissa': action['N']['mantissa'], 'N_exponent': action['N']['exponent'],
        'beta_mantissa': action['beta']['mantissa'], 'beta_exponent': action['beta']['exponent'],
        'subgap_tree_sha256': coverage['subgap_tree_sha256'],
        'original_raw_source_digest': equivalence['original_raw_source_digest'],
        'local_raw_source_digest': equivalence['local_raw_source_digest'],
        'v7_sha256': hashlib.sha256(encoded).hexdigest(),
        'physical_local_gate': result['physical_local_gate'],
        'physical_upstream_budget_component': result['physical_upstream_budget_component']}, indent=2))
