"""Add the shared cubic leading coefficient to the saved value residual.

The cubic inventory returns formal moments and leaves the physical residual
untouched. This owner is the executable successor: the leading term, not the
raw coefficient, is added on the exact target nodes of the saved finite-cutoff
value. A center scalar is never copied onto the other nodes. C_M stays
unconstructed, so the successor cannot be certified.
"""
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from flint import arb

from .nsc_ks_ball_trajectory import exact_upper
from .nsc_ks_cubic_inventory import (
    LARGE_CUTOFF_GROUPS, channel_moments, combine_basis, load_inventory)
from .nsc_ks_cutoff_bridge_evaluator import signed_family_angular_square_weights
from .nsc_ks_evaluation_binding import authenticate_payload
from .nsc_ks_finite_history_uv_coefficients import ARCHIVE_RHO_UP
from .nsc_ks_massive_cubic_uv import (
    CURRENT_HISTORY_IDENTITY, require_original_history)
from .nsc_ks_value_evaluator import (
    TARGET_COUNTS, assemble_residual, edge_change, geometry_change, merit,
    target_nodes)

import derive_nsc_evolved_incoming_constraints as evolved

VALUE_LABEL = 'v2-current-dop-g1024-d48-metadata2'
VALUE_RECORD = 'results/development/nsc-ks-gate-value-%s.json' % VALUE_LABEL
RETAINED_RHO_PAYLOAD = (
    'results/development/artifacts/nsc-ks-retained-response-sum/family-01-+1.npz')
SURFACE_COEFFICIENTS = 'results/development/nsc-incoming-surface-coefficients.json'
SURFACE_BRANCH = 'results/development/nsc-incoming-surface-regular-branch.json'
SOURCE_UPDATE = 'results/development/nsc-incoming-source-update-v5.json'
HISTORY_RECORD = 'results/development/nsc-ks-gate-history-lm-broyden.json'
CHANNEL_NAMES = ('N2', 'N4', 'NM', 'J2', 'J4', 'JM')
N_CHANNELS = ('N2', 'N4', 'NM')
BETA_CHANNELS = ('J2', 'J4', 'JM')
EXTRA_MASS_KEYS = ('surface_mass_correction', 'surface_mass', 'raw_coefficient')
EXPECTED_ANGULAR_SQUARE_WEIGHT = 107952.0


def _lab_root(root):
    root = Path(root).resolve()
    if not (root / 'src' / 'recursive_horizons').is_dir():
        raise ValueError('lab root required')
    return root


def _digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def _solver_token(value):
    if isinstance(value, str) and (value.startswith('0x') or value.startswith('-0x')):
        return float.fromhex(value)
    if isinstance(value, float):
        return float(value)
    return value


def _same_solver(left, right):
    """Family records store binary hex; the top record stores the same floats."""
    if not isinstance(left, dict) or not isinstance(right, dict) or set(left) != set(right):
        return False
    for key, item in left.items():
        if _solver_token(item) != _solver_token(right[key]):
            return False
    return True


def _endpoints(value):
    ball = value if isinstance(value, arb) else arb(float(value))
    if not ball.is_finite():
        raise ValueError('finite node required')
    return ball.lower().man_exp(), ball.upper().man_exp()


def _interval(ball):
    return {
        'lower': exact_upper(ball.lower()),
        'upper': exact_upper(ball.upper()),
        'radius_upper_float': float(ball.rad()),
    }


def midpoint_float(ball):
    """Binary64 image of the ball midpoint. This is not an enclosure."""
    return float(ball.mid())


def require_moments(moments):
    """Inventory moments only: both cuts already applied, no second energy fold."""
    if not isinstance(moments, dict):
        raise ValueError('cubic moments required')
    if moments.get('energy_folding_factor') != 1:
        raise ValueError('energy folding factor must stay 1')
    if moments.get('angular_families') != 60:
        raise ValueError('omitted angular family')
    if moments.get('C_M') is not None:
        raise ValueError('higher remainder was not constructed')
    if moments.get('physical_local_gate') != 'OPEN':
        raise ValueError('local gate remains OPEN')
    leading, raw = moments.get('leading'), moments.get('raw')
    if (not isinstance(leading, tuple) or not isinstance(raw, tuple)
            or len(leading) != 3 or len(raw) != 3):
        raise ValueError('three cubic moments required')
    for lead, original in zip(leading, raw):
        if not isinstance(lead, arb) or not isinstance(original, arb):
            raise ValueError('moment enclosures required')
        if not lead.is_finite() or not original.is_finite() or not lead >= 0 or not original >= 0:
            raise ValueError('nonnegative cubic moments required')
        if lead == original:
            raise ValueError('leading moment collapsed onto the raw moment')
    return moments


def _original_moments(root):
    rows = load_inventory(root)
    large = {int(row['group']) for row in rows if float(row['cutoff']) == 320.0}
    small = {int(row['group']) for row in rows if float(row['cutoff']) == 160.0}
    if large != set(LARGE_CUTOFF_GROUPS) or large & small or len(rows) != 60:
        raise ValueError('original family cutoff changed')
    return rows, require_moments(channel_moments(rows))


def _node_index(value, count):
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
        raise ValueError('node-profile mismatch')
    index = int(value)
    if index < 0 or index >= count:
        raise ValueError('node-profile mismatch')
    return index


def _channel_sequence(value, count, name):
    if isinstance(value, (arb, float, int, np.floating, np.integer)) or isinstance(value, bool):
        raise ValueError('center scalar must not be broadcast')
    if isinstance(value, np.ndarray) and value.ndim == 0:
        raise ValueError('center scalar must not be broadcast')
    try:
        sequence = list(value)
    except TypeError as error:
        raise ValueError('center scalar must not be broadcast') from error
    if len(sequence) == 1 and count != 1:
        raise ValueError('center scalar must not be broadcast')
    if len(sequence) != count:
        raise ValueError('node-profile mismatch')
    balls = []
    for item in sequence:
        if isinstance(item, bool) or item is None:
            raise ValueError('explicit ' + name + ' required')
        ball = item if isinstance(item, arb) else arb(item)
        if not ball.is_finite():
            raise ArithmeticError('finite ' + name + ' required')
        balls.append(ball)
    return balls


def _characteristic(basis, sign, count, identity):
    if not isinstance(basis, dict):
        raise ValueError('characteristic basis required')
    if type(basis.get('sign')) is not int or basis['sign'] != sign:
        raise ValueError('both characteristic signs required')
    if basis.get('massless_A') is not True:
        raise ValueError('massless unit-angular basis required')
    if basis.get('profile_identity') != identity:
        raise ValueError('original history 0b0e4ced required')
    if basis.get('C_M') is not None or basis.get('higher_uv_remainder_C_M') is not None:
        raise ValueError('higher remainder was not constructed')
    extra = set(basis.get('channels', {})) & set(EXTRA_MASS_KEYS)
    if extra or any(key in basis for key in EXTRA_MASS_KEYS):
        raise ValueError('NM already contains the surface mass term')
    channels = basis.get('channels')
    if not isinstance(channels, dict):
        raise ValueError('omitted cubic channel')
    unknown = set(channels) - set(CHANNEL_NAMES)
    if unknown:
        raise ValueError('NM already contains the surface mass term')
    missing = [name for name in CHANNEL_NAMES if name not in channels]
    if missing:
        raise ValueError('omitted cubic channel ' + missing[0])
    z_domain = _channel_sequence(basis.get('z_domain'), count, 'z_domain')
    rho = basis.get('rho_up')
    if isinstance(rho, bool) or rho is None:
        raise ValueError('explicit rho_up required')
    rho_ball = rho if isinstance(rho, arb) else arb(float(rho))
    evaluated = {}
    for name in CHANNEL_NAMES:
        evaluated[name] = _channel_sequence(channels[name], count, name)
    ends = [_endpoints(item) for item in z_domain]
    if len(set(ends)) != count:
        raise ValueError('center scalar must not be broadcast')
    return {
        'sign': sign, 'massless_A': True, 'profile_identity': identity,
        'C_M': None, 'rho_up': rho_ball, 'z_domain': z_domain, 'channels': evaluated,
        'endpoints': ends,
    }


def _subset(expected_nodes, node_indices):
    expected = np.asarray(expected_nodes, dtype=float)
    if expected.ndim != 1 or expected.size < 1 or not np.isfinite(expected).all():
        raise ValueError('exact target nodes required')
    if len({float(item).hex() for item in expected}) != expected.size:
        raise ValueError('node-profile mismatch')
    if node_indices is None:
        indices = list(range(expected.size))
    else:
        indices = [_node_index(item, expected.size) for item in node_indices]
        if len(set(indices)) != len(indices):
            raise ValueError('node-profile mismatch')
    return expected, indices, sorted(indices) == list(range(expected.size))


def apply_cubic_leading(value_residual, characteristic_bases, root, expected_nodes, *,
                         node_indices=None, coefficient_source='supplied-nodal-basis'):
    """Add leading N (minus) and beta (plus) once, on the supplied exact nodes.

    ``characteristic_bases`` maps each characteristic sign to nodal channels.
    Sequences shorter than the node list, a bare scalar, or one repeated
    ``z_domain`` are rejected instead of being broadcast. Raw moments and a
    second copy of the NM surface term are not added. The returned floats are
    midpoints; certification stays blocked while C_M is None.
    """
    root = _lab_root(root)
    _, identity = require_original_history(root)
    if identity != CURRENT_HISTORY_IDENTITY or not identity.startswith('0b0e4ced'):
        raise ValueError('original history 0b0e4ced required')
    _, moments = _original_moments(root)
    expected, indices, complete = _subset(expected_nodes, node_indices)
    count = len(indices)
    if not isinstance(characteristic_bases, dict) or set(characteristic_bases) != {1, -1}:
        raise ValueError('both characteristic signs required')
    prepared = {
        sign: _characteristic(characteristic_bases[sign], sign, count, identity)
        for sign in (1, -1)}
    if _endpoints(prepared[1]['rho_up']) != _endpoints(prepared[-1]['rho_up']):
        raise ValueError('node-profile mismatch')
    selected = expected[indices]
    for sign in (1, -1):
        if prepared[sign]['endpoints'] != [_endpoints(node) for node in selected]:
            raise ValueError('node-profile mismatch')
    value = np.asarray(value_residual, dtype=float)
    if value.shape == (2,):
        raise ValueError('center scalar must not be broadcast')
    if value.ndim != 2 or value.shape != (count, 2) or not np.isfinite(value).all():
        raise ValueError('node-profile mismatch')
    weights = moments['leading']
    leading_n = []
    leading_beta = []
    enclosures = []
    for position, node in enumerate(selected):
        pair = {}
        for sign in (1, -1):
            item = prepared[sign]
            pair[sign] = {
                'sign': sign, 'massless_A': True, 'profile_identity': identity,
                'C_M': None, 'rho_up': item['rho_up'], 'z_domain': item['z_domain'][position],
            }
            for name in CHANNEL_NAMES:
                pair[sign][name] = item['channels'][name][position]
        combined = combine_basis(pair, moments)
        if combined.get('C_M') is not None or combined.get('physical_residual_modified') is not False:
            raise ValueError('inventory assembly must leave the residual unmodified')
        manual_n = arb(0)
        manual_beta = arb(0)
        for sign in (1, -1):
            manual_n -= sum(
                (weight * pair[sign][name] for weight, name in zip(weights, N_CHANNELS)), arb(0))
            manual_beta += sum(
                (weight * pair[sign][name] for weight, name in zip(weights, BETA_CHANNELS)), arb(0))
        if not (combined['leading']['N'] - manual_n).contains(0):
            raise ValueError('N sign drifted from minus')
        if not (combined['leading']['beta'] - manual_beta).contains(0):
            raise ValueError('beta sign drifted from plus')
        if combined['leading']['N'] == combined['raw']['N'] or combined['leading']['beta'] == combined['raw']['beta']:
            raise ValueError('raw cubic coefficient was added instead of the leading term')
        leading_n.append(combined['leading']['N'])
        leading_beta.append(combined['leading']['beta'])
        enclosures.append({
            'node_index': indices[position],
            'z_hex': float(node).hex(),
            'N': _interval(combined['leading']['N']),
            'beta': _interval(combined['leading']['beta']),
        })
    if all(item == 0 for item in leading_n) and all(item == 0 for item in leading_beta):
        raise ValueError('cubic leading vanished')
    approximate = np.column_stack((
        [midpoint_float(item) for item in leading_n],
        [midpoint_float(item) for item in leading_beta])).astype(float)
    if approximate.shape != value.shape or not np.isfinite(approximate).all():
        raise ArithmeticError('approximate cubic leading is not finite')
    if np.array_equal(approximate, 0):
        raise ValueError('cubic leading vanished')
    residual = value + approximate
    report = {
        'residual': residual,
        'residual_before': value.copy(),
        'approximate_leading': approximate,
        'approximate_is_enclosure': False,
        'leading_N': tuple(leading_n),
        'leading_beta': tuple(leading_beta),
        'leading_enclosures': enclosures,
        'node_indices': indices,
        'node_count': int(expected.size),
        'completed_node_count': count,
        'all_exact_target_nodes': bool(complete),
        'center_was_broadcast': False,
        'value_residual_array_modified': not np.array_equal(residual, value),
        'physical_gate_modified': False,
        'physical_residual_modified_by_inventory_metadata': False,
        'cubic_leading_applications': 1,
        'used_raw_coefficient': False,
        'surface_mass_added_again': False,
        'energy_folding_factor': 1,
        'large_cutoff_groups': tuple(sorted(LARGE_CUTOFF_GROUPS)),
        'higher_uv_remainder_C_M': None,
        'C_M': None,
        'uniform_C4_on_I': None,
        'complete_UV_tail': None,
        'higher_remainder_constructed': False,
        'coefficient_source': coefficient_source,
        'profile_identity': identity,
        'physical_local_gate': 'OPEN',
        'physical_EXISTENCE_certificate': False,
        'physical_NONEXISTENCE_certificate': False,
        'successor_is_certified': False,
        'z_hex': [float(node).hex() for node in selected],
    }
    if not report['value_residual_array_modified']:
        raise ValueError('omitted cubic leading')
    return report


def certification_blockers(report):
    if not isinstance(report, dict):
        raise ValueError('cubic residual report required')
    blockers = []
    if report.get('higher_uv_remainder_C_M') is None or report.get('C_M') is None:
        blockers.append('higher remainder C_M is None')
    if report.get('higher_remainder_constructed') is not True:
        blockers.append('higher remainder was not constructed')
    if report.get('uniform_C4_on_I') is None or report.get('complete_UV_tail') is None:
        blockers.append('uniform C4 on I is None')
    if not report.get('all_exact_target_nodes'):
        blockers.append('cubic leading does not cover every exact target node')
    if report.get('approximate_is_enclosure') is not False:
        blockers.append('approximate values are not enclosures')
    if report.get('center_was_broadcast') is not False:
        blockers.append('center scalar was broadcast')
    if report.get('physical_local_gate') != 'OPEN':
        blockers.append('local gate remains OPEN')
    if report.get('successor_is_certified') is not False:
        blockers.append('successor is not certified')
    return blockers


def certify(report):
    """Refuse a certificate while the higher remainder is unconstructed."""
    blockers = certification_blockers(report)
    if blockers:
        raise ValueError('certification blocked: ' + '; '.join(blockers))
    raise ValueError('certification blocked: higher remainder C_M is None')


def _retained_rho_up(root):
    path = root / RETAINED_RHO_PAYLOAD
    with np.load(path, allow_pickle=False) as handle:
        batches = json.loads(np.asarray(handle['metadata_json']).tobytes())['upstream_batches']
    rhos = {float(batch['rho_up']) for batch in batches}
    if rhos != {float(ARCHIVE_RHO_UP)}:
        raise ValueError('original upstream rho changed')
    return float(ARCHIVE_RHO_UP)


def load_saved_finite_cutoff_value(root, label=VALUE_LABEL):
    """Replay the saved nodal value residual. No propagator is evolved."""
    if label != VALUE_LABEL:
        raise ValueError('saved finite-cutoff value label changed')
    root = _lab_root(root)
    family, identity = require_original_history(root)
    if identity != CURRENT_HISTORY_IDENTITY:
        raise ValueError('original history 0b0e4ced required')
    record_path = root / VALUE_RECORD
    record = json.loads(record_path.read_text())
    if record.get('profile_identity') != identity or record.get('history_path') != HISTORY_RECORD:
        raise ValueError('original history 0b0e4ced required')
    if int(record.get('node_count', -1)) not in TARGET_COUNTS:
        raise ValueError('exact target nodes required')
    if record.get('tangents') not in (None, 'zero') and record.get('schema') != 'NSC-KS-GATE-VALUE-v2':
        raise ValueError('saved value tangents are not zero')
    nodes = np.asarray(target_nodes(family, int(record['node_count'])), dtype=float)
    if nodes.size != int(record['node_count']) or float(nodes[len(nodes) // 2]) != float(family.center):
        raise ValueError('exact target nodes required')
    directory = root / 'results/development/artifacts' / ('nsc-ks-gate-value-' + label)
    entries = []
    for path in directory.glob('family-*.json'):
        meta = json.loads(path.read_text())
        entries.append((tuple(meta['positive_family']), path, meta))
    entries.sort()
    if len(entries) != 60:
        raise ValueError('sixty saved finite-cutoff families required')
    matter = np.zeros((nodes.size, 2), dtype=float)
    family_records = {}
    source_identities = []
    for key, path, meta in entries:
        if meta.get('profile_identity') != identity:
            raise ValueError('original history 0b0e4ced required')
        if meta.get('evaluation_identity') != record.get('evaluation_identity'):
            raise ValueError('saved value family is not this evaluation')
        if meta.get('tangents') != 'zero' or not _same_solver(meta.get('solver'), record.get('solver')):
            raise ValueError('saved value settings changed')
        if list(meta.get('positive_family')) != [int(key[0]), int(key[1])]:
            raise ValueError('node-profile mismatch')
        arrays = authenticate_payload(meta, path.with_suffix('.npz'), root, target=nodes)
        if not np.array_equal(arrays['z'], nodes):
            raise ValueError('node-profile mismatch')
        matter += arrays['matter_change']
        stored = arrays.get('metadata') or {}
        overlap = set(family_records) & set(stored.get('family_records', ()))
        if overlap:
            raise ValueError('signed family counted twice')
        family_records.update(stored.get('family_records', {}))
        source_identities.append(meta['source_identity'])
    if len(family_records) != 120:
        raise ValueError('one hundred twenty explicit signed families required')
    weights = signed_family_angular_square_weights(family_records)
    if not np.allclose(weights, EXPECTED_ANGULAR_SQUARE_WEIGHT, rtol=0, atol=1e-6):
        raise ValueError('signed angular-square weights must be 107952 for each energy sign')
    coeff = evolved.coefficients_from_records(
        json.loads((root / SURFACE_COEFFICIENTS).read_text()),
        json.loads((root / SURFACE_BRANCH).read_text()))
    baseline = np.asarray(
        json.loads((root / SOURCE_UPDATE).read_text())['baseline']['action_gradient_approximant'],
        dtype=float)
    if baseline.shape != (2,):
        raise ValueError('saved baseline is not the original center pair')
    geometry = geometry_change(family, nodes, coeff)
    edge = edge_change(
        family, nodes, weights, coeff['a'], _retained_rho_up(root),
        int(record['solver']['phase_nodes']))
    # Baseline broadcast belongs to the saved value owner. The cubic term is
    # added later, only as a nodal array, and never through this path.
    residual = assemble_residual(baseline, geometry, matter, edge)
    if merit(residual) != record['merit'] or np.max(np.abs(residual), axis=0).tolist() != record['constraint_maxima']:
        raise ValueError('saved finite-cutoff value replay changed')
    digest = sha256()
    for item in sorted(source_identities):
        digest.update(str(item).encode())
    return {
        'residual': residual,
        'matter': matter,
        'nodes': nodes,
        'merit': merit(residual),
        'register_merit': record['merit'],
        'constraint_maxima': list(record['constraint_maxima']),
        'profile_identity': identity,
        'node_count': int(record['node_count']),
        'label': label,
        'rho_up': _retained_rho_up(root),
        'history_path': HISTORY_RECORD,
        'evaluation_identity': record['evaluation_identity'],
        'source_identity_digest': digest.hexdigest(),
        'source_identity_count': len(set(source_identities)),
        'family_count': 60,
        'signed_family_count': 120,
        'weights': np.asarray(weights, dtype=float),
        'solver': dict(record['solver']),
        'center_index': int(nodes.size // 2),
        'record_sha256': _digest(record_path),
    }
