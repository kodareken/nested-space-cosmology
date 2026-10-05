"""Physical rho=1 source-column successor: split windows, traced construction, directed maps.

Successor of the v1 coverage owner. This module loads the stored three-column
fibers, splits the 59 straddling low/subgap panels at each authenticated V5
joined threshold, and builds row/panel manifests for the complete 164-panel /
49372-row column-accuracy universe. Axial values at rho=1 are the PDE
restriction G(1) F of those stored columns. Directed 2-norm Lipschitz maps
convert a future coherent fiber error into N,beta majorants for
baseline_low_subgap and the physical upstream slot.

V5 action certificates are not column reconstruction. Stored SVD
amplification, ODE tolerances and the historical 4.998e-9 aggregate are not
enclosures. A finite N,beta component is emitted only when every inventory
identity carries an enclosure payload. The missing primitive is a directed
a-posteriori enclosure of the stored rho=1 fibers against the exact physical
incoming modes for each construction method. Historical v1 bytes stay
immutable. Gate 2 is not closed here.
"""
from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

from flint import arb, ctx, fmpq
import numpy as np

from .nsc_ks_residual_error import nonnegative
from .nsc_ks_rho1_source_bound import (
    AT_OR_ABOVE_JOINED_ROWS_IN_REMAINING_PANELS,
    BELOW_JOINED_ROWS_IN_REMAINING_PANELS,
    FORBIDDEN_BOUND_METHODS as V1_FORBIDDEN_BOUND_METHODS,
    HISTORICAL_COARSE_AGGREGATE,
    INVENTORY_POSITIVE_PANELS,
    INVENTORY_POSITIVE_ROWS,
    REMAINING_RECONSTRUCTION_PANELS,
    REMAINING_RECONSTRUCTION_ROWS,
    SCHEMA as SCHEMA_V1,
    STRADDLING_REMAINING_PANELS,
    V5_ACTION_COVERED_PANELS,
    V5_ACTION_COVERED_ROWS,
    assert_ell0_baseline_retained,
    audit_row_threshold_partition,
    coherent_covariance_required,
    coverage_digest,
    coverage_identity,
    reject_historical_coarse,
    signed_channel_multiplicity,
)
from .nsc_ks_upstream_error import (
    ELL0_GROUPS,
    IDENTITIES_NOT_ERROR_BOUNDS,
    REMAINING_N_AFTER_KNOWN,
    combine_upstream_row_error,
    nine_component_upstream_feed,
    pack_upper,
    packed_upper_fraction,
    restore_packed_upper,
)
from .nsc_mode_resolved_cauchy_state import deterministic_npz_bytes
from .nsc_retarded_radial_response import ks_generator
from .nsc_transmitting_dirac_domain import S3


SCHEMA = 'NSC-KS-RHO1-SOURCE-BOUND-v2'
V1_RECORD = 'results/development/nsc-ks-rho1-source-bound-v1.json'
SPLIT_WINDOW_COUNT = 223
STRADDLING_BELOW_WINDOWS = 59
STRADDLING_ABOVE_WINDOWS = 59
NONSTRADDLING_REMAINING_WINDOWS = 43
WORKPREC = 80
HERMITIAN_TOLERANCE = 3e-13
SIDE_BELOW = 'below_joined'
SIDE_ABOVE_UNCOVERED = 'at_or_above_joined_uncovered'
SIDE_V5_ACTION = 'v5_action_covered'
SIDE_CODES = {SIDE_BELOW: 0, SIDE_ABOVE_UNCOVERED: 1, SIDE_V5_ACTION: 2}
KIND_CODES = {
    'low_mode_reconstruction': 0,
    'subgap_rows': 1,
    'middle_mode_columns': 2,
}
CONSTRUCTION_OWNERS = {
    'low_mode_reconstruction': {
        'method': 'accepted physical low-mode reconstruction',
        'stored_field': 'recovery mode_at_one, then restrict_resolved_modes at rho=1',
        'owner': 'src/recursive_horizons/nsc_pg_mode_recovery.py',
        'restriction_owner': 'src/recursive_horizons/nsc_common_ks_trace.py',
        'reconstruction_enclosure': None,
        'why_open': (
            'packet-observation SVD amplification is a numerical sensitivity '
            'diagnostic, not a directed a-posteriori enclosure of the stored '
            'three-column fiber against the exact physical incoming mode'),
    },
    'middle_mode_columns': {
        'method': 'same original middle producer as the incoming baseline',
        'stored_field': 'archived_middle_boundary_modes.mode_at_one, then restrict_resolved_modes at rho=1',
        'owner': 'src/recursive_horizons/nsc_pg_archived_high_energy_modes.py',
        'restriction_owner': 'src/recursive_horizons/nsc_common_ks_trace.py',
        'reconstruction_enclosure': None,
        'why_open': (
            'the archived order-16 middle producer matches the incoming baseline '
            'and carries a V5 action certificate; that action bound is not a '
            'column-reconstruction enclosure, and no middle-Riccati defect '
            'majorant is bound here'),
    },
    'subgap_rows': {
        'method': 'eight inherited real subgap rows; accuracy OPEN',
        'stored_field': 'mode-resolved Cauchy interior at rho=1, eight rows per signed channel',
        'owner': 'src/recursive_horizons/nsc_mode_resolved_cauchy_state.py',
        'restriction_owner': 'src/recursive_horizons/nsc_common_ks_trace.py',
        'reconstruction_enclosure': None,
        'why_open': (
            'the eight inherited real subgap rows are stored without a directed '
            'construction enclosure; Pauli covariance reconstruction is an '
            'identity, not that enclosure'),
    },
}
MISSING_PRIMITIVE = (
    'directed a-posteriori enclosure of stored rho=1 three-column fibers and '
    'their PDE axial restriction G(1)F against the exact physical incoming '
    'modes, for each construction method: accepted physical low-mode SVD '
    'recovery from packet observations, the archived middle Riccati producer, '
    'and the eight inherited real subgap rows. Stored SVD amplification, ODE '
    'tolerances and V5 action certificates are not that enclosure.'
)
SUCCESSOR_STATUS = (
    'OPEN: authenticated 164-panel/49372-row manifests with 59 explicit V5-'
    'threshold splits; directed generator/axial Lipschitz and N,beta '
    'propagation represented; reconstruction enclosure of stored rho=1 fibers '
    'against exact physical incoming modes remains missing; V5 action is not '
    'column error; directed N,beta budget not filled; Gate 2 not closed'
)
FORBIDDEN_BOUND_METHODS = V1_FORBIDDEN_BOUND_METHODS + (
    'v5_action_as_column',
    'action_bound_as_column_reconstruction',
    'ode_tolerance_as_bound',
    'svd_amplification_as_bound',
    'count_only_closure',
    'zero_filled_rows',
)
AXIAL_PDE_RESTRICTION = 'F_z = G(rho=1) F with homogeneous ks_generator; source energy is a label'


def _sha_file(root, path):
    return sha256((Path(root) / path).read_bytes()).hexdigest()


def reject_forbidden_bound_method(method):
    if method in FORBIDDEN_BOUND_METHODS:
        raise ValueError('forbidden bound method: ' + str(method))
    if method is None or not isinstance(method, str) or not method:
        raise ValueError('explicit reconstruction method required')
    return method


def reject_action_as_column_substitution(bound, v5_action_bound):
    """V5 action numbers cannot be relabelled as column reconstruction error."""
    if bound is None:
        return None
    number = nonnegative(bound, 'column reconstruction bound')
    if v5_action_bound is None:
        return number
    values = v5_action_bound
    if isinstance(values, dict):
        values = [values.get('N'), values.get('beta')]
    for item in values if isinstance(values, (list, tuple)) else (values,):
        if item is None:
            continue
        action = nonnegative(item, 'V5 action bound')
        if number == action or abs(float(number) - float(action)) <= 5e-16 * max(1.0, float(action)):
            raise ValueError('V5 action bound is not a rho=1 column reconstruction bound')
    return number


def attach_reconstruction_bound(entry, bound, method, *, v5_action_bound=None):
    reject_forbidden_bound_method(method)
    if bound is None:
        updated = dict(entry)
        updated['reconstruction_error_upper'] = None
        updated['method'] = method
        return updated
    reject_historical_coarse(bound)
    reject_action_as_column_substitution(bound, v5_action_bound)
    updated = dict(entry)
    updated['reconstruction_error_upper'] = nonnegative(bound, 'reconstruction error')
    updated['method'] = method
    return updated


def window_identity(group, sign, panel, side, n_rows=None):
    if side not in SIDE_CODES:
        raise ValueError('split-window side must be below, uncovered-above, or V5-action')
    key = coverage_identity(group, sign, panel) + (str(side),)
    if n_rows is None:
        return key
    if isinstance(n_rows, bool) or int(n_rows) != n_rows or int(n_rows) < 1:
        raise ValueError('positive selected row count required')
    return key + (int(n_rows),)


def assert_v1_immutable(root):
    """Historical v1 record and its hashed owners/inputs are not rewritten."""
    path = Path(root) / V1_RECORD
    record = json.loads(path.read_text())
    if record.get('schema') != SCHEMA_V1:
        raise ValueError('historical v1 source-bound schema is immutable')
    if not str(record.get('status', '')).startswith('OPEN'):
        raise ValueError('historical v1 source-bound core must remain OPEN')
    if record.get('remaining_physical_rho1_source_error_bound') is not None:
        raise ValueError('historical v1 may not grow a physical rho=1 bound')
    for rel, expected in {**record['source_hashes'], **record['input_hashes']}.items():
        digest = _sha_file(root, rel)
        if digest != expected:
            raise ValueError('historical v1 source-bound bytes changed: ' + rel)
    return sha256(path.read_bytes()).hexdigest()


def rho1_chart(bits=WORKPREC):
    """Exact rho=1 KS chart: r^2=2, a^2=3*pi/2-4."""
    if isinstance(bits, bool) or int(bits) != bits or int(bits) < 16:
        raise ValueError('positive enclosure precision required')
    with ctx.workprec(int(bits)):
        a_sq = (3 * arb.pi()) / 2 - 4
        r_sq = arb(2)
        if not a_sq > 0:
            raise ArithmeticError('positive KS axial scale squared at rho=1 required')
        a = a_sq.sqrt()
        r = r_sq.sqrt()
        if not a.is_finite() or not r.is_finite() or not a > 0 or not r > 0:
            raise ArithmeticError('finite positive rho=1 chart enclosure required')
        return {'a': a, 'r': r, 'a_sq': a_sq, 'r_sq': r_sq, 'rho': arb(1)}


def pack_lower(value):
    """Exact dyadic lower: negate the directed upper of -value."""
    packed = pack_upper(-value)
    return {'mantissa': str(-int(packed['mantissa'])), 'exponent': int(packed['exponent'])}


def arb_nonneg(value, name):
    number = nonnegative(value, name)
    return arb(fmpq(number.numerator, number.denominator))


def frobenius_enclosure(matrix, bits=WORKPREC):
    """Directed enclosure of the Frobenius norm of a stored binary64 fiber."""
    block = np.asarray(matrix, complex)
    if block.size == 0 or not np.isfinite(block.real).all() or not np.isfinite(block.imag).all():
        raise ValueError('finite source fiber required')
    with ctx.workprec(int(bits)):
        acc = arb(0)
        for entry in block.ravel():
            real = arb(float(entry.real))
            imag = arb(float(entry.imag))
            acc += real * real + imag * imag
        norm = acc.sqrt()
        if not norm.is_finite():
            raise ArithmeticError('finite fiber Frobenius enclosure required')
        return norm


def generator_two_norm_enclosure(energy, mass, angular, chart, bits=WORKPREC):
    """||G(1)||_2 = sqrt(m^2+(ell/r)^2+(E/a)^2)/a for the Pauli KS generator."""
    with ctx.workprec(int(bits)):
        mass_b = arb(float(mass))
        ell = arb(abs(float(angular)))
        energy_b = arb(abs(float(energy)))
        if not mass_b >= 0 or not ell >= 0:
            raise ValueError('nonnegative mass and absolute angular label required')
        a, r = chart['a'], chart['r']
        norm = ((mass_b * mass_b + (ell / r) ** 2 + (energy_b / a) ** 2).sqrt()) / a
        if not norm.is_finite() or not norm >= 0:
            raise ArithmeticError('finite generator 2-norm enclosure required')
        return norm


def weight_factor_enclosure(weight, bits=WORKPREC):
    """sqrt(w/(2 pi)) from the original quadrature weight, used once."""
    number = float(weight)
    if not np.isfinite(number) or number <= 0:
        raise ValueError('positive original quadrature weight required')
    with ctx.workprec(int(bits)):
        factor = (arb(number) / (2 * arb.pi())).sqrt()
        if not factor.is_finite() or not factor > 0:
            raise ArithmeticError('finite quadrature weight factor required')
        return factor


def pauli_generator_matches_ks(energy, mass, angular, *, bits=WORKPREC):
    """Numeric 2-norm of ks_generator must lie in the Pauli enclosure."""
    matrix = ks_generator(1.0, np.array([float(energy)]), float(mass), float(angular))[0]
    numeric = float(np.linalg.norm(matrix, 2))
    chart = rho1_chart(bits=bits)
    enclosure = generator_two_norm_enclosure(energy, mass, angular, chart, bits=bits)
    lower = float(restore_packed_upper(pack_lower(enclosure)))
    upper = float(restore_packed_upper(pack_upper(enclosure)))
    if numeric < lower or numeric > upper:
        raise ValueError('Pauli generator 2-norm enclosure missed ks_generator')
    return numeric, lower, upper


def _row_kind(method):
    if method.startswith('accepted physical low-mode'):
        return 'low_mode_reconstruction'
    if method.startswith('same original middle producer'):
        return 'middle_mode_columns'
    if method.startswith('eight inherited real subgap'):
        return 'subgap_rows'
    raise ValueError('untraced rho=1 source-column construction method')


def split_straddling_windows(arrays, coverage):
    """Split remaining panels at the authenticated V5 joined threshold."""
    if not isinstance(coverage, dict):
        raise TypeError('rho=1 coverage mapping required')
    windows = []
    straddling = 0
    below_rows = 0
    above_uncovered_rows = 0
    for item in coverage['remaining_windows']:
        energies = np.asarray(arrays[item['panel'] + '/energies'], float)
        if (energies.ndim != 1 or len(energies) != int(item['n_rows'])
                or not np.isfinite(energies).all()):
            raise ValueError('authenticated source-panel energy rows changed')
        joined = float(item['v5_joined'][0])
        below = energies < joined
        n_below = int(np.count_nonzero(below))
        n_above = int(len(energies) - n_below)
        kind = item['kind']
        if kind not in KIND_CODES:
            raise ValueError('remaining window has an untraced construction kind')
        base = {
            'group': int(item['group']),
            'sign': int(item['sign']),
            'panel': item['panel'],
            'kind': kind,
            'v5_joined': list(item['v5_joined']),
            'radius_history_response': item['radius_history_response'],
            'baseline_retained': True,
            'v5_action_is_not_column_reconstruction': True,
            'source_reconstruction': 'unbounded',
            'error_bound': None,
        }
        if n_below and n_above:
            straddling += 1
            below_idx = np.flatnonzero(below)
            above_idx = np.flatnonzero(~below)
            windows.append({
                **base, 'side': SIDE_BELOW, 'n_rows': n_below,
                'row_indices': [int(i) for i in below_idx],
                'energy_interval': [float(energies[below_idx].min()),
                                    float(energies[below_idx].max())],
                'straddling_parent': True,
            })
            windows.append({
                **base, 'side': SIDE_ABOVE_UNCOVERED, 'n_rows': n_above,
                'row_indices': [int(i) for i in above_idx],
                'energy_interval': [float(energies[above_idx].min()),
                                    float(energies[above_idx].max())],
                'straddling_parent': True,
            })
            below_rows += n_below
            above_uncovered_rows += n_above
        elif n_below:
            idx = np.arange(len(energies), dtype=int)
            windows.append({
                **base, 'side': SIDE_BELOW, 'n_rows': n_below,
                'row_indices': [int(i) for i in idx],
                'energy_interval': [float(energies.min()), float(energies.max())],
                'straddling_parent': False,
            })
            below_rows += n_below
        else:
            raise ValueError(
                'remaining reconstruction panel has no rows below its V5 joined threshold')
    if straddling != STRADDLING_REMAINING_PANELS:
        raise ValueError('straddling remaining panel count changed')
    if below_rows != BELOW_JOINED_ROWS_IN_REMAINING_PANELS:
        raise ValueError('below-joined row count changed after the V5 split')
    if above_uncovered_rows != AT_OR_ABOVE_JOINED_ROWS_IN_REMAINING_PANELS:
        raise ValueError('uncovered at-or-above row count changed after the V5 split')
    v5_rows = 0
    for item in coverage['v5_action_windows']:
        energies = np.asarray(arrays[item['panel'] + '/energies'], float)
        if (energies.ndim != 1 or len(energies) != int(item['n_rows'])
                or not np.isfinite(energies).all()):
            raise ValueError('authenticated V5 action energy rows changed')
        joined = float(item['v5_joined'][0])
        if np.any(energies < joined):
            raise ValueError('V5 action panel contains rows below its joined threshold')
        kind = item['kind']
        if kind not in KIND_CODES:
            raise ValueError('V5 action window has an untraced construction kind')
        idx = np.arange(len(energies), dtype=int)
        windows.append({
            'group': int(item['group']),
            'sign': int(item['sign']),
            'panel': item['panel'],
            'kind': kind,
            'side': SIDE_V5_ACTION,
            'n_rows': int(len(energies)),
            'row_indices': [int(i) for i in idx],
            'energy_interval': [float(energies.min()), float(energies.max())],
            'v5_joined': list(item['v5_joined']),
            'radius_history_response': item['radius_history_response'],
            'baseline_retained': True,
            'straddling_parent': False,
            'v5_action_is_not_column_reconstruction': True,
            'source_reconstruction': 'unbounded',
            'error_bound': None,
            'v5_action_covered_panel': True,
        })
        v5_rows += int(len(energies))
    if v5_rows != V5_ACTION_COVERED_ROWS:
        raise ValueError('V5 action-covered row count changed')
    if len(windows) != SPLIT_WINDOW_COUNT:
        raise ValueError('split-window count changed')
    identities = [
        list(window_identity(w['group'], w['sign'], w['panel'], w['side'], w['n_rows']))
        for w in windows]
    if len(identities) != len({tuple(item) for item in identities}):
        raise ValueError('duplicate split-window identity')
    by_side = {}
    for window in windows:
        by_side[window['side']] = by_side.get(window['side'], 0) + 1
    if (by_side.get(SIDE_BELOW) != STRADDLING_BELOW_WINDOWS + NONSTRADDLING_REMAINING_WINDOWS
            or by_side.get(SIDE_ABOVE_UNCOVERED) != STRADDLING_ABOVE_WINDOWS
            or by_side.get(SIDE_V5_ACTION) != V5_ACTION_COVERED_PANELS):
        raise ValueError('split-window side partition changed')
    return {
        'windows': windows,
        'split_window_count': len(windows),
        'straddling_remaining_panels': straddling,
        'straddling_below_windows': STRADDLING_BELOW_WINDOWS,
        'straddling_above_uncovered_windows': STRADDLING_ABOVE_WINDOWS,
        'nonstraddling_remaining_windows': NONSTRADDLING_REMAINING_WINDOWS,
        'v5_action_windows': V5_ACTION_COVERED_PANELS,
        'below_joined_rows': below_rows,
        'at_or_above_joined_uncovered_rows': above_uncovered_rows,
        'v5_action_rows': v5_rows,
        'window_identities': identities,
        'window_identity_digest': coverage_digest(identities),
        'v5_action_is_not_column_reconstruction': True,
        'straddling_panels_were_split': True,
    }


def reject_straddling_misclassification(arrays, coverage, split):
    """A straddling parent must not be kept as a single reconstruction window."""
    remaining = {item['panel']: item for item in coverage['remaining_windows']}
    split_by_panel = {}
    for window in split['windows']:
        split_by_panel.setdefault(window['panel'], []).append(window)
    n_straddle = 0
    for name, item in remaining.items():
        energies = np.asarray(arrays[name + '/energies'], float)
        joined = float(item['v5_joined'][0])
        n_below = int(np.count_nonzero(energies < joined))
        n_above = int(len(energies) - n_below)
        pieces = split_by_panel[name]
        if any(piece['side'] == SIDE_V5_ACTION for piece in pieces):
            raise ValueError('below-threshold rows cannot inherit V5 action coverage')
        if n_below and n_above:
            n_straddle += 1
            sides = {piece['side'] for piece in pieces}
            if sides != {SIDE_BELOW, SIDE_ABOVE_UNCOVERED}:
                raise ValueError('straddling panel was not split at its V5 joined threshold')
            rows = {piece['side']: piece['n_rows'] for piece in pieces}
            if rows[SIDE_BELOW] != n_below or rows[SIDE_ABOVE_UNCOVERED] != n_above:
                raise ValueError('straddling split row counts do not match the joined threshold')
        elif len(pieces) != 1 or pieces[0]['side'] != SIDE_BELOW:
            raise ValueError('non-straddling remaining panel must stay a below-joined window')
    if n_straddle != STRADDLING_REMAINING_PANELS:
        raise ValueError('straddling remaining panel count changed')
    for window in split['windows']:
        if window['side'] == SIDE_V5_ACTION and window['panel'] in remaining:
            raise ValueError('remaining reconstruction panel cannot be classified as V5 action')
    return True


def _minus_energy_isometry(amplitudes):
    mapped = np.einsum('ab,nbc->nac', S3, np.asarray(amplitudes, complex).conj())
    original = np.linalg.norm(amplitudes, axis=(1, 2))
    partner = np.linalg.norm(mapped, axis=(1, 2))
    if float(np.max(np.abs(original - partner))) > 0:
        raise ValueError('E/-E partner must preserve the Frobenius fiber norm')
    return mapped


def _multiplicity_arb(channel):
    per_sign, signs = signed_channel_multiplicity(channel)
    number = Fraction(per_sign).limit_denominator()
    nonnegative(number, 'signed-family multiplicity')
    if signs not in (1, 2):
        raise ValueError('one or two actual angular signs required')
    return arb_nonneg(number, 'signed-family multiplicity'), int(signs)


def row_propagation_coefficients(energy, weight, amplitudes, *, mass, angular,
                                 multiplicity, chart, bits=WORKPREC):
    """Directed linear/quadratic map from weighted fiber error to N,beta.

    The three columns stay a coherent Frobenius fiber. The minus-energy
    partner is the owned antiunitary isometry with original weights, counted
    as a second sector rather than a folded multiplicity. Covariance is the
    exact source law at the energy label (operator norm <= 1, ec=0).
    """
    f_unweighted = frobenius_enclosure(amplitudes, bits=bits)
    g = generator_two_norm_enclosure(energy, mass, angular, chart, bits=bits)
    factor = weight_factor_enclosure(weight, bits=bits)
    mu, _signs = (multiplicity if isinstance(multiplicity, tuple)
                  else (multiplicity, None))
    if not isinstance(mu, arb):
        mu = arb_nonneg(mu, 'signed-family multiplicity')
    with ctx.workprec(int(bits)):
        f = factor * f_unweighted
        mass_b = arb(float(mass))
        ell = arb(abs(float(angular)))
        a, r = chart['a'], chart['r']
        geom = mass_b + ell / r + g / a
        # Two isometric energy sectors, original weights, coherent Frobenius.
        linear_n = 4 * mu * f * geom
        quadratic_n = 6 * mu * geom
        linear_beta = 4 * mu * f * g
        quadratic_beta = 6 * mu * g
        values = (g, f, linear_n, linear_beta, quadratic_n, quadratic_beta)
        if any(not item.is_finite() or not item >= 0 for item in values):
            raise ArithmeticError('finite nonnegative N,beta Lipschitz enclosure required')
        return {
            'G_two_norm_upper': pack_upper(g),
            'weighted_fiber_norm_upper': pack_upper(f),
            'linear_N_upper': pack_upper(linear_n),
            'linear_beta_upper': pack_upper(linear_beta),
            'quadratic_N_upper': pack_upper(quadratic_n),
            'quadratic_beta_upper': pack_upper(quadratic_beta),
            'minus_energy_factor': 2,
            'minus_energy_retained': True,
            'weights_retained': True,
            'coherence_retained': True,
            'degeneracy_retained': True,
            'independent_column_max_used': False,
            'axial_restriction': AXIAL_PDE_RESTRICTION,
        }


def _zeros_packed(n):
    packed = pack_upper(arb(0))
    mantissa = np.full(n, packed['mantissa'].encode(), dtype='S48')
    exponent = np.full(n, int(packed['exponent']), dtype=np.int32)
    return mantissa, exponent


def trace_physical_rho1_universe(arrays, meta, split, channels, *, bits=WORKPREC):
    """One pass over every stored fiber: construction, split side, Lipschitz."""
    panels = meta['panels']
    names = sorted(panels, key=lambda name: (
        int(panels[name]['group']), int(panels[name]['angular_sign']), name))
    if len(names) != INVENTORY_POSITIVE_PANELS:
        raise ValueError('inventory panel count changed')
    name_index = {name: index for index, name in enumerate(names)}
    n_rows = INVENTORY_POSITIVE_ROWS
    row_panel = np.empty(n_rows, np.int32)
    row_index = np.empty(n_rows, np.int32)
    row_group = np.empty(n_rows, np.int16)
    row_sign = np.empty(n_rows, np.int8)
    row_side = np.empty(n_rows, np.int8)
    row_kind = np.empty(n_rows, np.int8)
    row_energy = np.empty(n_rows, np.float64)
    row_weight = np.empty(n_rows, np.float64)
    g_man, g_exp = _zeros_packed(n_rows)
    f_man, f_exp = _zeros_packed(n_rows)
    ln_man, ln_exp = _zeros_packed(n_rows)
    lb_man, lb_exp = _zeros_packed(n_rows)
    qn_man, qn_exp = _zeros_packed(n_rows)
    qb_man, qb_exp = _zeros_packed(n_rows)
    identities = []
    census = {kind: {'panels': 0, 'rows': 0, 'straddling_panels': 0}
              for kind in KIND_CODES}
    seen_panels = {kind: set() for kind in KIND_CODES}
    chart = rho1_chart(bits=bits)
    g_cache = {}
    cursor = 0
    with ctx.workprec(int(bits)):
        sum_lin_n = arb(0)
        sum_lin_b = arb(0)
        sum_quad_n = arb(0)
        sum_quad_b = arb(0)
        sum_f = arb(0)
        side_linear_n = {side: arb(0) for side in SIDE_CODES}
        for window in split['windows']:
            name = window['panel']
            desc = panels[name]
            group = int(window['group'])
            sign = int(window['sign'])
            if int(desc['group']) != group or int(desc['angular_sign']) != sign:
                raise ValueError('panel family/channel mismatch')
            kind = window['kind']
            traced = _row_kind(desc['provenance']['method'])
            if traced != kind:
                raise ValueError('stored construction method does not match the coverage kind')
            assert_ell0_baseline_retained(
                group, baseline_retained=True, panel=name,
                radius_response=window['radius_history_response'])
            energies = np.asarray(arrays[name + '/energies'], float)
            weights = np.asarray(arrays[name + '/weights'], float)
            amplitudes = np.asarray(arrays[name + '/amplitudes_at_one'], complex)
            covariance = np.asarray(arrays[name + '/covariance'], complex)
            coherent_covariance_required(covariance)
            if np.any(weights <= 0) or np.any(energies == 0):
                raise ValueError('complete nonzero real quadrature rows required')
            _minus_energy_isometry(amplitudes)
            channel = channels[group]
            mu, n_signs = _multiplicity_arb(channel)
            expected_signs = (1,) if float(channel['angular_eigenvalue']) == 0 else (-1, 1)
            if sign not in expected_signs or (n_signs == 1 and sign != 1):
                raise ValueError('panel sign is outside the retained channel inventory')
            mass = float(desc['mass'])
            angular = float(desc['angular_sign']) * float(desc['angular_magnitude'])
            for local in window['row_indices']:
                energy = float(energies[local])
                weight = float(weights[local])
                fiber = amplitudes[local]
                key = (mass, abs(angular), abs(energy))
                if key not in g_cache:
                    g_cache[key] = generator_two_norm_enclosure(
                        energy, mass, angular, chart, bits=bits)
                coeffs = row_propagation_coefficients(
                    energy, weight, fiber, mass=mass, angular=angular,
                    multiplicity=mu, chart=chart, bits=bits)
                if cursor >= n_rows:
                    raise ValueError('row manifest exceeded the inventory universe')
                row_panel[cursor] = name_index[name]
                row_index[cursor] = int(local)
                row_group[cursor] = group
                row_sign[cursor] = sign
                row_side[cursor] = SIDE_CODES[window['side']]
                row_kind[cursor] = KIND_CODES[kind]
                row_energy[cursor] = energy
                row_weight[cursor] = weight
                g_man[cursor] = coeffs['G_two_norm_upper']['mantissa'].encode()
                g_exp[cursor] = coeffs['G_two_norm_upper']['exponent']
                f_man[cursor] = coeffs['weighted_fiber_norm_upper']['mantissa'].encode()
                f_exp[cursor] = coeffs['weighted_fiber_norm_upper']['exponent']
                ln_man[cursor] = coeffs['linear_N_upper']['mantissa'].encode()
                ln_exp[cursor] = coeffs['linear_N_upper']['exponent']
                lb_man[cursor] = coeffs['linear_beta_upper']['mantissa'].encode()
                lb_exp[cursor] = coeffs['linear_beta_upper']['exponent']
                qn_man[cursor] = coeffs['quadratic_N_upper']['mantissa'].encode()
                qn_exp[cursor] = coeffs['quadratic_N_upper']['exponent']
                qb_man[cursor] = coeffs['quadratic_beta_upper']['mantissa'].encode()
                qb_exp[cursor] = coeffs['quadratic_beta_upper']['exponent']
                lin_n = restore_packed_upper(coeffs['linear_N_upper'])
                lin_b = restore_packed_upper(coeffs['linear_beta_upper'])
                quad_n = restore_packed_upper(coeffs['quadratic_N_upper'])
                quad_b = restore_packed_upper(coeffs['quadratic_beta_upper'])
                f_w = restore_packed_upper(coeffs['weighted_fiber_norm_upper'])
                sum_lin_n += lin_n
                sum_lin_b += lin_b
                sum_quad_n += quad_n
                sum_quad_b += quad_b
                sum_f += f_w
                side_linear_n[window['side']] += lin_n
                identities.append(list(coverage_identity(group, sign, name, int(local))))
                cursor += 1
            seen_panels[kind].add(name)
            census[kind]['rows'] += int(window['n_rows'])
        if cursor != n_rows:
            raise ValueError('row manifest does not cover the inventory universe')
        if len(identities) != len({tuple(item) for item in identities}):
            raise ValueError('duplicate rho=1 source row identity')
        for kind in KIND_CODES:
            census[kind]['panels'] = len(seen_panels[kind])
        census['low_mode_reconstruction']['straddling_panels'] = STRADDLING_REMAINING_PANELS
        census['subgap_rows']['straddling_panels'] = 0
        census['middle_mode_columns']['straddling_panels'] = 0
        if (census['low_mode_reconstruction']['panels'] != 64
                or census['subgap_rows']['panels'] != 38
                or census['middle_mode_columns']['panels'] != V5_ACTION_COVERED_PANELS):
            raise ValueError('construction-method panel census changed')
        if (census['low_mode_reconstruction']['rows']
                + census['subgap_rows']['rows']
                + census['middle_mode_columns']['rows'] != n_rows):
            raise ValueError('construction-method row census does not partition the universe')
        remaining = arb(float(REMAINING_N_AFTER_KNOWN))
        threshold = remaining / sum_lin_n
        arrays_out = {
            'panel_names': np.array(names),
            'row_panel_index': row_panel,
            'row_index': row_index,
            'row_group': row_group,
            'row_sign': row_sign,
            'row_side': row_side,
            'row_kind': row_kind,
            'row_energy': row_energy,
            'row_weight': row_weight,
            'G_two_norm_mantissa': g_man,
            'G_two_norm_exponent': g_exp,
            'weighted_fiber_norm_mantissa': f_man,
            'weighted_fiber_norm_exponent': f_exp,
            'linear_N_mantissa': ln_man,
            'linear_N_exponent': ln_exp,
            'linear_beta_mantissa': lb_man,
            'linear_beta_exponent': lb_exp,
            'quadratic_N_mantissa': qn_man,
            'quadratic_N_exponent': qn_exp,
            'quadratic_beta_mantissa': qb_man,
            'quadratic_beta_exponent': qb_exp,
        }
        side_packed = {side: pack_upper(value) for side, value in side_linear_n.items()}
        return {
            'arrays': arrays_out,
            'identities': identities,
            'identity_digest': coverage_digest(identities),
            'n_rows': n_rows,
            'n_panels': len(names),
            'unique_generator_keys': len(g_cache),
            'construction_census': census,
            'propagation_sums': {
                'linear_N_upper': pack_upper(sum_lin_n),
                'linear_beta_upper': pack_upper(sum_lin_b),
                'quadratic_N_upper': pack_upper(sum_quad_n),
                'quadratic_beta_upper': pack_upper(sum_quad_b),
                'weighted_fiber_norm_sum_upper': pack_upper(sum_f),
                'linear_N_by_side': side_packed,
            },
            'linear_indicator_uniform_weighted_fiber_threshold': pack_upper(threshold),
            'linear_indicator_is_not_a_certificate': True,
            'minus_energy_factor': 2,
            'reconstruction_enclosures_evaluated_rows': 0,
            'propagation_coefficients_evaluated_rows': n_rows,
        }


def manifest_payload_bytes(trace, split, *, extra=None):
    arrays = dict(trace['arrays'])
    metadata = {
        'schema': SCHEMA,
        'identity_digest': trace['identity_digest'],
        'window_identity_digest': split['window_identity_digest'],
        'n_rows': trace['n_rows'],
        'n_panels': trace['n_panels'],
        'split_window_count': split['split_window_count'],
        'straddling_remaining_panels': split['straddling_remaining_panels'],
        'construction_census': trace['construction_census'],
        'propagation_sums': trace['propagation_sums'],
        'missing_primitive': MISSING_PRIMITIVE,
        'v5_action_is_not_column_reconstruction': True,
        'reconstruction_enclosures_evaluated_rows': 0,
    }
    if extra:
        metadata.update(extra)
    arrays['metadata_json'] = np.frombuffer(
        json.dumps(metadata, sort_keys=True, separators=(',', ':')).encode(), np.uint8)
    return deterministic_npz_bytes(arrays)


def combine_directed_n_beta(contributions, *, declared_universe_size,
                            universe_identities=None, universe_digest=None):
    """Sum authenticated row enclosures. Missing payloads keep the component OPEN."""
    if (isinstance(declared_universe_size, bool)
            or int(declared_universe_size) != declared_universe_size
            or int(declared_universe_size) < 1):
        raise ValueError('positive declared universe size required')
    universe = int(declared_universe_size)
    rows = list(contributions)
    if len(rows) != universe:
        raise ValueError('contribution count does not match the declared universe')
    identities = []
    finite = True
    with ctx.workprec(WORKPREC):
        total_n = arb(0)
        total_beta = arb(0)
        for row in rows:
            reject_forbidden_bound_method(row.get('method', ''))
            if not row.get('coherence_retained') or not row.get('weights_retained'):
                raise ValueError('covariance, coherence and weights must be retained')
            if not row.get('minus_energy_retained'):
                raise ValueError('minus-energy sector must be retained')
            if not row.get('degeneracy_retained', True):
                raise ValueError('channel degeneracy must be retained')
            if row.get('independent_column_max_used'):
                raise ValueError('independent maxima of correlated columns are forbidden')
            identity = tuple(row['identity'])
            identities.append(identity)
            payload = row.get('enclosure_payload_sha256')
            n_value, beta_value = row.get('N'), row.get('beta')
            if n_value is None or beta_value is None or not payload:
                finite = False
                if n_value == 0 and beta_value == 0 and not payload:
                    raise ValueError('zero-filled rows cannot close a physical source bound')
                continue
            if not isinstance(payload, str) or len(payload) != 64:
                raise ValueError('enclosure payload must be a sha256 hex digest')
            n_num = nonnegative(n_value, 'N contribution')
            b_num = nonnegative(beta_value, 'beta contribution')
            total_n += arb_nonneg(n_num, 'N contribution')
            total_beta += arb_nonneg(b_num, 'beta contribution')
    if len(identities) != len(set(identities)):
        raise ValueError('double-counted rho=1 source contribution')
    if universe_identities is not None:
        expected = tuple(tuple(item) for item in universe_identities)
        if tuple(identities) != expected:
            raise ValueError('contribution identities do not match the authenticated universe')
    digest = coverage_digest([list(item) for item in identities])
    if universe_digest is not None and digest != universe_digest:
        raise ValueError('contribution identity digest does not match the universe')
    if not finite:
        return {
            'N': None,
            'beta': None,
            'status': (
                'OPEN: directed N,beta requires an enclosure payload on every '
                'authenticated inventory identity'),
            'declared_universe_size': universe,
            'identity_digest': digest,
            'independent_column_max_used': False,
            'reconstruction_enclosures_bound': 0,
        }
    if not total_n.is_finite() or not total_beta.is_finite():
        raise ArithmeticError('finite combined N,beta enclosure required')
    return {
        'N': pack_upper(total_n),
        'beta': pack_upper(total_beta),
        'status': 'conditional: every identity carries an enclosure payload',
        'declared_universe_size': universe,
        'identity_digest': digest,
        'independent_column_max_used': False,
        'reconstruction_enclosures_bound': universe,
        'physical_EXISTENCE_certificate': False,
        'gate_2_closed': False,
    }


def correlated_budget_feeds(*, below_joined_n_beta=None,
                            at_or_above_uncovered_n_beta=None,
                            v5_action_column_n_beta=None,
                            continuation_n_beta=None,
                            enclosure_payload_sha256=None,
                            universe_digest=None):
    """Map column error onto baseline_low_subgap and the physical upstream slot.

    V5 action-covered panels remain column-accuracy OPEN. Their action bound
    is never copied into covered_regions from this owner.
    """
    pieces = (below_joined_n_beta, at_or_above_uncovered_n_beta, v5_action_column_n_beta)
    authenticated = bool(enclosure_payload_sha256) and bool(universe_digest)
    if any(piece is not None for piece in pieces) and not authenticated:
        raise ValueError(
            'physical rho=1 N,beta requires an enclosure payload bound to the universe digest')

    def _pair(name, value):
        if value is None:
            return None
        if not isinstance(value, (tuple, list)) or len(value) != 2:
            raise ValueError('N,beta pair required for ' + name)
        return (nonnegative(value[0], name + ' N'), nonnegative(value[1], name + ' beta'))

    below = _pair('below_joined', below_joined_n_beta)
    above = _pair('uncovered_above', at_or_above_uncovered_n_beta)
    v5_column = _pair('v5_action_column', v5_action_column_n_beta)
    physical = None
    if below is not None and above is not None and v5_column is not None:
        physical = (below[0] + above[0] + v5_column[0],
                    below[1] + above[1] + v5_column[1])
    upstream = combine_upstream_row_error(
        None if continuation_n_beta is None else continuation_n_beta[0],
        None if physical is None else physical[0])
    feed = nine_component_upstream_feed(
        continuation_n_beta,
        None if physical is None else [float(physical[0]), float(physical[1])])
    return {
        'baseline_low_subgap_N_beta': (
            None if below is None else [str(below[0]), str(below[1])]),
        'baseline_low_subgap_status': (
            'OPEN: below-joined reconstruction enclosure missing'
            if below is None else
            'conditional: below-joined enclosure payload bound to the universe'),
        'at_or_above_joined_uncovered_N_beta': (
            None if above is None else [str(above[0]), str(above[1])]),
        'v5_action_column_N_beta': (
            None if v5_column is None else [str(v5_column[0]), str(v5_column[1])]),
        'v5_action_is_not_column_reconstruction': True,
        'covered_regions_not_fed_from_column_error': True,
        'physical_upstream_N_beta': (
            None if physical is None else [str(physical[0]), str(physical[1])]),
        'upstream_row_error': upstream,
        'nine_component_feed': feed,
        'enclosure_payload_sha256': enclosure_payload_sha256,
        'universe_digest': universe_digest,
        'gate_2_closed': False,
    }


def physical_source_error_for_upstream(physical_n_beta, continuation_n_beta=None, *,
                                       enclosure_payload_sha256=None,
                                       universe_digest=None):
    """Feed the upstream owner. A missing enclosure keeps the physical slot OPEN."""
    if physical_n_beta is None:
        return correlated_budget_feeds(continuation_n_beta=continuation_n_beta)
    if not enclosure_payload_sha256 or not universe_digest:
        raise ValueError(
            'physical rho=1 N,beta requires an enclosure payload bound to the universe digest')
    pair = (nonnegative(physical_n_beta[0], 'physical N'),
            nonnegative(physical_n_beta[1], 'physical beta'))
    return correlated_budget_feeds(
        below_joined_n_beta=pair,
        at_or_above_uncovered_n_beta=(0, 0),
        v5_action_column_n_beta=(0, 0),
        continuation_n_beta=continuation_n_beta,
        enclosure_payload_sha256=enclosure_payload_sha256,
        universe_digest=universe_digest)


def resource_forecast(trace, *, generators_launched=0):
    census = trace['construction_census']
    return {
        'measured_rows_processed': trace['n_rows'],
        'measured_unique_generator_keys': trace['unique_generator_keys'],
        'generators_launched': int(generators_launched),
        'missing_primitive': MISSING_PRIMITIVE,
        'construction_row_counts': {
            kind: {'panels': item['panels'], 'rows': item['rows']}
            for kind, item in census.items()},
        'forecast': {
            'low_mode_observation_map_defect_enclosure': (
                'one 6-packet observation ODE on [-1,1] plus a directed SVD/'
                'defect enclosure per stored low-mode energy; '
                + str(census['low_mode_reconstruction']['rows'])
                + ' rows; not launched'),
            'middle_riccati_defect_enclosure': (
                'a-posteriori defect majorant of the archived order-16 middle '
                'producer at rho=1; '
                + str(census['middle_mode_columns']['rows'])
                + ' rows; not launched'),
            'subgap_construction_enclosure': (
                'directed construction enclosure of the eight inherited real '
                'subgap rows; '
                + str(census['subgap_rows']['rows'])
                + ' rows; not launched'),
        },
        'analogous_existing_tool': (
            'continue_columns Bernstein/defect majorant in nsc_ks_upstream_error, '
            'which treats stored binary64 rho=1 columns as exact only for the '
            'numerical-continuation sub-bound'),
        'partial_or_resource_stop_cannot_close_gate_2': True,
    }


def _jsonable_windows(split):
    out = []
    for window in split['windows']:
        item = {
            'group': window['group'],
            'sign': window['sign'],
            'panel': window['panel'],
            'kind': window['kind'],
            'side': window['side'],
            'n_rows': window['n_rows'],
            'energy_interval': list(window['energy_interval']),
            'v5_joined': list(window['v5_joined']),
            'radius_history_response': window['radius_history_response'],
            'straddling_parent': bool(window.get('straddling_parent')),
            'source_reconstruction': 'unbounded',
            'error_bound': None,
            'v5_action_is_not_column_reconstruction': True,
            'baseline_retained': True,
        }
        out.append(item)
    return out


def validate_rho1_source_bound_v2_record(record):
    if record.get('schema') != SCHEMA:
        raise ValueError('unexpected rho=1 source-bound v2 schema')
    status = record.get('status')
    if not isinstance(status, str) or not status.startswith('OPEN'):
        raise ValueError('rho=1 source-bound v2 must remain OPEN')
    if record.get('missing_primitive') != MISSING_PRIMITIVE:
        raise ValueError('v2 must name the missing reconstruction-enclosure primitive')
    coverage = record['coverage']
    if coverage['inventory_positive_panels'] != INVENTORY_POSITIVE_PANELS:
        raise ValueError('inventory panel count changed')
    if coverage['inventory_positive_rows'] != INVENTORY_POSITIVE_ROWS:
        raise ValueError('inventory row count changed')
    if coverage['v5_action_covered_panels'] != V5_ACTION_COVERED_PANELS:
        raise ValueError('V5 action-covered panel count changed')
    if coverage['straddling_remaining_panels'] != STRADDLING_REMAINING_PANELS:
        raise ValueError('straddling remaining panel count changed')
    if coverage['split_window_count'] != SPLIT_WINDOW_COUNT:
        raise ValueError('split-window count changed')
    if not coverage.get('straddling_panels_were_split'):
        raise ValueError('v2 must split the 59 straddling panels')
    if not coverage.get('v5_action_is_not_column_reconstruction'):
        raise ValueError('V5 action coverage is not column reconstruction')
    if coverage.get('double_counted'):
        raise ValueError('coverage claims a double count')
    partition = record['row_threshold_partition']
    if partition['straddling_remaining_panels'] != STRADDLING_REMAINING_PANELS:
        raise ValueError('row-level V5 threshold partition changed')
    if partition['physical_column_accuracy_universe_rows'] != INVENTORY_POSITIVE_ROWS:
        raise ValueError('physical column-accuracy universe changed')
    split = record['split']
    if split['straddling_remaining_panels'] != STRADDLING_REMAINING_PANELS:
        raise ValueError('split straddling count changed')
    if coverage_digest(split['window_identities']) != split['window_identity_digest']:
        raise ValueError('split-window identity digest changed')
    if record.get('remaining_physical_rho1_source_error_bound') is not None:
        raise ValueError('v2 may not invent a missing physical rho=1 reconstruction bound')
    directed = record['directed_N_beta_contribution_bound']
    if directed.get('N') is not None or directed.get('beta') is not None:
        raise ValueError('v2 production record has no finite directed N,beta reconstruction budget')
    if directed.get('reconstruction_enclosures_bound', 0) != 0:
        raise ValueError('v2 production record evaluated no reconstruction enclosures')
    if record.get('historical_coarse_4p998e_minus9_is_current_remaining_bound'):
        raise ValueError('historical 4.998e-9 coarse aggregate is not the current remaining bound')
    if record.get('physical_EXISTENCE_certificate') or record.get('physical_NONEXISTENCE_certificate'):
        raise ValueError('OPEN source-bound successor cannot issue a physical certificate')
    if record.get('gate_2_closed'):
        raise ValueError('a partial source-accuracy successor cannot close Gate 2')
    if record.get('v1_rewritten'):
        raise ValueError('historical v1 bytes are immutable')
    sample = record['sample']
    if sample.get('reconstruction_bounds_evaluated_rows') != 0:
        raise ValueError('v2 sample may not claim evaluated reconstruction bounds')
    if sample.get('full_reconstruction_enclosure_run'):
        raise ValueError('v2 must not claim a complete reconstruction-enclosure campaign')
    if record['propagation']['propagation_coefficients_evaluated_rows'] != INVENTORY_POSITIVE_ROWS:
        raise ValueError('v2 must bind Lipschitz coefficients on the complete universe')
    if record['propagation']['reconstruction_enclosures_evaluated_rows'] != 0:
        raise ValueError('Lipschitz coefficients are not reconstruction enclosures')
    upstream = record['upstream_physical_slot']
    if (upstream.get('physical_upstream_N_beta') is not None
            or not str(upstream['nine_component_feed']['combined_status']).startswith('OPEN')):
        raise ValueError('v2 upstream physical source slot must remain OPEN')
    if not upstream.get('covered_regions_not_fed_from_column_error'):
        raise ValueError('column error must not be copied into the covered-regions action slot')
    if not record['propagation'].get('linear_indicator_is_not_a_certificate'):
        raise ValueError('uniform Lipschitz threshold is not a certificate')
    forecast = record['resource_forecast']
    if forecast.get('generators_launched', 1) != 0:
        raise ValueError('v2 must not launch source-construction generators')
    if not forecast.get('partial_or_resource_stop_cannot_close_gate_2'):
        raise ValueError('a resource stop cannot close Gate 2')
    return True
