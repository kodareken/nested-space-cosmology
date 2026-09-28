#!/usr/bin/env python3
"""Current-history profile-jet L1 comparison and optional one-cell residual.

No Dirac or source evolution. Parent supplies the saved DOP853 cell. This
adapter replaces the v2 interval-Clenshaw L1 with midpoint Chebyshev jets
plus endpoint remainders. It does not change production evaluation and it
is not a whole-field or physical-gate certificate.

Declared proof parameters (printed before any work):
  comparison CPU cap 60s; one-cell CPU cap 180s; bits 120; derivative
  order 8; transition panels 16; interior panels 4; Fourier M=1024, K=64;
  actual grid period, not nominal 0.4.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import argparse
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import signal
import sys
import time

import numpy as np
from flint import arb, ctx, fmpq

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))

from recursive_horizons.nsc_ks_ball_operator import AnalyticRadiusFamily, time_operator_enclosure
from recursive_horizons.nsc_ks_ball_trajectory import exact_upper, restored_upper
from recursive_horizons.nsc_ks_chebyshev_jet_bound import (
    enclose_profile_fourier_jet, monomial_derivative_l1,
)
from recursive_horizons.nsc_ks_current_history_bounds import radius_bounds, rational_record
from recursive_horizons.nsc_ks_difference_residual_polynomial import (
    DifferenceResidualPolynomial, ball_split, difference_operator_remainder_bounds,
    difference_polynomial_bounds, monomial_keys,
)
from recursive_horizons.nsc_ks_fourier_residual_bound import polynomial_residual_bounds, profile_sup_bounds
from recursive_horizons.nsc_ks_profile_fourier_bound import alias_and_tail_bounds
from recursive_horizons.nsc_ks_profile_identity import profile_identity
from recursive_horizons.nsc_ks_radius_enclosure import axial_profile_bounds, reciprocal_radius_tail
from recursive_horizons.nsc_ks_trajectory import TrajectorySegment
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily

CAPTURE = ROOT / 'results/development/nsc-ks-current-trajectory-pilot-v2.json'
V2_FIELD = ROOT / 'results/development/nsc-ks-current-field-pilot-v2.json'
OUTPUT = ROOT / 'results/development/nsc-ks-current-field-pilot-v3.json'
PROFILE_OUTPUT = ROOT / 'results/development/nsc-ks-current-field-pilot-v3-profile-bounds.json'
CELL_PROOF = ROOT / 'results/development/artifacts/nsc-ks-current-field-pilot-v3-cell.json'
OWNED = (
    'scripts/validate_nsc_ks_current_field_pilot_v3.py',
    'src/recursive_horizons/nsc_ks_chebyshev_jet_bound.py',
    'docs/nsc-ks-current-field-pilot-v3.md',
    'tests/test_nsc_ks_chebyshev_jet_bound.py',
)
INPUTS = (
    'results/development/nsc-ks-current-trajectory-pilot-v2.json',
    'results/development/nsc-ks-current-field-pilot-v2.json',
    'results/development/nsc-ks-gate-history-lm-broyden.json',
    'src/recursive_horizons/nsc_ks_current_history_bounds.py',
    'src/recursive_horizons/nsc_ks_difference_residual_polynomial.py',
    'src/recursive_horizons/nsc_ks_difference_error.py',
    'src/recursive_horizons/nsc_ks_residual_polynomial.py',
    'src/recursive_horizons/nsc_ks_fourier_residual_bound.py',
    'src/recursive_horizons/nsc_ks_profile_fourier_bound.py',
    'src/recursive_horizons/nsc_ks_ball_geometry.py',
    'src/recursive_horizons/nsc_ks_ball_operator.py',
    'src/recursive_horizons/nsc_ks_ball_trajectory.py',
    'src/recursive_horizons/nsc_ks_radius_enclosure.py',
    'src/recursive_horizons/nsc_ks_value_evaluator.py',
)
EXPECTED_IDENTITY = '0b0e4cedfdb695d342c7f7347740631e916bdaee3a977df56808c725b25dacb0'
EFFECTIVE_PERIOD = Fraction(1759218604441, 4398046511104)
BITS, TIME_DEGREE, RECIPROCAL_ORDER = 120, 8, 4
PILOT_CPU_CAP, CELL_CPU_CAP = 60.0, 180.0
FOURIER = {
    'derivative_order': 8,
    'transition_panels': 16,
    'interior_panels': 4,
    'quadrature_points': 1024,
    'retained_index': 64,
    'max_derivative': 2,
    'bits': BITS,
}
PROOF_PARAMETERS = {
    'precision_bits': BITS,
    'time_Taylor_degree': TIME_DEGREE,
    'reciprocal_order': RECIPROCAL_ORDER,
    'fourier': FOURIER,
    'pilot_CPU_cap': PILOT_CPU_CAP,
    'cell_CPU_cap': CELL_CPU_CAP,
    'comparison_uses_v2_saved_old_B8': True,
    'one_direction_metric': True,
    'actual_grid_period': True,
    'nominal_period_length_not_used': 0.4,
    'effective_period_rational': [EFFECTIVE_PERIOD.numerator, EFFECTIVE_PERIOD.denominator],
    'free_contains_zero_used_as_proof': False,
    'production_evaluation_unchanged': True,
    'hash_bound_module_unchanged': True,
}
TAIL_OBSTRUCTION = 'cheb32_mixed_profile_fourier_omitted_band_from_midpoint_jet_l1'
NEXT_GAP = (
    'transition_cutoff_jet_L1_still_dominates_omitted_Fourier_band; '
    'one_cell_unknown_global_field_and_source_error'
)
IMPROVEMENT_RATIO = 1e-6
CONDITIONAL_ERROR_MAPPING = {
    'production_split': "A'=G_ref A, D'=L_g D+M A, A z-independent",
    'error_equations_match': True,
    'R_D_formed_before_norms': True,
    'F_z_split_matches_envelope': (
        'weighted_reference_F_z uses E|e_A| because envelope A_z=0; '
        'weighted_difference_F_z is |e_{D,z}|+E|e_D| for D_z-i E D, not D_z-i E(A+D)'
    ),
    'matter_excludes_reference_density': 'difference_matter_error bounds D C A*+A C D*+D C D* only',
    'issues': [
        'whole-interval integrals are conservative Gronwall majorants, not exact integrating factors',
        'M e_A uses a uniform |e_A| bound times integral ||M||; valid only if those are cone integrals',
        'the module does not compute M_integral, residuals or source error; sampled DOP853 defects are not inputs',
        'Hermitian density uses a factor two; anti-Hermitian current does not, matching finite_matter_error',
        'one saved cell is not a backward-cone field certificate even after propagate_difference_error',
        'free_contains_zero only tests ball-contains-zero and is not used as an identity proof',
    ],
    'used_for_whole_field_certificate': False,
    'free_contains_zero_used_as_proof': False,
}


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT / path).read_bytes()).hexdigest()


def signatures():
    return {path: digest(path) for path in (*OWNED, *INPUTS)}


def write_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(data, indent=2, sort_keys=True, allow_nan=False) + '\n')
    os.replace(temporary, path)


def _scalar(value):
    return float(np.asarray(value).reshape(()))


def _arb_fraction(value):
    return arb(fmpq(value.numerator, value.denominator))


def pack_fraction(value):
    return exact_upper(_arb_fraction(value).upper())


def pack_tail(tail):
    return {
        'radius_ratio_upper': pack_fraction(tail['radius_ratio_upper']),
        'radius_lower': pack_fraction(tail['radius_lower']),
        'potential_derivative_tail_bounds': [pack_fraction(v) for v in tail['potential_derivative_tail_bounds']],
        'order': tail['order'],
        'scope': tail['scope'],
    }


def declare_parameters(kind):
    cap = PILOT_CPU_CAP if kind == 'pilot' else CELL_CPU_CAP
    print(json.dumps({
        'declared_before_run': True,
        'mode': kind,
        'CPU_cap': cap,
        'proof_parameters': PROOF_PARAMETERS,
        'no_Dirac_or_source_evolution': True,
        'physical_local_gate': 'OPEN',
    }, indent=2, sort_keys=True), flush=True)


def load_capture():
    if not CAPTURE.exists():
        raise FileNotFoundError('parent current-trajectory capture is required')
    record = json.loads(CAPTURE.read_text())
    if record.get('schema') != 'NSC-KS-CURRENT-TRAJECTORY-PILOT-v2':
        raise ValueError('unexpected current-trajectory capture schema')
    for name, expected in {**record['source_hashes'], **record['input_hashes']}.items():
        if digest(name) != expected:
            raise ValueError('pilot capture dependency changed: ' + name)
    payload = ROOT / record['payload']['path']
    if digest(payload) != record['payload']['sha256']:
        raise ValueError('pilot capture payload changed')
    if (record['cell_index'] != 122 or record['family'] != [14, 1]
            or record['grid_nodes'] != 1024 or record['source_columns'] != 12
            or len(record['source_rows']) != 4
            or record['profile_identity'] != EXPECTED_IDENTITY):
        raise ValueError('capture is not the declared current Cheb32 family14_1 cell122 pilot')
    with np.load(payload, allow_pickle=False) as arrays:
        data = {key: arrays[key] for key in arrays.files}
    return record, data


def load_family(capture):
    history = json.loads((ROOT / capture['history_path']).read_text())
    family = LocalIncomingFamily(np.array(history['history']['coefficients']))
    identity = profile_identity(family, include_normal_window=True)
    if identity != capture['profile_identity'] or identity != EXPECTED_IDENTITY:
        raise ValueError('current history identity differs from the captured cell')
    return family


def geometry_bounds(family, capture):
    bounds = radius_bounds(family, capture['rho_up'])
    rejected = (axial_profile_bounds(family.functions[0])[0],
                axial_profile_bounds(family.functions[1])[0])
    tail = reciprocal_radius_tail(
        bounds['w']['profile_bounds'], bounds['U']['profile_bounds'],
        normal_support=bounds['normal_support'],
        reference_radius_lower=bounds['reference_radius_lower'],
        axial_lower=bounds['axial_lower'],
        absolute_angular=abs(capture['angular']), order=RECIPROCAL_ORDER)
    return bounds, rejected, tail


def segment_from_capture(capture, arrays):
    grid = np.asarray(arrays['computational_z'], float)
    if len(grid) != capture['grid_nodes']:
        raise ValueError('captured spatial grid changed')
    spacing = Fraction(grid[1]) - Fraction(grid[0])
    period_q = spacing * len(grid)
    period = float(period_q)
    origin = float(grid[0])
    if abs(period - float(EFFECTIVE_PERIOD)) > 1e-18 and period_q != EFFECTIVE_PERIOD:
        # The extractor uses actual grid spacing; the dyadic 1759218604441/2^42
        # is the expected binary64 period of this 1024-node source Fourier grid.
        if abs(period - float(EFFECTIVE_PERIOD)) > 1e-12:
            raise ValueError('captured spatial period is not the measured source-Fourier period')
    segment = TrajectorySegment(
        _scalar(arrays['rho_start']), _scalar(arrays['rho_end']),
        np.asarray(arrays['start']), np.asarray(arrays['end']),
        np.asarray(arrays['dense_corrections']))
    expected = capture['cell_rho_interval']
    if [segment.rho_end, segment.rho_start] != expected:
        raise ValueError('saved cell covers another rho interval')
    weights = np.asarray(arrays['source_weights'], float)
    energies = np.asarray(arrays['source_energies'], float)
    if len(weights) != capture['source_columns'] or len(energies) != len(weights):
        raise ValueError('captured source columns changed')
    return segment, grid, origin, period, period_q, weights, energies


def load_v2_profiles():
    if not V2_FIELD.exists():
        raise FileNotFoundError('v2 field-pilot record is required for the old B8 comparison')
    saved = json.loads(V2_FIELD.read_text())
    if saved.get('schema') != 'NSC-KS-CURRENT-FIELD-PILOT-v2':
        raise ValueError('unexpected v2 field-pilot schema')
    if saved.get('profile_identity') != EXPECTED_IDENTITY:
        raise ValueError('v2 field-pilot history identity changed')
    return saved


def compare_profile_bounds(family, period, v2):
    keys = monomial_keys(RECIPROCAL_ORDER)
    l1, cells = monomial_derivative_l1(
        family, keys, order=FOURIER['derivative_order'],
        transition_panels=FOURIER['transition_panels'],
        interior_panels=FOURIER['interior_panels'], bits=BITS)
    rows = {}
    max_old_b8 = arb(0)
    max_new_b8 = arb(0)
    max_new_tail = arb(0)
    max_old_tail = arb(0)
    for key in keys:
        name = f'{key[0]}_{key[1]}'
        old = v2['profile_enclosure'][name]
        old_b8 = restored_upper(old['derivative_l1_upper'])
        new_b8 = l1[key]
        alias, tails = alias_and_tail_bounds(
            new_b8, period, FOURIER['derivative_order'], FOURIER['quadrature_points'],
            FOURIER['retained_index'], FOURIER['max_derivative'])
        old_alias = restored_upper(old['alias_error_upper'])
        old_tails = [restored_upper(v) for v in old['tail_upper']]
        rows[name] = {
            'powers': [key[0], key[1]],
            'old_derivative_l1_upper': exact_upper(old_b8),
            'new_derivative_l1_upper': exact_upper(new_b8),
            'old_alias_error_upper': exact_upper(old_alias),
            'new_alias_error_upper': exact_upper(alias),
            'old_tail_upper': [exact_upper(v) for v in old_tails],
            'new_tail_upper': [exact_upper(v) for v in tails],
            'B8_ratio_new_over_old_display': float(new_b8 / old_b8),
        }
        if old_b8 > max_old_b8:
            max_old_b8 = old_b8
        if new_b8 > max_new_b8:
            max_new_b8 = new_b8
        if tails[0] > max_new_tail:
            max_new_tail = tails[0]
        if old_tails[0] > max_old_tail:
            max_old_tail = old_tails[0]
    ratio = float(max_new_b8 / max_old_b8)
    improved = ratio <= IMPROVEMENT_RATIO and max_new_b8 < max_old_b8
    return {
        'settings': FOURIER,
        'derivative_cells': cells,
        'monomials': rows,
        'max_old_B8_display': float(max_old_b8),
        'max_new_B8_display': float(max_new_b8),
        'max_old_tail_order0_display': float(max_old_tail),
        'max_new_tail_order0_display': float(max_new_tail),
        'max_B8_ratio_new_over_old_display': ratio,
        'meaningful_improvement': improved,
        'improvement_ratio_threshold': IMPROVEMENT_RATIO,
        'dominating_key': '0_4',
        'method': (
            'midpoint Chebyshev jet plus endpoint/Lagrange remainder, '
            'owned plateau_series cutoff, actual 32-coefficient one-direction profiles'
        ),
    }


def enclose_profiles(family, origin, period):
    keys = monomial_keys(RECIPROCAL_ORDER)
    enclosed = enclose_profile_fourier_jet(family, keys, origin, period, **FOURIER)
    profiles = {key: {'coefficients': list(item.coefficients),
                      'tail': list(item.uniform_tail_bounds)}
                for key, item in enclosed.items()}
    inventory = {
        f'{p}_{q}': {
            'powers': [p, q],
            'derivative_l1_upper': exact_upper(item.derivative_l1_upper),
            'alias_error_upper': exact_upper(item.alias_error_upper),
            'tail_upper': [exact_upper(v) for v in item.uniform_tail_bounds],
            'derivative_cells': item.derivative_cells,
        } for (p, q), item in enclosed.items()}
    return profiles, inventory


def evaluate_cell(capture, arrays, family, bounds, tail):
    segment, grid, origin, period, period_q, weights, energies = segment_from_capture(capture, arrays)
    model = AnalyticRadiusFamily(family)
    keys = monomial_keys(RECIPROCAL_ORDER)
    with ctx.workprec(BITS):
        profiles, inventory = enclose_profiles(family, origin, period)
        max_tail = max(float(restored_upper(row['tail_upper'][0])) for row in inventory.values())
        obstruction = TAIL_OBSTRUCTION if max_tail > 1e3 else None
        field_x, field_a, field_d = ball_split(
            segment, weights, len(grid), period, bits=BITS)
        polynomials, errors = time_operator_enclosure(
            model, field_x.rho_start, field_x.rho_end, capture['angular'],
            degree=TIME_DEGREE, reciprocal_order=RECIPROCAL_ORDER, bits=BITS)
        present = {key: value for key, value in polynomials.items() if isinstance(key, tuple)}
        if set(present) != set(keys):
            raise ValueError('time-series monomials are not the fourteen reciprocal keys')
        residual = DifferenceResidualPolynomial(
            field_x, field_a, polynomials['inv_a2'], polynomials['inv_a'],
            polynomials['inv_ar'], present, capture['mass'], capture['angular'], energies)
        finite_profiles = {key: {**row, 'tail': (arb(0), arb(0))} for key, row in profiles.items()}
        polynomial = difference_polynomial_bounds(residual, profiles)
        finite = polynomial_residual_bounds(residual, finite_profiles)
        radius_tail = [_arb_fraction(v) for v in tail['potential_derivative_tail_bounds'][:2]]
        remainder = difference_operator_remainder_bounds(
            field_d, field_a, errors, profile_sup_bounds(profiles, field_x.length),
            radius_tail, capture['mass'], capture['angular'], energies)
        total = tuple((left + right).upper() for left, right in zip(polynomial, remainder))
        finite_vals = [float(v) for v in finite]
        rem_vals = [float(v) for v in remainder]
        poly_vals = [float(v) for v in polynomial]
        dominating = 'omitted_profile_Fourier_band'
        if poly_vals[0] <= 2 * finite_vals[0] and rem_vals[0] >= finite_vals[0]:
            dominating = 'time_or_radius_remainder'
        elif poly_vals[0] <= 2 * finite_vals[0]:
            dominating = 'polynomial_finite_Fourier_band'
    return {
        'profiles': inventory,
        'max_profile_tail_order0_display': max_tail,
        'named_obstruction': obstruction,
        'dominating_continuous_contribution': dominating,
        'effective_period_rational': [period_q.numerator, period_q.denominator],
        'bounds': {
            'polynomial_and_profile_tail': [exact_upper(v) for v in polynomial],
            'polynomial_finite_band': [exact_upper(v) for v in finite],
            'time_coefficient_and_radius_remainder': [exact_upper(v) for v in remainder],
            'total_continuous_normalized_residual': [exact_upper(v) for v in total],
        },
        'coverage': {
            'segment_index': capture['cell_index'],
            'accepted_steps': capture['accepted_steps'],
            'rho_interval': [segment.rho_end, segment.rho_start],
            'spatial_period_origin': origin,
            'spatial_period_length': period,
            'declared_period_length': capture['period_length'],
            'effective_period_rational': [EFFECTIVE_PERIOD.numerator, EFFECTIVE_PERIOD.denominator],
            'whole_spatial_period': True,
            'whole_time_cell': True,
            'all_history_segments': False,
            'source_columns': int(len(weights)),
            'reciprocal_keys': [list(key) for key in keys],
        },
        'background_residual_ball_contains_zero_is_not_an_identity_certificate': True,
    }


def _timed(cap, callback):
    def stop(*_):
        raise TimeoutError(f'current-field v3 exceeded {cap} CPU seconds')

    previous = signal.signal(signal.SIGPROF, stop)
    signal.setitimer(signal.ITIMER_PROF, cap)
    try:
        return callback()
    finally:
        signal.setitimer(signal.ITIMER_PROF, 0)
        signal.signal(signal.SIGPROF, previous)


def profile_record(comparison, runtime, sig):
    return {
        'schema': 'NSC-KS-CURRENT-FIELD-PILOT-PROFILE-BOUNDS-v3',
        'accountable_author': 'Douglas Ek',
        'status': (
            'OPEN: midpoint-jet L1 comparison; omitted Fourier band still dominates; '
            'physical local gate OPEN'
        ),
        'physical_local_gate': 'OPEN',
        'physical_EXISTENCE_certificate': False,
        'physical_NON_EXISTENCE_certificate': False,
        'profile_identity': EXPECTED_IDENTITY,
        'proof_parameters': PROOF_PARAMETERS,
        'comparison': comparison,
        'named_obstruction': TAIL_OBSTRUCTION if comparison['max_new_tail_order0_display'] > 1e3 else None,
        'named_next_gap': NEXT_GAP,
        'runtime': runtime,
        'source_hashes': {path: sig[path] for path in OWNED},
        'input_hashes': {path: sig[path] for path in INPUTS},
        'reproducer': 'PYTHONPATH=src .venv/validation/bin/python scripts/validate_nsc_ks_current_field_pilot_v3.py --check',
    }


def field_record(capture, family, bounds, rejected, tail, comparison, evaluation, runtime, sig):
    w_old, u_old = rejected
    obstruction = None if evaluation is None else evaluation['named_obstruction']
    status = (
        'OPEN: one-cell continuous residual bound with midpoint-jet L1; '
        'omitted Fourier band of Cheb32 U-powers still dominates; '
        'unknown global field/source error; physical local gate OPEN'
        if obstruction else
        'OPEN: one-cell continuous residual bound; whole-field and physical gates remain open'
        if evaluation is not None else
        'OPEN: midpoint-jet L1 comparison only; one-cell not run'
    )
    with ctx.workprec(BITS):
        packed_tail = pack_tail(tail)
    result = {
        'schema': 'NSC-KS-CURRENT-FIELD-PILOT-v3',
        'accountable_author': 'Douglas Ek',
        'status': status,
        'physical_local_gate': 'OPEN',
        'physical_EXISTENCE_certificate': False,
        'physical_NON_EXISTENCE_certificate': False,
        'profile_identity': capture['profile_identity'],
        'history_path': capture['history_path'],
        'family': capture['family'],
        'cell_index': capture['cell_index'],
        'source_rows': capture['source_rows'],
        'source_preparation_digest': capture['source_preparation_digest'],
        'state_law': 'C_Sigma[g]=U_g C_up U_gdagger; fixed upstream source; D=X-A',
        'proof_parameters': PROOF_PARAMETERS,
        'geometry': {
            'radius_bounds': rational_record(bounds),
            'reciprocal_radius_tail': packed_tail,
            'power_coefficient_absolute_sum_rejected': {
                'w0_binary64_upper': float(w_old),
                'U0_binary64_upper': float(u_old),
                'reason': 'Chebyshev endpoint bounds retain the represented basis',
            },
        },
        'method': {
            'precision_bits': BITS,
            'time_Taylor_degree': TIME_DEGREE,
            'reciprocal_order': RECIPROCAL_ORDER,
            'fourier': FOURIER,
            'profile_jet': (
                'midpoint Clenshaw of the represented Chebyshev polynomial plus '
                'mean-value/Lagrange remainder from global endpoint derivatives; '
                'owned plateau_series cutoff jets; one_direction_metric 32-coeff profiles'
            ),
            'old_fine_proof_transferred': False,
            'v2_interval_Clenshaw_L1_replaced': True,
            'free_contains_zero_used_as_proof': False,
        },
        'comparison': comparison,
        'named_next_gap': NEXT_GAP,
        'conditional_error_mapping': CONDITIONAL_ERROR_MAPPING,
        'scope': {
            'saved_cells': 0 if evaluation is None else 1,
            'full_source_covered': False,
            'continuous_full_history_field_error': None,
            'source_error_bound': None,
            'physical_local_gate': 'OPEN',
            'field_runs': 0,
            'source_runs': 0,
            'metric_timestep': False,
            'sample_max_used_as_bound': False,
            'propagate_difference_error_applied': False,
            'one_cell_unknown_global_field_and_source_error': True,
        },
        'runtime': runtime,
        'capture': {
            'path': str(CAPTURE.relative_to(ROOT)),
            'sha256': digest(CAPTURE),
            'payload': capture['payload'],
            'capture_CPU_seconds': capture['runtime']['CPU_seconds'],
        },
        'source_hashes': {path: sig[path] for path in OWNED},
        'input_hashes': {path: sig[path] for path in INPUTS},
        'reproducer': 'PYTHONPATH=src .venv/validation/bin/python scripts/validate_nsc_ks_current_field_pilot_v3.py --check',
    }
    if evaluation is not None:
        result['bounds'] = evaluation['bounds']
        result['coverage'] = evaluation['coverage']
        result['profile_enclosure'] = evaluation['profiles']
        result['named_obstruction'] = evaluation['named_obstruction']
        result['max_profile_tail_order0_display'] = evaluation['max_profile_tail_order0_display']
        result['dominating_continuous_contribution'] = evaluation['dominating_continuous_contribution']
        result['scope']['saved_cells'] = 1
    else:
        result['named_obstruction'] = TAIL_OBSTRUCTION
        result['bounds'] = None
        result['coverage'] = None
        result['profile_enclosure'] = None
    return result


def run_pilot():
    if PROFILE_OUTPUT.exists():
        raise FileExistsError('v3 profile-bound comparison exists; use --check')
    declare_parameters('pilot')
    sig = signatures()
    capture, arrays = load_capture()
    family = load_family(capture)
    v2 = load_v2_profiles()
    _segment, _grid, _origin, period, _period_q, _weights, _energies = segment_from_capture(capture, arrays)
    cpu, wall = time.process_time(), time.perf_counter()
    with ctx.workprec(BITS):
        comparison = _timed(PILOT_CPU_CAP, lambda: compare_profile_bounds(family, period, v2))
    runtime = {
        'CPU_seconds': time.process_time() - cpu,
        'wall_seconds': time.perf_counter() - wall,
        'CPU_cap': PILOT_CPU_CAP,
    }
    if sig != signatures():
        raise ValueError('v3 profile comparison source/input changed during evaluation')
    result = profile_record(comparison, runtime, sig)
    write_json(PROFILE_OUTPUT, result)
    runtime['storage_bytes'] = {'profile_record': PROFILE_OUTPUT.stat().st_size}
    result['runtime'] = runtime
    write_json(PROFILE_OUTPUT, result)
    return result


def run_cell():
    if OUTPUT.exists() or CELL_PROOF.exists():
        raise FileExistsError('v3 field-pilot record exists; use --check')
    declare_parameters('cell')
    sig = signatures()
    capture, arrays = load_capture()
    family = load_family(capture)
    bounds, rejected, tail = geometry_bounds(family, capture)
    v2 = load_v2_profiles()
    cpu, wall = time.process_time(), time.perf_counter()

    def work():
        _segment, _grid, origin, period, _period_q, _weights, _energies = segment_from_capture(capture, arrays)
        comparison = compare_profile_bounds(family, period, v2)
        evaluation = None
        if comparison['meaningful_improvement']:
            evaluation = evaluate_cell(capture, arrays, family, bounds, tail)
        return comparison, evaluation

    with ctx.workprec(BITS):
        comparison, evaluation = _timed(CELL_CPU_CAP, work)
    runtime = {
        'CPU_seconds': time.process_time() - cpu,
        'wall_seconds': time.perf_counter() - wall,
        'CPU_cap': CELL_CPU_CAP,
        'storage_bytes': {
            'capture_payload': capture['payload']['bytes'],
            'capture_record': CAPTURE.stat().st_size,
        },
        'one_cell_evaluated': evaluation is not None,
    }
    if sig != signatures():
        raise ValueError('v3 field-pilot source/input changed during evaluation')
    result = field_record(capture, family, bounds, rejected, tail, comparison, evaluation, runtime, sig)
    if evaluation is not None:
        write_json(CELL_PROOF, {
            'schema': 'NSC-KS-CURRENT-FIELD-PILOT-CELL-v3',
            'cell_index': capture['cell_index'],
            'bounds': evaluation['bounds'],
            'named_obstruction': evaluation['named_obstruction'],
            'dominating_continuous_contribution': evaluation['dominating_continuous_contribution'],
            'CPU_seconds': runtime['CPU_seconds'],
            'signature_digest': sha256(json.dumps(sig, sort_keys=True).encode()).hexdigest(),
            'physical_local_gate': 'OPEN',
            'free_contains_zero_used_as_proof': False,
        })
        runtime['storage_bytes']['cell_proof'] = CELL_PROOF.stat().st_size
    result['runtime'] = runtime
    write_json(OUTPUT, result)
    runtime['storage_bytes']['field_record'] = OUTPUT.stat().st_size
    result['runtime'] = runtime
    write_json(OUTPUT, result)
    return result


def check():
    target = OUTPUT if OUTPUT.exists() else PROFILE_OUTPUT
    if not target.exists():
        raise FileNotFoundError('v3 field-pilot record is missing; run --pilot or --run first')
    saved = json.loads(target.read_text())
    sig = signatures()
    for path, expected in {**saved['source_hashes'], **saved['input_hashes']}.items():
        if path.endswith('.md') and target == PROFILE_OUTPUT and path not in sig:
            continue
        if sig[path] != expected:
            raise ValueError('v3 owner/input changed: ' + path)
    if saved['physical_local_gate'] != 'OPEN' or saved['physical_EXISTENCE_certificate']:
        raise ValueError('the v3 pilot must not claim a physical gate')
    if saved['proof_parameters']['free_contains_zero_used_as_proof']:
        raise ValueError('ball-contains-zero must not be used as an identity proof')
    comparison = saved['comparison']
    if comparison['max_new_B8_display'] >= comparison['max_old_B8_display']:
        raise ValueError('recorded new B8 is not an improvement of the v2 L1')
    capture, _arrays = load_capture()
    if OUTPUT.exists():
        if digest(CAPTURE) != saved['capture']['sha256']:
            raise ValueError('pilot capture record changed')
        family = load_family(capture)
        bounds, _rejected, tail = geometry_bounds(family, capture)
        if rational_record(bounds)['radius_lower'] != saved['geometry']['radius_bounds']['radius_lower']:
            raise ValueError('current Chebyshev radius lower bound changed')
        with ctx.workprec(BITS):
            if pack_fraction(tail['radius_ratio_upper']) != saved['geometry']['reciprocal_radius_tail']['radius_ratio_upper']:
                raise ValueError('current reciprocal-radius remainder changed')
            if saved.get('bounds'):
                poly = saved['bounds']['polynomial_and_profile_tail']
                rem = saved['bounds']['time_coefficient_and_radius_remainder']
                total = saved['bounds']['total_continuous_normalized_residual']
                for a, b, t in zip(poly, rem, total):
                    if not restored_upper(t) >= (restored_upper(a) + restored_upper(b)).upper():
                        raise ValueError('one-cell residual contributions are not enclosed')
        if saved['coverage'] and saved['coverage']['all_history_segments']:
            raise ValueError('a one-cell pilot cannot certify the entire field evolution')
        if saved['conditional_error_mapping']['used_for_whole_field_certificate']:
            raise ValueError('conditional difference-error majorants were promoted to a field certificate')
        if saved['scope']['continuous_full_history_field_error'] is not None:
            raise ValueError('a one-cell pilot cannot certify the entire field evolution')
    return saved


def display(result):
    payload = {
        'status': result['status'],
        'physical_local_gate': result['physical_local_gate'],
        'named_obstruction': result.get('named_obstruction'),
        'named_next_gap': result.get('named_next_gap'),
        'comparison': {
            'max_old_B8_display': result['comparison']['max_old_B8_display'],
            'max_new_B8_display': result['comparison']['max_new_B8_display'],
            'max_old_tail_order0_display': result['comparison']['max_old_tail_order0_display'],
            'max_new_tail_order0_display': result['comparison']['max_new_tail_order0_display'],
            'max_B8_ratio_new_over_old_display': result['comparison']['max_B8_ratio_new_over_old_display'],
            'meaningful_improvement': result['comparison']['meaningful_improvement'],
            'U4': result['comparison']['monomials']['0_4'],
        },
        'runtime': result['runtime'],
    }
    if result.get('bounds'):
        with ctx.workprec(BITS):
            payload['total_continuous_residual_display'] = [
                float(restored_upper(v)) for v in result['bounds']['total_continuous_normalized_residual']]
            payload['polynomial_finite_band_display'] = [
                float(restored_upper(v)) for v in result['bounds']['polynomial_finite_band']]
            payload['time_radius_remainder_display'] = [
                float(restored_upper(v))
                for v in result['bounds']['time_coefficient_and_radius_remainder']]
        payload['dominating_continuous_contribution'] = result.get('dominating_continuous_contribution')
        payload['cell_index'] = result.get('cell_index')
    return payload


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument('--run', action='store_true', help='comparison then one-cell if improved, 180 CPU s')
    modes.add_argument('--pilot', action='store_true', help='14-monomial L1 comparison only, 60 CPU s')
    modes.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.check:
        result = check()
    elif args.pilot:
        result = run_pilot()
    else:
        result = run_cell()
    print(json.dumps(display(result), indent=2, sort_keys=True))
