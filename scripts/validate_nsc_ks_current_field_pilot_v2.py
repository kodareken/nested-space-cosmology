#!/usr/bin/env python3
"""One-cell continuous residual adapter for the current-history DOP853 capture.

No Dirac or source evolution. Parent supplies the saved cell. This adapter
bounds that cell's difference residual, or names a specific analytic
obstruction. It is not a whole-field or physical-gate certificate.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import argparse
from hashlib import sha256
import json
import os
from pathlib import Path
import signal
import sys
import time

import numpy as np
from flint import arb, ctx, fmpq

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
sys.path.insert(0, str(ROOT / 'scripts'))

from recursive_horizons.nsc_ks_ball_operator import AnalyticRadiusFamily, time_operator_enclosure
from recursive_horizons.nsc_ks_ball_trajectory import exact_upper, restored_upper
from recursive_horizons.nsc_ks_current_history_bounds import radius_bounds, rational_record
from recursive_horizons.nsc_ks_difference_residual_polynomial import (
    DifferenceResidualPolynomial, ball_split, difference_operator_remainder_bounds,
    difference_polynomial_bounds, monomial_keys,
)
from recursive_horizons.nsc_ks_fourier_residual_bound import polynomial_residual_bounds, profile_sup_bounds
from recursive_horizons.nsc_ks_profile_fourier_bound import enclose_profile_fourier
from recursive_horizons.nsc_ks_profile_identity import profile_identity
from recursive_horizons.nsc_ks_radius_enclosure import axial_profile_bounds, reciprocal_radius_tail
from recursive_horizons.nsc_ks_trajectory import TrajectorySegment
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily

CAPTURE = ROOT / 'results/development/nsc-ks-current-trajectory-pilot-v2.json'
OUTPUT = ROOT / 'results/development/nsc-ks-current-field-pilot-v2.json'
CELL_PROOF = ROOT / 'results/development/artifacts/nsc-ks-current-field-pilot-v2-cell.json'
OWNED = (
    'scripts/validate_nsc_ks_current_field_pilot_v2.py',
    'src/recursive_horizons/nsc_ks_difference_residual_polynomial.py',
    'docs/nsc-ks-current-field-pilot-v2.md',
)
INPUTS = (
    'results/development/nsc-ks-current-trajectory-pilot-v2.json',
    'src/recursive_horizons/nsc_ks_current_history_bounds.py',
    'src/recursive_horizons/nsc_ks_difference_error.py',
    'src/recursive_horizons/nsc_ks_residual_polynomial.py',
    'src/recursive_horizons/nsc_ks_fourier_residual_bound.py',
    'src/recursive_horizons/nsc_ks_profile_fourier_bound.py',
    'src/recursive_horizons/nsc_ks_ball_operator.py',
    'src/recursive_horizons/nsc_ks_ball_trajectory.py',
    'src/recursive_horizons/nsc_ks_radius_enclosure.py',
    'results/development/nsc-ks-gate-history-lm-broyden.json',
)
EXPECTED_IDENTITY = '0b0e4cedfdb695d342c7f7347740631e916bdaee3a977df56808c725b25dacb0'
BITS, TIME_DEGREE, RECIPROCAL_ORDER, CPU_CAP = 120, 8, 4, 180.0
FOURIER = {
    'derivative_order': 8,
    'transition_panels': 16,
    'quadrature_points': 1024,
    'retained_index': 64,
    'max_derivative': 2,
    'bits': BITS,
}
TAIL_OBSTRUCTION = 'cheb32_mixed_profile_fourier_omitted_band_from_global_derivative_l1'
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
    ],
    'used_for_whole_field_certificate': False,
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
    """Dyadic upper bound of an exact rational, without decimal digit explosion."""
    return exact_upper(_arb_fraction(value).upper())


def pack_tail(tail):
    return {
        'radius_ratio_upper': pack_fraction(tail['radius_ratio_upper']),
        'radius_lower': pack_fraction(tail['radius_lower']),
        'potential_derivative_tail_bounds': [pack_fraction(v) for v in tail['potential_derivative_tail_bounds']],
        'order': tail['order'],
        'scope': tail['scope'],
    }


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
    period = float((grid[1] - grid[0]) * len(grid))
    origin = float(grid[0])
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
    return segment, grid, origin, period, weights, energies


def enclose_profiles(model, origin, period):
    keys = monomial_keys(RECIPROCAL_ORDER)
    enclosed = enclose_profile_fourier(model, keys, origin, period, **FOURIER)
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
    segment, grid, origin, period, weights, energies = segment_from_capture(capture, arrays)
    model = AnalyticRadiusFamily(family)
    keys = monomial_keys(RECIPROCAL_ORDER)
    with ctx.workprec(BITS):
        profiles, inventory = enclose_profiles(model, origin, period)
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
    return {
        'profiles': inventory,
        'max_profile_tail_order0_display': max_tail,
        'named_obstruction': obstruction,
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
            'whole_spatial_period': True,
            'whole_time_cell': True,
            'all_history_segments': False,
            'source_columns': int(len(weights)),
            'reciprocal_keys': [list(key) for key in keys],
        },
        'zero_difference_background_cancelled': residual.free_contains_zero(),
    }


def record(capture, family, bounds, rejected, tail, evaluation, runtime, sig):
    w_old, u_old = rejected
    status = ('OPEN: one-cell continuous residual bound; Fourier omitted-band of Cheb32 U-powers '
              'blocks a deciding field certificate'
              if evaluation['named_obstruction'] else
              'OPEN: one-cell continuous residual bound; whole-field and physical gates remain open')
    with ctx.workprec(BITS):
        packed_tail = pack_tail(tail)
    return {
        'schema': 'NSC-KS-CURRENT-FIELD-PILOT-v2',
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
            'old_fine_proof_transferred': False,
            'old_fine_proof': 'ALPHA=0.001, U=0, 246 cells, four (n,0) profiles',
        },
        'bounds': evaluation['bounds'],
        'coverage': evaluation['coverage'],
        'profile_enclosure': evaluation['profiles'],
        'named_obstruction': evaluation['named_obstruction'],
        'max_profile_tail_order0_display': evaluation['max_profile_tail_order0_display'],
        'conditional_error_mapping': CONDITIONAL_ERROR_MAPPING,
        'scope': {
            'saved_cells': 1,
            'full_source_covered': False,
            'continuous_full_history_field_error': None,
            'source_error_bound': None,
            'physical_local_gate': 'OPEN',
            'field_runs': 0,
            'source_runs': 0,
            'metric_timestep': False,
            'sample_max_used_as_bound': False,
            'propagate_difference_error_applied': False,
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
        'reproducer': 'PYTHONPATH=src .venv/validation/bin/python scripts/validate_nsc_ks_current_field_pilot_v2.py --check',
    }


def run():
    if OUTPUT.exists() or CELL_PROOF.exists():
        raise FileExistsError('current-field pilot record exists; use --check')
    sig = signatures()
    capture, arrays = load_capture()
    family = load_family(capture)
    bounds, rejected, tail = geometry_bounds(family, capture)
    cpu, wall = time.process_time(), time.perf_counter()

    def stop(*_):
        raise TimeoutError('current-field pilot exceeded 180 CPU seconds')

    previous = signal.signal(signal.SIGPROF, stop)
    signal.setitimer(signal.ITIMER_PROF, CPU_CAP)
    try:
        with ctx.workprec(BITS):
            evaluation = evaluate_cell(capture, arrays, family, bounds, tail)
    finally:
        signal.setitimer(signal.ITIMER_PROF, 0)
        signal.signal(signal.SIGPROF, previous)
    runtime = {
        'CPU_seconds': time.process_time() - cpu,
        'wall_seconds': time.perf_counter() - wall,
        'CPU_cap': CPU_CAP,
        'storage_bytes': {
            'capture_payload': capture['payload']['bytes'],
            'capture_record': CAPTURE.stat().st_size,
        },
    }
    if sig != signatures():
        raise ValueError('current-field pilot source/input changed during evaluation')
    result = record(capture, family, bounds, rejected, tail, evaluation, runtime, sig)
    write_json(CELL_PROOF, {
        'schema': 'NSC-KS-CURRENT-FIELD-PILOT-CELL-v2',
        'cell_index': capture['cell_index'],
        'bounds': evaluation['bounds'],
        'named_obstruction': evaluation['named_obstruction'],
        'CPU_seconds': runtime['CPU_seconds'],
        'signature_digest': sha256(json.dumps(sig, sort_keys=True).encode()).hexdigest(),
        'physical_local_gate': 'OPEN',
    })
    runtime['storage_bytes']['cell_proof'] = CELL_PROOF.stat().st_size
    result['runtime'] = runtime
    write_json(OUTPUT, result)
    runtime['storage_bytes']['field_record'] = OUTPUT.stat().st_size
    write_json(OUTPUT, result)
    return result


def check():
    if not OUTPUT.exists():
        raise FileNotFoundError('current-field pilot record is missing; run the one-cell proof first')
    saved = json.loads(OUTPUT.read_text())
    sig = signatures()
    for path, expected in {**saved['source_hashes'], **saved['input_hashes']}.items():
        if sig[path] != expected:
            raise ValueError('current-field pilot owner/input changed: ' + path)
    capture, _arrays = load_capture()
    if digest(CAPTURE) != saved['capture']['sha256']:
        raise ValueError('pilot capture record changed')
    if capture['payload'] != saved['capture']['payload']:
        raise ValueError('pilot capture payload binding changed')
    family = load_family(capture)
    bounds, _rejected, tail = geometry_bounds(family, capture)
    if rational_record(bounds)['radius_lower'] != saved['geometry']['radius_bounds']['radius_lower']:
        raise ValueError('current Chebyshev radius lower bound changed')
    with ctx.workprec(BITS):
        if pack_fraction(tail['radius_ratio_upper']) != saved['geometry']['reciprocal_radius_tail']['radius_ratio_upper']:
            raise ValueError('current reciprocal-radius remainder changed')
    with ctx.workprec(BITS):
        poly = saved['bounds']['polynomial_and_profile_tail']
        rem = saved['bounds']['time_coefficient_and_radius_remainder']
        total = saved['bounds']['total_continuous_normalized_residual']
        for a, b, t in zip(poly, rem, total):
            if not restored_upper(t) >= (restored_upper(a) + restored_upper(b)).upper():
                raise ValueError('one-cell residual contributions are not enclosed')
    if saved['coverage']['all_history_segments'] or saved['scope']['continuous_full_history_field_error'] is not None:
        raise ValueError('a one-cell pilot cannot certify the entire field evolution')
    if saved['physical_local_gate'] != 'OPEN' or saved['physical_EXISTENCE_certificate']:
        raise ValueError('the pilot must not claim a physical gate')
    if saved['conditional_error_mapping']['used_for_whole_field_certificate']:
        raise ValueError('conditional difference-error majorants were promoted to a field certificate')
    return saved


def display(result):
    with ctx.workprec(BITS):
        total = [float(restored_upper(v))
                 for v in result['bounds']['total_continuous_normalized_residual']]
        finite = [float(restored_upper(v))
                  for v in result['bounds']['polynomial_finite_band']]
        eta = float(restored_upper(result['geometry']['reciprocal_radius_tail']['radius_ratio_upper']))
    return {
        'status': result['status'],
        'physical_local_gate': result['physical_local_gate'],
        'named_obstruction': result['named_obstruction'],
        'total_continuous_residual_display': total,
        'polynomial_finite_band_display': finite,
        'radius_lower_display': float(result['geometry']['radius_bounds']['radius_lower']['binary64_upper']),
        'eta_display': eta,
        'runtime': result['runtime'],
        'cell_index': result['cell_index'],
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument('--run', action='store_true')
    modes.add_argument('--pilot', action='store_true')
    modes.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = check() if args.check else run()
    print(json.dumps(display(result), indent=2, sort_keys=True))
