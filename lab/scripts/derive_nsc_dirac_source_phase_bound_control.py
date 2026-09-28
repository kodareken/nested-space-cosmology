#!/usr/bin/env python3
"""One-point directed f_s,z enclosure; no field or source campaign."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
import signal
import sys
import time

from flint import ctx

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
sys.path.insert(0, str(ROOT / 'scripts'))
import derive_nsc_local_prepared_response as P
from recursive_horizons.nsc_dirac_source_phase_bound import (
    ANGULAR_SQUARE_WEIGHT, BITS, CPU_CAP, MAX_CELLS, MAX_DEPTH, SOURCE_SIGNS,
    aggregate_phase_error_uppers, directed_source_phase_z_enclosure,
    pack_real_ball, pack_upper, unit_angular_target,
)
from recursive_horizons.nsc_ks_local_embedding import continuum_local_embedding
from recursive_horizons.nsc_ks_profile_identity import profile_identity
from recursive_horizons.nsc_ks_ball_geometry import AXIAL_LOWER, RHO_MAX

OUTPUT = ROOT / 'results/development/nsc-dirac-source-phase-bound-control.json'
CONTROL = ROOT / 'results/development/nsc-ks-cutoff-bridge-control.json'
OWNED = (
    'scripts/derive_nsc_dirac_source_phase_bound_control.py',
    'src/recursive_horizons/nsc_dirac_source_phase_bound.py',
    'docs/nsc-dirac-source-phase-bound.md',
    'tests/test_nsc_dirac_source_phase_bound.py',
)
INPUTS = (
    'results/development/nsc-ks-cutoff-bridge-control.json',
    'src/recursive_horizons/nsc_ks_ball_geometry.py',
    'src/recursive_horizons/nsc_ks_ball_operator.py',
    'src/recursive_horizons/nsc_ks_profile_identity.py',
    'src/recursive_horizons/nsc_ks_local_embedding.py',
    'src/recursive_horizons/nsc_local_incoming_family.py',
    'src/recursive_horizons/nsc_compatible_history_geometry.py',
    'scripts/derive_nsc_local_prepared_response.py',
)


def digest(path):
    return sha256((ROOT / path).read_bytes()).hexdigest()


def compute():
    control = json.loads(CONTROL.read_text())
    family = P.family(P.ALPHA)
    identity = profile_identity(family, include_normal_window=True)
    if identity != control['history_identity']:
        raise ValueError('enclosure must use the saved analytic history')
    rho_up = control['rho_up']
    z = control['z'][0]
    weights = control['signed_angular_square_weights']
    if list(weights) != [ANGULAR_SQUARE_WEIGHT, ANGULAR_SQUARE_WEIGHT]:
        raise ValueError('original aggregate angular-square weight per sign required')
    if not (1.0 < rho_up <= float(RHO_MAX)):
        raise ValueError('control rho_up must lie in the owned slab')
    embedding = continuum_local_embedding(
        rho_up=Q(*float(rho_up).as_integer_ratio()), rho_sigma=1,
        axial_lower=AXIAL_LOWER, period_left=Q(-205, 1024),
        period_length=Q(205, 512),
        physical_interval=(Q(-3, 100), Q(3, 100)),
        axial_support=(Q(-3, 50), Q(3, 50)),
        radius_positive=True)
    if not embedding['continuum_local_embedding']:
        raise ValueError('local embedding of this incoming interval failed')
    result = directed_source_phase_z_enclosure(
        family, z, rho_up=rho_up, bits=BITS, cpu_limit=CPU_CAP,
        max_depth=MAX_DEPTH, max_cells=MAX_CELLS)
    with ctx.workprec(BITS):
        target = unit_angular_target()
        balls = result.balls
        if result.obstruction is None and len(balls) == 2:
            newton, beta = aggregate_phase_error_uppers(*balls)
            radii = [pack_real_ball(v) for v in balls]
            newton_u, beta_u = pack_upper(newton), pack_upper(beta)
            pass_radius = all(v.rad() <= target for v in balls)
        else:
            newton_u = beta_u = None
            radii = [pack_real_ball(v) for v in balls] if balls else []
            pass_radius = False
    status = result.status
    record = {
        'schema': 'NSC-DIRAC-SOURCE-PHASE-BOUND-CONTROL-v1',
        'accountable_author': 'Douglas Ek',
        'status': status,
        'obstruction': result.obstruction,
        'history_identity': identity,
        'z': z,
        'rho_up': rho_up,
        'source_signs': list(SOURCE_SIGNS),
        'angular': 1.0,
        'angular_square_weight_per_energy_sign': ANGULAR_SQUARE_WEIGHT,
        'taylor_degree': result.taylor_degree,
        'remainder_order': result.remainder_order,
        'dyadic_scale': result.dyadic_scale,
        'cover_cells': result.cover_cells,
        'cells': list(result.cells),
        'max_depth': list(result.max_depth),
        'splits': list(result.splits),
        'limits': {'CPU_cap': CPU_CAP, 'max_cells': MAX_CELLS, 'max_depth': MAX_DEPTH,
                   'bits': BITS, 'unit_angular_radius_target': '2e-16'},
        'f_z_unit_angular': radii,
        'aggregate_phase_value_error_upper_N_beta': [newton_u, beta_u],
        'unit_angular_radius_pass': bool(pass_radius),
        'coefficient_error_bound': result.coefficient_error_bound,
        'remainder_error_bound': result.remainder_error_bound,
        'physical_error_bound': result.physical_error_bound,
        'source_tail_error_bound': result.source_tail_error_bound,
        'diagnostics': result.diagnostics,
        'local_embedding': {
            'reused': True,
            'continuum_local_embedding': embedding['continuum_local_embedding'],
            'global_matching_required': embedding['global_matching_required'],
        },
        'scope': {
            'pilot_points': 1,
            'all_21_incoming_nodes': False,
            'physical_local_gate': 'OPEN',
            'gauss_nodes_are_bounds': False,
            'new_field_or_source_runs': 0,
            'metric_timestep': False,
            'assembly_correction': False,
            'source_tail_error_bound': None,
            'field_error_bound': None,
            'between_node_error_bound': None,
        },
        'source_hashes': {p: digest(p) for p in OWNED},
        'input_hashes': {p: digest(p) for p in INPUTS},
        'reproducer': 'python scripts/derive_nsc_dirac_source_phase_bound_control.py --check',
    }
    canonical = json.loads(json.dumps(record, sort_keys=True, allow_nan=False))
    return canonical, result


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--record', action='store_true')
    group.add_argument('--check', action='store_true')
    args = parser.parse_args()

    def stop(*_):
        raise TimeoutError('phase-bound control exceeded 60 CPU seconds')

    previous = signal.signal(signal.SIGPROF, stop)
    signal.setitimer(signal.ITIMER_PROF, CPU_CAP)
    started = time.process_time()
    try:
        payload, result = compute()
    except TimeoutError:
        payload, result = {
            'schema': 'NSC-DIRAC-SOURCE-PHASE-BOUND-CONTROL-v1',
            'accountable_author': 'Douglas Ek',
            'status': 'OPEN: cpu_budget; physical local gate OPEN',
            'obstruction': 'cpu_budget',
            'source_tail_error_bound': None,
            'physical_error_bound': None,
            'coefficient_error_bound': None,
            'scope': {'physical_local_gate': 'OPEN', 'new_field_or_source_runs': 0},
            'reproducer': 'python scripts/derive_nsc_dirac_source_phase_bound_control.py --check',
        }, None
    finally:
        signal.setitimer(signal.ITIMER_PROF, 0)
        signal.signal(signal.SIGPROF, previous)
    cpu = time.process_time() - started
    payload['runtime'] = {'CPU_seconds': cpu if result is None else result.cpu_seconds,
                          'CPU_cap': CPU_CAP, 'wall_process_time': cpu}
    if args.record:
        if OUTPUT.exists():
            raise FileExistsError('phase-bound record exists; use --check')
        OUTPUT.write_text(json.dumps(payload, sort_keys=True, indent=2, allow_nan=False) + '\n')
    else:
        saved = json.loads(OUTPUT.read_text())
        saved.pop('runtime', None)
        compare = dict(payload)
        compare.pop('runtime', None)
        if compare != saved:
            raise ValueError('phase-bound replay differs')
    print(json.dumps({
        'status': payload['status'],
        'obstruction': payload.get('obstruction'),
        'unit_angular_radius_pass': payload.get('unit_angular_radius_pass'),
        'cells': payload.get('cells'),
        'max_depth': payload.get('max_depth'),
        'CPU_seconds': payload['runtime']['CPU_seconds'],
        'f_z_unit_angular': payload.get('f_z_unit_angular'),
        'aggregate_phase_value_error_upper_N_beta': payload.get(
            'aggregate_phase_value_error_upper_N_beta'),
    }, indent=2, sort_keys=True))
