#!/usr/bin/env python3
"""Directed phase-value accuracy at every node of the new coupled-history pilot."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import argparse
from hashlib import sha256
import json
from pathlib import Path
import signal
import sys
import time

import numpy as np
from flint import arb, ctx

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from recursive_horizons.nsc_dirac_source_phase import formal_source_phase_coefficient
from recursive_horizons.nsc_dirac_source_phase_bound import directed_source_phase_z_enclosure, pack_real_ball, pack_upper
from recursive_horizons.nsc_ks_cutoff_bridge_evaluator import source_edge_from_phase
from recursive_horizons.nsc_ks_ball_geometry import background_series
from recursive_horizons.nsc_ks_profile_identity import profile_identity
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily

OUTPUT = ROOT/'results/development/nsc-ks-coupled-phase-accuracy.json'
PILOT = 'results/development/nsc-ks-coupled-history-pilot.json'
BRANCH = 'results/development/nsc-incoming-surface-regular-branch.json'
LEDGER = 'results/development/nsc-ks-cutoff-bridge-control.json'
BITS, CAP = 120, 60.
OWNERS = ('scripts/derive_nsc_ks_coupled_phase_accuracy.py', 'docs/nsc-ks-coupled-phase-accuracy.md',
    'src/recursive_horizons/nsc_dirac_source_phase_bound.py',
    'src/recursive_horizons/nsc_dirac_source_phase.py',
    'src/recursive_horizons/nsc_ks_cutoff_bridge_evaluator.py',
    'src/recursive_horizons/nsc_ks_ball_geometry.py',
    'src/recursive_horizons/nsc_ks_ball_operator.py')


def digest(path):
    return sha256((ROOT/path).read_bytes()).hexdigest()


def compute():
    pilot = json.loads((ROOT/PILOT).read_text())
    branch = json.loads((ROOT/BRANCH).read_text())
    ledger = json.loads((ROOT/LEDGER).read_text())
    family = LocalIncomingFamily(pilot['history']['coefficients'])
    if profile_identity(family) != pilot['profile_identity'] or family.radius_lower_bound() <= 0:
        raise ValueError('same positive-radius coupled history required')
    rho_values = {p['rho_up'] for p in pilot['operators']}
    if len(rho_values) != 1:
        raise ValueError('one original upstream slice required')
    rho_up = rho_values.pop()
    z = np.asarray(pilot['target_nodes'])
    if not np.isfinite(z).all() or not 1.03 <= rho_up <= 33/32:
        raise ValueError('finite incoming nodes and owned preparation slab required')
    weights = ledger['signed_angular_square_weights']
    if weights != [107952., 107952.]:
        raise ValueError('same original full-source angular ledger required')
    a_numeric = float(np.mean(branch['branch']['intrinsic_intervals']['a']))
    started = time.process_time()
    rows, true_gradients = [], []
    with ctx.workprec(BITS):
        a_true = background_series(arb(1), 0).a[0]
        factor = arb(107952)/(2*arb.pi())
        for index, point in enumerate(z):
            remaining = CAP-(time.process_time()-started)
            if remaining <= 0:
                raise TimeoutError('coupled phase enclosure exceeded total60 CPU seconds')
            enclosure = directed_source_phase_z_enclosure(family, float(point),
                rho_up=rho_up, bits=BITS, cpu_limit=remaining)
            if enclosure.obstruction is not None or len(enclosure.balls) != 2:
                raise ArithmeticError(f'phase bound remains OPEN at node{index}: {enclosure.obstruction}')
            plus, minus = enclosure.balls
            truth = (factor*(plus-minus)/a_true, -factor*(plus+minus))
            true_gradients.append(truth)
            rows.append({'node': index, 'z': float(point),
                'f_z_unit_angular': [pack_real_ball(v) for v in enclosure.balls],
                'edge_gradient_balls_N_beta': [pack_real_ball(v) for v in truth],
                'cells': list(enclosure.cells), 'max_depth': list(enclosure.max_depth)})
        attempts = []
        selected = None
        # A value claim is compared directly to the true-gradient ball.
        # Its error is not merely the radius around another approximation.
        for count in (32, 192):
            phase = formal_source_phase_coefficient(family, z, angular=1., rho_up=rho_up, gauss_nodes=count)
            gradient, tangent = source_edge_from_phase(phase, weights, a_numeric)
            errors = [[abs(arb(float(value))-truth).upper() for value, truth in zip(row, exact)]
                      for row, exact in zip(gradient, true_gradients)]
            maximum = [arb(0), arb(0)]
            for row in errors:
                for k in (0, 1):
                    maximum[k] = maximum[k].max(row[k]).upper()
            passed = all(value <= arb(1)/10**11 for value in maximum)
            attempts.append({'gauss_nodes_per_normal_panel': count,
                'value_error_upper_N_beta': [pack_upper(v) for v in maximum],
                'phase_value_control_pass': passed})
            if passed:
                selected = {'gauss_nodes_per_normal_panel': count, 'phase_binding': phase.binding,
                    'edge_gradient': gradient.tolist(), 'edge_history_tangent': tangent.tolist(),
                    'pointwise_value_error_upper_N_beta': [[pack_upper(v) for v in row] for row in errors],
                    'value_error_upper_N_beta': [pack_upper(v) for v in maximum]}
                break
    return {
        'schema': 'NSC-KS-COUPLED-PHASE-ACCURACY-v1', 'accountable_author': 'Douglas Ek',
        'status': ('PASS' if selected is not None else 'OPEN')+': phase value on all coupled pilot nodes; other errors and local gate OPEN',
        'history_identity': pilot['profile_identity'], 'target_nodes': z.tolist(), 'rho_up': rho_up,
        'bits': BITS, 'phase_component_tolerance': 1e-11,
        'signed_angular_square_weights': weights, 'rows': rows, 'quadrature_attempts': attempts,
        'selected': selected,
        'scope': {'physical_local_gate': 'OPEN', 'all_declared_nodes_covered': True,
            'between_node_phase_error_bound': None, 'phase_tangent_error_bound': None,
            'finite_source_tail_bound': None, 'source_and_field_error_bound': None,
            'value_error_is_direct_distance_to_true_ball': True,
            'new_field_or_source_evolutions': 0, 'physical_optimizer_run': False, 'metric_timestep': False},
        'source_hashes': {p: digest(p) for p in OWNERS},
        'input_hashes': {p: digest(p) for p in (PILOT, BRANCH, LEDGER)},
    }


if __name__ == '__main__':
    p = argparse.ArgumentParser(); mode = p.add_mutually_exclusive_group(required=True)
    mode.add_argument('--record', action='store_true'); mode.add_argument('--check', action='store_true')
    args = p.parse_args()
    if args.record and OUTPUT.exists():
        raise FileExistsError('coupled phase accuracy exists; use --check')
    def stop(*_):
        raise TimeoutError('coupled phase accuracy exceeded60 CPU seconds')
    signal.signal(signal.SIGPROF, stop); signal.setitimer(signal.ITIMER_PROF, CAP)
    start = time.process_time()
    try:
        result = compute()
    finally:
        signal.setitimer(signal.ITIMER_PROF, 0)
    cpu = time.process_time()-start
    if args.record:
        result['runtime'] = {'CPU_seconds': cpu, 'CPU_cap': CAP}
        OUTPUT.write_text(json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+'\n')
    else:
        previous = json.loads(OUTPUT.read_text()); previous.pop('runtime')
        if result != previous:
            raise ValueError('coupled phase value-bound replay differs')
    print(json.dumps({'status': result['status'], 'nodes': len(result['rows']),
        'quadrature_attempts': result['quadrature_attempts'], 'CPU_seconds': cpu}, indent=2))
