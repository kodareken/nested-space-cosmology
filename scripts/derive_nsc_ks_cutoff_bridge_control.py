#!/usr/bin/env python3
"""Join the derived source-cutoff edge to the saved all-family raw control.

Geometry quadrature only. Original Dirac columns, source preparation,
weights and all raw result bytes are reused.
"""
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

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
sys.path.insert(0, str(ROOT/'scripts'))
import derive_nsc_local_prepared_response as P
from recursive_horizons.nsc_dirac_source_phase import formal_source_phase_coefficient
from recursive_horizons.nsc_ks_cutoff_bridge_evaluator import source_edge_from_phase
from recursive_horizons.nsc_ks_profile_identity import profile_identity
from recursive_horizons.nsc_pg_ks_metric_pullback import reference_chart

OUTPUT = ROOT/'results/development/nsc-ks-cutoff-bridge-control.json'
RAW = 'results/development/nsc-ks-retained-response-sum.json'
INVENTORY = 'results/development/nsc-ks-source-inventory.json'
FINE = 'results/development/nsc-ks-fine-trajectory.json'
BRIDGE = 'results/development/nsc-source-cutoff-bridge.json'
CAP = 30.
ORDERS = (96, 192)
OWNERS = ('scripts/derive_nsc_ks_cutoff_bridge_control.py',
    'docs/nsc-ks-cutoff-bridge-control.md',
    'src/recursive_horizons/nsc_dirac_source_phase.py',
    'src/recursive_horizons/nsc_ks_cutoff_bridge_evaluator.py',
    'src/recursive_horizons/nsc_compatible_history_geometry.py',
    'src/recursive_horizons/nsc_ks_profile_identity.py',
    'src/recursive_horizons/nsc_ks_spacetime_variation.py',
    'src/recursive_horizons/nsc_local_incoming_family.py',
    'scripts/derive_nsc_local_prepared_response.py')


def digest(path):
    return sha256((ROOT/path).read_bytes()).hexdigest()


def compute():
    input_hashes = {}
    def record(path, expected=None):
        actual = digest(path)
        if expected is not None and actual != expected:
            raise ValueError('changed bound record: '+path)
        input_hashes[path] = actual
        return json.loads((ROOT/path).read_text())
    def payload(description):
        path = description['path']
        if digest(path) != description['sha256']:
            raise ValueError('changed bound payload: '+path)
        input_hashes[path] = description['sha256']
        return ROOT/path

    raw = record(RAW)
    inventory = record(INVENTORY)
    fine = record(FINE)
    bridge = record(BRIDGE)
    if not bridge['status'].startswith('PASS:'):
        raise ValueError('the analytic cutoff bridge must be established first')
    for path, expected in {**bridge['source_hashes'], **bridge['input_hashes']}.items():
        if digest(path) != expected:
            raise ValueError('analytic bridge owner changed: '+path)
    with np.load(payload(inventory['payload']), allow_pickle=False) as f:
        meta = json.loads(f['metadata_json'].tobytes())
    with np.load(payload(fine['payload']), allow_pickle=False) as f:
        fine_meta = json.loads(f['metadata_json'].tobytes())
    required = {(p['group'], p['angular_sign']) for p in meta['panels'].values()
                if p['angular_magnitude'] != 0}
    if (len(required) != 60 or not raw['all_retained_nonzero_angular_families_included']
            or raw['required_positive_families'] != len(required)
            or raw['completed_positive_families'] != len(required)):
        raise ValueError('the complete saved60-family control is required')
    seen, signed_ids, rho_values = set(), set(), set()
    ledger = []
    for ref in raw['family_records']:
        family_record = record(ref['path'], ref['sha256'])
        group, sign = family_record['positive_family']
        key = (group, sign)
        if key in seen or key not in required:
            raise ValueError('original positive family must enter once')
        seen.add(key)
        expected = {f'{group}_{sign}_E+1', f'{group}_{-sign}_E-1'}
        if set(family_record['signed_families']) != expected:
            raise ValueError('explicit original signed partners required')
        signed_ids.update(expected)
        with np.load(payload(family_record['payload']), allow_pickle=False) as f:
            source_meta = json.loads(f['metadata_json'].tobytes())
        for batch in source_meta['upstream_batches']:
            if batch['group'] != group or batch['angular_sign'] != sign:
                raise ValueError('saved upstream channel changed')
            rho_values.add(batch['rho_up'])
        channel = meta['channels'][group]
        sign_count = len({sg for gp, sg in required if gp == group})
        multiplicity = channel['copy_count']*channel['degeneracy']/sign_count
        ell = channel['angular_eigenvalue']
        ledger.append({'group': group, 'positive_angular_sign': sign,
            'multiplicity_per_signed_family': multiplicity, 'absolute_angular': ell,
            'angular_square_weight_per_energy_sign': multiplicity*ell*ell})
    if seen != required or signed_ids != set(raw['signed_family_ids']) or len(rho_values) != 1:
        raise ValueError('same complete signed inventory and upstream slice required')
    rho_up = rho_values.pop()
    if rho_up != fine_meta['rho_up']:
        raise ValueError('saved histories must use the same upstream coordinate')
    family = P.family(P.ALPHA)
    identity = profile_identity(family)
    if identity != fine['history_identity']:
        raise ValueError('phase must use the actual saved analytic history')
    z = np.asarray(raw['z'])
    a = reference_chart(1.)[2]
    weight = sum(row['angular_square_weight_per_energy_sign'] for row in ledger)
    weights = [weight, weight]
    phases = [formal_source_phase_coefficient(family, z, angular=1., rho_up=rho_up,
               gauss_nodes=order) for order in ORDERS]
    pairs = [source_edge_from_phase(phase, weights, a) for phase in phases]
    edge, edge_tangent = pairs[-1]
    gradient = np.asarray(raw['raw_constraint_approximant'])
    tangent = np.asarray(raw['raw_history_jacobian'])
    if edge.shape != gradient.shape or edge_tangent.shape != tangent.shape:
        raise ValueError('same incoming nodes and retarded direction order required')
    step = 1e-6
    nearby = [formal_source_phase_coefficient(P.family(P.ALPHA+sign*step), z,
        angular=1., rho_up=rho_up, gauss_nodes=ORDERS[-1]) for sign in (1, -1)]
    fd = (source_edge_from_phase(nearby[0], weights, a)[0]
          - source_edge_from_phase(nearby[1], weights, a)[0])/(2*step)
    fd_error = np.max(abs(fd-edge_tangent[0]), axis=0)
    phase = phases[-1]
    return {
        'schema': 'NSC-KS-CUTOFF-BRIDGE-CONTROL-v1', 'accountable_author': 'Douglas Ek',
        'status': 'OPEN: all retained raw families plus derived edge; finite-tail and total error unbounded',
        'history_identity': identity, 'state_law': 'C_Sigma[g]=U_g C_up U_g†',
        'z': z.tolist(), 'rho_up': rho_up, 'physical_interval': [float(z[0]), float(z[-1])],
        'raw_action_gradient': gradient.tolist(), 'raw_history_jacobian': tangent.tolist(),
        'source_cutoff_edge_gradient': edge.tolist(), 'source_cutoff_edge_tangent': edge_tangent.tolist(),
        'joined_action_gradient': (gradient+edge).tolist(), 'joined_history_jacobian': (tangent+edge_tangent).tolist(),
        'raw_constraint_maxima': np.max(abs(gradient), axis=0).tolist(),
        'joined_constraint_maxima': np.max(abs(gradient+edge), axis=0).tolist(),
        'edge_maxima': np.max(abs(edge), axis=0).tolist(),
        'joined_matter_change_maxima': np.max(abs(np.asarray(raw['summed_matter_change'])+edge), axis=0).tolist(),
        'phase': {'binding': phase.binding, 'f_unit_angular': phase.f.tolist(),
            'f_z_unit_angular': phase.f_z.tolist(), 'delta_f_unit_angular': phase.delta_f.tolist(),
            'delta_f_z_unit_angular': phase.delta_f_z.tolist(), 'source_signs': list(phase.signs),
            'gauss_nodes_per_normal_panel': list(ORDERS)},
        'ledger': ledger, 'signed_angular_square_weights': weights,
        'control': {'quadrature_difference_N_beta': np.max(abs(pairs[1][0]-pairs[0][0]), axis=0).tolist(),
            'quadrature_tangent_difference_N_beta': np.max(abs(pairs[1][1]-pairs[0][1]), axis=(0, 1)).tolist(),
            'finite_difference_step': step, 'edge_tangent_residual_N_beta': fd_error.tolist(),
            'derivative_control_tolerance': 3e-8, 'derivative_control_pass': bool(max(fd_error) <= 3e-8),
            'convergence_indicators_are_bounds': False},
        'scope': {'physical_local_gate': 'OPEN', 'physical_root_claimed': False,
            'finite_source_tail_bound': None, 'phase_quadrature_bound': None,
            'full_source_and_field_error_bound': None, 'between_node_error_bound': None,
            'existing_raw_record_changed': False, 'source_reprepared': False,
            'new_field_or_source_evolutions': 0, 'new_action_term': False,
            'local_reference_coefficients_changed': False, 'Gamma_rest_assigned': False,
            'zero_angular_baseline_retained': True, 'metric_timestep': False},
        'source_hashes': {p: digest(p) for p in OWNERS}, 'input_hashes': input_hashes,
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--record', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.record and OUTPUT.exists():
        raise FileExistsError('bridge control already exists; use --check')
    def stop(*_):
        raise TimeoutError('source-cutoff geometry control exceeded30 CPU seconds')
    signal.signal(signal.SIGPROF, stop)
    signal.setitimer(signal.ITIMER_PROF, CAP)
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
        old = json.loads(OUTPUT.read_text()); old.pop('runtime')
        if result != old:
            raise ValueError('same-source cutoff control replay differs')
    print(json.dumps({k: result[k] for k in ('status', 'raw_constraint_maxima',
        'joined_constraint_maxima', 'edge_maxima', 'joined_matter_change_maxima', 'control')}, indent=2))
