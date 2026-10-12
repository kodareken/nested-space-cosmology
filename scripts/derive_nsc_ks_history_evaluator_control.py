#!/usr/bin/env python3
"""Bound the evaluator wiring check to one saved original source batch.

No source preparation is repeated. The original positive and negative
covariances and saved A_up are the controls, not fitted source parameters.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
sys.path.insert(0, str(ROOT/'scripts'))
from recursive_horizons.nsc_evolved_incoming_state import FixedSourcePreparation
from recursive_horizons.nsc_ks_batched_constraints import KSUpstreamBatch, map_negative_batch
from recursive_horizons.nsc_ks_source_inventory import ReferenceSourcePanel
from derive_nsc_evolved_incoming_constraints import coefficients_from_records

OUTPUT = ROOT/'results/development/nsc-ks-history-evaluator-control.json'
PAYLOAD = ROOT/'results/development/artifacts/nsc-ks-history-evaluator-control.npz'
OPTIONS = {'rtol': 2e-11, 'atol': 2e-13, 'max_step': .002}
CAP = 30.
OP_KEYS = ('energies', 'z', 'reference', 'difference', 'difference_z', 'tangent', 'tangent_z')
OWNERS = ('scripts/derive_nsc_ks_history_evaluator_control.py',
    'docs/nsc-ks-history-evaluator-control.md',
    'src/recursive_horizons/nsc_ks_history_evaluator.py',
    'src/recursive_horizons/nsc_ks_energy_propagator.py',
    'src/recursive_horizons/nsc_ks_batched_constraints.py',
    'src/recursive_horizons/nsc_ks_difference_envelope.py',
    'src/recursive_horizons/nsc_ks_source_envelope.py',
    'src/recursive_horizons/nsc_ks_profile_identity.py',
    'src/recursive_horizons/nsc_local_history_newton.py',
    'src/recursive_horizons/nsc_local_incoming_family.py',
    'src/recursive_horizons/nsc_evolved_incoming_constraints.py',
    'src/recursive_horizons/nsc_ks_signed_state.py')


def load_inputs():
    """Restore original low16_1 rows0:16 and their signed partner, without ODEs."""
    hashes = {}
    def bind(path, expected=None):
        path = Path(path)
        full = path if path.is_absolute() else ROOT/path
        digest = sha256(full.read_bytes()).hexdigest()
        if expected is not None and digest != expected:
            raise ValueError('bound source input changed: '+str(path))
        hashes[str(full.relative_to(ROOT))] = digest
        return full
    def record(path, expected=None):
        return json.loads(bind(path, expected).read_text())

    aggregate = record('results/development/nsc-ks-retained-response-sum.json')
    family_path = 'results/development/artifacts/nsc-ks-retained-response-sum/family-14-+1.json'
    selected = [r for r in aggregate['family_records'] if r['path'] == family_path]
    if len(selected) != 1:
        raise ValueError('one original group14 positive source receipt required')
    saved = record(family_path, selected[0]['sha256'])
    with np.load(bind(saved['payload']['path'], saved['payload']['sha256']), allow_pickle=False) as f:
        receipt = json.loads(f['metadata_json'].tobytes())['upstream_batches'][0]
        initial = np.array(f['upstream/'+receipt['panel_name']], copy=True)
    if receipt['original_panel'] != 'group14/low16_1' or receipt['rows'] != [0, 16]:
        raise ValueError('this bounded control owns exactly the saved first16 low-energy rows')

    inventory = record('results/development/nsc-ks-source-inventory.json')
    with np.load(bind(inventory['payload']['path'], inventory['payload']['sha256']), allow_pickle=False) as f:
        meta = json.loads(f['metadata_json'].tobytes())
        name = receipt['original_panel']
        p = meta['panels'][name]
        panel = ReferenceSourcePanel(name, p['group'], p['angular_sign'], p['mass'],
            p['angular_magnitude'], *(np.array(f[name+'/'+k], copy=True) for k in
            ('energies', 'weights', 'covariance', 'amplitudes_at_one')), p['provenance'])
    part = next(iter(panel.batches(16)))
    source = FixedSourcePreparation.from_signed_blocks(part.energies, part.weights, part.covariance)
    positive = KSUpstreamBatch(**receipt, source=source, initial_columns=initial)
    partner = part.negative_partner(meta['config'])
    negative_source = FixedSourcePreparation.from_signed_blocks(
        partner.energies, partner.weights, partner.covariance)
    negative = map_negative_batch(positive, negative_source)
    coefficient = record('results/development/nsc-incoming-surface-coefficients.json')
    branch = record('results/development/nsc-incoming-surface-regular-branch.json')
    baseline = record('results/development/nsc-incoming-source-update-v5.json')
    coefficients = coefficients_from_records(coefficient, branch)
    channel = meta['channels'][14]
    if channel['index'] != 14:
        raise ValueError('original channel ledger changed')
    return positive, negative, channel, coefficients, np.asarray(
        baseline['baseline']['action_gradient_approximant']), hashes


def evaluate(saved=None):
    from unittest.mock import patch
    import recursive_horizons.nsc_ks_history_evaluator as adapter
    from recursive_horizons.nsc_ks_energy_propagator import KSEnergyPropagator
    from recursive_horizons.nsc_ks_source_envelope import (
        computational_z_grid, KSEnvelopeBinding, _sample_axial_profiles, usual_axial_support)
    from recursive_horizons.nsc_ks_profile_identity import profile_identity
    from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily
    from recursive_horizons.nsc_local_history_newton import bind_history_evaluator_receipt

    positive, negative, channel, coefficients, baseline, hashes = load_inputs()
    central = np.zeros((2, 8)); central[0, 1] = .001
    direction = np.zeros((2, 8)); direction[0, :2] = (.25, 1.)
    direction[1, :2] = (.75, -.5)
    step = 1e-5
    families = [LocalIncomingFamily(central+x*step*direction) for x in (0, 1, -1)]
    solve = families[0].collocation_nodes(8)
    verify = families[0].collocation_nodes(17)
    grid = computational_z_grid(64)
    evaluator = adapter.KSHistoryEvaluator(
        [(positive, {**channel, 'angular_sign': 1}),
         (negative, {**channel, 'angular_sign': -1})],
        baseline_gradient=baseline, coefficients=coefficients,
        energy_interval=(0., .25), interpolant_degree=8, z_grid=grid,
        solve_nodes=solve, verification_nodes=verify, **OPTIONS)
    original_evolve = adapter.evolve_energy_propagator
    arrays, calls = {}, []
    def evolve(interval, degree, family, computational_z, target_z, mass, angular, rho_up, **options):
        index = len(calls)
        if index >= 3:
            raise AssertionError('at most three bounded candidate evolutions')
        identity = profile_identity(family)
        if saved is None:
            op = original_evolve(interval, degree, family, computational_z, target_z,
                                 mass, angular, rho_up, **options)
        else:
            expected = saved[1]['operators'][index]
            if identity != expected['history_identity']:
                raise ValueError('replay requested another analytic history')
            metric = family.metric()
            w, U = _sample_axial_profiles(metric.directions, computational_z)
            binding = KSEnvelopeBinding(computational_z, w, U, tuple(metric.amplitudes),
                tuple(d.inner_radius for d in metric.directions),
                tuple(d.outer_radius for d in metric.directions), usual_axial_support(),
                rho_up, 1., **OPTIONS)
            op = KSEnergyPropagator(interval, *(saved[0][f'operator{index}/{key}'] for key in OP_KEYS),
                mass, angular, rho_up, binding, expected['diagnostics'])
            if op.digest != expected['digest'] or not np.array_equal(op.z, target_z):
                raise ValueError('saved operator/profile/source restriction changed')
        for key in OP_KEYS:
            arrays[f'operator{index}/{key}'] = np.asarray(getattr(op, key))
        calls.append({'history_identity': identity, 'digest': op.digest,
                      'diagnostics': dict(op.diagnostics)})
        return op

    results = []
    expected_identity = None
    with patch.object(adapter, 'evolve_energy_propagator', side_effect=evolve):
        for family in families:
            raw = evaluator.evaluate(family, evaluator.target_nodes)
            receipt = bind_history_evaluator_receipt(raw, family,
                expected_nodes=evaluator.target_nodes,
                expected_source_identity=expected_identity)
            expected_identity = receipt.source_identity
            if receipt.history_jacobian.shape != (16, len(evaluator.target_nodes), 2):
                raise AssertionError('all16 retarded directions required')
            # These are exact selections from the same production evolution.
            for nodes in (solve, verify):
                subset = evaluator.evaluate(family, nodes)
                indices = np.searchsorted(evaluator.target_nodes, nodes)
                if not np.array_equal(subset['action_gradient'], raw['action_gradient'][indices]):
                    raise AssertionError('held-out node selection changed the assembly')
            results.append(raw)
    if len(calls) != 3 or evaluator.operator_evolution_count != 3:
        raise AssertionError('signed partners and node subsets must reuse each candidate operator')
    finite_difference = (results[1]['action_gradient']-results[2]['action_gradient'])/(2*step)
    tangent = np.tensordot(direction.ravel(), results[0]['history_jacobian'], axes=(0, 0))
    residual = np.max(abs(finite_difference-tangent), axis=0)
    geometry = np.tensordot(direction.ravel(), results[0]['local_reference_gradient_tangent'], axes=(0, 0))
    arrays['finite_difference'] = finite_difference
    arrays['retarded_tangent'] = tangent
    arrays['z'] = evaluator.target_nodes
    for index, value in enumerate(results):
        arrays[f'case{index}/action_gradient'] = value['action_gradient']
        arrays[f'case{index}/history_jacobian'] = value['history_jacobian']
    meta = {'operators': calls}
    arrays['metadata_json'] = np.frombuffer(json.dumps(meta, sort_keys=True).encode(), np.uint8)
    output = {
        'schema': 'NSC-KS-HISTORY-EVALUATOR-CONTROL-v1', 'accountable_author': 'Douglas Ek',
        'status': ('PASS' if max(residual) <= 3e-8 else 'OPEN')+': fixed-source evaluator derivative control; physical gate OPEN',
        'derivative_residual_N_beta': residual.tolist(), 'derivative_control_tolerance': 3e-8,
        'retarded_matter_tangent_maxima': np.max(abs(tangent-geometry), axis=0).tolist(),
        'raw_control_constraint_maxima': np.max(abs(results[0]['action_gradient']), axis=0).tolist(),
        'control': {'history': families[0].description(), 'direction': direction.tolist(),
            'step': step, 'source_rows_per_energy_sign': 16, 'group': 14,
            'source_panel': positive.panel_name, 'grid_nodes': 64,
            'source_energies_positive': np.unique(positive.source.energies).tolist(),
            'energy_interpolant_interval': [0., .25], 'energy_interpolant_degree': 8,
            'solve_nodes': solve.tolist(), 'verification_nodes': verify.tolist(),
            'target_union': evaluator.target_nodes.tolist(), 'retarded_directions': 16,
            'operators': calls, 'options': OPTIONS},
        'source_identity': {'positive_source': positive.source.digest,
            'negative_source': negative.source.digest, 'positive_preparation': positive.preparation_digest,
            'negative_preparation': negative.preparation_digest},
        'scope': {'physical_local_gate': 'OPEN', 'source_preparation_repeated': False,
            'positive_energy_result_doubled': False, 'seeds_are_solver_controls': True,
            'removed_state_derivative': False, 'source_cutoff_coincidence_bridge': 'OPEN',
            'full_error_budget': None, 'metric_timestep': False, 'physical_optimizer_run': False},
        'source_hashes': {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in OWNERS},
        'input_hashes': hashes,
    }
    return output, arrays


if __name__ == '__main__':
    import argparse
    import signal
    import time
    from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--record', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.record and (OUTPUT.exists() or PAYLOAD.exists()):
        raise FileExistsError('control already exists; use --check')
    saved = None
    if args.check:
        previous = json.loads(OUTPUT.read_text())
        for path, expected in {**previous['source_hashes'], **previous['input_hashes']}.items():
            if sha256((ROOT/path).read_bytes()).hexdigest() != expected:
                raise ValueError('control dependency changed: '+path)
        if sha256(PAYLOAD.read_bytes()).hexdigest() != previous['payload']['sha256']:
            raise ValueError('control payload changed')
        with np.load(PAYLOAD, allow_pickle=False) as f:
            values = {k: np.array(f[k], copy=True) for k in f.files}
        saved = values, json.loads(values['metadata_json'].tobytes())
    def stop(*_):
        raise TimeoutError('bounded evaluator control exceeded30 CPU seconds')
    signal.signal(signal.SIGPROF, stop)
    signal.setitimer(signal.ITIMER_PROF, CAP)
    start = time.process_time()
    try:
        result, arrays = evaluate(saved)
    finally:
        signal.setitimer(signal.ITIMER_PROF, 0)
    cpu = time.process_time()-start
    payload = deterministic_npz_bytes(arrays)
    result['payload'] = {'path': str(PAYLOAD.relative_to(ROOT)), 'sha256': sha256(payload).hexdigest(), 'bytes': len(payload)}
    if args.record:
        result['runtime'] = {'CPU_seconds': cpu, 'CPU_cap': CAP, 'new_operator_solves': 3}
        PAYLOAD.write_bytes(payload)
        OUTPUT.write_text(json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+'\n')
    else:
        previous.pop('runtime')
        if result != previous:
            raise ValueError('saved-operator evaluator replay differs')
    print(json.dumps({k: result[k] for k in ('status', 'derivative_residual_N_beta',
        'retarded_matter_tangent_maxima', 'raw_control_constraint_maxima')}, indent=2))
