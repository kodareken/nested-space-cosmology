#!/usr/bin/env python3
"""One bounded coupled-history pilot using both original angular signs of group14."""
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
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src')); sys.path.insert(0, str(ROOT/'scripts'))
import recursive_horizons.nsc_ks_history_evaluator as H
from recursive_horizons.nsc_ks_cutoff_bridge_evaluator import KSCutoffBridgeEvaluator
from recursive_horizons.nsc_ks_retained_upstream_archive import RetainedUpstreamArchive
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily
from recursive_horizons.nsc_local_history_newton import bind_history_evaluator_receipt
from recursive_horizons.nsc_ks_energy_propagator import KSEnergyPropagator
from recursive_horizons.nsc_ks_source_envelope import computational_z_grid, KSEnvelopeBinding, _sample_axial_profiles, usual_axial_support
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes
from derive_nsc_evolved_incoming_constraints import coefficients_from_records

OUTPUT = ROOT/'results/development/nsc-ks-coupled-history-pilot.json'
PAYLOAD = ROOT/'results/development/artifacts/nsc-ks-coupled-history-pilot.npz'
OPTIONS = {'rtol': 2e-11, 'atol': 2e-13, 'max_step': .002}
CAP = 30.
KEYS = ('energies', 'z', 'reference', 'difference', 'difference_z', 'tangent', 'tangent_z')
OWNERS = ('scripts/derive_nsc_ks_coupled_history_pilot.py', 'docs/nsc-ks-coupled-history-pilot.md',
    'src/recursive_horizons/nsc_ks_retained_upstream_archive.py',
    'src/recursive_horizons/nsc_ks_history_evaluator.py',
    'src/recursive_horizons/nsc_ks_cutoff_bridge_evaluator.py',
    'src/recursive_horizons/nsc_dirac_source_phase.py',
    'src/recursive_horizons/nsc_ks_energy_propagator.py',
    'src/recursive_horizons/nsc_ks_difference_envelope.py',
    'src/recursive_horizons/nsc_ks_batched_constraints.py')


def evaluate(saved=None):
    archive = RetainedUpstreamArchive(ROOT)
    coefficient = archive._record('results/development/nsc-incoming-surface-coefficients.json')
    branch = archive._record('results/development/nsc-incoming-surface-regular-branch.json')
    source = archive._record('results/development/nsc-incoming-source-update-v5.json')
    coeff = coefficients_from_records(coefficient, branch)
    baseline = np.asarray(source['baseline']['action_gradient_approximant'])
    # A geometric initial guess only. The actual source is evolved below;
    # the root of the frozen-C0 polynomial is not promoted to a physical state.
    slope = -baseline[1]/coeff['F'][0]
    U = -(baseline[0]+coeff['C']*slope*slope)/coeff['A0']
    control = np.zeros((2, 8)); control[0, 1] = .06*slope; control[1, 0] = U
    family = LocalIncomingFamily(control)
    entries = tuple(entry for key in ((14, -1), (14, 1)) for entry in archive.family_entries(key))
    raw = H.KSHistoryEvaluator(entries, baseline_gradient=baseline, coefficients=coeff,
        energy_interval=(0., 160.), interpolant_degree=32, z_grid=computational_z_grid(64),
        solve_nodes=family.collocation_nodes(8), verification_nodes=family.collocation_nodes(17),
        ell0_channels=archive.analytic_zero_channels(), **OPTIONS)
    joined = KSCutoffBridgeEvaluator(raw, phase_gauss_nodes=32)
    arrays, operators = {}, []
    actual_evolve = H.evolve_energy_propagator
    def evolve(interval, degree, candidate, grid, target, mass, angular, rho_up, **options):
        index = len(operators)
        if index >= 2:
            raise AssertionError('pilot permits exactly two signed angular operator solves')
        if saved is None:
            op = actual_evolve(interval, degree, candidate, grid, target, mass, angular, rho_up, **options)
        else:
            expected = saved[1]['operators'][index]
            if (mass, angular, rho_up) != tuple(expected[k] for k in ('mass', 'angular', 'rho_up')):
                raise ValueError('saved pilot channel changed')
            metric = candidate.metric(); w, u = _sample_axial_profiles(metric.directions, grid)
            binding = KSEnvelopeBinding(grid, w, u, tuple(metric.amplitudes),
                tuple(d.inner_radius for d in metric.directions), tuple(d.outer_radius for d in metric.directions),
                usual_axial_support(), rho_up, 1., **OPTIONS)
            op = KSEnergyPropagator(interval, *(saved[0][f'operator{index}/{key}'] for key in KEYS),
                mass, angular, rho_up, binding, expected['diagnostics'])
            if op.digest != expected['digest'] or not np.array_equal(op.z, target):
                raise ValueError('saved pilot operator/history binding changed')
        operators.append({'mass': mass, 'angular': angular, 'rho_up': rho_up,
            'digest': op.digest, 'diagnostics': dict(op.diagnostics)})
        for key in KEYS:
            arrays[f'operator{index}/{key}'] = np.asarray(getattr(op, key))
        return op
    with patch.object(H, 'evolve_energy_propagator', side_effect=evolve):
        result = joined.evaluate(family)
    receipt = bind_history_evaluator_receipt(result, family, expected_nodes=joined.target_nodes)
    if receipt.history_jacobian.shape != (16, len(joined.target_nodes), 2) or len(operators) != 2:
        raise AssertionError('all16 directions and both original angular signs required')
    for key in ('z', 'action_gradient', 'history_jacobian', 'raw_action_gradient', 'raw_history_jacobian',
                'source_cutoff_edge_gradient', 'source_cutoff_edge_history_tangent'):
        arrays[key] = result[key]
    arrays['metadata_json'] = np.frombuffer(json.dumps({'operators': operators}, sort_keys=True).encode(), np.uint8)
    matter_tangent = result['history_jacobian']-result['local_reference_gradient_tangent']
    report = {
        'schema': 'NSC-KS-COUPLED-HISTORY-PILOT-v1', 'accountable_author': 'Douglas Ek',
        'status': 'OPEN: coupled-history throughput control; only original group14 source responses included',
        'history': family.description(), 'profile_identity': result['profile_identity'],
        'seed': {'purpose': 'geometric solver control only', 'w_at_center': 0.,
            'w_z_at_center': float(slope), 'w_zz_at_center': 0., 'U_constant': float(U),
            'source_evolved_from_original_upstream': True, 'physical_root_claimed': False},
        'source_rows_per_energy_sign': sum(len(b.source.energies)//3 for b, _ in entries if b.energy_sign > 0),
        'signed_families': sorted(result['family_records']), 'signed_batches': len(entries),
        'grid_nodes': 64, 'energy_interpolation_degree': 32, 'energy_interval': [0., 160.],
        'target_nodes': joined.target_nodes.tolist(), 'retarded_directions': 16, 'options': OPTIONS,
        'constraint_maxima': np.max(abs(result['action_gradient']), axis=0).tolist(),
        'retarded_matter_plus_edge_jacobian_maxima': np.max(abs(matter_tangent), axis=(0, 1)).tolist(),
        'operators': operators,
        'scope': {'physical_local_gate': 'OPEN', 'full_source_coverage': False,
            'source_preparation_repeated': False, 'source_cutoff_edge_and_tangent_included': True,
            'finite_tail_and_full_error_bound': None, 'between_node_error_bound': None,
            'physical_optimizer_run': False, 'metric_timestep': False, 'new_action_term': False},
        'source_hashes': {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in OWNERS},
        'input_hashes': archive.input_hashes,
    }
    return report, arrays


if __name__ == '__main__':
    p = argparse.ArgumentParser(); mode = p.add_mutually_exclusive_group(required=True)
    mode.add_argument('--record', action='store_true'); mode.add_argument('--check', action='store_true')
    args = p.parse_args()
    if args.record and (OUTPUT.exists() or PAYLOAD.exists()):
        raise FileExistsError('coupled pilot exists; use --check')
    saved = None
    if args.check:
        previous = json.loads(OUTPUT.read_text())
        for path, expected in {**previous['source_hashes'], **previous['input_hashes']}.items():
            if sha256((ROOT/path).read_bytes()).hexdigest() != expected:
                raise ValueError('pilot dependency changed: '+path)
        if sha256(PAYLOAD.read_bytes()).hexdigest() != previous['payload']['sha256']:
            raise ValueError('pilot payload changed')
        with np.load(PAYLOAD, allow_pickle=False) as f:
            a = {k: np.array(f[k], copy=True) for k in f.files}
        saved = a, json.loads(a['metadata_json'].tobytes())
    def stop(*_):
        raise TimeoutError('two-operator coupled pilot exceeded30 CPU seconds')
    signal.signal(signal.SIGPROF, stop); signal.setitimer(signal.ITIMER_PROF, CAP)
    start = time.process_time()
    try:
        result, arrays = evaluate(saved)
    finally:
        signal.setitimer(signal.ITIMER_PROF, 0)
    cpu = time.process_time()-start
    data = deterministic_npz_bytes(arrays)
    result['payload'] = {'path': str(PAYLOAD.relative_to(ROOT)), 'sha256': sha256(data).hexdigest(), 'bytes': len(data)}
    if args.record:
        result['runtime'] = {'CPU_seconds': cpu, 'CPU_cap': CAP, 'new_operator_solves': 2}
        PAYLOAD.write_bytes(data); OUTPUT.write_text(json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+'\n')
    else:
        previous.pop('runtime')
        if result != previous:
            raise ValueError('coupled pilot replay differs')
    print(json.dumps({k: result[k] for k in ('status', 'seed', 'constraint_maxima',
        'retarded_matter_plus_edge_jacobian_maxima')}, indent=2))
    print(json.dumps({'CPU_seconds': cpu, 'new_operator_solves': 2 if args.record else 0}))
