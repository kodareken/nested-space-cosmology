#!/usr/bin/env python3
"""Refine only the failed spatial history controls; preserve the coarse record."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy.sparse import diags

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from recursive_horizons.nsc_transmitting_history_modes import (
    RelativeModePropagator, SuppliedKSHarmonicMetric, reference_fields_on_grid,
)
from recursive_horizons.nsc_transmitting_dirac_domain import S3
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes
from derive_nsc_transmitting_history_modes import MODE, TIME, AMPLITUDES, OMEGA, load, payload

BASE = 'results/development/nsc-transmitting-history-modes.json'
OUTPUT = 'results/development/nsc-transmitting-history-refinement.json'
SOURCES = ('scripts/refine_nsc_transmitting_history_modes.py',
           'docs/nsc-transmitting-history-refinement.md')
INPUTS = (BASE, MODE, 'scripts/derive_nsc_transmitting_history_modes.py',
          'src/recursive_horizons/nsc_transmitting_history_modes.py',
          'src/recursive_horizons/nsc_ks_spacetime_variation.py',
          'results/development/nsc-massive-signed-preparation.json')
TOL = 3e-8


def sha(path):
    return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()


def read(path):
    return json.loads((ROOT/path).read_text())


def reverse_energy_projection(value):
    n = len(value)//3
    permutation = (3*np.arange(n)[::-1, None]+np.arange(3)[None, :]).ravel()
    return value[np.ix_(permutation, permutation)].conj()


def one(job):
    name, positive, negative, coarse = job
    mass, angular = float(positive['mass'].item()), float(positive['angular'].item())
    if angular <= 0 or float(negative['angular'].item()) != -angular or float(negative['mass'].item()) != mass:
        raise ValueError('the imported opposite-angular pair is required')
    mapped = np.einsum('ab,ebs->eas', S3, positive['at_zero'][::-1].conj())
    map_error = float(np.max(abs(mapped-negative['at_zero'])))
    energy_error = float(np.max(abs(negative['energies']+positive['energies'][::-1])))
    if max(map_error, energy_error) > 3e-11:
        raise ValueError('signed map does not bind these source columns; do not reuse the paired solve')
    metric = SuppliedKSHarmonicMetric(tuple(AMPLITUDES), OMEGA, TIME)
    previous = coarse
    arrays, errors, diagnostics = {}, {}, {}
    for points in (1601, 3201):
        x = np.linspace(-2., 1.2, points)
        fields = reference_fields_on_grid(positive['energies'], mass, angular, positive['at_zero'], x)
        owner = RelativeModePropagator(x, mass, angular)
        result = owner.evolve(positive['energies'], fields, np.linspace(*TIME, 65), metric)
        projection = result['mode_projection']
        key = f'{points}_64'
        arrays[key+'/projection'] = projection
        errors[key] = float(np.max(abs(projection-previous)))
        diagnostics[key] = {k: v for k, v in result.items() if not isinstance(v, np.ndarray)}
        if errors[key] <= TOL:
            break
        previous = projection
    fine = owner.evolve(positive['energies'], fields, np.linspace(*TIME, 129), metric)
    arrays['final_projection'] = fine['mode_projection']
    arrays['difference_field'] = fine['difference_characteristic_field']
    arrays['reference_field'] = fine['reference_characteristic_field']
    arrays['x'] = x
    negative_owner = RelativeModePropagator(x, mass, -angular)
    coefficients = metric.values(sum(TIME)/2, x)
    Lp, _ = owner.generator(coefficients)
    Lm, _ = negative_owner.generator(coefficients)
    C = diags(np.r_[np.ones(points), -np.ones(points)])
    remainder = Lm-C@Lp.conj()@C
    operator_map_error = float(np.max(abs(remainder.data))) if remainder.nnz else 0.
    arrays['metadata_json'] = np.frombuffer(json.dumps({
        'points': points, 'source_map_residual': map_error, 'energy_map_residual': energy_error,
        'nonlinear_operator_map_residual': operator_map_error,
        'spatial_refinements': errors,
        'temporal_refinement': float(np.max(abs(fine['mode_projection']-projection))),
        'diagnostics': diagnostics,
        'fine': {k: v for k, v in fine.items() if not isinstance(v, np.ndarray)},
    }, sort_keys=True).encode(), np.uint8)
    return name, arrays


def prepare(workers):
    base = read(BASE)
    coarse, modes = payload(base), payload(read(MODE))
    failed = {k for k, r in base['per_family'].items() if r['spatial_response_refinement'] > TOL}
    positive = sorted(k for k in failed if k.endswith('_1'))
    if len(failed) != 2*len(positive):
        raise ValueError('this refinement reuses complete failed opposite-angular pairs only')
    def data(name):
        return {k[len(name)+1:]: v for k, v in modes.items() if k.startswith(name+'/')}
    jobs = []
    for name in positive:
        partner = name[:-1]+'-1'
        if partner not in failed:
            raise ValueError('failed signed partner absent')
        # Existing independently computed nonlinear responses check the map.
        error = np.max(abs(coarse[partner+'/801_128/projection']-
                           reverse_energy_projection(coarse[name+'/801_128/projection'])))
        if error > 3e-11:
            raise ValueError('coarse nonlinear response violates the imported signed map')
        jobs.append((name, data(name), data(partner), coarse[name+'/801_64/projection']))
    arrays = {}
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for count, (name, values) in enumerate(pool.map(one, jobs), 1):
            arrays.update({name+'/'+k: v for k, v in values.items()})
            meta = json.loads(values['metadata_json'].tobytes())
            print(f'{count}/{len(jobs)} failed pairs refined: {name}, N={meta["points"]}', flush=True)
    arrays['metadata_json'] = np.frombuffer(json.dumps({
        'signature': {p: sha(p) for p in (*SOURCES, *INPUTS)},
        'base_payload': base['payload'], 'failed_families': sorted(failed),
    }, sort_keys=True).encode(), np.uint8)
    raw = deterministic_npz_bytes(arrays)
    digest = hashlib.sha256(raw).hexdigest()
    out = ROOT/f'results/development/artifacts/nsc-transmitting-history-refinement.{digest}.npz'
    out.write_bytes(raw)
    return out


def make_record(path):
    arrays = load(path)
    meta = json.loads(arrays['metadata_json'].tobytes())
    for p, expected in meta['signature'].items():
        if sha(p) != expected:
            raise ValueError('refinement dependency changed: '+p)
    base = read(BASE)
    rows = {}
    for name, old in base['per_family'].items():
        if name not in meta['failed_families']:
            rows[name] = {'source': 'unchanged authenticated coarse record; not rerun',
                          'spatial_response_refinement': old['spatial_response_refinement'],
                          'temporal_response_refinement': old['temporal_response_refinement'],
                          'norm_flux_algebra': old['norm_flux_algebra'],
                          'forcing_phase_residual': old['forcing_phase_residual'],
                          'fine_left_trace': old['fine_left_trace'], 'fine_right_trace': old['fine_right_trace'],
                          'points': 801}
            continue
        positive = name.split('_')[0]+'_1'
        d = json.loads(arrays[positive+'/metadata_json'].tobytes())
        key = f'{d["points"]}_64'
        x = arrays[positive+'/x']
        weights = np.full(len(x), (x[-1]-x[0])/(len(x)-1))
        weights[[0, -1]] *= .5
        reconstructed = arrays[positive+'/reference_field'].conj().T@(
            np.tile(weights, 2)[:, None]*arrays[positive+'/difference_field'])
        row = {
            'source': 'refined spatial PDE' if name == positive else 'imported signed intertwiner; source and operator identities checked',
            'points': d['points'], 'spatial_response_refinement': d['spatial_refinements'][key],
            'temporal_response_refinement': d['temporal_refinement'],
            'all_new_spatial_comparisons': d['spatial_refinements'],
            'source_map_residual': d['source_map_residual'], 'energy_map_residual': d['energy_map_residual'],
            'nonlinear_operator_map_residual': d['nonlinear_operator_map_residual'],
            'projection_reconstruction': float(np.max(abs(reconstructed-arrays[positive+'/final_projection']))),
            'norm_flux_algebra': d['fine']['norm_flux_algebra'],
            'forcing_phase_residual': d['fine']['forcing_phase_residual'],
            'fine_left_trace': d['fine']['left_difference_trace'],
            'fine_right_trace': d['fine']['right_difference_trace'],
        }
        rows[name] = row
    keys = ('spatial_response_refinement', 'temporal_response_refinement', 'source_map_residual',
            'energy_map_residual', 'nonlinear_operator_map_residual', 'projection_reconstruction',
            'norm_flux_algebra', 'forcing_phase_residual', 'fine_left_trace', 'fine_right_trace')
    maxima = {k: max(r.get(k, 0.) for r in rows.values()) for k in keys}
    tolerances = {k: TOL if k in keys[:2] else 3e-11 for k in keys}
    failures = {k: v for k, v in maxima.items() if v > tolerances[k]}
    return {
        'schema': 'NSC-TRANSMITTING-HISTORY-REFINEMENT-v1', 'accountable_author': 'Douglas Ek',
        'status': 'PASS: composed conditional spatial history controls' if not failures else 'OPEN',
        'source_hashes': {p: sha(p) for p in SOURCES}, 'input_hashes': {p: sha(p) for p in INPUTS},
        'payload': {'path': str(path.relative_to(ROOT)), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
                    'bytes': path.stat().st_size},
        'locked_inputs': base['locked_inputs'],
        'changed_dependency': 'spatial resolution only for previously failed controls; tolerance unchanged',
        'failed_coarse_families': meta['failed_families'],
        'preserved_coarse_record': BASE,
        'domain_and_control': base['domain'],
        'per_family': rows, 'maxima': maxima, 'tolerances': tolerances, 'failures': failures,
        'gate': {**base['gate'], 'conditional_mode_propagation': 'PASS' if not failures else 'OPEN',
                 'spatial_history_controls': 'all 63 signed families covered by inherited or refined evidence',
                 'spectral_history_convergence': 'not established by this spatial refinement',
                 'full_EndpointBranchJets': 'OPEN', 'extended_stationarity': 'OPEN'},
        'comparison': {'fields': 'all', 'float_atol': 3e-12, 'float_rtol': 3e-10,
                       'exact': 'source/input hashes, structure and scope labels'},
        'reproducer': 'python3 scripts/refine_nsc_transmitting_history_modes.py --check',
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--prepare', action='store_true')
    p.add_argument('--check', action='store_true')
    p.add_argument('--workers', type=int, default=2)
    args = p.parse_args()
    if args.prepare == args.check:
        p.error('select --prepare or --check')
    path = prepare(args.workers) if args.prepare else ROOT/read(OUTPUT)['payload']['path']
    record = make_record(path)
    if args.check:
        from derive_nsc_transmitting_boundary_binding import compare
        compare(read(OUTPUT), record)
    else:
        (ROOT/OUTPUT).write_text(json.dumps(record, indent=2, sort_keys=True)+'\n')
    print(json.dumps({'status': record['status'], 'maxima': record['maxima'], 'failures': record['failures']}, indent=2), flush=True)
    if record['failures']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
