#!/usr/bin/env python3
"""Compute the new nonlinear common-representation Gaussian source."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT/'src'))
from recursive_horizons.nsc_nonlinear_ks_source import record_common_history
from recursive_horizons.nsc_transmitting_history_modes import SuppliedKSHarmonicMetric
from recursive_horizons.nsc_transmitting_history_jets import FourthOrderModePropagator, relative_branch_ctp_control
from recursive_horizons.nsc_transmitting_dirac_domain import MODE_TO_CURRENT, S3
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes
from derive_nsc_transmitting_history_jets import family, reflected_overlap, HISTORY
from derive_nsc_transmitting_history_modes import MODE, TIME, AMPLITUDES, OMEGA, payload, read, load

OUTPUT = 'results/development/nsc-nonlinear-ks-source.json'
JETS = 'results/development/nsc-transmitting-history-jets.json'
SOURCES = ('src/recursive_horizons/nsc_nonlinear_ks_source.py', 'scripts/derive_nsc_nonlinear_ks_source.py',
           'tests/test_nsc_nonlinear_ks_source.py', 'docs/nsc-nonlinear-ks-source.md')
INPUTS = (MODE, HISTORY, JETS, 'results/development/nsc-common-ks-trace.json',
          'src/recursive_horizons/nsc_transmitting_history_modes.py',
          'src/recursive_horizons/nsc_transmitting_history_jets.py',
          'src/recursive_horizons/nsc_ks_spacetime_variation.py',
          'src/recursive_horizons/nsc_pg_ks_metric_pullback.py',
          'src/recursive_horizons/nsc_transmitting_resolvent.py',
          'src/recursive_horizons/nsc_lorentzian.py',
          'src/recursive_horizons/nsc_transmitting_dirac_domain.py',
          'docs/nsc-common-ks-trace.md')
CONFIGS = ((401, 64), (801, 64), (801, 128))


def sha(p): return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()


def one(job):
    name, data, x0, reference, endpoint = job
    E = data['energies']; char = reference.reshape(2, len(x0), len(E), 3).transpose(1, 2, 0, 3)
    fields = np.einsum('ab,nebc->neac', MODE_TO_CURRENT, char)
    arrays = {}; diagnostics = {}
    for points, steps in CONFIGS:
        stride = (len(x0)-1)//(points-1); x = x0[::stride]
        owner = FourthOrderModePropagator(x, data['mass'].item(), data['angular'].item())
        result = record_common_history(owner, E, fields[::stride], np.linspace(*TIME, steps+1),
                                       SuppliedKSHarmonicMetric(tuple(AMPLITUDES), OMEGA, TIME))
        key = f'{points}_{steps}'
        for field in ('weak_KS', 'weak_PG'): arrays[key+'/'+field] = result[field]
        diagnostics[key] = {k: result[k] for k in ('flux_residual', 'causal_buffer', 'phase_residual')}
        if endpoint is not None and (points, steps) == CONFIGS[0]:
            diagnostics[key]['old_endpoint_field'] = float(np.max(abs(result['evolved_field']-endpoint)))
        if (points, steps) == CONFIGS[-1]:
            for field in ('record_times', 'KS_slice_fields', 'KS_slice_z_derivatives'):
                arrays[field] = result[field]
            if name == '14_1': arrays['time_weak_KS'] = result['time_weak_KS']
            energy = np.repeat(E, 3)
            diagnostics[key]['local_derivative_minus_incoming_label'] = float(np.max(abs(
                result['KS_slice_z_derivatives']+1j*energy[None, None, :]*result['KS_slice_fields'])))
    arrays['diagnostics_json'] = np.frombuffer(json.dumps(diagnostics, sort_keys=True).encode(), np.uint8)
    return name, arrays


def prepare(workers):
    modes, hist, old = (payload(read(p)) for p in (MODE, HISTORY, JETS))
    names = sorted({k.split('/')[0] for k in modes if '/' in k and k.split('/')[0].endswith('_1')})
    jobs = [(name, family(modes, name), hist[name+'/x'], hist[name+'/reference_field'],
             old['14_1/FD_control/evolved_field'] if name == '14_1' else None) for name in names]
    arrays = {}
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for i, (name, rows) in enumerate(pool.map(one, jobs), 1):
            arrays.update({name+'/'+k: v for k, v in rows.items()})
            if i % 4 == 0 or i == len(names): print(f'{i}/{len(names)} nonlinear KS source families; signed partners reused', flush=True)
    meta = {'signature': {p: sha(p) for p in (*SOURCES, *INPUTS)},
            'mode_payload': read(MODE)['payload'], 'jet_payload': read(JETS)['payload'],
            'history_payload': read(HISTORY)['payload']}
    arrays['metadata_json'] = np.frombuffer(json.dumps(meta, sort_keys=True).encode(), np.uint8)
    raw = deterministic_npz_bytes(arrays); digest = hashlib.sha256(raw).hexdigest()
    path = ROOT/f'results/development/artifacts/nsc-nonlinear-ks-source.{digest}.npz'; path.write_bytes(raw)
    return path


def make_record(path):
    a = load(path); meta = json.loads(a['metadata_json'].tobytes())
    for p, expected in meta['signature'].items():
        if sha(p) != expected: raise ValueError('nonlinear KS dependency changed: '+p)
    modes = payload({'payload': meta['mode_payload']}); old = payload({'payload': meta['jet_payload']})
    names = sorted({k.split('/')[0] for k in modes if '/' in k})
    rows, maxima = {}, {}
    for name in names:
        d = family(modes, name); positive = name.split('_')[0]+'_1'
        diagnostics = json.loads(a[positive+'/diagnostics_json'].tobytes())
        source, matrix = {}, {}; case = {}
        if name != positive:
            pd = family(modes, positive)
            case['signed_input_identity'] = float(np.max(abs(d['at_zero']-np.einsum('ab,ebs->eas', S3, pd['at_zero'][::-1].conj()))))
        for points, steps in CONFIGS:
            key = f'{points}_{steps}'
            ks, pg, target = a[positive+'/'+key+'/weak_KS'], a[positive+'/'+key+'/weak_PG'], 1j*old[positive+'/'+key+'/overlap']
            if name != positive:
                ks, pg, target = (-reflected_overlap(v) for v in (ks, pg, target))
            result = relative_branch_ctp_control(d['energies'], d['weights'], d['source'], d['projector'], -1j*ks)
            target_result = relative_branch_ctp_control(d['energies'], d['weights'], d['source'], d['projector'], -1j*target)
            source[key] = [float(v['derivative'].real) for v in result]
            matrix[key] = {'KS_PG': float(np.max(abs(ks-pg))), 'KS_endpoint': float(np.max(abs(ks-target))),
                           'source_endpoint': max(float(abs(v['derivative']-w['derivative'])) for v, w in zip(result, target_result))}
        fine = diagnostics['801_128']
        case.update({'source_by_configuration': source, 'matrix_by_configuration': matrix,
                     'KS_PG_matrix': matrix['801_128']['KS_PG'], 'KS_endpoint_matrix': matrix['801_128']['KS_endpoint'],
                     'Gaussian_source_endpoint': matrix['801_128']['source_endpoint'],
                     'spatial_source_refinement': float(np.max(abs(np.array(source['801_64'])-source['401_64']))),
                     'time_source_refinement': float(np.max(abs(np.array(source['801_128'])-source['801_64']))),
                     'flux': fine['flux_residual'], 'phase': fine['phase_residual'], 'causal_buffer': fine['causal_buffer'],
                     'local_derivative_minus_incoming_label': fine['local_derivative_minus_incoming_label']})
        if 'old_endpoint_field' in diagnostics['401_64']: case['old_endpoint_field'] = diagnostics['401_64']['old_endpoint_field']
        rows[name] = case
        for key in ('KS_PG_matrix', 'KS_endpoint_matrix', 'Gaussian_source_endpoint', 'spatial_source_refinement',
                    'time_source_refinement', 'flux', 'phase', 'signed_input_identity', 'old_endpoint_field'):
            if key in case: maxima[key] = max(maxima.get(key, 0.), case[key])
    tolerance = {k: 3e-11 if k in ('flux', 'phase', 'signed_input_identity', 'old_endpoint_field') else 3e-8 for k in maxima}
    failures = {k: v for k, v in maxima.items() if v > tolerance[k]}
    return {
        'schema': 'NSC-NONLINEAR-KS-SOURCE-v1', 'accountable_author': 'Douglas Ek',
        'status': 'PASS: nonlinear common-representation Gaussian source on inherited spectral sample' if not failures else 'OPEN',
        'source_hashes': {p: sha(p) for p in SOURCES}, 'input_hashes': {p: sha(p) for p in INPUTS},
        'payload': {'path': str(path.relative_to(ROOT)), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(), 'bytes': path.stat().st_size},
        'locked_inputs': read(MODE)['locked_inputs'], 'per_family': rows, 'maxima': maxima, 'tolerances': tolerance, 'failures': failures,
        'domain': {'history': 'same nonzero supplied raw-KS amplitudes and PG pulse as the locked nonlinear jets',
            'time_coordinates': list(TIME), 'amplitudes': AMPLITUDES.tolist(), 'omega': OMEGA,
            'frame': 'Psi_KS=diag(sqrt(-v_plus),sqrt(-v_minus))*chi_PG, actual metric',
            'derivative': 'd_z Psi=d_tau(frame)*chi+frame*d_tau(chi), actual PG equation and inflow',
            'trace_measure': 'dT dz=d tau d rho/a0 of the fixed chart',
            'source_sign': 'delta Gamma_G=-Tr(C_source*integrated weak H vertex); same dE/(2*pi) sample',
            'initial_frequency_is_local_momentum': False, 'phases_multiplied_twice': False,
            'field_recording': 'rho0 slice samples after full spatial evolution; not a closed few-mode bulk',
            'signed_families': len(rows), 'propagated_primary_families': 33,
            'full_spectral_integral': False, 'angular_copy_sum': False, 'old_generators_rerun': False},
        'gate': {'nonlinear_common_Gaussian_trace': 'PASS' if not failures else 'OPEN',
            'common_regulated_reference_trace': 'OPEN: same kernel regulator, boundary terms and spectral completion required',
            'full_raw_heat_equivalence': 'not an execution prerequisite of the selected canonical realization',
            'complete_source': 'OPEN', 'physical_EndpointBranchJets': 'OPEN', 'B2': 'OPEN', 'V_c': None,
            'extended_stationarity': 'OPEN', 'absolute_stress': None, 'nulls': None, 'updated_constraints': None,
            'metric_timestep': False, 'physical_history_selected': False, 'Z3': 'OUT OF SCOPE',
            'Gamma_rest_extra_independent_term': 'none introduced', 'homogeneous_nonexistence': 'preserved, not rerun'},
        'comparison': {'fields': 'all', 'float_atol': 3e-12, 'float_rtol': 3e-10, 'exact': 'hashes, structure, scope'},
        'reproducer': 'python3 scripts/derive_nsc_nonlinear_ks_source.py --check',
    }


def main():
    p = argparse.ArgumentParser(); p.add_argument('--prepare', action='store_true'); p.add_argument('--check', action='store_true')
    p.add_argument('--workers', type=int, default=3); args = p.parse_args()
    if args.prepare == args.check: p.error('choose prepare or check')
    path = prepare(args.workers) if args.prepare else ROOT/read(OUTPUT)['payload']['path']
    result = make_record(path)
    if args.check:
        from derive_nsc_transmitting_boundary_binding import compare
        compare(read(OUTPUT), result)
    else: (ROOT/OUTPUT).write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print(json.dumps({k: result[k] for k in ('status', 'maxima', 'failures')}, indent=2), flush=True)
    if result['failures']: raise SystemExit(1)


if __name__ == '__main__': main()
