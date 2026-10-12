#!/usr/bin/env python3
"""Nonlinear spatial metric jets and their existing Gaussian CTP contraction."""
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
from recursive_horizons.nsc_transmitting_history_modes import SuppliedKSHarmonicMetric
from recursive_horizons.nsc_transmitting_history_jets import (
    FourthOrderModePropagator, evolve_field_jets, relative_branch_ctp_control,
    raw_KS_amplitude_directions, generator_direction, SBP4_COEFFICIENT_SOURCE,
)
from recursive_horizons.nsc_transmitting_dirac_domain import MODE_TO_CURRENT, S3
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes
from derive_nsc_transmitting_history_modes import MODE, TIME, AMPLITUDES, OMEGA, load, payload

HISTORY = 'results/development/nsc-transmitting-history-modes.json'
OUTPUT = 'results/development/nsc-transmitting-history-jets.json'
SOURCES = ('src/recursive_horizons/nsc_transmitting_history_jets.py',
           'scripts/derive_nsc_transmitting_history_jets.py',
           'tests/test_nsc_transmitting_history_jets.py',
           'docs/nsc-transmitting-history-jets.md')
INPUTS = (MODE, HISTORY,
          'src/recursive_horizons/nsc_transmitting_history_modes.py',
          'src/recursive_horizons/nsc_pg_ks_metric_pullback.py',
          'src/recursive_horizons/nsc_ks_spacetime_variation.py',
          'src/recursive_horizons/nsc_transmitting_ctp_variation.py',
          'src/recursive_horizons/nsc_transmitting_dirac_domain.py',
          'scripts/derive_nsc_transmitting_history_modes.py',
          'scripts/refine_nsc_transmitting_history_modes.py',
          'docs/nsc-declared-action-scope.md')
CONFIGS = ((401, 64), (801, 64), (801, 128))


def sha(path):
    return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()


def read(path):
    return json.loads((ROOT/path).read_text())


def family(data, name):
    return {k[len(name)+1:]: v for k, v in data.items() if k.startswith(name+'/')}


def reflected_overlap(value):
    n = value.shape[-1]//3
    p = (3*np.arange(n)[::-1, None]+np.arange(3)[None, :]).ravel()
    return value[:, p][:, :, p].conj()


def one(job):
    name, d, partner, x0, reference = job
    E = d['energies']; mass, angular = float(d['mass'].item()), float(d['angular'].item())
    char = reference.reshape(2, len(x0), len(E), 3).transpose(1, 2, 0, 3)
    fields0 = np.einsum('ab,nebc->neac', MODE_TO_CURRENT, char)
    provider = SuppliedKSHarmonicMetric(tuple(AMPLITUDES), OMEGA, TIME)
    arrays, diagnostics = {}, {}
    if partner is not None:
        residual = float(np.max(abs(partner['at_zero']-
                          np.einsum('ab,ebs->eas', S3, d['at_zero'][::-1].conj()))))
        if residual > 3e-11 or np.max(abs(partner['energies']+E[::-1])) > 3e-11:
            raise ValueError('signed source columns cannot be identified')
        diagnostics['signed_source_residual'] = residual
    for points, steps in CONFIGS:
        stride = (len(x0)-1)//(points-1)
        if stride < 1 or (len(x0)-1) % (points-1):
            raise ValueError('reuse a nested authenticated reference grid')
        x, fields = x0[::stride], fields0[::stride]
        owner = FourthOrderModePropagator(x, mass, angular)
        result = evolve_field_jets(owner, E, fields, np.linspace(*TIME, steps+1), provider)
        key = f'{points}_{steps}'
        arrays[key+'/overlap'] = result['relative_branch_overlap_jets']
        diagnostics[key] = {k: v for k, v in result.items() if not isinstance(v, np.ndarray)}
        diagnostics[key]['SBP_residual'] = owner.sbp_residual
        if partner is not None and (points, steps) == CONFIGS[-1]:
            metric = provider.values(sum(TIME)/2, x)
            negative = FourthOrderModePropagator(x, mass, -angular)
            sign = diags(np.r_[np.ones(len(x)), -np.ones(len(x))])
            Lp, _ = owner.generator(metric); Lm, _ = negative.generator(metric)
            diff = Lm-sign@Lp.conj()@sign
            errors = [float(np.max(abs(diff.data))) if diff.nnz else 0.]
            for direction in raw_KS_amplitude_directions(provider, sum(TIME)/2, x, metric):
                dp, _ = generator_direction(owner, metric, direction)
                dm, _ = generator_direction(negative, metric, direction)
                diff = dm-sign@dp.conj()@sign
                errors.append(float(np.max(abs(diff.data))) if diff.nnz else 0.)
            diagnostics['signed_operator_and_jet_residual'] = max(errors)
        if name == '14_1' and (points, steps) == CONFIGS[0]:
            checks = []
            for B in range(4):
                h = 1e-5; shift = h*np.eye(4)[B]
                parts = [evolve_field_jets(owner, E, fields, np.linspace(*TIME, steps+1),
                    SuppliedKSHarmonicMetric(tuple(AMPLITUDES+sign*shift), OMEGA, TIME),
                    with_jets=False)['evolved_field'] for sign in (1., -1.)]
                checks.append(float(np.max(abs((parts[0]-parts[1])/(2*h)-result['tangent_fields'][B]))))
            diagnostics['independent_history_finite_difference'] = checks
            arrays['FD_control/evolved_field'] = result['evolved_field']
            arrays['FD_control/tangent_fields'] = result['tangent_fields']
    arrays['diagnostics_json'] = np.frombuffer(json.dumps(diagnostics, sort_keys=True).encode(), np.uint8)
    return name, arrays


def prepare(workers):
    modes, history = payload(read(MODE)), payload(read(HISTORY))
    names = sorted({k.split('/')[0] for k in modes if '/' in k and k.split('/')[0].endswith('_1')})
    jobs = []
    for name in names:
        d = family(modes, name)
        negative = name[:-1]+'-1'
        partner = family(modes, negative) if negative+'/energies' in modes else None
        jobs.append((name, d, partner, history[name+'/x'], history[name+'/reference_field']))
    arrays = {}
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for count, (name, values) in enumerate(pool.map(one, jobs), 1):
            arrays.update({name+'/'+k: v for k, v in values.items()})
            if count % 4 == 0 or count == len(jobs):
                print(f'{count}/{len(jobs)} nonlinear metric-jet families; signed partners reused', flush=True)
    arrays['metadata_json'] = np.frombuffer(json.dumps({
        'signature': {p: sha(p) for p in (*SOURCES, *INPUTS)},
        'mode_payload': read(MODE)['payload'], 'history_payload': read(HISTORY)['payload'],
        'external_coefficient_source': SBP4_COEFFICIENT_SOURCE,
    }, sort_keys=True).encode(), np.uint8)
    raw = deterministic_npz_bytes(arrays)
    digest = hashlib.sha256(raw).hexdigest()
    out = ROOT/f'results/development/artifacts/nsc-transmitting-history-jets.{digest}.npz'
    out.write_bytes(raw)
    return out


def make_record(path):
    a = load(path); metadata = json.loads(a['metadata_json'].tobytes())
    for p, expected in metadata['signature'].items():
        if sha(p) != expected:
            raise ValueError('nonlinear jet dependency changed: '+p)
    modes = payload({'payload': metadata['mode_payload']})
    rows, maxima = {}, {}
    for name in sorted({k.split('/')[0] for k in modes if '/' in k}):
        d = family(modes, name)
        positive = name.split('_')[0]+'_1'
        diagnostic = json.loads(a[positive+'/diagnostics_json'].tobytes())
        source_values, controls = {}, []
        for points, steps in CONFIGS:
            key = f'{points}_{steps}'
            A = a[positive+'/'+key+'/overlap']
            if name != positive:
                A = reflected_overlap(A)
            rows_ctp = relative_branch_ctp_control(d['energies'], d['weights'], d['source'], d['projector'], A)
            source_values[key] = [float(v['derivative'].real) for v in rows_ctp]
            controls.extend(rows_ctp)
        fine = diagnostic['801_128']
        row = {
            'group': int(d['group'].item()), 'angular_sign': int(d['sign'].item()),
            'source': 'spatial state/tangent solve' if name == positive else 'imported signed source/operator/jet map',
            'sampled_action_derivatives': source_values,
            'spatial_source_refinement': float(np.max(abs(np.array(source_values['801_64'])-source_values['401_64']))),
            'temporal_source_refinement': float(np.max(abs(np.array(source_values['801_128'])-source_values['801_64']))),
            'norm_flux': fine['norm_flux_residual'],
            'metric_flux_tangent': fine['metric_flux_tangent_residual'],
            'mode_overlap_tangent': fine['relative_overlap_tangent_residual'],
            'tangent_boundary_trace': fine['tangent_boundary_trace'],
            'CTP_tangent': max(v['tangent_residual'] for v in controls),
            'CTP_trace': max(v['trace_residual'] for v in controls),
            'imaginary_real_source': max(abs(v['derivative'].imag) for v in controls),
            'signed_source_residual': diagnostic.get('signed_source_residual', 0.),
            'signed_operator_and_jet_residual': diagnostic.get('signed_operator_and_jet_residual', 0.),
            'SBP_residual': fine['SBP_residual'], 'phase_residual': fine['phase_residual'],
            'sampled_causal_buffer': fine['sampled_causal_buffer'],
            'reference_grid_drift_in_collar': fine['reference_grid_drift_in_collar'],
        }
        if 'independent_history_finite_difference' in diagnostic:
            row['history_finite_difference_by_field'] = diagnostic['independent_history_finite_difference']
            row['history_finite_difference'] = max(row['history_finite_difference_by_field'])
        rows[name] = row
        for k in ('spatial_source_refinement', 'temporal_source_refinement', 'history_finite_difference',
                  'norm_flux', 'metric_flux_tangent', 'mode_overlap_tangent', 'tangent_boundary_trace',
                  'CTP_tangent', 'CTP_trace', 'imaginary_real_source', 'signed_source_residual',
                  'signed_operator_and_jet_residual', 'SBP_residual', 'phase_residual'):
            if k in row:
                maxima[k] = max(maxima.get(k, 0.), row[k])
    if len(rows) != 63 or {r['group'] for r in rows.values()} != set(range(33)):
        raise ValueError('all retained signed families required')
    tolerances = {k: 3e-8 if k in ('spatial_source_refinement', 'temporal_source_refinement',
                                  'history_finite_difference') else 3e-11 for k in maxima}
    failures = {k: v for k, v in maxima.items() if v > tolerances[k]}
    return {
        'schema': 'NSC-NONLINEAR-TRANSMITTING-HISTORY-JETS-v1', 'accountable_author': 'Douglas Ek',
        'status': 'PASS: nonlinear conditional field jets and sampled Gaussian variation' if not failures else 'OPEN',
        'source_hashes': {p: sha(p) for p in SOURCES}, 'input_hashes': {p: sha(p) for p in INPUTS},
        'payload': {'path': str(path.relative_to(ROOT)), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
                    'bytes': path.stat().st_size},
        'locked_inputs': read(MODE)['locked_inputs'],
        'domain': {
            'state': 'same authenticated reference mode fields and source covariance; initial preparation held fixed',
            'variation': 'four raw KS amplitudes of the supplied C-infinity response history',
            'history_coordinates': list(TIME), 'amplitudes': AMPLITUDES.tolist(), 'omega': OMEGA,
            'history_selection': 'none; response coordinates are not a physical duration or stationary solution',
            'tangent_equation': 'Ydot=L[g]Y+dL[g]Phi_g, Y(initial)=0',
            'metric_Jacobian': 'pg_log_jacobian at pg_to_ks(actual PG metric)',
            'overlap_kernel': 'A_B(Eo,Ei)=integral Phi_g,Eo dagger Y_B,Ei d rho = (U[g] dagger dU[g]) kernel',
            'Gaussian_contraction': 'relative current-history frame W=I; dWplus=A/2, dWminus=-A/2; existing ctp_first_variation',
            'measure': 'source-kernel weights sqrt(wo*wi)/(2*pi); covariance is source-fiber multiplication',
            'reported_source': 'bare Gaussian action derivatives on the inherited real-energy sample, per reduced family; not a full spectral or renormalized source',
            'new_boundary_term': False, 'closed_sampled_frequency_evolution': False,
        },
        'numerics': {
            'spatial_operator': 'imported SBP(4,2); existing NSC split/SAT expression',
            'coefficient_source': SBP4_COEFFICIENT_SOURCE,
            'configurations': [list(c) for c in CONFIGS],
            'reference_sampling': 'reused authenticated 801-point reference fields; nested spatial restriction only',
            'reference_grid_drift': 'reported separately; no numerical defect assigned as physical stress',
            'temporal_method': 'same midpoint augmented exponential for state and all tangents',
            'finite_difference_control': 'four fields on mixed family14, 401 points and64 steps; independent nonzero-history perturbations',
        },
        'per_family': rows, 'maxima': maxima, 'tolerances': tolerances, 'failures': failures,
        'gate': {
            'nonlinear_conditional_metric_jets': 'PASS' if not failures else 'OPEN',
            'sampled_Gaussian_contraction': 'PASS' if not failures else 'OPEN',
            'physical_EndpointBranchJets': 'OPEN: no physical two-endpoint variation family assigned',
            'complete_state_and_reference_integral': 'OPEN: real panels, analytic subgap and tails remain to be composed on the history',
            'B2': 'OPEN', 'V_c': None, 'extended_stationarity': 'OPEN',
            'Gamma_rest_extra_independent_term': 'none introduced; declared action inventory retained',
            'Weyl_time_node_diagnostic': 93.54264532195464,
            'absolute_stress': None, 'nulls': None, 'updated_constraints': None,
            'metric_timestep': False, 'physical_history_selected': False, 'Z3': 'OUT OF SCOPE',
            'homogeneous_nonexistence': 'imported in its original scope, not rerun',
        },
        'comparison': {'fields': 'all', 'float_atol': 3e-12, 'float_rtol': 3e-10,
                       'exact': 'source/input hashes, structure and scope labels'},
        'reproducer': 'python3 scripts/derive_nsc_transmitting_history_jets.py --check',
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--prepare', action='store_true'); p.add_argument('--check', action='store_true')
    p.add_argument('--workers', type=int, default=2)
    args = p.parse_args()
    if args.prepare == args.check:
        p.error('select --prepare or --check')
    path = prepare(args.workers) if args.prepare else ROOT/read(OUTPUT)['payload']['path']
    result = make_record(path)
    if args.check:
        from derive_nsc_transmitting_boundary_binding import compare
        compare(read(OUTPUT), result)
    else:
        (ROOT/OUTPUT).write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print(json.dumps({'status': result['status'], 'maxima': result['maxima'], 'failures': result['failures']}, indent=2), flush=True)
    if result['failures']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
