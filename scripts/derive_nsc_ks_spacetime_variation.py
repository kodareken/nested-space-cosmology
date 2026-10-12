#!/usr/bin/env python3
"""Check the new KS spacetime/CTP chain without rerunning mode generators."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import block_diag

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from recursive_horizons.nsc_ks_spacetime_variation import harmonic_mode_vertices, harmonic_dyson_derivative
from recursive_horizons.nsc_transmitting_ctp_variation import ctp_first_variation
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes

OUTPUT = 'results/development/nsc-ks-spacetime-variation.json'
MODE_RECORD = 'results/development/nsc-pg-ctp-mode-jets.json'
KS_RECORD = 'results/development/nsc-pg-ks-metric-pullback.json'
SOURCES = ('src/recursive_horizons/nsc_ks_spacetime_variation.py',
           'scripts/derive_nsc_ks_spacetime_variation.py',
           'tests/test_nsc_ks_spacetime_variation.py',
           'docs/nsc-ks-spacetime-variation.md')
INPUTS = (MODE_RECORD, KS_RECORD,
          'src/recursive_horizons/nsc_pg_ctp_mode_jets.py',
          'src/recursive_horizons/nsc_pg_ks_metric_pullback.py',
          'src/recursive_horizons/nsc_transmitting_resolvent.py',
          'src/recursive_horizons/nsc_common_time_bulk_split.py',
          'src/recursive_horizons/nsc_transmitting_ctp_variation.py',
          'src/recursive_horizons/nsc_lorentzian.py',
          'docs/nsc-adm-source-constraints.md',
          'docs/nsc-declared-action-scope.md',
          'results/development/nsc-smooth-seam-variation.json')
OMEGA = .4  # Existing response-transfer label, not a model parameter.


def digest(path):
    return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()


def read(path):
    return json.loads((ROOT/path).read_text())


def load(path):
    with np.load(path, allow_pickle=False) as a:
        return {k: a[k].copy() for k in a.files}


def verified_payload(record):
    p = record['payload']
    if digest(p['path']) != p['sha256']:
        raise ValueError('authenticated input changed: '+p['path'])
    return load(ROOT/p['path'])


def family_inputs(mode):
    for name in sorted({k.split('/')[0] for k in mode if '/' in k}):
        yield name, {k[len(name)+1:]: v for k, v in mode.items() if k.startswith(name+'/')}


def one_family(job):
    name, data = job
    fields = {k[7:]: v for k, v in data.items() if k.startswith('fields/')}
    mass, angular = float(data['mass'].item()), float(data['angular'].item())
    M, strong, finite = harmonic_mode_vertices(fields, mass, angular, OMEGA)
    negative, _, _ = harmonic_mode_vertices(fields, mass, angular, -OMEGA, difference_step=None)
    static, _, _ = harmonic_mode_vertices(fields, mass, angular, 0., difference_step=None)
    return name, {'vertex': M, 'strong': strong, 'finite': finite,
                  'negative': negative, 'static': static}


def prepare(workers):
    mode_record = read(MODE_RECORD)
    mode = verified_payload(mode_record)
    arrays = {}
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for count, (name, values) in enumerate(pool.map(one_family, family_inputs(mode)), 1):
            arrays.update({name+'/'+k: v for k, v in values.items()})
            if count % 8 == 0 or count == 63:
                print(f'{count}/63 new spacetime vertices; archived modes reused', flush=True)
    signature = {p: digest(p) for p in (*SOURCES, *INPUTS)}
    arrays['metadata_json'] = np.frombuffer(json.dumps({
        'signature': signature, 'mode_payload': mode_record['payload'],
        'KS_payload': read(KS_RECORD)['payload'], 'omega': OMEGA,
    }, sort_keys=True).encode(), np.uint8)
    raw = deterministic_npz_bytes(arrays)
    sha = hashlib.sha256(raw).hexdigest()
    path = ROOT/f'results/development/artifacts/nsc-ks-spacetime-variation.{sha}.npz'
    path.write_bytes(raw)
    return path


def evaluate(arrays, mode, archived_KS):
    rows, maxima = {}, {}
    for name, d in family_inputs(mode):
        M = arrays[name+'/vertex']
        P = block_diag(*d['projector'])
        active = np.flatnonzero(np.diag(P).real > .5)
        C = block_diag(*d['source'])[np.ix_(active, active)]
        E = np.repeat(d['energies'], 3)[active]
        rootw = np.repeat(np.sqrt(d['weights']/(2*np.pi)), 3)
        vertex = (M*rootw[None, :, None]*rootw[None, None, :])[:, active][:, :, active]
        # Conditional response window only. This is NOT a metric history.
        U, dU, generator = harmonic_dyson_derivative(E, vertex, OMEGA, .2, .6)
        controls = [ctp_first_variation(C, U, U, .5*du, -.5*du) for du in dU]
        g, sign = map(int, name.split('_'))
        row = {
            'group': g, 'angular_sign': sign,
            'weak_vs_strong': float(np.max(abs(M-arrays[name+'/strong']))),
            'weak_vs_perturbed_H': float(np.max(abs(M-arrays[name+'/finite']))),
            'conjugate_frequency': float(np.max(abs(arrays[name+'/negative']-M.swapaxes(-1, -2).conj()))),
            'zero_transfer_vs_archived_KS': float(np.max(abs(arrays[name+'/static']-archived_KS[name+'/raw_KS_vertex']))),
            'CTP_tangent': max(x['tangent_residual'] for x in controls),
            'CTP_unitarity': max(x['unitarity_residual'] for x in controls),
            'CTP_trace': max(float(abs(x['derivative']+np.trace(C@A))) for x, A in zip(controls, generator)),
            'imaginary_real_source': max(abs(x['derivative'].imag) for x in controls),
            'seam': float(d['fields/seam_residual'].item()),
            'omitted_spatial_phase_error': float(np.linalg.norm(M-arrays[name+'/static'])),
            'conditional_source_differential': [float(x['derivative'].real) for x in controls],
            'physical_history_selected': False,
        }
        rows[name] = row
        for key in ('weak_vs_strong', 'weak_vs_perturbed_H', 'conjugate_frequency',
                    'zero_transfer_vs_archived_KS', 'CTP_tangent', 'CTP_unitarity',
                    'CTP_trace', 'imaginary_real_source', 'seam'):
            maxima[key] = max(maxima.get(key, 0.), row[key])
    if len(rows) != 63 or {x['group'] for x in rows.values()} != set(range(33)):
        raise ValueError('all retained groups and signed families required')
    return rows, maxima


def make_record(path):
    arrays = load(path)
    meta = json.loads(arrays['metadata_json'].tobytes())
    for p, sha in meta['signature'].items():
        if digest(p) != sha:
            raise ValueError('new variation dependency changed: '+p)
    mode = verified_payload({'payload': meta['mode_payload']})
    KS = verified_payload({'payload': meta['KS_payload']})
    rows, maxima = evaluate(arrays, mode, KS)
    tolerances = {k: (3e-8 if k in ('weak_vs_strong', 'weak_vs_perturbed_H') else 3e-11) for k in maxima}
    failures = {k: v for k, v in maxima.items() if v > tolerances[k]}
    return {
        'schema': 'NSC-KS-SPACETIME-VARIATION-v1', 'accountable_author': 'Douglas Ek',
        'status': 'PASS: reference KS-spacetime directional CTP derivative' if not failures else 'FAIL',
        'source_hashes': {p: digest(p) for p in SOURCES},
        'input_hashes': {p: digest(p) for p in INPUTS},
        'payload': {'path': str(path.relative_to(ROOT)), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
                    'bytes': path.stat().st_size},
        'locked_inputs': read(MODE_RECORD)['locked_inputs'],
        'domain': {
            'operator': 'same transmitting whole-line canonical PG Dirac operator',
            'chart': 'T=T(rho), z=tau+S(rho); reference chart held fixed',
            'Cauchy_representation': 'fixed PG canonical half-density and physical source C0',
            'metric_direction': 'raw KS xi(T,z); evaluated f(T(rho))=s(rho), Fourier in z',
            'radial_chain_rule': 'd_rho(J xi)=J_prime xi+J(-xi_T/a0+beta0 xi_z/a0^2)',
            'canonical_measure': 'already included in owned Hamiltonian; no second half-density force',
            'source_measure': 'dE/(2*pi); finite samples are kernel quadrature controls, not closed bulk',
            'evaluated_frequencies': 'same two positive and paired negative derivative-control nodes per signed family',
            'omega_response_label': OMEGA, 'conditional_time_arguments': [.2, .6],
            'time_argument_role': 'algebraic Duhamel/CTP check only; no selected duration or metric history',
            'spatial_seam': 'fixed rho=0; continuous fields and common smooth metric variations',
        },
        'per_family': rows, 'maxima': maxima, 'tolerances': tolerances, 'failures': failures,
        'jet_inventory': {
            'raw_KS_spacetime_metric_jet': 'computed including J_prime, T_prime and S_prime',
            'reference_unitary_directional_derivative': 'computed for real harmonic metric direction',
            'Gaussian_CTP_contraction': 'computed on authenticated source fibers',
            'physical_KS_endpoint_variation_family': None,
            'full_transmitting_EndpointBranchJets': None,
            'nonlinear_history_U_g': None,
        },
        'gate': {
            'B1_spacetime_direction': 'PASS' if not failures else 'FAIL',
            'B1_full': 'OPEN: physical endpoint/history family and nonlinear U[g] not evaluated',
            'B2': 'OPEN: no full renormalized metric variation assembled',
            'Gamma_rest_extra_independent_term': 'not declared; none introduced',
            'full_Gaussian_Schur_count': 1,
            'Weyl_time_node_diagnostic': 93.54264532195464,
            'physical_interface_mismatch': None, 'V_c': None, 'extended_stationarity': 'OPEN',
            'absolute_stress': None, 'new_nulls': None, 'updated_constraints': None,
            'history_selected': False, 'metric_timestep': False, 'Z3': 'OUT OF SCOPE',
            'homogeneous_nonexistence': 'imported unchanged, not re-evaluated or extended',
        },
        'comparison': {'fields': 'all', 'float_atol': 3e-12, 'float_rtol': 3e-10,
                       'exact': 'source/input hashes, structure and scope labels'},
        'reproducer': 'python3 scripts/derive_nsc_ks_spacetime_variation.py --check',
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--prepare', action='store_true')
    p.add_argument('--workers', type=int, default=2)
    p.add_argument('--check', action='store_true')
    args = p.parse_args()
    if args.prepare == args.check:
        p.error('select exactly one of --prepare or --check')
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
