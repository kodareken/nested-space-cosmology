#!/usr/bin/env python3
"""Exact full-current angular parity with the original paired source ledger."""
import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT/'results/development/nsc-source-tail-angular-parity.json'
INVENTORY = 'results/development/nsc-ks-source-inventory.json'
JOINED = 'results/development/nsc-ks-cutoff-bridge-control.json'
BRIDGE = 'results/development/nsc-source-cutoff-bridge.json'
OWNERS = ('scripts/check_nsc_source_tail_angular_parity.py',
          'docs/nsc-source-tail-angular-parity.md')


def compute():
    a, m, ell, t = sp.symbols('a m ell t', real=True, nonzero=True)
    r = sp.Symbol('r', positive=True)
    S3, I2 = sp.diag(1, -1), sp.eye(2)
    V = sp.Matrix([[0, -m-sp.I*ell/r], [-m+sp.I*ell/r, 0]])
    identities = {}
    scalar = lambda value: value[0]
    for s in (1, -1):
        D = s*S3
        major = 0 if s == 1 else 1
        minor = 1-major
        A = [sp.eye(2)[:, major]]+[sp.Matrix(sp.symbols(f'a{s}_{j}_0:2')) for j in range(1, 5)]
        Az = [sp.zeros(2, 1)]+[sp.Matrix(sp.symbols(f'z{s}_{j}_0:2')) for j in range(1, 5)]
        Ar = sp.Matrix(sp.symbols(f'r{s}_0:2'))
        b = sp.Matrix(sp.symbols(f'b{s}_0:2'))
        bz = sp.Matrix(sp.symbols(f'c{s}_0:2'))
        L = Ar-S3*bz/a**2-sp.I*V*b/a
        transformed = D*Ar.conjugate()-S3*D*bz.conjugate()/a**2-sp.I*V.conjugate()*D*b.conjugate()/a
        identities[f's{s:+d}/generator_involution'] = [str(sp.expand(x)) for x in transformed-D*L.conjugate()]
        q = sp.Symbol('q')
        def riccati(v, potential, inv_e):
            return sp.I*potential[minor, major]/a+2*sp.I*v/(a*a*inv_e)-sp.I*potential[major, minor]*v*v/a
        identities[f's{s:+d}/initial_Riccati_involution'] = str(sp.expand(
            riccati(-sp.conjugate(q), V.conjugate(), t)+sp.conjugate(riccati(q, V, -t))))
        B = [(-1)**j*D*v.conjugate() for j, v in enumerate(A)]
        Bz = [(-1)**j*D*v.conjugate() for j, v in enumerate(Az)]
        def current(coeff, derivative, W, order):
            env = sum((scalar(coeff[j].conjugate().T*W*derivative[order-j])
                       -scalar(derivative[j].conjugate().T*W*coeff[order-j])
                       for j in range(order+1)), sp.Integer(0))/(2*sp.I)
            density = sum((scalar(coeff[j].conjugate().T*W*coeff[order+1-j])
                           for j in range(order+2)), sp.Integer(0))
            return env-s*density
        def potential(coeff, potential, order):
            return sum((scalar(coeff[j].conjugate().T*potential*coeff[order-j])
                        for j in range(order+1)), sp.Integer(0))
        for order in (1, 2, 3):
            sign = (-1)**(order+1)
            old_beta, new_beta = current(A, Az, I2, order), current(B, Bz, I2, order)
            old_N = -potential(A, V, order)-current(A, Az, S3, order)/a
            new_N = -potential(B, V.conjugate(), order)-current(B, Bz, S3, order)/a
            identities[f's{s:+d}/beta_p{order}'] = str(sp.expand(new_beta-sign*old_beta))
            identities[f's{s:+d}/N_p{order}'] = str(sp.expand(new_N-sign*old_N))
    leaves = [v for values in identities.values() for v in (values if isinstance(values, list) else [values])]
    if set(leaves) != {'0'}:
        raise AssertionError(identities)

    inventory = json.loads((ROOT/INVENTORY).read_text())
    joined = json.loads((ROOT/JOINED).read_text())
    payload = inventory['payload']
    if sha256((ROOT/payload['path']).read_bytes()).hexdigest() != payload['sha256']:
        raise ValueError('original source inventory changed')
    with np.load(ROOT/payload['path'], allow_pickle=False) as f:
        meta = json.loads(f['metadata_json'].tobytes())
        rows = {}
        for name, panel in meta['panels'].items():
            if panel['angular_magnitude'] == 0:
                continue
            key = (panel['group'], panel['angular_sign'])
            rows.setdefault(key, []).append(np.column_stack((f[name+'/energies'], f[name+'/weights'])))
    groups = sorted({key[0] for key in rows})
    checks = []
    for group in groups:
        pairs = []
        for sign in (1, -1):
            pair = np.concatenate(rows[(group, sign)])
            pairs.append(pair[np.lexsort((pair[:, 1], pair[:, 0]))])
        if not np.array_equal(pairs[0], pairs[1]):
            raise ValueError('angular partners must retain the same source labels and quadrature weights')
        ledger = [row for row in joined['ledger'] if row['group'] == group]
        if (len(ledger) != 2 or {r['positive_angular_sign'] for r in ledger} != {-1, 1}
                or ledger[0]['multiplicity_per_signed_family'] != ledger[1]['multiplicity_per_signed_family']
                or ledger[0]['absolute_angular'] != ledger[1]['absolute_angular']):
            raise ValueError('equal original angular-pair multiplicity required')
        checks.append({'group': group, 'positive_rows_per_angular_sign': len(pairs[0]),
            'source_grid_and_weight_residual': 0.,
            'multiplicity_per_signed_family': ledger[0]['multiplicity_per_signed_family']})
    if len(checks) != 30:
        raise ValueError('all30 original nonzero-angular pairs required')
    return {
        'schema': 'NSC-SOURCE-TAIL-ANGULAR-PARITY-v1', 'accountable_author': 'Douglas Ek',
        'status': 'PASS: paired raw source tail is O(Lambda^-2); numerical bound and local gate OPEN',
        'exact_identity_residuals': identities, 'paired_source_checks': checks,
        'coefficient_involution': 'A_j(-ell)=(-1)^j*(s*S3)*conjugate(A_j(ell))',
        'full_N_beta_coefficient_parity': 'at e^-p: (-1)^(p+1)',
        'leading_changed_current_cancellation': 'owned intrinsic match and delta n2_s=s*f_z_s',
        'first_allowed_raw_paired_power': 3,
        'full_inverse_cubic_coefficient_proven_nonzero': False,
        'sufficient_uniform_remainder_expansion_order': 4,
        'physical_input_energy_replaced_by_outgoing_momentum': False,
        'numerical_tail_constant': None, 'numerical_tail_tangent_constant': None,
        'exponential_coherent_source_remainder_retained': True,
        'physical_local_gate': 'OPEN', 'new_field_or_source_evolutions': 0,
        'scope': 'same smooth compact radius history for each original angular family; unchanged source and intrinsic Sigma data',
        'source_hashes': {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in OWNERS},
        'input_hashes': {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in
            (INVENTORY, payload['path'], JOINED, BRIDGE)},
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--record', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = compute()
    if args.record:
        if OUTPUT.exists():
            raise FileExistsError('angular parity record exists; use --check')
        OUTPUT.write_text(json.dumps(result, sort_keys=True, indent=2)+'\n')
    elif result != json.loads(OUTPUT.read_text()):
        raise ValueError('angular-pair source-tail replay differs')
    print(json.dumps({'status': result['status'], 'angular_pairs': len(result['paired_source_checks']),
                      'exact_residuals': result['exact_identity_residuals']}, indent=2))
