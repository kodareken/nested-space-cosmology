#!/usr/bin/env python3
"""Exact unitary cutoff control; this supplies no NSC stress correction.

Reuse the common-kernel owner's independent source cutoff and point split.
The question is whether per-column phase reasoning can exclude a finite
cutoff/coincidence exchange. A bounded Hermitian Jacobi operator gives a
globally unitary counterexample, checked without any physical source run.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / 'results/development/nsc-source-cutoff-unitary-control.json'
OWNERS = ('scripts/check_nsc_source_cutoff_unitary_control.py',
          'docs/nsc-source-cutoff-unitary-control.md')


def compute():
    N = sp.Symbol('N', integer=True, positive=True)
    z, eta = sp.symbols('z eta', real=True)
    kbar = -N - sp.Rational(1, 2)
    off = sp.I / (2 * (N + sp.Rational(1, 2)))
    # e_n(z)=exp(-i*n*z)/sqrt(2*pi). Only the two cutoff-edge entries
    # of i[B,P_N] are nonzero. B is real symmetric, not a per-fibre phase.
    kernel = sp.simplify((off * sp.exp(-sp.I*z)
                         - off * sp.exp(sp.I*z))
                        * sp.exp(sp.I*kbar*eta) / (2*sp.pi))
    expected = sp.sin(z) * sp.exp(sp.I*kbar*eta) / (2*sp.pi*(N + sp.Rational(1, 2)))
    current = sp.simplify(-sp.I * sp.diff(kernel, eta))
    expected_current = -sp.sin(z) * sp.exp(sp.I*kbar*eta) / (2*sp.pi)
    exact = {
        'kernel_identity': str(sp.simplify(sp.expand_complex(kernel/expected)-1)),
        'current_identity': str(sp.simplify(sp.expand((current-expected_current).rewrite(sp.exp)))),
        'zero_split_current': str(sp.simplify(sp.expand((current.subs(eta, 0) + sp.sin(z)/(2*sp.pi)).rewrite(sp.exp)))),
    }
    if set(exact.values()) != {'0'}:
        raise AssertionError(exact)

    # Independent finite-matrix exponentiation. Leave several empty modes
    # beyond the cutoff; no physical history or source is used here.
    cutoff, size, step = 8, 13, 1e-4
    B = np.zeros((size, size), complex)
    for n in range(1, size):
        B[n-1, n] = B[n, n-1] = 1/(2*(n+.5))
    P = np.diag(np.arange(1, size+1) <= cutoff).astype(complex)
    derivative = 1j*(B@P-P@B)
    plus, minus = expm(1j*step*B), expm(-1j*step*B)
    difference = (plus@P@plus.conj().T-minus@P@minus.conj().T)/(2*step)
    residual = float(np.max(np.abs(difference-derivative)))
    unitary = float(np.max(np.abs(plus@plus.conj().T-np.eye(size))))
    point, split = .7, .13
    labels = -np.arange(1, size+1)
    left = np.exp(1j*labels*(point+split/2))/np.sqrt(2*np.pi)
    right = np.exp(1j*labels*(point-split/2))/np.sqrt(2*np.pi)
    means = (labels[:, None]+labels[None, :])/2
    matrix_current = np.einsum('p,pq,pq,q->', left, derivative, means, right.conj())
    formula = -np.sin(point)*np.exp(-1j*(cutoff+.5)*split)/(2*np.pi)
    contraction = float(abs(matrix_current-formula))
    if residual > 1e-10 or unitary > 1e-13 or contraction > 1e-13:
        raise AssertionError((residual, unitary, contraction))

    return {
        'schema': 'NSC-SOURCE-CUTOFF-UNITARY-CONTROL-v1',
        'accountable_author': 'Douglas Ek',
        'status': 'PASS: exact globally unitary control; NSC exchange remains OPEN',
        'exact_identity_residuals': exact,
        'generator': 'B[n+1,n]=B[n,n+1]=1/(2*(n+1/2)), n>=1; negative Fourier sector',
        'unitary': 'U(t)=exp(i*t*B); B bounded self-adjoint; U commutes with the full negative-sector projector',
        'finite_cutoff': 'P_N projects onto n=1,...,N; derivative at t=0 is i[B,P_N]',
        'kernel_derivative': 'sin(z)*exp(-i*(N+1/2)*eta)/(2*pi*(N+1/2))',
        'current_derivative': '-sin(z)*exp(-i*(N+1/2)*eta)/(2*pi)',
        'kernel_uniform_bound': '1/(2*pi*(N+1/2)) -> 0',
        'zero_split_then_N_limit': '-sin(z)/(2*pi)',
        'full_projector_current_change': '0 exactly',
        'Abel_current_factor': '(1-r)*exp(-3*i*eta/2)/(1-r*exp(-i*eta)), 0<r<1',
        'Abel_then_coincidence_limit': '0; take r->1 at fixed eta not in 2*pi*Z, then eta->0',
        'numerical_control': {'cutoff': cutoff, 'matrix_size': size, 'tangent_step': step,
            'tangent_residual': residual, 'tangent_tolerance': 1e-10,
            'unitarity_residual': unitary, 'matrix_current_residual': contraction,
            'matrix_tolerance': 1e-13},
        'scope': {'NSC_Dirac_evolution_used': False, 'new_physical_source_runs': 0,
            'NSC_cutoff_exchange_determined': False, 'stress_correction_inserted': False,
            'physical_local_gate': 'OPEN', 'metric_timestep': False},
        'source_hashes': {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in OWNERS},
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
            raise FileExistsError('unitary control already recorded; use --check')
        OUTPUT.write_text(json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+'\n')
    elif result != json.loads(OUTPUT.read_text()):
        raise ValueError('unitary control replay differs')
    print(json.dumps({k: result[k] for k in ('status', 'exact_identity_residuals', 'numerical_control')}, indent=2))
