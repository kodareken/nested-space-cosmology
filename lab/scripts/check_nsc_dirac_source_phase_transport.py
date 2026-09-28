#!/usr/bin/env python3
"""Exact leading density transport in the owned source-frequency Dirac PDE."""
import argparse
from hashlib import sha256
import json
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT/'results/development/nsc-dirac-source-phase-transport.json'
OWNERS = ('scripts/check_nsc_dirac_source_phase_transport.py',
          'docs/nsc-dirac-source-phase-transport.md',
          'src/recursive_horizons/nsc_ks_signed_state.py',
          'src/recursive_horizons/nsc_ks_difference_envelope.py')


def compute():
    rho, z = sp.symbols('rho z', real=True)
    e = sp.Symbol('e', positive=True)
    m, ell = sp.symbols('m ell', real=True)
    a = sp.Function('a', real=True)(rho)
    r = sp.Function('r', real=True)(rho, z)
    phase = sp.Function('F', real=True)(rho, z)
    S3 = sp.diag(1, -1)
    V = sp.Matrix([[0, -m-sp.I*ell/r], [-m+sp.I*ell/r, 0]])
    Q = m*m+ell*ell/(r*r)
    exact = {}
    for s in (1, -1):
        major = 0 if s == 1 else 1
        minor = 1-major
        G = sp.I*V/a-sp.I*s*e*S3/a**2
        skew = G+G.conjugate().T
        exact[f's{s:+d}/antihermitian_generator'] = [str(sp.simplify(v)) for v in skew]
        q1 = -a*V[minor, major]/2
        minor_equation = 2*sp.I*q1/a**2+sp.I*V[minor, major]/a
        exact[f's{s:+d}/minor_coefficient'] = str(sp.simplify(minor_equation))
        # D_{-s} q1 enters the next minor coefficient. Major b1=i*F
        # obeys D_s F=-Q/2; this relation is used only for its derivative.
        q2 = a*a/(2*sp.I)*(sp.diff(q1, rho)+s*sp.diff(q1, z)/a**2
                              -sp.I*V[minor, major]*(sp.I*phase)/a)
        major_second_rhs = sp.I*V[major, minor]*q2/a
        minor_density = a*a*Q/4
        density_transport = (2*sp.re(sp.expand_complex(major_second_rhs))
            -phase*Q+sp.diff(minor_density, rho)-s*sp.diff(minor_density, z)/a**2)
        expression = density_transport+s*sp.diff(Q, z)/2
        # Derivatives of the declared real smooth geometry are real. SymPy
        # does not inherit that assumption from a real-valued Function.
        real_jets = {jet: sp.Symbol('real_jet_'+str(i), real=True)
                     for i, jet in enumerate(sorted(expression.atoms(sp.Derivative), key=str))}
        residual = sp.simplify(sp.expand_complex(expression.xreplace(real_jets)))
        exact[f's{s:+d}/density_transport'] = str(residual)
        fp = sp.Symbol('f_z', real=True)
        exact[f's{s:+d}/coincidence_inverse_energy_current'] = str(sp.simplify(fp-s*(s*fp)))
    leaves = [entry for value in exact.values() for entry in (value if isinstance(value, list) else [value])]
    if set(leaves) != {'0'}:
        raise AssertionError(exact)
    return {
        'schema': 'NSC-DIRAC-SOURCE-PHASE-TRANSPORT-v1', 'accountable_author': 'Douglas Ek',
        'status': 'PASS: exact formal source-column coefficient; uniform remainder OPEN',
        'exact_residuals': exact,
        'source_label': 'E=s*e, e>0; principal occupied S3=s; carrier exp(-i*s*e*z)',
        'exact_continuity': 'partial_rho n = partial_z(X† S3 X)/a²',
        'minor_coefficient': 'A1[-s]=-a*V[-s,s]/2',
        'density_transport': '(partial_rho-s/a² partial_z) n2 = -s/2 partial_z(V²)',
        'relative_phase': 'f_s=ell²/2 integral_1^rho_up (r_g^-2-r_ref^-2)(rho,z-s*D(rho)) d rho',
        'relative_density_coefficient': 'delta n2=s*partial_z f_s',
        'source_scope': 'occupied affine-vacuum asymptotics plus owned exponential coherent source remainder; no replacement of C_up',
        'scope': {'formal_coefficient_only': True, 'uniform_remainder_bound': None,
            'finite_source_cutoff_error': None, 'stress_term_inserted': False,
            'physical_local_gate': 'OPEN', 'new_field_or_source_runs': 0, 'metric_timestep': False},
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
            raise FileExistsError('source phase identity exists; use --check')
        OUTPUT.write_text(json.dumps(result, sort_keys=True, indent=2)+'\n')
    elif result != json.loads(OUTPUT.read_text()):
        raise ValueError('source phase identity replay differs')
    print(json.dumps({k: result[k] for k in ('status', 'exact_residuals')}, indent=2))
