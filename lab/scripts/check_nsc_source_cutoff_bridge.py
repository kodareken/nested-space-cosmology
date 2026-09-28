#!/usr/bin/env python3
"""Exact algebra behind the fixed-source differentiated-kernel bridge."""
import argparse
from hashlib import sha256
import json
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT/'results/development/nsc-source-cutoff-bridge.json'
OWNERS = ('scripts/check_nsc_source_cutoff_bridge.py', 'docs/nsc-source-cutoff-bridge.md',
          'src/recursive_horizons/nsc_ks_residual_error.py',
          'src/recursive_horizons/nsc_evolved_incoming_constraints.py',
          'src/recursive_horizons/nsc_common_subtracted_ks_source.py')
INPUTS = ('results/development/nsc-dirac-source-phase-transport.json',
          'results/development/nsc-ks-local-embedding.json',
          'docs/nsc-incoming-fixed-transfer.md', 'docs/nsc-reference-band-bulk.md')


def compute():
    e, K, a, mu = sp.symbols('e K a mu', positive=True)
    eta = sp.Symbol('eta', real=True)
    f = sp.Symbol('f_z', real=True)
    S3 = sp.diag(1, -1)
    identities = {}
    for s in (1, -1):
        occupied = (sp.eye(2)+s*S3)/2
        other = sp.eye(2)-occupied
        for order in (3, 4):
            coeff = [sp.Matrix(sp.symbols(f'A{j}_0:2')) for j in range(order+1)]
            coeff[0] = occupied*coeff[0]
            last = other*sp.Matrix(sp.symbols('R_0:2'))
            L = [2*sp.I/a**2*other*coeff[j+1] for j in range(order)]+[last]
            defect = sum((L[j]/e**j for j in range(order+1)), sp.zeros(2, 1))
            defect -= 2*sp.I*e/a**2*other*sum((coeff[j]/e**j for j in range(order+1)), sp.zeros(2, 1))
            identities[f's{s:+d}/finite_recurrence_M{order}'] = [str(sp.simplify(v)) for v in defect-last/e**order]
        I1 = sp.Function('I1')(eta)
        I2 = sp.exp(-sp.I*s*K*eta)/K-sp.I*s*eta*I1
        paired = sp.simplify(sp.I*eta*f*I1+s*f*I2)
        edge = sp.simplify(-sp.I*sp.diff(paired, eta).subs(eta, 0))
        integrand = sp.exp(-sp.I*s*e*eta)*(s*f/e**2+sp.I*eta*f/e)
        direct = sp.simplify((-sp.I*sp.diff(integrand, eta)).subs(eta, 0))
        identities[f's{s:+d}/kernel_primitive'] = str(sp.simplify(paired-s*f*sp.exp(-sp.I*s*K*eta)/K))
        identities[f's{s:+d}/ordered_current_edge'] = str(sp.simplify(edge+f))
        identities[f's{s:+d}/coincidence_first_current'] = str(direct)
        J = edge*occupied/(2*sp.pi)
        N, beta = -mu*sp.trace(S3*J)/a, mu*sp.trace(J)
        identities[f's{s:+d}/N_vertex_sign'] = str(sp.simplify(N-mu*s*f/(2*sp.pi*a)))
        identities[f's{s:+d}/beta_vertex_sign'] = str(sp.simplify(beta+mu*f/(2*sp.pi)))
    for order in (3, 4):
        tail = sp.integrate(e**(1-order), (e, K, sp.oo))
        identities[f'integrable_current_error_M{order}'] = str(sp.simplify(tail-K**(2-order)/(order-2)))
    leaves = [x for value in identities.values() for x in (value if isinstance(value, list) else [value])]
    if set(leaves) != {'0'}:
        raise AssertionError(identities)
    return {
        'schema': 'NSC-SOURCE-CUTOFF-BRIDGE-v1', 'accountable_author': 'Douglas Ek',
        'status': 'PASS: scoped analytic source-cutoff identity; finite-cutoff error and local gate OPEN',
        'exact_algebra_residuals': identities,
        'state_law': 'C_Sigma[g]=U_g C_up U_g†; C_up unchanged',
        'history_class': 'smooth compact pure-radius histories, a=a_ref>0, r_g>0, intrinsic Sigma fields fixed, reference near rho_up, owned local embedding',
        'source_scope': 'original coherent source; affine-vacuum initial expansion with exponential coherent remainder',
        'uniform_remainder': {
            'finite_expansion_order': 'M>=3',
            'envelope_H2_error': 'C_M/e^M, finite uniform C_M from initial expansion and smooth compact coefficient transport',
            'differentiated_kernel_error': 'O(e^(1-M)), integrable',
            'numerical_C_M': None},
        'bridge': {
            'N_edge': 'mu*(f_z_plus-f_z_minus)/(2*pi*a_Sigma)',
            'beta_edge': '-mu*(f_z_plus+f_z_minus)/(2*pi)',
            'raw_plus_edge_plus_tail': True,
            'omitted_raw_tail_order': 'O(1/Lambda)',
            'finite_panel_tail_error_bound': None,
            'outgoing_Weyl_momentum_identified_with_source_E': False,
            'kernel_derivative': 'off-diagonal distributional derivative; equivalently Abel-regularized oscillatory integral'},
        'scope': {'physical_local_gate': 'OPEN', 'physical_root_selected': False,
            'constraint_assembly_changed_by_this_record': False, 'new_action_term': False,
            'Gamma_rest_assigned': False, 'band_bulk_reopened': False,
            'local_reference_coefficients_changed': False,
            'new_field_or_source_evolutions': 0, 'metric_timestep': False},
        'source_hashes': {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in OWNERS},
        'input_hashes': {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in INPUTS},
    }


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    mode = p.add_mutually_exclusive_group(required=True)
    mode.add_argument('--record', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = p.parse_args()
    result = compute()
    if args.record:
        if OUTPUT.exists():
            raise FileExistsError('source-cutoff bridge already recorded; use --check')
        OUTPUT.write_text(json.dumps(result, sort_keys=True, indent=2)+'\n')
    elif result != json.loads(OUTPUT.read_text()):
        raise ValueError('source-cutoff bridge replay differs')
    print(json.dumps({'status': result['status'], 'exact_residuals': result['exact_algebra_residuals']}, indent=2))
