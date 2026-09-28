#!/usr/bin/env python3
"""Bind the inherited LLL geometric allocation to the existing KS history."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src')); sys.path.insert(0, str(ROOT/'tests'))
from recursive_horizons.nsc_lll_geometric_history import LLLGeometricHistoryAllocation
from recursive_horizons.nsc_spherical_local_history import SphericalHistoryGrid, OwnedCompactKSHistory
from recursive_horizons.nsc_horizon_source import conformal_stress
from test_nsc_lll_geometric_history import stationary_control

OUTPUT = 'results/development/nsc-lll-geometric-history.json'
SOURCES = ('src/recursive_horizons/nsc_lll_geometric_history.py',
           'scripts/derive_nsc_lll_geometric_history.py',
           'tests/test_nsc_lll_geometric_history.py', 'docs/nsc-lll-geometric-history.md')
INPUTS = ('results/development/nsc-spherical-local-history.json',
          'results/development/nsc-nonlinear-ks-source.json',
          'results/development/nsc-reference-band-action.json',
          'src/recursive_horizons/nsc_spherical_local_history.py',
          'src/recursive_horizons/nsc_ks_spacetime_variation.py',
          'src/recursive_horizons/nsc_lorentzian.py',
          'src/recursive_horizons/nsc_transmitting_resolvent.py',
          'src/recursive_horizons/nsc_horizon_source.py',
          'src/recursive_horizons/nsc_mode_resolved_cauchy_state.py',
          'src/recursive_horizons/nsc_charged_ctp_neck.py',
          'src/recursive_horizons/nsc_reference_band_action.py',
          'docs/nsc-horizon-source.md', 'docs/nsc-mode-resolved-cauchy-state.md',
          'docs/nsc-declared-action-scope.md', 'docs/nsc-common-ks-trace.md')
CONFIGS = ((64, 64), (96, 64), (96, 96))
AMPLITUDES = np.array([.002, .001, .003, .001])


def sha(path): return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
def read(path): return json.loads((ROOT/path).read_text())


def exact_nullmap_cancellation():
    """Apply imported stress/Schwarzian identities to the NSC normal order."""
    a, Lz, Lzz, H, Hz, p, pp, q, qp = sp.symbols('a Lz Lzz H Hz p pp q qp')
    # p=U_z, q=V_z; u_z sigma_u=P, v_z sigma_v=Q.
    P, Q = (Lz-a*H-p)/2, (Lz+a*H-q)/2
    Pz, Qz = (Lzz-a*Lz*H-a*Hz-pp)/2, (Lzz+a*Lz*H+a*Hz-qp)/2
    sigmazz = Lzz-(pp+qp)/2
    density = sigmazz-P*p-Q*q-P*P-Q*Q+(pp+qp)/2-(p*p+q*q)/4
    momentum = Qz-Q*q-Pz+P*p-Q*Q+P*P-(pp-qp)/2+(p*p-q*q)/4
    # Units in these two expressions are c/(12pi a^2).
    return {'density': str(sp.simplify(density-(Lzz-(Lz*Lz+a*a*H*H)/2))),
            'momentum': str(sp.simplify(momentum-a*Hz))}


def make_record():
    ledger = read(INPUTS[0])['locked_inputs']
    owner = LLLGeometricHistoryAllocation(ledger['magnetic_flux'])
    evaluations = {}; maxima = {'action_vs_stress_pairing': 0., 'trace_identity': 0.,
                                'inherited_stationary_allocation': 0., 'raw_metric_reconstruction': 0.}
    for nt, nx in CONFIGS:
        h = OwnedCompactKSHistory(SphericalHistoryGrid.gaussian(nt, nx))
        result = owner.evaluate(h, AMPLITUDES)
        gradient = owner.action_gradient(h, AMPLITUDES)
        q = result['geometry']; rho, flux, pressure, _ = np.moveaxis(result['stress_2D_allocation'], -1, 0)
        raw_expected = (1+AMPLITUDES[0]*h.envelope, AMPLITUDES[1]*h.envelope,
                        h.axial+AMPLITUDES[2]*h.envelope, h.radius+AMPLITUDES[3]*h.envelope)
        residuals = {
            'action_vs_stress_pairing': float(np.max(abs(gradient-result['action_gradient_from_stress']))),
            'trace_identity': float(np.max(abs(rho-pressure-abs(owner.magnetic_flux)*q['R2']/(24*np.pi)))),
            'inherited_stationary_allocation': stationary_control(h),
            'raw_metric_reconstruction': max(float(np.max(abs(q[k]-v))) for k, v in zip(('N', 'beta', 'a', 'r'), raw_expected)),
        }
        for k, v in residuals.items(): maxima[k] = max(maxima[k], v)
        baseline = owner.evaluate(h, np.zeros(4))
        evaluations[f'{nt}_{nx}'] = {
            'action': float(result['action']),
            'action_minus_fixed_baseline': float(result['action']-baseline['action']),
            'action_gradient': gradient.tolist(),
            'stress_pairing_gradient': result['action_gradient_from_stress'].tolist(),
            'action_force': (-gradient).tolist(), 'residuals': residuals,
            'allocation_2D_ranges': {k: [float(v.min()), float(v.max())] for k, v in
                zip(('delta_rho', 'delta_T01', 'delta_p_parallel'), (rho, flux, pressure))},
        }
    vector = lambda key: np.array([evaluations[key]['action'], *evaluations[key]['action_gradient'],
                                  *evaluations[key]['stress_pairing_gradient']])
    maxima['time_quadrature'] = float(np.max(abs(vector('96_64')-vector('64_64'))))
    maxima['radial_quadrature'] = float(np.max(abs(vector('96_96')-vector('96_64'))))
    # Independent fourth-order real differences of the new action on this same history.
    differences = []; step = 1e-5
    for B in range(4):
        vals = [owner.evaluate(h, AMPLITUDES+k*step*np.eye(4)[B])['action'] for k in (-2, -1, 1, 2)]
        fd = (vals[0]-8*vals[1]+8*vals[2]-vals[3])/(12*step)
        differences.append(float(abs(fd-gradient[B])))
    maxima['complex_step_vs_real_action_difference'] = max(differences)
    symbolic = exact_nullmap_cancellation()
    tol = {k: 3e-11 for k in maxima}
    failures = {k: v for k, v in maxima.items() if not np.isfinite(v) or v > tol[k]}
    failures.update({k: v for k, v in symbolic.items() if v != '0'})
    A, Ap, App = 1-3*np.pi/2, 6., -3*np.pi
    uu, uv, vv = conformal_stress(A, Ap, App, 0., 0., central_charge=abs(owner.magnetic_flux))
    neck = [(uu-2*uv+vv)/(-A), (-uu+vv)/(-A), (uu+2*uv+vv)/(-A), 0.]
    return {
        'schema': 'NSC-LLL-GEOMETRIC-HISTORY-v1', 'accountable_author': 'Douglas Ek',
        'status': 'PASS: inherited LLL geometric allocation on the supplied KS history' if not failures else 'OPEN',
        'source_hashes': {p: sha(p) for p in SOURCES}, 'input_hashes': {p: sha(p) for p in INPUTS},
        'locked_inputs': ledger, 'evaluations': evaluations,
        'maxima': maxima, 'tolerances': tol, 'failures': failures,
        'exact_nullmap_cancellation': symbolic,
        'inherited_neck_allocation_2D': dict(zip(('delta_rho', 'delta_T01', 'delta_p_parallel', 'delta_p_perp'), neck)),
        'real_action_difference': {'step': step, 'stencil_order': 4, 'residuals': differences},
        'domain': {
            'metric': 'physical radial g2; canonical half-density chi=sqrt(a*r^2)*psi',
            'history': 'existing supplied compact harmonic raw-KS metric, fixed PG chart',
            'amplitudes': AMPLITUDES.tolist(), 'omega': .4, 'PG_time_support': [0., .05], 'radial_support': [-1., 1.],
            'central_charge': abs(owner.magnetic_flux), 'multiplicity': 'one factor abs(q_mag); no further copy/partner factor',
            'fields': ['N_K', 'beta_K', 'a_K=q_ADM', 'r_K'],
            'stress_order': ['rho', 'covariant_orthonormal_T01', 'p_parallel', 'p_perp'],
            'trace_measure': 'dT dz = dPGtime drho/a0; fixed chart, not 1/a_actual',
            'normal_order': 'spatial canonical LLL chiral vacuum, with Schwarzian state transformation',
            'allocation': 'T_covariant - T_canonical; state functions and null maps cancel',
            'action': 'c/(24*pi) integral(-N*a*Ha^2 + 2*N_z*L_z/a - N*L_z^2/a) dT dz',
            'force_convention': 'minus action gradient',
            'CTP': 'same allocation S[g_plus]-S[g_minus]; relative-branch variation',
            'sphere_projection': 'radial allocation divided by 4*pi*r^2; zero sphere pressure in this sector',
            'boundary': 'compact boundary-identical variations; smooth fixed rho=0 cut, no new endpoint force',
            'double_count_guard': 'replaces stationary geometric LLL evaluation; do not add its old tensor again',
            'baseline': 'fixed baseline used only to report action differences, not a vacuum subtraction',
            'history_selected': False, 'state_changed': False},
        'gate': {'LLL_geometric_allocation': 'PASS' if not failures else 'OPEN',
            'common_subtracted_spectral_trace': 'OPEN: state and ad4/band connection terms require common-kernel pairing and convergence',
            'full_renormalized_stress': None, 'nulls': None, 'updated_constraints': None,
            'B2': 'OPEN', 'V_c': None, 'extended_stationarity': 'OPEN', 'metric_timestep': False,
            'Z3': 'OUT OF SCOPE', 'extra_Gamma_rest': 'none introduced',
            'Weyl_time_node_diagnostic': 93.54264532195464,
            'homogeneous_nonexistence': 'preserved in original scope; not rerun', 'PDF_bumped': False},
        'comparison': {'fields': 'all', 'float_atol': 3e-12, 'float_rtol': 3e-11, 'exact': 'hashes, structure, symbolic identities and scope'},
        'reproducer': 'python3 scripts/derive_nsc_lll_geometric_history.py --check',
    }


def compare(expected, actual, path='$'):
    if type(expected) is not type(actual): raise AssertionError('type mismatch '+path)
    if isinstance(expected, dict):
        if expected.keys() != actual.keys(): raise AssertionError('keys '+path)
        for k in expected: compare(expected[k], actual[k], path+'.'+k)
    elif isinstance(expected, list):
        if len(expected) != len(actual): raise AssertionError('length '+path)
        for i, (a, b) in enumerate(zip(expected, actual)): compare(a, b, path+f'[{i}]')
    elif isinstance(expected, float):
        if not np.isclose(expected, actual, atol=3e-12, rtol=3e-11): raise AssertionError('numeric mismatch '+path)
    elif expected != actual: raise AssertionError('mismatch '+path)


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--check', action='store_true'); args = parser.parse_args()
    result = make_record()
    if args.check: compare(read(OUTPUT), result)
    else: (ROOT/OUTPUT).write_text(json.dumps(result, sort_keys=True, indent=2)+'\n')
    print(json.dumps({k: result[k] for k in ('status', 'maxima', 'failures')}, indent=2), flush=True)
    if result['failures']: raise SystemExit(1)


if __name__ == '__main__': main()
