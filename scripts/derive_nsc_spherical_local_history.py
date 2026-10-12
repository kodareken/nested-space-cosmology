#!/usr/bin/env python3
"""Evaluate only the new spatially varying local-history contribution."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src')); sys.path.insert(0, str(ROOT/'tests'))
from recursive_horizons.nsc_spherical_local_history import (
    SphericalHistoryGrid, LockedSphericalLocalAction, OwnedCompactKSHistory,
)
from recursive_horizons.nsc_transmitting_history_modes import SuppliedKSHarmonicMetric
from test_nsc_spherical_local_history import limit_control

OUTPUT = 'results/development/nsc-spherical-local-history.json'
SOURCES = ('src/recursive_horizons/nsc_spherical_local_history.py',
           'scripts/derive_nsc_spherical_local_history.py', 'tests/test_nsc_spherical_local_history.py',
           'docs/nsc-spherical-local-history.md')
INPUTS = ('results/development/nsc-subgap-history-response.json',
          'src/recursive_horizons/nsc_general_ks_local_history.py',
          'src/recursive_horizons/nsc_finite_terms.py',
          'src/recursive_horizons/nsc_pg_ks_metric_pullback.py',
          'src/recursive_horizons/nsc_ks_spacetime_variation.py',
          'src/recursive_horizons/nsc_transmitting_history_modes.py',
          'src/recursive_horizons/nsc_lorentzian.py',
          'src/recursive_horizons/nsc_transmitting_resolvent.py',
          'docs/nsc-declared-action-scope.md', 'docs/nsc-smooth-seam-variation.md')
CONFIGS = ((128, 128), (192, 128), (192, 160))
AMPLITUDES = np.array([.002, .001, .003, .001])


def sha(p): return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def read(p): return json.loads((ROOT/p).read_text())


def actual_metric(t, x, amplitudes):
    N, b, q, r = SuppliedKSHarmonicMetric(tuple(amplitudes), .4, (0., .05)).values(t, np.array([x]))[:, 0]
    return np.array([[N*N-q*q*b*b, -q*q*b], [-q*q*b, -q*q]]), r


def independent_coordinate_jets():
    t, x = .019, .23
    # Only values from the pre-existing metric owner enter the FD side.
    grid = SphericalHistoryGrid(np.array([t]), np.array([x]), np.array([.05]), np.array([2.]), np.zeros((1, 1)), np.zeros((1, 1)))
    g, r, (dg, ddg, dr, ddr) = OwnedCompactKSHistory(grid).fields_and_jets(AMPLITUDES)
    # Second differences of order-one metric values amplify roundoff as h^-2.
    # The smaller 2e-5 time step gave 1.41e-7 relative cancellation error;
    # this wider fourth-order stencil resolves the same fixed derivative.
    steps = (8e-5, 1e-3)
    errors = {'value': max(float(np.max(abs(g[0, 0]-actual_metric(t, x, AMPLITUDES)[0]))),
                           float(abs(r[0, 0]-actual_metric(t, x, AMPLITUDES)[1]))),
              'first_relative': 0., 'second_relative': 0.}
    def fields(tt, xx):
        G, R = actual_metric(tt, xx, AMPLITUDES)
        return np.r_[G.ravel(), R]
    initial = fields(t, x)
    for axis, h in enumerate(steps):
        def f(k): return fields(t+k*h if axis == 0 else t, x+k*h if axis == 1 else x)
        fd1 = (f(-2)-8*f(-1)+8*f(1)-f(2))/(12*h)
        fd2 = (-f(2)+16*f(1)-30*initial+16*f(-1)-f(-2))/(12*h*h)
        target1 = np.r_[dg[0, 0, axis].ravel(), dr[0, 0, axis]]
        target2 = np.r_[ddg[0, 0, axis, axis].ravel(), ddr[0, 0, axis, axis]]
        errors['first_relative'] = max(errors['first_relative'], float(np.max(abs(fd1-target1))/max(1., np.max(abs(target1)))))
        errors['second_relative'] = max(errors['second_relative'], float(np.max(abs(fd2-target2))/max(1., np.max(abs(target2)))))
    coefficients = {-2: 1., -1: -8., 1: 8., 2: -1.}
    mixed = sum(wi*wj*fields(t+i*steps[0], x+j*steps[1]) for i, wi in coefficients.items() for j, wj in coefficients.items())/(144*np.prod(steps))
    target = np.r_[ddg[0, 0, 0, 1].ravel(), ddr[0, 0, 0, 1]]
    errors['mixed_relative'] = float(np.max(abs(mixed-target))/max(1., np.max(abs(target))))
    return errors


def make_record():
    ledger = read(INPUTS[0])['locked_inputs']; action = LockedSphericalLocalAction(ledger)
    rows = {}; maximum_identity = 0.
    for nt, nx in CONFIGS:
        grid = SphericalHistoryGrid.gaussian(nt, nx); history = OwnedCompactKSHistory(grid)
        values, curv = action.actions(grid, *history.fields_and_jets(AMPLITUDES))
        baseline, _ = action.actions(grid, *history.fields_and_jets(np.zeros(4)))
        variation = history.compact_variation(action, AMPLITUDES)
        maximum_identity = max(maximum_identity, float(np.max(abs(curv['base_ricci_identity']))))
        rows[f'{nt}_{nx}'] = {
            'fixed_baseline_action_difference': {k: float((v-baseline[k]).real) for k, v in values.items()},
            'channel_action_gradients': {k: v.tolist() for k, v in variation['channel_gradients'].items()},
            'local_action_gradient': variation['local_action_gradient'].tolist(),
            'local_action_force': variation['local_action_force'].tolist(),
            'Euler_bulk_residual': variation['Euler_bulk_residual'],
            'fixed_boundary_derivatives': variation['fixed_boundary_derivatives'],
        }
    # Same finite local functional perturbed directly; no field equation solve.
    checks = {}; h = 1e-5
    for B in range(4):
        values = []
        for offset in (-2., -1., 1., 2.):
            shifted = AMPLITUDES+offset*h*np.eye(4)[B]
            raw, _ = action.actions(grid, *history.fields_and_jets(shifted)); values.append(raw)
        for name in rows['192_160']['channel_action_gradients']:
            fd = (values[0][name]-8*values[1][name]+8*values[2][name]-values[3][name])/(12*h)
            exact = rows['192_160']['channel_action_gradients'][name][B]
            checks.setdefault(name, []).append(float(abs(fd-exact)))
    limits = {k: limit_control(k) for k in ('static', 'homogeneous')}
    coordinate = independent_coordinate_jets()
    get = lambda key: np.array([rows[key]['channel_action_gradients'][c] for c in
        ('einstein_bulk', 'maxwell_bulk', 'weyl_bulk')])
    maxima = {
        'time_quadrature': float(np.max(abs(get('192_128')-get('128_128')))),
        'spatial_quadrature': float(np.max(abs(get('192_160')-get('192_128')))),
        'Euler_bulk_variation': rows['192_160']['Euler_bulk_residual'],
        'base_Ricci_identity': maximum_identity,
        'existing_static_homogeneous_limits': max(v for c in limits.values() for v in c.values()),
        'action_complex_step_vs_independent_variation': max(v for c in checks.values() for v in c),
        **{'metric_'+k: v for k, v in coordinate.items()},
    }
    tolerance = {k: (3e-11 if k in ('Euler_bulk_variation', 'base_Ricci_identity',
                    'existing_static_homogeneous_limits', 'metric_value') else 3e-8) for k in maxima}
    failures = {k: v for k, v in maxima.items() if v > tolerance[k]}
    return {
        'schema': 'NSC-SPHERICAL-LOCAL-HISTORY-v1', 'accountable_author': 'Douglas Ek',
        'status': 'PASS: same-action local variation on the supplied spherical spacetime history' if not failures else 'OPEN',
        'source_hashes': {p: sha(p) for p in SOURCES}, 'input_hashes': {p: sha(p) for p in INPUTS},
        'locked_inputs': ledger, 'evaluations': rows, 'maxima': maxima, 'tolerances': tolerance, 'failures': failures,
        'independent_action_variations': checks, 'inherited_curvature_limits': limits,
        'coordinate_verification': {'time_difference_step': 8e-5, 'radial_difference_step': 1e-3,
            'stencil': 'fourth-order centered; errors normalized by max(1,max_abs(exact_jet))',
            'smaller_time_step_cancellation_diagnostic': {'step': 2e-5, 'second_relative_error': 1.4142211456756512e-7}},
        'domain': {
            'history': 'same supplied raw-KS compact harmonic family, pulled into fixed PG coordinates',
            'amplitudes': AMPLITUDES.tolist(), 'omega': .4, 'time_support': [0., .05], 'radial_support': [-1., 1.],
            'history_selected': False, 'metric_coordinates': ['N_K', 'beta_K', 'a_K=q_ADM', 'r_K'],
            'signature': '+---; same Riemann sign as nsc_finite_terms',
            'action': '-4*pi*integral sqrt(-det(g2))*r^2*(A*R+C_gauge*q_mag^2/(2*r^4)+C_Weyl*C^2+C_Euler*E4+C_boxR*BoxR)',
            'GHY': 'owned geometric boundary term retained; compact family has zero boundary-jet variation',
            'Euler_BoxR': 'same constant coefficients; fixed boundary derivatives zero; Euler bulk integration is a diagnostic only',
            'baseline': 'S[g(alpha)]-S[g(0)] only for reporting; fixed g(0), not a vacuum normalization or subtraction',
            'metric_coordinate_jets': 'exact second-order Taylor products of the owned raw-KS metric and fixed chart',
            'force_convention': 'force is minus action gradient; both are stored',
            'multiplicity': '4*pi sphere integration only; coefficients already include their field content',
            'seam': 'one smooth geometry across fixed rho=0; no additional oriented term'},
        'gate': {'spatially_varying_local_action': 'PASS' if not failures else 'OPEN',
            'quantum_reference_history': 'OPEN: homogeneous Bloch subtraction is not the general spatially varying reference',
            'full_renormalized_source': 'OPEN', 'B2': 'OPEN', 'V_c': None, 'extended_stationarity': 'OPEN',
            'absolute_stress': None, 'nulls': None, 'updated_constraints': None, 'metric_timestep': False,
            'Gamma_rest_extra_independent_term': 'none introduced', 'Z3': 'OUT OF SCOPE',
            'Weyl_time_node_diagnostic': 93.54264532195464,
            'homogeneous_nonexistence': 'preserved in its original scope; no rerun'},
        'comparison': {'fields': 'all', 'float_atol': 3e-10, 'float_rtol': 3e-10, 'exact': 'hashes, structure, scope'},
        'reproducer': 'python3 scripts/derive_nsc_spherical_local_history.py --check',
    }


def main():
    p = argparse.ArgumentParser(); p.add_argument('--check', action='store_true'); args = p.parse_args()
    result = make_record()
    if args.check:
        from derive_nsc_transmitting_boundary_binding import compare
        compare(read(OUTPUT), result)
    else:
        (ROOT/OUTPUT).write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print(json.dumps({k: result[k] for k in ('status', 'maxima', 'failures')}, indent=2), flush=True)
    if result['failures']: raise SystemExit(1)


if __name__ == '__main__': main()
