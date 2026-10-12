#!/usr/bin/env python3
"""Verify only the charged light allocation and its canonical-KS action."""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys

import mpmath as mp
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/'src'), str(ROOT/'tests')]
from recursive_horizons.nsc_magnetic_light_reference import magnetic_light_spectrum, homogeneous_magnetic_restoration
from recursive_horizons.nsc_light_restoration_action import LightRestorationAction, conformal_metric_jets, weyl_euler_average
from recursive_horizons.nsc_angular_stress import static_cylinder_density, profile_jets
from recursive_horizons.nsc_spherical_local_history import SphericalHistoryGrid
from test_nsc_magnetic_light_reference import continued_zeta, physical_lll_geometry
from test_nsc_light_restoration_action import homogeneous_variation, full_geometric_tensor


OUTPUT = 'results/development/nsc-magnetic-light-restoration.json'
SOURCES = (
    'src/recursive_horizons/nsc_magnetic_light_reference.py',
    'src/recursive_horizons/nsc_light_restoration_action.py',
    'tests/test_nsc_magnetic_light_reference.py',
    'tests/test_nsc_light_restoration_action.py',
    'scripts/check_nsc_magnetic_light_restoration.py',
    'docs/nsc-magnetic-light-restoration.md',
)
OWNERS = (
    'src/recursive_horizons/nsc_angular_stress.py',
    'src/recursive_horizons/nsc_spherical_local_history.py',
    'src/recursive_horizons/nsc_horizon_source.py',
    'src/recursive_horizons/nsc_lll_geometric_history.py',
    'src/recursive_horizons/nsc_compact_matching.py',
    'docs/nsc-angular-stress.md', 'docs/nsc-horizon-source.md',
    'docs/nsc-lll-geometric-history.md', 'docs/nsc-compact-matching.md',
    'docs/nsc-declared-action-scope.md', 'docs/nsc-common-ks-trace.md',
)
POINTERS = (
    'results/development/angular-stress.json',
    'results/development/nsc-compact-matched-restart.json',
    'results/development/charged-ctp-neck-source.json',
    'results/development/nsc-mode-resolved-cauchy-state.json',
    'results/development/nsc-lll-geometric-history.json',
    'results/development/nsc-spherical-local-history.json',
)


def sha(path): return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
def read(path): return json.loads((ROOT/path).read_text())


def spectral_row(spectrum):
    return {
        'charge': spectrum.charge, 'Z_minus_one': str(spectrum.z_minus_one),
        'Z_nonzero_zero': str(spectrum.z_nonzero_zero), 'H_q': str(spectrum.harmonic_number),
        'FP_Z_one': float(spectrum.finite_part_one),
        'Zprime_minus_one': float(spectrum.zprime_minus_one),
        'cylinder_density_R1_mu1': float(spectrum.cylinder_density()),
        'fourth_order_harmonic_mu1': float(spectrum.fourth_order_harmonic),
        'series_terms': spectrum.series_terms, 'last_series_index': spectrum.last_series_index,
        'Zprime_series_tail_bound_exact': str(spectrum.derivative_tail_bound),
        'Zprime_series_tail_bound': float(spectrum.derivative_tail_bound),
        'cylinder_series_error_bound': float(spectrum.cylinder_series_error_bound()),
        'decimal_precision': spectrum.decimal_precision,
    }


def action_control(spectrum, nt, nz):
    grid = SphericalHistoryGrid.gaussian(nt, nz, (-2.4, -1.3), (-.7, .7))
    owner = LightRestorationAction(spectrum)
    g, rgrid, jets, envelope, (N, a, r) = homogeneous_variation(grid, np.zeros(4))
    base_actions, _ = owner.actions(grid, g, rgrid, jets, canonical_coordinates='KS_time_reparametrization')
    rho, parallel, sphere = np.array([full_geometric_tensor(theta, spectrum) for theta in -grid.time]).T
    density = np.stack((-4*np.pi*a*r*r*rho, np.zeros_like(a),
                        4*np.pi*N*r*r*parallel, 8*np.pi*N*a*r*sphere), axis=-1)
    expected = grid.integral(density[:, None, :]*envelope[..., None])
    per_channel = {name: np.zeros(4) for name in base_actions}; step = 1e-28
    for B in range(4):
        amplitudes = np.zeros(4, complex); amplitudes[B] = 1j*step
        G, R, J, _, _ = homogeneous_variation(grid, amplitudes)
        actions, _ = owner.actions(grid, G, R, J, canonical_coordinates='KS_time_reparametrization')
        for name, value in actions.items(): per_channel[name][B] = value.imag/step
    actual = sum(per_channel.values())
    return {
        'grid': [nt, nz], 'reference_actions': {k: float(v) for k, v in base_actions.items()},
        'action_gradient_by_channel': {k: v.tolist() for k, v in per_channel.items()},
        'action_gradient': actual.tolist(), 'independent_tensor_pairing': expected.tolist(),
        'absolute_pairing_residual': abs(actual-expected).tolist(),
    }


def homogeneous_controls(spectrum):
    rows = []; h = 8e-5
    for theta in (1.4, np.pi/2, 2., 3*np.pi/4):
        row = homogeneous_magnetic_restoration(theta, spectrum)
        lll = physical_lll_geometry(theta, spectrum.charge)
        trace = row['rho']+lll[0]-row['p_parallel']-lll[1]-2*row['p_sphere']
        target = row['neutral_gravitational_trace']+row['gauge_trace']
        ward = {}
        for include_lll in (False, True):
            def tensor(x):
                value = homogeneous_magnetic_restoration(x, spectrum)
                base = np.array([value['rho'], value['p_parallel'], value['p_sphere']])
                return base+(physical_lll_geometry(x, spectrum.charge) if include_lll else 0)
            derivative = (tensor(theta-2*h)[0]-8*tensor(theta-h)[0]
                          +8*tensor(theta+h)[0]-tensor(theta+2*h)[0])/(12*h)
            rho, parallel, sphere = tensor(theta)
            W, W1, _, _, _ = profile_jets(theta); sigma1 = -np.cos(theta)/np.sin(theta)
            connection = (sigma1+W1/(2*W))*(rho+parallel)+2*sigma1*(rho+sphere)
            ward['with_physical_LLL' if include_lll else 'nonzero_angular_only'] = {
                'absolute': float(abs(derivative+connection)),
                'scaled': float(abs(derivative+connection)/(1+abs(derivative)+abs(connection))),
            }
        rows.append({
            'theta': float(theta), 'rho_coordinate': row['rho_coordinate'], 'sphere_radius': row['radius'],
            'nonzero_angular_restoration_trace': row['trace'],
            'physical_LLL_geometric_trace': float(lll[0]-lll[1]),
            'total_light_restoration_trace': float(trace), 'charged_light_anomaly': target,
            'trace_residual': float(abs(trace-target)),
            'LLL_trace_allocation_residual': float(abs(lll[0]-lll[1]-row['excluded_physical_LLL_trace'])),
            'Ward_residuals': ward,
        })
    return rows


def make_record():
    locked = read('results/development/nsc-compact-matched-restart.json')['locked_inputs']
    if locked['magnetic_flux'] != 4: raise ValueError('retained q=4 branch required')
    before = {p: sha(p) for p in (*SOURCES, *OWNERS, *POINTERS)}
    charged = magnetic_light_spectrum(4, series_tolerance='1e-35', decimal_precision=80)
    neutral = magnetic_light_spectrum(0)
    coarse = magnetic_light_spectrum(4, series_tolerance='1e-22', decimal_precision=60)
    with mp.workdps(80):
        step = mp.mpf('1e-15'); terms = charged.last_series_index+12
        zm = continued_zeta(-1+1j*step, 4, terms)
        z0 = continued_zeta(1j*step, 4, terms)
        fp = continued_zeta(1+1j*step, 4, terms)-2/(1j*step)
        spectral_residuals = {
            'Z_minus_one': float(abs(zm.real+mp.mpf(79)/30)),
            'Z_nonzero_zero': float(abs(z0.real+mp.mpf(13)/3)),
            'FP_Z_one': float(abs(fp.real-charged.finite_part_one)),
            'Zprime_minus_one': float(abs(zm.imag/step-charged.zprime_minus_one)),
        }
        derivative_change = coarse.zprime_minus_one-charged.zprime_minus_one
        bound = mp.mpf(coarse.derivative_tail_bound.numerator)/coarse.derivative_tail_bound.denominator
        convergence = {'coarse': spectral_row(coarse), 'Zprime_coarse_minus_resolved': float(derivative_change),
                       'positive_difference_within_analytic_tail_bound': bool(0 < derivative_change < bound),
                       'bound_scope': 'exact spectral-series truncation bound; not a floating-point interval certificate'}
    q = charged.charge
    exact = {
        'gauge_a4_trace_matching': str(-(charged.z_minus_one-neutral.z_minus_one)/8-Fraction(q*q, 48)),
        'full_angular_zero_coefficient_after_LLL': str(charged.z_nonzero_zero+q+Fraction(1, 3)),
        'gauge_a4_F_squared_coefficient': '2/3',
        'monopole_F_squared_times_r4': str(Fraction(q*q, 2)),
        'gauge_trace_times_pi2_r4': str(Fraction(q*q, 48)),
        'FP_Z_one_formula': '4*EulerGamma-2*H_q',
        'fourth_harmonic_formula_mu1': 'EulerGamma/2-H_q/4',
    }
    neutral_residual = abs(float(neutral.cylinder_density())-static_cylinder_density())
    neutral_harmonic = abs(float(neutral.fourth_order_harmonic)-float(mp.euler)/2)
    action = {'q0/128_128': action_control(neutral, 128, 128),
              'q4/96_96': action_control(charged, 96, 96),
              'q4/128_128': action_control(charged, 128, 128)}
    controls = homogeneous_controls(charged)
    grid = SphericalHistoryGrid.gaussian(12, 12, (-2.4, -1.3), (-.7, .7))
    g, radius, jets, _, _ = homogeneous_variation(grid, np.array([.001, .002, .003, .001]))
    three = weyl_euler_average(grid, g, radius, jets)
    five = weyl_euler_average(grid, g, radius, jets, quadrature_order=5)
    G, R, J = conformal_metric_jets(g, radius, jets, 1.)
    physical_endpoint = max(float(np.max(abs(G-g))), float(np.max(abs(R-radius))),
                            *(float(np.max(abs(a-b))) for a, b in zip(J, jets)))
    G, R, J = conformal_metric_jets(g, radius, jets, 0.)
    barred_endpoint = max(float(np.max(abs(G*radius[..., None, None]**2-g))),
                          float(np.max(abs(R-1))), float(np.max(abs(J[2]))), float(np.max(abs(J[3]))))
    maxima = {
        'spectral_off_pole_continuation': max(spectral_residuals.values()),
        'neutral_cylinder_normalization': float(neutral_residual),
        'neutral_fourth_harmonic': float(neutral_harmonic),
        'action_vs_tensor_resolved': max(max(action[k]['absolute_pairing_residual']) for k in ('q0/128_128', 'q4/128_128')),
        'action_quadrature_96_to_128': float(np.max(abs(np.array(action['q4/96_96']['action_gradient'])-action['q4/128_128']['action_gradient']))),
        'charged_trace': max(row['trace_residual'] for row in controls),
        'LLL_trace_allocation': max(row['LLL_trace_allocation_residual'] for row in controls),
        'source_conservation_scaled': max(item['scaled'] for row in controls for item in row['Ward_residuals'].values()),
        'formal_Weyl_quadrature_relative': float(np.max(abs(three-five)/(1+abs(five)))),
        'formal_Weyl_path_endpoints': max(physical_endpoint, barred_endpoint),
    }
    tolerances = {name: 3e-14 for name in maxima}
    tolerances.update(spectral_off_pole_continuation=1e-25, action_vs_tensor_resolved=3e-11,
                      action_quadrature_96_to_128=3e-9, source_conservation_scaled=3e-9,
                      formal_Weyl_quadrature_relative=3e-12)
    failures = {k: v for k, v in maxima.items() if not np.isfinite(v) or v > tolerances[k]}
    for name in ('gauge_a4_trace_matching', 'full_angular_zero_coefficient_after_LLL'):
        if exact[name] != '0': failures[name] = exact[name]
    if not convergence['positive_difference_within_analytic_tail_bound']: failures['spectral_series_convergence'] = False
    if before != {p: sha(p) for p in before}: raise ValueError('restoration dependency changed during the focused check')
    return {
        'schema': 'NSC-MAGNETIC-LIGHT-RESTORATION-v1', 'accountable_author': 'Douglas Ek',
        'status': 'PASS: charged light allocation and canonical-KS restoration action' if not failures else 'OPEN',
        'source_hashes': {p: before[p] for p in SOURCES}, 'dependency_hashes': {p: before[p] for p in OWNERS},
        'frozen_result_pointers': {p: {'sha256': before[p], 'generator_rerun': False} for p in POINTERS},
        'locked_inputs': {'magnetic_flux': 4, 'normalization_mu': 1., 'physical_charge_scan': False},
        'charged_spectrum': spectral_row(charged), 'neutral_control': spectral_row(neutral),
        'exact_matching_identities': exact, 'off_pole_spectral_residuals': spectral_residuals,
        'series_convergence': convergence, 'homogeneous_controls': controls,
        'action_variations': action, 'maxima': maxima, 'tolerances': tolerances, 'failures': failures,
        'domain': {
            'source': 'existing charged massless 4D Dirac light field; existing normalization and magnetic spectrum',
            'spectral_continuation': 'Hurwitz expansion with the j=2 pole-times-zero term retained',
            'canonical_foliation': 'same KS spatial slices as C0; KS time reparameterization allowed; PG inputs rejected',
            'metric_input': 'arbitrary supplied spherical metric with analytic first/second jets',
            'formal_Weyl_parameter': 's in [0,1], integrated with three Gauss nodes; not physical time or a field',
            'action_variation_fields': ['N', 'beta', 'a', 'r'],
            'control_coordinates': {'time': 't=-theta on the old geometry', 'interval': [-2.4, -1.3],
                                    'spatial_interval': [-.7, .7], 'physical_duration_selected': False},
            'action_channels': ['cylinder', 'bar_radial_R2_squared', 'LLL_geometry_bar', 'WZ_Weyl', 'WZ_Euler', 'WZ_boxR', 'WZ_gauge'],
            'tensor_action_identity': 'nonzero angular restoration plus physical LLL geometry = full light restoration action variation',
            'LLL_counting': 'action includes geometric LLL once; add LLL state only, not another geometric allocation',
            'compact_induced_counting': 'excluded; locked compact local action is added separately once',
            'boundary_scope': 'compact boundary-identical variation controls; general boundary terms are not discarded by this record',
            'new_independently_weighted_term': False, 'old_scientific_outputs_modified': False,
        },
        'gate': {'matched_light_allocation_action': 'PASS' if not failures else 'OPEN',
                 'complete_incoming_stress': 'OPEN: state-minus-reference source, remaining spectral error and complete same-slice ledger',
                 'complete_stress_tensor': None, 'constraints_solved': False, 'physical_metric_evolution': False,
                 'state_or_occupations_changed': False, 'physical_Cauchy_data_selected': False,
                 'PDF_or_public_release_changed': False},
        'comparison': {'fields': 'all', 'float_atol': 3e-13, 'float_rtol': 3e-12,
                       'exact': 'fractions, hashes, structure, scope and counting declarations'},
        'reproducer': 'python3 scripts/check_nsc_magnetic_light_restoration.py --check',
    }


def compare(expected, actual, path='$'):
    if type(expected) is not type(actual): raise AssertionError('type mismatch '+path)
    if isinstance(expected, dict):
        if expected.keys() != actual.keys(): raise AssertionError('keys differ '+path)
        for key in expected: compare(expected[key], actual[key], path+'.'+key)
    elif isinstance(expected, list):
        if len(expected) != len(actual): raise AssertionError('list length differs '+path)
        for i, (a, b) in enumerate(zip(expected, actual)): compare(a, b, path+f'[{i}]')
    elif isinstance(expected, float):
        if not np.isfinite(actual) or not np.isclose(expected, actual, atol=3e-13, rtol=3e-12):
            raise AssertionError('numeric mismatch '+path)
    elif expected != actual: raise AssertionError('value mismatch '+path)


def main():
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', action='store_true'); mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.write and (ROOT/OUTPUT).exists(): raise ValueError('new record already exists; use the all-field --check replay')
    result = make_record()
    if result['failures']:
        print(json.dumps({'status': result['status'], 'failures': result['failures']}, indent=2))
        raise SystemExit(1)
    if args.check: compare(read(OUTPUT), result)
    else: (ROOT/OUTPUT).write_text(json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+'\n')
    print(json.dumps({k: result[k] for k in ('status', 'maxima', 'failures')}, indent=2), flush=True)


if __name__ == '__main__': main()
