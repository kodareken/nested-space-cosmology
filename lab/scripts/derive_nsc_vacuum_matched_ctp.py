#!/usr/bin/env python3
"""Construct one vacuum-matched real-time Dirac/geometry prescription."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from recursive_horizons.nsc_covariant_operator import smooth_metric, SIGMA3
from recursive_horizons.nsc_influence import canonical_hamiltonian, geometric_vertex
from recursive_horizons.nsc_response_matching import euclidean_kernel, retarded_susceptibility
from recursive_horizons.nsc_vacuum_matched_ctp import VacuumMatchedDiracAction, vacuum_matched_influence
from check_nsc_compact_boundary_action import authenticated_record
from check_nsc_compact_casimir import compare
from check_nsc_vacuum_charge_matching import hashes

OUTPUT = ROOT/'results/development/vacuum-matched-ctp.json'
INPUTS = ('results/development/gaussian-cutoff-measure.json',
          'results/nsc-11-response-matching.json')
SOURCES = ('scripts/derive_nsc_vacuum_matched_ctp.py',
           'src/recursive_horizons/nsc_vacuum_matched_ctp.py',
           'src/recursive_horizons/nsc_covariant_operator.py',
           'src/recursive_horizons/nsc_influence.py',
           'src/recursive_horizons/nsc_response_matching.py',
           'src/recursive_horizons/nsc_regulated.py',
           'docs/nsc-vacuum-matched-ctp.md',
           'scripts/check_nsc_compact_boundary_action.py',
           'scripts/check_nsc_compact_casimir.py',
           'scripts/check_nsc_vacuum_charge_matching.py')


def calculate():
    for path in INPUTS:
        authenticated_record(path)
    metric = smooth_metric(12)
    h, v = canonical_hamiltonian(metric), geometric_vertex(metric)
    action = VacuumMatchedDiracAction(h, v, np.kron(SIGMA3, np.eye(metric.points)), 4.)
    nu = .7
    counter = action.euclidean_counter_kernel(nu*nu)
    canonical = retarded_susceptibility(action.data, 1j*nu)
    raw = euclidean_kernel(metric, nu, 4., frequency_points=128)
    binding = counter+canonical-raw
    assert abs(binding) < 5e-10
    real_rows = []
    for w in (0., .4, 1.1, 2.2):
        value = action.continued_force_kernel(w)
        real_rows.append({'frequency': w, 'continued_force_kernel': value,
                          'geometric_branch_action_kernel': -value})
    period = 2*np.pi
    amplitudes = (.02, -.015)
    phases = [action.relative_branch_action(period, {0: a}) for a in amplitudes]
    plus, minus = (expm(-1j*period*(h+a*v)) for a in amplitudes)
    state = vacuum_matched_influence(action.reference_covariance, plus, minus, *phases)
    unchanged_modulus = abs(state['amplitude'])-abs(state['canonical_amplitude'])
    phase_identity = abs(np.exp(1j*state['principal_action'])-state['amplitude'])
    assert abs(unchanged_modulus) < 1e-13 and phase_identity < 1e-13
    pair = lambda z: [float(z.real), float(z.imag)]
    return {
        'schema': 'NSC-VACUUM-MATCHED-CTP-v1',
        'status': 'explicit working prescription; calculated full-frequency quadratic geometric branch in a finite Dirac channel',
        'source_hashes': hashes(SOURCES), 'input_hashes': hashes(INPUTS),
        'prescription': {
            'state': 'retain canonical Gaussian C0 and its unitary histories',
            'vacuum_matching': 'B_E=Gamma_heat_E-Gamma_canonical_vac_E, same normalization and domain',
            'continuation': 'continue finite short-proper-time vertex in nu^2 to -omega^2; no loop-contour rotation',
            'CTP': 'Z_phys=exp[i(B_L[+]-B_L[-])] Z_canonical[C0]',
            'state_differences': 'canonical at fixed geometry; B is fixed by zero-temperature vacuum matching',
            'extra_physical_rule': True,
            'uniquely_forced_by_Euclidean_determinant': False,
            'raw_thermal_determinant_retained_as_physical_partition_function': False,
            'new_fitted_coefficients': False},
        'domain': {'spatial_points': 12, 'geometry': 'existing smooth AP axial cell, R=2, a=1, N=q=1',
                   'field': 'one two-component radial block, kappa=1; full four-dimensional multiplicity 4*kappa is not applied',
                   'probe': 'r(t,x)=r0(x)/(1+J(t)s(x)); H=H0+J V',
                   'cutoff': 4., 'proper_time_and_Feynman_nodes': 40,
                   'approximation': 'quadratic in probe amplitude, full frequency dependence',
                   'units': 'hbar=c=1; kernel and static conversion have energy units'},
        'static_conversion': {'constant_energy': action.constant, 'linear_force': action.linear},
        'Euclidean_binding': {'frequency': nu, 'counter_kernel': pair(counter),
                              'canonical_kernel': pair(canonical), 'old_covariant_heat_kernel': raw,
                              'reconstruction_residual': pair(binding), 'old_frequency_quadrature_nodes': 128},
        'real_time_vertices': real_rows,
        'finite_history': {'period': float(period), 'constant_probe_amplitudes': list(amplitudes),
                           'relative_geometric_branch_actions': phases,
                           'canonical_amplitude': pair(state['canonical_amplitude']),
                           'prescription_amplitude': pair(state['amplitude']),
                           'geometric_phase_difference': state['geometric_phase_difference'],
                           'modulus_change': float(unchanged_modulus),
                           'principal_action': pair(state['principal_action']),
                           'amplitude_action_residual': float(phase_identity)},
        'sinusoidal_geometric_phase': {'period': float(period), 'J1_real': 0., 'J1_imag': .01,
                                       'relative_branch_action': action.relative_branch_action(period, {1: .01j})},
        'scope': {'normalized_positive_overlap_kernel_at_finite_regulator': True,
                  'canonical_matter_dynamics_changed': False,
                  'bare_geometric_vertex_is_a_retarded_correlator': False,
                  'coupled_metric_causality_or_stability_proved': False,
                  'all_orders_interacting_nonlocal_unitarity_imported': False,
                  'full_curved_nonperturbative_continuation_evaluated': False,
                  'absolute_continuum_source_or_self_sourced_geometry_solved': False},
        'comparison': {'fields': 'all', 'float_atol': 3e-9, 'float_rtol': 3e-8,
                       'exact': 'structure, strings, non-float values, source and input hashes',
                       'exceptions': []},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true')
    parser.add_argument('--output', type=Path, default=OUTPUT)
    args = parser.parse_args()
    result = calculate()
    if args.check:
        compare(json.loads(args.output.read_text()), result)
        print('Vacuum-matched CTP construction reproduced: all fields agree; no old generator ran.')
    else:
        with args.output.open('x') as handle:
            json.dump(result, handle, indent=2, sort_keys=True, allow_nan=False)
            handle.write('\n')
        print(f'Wrote {args.output}')


if __name__ == '__main__':
    main()
