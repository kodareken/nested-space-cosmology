#!/usr/bin/env python3
"""Bounded derivative control: supplied geometry, initial columns and inflow.

Synthetic finite columns here are numerical controls only. No spectral source,
physical parent, metric solution or incoming C0 match is supplied by this run.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from recursive_horizons.nsc_compatible_history_geometry import (
    CompatibleRadiusDirection, CompatibleIncomingMetric,
)
from recursive_horizons.nsc_prepared_history_jets import (
    evolve_prepared_field_jets, restriction_covariance_tangent,
)
from recursive_horizons.nsc_transmitting_history_jets import FourthOrderModePropagator
from recursive_horizons.nsc_transmitting_dirac_domain import TransmittingDiracSeamDomain, MODE_TO_CURRENT
from recursive_horizons.nsc_ks_spacetime_variation import chart_coordinates

OUTPUT = ROOT/'results/development/nsc-compatible-prepared-history.json'
SOURCES = (
    'src/recursive_horizons/nsc_compatible_history_geometry.py',
    'src/recursive_horizons/nsc_prepared_history_jets.py',
    'src/recursive_horizons/nsc_transmitting_history_jets.py',
    'src/recursive_horizons/nsc_transmitting_history_modes.py',
    'src/recursive_horizons/nsc_ks_spacetime_variation.py',
    'src/recursive_horizons/nsc_pg_ks_metric_pullback.py',
    'src/recursive_horizons/nsc_transmitting_dirac_domain.py',
    'scripts/derive_nsc_compatible_prepared_history.py',
    'tests/test_nsc_compatible_prepared_history.py',
    'docs/nsc-compatible-prepared-history.md',
)


def wave(frequency, phase=0., scale=1.):
    return lambda z, order: float(scale*np.real((1j*frequency)**order*np.exp(1j*(frequency*z+phase))))


def control_problem():
    x = np.linspace(-2., 1.2, 65)
    owner = FourthOrderModePropagator(x, 1.3, 2.)
    directions = (
        CompatibleRadiusDirection(wave(.7), wave(.4, .3), .015, .08),
        CompatibleRadiusDirection(wave(.5, .6), wave(.3, .9), .015, .08),
    )
    theta = np.array([.003, -.002])
    times = np.linspace(0., .04, 5)  # caller's numerical control window only
    base = np.vstack((np.column_stack((np.exp(-x*x), np.exp(.3j*x))),
                      np.column_stack((.2*np.exp(.2j*x), np.exp(-.3*x*x)))))
    initial_tangents = np.array([(.1+.03j)*base, (.04-.08j)*base[:, ::-1]])
    n, sources = len(x), base.shape[1]
    boundary = np.array([n-1, 2*n-1])
    incident = np.array([[.4+.1j, -.2j], [.2, .3-.1j]])
    incident_tangents = np.array([.15j*incident, .08*incident[:, ::-1]])
    N, beta, q, _ = owner.reference
    speed = np.array([N[-1]/q[-1]-beta[-1], -N[-1]/q[-1]-beta[-1]])
    sat = -speed/owner.weights[-1]
    frequencies = np.array([.2, .6])

    def inputs(at_theta, *, omit_initial=False, omit_inflow=False, omit_geometry=False):
        provider = CompatibleIncomingMetric(tuple(at_theta), directions)
        if omit_geometry:
            provider = SimpleNamespace(values=provider.values,
                log_directions=lambda time, x, metric: np.zeros((2, 4, len(x))))
        initial = base+np.einsum('d,dij->ij', at_theta, initial_tangents)
        def incoming(time):
            phase = np.exp(-1j*frequencies*time)
            values = incident+np.einsum('d,dij->ij', at_theta, incident_tangents)
            forcing = np.zeros((2*n, sources), complex)
            forcing[boundary] = sat[:, None]*values*phase
            tangents = np.zeros((2, 2*n, sources), complex)
            if not omit_inflow:
                tangents[:, boundary] = sat[None, :, None]*incident_tangents*phase
            return forcing, tangents
        initial_jets = np.zeros_like(initial_tangents) if omit_initial else initial_tangents
        return provider, initial, initial_jets, incoming

    incoming_index = np.flatnonzero(x == 1.)
    if len(incoming_index) != 1:
        raise ValueError('an exact rho1 grid node is required; no implicit interpolation')
    sample_rows = np.array([incoming_index[0], n+incoming_index[0]])
    def solve(at_theta, **kwargs):
        return evolve_prepared_field_jets(owner, times, *inputs(at_theta, **kwargs),
                                          sample_rows=sample_rows)
    return owner, theta, times, solve, directions


def run_control():
    owner, theta, times, solve, directions = control_problem()
    result = solve(theta)
    step = 2e-5
    runs = [[solve(theta+sign*step*np.eye(2)[direction])
              for sign in (1., -1.)] for direction in range(2)]
    parts = [[run['evolved_field'] for run in pair] for pair in runs]
    finite = np.array([(plus-minus)/(2*step) for plus, minus in parts])
    residual = float(np.max(abs(finite-result['tangent_fields'])))
    no_geometry = solve(theta, omit_geometry=True)
    no_initial = solve(theta, omit_initial=True)
    no_inflow = solve(theta, omit_inflow=True)
    omission = {
        'metric_geometry': float(np.max(abs(result['tangent_fields']-no_geometry['tangent_fields']))),
        'initial_preparation': float(np.max(abs(result['tangent_fields']-no_initial['tangent_fields']))),
        'incoming_preparation': float(np.max(abs(result['tangent_fields']-no_inflow['tangent_fields']))),
    }
    # Restrict the actual time-dependent characteristic columns at rho1.
    # Intrinsic metric and normal are fixed there for this radius family.
    index = int(np.flatnonzero(owner.x == 1.)[0])
    N, beta, q, radius = owner.reference[:, index]
    domain = TransmittingDiracSeamDomain(N, q, beta, radius,
                                        surface_id='fixed-incoming-rho1-control')
    _, inverse, _ = domain.trace_map(np.ones(1))
    restriction = inverse[0] @ MODE_TO_CURRENT/(radius*np.sqrt(q))
    def trace(run):
        return np.einsum('ab,tbc->tac', restriction, run['sampled_fields']).reshape(-1, 2)
    F = trace(result)
    dF = np.einsum('ab,dtbc->dtac', restriction,
                   result['sampled_tangent_fields']).reshape(2, -1, 2)
    trace_finite = np.array([(trace(plus)-trace(minus))/(2*step) for plus, minus in runs])
    trace_error = float(np.max(abs(trace_finite-dF)))
    Cbase = np.array([[.45, .07+.04j], [.07-.04j, .6]])
    dC = np.array([[[.02, .01j], [-.01j, -.03]],
                   [[-.01, .015+.005j], [.015-.005j, .025]]])
    C = Cbase+np.einsum('d,dij->ij', theta, dC)
    analytic = restriction_covariance_tangent(F, dF, C, dC)
    covfinite = []
    for direction, (plus, minus) in enumerate(runs):
        Cp, Cm = C+step*dC[direction], C-step*dC[direction]
        Fp, Fm = trace(plus), trace(minus)
        covfinite.append((Fp@Cp@Fp.conj().T-Fm@Cm@Fm.conj().T)/(2*step))
    covariance_error = float(np.max(abs(analytic-np.asarray(covfinite))))
    fixed_source = restriction_covariance_tangent(F, dF, C, np.zeros_like(dC))
    source_term = float(np.max(abs(analytic-fixed_source)))
    flux_keys = [k for k in result if 'residual' in k and np.isscalar(result[k]) and result[k] is not None]
    algebra = {k: float(result[k]) for k in flux_keys}
    tolerance = 3e-8
    if max(residual, trace_error, covariance_error) >= tolerance:
        raise ArithmeticError('complete supplied-family derivative fails finite-difference control')
    if min(omission.values()) <= 1e-6 or source_term <= 1e-6:
        raise ArithmeticError('preparation omission controls are not discriminating')
    if any(value >= 3e-11 for key, value in algebra.items() if 'flux' in key):
        raise ArithmeticError('same-generator flux algebra failed')
    return {
        'schema': 'NSC-COMPATIBLE-PREPARED-HISTORY-v1',
        'accountable_author': 'Douglas Ek',
        'status': 'PASS: supplied-family derivative control; physical preparation and extended stationarity OPEN',
        'control': {'synthetic_columns_only': True, 'grid_points': len(owner.x),
                    'source_columns': 2, 'directions': 2, 'midpoint_steps': len(times)-1,
                    'PG_time_window': [float(times[0]), float(times[-1])],
                    'selected_physical_duration': None,
                    'incoming_center': 's=T(rho)-T(1); z=tau+S(rho)',
                    'numerical_normal_window_radii': [.015, .08],
                    'finite_difference_step': step,
                    'field_derivative_tolerance': tolerance,
                    'covariance_derivative_tolerance': tolerance,
                    'flux_algebra_tolerance': 3e-11,
                    'covariance_projection': 'actual rho1 normal/half-density trace of synthetic time-dependent columns',
                    'incoming_z_nodes': (times+chart_coordinates(1.)[1]).tolist(),
                    'stationary_phase_reapplied': False,
                    'incoming_trace_interpolation': None},
        'residuals': {'complete_field_derivative': residual,
                      'incoming_trace_derivative': trace_error,
                      'complete_coherent_covariance_derivative': covariance_error, **algebra},
        'omitted_term_effects': {**omission, 'source_covariance': source_term},
        'scope': {'physical_parent_preparation': 'OPEN',
                  'incoming_C0_match': 'OPEN', 'continuous_spectral_completion': 'OPEN',
                  'continuum_error_bound': None, 'physical_endpoint_family': None,
                  'global_initial_data': 'OPEN', 'extended_stationarity': 'OPEN',
                  'metric_timestep': False, 'source_stress_computed': False,
                  'fit_or_added_action_term': False},
    }


def provenance():
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCES}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true', help='authenticate the completed control without rerunning it')
    args = parser.parse_args()
    if args.check:
        saved = json.loads(OUTPUT.read_text())
        if saved['source_sha256'] != provenance():
            raise ValueError('control provenance changed; inspect the affected dependency before rerunning')
        print(json.dumps({'status': saved['status'], 'residuals': saved['residuals']}, indent=2))
    else:
        record = run_control()
        record['source_sha256'] = provenance()
        OUTPUT.write_text(json.dumps(record, indent=2, allow_nan=False)+'\n')
        print(json.dumps({'status': record['status'], 'residuals': record['residuals'],
                          'omitted_term_effects': record['omitted_term_effects']}, indent=2))


if __name__ == '__main__':
    main()
