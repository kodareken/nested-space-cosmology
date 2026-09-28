#!/usr/bin/env python3
"""Bind the finite proper-time profile to its Gaussian measure exactly."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
from scipy.special import exp1
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from recursive_horizons.nsc_gaussian_measure import warped_gaussian_measure
from recursive_horizons.nsc_regulated import OperatorConventions, RegulatedOperator
from check_nsc_compact_boundary_action import authenticated_record
from check_nsc_compact_casimir import compare
from check_nsc_vacuum_charge_matching import hashes

OUTPUT = ROOT/'results/development/gaussian-cutoff-measure.json'
INPUTS = ('results/development/canonical-spectral-bridge.json',
          'results/development/compact-anomaly-bridge.json')
SOURCES = ('scripts/derive_nsc_gaussian_measure.py',
           'src/recursive_horizons/nsc_gaussian_measure.py',
           'src/recursive_horizons/nsc_regulated.py',
           'docs/nsc-gaussian-cutoff-measure.md',
           'scripts/check_nsc_compact_boundary_action.py',
           'scripts/check_nsc_compact_casimir.py',
           'scripts/check_nsc_vacuum_charge_matching.py')


def calculate():
    for path in INPUTS:
        authenticated_record(path)
    d, cutoff, m = sp.symbols('d Lambda M', positive=True)
    z = d*d/cutoff**2
    ein = sp.EulerGamma+sp.log(z)+sp.expint(1, z)
    logkin = sp.log(d)-ein/2
    target = sp.expint(1, z)/2+sp.log(m/cutoff)+sp.EulerGamma/2
    residuals = {
        'profile_reconstruction': -logkin+sp.log(m)-target,
        'fixed_cutoff_logarithmic_derivative': sp.diff(logkin, d)-sp.exp(-z)/d,
        'cutoff_derivative_with_normalization':
            sp.diff(-logkin, cutoff)-(sp.exp(-z)-1)/cutoff,
    }
    exact = {}
    for name, value in residuals.items():
        assert sp.simplify(value) == 0, name
        exact[name] = '0'
    # One finite noncommuting congruence checks ordering and measure signs.
    # It is an algebraic integration example, not a spacetime source model.
    d0 = np.array([[1.2, .4+.2j, -.1, .3], [.4-.2j, -1.7, .2j, .1],
                   [-.1, -.2j, 2.1, -.25j], [.3, .1, .25j, -.8]], complex)
    sigma = np.diag([.13, -.07, .04, -.11])
    lam, mu = 1.6, .9
    data = warped_gaussian_measure(d0, sigma, lam, mu)
    a, k, ds = (data[x] for x in ('canonical_transform', 'kinetic', 'covariant_operator'))
    previous = RegulatedOperator(ds, OperatorConventions(cutoff=lam, normalization=mu))
    heat = float(.5*np.sum(exp1((previous.values/lam)**2)))
    action_from_jacobian = data['canonical_action']+data['measure_action']
    numerical = {
        'ordered_action_error': float(np.max(abs(a.conj().T@d0@a-k))),
        'normalized_determinant_error': data['kinetic_action']-previous.action(),
        'jacobian_action_error': action_from_jacobian-data['kinetic_action'],
        'raw_heat_reconstruction_error': data['raw_heat_action_from_measure']-heat,
        'noncommuting_warp_factor_norm': float(np.linalg.norm(data['warp']@data['factor_root']-data['factor_root']@data['warp'])),
        'normalized_action': previous.action(), 'raw_heat_action': heat,
        'canonical_action': data['canonical_action'], 'measure_action': data['measure_action'],
        'normalization_shift': data['normalization_shift'],
        'euclidean_kinetic_bound': data['euclidean_kinetic_bound'],
    }
    for name in ('ordered_action_error', 'normalized_determinant_error',
                 'jacobian_action_error', 'raw_heat_reconstruction_error'):
        assert abs(numerical[name]) < 2e-12, name
    assert numerical['noncommuting_warp_factor_norm'] > 1e-3
    return {
        'schema': 'NSC-GAUSSIAN-CUTOFF-MEASURE-v1',
        'status': 'exact finite Gaussian kinetic representation and ordered measure Jacobian',
        'source_hashes': hashes(SOURCES), 'input_hashes': hashes(INPUTS),
        'primary_identity': 'DLMF 6.2.3-6.2.4: Ein(z)=E1(z)+log(z)+EulerGamma',
        'primary_url': 'https://dlmf.nist.gov/6.2', 'exact_identities': exact,
        'construction': {
            'kinetic': 'K_Lambda(D)=D exp[-Ein(D^2/Lambda^2)/2]',
            'raw_action': 'Gamma_heat=-log|det(K_Lambda/M)|-N[log(M/Lambda)+EulerGamma/2]',
            'warp_order': 'D_s=W D0 W; W=exp(-sigma/2); chi=W F_s^(1/2) psi',
            'negative_log_Grassmann_Jacobian': 'Tr sigma + (1/2)Tr Ein(D_s^2/Lambda^2)',
            'source_variation_fixed_rank_cutoff_M': '-Tr[exp(-D^2/Lambda^2) D^-1 delta D]',
            'scale_covariance': 'K_(Omega Lambda)(Omega U D U†)=Omega U K_Lambda(D) U† for unitary U',
            'normalized_positive_operator': 'K_Lambda^2/(Lambda^2 exp(-EulerGamma))=exp[-E1(D^2/Lambda^2)]',
            'Fredholm_modulus': 'Gamma_heat=-(1/2)log det_F exp[-E1(D^2/Lambda^2)] when E1(D^2/Lambda^2) is trace class',
        },
        'example': {'dimension': 4, 'cutoff': lam, 'normalization': mu,
                    'reference_real': d0.real.tolist(), 'reference_imag': d0.imag.tolist(),
                    'sigma_diagonal': np.diag(sigma).tolist(), **numerical},
        'scope': {'new_parameters_or_fields': False,
                  'canonical_real_time_Hamiltonian_replaced': False,
                  'continuum_trace_terms_separately_finite': False,
                  'continuum_normalized_modulus_has_trace_class_condition': True,
                  'extra_zero_of_entire_factor': False,
                  'pole_agreement_proves_state_equivalence': False,
                  'real_time_contour_initial_state_and_observables_fixed': False,
                  'full_physical_closure_derived': False},
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
        print('Gaussian cutoff measure reproduced: all identities and record fields agree.')
    else:
        with args.output.open('x') as handle:
            json.dump(result, handle, indent=2, sort_keys=True, allow_nan=False)
            handle.write('\n')
        print(f'Wrote {args.output}')


if __name__ == '__main__':
    main()
