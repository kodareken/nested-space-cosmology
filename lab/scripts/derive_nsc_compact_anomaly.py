#!/usr/bin/env python3
"""Apply the published doubled-Dirac anomaly to the existing compact domain.

Only the new domain/import algebra and invariant substitution are evaluated.
No spectrum, heat-coefficient derivation or old source generator is run.
"""
import argparse
import json
from pathlib import Path
import sys

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from recursive_horizons.nsc_spinor_bridge import weyl_matrices, pauli
from check_nsc_compact_boundary_action import authenticated_record
from check_nsc_compact_casimir import compare
from check_nsc_vacuum_charge_matching import hashes

OUTPUT = ROOT/'results/development/compact-anomaly-bridge.json'
INPUTS = ('results/development/compact-boundary-action.json',
          'results/development/canonical-spectral-bridge.json')
SOURCES = ('scripts/derive_nsc_compact_anomaly.py',
           'docs/nsc-compact-anomaly-bridge.md',
           'src/recursive_horizons/nsc_spinor_bridge.py',
           'scripts/check_nsc_compact_boundary_action.py',
           'scripts/check_nsc_compact_casimir.py',
           'scripts/check_nsc_vacuum_charge_matching.py')


def calculate():
    domain, bridge = (authenticated_record(x) for x in INPUTS)
    frame = weyl_matrices()
    t1, t2, t3 = pauli()
    ident = sp.eye(8)
    u = sp.kronecker_product(t1, sp.eye(4))
    star = sp.kronecker_product(t2, sp.eye(4))
    grading = sp.kronecker_product(t3, sp.eye(4))
    gamma = frame['gamma']
    tangent = [gamma[0], *(sp.I*x for x in gamma[1:])]
    normals = sp.kronecker_product(sp.eye(2), frame['gamma5'])
    rows = []
    for eta in (-1, 1):
        chi = -eta*sp.kronecker_product(t3, frame['gamma5'])
        p = (ident+chi)/2
        for epsilon in (-1, 1):
            gn = epsilon*u*normals
            residuals = {
                'copy_exchange_maps_allowed_to_adjoint': u*p*u-(ident-p),
                'new_normal_current_isotropic': p*gn*p,
                'new_normal_anticommutes_with_chi': chi*gn+gn*chi,
                'published_chi_with_component_orientation': chi-(-eta*epsilon)*sp.I*star*gn,
                'spectral_grading_preserves_domain': grading*p-p*grading,
                'spectral_grading_anticommutes_with_normal': grading*gn+gn*grading,
            }
            for a, g in enumerate(tangent):
                new = u*sp.kronecker_product(sp.eye(2), g)
                residuals[f'chi_commutes_with_new_tangent_{a}'] = chi*new-new*chi
                residuals[f'grading_anticommutes_with_new_tangent_{a}'] = grading*new+new*grading
            for value in residuals.values():
                assert value == sp.zeros(8)
            rows.append({'eta': eta, 'outward_normal': epsilon,
                         'bag_component_sign': -eta*epsilon,
                         'exact_residuals': {name: '0' for name in residuals}})

    # Published charges: PRD108 085015, Sec. VI.A.3; no rederivation.
    charges = [0, -sp.Rational(429, 4), 621, 0, 0, 270, 990, -210, -360]
    invariants = sp.symbols('E4 I1 I2 I3 I4 I5 I6 I7 I8')
    polynomial = sum(c*i for c, i in zip(charges, invariants))
    vanishing = {invariants[k]: 0 for k in (1, 2, 5, 6, 7, 8)}
    reduced = sp.expand(polynomial.subs(vanishing))
    assert reduced == 0
    return {
        'schema': 'NSC-COMPACT-ANOMALY-BRIDGE-v1',
        'status': 'published universal conformal anomaly applied to the paired compact modulus',
        'source_hashes': hashes(SOURCES), 'input_hashes': hashes(INPUTS),
        'primary_source': {
            'title': 'Fermions, boundaries, and conformal and chiral anomalies in d=3, 4 and 5 dimensions',
            'authors': 'Amin Faraji Astaneh and Sergey N. Solodukhin',
            'doi': '10.1103/PhysRevD.108.085015',
            'url': 'https://scoap3-prod-backend.s3.cern.ch/media/files/81024/10.1103/PhysRevD.108.085015.pdf',
            'equations': '80-89; charges Sec VI.A.3 p15; gauge Eq96',
        },
        'input_domain': domain['domain'],
        'domain_map': {'left_unitary': 'U=tau1 tensor I4; not a similarity',
                       'self_adjoint_representative': 'D_tilde=U D on the original allowed domain',
                       'closed_operator_square': 'D_tilde^2=D_dagger D',
                       'assumptions': 'identical torsionless bulk copies, common gauge representation commuting with U, original smooth mixed domains',
                       'pointwise_import_checks': rows},
        'published_charges': {str(k): str(v) for k, v in zip(invariants, charges)},
        'gravitational_anomaly_normalization': '1/[5760 (4 pi)^2]',
        'geometric_application': {
            'metric': 'exp(2 sigma(Y)) (g4+dY^2), g4 independent of Y',
            'boundary_shape': 'K_ab=(K/4)h_ab, hence K_hat=0',
            'mixed_Weyl_tensor': 'W_nabc=0, by the product geometry and conformal covariance',
            'vanishing_invariants': [str(k) for k in vanishing],
            'unrestricted_intrinsic_invariants': ['E4', 'I3', 'I4'],
            'universal_gravitational_anomaly': str(reduced),
            'gauge_condition': 'F_aY=0 (includes tangential Y-independent gauge backgrounds)',
            'universal_gauge_anomaly': '0: Eq96 is proportional to F_an F^an',
            'first_variation_argument': 'On the unwarped product K_hat=0, W_nabc=0, nabla_n W_anbn=0; all surviving invariants start at quadratic or higher order. Boundary-constant sigma preserves this under the conformal map.',
            'universal_Wess_Zumino_value_and_first_source': '0 in the stated representative and conformal-product class'},
        'phase_scope': {
            'spectral_grading': 'tau3 tensor I4',
            'eta_invariant': '0 for the self-adjoint representative with a discrete eta-regularizable spectrum',
            'zeta_phase_without_zero_modes': '1 when a5=0, closed tangential sections/no extra boundary terms, fixed spectral cut',
            'phase_identity': 'det_zeta D_tilde=det_zeta(D_tilde^2)^(1/2) exp[i*pi*zeta_(D_tilde^2)(0)/2]',
            'original_D_phase_fixed_by_left_unitarity_alone': False,
            'zero_modes_and_global_phase_on_noncompact_transition_resolved': False},
        'retained_finite_cutoff_warp_source': bridge['allocation']['warp_path_moment_derivatives'],
        'scope': {'universal_parity_even_conformal_piece_evaluated': True,
                  'finite_cutoff_warp_functional_zero': False,
                  'scheme_dependent_trivial_anomalies_set_to_zero': False,
                  'full_chiral_or_gauge_anomaly_cancellation_claimed': False,
                  'full_C5_plus_J_real_time_prescription_derived': False,
                  'self_sourced_solution_derived': False},
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
        print('Compact anomaly application reproduced: exact domain algebra and all fields agree.')
    else:
        with args.output.open('x') as handle:
            json.dump(result, handle, sort_keys=True, indent=2, allow_nan=False)
            handle.write('\n')
        print(f'Wrote {args.output}')


if __name__ == '__main__':
    main()
