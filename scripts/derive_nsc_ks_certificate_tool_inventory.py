#!/usr/bin/env python3
"""Inventory of certificate tools other than the four spent currents.

No Dirac evolution. L(E), Pi-prime plus E_beta, the integrated defect and
the slot divergence are not recomputed. The campaign tests the Fourier
leftover-null gap and the n=128 whitened-image shortfall against the
geometry-remainder survival threshold. Neither is a class-wide certificate.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['VECLIB_MAXIMUM_THREADS'] = '1'
import argparse
from hashlib import sha256
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT/'results/development/nsc-ks-certificate-tool-inventory.json'
OWNERS = (
    'scripts/derive_nsc_ks_certificate_tool_inventory.py',
    'docs/nsc-ks-certificate-tool-inventory.md',
)
REGISTERS = {
    'fourier_leftover_null': ROOT/'results/development/nsc-ks-fourier-leftover-null.json',
    'leftover_identification': ROOT/'results/development/nsc-ks-n16-leftover-identification.json',
    'freed_jet': ROOT/'results/development/nsc-ks-n16-freed-jet-proposal.json',
    'n32_geometry_leftover': ROOT/'results/development/nsc-ks-n32-geometry-leftover.json',
    'fourier_wu_leftover': ROOT/'results/development/nsc-ks-fourier-wu-leftover.json',
    'assembled_n32': ROOT/'results/development/nsc-ks-assembled-chebyshev-n32-stall.json',
    'cokernel': ROOT/'results/development/nsc-ks-declared-wu-cokernel.json',
    'image': ROOT/'results/development/nsc-ks-n128-image-channels.json',
    'class_currents': ROOT/'results/development/nsc-ks-class-current-neighborhood.json',
}
EXISTENCE = 3e-11


def digest(path):
    path = Path(path)
    return sha256((path if path.is_absolute() else ROOT/path).read_bytes()).hexdigest()


def load(path):
    return json.loads(path.read_text())


def pair_above(values, level):
    values = [float(item) for item in values]
    return bool(values[0] > level and values[1] > level)


def compute():
    saved = {name: load(path) for name, path in REGISTERS.items()}
    for name, record in saved.items():
        if record.get('physical_NONEXISTENCE_certificate') or record.get('physical_EXISTENCE_certificate'):
            raise ValueError(name + ' already carries a physical gate certificate')
    currents = saved['class_currents']
    image = saved['image']
    null = saved['fourier_leftover_null']
    if currents['search_result'] != 'class_current_neighborhood_not_sign_stable':
        raise ValueError('spent-current search result changed')
    if image['stall'] != 'n128_sup_norm_image_outside_existence_ball':
        raise ValueError('whitened-image stall changed')
    if image['joint_linear_gain_positive'] is not True:
        raise ValueError('whitened image no longer has positive joint gain')
    survival = float(currents['integral_survival_threshold'])
    if not survival > EXISTENCE:
        raise ValueError('survival threshold must sit above the existence ball')
    null_gap = float(null['neighborhood_min'])
    geometry = float(null['pieces']['geometry'])
    assembled = float(null['pieces']['assembled'])
    edge = float(null['pieces']['edge'])
    identification = saved['leftover_identification']
    fourier = saved['fourier_wu_leftover']
    n32 = saved['n32_geometry_leftover']
    assembled_n32 = saved['assembled_n32']
    freed = saved['freed_jet']
    short = image['image_short_of_existence_ball']
    joint = float(image['joint_n_gain_at_beta_drop_1e-8'][-1]['n_gain'])
    leftover_null_test = {
        'neighborhood_min': null_gap,
        'neighborhood_max': float(null['neighborhood_max']),
        'walk_min': float(null['walk_min']),
        'walk_max': float(null['walk_max']),
        'sign_stable_on_its_probe': bool(
            null['walk_sign_stable'] and null['neighborhood_sign_stable']),
        'survives_geometry_remainder_threshold': bool(null_gap > survival),
        'geometry_piece': geometry,
        'assembled_piece': assembled,
        'edge_piece': edge,
        'geometry_abs_exceeds_assembled': bool(abs(geometry) > abs(assembled)),
        'finite_span_is_not_the_class': bool(
            null['scope']['finite_span_is_not_the_class']),
        'unenclosed_edge_tail': null['unenclosed_edge_tail'],
        'edge_can_eat_gap': bool(null['edge_can_eat_gap']),
        'class_certificate': False,
    }
    if not leftover_null_test['sign_stable_on_its_probe']:
        raise ValueError('Fourier leftover-null lost the recorded sign stability')
    if not leftover_null_test['survives_geometry_remainder_threshold']:
        raise ValueError('Fourier leftover-null gap fell through the survival threshold')
    if not leftover_null_test['geometry_abs_exceeds_assembled']:
        raise ValueError('Fourier geometry piece no longer cancels the assembled current')
    if leftover_null_test['unenclosed_edge_tail'] is not None:
        raise ValueError('edge tail became enclosed; this inventory must be rewritten')
    image_test = {
        'short_N': float(short['N']),
        'short_beta': float(short['beta']),
        'joint_n_gain_at_beta_drop_1e-8': joint,
        'joint_linear_gain_positive': True,
        'continuous_cokernel_trivial': bool(
            image['cokernel']['continuous_cokernel_trivial']),
        'both_components_short_of_ball': bool(
            short['N'] > 0.0 and short['beta'] > 0.0),
        'class_certificate': False,
    }
    if not image_test['both_components_short_of_ball'] or not image_test['continuous_cokernel_trivial']:
        raise ValueError('whitened-image obstruction changed')
    if joint <= 0.0:
        raise ValueError('joint gain is no longer positive; refuse a silent certificate')
    inventory = [
        {
            'tool': 'fourier_leftover_null',
            'certifies': (
                'Unit left-null of the iterate4 Fourier geometry Jacobian. '
                'Sign-stable on that walk and Fourier neighborhood. Finite span, '
                'geometry piece larger than the assembled current, edge tail unenclosed.'),
            'class_certificate': False,
        },
        {
            'tool': 'n16_leftover_direction',
            'certifies': (
                'Projection of the five-walk residuals on the n=16 leftover direction. '
                'The walk projection changes sign.'),
            'walk_projection_min': float(identification['walk_projection_min']),
            'walk_projection_max': float(identification['walk_projection_max']),
            'class_certificate': False,
        },
        {
            'tool': 'n16_freed_jet',
            'certifies': (
                'Dropping the solver jet kills the solve-node leftover. The all-node '
                'leftover stays above 3e-11.'),
            'unclipped_all_node_leftover': identification['dropped_jet_unclipped_all_node_leftover'],
            'reaches_existence': bool(identification['dropped_jet_unclipped_reaches_existence']),
            'class_certificate': False,
        },
        {
            'tool': 'n32_geometry_leftover',
            'certifies': (
                'Extra Chebyshev geometry columns at iterate4. The all-node leftover '
                'equals the measured residual and stays above 3e-11. High-mode matter '
                'is not invented.'),
            'all_node_unclipped_leftover': n32['all_node_unclipped_leftover'],
            'above_existence': pair_above(n32['all_node_unclipped_leftover'], EXISTENCE),
            'class_certificate': False,
        },
        {
            'tool': 'fourier_wu_leftover',
            'certifies': (
                'Dirac-free Fourier-on-I geometry image at iterate4. The unclipped '
                'all-node leftover stays above 3e-11.'),
            'all_node_unclipped_leftover': fourier['all_node_unclipped_leftover'],
            'above_existence': pair_above(fourier['all_node_unclipped_leftover'], EXISTENCE),
            'class_certificate': False,
        },
        {
            'tool': 'assembled_chebyshev_n32',
            'certifies': (
                'Assembled n=32 image with reused matter columns. The linear leftover '
                'equals the measured residual and stays above 3e-11.'),
            'assembled_unclipped_leftover': assembled_n32['assembled_unclipped_leftover'],
            'above_existence': pair_above(
                assembled_n32['assembled_unclipped_leftover'], EXISTENCE),
            'class_certificate': False,
        },
        {
            'tool': 'declared_wu_cokernel',
            'certifies': (
                'Principal determinant of compact (w,U) geometry excludes zero. '
                'That is invertibility, not an obstruction.'),
            'status': saved['cokernel']['status'],
            'class_certificate': False,
        },
        {
            'tool': 'n128_whitened_image',
            'certifies': (
                'Stable-rank sup-norm image at ad759424 misses the 3e-11 ball. '
                'Joint linear gain stays positive, and the cokernel is trivial. '
                'One truncation is not every declared history.'),
            'class_certificate': False,
        },
        {
            'tool': 'deciding_error_channels',
            'certifies': (
                'UV tail, full between-node remainder, field error and low/subgap '
                'stay None at ad759424. Missing bounds block a certificate.'),
            'missing': sorted(image['channels']['missing_paths']),
            'class_certificate': False,
        },
        {
            'tool': 'frozen_C0_P_K',
            'certifies': (
                'Frozen incoming C0 has P_K < -0.022178625. That exclusion is a '
                'different history class from evolved C_Sigma[g].'),
            'applies_to_evolved_class': False,
            'class_certificate': False,
        },
        {
            'tool': 'n16_freed_jet_owner_stall',
            'certifies': freed['named_search_stall'],
            'class_certificate': False,
        },
    ]
    if any(item['class_certificate'] for item in inventory):
        raise ValueError('inventory item became a class certificate')
    if freed['evolve_authorized'] or fourier['evolve_authorized'] or n32['evolve_authorized']:
        raise ValueError('an old leftover owner now authorizes evolution')
    certificate = False
    return {
        'schema': 'NSC-KS-CERTIFICATE-TOOL-INVENTORY-v1',
        'accountable_author': 'Douglas Ek',
        'status': (
            'OPEN: declared-class certificate tools are exhausted; '
            'no scoped NON-EXISTENCE'),
        'stall': 'declared_class_certificate_tools_exhausted',
        'prior_stall': image['stall'],
        'prior_search_result': currents['search_result'],
        'profile_identity': image['profile_identity'],
        'lifted_profile_identity': image['lifted_profile_identity'],
        'best_constraint_maxima': image['best_constraint_maxima'],
        'existence_tolerance': EXISTENCE,
        'survival_threshold': survival,
        'families_evolved': 0,
        'probe_authorized': False,
        'beta_primary_evolved': False,
        'higher_n_authorized': False,
        'spent_currents_not_rescored': True,
        'physical_EXISTENCE_certificate': False,
        'physical_NONEXISTENCE_certificate': certificate,
        'campaign': {
            'candidate': 'fourier_leftover_null_and_n128_whitened_image',
            'leftover_null': leftover_null_test,
            'whitened_image': image_test,
            'result': 'both_fail_as_class_certificates',
        },
        'inventory': inventory,
        'missing_primitive': (
            'A necessary relation on the evolved state, or a residual-image '
            'lower bound, that holds for every history in the declared (w,U) '
            'class on I and stays strictly above the geometry-remainder survival '
            'threshold after the changed-history edge and UV tail is enclosed. '
            'Finite-span nulls, a trivial principal cokernel, unenclosed tails, '
            'and one-history whitened images do not supply it.'),
        'next_owner_action': (
            'Do not evolve the beta-primary ray, do not open n=256, and do not '
            're-score L(E), Pi-prime plus E_beta, the integrated defect, or the '
            'slot divergence. The missing object is the class-wide enclosed '
            'relation named in missing_primitive. It is not among these tools.'),
        'source_hashes': {path: digest(path) for path in OWNERS},
        'input_hashes': {
            str(path.relative_to(ROOT)): digest(path) for path in REGISTERS.values()},
    }


def check():
    saved = json.loads(OUTPUT.read_text())
    fresh = json.loads(json.dumps(compute(), sort_keys=True))
    saved_norm = json.loads(json.dumps(saved, sort_keys=True))
    if saved_norm != fresh:
        raise ValueError('certificate-tool inventory changed')
    if saved['physical_NONEXISTENCE_certificate'] or saved['families_evolved']:
        raise ValueError('exhausted-tool register must stay OPEN and unevolved')
    if saved['higher_n_authorized'] or saved['beta_primary_evolved']:
        raise ValueError('n=256 and the beta-primary ray stay unauthorized')
    if saved['stall'] != 'declared_class_certificate_tools_exhausted':
        raise ValueError('exhaustion stall changed')
    return saved


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.write:
        if OUTPUT.exists():
            raise SystemExit('certificate-tool inventory already exists')
        payload = compute()
        OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True)+'\n')
        print(json.dumps({
            'stall': payload['stall'],
            'campaign_result': payload['campaign']['result'],
            'leftover_null_gap': payload['campaign']['leftover_null']['neighborhood_min'],
            'geometry_piece': payload['campaign']['leftover_null']['geometry_piece'],
            'image_short': [
                payload['campaign']['whitened_image']['short_N'],
                payload['campaign']['whitened_image']['short_beta'],
            ],
            'joint_n_gain': payload['campaign']['whitened_image']['joint_n_gain_at_beta_drop_1e-8'],
            'physical_NONEXISTENCE_certificate': payload['physical_NONEXISTENCE_certificate'],
        }, indent=2))
    else:
        saved = check()
        print(json.dumps({
            'stall': saved['stall'],
            'campaign_result': saved['campaign']['result'],
            'physical_NONEXISTENCE_certificate': saved['physical_NONEXISTENCE_certificate'],
        }, indent=2))
