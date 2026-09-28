#!/usr/bin/env python3
"""Persist the finite incoming preparation theorem; exact algebra only."""
import argparse
from fractions import Fraction
from hashlib import sha256
import json
import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from recursive_horizons.nsc_incoming_cauchy_jets import NORMAL_ENTRIES
from recursive_horizons.nsc_incoming_finite_preparation import symbolic_proof
from recursive_horizons.nsc_spatial_reference_symbol import SymbolJet

OUTPUT = 'results/development/nsc-incoming-finite-preparation.json'
SOURCES = ('src/recursive_horizons/nsc_incoming_finite_preparation.py',
           'scripts/derive_nsc_incoming_finite_preparation.py',
           'tests/test_nsc_incoming_finite_preparation.py',
           'docs/nsc-incoming-finite-preparation.md')
RECORDS = ('results/development/nsc-incoming-constraint-gate.json',
           'results/development/nsc-incoming-jet-rigidity.json',
           'results/development/nsc-incoming-fixed-transfer.json',
           'results/development/nsc-incoming-vacuum-tail-bound.json')
OWNERS = tuple('src/recursive_horizons/'+name+'.py' for name in (
    'nsc_incoming_cauchy_jets', 'nsc_spatial_reference_symbol',
    'nsc_paired_horizon_preparation', 'nsc_pg_massive_modes', 'nsc_pg_high_energy',
    'nsc_incoming_vacuum_tail_bound', 'nsc_incoming_constraint_gate',
    'nsc_incoming_joint_constraints', 'nsc_incoming_local_constraints',
    'nsc_pg_ks_metric_pullback', 'nsc_ks_spacetime_variation', 'nsc_lorentzian',
    'nsc_transmitting_dirac_domain'))
MARGIN = Fraction(5887634380411363, 265464355000000000)


def digest(path):
    return sha256((ROOT/path).read_bytes()).hexdigest()


def authenticate(record, label):
    for section in ('source_hashes', 'input_hashes'):
        for path, expected in record.get(section, {}).items():
            if digest(path) != expected:
                raise ValueError(label+' '+section+' changed: '+path)
    payloads = list(record.get('reused_artifacts', []))
    if 'payload' in record:
        payloads.append(record['payload'])
    for item in payloads:
        if digest(item['path']) != item['sha256']:
            raise ValueError(label+' artifact changed: '+item['path'])


def inputs():
    gate, jet, transfer, tail = [json.loads((ROOT/path).read_text()) for path in RECORDS]
    authenticate(gate, 'constraint-gate')
    authenticate(jet, 'jet-rigidity')
    authenticate(transfer, 'fixed-transfer')
    authenticate(tail, 'vacuum-tail-bound')
    if not gate['status'].startswith('NON-EXISTENCE PASS'):
        raise ValueError('authenticated raw-current witness required')
    witness = gate['source_bound']['rational_witness']
    margin = Fraction(witness['strict_negative_power_margin'])
    if not witness['strict_sign'] or margin != MARGIN or margin <= 0:
        raise ValueError('exact constraint-gate Killing-power margin required')
    if (not jet['status'].startswith('PASS:') or jet['result']['slot_count'] != 20
            or jet['residuals']['exact_symbolic_identities'] != 0
            or jet['result']['sufficient_for_full_matching']
            or jet['result']['finite_energy_threshold'] is not None):
        raise ValueError('completed twenty-slot rigidity receipt required')
    if (not transfer['status'].startswith('PASS:')
            or transfer['residuals']['exact_symbolic_identities'] != 0
            or transfer['scope']['field_radial_source_old_producer_or_energy_scan_runs'] != 0):
        raise ValueError('completed physical fixed-transfer receipt required')
    group = next(row for row in tail['groups'] if row['group'] == 14)
    if (group['initial_projector_difference'] != 0. or group['finite_offset_archive_accuracy_certified']
            or group['finite_offset_initial_projector_term_bound'] is not None
            or group['initial_scope'] != 'declared affine-horizon limit; finite-delta input artifacts unchanged'):
        raise ValueError('same exact affine-horizon P16 bound required')
    signs = sorted(row['sign'] for row in group['per_sign'])
    if signs != [-1, 1]:
        raise ValueError('both angular P16 receipts required')
    for row in group['per_sign']:
        bounds = row['I16_to_I32_upper']
        if len(bounds) != 17 or not all(math.isfinite(value) and value >= 0 for value in bounds):
            raise ValueError('finite nonnegative full-collar angular bounds required')
    if jet['reused_artifacts'][0]['sha256'] != tail['payload']['sha256']:
        raise ValueError('jet-rigidity and vacuum-tail P16 artifacts disagree')
    channel = jet['reused_physical_evidence']['source_bindings']['channel']
    if channel['index'] != 14 or channel['compact_mass']*channel['angular_eigenvalue'] == 0:
        raise ValueError('authenticated nondegenerate retained block required')
    if len(NORMAL_ENTRIES)*2 != 20 or SymbolJet.constant(1.).physical != 4:
        raise ValueError('owned twenty a/r slots and derivative-order-four SymbolJet required')
    return gate, jet, transfer, tail, group, witness, channel


def make_record():
    gate, jet, transfer, tail, group, witness, channel = inputs()
    proof = symbolic_proof()
    rows = proof['finite_endpoint_maps']
    source = proof['signed_source']
    if (proof['exact_scalar_residual_count'] != 33
            or any(set(row['exact_residuals']) != {'0'} or row['amplitude_derivative_used']
                   or not row['nonlinear_lower_jets_retained'] for row in rows)
            or set(source['exact_residuals']) != {'0'}
            or source['tail_receipt_sign_meaning'] != 'angular sign, not energy sign'):
        raise ArithmeticError('finite preparation identities failed')
    return {
        'schema': 'NSC-INCOMING-FINITE-PREPARATION-v1',
        'accountable_author': 'Douglas Ek',
        'status': ('NON-EXISTENCE PASS: declared incoming finite compact prepared fixed-C0 class; '
                   'extended transmitting stationarity OPEN'),
        'result': {
            'finite_normal_identities': [
                'partial_T^k a_g(z)=partial_T^k a_ref(z)',
                'partial_T^k r_g(z)=partial_T^k r_ref(z)',
            ],
            'normal_orders': [1, 2, 3, 4],
            'normal_mixed_slot_count': 20,
            'finite_differences_not_amplitude_derivatives': True,
            'physical_remainder': 'C_g-P_g : H^s -> H^{s+6}',
            'physical_remainder_constant': None,
            'finite_energy_threshold': None,
            'uniform_pointwise_energy_error_claimed': False,
            'sufficient_for_C0_matching': False,
            'included_shift': 'E_beta[g,C0]=E_beta[g_ref,C0]=P_K<-M',
            'M': str(MARGIN),
            'M_numerator': MARGIN.numerator,
            'M_denominator': MARGIN.denominator,
            'classical_projector_extension': (
                'analytic finite symbol through degree -6; distinct from implemented derivative-order-4 SymbolJet'),
            'implemented_SymbolJet_derivative_order': 4,
            'executed_SymbolJet_at_order_six': False,
        },
        'assumptions': {
            'history_class': (
                'smooth finite real metric histories on a compact KS slab between an unchanged upstream slice '
                'and incoming Sigma, inside the regular trapped chart, compact axial support, N,a,r positive'),
            'incoming_slice': 'spacelike; g^(TT)=1/N^2>0; reference radial normal a0^2=beta_PG^2-1>0',
            'connected_PG_chart': 'F>0, s_K>0 from the metric-pullback owner',
            'intrinsic_and_N_beta': (
                'N,beta,a,r fixed as functions along Sigma; N/beta normal jets through four also fixed'),
            'free_slots': 'the twenty a/r normal/mixed slots of IncomingNormalJetChange',
            'preparation': 'same upstream state and original source law; C_g(Sigma)=C0 exactly',
            'energy_spectrum': (
                'existing unbounded signed-real-energy continuum; no finite physical energy cutoff '
                'is introduced or covered; bounded-frequency partitions are estimate devices only'),
            'amplitudes': 'finite, not infinitesimal; no selected numerical duration',
        },
        'checked_finite_algebra': proof,
        'analytic_physical_bridge': {
            'upstream_P16': (
                'translation-invariant Fourier multiplier O(E^-16) from the positive full-collar bound; '
                'Hs to Hs+16. No energy-differentiability of that bound is used.'),
            'P16_per_sign': 'angular',
            'negative_energy_cone': (
                'charge conjugation of the existing source/mode owner: '
                'C(-E,ell)=I-sigma3*conj(C(E,-ell))*sigma3; not the tail per_sign label'),
            'finite_classical_projector': (
                'existence through degree -6 using owned diagonal/off-diagonal compatibility; '
                'R_g in Psi^{-6}; not an evaluation of stored SymbolJets'),
            'Duhamel': (
                'C_g(T)-P_g(T)=U D_u U^dagger - int U R_g U^dagger; both terms map Hs to Hs+6'),
            'endpoint_coefficient_bridge': (
                'C_g(Sigma)=C0 implies P_g(Sigma)-P_ref(Sigma) also Hs to Hs+6, so classical '
                'coefficients through degree -5 agree'),
            'remainder_constant': None,
            'uniform_pointwise_energy_error': False,
        },
        'independent_reviews': {
            'finite_jet_triangularity': '69ac44c5-8b55-41bf-ba3c-200822a36b06',
            'physical_projector_remainder_and_spacelike_KS': 'a4589cc3-11b0-4c92-9dea-a7b98d37649b',
        },
        'reused_physical_evidence': {
            'constraint_gate_record': RECORDS[0],
            'jet_rigidity_record': RECORDS[1],
            'fixed_transfer_record': RECORDS[2],
            'vacuum_tail_record': RECORDS[3],
            'raw_current_witness': witness['strict_negative_power_margin'],
            'background': {
                'group': 14,
                'angular_signs': [-1, 1],
                'starting_tail_domain': group['lower'],
                'full_collar_I16_to_I32_upper_by_angular_sign': {
                    str(row['sign']): row['I16_to_I32_upper'] for row in group['per_sign']
                },
            },
            'source_bindings': jet['reused_physical_evidence']['source_bindings'],
            'included_assembly_owner': 'src/recursive_horizons/nsc_incoming_joint_constraints.py',
            'lapse_source_numerical_accuracy_used': False,
            'physical_generators_rerun': False,
        },
        'residuals': {'exact_symbolic_identities': 0},
        'verification_tolerances': {'exact_symbolic_identities': 0},
        'scope': {
            'excluded_class': (
                'smooth finite compact same-upstream/source histories with C_g(Sigma)=C0, '
                'fixed intrinsic N,a,r,beta on Sigma, and fixed N/beta normal jets through four'),
            'finite_compact_prepared_fixed_C0_class': 'NON-EXISTENCE',
            'extended_transmitting_stationarity': 'OPEN',
            'broader_NSC_architecture_excluded': False,
            'other_parent_data_or_evolved_incoming_covariance_excluded': False,
            'higher_normal_jets_classified': False,
            'sufficient_C0_matching': False,
            'uniform_pointwise_energy_error': False,
            'computed_Sobolev_constant': None,
            'finite_energy_threshold': None,
            'classical_projector_is_implemented_SymbolJet': False,
            'P16_per_sign_meaning': 'angular',
            'source_law_or_boundary_prescription_changed': False,
            'finite_physical_energy_cutoff_introduced': False,
            'new_Hadamard_gate_or_cutoff': False,
            'physical_field_source_radial_producer_or_energy_scan_runs': 0,
            'metric_evolution_or_publication': False,
        },
        'source_hashes': {path: digest(path) for path in SOURCES},
        'input_hashes': {path: digest(path) for path in (*RECORDS, *OWNERS)},
        'reused_artifacts': [tail['payload']],
        'reproducer': 'python3 scripts/derive_nsc_incoming_finite_preparation.py --check',
    }


def main():
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    record = make_record()
    path = ROOT/OUTPUT
    if args.write:
        path.write_text(json.dumps(record, indent=2, sort_keys=True, allow_nan=False)+'\n')
    elif json.loads(path.read_text()) != record:
        raise ValueError('finite-preparation proof receipt differs')
    print(json.dumps({
        'status': record['status'],
        'normal_orders': record['result']['normal_orders'],
        'slot_count': record['result']['normal_mixed_slot_count'],
        'M': record['result']['M'],
        'physical_remainder': record['result']['physical_remainder'],
        'exact_scalar_residual_count': record['checked_finite_algebra']['exact_scalar_residual_count'],
        'residuals': record['residuals'],
        'excluded_class': record['scope']['excluded_class'],
        'extended_transmitting_stationarity': record['scope']['extended_transmitting_stationarity'],
    }, indent=2))


if __name__ == '__main__':
    main()
