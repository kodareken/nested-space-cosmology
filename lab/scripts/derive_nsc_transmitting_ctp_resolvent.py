#!/usr/bin/env python3
"""Evaluate only the new NSC CTP memory/frequency-transfer connection."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from recursive_horizons.nsc_transmitting_ctp_resolvent import (
    ArchivedResolventField, CTPRetardedMemory, TransmittingCrossResolvent,
    ctp_completion_inventory, frequency_transfer_driven, frequency_transfer_weak,
)
from recursive_horizons.nsc_transmitting_ctp_variation import endpoint_ctp_pullback
from recursive_horizons.nsc_transmitting_boundary_history import EndpointVariation, match_endpoint_variations
from recursive_horizons.nsc_mode_resolved_cauchy_state import deterministic_npz_bytes
from derive_nsc_transmitting_boundary_binding import compare

OUTPUT = 'results/development/nsc-transmitting-ctp-resolvent.json'
SOURCES = (
    'scripts/derive_nsc_transmitting_ctp_resolvent.py',
    'src/recursive_horizons/nsc_transmitting_ctp_resolvent.py',
    'tests/test_nsc_transmitting_ctp_resolvent.py',
    'docs/nsc-transmitting-ctp-resolvent.md',
)
INPUTS = (
    'results/development/nsc-transmitting-cross-resolvent.json',
    'results/development/nsc-transmitting-boundary-remainder.json',
    'results/development/nsc-mode-resolved-cauchy-state.json',
    'results/development/nsc-general-ks-same-action-history.json',
    'src/recursive_horizons/nsc_transmitting_resolvent.py',
    'src/recursive_horizons/nsc_transmitting_ctp_variation.py',
    'src/recursive_horizons/nsc_transmitting_boundary_history.py',
    'src/recursive_horizons/nsc_boundary_state.py',
    'src/recursive_horizons/nsc_influence.py',
    'src/recursive_horizons/nsc_causal_common.py',
    'src/recursive_horizons/nsc_mode_resolved_cauchy_state.py',
    'scripts/derive_nsc_transmitting_boundary_binding.py',
)


def digest(p): return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def residual(a, reference): return float(np.linalg.norm(a)/max(1., np.linalg.norm(reference)))


def calculate():
    old, remainder, state, homogeneous = [json.loads((ROOT/p).read_text()) for p in INPUTS[:4]]
    payload = old['payload']
    if digest(payload['path']) != payload['sha256']: raise ValueError('locked response artifact changed')
    arrays = np.load(ROOT/payload['path'], allow_pickle=False)
    x, weights = arrays['rho_quadrature'], arrays['rho_weights']
    zi = complex(*old['parameters']['z'])
    zo = complex(old['parameters']['identity_check_energy'][0], zi.imag)
    if not homogeneous['gate']['nonexistence_in_declared_class_proved']:
        raise ValueError('homogeneous certificate changed')
    if old['locked_inputs'] != remainder['locked_inputs']:
        raise ValueError('locked action ledger differs')
    rows = []
    saved = {name: [] for name in (
        'outgoing_response', 'outgoing_advanced_response', 'delta_retarded_response',
        'delta_inverse_retarded_memory', 'outgoing_fields', 'outgoing_adjoint_fields',
        'outgoing_field_derivatives', 'outgoing_adjoint_field_derivatives',
    )}
    for index, item in enumerate(old['channels']):
        channel = item['channel']
        owner = TransmittingCrossResolvent(channel['compact_mass'], channel['angular_eigenvalue'])
        incoming = ArchivedResolventField(owner, zi, x, weights, arrays['upper_fields'][index])
        incoming_adj = ArchivedResolventField(owner, zi.conjugate(), x, weights, arrays['lower_fields'][index])
        outgoing = owner.solve(zo, with_jets=False)
        outgoing_adj = owner.solve(zo.conjugate(), with_jets=False)
        driven, seam = frequency_transfer_driven(owner, zo, incoming)
        weak = frequency_transfer_weak(owner, outgoing_adj, incoming, x, weights)
        advanced, advanced_seam = frequency_transfer_driven(owner, zi.conjugate(), outgoing_adj)
        # Zero-transfer check reconstructs the old static kernel from its
        # stored fields. It does not re-run the old source or static generator.
        diagonal = frequency_transfer_weak(owner, incoming_adj, incoming, x, weights)
        memory_in = CTPRetardedMemory(arrays['upper_response'][index], arrays['lower_response'][index])
        memory_out = CTPRetardedMemory(outgoing.response, outgoing_adj.response)
        delta_d = memory_out.relative_vertex(memory_in, driven)
        reconstruction = -memory_out.G_retarded@delta_d@memory_in.G_retarded
        checks = {
            'outgoing_adjoint': residual(outgoing_adj.response-outgoing.response.conj().T, outgoing.response),
            'driven_vs_weak_two_frequency': residual(driven-weak, driven),
            'advanced_frequency_exchange': residual(advanced-driven.swapaxes(-1,-2).conj(), driven),
            'inverse_memory_Dyson_reconstruction': residual(reconstruction+driven, driven),
            'zero_transfer_static_recovery': residual(diagonal-arrays['upper_metric_jets'][index], arrays['upper_metric_jets'][index]),
            'archived_field_seam': max(incoming.seam_residual, incoming_adj.seam_residual),
            'driven_seam': max(seam, advanced_seam),
            'retarded_reverse_response': float(np.linalg.norm(outgoing.response[:2,2:])),
            'retarded_reverse_metric_transfer': float(np.linalg.norm(driven[:,:2,2:])),
            'retarded_inverse_identity': float(np.linalg.norm(memory_out.G_retarded@memory_out.D_retarded-np.eye(4))),
        }
        rows.append({
            'channel_index': channel['index'], 'family': channel['family'],
            'residuals': checks,
            'two_frequency_metric_response_norm': float(np.linalg.norm(driven)),
            'two_frequency_inverse_memory_vertex_norm': float(np.linalg.norm(delta_d)),
            'frequency_transfer_difference_from_static': float(np.linalg.norm(driven-arrays['upper_metric_jets'][index])),
            'forward_response_norm': float(np.linalg.norm(outgoing.response[2:,:2])),
        })
        for name, value in (
            ('outgoing_response', outgoing.response), ('outgoing_advanced_response', outgoing_adj.response),
            ('delta_retarded_response', -driven), ('delta_inverse_retarded_memory', delta_d),
            ('outgoing_fields', np.array([outgoing.field(t) for t in x])),
            ('outgoing_adjoint_fields', np.array([outgoing_adj.field(t) for t in x])),
            ('outgoing_field_derivatives', np.array([outgoing.field_derivative(t) for t in x])),
            ('outgoing_adjoint_field_derivatives', np.array([outgoing_adj.field_derivative(t) for t in x])),
        ): saved[name].append(value)
    maxima = {key: max(row['residuals'][key] for row in rows) for key in rows[0]['residuals']}
    tolerance = 3e-11
    if max(maxima.values()) > tolerance: raise ArithmeticError(json.dumps(maxima, indent=2))
    if not all(row['frequency_transfer_difference_from_static'] > 1e-10 for row in rows):
        raise ArithmeticError('temporal transfer indistinguishable from the static vertex')
    covector = EndpointVariation(remainder['endpoint_basis']['id'], tuple(
        tuple(remainder['locked_local_components'][key]) for key in ('N','beta','q_ADM','r')))
    c2 = endpoint_ctp_pullback(None, None, None, covector, None)
    c3 = match_endpoint_variations(covector, None)
    saved = {k: np.array(v) for k,v in saved.items()}
    saved.update(rho_quadrature=x, rho_weights=weights, channel_indices=arrays['channel_indices'])
    data = deterministic_npz_bytes(saved); sha = hashlib.sha256(data).hexdigest()
    return {
        'schema': 'NSC-TRANSMITTING-CTP-RESOLVENT-v1',
        'status': 'C1 retarded-memory and harmonic metric-transfer PASS; state-dependent C1 and C2-C4 OPEN',
        'source_hashes': {p: digest(p) for p in SOURCES},
        'input_hashes': {p: digest(p) for p in (*INPUTS, payload['path'])},
        'locked_inputs': old['locked_inputs'],
        'definition': {
            'domain': old['definition']['domain'],
            'projection': 'same canonical packets J as Z2a; fixed rho=0; no new independent room variables',
            'retarded_conversion': 'G_R=-R_packet; D_R=-K; G_A(z*)=G_R(z)^dagger',
            'memory': 'full energy dependence retained in D_R and its cross block; no constant B approximation',
            'CTP_triangular_blocks': 'G=[[G_R,G_K],[0,G_A]]; D=[[D_R,-D_R o G_K o D_A],[0,D_A]]',
            'state_lesser': 'G_less(t,s)=i J^dagger U(t,t0) C0 U(s,t0)^dagger J',
            'state_greater': 'G_greater(t,s)=-i J^dagger U(t,t0)(I-C0)U(s,t0)^dagger J',
            'preparation': 'full PG covariance including packet/complement cross blocks; not inferred from retarded data',
            'existing_Gaussian_owner': 'GaussianBoundaryState.reconstruct and overlap_schur; nsc_influence.influence',
            'full_action': 'Q=I-C0+C0 Vminus^dagger Vplus; retain complementary determinant and all preparation blocks',
            'physical_preparation_available': False,
            'stored_state_basis': state['basis'],
            'stored_seed_covariance_assigned_to_PG_packets': False,
        },
        'temporal_metric_vertex': {
            'input_energy': [zi.real, zi.imag], 'output_energy': [zo.real, zo.imag],
            'omega': float(zo.real-zi.real),
            'frequency_provenance': 'difference of existing real energy ports; shared imaginary part; response probe, not physical frequency selection',
            'fields': ['log_N','beta','log_q_PG','log_r'],
            'source': 'existing compact spatial metric vertex times exp(-i omega t); real source includes conjugate harmonic',
            'delta_R': '-J^dagger R(zout) deltaH R(zin) J',
            'delta_D_R': 'K(zout) delta_R(zout,zin) K(zin)',
            'adjoint': 'delta_R_A(zin*,zout*)=delta_R_R(zout,zin)^dagger',
            'finite_history_duration': None,
            'full_time_history_reconstructed': False,
        },
        'verification': {'maximum_residuals': maxima, 'tolerance': tolerance, 'retained_channels': len(rows),
            'earlier_static_residuals_reused': old['verification']['maximum_residuals']},
        'channels': rows,
        'payload': {'path': f'results/development/artifacts/nsc-transmitting-ctp-resolvent.{sha}.npz',
            'sha256': sha, 'bytes': len(data),
            'description': 'outgoing resolvent fields and new two-frequency retarded/ inverse-memory metric vertices; no state or stress'},
        'completion_inventory': ctp_completion_inventory(),
        'C2_existing_endpoint_consumer': c2,
        'C3_existing_Weyl_match_consumer': c3,
        'locked_endpoint_basis': remainder['endpoint_basis'],
        'locked_eight_Weyl_coefficients': remainder['locked_local_components'],
        'locked_Weyl_extraction_residual': remainder['locked_extraction_residual'],
        'preserved_Weyl_value': remainder['preserved_Weyl_value'],
        'C4': {'decision': 'OPEN', 'reason': 'full PG preparation, physical time-history jets, endpoint pullback and Gamma_rest boundary derivative are unevaluated',
            'homogeneous_certificate_preserved': homogeneous['gate']['nonexistence_in_declared_class_proved'],
            'physical_Vc': None, 'metric_stationarity_residual': None},
        'scope': {'parameters_refitted': False, 'prior_scientific_files_modified': False,
            'old_generators_rerun': False, 'new_physical_terms': False, 'given_B_used': False,
            'state_or_stress_fabricated': False, 'optimizer_started': False, 'metric_timestep_started': False,
            'Z3': 'OUT OF SCOPE: fixed coordinate cut'},
        'comparison': old['comparison'],
    }, data


def main():
    p=argparse.ArgumentParser(description=__doc__)
    g=p.add_mutually_exclusive_group();g.add_argument('--check',action='store_true');g.add_argument('--write',action='store_true')
    args=p.parse_args();record,data=calculate();out=ROOT/OUTPUT;artifact=ROOT/record['payload']['path']
    if args.check:
        compare(json.loads(out.read_text()), record)
        if artifact.read_bytes()!=data: raise AssertionError('CTP response artifact differs')
        print('New CTP memory/temporal vertices verified; state-dependent C1 and C2-C4 OPEN')
    elif args.write:
        if out.exists() or artifact.exists(): raise FileExistsError('refusing immutable record overwrite')
        artifact.write_bytes(data);out.write_text(json.dumps(record,indent=2,sort_keys=True)+'\n')
        print('Wrote',OUTPUT,record['verification'])
    else: print(json.dumps(record,indent=2,sort_keys=True))
    return 0


if __name__=='__main__': raise SystemExit(main())
