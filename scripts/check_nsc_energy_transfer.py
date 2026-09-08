#!/usr/bin/env python3
"""Compute state-defined Dirac energy transfer on the smooth compact throat."""
import argparse
import hashlib
import inspect
import json
from pathlib import Path
import sys

import numpy as np
from scipy.integrate import simpson
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import expm_multiply

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from recursive_horizons.nsc_energy_transfer import (
    dirac_hamiltonian, positions, regional_operators, split_cut_currents,
    vacuum, expectation, positive_packet, evolve_packet, packet_expectations,
    integrated_current, schur_state_data,
)
from recursive_horizons.nsc_shape_response import smooth_metric
from check_nsc_scale_closure import compare
OUTPUT = ROOT/'results/nsc-6-energy-transfer.json'


def calculate():
    stationary = []
    for eta in (0., .5):
        for kappa in (1, 2, 4):
            metric = smooth_metric(2., 64)
            h = dirac_hamiltonian(metric, kappa, eta)
            e, v, c = vacuum(h)
            ops = regional_operators(h, positions(metric) >= 0)
            cuts = split_cut_currents(ops['current_a'], positions(metric), metric.length)
            # Generic static spectral weighting is a stationarity control, not
            # a definition of the covariantly renormalized vacuum stress.
            weighted = (v * ((e < 0)*np.exp(-e*e/16))[None, :]) @ v.T.conj()
            current = expectation(c, ops['current_a'])
            assert abs(current) < 1e-11
            stationary.append({'eta': eta, 'kappa': kappa,
                'occupied_modes': int(np.sum(e < 0)),
                'commutator_residual': float(np.max(np.abs(h@c-c@h))),
                'net_energy_transfer': current,
                'cut_energy_transfers': {key: expectation(c, op) for key, op in cuts.items()},
                'weighted_spectral_control_transfer': expectation(weighted, ops['current_a']),
                'response': schur_state_data(h, positions(metric) >= 0, .3+.4j, c)})

    times = np.linspace(0, 1.2, 49)
    rows = []
    for points in (32, 64, 128, 256, 512):
        metric = smooth_metric(2., points)
        x = positions(metric)
        h = dirac_hamiltonian(metric)
        e, v, c = vacuum(h)
        ops = regional_operators(h, x >= 0)
        cuts = split_cut_currents(ops['current_a'], x, metric.length)
        psi = positive_packet(metric, e, v)
        states = evolve_packet(e, v, psi, times)
        history = {name: packet_expectations(states, op).tolist() for name, op in ops.items()}
        history.update({name+'_current': packet_expectations(states, op).tolist() for name, op in cuts.items()})
        energy = packet_expectations(states, h)
        integrated = {name: integrated_current(e, v, psi, op, times[-1]) for name, op in cuts.items()}
        change = history['energy_a'][-1] - history['energy_a'][0]
        closure = abs(change - sum(integrated.values()))
        assert closure < 1e-10
        assert np.ptp(energy) < 1e-10
        row = {'points': points, 'packet_energy': float(energy[0]),
               'region_energy_change': change, 'integrated_cut_transfers': integrated,
               'integrated_balance_residual': closure,
               'total_energy_drift': float(np.ptp(energy)),
               'norm_drift': float(np.max(np.abs(np.sum(abs(states)**2, axis=0)-1))),
               'vacuum_overlap': float(np.linalg.norm(c@psi)),
               'times': times.tolist(), 'excitation_history': history}
        if points == 64:
            independent = expm_multiply(-1j*csr_matrix(h), psi, start=0, stop=times[-1], num=len(times)).T
            row['sparse_exponential_state_error'] = float(np.max(np.abs(independent-states)))
            assert row['sparse_exponential_state_error'] < 1e-11
            dense_times = np.linspace(0, times[-1], 769)
            dense_states = evolve_packet(e, v, psi, dense_times)
            quadrature = simpson(packet_expectations(dense_states, ops['current_a']), x=dense_times)
            row['independent_quadrature_balance_residual'] = float(abs(quadrature-change))
            assert row['independent_quadrature_balance_residual'] < 1e-7
            # The same retarded Schur response has zero vacuum flow and
            # nonzero excited-state flow; the response alone cannot choose Q.
            excited = c + np.outer(states[:, 24], states[:, 24].conj())
            row['excited_state_response'] = schur_state_data(h, x>=0, .3+.4j, excited)
            row['excited_state_transfer_at_t_0_6'] = expectation(excited, ops['current_a'])
        rows.append(row)
    refinements = [{'coarse_points': a['points'], 'fine_points': b['points'],
                    'energy_difference': abs(a['packet_energy']-b['packet_energy']),
                    'transferred_energy_difference': abs(a['region_energy_change']-b['region_energy_change'])}
                   for a,b in zip(rows,rows[1:])]
    assert refinements[-1]['energy_difference'] < refinements[-2]['energy_difference']
    return {'schema': 'nsc-energy-transfer-v1', 'artifact_id': 'NSC-6-ENERGY-TRANSFER',
        'classification': 'exact_static_equilibrium_zero_transfer_and_finite_energy_excitation_transport_on_smooth_compact_geometry',
        'source_hashes': {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in (
            'src/recursive_horizons/nsc_energy_transfer.py', 'src/recursive_horizons/nsc_shape_response.py',
            'scripts/check_nsc_energy_transfer.py')},
        'comparison_function_sha256': hashlib.sha256(inspect.getsource(compare).encode()).hexdigest(),
        'input_hashes': {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in (
            'results/nsc-4-smooth-geometry.json', 'results/nsc-4-shape-response.json',
            'results/nsc-5-plateau-conditions.json')},
        'conventions': {'metric': 'static smooth compact axial S1, a=1,R=2,N=q=1',
            'operator': 'existing node/edge H_kappa, no added scalar field or free link',
            'units': 'hbar=c=a_throat=1; energies scale as 1/a, integrated transfer as 1/a, rates as 1/a^2',
            'region': 'x>=0, with throat and compact seam cuts; lattice cut approaches x=0 as h->0',
            'covariance': 'C_ij=<c_j^dagger c_i>; ground state C=theta(-H)',
            'energy_split': 'h_A={P_A,H}/2, h_B=H-h_A; half interaction energy per region',
            'angular_counting': 'one representative 2-component radial channel per kappa; full vacuum multiplier 4*kappa still gives zero; packet excites one mode',
            'packet': 'positive-energy projection of C-infinity bump centered -0.6, width 0.35, momentum 4, +sigma2 spinor; AP topology',
            'packet_origin': 'prepared excitation, not created by a computed collapse or vacuum process',
            'local_energy_caution': 'regional Dirac energy is not a probability and can be negative even in a positive-energy state'},
        'exact_identities': {'covariance_evolution': 'dot C=-i[H,C]',
            'energy_balance': 'dot E_A=Tr(C i[H,h_A])+Tr(C dot h_A)',
            'static_current': 'J_A=i[H,h_A]=i[H^2,P_A]/2; J_A+J_B=0',
            'equilibrium': '[H,C]=0 implies Tr(C J_A)=0 for every fixed regional energy split',
            'finite_system_time_average': '(E_A(T)-E_A(0))/T -> 0 for fixed finite H and state; no uniform infinite-volume/time limit asserted',
            'state_memory': 'eliminating B yields retarded memory AND B exp(-i H_BB t) psi_B(0); its covariance and initial cross-correlations are separate data',
            'plateau_requirement': 'Q/(H rho_b)=1-2 q_dec at a constant f<1; acceleration requires Q/(H rho_b)>1',
            'present_value_control': 'if f=0.95,q_dec=-0.55 were already a plateau, epsilon_required=0.105, Q/(H rho_b)=2.1; these are benchmark inputs, not predictions'},
        'stationary_vacuum': stationary, 'prepared_excitation': rows, 'spatial_refinement': refinements,
        'unresolved': ['absolute renormalized stress on varying throat', 'time-dependent metric and its work/backreaction',
            'physical parent/child state or in-state and causal continuation', 'baryon/dark subsystem definition',
            'coarse-graining to proper-volume cosmological Q and H', 'stationary scale and inheritance ratio'],
        'nonclaims': {'cosmological_Q_derived': False, 'accelerating_plateau_sourced': False,
            'negative_null_stress_implies_energy_flow': False, 'all_time_dependent_vacua_have_zero_flux': False,
            'full_physical_closure': False, 'new_to_world_theorem': False}, 'terminal': True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--check', action='store_true')
    group.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.output and args.output.exists():
        raise FileExistsError('refusing to overwrite result')
    record = calculate()
    if args.check:
        compare(json.loads(OUTPUT.read_text()), record)
        print('Static vacuum zero transfer and excited Dirac energy transport reproduced; every field checked.')
    elif args.output:
        with args.output.open('x') as stream:
            json.dump(record, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write('\n')
        print(args.output)
    else:
        print(json.dumps(record, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
