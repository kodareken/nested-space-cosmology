#!/usr/bin/env python3
"""Define the actual common-time bulk projectors and check the Dirac domain."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from recursive_horizons.nsc_common_time_bulk_split import (
    CommonTimeBulkSplit, BulkOperatorDomainError, require_common_time_bulk_split,
)
from recursive_horizons.nsc_causal_common import LinkHistory
from recursive_horizons.nsc_transmitting_dirac_domain import I2, S2
from derive_nsc_transmitting_boundary_binding import compare

OUTPUT = 'results/development/nsc-common-time-bulk-split.json'
INPUTS = {
    'T': 'results/development/nsc-transmitting-dirac-domain.json',
    'negative': 'results/development/nsc-hamiltonian-trace-representation.json',
    'state': 'results/development/nsc-mode-resolved-cauchy-state.json',
    'dirac': 'results/nsc-4-dirac-tetrad.json',
    'spatial': 'results/nsc-8-chiral-boundary.json',
}
SOURCES = (
    'src/recursive_horizons/nsc_common_time_bulk_split.py',
    'scripts/define_nsc_common_time_bulk_split.py',
    'tests/test_nsc_common_time_bulk_split.py',
    'docs/nsc-common-time-bulk-split.md',
)
DEPENDENCIES = (
    'src/recursive_horizons/nsc_causal_common.py',
    'src/recursive_horizons/nsc_lorentzian.py',
    'src/recursive_horizons/nsc_adm_source.py',
    'src/recursive_horizons/nsc_transmitting_dirac_domain.py',
    'scripts/derive_nsc_transmitting_boundary_binding.py',
)


def sha(path):return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()


def calculate():
    inputs = {key: json.loads((ROOT/path).read_text()) for key, path in INPUTS.items()}
    t, state, negative = (inputs[k] for k in ('T', 'state', 'negative'))
    for key in ('A', 'magnetic_flux', 'Omega', 'zeta', 'V_full'):
        if t['locked_inputs'][key] != state['locked_inputs'][key] or t['locked_inputs'][key] != negative['locked_inputs'][key]:
            raise ValueError('locked coefficient or scale changed')
    payload = t['payload']
    if sha(payload['path']) != payload['sha256']:raise ValueError('T artifact changed')
    with np.load(ROOT/payload['path'], allow_pickle=False) as f:
        weights_T = np.array(f['quadrature_weight'])
        metric_T = np.array(f['trace_metric'])
    domain = require_common_time_bulk_split(CommonTimeBulkSplit())
    radius = inputs['spatial']['conventions']['radius']
    rows = []
    max_residuals = {name: 0. for name in ('idempotence', 'orthogonality', 'completeness', 'canonical_CAR', 'spatial_half_density')}
    # Spatial midpoint grids are rank/norm probes only. Counts reuse the
    # retained resolution budget, NOT a frequency-to-rho basis transformation.
    for channel in state['channels']:
        count = channel['sample_count']
        rho = -radius+(np.arange(count)+0.5)*2*radius/count
        weights = np.full(count, 2*radius/count)
        spheres = np.sqrt(1+rho*rho)
        maps = domain.projectors(rho)
        pp, pc, s = maps['P_parent'], maps['P_child'], maps['canonical_reorder']
        identity = np.eye(2*count)
        inverse_half, gram = domain.spatial_half_density_map(weights, np.ones(count), spheres)
        residuals = {
            'idempotence': float(max(np.max(abs(pp@pp-pp)), np.max(abs(pc@pc-pc)))),
            'orthogonality': float(np.max(abs(pp@pc))),
            'completeness': float(np.max(abs(pp+pc-identity))),
            'canonical_CAR': float(np.max(abs(s@s.T-identity))),
            'spatial_half_density': float(np.max(abs(inverse_half*gram*inverse_half-1))),
        }
        for key, value in residuals.items():max_residuals[key] = max(max_residuals[key], value)
        if len(np.unique(np.argmax(s, axis=1))) != 2*count:raise ArithmeticError('bulk variables duplicated or lost')
        rows.append({'channel': channel['index'], 'compact_level': channel['compact_level'],
                     'angular_level': channel['angular_level'], 'spatial_probe_nodes': count,
                     'parent_rank': int(np.trace(pp)), 'child_rank': int(np.trace(pc)),
                     'independent_two_space_rank': 2*count, 'residuals': residuals,
                     'bounded_local_potential': float(np.hypot(channel['compact_mass'], channel['angular_eigenvalue']))})
    geometry = t['geometry']
    N, q, beta, r = (geometry[k] for k in ('N', 'q_PG', 'beta', 'r'))
    v0 = domain.principal(N, q, beta)
    # Same geometric cut, different integrations: T is the full spacelike
    # rho=0 seam; on each common-time slice its intersection is a two-sphere.
    compatibility = float(np.max(abs(metric_T+q*r*r*weights_T[:, None, None]*v0)))
    orbit = np.array([[N*N-q*q*beta*beta, -q*q*beta], [-q*q*beta, -q*q]])
    time_norm = float(np.linalg.inv(orbit)[0, 0])
    time_residual = abs(time_norm-1/(N*N))
    defect = domain.projection_domain_defect(N=N, q_PG=q, beta=beta)
    coefficient = defect['delta_coefficient']
    delta_antihermitian = float(np.max(abs(coefficient+coefficient.conj().T)))
    expected_min = beta-N/q;expected_max = beta+N/q
    singular_residual = max(abs(defect['smallest_singular_value']-expected_min), abs(defect['largest_singular_value']-expected_max))
    tol = 3e-11
    if max(*max_residuals.values(), compatibility, time_residual, delta_antihermitian, singular_residual) > tol:
        raise ArithmeticError('common-time domain definition check failed')
    if defect['smallest_singular_value'] <= tol:raise ArithmeticError('noncharacteristic seam hypothesis changed')
    rejected = {}
    for label, call in (
        ('given_B', lambda: require_common_time_bulk_split(LinkHistory(np.zeros((1,1,1))))),
        ('trace_graph_as_bulk_projectors', lambda: require_common_time_bulk_split(t)),
        ('ordinary_B_from_domain_unsafe_product', domain.ordinary_hamiltonian_link),
        ('copy_seed_frequency_covariance_to_PG_slice', lambda: domain.attach_seed_trace_covariance(np.eye(2))),
    ):
        try:call()
        except BulkOperatorDomainError as error:rejected[label] = str(error)
        else:raise AssertionError('invalid path accepted: '+label)
    return {
        'schema': 'NSC-COMMON-TIME-BULK-SPLIT-v1',
        'status': 'Z1 DEFINITION PASS: independent common-PG-time bulk spaces; Z2 ordinary link product OPEN for operator-domain reasons',
        'source_hashes': {p:sha(p) for p in SOURCES},
        'input_hashes': {p:sha(p) for p in (*INPUTS.values(), *DEPENDENCIES)},
        'authenticated_T_payload': payload,
        'locked_inputs': t['locked_inputs'],
        'definition': {
            'common_time': 'inherited PG tau; tau=constant is a whole-line Cauchy slice',
            'canonical_field': 'chi=r sqrt(q_PG) psi; L2(d rho) per reduced channel',
            'Hilbert_space': 'directsum of 33 retained channel copies of L2(R_rho;C2), with existing degeneracy/compact labels unchanged',
            'parent': 'L2((0,infinity);C2) in each channel',
            'child': 'L2((-infinity,0);C2) in each channel',
            'P_parent': 'multiplication by 1_(rho>0)',
            'P_child': 'multiplication by 1_(rho<0)',
            'projector_domain': 'entire L2 Hilbert space; value at rho=0 irrelevant to L2 classes',
            'transmitting_domain': 'graph closure of smooth compact-support whole-line spinors under the inherited Dirac operator; near noncharacteristic rho=0, piecewise H1 representatives have equal one-sided traces',
            'Hilbert_split_not_operator_domain_product': True,
            'no_new_outer_boundary_condition': True,
            'finite_probe_is_not_a_reflecting_or_periodic_Hamiltonian': True,
        },
        'operator': {
            'expression': 'H_c[g]=-i(v partial_rho+v_prime/2)+M_c; v=(N/q_PG)sigma2-beta I',
            'retained_local_term': 'M_c=-N*m_c*sigma1-N*lambda_c/r*sigma3 in the current basis fixed by T; a constant current-preserving phase converts the massless term to the old +lambda/r sigma1 convention',
            'realization': 'inherited self-adjoint whole-line closure on the fixed unwarped PG background; bounded retained multiplication terms preserve that realization',
            'imported_completion': inputs['dirac']['complete_line'],
            'maximum_retained_potential_bound': max(row['bounded_local_potential'] for row in rows),
            'off_shell_class': 'smooth compactly supported coefficient variations of the inherited background retaining positive N,q,r and noncharacteristic seam; no geometry solution selected',
            'full_state_on_common_time_slice': None,
            'seed_frequency_nodes_identified_with_spatial_nodes': False,
        },
        'finite_verification': {
            'purpose': 'orthogonal-support/CAR witnesses, not time evolution or seed-state reconstruction',
            'window_radius_reused': radius,
            'spatial_grid': 'midpoints in [-R,R]; zero is excluded; only the per-channel sample count is reused as a numerical budget',
            'rows': rows, 'maximum_residuals': max_residuals, 'tolerance': tol,
            'total_parent_rank': sum(row['parent_rank'] for row in rows),
            'total_child_rank': sum(row['child_rank'] for row in rows),
            'total_independent_rank': sum(row['independent_two_space_rank'] for row in rows),
            'rank_defect': 0,
            'old_duplicated_trace_CAR_mismatch': negative['operator_CAR_certificate']['measured_maximum_spectral_norm_mismatch'],
            'old_negative_certificate_reused_not_rerun': True,
        },
        'T_compatibility': {
            'condition': 'psi_parent(tau,0)=psi_child(tau,0) is an operator-domain trace condition for every tau; not equality of bulk canonical variables',
            'current_relation': 'G_T,j=-q_PG*r^2*w_j*v(0) at the inherited seam; opposite orientations cancel',
            'current_relation_residual': compatibility,
            'g_inverse_tautau': time_norm, 'time_normal_residual': time_residual,
            'T_full_seam_distinct_from_slice_intersection': True,
        },
        'Z2_domain_check': {
            'commutator': defect['identity'],
            'delta_coefficient': [[[float(z.real),float(z.imag)] for z in row] for row in coefficient],
            'coefficient_smallest_singular_value': defect['smallest_singular_value'],
            'coefficient_largest_singular_value': defect['largest_singular_value'],
            'coefficient_singular_value_residual': singular_residual,
            'coefficient_antihermiticity_residual': delta_antihermitian,
            'P_child_preserves_transmitting_operator_domain': False,
            'ordinary_product_P_parent_H_P_child': None,
            'physical_B': None, 'physical_dB_dgDelta': None,
            'zero_trace_core_B_zero_promoted_to_transmitting_link': False,
            'delta_coefficient_renamed_as_B': False,
            'exact_remaining_owner': 'a domain-respecting weak/Galerkin or boundary-resolvent realization of the transmitted Dirac operator before identifying a finite B block',
            'bounded_response_already_defined': 'P_parent (H_D-z)^-1 P_child for Im z nonzero, using the inherited self-adjoint H_D',
            'bounded_response_evaluated': False,
        },
        'moving_surface': {'status': 'OUT OF SCOPE', 'reason': 'this definition keeps the inherited coordinate cut rho=0; it selects no moving/dynamical neck', 'residual': None},
        'downstream': {'Z1': 'PASS', 'Z2': 'OPEN', 'Z3': 'fixed cut; moving jets out of scope', 'Z4': 'OPEN',
                       'Gamma_rest_boundary_derivative': None, 'two_sided_Weyl_mismatch': None,
                       'physical_Vc': None, 'stationary_history': None, 'extended_nonexistence_claimed': False},
        'rejected_paths': rejected,
        'scope': {'new_physical_terms_or_particle_copies': False, 'parameters_refitted': False,
                  'given_B_used': False, 'old_generators_rerun': False, 'optimizer_started': False,
                  'metric_timestep_started': False, 'coupled_evolution_reopened': False,
                  'finite_stress_fabricated': False, 'physical_history_selected': False},
        'comparison': t['comparison'],
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    group = p.add_mutually_exclusive_group();group.add_argument('--check',action='store_true');group.add_argument('--write',action='store_true')
    args=p.parse_args();record=calculate();out=ROOT/OUTPUT
    if args.check:compare(json.loads(out.read_text()),record);print('Common-time bulk split and sharp-projector domain record verified')
    elif args.write:
        if out.exists():raise FileExistsError('refusing overwrite')
        out.write_text(json.dumps(record,indent=2,sort_keys=True,allow_nan=False)+'\n');print(record['status'])
    else:print(json.dumps(record,indent=2,sort_keys=True,allow_nan=False))


if __name__=='__main__':main()
