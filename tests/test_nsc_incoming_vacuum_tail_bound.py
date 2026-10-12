"""Directed interval enclosure and the owned scalar/projector identity."""
import json
from copy import deepcopy
from pathlib import Path

import mpmath as mp
import numpy as np
import pytest

from recursive_horizons.nsc_incoming_vacuum_tail_bound import (
    IncomingVacuumTailBoundGeometry, recurrence_cancellation_identity,
    _precision, _lo, _hi, _A, _Aprime, _defect_coefficients,
    vacuum_tail_from_coefficients, aggregate_vacuum_tail_bounds,
)
from recursive_horizons.nsc_pg_high_energy import riccati_coefficients
from recursive_horizons.nsc_transmitting_dirac_domain import S1, S2, S3


ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope='module')
def owned():
    channels = json.loads((ROOT/'results/development/nsc-mode-resolved-cauchy-state.json').read_text())['channels']
    config = json.loads((ROOT/'results/development/nsc-compact-matched-restart.json').read_text())['scattering_provenance']['config']
    owner = IncomingVacuumTailBoundGeometry(config, intervals=8, precision=40)
    return owner, channels, config


def test_worst_low_split_group_has_actual_directed_upper_bound(owned):
    owner, channels, config = owned
    result = owner.bound_channel(channels[22], 160.)
    assert result['status'] == 'PASS'
    assert max(result['vacuum_tail_error_upper']) < 3e-12
    assert result['vacuum_tail_error_upper'][2] == 0.
    assert result['lower_order_cancellation'] == ['0']*15
    assert result['initial_projector_difference'] == 0.
    assert 'finite-delta input artifacts unchanged' in result['initial_scope']
    with _precision(40):
        left, right = map(mp.mpf, result['horizon_enclosure'])
        assert _hi(_A(mp.iv.mpf(left))) < 0
        assert _lo(_A(mp.iv.mpf(right))) > 0
        assert _lo(owner.horizon_derivative) > 0
        kl, ku = map(mp.mpf, result['kappa_rounding_enclosure'])
        assert kl <= mp.mpf(config['surface_gravity']) <= ku
    print(json.dumps({'group': 22, 'intervals': 8,
                      'vacuum_tail_error_upper': result['vacuum_tail_error_upper'],
                      'horizon_enclosure': result['horizon_enclosure'],
                      'kappa_rounding_enclosure': result['kappa_rounding_enclosure'],
                      'endpoint_weight_upper': result['endpoint_weight_upper']}))


def test_scalar_Riccati_defect_equals_normalized_projector_defect():
    # This finite scalar/Pauli identity fixes the factor and orientation used
    # by the enclosure; the values are algebraic controls, not a state/history.
    rho, E, mass, ell = 1.3, 2.4, np.pi/2, np.sqrt(5)
    S, derivative = .3+.2j, .07-.04j
    A = 1+3*rho+3*(1+rho*rho)*(np.arctan(rho)-np.pi/2)
    Aprime = 6+6*rho*(np.arctan(rho)-np.pi/2)
    a, radius = np.sqrt(-A), np.sqrt(1+rho*rho)
    aprime = -Aprime/(2*a)
    ratio = -1j*a*S
    dratio = -1j*(aprime*S+a*derivative)
    v, dv = np.array([1., ratio]), np.array([0., dratio])
    denominator = 1+abs(ratio)**2
    P = np.outer(v, v.conj())/denominator
    dP = (np.outer(dv, v.conj())+np.outer(v, dv.conj()))/denominator
    dP -= np.outer(v, v.conj())*2*np.real(ratio.conjugate()*dratio)/denominator**2
    K = mass/a*S1-ell/(radius*a)*S2+E/a**2*S3
    matrix_defect = np.linalg.norm(dP+1j*(K@P-P@K), 2)
    U = ell/radius+1j*mass
    polynomial = 2*E*S-U-1j*(A*derivative+Aprime*S/2)-A*U.conjugate()*S*S
    scalar_defect = abs(polynomial)/(a*(1+a*a*abs(S)**2))
    assert abs(matrix_defect-scalar_defect) < 3e-14
    assert recurrence_cancellation_identity() == ['0']*15


def test_interval_residual_coefficients_enclose_existing_point_recursion(owned):
    owner, channels, _ = owned
    channel = channels[22]
    # A sanity check of coefficient ownership, not the proof of enclosure:
    # directed arithmetic over each whole cell supplies that enclosure.
    with _precision(40):
        rho, A, sphere = owner.cells[4]
        mass = mp.iv.mpf(channel['compact_mass'])
        angular = mp.iv.mpf(channel['angular_eigenvalue'])
        bounds = _defect_coefficients(A, sphere, mass, angular)
        point = float((_lo(rho)+_hi(rho))/2)
        values, derivatives = riccati_coefficients([point], channel['compact_mass'], channel['angular_eigenvalue'], 16)
        values, derivatives = values[0], derivatives[0]
        metric = 1+3*point+3*(1+point*point)*(np.arctan(point)-np.pi/2)
        ap = 6+6*point*(np.arctan(point)-np.pi/2)
        Uminus = channel['angular_eigenvalue']/np.sqrt(1+point*point)-1j*channel['compact_mass']
        for power, bound in zip(range(16, 33), bounds):
            coefficient = -metric*Uminus*sum(values[j-1]*values[power-j-1]
                for j in range(1, 17) if 1 <= power-j <= 16)
            if power == 16:
                coefficient -= 1j*(metric*derivatives[-1]+ap*values[-1]/2)
            assert _lo(bound.real) <= coefficient.real <= _hi(bound.real)
            assert _lo(bound.imag) <= coefficient.imag <= _hi(bound.imag)


def test_no_new_split_or_unowned_partition(owned):
    owner, channels, config = owned
    with pytest.raises(ValueError, match='archived energy split'):
        owner.bound_channel(channels[22], 300.)
    with pytest.raises(ValueError, match='bounded partition'):
        IncomingVacuumTailBoundGeometry(config, intervals=1000)


def test_directed_cells_cover_entire_collar_without_gaps(owned):
    owner, _, _ = owned
    with _precision(40):
        assert _hi(owner.cells[0][0]) >= _hi(owner.horizon)
        assert _lo(owner.cells[-1][0]) <= 1
        for previous, following in zip(owner.cells, owner.cells[1:]):
            assert _hi(following[0]) >= _lo(previous[0])


def test_all_retained_bound_replay_and_no_finite_start_claim(owned):
    _, channels, _ = owned
    record = json.loads((ROOT/'results/development/nsc-incoming-vacuum-tail-bound.json').read_text())
    groups = record['groups']
    assert len(groups) == 32
    for row in groups:
        replay = vacuum_tail_from_coefficients(channels[row['group']], row['lower'], row['per_sign'])
        assert replay == row['vacuum_tail_error_upper']
        assert max(replay) <= record['tolerances']['per_group_component_upper']
        assert replay[2] == 0.
    result = aggregate_vacuum_tail_bounds(groups)
    assert result == record['aggregate']
    assert max(result['vacuum_tail_error_upper']) <= result['component_budget']
    assert result['finite_offset_initial_projector_term_bound'] is None
    assert not result['finite_offset_archive_accuracy_certified']
    assert record['remaining_initial_term']['value_or_bound'] is None
    assert not record['scope']['full_source_convergence_claimed']
    assert not record['scope']['changed_normal_jets_covered']
    assert not record['scope']['Gamma_rest_assigned']
    with pytest.raises(ValueError, match='all32'):
        aggregate_vacuum_tail_bounds(groups[:-1])
    wrong = deepcopy(groups); wrong[0]['finite_offset_initial_projector_term_bound'] = 0.
    with pytest.raises(ValueError, match='finite-offset'):
        aggregate_vacuum_tail_bounds(wrong)
    wrong = deepcopy(groups[21]['per_sign']); wrong.append(wrong[0])
    with pytest.raises(ValueError, match='angular sign'):
        vacuum_tail_from_coefficients(channels[22], 160., wrong)
