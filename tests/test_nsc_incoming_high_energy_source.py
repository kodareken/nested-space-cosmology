"""Stable local observable of archived massive fields; no ODE/tail rerun."""
import hashlib
import json
from pathlib import Path

import mpmath as mp
import numpy as np
import pytest

from recursive_horizons.nsc_incoming_high_energy_source import incoming_high_energy_source
from recursive_horizons.nsc_incoming_state_moments import incoming_adiabatic_reference, AXIAL, RADIUS
from recursive_horizons.nsc_paired_horizon_preparation import source_covariance
from recursive_horizons.nsc_pg_high_energy import riccati_coefficients


ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope='module')
def retained():
    record = json.loads((ROOT/'results/development/nsc-pg-retained-covariance.json').read_text())
    path = ROOT/record['payload']['path']
    assert hashlib.sha256(path.read_bytes()).hexdigest() == record['payload']['sha256']
    with np.load(path, allow_pickle=False) as payload:
        yield payload


def _inputs(payload, panel):
    metadata = json.loads(payload[panel+'/metadata_json'].tobytes())
    all_E = payload[panel+'/energies']
    targets = (40., 100., metadata['right']-1.)
    E = all_E[[np.argmin(abs(all_E-value)) for value in targets]]
    mass = metadata['channel']['compact_mass']
    config = metadata['config']
    C = np.array([source_covariance(e, config['surface_gravity'], config['omega'], mass) for e in E])
    return metadata, E, C


def _high_precision_unrationalized(E, mass, angular, source, higher):
    coefficients, _ = riccati_coefficients([1.], mass, angular, 16)
    answer = []
    with mp.workdps(60):
        a, r, m, ell = map(mp.mpf, (float(AXIAL), float(RADIUS), float(mass), float(angular)))
        for i, energy in enumerate(E):
            e = mp.mpf(float(energy))
            series = sum(mp.mpc(float(c.real), float(c.imag))/(2*e)**(j+1)
                         for j, c in enumerate(coefficients[0]))
            u = a*a*abs(series)**2
            b = [2*a*series.imag/(1+u), -2*a*series.real/(1+u), (1-u)/(1+u)]
            h = [-m, ell/r, -e/a]
            gap = mp.sqrt(sum(v*v for v in h))
            reference = [-v/gap for v in h]
            for k in range(3):
                reference[k] += sum(mp.mpf(float(higher[j, i, k])) for j in range(4))
            f, n = map(mp.mpf, (float(source[i, 0, 0].real), float(source[i, 2, 2].real)))
            delta = [(1-f-n)*b[k]-reference[k] for k in range(3)]
            answer.append([sum(h[k]*delta[k] for k in range(3)),
                           -e/a*delta[2], e/a*(n-f), ell/(2*r)*delta[1]])
    return np.array(answer, float)


@pytest.mark.parametrize('panel', ['mid/12_1', 'mid/23_1', 'mid/32_-1'])
def test_same_producer_stable_subtraction_matches_independent_arithmetic(retained, panel):
    metadata, E, C = _inputs(retained, panel)
    before = C.copy()
    result = incoming_high_energy_source(E, metadata, C)
    mass = metadata['channel']['compact_mass']
    angular = metadata['angular_sign']*metadata['channel']['angular_eigenvalue']
    _, orders = incoming_adiabatic_reference(E, mass, angular)
    independent = _high_precision_unrationalized(E, mass, angular, C, orders[1:])
    arithmetic_error = float(np.max(abs(result['kernels']-independent)))
    assert arithmetic_error < 3e-16
    assert max(result['diagnostics']['canonical_partner_bloch_residual']) < 3e-14
    assert max(result['diagnostics']['zeroth_reference_residual']) < 3e-14
    assert max(result['diagnostics']['direct_kernel_difference']) < 3e-11
    assert np.array_equal(C, before)
    assert np.array_equal(result['source_covariance'], before)
    assert np.array_equal(result['source_horizon_coherence'], C[:, 0, 1])
    assert np.count_nonzero(result['projected_horizon_coherence']) == 0
    assert not result['new_cutoff_or_filter']
    assert C[0, 2, 2].real > 0
    assert result['thermal_kernels'][0, 2] > 0
    assert result['kernels'][0, 2] == E[0]/AXIAL*(C[0, 2, 2].real-C[0, 0, 0].real)
    print(json.dumps({'panel': panel, 'energies': E.tolist(),
                      'unrationalized_60_digit_residual': arithmetic_error,
                      'direct_double_subtraction_difference': float(max(result['diagnostics']['direct_kernel_difference'])),
                      'stable_roundoff_indicator': float(result['roundoff_indicator'].max()),
                      'tiny_thermal_T01_retained': float(result['thermal_kernels'][0, 2])}))


def test_changed_tiny_source_is_not_silently_replaced(retained):
    metadata, E, C = _inputs(retained, 'mid/12_1')
    C[0, 2, 2] *= 2
    with pytest.raises(ValueError, match='exact unchanged horizon/incoming source'):
        incoming_high_energy_source(E, metadata, C)
