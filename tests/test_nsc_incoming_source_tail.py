"""Tail arithmetic of the owned approximation; no field ODE or scattering."""
import hashlib
import json
from pathlib import Path

import mpmath as mp
import numpy as np
import pytest

from recursive_horizons.nsc_incoming_source_tail import (
    _riccati_at_one, incoming_source_tail, evaluate_tail_series,
    transformed_tail_quadrature,
)
from recursive_horizons.nsc_pg_high_energy import riccati_coefficients


ROOT = Path(__file__).resolve().parents[1]
SOURCE_PAYLOAD = 'results/development/artifacts/nsc-incoming-spectral-panels.efbce40b2f959c7d0294689661d4e2e7193ade965ece0918c72c5a072ec8ad73.npz'


@pytest.fixture(scope='module')
def inputs():
    channels = json.loads((ROOT/'results/development/nsc-mode-resolved-cauchy-state.json').read_text())['channels']
    path = ROOT/SOURCE_PAYLOAD
    assert hashlib.sha256(path.read_bytes()).hexdigest() == 'efbce40b2f959c7d0294689661d4e2e7193ade965ece0918c72c5a072ec8ad73'
    with np.load(path, allow_pickle=False) as data:
        yield channels, data


@pytest.mark.parametrize('group,lower', [(13, 160.), (22, 160.), (32, 320.)])
def test_same_approximation_tail_and_stored_source_band(inputs, group, lower):
    channels, data = inputs
    channel = channels[group]
    result = incoming_source_tail(channel, lower, maximum_order=13, precision=50)
    q16 = transformed_tail_quadrature(channel, lower, points=16, precision=50)
    q24 = transformed_tail_quadrature(channel, lower, points=24, precision=60)
    series = result['corrections_by_order'][13]
    assert max(abs(q24-q16)) < 3e-20
    assert max(abs(series-q24)) < 3e-18
    assert result['nonintegrable_coefficient_residual'] < 1e-25
    assert result['nonintegrable_coefficients_exactly_zero']
    assert np.count_nonzero(result['coefficients_rho_parallel_sphere'][:2]) == 0
    assert result['regular_endpoint_identity']['vacuum_minus_ad1'] == ['0', '0', '0']
    assert result['numerical_differentiation_chop'] is False
    assert result['rigorous_physical_remainder_bound'] is None
    names = ['group13/mid'] if group == 13 else [f'mid/{group}_1', f'mid/{group}_-1']
    energies = data[names[0]+'/energies']
    sample = (energies > .75*lower) & (energies < lower)
    expected = result['multiplicity_factor']*sum(data[name+'/kernels'][sample] for name in names)
    indicator = result['multiplicity_factor']*sum(data[name+'/arithmetic_indicator'][sample] for name in names)
    actual = evaluate_tail_series(result, energies[sample])
    lower_order = evaluate_tail_series(result, energies[sample], order=11)
    error = abs(actual-expected)
    assert np.max(error) < 3e-13
    assert np.all(error <= 16*indicator+2*abs(actual-lower_order)+3e-20)
    with mp.workdps(50):
        precise = np.array([complex(v) for v in _riccati_at_one(
            mp.mpf(float(channel['compact_mass'])), mp.mpf(float(channel['angular_eigenvalue'])))])
    ordinary = riccati_coefficients([1.], channel['compact_mass'], channel['angular_eigenvalue'], 16)[0][0]
    recurrence_error = float(np.max(abs(precise-ordinary)/(1+abs(ordinary))))
    assert recurrence_error < 3e-11
    print(json.dumps({'group': group, 'lower': lower, 'tail13': series.tolist(),
                      'order11_13': result['series_order_difference'].tolist(),
                      'quadrature16_24': float(np.max(abs(q24-q16))),
                      'series_vs_quadrature': float(np.max(abs(series-q24))),
                      'stored_source_band_error': float(np.max(error)),
                      'Riccati_precision_conversion': recurrence_error}))


def test_wrong_existing_split_is_rejected(inputs):
    channels, _ = inputs
    with pytest.raises(ValueError, match='archived middle/tail split'):
        incoming_source_tail(channels[32], 160.)
    with pytest.raises(ValueError, match='archived middle/tail split'):
        transformed_tail_quadrature(channels[32], 160.)
    with pytest.raises(ValueError, match='owned lower split'):
        incoming_source_tail(channels[13], 500.)
