"""Metric interpolation and the observer's actual radar clock, without source data."""
import numpy as np
import pytest

from measure_nsc_conformal_lightcone import FourierHermiteMetric, observer_radar_distances
from recursive_horizons.nsc_atlas_characteristic_initialization import Budget


def test_saved_rate_hermite_and_fourier_recover_an_analytic_metric():
    times = np.array([-.1, .05, .3])
    length, count = 8., 9
    nodes = length*np.arange(count)/count
    wave = 2*np.pi/length

    def fields(t, x):
        return (3+.1*np.cos(wave*x)+.02*t+.03*t**3*np.sin(wave*x),
                .7+.01*t**2*np.cos(2*wave*x))

    radius, Q = fields(times[:, None], nodes[None, :])
    rt = .02+.09*times[:, None]**2*np.sin(wave*nodes[None, :])
    Qt = .02*times[:, None]*np.cos(2*wave*nodes[None, :])
    budget = Budget(maximum_evaluations=10000, cpu_seconds=10.)
    metric = FourierHermiteMetric(times, radius, Q, rt, Qt, length=length, budget=budget)
    t, x = np.array([-.05, .05, .2, .3]), np.array([.2, 1., -1., 9.3])
    measured = metric(t, x)
    exact_r, exact_Q = fields(t, x)
    assert measured['radius'] == pytest.approx(exact_r, abs=3e-14)
    assert measured['Q'] == pytest.approx(exact_Q, abs=3e-14)
    assert budget.evaluations == 4*count*len(t)
    with pytest.raises(ValueError, match='window'):
        metric(.5, 1.)


def test_radar_uses_round_trip_proper_clock_in_flat_and_de_sitter_metrics():
    lookbacks = (.025, .075, .125)

    def minkowski(t, x):
        t, x = np.broadcast_arrays(t, x)
        radius = 2+x
        return {'radius': radius, 'Q': 1/radius}

    result = observer_radar_distances(minkowski, np.array([.05, .175, .3]),
        reception_time=.3, reception_x=1., lookbacks=lookbacks, quadrature_order=32)
    assert result == pytest.approx(lookbacks, abs=2e-15)
    H = .8

    def de_sitter(t, x):
        t, x = np.broadcast_arrays(t, x)
        a = -1/(H*t)
        return {'radius': a*(2+x), 'Q': 1/(2+x)}

    reception = -1/H
    times = np.array([reception-.25, reception-.12, reception])
    exact = np.log1p(2*H*np.asarray(lookbacks))/(2*H)
    for order in (32, 64):
        result = observer_radar_distances(de_sitter, times, reception_time=reception,
            reception_x=1., lookbacks=lookbacks, quadrature_order=order)
        assert result == pytest.approx(exact, abs=3e-15)
    with pytest.raises(ValueError, match='interval'):
        observer_radar_distances(de_sitter, times, reception_time=reception,
            reception_x=1., lookbacks=(.2,), quadrature_order=32)
