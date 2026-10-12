"""Narrow checks of the new actual-metric source representation."""
import numpy as np

from recursive_horizons.nsc_nonlinear_ks_source import OwnedMetricSamples
from recursive_horizons.nsc_transmitting_history_modes import SuppliedKSHarmonicMetric
from recursive_horizons.nsc_transmitting_history_jets import raw_KS_amplitude_directions


def test_vectorized_history_and_vertices_are_the_owned_history():
    x = np.linspace(-2., 1.2, 17)
    p = SuppliedKSHarmonicMetric((.002, .001, .003, .001), .4, (0., .05))
    owner = OwnedMetricSamples(x, p)
    for t in (0., .019, .05):
        s = owner.at(t); expected = p.values(t, x)
        assert np.max(abs(s['metric']-expected)) < 3e-14
        J = raw_KS_amplitude_directions(p, t, x, expected)
        assert np.max(abs(s['directions']-J)) < 3e-14


def test_actual_slice_frame_time_derivative():
    x = np.linspace(-.8, .8, 11)
    p = SuppliedKSHarmonicMetric((.002, .001, .003, .001), .4, (0., .05))
    owner = OwnedMetricSamples(x, p); t, h = .019, 1e-5
    fd = (owner.at(t-2*h)['frame']-8*owner.at(t-h)['frame']+
          8*owner.at(t+h)['frame']-owner.at(t+2*h)['frame'])/(12*h)
    assert np.max(abs(fd-owner.at(t)['frame_t'])) < 3e-8
