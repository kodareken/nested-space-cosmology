"""Focused checks of the action-level reference construction."""
import numpy as np
import pytest

from recursive_horizons.nsc_spatial_reference_symbol import supplied_KS_metric_jets, reference_projector
from recursive_horizons.nsc_reference_band_action import band_frame, values, exact_massless_reference


def test_band_frame_and_time_connection_are_both_required():
    fields = supplied_KS_metric_jets(.019, .23)
    P = reference_projector(*fields, [.2], [np.pi/2], [np.sqrt(5)])
    Pi = np.array([np.diag([1., 0.])]); frame = band_frame(P, Pi)
    assert max(v.max() for v in frame['scaled_residuals'].values()) < 3e-11
    without_connection = values(frame['effective'])-values(frame['connection'])
    defect = without_connection@Pi-Pi@without_connection
    assert np.max(abs(defect)) > 1e-3
    N, beta, a, r = fields
    expected = -P['gap'][0]-.2*beta.value[0, 0, 0]
    assert abs(frame['band_energy'][0, 0]-expected) < 3e-14
    with pytest.raises(ValueError, match='chiral'):
        exact_massless_reference(P, [0.], [0.], [.2])


def test_massless_band_chart_cannot_be_normalized_through_zero():
    P = reference_projector(*supplied_KS_metric_jets(.019, .23), [.2], [0.], [0.])
    P = exact_massless_reference(P, [0.], [0.], [.2])
    with pytest.raises(ValueError, match='zero-overlap'):
        band_frame(P, np.array([np.diag([1., 0.])]))
    frame = band_frame(P, np.array([np.diag([0., 1.])]))
    assert max(v.max() for v in frame['scaled_residuals'].values()) < 3e-11
    assert np.max(abs(P['orders'][1:])) == 0.
