"""Reject the shortcuts that the common Cauchy construction replaces."""
import numpy as np
import pytest

from recursive_horizons.nsc_common_ks_trace import restrict_resolved_modes, incoming_KS_state


def test_seed_covariance_is_not_a_resolved_mode_map():
    with pytest.raises(ValueError, match='source mode columns'):
        restrict_resolved_modes(1., [1., 2.], np.tile(np.eye(2), (2, 1, 1)))


def test_slice_inside_perturbation_cannot_be_declared_unchanged_input():
    with pytest.raises(ValueError, match='upstream'):
        incoming_KS_state([1.], np.zeros((1, 2, 3)), np.zeros((1, 3, 3)), np.zeros((1, 3, 3)), rho=.5)


def test_contour_frequency_is_not_a_physical_covariance_state():
    with pytest.raises(ValueError, match='real physical frequencies'):
        restrict_resolved_modes(1., [1.+.01j], np.zeros((1, 2, 3)))
