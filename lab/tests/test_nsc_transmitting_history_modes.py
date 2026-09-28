"""New spatial history owner: flux, reference limit and domain guards."""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'src'))
from recursive_horizons.nsc_transmitting_history_modes import (
    RelativeModePropagator, SuppliedKSHarmonicMetric, reference_fields_on_grid,
)


def test_supplied_smooth_history_keeps_the_preparation_and_flux_domain():
    x = np.linspace(-2., 1.2, 81)
    owner = RelativeModePropagator(x, 1.3, 2.)  # Algebra control, no sector assignment.
    metric = SuppliedKSHarmonicMetric((.002, .001, .003, .001), support_time=(0., .05))
    assert np.max(abs(metric.values(0., x)-owner.reference)) < 3e-11
    assert np.max(abs(metric.values(.05, x)-owner.reference)) < 3e-11
    _, residual = owner.generator(metric.values(.025, x))
    assert residual['norm_flux_algebra'] < 3e-11
    assert residual['maximum_characteristic_speed'] < 0


def test_relative_reference_limit_does_not_fake_a_closed_source_space():
    x = np.linspace(-2., 1.2, 33)
    owner = RelativeModePropagator(x, 1.3, 2.)
    # Only the zero forcing limit is being tested; this array is not C_PG.
    fields = np.zeros((len(x), 2, 2, 3), complex)
    out = owner.evolve(np.array([-.2, .2]), fields, np.array([0., .01]),
                       SuppliedKSHarmonicMetric((0., 0., 0., 0.), support_time=(0., .01)))
    assert np.max(abs(out['difference_characteristic_field'])) == 0
    assert out['complete_EndpointBranchJets'] is None
    assert out['full_covariance_integral'] is None
    assert out['physical_history_selected'] is False


def test_contour_and_causal_domain_misuse_are_rejected():
    x = np.linspace(-2., 1.2, 33)
    with pytest.raises(ValueError, match='real source energies'):
        reference_fields_on_grid(np.array([.2+.1j]), 1., 2., np.zeros((1, 2, 3)), x)
    owner = RelativeModePropagator(x, 1., 2.)
    with pytest.raises(ValueError, match='causal buffer'):
        owner.evolve(np.array([.2]), np.zeros((len(x), 1, 2, 3)), np.array([0., 1.]),
                     SuppliedKSHarmonicMetric((0., 0., 0., 0.)))
    with pytest.raises(ValueError, match='reference exterior'):
        metric = owner.reference.copy()
        metric[0, 0] += .01
        owner.generator(metric)
