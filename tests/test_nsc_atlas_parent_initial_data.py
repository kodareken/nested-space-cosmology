"""Independent integration controls for published initial-data mathematics."""
import json

import numpy as np
import pytest
from numpy.testing import assert_allclose

from recursive_horizons.nsc_atlas_parent_initial_data import (
    BrillLindquistData, adm_masses, conformal_factor, conformal_fields, parent_record,
    sphere_distortion, two_punctures,
)


def test_lightweight_plot_evaluator_supports_grid_shape_and_masked_domain():
    data = two_punctures()
    x, z = np.meshgrid(np.linspace(-3, 3, 11), np.linspace(-4, 4, 13))
    points = np.stack([x, np.zeros_like(x), z], axis=-1)
    psi = conformal_factor(points, data.masses, data.positions)
    assert psi.shape == (13, 11)
    assert_allclose(psi, conformal_fields(points, data)["psi"], rtol=1e-15)
    with pytest.raises(ValueError):
        conformal_factor([[0, 0, -6]], data.masses, data.positions)


def test_analytic_derivatives_against_finite_differences():
    data = two_punctures()
    point = np.array([2.1, -1.3, 0.8])
    step = 2e-5
    fields = conformal_fields(point, data)
    directions = np.eye(3) * step
    finite_gradient = np.array([(conformal_fields(point + e, data)["psi"]
                                 - conformal_fields(point - e, data)["psi"]) / (2 * step)
                                for e in directions])
    finite_hessian = np.stack([(conformal_fields(point + e, data)["gradient"]
                                - conformal_fields(point - e, data)["gradient"]) / (2 * step)
                               for e in directions], axis=-1)
    assert_allclose(fields["gradient"], finite_gradient, atol=4e-12, rtol=2e-8)
    assert_allclose(fields["hessian"], finite_hessian, atol=4e-13, rtol=2e-8)


def test_vacuum_scalar_constraint_and_nonflat_spatial_metric():
    data = two_punctures()
    points = [[2.0, 1.0, -6.0], [3.0, 4.0, 2.0], [0.0, 0.0, 0.0]]
    fields = conformal_fields(points, data)
    assert_allclose(fields["scalar_curvature"], 0.0, atol=2e-17)
    assert np.all(fields["psi"] > 1)
    # A separate second-difference Laplacian of psi checks harmonicity.
    point = np.array(points[1])
    step = 1e-3
    laplace = sum(float(conformal_fields(point + e, data)["psi"])
                  + float(conformal_fields(point - e, data)["psi"])
                  - 2 * float(conformal_fields(point, data)["psi"])
                  for e in np.eye(3) * step) / step**2
    assert abs(laplace) < 3e-9


def test_adm_charges_against_outer_flux_and_inverted_end_asymptotics():
    data = two_punctures()
    charges = adm_masses(data)
    assert charges["outer_end"] == 2.5
    assert charges["puncture_ends"] == [1.0625, 1.5625]
    assert charges["outer_minus_sum_puncture_ends"] == -0.125
    # ADM conformal flux -1/(2pi) integral d_r psi dS; a distant
    # spherical quadrature is independent of the explicit mass-sum helper.
    mu, weights = np.polynomial.legendre.leggauss(24)
    radius = 40.0
    directions = np.stack([np.sqrt(1 - mu**2), np.zeros_like(mu), mu], axis=-1)
    gradient = conformal_fields(radius * directions, data)["gradient"]
    flux_mass = -radius**2 * np.dot(weights, np.sum(gradient * directions, axis=-1))
    assert_allclose(flux_mass, charges["outer_end"], atol=2e-14)
    for i, mass in enumerate(data.masses):
        inverse_radius = 1e5
        point = np.asarray(data.positions[i]) + np.array([mass**2 / (4 * inverse_radius), 0, 0])
        inverse_psi = mass / (2 * inverse_radius) * conformal_fields(point, data)["psi"]
        inferred = 2 * inverse_radius * (inverse_psi - 1)
        assert_allclose(inferred, charges["puncture_ends"][i], atol=1e-9)


def test_scale_covariance_keeps_dimensionless_geometry():
    data = two_punctures()
    factor = 7.25
    scaled = BrillLindquistData(tuple(factor * m for m in data.masses),
                               tuple(tuple(factor * np.asarray(p)) for p in data.positions))
    points = np.array([[2.0, 1.0, -6.0], [3.0, 4.0, 2.0]])
    original = conformal_fields(points, data)
    changed = conformal_fields(points * factor, scaled)
    for name, power in [("psi", 0), ("metric", 0), ("gradient", 1), ("hessian", 2)]:
        assert_allclose(changed[name] * factor**power, original[name], rtol=2e-14, atol=1e-16)
    assert_allclose(adm_masses(scaled)["puncture_ends"], np.array(adm_masses(data)["puncture_ends"]) * factor)
    first, second = sphere_distortion(data, 0, 1), sphere_distortion(scaled, 0, factor)
    for key in ("psi", "relative_psi_error", "relative_angular_metric_error",
                "gaussian_curvature_times_coordinate_radius_squared", "quadrupole_over_monopole"):
        assert_allclose(first[key], second[key], atol=2e-15)


def test_local_spheres_retain_angular_information_and_mean_value():
    data = two_punctures()
    profile = sphere_distortion(data, 0, 1)
    assert profile["psi"][-1] > profile["psi"][0]
    assert profile["is_exactly_round"] is False
    assert np.ptp(profile["gaussian_curvature_times_coordinate_radius_squared"]) > 0
    # The harmonic companion's average over a sphere equals its centre
    # value; evaluating the full 3D adapter gives an independent control.
    mu, weights = np.polynomial.legendre.leggauss(24)
    directions = np.stack([np.sqrt(1 - mu**2), np.zeros_like(mu), mu], axis=-1)
    points = np.array(data.positions[0]) + directions
    average = np.dot(weights, conformal_fields(points, data)["psi"]) / 2
    assert_allclose(average, profile["monopole_psi"], atol=1e-14)
    expected_max = (1.5 / 2) * (1 / 11 - 1 / 12) / profile["monopole_psi"]
    assert_allclose(profile["maximum_absolute_relative_psi_error"], expected_max, atol=2e-16)
    assert sphere_distortion(data, 1, 1.5)["maximum_absolute_relative_psi_error"] != profile["maximum_absolute_relative_psi_error"]


@pytest.mark.parametrize("points,exclusion", [([0, 0, -6], 0), ([0, 0, -5], 1), ([np.nan, 0, 0], 0), ([1, 2], 0), ([1, 2, 3], -1)])
def test_invalid_chart_samples_rejected(points, exclusion):
    with pytest.raises(ValueError):
        conformal_fields(points, two_punctures(), exclusion_radius=exclusion)


def test_invalid_punctures_and_local_sphere_domains_rejected():
    for masses, positions in [((0,), ((0, 0, 0),)), ((True,), ((0, 0, 0),)),
                              ((1, 2), ((0, 0, 0), (0, 0, 0))), ((1,), ((0, np.inf, 0),))]:
        with pytest.raises(ValueError):
            BrillLindquistData(masses, positions)
    with pytest.raises(ValueError):
        sphere_distortion(two_punctures(), 0, 12)
    with pytest.raises(ValueError):
        sphere_distortion(two_punctures(), 0, 1, samples=4098)


def test_parent_record_is_deterministic_and_retains_missing_physics():
    record = parent_record()
    assert json.dumps(record, sort_keys=True, allow_nan=False) == json.dumps(parent_record(), sort_keys=True, allow_nan=False)
    assert not record["globally_spherical"]
    assert record["causal_channels"] == []
    assert not record["child_formation_established"]
    assert not record["direct_spherical_child_attachment_established"]
    assert not record["horizons"]["computed"]
    assert not record["horizons"]["common_horizon_excluded"]
    assert record["clock"]["lapse"] is None
    assert record["admission"]["novelty_claim"] is False
