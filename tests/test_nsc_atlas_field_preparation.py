"""Independent geometry, characteristic data and operator-bound checks."""
import json
import math

import numpy as np
import pytest

from recursive_horizons.nsc_atlas_field_preparation import (
    coefficient_from_radius_derivatives, exterior_coefficient, preparation_record,
    radius_from_null, radius_product, require_observer_time, source_profile,
    validate_preparation_domain, validate_preparation_inputs,
)


@pytest.fixture(scope="module")
def record():
    return preparation_record()


def test_retained_side_and_radius_bounds_along_incoming_rays(record):
    geometry = record["geometry"]
    wall = record["boundary_samples"]["wall"]
    U, V = np.asarray(wall["U"]), np.asarray(wall["V"])
    assert np.all(np.diff(U) > 0) and np.all(np.diff(V) > 0)
    for fraction in (0., .25, .5, .75, 1.):
        ray_U = geometry["U0"]+fraction*(U-geometry["U0"])
        assert np.all(ray_U <= U)  # Actual retained side, not quadrant alone.
        radius = radius_from_null(ray_U, V)
        assert np.min(radius) >= .4-1e-14
        assert np.max(radius) <= .55+1e-14
        assert np.allclose(radius_product(radius), ray_U*V, atol=3e-16, rtol=2e-15)
    assert np.allclose(radius_from_null(U, V), wall["R_over_ell"], atol=2e-15)


def test_boundaries_are_initial_data_before_the_old_slice(record):
    boundary = record["boundary_samples"]
    g = record["geometry"]
    assert np.all(np.asarray(boundary["source"]["U"]) == g["U0"])
    assert np.all(np.asarray(boundary["bottom"]["V"]) == g["Va"])
    assert boundary["bottom"]["U"][-1] == g["Ua"]
    assert record["preparation"]["cap_initial_t_over_ell"] < -1
    assert record["preparation"]["remote_generation_solved"] is False


def test_kernel_from_implicit_geometry_and_schwarzschild_operator():
    x = np.array([.4, .425, .475, .5, .55])
    direct = exterior_coefficient(x)
    assert np.allclose(coefficient_from_radius_derivatives(x), direct, rtol=3e-15)
    # Independently transform the borrowed tortoise-coordinate scalar potential.
    M = .4
    potential = (1-2*M/x)*2*M/x**3
    assert np.allclose(-4*M*M*potential/radius_product(x), direct, rtol=5e-16)
    # Mixed finite differences of the radius itself do not use either kernel formula.
    U = -np.sqrt(radius_product(x))
    V = U.copy()
    h = 2e-5
    mixed = (radius_from_null(U+h, V+h)-radius_from_null(U+h, V-h)
             -radius_from_null(U-h, V+h)+radius_from_null(U-h, V-h))/(4*h*h)
    assert np.allclose(-mixed/x, direct, rtol=2e-6)


def test_compact_input_is_independent_and_corner_compatible(record):
    Va, Vb = record["geometry"]["Va"], record["geometry"]["Vb"]
    value, derivative = source_profile(np.array([Va, (Va+Vb)/2, Vb]), Va, Vb)
    assert np.array_equal(value, [0., 1., 0.])
    assert np.array_equal(derivative, [0., 0., 0.])
    assert record["preparation"]["independent_of_target_waveform"]
    validate_preparation_inputs(source_corner=float(value[0]))
    with pytest.raises(ValueError, match="admitted"):
        source_profile([Va-.01], Va, Vb)


def test_polynomial_has_exactly_c3_edge_matching():
    polynomial = np.polynomial.Polynomial([1., 0., -1.])**4
    for degree in range(4):
        assert np.array_equal(polynomial.deriv(degree)([-1., 1.]), [0., 0.])
    assert np.all(polynomial.deriv(4)([-1., 1.]) != 0)


@pytest.mark.parametrize("changed", [
    {"source_corner": 1.},
    {"source_corner": 1., "bottom_corner": 1.},
    {"initial_cap_field": 1.},
    {"initial_cap_time_derivative": 1.},
    {"xi": 0.},
    {"surface_scalar_coupling": .1},
    {"source_corner": math.nan},
])
def test_incompatible_data_or_changed_theory_requires_new_candidate(changed):
    with pytest.raises(ValueError):
        validate_preparation_inputs(**changed)


def test_same_quadrant_does_not_certify_retained_preparation(record):
    g = record["geometry"]
    with pytest.raises(ValueError, match="retained"):
        validate_preparation_domain(g["Ua"]+.001, g["Ua"], g["Ub"], g["Va"], g["Vb"])
    with pytest.raises(ValueError, match="areal radius"):
        validate_preparation_domain(g["U0"]-.001, g["Ua"], g["Ub"], g["Va"], g["Vb"])
    with pytest.raises(ValueError, match="sector IV"):
        radius_from_null(-1., -1.)


def test_analytic_contraction_bound_is_independent_of_clock_quadrature(record):
    operator = record["operator"]
    bound = operator["analytic_norm_bound"]
    assert 0 < operator["rectangle_norm_bound"] <= bound < .14
    # AM-GM bounds the coordinate area for every positive ratio, without a wall integration.
    Fc, Fb = float(radius_product(.4)), float(radius_product(.55))
    ratios = np.exp(np.linspace(-3, 3, 51))
    areas = Fc+Fb-Fc*ratios-Fb/ratios
    assert np.all(areas*float(exterior_coefficient(.4)) <= bound+2e-14)
    assert operator["midpoint_wall_reduced_field_min_over_epsilon"] > .83
    assert operator["solution_iteration_performed"] is False
    assert operator["new_waveform_evaluated"] is False


def test_full_wall_window_has_no_return_and_regular_central_future(record):
    observer = record["observer"]
    wall = record["boundary_samples"]["wall"]
    assert observer["no_return_margin_over_ell"] > 2.6
    assert observer["analytic_no_return_log_margin_lower_bound"] > .2
    delta_t = wall["t_flat_over_ell"][-1]-wall["t_flat_over_ell"][0]
    assert delta_t < observer["flat_time_increment_upper_bound_over_ell"]
    assert observer["initial_cap_t_over_ell"] < observer["arrival_window_over_ell"][0]
    assert np.all(np.diff(wall["arrival_t_C_over_ell"]) > 0)
    assert np.all(np.asarray(wall["dT_dtau"]) > 0)
    assert record["checks"]["normal_orientation_max_error"] < 2e-14
    assert record["checks"]["reduced_field_jump_coefficient"] == 0.


def test_observer_clock_domain_does_not_extend_to_zero(record):
    window = record["observer"]["arrival_window_over_ell"]
    for endpoint in window:
        require_observer_time(endpoint, window)
        require_observer_time(3*endpoint, window, 3)
    for outside in (window[0]-.001, window[1]+.001, 0., math.inf):
        with pytest.raises(ValueError, match="outside"):
            require_observer_time(outside, window)


def test_unit_rescaling_preserves_preparation_and_operator(record):
    scaled = preparation_record(2.75)
    assert scaled["inputs"]["length_scale"] == 2.75
    assert scaled["geometry"] == record["geometry"]
    assert scaled["operator"] == record["operator"]
    for key in ("arrival_window_over_ell", "initial_cap_t_over_ell", "no_return_margin_over_ell"):
        assert np.allclose(scaled["observer"][key], record["observer"][key], rtol=2e-14, atol=2e-14)
    U, V = record["geometry"]["U0"], record["geometry"]["Va"]
    assert float(radius_from_null(U, V, 2.75)) == pytest.approx(2.75*.4)
    with pytest.raises(ValueError, match="scale"):
        preparation_record(0)


def test_record_preserves_causal_exclusion_and_finite_budget(record):
    assert record["access"][0]["verdict"] == "excluded_in_reference"
    assert record["access"][1]["verdict"] == "finite_conditional_preparation_admitted"
    assert record["budget"]["boundary_samples"] == 163
    assert record["budget"]["boundary_samples"] <= record["budget"]["maximum_samples"]
    assert record["checks"]["quadrature_32_64_difference"] < 1e-9
    assert not record["checks"]["quadrature_is_rigorous_enclosure"]
    json.dumps(record, allow_nan=False)
