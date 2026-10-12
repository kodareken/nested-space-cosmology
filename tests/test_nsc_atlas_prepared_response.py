import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import j0, j1

from recursive_horizons import nsc_atlas_prepared_response as response
from recursive_horizons.nsc_atlas_field_preparation import (
    preparation_record, exterior_coefficient, radius_from_null,
)


@pytest.fixture(scope="module")
def record():
    return response.prepared_response_record()


def test_exact_transverse_kernel_primitive_has_independent_quadrature_check():
    a, b, v = -.315, -.29, -2.45
    integral = quad(lambda u: float(exterior_coefficient(radius_from_null(u, v))), a, b,
                    epsabs=1e-13, epsrel=1e-13)[0]
    assert float(response.transverse_kernel_integral(a, b, v)) == pytest.approx(integral, rel=2e-12)
    with pytest.raises(ValueError, match="reversed"):
        response.transverse_kernel_integral(b, a, v)


def test_prepared_input_is_retained_and_not_selected_from_response(record):
    prior = preparation_record()
    assert record["preparation"] == prior["preparation"]
    assert record["inputs"] == prior["inputs"]
    assert not record["method"]["exterior_kernel_evaluated_in_excised_cap"]
    assert record["preparation"]["remote_generation_solved"] is False


def test_post_input_tail_obeys_independent_positive_analytic_brackets(record):
    s = np.array(record["samples"]["source_coordinate"])
    g = np.array(record["samples"]["source_chi_over_epsilon"])
    J = np.array(record["samples"]["wall_dchi_dV_over_epsilon"])
    phi = np.array(record["samples"]["ell_phi_C_over_epsilon"])
    assert np.all(g[s > .75] == 0)
    assert np.all(J[s > .75] >= record["analytic_tail"]["uniform_lower_J_over_epsilon"])
    assert np.all(phi[s > .75] > 0)
    for row in record["analytic_tail"]["independent_controls"]:
        assert row["lower_from_B"] <= row["J"] <= row["upper_from_B"]
        assert row["B_quadrature_difference"] < 1e-8*row["B_64"]


def test_refinement_is_small_relative_to_the_tail_not_just_the_large_pulse(record):
    checks = record["checks"]
    assert checks["normalized_mixed_refinement_error"] <= 1
    assert checks["normalized_geometry_refinement_error"] <= 1
    assert checks["tail_J_refinement_difference"] < .1*record["analytic_tail"]["uniform_lower_J_over_epsilon"]
    assert checks["geometry_129_257_tail_J_difference"] < .1*record["analytic_tail"]["uniform_lower_J_over_epsilon"]
    assert checks["refinement_is_rigorous_bound"] is False
    assert checks["observer_area_identity_difference"] < 5e-5


def test_no_kernel_has_exact_zero_correction_and_no_spurious_tail(monkeypatch):
    prepared = preparation_record()
    geometry = response._geometry(129)
    monkeypatch.setattr(response, "exterior_coefficient", lambda x: np.zeros_like(x))
    extract, _ = response._solve(129, geometry, prepared, response._Budget())
    sample = extract(np.linspace(0, 1, 129))
    assert np.array_equal(sample["source_chi_over_epsilon"], sample["wall_chi_over_epsilon"])
    assert np.array_equal(sample["ell_phi_C_over_epsilon"], sample["zero_kernel_phi_C_over_epsilon"])
    assert np.all(sample["ell_phi_C_over_epsilon"][97:] == 0)
    with pytest.raises(ValueError, match="exits"):
        extract([1.001])


def test_dimensionless_response_survives_regional_unit_rescaling(record):
    scaled = response.prepared_response_record(7.25)
    assert scaled["inputs"]["length_scale"] == 7.25
    assert scaled["samples"] == record["samples"]
    assert scaled["analytic_tail"] == record["analytic_tail"]


@pytest.mark.parametrize("value", [0, -1, float("nan"), float("inf")])
def test_invalid_units_do_not_enter_evolution(value):
    with pytest.raises(ValueError, match="length scale"):
        response.prepared_response_record(value)


def test_exhausted_budget_and_unadmitted_discretization_fail():
    budget = response._Budget()
    with pytest.raises(RuntimeError, match="budget"):
        budget.check(response.MAX_KERNEL_VALUES+1)
    with pytest.raises(ValueError, match="field grids"):
        response._solve(1025, {}, {}, response._Budget())
    with pytest.raises(ValueError, match="geometry"):
        response._geometry(513)


def test_budget_is_separate_from_historical_boundary_sample_limit(record):
    budget = record["budget"]
    assert budget["field_nodes"] == [129, 257, 513]
    assert max(budget["rhs_calls"]) <= budget["rhs_calls_max_per_solve"]
    assert budget["kernel_values"] <= budget["kernel_values_max"]
    assert budget["output_samples"] == 257
    assert budget["cpu_seconds_max"] == 60


def test_constant_kernel_rectangle_matches_borrowed_bessel_solution(monkeypatch):
    # A manufactured numerical benchmark, not another admitted spacetime.
    prepared = {"geometry": {"U0": -.3, "Ua": -.2, "Ub": -.2, "Va": -2.5, "Vb": -2.25}}
    geometry = {
        "wall": lambda v: -.2,
        "arrival": lambda v, derivative=0: np.ones_like(v) if derivative else np.asarray(v),
        "Vdot": lambda v: np.ones_like(v),
    }
    monkeypatch.setattr(response, "exterior_coefficient", lambda x: np.ones_like(x))
    monkeypatch.setattr(response, "source_profile", lambda v, a, b: (np.asarray(v)-a, np.ones_like(v)))
    extract, _ = response._solve(257, geometry, prepared, response._Budget())
    t = np.linspace(.01, 1, 65)
    sample = extract(t)
    y = .25*t
    z = .1*y
    expected = y*j1(2*np.sqrt(z))/np.sqrt(z)
    assert np.max(abs(sample["wall_chi_over_epsilon"]-expected)) < 2e-8
    assert np.max(abs(sample["wall_dchi_dV_over_epsilon"]-j0(2*np.sqrt(z)))) < 2e-8


def test_fixed_linear_field_scales_without_parameter_retuning(monkeypatch):
    prepared, geometry = preparation_record(), response._geometry(129)
    first, _ = response._solve(129, geometry, prepared, response._Budget())
    t = np.array([.375, .5, .625, .8125, 1.])
    baseline = first(t)
    source = response.source_profile
    monkeypatch.setattr(response, "source_profile", lambda *args: tuple(3*x for x in source(*args)))
    second, _ = response._solve(129, geometry, prepared, response._Budget())
    scaled = second(t)
    assert np.allclose(scaled["ell_phi_C_over_epsilon"], 3*baseline["ell_phi_C_over_epsilon"],
                       rtol=1e-6, atol=1e-8)


def test_interpolation_rejects_domain_exits_and_interior_derivative_reversal():
    geometry = response._geometry(129)
    assert np.isnan(geometry["wall"](geometry["Va"]-.01))
    assert np.isnan(geometry["arrival"](geometry["Vb"]+.01))
    # Positive endpoint slopes alone do not certify a monotone cubic.
    curve = response.CubicHermiteSpline([0., 1.], [0., .1], [1., 1.])
    assert response._minimum_derivative(curve) < 0
