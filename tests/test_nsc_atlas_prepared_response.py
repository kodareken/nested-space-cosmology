import numpy as np
import pytest
from scipy.integrate import quad
from scipy.special import j0, j1, jv

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


def test_bulk_constant_kernel_matches_independent_field_and_two_derivatives(monkeypatch):
    prepared = {"geometry": {"U0": -.3, "Ua": -.2, "Ub": -.2, "Va": -2.5, "Vb": -2.25}}
    geometry = {
        "wall": lambda v: -.2,
        "arrival": lambda v, derivative=0: np.ones_like(v) if derivative else np.asarray(v),
        "Vdot": lambda v: np.ones_like(v),
    }
    monkeypatch.setattr(response, "exterior_coefficient", lambda x: np.ones_like(x))
    monkeypatch.setattr(response, "source_profile", lambda v, a, b: (np.asarray(v)-a, np.ones_like(v)))
    wall, calls, bulk = response._solve(
        129, geometry, prepared, response.PreparedFieldBudget(), include_bulk=True)
    # Off-grid points test a partial integration cell, not only stored nodes.
    points = np.array([-.3, -.284321, -.231234, -.2])
    t, y = .83, .25*.83
    z = (points+.3)*y
    argument = 2*np.sqrt(z[1:])
    expected_chi = np.r_[y, y*j1(argument)/np.sqrt(z[1:])]
    expected_U = np.r_[-y*y/2, -y*y*jv(2, argument)/z[1:]]
    for order in (32, 64):
        sample = bulk(points, t, order=order)
        assert np.max(abs(sample["chi"]-expected_chi)) < 2e-9
        assert np.max(abs(sample["chi_U"]-expected_U)) < 2e-9
        assert np.max(abs(sample["chi_V"]-j0(2*np.sqrt(z)))) < 2e-9
        assert sample["chi"][-1] == pytest.approx(wall([t])["wall_chi_over_epsilon"][0], abs=1e-14)
    assert calls < response.MAX_RHS_PER_SOLVE


def test_bulk_zero_kernel_preserves_source_and_canonical_radius_factors(monkeypatch):
    monkeypatch.setattr(response, "exterior_coefficient", lambda x: np.zeros_like(x))
    field = response.solve_prepared_field(nodes=129, geometry_knots=129)
    points = np.linspace(field.domain["U0"], field.domain["Ua"], 7)
    sample = field.bulk(points, .4)
    V = field.domain["Va"]+.4*(field.domain["Vb"]-field.domain["Va"])
    g, gp = response.source_profile(V, field.domain["Va"], field.domain["Vb"])
    assert np.all(sample["chi"] == g)
    assert np.all(sample["chi_U"] == 0)
    assert np.all(sample["chi_V"] == gp)
    # Independent implicit-radius identity UV=F(x), not a finite difference.
    Fprime = -sample["x"]*np.exp(sample["x"]/.8)/.64
    assert np.allclose(Fprime*sample["x_U"], V, rtol=2e-15)
    assert np.allclose(Fprime*sample["x_V"], points, rtol=2e-15)
    assert np.allclose(sample["p"], g/sample["x"], rtol=2e-15)
    assert np.allclose(sample["p_U"], -g*sample["x_U"]/sample["x"]**2, rtol=2e-15)
    assert np.allclose(sample["p_V"], gp/sample["x"]-g*sample["x_V"]/sample["x"]**2, rtol=2e-15)
    assert np.all(field.bulk(points, .9)["chi"] == 0)


def test_bulk_domain_and_query_budget_guards_precede_evaluation(monkeypatch):
    monkeypatch.setattr(response, "exterior_coefficient", lambda x: np.zeros_like(x))
    budget = response.PreparedFieldBudget()
    field = response.solve_prepared_field(129, 129, budget)
    U0, Ua = field.domain["U0"], field.domain["Ua"]
    for points in ([U0-1e-12], [Ua+1e-12], [np.nan], [], [[U0]], np.full(514, U0)):
        with pytest.raises(ValueError, match="rectangle"):
            field.bulk(points, .5)
    for t in (-1e-12, 1+1e-12, np.nan, np.inf, [.5], True):
        with pytest.raises(ValueError, match="window"):
            field.bulk([U0], t)
    for order in (16, 32.0, True):
        with pytest.raises(ValueError, match="Gauss"):
            field.bulk([U0], .5, order=order)
    for times in ([], [[.5]], np.full(514, .5)):
        with pytest.raises(ValueError, match="wall query"):
            field.wall(times)
    with pytest.raises(ValueError, match="window"):
        field.wall([1.01])
    before = budget.kernel_values
    field.bulk([U0, Ua], .9, order=32)
    assert budget.kernel_values-before == 129+2*32*3
    budget.kernel_values = budget.kernel_values_max
    with pytest.raises(RuntimeError, match="budget"):
        field.bulk([U0], .5)


@pytest.mark.parametrize("kwargs", [
    {"cpu_seconds": 0}, {"cpu_seconds": 121}, {"cpu_seconds": np.nan},
    {"cpu_seconds": "120"}, {"cpu_seconds": True}, {"kernel_values_max": 0},
    {"kernel_values_max": 30_000_001}, {"kernel_values_max": 100.0},
])
def test_bulk_budget_cannot_extend_admitted_ceilings(kwargs):
    with pytest.raises(ValueError, match="budget"):
        response.PreparedFieldBudget(**kwargs)


def test_bulk_budget_clock_and_invalid_factory_arguments(monkeypatch):
    budget = response.PreparedFieldBudget(cpu_seconds=1, kernel_values_max=100)
    for count in (-1, .5, True):
        with pytest.raises(ValueError, match="nonnegative integer"):
            budget.check(count)
    monkeypatch.setattr(response.time, "process_time", lambda: budget.started+2)
    with pytest.raises(RuntimeError, match="budget"):
        budget.check()
    for kwargs in ({"nodes": 129.0}, {"nodes": True}, {"nodes": 1025},
                   {"geometry_knots": 129.0}, {"geometry_knots": 513}, {"budget": {}}):
        with pytest.raises(ValueError):
            response.solve_prepared_field(**kwargs)
