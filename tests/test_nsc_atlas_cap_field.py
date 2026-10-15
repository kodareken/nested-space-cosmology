from types import SimpleNamespace

import numpy as np
import pytest
from numpy.polynomial import Polynomial

from recursive_horizons import nsc_atlas_cap_field as cap_module
from recursive_horizons.nsc_atlas_cap_field import CapLeadingField, source_second_derivative
from recursive_horizons.nsc_atlas_cap_wall_response import field_mass_ingredients, reference_coefficients
from recursive_horizons.nsc_atlas_field_preparation import preparation_record, radius_from_null, source_profile
from recursive_horizons.nsc_atlas_prepared_response import PreparedFieldBudget, _geometry, solve_prepared_field


class PolynomialCap(CapLeadingField):
    """Independent analytic profile; no prepared-source physical claim."""
    def __init__(self, degree, order=8, centre_radius=1e-4):
        self.degree = degree
        self.centre_radius = centre_radius
        self.domain = dict(Ta=-1., Tb=-.6, qa=-.4, qb=.6, qmin=-1.5)
        self.edges = np.array([-1.5, .6])
        self.budget = PreparedFieldBudget()
        self._centre_gauss = np.polynomial.legendre.leggauss(order)

    def profile(self, q, *, second=False):
        q = np.asarray(q)
        result = {"f": q**self.degree, "f_prime": self.degree*q**(self.degree-1)}
        if second:
            result["f_second"] = self.degree*(self.degree-1)*q**(self.degree-2)
        return result

    def radial_limit(self, T):
        return np.full_like(np.asarray(T), .5)


@pytest.mark.parametrize("order", [8, 16])
@pytest.mark.parametrize("degree", [2, 3])
@pytest.mark.parametrize("centre_radius", [1e-4, 1e-3])
def test_independent_polynomial_jets_at_radius_one_e_minus_eight(degree, order, centre_radius):
    cap = PolynomialCap(degree, order, centre_radius)
    T, r = .13, 1e-8
    row = cap.jet(T, r)
    if degree == 2:
        expected_p, expected_T, expected_r = 4*T, 4., 8*T*r/3
    else:
        expected_p, expected_T, expected_r = 6*T*T, 12*T, 4*(T*T+1)*r
    assert float(row["p"]) == pytest.approx(expected_p, rel=3e-15, abs=3e-15)
    assert float(row["p_T"]) == pytest.approx(expected_T, rel=3e-15, abs=3e-15)
    assert float(row["p_r"]) == pytest.approx(expected_r, rel=2e-7, abs=2e-15)


def test_polynomial_reflected_field_and_exact_centre():
    cap = PolynomialCap(3)
    T = np.array([0., .1, .13])
    r = np.array([.2, .12, 0.])
    row = cap.jet(T, r)
    d = np.arctanh(r[:2])
    expected = (6*T[:2]**2*d+2*d**3)/r[:2]
    np.testing.assert_allclose(row["p"][:2], expected, rtol=2e-14)
    np.testing.assert_allclose(row["p_T"][:2], 12*T[:2]*d/r[:2], rtol=2e-14)
    assert row["p"][0] == pytest.approx(2*np.arctanh(.2)**3/.2)
    assert row["p"][2] == pytest.approx(6*T[2]**2)
    assert row["p_T"][2] == pytest.approx(12*T[2])
    assert row["p_r"][2] == 0


def test_small_radius_quadrature_splits_profile_edges():
    cap = PolynomialCap(3, centre_radius=1e-3)
    cap.edges = np.array([0.])

    def profile(q, *, second=False):
        q = np.asarray(q)
        result = {"f": np.maximum(q, 0.)**4, "f_prime": 4*np.maximum(q, 0.)**3}
        if second:
            result["f_second"] = 12*np.maximum(q, 0.)**2
        return result

    cap.profile = profile
    r = 3e-4
    d = np.arctanh(r)
    row = cap.jet(0., r)
    assert float(row["p"]) == pytest.approx(d**4/r, rel=3e-14)
    assert float(row["p_T"]) == pytest.approx(4*d**3/r, rel=3e-14)
    assert float(row["p_r"]) == pytest.approx(4*d**3/(r*(1-r*r))-d**4/r**2, rel=3e-14)
    assert float(row["p_center"]) == 0


@pytest.mark.parametrize("threshold", [1e-4, 1e-3])
def test_polynomial_jets_on_both_sides_of_centre_threshold(threshold):
    cap = PolynomialCap(3, centre_radius=threshold)
    r = np.array([np.nextafter(threshold, 0.), threshold, np.nextafter(threshold, np.inf)])
    T = .13
    row = cap.jet(T, r)
    ratio = np.arctanh(r)/r
    # Independent closed form for f(q)=q^3, evaluated without cancellation.
    radius_derivative = r*(2/3+r*r*(4/5+r*r*(6/7+r*r*8/9)))
    expected_p = 6*T*T*ratio+2*r*r*ratio**3
    expected_T = 12*T*ratio
    expected_r = 6*T*T*radius_derivative+2*r*ratio**2*(3/(1-r*r)-ratio)
    np.testing.assert_allclose(row["p"], expected_p, rtol=4e-13, atol=4e-15)
    np.testing.assert_allclose(row["p_T"], expected_T, rtol=4e-13, atol=4e-13)
    np.testing.assert_allclose(row["p_r"], expected_r, rtol=1e-7, atol=2e-11)


@pytest.fixture
def fake_prepared():
    prepared, geometry = preparation_record(), _geometry(129)
    domain = {key: prepared["geometry"][key] for key in ("U0", "Ua", "Ub", "Va", "Vb")}
    grid = np.linspace(domain["U0"], domain["Ub"], 129)
    budget = PreparedFieldBudget()

    def wall(t):
        t = np.asarray(t)
        V = domain["Va"]+(domain["Vb"]-domain["Va"])*t
        g, gp = source_profile(V, domain["Va"], domain["Vb"])
        budget.check(len(grid)*t.size)
        return {"wall_chi_over_epsilon": g, "wall_dchi_dV_over_epsilon": gp}

    return SimpleNamespace(geometry=geometry, domain=domain, budget=budget, wall=wall,
        dense_data={"grid": grid, "cuts": np.array([0., 1.]), "prepared": prepared,
                    "segments": (lambda t: np.zeros((len(grid),)+np.asarray(t).shape),)})


def test_analytic_compact_source_second_derivative():
    B = Polynomial([1., 0., -4., 0., 6., 0., -4., 0., 1.])
    Va, Vb = -2.5, -2.2
    z = np.array([-.9, -.2, 0., .37, .9])
    V = (Va+Vb)/2+z*(Vb-Va)/4
    np.testing.assert_allclose(source_second_derivative(V, Va, Vb), B.deriv(2)(z)*(4/(Vb-Va))**2,
                               rtol=3e-13, atol=1e-10)
    assert np.all(source_second_derivative(np.array([Va, Vb]), Va, Vb) == 0)


def test_profile_second_derivative_zero_kernel_control(fake_prepared, monkeypatch):
    monkeypatch.setattr(cap_module, "exterior_coefficient", lambda x: np.zeros_like(x))
    cap = CapLeadingField(fake_prepared, derivative_order=4)
    V = cap.Va+cap.dV*np.array([.31, .44, .62])
    q = cap.geometry["arrival"](V)
    result = cap.profile(q, second=True)
    _, gp = source_profile(V, cap.Va, cap.Vb)
    qV, qVV = cap.geometry["arrival"](V, 1), cap.geometry["arrival"](V, 2)
    expected = (source_second_derivative(V, cap.Va, cap.Vb)*qV-gp*qVV)/qV**3
    np.testing.assert_allclose(result["f_second"], expected, rtol=2e-11, atol=2e-10)
    calls = cap.evaluation_counts["wall_points"]
    count = cap.budget.kernel_values
    again = cap.profile(np.r_[q, q], second=True)
    np.testing.assert_array_equal(again["f_second"], np.tile(result["f_second"], 2))
    assert cap.evaluation_counts["wall_points"] == calls
    assert cap.budget.kernel_values == count


def test_quiet_interval_and_finite_cap_domain(fake_prepared):
    cap = CapLeadingField(fake_prepared)
    q = np.linspace(cap.domain["qmin"], cap.domain["qa"], 17)
    result = cap.profile(q, second=True)
    assert all(np.all(value == 0) for value in result.values())
    row = cap.jet(cap.domain["Ta"], np.linspace(0., cap.domain["xa"], 9))
    assert all(np.all(value == 0) for value in row.values())
    assert float(cap.radial_limit(cap.domain["qb"])) == 0
    assert float(cap.radial_limit(cap.domain["Tb"])) == pytest.approx(cap.domain["xb"], abs=2e-15)
    for qbad in (cap.domain["qmin"]-1e-10, cap.domain["qb"]+1e-10):
        with pytest.raises(ValueError, match="supplied"):
            cap.profile(qbad)
    for T, r in ((cap.domain["Ta"]-1e-8, 0.), (cap.domain["qb"]+1e-8, 0.),
                 (cap.domain["qb"], .01), (cap.domain["Ta"], cap.domain["xa"]+.01),
                 (cap.domain["Ta"], -.1)):
        with pytest.raises(ValueError, match="domain"):
            cap.jet(T, r)


def test_cap_time_allows_only_endpoint_roundoff(fake_prepared):
    cap = CapLeadingField(fake_prepared)
    Ta = cap.domain["Ta"]
    rounded = np.nextafter(Ta, -np.inf)
    assert float(cap.radial_limit(rounded)) == float(cap.radial_limit(Ta))
    row = cap.jet(rounded, cap.domain["xa"])
    assert all(float(value) == 0 for value in row.values())


@pytest.mark.parametrize("value", [[], [[-.4]], [np.nan], [np.inf], True, np.full(514, -.4)])
def test_bad_profile_inputs(fake_prepared, value):
    with pytest.raises(ValueError):
        CapLeadingField(fake_prepared).profile(value)


@pytest.mark.parametrize("kwargs", [{"derivative_order": 3}, {"derivative_order": True},
                                    {"small_radius_order": 4}, {"small_radius_order": 8.},
                                    {"centre_radius": True}, {"centre_radius": 0.},
                                    {"centre_radius": 1e-5}, {"centre_radius": float("nan")},
                                    {"centre_radius": "1e-4"}])
def test_unadmitted_orders(fake_prepared, kwargs):
    with pytest.raises(ValueError):
        CapLeadingField(fake_prepared, **kwargs)


def test_shared_budget_clock_still_checked_on_cache_hits(fake_prepared):
    cap = CapLeadingField(fake_prepared)
    q = np.mean([cap.domain["qa"], cap.domain["qb"]])
    cap.profile(q)
    cap.budget.kernel_values = cap.budget.kernel_values_max+1
    with pytest.raises(RuntimeError, match="budget"):
        cap.profile(q)


@pytest.mark.parametrize("kind", ["margin", "quiet"])
def test_rejects_missing_quiet_or_no_return_certificate(fake_prepared, kind):
    certificate = fake_prepared.dense_data["prepared"]
    if kind == "margin":
        certificate["observer"]["analytic_no_return_log_margin_lower_bound"] = -1.
    else:
        certificate["preparation"]["cap_data"] = "unspecified"
    with pytest.raises(ValueError, match="certificate"):
        CapLeadingField(fake_prepared)


def test_existing_transmitted_wall_first_jets_and_regular_centre():
    budget = PreparedFieldBudget()
    prepared = solve_prepared_field(129, 129, budget)
    cap = CapLeadingField(prepared)
    t = np.array([0., .25, .37, .5, .64, .75, 1.])
    wall = prepared.wall(t)
    V = wall["V"]
    x = radius_from_null(prepared.geometry["wall"](V), V)
    q = prepared.geometry["arrival"](V)
    row = cap.jet(q-np.arctanh(x), x)
    admitted_x = np.clip(x, .45, .55)
    g = reference_coefficients(admitted_x)
    w = wall["wall_chi_over_epsilon"]
    j = prepared.geometry["Vdot"](V)*wall["wall_dchi_dV_over_epsilon"]
    expected = field_mass_ingredients(admitted_x, w, j)
    A = g["beta_minus"]/(1-x*x)*row["p_T"]+g["v"]*row["p_r"]
    B = g["v"]/(1-x*x)*row["p_T"]+g["beta_minus"]*row["p_r"]
    np.testing.assert_allclose(row["p"], expected["p"], rtol=2e-12, atol=1e-12)
    np.testing.assert_allclose(A, expected["A"], rtol=3e-7, atol=3e-7)
    np.testing.assert_allclose(B, expected["B_minus"], rtol=3e-7, atol=3e-7)
    T = float(np.mean([cap.domain["qa"], cap.domain["qb"]]))
    tiny = cap.jet(T, np.array([0., 1e-8]))
    np.testing.assert_allclose(tiny["p"], tiny["p_center"], rtol=3e-13, atol=3e-13)
    assert np.all(np.isfinite(tiny["p_T"])) and abs(tiny["p_r"][-1]) < 1e-3
    assert cap.evaluation_counts["second_points"] > 0
