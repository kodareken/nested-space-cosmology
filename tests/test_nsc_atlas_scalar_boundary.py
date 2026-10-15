"""Manufactured checks for the unevolved cubic boundary input; no physical PDE."""
import time
from types import SimpleNamespace

import numpy as np
import pytest

from recursive_horizons import nsc_atlas_scalar_boundary as model
from recursive_horizons.nsc_atlas_bulk_response import BulkBudget
from recursive_horizons.nsc_atlas_cap_wall_response import reference_coefficients
from recursive_horizons.nsc_atlas_prepared_response import PreparedFieldBudget


@pytest.fixture(autouse=True)
def manufactured_cpu_budget():
    started = time.process_time()
    yield
    assert time.process_time()-started < 30


class ManufacturedCap:
    """Regular C3 incoming profile on a pre-arrival cap; no dense field solve."""
    def __init__(self, scale=1.):
        self.scale = scale
        self.budget = PreparedFieldBudget()
        self.prepared_field = object()
        self.domain = dict(U0=-3., Ua=-2., Ub=-1., Va=-2., Vb=-1., qa=0., qb=.33,
                           Ta=-np.arctanh(.45), Tb=.33-np.arctanh(.55),
                           qmin=-2*np.arctanh(.45), xa=.45, xb=.55)
        self.geometry = {"arrival": lambda V: .33*(np.asarray(V)+2)}
        self.edges = np.array([0., .11, .22, .33])
        self.quiet_certificate = {"quiet_cap": True, "no_return_margin_over_ell": .5,
                                  "analytic_no_return_log_margin_lower_bound": .5}
        self.evaluation_counts = dict(profile_requested=0, wall_points=0, second_points=0, cache_hits=0)

    def profile(self, q):
        q = np.asarray(q, dtype=float)
        if np.any(q < self.domain["qmin"]-1e-14) or np.any(q > self.domain["qb"]+1e-14):
            raise ValueError("manufactured profile domain")
        self.evaluation_counts["profile_requested"] += q.size
        qp = np.maximum(q, 0.)
        return {"f": self.scale*qp**4, "f_prime": self.scale*4*qp**3}

    def radial_limit(self, T):
        T = np.asarray(T, dtype=float)
        if np.any(T < self.domain["Ta"]-1e-14) or np.any(T > self.domain["Tb"]+1e-14):
            raise ValueError("manufactured cap time domain")
        return np.full_like(T, .55)

    def jet(self, T, r):
        T, r = np.broadcast_arrays(T, r)
        limit = self.radial_limit(T)
        if np.any(r < 0) or np.any(r > limit+1e-14):
            raise ValueError("manufactured cap radius domain")
        profile = self.profile(T+np.arctanh(r))
        f, fp = profile["f"], profile["f_prime"]
        p = np.divide(f, r, out=np.zeros_like(r), where=r != 0)
        pT = np.divide(fp, r, out=np.zeros_like(r), where=r != 0)
        pr = np.divide(fp/(1-r*r)-p, r, out=np.zeros_like(r), where=r != 0)
        return {"p": p, "p_T": pT, "p_r": pr, "p_center": np.zeros_like(r)}


class ManufacturedWall:
    def __init__(self, cap):
        self.cap = cap
        self.geometry = {key: cap.domain[key] for key in model.DOMAIN_KEYS[:5]}
        self.budget = BulkBudget(started=cap.budget.started, external=cap.budget)
        self.wedge = SimpleNamespace(prepared_field=cap.prepared_field)
        self.calls = 0

    def sample(self, t):
        self.calls += 1
        x, V = .45+.1*t, -2+t
        g = {key: float(value) for key, value in reference_coefficients(x).items()}
        w = float(self.cap.profile(np.array([.33*t]))["f"][0])
        return {**g, "V": V, "w": w, "p": w/x,
                "rho": .1*t*self.cap.scale**2, "V1": .02*t*self.cap.scale**2,
                "Vdot": -V/(1.6*(g["v"]-g["beta_plus"]))}


def adapter(scale=1., **kwargs):
    cap = ManufacturedCap(scale)
    wall = ManufacturedWall(cap)
    return model.ScalarBoundary(cap, wall, **kwargs)


def test_zero_source_is_exactly_zero_and_default_summary_has_nine_cuts():
    boundary = adapter(0., radial_order=4, ray_order=4)
    result = boundary.summary()
    assert result["schema"] == model.SCHEMA
    assert len(result["boundary_samples"]) == result["boundary_evaluations"] == 9
    for row in result["boundary_samples"]:
        for key in ("J_minus", "J_minus_r", "c1_r", "backscatter_term", "displacement_term",
                    "positive_integral", "negative_integral", "absolute_integral"):
            assert row[key] == 0
    assert all(value == 0 for value in result["maximum_absolute_defects"].values())
    assert result["minimum_domain_margins"]["active_radius_lower"] > .4
    assert result["scalar_kernel_values"] == 0
    assert result["profile_counts"]["second_points"] == 0


def test_cubic_boundary_scaling_includes_metric_and_wall_pullback():
    first, second = adapter(1., radial_order=16, ray_order=16), adapter(-2., radial_order=16, ray_order=16)
    a, b = first.sample(.8), second.sample(.8)
    assert abs(a["J_minus"]) > 1e-12
    assert abs(a["displacement_term"]) > 1e-8
    for key in ("J_minus", "J_minus_r", "c1_r", "backscatter_term", "displacement_term"):
        assert b[key] == pytest.approx(-8*a[key], rel=3e-13, abs=1e-18)
    potential_scale = max(abs(a["K1_sample_min"]), abs(a["K1_sample_max"]))
    assert a["potential_identity_defect"] <= 64*np.finfo(float).eps*potential_scale
    assert a["wall_flux_defect"] < 1e-18


def test_q_and_r_quadratures_agree_on_independently_placed_nodes():
    boundary = adapter(radial_order=16, ray_order=16)
    row = boundary.sample(1.)
    assert row["J_minus"] == pytest.approx(row["J_minus_r"], rel=3e-12, abs=1e-18)
    assert row["positive_integral"]+row["negative_integral"] == row["J_minus"]
    assert row["positive_integral"]-row["negative_integral"] == row["absolute_integral"]
    assert row["positive_integral"] >= 0 >= row["negative_integral"]
    assert row["ray_nodes"] == row["ray_nodes_r"] == 48


def test_forced_null_polynomial_integral_has_the_correct_sign_and_jacobian(monkeypatch):
    cap = ManufacturedCap()
    budget = ManufacturedWall(cap).budget
    # Manufacture a metric coefficient giving K1=1+q+q^2.  This is an
    # independent quadrature test, not a claimed solution of Einstein's law.
    def polynomial_metric(field, T, r, *, radial_order, budget):
        q, A = T+np.arctanh(r), 1-r*r
        jets = field.jet(T, r)
        p, pT, pr = jets["p"], jets["p_T"], jets["p_r"]
        K = 1+q+q*q
        nu = r**3*(4*K/A+p*p/3+pT*pT/(6*A)-A*pr*pr/6)/2
        budget.check(metric=len(r))
        return {**jets, "nu_E": nu}
    monkeypatch.setattr(model, "radial_metric", polynomial_metric)
    upper, u = .27, -.9
    expected = -(upper**5/5+upper**6/6+upper**7/7)
    for coordinate in ("q", "r"):
        result = model.outgoing_cap_flux(cap, u, upper, budget=budget, ray_order=16,
                                         radial_order=4, coordinate=coordinate)
        assert result["J_minus"] == pytest.approx(expected, rel=3e-14, abs=1e-18)
        assert result["positive_integral"] == 0
        assert result["negative_integral"] < 0


def test_beta_plus_zero_uses_a_regular_signed_boundary_factor():
    cap = ManufacturedCap()
    wall = ManufacturedWall(cap)
    x = .16**(1/3)
    row = wall.sample((x-.45)/.1)
    assert abs(row["beta_plus"]) < 2e-15
    response = model.wall_boundary_coefficients(row, .07)
    assert all(np.isfinite(value) for value in response.values())
    assert response["wall_flux_defect"] == pytest.approx(0, abs=3e-17)
    expected = (2*(row["beta_minus"]-row["v"])*.07/((1-x*x)*row["v"])
                -(row["rho"]-x*row["p"]**2/12-row["v"]*row["V1"]/row["Vdot"])
                *.8*row["w"]/(x**3*row["v"]**2))
    assert response["c1_r"] == pytest.approx(expected, abs=3e-16)


def test_profile_conversion_error_is_exposed_by_potential_identity():
    r, nu, f, fp = .5, .2, .1, .3
    p, pT, pr = f/r, fp/r, fp/(r*(1-r*r))-f/r**2
    correct = model.cap_potential(r, nu, p, pT, pr, f, fp)
    wrong = model.cap_potential(r, nu, p, pT, pr+1, f, fp)
    assert float(correct["potential_identity_defect"]) == pytest.approx(0, abs=1e-16)
    assert abs(float(wrong["potential_identity_defect"])) > .01


def test_cache_avoids_repeated_queries_and_returns_an_independent_mapping():
    boundary = adapter(radial_order=4, ray_order=4)
    row = boundary.sample(.5)
    counts = (boundary.budget.metric_values, boundary.cap_field.evaluation_counts.copy(), boundary.wall_solution.calls)
    row["c1_r"] = 123
    assert boundary.sample(.5)["c1_r"] != 123
    assert counts == (boundary.budget.metric_values, boundary.cap_field.evaluation_counts, boundary.wall_solution.calls)
    assert boundary.summary([.5])["boundary_evaluations"] == 1


def test_different_preparations_domains_and_budgets_are_rejected():
    cap, other = ManufacturedCap(), ManufacturedCap()
    wall = ManufacturedWall(cap)
    wall.wedge.prepared_field = other.prepared_field
    with pytest.raises(ValueError, match="same prepared"):
        model.ScalarBoundary(cap, wall)
    wall.wedge.prepared_field = cap.prepared_field
    with pytest.raises(ValueError, match="shared"):
        model.ScalarBoundary(cap, wall, budget=BulkBudget())
    with pytest.raises(ValueError, match="cumulative metric"):
        model.ScalarBoundary(cap, wall, budget=BulkBudget(started=cap.budget.started, external=cap.budget))
    wall.geometry["Va"] += .01
    with pytest.raises(ValueError, match="domains differ"):
        model.ScalarBoundary(cap, wall)


@pytest.mark.parametrize("bad", [True, np.nan, -1., 1.1])
def test_bad_wall_times_fail_before_computation(bad):
    boundary = adapter(radial_order=4, ray_order=4)
    with pytest.raises(ValueError):
        boundary.sample(bad)
    assert boundary.wall_solution.calls == 0


@pytest.mark.parametrize("kwargs", [dict(radial_order=16.), dict(ray_order=True), dict(ray_order=64)])
def test_unadmitted_orders_are_rejected(kwargs):
    with pytest.raises(ValueError, match="order"):
        adapter(**kwargs)


def test_ray_domain_and_exhausted_cumulative_budget_fail_closed():
    cap = ManufacturedCap()
    budget = ManufacturedWall(cap).budget
    for u, q in ((.01, .2), (-1.1, .2), (-.9, .4), (-.9, -.01)):
        with pytest.raises(ValueError, match="domain"):
            model.outgoing_cap_flux(cap, u, q, budget=budget)
    with pytest.raises(ValueError, match="coordinate"):
        model.outgoing_cap_flux(cap, -.9, .2, budget=budget, coordinate="T")
    budget.started -= 121
    cap.budget.started = budget.started
    with pytest.raises(RuntimeError, match="budget"):
        model.outgoing_cap_flux(cap, -.9, .2, budget=budget)
