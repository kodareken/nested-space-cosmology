"""Analytic profile moments and manufactured boundary jets; no field solve."""
import time
from types import SimpleNamespace

import numpy as np
import pytest
import sympy as sp
from scipy.interpolate import CubicHermiteSpline

from recursive_horizons import nsc_atlas_profile_moment as model
from recursive_horizons.nsc_atlas_bulk_response import BulkBudget
from recursive_horizons.nsc_atlas_cap_response import radial_metric
from recursive_horizons.nsc_atlas_cap_wall_response import reference_coefficients
from recursive_horizons.nsc_atlas_prepared_response import PreparedFieldBudget
from recursive_horizons.nsc_atlas_scalar_boundary import ScalarBoundary
from recursive_horizons.nsc_atlas_scalar_lift import BoundaryJet, BoundaryLift


@pytest.fixture(autouse=True)
def analytic_budget():
    started = time.process_time()
    yield
    assert time.process_time()-started < 30


class QuarticCap:
    def __init__(self, scale=1., arrival_knots=5):
        self.scale = scale
        self.budget = PreparedFieldBudget()
        self.domain = dict(U0=-3., Ua=-2., Ub=-1., Va=-2., Vb=-1., qa=0., qb=.33,
                           Ta=-np.arctanh(.45), Tb=.33-np.arctanh(.55),
                           qmin=-2*np.arctanh(.45), xa=.45, xb=.55)
        V = np.linspace(-2., -1., arrival_knots)
        self.geometry = {"arrival": CubicHermiteSpline(V, .33*(V+2), np.full(V.size, .33)),
                         "wall": CubicHermiteSpline([-2., -1.], [-2., -1.], [1., 1.])}
        self.prepared_field = SimpleNamespace(budget=self.budget,
            dense_data={"cuts": np.array([0., .25, .43, .75, 1.])})
        self.edges = np.array([0., .11, .22, .33])
        self.quiet_certificate = {"quiet_cap": True}
        self.calls = []

    def profile(self, q, *, second=False):
        q = np.asarray(q)
        assert second is False
        assert 1 <= q.size <= 513
        assert np.all(q >= self.domain["qmin"]) and np.all(q <= self.domain["qb"])
        self.calls.append(q.copy())
        self.budget.check(17*q.size)
        active = np.maximum(q, 0.)
        return {"f": self.scale*active**4, "f_prime": self.scale*4*active**3}

    def radial_limit(self, T):
        T = np.asarray(T)
        if np.any(T < self.domain["Ta"]-2e-15) or np.any(T >= self.domain["qa"]):
            raise ValueError("manufactured cap time")
        return np.full_like(T, .55)

    def jet(self, T, r):
        T, r = np.broadcast_arrays(T, r)
        q = T+np.arctanh(r)
        assert np.all(q <= self.domain["qb"]+8*np.finfo(float).eps)
        q = np.clip(q, self.domain["qmin"], self.domain["qb"])
        p = self.profile(q)
        f, fp = p["f"], p["f_prime"]
        value = np.divide(f, r, out=np.zeros_like(r), where=r != 0)
        return {"p": value,
                "p_T": np.divide(fp, r, out=np.zeros_like(r), where=r != 0),
                "p_r": np.divide(fp/(1-r*r)-value, r, out=np.zeros_like(r), where=r != 0),
                "p_center": np.zeros_like(r)}


def prepared_moment(order=4, scale=1., arrival_knots=5):
    cap = QuarticCap(scale, arrival_knots)
    budget = BulkBudget(started=cap.budget.started, external=cap.budget)
    return model.ProfileMoment(cap, order, budget=budget)


def test_quartic_profile_moment_and_all_required_cell_edges_are_exact():
    moment = prepared_moment()
    q = np.r_[-.5, 0., np.linspace(.001, .33, 37)]
    assert moment.value(q) == pytest.approx(8*np.maximum(q, 0.)**7/7, abs=3e-19)
    assert np.all(np.diff(moment.prefix) >= 0)
    cap = moment.cap_field
    mapped = cap.geometry["arrival"](-2+cap.prepared_field.dense_data["cuts"])
    for value in np.r_[cap.edges, mapped, cap.geometry["arrival"](cap.geometry["arrival"].x)]:
        assert np.min(abs(moment.cuts-value)) < 1e-16
    fine = prepared_moment(8)
    assert fine.value(q) == pytest.approx(moment.value(q), abs=3e-19)


def test_mass_and_first_jets_match_analytic_profile_and_independent_radial_metric():
    moment = prepared_moment()
    T = np.array([-.45, -.4, -.35])
    r = np.array([.5, .48, .49])
    row = moment.mass(T, r)
    q, A = T+np.arctanh(r), 1-r*r
    expected = 8*q**7/7-A*q**8/(4*r)
    assert row["nu_E"] == pytest.approx(expected, abs=3e-19)
    assert row["p"] == pytest.approx(q**4/r, abs=2e-18)
    assert row["p_T"] == pytest.approx(4*q**3/r, abs=2e-17)
    assert row["p_r"] == pytest.approx(4*q**3/(r*A)-q**4/r**2, abs=2e-17)
    direct = radial_metric(moment.cap_field, T, r, radial_order=16, budget=moment.budget)
    assert row["nu_E"] == pytest.approx(direct["nu_E"], abs=5e-18)
    assert np.all(row["p_center"] == 0)
    # Potential energy is part of this equality; the kinetic term alone fails.
    assert np.max(abs(row["nu_E"]-direct["mass_kinetic"])) > 1e-9


def test_profile_batches_partial_cache_quiet_extension_and_counts():
    moment = prepared_moment(arrival_knots=257)
    assert moment.summary()["prefix_nodes"] > 513
    assert max(len(values) for values in moment.cap_field.calls) <= 513
    queries = np.array([.031234, .031234, .187654, .291234])
    before = moment.partial_nodes
    result = moment.value(queries)
    assert moment.partial_nodes-before == 3*moment.order
    profile_count = moment.profile_points
    assert moment.value(queries) == pytest.approx(result, abs=0)
    assert moment.profile_points == profile_count
    assert moment.cache_hits >= 3
    assert np.all(moment.value([-.7, -.1, 0.]) == 0)
    assert moment.profile_points == profile_count
    report = moment.summary()
    assert report["settings"] == {"moment_order": 4}
    assert report["cached_values"] == 3
    assert report["terminal_C"] == pytest.approx(8*.33**7/7, abs=3e-19)
    assert report["schema"] == model.SCHEMA


def test_partial_cache_has_a_fixed_storage_bound(monkeypatch):
    moment = prepared_moment()
    monkeypatch.setattr(model, "MAX_CACHE", 3)
    moment.value([.001, .002, .003, .004, .005])
    assert len(moment.cache) == 3
    assert list(moment.cache) == [.003, .004, .005]


def test_arbitrarily_small_positive_radius_is_exactly_quiet_without_cancellation():
    moment = prepared_moment()
    row = moment.mass(-.4, np.array([1e-300, 1e-150, 1e-12]))
    for key in ("nu_E", "p", "p_T", "p_r", "C"):
        assert np.all(row[key] == 0)


def test_joint_moving_ray_derivative_matches_differentiated_fixed_interval():
    moment = prepared_moment()
    s, z = sp.symbols("s z", real=True)
    u, upper = -sp.Rational(9, 10)+3*s/100, sp.Rational(1, 5)+s/25
    q = upper*z
    r = sp.tanh((q-u)/2)
    A, f, fp, C = 1-r*r, q**4, 4*q**3, 8*q**7/7
    K = A*C/(2*r**3)-A*f*fp/(12*r**3)-A*f*f/(12*r**4)
    fixed = -K*f*upper
    value = sp.lambdify((s, z), fixed, "numpy")
    derivative = sp.lambdify((s, z), sp.diff(fixed, s), "numpy")
    nodes, weights = np.polynomial.legendre.leggauss(64)
    points, weights = (nodes+1)/2, weights/2
    location = .4
    row = model.outgoing_moment_jet(moment, float(u.subs(s, location)), float(upper.subs(s, location)))
    total = -.04*row["K1_wall"]*row["wall_profile"]-.03*row["integral_Ku_f"]
    assert row["J_minus"] == pytest.approx(float(weights @ value(location, points)), rel=5e-13, abs=1e-20)
    assert total == pytest.approx(float(weights @ derivative(location, points)), rel=8e-13, abs=1e-20)
    assert row["J_minus_u"] == -row["integral_Ku_f"]


def make_boundary(scale=0.):
    cap = QuarticCap(scale)
    budget = BulkBudget(started=cap.budget.started, external=cap.budget)
    def wall_sample(t):
        x, V = .45+.1*t, -2+t
        ref = {key: float(value) for key, value in reference_coefficients(x).items()}
        Vs = -V/(1.6*(ref["v"]-ref["beta_plus"]))
        q = .33*t
        w, j = scale*q**4, scale*4*q**3*.33*Vs
        p = w/x
        return {**ref, "V": V, "Vdot": Vs, "w": w, "j": j, "p": p,
                "A": (j-ref["v"]*p)/x,
                **dict.fromkeys(("rho", "rho_s", "V1", "V1_s"), 0.)}
    wall = SimpleNamespace(budget=budget, sample=wall_sample, geometry=cap.domain,
                           wedge=SimpleNamespace(prepared_field=cap.prepared_field))
    return ScalarBoundary(cap, wall, radial_order=16, ray_order=16, budget=budget)


def test_quiet_boundary_lift_and_shared_cache_are_compatible():
    boundary = make_boundary()
    jets = model.MomentBoundaryJet(boundary)
    assert jets.settings["moment_order"] == 4
    assert "moment_order" not in boundary.settings
    lift = BoundaryLift(jets, knots=5, validation_times=(.125, .375, .625, .875))
    for row in jets.cache.values():
        for key in ("J_minus", "J_minus_s", "b", "b_U", "b_V", "c1_r", "c1_r_s"):
            assert row[key] == 0
    assert lift.summary()["terminal_lift"] == 0
    assert len(jets.cache) == 9
    before = jets.moment.summary()
    jets.sample(.5)
    assert jets.moment.summary() == before


def test_moment_boundary_matches_original_radial_jet_on_nonzero_manufactured_profile():
    boundary = make_boundary(1.)
    old = BoundaryJet(boundary)
    new = model.MomentBoundaryJet(boundary)
    for t in (.3, .65, 1.):
        reference, row = old.sample(t), new.sample(t)
        assert set(row) == set(reference)
        for key in ("J_minus", "integral_Ku_f", "J_minus_s", "b", "b_U", "b_V", "c1_r"):
            assert row[key] == pytest.approx(reference[key], rel=5e-12, abs=1e-18)


@pytest.mark.parametrize("order", [True, 4., 2, 16])
def test_moment_order_rejected(order):
    with pytest.raises(ValueError, match="order"):
        prepared_moment(order)


def test_domains_shared_budget_and_quiet_certificate_fail_closed():
    moment = prepared_moment()
    for labels in (True, np.nan, -.99, .331, np.zeros((2, 2)), np.zeros(514)):
        with pytest.raises(ValueError):
            moment.value(labels)
    for T, r in ((0., .5), (-1., .5), (-.4, 0.), (-.4, 1.), (-.4, .56)):
        with pytest.raises(ValueError):
            moment.mass(T, r)
    cap = QuarticCap()
    with pytest.raises(ValueError, match="shared"):
        model.ProfileMoment(cap, budget=BulkBudget())
    changed_start = BulkBudget(started=cap.budget.started-1, external=cap.budget)
    with pytest.raises(ValueError, match="clock"):
        model.ProfileMoment(cap, budget=changed_start)
    cap.quiet_certificate["quiet_cap"] = False
    with pytest.raises(ValueError, match="quiet"):
        model.ProfileMoment(cap)


def test_invalid_profile_and_nonmonotone_cuts_are_rejected():
    cap = QuarticCap()
    cap.prepared_field.dense_data["cuts"] = np.array([0., .5, .4, 1.])
    with pytest.raises(ValueError, match="partition"):
        model.ProfileMoment(cap)
    cap = QuarticCap()
    cap.profile = lambda q, **kwargs: {"f": np.zeros_like(q), "f_prime": np.full_like(q, np.nan)}
    with pytest.raises(ValueError, match="profile"):
        model.ProfileMoment(cap)


def test_moment_boundary_rejects_another_profile_or_budget():
    boundary = make_boundary()
    with pytest.raises(ValueError, match="same profile"):
        model.MomentBoundaryJet(boundary, moment=prepared_moment())
    moment = model.ProfileMoment(boundary.cap_field, budget=boundary.budget)
    with pytest.raises(ValueError, match="same profile"):
        model.MomentBoundaryJet(boundary, moment=moment, moment_order=8)


def test_expired_budget_rejects_cached_queries():
    moment = prepared_moment()
    moment.value(.123)
    moment.budget.started -= 121
    moment.cap_field.budget.started = moment.budget.started
    with pytest.raises(RuntimeError, match="budget"):
        moment.value(.123)
