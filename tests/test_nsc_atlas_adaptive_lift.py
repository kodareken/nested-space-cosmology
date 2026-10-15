"""Analytic fixed-stencil and provenance controls; no physical preparation."""
from types import SimpleNamespace
import time

import numpy as np
from numpy.polynomial import Polynomial
import pytest

from recursive_horizons.nsc_atlas_adaptive_lift import AdaptiveBoundaryLift, BLOCKS


@pytest.fixture(autouse=True)
def analytic_budget():
    started = time.process_time()
    yield
    assert time.process_time()-started < 30


class PolynomialJets:
    """Ducktyped provider with analytic jets and a nonlinear reference wall."""
    def __init__(self, degree=9):
        self.domain = dict(U0=-.32, Ua=-.30, Ub=-.27, Va=-2., Vb=-1.)
        self.cache = {}
        self.budget = SimpleNamespace(check=lambda: None)
        self.polynomial = Polynomial([0., 0., .3, -.2, .1, .4, -.7, .2, .1, -.15][:degree+1])
        def wall(V):
            t = np.asarray(V)+2
            return self.domain["Ua"]+.03*(t+t*t)/2
        self.geometry = {"wall": wall}

    def sample(self, t):
        if t not in self.cache:
            V = -2+t
            U = float(self.geometry["wall"](V))
            z, Uv = (U+.30)/.03, .03*(1+2*t)/2
            b, slope = float(self.polynomial(z)), float(self.polynomial.deriv()(z)/.03)
            self.cache[t] = dict(t=t, U=U, V=V, b=b, b_U=slope, b_V=slope*Uv, U_V=Uv)
        return dict(self.cache[t])


def knot_times():
    # Refine only selected cells while keeping every original source edge.
    return sorted(set(np.linspace(0., 1., 17)) | {9/32, 23/32})


@pytest.mark.parametrize("nodes,degree", ((4, 7), (5, 9)))
def test_polynomial_reproduction_and_exact_integrated_lift(nodes, degree):
    jets = PolynomialJets(degree)
    knots = knot_times()
    withheld = [(a+b)/2 for a, b in zip(knots[:-1], knots[1:])]
    lift = AdaptiveBoundaryLift(jets, knot_times=knots, nodes=nodes, validation_times=withheld)
    U = np.linspace(-.32, -.27, 101)
    z = np.maximum((U+.30)/.03, 0.)
    assert lift.b_at_U(U) == pytest.approx(jets.polynomial(z), abs=2e-13)
    assert .03*lift.b_U(U) == pytest.approx(jets.polynomial.deriv()(z), abs=2e-12)
    assert lift.B(U) == pytest.approx(.03*jets.polynomial.integ()(z), abs=3e-15)
    V = np.linspace(-2., -1., 41)
    assert lift.b(V) == pytest.approx(lift.b_at_U(jets.geometry["wall"](V)), abs=1e-14)
    assert lift.B_at_wall(V) == pytest.approx(lift.B(jets.geometry["wall"](V)), abs=1e-15)
    assert max(lift.summary()["maximum_absolute_defects"].values()) < 2e-12
    for index, piece in enumerate(lift._pieces):
        width = lift.cells[index]["width"]
        eta = np.linspace(0., 1., 11)
        assert lift._primitives[index].deriv()(eta)/width == pytest.approx(piece(eta), abs=2e-15)


def test_c1_boundary_and_c2_lift_across_changing_stencils_and_source_edges():
    jets = PolynomialJets()
    # Smooth nonpolynomial data tests continuity without relying on reproduction.
    original_sample = jets.sample
    def sample(t):
        row = original_sample(t)
        z = (row["U"]+.30)/.03
        row["b"] = z*z*np.sin(4*z)
        row["b_U"] = (2*z*np.sin(4*z)+4*z*z*np.cos(4*z))/.03
        row["b_V"] = row["b_U"]*row["U_V"]
        jets.cache[t] = dict(row)
        return row
    jets.sample = sample
    lift = AdaptiveBoundaryLift(jets, knot_times=knot_times(), nodes=5)
    changes = 0
    for i in range(1, len(lift._pieces)):
        before, after = lift._pieces[i-1], lift._pieces[i]
        wl, wr = lift.cells[i-1]["width"], lift.cells[i]["width"]
        assert before(1.) == pytest.approx(after(0.), abs=4e-13)
        assert before.deriv()(1.)/wl == pytest.approx(after.deriv()(0.)/wr, abs=2e-10)
        assert lift.prefix[i-1]+lift._primitives[i-1](1.) == pytest.approx(lift.prefix[i], abs=3e-17)
        assert lift._primitives[i-1].deriv()(1.)/wl == pytest.approx(lift._primitives[i].deriv()(0.)/wr, abs=4e-13)
        assert lift._primitives[i-1].deriv(2)(1.)/wl**2 == pytest.approx(lift._primitives[i].deriv(2)(0.)/wr**2, abs=2e-10)
        changes += lift.cells[i-1]["training_indices"] != lift.cells[i]["training_indices"]
    assert changes >= 4
    assert lift.summary()["quiet_anchor"] == {"b": 0., "b_U": 0., "B": 0.}
    assert np.array_equal(lift.B([-.32, -.31, -.30]), np.zeros(3))
    assert np.array_equal(lift.b_at_U([-.32, -.31, -.30]), np.zeros(3))
    assert np.array_equal(lift.b_U([-.32, -.31, -.30]), np.zeros(3))


def test_fixed_stencils_source_edges_and_heldout_values_never_enter_fit():
    jets = PolynomialJets()
    knots = knot_times()
    lift = AdaptiveBoundaryLift(jets, knot_times=knots)
    coefficients = [list(cell["coefficients_in_cell_coordinate"]) for cell in lift.cells]
    for cell in lift.cells:
        lo, hi = BLOCKS[cell["block"]]
        assert all(lo <= t <= hi for t in cell["training_times"])
        assert cell["index"] in cell["training_indices"]
        assert cell["index"]+1 in cell["training_indices"]
    withheld = [(a+b)/2 for a, b in zip(knots[:-1], knots[1:])]
    for t in withheld:
        jets.sample(t)
        jets.cache[t]["b"] += 100.
        jets.cache[t]["b_U"] -= 1000.
    lift.validate(withheld)
    assert coefficients == [cell["coefficients_in_cell_coordinate"] for cell in lift.cells]
    assert all(not set(cell["training_times"]) & set(withheld) for cell in lift.cells)
    assert lift.summary()["maximum_absolute_defects"]["b_defect"] > 99.
    for cell in lift.cells:
        queries = np.linspace(cell["U_left"], cell["U_right"], 31)[1:-1]
        eta = (queries-cell["U_left"])/cell["width"]
        assert lift.b_at_U(queries) == pytest.approx(Polynomial(cell["coefficients_in_cell_coordinate"])(eta), abs=1e-14)
    count = len(jets.cache)
    lift.validate(withheld)
    assert len(jets.cache) == count


@pytest.mark.parametrize("knots,nodes", (([0., .25, .75, 1.], 4),
    ([t for t in np.linspace(0., 1., 17) if t != .25], 4),
    (list(reversed(np.linspace(0., 1., 17))), 4),
    (list(np.linspace(0., 1., 17))+[1.], 4),
    (np.linspace(0., 1., 17), True), (np.linspace(0., 1., 17), 4.),
    (np.linspace(0., 1., 17), 6), (None, 4),
    ([*np.linspace(0., 1., 17)[:-1], True], 4)))
def test_bad_training_inputs_rejected_before_direct_queries(knots, nodes):
    jets = PolynomialJets()
    with pytest.raises(ValueError):
        AdaptiveBoundaryLift(jets, knot_times=knots, nodes=nodes)
    assert not jets.cache


def test_shared_cache_limit_is_preflighted_and_validation_roles_checked():
    jets = PolynomialJets()
    knots = np.linspace(0., 1., 65)
    withheld = (knots[:-1]+knots[1:])/2
    lift = AdaptiveBoundaryLift(jets, knot_times=knots, validation_times=withheld)
    assert len(jets.cache) == 129
    with pytest.raises(ValueError, match="129"):
        lift.validate([1/257])
    assert len(jets.cache) == 129
    with pytest.raises(ValueError, match="fitted knot"):
        lift.validate([.25])
    fresh = PolynomialJets()
    with pytest.raises(ValueError, match="129"):
        AdaptiveBoundaryLift(fresh, knot_times=np.linspace(0., 1., 129), validation_times=[1/257])
    assert not fresh.cache


def test_invalid_domains_anchors_queries_and_budget():
    jets = PolynomialJets()
    lift = AdaptiveBoundaryLift(jets, knot_times=np.linspace(0., 1., 17))
    for bad in (True, np.nan, np.zeros((2, 2)), [], [-.33], [-.269], np.ones(514)):
        with pytest.raises(ValueError):
            lift.B(bad)
    for bad in (-2.01, -.99, np.inf):
        with pytest.raises(ValueError):
            lift.b(bad)
    jets.geometry["wall"] = lambda V: np.full_like(np.asarray(V), 1., dtype=float)
    with pytest.raises(ValueError, match="wall map"):
        lift.B_at_wall(-1.5)
    bad = PolynomialJets()
    bad.sample(0.)
    bad.cache[0.]["b_U"] = 1.
    with pytest.raises(ValueError, match="quiet"):
        AdaptiveBoundaryLift(bad, knot_times=np.linspace(0., 1., 17))
    def expired():
        raise RuntimeError("shared budget expired")
    jets.budget.check = expired
    with pytest.raises(RuntimeError, match="budget"):
        lift.B(-.29)
