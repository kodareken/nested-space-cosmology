"""Shared bounded stage/driver controls; no prepared physical trajectory."""
import copy
import time

import numpy as np
import pytest

from recursive_horizons.nsc_atlas_characteristic_initialization import Budget
from recursive_horizons.nsc_atlas_characteristic_background import CharacteristicProfile
from recursive_horizons.nsc_atlas_reference_wall_gauge import ReferenceWallGauge
from recursive_horizons.nsc_atlas_compatible_evolution import quiet_departures
from recursive_horizons.nsc_atlas_background_evolution import background_stage
from recursive_horizons.nsc_atlas_departure_gradient import first_order_departure_gradient
from recursive_horizons.nsc_atlas_gradient_evolution import coupled_stage
from test_nsc_atlas_reference_wall_gauge import GEOMETRY


CONTROL = Budget(maximum_evaluations=300000, cpu_seconds=30.)
TUBE = {"r": [.2, .8], "L": [-2., 2.], "c": [-4., 4.]}


@pytest.fixture(scope="session")
def gradient_control_budget():
    yield CONTROL
    CONTROL.check()
    print("gradient stage/driver controls:", CONTROL.evaluations, "points;",
          time.process_time()-CONTROL.started, "process CPU seconds")


def profile(scale=0.):
    q = np.array([0., .3, .6, 1.3])
    return CharacteristicProfile(q, scale*q**3, 3*scale*q**2, 6*scale*q)


def source(x, slope=0.):
    return dict(r=x, log_C=0., c=slope*(x-.4), r_x=1., log_C_x=0., c_x=slope)


def stage_options(budget, *, n=4, active=False, alpha=0.):
    eta, y = np.linspace(0., 1., n+1), np.linspace(0., 1., 4)
    ref = ReferenceWallGauge(GEOMETRY, budget=budget)
    return dict(x=.4, eta=eta, y_faces=y, source=source(.4, .07 if active else 0.),
        gauge=ref.sample(.4), alpha=alpha, tube=TUBE, budget=budget,
        profile=profile(.025 if active else 0.))


def initial_state(n=4, active=False):
    value = quiet_departures(exterior_cells=n, cap_cells=3)
    if active:
        eta = np.linspace(0., 1., n+1)[1:-1]
        value["departure"]["r"] = .0002*eta*(1-eta)
        value["departure"]["L"] = -.0004*eta*(1-eta)
        value["departure"]["c"] = .001*eta**2
        r = .45*np.linspace(0., 1., 4)
        value["cap"]["Ip"] = .007*np.diff(r*r)/2
        value["cap"]["Iq"] = .004*np.diff(r)
    return value


def assert_tree_equal(actual, expected, budget):
    """Compare every retained field, including arrays inside physical helpers."""
    if isinstance(expected, dict):
        assert set(expected) <= set(actual)
        for key in expected:
            assert_tree_equal(actual[key], expected[key], budget)
    elif isinstance(expected, np.ndarray):
        budget.check(expected.size)
        np.testing.assert_array_equal(actual, expected)
    elif isinstance(expected, (list, tuple)):
        assert len(actual) == len(expected)
        for a, e in zip(actual, expected):
            assert_tree_equal(a, e, budget)
    elif isinstance(expected, (str, float, int, bool, type(None), np.number)):
        assert actual == expected
    else:
        assert actual is expected


@pytest.mark.parametrize("active,alpha,n", [(False, 0., 4), (False, .03, 8),
                                          (True, 0., 4), (True, .03, 8)])
def test_first_order_callback_reproduces_every_frozen_stage_field(gradient_control_budget, active, alpha, n):
    budget = gradient_control_budget
    kw = stage_options(budget, n=n, active=active, alpha=alpha)
    state = initial_state(n, active)
    untouched = copy.deepcopy(state)
    old = background_stage(state, **kw)
    before = budget.evaluations
    new = coupled_stage(state, departure_gradient=first_order_departure_gradient, **kw)
    assert new["charged_evaluations"] == budget.evaluations-before
    assert_tree_equal(new, old, budget)
    assert_tree_equal(state, untouched, budget)


@pytest.mark.parametrize("alpha", [0., .03])
def test_quadratic_default_preserves_quiet_vacuum_and_charges_gradient(gradient_control_budget, alpha):
    budget = gradient_control_budget
    kw = stage_options(budget, n=8, alpha=alpha)
    state = initial_state(8)
    before = budget.evaluations
    row = coupled_stage(state, **kw)
    assert row["charged_evaluations"] == budget.evaluations-before
    for key in ("r", "L"):
        np.testing.assert_allclose(row["rates"]["departure"][key], 0., atol=4e-13)
    np.testing.assert_array_equal(row["rates"]["departure"]["c"], 0.)
    for key in ("Ip", "Iq"):
        np.testing.assert_array_equal(row["rates"]["cap"][key], 0.)
    assert max(abs(v) for v in row["rates"]["wall"].values()) < 4e-13
    assert row["cap_metric"]["budget"] is budget


def test_one_returned_gradient_controls_prefix_ale_and_wall_defect(gradient_control_budget):
    budget = gradient_control_budget
    kw = stage_options(budget, active=True, alpha=.01)
    seen = []
    def chosen(eta, values, *, budget):
        assert not eta.flags.writeable and not values.flags.writeable
        np.testing.assert_array_equal(values[:, 0], 0.)
        seen.append(values.copy())
        budget.check(values.size)
        return np.stack((.001*(1+eta), -.002*eta, .003*(1-eta)))
    row = coupled_stage(initial_state(active=True), departure_gradient=chosen, **kw)
    assert len(seen) == 1
    eta, W, Wx = kw["eta"], float(kw["gauge"]["W"]), float(kw["gauge"]["W_x"])
    ext = row["exterior_prefix"]
    for key, jet, delta, source_rate in (("r", "r_Ubar", "delta_r_x", 0.),
        ("L", "log_C_Ubar", "delta_log_C_x", 0.), ("c", "c_Ubar", "delta_c_x", .07)):
        G = row["departure_gradients"][key]
        base = 0. if key == "c" else ext["vacuum"][jet]
        np.testing.assert_allclose(row["exterior_first_jets"][jet], base+G/W, rtol=2e-15, atol=2e-15)
        expected = (ext[delta]-source_rate+eta*Wx/W*G)[1:-1]
        np.testing.assert_allclose(row["rates"]["departure"][key], expected, rtol=2e-15, atol=2e-15)
        assert row["wall_kinematic_defect"][key] == pytest.approx(
            row["wall_rates"][key]-row["exterior_longitudinal"][key][-1]
            -Wx*row["exterior_first_jets"][jet][-1], abs=2e-15)
        assert row["exterior"][key][0] == kw["source"]["log_C" if key == "L" else key]
        assert row["exterior"][key][-1] == row["physical_wall"][key]


def test_active_scalar_has_one_canonical_and_conformal_wall_flux(gradient_control_budget):
    budget = gradient_control_budget
    alpha, amplitude = .01, .025
    state = initial_state(active=True)
    row = coupled_stage(state, **stage_options(budget, active=True, alpha=alpha))
    wall, cap, hybrid = row["physical_wall"], row["cap_rhs"], row["hybrid"]
    R, T = wall["r"], wall["T"]
    qp, qm = max(T+np.arctanh(R), 0.), max(T-np.arctanh(R), 0.)
    chi = amplitude*(qp**3-qm**3)+np.cumsum(state["cap"]["Iq"])[-1]
    scale = np.sqrt(alpha/6)
    assert wall["c"] == pytest.approx(R*np.sinh(scale*chi/R)/scale, abs=2e-15)
    p0 = 3*amplitude*(qp*qp-qm*qm)/(1-R*R)
    q0 = 3*amplitude*(qp*qp+qm*qm)/(1-R*R)
    background_rate = row["wall_rates"]["T"]*(1-R*R)*p0+row["wall_rates"]["r"]*q0
    assert abs(background_rate) > 1e-4
    assert cap["wall_chi_x"] == pytest.approx(background_rate+sum(cap["rates"]["Iq"]), abs=3e-15)
    assert abs(row["canonical_wall_flux_defect"]) < 3e-14
    assert abs(row["conformal_wall_flux_defect"]) < 3e-14
    assert cap["incoming_departure"] == pytest.approx(R*hybrid["w_plus"]+chi/R-p0-q0, abs=3e-15)
    assert abs(cap["incoming_departure"]) > 1e-4


@pytest.mark.parametrize("kind", ["not_callable", "wrong_shape", "nonfinite", "mutation"])
def test_invalid_gradient_stops_before_dependent_prefix(gradient_control_budget, kind):
    kw = stage_options(gradient_control_budget)
    def invalid(eta, values, *, budget):
        budget.check(values.size)
        if kind == "mutation":
            values[0, 0] = 1.
        return np.zeros(eta.shape) if kind == "wrong_shape" else np.full(values.shape, np.nan)
    with pytest.raises(ValueError):
        coupled_stage(initial_state(), departure_gradient=None if kind == "not_callable" else invalid, **kw)


class ExhaustionBudget(Budget):
    """Deliberate local exhaustion still debits the one control allowance."""
    def __init__(self, parent, limit=1):
        super().__init__(maximum_evaluations=limit, cpu_seconds=parent.cpu_seconds, started=parent.started)
        self.parent = parent

    def check(self, count=0):
        self.parent.check(count)
        super().check(count)


def test_stage_budget_exhaustion_keeps_same_clock_and_counts(gradient_control_budget):
    kw = stage_options(gradient_control_budget)
    child = ExhaustionBudget(gradient_control_budget)
    kw["budget"] = child
    with pytest.raises(RuntimeError, match="budget exhausted"):
        coupled_stage(initial_state(), **kw)
    assert child.evaluations > child.maximum_evaluations
    assert child.started == gradient_control_budget.started


@pytest.mark.parametrize("started", [True, -1., float("nan")])
def test_standalone_stage_checks_budget_origin_even_with_custom_callback(gradient_control_budget, started):
    kw = stage_options(gradient_control_budget)
    child = ExhaustionBudget(gradient_control_budget)
    child.started = started
    kw["budget"] = child
    with pytest.raises(ValueError, match="budget"):
        coupled_stage(initial_state(), departure_gradient=first_order_departure_gradient, **kw)
    assert child.evaluations == 0
