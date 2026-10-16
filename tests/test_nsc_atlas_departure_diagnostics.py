"""Analytic vacuum and independent non-solution controls for saved departures."""
import copy
import time

import numpy as np
import pytest
from scipy.special import lambertw

from recursive_horizons.nsc_atlas_characteristic_initialization import Budget
from recursive_horizons.nsc_atlas_departure_diagnostics import departure_history_diagnostics

TUBE = {"r": (.2, .8), "L": (-2., 2.), "c": (-3., 3.)}


@pytest.fixture(scope="module")
def budget():
    shared = Budget(maximum_evaluations=20000, cpu_seconds=30.)
    yield shared
    shared.check()
    print("departure diagnostic controls:", shared.evaluations, "points;",
          time.process_time()-shared.started, "process CPU seconds")


def analytic_vacuum(U, x):
    """Independent closed Schwarzschild formula; no diagnostic/solver call."""
    r = .8*(1+lambertw(-(1-U/3.2)*(1-x/.8)*np.exp(x/.8)/np.e).real)
    r = np.where(U == 0., x, r)
    L = np.log(x/r)+(x-r)/.8
    C = np.exp(L)
    ru, rx = C*(.8-x)/(4*x), (1-U/3.2)*C
    lu, lx = -(1/r+1/.8)*ru, 1/x+1/.8-(1/r+1/.8)*rx
    return {"r": r, "L": L, "r_Ubar": ru, "r_x": rx,
            "r_Ubar_Ubar": lu*ru, "r_xx": lx*rx,
            "L_Ubar": lu, "L_x": lx}


def history(*, nt=7, ne=9, nonzero=False, quadratic=False):
    x, eta = np.linspace(.4, .47, nt), np.linspace(0., 1., ne)
    W = .2+.3*x+(.04*x*x if quadratic else 0.)
    Wx = .3+(.08*x if quadratic else np.zeros(nt))
    gauge = {"W": W, "W_x": Wx, "W_xx": np.full(nt, .08 if quadratic else 0.)}
    U, t = W[:, None]*eta, x[:, None]
    vac = analytic_vacuum(U, t)
    source = {"r": x.copy(), "r_x": np.ones(nt), "r_xx": np.zeros(nt),
              "log_C": np.zeros(nt), "log_C_x": np.zeros(nt),
              "c": np.zeros(nt), "c_x": np.zeros(nt)}
    Z = {key: np.zeros_like(U) for key in ("r", "L", "c")}
    expected = {key: vac[key].copy() for key in vac if key not in ("r", "L")}
    expected.update(c_Ubar=np.zeros_like(U), c_x=np.zeros_like(U))
    if nonzero:
        # A steep known source offset must retain its analytic curvature,
        # rather than the curvature of these few sampled source values.
        offset = .001*np.sin(70*x)
        source.update(r=x+offset, r_x=1+.07*np.cos(70*x), r_xx=-4.9*np.sin(70*x),
                      log_C=.02*x*x, log_C_x=.04*x, c=.1+.03*x*x, c_x=.06*x)
        b, cross = (0., 0.) if quadratic else (.03, .04)
        Z = {"r": .08*U+b*U*U+cross*U*t, "L": .07*U, "c": .12*U}
        expected["r_Ubar"] += .08+2*b*U+cross*t
        expected["r_Ubar_Ubar"] += 2*b
        expected["r_x"] += .07*np.cos(70*t)+cross*U
        expected["r_xx"] += -4.9*np.sin(70*t)
        expected["L_Ubar"] += .07
        expected["L_x"] += .04*t
        expected["c_Ubar"] += .12
        expected["c_x"] += .06*t
    X = {"r": vac["r"]+(source["r"]-x)[:, None]+Z["r"],
         "L": vac["L"]+source["log_C"][:, None]+Z["L"],
         "c": source["c"][:, None]+Z["c"]}
    snapshots = [{"x": float(at), "state": {
        "departure": {key: Z[key][i, 1:-1].copy() for key in Z}},
        "physical": {"wall": {key: float(X[key][i, -1]) for key in X}}}
        for i, at in enumerate(x)]
    return {"snapshots": snapshots, "eta": eta, "gauge": gauge, "source": source}, X, expected


def diagnose(data, budget, *, alpha=.2, fronts=(), tube=None):
    return departure_history_diagnostics(**data, alpha=alpha, tube=TUBE if tube is None else tube,
                                          source_fronts=fronts, budget=budget)


@pytest.mark.parametrize("alpha", [0., .2])
def test_exact_vacuum_constraints_do_not_require_grid_refinement(alpha, budget):
    data, _, _ = history(quadratic=True)
    row = diagnose(data, budget, alpha=alpha)
    np.testing.assert_array_equal(row["E_Ubar"], 0.)
    np.testing.assert_array_equal(row["E_x"], 0.)
    np.testing.assert_allclose(row["mass"], .4, atol=3e-16, rtol=0.)
    np.testing.assert_array_equal(row["source_constraint"], 0.)
    assert row["scope"]["mixed_equations_used"] is False
    assert row["scope"]["vacuum_background_differenced"] is False


@pytest.mark.parametrize("quadratic", [False, True])
def test_manufactured_departures_retain_nonzero_constraint_violations(quadratic, budget):
    data, X, exact = history(nonzero=True, quadratic=quadratic)
    row = diagnose(data, budget)
    take = np.s_[1:-1, 1:-1]
    for key, values in exact.items():
        np.testing.assert_allclose(row["derivatives"][key], values[take], atol=3e-11, rtol=1e-12)
    r, L, c = (X[key][take] for key in ("r", "L", "c"))
    u = c/r
    for name, second, result in (("Ubar", "r_Ubar_Ubar", "E_Ubar"), ("x", "r_xx", "E_x")):
        ra, la, ca = (exact[key+"_"+name][take] for key in ("r", "L", "c"))
        kinetic = .2*r*((ca-u*ra)/r)**2/(2*(1+.2*u*u/6))
        target = exact[second][take]-la*ra+kinetic
        np.testing.assert_allclose(row[result], target, atol=3e-11, rtol=1e-12)
        assert np.max(abs(row[result])) > .001
    mass = r*(1+4*exact["r_Ubar"][take]*exact["r_x"][take]/np.exp(L))/2
    np.testing.assert_allclose(row["mass"], mass, atol=2e-13, rtol=0.)
    source = data["source"]
    us = source["c"]/source["r"]
    ux = (source["c_x"]-us*source["r_x"])/source["r"]
    target = source["r_xx"]-source["log_C_x"]*source["r_x"]+.1*source["r"]*ux*ux/(1+.2*us*us/6)
    np.testing.assert_allclose(row["source_constraint"], target, atol=1e-15)
    x = np.array([s["x"] for s in data["snapshots"]])
    naive = np.diff(source["r"], n=2)/(x[1]-x[0])**2
    assert np.max(abs(naive-source["r_xx"][1:-1])) > .1
    assert row["scope"]["source_error_removed"] is False


def test_physical_wall_supplies_terminal_departure_and_front_partition(budget):
    data, X, _ = history(nonzero=True)
    before = copy.deepcopy(data)
    count = budget.evaluations
    row = diagnose(data, budget, fronts=(.435,))
    expected = 4*7*9+9*5*7+7
    assert budget.evaluations-count == row["settings"]["charged_evaluations"] == expected
    for key in X:
        np.testing.assert_allclose(row["states"][key], X[key], atol=2e-16)
        np.testing.assert_array_equal(row["states"][key][:, -1], [s["physical"]["wall"][key] for s in data["snapshots"]])
        np.testing.assert_array_equal(data["snapshots"][0]["state"]["departure"][key], before["snapshots"][0]["state"]["departure"][key])
    assert row["E_x"].shape == (5, 7)
    near = row["masks"]["front_near"]
    x = row["x"]
    np.testing.assert_array_equal(near[:, 0], (x[:-2] <= .435) & (.435 <= x[2:]))
    assert row["summaries"]["E_x"]["all_interior"]["count"] == 35
    assert row["tube_margins"]["r"]["lower_margin"] > 0
    row["states"]["r"][:] = 100
    assert data["snapshots"][-1]["physical"]["wall"]["r"] < 1


@pytest.mark.parametrize("change", ["many_times", "few_eta", "bad_time", "bad_eta", "alpha", "bool",
                                    "missing_wall", "bad_departure", "gauge", "source", "fronts", "tube"])
def test_malformed_inputs_and_domain_exits_are_rejected(change, budget):
    data, _, _ = history(nt=5, ne=5)
    kwargs = {"alpha": 0., "fronts": (), "tube": copy.deepcopy(TUBE)}
    if change == "many_times": data["snapshots"] = data["snapshots"]*26
    if change == "few_eta": data["eta"] = [0., .5, 1.]
    if change == "bad_time": data["snapshots"][1]["x"] += .001
    if change == "bad_eta": data["eta"][1] += .001
    if change == "alpha": kwargs["alpha"] = -1.
    if change == "bool": data["source"]["c_x"] = [True, 0., 0., 0., 0.]
    if change == "missing_wall": del data["snapshots"][1]["physical"]["wall"]
    if change == "bad_departure": data["snapshots"][1]["state"]["departure"]["r"][0] = np.nan
    if change == "gauge": data["gauge"]["W"][:] = 4.
    if change == "source": data["source"]["r_xx"] = [0.]*4
    if change == "fronts": kwargs["fronts"] = [.45, .42]
    if change == "tube": kwargs["tube"]["r"] = (.6, .8)
    with pytest.raises(ValueError):
        diagnose(data, budget, **kwargs)


def test_budget_exhaustion_stops_algebra_before_derivatives():
    data, _, _ = history(nt=5, ne=5)
    budget = Budget(maximum_evaluations=4, cpu_seconds=30.)
    with pytest.raises(RuntimeError, match="budget exhausted"):
        diagnose(data, budget)
    assert budget.evaluations == 5
