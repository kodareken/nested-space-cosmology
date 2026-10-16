"""Independent fixed-Ubar transformation and exact-vacuum FD controls."""
import copy
import time

import numpy as np
import pytest
from scipy.special import lambertw

from recursive_horizons.nsc_atlas_characteristic_initialization import Budget
from recursive_horizons.nsc_atlas_coupled_diagnostics import (
    coupled_history_diagnostics, compare_rectangle_overlap,
)

TUBE = {"r": (.2, 1.), "L": (-2., 2.), "c": (-3., 3.)}


@pytest.fixture(scope="module")
def budget():
    result = Budget(maximum_evaluations=20000, cpu_seconds=30)
    yield result
    result.check()
    print("coupled diagnostic evaluations:", result.evaluations,
          "process CPU:", time.process_time()-result.started)


def snapshots(x, X, source):
    source_values = {"r": source["r"], "L": source["log_C"], "c": source["c"]}
    return [{"x": float(at), "state": {
        "increment": {key: X[key][i, 1:-1]-source_values[key][i] for key in X},
        "wall": {**{key: X[key][i, -1] for key in X}, "T": 0., "s": 0.},
        "cap": {"f": [0.]*3, "Pi": [0.]*3, "Q": [0.]*3}}}
        for i, at in enumerate(x)]


def manufactured(*, quadratic_gauge=False, nt=7, ne=9):
    x, eta = np.linspace(.2, .4, nt), np.linspace(0, 1, ne)
    W = .4+.2*x+(.07*x*x if quadratic_gauge else 0)
    Wx = .2+(.14*x if quadratic_gauge else 0)
    if np.ndim(Wx) == 0:
        Wx = np.full_like(x, Wx)
    gauge = {"W": W, "W_x": Wx, "W_xx": np.full_like(x, .14 if quadratic_gauge else 0.)}
    U, t = W[:, None]*eta, x[:, None]
    # For quadratic W use physical fields linear in U, so every numerical
    # partial needed for the coordinate cancellation is polynomial-exact.
    a2, cross = (0., 0.) if quadratic_gauge else (.03, .05)
    X = {"r": .6+.1*U+a2*U*U+.04*t+.02*t*t+cross*U*t,
         "L": .1+.2*U-.1*t, "c": .2+.1*U+.2*t}
    source = {"r": .6+.04*x+.02*x*x, "r_x": .04+.04*x,
              "r_xx": np.full_like(x, .04), "log_C": .1-.1*x,
              "log_C_x": np.full_like(x, -.1), "c": .2+.2*x,
              "c_x": np.full_like(x, .2)}
    exact = {"r_Ubar": .1+2*a2*U+cross*t,
             "r_Ubar_Ubar": np.full_like(U, 2*a2),
             "r_x": .04+.04*t+cross*U, "r_xx": np.full_like(U, .04),
             "L_Ubar": np.full_like(U, .2), "L_x": np.full_like(U, -.1),
             "c_Ubar": np.full_like(U, .1), "c_x": np.full_like(U, .2)}
    # Broadcast a constant-in-U exact derivative to the full rectangular shape.
    exact = {key: np.broadcast_to(value, U.shape) for key, value in exact.items()}
    return {"snapshots": snapshots(x, X, source), "eta": eta,
            "gauge": gauge, "source": source}, X, exact


def diagnose(data, budget, *, alpha=.3, fronts=()):
    return coupled_history_diagnostics(**data, alpha=alpha, tube=TUBE,
                                       budget=budget, source_fronts=fronts)


@pytest.mark.parametrize("quadratic_gauge", [False, True])
def test_polynomial_fixed_coordinate_jets_constraints_and_mass(quadratic_gauge, budget):
    data, X, exact = manufactured(quadratic_gauge=quadratic_gauge)
    row = diagnose(data, budget)
    take = np.s_[1:-1, 1:-1]
    for key, values in exact.items():
        np.testing.assert_allclose(row["derivatives"][key], values[take], atol=4e-12, rtol=4e-12)
    r, L, c = (X[key][take] for key in ("r", "L", "c"))
    u = c/r
    D = 1+.3*u*u/6
    for coordinate, second, result in (("Ubar", "r_Ubar_Ubar", "E_Ubar"), ("x", "r_xx", "E_x")):
        ra, La, ca = (exact[key+"_"+coordinate][take] for key in ("r", "L", "c"))
        ua = (ca-u*ra)/r
        expected = exact[second][take]-La*ra+.3*r*ua*ua/(2*D)
        np.testing.assert_allclose(row[result], expected, atol=4e-12, rtol=4e-12)
        assert np.max(abs(row[result])) > .005
    mass = r*(1+4*exact["r_Ubar"][take]*exact["r_x"][take]/np.exp(L))/2
    np.testing.assert_allclose(row["mass"], mass, atol=2e-14)
    # Omitting the W_xx term for the quadratic map gives an O(1) error,
    # even though every saved field is a simple analytic polynomial.
    if quadratic_gauge:
        omitted = row["interior_eta"][None, :]*data["gauge"]["W_xx"][1:-1, None]*.1
        assert np.max(abs(omitted)) > .01
    assert row["scope"]["mixed_equations_used"] is False
    assert row["scope"]["cap_constraint_evaluated"] is False


def schwarzschild(n):
    x, eta = np.linspace(.4, .47, n), np.linspace(0, 1, n)
    W, Wx = .22+.4*(x-.4)+.2*(x-.4)**2, .4+.4*(x-.4)
    U, t = W[:, None]*eta, x[:, None]
    product = (1-U/3.2)*(1-t/.8)*np.exp(t/.8)
    r = .8*(1+lambertw(-product/np.e).real)
    L = np.log(t/r)+(t-r)/.8
    # A constant normalized physical scalar has zero stress.  Its invariant
    # c remains nonzero and exercises the diagnostic's frame conversion.
    X = {"r": r, "L": L, "c": .3*r}
    source = {"r": x, "r_x": np.ones(n), "r_xx": np.zeros(n),
              "log_C": np.zeros(n), "log_C_x": np.zeros(n),
              "c": .3*x, "c_x": np.full(n, .3)}
    return {"snapshots": snapshots(x, X, source), "eta": eta,
            "gauge": {"W": W, "W_x": Wx, "W_xx": np.full(n, .4)}, "source": source}


def test_exact_schwarzschild_constraint_and_geometric_mass_differences_refine(budget):
    errors, common_errors = [], []
    for n in (17, 33):
        row = diagnose(schwarzschild(n), budget, alpha=0.)
        errors.append(np.array([np.max(abs(row["E_Ubar"])), np.max(abs(row["E_x"])),
                                np.max(abs(row["mass"]-.4))]))
        # The finer full interior includes new points closer to each excluded
        # edge.  Estimate the observed ratio on identical physical nodes.
        common = np.s_[:, :] if n == 17 else np.s_[1::2, 1::2]
        common_errors.append(np.array([np.max(abs(row["E_Ubar"][common])),
                                       np.max(abs(row["E_x"][common])),
                                       np.max(abs(row["mass"][common]-.4))]))
        np.testing.assert_array_equal(row["source_constraint"], np.zeros(n))
        assert np.all(row["derivatives"]["r_Ubar"] > 0)
        assert np.all(row["derivatives"]["r_x"] > 0)
    ratios = common_errors[0]/common_errors[1]
    print("Schwarzschild diagnostic errors", [value.tolist() for value in errors],
          "common-node ratios", ratios.tolist())
    assert np.all(ratios > 3.7)
    assert np.all(ratios < 4.3)


def test_analytic_source_curvature_is_retained_separately(budget):
    data, X, _ = manufactured(quadratic_gauge=True, nt=7, ne=7)
    x = np.array([row["x"] for row in data["snapshots"]])
    W, eta = data["gauge"]["W"], data["eta"]
    source = data["source"]
    source.update(r=.6+.005*np.sin(70*x), r_x=.35*np.cos(70*x),
                  r_xx=-24.5*np.sin(70*x), log_C=np.zeros(7), log_C_x=np.zeros(7),
                  c=np.zeros(7), c_x=np.zeros(7))
    X = {"r": source["r"][:, None]+.1*W[:, None]*eta,
         "L": np.zeros((7, 7)), "c": np.zeros((7, 7))}
    data["snapshots"] = snapshots(x, X, source)
    row = diagnose(data, budget, alpha=0.)
    np.testing.assert_allclose(row["E_x"], source["r_xx"][1:-1, None]+np.zeros((5, 5)), atol=1e-12)
    np.testing.assert_array_equal(row["source_constraint"], source["r_xx"])
    naive_source_second = np.diff(source["r"], n=2)/(x[1]-x[0])**2
    assert np.max(abs(naive_source_second-source["r_xx"][1:-1])) > 1
    assert row["scope"]["source_error_removed"] is False


def test_front_masks_exclude_edges_and_summary_partitions(budget):
    data, _, _ = manufactured(nt=7, ne=7)
    row = diagnose(data, budget, fronts=(.31, .6))
    near = row["masks"]["front_near"]
    np.testing.assert_array_equal(near[:, 0], [False, False, True, True, False])
    assert near.shape == (5, 5)
    assert row["summaries"]["E_x"]["front_near"]["count"] == 10
    assert row["summaries"]["E_x"]["away_from_front"]["count"] == 15
    assert row["summaries"]["E_x"]["all_interior"]["count"] == 25
    assert np.min(row["interior_x"]) > data["snapshots"][0]["x"]
    assert np.max(row["interior_eta"]) < 1
    assert row["tube_margins"]["r"]["lower_margin"] > 0


def test_counts_input_immutability_and_no_boundary_derivative_claim(budget):
    data, _, _ = manufactured(nt=6, ne=6)
    original = copy.deepcopy(data)
    before = budget.evaluations
    row = diagnose(data, budget)
    assert budget.evaluations-before == row["settings"]["charged_evaluations"] == 258
    for key in data["source"]:
        np.testing.assert_array_equal(data["source"][key], original["source"][key])
    assert row["settings"]["boundary_stencils"] == "excluded"
    assert row["summaries"]["E_x"]["front_near"]["max_abs"] is None
    row["states"]["r"][:] = 99
    assert data["snapshots"][0]["state"]["wall"]["r"] < 1


@pytest.mark.parametrize("kind", ["few_times", "few_eta", "nonuniform_x", "nonuniform_eta", "wrong_endpoint",
                                  "negative_alpha", "bad_source", "missing_gauge", "nonpositive_gauge",
                                  "bad_increment", "tube", "fronts", "boolean", "budget"])
def test_malformed_inputs_stop_before_derivative_charges(kind, budget):
    data, _, _ = manufactured(nt=5, ne=5)
    kwargs = dict(alpha=.3, tube=copy.deepcopy(TUBE), budget=budget, source_fronts=())
    if kind == "few_times": data["snapshots"] = data["snapshots"][:3]
    if kind == "few_eta": data["eta"] = np.linspace(0, 1, 3)
    if kind == "nonuniform_x": data["snapshots"][1]["x"] += .002
    if kind == "nonuniform_eta": data["eta"][1] += .002
    if kind == "wrong_endpoint": data["eta"] += .1
    if kind == "negative_alpha": kwargs["alpha"] = -.1
    if kind == "bad_source": data["source"]["r_xx"] = np.zeros(4)
    if kind == "missing_gauge": del data["gauge"]["W_xx"]
    if kind == "nonpositive_gauge": data["gauge"]["W_x"][2] = 0
    if kind == "bad_increment": data["snapshots"][1]["state"]["increment"]["r"][0] = np.nan
    if kind == "tube": kwargs["tube"]["r"] = (.8, 1.)
    if kind == "fronts": kwargs["source_fronts"] = [.35, .25]
    if kind == "boolean": data["source"]["c_x"] = [True, .2, .2, .2, .2]
    if kind == "budget": kwargs["budget"] = None
    before = budget.evaluations
    with pytest.raises(ValueError):
        coupled_history_diagnostics(**data, **kwargs)
    assert budget.evaluations == before


def test_budget_exhaustion_keeps_original_counts_and_clock():
    data, _, _ = manufactured(nt=5, ne=5)
    limited = Budget(maximum_evaluations=10, cpu_seconds=30, evaluations=2)
    started = limited.started
    with pytest.raises(RuntimeError, match="budget exhausted"):
        diagnose(data, limited)
    assert limited.evaluations == 77
    assert limited.started == started


def overlap_data(quadratic=True):
    x, eta = np.linspace(.1, .3, 5), np.linspace(0., 1., 5)
    W = .4+.1*x
    old_x, old_U = x[::2].copy(), np.linspace(0., .3, 4)
    q = 1. if quadratic else 0.
    def fields(U, t):
        return {"r": .6+.1*U+q*.2*U*U+.1*t,
                "L": .1-.3*U+q*.4*U*U-.2*t,
                "c": .3+.2*U+q*.6*U*U+.3*t}
    current = fields(W[:, None]*eta, x[:, None])
    old = fields(old_U[:, None]+np.zeros_like(old_x), old_x[None, :])
    return [current, x, eta, W, old_x, old_U, old]


@pytest.mark.parametrize("quadratic", [False, True])
def test_rectangle_overlap_polynomial_rules_and_orientation(quadratic, budget):
    data = overlap_data(quadratic)
    before = budget.evaluations
    row = compare_rectangle_overlap(*data, budget=budget)
    assert budget.evaluations-before == row["settings"]["charged_evaluations"] == 180
    np.testing.assert_array_equal(row["matched_time_indices"], [0, 2, 4])
    np.testing.assert_array_equal(row["time_coordinate_difference"], np.zeros(3))
    assert row["difference_linear"]["r"].shape == (4, 3)
    for key in ("r", "L", "c"):
        np.testing.assert_allclose(row["difference_quadratic"][key], 0., atol=3e-16)
        np.testing.assert_allclose(row["difference_quadratic"][key]-row["difference_linear"][key],
                                   row["interpolation_sensitivity"][key], atol=3e-16)
        if quadratic:
            assert row["summaries"][key]["linear"]["full"]["max_abs"] > 1e-5
        else:
            np.testing.assert_allclose(row["difference_linear"][key], 0., atol=3e-16)
        assert row["summaries"][key]["linear"]["full"]["count"] == 12
        assert row["summaries"][key]["linear"]["interior"]["count"] == 2
    assert row["scope"]["interpolation_error_enclosed"] is False
    assert row["scope"]["time_interpolations"] == 0


def test_rectangle_overlap_accepts_only_reported_coordinate_roundoff(budget):
    data = overlap_data(False)
    data[4] = np.nextafter(data[4], np.inf)
    data[5] = np.linspace(0, np.nextafter(data[3][0], np.inf), 4)
    # New old-grid coordinates need their own exact linear reference values.
    U, t = data[5][:, None], data[4][None, :]
    data[6] = {"r": .6+.1*U+.1*t, "L": .1-.3*U-.2*t, "c": .3+.2*U+.3*t}
    row = compare_rectangle_overlap(*data, budget=budget)
    assert 0 < row["max_abs_time_coordinate_difference"] < 1e-14
    assert 0 < row["max_eta_endpoint_roundoff_adjustment"] < 8*np.finfo(float).eps
    assert row["eta_used"][-1, 0] == 1
    for values in row["difference_linear"].values():
        np.testing.assert_allclose(values, 0., atol=3e-16)


@pytest.mark.parametrize("kind", ["unmatched_time", "out_of_domain", "bad_shape", "bad_W", "too_few_eta", "boolean"])
def test_rectangle_overlap_never_extrapolates_or_interpolates_time(kind, budget):
    data = overlap_data()
    if kind == "unmatched_time": data[4] += 1e-8
    if kind == "out_of_domain": data[5] = np.linspace(0., .5, 4)
    if kind == "bad_shape": data[6]["r"] = data[6]["r"].T
    if kind == "bad_W": data[3][1] = 0.
    if kind == "too_few_eta": data[2] = [0., 1.]
    if kind == "boolean": data[3] = [True, .4, .4, .4, .4]
    before = budget.evaluations
    with pytest.raises(ValueError):
        compare_rectangle_overlap(*data, budget=budget)
    assert budget.evaluations == before
