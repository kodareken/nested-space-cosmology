"""Local algebra and analytic vacuum controls, with no prepared-source march."""
import copy
import time

import numpy as np
import pytest
from scipy.special import lambertw

from recursive_horizons.nsc_atlas_characteristic_cell import (
    CellSolveError, cell_evaluation, solve_characteristic_cell,
)
from recursive_horizons.nsc_atlas_characteristic_initialization import Budget


@pytest.fixture(autouse=True)
def bounded_test_cpu():
    start = time.process_time()
    yield
    assert time.process_time()-start < 30


TUBE = dict(r=(.1, 1.5), L=(-3., 3.), c=(-3., 3.))


def args():
    return (dict(r=.5, L=.1, c=.2), dict(r=.51, L=.13, c=.24),
            dict(r=.515, L=.08, c=.21))


def vector(row):
    return np.array([row[key] for key in ("r", "L", "c")])


def corners(function, u, v, du, dv):
    return function(u, v), function(u, v+dv), function(u+du, v), function(u+du, v+dv)


def minkowski(u, v):
    r = .6+(v-u)/2
    return dict(r=r, L=0., c=.37*r)


def trapped_affine(u, x):
    # Ubar=0 is the original source. This is the same shifted affine gauge,
    # derived by transforming Kruskal U,V, not a call to the cell equations.
    product = (1-u/3.2)*(1-x/.8)*np.exp(x/.8)
    r = float(.8*(1+lambertw(-product/np.e).real))
    C = x/r*np.exp((x-r)/.8)
    return dict(r=r, L=np.log(C), c=.37*r)


def de_sitter(u, v):
    r = np.tanh((v-u)/2)
    return dict(r=r, L=np.log(1-r*r), c=.2+.17*u-.11*v)


def test_exact_minkowski_cell_with_finite_alpha_constant_canonical_scalar():
    du, dv = .02, .03
    sw, nw, se, exact = corners(minkowski, 0., 0., du, dv)
    budget = Budget(cpu_seconds=30.)
    row = solve_characteristic_cell(sw, nw, se, du=du, dv=dv, alpha=.7,
                                    region="exterior", tube=TUBE, budget=budget)
    np.testing.assert_allclose(vector(row["ne"]), vector(exact), rtol=0, atol=3e-16)
    assert row["updates"] == row["backtrack_trials"] == 0
    assert row["seed_kind"] == "affine_predictor"
    assert budget.evaluations == row["rhs_evaluations"]+row["jacobian_evaluations"] == 2
    assert row["residual_norm"] <= 1 and row["correction_norm"] <= 1
    assert row["first_jets"]["r_U"] == pytest.approx(-.5)
    assert row["first_jets"]["r_V"] == pytest.approx(.5)


@pytest.mark.parametrize("analytic,origin,widths,alpha,region", [
    (trapped_affine, (0., .4), (.08, .02), .7, "exterior"),
    (de_sitter, (0., 1.), (.08, .06), 0., "cap"),
])
def test_curved_vacuum_one_cell_has_fourth_order_local_error(analytic, origin, widths, alpha, region):
    errors, defects = [], []
    budget = Budget(cpu_seconds=30.)
    for factor in (1., .5, .25):
        du, dv = (factor*h for h in widths)
        sw, nw, se, exact = corners(analytic, *origin, du, dv)
        defect = cell_evaluation(exact, sw, nw, se, du=du, dv=dv, alpha=alpha,
                                 region=region, tube=TUBE, budget=budget)
        row = solve_characteristic_cell(sw, nw, se, du=du, dv=dv, alpha=alpha,
            region=region, tube=TUBE, budget=budget, atol=(2e-15,)*3, rtol=0.)
        errors.append(np.max(abs(vector(row["ne"])-vector(exact))))
        defects.append(np.max(abs(defect["residual"])))
        assert row["updates"] > 0
        assert row["residual_norm"] <= 1 and row["correction_norm"] <= 1
    # Exact corners on three sides leave local O(du*dv*(du²+dv²)) error.
    # This is not a fourth-order assertion for a propagated grid.
    for values in (errors, defects):
        assert values[0] > 1e-10
        assert all(13 < values[i]/values[i+1] < 19 for i in range(2)), values


def test_affine_vacuum_source_and_orientation_are_the_original_ones():
    for x in (.4, .45, .5):
        row = trapped_affine(0., x)
        assert row["r"] == pytest.approx(x, abs=3e-16)
        assert row["L"] == pytest.approx(0., abs=2e-15)
    sw, nw, se, ne = corners(trapped_affine, 0., .4, .04, .01)
    row = cell_evaluation(ne, sw, nw, se, du=.04, dv=.01, alpha=.7,
                          region="exterior", tube=TUBE)
    assert row["first_jets"]["r_U"] > 0 and row["first_jets"]["r_V"] > 0


@pytest.mark.parametrize("region,alpha", [("cap", .6), ("exterior", .4), ("cap", 0.)])
def test_analytic_newton_jacobian_against_independent_algebra_perturbation(region, alpha):
    sw, nw, se = args()
    ne = dict(r=.529, L=.108, c=.267)
    kwargs = dict(du=.03, dv=.02, alpha=alpha, region=region, tube=TUBE)
    row = cell_evaluation(ne, sw, nw, se, jacobian=True, **kwargs)
    # Differentiate the algebraic residual only, never sampled physical fields.
    for direction in (np.array([.3, -.4, .7]), np.array([-.5, .2, .1])):
        step = 2e-6
        plus = dict(zip(("r", "L", "c"), vector(ne)+step*direction))
        minus = dict(zip(("r", "L", "c"), vector(ne)-step*direction))
        rp = cell_evaluation(plus, sw, nw, se, **kwargs)["residual"]
        rm = cell_evaluation(minus, sw, nw, se, **kwargs)["residual"]
        np.testing.assert_allclose(row["jacobian"]@direction, (rp-rm)/(2*step), rtol=2e-9, atol=2e-10)


def test_arbitrary_nonlinear_algebra_root_and_fixed_scales():
    known = args()
    budget = Budget(cpu_seconds=30.)
    row = solve_characteristic_cell(*known, du=.03, dv=.02, alpha=.6, region="cap",
        tube=TUBE, budget=budget, atol=(1e-13, 2e-13, 3e-13), rtol=1e-12,
        reference_scales=(.4, 1., 1.), initial_guess=dict(r=.9, L=.7, c=-.2))
    independent = cell_evaluation(row["ne"], *known, du=.03, dv=.02, alpha=.6,
                                  region="cap", tube=TUBE)
    np.testing.assert_array_equal(row["residual"], independent["residual"])
    expected = np.array([1e-13, 2e-13, 3e-13])+1e-12*np.array([.515, 1., 1.])
    np.testing.assert_array_equal(row["scales"], expected)
    assert row["seed_kind"] == "supplied" and 0 < row["updates"] <= 12
    assert row["backtrack_trials"] <= 8*row["updates"]
    assert row["residual_norm"] <= 1 and row["correction_norm"] <= 1
    assert budget.evaluations == row["rhs_evaluations"]+row["jacobian_evaluations"]


def test_corner_values_are_not_mutated_and_midpoint_orientation_is_explicit():
    known = args()
    before = copy.deepcopy(known)
    ne = dict(r=.529, L=.108, c=.267)
    row = cell_evaluation(ne, *known, du=.03, dv=.02, alpha=.6, region="cap", tube=TUBE)
    sw, nw, se = map(vector, known)
    np.testing.assert_allclose(vector(row["midpoint"]), (vector(ne)+sw+nw+se)/4)
    assert row["first_jets"]["c_U"] == pytest.approx((ne["c"]+se[2]-nw[2]-sw[2])/.06)
    assert row["first_jets"]["c_V"] == pytest.approx((ne["c"]+nw[2]-se[2]-sw[2])/.04)
    assert known == before


def test_update_limit_is_an_explicit_failure_without_hidden_retry():
    budget = Budget(cpu_seconds=30.)
    with pytest.raises(CellSolveError) as caught:
        solve_characteristic_cell(*args(), du=.03, dv=.02, alpha=.6, region="cap", tube=TUBE,
            budget=budget, max_updates=1, initial_guess=dict(r=.9, L=.7, c=-.2))
    failure = caught.value
    assert failure.reason == "newton_update_limit"
    assert failure.diagnostics["updates"] == 1
    assert failure.diagnostics["backtrack_trials"] <= 8
    assert budget.evaluations == failure.diagnostics["rhs_evaluations"]+failure.diagnostics["jacobian_evaluations"]


def test_budget_is_shared_counted_and_not_reset():
    budget = Budget(maximum_evaluations=2, cpu_seconds=30.)
    row = cell_evaluation(dict(r=.529, L=.108, c=.267), *args(), du=.03, dv=.02,
        alpha=.6, region="cap", tube=TUBE, budget=budget, jacobian=True)
    assert budget.evaluations == row["rhs_evaluations"]+row["jacobian_evaluations"] == 2
    with pytest.raises(CellSolveError, match="budget_exhausted"):
        solve_characteristic_cell(*args(), du=.03, dv=.02, alpha=.6, region="cap", tube=TUBE, budget=budget)
    assert budget.evaluations == 3


def test_ill_conditioned_limit_is_reported():
    with pytest.raises(CellSolveError, match="ill_conditioned_jacobian"):
        solve_characteristic_cell(*args(), du=.03, dv=.02, alpha=.6, region="cap", tube=TUBE,
                                  maximum_condition=1.)


def test_tube_rejection_has_finite_backtracking_and_a_recorded_seed():
    known = (dict(r=.5, L=0., c=0.), dict(r=.4, L=0., c=0.), dict(r=.4, L=0., c=0.))
    tube = dict(r=(.35, .6), L=(-1., 1.), c=(-1., 1.))
    with pytest.raises(CellSolveError) as caught:
        solve_characteristic_cell(*known, du=.01, dv=.01, alpha=0., region="exterior",
                                  tube=tube, max_backtracks=1)
    failure = caught.value
    assert failure.reason == "backtracking_limit"
    assert failure.diagnostics["seed_kind"] == "convex_corner_seed"
    assert failure.diagnostics["updates"] == 0
    assert failure.diagnostics["backtrack_trials"] == 1
    assert failure.diagnostics["rhs_evaluations"] == failure.diagnostics["jacobian_evaluations"] == 1
    assert failure.diagnostics["ne"]["r"] == .4


def test_expired_shared_budget_is_a_diagnostic_failure_before_evaluation():
    budget = Budget(maximum_evaluations=0, evaluations=1, cpu_seconds=30.)
    with pytest.raises(CellSolveError) as caught:
        solve_characteristic_cell(*args(), du=.03, dv=.02, alpha=.6,
                                  region="cap", tube=TUBE, budget=budget)
    assert caught.value.reason == "budget_exhausted"
    assert caught.value.diagnostics["rhs_evaluations"] == 0


def test_nonrepresentable_first_jets_are_rejected_before_bulk_evaluation():
    budget = Budget(cpu_seconds=30.)
    with pytest.raises(ValueError, match="first jets"):
        cell_evaluation(dict(r=.529, L=.108, c=.267), *args(), du=1e-320, dv=.02,
                        alpha=.6, region="cap", tube=TUBE, budget=budget)
    assert budget.evaluations == 0


@pytest.mark.parametrize("setting,value", [("du", 0.), ("dv", -1.), ("alpha", -.1),
    ("du", True), ("dv", np.nan), ("alpha", np.inf), ("region", "vacuum"),
    ("max_updates", 13), ("max_backtracks", 9), ("max_updates", True),
    ("rtol", -1.), ("atol", (0., 1e-12, 1e-12)), ("reference_scales", (1., 0., 1.))])
def test_invalid_parameters_do_not_start_iteration(setting, value):
    options = dict(du=.03, dv=.02, alpha=.6, region="cap", tube=TUBE)
    options[setting] = value
    with pytest.raises(ValueError):
        solve_characteristic_cell(*args(), **options)


@pytest.mark.parametrize("tube", [dict(r=(0., 1.), L=(-1., 1.), c=(-1., 1.)),
    dict(r=(.1, 1.), L=(-1000., 1.), c=(-1., 1.)),
    dict(r=(.1, 1.), L=(-1., 1000.), c=(-1., 1.)),
    dict(r=(.6, 1.), L=(-1., 1.), c=(-1., 1.)),
    dict(r=(1., .1), L=(-1., 1.), c=(-1., 1.)), dict(r=(.1, 1.))])
def test_invalid_or_incompatible_tube(tube):
    with pytest.raises(ValueError):
        solve_characteristic_cell(*args(), du=.03, dv=.02, alpha=.6, region="cap", tube=tube)


def test_midpoint_and_ne_tube_checks_and_no_second_derivative_claim():
    with pytest.raises(ValueError, match="NE"):
        cell_evaluation(dict(r=2., L=.1, c=.2), *args(), du=.03, dv=.02,
                        alpha=.6, region="cap", tube=TUBE)
    row = cell_evaluation(dict(r=.529, L=.108, c=.267), *args(), du=.03, dv=.02,
                          alpha=.6, region="cap", tube=TUBE)
    assert not any("constraint" in key for key in row)
