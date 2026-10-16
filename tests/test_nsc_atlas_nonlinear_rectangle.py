"""Manufactured grid controls only; no original prepared-source construction."""
import copy
import time

import numpy as np
import pytest
from scipy.special import lambertw

from recursive_horizons.nsc_atlas_characteristic_initialization import Budget
from recursive_horizons.nsc_atlas_nonlinear_rectangle import (
    RectangleSolveError, compare_common_grids, grid_diagnostics, solve_rectangle,
)

TUBE = dict(r=(.1, 1.5), L=(-3., 3.), c=(-3., 3.))


@pytest.fixture(scope="module")
def aggregate_budget():
    budget = Budget(cpu_seconds=30., maximum_evaluations=20000)
    yield budget
    budget.check()
    assert budget.evaluations <= 20000


@pytest.fixture(autouse=True)
def bounded_test_cpu(aggregate_budget):
    start = time.process_time()
    yield
    assert time.process_time()-start < 30
    aggregate_budget.check()


def axes(cells=8):
    return np.linspace(0., .16, cells+1), np.linspace(.4, .48, cells+1)


def boundaries(states):
    return ({key: value[0, :].copy() for key, value in states.items()},
            {key: value[:, 0].copy() for key, value in states.items()})


def vacuum(Ubar, x):
    U, X = Ubar[:, None], x[None, :]
    q = 1-U/3.2
    product = q*(1-X/.8)*np.exp(X/.8)
    r = .8*(1+lambertw(-product/np.e).real)
    L = np.log(X/r)+(X-r)/.8
    C = np.exp(L)
    rU, rx = C*(.8-X)/(4*X), q*C
    return dict(r=r, L=L, c=np.zeros_like(r)), dict(
        r_Ubar=rU, r_x=rx, log_C_Ubar=-(1/r+1/.8)*rU,
        log_C_x=1/X+1/.8-(1/r+1/.8)*rx,
        c_Ubar=np.zeros_like(r), c_x=np.zeros_like(r))


def polynomial(Ubar, x):
    U, X = Ubar[:, None], x[None, :]
    r = .6+.2*U+.3*X
    L = .1*U+.2*X
    p = .4+.7*U-.3*X
    return dict(r=r, L=L, c=r*p), p


def test_nonzero_manufactured_constraints_are_detected_independently():
    U, x = axes()
    states, p = polynomial(U, x)
    alpha = .3
    result = grid_diagnostics(U, x, states, alpha=alpha)
    D = 1+alpha*p*p/6
    expected_U = -.1*.2+alpha*states["r"]*.7**2/(2*D)
    expected_x = -.2*.3+alpha*states["r"]*.3**2/(2*D)
    np.testing.assert_allclose(result["E_U"], expected_U, atol=2e-11, rtol=0)
    np.testing.assert_allclose(result["E_x"], expected_x, atol=2e-11, rtol=0)
    expected_mass = states["r"]*(1+4*.2*.3*np.exp(-states["L"]))/2
    np.testing.assert_allclose(result["m"], expected_mass, atol=2e-13, rtol=0)
    assert result["summary"]["interior"]["E_x"]["max_abs"] > .04


def test_pure_second_derivatives_and_one_sided_boundaries_are_not_omitted():
    Ubar, x = axes()
    U, X = Ubar[:, None], x[None, :]
    r = .6+.2*U+.3*X+.4*U*U-.1*X*X
    L = .1*U+.2*X
    row = grid_diagnostics(Ubar, x, dict(r=r, L=L, c=.37*r), alpha=.7)
    np.testing.assert_allclose(row["E_U"], .8-.1*(.2+.8*U)+np.zeros_like(X), atol=2e-11, rtol=0)
    np.testing.assert_allclose(row["E_x"], -.2-.2*(.3-.2*X)+np.zeros_like(U), atol=2e-11, rtol=0)
    np.testing.assert_allclose(row["derivatives"]["r_Ubar_Ubar"], .8, atol=2e-11, rtol=0)
    np.testing.assert_allclose(row["derivatives"]["r_x_x"], -.2, atol=2e-11, rtol=0)


def test_exact_schwarzschild_nodal_diagnostics_refine_independently_of_any_march():
    diagnostics = []
    for cells in (16, 32, 64):
        U, x = axes(cells)
        states, _ = vacuum(U, x)
        diagnostics.append(grid_diagnostics(U, x, states, alpha=0., source_fronts=(.421, .457)))
    for coarse, fine in zip(diagnostics[:-1], diagnostics[1:]):
        common = compare_common_grids(coarse, fine)
        for key in ("E_U", "E_x"):
            coarse_error = np.max(abs(common["coarse"][key]))
            fine_error = np.max(abs(common["fine"][key]))
            assert 3.5 < coarse_error/fine_error < 4.5
        c_mass = np.max(abs(common["coarse"]["m"]-.4))
        f_mass = np.max(abs(common["fine"]["m"]-.4))
        assert 3.5 < c_mass/f_mass < 4.5


def test_shared_grid_masks_and_paired_baselines_do_not_hide_raw_constraints():
    coarse, fine, coarse_zero, fine_zero = [] , [], [], []
    results = []
    for cells in (8, 16):
        U, x = axes(cells)
        states, _ = polynomial(U, x)
        results.append([grid_diagnostics(U, x, states, alpha=a, source_fronts=(.421, .457))
                        for a in (.3, 0.)])
    coarse, coarse_zero = results[0]
    fine, fine_zero = results[1]
    result = compare_common_grids(coarse, fine, coarse_zero=coarse_zero, fine_zero=fine_zero)
    expected_near = coarse["masks"]["front_near"][1:-1, 1:-1] | fine["masks"]["front_near"][2:-2:2, 2:-2:2]
    np.testing.assert_array_equal(result["masks"]["front_near"], expected_near)
    assert result["coarse_summary"]["interior"]["E_x"]["max_abs"] > .04
    assert np.max(abs(result["paired"]["differences"]["E_x"])) < 3e-11
    np.testing.assert_allclose(result["paired"]["coarse"]["E_U"],
        coarse["E_U"][1:-1, 1:-1]-coarse_zero["E_U"][1:-1, 1:-1])
    assert result["masks"]["interior"].shape == (7, 7)
    assert np.any(expected_near) and np.any(result["masks"]["away"])


def test_boundary_first_jet_defects_keep_signed_postpulse_scalar_memory():
    U, x = axes()
    states, _ = polynomial(U, x)
    # This independent polynomial has a nonzero transverse scalar gradient
    # beyond the declared front labels, testing bookkeeping without a pulse.
    source = dict(r_Ubar=np.full(len(x), .2), r_x=np.full(len(x), .3),
                  log_C_Ubar=np.full(len(x), .1), log_C_x=np.full(len(x), .2),
                  c_Ubar=.2*(.4-.3*x)+.7*(.6+.3*x),
                  c_x=.3*(.4-.3*x)-.3*(.6+.3*x))
    p = .4+.7*U-.3*x[0]
    r = .6+.2*U+.3*x[0]
    bottom = dict(r_Ubar=np.full(len(U), .2), r_x=np.full(len(U), .3),
                  log_C_Ubar=np.full(len(U), .1), log_C_x=np.full(len(U), .2),
                  c_Ubar=.2*p+.7*r, c_x=.3*p-.3*r)
    exact = grid_diagnostics(U, x, states, alpha=.3, source_fronts=(.421, .457),
                             source_jets=source, bottom_jets=bottom)
    for boundary in ("source", "bottom"):
        for values in exact["boundary_defects"][boundary].values():
            np.testing.assert_allclose(values, 0., atol=2e-13, rtol=0)
    assert np.all(source["c_Ubar"][exact["source_masks"]["postpulse"]] > .4)
    source["c_Ubar"] = source["c_Ubar"]+.012
    shifted = grid_diagnostics(U, x, states, alpha=.3, source_fronts=(.421, .457), source_jets=source)
    np.testing.assert_allclose(shifted["boundary_defects"]["source"]["c_Ubar"], -.012, atol=2e-13)
    assert shifted["boundary_summary"]["source"]["postpulse"]["c_Ubar"]["max_abs"] == pytest.approx(.012)


def test_exact_finite_alpha_minkowski_march_and_immutable_boundaries(aggregate_budget):
    U, x = np.linspace(0., .1, 5), np.linspace(0., .12, 5)
    r = .6+(x[None, :]-U[:, None])/2
    exact = dict(r=r, L=np.zeros_like(r), c=.37*r)
    source, bottom = boundaries(exact)
    before = copy.deepcopy((source, bottom))
    for values in (source, bottom):
        for value in values.values():
            value.flags.writeable = False
    row = solve_rectangle(U, x, source, bottom, alpha=.7, tube=TUBE, budget=aggregate_budget)
    for name in exact:
        np.testing.assert_allclose(row["states"][name], exact[name], atol=3e-15, rtol=0)
        np.testing.assert_array_equal(source[name], before[0][name])
        np.testing.assert_array_equal(bottom[name], before[1][name])
    assert row["summary"]["cell_count"] == 16
    assert row["summary"]["budget_end"]-row["summary"]["budget_start"] == (
        row["summary"]["rhs_evaluations"]+row["summary"]["jacobian_evaluations"])
    assert row["summary"]["max_residual_norm"] <= 1
    assert row["summary"]["max_correction_norm"] <= 1
    assert sum(row["summary"]["seed_counts"].values()) == 16
    assert row["settings"]["axis_order"] == ["Ubar", "x"]


def test_marched_zero_field_schwarzschild_geometry_has_second_order_global_refinement(aggregate_budget):
    errors = []
    for cells in (2, 4, 8):
        U, x = axes(cells)
        exact, _ = vacuum(U, x)
        source, bottom = boundaries(exact)
        row = solve_rectangle(U, x, source, bottom, alpha=0., tube=TUBE, budget=aggregate_budget)
        errors.append(max(np.max(abs(row["states"][key]-exact[key])) for key in ("r", "L")))
        np.testing.assert_array_equal(row["states"]["c"], np.zeros_like(exact["c"]))
        assert row["summary"]["max_residual_norm"] <= 1
    assert errors[0] > 1e-7
    assert all(3.5 < errors[i]/errors[i+1] < 4.5 for i in range(2)), errors


def test_cell_failure_retains_coordinates_corners_and_cumulative_counts():
    U, x = axes(2)
    r = .6+(x[None, :]-U[:, None])/2
    source, bottom = boundaries(dict(r=r, L=np.zeros_like(r), c=.37*r))
    budget = Budget(maximum_evaluations=2, cpu_seconds=30.)
    with pytest.raises(RectangleSolveError) as caught:
        solve_rectangle(U, x, source, bottom, alpha=.7, tube=TUBE, budget=budget)
    data = caught.value.diagnostics
    assert data["reason"] == "budget_exhausted" and data["indices"] == [1, 2]
    assert data["completed_cells"] == 1 and data["budget_end"] == 3
    assert data["coordinates"] == dict(Ubar=U[:2].tolist(), x=x[1:].tolist())
    assert set(data["input_corners"]) == {"sw", "nw", "se"}


def test_newton_failure_is_not_retried_or_subdivided(aggregate_budget):
    source = dict(r=[.5, .51], L=[.1, .13], c=[.2, .24])
    bottom = dict(r=[.5, .515], L=[.1, .08], c=[.2, .21])
    with pytest.raises(RectangleSolveError) as caught:
        solve_rectangle([0., .03], [.4, .42], source, bottom, alpha=.6, tube=TUBE,
                        budget=aggregate_budget, cell_options=dict(max_updates=1, atol=(1e-15,)*3, rtol=0.))
    assert caught.value.reason == "newton_update_limit"
    assert caught.value.diagnostics["completed_cells"] == 0
    assert caught.value.diagnostics["cell"]["updates"] == 1


@pytest.mark.parametrize("bad_axis", [[0., 0.], [0., .1, .21], [0.], [0., np.nan],
                                      [False, .1], np.linspace(0., 1., 66)])
def test_bad_axes_rejected_before_march(aggregate_budget, bad_axis):
    source = dict(r=[.5, .6], L=[0., 0.], c=[0., 0.])
    with pytest.raises(ValueError):
        solve_rectangle(bad_axis, [.4, .5], source, source, alpha=0., tube=TUBE, budget=aggregate_budget)


def test_invalid_boundary_corner_tube_and_options(aggregate_budget):
    U, x = axes(4)
    exact, _ = vacuum(U, x)
    source, bottom = boundaries(exact)
    for bad in ({**source, "extra": source["r"]}, {**source, "c": [0., 0.]},
                {**source, "r": np.zeros(len(x))}, {**source, "L": np.full(len(x), 1000.)}):
        with pytest.raises(ValueError):
            solve_rectangle(U, x, bad, bottom, alpha=0., tube=TUBE, budget=aggregate_budget)
    changed = copy.deepcopy(bottom)
    changed["r"][0] += 1e-6
    with pytest.raises(ValueError, match="corner"):
        solve_rectangle(U, x, source, changed, alpha=0., tube=TUBE, budget=aggregate_budget)
    with pytest.raises(ValueError, match="tube"):
        solve_rectangle(U, x, source, bottom, alpha=0., tube=dict(TUBE, r=(.5, 1.)), budget=aggregate_budget)
    with pytest.raises(ValueError, match="cell_options"):
        solve_rectangle(U, x, source, bottom, alpha=0., tube=TUBE, budget=aggregate_budget, cell_options=dict(region="cap"))
    with pytest.raises(ValueError, match="alpha"):
        solve_rectangle(U, x, source, bottom, alpha=-.1, tube=TUBE, budget=aggregate_budget)


def test_diagnostic_and_comparison_domains_are_checked():
    U, x = axes(8)
    states, _ = polynomial(U, x)
    row = grid_diagnostics(U, x, states, alpha=.3)
    with pytest.raises(ValueError, match="4..65"):
        grid_diagnostics(U[:3], x, {key: value[:3] for key, value in states.items()}, alpha=.3)
    with pytest.raises(ValueError, match="source fronts"):
        grid_diagnostics(U, x, states, alpha=.3, source_fronts=(.4, .45))
    with pytest.raises(ValueError, match="six"):
        grid_diagnostics(U, x, states, alpha=.3, source_jets={})
    with pytest.raises(ValueError, match="two-to-one"):
        compare_common_grids(row, row)
    U2, x2 = axes(16)
    states2, _ = polynomial(U2, x2)
    fine = grid_diagnostics(U2, x2, states2, alpha=.3)
    with pytest.raises(ValueError, match="both zero"):
        compare_common_grids(row, fine, coarse_zero=row)
    with pytest.raises(ValueError, match="zero baselines"):
        compare_common_grids(row, fine, coarse_zero=row, fine_zero=fine)
