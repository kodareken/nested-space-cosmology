"""Manufactured integral/interpolation checks; the physical pilot is separate."""
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.interpolate import CubicHermiteSpline

from recursive_horizons import nsc_atlas_prepared_response as prepared_module
from recursive_horizons import nsc_atlas_wedge_field as wedge_module
from recursive_horizons.nsc_atlas_prepared_response import PreparedFieldBudget, _bulk_evaluator
from recursive_horizons.nsc_atlas_wedge_field import WedgeField


@pytest.fixture
def manufactured(monkeypatch):
    domain = {"U0": -.32130105902052203, "Ua": -.299270297350313,
              "Ub": -.274746803092629, "Va": -2.565695357068247,
              "Vb": -2.262011613052105}
    Va, Vb, Ua, Ub = (domain[k] for k in ("Va", "Vb", "Ua", "Ub"))
    slope = (Ub-Ua)/(Vb-Va)
    wall = CubicHermiteSpline([Va, Vb], [Ua, Ub], [slope, slope], extrapolate=False)
    geometry = {**domain, "wall": wall}
    grid = np.linspace(domain["U0"], Ub, 17)
    cuts = np.array([0., .25, .75, 1.])
    segments = tuple(lambda t: np.zeros((len(grid), np.size(t))) for _ in cuts[:-1])
    prepared = {"geometry": domain, "operator": {"analytic_norm_bound": .14}}
    budget = PreparedFieldBudget()
    field = SimpleNamespace(domain=domain, geometry=geometry, budget=budget,
                            dense_data={"grid": grid, "cuts": cuts, "segments": segments,
                                        "prepared": prepared})
    source = lambda V, Va, Vb: (np.ones_like(np.asarray(V)), np.zeros_like(np.asarray(V)))
    kernel = lambda x: np.full_like(np.asarray(x), 2.)
    for module in (prepared_module, wedge_module):
        monkeypatch.setattr(module, "source_profile", source)
        monkeypatch.setattr(module, "exterior_coefficient", kernel)
    field.bulk = _bulk_evaluator(grid, cuts, segments, geometry, prepared, budget)
    return field


def test_rectangle_overlap(manufactured):
    wedge = WedgeField(manufactured)
    U = np.linspace(wedge.U0, wedge.Ua, 11)[::-1]
    for t in (0., .25, .47, .75, 1.):
        old = manufactured.bulk(U, t, order=32)
        new = wedge.bulk(U, t, order=32)
        for key in old:
            np.testing.assert_allclose(new[key], old[key], rtol=2e-14, atol=2e-14)
        assert np.all(new["lower_V"] == wedge.Va)
        assert np.all(new["chi_U_lower_limit_bound"] == 0)


@pytest.mark.parametrize("order", [32, 64])
def test_birth_integral_and_wall_condition(manufactured, order):
    wedge = WedgeField(manufactured)
    t = .63
    V = wedge.Va+t*wedge.dV
    wall = float(wedge.geometry["wall"](V))
    U = np.array([wedge.U0, wedge.Ua, (wedge.Ua+wall)/2, wall])
    result = wedge.bulk(U, t, order=order)
    analytic_lower = np.where(U <= wedge.Ua, wedge.Va,
                             wedge.Va+(U-wedge.Ua)*wedge.dV/(wedge.Ub-wedge.Ua))
    np.testing.assert_allclose(result["chi_U"], -2*(V-analytic_lower), atol=3e-15)
    np.testing.assert_allclose(result["chi_V"], -2*(U-wedge.U0), atol=3e-16)
    assert result["chi_U"][-1] == 0
    assert result["p_U"][-1] != 0
    np.testing.assert_allclose(result["p_U"][-1], -result["p"][-1]*result["x_U"][-1]/result["x"][-1])
    assert np.all(result["lower_V_bracket_width"] >= 0)
    assert np.max(result["lower_V_bracket_width"]) < 3e-15
    assert np.all(result["chi_U_lower_limit_bound"] >= 0)


def test_paired_interpolation_matches_linear_field_across_partial_cells(manufactured):
    grid = manufactured.dense_data["grid"]
    Va, Vb = manufactured.domain["Va"], manufactured.domain["Vb"]
    wall = manufactured.geometry["wall"]

    def dense(t):
        t = np.atleast_1d(t)
        return .7*np.minimum(grid[:, None], wall(Va+(Vb-Va)*t)[None, :])*t[None, :]

    manufactured.dense_data["segments"] = (dense, dense, dense)
    wedge = WedgeField(manufactured)
    V = np.linspace(Va, Vb, 1500)
    U = wedge.U0+(wall(V)-wedge.U0)*np.linspace(0., 1., V.size)
    values = wedge._paired_values(U, V)
    np.testing.assert_allclose(values, 1+.7*U*(V-Va)/(Vb-Va), atol=4e-16)


def test_kernel_never_evaluated_across_wall(manufactured, monkeypatch):
    wedge = WedgeField(manufactured)
    original = wedge_module.radius_from_null
    seen = []

    def guarded(U, V):
        U, V = np.broadcast_arrays(U, V)
        assert np.all(U <= wedge.geometry["wall"](V))
        assert np.all(V >= wedge.Va) and np.all(V <= wedge.Vb)
        seen.append(U.size)
        return original(U, V)

    monkeypatch.setattr(wedge_module, "radius_from_null", guarded)
    t = .77
    wall = float(wedge.geometry["wall"](wedge.Va+t*wedge.dV))
    U = np.linspace(wedge.U0, wall, 43)
    result = wedge.bulk(U, t, order=64)
    assert seen and result["chi_U"][-1] == 0
    assert manufactured.budget.kernel_values <= 1+len(wedge.grid)+4*len(U)*64


def test_quiet_zero_field(manufactured, monkeypatch):
    monkeypatch.setattr(wedge_module, "source_profile", lambda V, *_: (np.zeros_like(V), np.zeros_like(V)))
    wedge = WedgeField(manufactured)
    result = wedge.bulk(np.linspace(wedge.U0, wedge.Ub, 9), 1.)
    for key in ("chi", "chi_U", "chi_V", "p", "p_U", "p_V"):
        assert np.all(result[key] == 0)


@pytest.mark.parametrize("points,t,order", [
    ([], .5, 32), ([[0]], .5, 32), ([np.nan], .5, 32),
    ([-.32], True, 32), ([-.32], -.1, 32), ([-.32], 1.1, 32),
    ([-.32], np.nan, 32), ([-.32], .5, True), ([-.32], .5, 48),
    ([-.32], .5, 32.), ([-.4], .5, 32), ([-.28], .1, 32),
])
def test_rejects_bad_queries(manufactured, points, t, order):
    with pytest.raises(ValueError):
        WedgeField(manufactured).bulk(points, t, order=order)


def test_rejects_spline_domain_and_dense_partition(manufactured):
    manufactured.dense_data["cuts"] = np.array([0., .25, .25, 1.])
    with pytest.raises(ValueError, match="partition"):
        WedgeField(manufactured)
    manufactured.dense_data["cuts"] = np.array([0., .25, .75, 1.])
    manufactured.geometry["wall"] = CubicHermiteSpline([-2.5, -2.3], [-.30, -.28], [.1, .1])
    with pytest.raises(ValueError, match="cover"):
        WedgeField(manufactured)


def test_shared_budget_limits_queries(manufactured):
    manufactured.budget = PreparedFieldBudget(kernel_values_max=1)
    wedge = WedgeField(manufactured)
    with pytest.raises(RuntimeError, match="budget"):
        wedge.bulk([wedge.U0], .5)


def test_solve_return_shapes_are_additive(manufactured, monkeypatch):
    def fake_solve(rate, interval, state, **kwargs):
        dense = lambda t: np.zeros((len(state), np.size(t)))
        return SimpleNamespace(success=True, t=np.array(interval),
                               y=np.zeros((len(state), 2)), sol=dense)

    monkeypatch.setattr(prepared_module, "solve_ivp", fake_solve)
    arguments = (129, manufactured.geometry, manufactured.dense_data["prepared"], manufactured.budget)
    assert len(prepared_module._solve(*arguments)) == 2
    assert len(prepared_module._solve(*arguments, include_bulk=True)) == 3
    result = prepared_module._solve(*arguments, include_dense=True)
    assert len(result) == 4 and callable(result[2])
    assert set(result[3]) == {"grid", "cuts", "segments", "prepared"}
