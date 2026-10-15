"""Manufactured cumulative forcing checks; no physical field preparation."""
import time
from types import SimpleNamespace

import numpy as np
import pytest

from recursive_horizons import nsc_atlas_scalar_forcing as model
from recursive_horizons.nsc_atlas_bulk_response import BulkBudget, background
from recursive_horizons.nsc_atlas_field_preparation import radius_from_null
from recursive_horizons.nsc_atlas_prepared_response import PreparedFieldBudget


@pytest.fixture(autouse=True)
def manufactured_cpu_budget():
    started = time.process_time()
    yield
    assert time.process_time()-started < 30


@pytest.fixture
def geometry():
    # Fixed saved endpoints; no wall integration or prepared-field solve.
    U0, Ua, Ub = -.32130105902052203, -.2992702973503126, -.2747468030926288
    Va, Vb = -2.565695357068247, -2.262011613052105
    support = [Va+(Vb-Va)/4, Va+3*(Vb-Va)/4]
    return dict(order=64, U0=U0, Ua=Ua, Ub=Ub, Va=Va, Vb=Vb,
                x_corner=.4, x_last=float(radius_from_null(U0, Vb)),
                support_V=support, support_x=radius_from_null(U0, support).tolist())


class ManufacturedWedge:
    def __init__(self, geometry, nodes=17, scale=1.):
        self.domain = dict(geometry)
        self.grid = np.linspace(geometry["U0"], geometry["Ub"], nodes)
        self.budget = PreparedFieldBudget()
        self.prepared_field = SimpleNamespace(budget=self.budget)
        self.geometry = {"wall": self.wall}
        self.scale, self.calls = scale, []

    def wall(self, V):
        d = self.domain
        V = np.asarray(V)
        t = (V-d["Va"])/(d["Vb"]-d["Va"])
        return np.where(V == d["Va"], d["Ua"],
                        np.where(V == d["Vb"], d["Ub"], d["Ua"]+t*(d["Ub"]-d["Ua"])))

    def bulk(self, U, t, order=32):
        V = self.domain["Va"]+t*(self.domain["Vb"]-self.domain["Va"])
        assert U.ndim == 1 and 1 <= U.size <= 513
        assert np.all(U >= self.domain["U0"]) and np.all(U <= self.wall(V))
        assert order in (32, 64)
        self.calls.append(U.copy())
        self.budget.check(U.size)
        bg = model.background(U, V)
        q = U-self.domain["U0"]
        p = self.scale*(.3+q+q*q/2)
        pU, pV = self.scale*(1+q), np.full_like(U, .2*self.scale)
        return {"p": p, "p_U": pU, "p_V": pV,
                "chi": bg["r"]*p,
                "chi_U": bg["r_U"]*p+bg["r"]*pU,
                "chi_V": bg["r_V"]*p+bg["r"]*pV,
                "chi_U_lower_limit_bound": np.zeros_like(U)}


def polynomial_background(U, V):
    U = np.asarray(U)
    return {key: np.full_like(U, value) for key, value in
            (("r", .5), ("r_U", 2.), ("r_V", 3.), ("C", 4.))}


def install_source(monkeypatch, scale=1., *, p=.3, pV=.2):
    calls = []
    def source(V, geometry, *, order, source_scale, budget):
        assert source_scale == scale
        calls.append(V)
        budget.check(source=order)
        return dict(p=p*scale, p_V=pV*scale, mu_E=.7*scale**2,
                    zeta_E=-.2*scale**2, H_E=.4*scale**2, H_E_V=.8*scale**2)
    monkeypatch.setattr(model, "source_traces", source)
    return calls


def test_cumulative_polynomial_matches_all_prefixes_and_conservative_forcing(geometry, monkeypatch):
    monkeypatch.setattr(model, "background", polynomial_background)
    source_calls = install_source(monkeypatch)
    field = ManufacturedWedge(geometry)
    grid = np.linspace(geometry["U0"], geometry["Ub"], 12)
    forcing = model.ScalarForcing(field, geometry, grid, radial_order=2)
    row = forcing.slice(.63)
    q = row["U"]-geometry["U0"]
    integral = q+q*q+q**3/3
    mu, zeta = .7-3*integral/16, -.2+integral/8
    assert row["mu_E"] == pytest.approx(mu, abs=2e-15)
    assert row["zeta_E"] == pytest.approx(zeta, abs=2e-15)
    p, pU, pV = .3+q+q*q/2, 1+q, .2
    chi, chiU = .5*p, 2*p+.5*pU
    Z0, A0 = -1.6/row["V"], 1-.8/.5
    expected_F = Z0*(A0*zeta-2*mu/.5)*chiU/4
    expected_G = (4*(.4*zeta+mu)/(2*.5**3)-pU*pV/6)*chi
    assert row["F"] == pytest.approx(expected_F, abs=3e-15)
    assert row["G"] == pytest.approx(expected_G, abs=3e-15)
    assert row["K0"] == pytest.approx(np.full(grid.size, 6.4))
    curvature_omitted = 4*(.4*zeta+mu)/(2*.5**3)*chi
    assert np.max(abs(row["G"]-curvature_omitted)) > 1e-3
    assert len(source_calls) == 1
    assert row["source_p_defect"] == row["source_p_V_defect"] == 0
    assert np.all(np.diff(row["mu_E"]) <= 0)
    assert np.all(np.diff(row["zeta_E"]) >= 0)


def test_union_cells_clipping_and_unique_bounded_queries(geometry, monkeypatch):
    monkeypatch.setattr(model, "background", polynomial_background)
    install_source(monkeypatch)
    field = ManufacturedWedge(geometry, nodes=257)
    forcing = model.ScalarForcing(field, geometry, field.grid, radial_order=4)
    row = forcing.slice(1.)
    queried = np.concatenate(field.calls)
    assert np.array_equal(queried, np.unique(queried))
    assert row["radial_panels"] == 256
    assert row["radial_nodes"] == 1024
    assert row["field_points"] == queried.size == 1281
    assert row["field_batches"] == len(field.calls) == 3
    assert max(len(values) for values in field.calls) <= 513
    assert forcing.budget.metric_values == queried.size
    assert forcing.field_points == queried.size
    middle = forcing.slice(.31)
    inactive = field.grid >= middle["U_wall"]
    assert inactive.any()
    for key in ("U", "F", "G", "mu_E", "zeta_E", "K0"):
        assert np.all(middle[key][inactive] == middle[key][-1])
    assert all(np.all(points <= middle["U_wall"]) for points in field.calls[3:])
    assert forcing.slice_evaluations == 2


def test_quiet_control_is_exactly_zero_and_source_anchors_preserved(geometry, monkeypatch):
    monkeypatch.setattr(model, "background", polynomial_background)
    install_source(monkeypatch, 0.)
    field = ManufacturedWedge(geometry, scale=0.)
    forcing = model.ScalarForcing(field, geometry, field.grid, source_scale=0.)
    for t in (0., .5, 1.):
        row = forcing.slice(t)
        for key in ("F", "G", "mu_E", "zeta_E", "g1", "g1_V"):
            assert np.all(np.asarray(row[key]) == 0)
        assert row["mu_E"][0] == row["source_mu_E"]
        assert row["zeta_E"][0] == row["source_zeta_E"]


def test_one_ulp_wall_endpoint_keeps_the_actual_boundary_target(geometry, monkeypatch):
    monkeypatch.setattr(model, "background", polynomial_background)
    install_source(monkeypatch)
    field = ManufacturedWedge(geometry)
    old_wall = field.wall
    terminal = np.nextafter(geometry["Ub"], np.inf)
    field.wall = lambda V: np.where(np.asarray(V) == geometry["Vb"], terminal, old_wall(V))
    field.geometry["wall"] = field.wall
    row = model.ScalarForcing(field, geometry, field.grid).slice(1.)
    assert row["U_wall"] == row["U"][-1] == terminal
    assert np.concatenate(field.calls)[-1] == terminal
    assert np.all(row["U"] <= terminal)


def test_source_derivative_with_actual_background_constant_physical_scalar(geometry, monkeypatch):
    field = ManufacturedWedge(geometry)
    def constant_scalar(U, t, order=32):
        V = geometry["Va"]+t*(geometry["Vb"]-geometry["Va"])
        bg = background(U, V)
        return {"p": np.ones_like(U), "p_U": np.zeros_like(U), "p_V": np.zeros_like(U),
                "chi": bg["r"], "chi_U": bg["r_U"], "chi_V": bg["r_V"],
                "chi_U_lower_limit_bound": np.zeros_like(U)}
    field.bulk = constant_scalar
    def source(V, geometry, **kwargs):
        return dict(p=1., p_V=0., mu_E=0., zeta_E=0., H_E=V*V+V, H_E_V=2*V+1)
    monkeypatch.setattr(model, "source_traces", source)
    row = model.ScalarForcing(field, geometry, field.grid).slice(.6)
    # chi=r solves the background reduced wave equation, so chi_r=1 and
    # d[-H(V)*chi_r]/dV=-H'(V) without numerical differencing.
    assert row["source_chi_r"] == 1
    assert row["g1"] == pytest.approx(-row["V"]**2-row["V"])
    assert row["g1_V"] == pytest.approx(-2*row["V"]-1, abs=3e-15)


def test_transverse_source_memory_is_retained_after_zero_incident_trace(geometry, monkeypatch):
    monkeypatch.setattr(model, "background", polynomial_background)
    install_source(monkeypatch, p=0., pV=0.)
    field = ManufacturedWedge(geometry)
    def transverse_only(U, t, order=32):
        return {"p": np.zeros_like(U), "p_U": np.ones_like(U), "p_V": np.zeros_like(U),
                "chi": np.zeros_like(U), "chi_U": np.full_like(U, .5),
                "chi_V": np.zeros_like(U), "chi_U_lower_limit_bound": np.zeros_like(U)}
    field.bulk = transverse_only
    row = model.ScalarForcing(field, geometry, field.grid).slice(1.)
    assert row["source_chi"] == 0
    assert row["source_chi_r"] == .25
    assert row["g1"] == -.1
    assert row["g1_V"] == pytest.approx(-(.8+.4*row["source_P0_r"]/2)*.25)
    assert row["F"][0] != 0


@pytest.mark.parametrize("kwargs", [dict(radial_order=3), dict(field_order=16),
    dict(source_order=True), dict(source_scale=True), dict(source_scale=.5)])
def test_orders_and_source_control_fail_closed(geometry, kwargs):
    field = ManufacturedWedge(geometry)
    with pytest.raises(ValueError):
        model.ScalarForcing(field, geometry, field.grid, **kwargs)


@pytest.mark.parametrize("bad_grid", [np.array([[0., 1.]]), np.array([0., 1.]),
    np.array([np.nan, np.nan]), np.linspace(-.4, -.2, 514)])
def test_grid_validation(geometry, bad_grid):
    field = ManufacturedWedge(geometry)
    with pytest.raises(ValueError, match="grid"):
        model.ScalarForcing(field, geometry, bad_grid)


def test_shared_budget_domain_and_original_grid_validation(geometry):
    field = ManufacturedWedge(geometry)
    with pytest.raises(ValueError, match="shared"):
        model.ScalarForcing(field, geometry, field.grid, budget=BulkBudget())
    shifted = BulkBudget(started=field.budget.started-1, external=field.budget)
    with pytest.raises(ValueError, match="start"):
        model.ScalarForcing(field, geometry, field.grid, budget=shifted)
    field.domain["Ua"] += 1e-4
    with pytest.raises(ValueError, match="domains"):
        model.ScalarForcing(field, geometry, field.grid)
    field.domain = dict(geometry)
    field.grid[0] += 1e-4
    with pytest.raises(ValueError, match="field grid"):
        model.ScalarForcing(field, geometry, np.linspace(geometry["U0"], geometry["Ub"], 17))
    field.grid = np.linspace(geometry["U0"], geometry["Ub"], 17)
    field.prepared_field.budget = PreparedFieldBudget()
    with pytest.raises(ValueError, match="same budget"):
        model.ScalarForcing(field, geometry, field.grid)


@pytest.mark.parametrize("t", [True, np.nan, -.1, 1.1, np.array([.5])])
def test_time_domain(geometry, t):
    field = ManufacturedWedge(geometry)
    with pytest.raises(ValueError):
        model.ScalarForcing(field, geometry, field.grid).slice(t)


def test_exhausted_shared_budget_fails_before_queries(geometry):
    field = ManufacturedWedge(geometry)
    field.budget.started = time.process_time()-121
    with pytest.raises(RuntimeError, match="budget"):
        model.ScalarForcing(field, geometry, field.grid)
    assert field.calls == []


@pytest.mark.parametrize("bad_key", ["p", "chi_U", "chi_U_lower_limit_bound"])
def test_invalid_field_jets_fail_closed(geometry, monkeypatch, bad_key):
    monkeypatch.setattr(model, "background", polynomial_background)
    field = ManufacturedWedge(geometry)
    original = field.bulk
    def bad_field(U, t, order=32):
        row = original(U, t, order)
        row[bad_key] = np.full_like(U, np.nan)
        return row
    field.bulk = bad_field
    with pytest.raises(ValueError, match="field jet"):
        model.ScalarForcing(field, geometry, field.grid).slice(.5)


def test_source_gauge_mismatch_and_wall_exit_fail_closed(geometry, monkeypatch):
    monkeypatch.setattr(model, "background", polynomial_background)
    install_source(monkeypatch, p=1.)
    field = ManufacturedWedge(geometry)
    forcing = model.ScalarForcing(field, geometry, field.grid)
    with pytest.raises(ValueError, match="incompatible"):
        forcing.slice(.5)
    field.geometry["wall"] = lambda V: np.array(.1)
    with pytest.raises(ValueError, match="retained"):
        forcing.slice(.5)
