"""Manufactured/saved-source snapshot controls without a physical solve."""
import copy
import json
from pathlib import Path
import time

import numpy as np
import pytest

from recursive_horizons.nsc_atlas_characteristic_initialization import Budget
from recursive_horizons.nsc_atlas_characteristic_trace import SourceSnapshot
from recursive_horizons.nsc_atlas_nonlinear_source import source_amplitude_bound


ROOT = Path(__file__).resolve().parents[1]
GEOMETRY = dict(U0=-.32130105902052203, Ua=-.2992702973503126,
    Ub=-.2747468030926288, Va=-2.565695357068247, Vb=-2.262011613052105,
    order=64, support_V=[-2.489774421064211, -2.3379325490561405],
    support_x=[.42269179059454887, .4631925883576189], x_corner=.4, x_last=.4815042082576364)


@pytest.fixture(autouse=True)
def bounded_test():
    started = time.process_time()
    yield
    assert time.process_time()-started < 30


def exact_pulse(x):
    """Direct derivative of the original compact polynomial in V(x)."""
    U0, Va, Vb = (GEOMETRY[key] for key in ("U0", "Va", "Vb"))
    V = (1-x/.8)*np.exp(x/.8)/U0
    z = 4*(V-(Va+Vb)/2)/(Vb-Va)
    zx = -4*x*np.exp(x/.8)/(.64*U0*(Vb-Va))
    zxx = -4*np.exp(x/.8)*(1+x/.8)/(.64*U0*(Vb-Va))
    active = abs(z) < 1
    g = np.where(active, (1-z*z)**4, 0.)
    gz = np.where(active, -8*z*(1-z*z)**3, 0.)
    gzz = np.where(active, 8*(1-z*z)**2*(7*z*z-1), 0.)
    return g, gz*zx, gzz*zx*zx+gz*zxx


def manufactured(nodes=17, degree=3, alpha=None):
    x = np.linspace(.4, GEOMETRY["x_last"], nodes)
    alpha = source_amplitude_bound(GEOMETRY)["test_alpha"] if alpha is None else alpha
    t = x-.4
    h, hx = 100*t**degree, 100*degree*t**(degree-1)
    g, gx, _ = exact_pulse(x)
    y, yx = x+alpha*h, 1+alpha*hx
    F = 1-alpha*g*g/(6*y*y)
    r = np.sqrt(y*y-alpha*g*g/6)
    rx = (y*yx-alpha*g*gx/6)/r
    p, px = g/y, (gx-g*yx/y)/y
    a = (.8/x-1)/4
    rows = dict(x=x, h=h, h_x=hx, rJ=y, rJ_x=yx, r=r, C=F,
        log_C=np.log1p(-alpha*g*g/(6*y*y)), c=g, r_x=rx, c_x=gx,
        log_C_x=-alpha*p*px/(3*F), r_Ubar=a, c_Ubar=-t*t,
        log_C_Ubar=-(1/.8+1/x)*a)
    return alpha, x, rows


def make(nodes=17, degree=3, alpha=None, budget=None):
    alpha, x, rows = manufactured(nodes, degree, alpha)
    return SourceSnapshot(GEOMETRY, alpha=alpha, x_nodes=x, source_rows=rows, budget=budget)


def test_cubic_radius_correction_and_exact_pulse_at_withheld_points():
    snapshot = make()
    x = np.linspace(.4, GEOMETRY["x_last"], 65)
    row = snapshot.sample(x)
    t = x-.4
    np.testing.assert_allclose(row["h"], 100*t**3, atol=3e-17)
    np.testing.assert_allclose(row["h_x"], 300*t*t, atol=8e-15)
    np.testing.assert_allclose(row["h_xx"], 600*t, atol=8e-12)
    for key, expected in zip(("c", "c_x", "c_xx"), exact_pulse(x)):
        np.testing.assert_allclose(row[key], expected, rtol=2e-12, atol=1e-10)
    # h is deliberately not an ODE solution: do not hide the actual defect.
    assert np.max(abs(row["source_constraint"])) > 1e-3
    np.testing.assert_allclose(row["source_constraint"], row["source_constraint_jordan"], atol=2e-14)
    np.testing.assert_allclose(row["source_constraint"], row["source_constraint_from_curvature_defect"], atol=2e-14)
    assert np.max(abs(row["h_xx"]-row["source_ode_h_xx"])) > 1


def test_quartic_manufacture_refines_value_and_actual_curvature_independently():
    errors = []
    for nodes in (9, 17):
        snapshot = make(nodes=nodes, degree=4)
        x = (snapshot.x_nodes[:-1]+snapshot.x_nodes[1:])/2
        row = snapshot.sample(x)
        errors.append((np.max(abs(row["h"]-100*(x-.4)**4)),
                       np.max(abs(row["h_xx"]-1200*(x-.4)**2))))
    assert errors[0][0]/errors[1][0] == pytest.approx(16., rel=2e-8)
    assert errors[0][1]/errors[1][1] == pytest.approx(4., rel=2e-9)


def test_alpha_zero_preserves_normalized_pulse_and_vacuum_source():
    snapshot = make(alpha=0.)
    x = np.linspace(.4, GEOMETRY["x_last"], 49)
    row = snapshot.sample(x)
    np.testing.assert_array_equal(row["r"], x)
    np.testing.assert_array_equal(row["r_x"], np.ones_like(x))
    np.testing.assert_array_equal(row["C"], np.ones_like(x))
    np.testing.assert_array_equal(row["log_C_x"], np.zeros_like(x))
    assert np.max(row["c"]) > .9
    assert abs(row["c_Ubar"][-1]) > 0
    assert row["c"][-1] == row["c_x"][-1] == 0


def test_inputs_are_copied_and_query_modes_preserve_shape():
    alpha, x, rows = manufactured()
    geometry = copy.deepcopy(GEOMETRY)
    snapshot = SourceSnapshot(geometry, alpha=alpha, x_nodes=x, source_rows=rows)
    before = snapshot.sample(float(x[5]))
    x[:] = 0
    rows["h"][:] = 999
    geometry["support_x"][0] = 0
    exposed = snapshot.x_nodes
    exposed[:] = 0
    after = snapshot.sample(float(np.linspace(.4, GEOMETRY["x_last"], 17)[5]))
    assert all(np.array_equal(before[key], after[key]) for key in before)
    assert after["r"].shape == ()
    assert snapshot.sample([.4])["r"].shape == (1,)


@pytest.mark.parametrize("query", [True, [.4, True], float("nan"), float("inf"),
    [[.4]], np.linspace(.4, GEOMETRY["x_last"], 514),
    np.nextafter(.4, -np.inf), np.nextafter(GEOMETRY["x_last"], np.inf)])
def test_no_extrapolation_nonfinite_boolean_or_oversized_query(query):
    snapshot = make()
    before = snapshot.evaluations
    with pytest.raises(ValueError):
        snapshot.sample(query)
    assert snapshot.evaluations == before


@pytest.mark.parametrize("change", ["knots", "frame", "corner", "pulse", "alpha"])
def test_incompatible_source_data_is_rejected(change):
    alpha, x, rows = manufactured()
    if change == "knots":
        x[4] = x[3]
    elif change == "frame":
        rows["r_x"][8] += .01
    elif change == "corner":
        rows["r_Ubar"][0] += .01
    elif change == "pulse":
        rows["c"][8] += .01
    else:
        alpha = 2*alpha
    with pytest.raises(ValueError):
        SourceSnapshot(GEOMETRY, alpha=alpha, x_nodes=x, source_rows=rows)


def test_live_capture_queries_once_then_survives_expired_source_child(monkeypatch):
    import recursive_horizons.nsc_atlas_characteristic_trace as module
    now = [100.]
    monkeypatch.setattr(module.time, "process_time", lambda: now[0])
    parent = Budget(maximum_evaluations=60000, cpu_seconds=120., started=100.)
    class Child(Budget):
        def check(self, count=0):
            parent.check(count)
            super().check(count)
    child = Child(maximum_evaluations=20000, cpu_seconds=30., started=100.)
    child.parent = parent
    alpha, x, rows = manufactured()
    class Live:
        geometry, budget = GEOMETRY, child
        calls = 0
        def sample(self, requested):
            self.calls += 1
            self.budget.check(len(requested))
            return copy.deepcopy(rows)
    live = Live()
    live.alpha = alpha
    snapshot = SourceSnapshot.from_live(live, x, budget=parent)
    assert live.calls == 1
    now[0] = 140.
    with pytest.raises(RuntimeError, match="budget exhausted"):
        child.check()
    snapshot.sample(x)
    assert live.calls == 1 and child.evaluations == len(x)
    assert parent.evaluations == len(x)+snapshot.evaluations
    now[0] = 221.
    with pytest.raises(RuntimeError, match="budget exhausted"):
        snapshot.sample(x[0])


def test_local_and_shared_query_ceilings_are_enforced():
    shared = Budget(maximum_evaluations=60000, cpu_seconds=30.)
    snapshot = make(budget=shared)
    snapshot.evaluations = 20000
    with pytest.raises(RuntimeError, match="snapshot evaluation"):
        snapshot.sample(.4)
    assert snapshot.evaluations == 20001
    fresh = make(budget=shared)
    shared.evaluations = shared.maximum_evaluations
    with pytest.raises(RuntimeError, match="budget exhausted"):
        fresh.sample(.4)


def test_saved_source_knots_recover_values_without_integration_or_march(monkeypatch):
    import record_nsc_atlas_nonlinear_rectangle as record_tool
    import recursive_horizons.nsc_atlas_nonlinear_source as source_module
    import recursive_horizons.nsc_atlas_nonlinear_rectangle as rectangle_module
    def forbidden(*args, **kwargs):
        pytest.fail("physical integration during saved source reuse")
    monkeypatch.setattr(source_module, "ExactSource", forbidden)
    monkeypatch.setattr(rectangle_module, "solve_rectangle", forbidden)
    monkeypatch.setattr(rectangle_module, "solve_characteristic_cell", forbidden)
    path = ROOT/record_tool.RELATIVE_OUTPUT
    before = (path/"rectangle.json").read_bytes()
    record_tool.check(path)
    record = json.loads(before)
    data = record["initial_data"]["finite"]
    snapshot = SourceSnapshot(record["geometry"], alpha=data["alpha"],
        x_nodes=data["source"]["x"], source_rows=data["source"])
    row = snapshot.sample(snapshot.x_nodes)
    for key in ("h", "h_x", "r", "C", "log_C", "c", "r_x", "log_C_x", "c_x",
                "r_Ubar", "log_C_Ubar", "c_Ubar"):
        np.testing.assert_allclose(row[key], data["source"][key], rtol=2e-13, atol=2e-13)
    # Its polynomial curvature is honestly separate from the integrated source.
    assert np.max(abs(row["source_constraint"])) > 1e-8
    np.testing.assert_allclose(row["source_constraint"], row["source_constraint_from_curvature_defect"], atol=2e-14)
    summary = snapshot.summary()
    assert summary["source_integrations"] == summary["regional_marches"] == 0
    assert summary["source_numerical_error_is_enclosed"] is False
    assert summary["interpolation_error_is_enclosed"] is False
    assert before == (path/"rectangle.json").read_bytes()
