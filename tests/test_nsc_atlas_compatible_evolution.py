"""Independent vacuum, common-field and fixed-step controls for the successor."""
import copy
import time

import numpy as np
import pytest

from recursive_horizons.nsc_atlas_characteristic_initialization import Budget
from recursive_horizons.nsc_atlas_reference_wall_gauge import ReferenceWallGauge
from recursive_horizons.nsc_atlas_compatible_evolution import (
    quiet_departures, compatible_stage, evolve_compatible, CoupledEvolutionError,
)
from test_nsc_atlas_reference_wall_gauge import GEOMETRY

TUBE = {"r": [.2, .8], "L": [-2., 2.], "c": [-4., 4.]}


@pytest.fixture(scope="module")
def budget():
    shared = Budget(maximum_evaluations=20000, cpu_seconds=30.)
    yield shared
    shared.check()
    print("compatible controls:", shared.evaluations, "points;",
          time.process_time()-shared.started, "process CPU seconds")


@pytest.fixture(scope="module")
def gauge(budget):
    return ReferenceWallGauge(GEOMETRY, budget=budget)


def source(x, scalar=False):
    return dict(r=x, log_C=0., c=0., r_x=1., log_C_x=0.,
                c_x=.7 if scalar else 0., source_constraint=0.)


def grids(n=8, nc=4):
    return np.linspace(0., 1., n+1), np.linspace(0., 1., nc+1)


def nonquiet(state, y, radius):
    state = copy.deepcopy(state)
    r = y*radius
    # Independently integrated regular polynomial p=.2*r, q=.3+.1*r^2.
    state["cap"]["Ip"] = .1*np.diff(r*r)
    state["cap"]["Iq"] = .3*np.diff(r)+.1*np.diff(r**3)/3
    n = len(state["departure"]["c"])+1
    state["departure"]["c"] = .15*np.linspace(0, 1, n+1)[1:-1]**2
    return state


def evaluate(state, gauge, budget, alpha=0., scalar=False, n=8, nc=4):
    eta, y = grids(n, nc)
    return compatible_stage(state, x=.4, eta=eta, y_faces=y, source=source(.4, scalar),
        gauge=gauge.sample(.4), alpha=alpha, tube=TUBE, budget=budget)


@pytest.mark.parametrize("n", [4, 8, 16])
def test_exact_vacuum_stage_has_no_geometric_departure(n, gauge, budget):
    result = evaluate(quiet_departures(exterior_cells=n, cap_cells=4), gauge, budget, n=n)
    for key in ("r", "L"):
        np.testing.assert_allclose(result["rates"]["departure"][key], 0., atol=3e-13)
    for rate in result["rates"]["wall"].values():
        assert abs(rate) < 3e-13
    for value in result["rates"]["cap"].values():
        np.testing.assert_array_equal(value, 0.)
    np.testing.assert_array_equal(result["rates"]["departure"]["c"], 0.)
    assert max(abs(v) for v in result["vacuum_wall_anchor_defect"].values()) < 3e-15


def test_finite_amplitude_quiet_limit(gauge, budget):
    state = quiet_departures(exterior_cells=8, cap_cells=4)
    zero, finite = evaluate(state, gauge, budget), evaluate(state, gauge, budget, alpha=.03)
    for group in zero["rates"]:
        for key in zero["rates"][group]:
            np.testing.assert_array_equal(zero["rates"][group][key], finite["rates"][group][key])


def test_zero_amplitude_geometry_independent_of_nonzero_scalar(gauge, budget):
    quiet = quiet_departures(exterior_cells=8, cap_cells=4)
    state = nonquiet(quiet, grids()[1], .45)
    a, b = evaluate(quiet, gauge, budget), evaluate(state, gauge, budget, scalar=True)
    for key in ("r", "L"):
        np.testing.assert_array_equal(a["rates"]["departure"][key], b["rates"]["departure"][key])
    for key in ("r", "L", "T", "s"):
        assert a["rates"]["wall"][key] == b["rates"]["wall"][key]
    assert np.max(abs(b["rates"]["cap"]["Ip"])) > 0


def test_single_wall_field_closes_both_frame_chain_rules(gauge, budget):
    state = nonquiet(quiet_departures(exterior_cells=8, cap_cells=4), grids()[1], .45)
    result = evaluate(state, gauge, budget, alpha=.01, scalar=True)
    wall, h = result["physical_wall"], result["hybrid"]
    chi = float(np.sum(state["cap"]["Iq"]))
    f = chi/wall["r"]
    expected_c = wall["r"]*np.sqrt(6/.01)*np.sinh(np.sqrt(.01/6)*f)
    assert wall["c"] == pytest.approx(expected_c, abs=2e-15)
    assert abs(result["canonical_wall_flux_defect"]) < 5e-14
    assert abs(result["conformal_wall_flux_defect"]) < 5e-14
    assert h["eta_plus_s"] == pytest.approx(float(gauge.sample(.4)["eta_plus_x"])*h["Vdot_plus"], abs=3e-14)
    assert max(abs(v) for v in result["rates"]["wall"].values()) > 1e-6


def test_numerical_wall_jet_defect_is_retained(gauge, budget):
    state = nonquiet(quiet_departures(exterior_cells=8, cap_cells=4), grids()[1], .45)
    row = evaluate(state, gauge, budget, alpha=.01, scalar=True)
    h, wx = row["hybrid"], float(gauge.sample(.4)["W_x"])
    required = {"r": h["r_U_plus"], "L": h["log_C_U_plus"], "c": h["c_U_plus"]}
    names = {"r": "r_Ubar", "L": "log_C_Ubar", "c": "c_Ubar"}
    for key in required:
        numerical = row["exterior_first_jets"][names[key]][-1]
        assert row["wall_kinematic_defect"][key] == pytest.approx(wx*(required[key]-numerical), abs=3e-14)
    assert abs(row["wall_kinematic_defect"]["c"]) > 1e-3


@pytest.mark.parametrize("scalar", [False, True])
def test_fixed_steps_preserve_moving_vacuum(gauge, budget, scalar):
    eta, y = grids()
    initial = quiet_departures(exterior_cells=8, cap_cells=4)
    if scalar:
        initial = nonquiet(initial, y, .45)
    x = np.linspace(.4, .42, 5)
    result = evolve_compatible(initial, x, eta=eta, y_faces=y,
        source=lambda at: source(at, scalar), gauge=gauge.sample, alpha=0.,
        tube=TUBE, budget=budget, save_every=2)
    assert result["steps"] == 4 and result["stages"] == 9 and result["retries"] == 0
    assert len(result["snapshots"]) == 3
    for snapshot in result["snapshots"]:
        for key in ("r", "L"):
            np.testing.assert_allclose(snapshot["state"]["departure"][key], 0., atol=3e-14)
        for value in snapshot["state"]["wall"].values():
            assert abs(value) < 3e-14
        reference = gauge.sample(snapshot["x"])
        for key in ("r", "L", "T", "s"):
            assert snapshot["physical"]["wall"][key] == pytest.approx(float(reference[key]), abs=3e-14)
    if not scalar:
        np.testing.assert_array_equal(result["snapshots"][-1]["state"]["cap"]["Iq"], 0.)


def test_tube_failure_retains_the_attempted_state(gauge, budget):
    eta, y = grids()
    initial = quiet_departures(exterior_cells=8, cap_cells=4)
    initial["wall"]["L"] = 3.
    with pytest.raises(CoupledEvolutionError) as caught:
        evolve_compatible(initial, [.4, .401], eta=eta, y_faces=y, source=source,
                          gauge=gauge.sample, alpha=0., tube=TUBE, budget=budget)
    diagnostic = caught.value.diagnostic
    assert diagnostic["phase"] == "first_stage" and diagnostic["completed_stages"] == 0
    assert diagnostic["attempt_state"]["wall"]["L"] == 3.


def test_cumulative_budget_cannot_be_reset(gauge):
    eta, y = grids()
    exhausted = Budget(maximum_evaluations=1, cpu_seconds=30.)
    with pytest.raises(CoupledEvolutionError):
        evolve_compatible(quiet_departures(exterior_cells=8, cap_cells=4), [.4, .401],
            eta=eta, y_faces=y, source=source, gauge=gauge.sample, alpha=0.,
            tube=TUBE, budget=exhausted)
