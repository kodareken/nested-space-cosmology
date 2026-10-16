"""Quiet-vacuum, simultaneous-interface and fixed-step controls."""
import copy
import time

import numpy as np
import pytest

from recursive_horizons.nsc_atlas_characteristic_initialization import Budget
from recursive_horizons.nsc_atlas_field_preparation import radius_from_null
from recursive_horizons.nsc_atlas_coupled_evolution import (
    canonical_field, quiet_initial_state, coupled_stage, evolve_coupled,
    CoupledEvolutionError,
)


TUBE = {"r": [.2, .8], "L": [-2., 2.], "c": [-4., 4.]}


@pytest.fixture(scope="module")
def budget():
    result = Budget(maximum_evaluations=20000, cpu_seconds=30.)
    yield result
    result.check()
    print("coupled manufactured evaluations:", result.evaluations,
          "process CPU:", time.process_time()-result.started)


def quiet(n=8, nc=4):
    # Independent exact Schwarzschild data at the fixed quiet wall event.
    x, R = .4, .45
    U0, Va, Ua = -.32130105902052203, -2.565695357068247, -.2992702973503126
    CJ = -3.2/U0
    W = CJ*(Ua-U0)
    eta = np.linspace(0, 1, n+1)
    U = U0+eta*W/CJ
    r = radius_from_null(U, Va)
    L = np.log(x/r)+(x-r)/.8
    Vx = -x*np.exp(x/.8)/(.64*U0)
    bp, bm = .2/R**2-1.25*R, .2/R**2+.75*R
    v = np.sqrt(bm*bm-(1-R*R))
    rU = -.64*Va*np.exp(-R/.8)/(R*CJ)
    rx = -.64*Ua*np.exp(-R/.8)*Vx/R
    C = np.exp(L[-1])
    xs = (v+bp)/(2*rx)
    Us = 1/(C*xs)
    Wx = Us/xs
    Rx = v/xs
    Lx = -(1/.8+1/R)*rx+1/x+1/.8
    totalLx = -(1/.8+1/R)*Rx+1/x+1/.8
    boostx = (-.4/R**3-1.25)/xs-Lx+totalLx/2
    source = dict(r=x, log_C=0., c=0., r_x=1., log_C_x=0., c_x=0., source_constraint=0.)
    bottom = dict(r=r, log_C=L, c=np.zeros_like(r))
    gauge = dict(W=W, W_x=Wx, eta_plus=-.5*np.log(Wx), eta_plus_x=boostx)
    state = quiet_initial_state(bottom, source, cap_cells=nc)
    expected = dict(r=Rx, L=totalLx, c=0., T=bm/((1-R*R)*xs), s=1/xs,
                    r_Ubar=rU, r_x=rx)
    return eta, np.linspace(0, 1, nc+1), source, gauge, state, expected


def stage(data, budget, alpha=0., diagnostics=False):
    eta, y, source, gauge, state, _ = data
    return coupled_stage(state, eta=eta, y_faces=y, source=source, gauge=gauge,
                         alpha=alpha, tube=TUBE, budget=budget, diagnostics=diagnostics)


def test_quiet_vacuum_wall_rates_and_spatial_refinement(budget):
    errors = []
    for n in (8, 16, 32):
        data = quiet(n)
        result = stage(data, budget, diagnostics=True)
        expected = data[-1]
        errors.append(abs(result["rates"]["wall"]["L"]-expected["L"]))
        for key in ("r", "T", "s"):
            assert result["rates"]["wall"][key] == pytest.approx(expected[key], abs=.0002)
        assert result["rates"]["wall"]["c"] == 0.
        for value in result["rates"]["cap"].values():
            np.testing.assert_array_equal(value, 0.)
        np.testing.assert_array_equal(result["rates"]["increment"]["c"], 0.)
        assert result["cap_diagnostics"]["max_abs_Phi_minus_fd_f_r"] == 0.
    assert 1.5 < errors[0]/errors[1] < 2.5
    assert 1.5 < errors[1]/errors[2] < 2.5


def nonquiet(data):
    data = copy.deepcopy(data)
    eta, _, source, _, state, _ = data
    state["wall"]["c"] = .15
    state["increment"]["c"] = .15*eta[1:-1]**2
    for key, value in (("f", .3), ("Pi", .2), ("Q", -.1)):
        state["cap"][key][:] = value
    source["c_x"] = .7
    return data


def test_zero_amplitude_geometry_independent_of_normalized_scalar(budget):
    a, b = stage(quiet(), budget), stage(nonquiet(quiet()), budget)
    for key in ("r", "L", "T", "s"):
        assert a["rates"]["wall"][key] == b["rates"]["wall"][key]
    for key in ("r", "L"):
        np.testing.assert_array_equal(a["rates"]["increment"][key], b["rates"]["increment"][key])
    assert b["rates"]["wall"]["c"] != 0
    assert np.max(abs(b["rates"]["cap"]["Pi"])) > 0


def test_actual_wall_lapse_rate_closes_acceleration_and_value_chain_rule(budget):
    data = nonquiet(quiet())
    row = stage(data, budget, alpha=.01)
    h = row["hybrid"]
    xs = h["Vdot_plus"]
    assert h["eta_plus_s"] == pytest.approx(data[3]["eta_plus_x"]*xs, abs=2e-14)
    assert h["Udot_plus"]/xs == pytest.approx(data[3]["W_x"], abs=2e-14)
    wall = data[4]["wall"]
    rate = row["rates"]["wall"]
    u = wall["c"]/wall["r"]
    field_rate = (rate["c"]-u*rate["r"])/(wall["r"]*np.sqrt(1+.01*u*u/6))
    assert field_rate == pytest.approx(h["f_s"]/xs, abs=1e-14)
    required = dict(r=h["r_U_plus"], L=h["log_C_U_plus"], c=h["c_U_plus"])
    numeric = dict(r=row["exterior_first_jets"]["r_Ubar"][-1],
                   L=row["exterior_first_jets"]["log_C_Ubar"][-1],
                   c=row["exterior_first_jets"]["c_Ubar"][-1])
    for key in ("r", "L", "c"):
        assert row["wall_kinematic_defect"][key] == pytest.approx(
            data[3]["W_x"]*(required[key]-numeric[key]), abs=2e-14)
    assert max(abs(v) for v in row["wall_kinematic_defect"].values()) > 1e-4


def test_cap_stage_uses_actual_common_radius_and_clock(budget):
    data = nonquiet(quiet())
    row = stage(data, budget, alpha=.01)
    assert row["cap_metric"]["radius"] == data[4]["wall"]["r"]
    np.testing.assert_array_equal(row["cap_metric"]["state"]["Pi"], data[4]["cap"]["Pi"])
    assert row["wall_f"] == pytest.approx(canonical_field(.15, .45, .01))
    assert row["wall_value_cell_difference"] == pytest.approx(row["wall_f"]-.3)


def test_fixed_step_integrator_stage_times_and_quiet_field(budget):
    data = quiet(4, 3)
    eta, y, source, gauge, state, _ = data
    calls = []
    def recorded(at):
        calls.append(at)
        return source
    xs = np.linspace(.4, .400002, 3)
    out = evolve_coupled(state, xs, eta=eta, y_faces=y, source=recorded,
        gauge=lambda _: gauge, alpha=0., tube=TUBE, budget=budget, save_every=2)
    # A manufactured frozen-coefficient input tests the driver, not a
    # physical source or reference trajectory.
    np.testing.assert_array_equal(calls, [xs[0], xs[1], xs[1], xs[2], xs[2]])
    assert out["stages"] == 5 and out["steps"] == 2 and out["retries"] == 0
    assert len(out["snapshots"]) == 2
    assert out["snapshots"][-1]["state"]["wall"]["T"] > 0
    for array in out["snapshots"][-1]["state"]["cap"].values():
        np.testing.assert_array_equal(array, 0.)
    np.testing.assert_array_equal(state["cap"]["f"], 0.)


def test_courant_stop_has_first_stage_context_without_retry(budget):
    eta, y, source, gauge, state, _ = quiet(4, 3)
    calls = []
    def sample(at):
        calls.append(at)
        return source
    with pytest.raises(CoupledEvolutionError) as failure:
        evolve_coupled(state, [.4, .48], eta=eta, y_faces=y, source=sample,
            gauge=lambda _: gauge, alpha=0., tube=TUBE, budget=budget, maximum_courant=.001)
    assert calls == [.4]
    assert failure.value.diagnostic["phase"] == "first_stage"
    assert failure.value.diagnostic["completed_stages"] == 1
    assert "Courant" in failure.value.diagnostic["reason"]


@pytest.mark.parametrize("change", ["boost", "W", "grid", "tube", "alpha", "shape"])
def test_incompatible_stage_inputs_rejected(change, budget):
    data = list(quiet())
    if change == "boost": data[3]["eta_plus"] += .2
    if change == "W": data[3]["W"] = 0.
    if change == "grid": data[0][2] += .01
    if change == "tube": data[4]["wall"]["r"] = 1.
    if change == "shape": data[4]["cap"]["f"] = np.zeros(2)
    with pytest.raises(ValueError):
        stage(data, budget, alpha=-1. if change == "alpha" else 0.)


def test_exhausted_budget_propagates_without_reset():
    b = Budget(maximum_evaluations=1, cpu_seconds=30.)
    with pytest.raises(RuntimeError):
        stage(quiet(), b)
    assert b.evaluations > b.maximum_evaluations


def test_common_corner_and_nonquiet_initial_data_rejected():
    _, _, source, _, _, _ = quiet()
    bottom = dict(r=np.linspace(.4, .45, 4), log_C=np.zeros(4), c=np.zeros(4))
    bottom["c"][-1] = .1
    with pytest.raises(ValueError, match="quiet"):
        quiet_initial_state(bottom, source, cap_cells=3)
    bottom["c"][-1] = 0.
    bottom["r"][0] = .41
    with pytest.raises(ValueError, match="corner"):
        quiet_initial_state(bottom, source, cap_cells=3)
