"""Explicit-callback time controls sharing the stage suite's single allowance."""
import copy
from functools import partial
import time

import numpy as np
import pytest

from recursive_horizons.nsc_atlas_coupled_evolution import CoupledEvolutionError
from recursive_horizons.nsc_atlas_background_ssprk3 import advance_background_ssprk3
from recursive_horizons.nsc_atlas_stage_ssprk3 import advance_stage_ssprk3, evolve_stage_ssprk3
from recursive_horizons.nsc_atlas_gradient_evolution import coupled_stage
from recursive_horizons.nsc_atlas_departure_gradient import first_order_departure_gradient
from recursive_horizons.nsc_atlas_reference_wall_gauge import ReferenceWallGauge
from test_nsc_atlas_reference_wall_gauge import GEOMETRY
from test_nsc_atlas_gradient_evolution import (
    gradient_control_budget, TUBE, profile, source, initial_state,
    assert_tree_equal, ExhaustionBudget,
)


ETA, Y = np.linspace(0., 1., 5), np.linspace(0., 1., 4)


def options(budget):
    return dict(eta=ETA, y_faces=Y, source=lambda x: dict(x=x),
        gauge=lambda x: dict(x=x, r=.4, T=.1, L=0., s=0.), alpha=0.,
        tube=TUBE, budget=budget, profile=None)


def ode_stage(rhs, calls):
    def evaluate(state, *, x, source, gauge, budget, **kwargs):
        calls.append((x, copy.deepcopy(state), source, gauge, budget))
        budget.check(1)
        rates = initial_state()
        rates["wall"]["s"] = rhs(x, float(state["wall"]["s"]))
        wall = {k: float(gauge[k])+float(state["wall"][k]) for k in ("r", "L", "T", "s")}
        wall["c"] = 0.
        return dict(rates=rates, courant_rate=1., physical_wall=wall,
            physical_cap=copy.deepcopy(state["cap"]), physical_increment=copy.deepcopy(state["departure"]))
    return evaluate


def test_exact_nonmonotone_cache_keys_and_quadratic_time_weights(gradient_control_budget):
    calls, lookups = [], []
    left, right = .2, .4
    mid = float(np.nextafter(left+(right-left)/2, 1.))
    keys = {left, right, mid}
    kw = options(gradient_control_budget)
    def src(x):
        assert x in keys
        lookups.append(("source", x))
        return dict(x=x)
    def ref(x):
        assert x in keys
        lookups.append(("gauge", x))
        return dict(x=x, r=.4, L=0., T=.1, s=0.)
    kw.update(source=src, gauge=ref)
    initial = initial_state(); initial["wall"]["s"] = .7
    row = advance_stage_ssprk3(initial, left, mid, right,
        stage=ode_stage(lambda x, s: x*x, calls), **kw)
    assert row["stage_times"] == [left, right, mid] == [c[0] for c in calls]
    assert lookups == [(name, x) for x in (left, right, mid) for name in ("source", "gauge")]
    assert row["accepted_state"]["wall"]["s"] == pytest.approx(.7+(right**3-left**3)/3, abs=2e-16)
    assert row["stage_calls"] == 3 and row["retries"] == 0
    assert all(c[-1] is gradient_control_budget for c in calls)
    assert initial["wall"]["s"] == .7


@pytest.mark.parametrize("kind", ["nonlinear", "mixed"])
def test_callback_rk_order_on_independent_nonlinear_odes(gradient_control_budget, kind):
    rhs = (lambda x, y: y*y) if kind == "nonlinear" else (lambda x, y: x*y)
    exact = .7/(1-.7*.3) if kind == "nonlinear" else .7*np.exp(.3**2/2)
    errors = []
    for n in (8, 16, 32):
        fine = np.linspace(0., .3, 2*n+1)
        calls = []
        value = initial_state(); value["wall"]["s"] = .7
        row = evolve_stage_ssprk3(value, fine[::2], mid_nodes=fine[1::2],
            stage=ode_stage(rhs, calls), save_every=n, **options(gradient_control_budget))
        errors.append(abs(row["snapshots"][-1]["state"]["wall"]["s"]-exact))
        assert row["stages"] == len(calls) == 3*n+1
        assert len(row["snapshots"]) == 2
    ratios = np.array(errors[:-1])/errors[1:]
    print("callback RK", kind, "errors", errors, "ratios", ratios.tolist())
    assert np.all((ratios > 7.3) & (ratios < 8.7))


def test_same_callback_supplies_final_snapshot_without_state_aliases(gradient_control_budget):
    calls = []
    row = evolve_stage_ssprk3(initial_state(), [0., .1, .2], mid_nodes=[.05, .15],
        stage=ode_stage(lambda x, y: 1., calls), **options(gradient_control_budget))
    assert [c[0] for c in calls] == [0., .1, .05, .1, .2, .15, .2]
    assert row["stages"] == 7 and row["steps"] == 2 and row["retries"] == 0
    assert [s["x"] for s in row["snapshots"]] == [0., .1, .2]
    for snapshot in row["snapshots"]:
        assert set(snapshot) == {"x", "state", "physical"}
        assert set(snapshot["physical"]) == {"increment", "wall", "cap"}
        assert snapshot["state"]["wall"]["s"] == pytest.approx(snapshot["x"], abs=1e-16)
    last = row["snapshots"][-1]
    last["physical"]["cap"]["Iq"][0] = 1.
    assert last["state"]["cap"]["Iq"][0] == 0.
    assert row["final_stage"]["physical_cap"]["Iq"][0] == 0.


def test_curved_reference_geometry_is_restored_at_actual_stage_time(gradient_control_budget):
    calls = []
    base = ode_stage(lambda x, y: 0., calls)
    def ref(x): return dict(r=.4+.1*x+.02*x**3, T=.1+.2*x+.05*x**3, L=0., s=0.)
    def stage(state, **kwargs):
        row = base(state, **kwargs)
        x = kwargs["x"]
        row["rates"]["wall"]["r"] = -.06*x*x
        row["rates"]["wall"]["T"] = -.15*x*x
        return row
    kw = options(gradient_control_budget); kw["gauge"] = ref
    row = advance_stage_ssprk3(initial_state(), 0., .02, .04, stage=stage, **kw)
    assert row["stages"][1]["physical_wall"]["r"] == pytest.approx(ref(.04)["r"], abs=1e-16)
    accepted = row["accepted_state"]["wall"]
    assert ref(.04)["r"]+accepted["r"] == pytest.approx(.404, abs=1e-16)
    assert ref(.04)["T"]+accepted["T"] == pytest.approx(.108, abs=1e-16)


def test_first_order_explicit_core_matches_frozen_three_stage_step(gradient_control_budget):
    budget = gradient_control_budget
    ref = ReferenceWallGauge(GEOMETRY, budget=budget)
    kw = dict(eta=ETA, y_faces=Y, source=source, gauge=ref.sample, alpha=.03,
        tube=TUBE, budget=budget, profile=profile())
    old = advance_background_ssprk3(initial_state(), .4, .4001, .4002, **kw)
    new = advance_stage_ssprk3(initial_state(), .4, .4001, .4002,
        stage=partial(coupled_stage, departure_gradient=first_order_departure_gradient), **kw)
    assert_tree_equal(new, old, budget)


def test_actual_quadratic_stage_callback_preserves_manufactured_quiet_vacuum(gradient_control_budget):
    budget = gradient_control_budget
    ref = ReferenceWallGauge(GEOMETRY, budget=budget)
    row = evolve_stage_ssprk3(initial_state(), [.4, .4002], mid_nodes=[.4001],
        stage=coupled_stage, eta=ETA, y_faces=Y, source=source, gauge=ref.sample,
        alpha=.03, tube=TUBE, budget=budget, profile=profile())
    assert row["stages"] == 4
    for snapshot in row["snapshots"]:
        for key in ("r", "L"):
            np.testing.assert_allclose(snapshot["state"]["departure"][key], 0., atol=4e-14)
        np.testing.assert_array_equal(snapshot["state"]["departure"]["c"], 0.)
        for key in ("Ip", "Iq"):
            np.testing.assert_array_equal(snapshot["state"]["cap"][key], 0.)
    assert "departure_gradients" in row["final_stage"]
    assert row["final_stage"]["charged_evaluations"] > 0


@pytest.mark.parametrize("failure,expected_calls,phase", [
    ("source", 0, "stage_0_inputs"), ("callback", 1, "stage_0_operator"),
    ("courant", 1, "stage_0_operator"), ("midpoint", 0, "validation"),
    ("not_callable", 0, "validation")])
def test_failure_counts_actual_attempts_and_keeps_initial_state(gradient_control_budget, failure, expected_calls, phase):
    calls = []
    callback = ode_stage(lambda x, y: 1., calls)
    kw = options(gradient_control_budget)
    mid = .05
    if failure == "source":
        def broken(x): raise KeyError("missing exact source key")
        kw["source"] = broken
    elif failure == "callback":
        def callback(state, **kwargs):
            kwargs["budget"].check(1)
            raise ValueError("deliberate stage failure")
    elif failure == "courant": kw["maximum_courant"] = .01
    elif failure == "midpoint": mid = .06
    elif failure == "not_callable": callback = None
    with pytest.raises(CoupledEvolutionError) as caught:
        advance_stage_ssprk3(initial_state(), 0., mid, .1, stage=callback, **kw)
    error = caught.value.diagnostic
    assert error["stage_calls"] == expected_calls and error["phase"] == phase
    assert error["charged_evaluations"] == gradient_control_budget.evaluations
    assert_tree_equal(error["accepted_state"], initial_state(), gradient_control_budget)


@pytest.mark.parametrize("kind", ["final_input", "final_callback", "final_nonfinite"])
def test_final_callback_failures_are_distinct_and_counted(gradient_control_budget, kind):
    calls, input_calls = [], []
    base = ode_stage(lambda x, y: 1., calls)
    def stage(state, **kwargs):
        if len(calls) == 3 and kind == "final_callback":
            kwargs["budget"].check(1)
            raise ValueError("final callback failure")
        result = base(state, **kwargs)
        if len(calls) == 4 and kind == "final_nonfinite": result["courant_rate"] = np.nan
        return result
    kw = options(gradient_control_budget)
    def src(x):
        input_calls.append(x)
        if len(input_calls) == 4 and kind == "final_input": raise KeyError("final cache key missing")
        return dict(x=x)
    kw["source"] = src
    with pytest.raises(CoupledEvolutionError) as caught:
        evolve_stage_ssprk3(initial_state(), [0., .1], mid_nodes=[.05], stage=stage, **kw)
    err = caught.value.diagnostic
    assert err["stage_calls"] == (3 if kind == "final_input" else 4)
    assert err["phase"] == ("final_inputs" if kind == "final_input" else "final_operator")
    assert err["accepted_state"]["wall"]["s"] == pytest.approx(.1, abs=1e-16)


def test_exhaustion_during_callback_counts_attempt_without_reset(gradient_control_budget):
    child = ExhaustionBudget(gradient_control_budget, limit=1)
    kw = options(child)
    with pytest.raises(CoupledEvolutionError) as caught:
        evolve_stage_ssprk3(initial_state(), [0., .1], mid_nodes=[.05],
            stage=ode_stage(lambda x, y: 1., []), **kw)
    err = caught.value.diagnostic
    assert err["stage_calls"] == 2
    assert err["local_failure"]["completed_stages"] == 1
    assert err["charged_evaluations"] == 2 == child.evaluations
    assert child.started == gradient_control_budget.started


@pytest.mark.parametrize("started", [True, -1., float("nan"), float("inf"), "future"])
def test_invalid_live_clock_rejected_before_stage_or_input_lookup(gradient_control_budget, started):
    child = ExhaustionBudget(gradient_control_budget)
    child.started = time.process_time()+100. if started == "future" else started
    calls = []
    with pytest.raises(CoupledEvolutionError, match="budget") as caught:
        advance_stage_ssprk3(initial_state(), 0., .05, .1,
            stage=ode_stage(lambda x, y: 1., calls), **options(child))
    assert calls == [] and caught.value.diagnostic["stage_calls"] == 0
