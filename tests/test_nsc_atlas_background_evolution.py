"""Manufactured coupled controls; no prepared source or retained-history solve."""
import copy
import time

import numpy as np
import pytest

from recursive_horizons.nsc_atlas_characteristic_initialization import Budget
from recursive_horizons.nsc_atlas_characteristic_background import CharacteristicProfile
from recursive_horizons.nsc_atlas_reference_wall_gauge import ReferenceWallGauge
from recursive_horizons.nsc_atlas_background_evolution import (
    background_stage, evolve_background, quiet_departures, CoupledEvolutionError,
)
from test_nsc_atlas_reference_wall_gauge import GEOMETRY

TUBE = {"r": [.2, .8], "L": [-2., 2.], "c": [-4., 4.]}


@pytest.fixture(scope="module")
def budget():
    shared = Budget(maximum_evaluations=100000, cpu_seconds=30.)
    yield shared
    shared.check()
    print("background coupled controls:", shared.evaluations, "points;",
          time.process_time()-shared.started, "process CPU seconds")


@pytest.fixture(scope="module")
def gauge(budget):
    return ReferenceWallGauge(GEOMETRY, budget=budget)


def profile(scale=0.):
    labels = np.array([0., .3, .6, 1.3])
    return CharacteristicProfile(labels, scale*labels**3,
        3*scale*labels**2, 6*scale*labels)


def source(x, slope=0.):
    return dict(r=x, log_C=0., c=slope*(x-.4), r_x=1., log_C_x=0.,
                c_x=slope, source_constraint=0.)


def grids(n=8, nc=3):
    return np.linspace(0., 1., n+1), np.linspace(0., 1., nc+1)


def stage(state, gauge, budget, *, prof=None, alpha=0., slope=0., n=8, nc=3):
    eta, y = grids(n, nc)
    return background_stage(state, x=.4, eta=eta, y_faces=y,
        source=source(.4, slope), gauge=gauge.sample(.4), alpha=alpha,
        tube=TUBE, budget=budget, profile=profile() if prof is None else prof)


@pytest.mark.parametrize("n", [4, 8])
def test_quiet_stage_retains_exact_vacuum_and_zero_cap(n, gauge, budget):
    row = stage(quiet_departures(exterior_cells=n, cap_cells=3),
                gauge, budget, n=n)
    for key in ("r", "L"):
        np.testing.assert_allclose(row["rates"]["departure"][key], 0., atol=3e-13)
    for value in row["rates"]["wall"].values():
        assert abs(value) < 3e-13
    np.testing.assert_array_equal(row["rates"]["departure"]["c"], 0.)
    for key in ("Ip", "Iq"):
        np.testing.assert_array_equal(row["rates"]["cap"][key], 0.)
        np.testing.assert_array_equal(row["physical_cap"][key], 0.)
    assert row["cap_metric"]["budget"] is budget


def test_finite_amplitude_quiet_limit_uses_same_coupled_equations(gauge, budget):
    state = quiet_departures(exterior_cells=8, cap_cells=3)
    zero, finite = (stage(state, gauge, budget, alpha=a) for a in (0., .03))
    for group in zero["rates"]:
        for key in zero["rates"][group]:
            np.testing.assert_array_equal(zero["rates"][group][key], finite["rates"][group][key])


def test_active_profile_wall_flux_includes_background_and_incoming_difference(gauge, budget):
    alpha, amp = .01, .025
    initial = quiet_departures(exterior_cells=8, cap_cells=3)
    r = .45*grids()[1]
    initial["cap"]["Ip"] = .007*np.diff(r*r)/2
    initial["cap"]["Iq"] = .004*np.diff(r)
    row = stage(initial, gauge, budget, prof=profile(amp), alpha=alpha, slope=.07)
    wall, cap, hybrid = row["physical_wall"], row["cap_rhs"], row["hybrid"]
    R, T = wall["r"], wall["T"]
    plus, minus = max(T+np.arctanh(R), 0.), max(T-np.arctanh(R), 0.)
    chi0 = amp*(plus**3-minus**3)
    chi = chi0+sum(initial["cap"]["Iq"])
    f = chi/R
    expected_c = R*np.sinh(np.sqrt(alpha/6)*f)/np.sqrt(alpha/6)
    assert wall["c"] == pytest.approx(expected_c, abs=2e-15)
    assert sum(row["physical_cap"]["Iq"]) == pytest.approx(chi, abs=2e-15)
    p0 = 3*amp*(plus**2-minus**2)/(1-R*R)
    q0 = 3*amp*(plus**2+minus**2)/(1-R*R)
    baseline_rate = row["wall_rates"]["T"]*(1-R*R)*p0+row["wall_rates"]["r"]*q0
    assert abs(baseline_rate) > 1e-4
    assert cap["wall_chi_x"] == pytest.approx(baseline_rate+sum(cap["rates"]["Iq"]), abs=3e-15)
    assert abs(row["canonical_wall_flux_defect"]) < 3e-14
    assert abs(row["conformal_wall_flux_defect"]) < 3e-14
    assert cap["incoming_departure"] == pytest.approx(R*hybrid["w_plus"]+f-p0-q0, abs=2e-15)
    assert abs(cap["incoming_departure"]) > 1e-4
    assert np.max(abs(cap["rates"]["Iq"])) > 1e-4
    # The exterior first-order endpoint defect is measured, never imposed away.
    wx = float(gauge.sample(.4)["W_x"])
    assert row["wall_kinematic_defect"]["c"] == pytest.approx(
        wx*(hybrid["c_U_plus"]-row["exterior_first_jets"]["c_Ubar"][-1]), abs=3e-14)


@pytest.mark.parametrize("alpha", [0., .03])
def test_fixed_steps_keep_quiet_moving_vacuum(alpha, gauge, budget):
    eta, y = grids()
    initial = quiet_departures(exterior_cells=8, cap_cells=3)
    row = evolve_background(initial, np.linspace(.4, .42, 5), eta=eta,
        y_faces=y, source=source, gauge=gauge.sample, alpha=alpha,
        tube=TUBE, budget=budget, profile=profile(), save_every=2)
    assert row["steps"] == 4 and row["stages"] == 9 and row["retries"] == 0
    assert row["max_courant"] <= .4 and len(row["snapshots"]) == 3
    for sample in row["snapshots"]:
        for key in ("r", "L"):
            np.testing.assert_allclose(sample["state"]["departure"][key], 0., atol=4e-14)
        np.testing.assert_array_equal(sample["state"]["departure"]["c"], 0.)
        reference = gauge.sample(sample["x"])
        for key in ("r", "L", "T", "s"):
            assert sample["physical"]["wall"][key] == pytest.approx(float(reference[key]), abs=4e-14)
        for key in ("Ip", "Iq"):
            np.testing.assert_array_equal(sample["physical"]["cap"][key], 0.)
            np.testing.assert_array_equal(sample["state"]["cap"][key], 0.)


def test_active_snapshots_recover_total_cap_from_actual_clock_and_radius(gauge, budget):
    eta, y = grids()
    amp = .015
    initial = quiet_departures(exterior_cells=8, cap_cells=3)
    original = copy.deepcopy(initial)
    row = evolve_background(initial, np.linspace(.4, .402, 3), eta=eta,
        y_faces=y, source=lambda x: source(x, .04), gauge=gauge.sample,
        alpha=0., tube=TUBE, budget=budget, profile=profile(amp))
    assert row["steps"] == 2 and row["stages"] == 5
    for sample in row["snapshots"]:
        wall = sample["physical"]["wall"]
        R, T = wall["r"], wall["T"]
        xi = np.arctanh(R*y)
        Fplus, Fminus = amp*np.maximum(T+xi, 0.)**3, amp*np.maximum(T-xi, 0.)**3
        for key, primitive in (("Ip", Fplus+Fminus), ("Iq", Fplus-Fminus)):
            expected = np.diff(primitive)+sample["state"]["cap"][key]
            np.testing.assert_allclose(sample["physical"]["cap"][key], expected, rtol=2e-14, atol=2e-17)
        assert wall["c"] == pytest.approx(sum(sample["physical"]["cap"]["Iq"]), abs=2e-16)
        for key in ("r", "L"):
            np.testing.assert_allclose(sample["state"]["departure"][key], 0., atol=4e-14)
    assert np.max(abs(row["snapshots"][-1]["state"]["cap"]["Iq"])) > 1e-7
    for group in original:
        for key in original[group]:
            np.testing.assert_array_equal(initial[group][key], original[group][key])
    # Total snapshots and accepted raw state must not share mutable arrays.
    sample = row["snapshots"][-1]
    raw = sample["state"]["cap"]["Iq"].copy()
    sample["physical"]["cap"]["Iq"][0] += 1.
    np.testing.assert_array_equal(sample["state"]["cap"]["Iq"], raw)


@pytest.mark.parametrize("failure", ["profile", "tube", "courant"])
def test_failed_stage_retains_accepted_and_attempted_state(failure, gauge, budget):
    eta, y = grids()
    initial = quiet_departures(exterior_cells=8, cap_cells=3)
    prof, courant = profile(), .4
    if failure == "profile":
        prof = CharacteristicProfile([0., .01], [0., 0.], [0., 0.], [0., 0.])
    elif failure == "tube":
        initial["wall"]["L"] = 3.
    else:
        courant = 1e-10
    with pytest.raises(CoupledEvolutionError) as caught:
        evolve_background(initial, [.4, .401], eta=eta, y_faces=y, source=source,
            gauge=gauge.sample, alpha=0., tube=TUBE, budget=budget, profile=prof,
            maximum_courant=courant)
    error = caught.value.diagnostic
    assert error["phase"] == "first_stage" and error["step"] == 0
    assert error["completed_stages"] == (1 if failure == "courant" else 0)
    assert error["charged_evaluations"] == budget.evaluations
    for key in ("accepted_state", "attempt_state"):
        for field in ("Ip", "Iq"):
            np.testing.assert_array_equal(error[key]["cap"][field], initial["cap"][field])


def test_exhausted_cumulative_budget_does_not_restart(gauge):
    eta, y = grids()
    exhausted = Budget(maximum_evaluations=1, cpu_seconds=30.)
    with pytest.raises(CoupledEvolutionError, match="budget exhausted") as caught:
        evolve_background(quiet_departures(exterior_cells=8, cap_cells=3),
            [.4, .401], eta=eta, y_faces=y, source=source, gauge=gauge.sample,
            alpha=0., tube=TUBE, budget=exhausted, profile=profile())
    assert caught.value.diagnostic["charged_evaluations"] == exhausted.evaluations
    assert exhausted.evaluations > exhausted.maximum_evaluations
    assert caught.value.diagnostic["completed_stages"] == 0
