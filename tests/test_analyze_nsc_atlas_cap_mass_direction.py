"""Independent directional-mass controls; no retained physical analysis."""
import copy
import time

import numpy as np
import pytest

import analyze_nsc_atlas_cap_mass_direction as tool
from recursive_horizons.nsc_atlas_characteristic_initialization import Budget

_CONTROL_BUDGET = None
_REPORTED = False


@pytest.fixture(scope="session")
def budget():
    # Imported fixture registrations share one accumulator across all suites.
    global _CONTROL_BUDGET, _REPORTED
    if _CONTROL_BUDGET is None:
        _CONTROL_BUDGET = Budget(maximum_evaluations=300000, cpu_seconds=30.)
    yield _CONTROL_BUDGET
    _CONTROL_BUDGET.check()
    if not _REPORTED:
        print("cap mass direction controls", _CONTROL_BUDGET.evaluations,
              "points", time.process_time()-_CONTROL_BUDGET.started, "process CPU seconds")
        _REPORTED = True


def test_backward_polynomial_coefficients_and_one_sided_kink(budget):
    h = 1/32
    offsets = np.array([0., -2*h, -h, -h/2, -h/4])
    a, b, c, d = .3, -.2, .7, -.4
    values = 1+a*offsets+b*offsets**2+c*offsets**3+d*offsets**4
    before = values.copy()
    rows = tool.backward_derivatives(values, h, budget)
    assert len(rows) == 3
    for i, row in enumerate(rows):
        width = h/2**i
        assert row["width"] == width
        assert row["D_n"] == pytest.approx(a-2*c*width**2+6*d*width**3, rel=0, abs=3e-14)
    np.testing.assert_array_equal(values, before)
    kink = .3+2*offsets+3*np.abs(offsets)
    rows = tool.backward_derivatives(kink, h, budget)
    assert all(row["D_n"] == pytest.approx(-1.) for row in rows)
    # A stable backward estimate is explicitly different from the forward
    # derivative at this branch point; refinement alone does not remove it.
    assert rows[-1]["D_n"] != 5.


def test_moving_multi_input_functional_includes_geometry_background_and_state(budget):
    R, T = .5, .125
    state = dict(Ip=np.linspace(.001, .002, 64), Iq=np.linspace(-.003, .002, 64))
    rates = dict(R=.125, T=.25, Ip=np.full(64, .0001), Iq=np.full(64, -.0002))
    original = copy.deepcopy((state, rates))
    h = 1/64
    P, Q = (sum(state[k]) for k in ("Ip", "Iq"))
    Px, Qx = (sum(rates[k]) for k in ("Ip", "Iq"))
    exact = 2*R*rates["R"]*T+R*R*rates["T"]+Px*Q+P*Qx+3*T*T*rates["T"]
    cubic = rates["R"]**2*rates["T"]+rates["T"]**3
    values = []
    for offset in (0., -2*h, -h, -h/2, -h/4):
        probe = tool.linear_probe(R, T, state, rates, offset, budget)
        pr, pt = probe["R"], probe["T"]
        values.append(pr*pr*pt+sum(probe["state"]["Ip"])*sum(probe["state"]["Iq"])+pt**3)
    rows = tool.backward_derivatives(values, h, budget)
    for row in rows:
        assert row["D_n"] == pytest.approx(exact-2*cubic*row["width"]**2, rel=0, abs=1e-14)
    for source, before in zip((state, rates), original):
        for key in source:
            np.testing.assert_array_equal(source[key], before[key])
    assert exact != pytest.approx(Px*Q+P*Qx)
    with pytest.raises(ValueError):
        tool.linear_probe(R, T, state, rates, h, budget)


def test_stored_proper_rates_convert_once_to_x_label(budget):
    original = dict(r_s=.7, T_s=1.4, Vdot_plus=.35, w_plus=-.4)
    direction = tool.stored_direction(original, budget)
    assert direction == dict(R_x=2., T_x=4., w_plus=-.4)
    assert original == dict(r_s=.7, T_s=1.4, Vdot_plus=.35, w_plus=-.4)
    for key, value in (("Vdot_plus", 0.), ("T_s", -1.), ("w_plus", float("nan"))):
        with pytest.raises(ValueError):
            tool.stored_direction(dict(original, **{key:value}), budget)


def test_global_balance_and_history_bridge_keep_changed_rates_separate(budget):
    D_n, D_history, n_r, rho, J = 3.9, 3.7, .6, .55, .7
    Rh, Th, Rm, Tm = .4, 1.2, .45, 1.3
    history_n_T = (D_history-Rm*n_r)/Tm
    row = tool.balance_quantities(D_n, R_x=Rh, T_x=Th, n_r=n_r, rho_n=rho,
        J=J, history_n_T=history_n_T, history_R_x=Rm, history_T_x=Tm, budget=budget)
    assert row["e_m_dir"] == pytest.approx((D_n-Rh*n_r)/Th-J)
    assert row["B_mass"] == pytest.approx(D_n-Rh*rho-Th*J)
    assert row["e_r"] == pytest.approx(n_r-rho)
    assert row["identity_defect"] == pytest.approx(0., abs=1e-15)
    assert row["history_D_n"] == pytest.approx(D_history)
    assert row["e_hist"] == pytest.approx(history_n_T-J)
    assert row["e_hist_hybrid"] == pytest.approx((D_history-Rh*n_r)/Th-J)
    assert row["history_rate_bridge"] == pytest.approx(row["e_hist_hybrid"]-row["e_hist"])
    assert row["history_rate_bridge"] != 0.
    assert row["direction_minus_history_same_rates"] == pytest.approx((D_n-D_history)/Th)
    assert row["derivative_bridge_defect"] == pytest.approx(0., abs=1e-15)
    shifted = tool.balance_quantities(D_n, R_x=Rh, T_x=Th, n_r=n_r, rho_n=rho,
        J=J+.08, history_n_T=history_n_T, history_R_x=Rm, history_T_x=Tm, budget=budget)
    assert shifted["B_mass"]-row["B_mass"] == pytest.approx(-Th*.08)
    assert shifted["e_m_dir"]-row["e_m_dir"] == pytest.approx(-.08)
    assert shifted["direction_minus_history_same_rates"] == pytest.approx(row["direction_minus_history_same_rates"])


def test_limiter_ties_axis_override_and_branch_change_are_visible(budget):
    R = .5
    r = np.linspace(0, R, 65)
    dr = np.diff(r)
    midpoint = (r[:-1]+r[1:])/2
    plus, minus = midpoint, np.zeros(64)
    state = dict(Ip=dr*(plus+minus)/2, Iq=dr*(plus-minus)/2)
    before = copy.deepcopy(state)
    signature = tool.limiter_signature(R, state, budget)
    np.testing.assert_array_equal(signature["plus"]["secants"], np.ones(63))
    np.testing.assert_array_equal(signature["minus"]["secants"], np.zeros(63))
    assert signature["plus"]["selection"] == ["axis_override"]+["tie"]*62+["endpoint"]
    assert signature["minus"]["selection"] == ["axis_override"]+["zero"]*62+["endpoint"]
    changed_plus = plus.copy()
    changed_plus[32] = 0.
    changed = dict(Ip=dr*changed_plus/2, Iq=dr*changed_plus/2)
    next_signature = tool.limiter_signature(R, changed, budget)
    delta = tool.limiter_changes(signature, next_signature, budget)
    assert delta["plus"]["selection_changed_cells"]
    assert delta["plus"]["secant_sign_changed_indices"]
    assert not delta["minus"]["selection_changed_cells"]
    assert 0 not in delta["plus"]["selection_changed_cells"]
    for key in state:
        np.testing.assert_array_equal(state[key], before[key])


def test_exhaustion_is_charged_to_the_common_allowance(budget):
    class Child(Budget):
        def check(self, count=0):
            budget.check(count)
            super().check(count)
    child = Child(maximum_evaluations=1, cpu_seconds=30., started=budget.started)
    start = budget.evaluations
    with pytest.raises(RuntimeError, match="exhausted"):
        tool.backward_derivatives([1., 1., 1., 1., 1.], .01, child)
    assert budget.evaluations-start == child.evaluations > 1
