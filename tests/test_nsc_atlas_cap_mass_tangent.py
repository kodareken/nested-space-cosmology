"""Synthetic mass oracles, branch transitions and finite-volume identities."""
import copy

import numpy as np
import pytest

from recursive_horizons.nsc_atlas_characteristic_background import CharacteristicProfile
from recursive_horizons.nsc_atlas_background_cap import background_cap_metric, background_cap_rhs
from recursive_horizons import nsc_atlas_cap_mass_tangent as model
from test_nsc_atlas_cap_mass_derivatives import budget, report_allowance


def profile(knots=(-1., 0., 1.), zero=False):
    x = np.asarray(knots)
    a = 0. if zero else 1.
    return CharacteristicProfile(x, a*x**3, 3*a*x*x, 6*a*x)


def metric(budget, *, R=.4, T=.1, alpha=.02, source=None, state=None):
    if state is None:
        state = dict(Ip=np.array([.002, .005, .013, .025]),
                     Iq=np.array([.003, .001, -.0005, -.002]))
    return background_cap_metric(np.linspace(0, 1, 5), R, state, T=T,
        profile=source or profile(), alpha=alpha, quadrature_order=4, budget=budget)


def backward_oracle(base, direction, h, budget):
    values = []
    for delta in (-h, -2*h):
        row = background_cap_metric(base["y_faces"], base["radius"]+delta*direction["R"],
            {key: base["state"][key]+delta*direction[key] for key in ("Ip", "Iq")},
            T=base["T"]+delta*direction["T"], profile=base["profile"], alpha=base["alpha"],
            quadrature_order=base["quadrature_order"], budget=budget)
        values.append(row["n_faces"][-1])
    return (3*base["n_faces"][-1]-4*values[0]+values[1])/(2*h)


@pytest.mark.parametrize("args,winner", [
    ((0., 0., 2., 1.), "right"), ((0., 0., -2., -1.), "right"),
    ((0., 0., 2., -1.), "zero"), ((0., 0., 0., 1.), "zero"),
    ((0., 1., -2., 0.), "left"), ((0., 1., 2., 0.), "zero"),
    ((1., 1., 2., 1.), "left"), ((-1., -1., 2., 1.), "right"),
    ((-2., -1., 0., 0.), "right"), ((-1., 2., 1., 1.), "zero"),
])
def test_common_backward_minmod_branch(args, winner, budget):
    assert model.limiter_selection(*args) == winner
    budget.check(4)


def test_exact_polynomial_mass_and_component_contractions(budget):
    R, p, q, Rx = .45, .3, .2, .07
    r = R*np.linspace(0, 1, 5)
    state = dict(Ip=p*np.diff(r*r)/2, Iq=q*np.diff(r))
    base = metric(budget, R=R, alpha=0., source=profile(zero=True), state=state)
    direction = dict(R=Rx, T=1., Ip=2*Rx*state["Ip"]/R, Iq=Rx*state["Iq"]/R)
    row = model.mass_tangent(base, direction, budget=budget)
    mass = q*q*R**3/6+p*p*(R**3/12-R**5/20)
    rate = Rx*(q*q*R*R/2+p*p*(R*R/4-R**4/4))
    assert row["N"] == pytest.approx(mass, rel=2e-13)
    assert row["forward"] == pytest.approx(rate, rel=2e-12)
    assert row["reverse"] == pytest.approx(rate, rel=2e-12)
    assert row["Ip_weights"]@(state["Ip"]/p) == pytest.approx(2*p*(R**3/12-R**5/20), rel=2e-12)
    assert row["Iq_weights"]@(state["Iq"]/q) == pytest.approx(2*q*R**3/6, rel=2e-12)


def test_smooth_direction_matches_unchanged_metric_without_mutation(budget):
    base = metric(budget)
    before = copy.deepcopy(base["state"])
    direction = dict(R=.02, T=.3, Ip=np.array([.01, .03, .02, .04]), Iq=np.array([-.01, .02, .04, .01]))
    row = model.mass_tangent(base, direction, budget=budget)
    expected = backward_oracle(base, direction, 1e-5, budget)
    assert row["forward"] == pytest.approx(row["reverse"], rel=2e-13, abs=2e-13)
    assert row["reverse"] == pytest.approx(expected, rel=2e-6, abs=3e-7)
    for key in before: np.testing.assert_array_equal(base["state"][key], before[key])


def test_zero_departure_limiter_opens_on_the_full_backward_direction(budget):
    base = metric(budget, state=dict(Ip=np.zeros(4), Iq=np.zeros(4)))
    direction = dict(R=0., T=.2, Ip=np.array([.1, .2, .4, .8]), Iq=np.array([.03, .08, .11, .19]))
    row = model.mass_tangent(base, direction, budget=budget)
    assert any(value in {"left", "right"} for values in row["limiter_branches"].values() for value in values)
    expected = backward_oracle(base, direction, 1e-6, budget)
    assert row["reverse"] == pytest.approx(expected, rel=2e-6, abs=3e-7)


@pytest.mark.parametrize("clock_rate,owner", [(.2, 2), (-.2, 1)])
def test_coincident_profile_and_cell_cut_keeps_the_born_panel(budget, clock_rate, owner):
    R = 2*np.tanh(.2)
    base = metric(budget, R=R, T=0., source=profile((-1., 0., .2, 1.)))
    direction = dict(R=0., T=clock_rate, Ip=np.zeros(4), Iq=np.zeros(4))
    row = model.mass_tangent(base, direction, budget=budget)
    assert row["born_panels"] >= 1
    born = [p for p in row["panels"] if p["born"] and p["left"] == R/2]
    assert len(born) == 1 and born[0]["owner"] == owner
    expected = backward_oracle(base, direction, 1e-5, budget)
    assert row["forward"] == pytest.approx(row["reverse"], rel=2e-13, abs=2e-13)
    assert row["reverse"] == pytest.approx(expected, rel=2e-6, abs=3e-7)


def test_departure_summation_by_parts_and_incoming_trace_sensitivity(budget):
    base = metric(budget)
    Tx, Rx, incoming, delta = 1., .02, .7, .017
    rhs = background_cap_rhs(base, T_x=Tx, radius_x=Rx, incoming_w_plus=incoming, budget=budget)
    direction = dict(R=Rx, T=Tx, **rhs["rates"])
    tangent = model.mass_tangent(base, direction, budget=budget)
    ledger = model.mass_balance_ledger(base, rhs, tangent, R_x=Rx, T_x=Tx, budget=budget)
    assert rhs["q_flux"][0] == 0 and rhs["p_flux"][0] != 0
    assert ledger["centre"] != 0
    assert ledger["total"] == pytest.approx(tangent["forward"], rel=3e-12, abs=3e-13)
    for key in ("telescoping_defect", "interior_split_defect", "maximum_face_split_defect"):
        assert abs(ledger[key]) < 5e-13, key
    other = background_cap_rhs(base, T_x=Tx, radius_x=Rx, incoming_w_plus=incoming+delta, budget=budget)
    change = sum(tangent[key+"_weights"]@(other["rates"][key]-rhs["rates"][key]) for key in ("Ip", "Iq"))
    assert change == pytest.approx(delta*ledger["incoming_mass_rate_derivative"], rel=3e-12, abs=3e-13)


def test_primal_reference_drift_is_rejected(budget):
    base = metric(budget)
    base["n_faces"][2] += .01
    with pytest.raises(ValueError, match="primal metric"):
        model.mass_tangent(base, dict(R=0., T=1., Ip=np.zeros(4), Iq=np.zeros(4)), budget=budget)


def test_outward_backward_direction_at_upper_domain_is_rejected(budget):
    R, T = .4, .1
    top = T+np.arctanh(R)
    base = metric(budget, R=R, T=T, source=profile((-1., 0., top)))
    with pytest.raises(ValueError, match="upper profile domain"):
        model.mass_tangent(base, dict(R=0., T=-1., Ip=np.zeros(4), Iq=np.zeros(4)), budget=budget)


def test_shared_budget_exhaustion_in_tangent(budget):
    from recursive_horizons.nsc_atlas_characteristic_initialization import Budget
    base = metric(budget)
    child = Budget(maximum_evaluations=1, cpu_seconds=30)
    base["budget"] = child
    try:
        with pytest.raises(RuntimeError, match="budget"):
            model.mass_tangent(base, dict(R=0., T=1., Ip=np.zeros(4), Iq=np.zeros(4)), budget=child)
    finally:
        budget.check(child.evaluations)


def test_exact_radius_formula_switch_stops_before_field_evaluation(budget):
    z = model._Jet(np.array([.01]), np.zeros((1, 2)))
    with pytest.raises(ValueError, match="small-radius formula switch"):
        model._total(None, None, None, 0, z, budget)


@pytest.mark.parametrize("sign", [-1., 1.])
def test_exact_potential_formula_switch_stops(monkeypatch, budget, sign):
    f = model._Jet(np.array([sign*1e-4]), np.zeros((1, 2)))
    monkeypatch.setattr(model, "_total", lambda *a: {"f":f})
    with pytest.raises(ValueError, match="potential formula switch"):
        model._integrands(None, None, None, 0, None, 6., budget)
