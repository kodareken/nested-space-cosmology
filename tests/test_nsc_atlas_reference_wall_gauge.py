"""Bounded reference controls; no physical source or regional evolution."""
import copy
import time

import numpy as np
import pytest
from scipy.special import lambertw

from recursive_horizons.nsc_atlas_characteristic_initialization import Budget
from recursive_horizons.nsc_atlas_causal_reference import reference_wall_events
from recursive_horizons.nsc_atlas_nonlinear_boundary import hybrid_wall_boundary
from recursive_horizons.nsc_atlas_reference_wall_gauge import ReferenceWallGauge


GEOMETRY = dict(U0=-.32130105902052203, Ua=-.2992702973503126,
    Ub=-.2747468030926288, Va=-2.565695357068247, Vb=-2.262011613052105,
    order=64, support_V=[-2.489774421064211, -2.3379325490561405],
    support_x=[.42269179059454887, .4631925883576189], x_corner=.4, x_last=.4815042082576364)


@pytest.fixture(scope="module")
def budget():
    shared = Budget(maximum_evaluations=20000, cpu_seconds=30.)
    yield shared
    shared.check()
    print(f"reference gauge controls: {shared.evaluations} charged evaluations, "
          f"{time.process_time()-shared.started:.6f} process CPU seconds")


@pytest.fixture(scope="module")
def gauge(budget):
    return ReferenceWallGauge(GEOMETRY, budget=budget)


@pytest.fixture(scope="module")
def regular_rows(gauge):
    return gauge.sample(np.linspace(.4, GEOMETRY["x_last"], 9))


@pytest.fixture(scope="module")
def independent_reference(budget, gauge):
    radii = np.array([.45, .5, .16**(1/3), .55])
    # Original reference: one anchor panel, 4+3+3+3 panels for
    # each of advanced time and proper time, plus four direct states.
    budget.check(27*64+4)
    reference = reference_wall_events(radii, order=64)
    source_radius = .8*(1+lambertw(-GEOMETRY["U0"]*reference["V"]/np.e).real)
    source_radius[[0, -1]] = [.4, GEOMETRY["x_last"]]
    return radii, reference, gauge.sample(source_radius)


def test_matches_original_reference_coordinates_and_proper_rates(independent_reference):
    radii, reference, row = independent_reference
    CJ0 = -3.2/GEOMETRY["U0"]
    np.testing.assert_allclose(row["r"], radii, rtol=0, atol=3e-15)
    for key in ("U", "V", "beta_minus", "beta_plus"):
        np.testing.assert_allclose(row[key], reference[key], rtol=2e-13, atol=3e-15)
    np.testing.assert_allclose(row["W"], CJ0*(reference["U"]-GEOMETRY["U0"]), atol=2e-14)
    np.testing.assert_allclose(row["Ubar_s"], CJ0*reference["ell_dU_dtau"], rtol=3e-13)
    V_x = -row["x"]*np.exp(row["x"]/.8)/(.64*GEOMETRY["U0"])
    np.testing.assert_allclose(row["x_s"], reference["ell_dV_dtau"]/V_x, rtol=3e-13)
    np.testing.assert_allclose(row["s"], reference["tau_wall_over_ell"]-reference["tau_wall_over_ell"][0], atol=3e-14)
    assert abs(row["beta_plus"][2]) < 3e-14
    assert np.all(np.isfinite(row["eta_plus_x"]))


def test_affine_vacuum_metric_and_total_derivative_identities(regular_rows):
    row = regular_rows
    q = 1-row["W"]/3.2
    product = q*(1-row["x"]/.8)*np.exp(row["x"]/.8)
    exact_r = .8*(1+lambertw(-product/np.e).real)
    np.testing.assert_allclose(row["r"], exact_r, atol=2e-15)
    np.testing.assert_allclose(row["r_x_total"], row["r_Ubar"]*row["W_x"]+row["r_x"], atol=5e-15)
    np.testing.assert_allclose(row["L_x_total"], row["log_C_Ubar"]*row["W_x"]+row["log_C_x"], atol=2e-14)
    mass = row["r"]*(1+4*row["r_Ubar"]*row["r_x"]/row["C"])/2
    np.testing.assert_allclose(mass, .4, atol=4e-16)
    np.testing.assert_allclose(row["C"]*row["W_x"], row["s_x"]**2, rtol=2e-14)
    assert np.max(abs(row["normalization_defect"])) < 2e-14
    np.testing.assert_allclose(row["Ubar_s"], np.exp(-row["eta_plus"])/np.sqrt(row["C"]), rtol=2e-14)
    np.testing.assert_allclose(row["x_s"], np.exp(row["eta_plus"])/np.sqrt(row["C"]), rtol=2e-14)
    np.testing.assert_allclose(row["r_s"], row["r_x_total"]/row["s_x"], rtol=2e-14)
    np.testing.assert_allclose(row["beta_plus"], -row["r_Ubar"]*row["Ubar_s"]+row["r_x"]*row["x_s"], atol=3e-15)


def test_hybrid_wall_law_independently_recovers_analytic_acceleration(independent_reference):
    _, _, row = independent_reference
    result = hybrid_wall_boundary(
        {"r": row["r"], "c": 0., "eta_plus": row["eta_plus"]},
        {"r_V": row["r_x"], "c_V": 0., "log_C_V": row["log_C_x"]},
        {"A": 1-row["r"]**2, "delta": 0., "w_minus": 0.},
        {"C_wall": row["C"], "log_C_s": row["L_x_total"]/row["s_x"]}, alpha=0.)
    for key in ("r_s", "T_s", "beta_minus", "beta_plus"):
        np.testing.assert_allclose(result[key], row[key], rtol=2e-13, atol=4e-15)
    np.testing.assert_allclose(result["eta_plus_s"], row["eta_plus_x"]/row["s_x"], rtol=2e-13)
    np.testing.assert_allclose(result["r_U_plus"], row["r_Ubar"], atol=2e-15)
    np.testing.assert_allclose(result["log_C_U_plus"], row["log_C_Ubar"], atol=2e-14)
    assert np.all(result["c_s"] == 0)


def test_clock_origin_and_quiet_corner(gauge):
    row = gauge.sample(.4)
    assert row["W"] == gauge.domain["W_a"]
    assert row["r"] == .45
    assert row["T"] == row["s"] == 0.
    assert row["c"] == row["c_Ubar"] == row["c_x"] == 0.
    np.testing.assert_allclose(row["T_x"], row["beta_minus"]*row["s_x"]/(1-row["r"]**2), rtol=0)
    assert "T=s=0" in gauge.summary()["clock_origin"]
    assert gauge.summary()["source_integrations"] == gauge.summary()["regional_marches"] == 0


def test_smooth_quadrature_refinement_preserves_derivatives(gauge, regular_rows, budget):
    coarse = ReferenceWallGauge(GEOMETRY, order=32, budget=budget)
    rows = coarse.sample(regular_rows["x"])
    for key in ("W", "W_x", "eta_plus", "eta_plus_x", "r", "r_x_total", "L_x_total", "T", "T_x", "s_x"):
        np.testing.assert_allclose(rows[key], regular_rows[key], rtol=3e-13, atol=2e-13)
    assert max(coarse.summary()["maximum_inverse_iterations"], gauge.summary()["maximum_inverse_iterations"]) <= 32
    assert abs(coarse.summary()["endpoint_integral_defect"]) < 3e-13
    assert gauge.summary()["comparison_is_rigorous_enclosure"] is False


def test_cache_reuses_quadrature_but_keeps_live_query_cost(gauge, regular_rows, budget):
    before = gauge.summary()
    parent_before = budget.evaluations
    result = gauge.sample(regular_rows["x"])
    result["W"][:] = -1  # Returned arrays do not mutate cached scalar rows.
    after = gauge.summary()
    assert after["quadrature_nodes"] == before["quadrature_nodes"]
    assert after["inverse_calls"] == before["inverse_calls"]
    assert after["cache_hits"]-before["cache_hits"] == 9
    assert budget.evaluations-parent_before == 9
    assert np.all(gauge.sample(regular_rows["x"])["W"] > 0)


def test_geometry_is_copied_and_not_reanchored(budget):
    supplied = copy.deepcopy(GEOMETRY)
    other = ReferenceWallGauge(supplied, budget=budget)
    supplied["Ua"] = 0.
    exposed = other.geometry
    exposed["support_x"][0] = 0.
    assert other.geometry == GEOMETRY
    assert other.sample(.4)["W"] == other.domain["W_a"]


@pytest.mark.parametrize("bad", [True, [True, .42], np.nan, np.inf, 1+0j, [],
    [[.42]], [.4]*514, np.nextafter(.4, -np.inf), np.nextafter(GEOMETRY["x_last"], np.inf)])
def test_rejects_nonreal_malformed_and_outside_queries(gauge, bad):
    with pytest.raises(ValueError):
        gauge.sample(bad)


@pytest.mark.parametrize("order", [True, 64., 16, 128])
def test_rejects_unadmitted_quadrature_orders(budget, order):
    with pytest.raises(ValueError, match="quadrature order"):
        ReferenceWallGauge(GEOMETRY, order=order, budget=budget)


def test_changed_wall_geometry_and_source_corner_rejected(budget):
    changed = copy.deepcopy(GEOMETRY)
    changed["Ua"] += 1e-5
    with pytest.raises(ValueError):
        ReferenceWallGauge(changed, budget=budget)
    changed = copy.deepcopy(GEOMETRY)
    changed["U0"] *= 1.001
    # Keep the old shape validator consistent to isolate the source corner.
    changed["x_last"] = float(.8*(1+lambertw(-changed["U0"]*changed["Vb"]/np.e).real))
    changed["support_x"] = (.8*(1+lambertw(-changed["U0"]*np.array(changed["support_V"])/np.e).real)).tolist()
    with pytest.raises(ValueError):
        ReferenceWallGauge(changed, budget=budget)


def test_shared_budget_has_no_private_lifetime_or_count_reset(budget, monkeypatch):
    # A larger caller ceiling is accepted; actual work still debits the same
    # aggregate test budget. Expiration applies even to an already cached row.
    original = budget.maximum_evaluations
    monkeypatch.setattr(budget, "maximum_evaluations", 60000)
    other = ReferenceWallGauge(GEOMETRY, budget=budget)
    assert other.budget is budget
    other.sample(.4)
    monkeypatch.setattr(budget, "maximum_evaluations", original)
    with monkeypatch.context() as patch:
        patch.setattr(budget, "started", time.process_time()-31.)
        with pytest.raises(RuntimeError, match="budget exhausted"):
            other.sample(.4)
    with monkeypatch.context() as patch:
        patch.setattr(budget, "maximum_evaluations", budget.evaluations)
        with pytest.raises(RuntimeError, match="budget exhausted"):
            other.sample(.4)


@pytest.mark.parametrize("field,value", [("cpu_seconds", float("inf")),
    ("maximum_evaluations", True), ("evaluations", -1), ("started", float("nan"))])
def test_rejects_invalid_budget_metadata(budget, monkeypatch, field, value):
    with monkeypatch.context() as patch:
        patch.setattr(budget, field, value)
        with pytest.raises(ValueError, match="live cumulative"):
            ReferenceWallGauge(GEOMETRY, budget=budget)
