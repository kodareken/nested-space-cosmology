"""Reference-backed observer checks. No PG trajectory is launched."""
import os

for _name in (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
):
    os.environ[_name] = "1"

import numpy as np
import pytest

from recursive_horizons import nsc_rn_observables as obs
from recursive_horizons import nsc_rn_reference as reference


def _model(ratio=1.04, r_m=1.0):
    return obs.reference_model(ratio * r_m, r_m)


def test_vacuum_speeds_match_the_reference_and_clocks_follow_g_tt():
    model = _model()
    horizons = model.horizons
    excision = 0.5 * (horizons["r_minus"] + horizons["r_plus"])
    radii = np.linspace(excision, 8.0, 6)
    lapse = np.asarray(model.vacuum_lapse(radii), dtype=float)
    shift = np.asarray(model.beta_rn(radii), dtype=float)
    outgoing, ingoing = obs.chart_speeds(lapse, shift)
    reference_out, reference_in = model.characteristic_speeds(radii)
    np.testing.assert_allclose(outgoing, reference_out, atol=1e-12)
    np.testing.assert_allclose(ingoing, reference_in, atol=1e-12)
    assert bool(np.asarray(model.pure_outflow(excision)))
    with pytest.raises(obs.ObservableError):
        obs.static_killing_redshift(model.f(excision))
    assert obs.static_killing_redshift(model.f(8.0)) > 0.0
    assert obs.pg_normal_redshift(float(np.asarray(model.vacuum_lapse(excision)))) == pytest.approx(1.0)
    charged = obs.charged_mass(radii, lapse, shift, model.magnetic_r2)
    np.testing.assert_allclose(charged, model.charged_mass_pg(radii), atol=1e-12)
    np.testing.assert_allclose(charged, model.mass, atol=1e-12)


def test_reference_curvature_rows_and_endpoint_jet_are_not_repaired():
    model = _model()
    radii = np.linspace(3.0, 8.0, 5)
    rows = obs.curvature_from_reference(model, radii)
    assert rows["rows"] == 5 and rows["naive_D_at_D"] is False
    assert np.max(np.abs(rows["R4"])) == pytest.approx(0.0)
    analytic = np.asarray(model.beta_second_derivative(radii), dtype=float)
    halved = analytic.copy()
    halved[0] *= 0.5
    halved[-1] *= 0.5
    bias = obs.endpoint_second_derivative_bias(halved, analytic)
    assert bias["endpoint_half_actual"] is True
    assert bias["artifact"] == "D_at_D_half_actual"
    assert bias["repaired"] is False and bias["used_as_metric_jet"] is False
    intact = obs.endpoint_second_derivative_bias(analytic, analytic)
    assert intact["endpoint_half_actual"] is False
    jets = model.pg_metric_jets(radii)
    measured = obs.curvature_from_pg_jets(
        radii, jets["N"], jets["beta"],
        lapse_r=jets["N_r"], shift_r=jets["beta_r"],
        lapse_rr=np.zeros(radii.size), shift_rr=jets["beta_rr"],
        lapse_t=np.zeros(radii.size), shift_t=np.zeros(radii.size),
    )
    assert measured["status"] == "stationary_radial_jets"
    assert measured["static_formula_applied"] is False
    closed = model.invariants(radii)
    assert np.max(np.abs(measured["R4"] - closed["R4"])) < 1e-6
    assert np.max(np.abs(measured["Ricci2"] - closed["Ricci2"])) / np.max(closed["Ricci2"]) < 1e-5
    evolving = obs.curvature_from_pg_jets(
        radii, jets["N"], jets["beta"],
        lapse_r=jets["N_r"], shift_r=jets["beta_r"],
        lapse_rr=np.zeros(radii.size), shift_rr=jets["beta_rr"],
        lapse_t=np.full(radii.size, 1e-3), shift_t=np.zeros(radii.size),
    )
    assert evolving["status"] == "unresolved_time_dependent_curvature"
    assert evolving["R4"] is None and evolving["static_formula_applied"] is False


def test_gram_and_car_do_not_apply_multiplicity_twice():
    shell = reference.manufactured_quadrature_control(
        center=6.0, width=0.4, r_inner=4.0, r_outer=10.0, points=32,
        occupations=np.array([0.5, 0.25]),
    )
    assert shell["phi"].shape == (2, 32, 2)
    assert shell["multiplicity"] == 4 and shell["multiplicity_counted_once"] is True
    assert shell["positivity_unestablished"] is True
    assert shell["prepared_positive_energy_source"] is False
    assert shell["exterior_killing_frequency"] is None
    accounts = obs.shell_accounts(
        shell["phi"], shell["weights"], shell["radius"],
        killing_frequency=0.3,
    )
    density = np.sum(np.abs(shell["phi"]) ** 2, axis=(0, 2))
    car = float(np.dot(shell["weights"], density / shell["radius"]))
    assert accounts["probability_remaining"] == pytest.approx(np.real(np.trace(shell["gram"])))
    assert accounts["car_angular"] == pytest.approx(car)
    assert accounts["car_angular"] != pytest.approx(4.0 * car)
    assert accounts["multiplicity_in_gram"] is False
    assert accounts["multiplicity_in_H"] is False
    assert accounts["multiplicity_applications_in_car"] == 1
    assert accounts["filled_sea"] is False and accounts["V4"] == 0.0
    assert accounts["occupations_normalized_to_survivors"] is False


def test_ledger_accounts_stay_distinct_on_a_reference_excision():
    model = _model()
    horizons = model.horizons
    excision = 0.5 * (horizons["r_minus"] + horizons["r_plus"])
    radii = np.linspace(excision, 8.0, 8)
    weights = reference.radial_sbp_weights(excision, 8.0, 8)
    lapse = np.asarray(model.vacuum_lapse(radii), dtype=float)
    shift = np.asarray(model.beta_rn(radii), dtype=float)
    phi = np.zeros((2, 8, 1), dtype=complex)
    phi[0, 0, 0] = 1.0
    phi[0, -1, 0] = 1.0
    outgoing, _ingoing = model.characteristic_speeds(radii)
    ledger = obs.boundary_ledger(
        phi=phi, weights=weights, radius=radii, lapse=lapse, shift=shift,
        killing_frequency=0.4, metric_tt=np.asarray(model.f(radii), dtype=float),
        sat_debit=0.25, sat_enabled=True, horizon_radius=float(radii[3]),
        horizon_radius_rate=0.1,
    )
    assert ledger["excision_outflow"] == pytest.approx(-float(outgoing[0]))
    assert ledger["outer_outflow"] == pytest.approx(float(outgoing[-1]))
    assert ledger["sat_debit"] == pytest.approx(0.25)
    assert ledger["sat_added_into_outflow"] is False
    assert ledger["pure_outflow_excision"] is True
    assert ledger["multiplicity_in_H"] is False
    assert "normalized_survivors" not in ledger
    bare = obs.boundary_ledger(
        phi=phi, weights=weights, radius=radii, lapse=lapse, shift=shift,
        killing_frequency=0.4, metric_tt=np.asarray(model.f(radii), dtype=float),
    )
    assert bare["sat_debit"] == 0.0
    assert bare["outer_outflow"] == ledger["outer_outflow"]
    assert bare["excision_outflow"] == ledger["excision_outflow"]
    moved = obs.horizon_transport(2.0, 0.5, 0.25)
    still = obs.horizon_transport(2.0, 4.0, 0.25)
    assert moved["horizon_reynolds"] == pytest.approx(0.5)
    assert still["horizon_reynolds"] == moved["horizon_reynolds"]
    assert still["horizon_fixed_flux"] == pytest.approx(4.0)


def test_one_enclosed_mass_uses_the_reference_and_keeps_the_sample_lapse():
    model = _model()
    radius = 6.0
    initial = model.pg_metric_jets(radius)
    sample = {
        "radius": radius,
        "lapse": 1.02,
        "shift": float(np.asarray(initial["beta"]).reshape(-1)[0]) + 0.01,
    }
    compared = obs.compare_station(sample, model, model.mass + 0.02)
    assert compared["per_radius_refit"] is False
    assert compared["coupled_lapse_replaced"] is False
    assert compared["sample_lapse"] == pytest.approx(1.02)
    assert compared["initial_rn_lapse"] == pytest.approx(1.0)
    assert compared["matched_rn_mass"] == pytest.approx(model.mass + 0.02)
    assert compared["magnetic_r2"] == pytest.approx(model.magnetic_r2)
    with pytest.raises(obs.ObservableError):
        obs.compare_station(sample, model, np.array([model.mass, model.mass + 0.01]))


def test_throat_grid_and_full_cfl():
    spacing = obs.throat_spacing(0.16)
    assert spacing["spacing"] == pytest.approx(0.01)
    assert spacing["derivative_ladder"] is False
    assert "lambda_min" in spacing["pending_scales"]
    tightened = obs.throat_spacing(0.16, lambda_min=0.12, sigma=0.48)
    assert tightened["spacing"] == pytest.approx(min(0.01, 0.005, 0.02))
    assert obs.route_points(1.0, 0.1, 32) == 32
    assert obs.route_points(10.0, 0.1, 64) >= 64
    assert obs.full_cfl_dt(1.0, 0.5, 0.1) == pytest.approx(0.1 / 1.5)
    with pytest.raises(obs.ObservableError):
        obs.require_phi(np.ones((2, 4)))
