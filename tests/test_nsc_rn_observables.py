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
    assert evolving["status"] == "sourced_time_dependent_jets"
    assert evolving["static_formula_applied"] is False
    assert evolving["trace_free_condition_imposed"] is False
    assert evolving["lapse_tr_used"] is False
    assert np.all(np.isfinite(evolving["R4"]))
    clock = 1.2
    reparameterized = obs.curvature_from_pg_jets(
        radii, clock * jets["N"], clock * jets["beta"],
        lapse_r=np.zeros(radii.size), shift_r=clock * jets["beta_r"],
        lapse_rr=np.zeros(radii.size), shift_rr=clock * jets["beta_rr"],
        lapse_t=np.full(radii.size, 0.3), shift_t=0.3 * jets["beta"],
        shift_tr=0.3 * jets["beta_r"],
    )
    assert np.max(np.abs(reparameterized["R4"] - closed["R4"])) < 1e-5
    assert np.max(np.abs(reparameterized["Ricci2"] - closed["Ricci2"])) / np.max(closed["Ricci2"]) < 1e-4


def test_nonstationary_metric_matches_an_independent_christoffel_scalar():
    radius = np.array([4.0])
    lapse = np.array([1.15])
    lapse_r = np.array([0.03])
    lapse_rr = np.array([-0.01])
    shift = np.array([0.4])
    shift_r = np.array([-0.02])
    shift_rr = np.array([0.004])
    lapse_t = np.array([0.05])
    shift_t = np.array([-0.01])
    shift_tr = np.array([0.008])
    measured = obs.curvature_from_pg_jets(
        radius, lapse, shift, lapse_r=lapse_r, shift_r=shift_r,
        lapse_rr=lapse_rr, shift_rr=shift_rr, lapse_t=lapse_t, shift_t=shift_t,
        shift_tr=shift_tr,
    )
    theta = 0.5 * np.pi
    step = 1.0e-5

    def metric(radius_value, theta_value, lapse_value, shift_value):
        chart = np.zeros((4, 4))
        chart[0, 0] = lapse_value ** 2 - shift_value ** 2
        chart[0, 1] = chart[1, 0] = -shift_value
        chart[1, 1] = -1.0
        chart[2, 2] = -radius_value ** 2
        chart[3, 3] = -radius_value ** 2 * np.sin(theta_value) ** 2
        return chart

    def state(radius_shift, time_shift):
        radius_value = float(radius[0] + radius_shift)
        time = time_shift
        lapse_value = float(
            lapse[0] + time * lapse_t[0] + radius_shift * lapse_r[0]
            + 0.5 * radius_shift ** 2 * lapse_rr[0]
        )
        shift_value = float(
            shift[0] + time * shift_t[0] + radius_shift * shift_r[0]
            + time * radius_shift * shift_tr[0] + 0.5 * radius_shift ** 2 * shift_rr[0]
        )
        return radius_value, lapse_value, shift_value

    def connection(radius_shift, time_shift, theta_value):
        radius_value, lapse_value, shift_value = state(radius_shift, time_shift)
        inverse = np.linalg.inv(metric(radius_value, theta_value, lapse_value, shift_value))
        samples = {}
        for label, item in (
            ("r", state(radius_shift + step, time_shift)),
            ("rm", state(radius_shift - step, time_shift)),
            ("t", state(radius_shift, time_shift + step)),
            ("tm", state(radius_shift, time_shift - step)),
        ):
            samples[label] = metric(item[0], theta_value, item[1], item[2])
        up = metric(radius_value, theta_value + step, lapse_value, shift_value)
        down = metric(radius_value, theta_value - step, lapse_value, shift_value)
        derivative = [
            (samples["t"] - samples["tm"]) / (2.0 * step),
            (samples["r"] - samples["rm"]) / (2.0 * step),
            (up - down) / (2.0 * step),
            np.zeros((4, 4)),
        ]
        gamma = np.zeros((4, 4, 4))
        for lam in range(4):
            for mu in range(4):
                for nu in range(4):
                    total = 0.0
                    for sig in range(4):
                        total += inverse[lam, sig] * (
                            derivative[mu][nu, sig] + derivative[nu][mu, sig]
                            - derivative[sig][mu, nu]
                        )
                    gamma[lam, mu, nu] = 0.5 * total
        return gamma

    def scalar_from(gamma, partial):
        riemann = np.zeros((4, 4, 4, 4))
        for a in range(4):
            for b in range(4):
                for c in range(4):
                    for d in range(4):
                        quadratic = 0.0
                        for mid in range(4):
                            quadratic += (
                                gamma[a, c, mid] * gamma[mid, d, b]
                                - gamma[a, d, mid] * gamma[mid, c, b]
                            )
                        riemann[a, b, c, d] = (
                            partial[c][a, d, b] - partial[d][a, c, b] + quadratic
                        )
        radius_value, lapse_value, shift_value = state(0.0, 0.0)
        inverse = np.linalg.inv(metric(radius_value, theta, lapse_value, shift_value))
        ricci = np.zeros((4, 4))
        for b in range(4):
            for d in range(4):
                ricci[b, d] = sum(riemann[a, b, a, d] for a in range(4))
        return float(sum(inverse[b, d] * ricci[b, d] for b in range(4) for d in range(4)))

    gamma0 = connection(0.0, 0.0, theta)
    partial = [
        (connection(0.0, step, theta) - connection(0.0, -step, theta)) / (2.0 * step),
        (connection(step, 0.0, theta) - connection(-step, 0.0, theta)) / (2.0 * step),
        (connection(0.0, 0.0, theta + step) - connection(0.0, 0.0, theta - step)) / (2.0 * step),
        np.zeros((4, 4, 4)),
    ]
    assert scalar_from(gamma0, partial) == pytest.approx(float(measured["R4"][0]), abs=1.0e-4)
    assert measured["trace_free_condition_imposed"] is False
    assert measured["lapse_tr_used"] is False


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
