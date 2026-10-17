"""Piecewise checks for the two-level nested slice. The full record stays frozen."""
import math
from fractions import Fraction

import numpy as np
import pytest

from recursive_horizons.nsc_atlas_connection_slice import schwarzschild_bridge
from recursive_horizons.nsc_atlas_nested_slice import (
    EvaluationBudget, bridge_amplitude, bridge_chart, flat_cap_state, gauss_integral,
    kottler_horizons, metric_lapse, metric_lapse_second, positive_sheet_state,
    positive_sheet_velocity, signed_slope, taper_velocity, tension_wall,
    transition_velocity, vacuum_slice_state,
)


def _scalar(value):
    return float(np.asarray(value).reshape(-1)[0])


def test_flat_cap_is_the_regular_de_sitter_sheet():
    hubble, radius = 1.7, 0.4
    state = flat_cap_state(radius, hubble)
    assert state["Rprime"] == 1
    assert state["velocity"] == pytest.approx(hubble*radius)
    assert state["theta"] == pytest.approx(3*hubble)
    assert state["scalar3"] == 0
    assert state["shear_squared"] == 0
    assert state["hamiltonian"] == 0
    assert state["momentum"] == 0
    assert state["misner_sharp"] == pytest.approx(hubble**2 * radius**3 / 2)
    centre = flat_cap_state(0.0, hubble)
    assert centre["radius"] == 0
    assert centre["Kll"] == -hubble
    assert centre["Ktheta"] == -hubble
    assert centre["misner_sharp"] == 0


def test_inner_wall_matches_rational_junction_arithmetic():
    radius, tension, mass_plus = Fraction(4, 5), Fraction(2), Fraction(1, 100)
    lapse_minus = 1-radius**2
    lapse_plus = 1-2*mass_plus/radius-radius**2
    beta = (lapse_minus-lapse_plus+(tension*radius)**2)/(2*tension*radius)
    assert beta == Fraction(517, 640)
    assert beta-tension*radius == Fraction(-507, 640)
    assert lapse_minus == Fraction(9, 25)
    motion = beta**2-lapse_minus
    predicted = ((1-1-tension**2)*radius**3/2 + tension*radius**2*beta)
    assert predicted == mass_plus
    wall = tension_wall(0.8, 0.0, 1.0, 0.01, 1.0, 2.0)
    assert wall["beta_minus"] == pytest.approx(float(Fraction(517, 640)), abs=1e-14)
    assert wall["beta_plus"] == pytest.approx(float(Fraction(-507, 640)), abs=1e-14)
    assert wall["Rdot"]**2 == pytest.approx(float(motion), abs=1e-14)
    assert wall["mean_trace"] == pytest.approx(0.0, abs=1e-12)
    assert wall["mass_map_residual"] == pytest.approx(0.0, abs=1e-12)
    assert wall["normal_force_residual"] == pytest.approx(0.0, abs=1e-12)
    assert wall["angular_jump"] == pytest.approx(-2.0, abs=1e-12)
    assert wall["temporal_jump"] == pytest.approx(-2.0, abs=1e-12)
    assert wall["trace_jump"] == pytest.approx(-6.0, abs=1e-12)
    for key in ("lorentz_residual", "minus_rate_residual", "minus_beta_residual",
                "plus_rate_residual", "plus_beta_residual"):
        assert wall[key] == pytest.approx(0.0, abs=1e-12)
    assert wall["Rprime_minus"] == 1
    assert wall["Rprime_plus"] < 0


def test_outer_wall_matches_rational_mass_and_mean_trace():
    radius, tension = Fraction(3, 2), Fraction(2)
    mass_minus, mass_plus = Fraction(1, 100), Fraction(2, 5)
    lapse_minus = 1-2*mass_minus/radius-radius**2
    lapse_plus = 1-2*mass_plus/radius
    beta = (lapse_minus-lapse_plus+(tension*radius)**2)/(2*tension*radius)
    assert beta == Fraction(727, 600)
    assert beta-tension*radius == Fraction(-1073, 600)
    assert lapse_minus == Fraction(-379, 300)
    predicted = (mass_minus + (1-0-tension**2)*radius**3/2 + tension*radius**2*beta)
    assert predicted == mass_plus
    wall = tension_wall(1.5, 0.01, 1.0, 0.4, 0.0, 2.0)
    assert wall["beta_minus"] == pytest.approx(float(Fraction(727, 600)), abs=1e-14)
    assert wall["beta_plus"] == pytest.approx(float(Fraction(-1073, 600)), abs=1e-14)
    assert wall["Rdot"]**2 == pytest.approx(float(beta**2-lapse_minus), abs=1e-12)
    assert wall["mean_trace"] == pytest.approx(-0.75, abs=1e-12)
    assert wall["mass_map_residual"] == pytest.approx(0.0, abs=1e-12)
    assert wall["normal_force_residual"] == pytest.approx(0.0, abs=1e-12)
    assert wall["temporal_jump"] == pytest.approx(-2.0, abs=1e-12)
    assert wall["trace_jump"] == pytest.approx(-6.0, abs=1e-12)
    assert wall["Rprime_plus"] < 0


def test_h0_bridge_matches_imported_schwarzschild_chart():
    mass = 0.4
    horizons = kottler_horizons(mass, 0.0)
    assert horizons["rb"] == 2*mass
    assert horizons["mass_identity_residual"] == 0
    parameters = np.array([-0.8, -0.25, 0.0, 0.55, 1.4])
    coordinate = math.sqrt(2*mass)*np.sinh(parameters)
    reference = schwarzschild_bridge(parameters, mass)
    chart = bridge_chart(coordinate, horizons["rb"], 0.0, mass)
    assert chart["radius"] == pytest.approx(reference["R"], abs=1e-12)
    assert chart["Rprime"] == pytest.approx(reference["Rprime"], abs=1e-12)
    end = float(coordinate[-1])

    def length_density(y):
        return bridge_chart(y, 2*mass, 0.0, mass)["dl_dy"]

    length = gauss_integral(length_density, 0.0, end, 64, EvaluationBudget())
    assert float(length[0]) == pytest.approx(float(reference["l"][-1]-reference["l"][2]), abs=1e-12)


def test_kottler_throat_is_regular_and_keeps_the_supplied_mass():
    mass, hubble = 0.01, 1.0
    horizons = kottler_horizons(mass, hubble)
    throat = horizons["rb"]
    polynomial = throat**3 - throat + 2*mass
    assert 2*mass < throat < 3*mass
    assert abs(polynomial) < 1e-12
    assert abs(horizons["mass_identity_residual"]) < 1e-12
    amplitude = _scalar(bridge_amplitude(throat, throat, hubble))
    assert amplitude > 0
    state = bridge_chart(0.0, throat, hubble, mass)
    assert state["Rprime"] == 0
    assert state["dl_dy"] == pytest.approx(2/math.sqrt(amplitude), rel=1e-12)
    assert state["Rsecond"] == pytest.approx(mass/throat**2 - hubble**2 * throat, abs=1e-12)
    assert abs(state["hamiltonian"]) < 1e-8
    assert abs(state["momentum"]) < 1e-8
    assert abs(state["Kll"]) < 1e-12 and abs(state["Ktheta"]) < 1e-12
    sample = 1e-6
    radius = throat+sample**2
    expected = math.sqrt((1-hubble**2 * (radius**2+radius*throat+throat**2))/radius)
    displaced = bridge_chart(sample, throat, hubble, mass)
    assert displaced["Rprime"]/sample == pytest.approx(expected, rel=1e-12)
    assert displaced["Rprime"]/sample == pytest.approx(math.sqrt(amplitude), rel=1e-8)
    with pytest.raises(ValueError):
        kottler_horizons(1.0, 1.0)


def test_expanding_sheet_is_regular_across_the_cosmological_root():
    mass, hubble = 0.01, 1.0
    root = kottler_horizons(mass, hubble)["rc"]
    assert 1/math.sqrt(3) < root < 1/hubble
    assert abs(root**3 - root + 2*mass) < 1e-12
    state = positive_sheet_state(root, mass, hubble)
    assert state["f"] == pytest.approx(0.0, abs=1e-12)
    assert state["Rprime"] == 1
    assert state["f"]+state["velocity"]**2 == pytest.approx(1.0, abs=1e-12)
    assert math.isfinite(state["Kll"]) and math.isfinite(state["Ktheta"])
    for shift in (-0.02, 0.02):
        nearby = positive_sheet_state(root+shift, mass, hubble)
        assert nearby["Rprime"] == 1
        assert abs(nearby["hamiltonian"]) < 1e-8
        assert math.isfinite(nearby["Kll"])


def test_length_rescaling_preserves_dimensionless_wall_and_roots():
    factor = 4.5
    inner = tension_wall(0.8, 0.0, 1.0, 0.01, 1.0, 2.0)
    scaled = tension_wall(0.8*factor, 0.0, 1/factor, 0.01*factor, 1/factor, 2/factor)
    for key in ("beta_minus", "beta_plus", "Rdot", "gamma", "v", "Rprime_plus", "V_minus"):
        assert scaled[key] == pytest.approx(inner[key], abs=1e-12)
    assert scaled["energy"]*factor == pytest.approx(inner["energy"], abs=1e-12)
    assert scaled["mean_trace"]*factor == pytest.approx(inner["mean_trace"], abs=1e-12)
    unit = kottler_horizons(0.01, 1.0)
    enlarged = kottler_horizons(0.01*factor, 1/factor)
    assert enlarged["rb"]/factor == pytest.approx(unit["rb"], rel=1e-12)
    assert enlarged["rc"]/factor == pytest.approx(unit["rc"], rel=1e-12)


def test_wall_deltas_use_one_over_gamma_and_slice_k_stays_bounded():
    wall = tension_wall(0.8, 0.0, 1.0, 0.01, 1.0, 2.0)
    jump = wall["Rprime_plus"]-1
    assert -16*math.pi*wall["radius"]*jump == pytest.approx(
        16*math.pi*wall["tension"]*wall["radius"]**2 * wall["gamma"], rel=1e-12)
    assert -4*jump/wall["radius"] == pytest.approx(16*math.pi*wall["energy"], rel=1e-12)
    assert 2*(wall["V_plus"]-wall["V_minus"])/wall["radius"] == pytest.approx(
        8*math.pi*wall["momentum"], rel=1e-12)
    metric = np.diag([-1.0, 1.0, 1.0, 1.0])
    normal = np.array([wall["gamma_v"], wall["gamma"], 0.0, 0.0])
    lowered = metric@normal
    tensor = -wall["sigma"]*(metric-np.outer(lowered, lowered))/wall["gamma"]
    assert tensor[0, 0] == pytest.approx(wall["energy"], rel=1e-12)
    assert -tensor[0, 1] == pytest.approx(wall["momentum"], rel=1e-12)
    assert abs(tensor[0, 0]-wall["sigma"]) > 1e-4
    velocity, velocity_radius = taper_velocity(wall["radius"], 0.1, wall["radius"], wall["V_plus"])
    exterior = vacuum_slice_state(
        wall["radius"], 0.01, 1.0, velocity, velocity_radius,
        signed_slope(wall["radius"], 0.01, 1.0, velocity, -1))
    interior = flat_cap_state(wall["radius"], 1.0)
    for state in (exterior, interior):
        assert abs(_scalar(state["Kll"])) < 20
        assert abs(_scalar(state["Ktheta"])) < 20
    assert _scalar(velocity) == pytest.approx(wall["V_plus"], abs=1e-12)
    assert _scalar(velocity_radius) == pytest.approx(0.0, abs=1e-12)


def test_transition_ends_match_the_static_bridge_and_expanding_sheet():
    low, high, mass, hubble = 0.5, 0.9, 0.01, 1.0
    stopped, stopped_radius = transition_velocity(low, low, high, mass, hubble)
    assert _scalar(stopped) == pytest.approx(0.0, abs=1e-15)
    assert _scalar(stopped_radius) == pytest.approx(0.0, abs=1e-15)
    matched, matched_radius = transition_velocity(high, low, high, mass, hubble)
    sheet, sheet_radius = positive_sheet_velocity(high, mass, hubble)
    assert _scalar(matched) == pytest.approx(_scalar(sheet), abs=1e-12)
    assert _scalar(matched_radius) == pytest.approx(_scalar(sheet_radius), abs=1e-12)
    state = vacuum_slice_state(high, mass, hubble, matched, matched_radius,
                               signed_slope(high, mass, hubble, matched, 1))
    assert abs(_scalar(state["hamiltonian"])) < 1e-8
    for radius, mass_value, hubble_value in ((0.1, 0.01, 1.0), (0.8, 0.01, 1.0),
                                             (1.0, 0.4, 0.0), (1.5, 0.4, 0.0)):
        assert _scalar(metric_lapse(radius, mass_value, hubble_value)) > 0
        assert _scalar(metric_lapse_second(radius, mass_value, hubble_value)) < 0


def test_gauss_ball_volume_agrees_at_the_two_admitted_orders():
    def area(radius):
        return 4*math.pi*np.asarray(radius, dtype=float)**2

    budget = EvaluationBudget()
    fine = gauss_integral(area, 0.0, 0.8, 64, budget)
    coarse = gauss_integral(area, 0.0, 0.8, 32, budget)
    exact = 4*math.pi*0.8**3/3
    assert float(fine[0]) == pytest.approx(exact, abs=1e-12)
    assert abs(float(coarse[0]-fine[0])) <= 1e-7*max(1.0, abs(float(fine[0])))
    assert budget.count == 96
    with pytest.raises(ValueError):
        gauss_integral(area, 0.0, 1.0, 16, EvaluationBudget())


def test_taper_distance_quadrature_stays_positive_and_close():
    wall = tension_wall(0.8, 0.0, 1.0, 0.01, 1.0, 2.0)

    def length_density(radius):
        velocity, velocity_radius = taper_velocity(radius, 0.1, 0.8, wall["V_plus"])
        slope = signed_slope(radius, 0.01, 1.0, velocity, -1)
        state = vacuum_slice_state(radius, 0.01, 1.0, velocity, velocity_radius, slope)
        assert abs(_scalar(np.max(np.abs(state["hamiltonian"])))) < 1e-8
        return 1/state["Rprime"]

    budget = EvaluationBudget()
    fine = gauss_integral(length_density, 0.8, 0.1, 64, budget)
    coarse = gauss_integral(length_density, 0.8, 0.1, 32, budget)
    assert float(fine[0]) > 0
    assert abs(float(coarse[0]-fine[0])) <= 1e-7*max(1.0, abs(float(fine[0])))
