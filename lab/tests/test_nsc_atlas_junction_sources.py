"""Independent controls for the finite source/geometry incompatibility screen."""

import copy
import json
import math

import pytest
import sympy as sp

from recursive_horizons.nsc_atlas_junction_sources import (
    constant_wall_polynomial,
    junction_record,
    radiation_junction,
)


def test_exact_polynomial_obstruction_in_both_cap_orientations():
    # Independently eliminate the imported angular equation without using
    # its implemented polynomial. k may be negative, zero, or positive.
    r, a, m, k, c = sp.symbols("r a m k c", real=True)
    s = sp.symbols("s", positive=True)
    beta_plus_squared = c*c + a*a*s**4/r**2 - 2*m/r
    necessary = sp.Poly(sp.expand(r*r*((c-k*r)**2-beta_plus_squared)), r)
    assert necessary.as_dict() == {(4,): k*k, (3,): -2*c*k, (1,): 2*m, (0,): -a*a*s**4}
    assert necessary.nth(0).subs(a, 1) == -s**4
    # Neither sign of the interior angular normal changes the obstruction.
    for orientation in (sp.Rational(1, 2), -sp.Rational(1, 2)):
        assert necessary.nth(0).subs(c, orientation) == -a*a*s**4


def test_bulk_friedmann_equations_and_no_false_bounce():
    for eta in (1.25, 1.35, math.pi/2, 1.8):
        q = radiation_junction(eta)
        adot = 1/math.tan(eta)
        addot = -1/q["a"]**3
        assert adot*adot+1 == pytest.approx(8*math.pi*q["rho_child"]*q["a"]**2/3)
        assert addot/q["a"] == pytest.approx(-4*math.pi*(q["rho_child"]+3*q["p_child"])/3)
    turning = radiation_junction(math.pi/2)
    assert turning["m_child_boundary"] == pytest.approx(math.sin(math.pi/3)**3/2)
    assert turning["m_child_boundary"] > 0


def test_required_stress_conserves_but_fixed_wall_fails():
    record = junction_record()
    assert record["status"] == "INCOMPATIBLE_FIXED_SOURCE"
    assert len(record["samples"]) == 33
    for row in record["samples"]:
        assert row["f_plus"] > 0
        assert row["A_sigma_required"] > 0
        assert row["A_pressure_required"] < 0
        assert abs(row["A2_codazzi_residual"]) < 1e-15
        assert row["A_eos_residual"] < -0.03
        assert abs(row["angular_residual"]) > 0.1
    assert len({row["A_sigma_wall"] for row in record["samples"]}) == 1
    assert len({row["A_sigma_required"] for row in record["samples"]}) == 33


def test_independent_metric_conservation_and_dust_controls():
    q = junction_record()["independent_checks"]
    assert q["metric_acceleration_max_abs_error"] < 2e-8
    assert q["finite_difference_codazzi_max_abs_error"] < 2e-10
    assert q["normal_einstein_jump_max_abs_error"] < 1e-14
    assert q["oppenheimer_snyder_max_abs_sigma"] < 1e-14
    assert q["oppenheimer_snyder_max_abs_pressure"] < 1e-14
    assert q["tested_cap_orientations"] == 2


def test_unit_rescaling_preserves_dimensionless_decision():
    ref = radiation_junction(1.35)
    for scale in (1e-6, 7.0, 1e6):
        q = radiation_junction(1.35, a_max=scale, mass=0.4*scale)
        for name in ("beta_minus", "beta_plus", "f_plus", "Rdot"):
            assert q[name] == pytest.approx(ref[name])
        assert q["R"]/scale == pytest.approx(ref["R"])
        assert q["tau"]/scale == pytest.approx(ref["tau"])
        for name in ("sigma", "pressure", "Ktau_plus", "Ktheta_plus"):
            assert q[name]*scale == pytest.approx(ref[name])
        p = constant_wall_polynomial(q["R"], a_max=scale, mass=.4*scale,
            chi=math.pi/3, sigma_wall=1/(64*math.pi*scale))/scale**2
        p_ref = constant_wall_polynomial(ref["R"], a_max=1, mass=.4,
            chi=math.pi/3, sigma_wall=1/(64*math.pi))
        assert p == pytest.approx(p_ref)


def test_large_child_is_a_distinct_normal_orientation():
    cap = radiation_junction(1.35)
    large = radiation_junction(1.35, chi=2*math.pi/3)
    assert large["R"] == pytest.approx(cap["R"])
    assert large["beta_minus"] == pytest.approx(-cap["beta_minus"])
    assert large["sigma"] < 0 < cap["sigma"]
    assert large["pressure"] != pytest.approx(cap["pressure"])
    assert abs(large["codazzi_residual"]) < 1e-15


def test_fixed_source_is_not_refitted_under_mass_perturbation():
    sigma_wall = 1/(64*math.pi)
    for mass in (.39, .4, .405):
        q = radiation_junction(1.35, mass=mass)
        residual = q["beta_minus"]-q["beta_plus"]-4*math.pi*sigma_wall*q["R"]
        assert abs(residual) > .05
        assert abs(q["pressure"]+q["sigma"]) > .03
    # These are source consistency checks, not perturbation stability.


@pytest.mark.parametrize("kwargs", [
    {"eta": 0.0}, {"eta": math.pi}, {"eta": float("nan")},
    {"eta": 1.35, "chi": 0}, {"eta": 1.35, "chi": math.pi},
    {"eta": 1.35, "a_max": 0}, {"eta": 1.35, "mass": -1},
    {"eta": .3},  # The comoving radiation clock exists but this exterior chart does not.
])
def test_domain_exit_is_rejected(kwargs):
    with pytest.raises(ValueError):
        radiation_junction(**kwargs)


def test_record_is_deterministic_and_cannot_be_mutated_across_calls():
    original = junction_record()
    assert json.dumps(original, sort_keys=True, allow_nan=False) == json.dumps(junction_record(), sort_keys=True, allow_nan=False)
    modified = copy.deepcopy(original)
    modified["theory"]["M_over_A"] = .1
    modified["samples"][0]["A_sigma_wall"] = 200
    assert junction_record() == original


def test_squared_condition_does_not_hide_an_orientation_choice():
    q = radiation_junction(1.35)
    # Angular matching at one point can be forced by defining its stress;
    # that does not supply an independent source or match the pressure law.
    at_one_point = q["sigma"]
    assert abs(constant_wall_polynomial(q["R"], a_max=1, mass=.4,
        chi=math.pi/3, sigma_wall=at_one_point)) < 1e-15
    assert abs(q["pressure"]+at_one_point) > .03
    elsewhere = radiation_junction(1.4)
    assert abs(constant_wall_polynomial(elsewhere["R"], a_max=1, mass=.4,
        chi=math.pi/3, sigma_wall=at_one_point)) > 1e-4
