"""Oracle comparisons for one finite channel-turnover witness.

References here are independent of the calculated record. A missing module
fails; it is not a skipped success. The physical gate is not closed.
"""
from __future__ import annotations

import json
import math
from fractions import Fraction

import pytest
import sympy as sp

from recursive_horizons.nsc_finite_turnover import compute_result
from recursive_horizons.nsc_nested_qualities import finite_window

SCHEMA = "NSC-FINITE-CHANNEL-TURNOVER-v1"
PASS = "PASS_FINITE_CHANNEL_TURNOVER"
EXACT = {
    "G_cut_exact": Fraction("8986745/950450651136"),
    "total_inventory_exact": Fraction("7246127239/2477260800"),
    "sign_control_value_exact": Fraction(1, 4),
}
NUMERIC = {
    "probe_response_change_fro": 6.87914636e-5,
    "omitted_memory_error_fro": 0.00409574031,
    "omitted_drive_error_fro": 0.00220041986,
    "dropped_cross_covariance_error_fro": 6.46603026e-5,
    "half_identity_response_change_fro": 8.26498528e-4,
}


def as_fraction(value):
    if isinstance(value, str):
        return Fraction(value)
    if type(value) is int:
        return Fraction(value)
    raise AssertionError(f"exact metric must be a calculated fraction, not {value!r}")


def as_float(value):
    if isinstance(value, bool) or not isinstance(value, (int, float, str)):
        raise AssertionError(f"numeric metric must be a real number, not {value!r}")
    number = float(value)
    if not math.isfinite(number):
        raise AssertionError(f"numeric metric is not finite: {value!r}")
    return number


@pytest.fixture(scope="module")
def result():
    value = compute_result()
    json.dumps(value, sort_keys=True, allow_nan=False)
    return value


def test_exact_fractions_and_nonzero_sign(result):
    metrics = result["metrics"]
    for key, expected in EXACT.items():
        assert as_fraction(metrics[key]) == expected
    sign = as_fraction(metrics["sign_control_value_exact"])
    assert sign != 0 and sign > 0
    assert as_fraction(metrics["G_cut_exact"]) > 0
    assert as_fraction(metrics["total_inventory_exact"]) > 0


def test_numeric_references_and_far_probe(result):
    metrics = result["metrics"]
    for key, expected in NUMERIC.items():
        assert as_float(metrics[key]) == pytest.approx(expected, rel=1e-8, abs=1e-15)
    assert as_float(metrics["far_probe_response_change_fro"]) > 0


def test_verdict_follows_checks_and_assumptions_stay_limited(result):
    checks = result["checks"]
    assert result["schema"] == SCHEMA
    assert result["physical_local_gate"] == "OPEN"
    assert isinstance(result["inputs"], (dict, list))
    assert isinstance(checks, dict) and checks
    assert all(type(item) is bool for item in checks.values())
    passed = all(checks.values())
    assert (result["verdict"] == PASS) is passed
    assert passed
    assumptions = result["assumptions"]
    assert isinstance(assumptions, (dict, list, str))
    if isinstance(assumptions, dict):
        for key in (
            "eternity_claimed",
            "cosmological_identification",
            "new_physical_laws_claimed",
        ):
            if key in assumptions:
                assert assumptions[key] is False
        if "physical_local_gate" in assumptions:
            assert assumptions["physical_local_gate"] == "OPEN"


def test_frozen_six_dimensional_window_control():
    local = sp.Matrix([[1, sp.I / 5], [-sp.I / 5, 2]])
    link = sp.Matrix([
        [sp.Rational(1, 4), sp.I / 7],
        [sp.Rational(1, 9), sp.Rational(1, 6)],
    ])
    window = finite_window(local, link, sp.Rational(3, 2), 0, 3)
    assert window.shape == (6, 6)
    assert window == window.H
    stationary = sp.eye(6) / 2
    state = stationary - (window ** 3) / 2048
    assert sp.expand(state - state.H) == sp.zeros(6)
    assert sp.expand(state - stationary) != sp.zeros(6)
    assert sp.expand(window * stationary - stationary * window) == sp.zeros(6)


def _independent_turnover_fixture():
    import sympy as s
    from recursive_horizons.nsc_nested_qualities import finite_window
    H = s.Matrix([[1, s.I/5], [-s.I/5, 2]])
    B = s.Matrix([[s.Rational(1,4), s.I/7], [s.Rational(1,9), s.Rational(1,6)]])
    J = finite_window(H, B, s.Rational(3,2), 0, 3)
    C = s.eye(6)/2 - J**3/s.Integer(2048)
    low = s.Rational(3,2) - s.sqrt(29)/10
    high = s.Rational(3,2) + s.sqrt(29)/10
    qlow = (high*s.eye(2)-H)/(high-low)
    ps = []
    for region in range(3):
        for q in (qlow, s.eye(2)-qlow):
            p = s.zeros(6)
            p[2*region:2*region+2, 2*region:2*region+2] = q
            ps.append(p)
    return J, C, ps


def test_exact_row_bound_supports_admissibility():
    """Hermiticity and ||J||2 <= ||J||inf <= 8 imply spec(C) in [1/4,3/4]."""
    import sympy as s
    J, C, _ = _independent_turnover_fixture()
    assert (J-J.H).applyfunc(s.simplify) == s.zeros(6)
    bound = max(sum(s.Abs(J[i,j]) for j in range(6)) for i in range(6))
    assert bound == s.Rational(379,70)
    assert bound <= 8
    assert (C-C.H).applyfunc(s.simplify) == s.zeros(6)
    assert (J*C-C*J).applyfunc(s.simplify) == s.zeros(6)
    assert C == s.eye(6)/2 - J**3/s.Integer(2048)


def test_joint_representation_change_preserves_cut_activity():
    import sympy as s
    from recursive_horizons.nsc_finite_turnover import compute_result
    J, C, ps = _independent_turnover_fixture()
    R = s.eye(6)
    R[:2,:2] = s.Matrix([[s.Rational(3,5), -s.Rational(4,5)],
                         [s.Rational(4,5), s.Rational(3,5)]])
    assert R.H*R == s.eye(6)
    Jr = (R*J*R.H).applyfunc(s.expand)
    Cr = (R*C*R.H).applyfunc(s.expand)
    prs = [(R*p*R.H).applyfunc(s.expand) for p in ps]
    currents = [
        s.simplify(2*s.im(s.trace(prs[a]*Jr*prs[b]*Cr)))
        for a in range(2) for b in range(2,6)
    ]
    activity = s.simplify(sum(s.Abs(j) for j in currents))
    expected = s.Rational(8986745,950450651136)
    assert s.simplify(activity-expected) == 0
    assert expected > 0
    assert expected == s.sympify(compute_result()["metrics"]["G_cut_exact"])

