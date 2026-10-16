"""Supplied wall chain, computational curvature and finite-domain controls."""
import time

import numpy as np
import pytest

from recursive_horizons.nsc_atlas_background_profile import (
    CURVATURE_METHOD, build_background_profile,
)
from recursive_horizons.nsc_atlas_characteristic_initialization import Budget


@pytest.fixture(scope="module")
def budget():
    result = Budget(maximum_evaluations=20000, cpu_seconds=30)
    yield result
    result.check()
    print("background-profile manufactured evaluations:", result.evaluations,
          "process CPU:", time.process_time()-result.started)


def cubic_data(n=9):
    z = np.linspace(0., .5, n)
    labels = 1+z
    radius = np.full(n, .75)
    rates = 1+labels*labels
    walls = {"T": labels-np.arctanh(radius), "r": radius, "c": z**3}
    wall_rates = {"T": rates, "r": np.zeros(n), "c": 3*z*z*rates}
    return walls, wall_rates, labels


def build(walls, rates, budget, alpha=0.):
    return build_background_profile(walls, rates, alpha=alpha, budget=budget)


def test_known_cubic_chain_and_same_interpolant_jets(budget):
    walls, rates, labels = cubic_data()
    result = build(walls, rates, budget)
    z = labels-1
    np.testing.assert_array_equal(result["labels"], labels)
    np.testing.assert_array_equal(result["values"], z**3)
    np.testing.assert_allclose(result["first"], 3*z*z, atol=1e-16)
    np.testing.assert_allclose(result["second"], 6*z, atol=1e-14)
    np.testing.assert_array_equal(result["arrival_rate"], rates["T"])
    query = np.linspace(1., 1.5, 19)
    jet = result["profile"].jet(query, budget=budget)
    for name, expected in (("F", (query-1)**3), ("F_prime", 3*(query-1)**2),
                           ("F_second", 6*(query-1)), ("F_third", np.full(19, 6.))):
        np.testing.assert_allclose(jet[name], expected, rtol=1e-11, atol=1e-12)
    assert result["curvature_method"] == CURVATURE_METHOD
    assert "not physical second-jet evidence" in result["scope"]
    diag = result["diagnostics"]
    assert diag["no_return_margin"] == pytest.approx(np.arctanh(.75)-.5)
    assert diag["quiet_anchor"] == {"value": 0., "first": 0., "second": 0.}
    assert diag["max_abs_material_chain_defect"] < 1e-15
    assert diag["max_abs_knot_value_defect"] < 1e-15
    assert diag["max_abs_knot_first_defect"] < 1e-14
    assert diag["max_abs_knot_second_defect"] < 1e-12
    assert result["charged_evaluations"] == 10*len(labels)


def test_nonconstant_radius_material_chain_uses_both_rates(budget):
    labels = 1+np.linspace(0., .5, 13)
    radii = np.tanh(.8+.1*labels)
    label_rate = .7+labels*labels
    z = np.maximum(labels-labels[3], 0.)
    walls = {"T": labels-np.arctanh(radii), "r": radii, "c": z**3}
    rates = {"T": .9*label_rate, "r": .1*(1-radii*radii)*label_rate,
             "c": 3*z*z*label_rate}
    result = build(walls, rates, budget)
    np.testing.assert_allclose(result["arrival_rate"], label_rate, atol=5e-16)
    np.testing.assert_allclose(result["first"], 3*z*z, atol=2e-16)
    # Away from the front the sampled first jet is quadratic, so the declared
    # nonuniform second-order curvature recovers its independent linear jet.
    np.testing.assert_allclose(result["second"][5:], 6*z[5:], atol=2e-14)
    changed_parameter = build(walls, {key: 2*value for key, value in rates.items()}, budget)
    np.testing.assert_array_equal(changed_parameter["first"], result["first"])
    np.testing.assert_array_equal(changed_parameter["second"], result["second"])


def test_first_jets_are_supplied_not_waveform_differences(budget):
    walls, rates, _ = cubic_data()
    baseline = build(walls, rates, budget)
    changed_values = {key: value.copy() for key, value in walls.items()}
    changed_values["c"][4] += .1
    changed = build(changed_values, rates, budget)
    np.testing.assert_array_equal(changed["first"], baseline["first"])
    np.testing.assert_array_equal(changed["second"], baseline["second"])
    assert changed["values"][4]-baseline["values"][4] == pytest.approx(.1)
    # This deliberately inconsistent manufacture demonstrates why provenance
    # of the supplied same-stage rates remains a caller obligation.
    walls["c"][:] = 20
    rates["c"][:] = 30
    assert baseline["walls"]["c"][0] == baseline["wall_rates"]["c"][0] == 0


@pytest.mark.parametrize("n", [65, 129, 513])
def test_caller_selected_rows_and_finite_coverage(n, budget):
    walls, rates, labels = cubic_data(n)
    result = build(walls, rates, budget)
    assert result["profile"].breakpoints.size == n
    assert result["label_range"] == [labels[0], labels[-1]]
    assert result["profile"].quiet_extension
    quiet = result["profile"].jet(np.array([-1., .5]), budget=budget)
    for value in quiet.values():
        np.testing.assert_array_equal(value, np.zeros(2))
    with pytest.raises(ValueError, match="finite labels"):
        result["profile"].jet(np.nextafter(labels[-1], np.inf), budget=budget)


@pytest.mark.parametrize("alpha", [True, .1, -1e-10, np.inf, np.nan, "0", [0.]])
def test_only_zero_case_is_accepted(alpha, budget):
    walls, rates, _ = cubic_data()
    with pytest.raises(ValueError, match="alpha=0"):
        build(walls, rates, budget, alpha)


@pytest.mark.parametrize("which", ["value", "first", "second"])
def test_each_quiet_anchor_condition_is_required(which, budget):
    walls, rates, _ = cubic_data()
    if which == "value":
        walls["c"][0] = .01
    elif which == "first":
        rates["c"][0] = .01
    else:
        rates["c"][1] += .01
    with pytest.raises(ValueError, match="exactly quiet"):
        build(walls, rates, budget)


def test_returning_or_nonmonotone_or_stationary_input_rejected(budget):
    walls, rates, _ = cubic_data(5)
    # A wider binary grid retains exact quiet cubic curvature, but reaches
    # the first arrival. It cannot use F=c as a no-return wall trace.
    z = np.arange(5)/4
    walls["T"] = 1+z-np.arctanh(.75)
    walls["c"] = z**3
    rates["T"][:] = 1
    rates["c"] = 3*z*z
    with pytest.raises(ValueError, match="first central arrival"):
        build(walls, rates, budget)
    walls, rates, _ = cubic_data()
    walls["T"][2] = walls["T"][1]
    with pytest.raises(ValueError, match="labels must increase"):
        build(walls, rates, budget)
    walls, rates, _ = cubic_data()
    rates["T"][2] = 0
    with pytest.raises(ValueError, match="rates must be positive"):
        build(walls, rates, budget)


@pytest.mark.parametrize("kind", ["missing", "scalar", "short", "large", "shape",
                                  "boolean", "mixed_boolean", "nan", "radius_zero", "radius_one"])
def test_input_guards(kind, budget):
    walls, rates, _ = cubic_data()
    if kind == "missing":
        del walls["r"]
    elif kind == "scalar":
        rates["r"] = 0.
    elif kind == "short":
        walls = {key: value[:2] for key, value in walls.items()}
    elif kind == "large":
        walls = {key: np.zeros(514) for key in walls}
    elif kind == "shape":
        rates["r"] = rates["r"][:-1]
    elif kind == "boolean":
        walls["T"] = np.zeros(9, dtype=bool)
    elif kind == "mixed_boolean":
        walls["T"] = [False, *walls["T"][1:]]
    elif kind == "nan":
        rates["r"][4] = np.nan
    elif kind == "radius_zero":
        walls["r"][1] = 0.
    else:
        walls["r"][1] = 1.
    with pytest.raises(ValueError):
        build(walls, rates, budget)


def test_cumulative_budget_and_no_internal_reset(budget):
    walls, rates, _ = cubic_data()
    with pytest.raises(ValueError, match="shared cumulative"):
        build(walls, rates, None)
    limited = Budget(maximum_evaluations=179, cpu_seconds=30)
    build(walls, rates, limited)
    assert limited.evaluations == 90
    with pytest.raises(RuntimeError, match="budget exhausted"):
        build(walls, rates, limited)
    assert limited.evaluations == 180
