"""Exact source controls against vacuum, frame constraints and linear moments."""
import copy
import time

import numpy as np
import pytest

from recursive_horizons.nsc_atlas_characteristic_initialization import (
    Budget, constraint_samples, source_geometry, source_jet,
)
from recursive_horizons.nsc_atlas_nonlinear_bulk import bulk_operator
from recursive_horizons.nsc_atlas_nonlinear_source import ExactSource, source_amplitude_bound


@pytest.fixture(autouse=True)
def bounded_test_cpu():
    started = time.process_time()
    yield
    assert time.process_time()-started < 30


def build(alpha=None, **kwargs):
    budget = Budget(cpu_seconds=30)
    geometry = source_geometry(64, budget)
    threshold = source_amplitude_bound(geometry, budget=budget)["test_alpha"]
    return ExactSource(geometry, alpha=threshold if alpha is None else alpha,
                       budget=budget, **kwargs)


@pytest.fixture(scope="module")
def zero():
    return build(alpha=0.)


@pytest.fixture(scope="module")
def finite():
    return build()


def test_continuous_pulse_bounds_and_strict_admitted_amplitude(finite):
    b = finite.bounds
    assert finite.alpha == b["alpha_star"]/2
    assert finite.alpha*b["Lbar"]*b["length"] <= .25
    assert .4-finite.alpha*b["Lbar"]*b["length"]**2/2 > b["r_min"]
    assert 1-finite.alpha/(6*b["r_min"]**2) >= .75
    # Exact extrema of B*B_z and |B*B_zz-2*B_z^2| occur here.
    z1, z2 = 1/np.sqrt(15), 1/np.sqrt(21)
    assert 8*z1*(1-z1*z1)**7 == pytest.approx(8/np.sqrt(15)*(14/15)**7)
    assert 8*(1-z2*z2)**6*(1+9*z2*z2) == pytest.approx((80/7)*(20/21)**6)


def test_quiet_corner_and_initial_segment_are_exact_schwarzschild(finite):
    x = np.linspace(.4, finite.geometry["support_x"][0], 9)
    row = finite.sample(x)
    a = (.8/x-1)/4
    np.testing.assert_allclose(row["r"], x, rtol=0, atol=2e-15)
    np.testing.assert_array_equal(row["rJ_x"], np.ones_like(x))
    np.testing.assert_allclose(row["r_Ubar"], a, rtol=2e-10, atol=2e-12)
    np.testing.assert_allclose(row["log_C_Ubar"], -(1/.8+1/x)*a, atol=2e-11)
    for key in ("c", "c_x", "c_Ubar", "h", "h_x", "log_C_x"):
        np.testing.assert_array_equal(row[key], np.zeros_like(x))
    np.testing.assert_allclose(row["C"], 1., atol=0.)
    np.testing.assert_allclose(row["m"], .4, atol=2e-11)


def test_alpha_zero_geometry_mass_and_normalized_transverse_field(zero):
    x = np.linspace(.4, zero.geometry["x_last"], 17)
    row = zero.sample(x)
    np.testing.assert_array_equal(row["r"], x)
    np.testing.assert_array_equal(row["r_x"], np.ones_like(x))
    np.testing.assert_allclose(row["r_Ubar"], (.8/x-1)/4, atol=2e-11)
    np.testing.assert_allclose(row["m"], .4, atol=2e-11)
    nodes, weights = np.polynomial.legendre.leggauss(64)
    lower, upper_support = zero.geometry["support_x"]
    expected = []
    for endpoint in x:
        upper = min(endpoint, upper_support)
        if upper <= lower:
            expected.append(0.)
            continue
        q = (lower+upper)/2+(upper-lower)*nodes/2
        zero.budget.check(len(q))
        p, _, _ = source_jet(q, zero.geometry)
        # At alpha=0, K=.2/x^3. This independent one-dimensional Green
        # integral must retain c_Ubar although the metric is exactly vacuum.
        expected.append(-.2*(upper-lower)/2*float(weights @ (p/q**2)))
    np.testing.assert_allclose(row["c_Ubar"], expected, rtol=2e-9, atol=2e-11)
    assert abs(row["c_Ubar"][-1]) > .01


def test_zero_amplitude_scaled_response_matches_independent_improved_stress_moments(zero):
    moments = constraint_samples(zero.geometry, samples=17, order=64, budget=zero.budget)
    row = zero.sample(moments["x"])
    np.testing.assert_allclose(row["h"], moments["h"], rtol=3e-9, atol=2e-9)
    np.testing.assert_allclose(row["h_x"], moments["h_prime"], rtol=3e-9, atol=2e-9)
    assert row["h"][-1] < 0 and row["h_x"][-1] < 0


def test_finite_amplitude_chain_rule_constraint_and_mixed_data(finite):
    x = np.linspace(.4, finite.geometry["x_last"], 33)
    row = finite.sample(x)
    finite.budget.check(len(x))
    p, px, pxx = source_jet(x, finite.geometry)
    g, gx, gxx = x*p, p+x*px, 2*px+x*pxx
    y, yx, alpha = row["rJ"], row["rJ_x"], finite.alpha
    yxx = alpha*((g*gxx-2*gx*gx)/(6*y)+g*gx*yx/(3*y*y))
    # Independently differentiate rE^2=y^2-alpha*g^2/6, and the
    # canonical Psi=sqrt6*atanh(sqrt(alpha/6)*g/y).
    r = np.sqrt(y*y-alpha*g*g/6)
    rx = (y*yx-alpha*g*gx/6)/r
    rxx = (yx*yx+y*yxx-alpha*(gx*gx+g*gxx)/6-rx*rx)/r
    z, zx = np.sqrt(alpha)*g/y, np.sqrt(alpha)*(gx/y-g*yx/y**2)
    F = 1-z*z/6
    Lx = -z*zx/(3*F)
    psi_x = zx/F
    constraint = rxx-Lx*rx+r*psi_x*psi_x/2
    np.testing.assert_allclose(row["r"], r, rtol=2e-15)
    np.testing.assert_allclose(np.exp(row["log_C"]), F, rtol=2e-15)
    np.testing.assert_array_equal(row["log_C"], np.log1p(-alpha*g*g/(6*y*y)))
    np.testing.assert_allclose(row["log_C_x"], Lx, atol=3e-15)
    np.testing.assert_allclose(constraint, 0., atol=3e-14)
    np.testing.assert_allclose(row["source_constraint"], constraint, atol=3e-14)
    np.testing.assert_allclose(row["source_constraint_jordan"], 0., atol=3e-14)
    data = bulk_operator(dict(r=row["r"], C=row["C"], c=row["c"]),
        dict(r_U=row["r_Ubar"], r_V=row["r_x"], c_U=row["c_Ubar"], c_V=row["c_x"]),
        alpha=alpha, region="exterior")
    np.testing.assert_allclose(data["m"], row["m"], atol=2e-15)
    np.testing.assert_allclose(row["r_xx"]-row["log_C_x"]*row["r_x"],
                               data["constraint_V_rhs"], atol=3e-14)
    assert np.min(row["C"]) > .75 and np.min(row["rJ_x"]) > .75


def test_support_edges_and_post_pulse_affine_jordan_radius(finite):
    low, high = finite.geometry["support_x"]
    end = finite.geometry["x_last"]
    x = np.array([low, high, (high+end)/2, end])
    row = finite.sample(x)
    np.testing.assert_array_equal(row["c"], np.zeros(4))
    np.testing.assert_array_equal(row["c_x"], np.zeros(4))
    np.testing.assert_allclose(row["rJ_x"][1:], row["rJ_x"][1], atol=2e-15)
    np.testing.assert_allclose(row["rJ"][1:], row["rJ"][1]+row["rJ_x"][1]*(x[1:]-high), atol=2e-15)
    # The pulse's transverse scalar memory survives its longitudinal support.
    assert row["c_Ubar"][-1] != 0.
    np.testing.assert_allclose(row["c_Ubar"][1:], row["c_Ubar"][1], atol=2e-13)


def test_quiet_bottom_reference_null_geometry_and_corner_compatibility(finite):
    Ubar = np.linspace(0., finite.domain["Ubar_b"], 13)
    row = finite.quiet_bottom(Ubar)
    U, V, r = row["U"], finite.geometry["Va"], row["r"]
    # Independent implicit Kruskal relation and vacuum mass in the retained
    # trapped sector; both radius derivatives are positive on this branch.
    np.testing.assert_allclose((1-r/.8)*np.exp(r/.8), U*V, atol=2e-15)
    np.testing.assert_allclose(row["m"], .4, atol=2e-15)
    assert np.all(row["r_Ubar"] > 0) and np.all(row["r_x"] > 0)
    assert np.all(row["C"] > 0)
    assert abs(row["log_C_x"][-1]) > 1e-4
    corner = finite.sample(.4)
    np.testing.assert_allclose(np.exp(row["log_C"]), row["C"], rtol=2e-15)
    for key in ("r", "C", "log_C", "c", "r_Ubar", "r_x", "c_Ubar", "c_x", "log_C_Ubar", "log_C_x", "m"):
        assert row[key][0] == pytest.approx(corner[key], rel=2e-14, abs=2e-14)


def test_tighter_tolerance_refines_source_without_resetting_shared_budget():
    budget = Budget(cpu_seconds=30)
    geometry = source_geometry(64, budget)
    alpha = source_amplitude_bound(geometry, budget=budget)["test_alpha"]
    coarse = ExactSource(geometry, alpha=alpha, budget=budget)
    count = budget.evaluations
    fine = ExactSource(geometry, alpha=alpha, budget=budget, rtol=2e-12, atol=1e-14)
    x = np.linspace(.4, geometry["x_last"], 33)
    a, b = coarse.sample(x), fine.sample(x)
    keys = ("h", "h_x", "r", "C", "log_C", "r_x", "r_Ubar", "c_Ubar", "log_C_Ubar", "log_C_x", "m")
    for key in keys:
        np.testing.assert_allclose(a[key], b[key], rtol=2e-8, atol=2e-9)
    assert budget.evaluations > count and coarse.budget is fine.budget
    assert budget.evaluations <= 20000 and coarse.nfev > 0 and fine.nfev > 0
    # Emitted only when pytest capture is disabled; this is the existing
    # bounded control's receipt, not a separate trajectory calculation.
    print({"source_refinement": {"alpha_star": coarse.bounds["alpha_star"], "alpha": alpha,
        "maximum_differences": {key: float(np.max(abs(a[key]-b[key]))) for key in keys},
        "constraint_max": max(float(np.max(abs(a["source_constraint"]))),
                              float(np.max(abs(b["source_constraint"])))),
        "nfev": [coarse.nfev, fine.nfev], "cumulative_evaluations": budget.evaluations}})


@pytest.mark.parametrize("invalid", [True, None, "0.4", [], [True, .4],
                                     np.nan, np.inf, 1j, [[.4]], np.zeros(514)])
def test_invalid_dense_queries_rejected(zero, invalid):
    with pytest.raises(ValueError):
        zero.sample(invalid)
    with pytest.raises(ValueError):
        zero.quiet_bottom(invalid)


def test_no_source_or_bottom_extrapolation(zero):
    for x in (np.nextafter(.4, -np.inf), np.nextafter(zero.domain["x_b"], np.inf)):
        with pytest.raises(ValueError, match="closed domain"):
            zero.sample(x)
    for U in (-1e-15, np.nextafter(zero.domain["Ubar_b"], np.inf)):
        with pytest.raises(ValueError, match="closed domain"):
            zero.quiet_bottom(U)


@pytest.mark.parametrize("key,value", [("alpha", -.1), ("alpha", True), ("alpha", np.inf),
    ("rtol", 0.), ("rtol", 1e-15), ("rtol", .1), ("atol", 0.),
    ("max_step", 0.), ("max_step", 1.)])
def test_invalid_constructor_controls(zero, key, value):
    kwargs = dict(alpha=0.)
    kwargs[key] = value
    with pytest.raises(ValueError):
        ExactSource(zero.geometry, **kwargs)


def test_source_threshold_and_changed_geometry_rejected(zero):
    with pytest.raises(ValueError, match="alpha_star"):
        ExactSource(zero.geometry, alpha=zero.bounds["alpha_star"])
    geometry = copy.deepcopy(zero.geometry)
    geometry["support_V"][0] += .01
    with pytest.raises(ValueError, match="geometry"):
        ExactSource(geometry, alpha=0.)


@pytest.mark.parametrize("kwargs", [dict(maximum_evaluations=20001), dict(cpu_seconds=31),
    dict(maximum_evaluations=-1), dict(evaluations=-1), dict(evaluations=True),
    dict(cpu_seconds=np.nan), dict(started=float("inf"))])
def test_malformed_budgets_rejected(zero, kwargs):
    settings = dict(cpu_seconds=30)
    settings.update(kwargs)
    with pytest.raises(ValueError):
        ExactSource(zero.geometry, alpha=0., budget=Budget(**settings))


def test_count_and_cpu_exhaustion_stop_construction(zero):
    with pytest.raises(RuntimeError, match="exhausted"):
        ExactSource(zero.geometry, alpha=0., budget=Budget(maximum_evaluations=5, cpu_seconds=30))
    with pytest.raises(RuntimeError, match="exhausted"):
        ExactSource(zero.geometry, alpha=0., budget=Budget(cpu_seconds=.01, started=time.process_time()-1))


def test_geometry_is_copied_and_public_queries_are_counted():
    source = build(alpha=0.)
    original = copy.deepcopy(source.geometry)
    count = source.budget.evaluations
    source.sample([.4, source.domain["x_b"]])
    source.quiet_bottom([0., source.domain["Ubar_b"]])
    assert source.budget.evaluations == count+4
    assert source.geometry == original
