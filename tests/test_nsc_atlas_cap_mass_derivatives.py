"""Independent high-precision and polynomial controls; no retained field solve."""
import numpy as np
import pytest
import mpmath as mp

from control_budget import ControlBudget
from recursive_horizons.nsc_atlas_characteristic_background import CharacteristicProfile
from recursive_horizons import nsc_atlas_cap_mass_derivatives as model


CONTROLS = ControlBudget(maximum_evaluations=1_000_000, cpu_seconds=30)


@pytest.fixture
def budget():
    with CONTROLS.active() as result:
        yield result


@pytest.fixture(scope="module", autouse=True)
def report_allowance():
    yield
    if CONTROLS.budget is not None:
        print("mass derivative controls:", CONTROLS.budget.evaluations,
              "points;", CONTROLS.active_cpu_seconds, "active CPU seconds")


def cubic():
    z = np.array([-1., 0., 1.])
    return CharacteristicProfile(z, z**3, 3*z*z, 6*z)


def positive_cubic():
    z = np.array([-.5, 0., .5])
    p = np.maximum(z, 0)
    return CharacteristicProfile(z, p**3, 3*p*p, 6*p)


def cubic_oracle(T, r, name, small):
    """Closed cubic field before any spline, evaluated at high precision."""
    d, A = mp.atanh(r), 1-r*r
    f = ((T+d)**3-(T-d)**3)/r
    Pi = (3*(T+d)**2-3*(T-d)**2)/(r*A)
    if small:
        C = sum(mp.mpf(2*k)/(2*k+1)*r**(2*k-1) for k in range(1, 6))
        Phi = C*(6*T*T+2*d*d)+4*d*d/(r*A)
    else:
        Phi = ((3*(T+d)**2+3*(T-d)**2)/A-f)/r
    return dict(f=f, Pi=Pi, Phi=Phi)[name]


@pytest.mark.parametrize("r", [1e-8, .003, .01, .0101, .2])
def test_cubic_field_and_both_partials_against_independent_oracle(r, budget):
    T = .17
    row = model.background_partials(cubic(), T, r, budget=budget)
    with mp.workdps(70):
        tt, rr = mp.mpf(T), mp.mpf(r)
        for name in ("f", "Pi", "Phi"):
            fn = lambda t, x: cubic_oracle(t, x, name, r <= .01)
            expected = {name: fn(tt, rr), name+"_T": mp.diff(lambda t: fn(t, rr), tt),
                        name+"_r": mp.diff(lambda x: fn(tt, x), rr)}
            for key, value in expected.items():
                assert float(row[key]) == pytest.approx(float(value), rel=2e-10, abs=2e-9 if r > .01 else 1e-11), key
    budget.check(9)


@pytest.mark.parametrize("T", [-.001, 0., .001])
def test_third_moments_cross_a_c2_knot_without_a_jump_term(T, budget):
    d = .004
    z = -T/d
    # F''' is exactly 6 on z>z0 and zero on z<z0.
    expected = np.array([6*(1-z), 3*(1-z*z), 2*(1-z**3)])
    result = model.third_moments(positive_cubic(), T, d, budget=budget)
    np.testing.assert_allclose(result, expected, rtol=3e-14, atol=3e-14)


def test_c2_axis_limit_has_nonzero_radial_Pi_derivative(budget):
    row = model.background_partials(positive_cubic(), 0., [1e-8, .003], budget=budget)
    # Half of the F''' jump is 3; C2 alone does not give Pi_r(0+)=0.
    assert row["Pi_r"][0] == pytest.approx(3., abs=2e-12)
    with mp.workdps(70):
        for r, actual in zip((1e-8, .003), row["Pi_r"]):
            expected = mp.diff(lambda x: 3*mp.atanh(x)**2/(x*(1-x*x)), mp.mpf(r))
            assert actual == pytest.approx(float(expected), rel=3e-12)


def test_quiet_extension_exactly_zero_and_inputs_are_unchanged(budget):
    profile = positive_cubic()
    r = np.array([1e-8, .003, .1]); before = r.copy()
    coefficients = profile._spline.c.copy()
    result = model.background_partials(profile, -.8, r, budget=budget)
    assert all(np.array_equal(value, np.zeros(3)) for value in result.values())
    np.testing.assert_array_equal(r, before)
    np.testing.assert_array_equal(profile._spline.c, coefficients)


@pytest.mark.parametrize("r", [0., -1., 1., np.nan, np.ones(514)*.1])
def test_positive_bounded_radial_domain(r, budget):
    with pytest.raises(ValueError):
        model.background_partials(cubic(), 0., r, budget=budget)


def test_no_upper_extrapolation_and_no_unjustified_lower_extension(budget):
    with pytest.raises(ValueError, match="finite labels"):
        model.background_partials(cubic(), 1., .01, budget=budget)
    with pytest.raises(ValueError, match="finite labels"):
        model.background_partials(cubic(), -1., .01, budget=budget)
    with pytest.raises(ValueError, match="finite labels"):
        model.third_moments(cubic(), 1., .01, budget=budget)


@pytest.mark.parametrize("alpha", [0., .1, 100.])
def test_potential_partial_against_high_precision_actual_branch(alpha, budget):
    f = np.array([-.3, -1e-5, 0., 1e-5, .3])
    row = model.potential_partial(f, alpha, budget=budget)
    with mp.workdps(70):
        for value, expected_diff, actual in zip(f, row["Pdiff"], row["Pdiff_f"]):
            k = mp.sqrt(mp.mpf(alpha)/6)
            small = abs(np.sqrt(alpha/6)*value) < 1e-4
            def potential(x):
                q = k*x
                s = 1+q*q/6+q**4/120 if small else mp.sinh(q)/q
                return x*x*s*s*(mp.cosh(q)**2+1)/2
            xx = mp.mpf(value)
            assert expected_diff == pytest.approx(float(potential(xx)), rel=8e-15, abs=1e-28)
            assert actual == pytest.approx(float(mp.diff(potential, xx)), rel=8e-15, abs=1e-28)
    budget.check(10)


def test_symbolic_series_distinguishes_actual_potential_from_continuum_shortcut(budget):
    # Their production difference is below binary64 roundoff. This exact
    # rational coefficient control verifies the consequential derived identity.
    import sympy as sp
    q = sp.Symbol("q")
    s = 1+q*q/6+q**4/120
    actual = sp.diff(q*q*s*s*(sp.cosh(q)**2+1)/2, q)/q
    shortcut = 2*s*sp.cosh(q)**3
    assert sp.series(actual, q, 0, 7).removeO().coeff(q, 6) == sp.Rational(37, 45)
    assert sp.series(shortcut, q, 0, 7).removeO().coeff(q, 6) == sp.Rational(33, 40)
    budget.check(2)


@pytest.mark.parametrize("f", [0., .3])
def test_scalar_potential_shape_is_preserved(f, budget):
    scalar = model.potential_partial(f, .1, budget=budget)
    vector = model.potential_partial([f], .1, budget=budget)
    for key in scalar:
        assert np.shape(scalar[key]) == ()
        assert scalar[key] == vector[key][0]


def test_shared_budget_exhaustion_precedes_partial_calculation(budget):
    from recursive_horizons.nsc_atlas_characteristic_initialization import Budget
    child = Budget(maximum_evaluations=1, cpu_seconds=30)
    try:
        with pytest.raises(RuntimeError, match="budget"):
            model.background_partials(cubic(), .1, [.001, .2], budget=child)
    finally:
        budget.check(child.evaluations)
