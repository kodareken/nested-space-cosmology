"""Independent polynomial, reflected-axis and characteristic-transport controls."""
import time

import numpy as np
import pytest
from numpy.polynomial import Polynomial
from numpy.polynomial.legendre import leggauss

from recursive_horizons.nsc_atlas_characteristic_background import CharacteristicProfile
from recursive_horizons.nsc_atlas_characteristic_initialization import Budget


@pytest.fixture(scope="module")
def budget():
    result = Budget(maximum_evaluations=20000, cpu_seconds=30)
    yield result
    result.check()
    print("characteristic-background manufactured evaluations:", result.evaluations,
          "process CPU:", time.process_time()-result.started)


def polynomial_profile(polynomial, knots=(-2., -.3, .4, 2.)):
    knots = np.asarray(knots)
    return CharacteristicProfile(knots, polynomial(knots), polynomial.deriv()(knots),
                                 polynomial.deriv(2)(knots))


def radial_coefficient(r):
    # Independent convergent series for (r/(1-r*r)-atanh(r))/r**2.
    return sum((2*k/(2*k+1))*r**(2*k-1) for k in range(1, 20))


def test_quintic_reproduction_all_jets_and_copied_data(budget):
    p = Polynomial([.2, -.4, .3, .1, -.2, .05])
    labels = np.array([-1., -.2, .3, 1.])
    values = p(labels)
    profile = CharacteristicProfile(labels, values, p.deriv()(labels), p.deriv(2)(labels))
    labels[:] = 10
    values[:] = 100
    queries = np.r_[np.linspace(-1, 1, 31), -.2, .3]
    row = profile.jet(queries, budget=budget)
    for order, name in enumerate(("F", "F_prime", "F_second", "F_third")):
        np.testing.assert_allclose(row[name], p.deriv(order)(queries), rtol=2e-12, atol=2e-13)
    copy = profile.breakpoints
    copy[:] = 20
    np.testing.assert_array_equal(profile.breakpoints, [-1, -.2, .3, 1])


@pytest.mark.parametrize("T", [-.3, .2, .4])
def test_quadratic_reflected_profile_and_stable_small_radius(T, budget):
    profile = polynomial_profile(Polynomial([0., 0., 1.]))
    r = np.array([0., 1e-8, 1e-6, .00999, .01, .01001, .1, .4])
    row = profile.background(T, r, budget=budget)
    d, A = np.arctanh(r), 1-r*r
    ratio = np.ones_like(r)
    ratio[1:] = d[1:]/r[1:]
    expected = {"chi": 4*T*d, "p": 4*d/A, "q": 4*T/A,
                "f": 4*T*ratio, "Pi": 4*ratio/A,
                "Phi": 4*T*radial_coefficient(r)}
    for key, value in expected.items():
        if key == "Phi":
            # The direct branch divides spline endpoint roundoff by r^2.
            # Keep the centre test separate, so this does not relax r=1e-8.
            small = r <= .01
            np.testing.assert_allclose(row[key][small], value[small], rtol=3e-9, atol=3e-12)
            np.testing.assert_allclose(row[key][~small], value[~small], rtol=3e-9, atol=1e-10)
        else:
            np.testing.assert_allclose(row[key], value, rtol=3e-12, atol=3e-14)
    assert row["Phi"][1] == pytest.approx(8*T*1e-8/3, abs=1e-14)
    assert row["chi"][0] == row["p"][0] == row["Phi"][0] == 0


def test_c2_knot_centre_has_finite_jets_without_false_smoothness(budget):
    left, right, b = .2, .7, .3
    x = np.array([-1., 0., 1.])
    a = np.where(x < 0, left, right)
    profile = CharacteristicProfile(x, x+b*x*x+a*x**3,
                                    1+2*b*x+3*a*x*x, 2*b+6*a*x)
    r = np.array([0., 1e-8, 1e-5, .0009, .03])
    d, A = np.arctanh(r), 1-r*r
    ratio = np.ones_like(r)
    ratio[1:] = d[1:]/r[1:]
    row = profile.background(0., r, budget=budget)
    expected_pi = 4*b*ratio/A+3*(right-left)*d*ratio/A
    np.testing.assert_allclose(row["Pi"], expected_pi, rtol=2e-11, atol=1e-14)
    np.testing.assert_allclose(row["f"], 2*ratio+(left+right)*d*d*ratio,
                               rtol=1e-12, atol=1e-14)
    assert row["Phi"][0] == 0
    assert abs(row["Phi"][1]) < 5e-8
    # C2 suffices for finite first jets. A third-derivative jump leaves an O(r)
    # term in Pi, so a general profile does not promise an even Taylor series.
    assert (row["Pi"][2]-row["Pi"][0])/r[2] == pytest.approx(3*(right-left), abs=2e-5)
    jet = profile.jet(np.array([-1e-7, 0., 1e-7]), budget=budget)
    assert jet["F_third"][0] == pytest.approx(6*left, abs=1e-13)
    assert jet["F_third"][1] == pytest.approx(6*right, abs=1e-13)


def test_quiet_extension_and_nonzero_reflection(budget):
    profile = polynomial_profile(Polynomial([0., 0., 0., 1.]), (0., .5, 1.5))
    row = profile.jet([-10., -1., 0., .2], budget=budget)
    np.testing.assert_array_equal(row["F"][:3], np.zeros(3))
    np.testing.assert_array_equal(row["F_second"][:3], np.zeros(3))
    quiet = profile.background(-2., [0., .2], budget=budget)
    for key in ("chi", "p", "q", "f", "Pi", "Phi"):
        np.testing.assert_array_equal(quiet[key], np.zeros(2))
    r, T = .1, .5
    reflected = profile.background(T, r, budget=budget)
    d = np.arctanh(r)
    assert reflected["chi"] == pytest.approx((T+d)**3-(T-d)**3)
    assert abs(float(reflected["chi"])-(T+d)**3) > .05
    assert profile.domain == {"first_label": 0., "last_label": 1.5, "quiet_below": True}


def test_endpoint_integrals_match_independent_quadrature_with_reflection(budget):
    polynomial = Polynomial([.1, .4, -.2, .1, .02, -.03])
    profile = polynomial_profile(polynomial)
    T = .2
    faces = np.array([0., .07, .2, .35, .55])
    row = profile.integrals(T, faces, budget=budget)
    nodes, weights = leggauss(32)
    expected_p, expected_q = [], []
    for a, b in zip(faces[:-1], faces[1:]):
        r = (a+b)/2+(b-a)*nodes/2
        budget.check(len(r))
        d, A = np.arctanh(r), 1-r*r
        plus, minus = polynomial.deriv()(T+d), polynomial.deriv()(T-d)
        expected_p.append((b-a)/2*np.dot(weights, (plus-minus)/A))
        expected_q.append((b-a)/2*np.dot(weights, (plus+minus)/A))
    np.testing.assert_allclose(row["Ip"], expected_p, atol=3e-15)
    np.testing.assert_allclose(row["Iq"], expected_q, atol=3e-15)
    assert np.sum(row["Iq"]) == pytest.approx(row["chi_faces"][-1], abs=2e-16)
    assert row["chi_faces"][0] == 0


def test_analytic_reduced_pde_and_moving_integral_flux(budget):
    poly = Polynomial([.3, .1, -.4, .2, .05])
    profile = polynomial_profile(poly)
    T, R, Tx, Rx = .2, .55, .7, -.1
    y = np.linspace(0, 1, 13)
    r = R*y
    data = profile.background(T, r, budget=budget)
    d, A = np.arctanh(r), 1-r*r
    # Independently integrate analytic time derivatives and compare spatial
    # fluxes of the evaluator: p_T=(Aq)_r and q_T=(Ap)_r.
    nodes, weights = leggauss(16)
    integrated_pT, integrated_qT = [], []
    for a, b in zip(r[:-1], r[1:]):
        query = (a+b)/2+(b-a)*nodes/2
        budget.check(len(query))
        distance = np.arctanh(query)
        plus2 = poly.deriv(2)(T+distance)
        minus2 = poly.deriv(2)(T-distance)
        integrated_pT.append((b-a)/2*np.dot(weights, (plus2-minus2)/(1-query*query)))
        integrated_qT.append((b-a)/2*np.dot(weights, (plus2+minus2)/(1-query*query)))
    np.testing.assert_allclose(np.diff(A*data["q"]), integrated_pT, atol=2e-15)
    np.testing.assert_allclose(np.diff(A*data["p"]), integrated_qT, atol=2e-15)
    plus1, minus1 = poly.deriv()(T+d), poly.deriv()(T-d)
    endpoint_p_rate = plus1*(Tx+Rx*y/A)+minus1*(Tx-Rx*y/A)
    endpoint_q_rate = plus1*(Tx+Rx*y/A)-minus1*(Tx-Rx*y/A)
    fluxp = Tx*A*data["q"]+Rx*y*data["p"]
    fluxq = Tx*A*data["p"]+Rx*y*data["q"]
    np.testing.assert_allclose(np.diff(fluxp), np.diff(endpoint_p_rate), atol=2e-15)
    np.testing.assert_allclose(np.diff(fluxq), np.diff(endpoint_q_rate), atol=2e-15)
    assert fluxq[0] == 0
    # Actual mapped-cell integral derivative from an independent polynomial
    # endpoint chain, not a difference of sampled waveform values.
    assert sum(np.diff(fluxq)) == pytest.approx(endpoint_q_rate[-1], abs=5e-15)


def test_exact_incoming_characteristic_translation_has_no_outgoing_wave(budget):
    # Quiet below zero, F(q)=q^3 on [0,1]. All retained times precede centre
    # arrival; this is an exact characteristic propagation oracle, not an ODE.
    profile = polynomial_profile(Polynomial([0., 0., 0., 1.]), (0., .3, 1.))
    for T in (-.2, -.1, 0.):
        R = .5+.2*T
        faces = R*np.linspace(0, 1, 17)
        row = profile.integrals(T, faces, budget=budget)
        incoming = np.maximum(T+np.arctanh(faces), 0.)**3
        np.testing.assert_allclose(row["Iq"], np.diff(incoming), atol=2e-15)
        np.testing.assert_allclose(row["Ip"], row["Iq"], atol=2e-15)
        fields = profile.background(T, faces, budget=budget)
        np.testing.assert_allclose(fields["p"], fields["q"], atol=2e-15)


def test_mapped_knots_split_both_null_labels(budget):
    profile = polynomial_profile(Polynomial([0., 0., 0., 1.]), (0., .2, .5, 1.))
    cuts = profile.radial_cuts(.4, .5, budget=budget)
    np.testing.assert_allclose(cuts, [0., np.tanh(.1), np.tanh(.2), np.tanh(.4), .5], atol=1e-16)
    for a, b in zip(cuts[:-1], cuts[1:]):
        for sign in (-1, 1):
            endpoints = sorted((.4+sign*np.arctanh(a), .4+sign*np.arctanh(b)))
            assert not any(endpoints[0]+1e-14 < k < endpoints[1]-1e-14 for k in profile.breakpoints)


@pytest.mark.parametrize("jets", [([1., 1.], [0., 0.], [0., 0.]),
                                 ([0., 1.], [1., 0.], [0., 0.]),
                                 ([0., 1.], [0., 0.], [1., 0.])])
def test_lower_extension_requires_all_three_quiet_jets(jets, budget):
    profile = CharacteristicProfile([0., 1.], *jets)
    assert not profile.quiet_extension
    with pytest.raises(ValueError, match="finite labels"):
        profile.jet(-1e-12, budget=budget)


@pytest.mark.parametrize("labels,values,first,second", [
    ([0], [0], [0], [0]), ([0, 0], [0, 0], [0, 0], [0, 0]),
    ([0, 1], [0], [0, 0], [0, 0]), ([0, 1], [False, 1.], [0, 0], [0, 0]),
    ([0, 1], [0, 1], [0, np.nan], [0, 0]),
    ([[0, 1]], [0, 1], [0, 0], [0, 0]),
    (np.arange(514), np.zeros(514), np.zeros(514), np.zeros(514)),
])
def test_constructor_guards(labels, values, first, second):
    with pytest.raises(ValueError):
        CharacteristicProfile(labels, values, first, second)


def test_query_guards_and_budget_are_shared(budget):
    profile = polynomial_profile(Polynomial([0., 0., 0., 1.]), (0., 1.))
    for labels in (1.+1e-15, [0., True], np.zeros(514), np.nan):
        with pytest.raises(ValueError):
            profile.jet(labels, budget=budget)
    for T, r in ((.9, .5), (0., -.1), (0., 1.), (True, .1), ([0., .1], [0., .1, .2])):
        with pytest.raises(ValueError):
            profile.background(T, r, budget=budget)
    for faces in ([.1, .2], [0., .2, .1], [0., 1.], [0.]):
        with pytest.raises(ValueError):
            profile.integrals(0., faces, budget=budget)
    for T, r in ((0., 0.), (.9, .5), ([0., .1], .1)):
        with pytest.raises(ValueError):
            profile.radial_cuts(T, r, budget=budget)
    with pytest.raises(ValueError, match="shared cumulative"):
        profile.jet(0., budget=None)
    tiny = Budget(maximum_evaluations=2, cpu_seconds=30)
    profile.jet([0., .1], budget=tiny)
    with pytest.raises(RuntimeError, match="budget exhausted"):
        profile.jet(.2, budget=tiny)
    assert tiny.evaluations == 3
