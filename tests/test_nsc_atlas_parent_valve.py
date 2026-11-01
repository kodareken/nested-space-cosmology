"""Portable algebra controls; no history caches, new states or evolution."""
import math

import numpy as np
import pytest
import sympy as sp

from scripts.check_nsc_atlas_parent_valve import (
    conformal_sections,
    curvature_invariants,
    mcrit_relation,
)


def _tensor_invariants(A, B, C, D, F):
    """Contract an explicit orthonormal Riemann tensor, signature -+++."""
    tensor = np.zeros((4, 4, 4, 4))
    for (a, b, c, d), value in (
        ((0, 1, 0, 1), A), ((0, 2, 0, 2), B), ((0, 3, 0, 3), B),
        ((1, 2, 1, 2), C), ((1, 3, 1, 3), C), ((2, 3, 2, 3), D),
        ((0, 2, 1, 2), F), ((0, 3, 1, 3), F),
    ):
        for i, j, sij in ((a, b, 1), (b, a, -1)):
            for k, l, skl in ((c, d, 1), (d, c, -1)):
                tensor[i, j, k, l] = sij * skl * value
                tensor[k, l, i, j] = sij * skl * value
    signs = np.array([-1., 1., 1., 1.])
    metric = np.diag(signs)
    ricci = np.einsum('a,abad->bd', signs, tensor)
    scalar = float(np.einsum('a,aa', signs, ricci))
    weyl = np.empty_like(tensor)
    for a, b, c, d in np.ndindex(tensor.shape):
        weyl[a, b, c, d] = tensor[a, b, c, d] - .5 * (
            metric[a, c] * ricci[b, d] - metric[a, d] * ricci[b, c]
            - metric[b, c] * ricci[a, d] + metric[b, d] * ricci[a, c]
        ) + scalar / 6 * (
            metric[a, c] * metric[b, d] - metric[a, d] * metric[b, c]
        )
    def norm4(value):
        return float(np.einsum('a,b,c,d,abcd,abcd', signs, signs, signs, signs,
                               value, value))
    return dict(Ricci_scalar=scalar,
        Ricci_squared=float(np.einsum('a,b,ab,ab', signs, signs, ricci, ricci)),
        Kretschmann=norm4(tensor), Weyl_squared=norm4(weyl))


@pytest.mark.parametrize('mass,radius', [(1., 2.), (.7, .3)])
def test_schwarzschild_curvature_is_vacuum_and_has_expected_weyl(mass, radius):
    section = mass / radius**3
    row = curvature_invariants(-2*section, section, -section, 2*section)
    assert row['Ricci_scalar'] == pytest.approx(0., abs=1e-13)
    assert row['Ricci_squared'] == pytest.approx(0., abs=1e-13)
    for key in ('Kretschmann', 'Weyl_squared'):
        assert row[key] == pytest.approx(48*mass**2/radius**6, rel=2e-14)


@pytest.mark.parametrize('H', [.3, 1.7])
def test_de_sitter_constant_curvature_control(H):
    row = curvature_invariants(-H**2, -H**2, H**2, H**2)
    assert row['Ricci_scalar'] == pytest.approx(12*H**2, rel=2e-14)
    assert row['Ricci_squared'] == pytest.approx(36*H**4, rel=2e-14)
    assert row['Kretschmann'] == pytest.approx(24*H**4, rel=2e-14)
    assert row['Weyl_squared'] == pytest.approx(0., abs=1e-13)


def test_general_curvature_matches_explicit_tensor_with_mixed_section():
    sections = (.7, -.4, .9, 1.1, .2)
    independent = _tensor_invariants(*sections)
    actual = curvature_invariants(*sections)
    for key, value in independent.items():
        assert actual[key] == pytest.approx(value, rel=2e-14, abs=2e-14)
    assert actual['Kretschmann'] - 2*actual['Ricci_squared'] + (
        actual['Ricci_scalar']**2/3
    ) == pytest.approx(actual['Weyl_squared'], rel=2e-14, abs=2e-14)
    assert actual['decomposition_residual'] == pytest.approx(0., abs=2e-14)


@pytest.mark.parametrize('power', [1., 1.5])
def test_conformal_sections_from_direct_rescaled_warped_metric(power):
    # A locally normalized conformal base chart fixes all physical curvature
    # sections at its origin, without invoking a source or Einstein equation.
    t, x = sp.symbols('t x', real=True)
    radius = sp.Rational(3, 5)
    velocity = -sp.Rational(7, 4)
    A, B, C = -sp.Rational(5, 3), sp.Rational(4, 7), -sp.Rational(2, 5)
    p = sp.Rational(str(power))
    a = -A*t*t/2
    r = radius + velocity*t - B*radius*t*t/2 - C*radius*x*x/2
    a_new = a - p*sp.log(r)
    r_new = r**(1-p)
    at_origin = {t: 0, x: 0}
    e = sp.exp(-2*a_new)
    direct = (
        -e*(sp.diff(a_new, t, 2)-sp.diff(a_new, x, 2)),
        -e*(sp.diff(r_new, t, 2)-sp.diff(a_new, t)*sp.diff(r_new, t)
            - sp.diff(a_new, x)*sp.diff(r_new, x))/r_new,
        -e*(sp.diff(r_new, x, 2)-sp.diff(a_new, t)*sp.diff(r_new, t)
            - sp.diff(a_new, x)*sp.diff(r_new, x))/r_new,
        (1+e*(sp.diff(r_new, t)**2-sp.diff(r_new, x)**2))/r_new**2,
    )
    D = (1+velocity**2)/radius**2
    old = dict(zip(('A', 'B', 'C', 'D'), map(float, (A, B, C, D))))
    row = conformal_sections(float(radius), float(velocity), old, power)
    expected = [float(sp.simplify(value.subs(at_origin))) for value in direct]
    np.testing.assert_allclose(list(row['sectional_components'].values()), expected,
                               rtol=3e-14, atol=3e-14)
    new_r = row['rescaled_areal_radius']
    new_u = row['rescaled_normal_radius_rate']
    assert row['sectional_components']['D'] == pytest.approx(
        (1+new_u**2)/new_r**2, rel=3e-14)
    old_weyl = curvature_invariants(**old)['Weyl_squared']
    assert row['invariants']['Weyl_squared'] == pytest.approx(
        float(radius)**(4*power)*old_weyl, rel=3e-14)
    if power == 1.:
        assert new_r == 1.
        assert row['sectional_components']['B'] == 0.
        assert row['sectional_components']['C'] == 0.
        assert row['sectional_components']['D'] == pytest.approx(1.)


def test_conformal_vacuum_dominant_cut_has_negative_angular_nec():
    # This is an algebraic metric diagnostic, not a prepared matter state or
    # a claim that finite conformal curvature supplies a gluing theorem.
    mass, radius = .7, .15
    U = -math.sqrt(2*mass/radius-1)
    f = mass/radius**3
    row = conformal_sections(radius, U, dict(A=-2*f, B=f, C=-f, D=2*f), 1.5)
    tensor = _tensor_invariants(*row['sectional_components'].values(), 0.)
    assert row['invariants']['Kretschmann'] == pytest.approx(tensor['Kretschmann'], rel=2e-14)
    a, b, c, d = row['sectional_components'].values()
    # Contract Einstein with orthonormal angular/radial null vectors.
    einstein_00, einstein_11, einstein_22 = 2*c+d, 2*b-d, a+b-c
    actual = row['effective_Einstein_stress']
    assert actual['angular_NEC'] == pytest.approx(einstein_00+einstein_22)
    assert actual['radial_NEC'] == pytest.approx(einstein_00+einstein_11)
    assert actual['angular_NEC'] < 0 < actual['radial_NEC']


def test_conformal_reflection_identity_is_required():
    with pytest.raises(ValueError, match='reflection radius-gradient identity'):
        conformal_sections(1., -1., dict(A=1., B=1., C=1., D=1.), 1.5)


@pytest.fixture(scope='module')
def loading_record():
    return mcrit_relation()


def test_loading_fold_and_nonidentity_by_exact_elimination(loading_record):
    x, s = sp.symbols('x s', positive=True)
    V = 7*sp.exp(-12*s)-sp.Rational(15, 2)*sp.exp(-14*s)+sp.Rational(25, 32)*sp.exp(-24*s)
    # Stationarity of V+M exp(-6s), at fixed oscillator-action loading M.
    M = sp.exp(6*s)*sp.diff(V, s)/6
    fold = sp.expand(4*sp.diff(M, s).subs(s, -sp.log(x)/2)/x**3)
    P = sp.Poly(fold, x)
    assert P.as_expr() == 225*x**6-560*x+336
    record_P = sp.sympify(loading_record['exact']['fold_polynomial'], locals={'x': x})
    assert sp.expand(record_P-P.as_expr()) == 0
    M_x = sp.expand(M.subs(s, -sp.log(x)/2))
    reduced = x**3*(175*x-168)/18
    assert sp.rem(M_x-reduced, P.as_expr(), x) == 0
    equality = sp.together(reduced-sp.Rational(3, 8)).as_numer_denom()[0]
    Q = sp.Poly(equality, x).primitive()[1]
    assert sp.gcd(P, Q).degree() == 0
    assert Q.as_expr() == 700*x**4-672*x**3-27
    record_Q = sp.sympify(loading_record['exact']['equality_polynomial'], locals={'x': x})
    assert sp.expand(record_Q-Q.as_expr()) == 0
    lo, hi = map(sp.Rational, loading_record['exact']['isolating_interval'])
    assert sp.Rational(99, 100) < lo < hi < 1
    assert P.eval(lo)*P.eval(hi) < 0
    assert P.count_roots(lo, hi) == 1
    xi = loading_record['fold_x']
    assert xi == pytest.approx(float((lo+hi)/2), abs=2e-15)
    critical = float(M_x.subs(x, xi))
    assert loading_record['Mcrit'] == pytest.approx(critical, abs=2e-14)
    assert loading_record['ratio'] == pytest.approx(critical/(3/8), abs=6e-14)
    assert 0 < loading_record['percent_above_E0'] < .022
    assert loading_record['actual_density_at_fold'] == pytest.approx(
        critical*math.exp(-6*loading_record['fold_s']), abs=2e-14)


def test_matching_fixes_M_at_zero_but_does_not_fix_fold_equality(loading_record):
    s, m = sp.symbols('s m', real=True)
    A, B, T = 3*m+sp.Rational(47, 8), 3*m+sp.Rational(51, 8), sp.Rational(25, 32)
    V = A*sp.exp(-12*s)-B*sp.exp(-14*s)+T*sp.exp(-24*s)
    M = sp.exp(6*s)*sp.diff(V, s)/6
    assert sp.simplify(M.subs(s, 0)-m) == 0
    derivative = sp.simplify(sp.diff(M, s).subs(s, 0))
    second = sp.simplify(sp.diff(M, s, 2).subs(s, 0))
    assert derivative == sp.Rational(31, 4)-20*m
    assert second == 232*m-sp.Rational(967, 2)
    exact_fold_m = sp.solve(derivative, m)[0]
    assert exact_fold_m == sp.Rational(31, 80)
    assert second.subs(m, exact_fold_m) < 0
    # An allowed symbolic member puts the maximum at the matched point
    # exactly; the fiducial m=3/8 instead has a nonzero slope there.
    assert derivative.subs(m, sp.Rational(3, 8)) == sp.Rational(1, 4)
    assert sp.Rational(loading_record['exact']['matched_fold_m']) == exact_fold_m
    assert sp.sympify(loading_record['exact']['matched_M_at0'], locals={'m': m}) == m
    E0 = (1-sp.Rational(1, 4))*(sp.Rational(3, 2)-sp.Rational(1, 2))**2/2
    assert E0 == sp.Rational(3, 8)
    assert loading_record['E0'] == float(E0)
