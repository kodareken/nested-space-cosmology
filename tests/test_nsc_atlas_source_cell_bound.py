"""Independent residual, frame-derivative and interval-input controls."""
from fractions import Fraction as Q

import mpmath as mp
import pytest
import sympy as s
from flint import arb, ctx

from control_budget import ControlBudget
from recursive_horizons.nsc_atlas_source_cell_bound import (
    _ball, _coordinate, _pulse_jets, frame_jacobian, radial_error_bounds,
    source_cell_error_bound,
)

CONTROLS = ControlBudget(maximum_evaluations=150_000, cpu_seconds=15.)
GEOMETRY = dict(U0=-1., Va=-.84, Vb=-.7, support_x=[.42, .6])


@pytest.fixture
def budget():
    with CONTROLS.active() as value:
        yield value


@pytest.fixture(scope='module', autouse=True)
def report():
    yield
    print('source cell controls:', CONTROLS.budget.evaluations, 'charged values;',
          CONTROLS.active_cpu_seconds, 'active CPU seconds')


def test_radial_rhs_partials_follow_from_original_constraint(budget):
    x, h, w, al, g, gx, gxx = s.symbols('x h w alpha g gx gxx')
    r, v = x+al*h, 1+al*w
    rhs = (g*gxx-2*gx**2)/(6*r)+g*gx*v/(3*r**2)
    budget.check(18000)
    assert s.cancel(s.diff(rhs, h)+al*(g*gxx-2*gx**2)/(6*r**2)+2*al*g*gx*v/(3*r**3)) == 0
    assert s.cancel(s.diff(rhs, w)-al*g*gx/(3*r**2)) == 0


def test_integral_bound_encloses_exact_forced_solution_and_zero_limit(budget):
    with ctx.workprec(160), mp.workdps(80):
        l, D, Lh, Lw = map(_ball, (Q(1, 10), Q(1, 5), 2, 3))
        theta, H, W = radial_error_bounds(l, D, Lh, Lw)
        r1, r2 = (3+mp.sqrt(17))/2, (3-mp.sqrt(17))/2
        t = mp.mpf(1)/10; drive = mp.mpf(1)/5
        exact_h = drive/(r1-r2)*((mp.exp(r1*t)-1)/r1-(mp.exp(r2*t)-1)/r2)
        exact_w = drive/(r1-r2)*(mp.exp(r1*t)-mp.exp(r2*t))
        assert theta < 1 and H > arb(str(exact_h)) and W > arb(str(exact_w))
        _, H0, W0 = radial_error_bounds(l, D, arb(0), arb(0))
        assert H0 >= l*l*D/2 and W0 >= l*D
        _, Hz, Wz = radial_error_bounds(l, arb(0), Lh, Lw)
        assert Hz == Wz == 0
    budget.check(4000)


def test_frame_derivatives_against_independent_high_precision_differentiation(budget):
    values = (Q(2, 5), Q(101, 100), Q(1, 100), Q(1, 10))
    with ctx.workprec(160), mp.workdps(85):
        al = mp.mpf(1)/50
        point = tuple(mp.mpf(v.numerator)/v.denominator for v in values)
        def frame(a, v, g, gx):
            F = 1-al*g*g/(6*a*a); r = mp.sqrt(a*a-al*g*g/6)
            return (r, mp.log(F), g, (a*v-al*g*gx/6)/r,
                    -al*(g*gx/a**2-g*g*v/a**3)/(3*F), gx)
        got = frame_jacobian(*(map(_ball, values)), _ball(Q(1, 50)))
        for index, name in enumerate(('r', 'log_C', 'c', 'r_x', 'log_C_x', 'c_x')):
            for variable in range(4):
                order = tuple(int(i == variable) for i in range(4))
                exact = mp.diff(lambda *args: frame(*args)[index], point, order)
                assert got[name][variable].contains(arb(str(exact))), (name, variable)
    budget.check(12000)


def test_compact_jet_hulls_enclose_independent_formulas_and_quiet_piece(budget):
    with ctx.workprec(160), mp.workdps(85):
        bounds = _pulse_jets(_ball(.4), _ball(.44), GEOMETRY)
        scale = mp.mpf(Q(.8).numerator)/Q(.8).denominator
        u, va, vb = (mp.mpf(Q(GEOMETRY[k]).numerator)/Q(GEOMETRY[k]).denominator for k in ('U0', 'Va', 'Vb'))
        def primitive(x):
            z = 4*((1-x/scale)*mp.exp(x/scale)/u-(va+vb)/2)/(vb-va)
            return (1-z*z)**4 if abs(z) < 1 else mp.mpf(0)
        for value in (Q(.4), Q(.42), Q(.44)):
            x = mp.mpf(value.numerator)/value.denominator
            for order in range(3):
                assert bounds[order].contains(arb(str(mp.diff(primitive, x, order))))
        assert all(v == 0 for v in _pulse_jets(_ball(.39), _ball(.4), GEOMETRY))
    budget.check(4000)


def test_original_input_certificate_covers_both_panels_and_mask(budget):
    result = source_cell_error_bound(GEOMETRY, alpha=.00001, left=.4, right=.44,
        power=(0, 0, Q(1, 10**7), Q(-2, 10**7)), budget=budget)
    assert result['panel_count'] == 2
    assert result['tube']['exact_continuation_proved']
    assert 0 < result['radial']['theta']['upper_float'] < 1
    assert result['radial']['residual']['upper_float'] > 0
    assert set(result['frame_error']) == {'r', 'log_C', 'c', 'r_x', 'log_C_x', 'c_x'}
    # A deliberately delayed front must not silently erase its nonzero sliver.
    assert result['front_mask_difference']['c']['upper_float'] > 0
    for field in ('c', 'c_x'):
        assert result['frame_error'][field]['upper_float'] >= result['front_mask_difference'][field]['upper_float']


def test_zero_coupling_has_exact_vacuum_metric_and_domain_failures(budget):
    result = source_cell_error_bound(GEOMETRY, alpha=0., left=.4, right=.44,
                                     power=(0, 0, 0, 0), budget=budget)
    for field in ('r', 'r_x', 'log_C', 'log_C_x'):
        assert result['frame_error'][field]['upper_float'] == 0
    with pytest.raises(ValueError, match='zero-left-jet'):
        source_cell_error_bound(GEOMETRY, alpha=0., left=.4, right=.44,
                                power=(1, 0, 0, 0), budget=budget)
    with pytest.raises(ValueError, match='continuation tube'):
        source_cell_error_bound(GEOMETRY, alpha=1., left=.4, right=.44,
                                power=(0, 0, 0, 0), budget=budget)
    with pytest.raises(ValueError, match='theta'):
        radial_error_bounds(arb(1), arb(1), arb(2), arb(0))
