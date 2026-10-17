"""Independent positive-kernel controls; no retained physical field solve."""
from fractions import Fraction as Q
from math import factorial

import pytest
import sympy as s

from control_budget import ControlBudget
from recursive_horizons.nsc_atlas_onset_endpoint_bound import (
    monotone_wall_enclosure, prefix_kernel_integral,
)

CONTROLS = ControlBudget(maximum_evaluations=20000, cpu_seconds=20.)


@pytest.fixture
def budget():
    with CONTROLS.active() as value:
        yield value


@pytest.fixture(scope='module', autouse=True)
def report():
    yield
    print('onset endpoint controls:', CONTROLS.budget.evaluations,
          'charged values;', CONTROLS.active_cpu_seconds, 'active CPU seconds')


def inputs():
    return dict(left=Q(1, 4), right=Q(1, 3), radius_lower=Q(2, 5),
                radius_upper=Q(1, 2), g=(Q(1, 12),)*2, g_x=(Q(1),)*2)


def test_kernel_change_of_variable_and_monotonicities(budget):
    budget.check(3000)
    a, x, z, r, C = s.symbols('a x z r C', positive=True)
    K = C*a/(4*r**3)
    radial_derivative = C*(a-x)/(4*x)
    independent = s.integrate(s.cancel(K/radial_derivative), (r, x, z))
    B = a*x/(2*(a-x))*(1/x**2-1/z**2)
    assert s.factor(independent-B) == 0
    assert s.factor(s.diff(B, z)-a*x/((a-x)*z**3)) == 0
    negative = -((a-x)**2/(a*x**2)+(a*a-z*z)/(a*z*z))
    assert s.factor(s.diff(B, x)-a*negative/(2*(a-x)**2)) == 0
    exact = independent.subs({a:s.Rational(4, 5), x:s.Rational(1, 3), z:s.Rational(1, 2)})
    assert prefix_kernel_integral(Q(1, 3), Q(1, 2)) == Q(int(exact.p), int(exact.q))


def test_independent_volterra_series(budget):
    """Constant positive kernel, linear quiet source; alternating exact tails."""
    budget.check(1000)
    row = monotone_wall_enclosure(**inputs(), budget=budget)
    t = Q(1, 12)
    M = sum(row['B_endpoint'].values())/2
    # The constant-kernel rectangular problem has chi(U,t)=sum
    # (-kU)^n t^(n+1)/(n!(n+1)!), obtained by repeated exact integrations.
    def partial(n, derivative=False):
        return sum((-M)**j*t**(j if derivative else j+1)
                   / (factorial(j)*factorial(j if derivative else j+1))
                   for j in range(n+1))
    for key, derivative in [('wall_value', False), ('wall_derivative', True)]:
        lo, hi = partial(9, derivative), partial(10, derivative)
        assert lo <= hi
        assert row[key]['lower'] <= lo <= hi <= row[key]['upper']
    assert 0 <= t-partial(10) <= t-partial(9) <= row['source_minus_wall']['upper']
    gap = row['source_derivative_minus_wall_derivative']
    assert gap['lower'] <= 1-partial(10, True) <= 1-partial(9, True) <= gap['upper']


def test_moving_clipped_prefix_has_zero_wall_transverse_derivative(budget):
    budget.check(2000)
    U, x, w0, speed, t = s.symbols('U x w0 speed t', positive=True)
    entry = (U-w0)/speed
    wall = w0+speed*t
    # Independently integrate f(v,t)=v+t with K=1, splitting where the
    # moving prefix first reaches U. These are wall-born columns, U>w0.
    integral = (s.integrate(wall**2/2+t*wall, (t, 0, entry))
                + s.integrate(U**2/2+t*U, (t, entry, x)))
    W = w0+speed*x
    assert s.factor(s.diff(integral, U).subs(U, W)) == 0
    along_wall = s.diff(integral.subs(U, W), x)
    assert s.factor(along_wall-(W**2/2+x*W)) == 0


def test_quiet_and_source_rescaling(budget):
    args = inputs()
    quiet = monotone_wall_enclosure(**{**args, 'g':(0, 0), 'g_x':(0, 0)}, budget=budget)
    output_keys = ('wall_value', 'wall_derivative', 'source_minus_wall',
                   'source_derivative_minus_wall_derivative')
    assert all(set(quiet[k].values()) == {0} for k in output_keys)
    base = monotone_wall_enclosure(**args, budget=budget)
    scaled = monotone_wall_enclosure(**{**args, 'g':(Q(1, 4),)*2, 'g_x':(3, 3)}, budget=budget)
    for key in output_keys:
        assert scaled[key] == {name:3*value for name, value in base[key].items()}
    assert scaled['q'] == base['q']


@pytest.mark.parametrize('change', [
    {'left':.25}, {'left':True}, {'g':(1, 0)}, {'g_x':(-1, 0)},
    {'radius_lower':Q(1, 3)}, {'radius_upper':Q(4, 5)},
    {'left':Q(1, 100), 'right':Q(1, 2), 'radius_lower':Q(3, 5), 'radius_upper':Q(7, 10)},
])
def test_invalid_input_or_failed_contraction_stops(change, budget):
    with pytest.raises(ValueError):
        monotone_wall_enclosure(**{**inputs(), **change}, budget=budget)


def test_budget_and_interval_inputs(budget):
    class Spent:
        def check(self, count=0):
            raise RuntimeError('spent')
    with pytest.raises(RuntimeError, match='spent'):
        monotone_wall_enclosure(**inputs(), budget=Spent())
    args = inputs()
    row = monotone_wall_enclosure(**{**args, 'g':(Q(1, 13), Q(1, 11)),
                                    'g_x':(Q(9, 10), Q(11, 10))}, budget=budget)
    exact = monotone_wall_enclosure(**args, budget=budget)
    for key in ('wall_value', 'wall_derivative'):
        assert row[key]['lower'] <= exact[key]['lower'] <= exact[key]['upper'] <= row[key]['upper']
