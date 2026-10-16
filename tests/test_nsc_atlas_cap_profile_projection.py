"""Independent polynomial and piecewise-integral projection controls."""
import numpy as np
from numpy.polynomial import Polynomial, Legendre
import pytest
import sympy as s

from control_budget import ControlBudget
from recursive_horizons.nsc_atlas_characteristic_background import CharacteristicProfile
from recursive_horizons.nsc_atlas_cap_profile_projection import project_characteristic_profile


CONTROLS = ControlBudget(maximum_evaluations=2_000_000,cpu_seconds=30.)


@pytest.fixture
def budget():
    with CONTROLS.active() as value:
        yield value


@pytest.fixture(scope='module',autouse=True)
def receipt():
    yield
    print('profile projection controls:',CONTROLS.budget.evaluations,'charged values;',
          CONTROLS.active_cpu_seconds,'active CPU seconds')


def profile(poly,knots):
    knots = np.asarray(knots,dtype=float)
    return CharacteristicProfile(knots,poly(knots),poly.deriv()(knots),poly.deriv(2)(knots))


@pytest.mark.parametrize('count',[3,16])
def test_global_quintic_reproduction_includes_centre_and_wall(budget,count):
    F = Polynomial([.5,-.25,.125,.5,-.125,.25])
    T,Z = .125,.5
    row = project_characteristic_profile(profile(F,[-1.,0.,.25,1.]),T,Z,count=count,budget=budget)
    plus,minus = Polynomial([T,Z]),Polynomial([T,-Z])
    for key,field in [('chi',F(plus)-F(minus)),('k',F.deriv()(plus)-F.deriv()(minus))]:
        expected = field.convert(kind=Legendre).coef[1::2]
        expected = np.pad(expected,(0,count-len(expected)))
        np.testing.assert_allclose(row['coefficients'][key],expected,rtol=0,atol=2e-14)
    tolerance = 2e-13 if count == 3 else 1e-11
    assert max(row['residual_norms'].values()) < tolerance
    assert max(abs(value) for value in row['centre_error'].values()) < tolerance
    assert max(abs(value) for value in row['wall_error'].values()) < 10*tolerance


def test_piecewise_projection_and_norms_match_exact_integrals(budget):
    row = project_characteristic_profile(profile(Polynomial([0,0,0,1]),[0.,.25,1.]),
                                         .125,.5,count=3,budget=budget)
    x = s.symbols('x',real=True)
    T,Z = s.Rational(1,8),s.Rational(1,2)
    plus,minus = T+Z*x,T-Z*x
    pieces = dict(chi=[plus**3-minus**3,plus**3],k=[3*(plus**2-minus**2),3*plus**2])
    intervals = [(0,s.Rational(1,4)),(s.Rational(1,4),1)]
    def integral(expressions):
        return sum(s.integrate(s.expand(expression),(x,*interval)) for expression,interval in zip(expressions,intervals))
    for key,values in pieces.items():
        expected = [(2*p+1)*integral([v*s.legendre(p,x) for v in values]) for p in [1,3,5]]
        np.testing.assert_allclose(row['coefficients'][key],np.array(expected,dtype=float),rtol=0,atol=2e-14)
        approximation = sum(c*s.legendre(p,x) for c,p in zip(expected,[1,3,5]))
        errors = [approximation-v for v in values]
        squared = Z*integral([e**2 for e in errors])
        derivative_squared = integral([s.diff(e,x)**2 for e in errors])/Z
        assert row['residual_norms'][key] == pytest.approx(float(s.sqrt(squared)),abs=2e-13)
        assert row['residual_norms'][key+'_xi'] == pytest.approx(float(s.sqrt(derivative_squared)),abs=2e-13)
    budget.check(6000)
    # The error in the imposed incoming trace is exactly the two field-jet errors.
    R = row['radius']
    assert row['wall_error']['incoming'] == pytest.approx(
        (row['wall_error']['p']+row['wall_error']['h'])/R,abs=2e-14)


def test_quiet_extension_is_exact_and_has_no_spurious_centre(budget):
    row = project_characteristic_profile(profile(Polynomial([0,0,0,1]),[0.,1.]),
                                         -1.,.5,count=16,budget=budget)
    assert row['panels'] == 1
    for key in ('residual_norms','reference_norms','wall_error','centre_error'):
        assert all(value == 0. for value in row[key].values())


@pytest.mark.parametrize('changes',[dict(T=.8),dict(Z=0.),dict(count=17)])
def test_invalid_domain_and_count_stop(changes,budget):
    data = dict(T=.125,Z=.5,count=3); data.update(changes)
    with pytest.raises(ValueError):
        project_characteristic_profile(profile(Polynomial([0,0,0,1]),[0.,1.]),**data,budget=budget)


def test_projection_uses_shared_budget(budget):
    maximum = budget.maximum_evaluations
    budget.maximum_evaluations = budget.evaluations+1
    try:
        with pytest.raises(RuntimeError,match='budget exhausted'):
            project_characteristic_profile(profile(Polynomial([0,0,0,1]),[0.,1.]),.125,.5,budget=budget)
    finally:
        budget.maximum_evaluations = maximum
