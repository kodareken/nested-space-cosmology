"""Same-field basis equivalence and one fixed sixteen-mode reference."""
import numpy as np
import pytest

from control_budget import ControlBudget
from recursive_horizons import nsc_atlas_cap_polynomial_reference as model


CONTROLS = ControlBudget(maximum_evaluations=2500000,cpu_seconds=60.)


@pytest.fixture
def budget():
    with CONTROLS.active() as value:
        yield value


@pytest.fixture(scope='module',autouse=True)
def report():
    yield
    print('Legendre cap controls:',CONTROLS.budget.evaluations,'charged values;',
          CONTROLS.active_cpu_seconds,'active CPU seconds')


def inputs(count=16,**changes):
    j = np.arange(1,count+1,dtype=float)
    data = dict(Z=.5,chi_coefficients=.001*(-1.)**np.arange(count)/j**4,
                density_coefficients=.0015/j**4,alpha=.002,T_x=.7,
                radius_rate=.04,incoming_w_plus=.4,basis='legendre')
    data.update(changes)
    return data


def test_same_three_mode_field_is_independent_of_coordinate_basis(budget):
    # Columns are P1, P3, P5 expressed in the odd monomial basis.
    transform = np.array([[1.,-1.5,15/8],[0.,2.5,-70/8],[0.,0.,63/8]])
    a,d = np.array([.03,.02,-.005]),np.array([.01,.015,-.002])
    leg = model.polynomial_cap_reference(**inputs(3,chi_coefficients=a,density_coefficients=d),budget=budget)
    mono = model.polynomial_cap_reference(**inputs(3,chi_coefficients=transform@a,
                                                   density_coefficients=transform@d,basis='monomial'),budget=budget)
    for name in ['mass','lapse','mass_rate','physical_weak_flux','moving_wall_density','interface_balance_defect']:
        assert leg[name] == pytest.approx(mono[name],abs=2e-9)
    for name in ['chi','canonical_density']:
        np.testing.assert_allclose(leg['gradient'][name],transform.T@mono['gradient'][name],atol=2e-9,rtol=0)
        np.testing.assert_allclose(transform@leg['rates'][name],mono['rates'][name],atol=2e-9,rtol=0)
        np.testing.assert_allclose(leg['adjoint_gradient'][name],leg['gradient'][name],atol=2e-8,rtol=0)
    for name in ['f','Pi','Phi','delta']:
        assert leg['centre'][name] == pytest.approx(mono['centre'][name],abs=2e-12)
    assert leg['gradient']['Z'] == pytest.approx(mono['gradient']['Z'],abs=2e-9)
    assert leg['gram_condition'] == pytest.approx(11/3)


def test_sixteen_mode_geometric_limit_matches_exact_quadratic_mass(budget):
    data = inputs(alpha=0.)
    row = model.polynomial_cap_reference(**data,budget=budget)
    budget.check(4096)
    p = np.arange(1,32,2,dtype=float)
    Z,a,d = data['Z'],data['chi_coefficients'],data['density_coefficients']
    A0,R = 1/np.cosh(Z)**2,np.tanh(Z)
    M = np.diag(Z/(2*p+1))
    lower = np.minimum(p[:,None],p[None,:])
    K = lower*(lower+1)/(2*Z)
    S = K-A0/R*np.ones((16,16))
    expected = (a@S@a+d@M@d)/4
    gc,gk = S@a/2,M@d/2
    shape = (d@M@d/Z-a@K@a/Z+(2*A0+(A0/R)**2)*sum(a)**2)/4
    assert row['mass'] == pytest.approx(expected,abs=2e-9)
    np.testing.assert_allclose(row['gradient']['chi'],gc,rtol=0,atol=2e-9)
    np.testing.assert_allclose(row['gradient']['canonical_density'],gk,rtol=0,atol=2e-9)
    assert row['gradient']['Z'] == pytest.approx(shape,abs=2e-9)
    assert abs(row['wall']['momentum_discrepancy']) < 2e-9
    assert abs(row['balance_defect']) < 1e-10
    assert row['gram_condition'] == pytest.approx(21.)


def test_sixteen_mode_finite_alpha_uses_checked_gradients_and_shape(budget):
    row = model.polynomial_cap_reference(**inputs(),budget=budget)
    assert row['gradient_discrepancy'] < 2e-8
    assert abs(row['shape_discrepancy']) < 2e-8
    assert abs(row['balance_defect']) < 1e-10
    assert row['gram_condition'] == pytest.approx(21.)
    assert row['sampled_minimum_A'] > row['metric_certificate']['A_lower_float_down']
    assert row['adjoint_gradient']['Z'] == pytest.approx(row['gradient']['Z'],abs=2e-8)
    assert not row['numerical_error_enclosure']
    print('sixteen-mode manufactured:',{name:row[name] for name in
          ['mass','gradient_discrepancy','shape_discrepancy','balance_defect','gram_condition','coefficient_evaluations']})


@pytest.mark.parametrize('changes',[
    dict(basis='chebyshev'),dict(chi_coefficients=[0.]*17,density_coefficients=[0.]*17),
])
def test_unsupported_basis_or_count_stops_before_radial_solve(changes,budget,monkeypatch):
    def forbidden(*args,**kwargs):
        raise AssertionError('invalid basis started radial integration')
    monkeypatch.setattr(model.radial,'_radial_variations',forbidden)
    with pytest.raises(ValueError):
        model.polynomial_cap_reference(**inputs(**changes),budget=budget)
