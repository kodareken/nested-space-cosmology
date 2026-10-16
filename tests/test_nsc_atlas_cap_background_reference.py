"""Chain-rule and incoming-limit checks for the same mass with a background."""
import numpy as np
from numpy.polynomial import Polynomial, Legendre
import pytest

from control_budget import ControlBudget
from recursive_horizons.nsc_atlas_cap_background_jets import RetainedCapBackground
from recursive_horizons.nsc_atlas_cap_polynomial_basis import odd_legendre_operators
from recursive_horizons.nsc_atlas_cap_polynomial_reference import polynomial_cap_reference


CONTROLS = ControlBudget(maximum_evaluations=4_000_000,cpu_seconds=60.)


@pytest.fixture
def budget():
    with CONTROLS.active() as value:
        yield value


@pytest.fixture(scope='module',autouse=True)
def receipt():
    yield
    print('background reference controls:',CONTROLS.budget.evaluations,'charged values;',
          CONTROLS.active_cpu_seconds,'active CPU seconds')


def background(poly):
    knots=np.array([-1.,0.,.5,1.])
    return RetainedCapBackground(knots,poly(knots),poly.deriv()(knots),poly.deriv(2)(knots))


def test_same_polynomial_field_checks_coefficient_time_and_shape_chains(budget):
    F=Polynomial([0,0,0,.01])
    bg=background(F)
    T,Z=.125,.5
    c,d=np.array([.001,-.002]),np.array([.003,.001])
    plus,minus=Polynomial([T,Z]),Polynomial([T,-Z])
    def coeff(poly):
        v=poly.convert(kind=Legendre).coef[1::2]
        return np.pad(v,(0,2-len(v)))
    cb,kb=coeff(F(plus)-F(minus)),coeff(F.deriv()(plus)-F.deriv()(minus))
    ab,bb=kb,coeff(F.deriv(2)(plus)-F.deriv(2)(minus))
    options=dict(alpha=.002,T_x=.7,radius_rate=.04,incoming_w_plus=.2,basis='legendre')
    affine=polynomial_cap_reference(Z,c,d,background=bg,background_time=T,budget=budget,**options)
    full=polynomial_cap_reference(Z,c+cb,d+kb,budget=budget,**options)
    E=odd_legendre_operators(2,Z,budget=budget)['dilation']
    for key in ('mass','lapse','physical_weak_flux','mass_rate','interface_balance_defect'):
        assert affine[key] == pytest.approx(full[key],abs=2e-10)
    for key in ('chi','canonical_density'):
        np.testing.assert_allclose(affine['gradient'][key],full['gradient'][key],rtol=0,atol=2e-10)
    expected_time=full['gradient']['chi']@ab+full['gradient']['canonical_density']@bb
    expected_shape=full['gradient']['Z']+(full['gradient']['chi']@E@cb+full['gradient']['canonical_density']@E@kb)/Z
    assert affine['gradient']['T']==pytest.approx(expected_time,abs=2e-10)
    assert affine['gradient']['Z']==pytest.approx(expected_shape,abs=2e-10)
    for key,base,rate in [('chi',cb,ab),('canonical_density',kb,bb)]:
        converted=affine['rates'][key]+options['T_x']*rate+affine['Z_rate']*E@base/Z
        np.testing.assert_allclose(converted,full['rates'][key],rtol=0,atol=2e-9)
    assert abs(affine['explicit_time_remainder'])<2e-10
    assert abs(affine['remainder_identity_defect'])<2e-12
    assert affine['adjoint_gradient']['T']==pytest.approx(affine['gradient']['T'],abs=2e-8)


def test_exact_quartic_background_retains_the_derived_incompatible_remainder(budget):
    Z=.5
    options=dict(alpha=0.,T_x=1.,radius_rate=0.,basis='legendre',
                 background=background(Polynomial([0,0,0,0,1.])),background_time=0.)
    row=polynomial_cap_reference(Z,[0.],[0.],incoming_w_plus=0.,budget=budget,**options)
    assert row['mass_rate']==pytest.approx(-96*Z**6/5,abs=2e-10)
    assert row['physical_weak_flux']==pytest.approx(-32*Z**6,abs=2e-10)
    assert row['explicit_time_remainder']==pytest.approx(64*Z**6/5,abs=2e-10)
    assert row['balance_defect']==pytest.approx(64*Z**6/5,abs=2e-10)
    assert abs(row['remainder_identity_defect'])<2e-12
    A0,R=1/np.cosh(Z)**2,np.tanh(Z)
    matched=polynomial_cap_reference(Z,[0.],[0.],incoming_w_plus=8*Z**3/(A0*R),budget=budget,**options)
    assert abs(matched['explicit_time_remainder'])<2e-10
    assert abs(matched['rates']['chi'][0])<2e-10
    assert abs(matched['rates']['canonical_density'][0])<2e-10


def test_quiet_background_is_zero_without_a_false_positive_mass_guard(budget):
    poly=Polynomial([0,0,0,1.])
    knots=np.array([0.,.5,1.])
    bg=RetainedCapBackground(knots,poly(knots),poly.deriv()(knots),poly.deriv(2)(knots))
    row=polynomial_cap_reference(.5,[0.],[0.],alpha=.002,T_x=1.,radius_rate=0.,
        incoming_w_plus=0.,basis='legendre',background=bg,background_time=-1.,budget=budget)
    assert row['mass']==row['gradient']['T']==row['explicit_time_remainder']==0.


def test_nonzero_background_mass_underflow_is_rejected(budget):
    bg=background(Polynomial([0.,1e-200]))
    with pytest.raises(ValueError,match='mass is unresolved'):
        polynomial_cap_reference(.5,[0.],[0.],alpha=0.,T_x=1.,radius_rate=0.,
            incoming_w_plus=0.,basis='legendre',background=bg,background_time=0.,budget=budget)


@pytest.mark.parametrize('changes',[dict(background_time=0.),dict(background=object(),background_time=0.),
                                  dict(background=background(Polynomial([0.,1.])),background_time=0.)])
def test_background_requires_its_explicit_legendre_contract(changes,budget):
    with pytest.raises(ValueError):
        polynomial_cap_reference(.5,[0.],[0.],alpha=0.,T_x=1.,radius_rate=0.,incoming_w_plus=0.,
                                 budget=budget,**changes)
