"""Independent high-precision field and derivative controls for the background."""
import mpmath as mp
import numpy as np
from numpy.polynomial import Polynomial
import pytest

from control_budget import ControlBudget
from recursive_horizons.nsc_atlas_cap_background_jets import RetainedCapBackground


CONTROLS = ControlBudget(maximum_evaluations=300000,cpu_seconds=20.)


@pytest.fixture
def budget():
    with CONTROLS.active() as value:
        yield value


@pytest.fixture(scope='module',autouse=True)
def receipt():
    yield
    print('background jet controls:',CONTROLS.budget.evaluations,'charged values;',
          CONTROLS.active_cpu_seconds,'active CPU seconds')


def background(poly,knots):
    knots=np.array(knots,dtype=float)
    return RetainedCapBackground(knots,poly(knots),poly.deriv()(knots),poly.deriv(2)(knots))


def exact_fields(F,T,z):
    r,A0=mp.tanh(z),1/mp.cosh(z)**2
    chi=F(T+z)-F(T-z)
    k=mp.diff(F,T+z)-mp.diff(F,T-z)
    q=(mp.diff(F,T+z)+mp.diff(F,T-z))/A0
    return dict(f=chi/r,p=k/A0,h=q-chi/r)


@pytest.mark.parametrize('point',[1e-12,.005,.01,.01000001,.4])
def test_quintic_fields_and_both_derivatives_against_85_digits(point,budget):
    coefficients=[.5,-.25,.125,.5,-.125,.25]
    poly=Polynomial(coefficients)
    row=background(poly,[-1.,0.,.5,1.]).reduced_jets(.125,point,budget=budget)
    budget.check(1000)
    with mp.workdps(85):
        F=lambda value:sum(mp.mpf(c)*value**i for i,c in enumerate(coefficients))
        T,z=mp.mpf(.125),mp.mpf(point)
        exact=exact_fields(F,T,z)
        for name in ('f','p','h'):
            references={name:exact[name],name+'_T':mp.diff(lambda t:exact_fields(F,t,z)[name],T),
                        name+'_z':mp.diff(lambda x:exact_fields(F,T,x)[name],z)}
            for key,value in references.items():
                error=abs(mp.mpf(float(row[key]))-value)
                assert error <= mp.mpf('2e-10')*max(abs(value),mp.mpf('1e-35')),(point,key,error,value)
        assert row['chi_z'] == pytest.approx(float(mp.diff(F,T+z)+mp.diff(F,T-z)),abs=1e-13)
        assert row['k_T'] == pytest.approx(float(mp.diff(F,T+z,2)-mp.diff(F,T-z,2)),abs=1e-13)


def test_c2_axis_knot_uses_both_polynomial_pieces(budget):
    row=background(Polynomial([0,0,0,1]),[0.,.5,1.]).reduced_jets(0.,1e-8,budget=budget)
    with mp.workdps(85):
        F=lambda value:max(value,mp.mpf(0))**3
        z=mp.mpf(1e-8)
        for name in ('f','p','h'):
            exact=exact_fields(F,mp.mpf(0),z)[name]
            dT=mp.diff(lambda t:exact_fields(F,t,z)[name],mp.mpf(0))
            dz=mp.diff(lambda x:exact_fields(F,mp.mpf(0),x)[name],z)
            for key,value in [(name,exact),(name+'_T',dT),(name+'_z',dz)]:
                assert abs(mp.mpf(float(row[key]))-value)<mp.mpf('1e-12')*abs(value)
    budget.check(1000)


@pytest.mark.parametrize('T,z',[(.001,.009),(.001,.011),(0.,1e-12)])
def test_interior_c2_knot_is_compared_to_exact_piecewise_polynomials(T,z,budget):
    knots=np.array([-.02,0.,.02])
    factor=np.where(knots<0,2.,1.)
    item=RetainedCapBackground(knots,factor*knots**3,3*factor*knots**2,6*factor*knots)
    row=item.reduced_jets(T,z,budget=budget)
    with mp.workdps(85):
        F=lambda value:(2 if value<0 else 1)*value**3
        tt,zz=mp.mpf(T),mp.mpf(z)
        for name in ('f','p','h'):
            refs={name:exact_fields(F,tt,zz)[name],
                  name+'_T':mp.diff(lambda t:exact_fields(F,t,zz)[name],tt),
                  name+'_z':mp.diff(lambda x:exact_fields(F,tt,x)[name],zz)}
            for key,value in refs.items():
                # Coefficient conversion leaves tiny stored knot jumps; they
                # are numerical error, not extra singular physical sources.
                assert abs(mp.mpf(float(row[key]))-value)<mp.mpf('2e-13')
    budget.check(1000)


def test_axis_and_quiet_values_and_original_inputs_are_retained(budget):
    poly=Polynomial([0.,0.,1.])
    item=background(poly,[-1.,0.,1.])
    row=item.reduced_jets(.125,0.,budget=budget)
    # Existing Hermite coefficient conversion is floating arithmetic.
    np.testing.assert_allclose([row['f'],row['f_T'],row['p_z']],[.5,4.,4.],rtol=0,atol=2e-14)
    assert all(row[key]==0 for key in ('chi','k','k_T','p','h','f_z','h_z','p_T','h_T'))
    copied=item.inputs
    copied['values'][:]=99.
    assert item.inputs['values'].tolist()==[1.,0.,1.]
    quiet=background(Polynomial([0,0,0,1]),[0.,1.]).reduced_jets(-1.,.4,budget=budget)
    assert all(value==0 for key,value in quiet.items() if key not in ('r','A0'))


@pytest.mark.parametrize('T,z',[(1.,.1),(0.,-1.),(np.nan,.1)])
def test_domain_is_not_extrapolated(T,z,budget):
    item=background(Polynomial([0,0,0,1]),[0.,1.])
    with pytest.raises(ValueError):
        item.reduced_jets(T,z,budget=budget)


def test_shared_budget_is_enforced(budget):
    item=background(Polynomial([0,0,0,1]),[0.,1.])
    maximum=budget.maximum_evaluations
    budget.maximum_evaluations=budget.evaluations+1
    try:
        with pytest.raises(RuntimeError,match='budget exhausted'):
            item.reduced_jets(.2,.1,budget=budget)
    finally:
        budget.maximum_evaluations=maximum
