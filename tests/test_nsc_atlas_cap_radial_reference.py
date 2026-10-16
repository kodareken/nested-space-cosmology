"""Manufactured radial maps, with independent exact and fixed-r variations."""
import mpmath as mp
import numpy as np
import pytest
from scipy.integrate import solve_ivp

from control_budget import ControlBudget
from recursive_horizons.nsc_atlas_cap_radial_reference import (
    affine_radial_map, reference_potential,
)

CONTROLS = ControlBudget(maximum_evaluations=200000, cpu_seconds=20.)


@pytest.fixture
def budget():
    with CONTROLS.active() as value:
        yield value


@pytest.fixture(scope="module", autouse=True)
def report():
    yield
    print("radial reference controls:",CONTROLS.budget.evaluations,
          "charged values;",CONTROLS.active_cpu_seconds,"active CPU seconds")


def test_smooth_potential_and_tiny_native_gap_against_high_precision(budget):
    values = np.array([0.,1e-8,-1e-5,np.nextafter(1e-4,0.),1e-4,.2,-.5])
    row = reference_potential(values,6.,budget=budget)
    with mp.workdps(120):
        for i,value in enumerate(values):
            f = mp.mpf(float(value))
            def potential(z):
                sinc = mp.sinh(z)/z if z else mp.mpf(1)
                return z*z*sinc*sinc*(mp.cosh(z)**2+1)/2
            exact, exact_u = potential(f),mp.diff(potential,f)
            np.testing.assert_allclose(row['Pdiff'][i],float(exact),rtol=1e-14,atol=0.)
            np.testing.assert_allclose(row['Pdiff_f'][i],float(exact_u),rtol=1e-14,atol=0.)
            if abs(value) < 1e-4:
                def native(z):
                    series = 1+z*z/6+z**4/120
                    return z*z*series**2*(mp.cosh(z)**2+1)/2
                gaps = (exact-native(f),exact_u-mp.diff(native,f))
            else:
                gaps = (mp.mpf(0),mp.mpf(0))
            for name,bound,expected in zip(('implementation_gap','implementation_gap_f'),
                    ('gap_series_bound','gap_f_series_bound'),gaps):
                error = abs(mp.mpf(float(row[name][i]))-expected)
                assert error <= 64*np.finfo(float).eps*abs(expected)+mp.mpf(float(row[bound][i]))
    assert row['implementation_gap'][3] > 0
    assert row['implementation_gap'][4] == 0
    assert row['implementation_gap_f'][2] < 0


def test_zero_coupling_and_invalid_reference_potential(budget):
    f=np.array([-2.,0.,3.])
    row=reference_potential(f,0.,budget=budget)
    np.testing.assert_array_equal(row['Pdiff'],f*f)
    np.testing.assert_array_equal(row['Pdiff_f'],2*f)
    for key in ('implementation_gap','implementation_gap_f','gap_series_bound','gap_f_series_bound'):
        np.testing.assert_array_equal(row[key],0.)
    for f,alpha in (([np.nan],0.),([1.],-1.),([1.]*514,0.),([1e200],6.)):
        with pytest.raises(ValueError): reference_potential(f,alpha,budget=budget)


@pytest.mark.parametrize('f',[5e-5,.2])
def test_potential_driven_radial_map_matches_exact_constant_field_integral(f,budget):
    p=reference_potential([f],6.,budget=budget)
    P,U=float(p['Pdiff'][0]),float(p['Pdiff_f'][0]); fx=.07
    l,r,lx,rx=.1,.6,.03,-.02
    def coefficients(at):
        return dict(a=0.,b=at*at*P/2,a_x=0.,b_x=at*at*U*fx/2,a_r=0.,b_r=at*P)
    row=affine_radial_map(coefficients,left=l,right=r,left_rate=lx,right_rate=rx,budget=budget)
    volume=(r**3-l**3)/6; volume_rate=(r*r*rx-l*l*lx)/2
    assert row['S']==1. and row['S_x']==0.
    np.testing.assert_allclose(row['T'],P*volume,rtol=1e-13,atol=0.)
    np.testing.assert_allclose(row['T_x'],U*fx*volume+P*volume_rate,rtol=1e-13,atol=0.)
    # For this manufactured Pi=Phi=0 case the native Gauss mass integrates
    # r^2 exactly. The potential implementation gap is its only analytic
    # difference from the reference, carried independently without a
    # cancellation-prone subtraction of the two rounded masses.
    if f < 1e-4:
        assert -p['implementation_gap'][0]*volume < 0.
        assert -p['implementation_gap_f'][0]*fx*volume-p['implementation_gap'][0]*volume_rate < 0.
    else:
        assert p['implementation_gap'][0]==p['implementation_gap_f'][0]==0.


@pytest.mark.parametrize('a,ax',[(0.,.3),(1.4,-.2)])
def test_constant_coefficients_with_independently_moving_endpoints(a,ax,budget):
    b,bx,l,r,lx,rx = .7,.2,.1,.6,.03,-.02
    coefficients=lambda at:dict(a=a,b=b,a_x=ax,b_x=bx,a_r=0.,b_r=0.)
    row=affine_radial_map(coefficients,left=l,right=r,left_rate=lx,right_rate=rx,budget=budget)
    length,rate = r-l,rx-lx
    if a == 0:
        expected=dict(S=1.,T=b*length,S_x=-ax*length,
                      T_x=bx*length+b*rate-b*ax*length**2/2)
    else:
        S=np.exp(-a*length); Sx=-S*(ax*length+a*rate)
        expected=dict(S=S,T=b/a*(1-S),S_x=Sx,
            T_x=(bx*a-b*ax)/(a*a)*(1-S)-b/a*Sx)
    for key,value in expected.items():
        np.testing.assert_allclose(row[key],value,rtol=3e-11,atol=1e-13)
    assert row['evaluations'] == 16*row['coefficient_calls']
    assert not row['born'] and not row['uniform_error_enclosure']


def polynomial_coefficients(r):
    a0,ax,u,ux,v,vx=3.,.2,.1,.3,.2,-.4
    p,px,pr=r*r+u*r+v,ux*r+vx,2*r+u
    return dict(a=a0*r,b=pr+a0*r*p,a_x=ax*r,b_x=ux+ax*r*p+a0*r*px,
                a_r=a0,b_r=2+a0*p+a0*r*pr)


def test_nonconstant_polynomial_particular_solution_and_full_map_derivative(budget):
    l,r,lx,rx=.08,.7,.04,-.03
    row=affine_radial_map(polynomial_coefficients,left=l,right=r,left_rate=lx,right_rate=rx,budget=budget)
    p=lambda z:z*z+.1*z+.2
    pd=lambda z,rate:.3*z-.4+(2*z+.1)*rate
    S=np.exp(-1.5*(r*r-l*l))
    Sx=-S*(.1*(r*r-l*l)+3*(r*rx-l*lx))
    expected=dict(S=S,S_x=Sx,T=p(r)-S*p(l),T_x=pd(r,rx)-Sx*p(l)-S*pd(l,lx))
    for key,value in expected.items():
        np.testing.assert_allclose(row[key],value,rtol=5e-11,atol=2e-13)


def test_independent_fixed_radius_variation_includes_input_and_boundary_motion(budget):
    l,r,lx,rx=.08,.7,.04,-.03
    incoming,incoming_rate=1.2,-.17
    row=affine_radial_map(polynomial_coefficients,left=l,right=r,left_rate=lx,right_rate=rx,budget=budget)
    at_left=polynomial_coefficients(l)
    nu0=incoming_rate-lx*(at_left['b']-at_left['a']*incoming)
    def variation(at,state):
        budget.check(8)
        fields=polynomial_coefficients(at); n,nu=state
        return [fields['b']-fields['a']*n,
                fields['b_x']-fields['a_x']*n-fields['a']*nu]
    direct=solve_ivp(variation,(l,r),[incoming,nu0],method='DOP853',rtol=2e-12,atol=2e-14,max_step=.08)
    assert direct.success
    n,nu=direct.y[:,-1]; end=polynomial_coefficients(r)
    material=nu+rx*(end['b']-end['a']*n)
    np.testing.assert_allclose(row['S']*incoming+row['T'],n,rtol=6e-11,atol=3e-13)
    mapped=row['S_x']*incoming+row['S']*incoming_rate+row['T_x']
    np.testing.assert_allclose(mapped,material,rtol=6e-11,atol=3e-13)
    assert abs(mapped-nu)>1e-4


def test_backward_born_panel_keeps_its_derivative_and_rejects_the_other_ray(budget):
    r=.3; fields=polynomial_coefficients(r)
    row=affine_radial_map(polynomial_coefficients,left=r,right=r,left_rate=.2,right_rate=.1,budget=budget)
    assert row['S']==1 and row['T']==0 and row['born']
    assert row['S_x']==pytest.approx(.1*fields['a'])
    assert row['T_x']==pytest.approx(-.1*fields['b'])
    assert row['coefficient_calls']==1
    with pytest.raises(ValueError,match='backward ray'):
        affine_radial_map(polynomial_coefficients,left=r,right=r,left_rate=.1,right_rate=.2,budget=budget)


def test_reference_panel_errors_remain_local_and_budgeted(budget):
    args=dict(left=.1,right=.3,left_rate=0.,right_rate=0.,budget=budget)
    for callback in (lambda r:{},lambda r:dict(a=np.nan,b=0.,a_x=0.,b_x=0.,a_r=0.,b_r=0.)):
        with pytest.raises(ValueError): affine_radial_map(callback,**args)
    for invalid in (dict(right=.05),dict(right=1.),dict(rtol=1e-16),dict(atol=0.)):
        with pytest.raises(ValueError): affine_radial_map(polynomial_coefficients,**(args|invalid))
    maximum=budget.maximum_evaluations
    budget.maximum_evaluations=budget.evaluations+15
    try:
        with pytest.raises(RuntimeError,match='budget exhausted'):
            affine_radial_map(polynomial_coefficients,**args)
    finally:
        budget.maximum_evaluations=maximum
