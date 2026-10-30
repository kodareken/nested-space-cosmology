import math
import cmath
import pytest
import sympy as sp
from recursive_horizons.nsc_atlas_three_targets import calculation,room_band,sphere_band,orientation_imbalance,outside_energy_share
from recursive_horizons.unified_action import RecursiveScaleDilation


def test_band_definition_and_grid_refinement_are_not_a_selected_physical_cutoff():
    n,b,r,h=sp.symbols('n b r h',positive=True)
    assert sp.simplify(((n+b)/h)/(n/r)-(r/h)*(1+b/n))==0
    x=room_band(1.,.2,3.,.01)
    assert room_band(1.,.2,3.,.005)==2*x
    assert room_band(1.,.2,30.,.1)==pytest.approx(x)
    assert sphere_band(1,1)==1 and sphere_band(1000,1)==1000
    assert sphere_band(1000,2)==2000
    # At the angular ground, the stated wavelength is a circumference.
    c,R=sp.symbols('c R',positive=True)
    assert sp.simplify((2*sp.pi*c/(c/R))/R)==2*sp.pi


def test_existing_frequency_dilation_is_not_a_curvature_dilation_of_the_same_factor():
    omega=3.973074368754331
    d=RecursiveScaleDilation(omega)
    assert d.child_cutoff(2)==pytest.approx(2*omega)
    assert d.zeta==pytest.approx(omega**2)
    assert 1/math.sqrt(d.zeta)==pytest.approx(1/omega)
    assert d.compose(d).omega==pytest.approx(omega**2)
    assert not math.isclose(1/omega,1/math.sqrt(omega))


def test_phase_does_not_select_conjugate_population_imbalance():
    for phase in (0.,.3,math.pi):
        a=math.sqrt(.5)*cmath.exp(1j*phase)
        b=math.sqrt(.5)
        assert orientation_imbalance(abs(a)**2,abs(b)**2)==pytest.approx(0.,abs=2e-16)
    assert orientation_imbalance(.75,.25)==.5
    A,n,ng=sp.symbols('A n ng',positive=True)
    # A supplied orientation fraction still needs a baryon/photon conversion.
    eta=A*n/ng
    assert sp.diff(eta,ng)!=0


def test_density_ratio_does_not_select_outside_energy_or_cosmological_components():
    assert outside_energy_share(1)==.75
    assert outside_energy_share(1/3)==.5
    v=sp.symbols('v',positive=True)
    assert sp.simplify(sp.diff(3*v/(1+3*v),v)-3/(1+3*v)**2)==0


def test_targets_are_comparisons_only_and_missing_predictions_are_not_zero():
    a=calculation(3.973074368754331,.220983024212607,1000)
    b=calculation(3.973074368754331,.220983024212607,1000,targets=(120,130))
    assert a['inputs']==b['inputs']
    assert a['band']['inheritance']['conditional_rows'][2]['log_band_if_endpoint_factor_is_one']==pytest.approx(137.9540195195991)
    assert a['band']['inheritance']['curvature_interpretation_control']['log_band_if_endpoint_factor_is_one']==pytest.approx(137.9540195195991)
    assert a['band']['prediction'] is None and a['asymmetry']['prediction'] is None
    assert all(v is None for v in a['split']['predictions'].values())
    assert not a['physical_targets_used_for_selection']
