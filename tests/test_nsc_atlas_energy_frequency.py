import math
import pytest
import sympy as sp
from recursive_horizons.nsc_atlas_energy_frequency import (
    flat_equilibrium_ratios,radiation_temperature_ratio,thermal_log_rate,clock_redshift)


def test_pressure_balance_size_and_occupied_frequency_bound():
    R,t,D,h,c,n,q=sp.symbols('R t D h c n q',positive=True)
    energy=4*sp.pi*t*R**2+D/R
    Dstar=8*sp.pi*t*R**3
    assert sp.simplify(sp.diff(energy,R).subs(D,Dstar))==0
    assert sp.simplify(energy.subs(D,Dstar)-12*sp.pi*t*R**2)==0
    assert sp.simplify(sp.diff(energy,R,2).subs(D,Dstar)-24*sp.pi*t)==0
    A=n*q*(q+1)*(2*q+1)/3
    omega3=c**3*q**3*8*sp.pi*t/(h*c*A)
    bound=12*sp.pi*t*c*c/(n*h)
    assert sp.simplify(omega3/bound-2*q*q/((q+1)*(2*q+1)))==0
    assert sp.expand((q+1)*(2*q+1)-2*q*q)==3*q+1
    r=flat_equilibrium_ratios(999)
    assert r['radius_over_reference']<1 and r['total_energy_over_reference']<1
    assert r['fixed_level_frequency_over_reference']>1
    assert r['highest_occupied_frequency_over_reference']<1


def test_first_law_and_cooling_threshold_do_not_confuse_total_energy_with_temperature():
    R,T,a=sp.symbols('R T a',positive=True)
    Rdot,Tdot,P=sp.symbols('Rdot Tdot P',real=True)
    V=4*sp.pi*R**3/3;E=a*T**4*V;p=E/(3*V)
    Edot=sp.diff(E,T)*Tdot+sp.diff(E,R)*Rdot
    solved=sp.solve(Edot+p*sp.diff(V,R)*Rdot+P,Tdot)[0]
    assert sp.simplify(solved/T+Rdot/R+P/(4*E))==0
    assert radiation_temperature_ratio(.5,.5)==pytest.approx(math.sqrt(2))
    assert radiation_temperature_ratio(.5,1/32)==pytest.approx(1/math.sqrt(2))
    assert thermal_log_rate(1,-1,1,4)==0
    assert thermal_log_rate(1,-1,1,5)<0
    assert thermal_log_rate(1,-1,1,3)>0


def test_clock_calibrated_redshift_uses_both_photon_and_reference_frequency():
    ae,ao,oe,oo,nu0=sp.symbols('ae ao oe oo nu0',positive=True)
    received=nu0*oe*ae/ao
    reference=nu0*oo
    assert sp.simplify(reference/received-(ao*oo)/(ae*oe))==0
    assert clock_redshift(1,.3)<1
    assert clock_redshift(-1,-.3)>1
    # A decreasing Einstein-frame scale can look expanding only because
    # the dimensionless clock-calibrated size a_E*Omega actually grows.
    assert clock_redshift(.2,0)==pytest.approx(math.exp(.2))


def test_energy_eigenprojector_stays_static_while_phases_rotate():
    t,w1,w2=sp.symbols('t w1 w2',real=True)
    U=sp.diag(sp.exp(-sp.I*w1*t),sp.exp(-sp.I*w2*t))
    rho=sp.diag(1,0)
    assert sp.simplify(U*rho*sp.conjugate(U.T)-rho)==sp.zeros(2)
