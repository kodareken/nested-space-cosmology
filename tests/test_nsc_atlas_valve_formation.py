"""Consequential pressure-release identities, independent of the finder."""
import importlib.util
import math
from pathlib import Path

import numpy as np
import sympy as s

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('valve_check',ROOT/'scripts/check_nsc_atlas_valve_formation.py')
producer=importlib.util.module_from_spec(spec);spec.loader.exec_module(producer)


def test_fixed_action_fold_from_independent_potential_derivatives():
    from recursive_horizons.nsc_atlas_parent_circulation import Parent
    variable,M=s.symbols('s M',real=True)
    V=7*s.exp(-12*variable)-s.Rational(15,2)*s.exp(-14*variable)+s.Rational(25,32)*s.exp(-24*variable)
    effective=V+M*s.exp(-6*variable)
    loading=s.solve(s.diff(effective,variable),M)[0]
    x=s.symbols('x',positive=True)
    polynomial=s.simplify(4*s.diff(loading,variable).subs(variable,-s.log(x)/2)/x**3)
    assert s.expand(polynomial)==225*x**6-560*x+336
    row=producer.stabilization_threshold(Parent.retained())
    point={variable:row['fold_s'],M:row['critical_adiabatic_loading']}
    assert abs(float(s.diff(effective,variable).subs(point)))<1e-10
    assert abs(float(s.diff(effective,variable,2).subs(point)))<1e-9
    assert row['critical_adiabatic_loading']>3/8
    assert 0<row['critical_epsilon_at_scalar_maximum']<1


def test_instantaneous_barrier_and_source_maximum_are_independent_conditions():
    from recursive_horizons.nsc_atlas_parent_circulation import Parent
    parent=Parent.retained();row=producer.stabilization_threshold(parent)
    x=(14/15)**.2;c=row['instantaneous_chi_squared_critical']
    P=lambda u,chi2:42+1200*chi2-52.5*u+9.375*u**6
    assert abs(P(x,c))<2e-14
    assert P(x,c*1.01)>0 and P(x,c*.99)<0
    z=parent.mass_squared/parent.H**2
    a=math.sqrt(12/(z+3));b=math.sqrt(72/(z+8))
    X=np.array([a/b,math.sqrt((1-(a/b)**2)/2),0.,math.sqrt((1-(a/b)**2)/2)])
    f=a*X[0]+b*X[1]*X[3]
    assert math.isclose(float(X@X),1,rel_tol=0,abs_tol=2e-16)
    assert math.isclose(f,row['maximum_initial_field_at_epsilon1'],rel_tol=2e-15)
    assert row['instantaneous_critical_epsilon_at_scalar_maximum']<row['critical_epsilon_at_scalar_maximum']
