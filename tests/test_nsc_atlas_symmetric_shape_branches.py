"""Stationary-family isolation and retained scalar indices, without PDE solves."""
from fractions import Fraction as F
import math

import numpy as np
import pytest
from scipy.linalg import eigh
import sympy as sp

from recursive_horizons.nsc_atlas_sphere_flux_completion import (
    symmetric_radius_branches,log_radius_potential,log_radius_gradient,
    log_radius_hessian,log_radius_kinetic_metric,all_volume_extrema,
)


@pytest.fixture(scope='module')
def record():return symmetric_radius_branches()


def test_exact_positive_isolations_are_exhaustive_only_within_the_named_families(record):
    t=sp.symbols('t');A=sp.Integer(7);a=sp.Rational(5,4);S=t*t+t+1
    expected={
        '1+3':5*4**5*A**5*S**6-12*a**5*(t**3+28*t*t+28*t+15)**5,
        '2+2':5*2**5*A**5*(S*(t**3+1))**6-24*a**5*(5*t**6+14*t**5+14*t**4+6*t**3+14*t*t+14*t+5)**5,
    }
    for row in record['isolations']:
        polynomial=sp.Poly.from_list(list(row['integer_polynomial_coefficients']),t)
        primitive=sp.Poly(expected[row['branch']],t).clear_denoms()[1].primitive()[1]
        assert polynomial.as_expr()==primitive.as_expr() and polynomial.degree()=={'1+3':15,'2+2':30}[row['branch']]
        intervals=row['positive_root_intervals']
        assert len(intervals)==polynomial.count_roots(0,sp.oo)=={'1+3':3,'2+2':2}[row['branch']]
        for item in intervals:
            lo,hi=item['interval']
            assert isinstance(lo,F) and isinstance(hi,F) and 0<lo<=hi
            assert hi-lo<=F(1,10**14)
            assert polynomial.count_roots(sp.Rational(lo),sp.Rational(hi))==item['multiplicity']==1


def test_all_reported_points_solve_full_gradient_and_canonical_mass_problem(record):
    G=log_radius_kinetic_metric()
    for row in record['branches']:
        u=np.array(row['log_radii']);H=log_radius_hessian(u)
        np.testing.assert_allclose(log_radius_gradient(u),0.,atol=8e-11,rtol=0)
        assert row['V']==pytest.approx(log_radius_potential(u),abs=3e-15)
        modes=np.array(row['canonical_mass_modes']);masses=np.array(row['canonical_mass_squared'])
        np.testing.assert_allclose(modes.T@G@modes,np.eye(4),atol=2e-14,rtol=0)
        np.testing.assert_allclose(H@modes,(G@modes)*masses[None,:],atol=3e-13,rtol=2e-13)
        np.testing.assert_allclose(masses,eigh(H,G,eigvals_only=True),atol=3e-14,rtol=2e-14)
        assert row['homogeneous_negative_count']==int(sum(masses<0))
        if row['branch']!='isotropic':
            x=np.exp(-2*u);t=row['ratio']
            assert x[0]/x[-1]==pytest.approx(t,rel=3e-15)
            S=t*t+t+1
            if row['branch']=='1+3':assert x[-1]**5==pytest.approx(12/(5*S),rel=3e-15)
            else:assert x[-1]**5==pytest.approx(24/(5*S*(t**3+1)),rel=3e-15)


def test_reference_false_and_AdS_radius_minima_and_unique_index_one_branch(record):
    rows=record['branches'];isotropic=[r for r in rows if r['branch']=='isotropic']
    for point in all_volume_extrema(1.5):
        row=next(r for r in isotropic if r['volume_kind']==point['kind'])
        assert row['V']==pytest.approx(point['vacuum'],abs=3e-14)
    false=next(r for r in isotropic if r['volume_kind']=='minimum')
    assert false['retained_radius_local_minimum'] and false['V']>0
    assert false['homogeneous_negative_count']==false['S4_scalar_harmonic_negative_count']==0
    ads=[r for r in rows if r['V']<0]
    assert len(ads)==1 and ads[0]['retained_radius_local_minimum'] and ads[0]['homogeneous_negative_count']==0
    assert ads[0]['S4_scalar_harmonic_negative_count'] is None
    index_one=[r for r in rows if r['branch']=='1+3' and r['homogeneous_negative_count']==1]
    assert len(index_one)==1 and index_one[0]['S4_scalar_harmonic_negative_count']==6
    barrier=next(r for r in isotropic if r['volume_kind']=='maximum')
    assert barrier['homogeneous_negative_count']==barrier['S4_scalar_harmonic_negative_count']==4
    assert sorted(r['homogeneous_negative_count'] for r in rows if r['branch']=='1+3')==[0,1,3]
    assert next(r for r in rows if r['branch']=='2+2')['homogeneous_negative_count']==3


def test_S4_counts_use_strict_harmonic_inequality_and_correct_degeneracies(record):
    for row in record['branches']:
        if row['V']<=0:continue
        count=0
        for mass,channel in zip(row['canonical_mass_squared'],row['S4_negative_harmonics']):
            cutoff=-3*mass/row['V'];expected=[];ell=0
            while ell*(ell+3)<cutoff:
                expected.append(ell);count+=(ell+1)*(ell+2)*(2*ell+3)//6;ell+=1
            assert channel['ell']==tuple(expected)
            assert channel['degeneracy']==tuple(math.comb(l+4,4)-math.comb(l+2,4) if l>=2 else (1 if l==0 else 5) for l in expected)
        assert row['S4_scalar_harmonic_negative_count']==count
    assert symmetric_radius_branches() is record
