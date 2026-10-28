import math

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.linalg import solve

from recursive_horizons.nsc_atlas_parent_supernova import (
    SupernovaShapeLikelihood, fit_shape, flrw_distance, magnitude_error_bound,
    parent_distance, parent_fractions, parent_ray_angle_upper,
)
from recursive_horizons.nsc_atlas_inherited_matter import exact_initial_geometry


def test_flrw_limits_and_independent_null_integral():
    z = np.array([.03, .2, 1., 2.3])
    hel = z+.001
    np.testing.assert_allclose(flrw_distance(z,hel,matter=0), (1+hel)*z, rtol=2e-15)
    np.testing.assert_allclose(flrw_distance(z,hel,matter=1),
                               (1+hel)*2*(1-1/np.sqrt(1+z)), rtol=3e-15)
    for matter, curvature in ((.3,0.), (.4,-.1), (.2,.2)):
        expected=[]
        for zz,hh in zip(z,hel):
            I = quad(lambda x: 1/math.sqrt(matter*(1+x)**3+curvature*(1+x)**2+1-matter-curvature),0,zz,
                     epsabs=2e-13, epsrel=2e-13)[0]
            radial = (math.sin(math.sqrt(-curvature)*I)/math.sqrt(-curvature) if curvature<0 else
                      math.sinh(math.sqrt(curvature)*I)/math.sqrt(curvature) if curvature>0 else I)
            expected.append((1+hh)*radial)
        np.testing.assert_allclose(flrw_distance(z,hel,matter=matter,curvature=curvature),expected,rtol=3e-14)


def test_parent_preserves_source_and_domain_and_density_relation():
    z=np.array([.03,.25])
    d=exact_initial_geometry()
    for s in (.2,.5,.8):
        f=parent_fractions(s)
        invariant=(-f['curvature'])**3/(f['matter']**2*f['vacuum'])
        expected=float(d['curvature']**3/((d['excitation']/3)**2*(d['vacuum']/3)))
        assert invariant==pytest.approx(expected,rel=8e-15)
        Ao=1/s
        a0=1/math.sqrt(float(d['curvature']))
        I=quad(lambda A:1/(A*A*math.sqrt(float(d['vacuum']/3)+float(d['excitation']/3)/A**3-float(d['curvature'])/A**2)),
               Ao/1.25,Ao,epsabs=2e-13)[0]
        direct=1.25*Ao*a0*math.sin(I/a0)*math.sqrt(f['H_squared'])
        assert parent_distance(z,z,s)[1]==pytest.approx(direct,rel=2e-14)
    with pytest.raises(ValueError,match='precedes'):
        parent_distance(z,z,.81)
    assert parent_ray_angle_upper(2.3,1/3.3)<math.pi/2
    assert 0<magnitude_error_bound()<8e-10


def test_profile_offset_matches_augmented_normal_equation_and_units():
    distance=np.array([.03,.07,.18,.3])
    L=np.diag([.1,.2,.15,.13])+np.tril(np.ones((4,4))*.012,-1)
    C=L@L.T
    y=5*np.log10(distance)+18+np.array([.03,-.02,.04,-.01])
    ll=SupernovaShapeLikelihood(y,C)
    row=ll.evaluate(distance)
    # Independent full normal equation, without the production Cholesky path.
    inv=solve(C,np.eye(4),assume_a='sym')
    residual=y-5*np.log10(distance)
    offset=np.sum(inv@residual)/np.sum(inv)
    assert row['magnitude_intercept']==pytest.approx(offset,abs=2e-14)
    assert row['chi_squared']==pytest.approx((residual-offset)@inv@(residual-offset),abs=1e-13)
    assert ll.evaluate(distance*1e12)['chi_squared']==pytest.approx(row['chi_squared'],abs=2e-13)


def test_synthetic_shape_recovery_and_boundaries():
    z=np.linspace(.025,.7,50)
    target=flrw_distance(z,z,matter=.37)
    ll=SupernovaShapeLikelihood(5*np.log10(target)+18,np.eye(50)*.02**2)
    result=fit_shape(ll,lambda m:flrw_distance(z,z,matter=m),1.)
    assert result['parameter']==pytest.approx(.37,abs=2e-8)
    assert result['chi_squared']<1e-10
    vacuum=SupernovaShapeLikelihood(5*np.log10((1+z)*z)+18,np.eye(50)*.02**2)
    assert fit_shape(vacuum,lambda m:flrw_distance(z,z,matter=m),1.)['parameter']==0
