"""Low-order quadrature, normalization, symmetry and independent checks."""
import math

import numpy as np
import pytest

from recursive_horizons.nsc_atlas_s3_spectral import (
    build_surface_basis,harmonic_values,sphere_quadrature)


@pytest.fixture(scope='module')
def low_basis():return build_surface_basis(4)


def test_full_sphere_quadrature_volume_and_polynomial_moments(low_basis):
    b=low_basis;X=b.points
    np.testing.assert_allclose(np.sum(X*X,axis=1),1,atol=2e-15,rtol=0)
    assert np.all(b.weights>0)
    assert sum(b.weights)==pytest.approx(2*math.pi**2,rel=2e-14)
    for coordinate in X.T:
        assert np.dot(b.weights,coordinate)==pytest.approx(0,abs=2e-14)
        assert np.dot(b.weights,coordinate**2)==pytest.approx(math.pi**2/2,rel=2e-14)
    assert np.dot(b.weights,X[:,0]**2*X[:,1]**2)==pytest.approx(math.pi**2/12,rel=3e-14)
    # The default L=4 quadrature reaches degree28, not only the basis norm.
    power=28
    numerator=math.factorial(power)//(2**(power//2)*math.factorial(power//2))
    denominator=2**(power//2)*math.factorial(power//2+1)
    expected=2*math.pi**2*numerator/denominator
    assert np.dot(b.weights,X[:,3]**power)==pytest.approx(expected,rel=2e-13)


def test_default_l4_and_l6_harmonics_are_orthonormal():
    for L,count in ((4,11),(6,24)):
        b=build_surface_basis(L)
        assert len(b.labels)==count
        G=b.values.T@(b.weights[:,None]*b.values)
        np.testing.assert_allclose(G,np.eye(count),atol=1.7e-13,rtol=0)
        assert all(l!=1 and j%2==0 and (l+m)%2==0 for l,j,m in b.labels)
        np.testing.assert_array_equal(b.eigenvalues,[l*(l+2) for l,_,_ in b.labels])


def test_low_modes_match_independent_ambient_polynomials(low_basis):
    b=low_basis;z,x,y,u=b.points.T
    expected={
        (0,0,0):np.ones(len(u))/(math.sqrt(2)*math.pi),
        (2,0,0):(4*u*u-1)/(math.sqrt(2)*math.pi),
        (2,2,0):(2*z*z-x*x-y*y)/math.pi,
        (2,2,2):math.sqrt(3)*(x*x-y*y)/math.pi,
        (3,2,1):-4*math.sqrt(6)*z*x*u/math.pi}
    for label,value in expected.items():
        np.testing.assert_allclose(b.values[:,b.labels.index(label)],value,atol=1e-14,rtol=2e-13)


def test_retained_harmonics_obey_all_three_source_symmetries(low_basis):
    b=low_basis;sample=b.points[::97]
    original=harmonic_values(sample,b.labels)
    for signs in ((-1,-1,1,1),(-1,1,1,-1),(1,1,-1,1)):
        transformed=harmonic_values(sample*np.array(signs),b.labels)
        np.testing.assert_allclose(transformed,original,atol=9e-14,rtol=0)


def test_projection_and_laplacian_work_for_vectors_and_batches(low_basis):
    b=low_basis
    c=np.arange(1,len(b.labels)+1,dtype=float)/100
    C=np.c_[c,-2*c]
    np.testing.assert_allclose(b.project(b.evaluate(c)),c,atol=3e-14,rtol=0)
    np.testing.assert_allclose(b.project(b.evaluate(C)),C,atol=6e-14,rtol=0)
    z,x,y,u=b.points.T
    f=u*u-sum((z*z,x*x,y*y))/3
    cf=b.project(f)
    np.testing.assert_allclose(b.laplacian_values(cf),-8*f,atol=3e-13,rtol=0)
    np.testing.assert_allclose(b.laplacian_values(C),b.values@(-b.eigenvalues[:,None]*C),atol=0,rtol=0)


def test_general_evaluation_agrees_with_tensor_grid_and_handles_poles(low_basis):
    b=low_basis
    # Chunked evaluation uses only the points, without tensor-grid metadata.
    actual=harmonic_values(b.points[::53],b.labels)
    np.testing.assert_allclose(actual,b.values[::53],atol=1e-14,rtol=2e-13)
    poles=np.array([[0,0,0,1.],[0,0,0,-1.],[1.,0,0,0],[-1.,0,0,0]])
    value=harmonic_values(poles,b.labels)
    assert np.all(np.isfinite(value))
    assert value[0,b.labels.index((0,0,0))]==pytest.approx(1/(math.sqrt(2)*math.pi))


def test_quintic_source_projection_is_stable_under_independent_higher_quadrature(low_basis):
    b=low_basis
    c=np.zeros(len(b.labels));c[0]=math.sqrt(2)*math.pi
    c[1:]=np.arange(1,len(c))/300
    psi=b.evaluate(c);X=b.points
    source=(1+.07*(X[:,0]+X[:,1]*X[:,3])**2)*psi**5
    projected=b.project(source)
    points,weights=sphere_quadrature(18,19,35)
    Y=harmonic_values(points,b.labels);other=Y@c
    source_other=(1+.07*(points[:,0]+points[:,1]*points[:,3])**2)*other**5
    np.testing.assert_allclose(Y.T@(weights*source_other),projected,atol=9e-13,rtol=0)


def test_validation_and_immutable_basis_arrays(low_basis):
    b=low_basis
    assert not b.values.flags.writeable and not b.weights.flags.writeable
    with pytest.raises(ValueError,match='coefficients'):b.evaluate(np.ones(2))
    with pytest.raises(ValueError,match='samples'):b.project(np.ones(2))
    with pytest.raises(ValueError,match='unit'):harmonic_values(np.array([[0.,0,0,2]]),b.labels)
    with pytest.raises(ValueError,match='nonnegative'):build_surface_basis(-1)
    with pytest.raises(ValueError,match='positive'):sphere_quadrature(0,3,5)
