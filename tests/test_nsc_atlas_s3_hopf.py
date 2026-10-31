"""FFT transforms against the independently evaluated regular Wigner basis."""
import math

import numpy as np
import pytest

from recursive_horizons.nsc_atlas_s3_frame import build_frame_basis,arbitrary_points_values
from recursive_horizons.nsc_atlas_s3_hopf import build_hopf_frame_basis


@pytest.fixture(scope='module')
def hopf():return build_hopf_frame_basis(4)


def test_hopf_points_and_weights_cover_full_s3(hopf):
    b=hopf;X=b.points
    np.testing.assert_allclose(np.sum(X*X,axis=1),1,atol=1e-15,rtol=0)
    assert np.all(b.weights>0)
    assert sum(b.weights)==pytest.approx(2*math.pi**2,rel=3e-14)
    for coordinate in X.T:
        assert np.dot(b.weights,coordinate)==pytest.approx(0,abs=2e-14)
        assert np.dot(b.weights,coordinate**2)==pytest.approx(math.pi**2/2,rel=3e-14)
    assert np.dot(b.weights,X[:,0]**2*X[:,3]**2)==pytest.approx(math.pi**2/12,rel=3e-14)


def test_fft_and_direct_wigner_evaluation_agree_for_vector_and_batch(hopf):
    b=hopf;M=len(b.labels)
    c=np.cos(np.arange(M))/10
    C=np.c_[c,-2*c,np.sin(np.arange(M))/13]
    direct=arbitrary_points_values(b.points,b.labels)
    np.testing.assert_allclose(b.evaluate(c),direct@c,atol=3e-15,rtol=2e-13)
    np.testing.assert_allclose(b.evaluate(C),direct@C,atol=8e-15,rtol=2e-13)


def test_l2_and_l4_fft_transforms_are_orthonormal_without_a_dense_runtime_matrix():
    for L in (2,4):
        b=build_hopf_frame_basis(L);I=np.eye(len(b.labels))
        values=b.evaluate(I)
        np.testing.assert_allclose(b.project(values),I,atol=5e-14,rtol=0)
        np.testing.assert_allclose(values.T@(b.weights[:,None]*values),I,atol=5e-14,rtol=0)
        assert not hasattr(b,'values')


def test_sine_projection_phase_and_known_ambient_polynomials(hopf):
    b=hopf;X=b.points
    x=b.project(X[:,0])
    np.testing.assert_allclose(b.evaluate(x),X[:,0],atol=4e-15,rtol=0)
    assert x[b.labels.index((1,0,1,'sin'))]==pytest.approx(-math.pi/math.sqrt(2),rel=2e-14)
    f=.2*X[:,0]+.7*X[:,1]*X[:,3]
    c=b.project(f)
    np.testing.assert_allclose(b.evaluate(c),f,atol=3e-15,rtol=0)
    constant=b.project(np.ones(len(X)))
    assert constant[0]==pytest.approx(math.sqrt(2)*math.pi,rel=2e-14)
    assert np.max(abs(constant[1:]))<3e-14


def test_generators_labels_and_characters_are_the_unchanged_frame_basis(hopf):
    b=hopf;reference=build_frame_basis(4)
    assert b.labels==reference.labels
    np.testing.assert_array_equal(b.characters,reference.characters)
    c=np.sin(np.arange(len(b.labels)))/12
    direct=arbitrary_points_values(b.points,b.labels)
    for D,known in zip(b.generators,reference.generators):
        np.testing.assert_allclose((D-known).toarray(),0,atol=0,rtol=0)
        np.testing.assert_allclose(b.evaluate(D@c),direct@(known@c),atol=8e-15,rtol=2e-13)


def test_projection_matches_direct_volume_integral_and_independent_quadrature(hopf):
    b=hopf;X=b.points
    f=1+.3*X[:,0]-.4*X[:,1]*X[:,3]+.07*X[:,0]**2*X[:,2]**2
    direct=arbitrary_points_values(X,b.labels)
    np.testing.assert_allclose(b.project(f),direct.T@(b.weights*f),atol=8e-15,rtol=0)
    high=build_hopf_frame_basis(4,quadrature_order=(8,19));Y=high.points
    fhigh=1+.3*Y[:,0]-.4*Y[:,1]*Y[:,3]+.07*Y[:,0]**2*Y[:,2]**2
    np.testing.assert_allclose(b.project(f),high.project(fhigh),atol=8e-15,rtol=0)


def test_quadrature_and_shape_validation(hopf):
    b=hopf
    assert b.quadrature_shape==(5,13,13)
    assert not b.points.flags.writeable and not b.weights.flags.writeable
    assert build_hopf_frame_basis(2,quadrature_order=5).quadrature_shape==(5,11,11)
    with pytest.raises(ValueError,match='two-product'):build_hopf_frame_basis(4,(2,13))
    with pytest.raises(ValueError,match='coefficients'):b.evaluate(np.zeros(3))
    with pytest.raises(ValueError,match='samples'):b.project(np.zeros(3))
