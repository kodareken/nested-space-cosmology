"""Independent exact-metric controls for the null-curvature diagnostic."""
import numpy as np
import pytest

from recursive_horizons.nsc_atlas_null_geometry import (
    first_connection, geometry_from_jets, optical_tidal,
)


def _zeros(batch=()):
    return (np.zeros(batch + (4, 4)), np.zeros(batch + (3, 4, 4)),
            np.zeros(batch + (3, 3, 4, 4)), np.zeros(batch + (3, 4, 4)),
            np.zeros(batch + (4, 4)))


def _symmetries(curvature, atol=2e-12):
    R = curvature['Riemann']
    np.testing.assert_allclose(R + R.swapaxes(-4, -3), 0, atol=atol)
    np.testing.assert_allclose(R + R.swapaxes(-2, -1), 0, atol=atol)
    np.testing.assert_allclose(R, R.swapaxes(-4, -2).swapaxes(-3, -1), atol=atol)
    # R_abcd + R_acdb + R_adbc = 0.
    np.testing.assert_allclose(R + np.einsum('...acdb->...abcd', R)
                               + np.einsum('...adbc->...abcd', R), 0, atol=atol)
    np.testing.assert_allclose(curvature['Ricci'], curvature['Ricci'].swapaxes(-1, -2), atol=atol)


def test_minkowski_commuting_frame_and_batch():
    batch = (2, 3)
    g = np.broadcast_to(np.diag([-1., 1., 1., 1.]), batch + (4, 4))
    jets = _zeros(batch)
    result = geometry_from_jets(g, *jets, frame_scale=0)
    np.testing.assert_array_equal(result['inverse'], g)
    for name in ('Gamma', 'Riemann', 'Ricci', 'commutator_error'):
        np.testing.assert_array_equal(result[name], 0)
    first = first_connection(g, jets[0], jets[1], frame_scale=0)
    np.testing.assert_array_equal(first['Gamma'], result['Gamma'])


def test_closed_de_sitter_full_riemann_and_zero_null_tidal():
    H = .37
    time = np.array([-.8, 0., .5, 2.])
    a = np.cosh(H*time)
    da = H*np.sinh(H*time)
    dda = H**2*a
    g = np.zeros((len(time), 4, 4))
    g[:, 0, 0] = -1
    v, dg, d2g, dv, accg = _zeros((len(time),))
    for i in range(1, 4):
        g[:, i, i] = a*a
        v[:, i, i] = 2*a*da
        accg[:, i, i] = 2*(da*da+a*dda)
    result = geometry_from_jets(g, v, dg, d2g, dv, accg, frame_scale=H)
    expected = H**2*(np.einsum('nac,nbd->nabcd', g, g)-np.einsum('nad,nbc->nabcd', g, g))
    np.testing.assert_allclose(result['Riemann'], expected, atol=7e-16)
    np.testing.assert_allclose(result['Ricci'], 3*H**2*g, atol=7e-16)
    k = np.zeros((len(time), 4)); k[:, 0] = 1; k[:, 1] = 1/a
    screen = np.zeros((len(time), 2, 4)); screen[:, 0, 2] = 1/a; screen[:, 1, 3] = 1/a
    np.testing.assert_allclose(optical_tidal(result['Riemann'], screen, k), 0, atol=4e-16)
    _symmetries(result)


def test_inhomogeneous_ordered_frame_jets_and_algebraic_bianchi():
    # A_i X are the exact smooth left-invariant S3 frame directions.
    A = np.array([[[0,0,0,1],[0,0,1,0],[0,-1,0,0],[-1,0,0,0]],
                  [[0,0,-1,0],[0,0,0,1],[1,0,0,0],[0,-1,0,0]],
                  [[0,1,0,0],[-1,0,0,0],[0,0,0,1],[0,0,-1,0]]], dtype=float)
    X = np.array([.2, -.3, .4, .5]); X /= np.linalg.norm(X)
    r = .43
    vector = np.einsum('iab,b->ia', A, X)
    P = X[0]*X[1]+.3*X[2]**2+.2*X[3]
    grad = np.array([X[1], X[0], .6*X[2], .2])
    hess = np.zeros((4,4)); hess[0,1] = hess[1,0] = 1; hess[2,2] = .6
    B = .04*np.array([[.2,.3,-.2,.1],[.3,.5,.2,-.4],[-.2,.2,-.1,.3],[.1,-.4,.3,.7]])
    g = np.diag([-1., 2., 1.5, .9])+P*B
    v, dg, d2g, dv, accg = _zeros()
    v[:] = .4*X[1]*B
    accg[:] = .2*B
    for i in range(3):
        dg[i] = r*np.dot(vector[i], grad)*B
        dv[i] = .4*r*vector[i,1]*B
        for j in range(3):
            d2g[i,j] = r*r*(vector[i]@hess@vector[j]+grad@(A[j]@A[i]@X))*B
    result = geometry_from_jets(g, v, dg, d2g, dv, accg, frame_scale=r)
    assert result['commutator_error'] < 1e-16
    _symmetries(result)
    # A coordinate-Hessian substitution breaks the frame jet identity.
    damaged = geometry_from_jets(g, v, dg, (d2g+d2g.swapaxes(0,1))/2, dv, accg, frame_scale=r)
    assert damaged['commutator_error'] > .001


def test_schwarzschild_pg_tangential_null_tidal_eigenvalues():
    M, radius = .37, 3.8
    x = np.array([radius, 0., 0.]); speed = np.sqrt(2*M/radius)
    g = np.eye(4); g[0,0] = -1+2*M/radius
    g[0,1:] = g[1:,0] = np.sqrt(2*M)*x/radius**1.5
    v, dg, d2g, dv, accg = _zeros()
    eye = np.eye(3); c = np.sqrt(2*M)
    for j in range(3):
        dg[j,0,0] = -2*M*x[j]/radius**3
        for i in range(3):
            dg[j,0,i+1] = dg[j,i+1,0] = c*(eye[i,j]/radius**1.5-1.5*x[i]*x[j]/radius**3.5)
        for k in range(3):
            d2g[j,k,0,0] = -2*M*eye[j,k]/radius**3+6*M*x[j]*x[k]/radius**5
            for i in range(3):
                value = c*(-1.5*(eye[i,j]*x[k]+eye[i,k]*x[j]+eye[j,k]*x[i])/radius**3.5
                           +5.25*x[i]*x[j]*x[k]/radius**5.5)
                d2g[j,k,0,i+1] = d2g[j,k,i+1,0] = value
    result = geometry_from_jets(g, v, dg, d2g, dv, accg, frame_scale=0)
    np.testing.assert_allclose(result['Ricci'], 0, atol=4e-17)
    k = np.array([1., -speed, 1., 0.])  # n + tangential unit spatial direction
    screen = np.array([[0.,1.,0.,0.], [0.,0.,0.,1.]])
    np.testing.assert_allclose(k@g@k, 0, atol=2e-16)
    np.testing.assert_allclose(screen@g@k, 0, atol=2e-16)
    np.testing.assert_allclose(screen@g@screen.T, np.eye(2), atol=2e-16)
    tidal = optical_tidal(result['Riemann'], screen, k)
    np.testing.assert_allclose(np.linalg.eigvalsh(tidal), np.array([-3.,3.])*M/radius**3, atol=4e-17)
    inverse = result['inverse']; R = result['Riemann']
    square = np.einsum('abcd,efgh,ae,bf,cg,dh', R,R,inverse,inverse,inverse,inverse)
    np.testing.assert_allclose(square, 48*M*M/radius**6, atol=2e-18)
    _symmetries(result)


def test_validation_and_tidal_screen_broadcast():
    g = np.diag([-1.,1.,1.,1.]); jets = _zeros()
    for scale in (-1, np.nan, True):
        with pytest.raises(ValueError, match='frame_scale'):
            geometry_from_jets(g, *jets, frame_scale=scale)
    with pytest.raises(ValueError, match='spatial metric'):
        geometry_from_jets(-g, *jets, frame_scale=0)
    wrong = np.zeros((3,3,4,4)); wrong[0,1,0,1] = 1
    with pytest.raises(ValueError, match='symmetric'):
        geometry_from_jets(g, jets[0], jets[1], wrong, jets[3], jets[4], frame_scale=1)
    matrix = optical_tidal(np.zeros((5,4,4,4,4)), np.eye(4)[1:3], np.array([1.,0.,0.,1.]))
    assert matrix.shape == (5,2,2)
    np.testing.assert_array_equal(matrix, 0)
