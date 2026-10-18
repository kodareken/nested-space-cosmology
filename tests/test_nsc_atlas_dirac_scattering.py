"""Independent characteristic, action-variation and CAR controls; no P001 grids."""
import numpy as np
import pytest

from recursive_horizons.nsc_atlas_dirac_scattering import (
    batch_null_transfer, inherit_boundary_state, null_transfer, occupation_response,
    ordered_metric_response, reset_inherited_mode)
from recursive_horizons.nsc_boundary_state import GaussianBoundaryState


def test_characteristic_flux_zero_and_constant_plane_wave():
    zero,_ = null_transfer(np.zeros((8,8)),.125,.125)
    assert np.array_equal(zero,np.eye(16))
    errors = []
    for n in (8,16,32):
        h = 1/n
        labels = (np.arange(n)+.5)*h
        w,p = 1.7,1.1
        q = w*w/(4*p)
        source = np.r_[np.exp(-1j*p*labels),2*q/w*np.exp(-1j*q*labels)]*np.sqrt(h)
        target = np.r_[np.exp(-1j*(p*labels+q)),
                       2*q/w*np.exp(-1j*(p+q*labels))]*np.sqrt(h)
        u,_ = null_transfer(np.full((n,n),w),h,h)
        assert np.linalg.norm(u.conj().T@u-np.eye(2*n),ord=2) < 4e-14
        errors.append(np.linalg.norm(u@source-target))
    assert 3.7 < errors[0]/errors[1] < 4.3
    assert 3.7 < errors[1]/errors[2] < 4.3


def test_exact_static_solution_across_C1_wall():
    # H xi=0, xi=(exp(-integral W dx),0). W and W_x continuous
    # at x=0, W_xx jumps, as in the actual radar wall potential.
    def amplitude(x):
        return np.exp(-.7*x-3.2*np.maximum(x,0)**3/3)
    errors = []
    for n in (8,16,32):
        h = 1/n
        labels = (np.arange(n)+.5)*h
        x = (labels[None,:]-labels[:,None])/2
        w = .7+3.2*np.maximum(x,0)**2
        incoming = np.r_[amplitude(-labels/2),-1j*amplitude(labels/2)]*np.sqrt(h/2)
        target = np.r_[amplitude((1-labels)/2),
                       -1j*amplitude((labels-1)/2)]*np.sqrt(h/2)
        u,_ = null_transfer(w,h,h)
        errors.append(np.linalg.norm(u@incoming-target))
    assert 3.6 < errors[0]/errors[1] < 4.5
    assert 3.6 < errors[1]/errors[2] < 4.5


def test_metric_source_derivative_and_causality():
    n,h = 8,.125
    w = 1+np.arange(n*n).reshape(n,n)/300
    profiles = np.zeros((2,n,n))
    profiles[0,6,6] = 1.
    profiles[1,2,2] = 1.
    u,j = null_transfer(w,h,h,profiles=profiles)
    response = ordered_metric_response(j,profiles)
    eps = 2e-5
    up,jp = null_transfer(w*np.exp(eps*profiles[1]),h,h,profiles=profiles)
    um,jm = null_transfer(w*np.exp(-eps*profiles[1]),h,h,profiles=profiles)
    # Independent derivative of the full scattering matrix fixes the action
    # vertex sign and normalization; the later insertion checks retardedness.
    assert np.max(np.abs(1j*u.conj().T@(up-um)/(2*eps)-j[1])) < 4e-11
    assert np.max(np.abs((jp[0]-jm[0])/(2*eps)-response)) < 4e-11
    assert np.linalg.norm(response) > 1e-7
    _,future = null_transfer(w*np.exp(eps*profiles[0]),h,h,profiles=profiles)
    assert np.array_equal(future[1],j[1])
    # Spacelike point insertions commute, including their common past.
    profiles[0] = 0.
    profiles[1] = 0.
    profiles[0,2,6],profiles[1,6,2] = 1.,1.
    _,js = null_transfer(w,h,h,profiles=profiles)
    assert np.max(np.abs(js[0]@js[1]-js[1]@js[0])) < 2e-17
    with pytest.raises(ValueError):
        ordered_metric_response(js,profiles)


def test_angular_sign_units_and_parallel_batch():
    w = np.arange(16).reshape(4,4)/20+.7
    profiles = np.zeros((2,4,4))
    profiles[0,3,3],profiles[1,0,0] = 1.,1.
    u,j = null_transfer(w,.2,.3,profiles=profiles)
    scaled,scaled_j = null_transfer(w/7,1.4,2.1,profiles=profiles)
    assert np.allclose(scaled,u,rtol=0,atol=4e-16)
    assert np.allclose(scaled_j,j,rtol=0,atol=4e-16)
    histories,vertices = batch_null_transfer(np.array([w,-w]),.2,.3,profiles=profiles)
    sign = np.diag(np.r_[np.ones(4),-np.ones(4)])
    assert np.allclose(histories[0],u,rtol=0,atol=4e-16)
    assert np.allclose(histories[1],sign@u@sign,rtol=0,atol=4e-16)
    assert np.allclose(vertices[1],sign@j@sign,rtol=0,atol=4e-16)


def test_inherited_CAR_complement_and_single_multiplet_factor():
    rng = np.random.default_rng(4026)
    q,_ = np.linalg.qr(rng.normal(size=(6,6))+1j*rng.normal(size=(6,6)))
    c = q@np.diag(np.linspace(0,1,6))@q.conj().T
    f = rng.normal(size=6)+1j*rng.normal(size=6)
    f /= np.linalg.norm(f)
    empty = reset_inherited_mode(c,f,0.)
    occupied = reset_inherited_mode(c,f,.4)
    assert np.min(np.linalg.eigvalsh(occupied)) > -1e-14
    assert np.max(np.linalg.eigvalsh(occupied)) < 1+1e-14
    assert np.allclose(occupied-empty,.4*np.outer(f,f.conj()),rtol=0,atol=2e-16)
    matrix = rng.normal(size=(6,6))+1j*rng.normal(size=(6,6))
    response = matrix+matrix.conj().T
    direct = 6*np.trace((occupied-empty)@response).real
    assert occupation_response(response,f,angular_k=3,occupation=.4) == pytest.approx(direct,abs=5e-15)
    for bad in (-.1,1.01,np.nan,1j,True):
        with pytest.raises(ValueError):
            reset_inherited_mode(c,f,bad)


def test_domains():
    for potential in (np.ones((2,2))*1j,np.full((2,2),np.nan),np.ones(4)):
        with pytest.raises(ValueError):
            null_transfer(potential,.1,.1)
    with pytest.raises(ValueError):
        null_transfer(np.ones((2,2)),0.,.1)
    with pytest.raises(ValueError):
        null_transfer(np.ones((2,2)),.1,.1,profiles=np.ones((2,3,3)))


def test_parent_child_inheritance_retains_initial_correlations():
    rng = np.random.default_rng(4026)
    q,_ = np.linalg.qr(rng.normal(size=(6,6))+1j*rng.normal(size=(6,6)))
    covariance = q@np.diag(np.linspace(0,1,6))@q.conj().T
    incoming = GaussianBoundaryState(covariance,np.r_[np.zeros(3),np.ones(3)])
    u,_ = null_transfer(np.ones((3,3)),.2,.2)
    inherited = inherit_boundary_state(incoming,u,
        outgoing_parent_mask=np.r_[np.ones(3),np.zeros(3)])
    # The existing four-term Gaussian elimination is an independent block
    # reconstruction of the outgoing cap covariance from the actual transfer.
    terms = incoming.reconstruct(u[3:,3:],u[3:,:3])
    assert np.max(np.abs(terms['total']-inherited.blocks[-1])) < 3e-16
    assert np.linalg.norm(terms['initial_parent_child']) > 1e-3
    assert np.linalg.norm(terms['initial_child_parent']) > 1e-3
    assert np.linalg.eigvalsh(inherited.covariance) == pytest.approx(
        np.linalg.eigvalsh(covariance),abs=8e-16)
