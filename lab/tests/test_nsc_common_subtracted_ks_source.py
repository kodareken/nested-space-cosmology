"""Focused common-kernel signs, endpoint ordering and exact band accounting."""
from types import SimpleNamespace

import numpy as np
import pytest

from recursive_horizons.nsc_common_subtracted_ks_source import (
    SpectrumDeclaration, physical_source_kernel, reference_symbol_kernel,
    raw_ks_vertex_coefficients, symmetric_pairing_density, reference_band_remainder,
)
from recursive_horizons.nsc_nonlinear_ks_source import weak_KS_source
from recursive_horizons.nsc_lorentzian import geometry
from recursive_horizons.nsc_spatial_reference_symbol import (
    SymbolJet, SIGMA, supplied_KS_metric_jets, reference_projector,
)
from recursive_horizons.nsc_reference_band_action import band_frame


CONTROL = SpectrumDeclaration('control_only', 'analytic resolved-column sign fixture')


def source_fixture(z):
    # Three-source columns with nonzero within-source coherence and a real
    # spatial amplitude gradient. This is declared control data, not a UV sea.
    E = np.array([-1.1, -.3, .4, 1.2]); weights = np.array([.2, .4, .3, .1])
    k = np.array([-.9, -.25, .5, 1.1])
    F = np.array([[1., .2j, .3], [.1j, .8, -.2j]])
    C = np.broadcast_to(np.array([[.65, .08j, 0.], [-.08j, .3, .04], [0., .04, .45]]), (4, 3, 3))
    z = np.asarray(z)
    phase = np.exp(1j*z[..., None]*k)
    envelope = 1+.07*z
    psi = (envelope[..., None]*phase)[..., None, None]*F
    derivative = ((.07+1j*envelope[..., None]*k)*phase)[..., None, None]*F
    return E, weights, C, psi, derivative


def reference_fixture(batch, separation):
    k = np.array([-.8, .35, 1.3]); weights = np.array([.4, .6, .3])
    P = np.zeros((5, *batch, len(k), 2, 2), complex)
    P[0] = np.array([np.diag([1., 0.]), np.diag([0., 1.]), np.diag([0., 1.])])
    return reference_symbol_kernel(P, k, weights, separation=separation, declaration=CONTROL), P, k, weights


def test_zero_split_matches_existing_weak_source_and_reference_sign():
    x = np.array([-.5, 0., .5]); E, weights, C, psi, psi_z = source_fixture(x)
    source = physical_source_kernel(psi, psi, psi_z, psi_z, C, E, weights,
                                    separation=0., declaration=CONTROL)
    reference, P, k, kw = reference_fixture(x.shape, 0.)
    metric = np.array([[1.2, .1, 1.4, 1.1], [1.1, -.2, 1.5, 1.2], [1.3, .05, 1.6, 1.3]])
    shape = np.array([.4, .7, .3]); mass, angular = .3, .8
    vertices = raw_ks_vertex_coefficients(metric, mass, angular, envelopes=shape[:, None])
    pairing = symmetric_pairing_density(source, reference, vertices, vertices, compact_spatial_boundary=True)

    # Direct call to the inherited un-subtracted weak-vertex owner. Its
    # canonical frame is the identity for this synthetic field comparison.
    owner = SimpleNamespace(x=x, mass=mass, angular=angular, weights=np.array([.2, .3, .25]))
    sample = {'frame': np.ones((2, 3)), 'frame_t': np.zeros((2, 3)), 'KS': metric.T, 'shape': shape}
    phi = psi.transpose(2, 0, 1, 3).reshape(2*len(x), -1)
    phi_z = psi_z.transpose(2, 0, 1, 3).reshape(2*len(x), -1)
    weak, _, _ = weak_KS_source(owner, phi, phi_z, sample)
    blocks = weak.reshape(4, len(E), 3, len(E), 3)
    gaussian = -sum(weights[e]/(2*np.pi)*np.einsum('ij,bji->b', C[e], blocks[:, e, :, e, :]).real
                    for e in range(len(E)))
    a0 = np.sqrt(geometry(x)[0]**2-1)
    measure = owner.weights/a0
    reference_vertex = (vertices.multiplication[..., None, :, :]
                        +vertices.momentum[..., None, :, :]*k[None, None, :, None, None])
    reference_gradient = np.einsum('x,b->xb', measure, np.ones(4))*np.einsum(
        'xkij,xbkji,k->xb', P.sum(axis=0), reference_vertex, kw/(2*np.pi)).real
    assert np.max(abs(pairing.integrate(measure)-gaussian-reference_gradient.sum(axis=0))) < 3e-14
    assert np.max(abs(gaussian)) > 1e-3
    with pytest.raises(ValueError, match='complete source'):
        pairing.require_complete_trace(coincidence_convergence_evidence='not measured', channel_inventory_evidence='fixture')
    with pytest.raises(ValueError, match='noncompact'):
        symmetric_pairing_density(source, reference, vertices, vertices, compact_spatial_boundary=False)


def test_finite_split_uses_endpoint_coefficients_and_kernel_adjoint():
    z = np.array([-.2, .15]); eta = .13
    E, w, C, plus, plus_z = source_fixture(z+eta/2)
    _, _, _, minus, minus_z = source_fixture(z-eta/2)
    forward = physical_source_kernel(plus, minus, plus_z, minus_z, C, E, w, separation=eta, declaration=CONTROL)
    backward = physical_source_kernel(minus, plus, minus_z, plus_z, C, E, w, separation=-eta, declaration=CONTROL)
    assert np.max(abs(backward.kernel-forward.kernel.swapaxes(-1, -2).conj())) < 3e-15
    assert np.max(abs(backward.separation_derivative+forward.separation_derivative.swapaxes(-1, -2).conj())) < 3e-15
    reference, _, _, _ = reference_fixture(z.shape, eta)
    gp = np.broadcast_to([1.2, .1, 1.4, 1.1], (2, 4))
    gm = np.broadcast_to([1.1, -.1, 1.7, 1.3], (2, 4))
    vp = raw_ks_vertex_coefficients(gp, .3, .8, envelopes=np.array([.6, .4])[:, None])
    vm = raw_ks_vertex_coefficients(gm, .3, .8, envelopes=np.array([.2, .3])[:, None])
    result = symmetric_pairing_density(forward, reference, vp, vm, compact_spatial_boundary=True)
    D, Deta = forward.kernel-reference.kernel, forward.separation_derivative-reference.separation_derivative
    expected = np.array([[-np.trace(.5*(vp.multiplication[x, B]+vm.multiplication[x, B])@D[x]
                                   -.5j*(vp.momentum[x, B]+vm.momentum[x, B])@Deta[x]).real
                          for B in range(4)] for x in range(2)])
    assert np.max(abs(result.action_gradient_density-expected)) < 3e-15
    with pytest.raises(ValueError, match='same KS'):
        symmetric_pairing_density(backward, reference, vp, vm, compact_spatial_boundary=True)
    varying = np.array([eta, 2*eta])
    _, _, _, p, pz = source_fixture(z+varying/2)
    _, _, _, m, mz = source_fixture(z-varying/2)
    varying_source = physical_source_kernel(p, m, pz, mz, C, E, w, separation=varying, declaration=CONTROL)
    varying_reference, _, _, _ = reference_fixture(z.shape, varying)
    with pytest.raises(ValueError, match='constant separation'):
        symmetric_pairing_density(varying_source, varying_reference, vp, vm, compact_spatial_boundary=True)
    with pytest.raises(ValueError, match='resolved'):
        physical_source_kernel(plus[..., :2], minus[..., :2], plus_z[..., :2], minus_z[..., :2],
                               C, E, w, separation=eta, declaration=CONTROL)
    with pytest.raises(ValueError, match='convergence evidence'):
        SpectrumDeclaration('complete_spectrum', 'only a finite array')


def test_band_remainder_retains_exact_time_and_star_exchange():
    # One existing local reference case, with an exact constant spin-frame
    # rotation tangent. H'=RHR^dagger, U'=UR^dagger gives delta h=0,
    # while the unsummed local projector/vertex and connection terms need
    # not vanish. This tests the action identity without numerical reruns.
    P = reference_projector(*supplied_KS_metric_jets(.019, .23), [.2], [np.pi/2], [np.sqrt(5)])
    frame = band_frame(P, np.array([np.diag([1., 0.])]))
    A = SymbolJet.constant(.5j*SIGMA[0])
    delta_H = A*P['hamiltonian']-P['hamiltonian']*A
    delta_U = [-u*A for u in frame['U']]
    delta_effective = [h*0 for h in frame['effective']]
    remainder = reference_band_remainder(P, frame, delta_H, delta_U, delta_effective)
    assert np.max(remainder.decomposition_residual) < 3e-11
    assert np.max(remainder.band_identity_residual) < 3e-11
    assert np.max(abs(remainder.symmetric_vertex_orders)) > 1e-5
    assert np.max(abs(remainder.orders+remainder.symmetric_vertex_orders-remainder.action_derivative_orders)) < 3e-14
    eta, k, w = .11, np.array([.2]), np.array([.7])
    symmetric_pair = np.sum(np.exp(1j*k*eta)*w/(2*np.pi)*remainder.symmetric_vertex_orders.sum(axis=0)).real
    assert abs(remainder.fourier_pairing(k, w, separation=eta)+symmetric_pair) < 3e-14
