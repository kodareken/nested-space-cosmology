"""Frozen-slice operator controls, without PDE evolution or saved histories."""
import math

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.linalg import eigh
from threadpoolctl import threadpool_limits

import check_nsc_atlas_sphere_flux_completion as producer
from recursive_horizons.nsc_atlas_valve_spherical import SphericalValve


@pytest.fixture(scope="module")
def round_parent():
    with threadpool_limits(limits=1):
        system = SphericalValve(16, epsilon=0., reflection_symmetry=False)
    state = np.zeros_like(system.initial)
    state[0] = 1 / system.H
    return system, state


def test_round_spectrum_has_full_parity_and_the_physical_radius(round_parent):
    system, state = round_parent
    stiffness, gram, info = producer.matrices(system, state, 8, 97)
    eigenvalues = eigh(stiffness, gram, eigvals_only=True)
    degrees = np.arange(9)
    mass_squared = -2.5 * math.exp(-14 * system.sstar) + 3.125 * math.exp(-24 * system.sstar)
    expected = mass_squared + system.H**2 * degrees * (degrees + 2)
    np.testing.assert_allclose(eigenvalues, expected, rtol=3e-12, atol=2e-12)
    np.testing.assert_allclose(gram, np.eye(9) / system.H**3, rtol=3e-12, atol=2e-11)
    # In particular l=1,3 are present: cap reflection must not remove odd modes.
    assert eigenvalues[1] == pytest.approx(mass_squared + 3 * system.H**2, abs=2e-12)
    assert eigenvalues[3] == pytest.approx(mass_squared + 15 * system.H**2, abs=2e-12)
    assert info["constant_trial_Rayleigh"] == pytest.approx(mass_squared, abs=2e-12)


def test_off_shell_warped_metric_bilinear_agrees_with_direct_laplacian(round_parent):
    system, _ = round_parent
    u = system.u
    state = np.zeros_like(system.initial)
    state[0] = (1 + .07*u + .04*u*u) / system.H
    state[1] = .03 * (1 - .2*u) / system.H
    state[4] = math.sqrt(84) * (.025*u + .006*u*u)
    stiffness, gram, _ = producer.matrices(system, state, 5, 151)
    normalization = math.sqrt(2 * math.pi**2)
    phi_coeff = normalization * np.array([0., .3, 1., 0., 0., 0.])
    psi_coeff = normalization * np.array([-.2, 0., 0., 1., 0., 0.])

    def integrands(angle):
        # Direct warped-metric Laplacian in chi, independent of the stiffness
        # integration by parts and of the producer's spectral field jets.
        x = math.cos(angle)
        sine = math.sin(angle)
        b = (1 + .07*x + .04*x*x) / system.H
        bu = (.07 + .08*x) / system.H
        c = .03 * (1 - .2*x) / system.H
        cu = -.006 / system.H
        q = b + sine*sine*c
        qu = bu - 2*x*c + sine*sine*cu
        radius = sine*b
        radius_chi = x*b - sine*sine*bu
        q_chi = -sine*qu
        phi = 4*x*x - 1 + .6*x
        psi = 8*x*x*x - 4*x - .2
        psi_u = 24*x*x - 4
        psi_chi = -sine*psi_u
        psi_chichi = -x*psi_u + sine*sine*48*x
        laplacian = (psi_chichi + (2*radius_chi/radius - q_chi/q)*psi_chi) / q**2
        local_s = system.sstar + .025*x + .006*x*x
        mass = -2.5*math.exp(-14*local_s) + 3.125*math.exp(-24*local_s)
        measure = 4*math.pi*sine*sine*q*b*b
        return measure*phi*(-laplacian + mass*psi), measure*phi*psi

    direct_operator = quad(lambda x: integrands(x)[0], 0., math.pi, epsabs=2e-10, epsrel=2e-12)[0]
    direct_inner = quad(lambda x: integrands(x)[1], 0., math.pi, epsabs=2e-10, epsrel=2e-12)[0]
    assert phi_coeff @ stiffness @ psi_coeff == pytest.approx(direct_operator, rel=2e-11, abs=3e-9)
    assert phi_coeff @ gram @ psi_coeff == pytest.approx(direct_inner, rel=2e-11, abs=3e-9)
    # Smooth at the geometric poles does not mean zero derivative in u.
    assert 8*(-1.) + .6 != 0 and 8*(1.) + .6 != 0
    assert 24*(-1.)**2 - 4 != 0 and 24*(1.)**2 - 4 != 0
    pole_fields = system.pieces(state, [-1., 1.])
    np.testing.assert_allclose(pole_fields["q"], pole_fields["b"], rtol=0., atol=0.)


def test_negative_constant_trial_is_a_witness_not_a_state_selection(round_parent, monkeypatch):
    system, state = round_parent
    negative_mass_squared = -.1
    monkeypatch.setattr(producer, "shape_mass_squared", lambda s: np.full_like(s, negative_mass_squared))
    stiffness, gram, info = producer.matrices(system, state, 5, 73)
    eigenvalues = eigh(stiffness, gram, eigvals_only=True)
    assert stiffness[0, 0] / gram[0, 0] == pytest.approx(negative_mass_squared, abs=2e-14)
    assert info["constant_trial_Rayleigh"] == pytest.approx(negative_mass_squared, abs=2e-14)
    assert eigenvalues[0] == pytest.approx(negative_mass_squared, abs=2e-13)
    assert eigenvalues[1] > 0
