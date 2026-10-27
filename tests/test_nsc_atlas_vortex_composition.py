"""A direct nonlinear check of finite-resolution composition on a parent vortex."""
import math

import numpy as np
import pytest

from recursive_horizons.nsc_atlas_vortex_transfer import (
    KAPPA, bump, exact_theta_rate, jets, profile, reynolds,
)


def test_prepared_modes_compose_on_an_actual_axisymmetric_parent():
    modes, amplitudes = (4, 8, 16), (1., .25, .0625)
    theta = np.arange(65)*2*math.pi/65
    for r in (.75, 1.25):
        z = 0.
        a, ap, _ = profile(r)
        c, cp, _ = bump(z)
        # A smooth compact parent B=(0,r*a^2*chi^2,0) is divergence-free.
        B = np.array([0., r*a*a*c*c, 0.], dtype=float)
        Br = np.array([0., (a*a+2*r*a*ap)*c*c, 0.])
        Bz = np.array([0., 2*r*a*a*c*cp, 0.])
        U = np.repeat(B[:, None], len(theta), axis=1)
        Ur, Uz = [np.repeat(j[:, None], len(theta), axis=1) for j in (Br, Bz)]
        Ut = np.zeros_like(U)
        expected_stress = np.zeros((3, 3))
        for mode, amplitude in zip(modes, amplitudes):
            phase = np.exp(1j*(mode*theta+KAPPA*mode*math.log(r)))
            W, dr, dt, dz = jets(r, z, mode)
            for total, value in zip((U, Ur, Ut, Uz), (W, dr, dt, dz)):
                total += amplitude*np.real(value[:, None]*phase)
            expected_stress += amplitude**2*reynolds(r, z, mode)
        assert np.mean(U, axis=1) == pytest.approx(B, abs=5e-15)
        added_stress = U @ U.T/len(theta)-np.outer(B, B)
        assert added_stress == pytest.approx(expected_stress, abs=3e-14)
        # The pure-swirl parent's tangential convection is zero. Linear
        # viscous increments have zero angular mean; mean pressure has no
        # angular gradient in either evolution, so this is the NS rate difference.
        convective = U[0]*Ur[1]+U[1]*Ut[1]/r+U[2]*Uz[1]+U[0]*U[1]/r
        expected_rate = sum(e*e for e in amplitudes)*exact_theta_rate(r)
        assert -np.mean(convective) == pytest.approx(expected_rate, abs=3e-13)
