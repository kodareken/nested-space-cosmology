"""Actual-ODE metric jets checked against analytic uniform and Fourier flows."""
import numpy as np
import pytest

from measure_nsc_conformal_geometric_average import snapshot_jets
from recursive_horizons.nsc_spherical_feedback_action import partial_F, feedback_Z
from recursive_horizons.nsc_spherical_galerkin_coupling import periodic_interpolation


def test_saved_rate_direction_recovers_uniform_exponential_metric_jets():
    A, C_W, radial_H, Q_H = .3, .4, -.1, .2
    radius, Q = 2.1, .7
    Fr, Fchi = partial_F(radius, A, C_W)
    Z = float(feedback_Z(A))
    pc = 2*float(Fchi)*Q_H
    pr = float(Fr)/float(Fchi)*pc+2*Z*radial_H*radius
    interpolation = periodic_interpolation(7, 32, 8.)

    def coarse(value):
        return np.full((1, 7), value)

    def fine(value):
        return np.full((1, 32), value)

    arrays = dict(r=coarse(radius), Q=coarse(Q), p_r=coarse(pr), p_chi=coarse(pc),
                  rate_r=fine(radial_H*radius), rate_Q=fine(Q_H*Q),
                  rate_p_r=fine(radial_H*pr), rate_p_chi=fine(0.),
                  L=fine(Q), L_dot=fine(Q_H*Q))
    jets, gaps = snapshot_jets(interpolation, arrays, 0, A=A, C_W=C_W)
    assert max(gaps.values()) < 2e-14
    assert jets['radius_second_rate'] == pytest.approx(np.full(32, radial_H**2*radius))
    assert jets['radial_scale_second_rate'] == pytest.approx(
        np.full(32, (radial_H+Q_H)**2*radius*Q))


def test_second_Q_direction_keeps_the_outer_galerkin_projection():
    # Retain modes 0,+/-1. Multiplication creates cos(2x), which is projected out.
    A, C_W, radius, Q0, amplitude, b = .3, .4, 2., .7, .2, .15
    Fr, Fchi = partial_F(radius, A, C_W)
    ng, nq = 3, 32
    interpolation = periodic_interpolation(ng, nq, 2*np.pi)
    coarse_x, fine_x = 2*np.pi*np.arange(ng)/ng, 2*np.pi*np.arange(nq)/nq
    pc = 2*float(Fchi)*b*np.cos(coarse_x)
    actual_Qt = Q0*b*(np.cos(fine_x)+amplitude/2)
    arrays = {
        'r': np.full((1, ng), radius), 'Q': (Q0*(1+amplitude*np.cos(coarse_x)))[None, :],
        'p_chi': pc[None, :], 'p_r': (float(Fr)/float(Fchi)*pc)[None, :],
        'rate_r': np.zeros((1, nq)), 'rate_Q': actual_Qt[None, :],
        'rate_p_r': np.zeros((1, nq)), 'rate_p_chi': np.zeros((1, nq)),
        'L': (Q0*(1+amplitude*np.cos(fine_x)))[None, :], 'L_dot': actual_Qt[None, :],
    }
    jets, gaps = snapshot_jets(interpolation, arrays, 0, A=A, C_W=C_W)
    expected_Qtt = Q0*b*b*(.5+amplitude/2*np.cos(fine_x))
    assert max(gaps.values()) < 2e-14
    assert jets['radial_scale_second_rate']/radius == pytest.approx(expected_Qtt, abs=3e-15)
    assert np.max(abs(jets['radius_second_rate'])) < 2e-14
