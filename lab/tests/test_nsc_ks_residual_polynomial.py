"""Independent constant-operator residual and continuous-cell enclosure."""
import pytest
pytest.importorskip('flint')
import numpy as np
from flint import arb, acb, ctx

from recursive_horizons.nsc_ks_trajectory import TrajectorySegment
from recursive_horizons.nsc_ks_ball_trajectory import BallFourierSegment
from recursive_horizons.nsc_ks_residual_polynomial import ResidualPolynomial


def case():
    start = np.r_[np.array([.5, .25]), np.zeros(16)]
    end = 2*start
    segment = TrajectorySegment(1.03, 1.02, start, end, np.zeros((6, len(start))))
    field = BallFourierSegment(segment, [1.], 8, 1., bits=120)
    residual = ResidualPolynomial(field, [1], [1], [arb(1)/2], {}, .4, .8, [.7])
    return field, residual


def test_constant_operator_has_the_correct_two_spin_rows():
    field, residual = case()
    with ctx.workprec(120):
        coefficients = residual.point_coefficients(arb(1)/3, {}, 0)[0]
        h = arb(field.rho_end)-arb(field.rho_start)
        for q, row in enumerate(coefficients):
            factor = 1+arb(q)/residual.degree
            X0, X1 = factor/2, factor/4
            E, m, ell = arb(.7), arb(.4), arb(.8)
            expected0 = arb(1)/2-h*(-acb(0,E)*X0+(ell/2-acb(0,m))*X1)
            expected1 = arb(1)/4-h*(acb(0,E)*X1+(-ell/2-acb(0,m))*X0)
            assert row[0][0].overlaps(expected0)
            assert row[1][0].overlaps(expected1)


def test_cell_bound_covers_the_continuous_bernstein_residual():
    _, residual = case()
    bounds = residual.spatial_cell_bounds(arb(1)/2, arb(1)/4, {}, {}, taylor_order=3)
    assert bounds[0] > arb(.5)
    assert bounds[1] == 0


def test_actual_axial_wave_and_radius_potential_both_enter_the_operator():
    harmonic = np.array([1, 1j, -1, -1j]*2, complex)
    start = np.r_[np.zeros(2), np.stack((harmonic, 2*harmonic)).ravel()]
    segment = TrajectorySegment(1.03, 1.02, start, 2*start, np.zeros((6, len(start))))
    field = BallFourierSegment(segment, [1.], 8, 1., bits=120)
    residual = ResidualPolynomial(field, [1], [1], [arb(1)/2], {(1, 0): [arb(.3)]}, .4, .8, [.7])
    with ctx.workprec(120):
        z = arb(1)/16
        values = residual.point_coefficients(z, {(1, 0): [arb(1), arb(0)]}, 1)
        wave, h = 4*arb.pi(), arb(field.rho_end)-arb(field.rho_start)
        phase = acb(0, wave*z).exp()
        for q in range(residual.degree+1):
            X0 = (1+arb(q)/residual.degree)*phase
            X1 = 2*X0
            diagonal = acb(0, wave-arb(.7))
            expected = phase-h*(diagonal*X0+(arb(.8)/2+arb(.3)-acb(0,arb(.4)))*X1)
            assert values[0][q][0][0].overlaps(expected)
            assert values[1][q][0][0].overlaps(acb(0,wave)*expected)
