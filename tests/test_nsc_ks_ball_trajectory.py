"""Ball containment of the exact stored Fourier/polynomial reconstruction."""
import math

import numpy as np
import pytest

flint = pytest.importorskip('flint')
from flint import arb, acb, ctx
from recursive_horizons.nsc_ks_ball_trajectory import (
    BallFourierSegment, power_to_bernstein, polynomial_at, affine_restriction,
    exact_upper, restored_upper,
)
from recursive_horizons.nsc_ks_trajectory import TrajectorySegment


def test_exact_input_transform_and_derivative_containment():
    n = 8
    # A single Fourier harmonic with dyadic samples, so no trigonometric
    # input rounding is hidden in the control.
    harmonic = np.array([1, 1j, -1, -1j]*2, complex)
    start = np.r_[np.array([.5, .25]), np.stack((harmonic, 2*harmonic)).ravel()]
    end = 2*start
    segment = TrajectorySegment(1.03, 1.02, start, end, np.zeros((6, len(start)), complex))
    before = ctx.prec
    series = BallFourierSegment(segment, [1.], n, arb(1), bits=120)
    assert ctx.prec == before
    with ctx.workprec(120):
        value = series.value(arb(1)/2, arb(1)/16)
        phase = acb(0, arb.pi()/4).exp()
        assert value[0][0].overlaps(arb(3)/2*(arb(1)/2+phase))
        assert value[1][0].overlaps(arb(3)/2*(arb(1)/4+2*phase))
        derivative = series.value(arb(1)/2, arb(1)/16, axial_derivative=1)
        assert derivative[0][0].overlaps(arb(3)/2*acb(0, 4*arb.pi())*phase)


def test_bernstein_conversion_and_restriction_preserve_the_polynomial():
    with ctx.workprec(120):
        coefficients = list(map(acb, [1, 2, -3, 4]))
        restricted = affine_restriction(coefficients, arb(1)/4, arb(3)/4)
        bernstein = power_to_bernstein(restricted, 6)
        for point in (arb(0), arb(1)/3, arb(1)):
            direct = polynomial_at(coefficients, arb(1)/4+point/2)
            observed = sum((b*math.comb(6, k)*point**k*(1-point)**(6-k) for k, b in enumerate(bernstein)), acb(0))
            assert direct.overlaps(observed)
        bound = arb(2).sqrt().upper()
        assert restored_upper(exact_upper(bound)) == bound
