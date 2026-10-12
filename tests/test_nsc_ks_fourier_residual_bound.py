"""Global convolution bound includes actual mode shifts and continuous tails."""
import pytest
pytest.importorskip('flint')
from flint import arb,acb,ctx
from recursive_horizons.nsc_ks_fourier_residual_bound import polynomial_residual_bounds,profile_sup_bounds,operator_remainder_bounds
from test_nsc_ks_residual_polynomial import case
import numpy as np
from recursive_horizons.nsc_ks_trajectory import TrajectorySegment
from recursive_horizons.nsc_ks_ball_trajectory import BallFourierSegment
from recursive_horizons.nsc_ks_residual_polynomial import ResidualPolynomial,row_norm


def test_global_bound_covers_the_exact_constant_operator_control():
    field,residual=case()
    with ctx.workprec(120):
        bound=polynomial_residual_bounds(residual,{})
        assert bound[0]>arb(.5) and bound[1]==0
        leftover=operator_remainder_bounds(field,{'inv_a2':arb(0),'inv_a':arb(0),'inv_ar':arb(0)},
            {},(arb(0),arb(0)),.4,.8,[.7])
        assert leftover==(arb(0),arb(0))


def test_profile_norm_includes_finite_modes_and_omitted_band():
    with ctx.workprec(120):
        profile={(1,0):{'coefficients':[acb(1),acb(2),acb(1)],'tail':(arb(1)/10,arb(1)/20)}}
        norm=profile_sup_bounds(profile,1)[(1,0)]
        assert norm[0]>4 and norm[1]>4*arb.pi()


def test_nonzero_profile_convolution_bounds_both_mode_shifts_and_derivatives():
    harmonic=np.array([1,1j,-1,-1j]*2,complex)
    start=np.r_[np.zeros(2),np.stack((harmonic,2*harmonic)).ravel()]
    segment=TrajectorySegment(1.03,1.02,start,2*start,np.zeros((6,len(start))))
    field=BallFourierSegment(segment,[1.],8,1.,bits=120)
    residual=ResidualPolynomial(field,[1],[1],[arb(1)/2],{(1,0):[arb(.3)]},.4,.8,[.7])
    with ctx.workprec(120):
        profiles={(1,0):{'coefficients':[acb(arb(1)/2),acb(0),acb(arb(1)/2)],'tail':(arb(0),arb(0))}}
        bound=polynomial_residual_bounds(residual,profiles)
        for j in range(17):
            z=arb(j)/17
            jets={(1,0):[(2*arb.pi()*z).cos(),-2*arb.pi()*(2*arb.pi()*z).sin()]}
            point=residual.point_coefficients(z,jets,1)
            for derivative in range(2):
                for coefficient in point[derivative]:
                    for spin in range(2):
                        assert bound[derivative]>=row_norm(coefficient[spin])
