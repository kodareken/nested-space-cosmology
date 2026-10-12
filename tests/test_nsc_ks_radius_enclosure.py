"""Actual cutoff bounds and differentiated geometric-series remainder."""
from fractions import Fraction as Q

import mpmath as mp
import numpy as np

from recursive_horizons.nsc_ks_radius_enclosure import axial_profile_bounds, reciprocal_radius_tail
from recursive_horizons.nsc_local_incoming_family import LocalAxialFunction


def test_axial_bound_covers_both_transitions_and_representation_map():
    profile = LocalAxialFunction((.001, -.002, .0001), 1.2068780671724793)
    bounds = axial_profile_bounds(profile)
    for z in np.linspace(profile.center-profile.outer, profile.center+profile.outer, 701):
        assert all(abs(profile(z, j)) <= float(bounds[j]) for j in range(3))


def test_remainder_covers_two_spatial_derivatives_without_radius_linearization():
    with mp.workdps(60):
        # x(z)=eta0+eta1 sin(z); all three absolute bounds are explicit.
        eta0, eta1 = mp.mpf('0.00002'), mp.mpf('0.00001')
        w = (Q(3, 100000), Q(1, 100000), Q(1, 100000))
        result = reciprocal_radius_tail(w, (0, 0, 0), normal_support=1,
            reference_radius_lower=1, axial_lower=1, absolute_angular=1, order=4)
        def remainder(z):
            x = eta0+eta1*mp.sin(z)
            return -x/(1+x)-sum((-x)**n for n in range(1, 5))
        for z in (mp.mpf('0'), mp.mpf('.3'), mp.mpf('1.7'), mp.mpf('4.2')):
            for j, bound in enumerate(result['potential_derivative_tail_bounds']):
                assert abs(mp.diff(remainder, z, j)) <= mp.mpf(bound.numerator)/bound.denominator
