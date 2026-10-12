"""Preserve tight Chebyshev bounds and include the cubic U history term."""
from fractions import Fraction as Q

import numpy as np
import pytest

from recursive_horizons.nsc_local_incoming_family import LocalAxialFunction, LocalIncomingFamily
from recursive_horizons.nsc_ks_current_history_bounds import (
    chebyshev_profile_bounds, radius_bounds, rational_record, value_integral_bounds)
from recursive_horizons.nsc_ks_radius_enclosure import axial_profile_bounds


def test_high_degree_bound_retains_chebyshev_basis():
    coefficients = [0.0]*33
    coefficients[-1] = 1e-4
    profile = LocalAxialFunction(tuple(coefficients), 0.0)
    result = chebyshev_profile_bounds(profile)
    assert float(result["profile_bounds"][0]) < 0.000100001
    assert axial_profile_bounds(profile)[0] > result["profile_bounds"][0] * 1e8
    for order, upper in enumerate(result["profile_bounds"]):
        samples = [abs(profile(float(z), order)) for z in np.linspace(-0.06, 0.06, 501)]
        assert max(samples) <= float(upper)


def test_constant_u_retains_its_axial_cutoff_and_positive_radius_margin():
    coefficients = np.zeros((2, 8))
    coefficients[1, 0] = 12
    family = LocalIncomingFamily(coefficients)
    bound = radius_bounds(family, 1.03)
    assert bound["radius_lower"] < Q(7, 5)
    assert bound["radius_lower"] > Q(139, 100)
    assert bound["delta_radius_bounds"][1] > 0  # U cutoff derivative matters.
    integrals = value_integral_bounds(bound, 1.5, 2.0)
    assert integrals["K1"] > 0
    assert integrals["history_tangent_bound"] is None
    assert integrals["D"] == (Q(1.03)-1)/Q(4, 5)**2


def test_unproved_slab_and_nonpositive_radius_are_rejected():
    family = LocalIncomingFamily(np.zeros((2, 8)))
    with pytest.raises(ValueError, match="background slab"):
        radius_bounds(family, 1.1)
    coefficients = np.zeros((2, 8))
    coefficients[0, 0] = 100
    with pytest.raises(ValueError, match="positive radius"):
        radius_bounds(LocalIncomingFamily(coefficients), 1.03)


def test_serialization_rounds_inspection_upper_outward():
    value = Q(1, 10)
    packed = rational_record({"value": value})
    assert Q(packed["value"]["exact_rational"]) == value
    assert Q(packed["value"]["binary64_upper"]) >= value
