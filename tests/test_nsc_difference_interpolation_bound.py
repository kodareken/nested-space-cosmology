"""The difference transfer bounds the nonnegative quadratic stress error."""
from fractions import Fraction as Q

import numpy as np
import pytest
flint = pytest.importorskip("flint")
from flint import arb, ctx

from recursive_horizons.nsc_ks_ball_trajectory import restored_upper
from recursive_horizons.nsc_ks_finite_matter_error import finite_matter_error


def test_doubling_field_errors_costs_at_most_four_in_the_owned_stress_bound():
    rng = np.random.default_rng(2510)
    with ctx.workprec(160):
        for _ in range(20):
            f, z, e, ez, gamma = [arb(float(v)) for v in rng.uniform(0.001, 2.0, 5)]
            params = dict(mass=1.5, absolute_angular=2.3, axial_lower=0.8,
                          radius_lower=1.3, multiplicity=12, bits=160)
            old = finite_matter_error(f, z, e, ez, gamma, **params)
            new = finite_matter_error(f, z, 2*e, 2*ez, gamma, **params)
            for name in ("N", "beta"):
                assert restored_upper(new[name]) <= 4*restored_upper(old[name])
