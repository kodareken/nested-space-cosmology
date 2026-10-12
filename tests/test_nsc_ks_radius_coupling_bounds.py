from fractions import Fraction as Q

import numpy as np
import pytest
import sympy as sp

from recursive_horizons.nsc_ks_radius_coupling_bounds import (
    radius_coupling_integrals, history_coupling_bounds, with_radius_coupling_bounds)
from recursive_horizons.nsc_ks_current_field_cone import whole_cone_continuous_inputs
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily


def test_coupling_derivative_and_integral_are_exact():
    # Independent symbolic differentiation detects a lost reference factor or
    # accidental differentiation of the homogeneous reference radius.
    r0, d, dz, a, ell = sp.symbols("r0 d dz a ell", positive=True)
    m = ell*d/(a*r0*(r0+d))
    assert sp.simplify(sp.diff(m, d)*dz-ell*dz/(a*(r0+d)**2)) == 0
    s = sp.symbols("s", nonnegative=True)
    integral = sp.integrate(s*3+s**3*5/6, (s, 0, sp.Rational(1, 10)))
    value = radius_coupling_integrals([3, 7], [5, 11], Q(1, 10), 2, 1, 4)
    assert value["M_integral"] == 2*Q(integral)
    assert value["Mz_integral"] == 4*(Q(7, 200)+Q(11, 240000))


def test_zero_history_and_zero_angular_channel_have_exact_zero_coupling():
    assert set(radius_coupling_integrals([0, 0], [0, 0], Q(1, 50), 1, 1, 3).values()) == {0}
    assert set(radius_coupling_integrals([1, 2], [3, 4], Q(1, 50), 1, 1, 0).values()) == {0}
    with pytest.raises(ValueError):
        radius_coupling_integrals([1, 2], [3, 4], Q(1, 50), 1, 0, 3)


def test_inputs_bind_geometry_and_keep_reference_and_source_gaps():
    coefficients = np.zeros((2, 8))
    coefficients[:, :2] = [[.01, .002], [.03, -.001]]
    family = LocalIncomingFamily(coefficients)
    args = dict(mass=.1, angular=2., rho_up=1.03, source_energies=[1.])
    previous = whole_cone_continuous_inputs(family, **args)
    filled = with_radius_coupling_bounds(previous, family, 1.03, 2.)
    assert previous["propagation"]["M_integral"]["value"] is None
    assert all(v > 0 for v in history_coupling_bounds(family, 1.03, 2.).values())
    for name in ("M_integral", "Mz_integral"):
        assert filled["propagation"][name]["value"] is not None
    for name in ("reference_initial", "reference_residual"):
        assert filled["propagation"][name]["value"] is None
    assert not filled["complete_for_propagate_difference_error"]
    assert filled["physical_rho1_source_error"] is None
    with pytest.raises(ValueError, match="geometry binding"):
        with_radius_coupling_bounds(previous, family, 1.02, 2.)
    with pytest.raises(ValueError, match="angular channel"):
        with_radius_coupling_bounds(previous, family, 1.03, 3.)
