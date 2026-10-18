"""Independent Clifford/heat contraction, central-charge and scale controls."""
import math

import pytest
import sympy as s

from recursive_horizons.nsc_atlas_dirac_trace import dirac_bulk_trace


def test_spin_connection_trace_and_fermion_heat_coefficients():
    sigma = [s.Matrix([[0, 1], [1, 0]]), s.Matrix([[0, -s.I], [s.I, 0]]), s.diag(1, -1)]
    gamma = [s.kronecker_product(sigma[0], p) for p in sigma]
    gamma.append(s.kronecker_product(sigma[1], s.eye(2)))
    I = s.eye(4)
    for a in range(4):
        for b in range(4):
            assert gamma[a]*gamma[b]+gamma[b]*gamma[a] == 2*int(a==b)*I
    # Constant sectional curvature K=1: Riemann^2=24. Direct spin connection.
    spin_trace = sum(s.trace((gamma[a]*gamma[b]/2)**2)
                     for a in range(4) for b in range(4) if a != b)
    assert spin_trace == -s.Rational(1, 2)*24
    R, Ric, Riem = s.symbols('R Ric Riem')
    trE, trE2, trOmega2 = -R, R*R/4, -Riem/2
    heat = s.expand(4*(5*R*R-2*Ric+2*Riem)+60*R*trE+180*trE2+30*trOmega2)
    assert s.expand(heat-(5*R*R-8*Ric-7*Riem)) == 0
    C2, Euler = Riem-2*Ric+R*R/3, Riem-4*Ric+R*R
    # Fermionic sign: -a4, not the bosonic determinant.
    assert s.expand(-heat/360-(C2/s.Integer(20)-s.Rational(11, 360)*Euler)) == 0


def test_actual_quiet_cap_and_exterior_from_curvature_not_flat_kernel():
    cap = dirac_bulk_trace(R=12., Ricci_squared=36., Riemann_squared=24.)
    assert cap == pytest.approx(-11/(240*math.pi**2))
    r, M = .45, .4
    ext = dirac_bulk_trace(R=0., Ricci_squared=0., Riemann_squared=48*M*M/r**6)
    assert ext == pytest.approx(7*M*M/(120*math.pi**2*r**6))
    assert cap < 0 < ext
    assert dirac_bulk_trace(R=0., Ricci_squared=0., Riemann_squared=0.) == 0.


def test_dimensional_curvature_rescaling():
    base = dirac_bulk_trace(R=12., Ricci_squared=36., Riemann_squared=24.)
    scaled = dirac_bulk_trace(R=12./49, Ricci_squared=36./49**2, Riemann_squared=24./49**2)
    assert scaled == pytest.approx(base/7**4)


def test_reject_nonfinite_bulk_invariants():
    with pytest.raises(ValueError):
        dirac_bulk_trace(R=math.nan, Ricci_squared=0., Riemann_squared=0.)
