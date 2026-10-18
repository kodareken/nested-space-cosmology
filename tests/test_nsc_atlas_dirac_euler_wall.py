"""Smooth Lorentz curvature contraction, Euler transgression and wall identities."""
import itertools
import math

import pytest
import sympy as s

from recursive_horizons.nsc_atlas_dirac_trace import (
    dirac_euler_surface_trace, euler_chern_scalar, quiet_euler_surface,
)


def test_direct_lorentz_curvature_contraction_and_pure_trace_limit():
    t, b, tn, bn, A, B, k, delta = s.symbols('t b tn bn A B k delta')
    signature = (-1, 1, 1, 1)
    sectional = {(0, 1): A-t*b, (0, 2): A-t*b, (1, 2): B-b*b,
                 (0, 3): -(tn+t*t), (1, 3): -(bn+b*b), (2, 3): -(bn+b*b)}
    riemann = s.MutableDenseNDimArray.zeros(4, 4, 4, 4)
    for (i, j), value in sectional.items():
        value *= signature[i]*signature[j]
        riemann[i,j,i,j] = riemann[j,i,j,i] = value
        riemann[i,j,j,i] = riemann[j,i,i,j] = -value
    ricci = s.Matrix(4, 4, lambda i,j: sum(signature[m]*riemann[m,i,m,j] for m in range(4)))
    R = s.expand(sum(signature[i]*ricci[i,i] for i in range(4)))
    Ric2 = s.expand(sum(signature[i]*signature[j]*ricci[i,j]**2
                        for i,j in itertools.product(range(4), repeat=2)))
    Riem2 = s.expand(sum(math.prod(signature[i] for i in indices)*riemann[indices]**2
                         for indices in itertools.product(range(4), repeat=4)))
    Euler = s.expand(Riem2-4*Ric2+R*R)
    chern = -8*t*b*b+8*B*t+16*A*b
    singular = s.diff(Euler, tn)*tn+s.diff(Euler, bn)*bn
    assert s.expand(singular+s.diff(chern,t)*tn+s.diff(chern,b)*bn) == 0
    assert s.diff(Euler, tn, 2) == s.diff(Euler, bn, 2) == s.diff(Euler, tn, bn) == 0
    substitute = {tn: -k*delta, bn: -k*delta}
    assert s.expand(Riem2.subs(substitute)).coeff(delta, 2) == 12*k*k
    assert s.expand(Ric2.subs(substitute)).coeff(delta, 2) == 12*k*k
    assert s.expand((R*R).subs(substitute)).coeff(delta, 2) == 36*k*k
    Weyl2 = s.expand(Riem2-2*Ric2+R*R/3)
    assert s.expand(Weyl2.subs(substitute)).coeff(delta, 2) == 0
    assert s.expand(Weyl2.subs(substitute)).coeff(delta, 1) == 0
    assert s.expand(Euler.subs(substitute)).coeff(delta, 2) == 0


def test_euler_flux_sign_from_flat_ball_and_static_schwarzschild():
    radius, mass = s.symbols('r m', positive=True)
    f = 1-2*mass/radius
    t, b = mass/(radius**2*s.sqrt(f)), s.sqrt(f)/radius
    chern = -8*t*b*b+8*t/radius**2
    flux = s.simplify(s.sqrt(f)*radius**2*chern)
    assert s.simplify(-s.diff(flux,radius)/radius**2-48*mass**2/radius**6) == 0
    # Constant-curvature three-sphere boundary of a flat Euclidean ball.
    ball = euler_chern_scalar(Ktau=.5, Ksphere=.5,
                              intrinsic_acceleration=.25, intrinsic_curvature=.25)
    assert ball*2*math.pi**2*2**3 == pytest.approx(32*math.pi**2)
    numeric = euler_chern_scalar(Ktau=float(t.subs({radius:3.,mass:.4})),
        Ksphere=float(b.subs({radius:3.,mass:.4})), intrinsic_acceleration=0.,
        intrinsic_curvature=1/9)
    assert numeric == pytest.approx(float(chern.subs({radius:3.,mass:.4})))


def test_actual_signed_wall_identity_removes_mass_and_radius():
    r, m, k, H = s.symbols('r m k H', nonzero=True)
    bm = m/(k*r*r)+(k*k-H*H)*r/(2*k)
    bp = bm-k*r
    acceleration = bm*s.diff(bm,r)+H*H*r
    A = acceleration/r
    B = (bm*bm+H*H*r*r)/r**2
    def chern(beta):
        return -8*s.diff(beta,r)*(beta/r)**2+8*B*s.diff(beta,r)+16*A*beta/r
    surface = s.factor(chern(bm)-chern(bp))
    assert s.factor(surface-4*k*(k*k+3*H*H)) == 0
    assert s.factor(chern(-bp)-chern(-bm)-surface) == 0
    assert quiet_euler_surface(.5)['euler_surface'] == pytest.approx(56., abs=2e-12)
    coefficient = dirac_euler_surface_trace(56.)
    assert coefficient == pytest.approx(-77/(720*math.pi**2))
    assert dirac_euler_surface_trace(56./7**3) == pytest.approx(coefficient/7**3)
    assert dirac_euler_surface_trace(0.) == 0.


def test_nonfinite_or_unadmitted_wall_data_are_rejected():
    with pytest.raises(ValueError): dirac_euler_surface_trace(float('nan'))
    with pytest.raises(ValueError): quiet_euler_surface(.7)
    with pytest.raises(ValueError):
        euler_chern_scalar(Ktau=1., Ksphere=1., intrinsic_acceleration=0.,
                           intrinsic_curvature=float('inf'))
