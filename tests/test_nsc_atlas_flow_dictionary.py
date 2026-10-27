"""Independent geometric contractions for the atlas's published flow blocks."""
import math

import numpy as np
import pytest
import sympy as sp

from recursive_horizons.nsc_atlas_nested_slice import tension_wall


def test_pg_null_directions_and_invariant_mass():
    v, r = sp.symbols("v r", real=True, nonzero=True)
    g = sp.Matrix([[v*v-1, -v], [-v, 1]])
    kp, km = sp.Matrix([1, v+1]), sp.Matrix([1, v-1])
    assert g.det() == -1
    assert sp.expand((kp.T*g*kp)[0]) == 0
    assert sp.expand((km.T*g*km)[0]) == 0
    assert sp.expand((kp.T*g*km)[0]) == -2
    mass = sp.simplify(r*(1-g.inv()[1, 1])/2)
    assert mass == r*v*v/2
    # Opposite flow directions share the mass, but reverse the trapping signs.
    for flow, signs in [(-2, (-1, -1)), (0, (1, -1)), (2, (1, 1))]:
        assert tuple(np.sign([flow+1, flow-1])) == signs


def test_flat_flrw_is_exactly_a_time_dependent_river():
    a, H, r = sp.symbols("a H r", positive=True)
    # chi=r/a(t), so dchi=dr/a-H*r*dt/a.
    jacobian = sp.Matrix([[1, 0], [-H*r/a, 1/a]])
    transformed = sp.simplify(jacobian.T*sp.diag(-1, a*a)*jacobian)
    assert transformed == sp.Matrix([[H*H*r*r-1, -H*r], [-H*r, 1]])
    assert sp.simplify(r*(1-transformed.inv()[1, 1])/2) == H*H*r**3/2


def test_acoustic_sink_satisfies_continuity_and_forced_euler():
    r, mu, cs, rho0, r0 = sp.symbols("r mu cs rho0 r0", positive=True)
    u = -sp.sqrt(2*mu/r)
    rho = rho0*(r/r0)**(-sp.Rational(3, 2))
    pressure = cs*cs*rho
    forcing = -mu/r+sp.Rational(3, 2)*cs*cs*sp.log(r/r0)
    assert sp.simplify(sp.diff(r*r*rho*u, r)) == 0
    assert sp.simplify(u*sp.diff(u, r)+sp.diff(pressure, r)/rho
                       +sp.diff(forcing, r)) == 0
    # Its conformal factor is not constant: cone matching does not remove it.
    assert sp.simplify(sp.diff(sp.log(rho/cs), r)) == -sp.Rational(3, 2)/r


@pytest.mark.parametrize("R", [.5, 1.5, 2.])
def test_retained_israel_wall_is_exactly_a_boundary_stress_jump(R):
    wall = tension_wall(R, .01, 1., .4, 0., 2.)
    curvatures = []
    for side in ("minus", "plus"):
        K = np.diag([wall[f"worldtube_Ktau_{side}"],
                     wall[f"worldtube_Ktheta_{side}"],
                     wall[f"worldtube_Ktheta_{side}"]])
        # Mixed-index BY stress with one common oriented normal, G=1.
        curvatures.append((np.trace(K)*np.eye(3)-K)/(8*math.pi))
    assert curvatures[1]-curvatures[0] == pytest.approx(
        -wall["sigma"]*np.eye(3), abs=3e-15)


def test_static_boundary_pressure_by_direct_covariant_contraction():
    # Compute Christoffel-derived K_ab, not a second call to a BY calculator.
    R, f, fp, G = sp.symbols("R f fp G", positive=True)
    induced = sp.diag(-f, R*R, R*R)  # equatorial point
    normal_r_cov = 1/sp.sqrt(f)
    gamma_r_tt, gamma_r_ang = f*fp/2, -f*R
    Kcov = -normal_r_cov*sp.diag(gamma_r_tt, gamma_r_ang, gamma_r_ang)
    Kmixed = sp.simplify(induced.inv()*Kcov)
    stress = sp.simplify((sp.trace(Kmixed)*sp.eye(3)-Kmixed)/(8*sp.pi*G))
    assert sp.simplify(-stress[0, 0]+sp.sqrt(f)/(4*sp.pi*G*R)) == 0
    expected_p = (fp/(2*sp.sqrt(f))+sp.sqrt(f)/R)/(8*sp.pi*G)
    assert sp.simplify(stress[1, 1]-expected_p) == 0
    M = sp.symbols("M", positive=True)
    # alpha*p is regular at a Schwarzschild horizon even though local p diverges.
    redshifted = sp.simplify((sp.sqrt(f)*stress[1, 1]).subs(
        {f: 1-2*M/R, fp: 2*M/R**2}))
    assert sp.simplify(sp.limit(redshifted, R, 2*M, dir="+")
                       -1/(32*sp.pi*G*M)) == 0


def test_rindler_shear_coefficient_from_the_published_perturbed_metric():
    # Linear-in-velocity part of Bredberg et al. (1101.2451), eq. 14.
    # A divergence-free planar shear isolates the first derivative response.
    tau, r, x, y = sp.symbols("tau r x y", real=True)
    rc, G = sp.symbols("rc G", positive=True)
    epsilon, a, b = sp.symbols("epsilon a b", real=True)
    coordinates = [tau, r, x, y]
    velocities = [a*y, b*x]
    g0 = sp.Matrix([[-r, 1, 0, 0], [1, 0, 0, 0],
                    [0, 0, 1, 0], [0, 0, 0, 1]])
    dg = sp.zeros(4)
    for i, vi in enumerate(velocities, 2):
        dg[0, i] = dg[i, 0] = -(1-r/rc)*vi
        dg[1, i] = dg[i, 1] = -vi/rc
    metric = g0+epsilon*dg
    inverse = g0.inv()-epsilon*g0.inv()*dg*g0.inv()
    connection_r_xy = sum(inverse[1, j]*(sp.diff(metric[j, 3], x)
                           +sp.diff(metric[j, 2], y)
                           -sp.diff(metric[2, 3], coordinates[j]))/2
                           for j in range(4))
    linear_connection = sp.expand(connection_r_xy).coeff(epsilon, 1)
    # n_r=1/sqrt(rc), h_xy=0 at the cutoff; no trace term contributes.
    stress_xy = sp.simplify(linear_connection.subs(r, rc)
                            /(8*sp.pi*G*sp.sqrt(rc)))
    # Proper boundary time is sqrt(rc)*tau and physical V_i=v_i/sqrt(rc).
    strain_xy = (a+b)/(2*sp.sqrt(rc))
    eta = sp.simplify(-stress_xy/(2*strain_xy))
    assert eta == 1/(16*sp.pi*G)
