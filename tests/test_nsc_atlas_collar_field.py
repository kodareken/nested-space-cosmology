"""Independent background, original-action and continuum weak-equation controls."""
import math

import numpy as np
import pytest

from recursive_horizons.nsc_atlas_collar_field import (
    ETA, MASS, QuietCollar, assemble, continuum_force, hermite_basis,
    manufactured, manufactured_coefficients, positivity_bounds,
)
from recursive_horizons.nsc_atlas_conformal_action import time_kinetic_matrix


def test_normal_geodesics_jacobi_clock_and_parameter_derivatives():
    bounds = positivity_bounds()
    assert bounds['radius_lower'] > .42
    assert bounds['lapse_lower'] > .88
    n = np.array([-ETA, -.01, 0., .01, ETA])
    for R in (.45, .5, .55):
        g = QuietCollar(R)
        r, rn, N, Nn, w, wn = g.at(n)
        f = np.where(n <= 0, 1-r*r, 1-2*MASS/r)
        rnn = np.where(n <= 0, -r, MASS/r**2)
        assert rn*rn-g.v*g.v == pytest.approx(f, abs=2e-12)
        assert rn*Nn-rnn*N == pytest.approx(g.wall['radial_acceleration'], abs=2e-11)
        # Direct two-metric conformal curvature includes the wall-time term.
        box_r = -g.wall['radial_acceleration']/N+rn*Nn/N+rnn
        Nnn = np.where(n <= 0, -N, -2*MASS*N/r**3)
        conformal_R4 = -2*r*r*Nnn/N+2*r*box_r-2*(rn*rn-g.v*g.v)+2
        assert conformal_R4 == pytest.approx(np.where(n <= 0, 0., 12*MASS/r), abs=2e-11)
        assert np.min(r) >= bounds['radius_lower']
        assert np.min(N) >= bounds['lapse_lower']
        wall = g.wall_sides()
        assert wall[1, 1]-wall[1, 0] == pytest.approx(-2*R)
        assert wall[3, 1]-wall[3, 0] == pytest.approx(-2)
    h = 1e-6
    g = QuietCollar(.5).at(n)
    fd = (QuietCollar(.5+h).at(n)-QuietCollar(.5-h).at(n))/(2*h)
    assert fd[0] == pytest.approx(g[2], abs=2e-9)
    assert fd[2] == pytest.approx(g[4], abs=2e-9)
    assert fd[3] == pytest.approx(g[5], abs=2e-8)


def test_shared_cubic_representation_retains_actual_transmission_and_radius_kink():
    n = np.array([-ETA, -.013, 0., .007, ETA])
    B, D = hermite_basis(n)
    q, u, _ = manufactured_coefficients(.49)
    exact = manufactured(.49, n)
    assert B @ q == pytest.approx(exact['chi'], abs=2e-17)
    assert D @ q == pytest.approx(exact['normal'], abs=2e-16)
    assert B @ u == pytest.approx(exact['time'], abs=2e-16)
    # Separate one-sided cubic endpoint formulas have identical physical traces.
    assert q[2] == pytest.approx(exact['chi'][2])
    assert q[3]/ETA == pytest.approx(exact['normal'][2])
    g = QuietCollar(.49)
    r, rn, N, Nn, _, _ = g.wall_sides()
    normal_phi = q[3]/(ETA*r)-q[2]*rn/r**2
    assert np.diff(normal_phi)[0] == pytest.approx(2*q[2]/r[0], abs=2e-16)
    assert rn[1] != rn[0]
    assert np.diff((Nn-N*rn/r)/3)[0] == pytest.approx(0., abs=2e-15)


def test_finite_legendre_map_is_original_ADM_momentum_not_wave_velocity():
    a = assemble(.49)
    q, u, _ = manufactured_coefficients(.49)
    chi, chi_t = a.B @ q, a.B @ u
    density = []
    for i, r in enumerate(a.profile[0]):
        N = a.profile[2, i]
        matrix = time_kinetic_matrix(radius=r, field=chi[i], kappa=8*math.pi,
                                    radial_length=1/r, lapse=N/r)
        velocity = np.array([-a.geometry.v*N/r**2, a.geometry.v*N, chi_t[i]])
        density.append((matrix @ velocity)[2])
    original = a.B.T @ (a.weights*np.array(density))
    p = a.M @ u+a.S @ q
    assert p == pytest.approx(original, abs=2e-17)
    # Scaled Hermite slopes still have a conditioned six-dimensional solve.
    roundoff = 8*np.finfo(float).eps*np.linalg.cond(a.M)*np.max(np.abs(u))
    assert a.velocity(q, p) == pytest.approx(u, abs=roundoff)
    assert np.max(np.abs(a.S @ q)) > 1e-6
    assert np.min(np.linalg.eigvalsh(a.M)) > 0


def test_continuum_forcing_GHY_and_metric_work_match_the_finite_action():
    R, h = .49, 1e-6
    a = assemble(R)
    q, u, utt = manufactured_coefficients(R)
    force = continuum_force(a)
    assert a.M @ utt+a.geometry.v*a.MR @ u+a.K @ q == pytest.approx(force, abs=8e-15)
    # Differentiate independently assembled geometry and actual known field.
    before, after = assemble(R-h), assemble(R+h)
    for name in ('M', 'K', 'S'):
        fd = (getattr(after, name)-getattr(before, name))/(2*h)
        assert fd == pytest.approx(getattr(a, name+'R'), abs=2e-7)
    qm, um, _ = manufactured_coefficients(R-h)
    qp, up, _ = manufactured_coefficients(R+h)
    finite_energy_R = (after.energy(qp, up)-before.energy(qm, um))/(2*h)
    work = (-u @ a.MR @ u+q @ a.KR @ q)/2+u @ force/a.geometry.v
    assert finite_energy_R == pytest.approx(work, abs=2e-9)
    # The time-end scalar generating function shifts the momentum, not the force.
    qR, pR = u/a.geometry.v, a.SR @ q+a.S @ (u/a.geometry.v)+(force-a.K @ q)/a.geometry.v
    pm = before.M @ um+before.S @ qm
    pp = after.M @ up+after.S @ qp
    assert (pp-pm)/(2*h) == pytest.approx(pR, abs=2e-9)


def test_admitted_domains_are_enforced():
    with pytest.raises(ValueError): QuietCollar(.6)
    with pytest.raises(ValueError): QuietCollar(.5).at([ETA+.001])
    with pytest.raises(ValueError): hermite_basis([float('nan')])
    with pytest.raises(ValueError): assemble(.5, 2)
