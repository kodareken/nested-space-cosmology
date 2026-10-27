"""Independent cylindrical, nonlinear and energy controls for the vortex family."""
import math

import numpy as np
import pytest
from scipy.integrate import quad
import sympy as sp

from recursive_horizons.nsc_atlas_vortex_transfer import (
    KAPPA, angular_advection_average, bump, direct_integrals, dissipation,
    energy, exact_theta_rate, jets, mean_theta_rate, moments, profile,
    radial_fluxes, reynolds, stress_divergence,
)


def test_vector_potential_curl_is_symbolically_divergence_free():
    r, z, m, k = sp.symbols("r z m k", positive=True, real=True)
    a, b, c = sp.Function("a")(r), sp.Function("b")(r), sp.Function("c")(z)
    # Curl of (0,b*c*exp(i phase)/m,a*c*exp(i phase)/m), with phase removed.
    Wr = -b*sp.diff(c, z)/m+sp.I*a*c/r
    Wt = -sp.diff(a, r)*c/m-sp.I*k*a*c/r
    Wz = (sp.diff(b, r)+b/r)*c/m+sp.I*k*b*c/r
    divergence = sp.diff(Wr, r)+sp.I*k*m*Wr/r+Wr/r+sp.I*m*Wt/r+sp.diff(Wz, z)
    assert sp.simplify(divergence) == 0


def test_midplane_acceleration_identity_before_any_numerical_integration():
    r, m, k = sp.symbols("r m k", positive=True, real=True)
    a = sp.Function("a")(r)
    # At the symmetry plane chi'=0, R_rt=-k*a^2/(2r^2), d_z R_zt=0.
    Rrt = -k*a*a/(2*r*r)
    rate = -sp.diff(Rrt, r)-2*Rrt/r
    assert sp.simplify(rate-k*a*sp.diff(a, r)/r**2) == 0
    # Actual profile derivatives give both exact observable rates.
    assert exact_theta_rate(.75) == pytest.approx(-1024/81*math.exp(-2/3), abs=3e-14)
    assert exact_theta_rate(1.25) == pytest.approx(1024/225*math.exp(-2/3), abs=3e-14)


@pytest.mark.parametrize("mode", [4, 8, 16])
def test_direct_nonlinear_advection_gives_the_predicted_initial_ns_rate(mode):
    for r, z in [(.75, 0.), (1.25, 0.), (.8, .3), (1.15, -.4)]:
        direct = angular_advection_average(r, z, mode)
        stress = -stress_divergence(r, z, mode)[1]
        assert direct == pytest.approx(stress, abs=3e-13)
    assert mean_theta_rate(.75, mode) == pytest.approx(exact_theta_rate(.75), abs=3e-14)


def test_zero_velocity_mean_does_not_remove_the_positive_momentum_flux():
    r, z, m = .8, .2, 8
    W = jets(r, z, m)[0]
    theta = np.arange(65)*2*math.pi/65
    U = np.real(W[:, None]*np.exp(1j*(m*theta+KAPPA*m*math.log(r))))
    assert np.mean(U, axis=1) == pytest.approx(np.zeros(3), abs=3e-15)
    covariance = U @ U.T/U.shape[1]
    R = reynolds(r, z, m)
    assert covariance == pytest.approx(R, abs=3e-14)
    assert np.linalg.eigvalsh(R).min() > -2e-14
    # Any positive covariance has a minimum energy cost for these two fluxes.
    assert np.trace(R)/2 >= math.hypot(R[0, 1], R[0, 2])
    assert R[0, 1] != 0


def test_stress_splits_exactly_into_a_retained_term_and_inverse_mode_squared():
    # This independently checks R_m = R_inf + R_2/m^2 via three finite modes.
    r, z = .8, .3
    A, B, C = [reynolds(r, z, m) for m in (4, 8, 16)]
    assert A-B == pytest.approx(4*(B-C), abs=3e-15)
    limiting = (4*B-A)/3
    assert np.linalg.norm(limiting) > 1
    assert np.linalg.eigvalsh(limiting).min() > -2e-14


@pytest.mark.parametrize("mode", [4, 8, 16])
def test_energy_and_viscous_loss_against_direct_curl_and_gradient_integrals(mode):
    direct = direct_integrals(mode, order=128)
    assert energy(mode, order=128) == pytest.approx(direct["energy"], abs=3e-12)
    assert dissipation(mode, order=128) == pytest.approx(direct["curl_loss"], abs=2e-11)
    assert direct["curl_loss"] == pytest.approx(direct["gradient_loss"], abs=2e-11)
    assert energy(mode) > moments()["E0"] > 0


def test_momentum_fluxes_and_global_initial_torque_conserve_momentum():
    I = quad(lambda z: float(bump(z)[0])**2, -1., 1., epsabs=1e-13)[0]
    for r in (.75, 1., 1.25):
        for m in (4, 16):
            angular = 2*math.pi*r*r*quad(lambda z: float(reynolds(r, z, m)[0, 1]), -1., 1., epsabs=1e-12)[0]
            axial = 2*math.pi*r*quad(lambda z: float(reynolds(r, z, m)[0, 2]), -1., 1., epsabs=1e-12)[0]
            a = float(profile(r)[0])
            assert angular == pytest.approx(2*math.pi*a*a*I, abs=2e-12)
            assert axial == pytest.approx(math.pi*a*a*I/r, abs=2e-12)
    # Height-integrated initial mean angular acceleration, weighted by volume
    # and moment arm. The compact flux gives opposite torques, with zero sum.
    torque = quad(lambda r: -4*math.pi*I*float(profile(r)[0]*profile(r)[1]),
                  .5, 1.5, epsabs=1e-12)[0]
    assert abs(torque) < 2e-12


def test_fine_resolution_retains_force_but_increases_the_viscous_cost():
    m = (4, 8, 16)
    assert energy(m[0]) > energy(m[1]) > energy(m[2]) > moments()["E0"]
    assert dissipation(m[0]) < dissipation(m[1]) < dissipation(m[2])
    assert [float(mean_theta_rate(.75, j)) for j in m] == pytest.approx([exact_theta_rate(.75)]*3)


def test_compact_support_and_invalid_domains():
    for r, z in [(.25, 0.), (2., 0.), (.8, 1.5), (.5, 0.), (1.5, 0.)]:
        assert jets(r, z, 4)[0] == pytest.approx(np.zeros(3), abs=0.)
    with pytest.raises(ValueError):
        jets(0., 0., 4)
    with pytest.raises(ValueError):
        jets(1., 0., 1.5)
    with pytest.raises(ValueError):
        dissipation(4, 0.)


def test_source_core_scalings_are_energy_integrable_but_critical_norm_diverges():
    h = sp.symbols("h", positive=True)
    volume = sp.Rational(3, 2)-h
    velocity = -sp.Rational(1, 2)-h
    E = sp.simplify(2*velocity+volume)
    D = sp.simplify(2*(velocity-sp.Rational(1, 2))+volume)
    L3 = sp.simplify(velocity+volume/3)
    assert E == sp.Rational(1, 2)-3*h
    assert D == -sp.Rational(1, 2)-3*h
    assert L3 == -4*h/3
    assert E.subs(h, sp.Rational(1, 100)) > 0
    assert D.subs(h, sp.Rational(1, 100)) > -1
    assert L3 < 0
