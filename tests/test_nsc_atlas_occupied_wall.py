"""State counting, signed junctions, fixed-action response and causal clocks."""
import math

import mpmath as mp
import numpy as np
import pytest
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from recursive_horizons.nsc_atlas_occupied_wall import (
    OccupiedWall, filled_shell, nested_branch_curvature,
)
from recursive_horizons.nsc_atlas_parent_observables import DustParent


@pytest.mark.parametrize("q,g", [(0, 2), (1, 2), (7, 1), (1000, 2)])
def test_exact_spectrum_count_not_a_continuum_density(q, g):
    s = filled_shell(q, g)
    assert s["number"] == sum(2*g*n for n in range(1, q+1))
    assert s["energy_coefficient"] == sum(2*g*n*n for n in range(1, q+1))
    nxt = filled_shell(q+1, g)
    assert nxt["energy_coefficient"]-s["energy_coefficient"] == (
        q+1)*(nxt["number"]-s["number"])


def test_conformal_evolution_preserves_the_occupied_projector():
    # For intrinsic ds^2=R(eta)^2(-deta^2+dOmega^2), the rescaled
    # Hamiltonian is the same sphere operator. No adiabatic step is used.
    eig = np.array([-3., -1., 1., 2., 3.])
    occupied = np.diag([0., 0., 1., 1., 0.])
    U = np.diag(np.exp(-1j*eig*.73))
    assert U @ occupied @ U.conj().T == pytest.approx(occupied, abs=2e-16)
    for R in (.8, 1., 1.4):
        assert np.trace(occupied @ np.diag(eig/R)) == pytest.approx(3/R)


@pytest.mark.parametrize("R", [.8, 1., 1.4])
def test_pressure_is_the_area_derivative_and_conserves_surface_energy(R):
    w = OccupiedWall()
    d, t = w.gas_coefficient(1000), w.bare_k
    s = w.matter(R, 1000)
    # Geometric shell mass = 4piG*sigma0 R^2 + ell_P^2 A/R.
    mp_ = 2*t*R-d/R**2
    assert mp_+2*R*s["four_pi_G_pressure"] == pytest.approx(0., abs=2e-15)
    assert (-s["four_pi_G_sigma"]+2*s["four_pi_G_pressure"]) == pytest.approx(-3*t)
    assert s["four_pi_G_sigma"]+s["four_pi_G_pressure"] > 0


def test_reference_has_exact_signed_junctions_and_positive_radial_mode():
    w, s = OccupiedWall(), OccupiedWall().static_state(1000)
    assert s["radius"] == pytest.approx(1., abs=1e-13)
    assert s["mass"] == pytest.approx(9/50, abs=2e-15)
    assert s["beta_minus"] == pytest.approx(24/25)
    assert s["beta_plus"] == pytest.approx(4/5)
    assert s["four_pi_G_sigma"] == pytest.approx(4/25)
    assert s["four_pi_G_pressure"] == pytest.approx(11/150)
    assert s["V_second"] == pytest.approx(35159/7200, abs=2e-12)
    assert s["mass_radial_second_derivative"] == pytest.approx(35159/86400, abs=2e-13)
    for k in ("angular_junction_residual", "temporal_junction_residual",
              "potential_residual", "mass_stationarity_residual"):
        assert abs(s[k]) < 3e-14
    assert w.mass_jet(1.01, 1000)["mass"] > s["mass"]
    assert w.mass_jet(.99, 1000)["mass"] > s["mass"]


def test_mass_hessian_and_potential_are_independently_differentiated():
    import sympy as sp
    R = sp.symbols("R", positive=True)
    h, t, d, M = map(sp.Rational, ("49/625", "1/225", "7/45", "9/50"))
    ms = t*R**2+d/R
    mass = h*R**3/2+sp.sqrt(1-h*R*R)*ms-ms**2/(2*R)
    sj = ms/R
    V = 1-h*R*R-((1-h*R*R)-(1-2*M/R)+sj**2)**2/(4*sj**2)
    assert sp.diff(mass, R).subs(R, 1) == 0
    assert sp.diff(mass, R, 2).subs(R, 1) == sp.Rational(35159, 86400)
    assert sp.simplify(sp.diff(V, R, 2).subs(R, 1)) == sp.Rational(35159, 7200)


def test_neighboring_states_keep_every_coupling_fixed():
    w = OccupiedWall()
    states = [w.static_state(q) for q in (999, 1000, 1001)]
    assert states[0]["radius"] < states[1]["radius"] < states[2]["radius"]
    assert states[0]["mass"] < states[1]["mass"] < states[2]["mass"]
    for s in states:
        assert s["V_second"] == pytest.approx(s["V_second_mass_check"], abs=3e-13)
        assert s["V_second"] > 4
        assert s["beta_plus"] > 0
        assert s["gas_coefficient"]/s["energy_coefficient"] == w.planck_length_squared
        assert abs(s["temporal_junction_residual"]) < 3e-14


def test_high_precision_neighbor_and_mass_first_law():
    w = OccupiedWall()
    with mp.workdps(50):
        h, t = mp.mpf(49)/625, mp.mpf(1)/225
        gh = (mp.mpf(7)/45)/filled_shell(1000)["energy_coefficient"]
        A = filled_shell(1001)["energy_coefficient"]
        def mass(R, spectral):
            ms = t*R*R+gh*spectral/R
            return h*R**3/2+mp.sqrt(1-h*R*R)*ms-ms*ms/(2*R)
        R = mp.findroot(lambda r: mp.diff(lambda x: mass(x, A), r), (.99, 1.02))
        s = w.static_state(1001)
        assert float(R) == pytest.approx(s["radius"], abs=3e-14)
        assert float(mass(R, A)) == pytest.approx(s["mass"], abs=3e-15)
        derivative = mp.diff(lambda coefficient: mass(R, coefficient), A)
        assert float(derivative) == pytest.approx(
            s["mass_derivative_wrt_spectral_coefficient"], rel=2e-14)


def test_wormhole_orientation_cannot_borrow_the_radial_stability():
    s = nested_branch_curvature(1., .18, 49/625, 7/45)
    assert s["F_RR"] < 0 < s["F_z"]
    # Every numerator term is nonpositive, and M>0 makes it strictly negative.
    assert s["necessary_V_second"] < 0
    with pytest.raises(ValueError):
        nested_branch_curvature(.3, .18, 49/625, 7/45)


def test_central_light_across_all_three_regions_against_piecewise_ode():
    w = OccupiedWall()
    s = w.static_state(1000)
    signal = w.central_parent_signal(s)
    M, R, a, b = (s[k] for k in ("mass", "radius", "beta_minus", "beta_plus"))
    H = math.sqrt(w.H_squared)
    # Inside: dR/dt_c=1-H^2R^2; outside: dR/dT=f.
    cap = solve_ivp(lambda t,y: [1-H*H*y[0]**2], (0., 2.), [0.],
                    dense_output=True, rtol=2e-12, atol=2e-13, method="DOP853")
    tc = brentq(lambda t: cap.sol(t)[0]-R, .5, 1.5)
    Tw = a/b*tc
    parent = DustParent(M, signal["initial_join_radius"])
    tauw = brentq(lambda tau: parent.static_time(tau)-Tw, 0., 2.)
    def contact(t, y):
        return y[0]-parent.state(t)["radius"]
    contact.terminal = True
    sol = solve_ivp(lambda t,y: [(1-2*M/y[0])/parent.state(t)["f"]],
                    (tauw, 10.), [R], events=contact, method="DOP853",
                    rtol=2e-12, atol=2e-13)
    tb = sol.t_events[0][0]
    assert tb == pytest.approx(signal["join_reception_time"], abs=5e-11)
    def receiver(t,y):
        return y[0]-8
    receiver.terminal = True
    cosmo = solve_ivp(lambda t,y: [1/parent.state(t)["scale_factor"]],
                      (tb, 20.), [signal["initial_join_radius"]], events=receiver,
                      method="DOP853", rtol=2e-12, atol=2e-13)
    assert cosmo.t_events[0][0] == pytest.approx(signal["reception_time"], abs=8e-11)


def test_redshift_equals_neighboring_proper_clock_arrival_derivative():
    w = OccupiedWall()
    s = w.static_state(1000)
    h = 1e-4
    a, b = [w.central_parent_signal(s, emission_time=t) for t in (-h, h)]
    ratio = (b["reception_time"]-a["reception_time"])/(2*h)
    assert ratio == pytest.approx(w.central_parent_signal(s)["one_plus_redshift"], abs=2e-9)


def test_two_unequal_rooms_fit_the_same_cosmology_without_overlap():
    w = OccupiedWall()
    left, right = [w.central_parent_signal(w.static_state(q)) for q in (999,1001)]
    assert left["parent_density"] == right["parent_density"]
    assert left["initial_join_radius"]+right["initial_join_radius"] < 12
    # Centres at (-6,0,0),(6,0,0); each beacon travels along its own +y ray.
    # The distance from either ray to the other centre is at least 12>R_b.
    assert max(left["initial_join_radius"],right["initial_join_radius"]) < 12
    assert left["one_plus_redshift"] < right["one_plus_redshift"]


def test_unadmitted_states_and_nonstatic_core_are_not_silently_retuned():
    w = OccupiedWall()
    with pytest.raises(ValueError):
        w.static_state(1003)
    with pytest.raises(ValueError):
        filled_shell(-1)
    with pytest.raises(ValueError):
        w.mass_jet(4., 1000)
    with pytest.raises(ValueError):
        w.central_parent_signal(w.static_state(1000), observer_chi=2.)
