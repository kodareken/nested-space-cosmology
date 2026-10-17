"""Independent checks for the parent-end monopole response."""
import math
from functools import lru_cache

import pytest
import sympy as sp

from recursive_horizons.nsc_atlas_gravitational_reduction import (
    KAPPA, boundary_stress_from_extrinsic, exact_tube_identities,
    homogeneous_action_variations, homogeneous_energy_stress,
    schwarzschild_boundary_response,
    schwarzschild_normal_weyl,
)


def _tube_tau(mu, radius):
    """Outward curvature route, written out rather than calling the API algebra."""
    redshift = math.sqrt(1 - 2*mu/radius)
    radial = 0.0 if mu == 0 else mu/(radius*radius*redshift)
    angular = redshift/radius
    trace = radial + 2*angular
    reference_trace = 2/radius

    def component(value, total, reference_value, reference_total):
        return -((value - total) - (reference_value - reference_total))/KAPPA

    return {
        'tau_time_time': component(radial, trace, 0.0, reference_trace),
        'tau_angular': component(angular, trace, 1/radius, reference_trace),
    }


def test_extrinsic_tensor_matches_stable_boundary_response():
    # Owned finite parent boundary, one point: R = 4 ell, mu = 0.4 ell.
    response = schwarzschild_boundary_response(0.4, 4.0)
    stress = boundary_stress_from_extrinsic(0.4, 4.0)
    direct = _tube_tau(0.4, 4.0)
    assert stress['time_time_plus_surface_density'] == pytest.approx(0.0, abs=1e-15)
    assert stress['angular_minus_surface_pressure'] == pytest.approx(0.0, abs=1e-15)
    assert direct['tau_time_time'] == pytest.approx(-response['surface_energy_density'])
    assert direct['tau_angular'] == pytest.approx(response['surface_pressure'])
    redshift = math.sqrt(1 - 0.8/4)
    unstable = ((1 - 0.4/4)/redshift - 1)/(8*math.pi*4)
    assert response['surface_pressure'] == pytest.approx(unstable, abs=1e-15)
    assert response['quasilocal_energy'] == pytest.approx(0.8/(1 + redshift))
    assert response['adm_mass'] == 0.4
    assert response['quasilocal_energy'] != response['adm_mass']


def test_first_law_mass_recovery_and_adm_correction():
    mu, radius = 0.4, 4.0
    base = schwarzschild_boundary_response(mu, radius)
    step = 1e-6
    radial = schwarzschild_boundary_response(mu, radius + step)['quasilocal_energy']
    radial -= schwarzschild_boundary_response(mu, radius - step)['quasilocal_energy']
    radial /= 2*step
    assert radial == pytest.approx(
        -base['surface_pressure']*8*math.pi*radius, abs=1e-10)
    mass_step = 1e-7
    mass_derivative = schwarzschild_boundary_response(mu + mass_step, radius)
    mass_derivative = mass_derivative['quasilocal_energy']
    mass_derivative -= schwarzschild_boundary_response(mu - mass_step, radius)['quasilocal_energy']
    assert mass_derivative/(2*mass_step) == pytest.approx(base['denergy_dmass'], rel=1e-8)
    energy = base['quasilocal_energy']
    assert energy - energy**2/(2*radius) == pytest.approx(mu, abs=1e-12)
    small_mass, small_radius = 0.02, 1.0
    small = schwarzschild_boundary_response(small_mass, small_radius)
    series = (small_mass + small_mass**2/(2*small_radius)
              + small_mass**3/(2*small_radius**2))
    quartic = 5*small_mass**4/(8*small_radius**3)
    assert (small['quasilocal_energy'] - series)/quartic == pytest.approx(1.0, rel=0.05)


def test_scaling_and_zero_mass_limit():
    base = schwarzschild_boundary_response(0.4, 4.0)
    scaled = schwarzschild_boundary_response(1.2, 12.0)
    assert scaled['lapse_factor'] == pytest.approx(base['lapse_factor'])
    assert scaled['quasilocal_energy'] == pytest.approx(3*base['quasilocal_energy'])
    assert scaled['surface_energy_density'] == pytest.approx(base['surface_energy_density']/3)
    assert scaled['surface_pressure'] == pytest.approx(base['surface_pressure']/3)
    flat = schwarzschild_boundary_response(0.0, 4.0)
    assert flat['lapse_factor'] == pytest.approx(1.0)
    assert flat['quasilocal_energy'] == pytest.approx(0.0, abs=0.0)
    assert flat['surface_energy_density'] == pytest.approx(0.0, abs=0.0)
    assert flat['surface_pressure'] == pytest.approx(0.0, abs=0.0)
    assert flat['denergy_dmass'] == pytest.approx(1.0)
    with pytest.raises(ValueError):
        schwarzschild_boundary_response(0.4, 0.8)


def test_symbolic_tube_and_action_variations_agree():
    identities = exact_tube_identities()
    for value in identities.values():
        assert value == 0
    derived = homogeneous_action_variations()
    assert derived['pressure_plus_dE'] == 0
    assert derived['fixed_pressure'] == 0
    lapse, scale, coordinate_volume = sp.symbols('N a Vc', positive=True)
    energy = sp.Function('E')
    volume = scale**3*coordinate_volume
    integrand = -lapse*energy(volume)
    independent_pressure = sp.simplify(
        sp.diff(integrand, scale)/(lapse*sp.diff(volume, scale)))
    independent_density = sp.simplify(-sp.diff(integrand, lapse)/volume)
    assert sp.simplify(independent_pressure - derived['pressure']) == 0
    assert sp.simplify(independent_density - derived['density']) == 0
    assert derived['extensive_vacuum_defect'] == 0
    assert derived['extensive_pressure'] == -derived['extensive_density']


def test_fixed_energy_is_dust_and_extensive_energy_is_vacuum_form():
    dust = homogeneous_energy_stress(2.5, 0.0, 4.0)
    assert dust['rho_geom'] == pytest.approx(2.5/4)
    assert dust['p_geom'] == pytest.approx(0.0, abs=0.0)
    assert dust['kappa_rho'] == pytest.approx(KAPPA*dust['rho_geom'])
    assert dust['lambda_trace'] == pytest.approx(2*math.pi*dust['rho_geom'])
    assert dust['vacuum_form_defect'] == pytest.approx(dust['kappa_rho'])
    density = 0.3
    volume = 5.0
    vacuum = homogeneous_energy_stress(density*volume, density, volume)
    assert vacuum['p_geom'] == pytest.approx(-vacuum['rho_geom'])
    assert vacuum['vacuum_form_defect'] == pytest.approx(0.0, abs=1e-15)
    assert vacuum['lambda_trace'] == pytest.approx(KAPPA*density)
    # A constant Brown-York energy is not turned into vacuum form.
    tube = schwarzschild_boundary_response(0.4, 4.0)
    misread = homogeneous_energy_stress(tube['quasilocal_energy'], 0.0, 10.0)
    assert misread['vacuum_form_defect'] > 0.1


def test_normal_weyl_from_coordinate_riemann_and_radial_boost():
    t, r, theta, phi = sp.symbols('t r theta phi', real=True)
    mass, velocity = sp.symbols('M v', real=True)
    coordinates = (t, r, theta, phi)
    f = 1-2*mass/r
    metric = sp.diag(-f, 1/f, r**2, r**2*sp.sin(theta)**2)
    inverse = metric.inv()

    @lru_cache(None)
    def gamma(a, b, c):
        return sp.simplify(sum(inverse[a, d]*(
            sp.diff(metric[d, c], coordinates[b])+sp.diff(metric[d, b], coordinates[c])
            -sp.diff(metric[b, c], coordinates[d])) for d in range(4))/2)

    @lru_cache(None)
    def riemann(a, b, c, d):
        return sp.simplify(sp.diff(gamma(a, d, b), coordinates[c])
                           -sp.diff(gamma(a, c, b), coordinates[d])
                           +sum(gamma(a, c, e)*gamma(e, d, b)
                                -gamma(a, d, e)*gamma(e, c, b) for e in range(4)))

    ricci = sp.Matrix(4, 4, lambda a, b: sp.simplify(sum(riemann(c, a, c, b) for c in range(4))))
    assert ricci == sp.zeros(4)  # Therefore these Riemann components are Weyl.
    radial = sp.simplify(metric[0, 0]*riemann(0, 1, 0, 1))
    angular_time = sp.simplify(metric[0, 0]*riemann(0, 2, 0, 2)/(f*r**2))
    angular_radial = sp.simplify(metric[1, 1]*riemann(1, 2, 1, 2)*f/r**2)
    assert radial == -2*mass/r**3
    assert angular_time == mass/r**3
    boosted = sp.simplify((angular_radial+velocity**2*angular_time)/(1-velocity**2))
    assert boosted == -mass/r**3
    assert sp.simplify(-radial+2*boosted) == 0
    row = schwarzschild_normal_weyl(.4, 1.5)
    substitution = {mass: sp.Rational(2, 5), r: sp.Rational(3, 2)}
    assert row['C_u_n_u_n'] == pytest.approx(float(radial.subs(substitution)))
    assert row['C_theta_n_theta_n'] == pytest.approx(float(boosted.subs(substitution)))
    assert row['C_phi_n_phi_n'] == row['C_theta_n_theta_n']
    scaled = schwarzschild_normal_weyl(1.2, 4.5)
    for key, value in row.items():
        assert scaled[key] == pytest.approx(value/9)
    assert all(value == 0. for value in schwarzschild_normal_weyl(0., 1.5).values())


def test_pure_trace_shell_has_no_singular_weyl_tensor():
    # R_{n a n b}^{sing}=k*h_ab from [K_ab]=-k*h_ab. Normal index3.
    k = sp.symbols('k', real=True)
    metric = sp.diag(-1, 1, 1, 1)
    tangential = sp.diag(-k, k, k, 0)
    normal = (0, 0, 0, 1)

    def curvature(a, b, c, d):
        return (normal[a]*normal[c]*tangential[b, d]
                -normal[a]*normal[d]*tangential[b, c]
                -normal[b]*normal[c]*tangential[a, d]
                +normal[b]*normal[d]*tangential[a, c])

    ricci = sp.Matrix(4, 4, lambda a, b: sum(metric[i, i]*curvature(i, a, i, b) for i in range(4)))
    scalar = sum(metric[i, i]*ricci[i, i] for i in range(4))
    assert ricci == sp.diag(-k, k, k, 3*k)
    assert scalar == 6*k
    for a in range(4):
        for b in range(4):
            for c in range(4):
                for d in range(4):
                    weyl = (curvature(a, b, c, d)
                            -(metric[a, c]*ricci[d, b]-metric[a, d]*ricci[c, b]
                              -metric[b, c]*ricci[d, a]+metric[b, d]*ricci[c, a])/2
                            +scalar*(metric[a, c]*metric[d, b]-metric[a, d]*metric[c, b])/6)
                    assert sp.simplify(weyl) == 0
