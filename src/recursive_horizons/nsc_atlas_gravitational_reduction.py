"""Parent-end response of one quiet cap-wall-Schwarzschild body.

G = c = 1. ``mass_length`` is mu = GM/c^2. Quasilocal energy and homogeneous
energy are geometric lengths. Boundary ``surface_energy_density`` and
``surface_pressure`` are inverse lengths. Bulk ``rho_geom`` and ``p_geom``
are inverse lengths squared, and ``kappa_rho = 8*pi*rho_geom``.

The tube is the owned static right end, outward toward increasing areal
radius, with proper boundary lapse and Minkowski subtraction. Brown and
York, gr-qc/9209012, Section III Eqs. 3.14-3.19 and Section VI Eqs.
6.2-6.20, use the opposite sign for the two-sphere extrinsic curvature.
The expressions here are derived with the positive outward trace
k = 2*sqrt(f)/R and flat reference k0 = 2/R. Their final epsilon, pressure
and energy agree with those equations. The TeX source was not executed.

This is the leading parent-side monopole. It is not the cap-centre scalar
map, not a four-dimensional vacuum density, and not a constructed parent
cosmology. Finite-size, tidal, interaction and velocity corrections are
omitted, with no proved numerical bound.
"""
from __future__ import annotations

import math
import time

KAPPA = 8*math.pi
SCOPE = (
    'leading conserved-monopole response at the parent end; '
    'no population, no prescribed expansion, no vacuum identification '
    'from the tube pressure or the dust trace'
)


def _real(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(name + ' must be a finite real number')
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(name + ' must be a finite real number')
    return number


def schwarzschild_boundary_response(mass_length, radius):
    """Proper-lapse Brown-York response for R > 2*mu >= 0.

    ADM mass is ``mass_length``, not the finite-radius energy. Surface
    pressure is the tube stress, not a bulk pressure. The existing parent
    boundary R = 4*ell is one allowed input, not a radius campaign.
    """
    mu = _real(mass_length, 'mass length')
    radius = _real(radius, 'radius')
    if mu < 0 or radius <= 2*mu:
        raise ValueError('need radius > 2*mass_length >= 0')
    discriminant = 1 - 2*mu/radius
    redshift = math.sqrt(discriminant)
    # E = R(1-s) = 2*mu/(1+s) avoids cancellation at large radius.
    energy = 2*mu/(1 + redshift) if mu > 0 else 0.0
    area = 4*math.pi*radius*radius
    epsilon = energy/area
    # ((1-mu/R)/s-1)/(8*pi*R), rewritten without subtracting nearly equal terms.
    pressure = (mu/radius)**2/(4*math.pi*radius*redshift*(1 + redshift)**2)
    return {
        'mass_length': mu,
        'adm_mass': mu,
        'radius': radius,
        'lapse_factor': redshift,
        'quasilocal_energy': energy,
        'surface_energy_density': epsilon,
        'surface_pressure': pressure,
        'denergy_dmass': 1/redshift,
        'scope': SCOPE,
    }


def boundary_stress_from_extrinsic(mass_length, radius):
    """tau from the outward tube curvature, independent of the stable energy form.

    Orthonormal components K^a_b = (mu/(R^2 s), s/R, s/R). The flat reference
    is (0, 1/R, 1/R). tau^a_b = -[(K^a_b - K delta^a_b) - reference]/(8*pi).
    """
    base = schwarzschild_boundary_response(mass_length, radius)
    mu, radius = base['mass_length'], base['radius']
    redshift = base['lapse_factor']
    radial = mu/(radius*radius*redshift) if mu > 0 else 0.0
    angular = redshift/radius
    trace = radial + 2*angular
    reference_trace = 2/radius
    reference_angular = 1/radius

    def component(value, total, reference_value, reference_total):
        return -((value - total) - (reference_value - reference_total))/KAPPA

    time_time = component(radial, trace, 0.0, reference_trace)
    sphere = component(angular, trace, reference_angular, reference_trace)
    return {
        'K_time_time': radial,
        'K_angular': angular,
        'tau_time_time': time_time,
        'tau_angular': sphere,
        'time_time_plus_surface_density': time_time + base['surface_energy_density'],
        'angular_minus_surface_pressure': sphere - base['surface_pressure'],
    }


def schwarzschild_normal_weyl(mass_length, radius):
    """C_{a n b n} in an orthonormal radially moving wall frame, outside 2M.

    Signature -+++ and R^a_bcd = partial_c Gamma^a_db - partial_d Gamma^a_cb
    + Gamma^a_ce Gamma^e_db - Gamma^a_de Gamma^e_cb. The normal is spacelike;
    tangential indices are (u, theta, phi). A radial boost leaves these
    components unchanged. Values have inverse-length-squared units.
    This regular curvature tensor is not a completed higher-derivative junction.
    """
    mu = _real(mass_length, 'mass length')
    radius = _real(radius, 'radius')
    if mu < 0 or radius <= 2*mu:
        raise ValueError('need radius > 2*mass_length >= 0')
    amplitude = mu/radius**3
    return {'C_u_n_u_n': -2*amplitude,
            'C_theta_n_theta_n': -amplitude,
            'C_phi_n_phi_n': -amplitude,
            'induced_trace': 0.0}


def homogeneous_energy_stress(energy_length, denergy_dvolume, physical_volume):
    """Bulk stress from an on-shell energy functional of physical volume.

    ``denergy_dvolume`` is the metric variation of that functional, not a
    prescribed time derivative and not the Brown-York tube energy. A fixed
    energy is dust. ``energy = density*volume`` with constant density is
    vacuum form. No tolerance classifies the result.
    """
    energy = _real(energy_length, 'energy length')
    slope = _real(denergy_dvolume, 'energy derivative')
    volume = _real(physical_volume, 'physical volume')
    if volume <= 0:
        raise ValueError('physical volume must be positive')
    rho = energy/volume
    pressure = -slope
    return {
        'energy_length': energy,
        'denergy_dvolume': slope,
        'physical_volume': volume,
        'rho_geom': rho,
        'p_geom': pressure,
        'kappa_rho': KAPPA*rho,
        'kappa_p': KAPPA*pressure,
        'lambda_trace': 2*math.pi*(rho - 3*pressure),
        'vacuum_form_defect': KAPPA*(rho + pressure),
        'scope': SCOPE,
    }


def homogeneous_action_variations():
    """Vary S = -∫ N E(a^3 Vc) dt under a homogeneous lapse and scale.

    The spatial variation identifies p = -dE/dV. The lapse variation
    identifies rho = E/V. Fixed E is independent of a. Extensive E is a
    constant density times the physical volume.
    """
    import sympy as sp

    lapse, scale, coordinate_volume = sp.symbols('N a Vc', positive=True)
    energy = sp.Function('E')
    volume = scale**3*coordinate_volume
    integrand = -lapse*energy(volume)
    lapse_derivative = sp.diff(integrand, lapse)
    scale_derivative = sp.diff(integrand, scale)
    volume_derivative = sp.diff(volume, scale)
    pressure = sp.simplify(scale_derivative/(lapse*volume_derivative))
    density = sp.simplify(-lapse_derivative/volume)
    volume_symbol = sp.symbols('V', positive=True)
    pressure_on_volume = pressure.subs(volume, volume_symbol)
    fixed = sp.Symbol('E_fixed', real=True)
    density0 = sp.Symbol('rho0', real=True)
    fixed_integrand = -lapse*fixed
    fixed_pressure = sp.simplify(
        sp.diff(fixed_integrand, scale)/(lapse*volume_derivative))
    extensive_integrand = -lapse*density0*volume
    extensive_pressure = sp.simplify(
        sp.diff(extensive_integrand, scale)/(lapse*volume_derivative))
    extensive_density = sp.simplify(-sp.diff(extensive_integrand, lapse)/volume)
    return {
        'density': density,
        'pressure': pressure,
        'pressure_plus_dE': sp.simplify(
            pressure_on_volume + sp.diff(energy(volume_symbol), volume_symbol)),
        'fixed_pressure': fixed_pressure,
        'extensive_pressure': extensive_pressure,
        'extensive_density': extensive_density,
        'extensive_vacuum_defect': sp.simplify(extensive_density + extensive_pressure),
    }


def exact_tube_identities():
    """Symbolic tube relations used by the checks, not a radius campaign."""
    import sympy as sp

    mu, radius = sp.symbols('mu R', positive=True)
    redshift = sp.sqrt(1 - 2*mu/radius)
    energy = 2*mu/(1 + redshift)
    direct_energy = radius*(1 - redshift)
    direct_pressure = ((1 - mu/radius)/redshift - 1)/(8*sp.pi*radius)
    stable_pressure = (mu/radius)**2/(4*sp.pi*radius*redshift*(1 + redshift)**2)
    recovered = sp.simplify(energy - energy**2/(2*radius))
    differential = sp.simplify(
        sp.diff(energy, radius) + stable_pressure*8*sp.pi*radius)
    mass_differential = sp.simplify(sp.diff(energy, mu) - 1/redshift)
    series = energy.series(mu, 0, 3).removeO()
    return {
        'energy_forms': sp.simplify(energy - direct_energy),
        'pressure_forms': sp.simplify(stable_pressure - direct_pressure),
        'recovered_mass': sp.simplify(recovered - mu),
        'fixed_mass_work': differential,
        'mass_work': mass_differential,
        'adm_series': sp.simplify(series - (mu + mu**2/(2*radius))),
    }


def spherical_einstein_response(*, radius, mass_length, hubble=0.):
    """Mixed Einstein eigenvalues from the quiet spherical metric's jets.

    f=1-2*mu/R-H^2*R^2. No wall distribution or quantum backreaction is
    included. The returned trace is R/4; identifying it with Lambda also
    requires the four Einstein eigenvalues to agree. This direct metric
    calculation is independent of the initial-slice scalar fitting map.
    """
    radius = _real(radius, 'radius')
    mu = _real(mass_length, 'mass length')
    H = _real(hubble, 'Hubble inverse length')
    if radius <= 0 or mu < 0 or H < 0:
        raise ValueError('positive radius and nonnegative mass and H required')
    f = 1-2*mu/radius-H*H*radius*radius
    first = 2*mu/radius**2-2*H*H*radius
    second = -4*mu/radius**3-2*H*H
    radial = (radius*first+f-1)/radius**2
    angular = second/2+first/radius
    eigenvalues = (radial, radial, angular, angular)
    projection = -math.fsum(eigenvalues)/4
    return dict(f=f, mixed_einstein_eigenvalues=eigenvalues,
                lambda_trace=projection,
                vacuum_form_defect=max(abs(value+projection) for value in eigenvalues),
                weyl_squared=48*mu*mu/radius**6)


def paired_nesting_vacuum_test(*, ell_m=10000.):
    """One fresh one-level/two-level comparison, with the same exact outer end.

    Reuses the existing fixed quadratures and initial normal-clock rule.
    Changes in scalar fitting sources are NOT a tensor coarse graining or
    a vacuum term. An independent exterior metric comparison tests the
    actual quiet parent-side curvature at the common Rcut=4*ell.
    No physical evolution, grid scan, source fit or new state is used.
    """
    from .nsc_atlas_geometric_average import geometric_slice_average
    from .nsc_atlas_nested_slice import nested_slice_record

    ell = _real(ell_m, 'supplied physical ruler')
    if ell <= 0:
        raise ValueError('supplied physical ruler must be positive')
    started = time.process_time()
    base = geometric_slice_average(quadrature='adaptive')
    base32 = geometric_slice_average(quadrature='gauss32')
    base64 = geometric_slice_average(quadrature='gauss64')
    nested = nested_slice_record()
    mass = base['domain']['M_over_ell']
    cut = base['domain']['Rcut_right_over_ell']
    if (mass != nested['assumptions']['exterior_mass_over_ell']
            or cut != nested['analytic_end']['retained_R_over_ell']
            or base['domain']['k_ell'] != nested['assumptions']['k_ell']
            or base['domain']['H_ell'] != nested['assumptions']['inner_hubble_ell']):
        raise ValueError('the two geometries do not share the declared controls')
    keys = ('Q', 'mean_intrinsic_curvature', 'kappa_geometric_density',
            'kappa_geometric_pressure', 'geometric_trace_projection',
            'geometric_vacuum_condition')
    left = base['normalized_sources']
    right = nested['average_diagnostic']['scalar_outputs_in_ell_units']
    differences = {key: right[key]-left[key] for key in keys}
    # Keep the two geometries' discrepancies separate: cancellation between
    # them must not make the numerical indicator look artificially small.
    # Their absolute sum is still only an indicator, not an error bound.
    gap = nested['quadrature']['scalar_diagnostic_64_minus_32']
    comparisons = {key: dict(
        nested_64_minus_32=abs(gap[key]),
        baseline_adaptive_minus_32=abs(left[key]-base32['normalized_sources'][key]),
        baseline_adaptive_minus_64=abs(left[key]-base64['normalized_sources'][key]))
        for key in keys}
    indicators = {key: values['nested_64_minus_32']
                  +max(values['baseline_adaptive_minus_32'],
                       values['baseline_adaptive_minus_64'])
                  for key, values in comparisons.items()}
    density = differences['kappa_geometric_density']
    pressure = differences['kappa_geometric_pressure']
    # A trace and the independent rho+p condition must both be retained.
    trace = (density-3*pressure)/4
    nonvacuum = density+pressure
    scale = max(abs(density), abs(pressure))
    if max(abs(trace-differences['geometric_trace_projection']),
           abs(nonvacuum+differences['geometric_vacuum_condition']/3)) > 2e-12:
        raise RuntimeError('independent scalar stress contractions disagree')
    outer_base = spherical_einstein_response(radius=cut, mass_length=mass)
    outer_nested = spherical_einstein_response(
        radius=nested['analytic_end']['retained_R_over_ell'],
        mass_length=nested['assumptions']['exterior_mass_over_ell'],
        hubble=nested['assumptions']['exterior_hubble_ell'])
    exterior_delta = outer_nested['lambda_trace']-outer_base['lambda_trace']
    elapsed = time.process_time()-started
    if elapsed > 30:
        raise RuntimeError('paired nesting test exceeded 30 process-CPU seconds')
    return dict(
        controls=dict(M_over_ell=mass, Rcut_over_ell=cut, k_ell=base['domain']['k_ell'],
                      cap_H_ell=base['domain']['H_ell'],
                      clock='each declared instantaneous Gaussian normal; common static outer tube',
                      reference_length_m=ell, reference_length_is_supplied=True),
        baseline_scalar_sources={key: left[key] for key in keys},
        nested_scalar_sources={key: right[key] for key in keys},
        nested_minus_baseline=differences,
        separate_quadrature_comparisons=comparisons,
        summed_quadrature_indicators=indicators,
        geometric_source_scope=('Only the Q and mean-3-curvature fitting correction; '
            'supplied cap vacuum and wall/matter mean energy and pressure are '
            'not included. One slice cannot separate a possible constant '
            'component from other stress.'),
        stress_test=dict(kappa_delta_rho_plus_delta_p_times_ell2=nonvacuum,
                         delta_pressure_over_delta_density=pressure/density if density else None,
                         relative_nonvacuum_residual=abs(nonvacuum)/scale if scale else 0.,
                         scalar_trace_per_m2=trace/ell**2),
        exact_quiet_exterior=dict(
            baseline=outer_base, nested=outer_nested,
            delta_Lambda_times_ell2=exterior_delta,
            delta_Lambda_per_m2=exterior_delta/ell**2,
            interpretation='identical spherical vacuum metrics with fixed outer mass and Lambda; no quantum backreaction or collective parent'),
        counts=dict(new_evolution_grids=0, field_solves=0, evolution_steps=0,
                    sampled_spatial_values=(base['spatial_evaluations']
                        +base32['spatial_evaluations']+base64['spatial_evaluations']
                        +nested['geometric_evaluations'])),
        process_cpu_seconds=elapsed,
        conclusion='Nonzero instantaneous scalar correlation change, failing vacuum form; exactly zero quiet exterior Lambda change. A collective parent gravitational reduction remains unspecified.')
