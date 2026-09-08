#!/usr/bin/env python3
"""Reproduce covariant-measure obstructions and full-angular vacuum controls."""
import argparse
import hashlib
import inspect
import json
from pathlib import Path
import sys

import numpy as np
import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from recursive_horizons.nsc_covariant_measure import (  # noqa: E402
    antiperiodic_fermi_integral, cylinder_energy, cylinder_stress,
    finite_cutoff_cylinder_stress, finite_weyl_average, product_geometry_identities,
    spatial_heat_weyl_ratio, spin_difference, spin_difference_heat_integral,
    winding_theta_sums,
)
from check_nsc_scale_closure import compare  # noqa: E402

OUTPUT = ROOT / "results/nsc-4-covariant-measure.json"


def calculate():
    exact = product_geometry_identities()
    for name in ("euler_density", "rigid_four_dimensional_weyl_residual",
                 "local_weyl_weight_residual", "local_C_squared_radial_null_projection",
                 "local_C_squared_stress_trace"):
        assert exact[name] == "0"
    assert exact["weyl_squared"] == "4/(3*a**4)"
    average = finite_weyl_average()
    asymptote = average["asymptotic_integrand"]
    average["successive_mean_slopes"] = [
        (right["negative_orbit_integral"] - left["negative_orbit_integral"])
        / (right["extent"] - left["extent"])
        for left, right in zip(average["integrals"], average["integrals"][1:])]
    assert max(abs(v - asymptote) for v in average["successive_mean_slopes"]) < 1e-12

    counting = []
    # Independent explicit Hamiltonian blocks test the particle/antiparticle
    # factor before summing the field determinant or using winding formulas.
    for k in (1, 2, 3):
        p, radius = .3, 1.7
        block = np.array([[k / radius, p], [p, -k / radius]])
        values = np.linalg.eigvalsh(np.kron(np.eye(4 * k), block))
        actual = -.5 * np.abs(values).sum()
        expected = -4 * k * np.sqrt(p**2 + (k / radius)**2)
        assert abs(actual - expected) < 1e-12
        assert np.sum(values > 0) == np.sum(values < 0) == 4 * k
        counting.append({"k": k, "momentum": p, "sphere_radius": radius,
                         "positive_energy_multiplicity": int(np.sum(values > 0)),
                         "negative_energy_multiplicity": int(np.sum(values < 0)),
                         "minus_half_trace_absolute_H": float(actual),
                         "independent_mode_energy": float(expected)})
    heat_counting = [{"proper_time": t, "angular_max": 1024,
                      "rank_four_volume_weyl_ratio": spatial_heat_weyl_ratio(t)}
                     for t in (.1, .01, .001)]
    assert abs(heat_counting[-1]["rank_four_volume_weyl_ratio"] - 1) < 2e-4

    controls = []
    stress_controls = []
    for length in (1., 2., 4.):
        periodic, periodic_gradient = cylinder_energy(length, eta=0.)
        antiperiodic, antiperiodic_gradient = cylinder_energy(length, eta=.5)
        difference, gradient = spin_difference(length)
        heat, heat_error = spin_difference_heat_integral(length)
        assert abs(periodic - antiperiodic - difference) < 1e-12
        assert abs(periodic_gradient - antiperiodic_gradient - gradient) < 1e-12
        assert abs(difference - heat) < 2e-11
        assert difference > 0 and gradient < 0
        derivative_checks = []
        for step in (.002, .001, .0005):
            minus = spin_difference_heat_integral(length - step)[0]
            plus = spin_difference_heat_integral(length + step)[0]
            finite_difference = (plus - minus) / (2 * step)
            derivative_checks.append({"step": step, "heat_integral_centered_derivative": finite_difference,
                                      "absolute_derivative_error": abs(finite_difference - gradient)})
        assert derivative_checks[-1]["absolute_derivative_error"] < 4e-5
        controls.append({"length": length, "sphere_radius": 1., "angular_max": 64, "winding_max": 64,
                         "periodic_energy": periodic, "antiperiodic_energy": antiperiodic,
                         "periodic_spacing_derivative": periodic_gradient,
                         "antiperiodic_spacing_derivative": antiperiodic_gradient,
                         "spin_energy_difference": difference, "bessel_spacing_derivative": gradient,
                         "spin_pressure_difference": -gradient / (4 * np.pi),
                         "independent_winding_heat_integral": heat,
                         "quadrature_error_estimate_only": heat_error,
                         "independent_heat_derivative_checks": derivative_checks})
        periodic_stress = cylinder_stress(length, eta=0.)
        ap_stress = cylinder_stress(length, eta=.5)
        fermi = antiperiodic_fermi_integral(length)
        for key in ("energy", "length_derivative", "radius_derivative"):
            assert abs(ap_stress[key] - fermi[key]) < 3e-11
        assert periodic_stress["radial_null_stress"] > 0 > ap_stress["radial_null_stress"]
        assert fermi["energy"] < 0 < fermi["length_derivative"]
        for stress in (periodic_stress, ap_stress):
            assert abs(stress["relative_stress_trace"]) < 1e-12
            assert abs(stress["homogeneity_residual"]) < 1e-12
        difference_stress = {key: periodic_stress[key] - ap_stress[key] for key in (
            "relative_energy_density", "relative_radial_pressure", "relative_transverse_pressure",
            "radial_null_stress", "relative_stress_trace")}
        step = .0005
        independent_radius_derivative = (
            antiperiodic_fermi_integral(length, radius=1. + step)["energy"]
            - antiperiodic_fermi_integral(length, radius=1. - step)["energy"]) / (2 * step)
        assert abs(independent_radius_derivative - ap_stress["radius_derivative"]) < 3e-8
        stress_controls.append({"length": length, "sphere_radius": 1.,
                                "angular_max": 64, "winding_max": 64,
                                "periodic": periodic_stress, "antiperiodic": ap_stress,
                                "periodic_minus_antiperiodic": difference_stress,
                                "independent_AP_fermi_log_integral": fermi,
                                "radius_finite_difference_step": step,
                                "independent_AP_fermi_radius_finite_difference": independent_radius_derivative})

    reference_energy, reference_gradient = spin_difference(1.)
    convergence = {"length": 1., "sphere_radius": 1.,
                   "reference_angular_max": 64, "reference_winding_max": 64,
                   "angular": [], "winding": [], "tail_bound_proved": False}
    for kmax in (8, 16, 32):
        energy, gradient = spin_difference(1., angular_max=kmax)
        convergence["angular"].append({"angular_max": kmax, "winding_max": 64,
                                       "energy": energy, "spacing_derivative": gradient,
                                       "energy_difference_from_reference": abs(energy - reference_energy),
                                       "gradient_difference_from_reference": abs(gradient - reference_gradient)})
    for nmax in (3, 7, 15, 31):
        energy, gradient = spin_difference(1., winding_max=nmax)
        convergence["winding"].append({"angular_max": 64, "winding_max": nmax,
                                       "energy": energy, "spacing_derivative": gradient,
                                       "energy_difference_from_reference": abs(energy - reference_energy),
                                       "gradient_difference_from_reference": abs(gradient - reference_gradient)})
    assert convergence["angular"][-1]["gradient_difference_from_reference"] < 1e-8
    assert convergence["winding"][-1]["gradient_difference_from_reference"] < 1e-10
    dilation = 2.3
    scaled_energy, scaled_gradient = spin_difference(dilation, radius=dilation)
    unit_control = {"factor": dilation,
                    "energy_scaling_residual": scaled_energy * dilation - reference_energy,
                    "gradient_scaling_residual": scaled_gradient * dilation**2 - reference_gradient}
    assert max(abs(unit_control[k]) for k in ("energy_scaling_residual", "gradient_scaling_residual")) < 1e-12

    theta_controls = []
    with mp.workdps(100):
        for b_text in (".02", ".1", "1", "10"):
            b = float(b_text)
            q = mp.exp(-mp.mpf(b_text))
            exact_second = q * mp.diff(lambda z: mp.jtheta(4, 0, z), q) / 2
            zeroth, second = winding_theta_sums(b, .5)
            relative_error = float(abs(mp.mpf(second) / exact_second - 1))
            assert second < 0 and relative_error < 1e-12
            theta_controls.append({"b": b, "branch": "direct" if b >= 1 else "Poisson",
                                   "zeroth_winding_sum": zeroth, "second_winding_sum": second,
                                   "second_winding_sum_sign": -1,
                                   "log10_absolute_second_sum": float(np.log10(abs(second))),
                                   "independent_100_digit_theta_derivative": float(exact_second),
                                   "relative_second_sum_error": relative_error})
    cutoff_controls = []
    for length in (1., 2., 4.):
        for eta in (0., .5):
            continuum = cylinder_stress(length, eta=eta)
            previous_error = np.inf
            for cutoff in (1., 2., 4., 8., 16.):
                finite = finite_cutoff_cylinder_stress(length, cutoff=cutoff, eta=eta)
                null = finite["axial_null_stress_integral"]
                assert null > 0 if eta == 0 else null < 0
                assert abs(null - finite["axial_null_stress_from_energy_variation"]) < 1e-12
                assert abs(finite["cutoff_homogeneity_residual"]) < 1e-11
                error = abs(null - continuum["radial_null_stress"])
                assert error <= previous_error + 1e-12
                previous_error = error
                step = .0005
                minus = finite_cutoff_cylinder_stress(length - step, cutoff=cutoff, eta=eta)["energy"]
                plus = finite_cutoff_cylinder_stress(length + step, cutoff=cutoff, eta=eta)["energy"]
                derivative_fd = (plus - minus) / (2 * step)
                null_fd = finite["relative_energy_density"] - derivative_fd / (4 * np.pi)
                assert abs(null_fd - null) < 3e-6
                cutoff_controls.append({"length": length, "sphere_radius": 1., "eta": eta,
                                        "cutoff": cutoff, "cutoff_times_radius": cutoff,
                                        "angular_max": 64, "response": finite,
                                        "continuum_energy": continuum["energy"],
                                        "continuum_axial_null_stress": continuum["radial_null_stress"],
                                        "absolute_energy_error_from_continuum": abs(finite["energy"] - continuum["energy"]),
                                        "absolute_null_error_from_continuum": error,
                                        "fixed_cutoff_length_difference_step": step,
                                        "independent_length_derivative": derivative_fd,
                                        "independent_null_stress_from_energy_difference": null_fd})
            assert previous_error < 1e-11
    cutoff_angular = []
    for eta in (0., .5):
        reference = finite_cutoff_cylinder_stress(1., cutoff=4., eta=eta)
        for kmax in (8, 16, 32):
            value = finite_cutoff_cylinder_stress(1., cutoff=4., eta=eta, angular_max=kmax)
            cutoff_angular.append({"eta": eta, "length": 1., "sphere_radius": 1., "cutoff": 4.,
                                   "angular_max": kmax, "reference_angular_max": 64,
                                   "energy_error": abs(value["energy"] - reference["energy"]),
                                   "null_error": abs(value["axial_null_stress_integral"] - reference["axial_null_stress_integral"])})
        assert cutoff_angular[-1]["null_error"] < 1e-12

    return {
        "schema": "nsc-covariant-measure-v1", "artifact_id": "NSC-4-COVARIANT-MEASURE",
        "classification": "anomaly_insufficiency_and_finite_scale_average_obstruction_with_full_angular_vacuum_response_control",
        "source_hashes": {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in (
            "src/recursive_horizons/nsc_covariant_measure.py", "src/recursive_horizons/nsc_regulated.py",
            "scripts/check_nsc_covariant_measure.py")},
        "comparison_function_sha256": hashlib.sha256(inspect.getsource(compare).encode()).hexdigest(),
        "reproduction_tolerance": {"every_numeric_field_rtol": 1e-8, "every_numeric_field_atol": 1e-8,
                                   "other_fields_and_keys": "exact", "excluded_scientific_fields": []},
        "conventions": {"units": "hbar=c=1; L and a are lengths; E has inverse-length units; pressure has inverse-length^4 units",
                        "geometry": "R_t x S1_L x S2_a; Euclidean modulus uses R_tau x S1_L x S2_a",
                        "coordinate_metric": "diag(1,(L/(2*pi))^2,a^2,a^2*sin(theta)^2); x has period 2*pi",
                        "field": "one complex four-component massless Dirac field, no gauge field or inserted mass",
                        "state": "zero-temperature positive-frequency static ground state for each declared spatial spin structure",
                        "cylinder_regulator_limit": "UV-finite continuum one-loop compactification contribution; no finite physical cutoff is imposed in these Bessel, Gaussian, or Fermi-log controls",
                        "spin_structures": {"periodic_eta": 0., "antiperiodic_eta": .5},
                        "topology_and_direction_scope": "axial x is compact S1; this differs from the unwrapped periodic rho-throat array and from compactification of transverse carrier y; negative null stress here is along compact x",
                        "energy_subtraction": "zero winding / decompactified cylinder removed; spin difference compares the same local metric",
                        "anomaly_coefficient_added_to_computation": False},
        "shape_ambiguity": {
            "exact_coordinate_curvature_and_variation": exact,
            "unresolved_term": "alpha integral sqrt(g) C_abcd C^abcd in four dimensions",
            "alpha_role": "example of a finite Weyl-invariant coefficient not determined by anomaly matching; no value selected or fitted",
            "root_sensitivity_if_E0_LL_nonzero": "dL_star/dalpha=-16*pi/(3*a^2*E0_LL(L0)) at alpha=0 and fixed a",
            "interpretation": "anomaly data alone do not determine a spacing root; properly matched observables are unchanged by consistent scheme and coupling transformations",
            "demonstrates_ambiguity_of_fully_specified_uv_functional": False,
        },
        "flat_weyl_average": {
            "normalization": "Gamma=Tr E1(D^2/Lambda^2)/2+N[log(M/Lambda)+gamma_E/2]",
            "orbit": "D_phi=exp(-phi)D at fixed nonzero finite spectrum, N, Lambda, M",
            "negative_infinite_phi_limit_Gamma": "N[log(M/Lambda)+gamma_E/2]",
            "proof": "E1(e^(-2phi)d_j^2/Lambda^2) tends to zero for every nonzero fixed eigenvalue as phi tends to -infinity; exp(-Gamma) tends to a positive constant",
            "flat_integral_over_real_phi_converges": False,
            "numerical_control": average,
            "scope": "finite fixed-N prescribed determinant and flat measure, not a continuum no-go result or a derived gauge-volume prescription",
        },
        "determinant_counting": {
            "sphere_eigenvalues": "+/-k/a, k=1,2,...; 2k sphere eigenspinors per sign",
            "full_spatial_H_eigenvalues": "+/-sqrt(p^2+k^2/a^2), each with multiplicity 4k",
            "spatial_H_squared_heat_multiplicity": "8k per k and circle momentum",
            "vacuum_energy_per_k_and_momentum": "-Tr|H|/2=-4k*sqrt(p^2+k^2/a^2)",
            "matrix_counting_controls": counting, "rank_four_spatial_heat_controls": heat_counting,
        },
        "spin_structure_response": {
            "energy_formula": "E_eta=8/(pi*a) sum_k k^2 sum_n cos(2*pi*n*eta)*K1(n*k*L/a)/n",
            "difference_heat_integral": "DeltaE=4*L/pi integral_0^infinity dt/t^2 [sum_k k exp(-t*k^2/a^2)] [sum_positive_odd_n exp(-n^2*L^2/(4*t))]",
            "pressure_definition": "Delta p_x=-(partial DeltaE/partial L)/(4*pi*a^2), radius held fixed",
            "exact_sign": "DeltaE>0 and partial_L DeltaE<0 for all L,a>0; only positive odd windings survive and K1'(z)=-(K0(z)+K2(z))/2<0",
            "controls": controls, "independent_truncation_convergence": convergence,
            "natural_unit_scaling": unit_control,
            "local_counterterms_cancel_in_same_metric_spin_difference": True,
            "absolute_vacuum_stress_or_gravitational_spacing_fixed": False,
        },
        "radial_null_response": {
            "null_vector_convention": "orthonormal k^a=(1,1,0,0); T_ab k^a k^b=rho+p_x",
            "relative_stress_definitions": "rho=E/(4*pi*a^2*L); p_x=-E_L/(4*pi*a^2); p_perp=-E_a/(8*pi*a*L)",
            "AP_fermi_log": "E_AP=-(8/pi) sum_k k integral_0^infinity dp log(1+exp(-L*sqrt(p^2+k^2/a^2)))",
            "AP_exact_sign": "the Fermi logarithm is positive and its L derivative negative, so E_AP<0, E_AP,L>0 and rho_AP+p_x_AP<0 for every L,a>0",
            "periodic_exact_sign": "all periodic Bessel terms are positive with negative L derivative, so rho_P+p_x_P>0 for every L,a>0",
            "local_counterterm_cancellation": "every locally covariant metric counterterm tensor respects boosts in the locally flat (t,x) factor of this constant-radius product, so T_AB is proportional to its two-dimensional Minkowski metric and its radial-null projection vanishes",
            "decompactified_vacuum_cancellation": "the reference ground state on R_t x R_x x S2_a is also boost invariant in (t,x), hence its radial-null projection vanishes",
            "trace_identity": "L*E_L+a*E_a=-E implies -rho+p_x+2*p_perp=0 for the finite massless compactification contribution; the identical local anomaly cancels in the relative subtraction",
            "full_absolute_stress_tensor_fixed": False,
            "absolute_radial_null_projection_fixed_for_this_product_vacuum": True,
            "same_values_proved_at_finite_physical_spectral_cutoff": False,
            "controls": stress_controls,
            "physical_connection": "the same ordinary quantum Dirac field supplies the negative radial-null stress sign needed by a flaring throat in an Einstein-equation realization; magnitude and varying-radius self-consistency remain unsolved",
            "varying_radius_neck_stress_obtained_by_substituting_cylinder_values": False,
            "negative_kinetic_energy_field_added": False,
        },
        "finite_proper_time_cutoff": {
            "prescription": "same proper-time modulus lower limit Lambda^-2 after Euclidean frequency integration, with zero winding subtracted at the same cutoff; no unresolved complete anomaly/compensator is inserted",
            "energy_formula": "E_eta=2*L/pi integral_Lambda^-2^infinity dt/t^2 K_s(t) S0_eta(L^2/(4t)); K_s=sum_k k exp(-t*k^2/a^2)",
            "null_formula": "rho+p_x=L^2/(4*pi^2*a^2) integral_Lambda^-2^infinity dt/t^3 K_s(t) S2_eta(L^2/(4t))",
            "theta_definitions": "S0=sum_n>=1 cos(2*pi*n*eta) exp(-b*n^2); S2=sum_n>=1 n^2 cos(2*pi*n*eta) exp(-b*n^2)",
            "AP_sign_proof": "S2_AP=q*dtheta4(0,q)/dq/2<0 for q=exp(-b) in (0,1): theta4=product_m>=1 (1-q^(2m))*(1-q^(2m-1))^2 is positive and strictly decreasing; K_s and integration measure are positive",
            "negative_AP_compact_axial_null_sign_survives_every_positive_cutoff": True,
            "continuum_AP_length_derivative_positive_claim_extended_to_finite_cutoff": False,
            "theta_evaluation": "12 direct winding terms for b>=1; 8 Poisson-dual terms for b<1, including differentiated dual formula for S2; prevents small-b AP cancellation",
            "theta_control_precision_decimal_digits": 100,
            "theta_controls": theta_controls,
            "controls": cutoff_controls,
            "angular_convergence_at_largest_primary_cutoff": cutoff_angular,
            "metric_variations_hold_dimensional_cutoff_fixed": True,
            "finite_cutoff_scale_identity": "E+L*E_L+a*E_a-Lambda*E_Lambda=0; finite compactification trace need not vanish when Lambda is held fixed",
            "limits": "AP S2 can underflow to zero for very small b, so analytic sign is not inferred from floating-point zeros; cutoff 4 need not approximate continuum when Lambda*L is small; auxiliary cutoff 8 and 16 runs check approach to continuum",
            "identification_with_original_throat_radial_source": False,
            "model_wide_finite_cutoff_stationarity_or_compensator_derived": False,
        },
        "next_equation": {
            "geometry": "ds^2=-N(x)^2dt^2+q(x)^2dx^2+a(x)^2dOmega_2^2 with smooth periodic x",
            "variation": "delta Gamma_one/delta N=delta Gamma_one/delta q=delta Gamma_one/delta a=delta Gamma_one/delta Phi=0 before gauge fixing",
            "conservation_control": "p_x'+(N'/N)(rho+p_x)+2(a'/a)(p_x-p_perp)=0",
            "conditional_weyl_ward_identity": "2g_mn deltaGamma/delta g_mn-deltaGamma/delta phi=0 if delta g_mn=2sigma g_mn, delta phi=-sigma is an exact gauge symmetry",
            "required_owners": ["common covariant determinant and finite invariant remainder or ultraviolet measure",
                                "quantum state and dilaton measure or gauge volume",
                                "all angular sectors and varying-radius metric response",
                                "coupled unequal-cell metric and link variation to determine Omega"],
        },
        "primary_sources": [
            "https://arxiv.org/html/1106.3263v1",
            "https://arxiv.org/abs/0904.0612",
            "https://arxiv.org/html/gr-qc/9810032v2",
            "https://arxiv.org/abs/gr-qc/9505009",
            "https://dlmf.nist.gov/20.5#E4",
        ],
        "nonclaims": {"covariant_stationary_geometry_solved": False,
                      "motif_spacing_or_inheritance_ratio_predicted": False,
                      "full_varying_radius_matter_stress_computed": False,
                      "trapped_five_dimensional_state_or_action_derived": False,
                      "particle_species_or_child_cosmological_constant_predicted": False,
                      "new_to_world_physics_discovery_established": False},
        "terminal": True,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--check", action="store_true")
    group.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = calculate()
    if args.check:
        compare(json.loads(OUTPUT.read_text()), result)
        print("Covariant-measure and full-angular vacuum controls reproduced; every scientific field checked.")
    elif args.output:
        with args.output.open("x") as stream:
            json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write("\n")
        print(args.output)
    else:
        print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
