#!/usr/bin/env python3
"""Compute the spin-structure stress difference on the same smooth axial metric."""
import argparse
from dataclasses import replace
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from recursive_horizons.nsc_shape_response import (  # noqa: E402
    StaticAxialMetric, radial_reduction_identities, smooth_metric,
    spin_shape_difference, stress_summary,
)
from recursive_horizons.nsc_covariant_measure import cylinder_stress  # noqa: E402
from recursive_horizons.nsc_smooth_geometry import areal_radius, eight_pi_G_stress  # noqa: E402

OUTPUT = ROOT / "results/nsc-4-shape-response.json"
ABS_TOL = 2e-7
REL_TOL = 1e-6


def compare_record(expected, observed, pointer=""):
    """Every field is checked; integer/scope/provenance fields remain exact.

    SVD vacuum differences subtract extensive UV contributions, so absolute
    tolerance reflects their BLAS roundoff rather than floating-point equality.
    Runtime Ward, sign, finite-variation and convergence assertions are separate.
    """
    if isinstance(expected, bool) or isinstance(observed, bool):
        if expected is not observed:
            raise RuntimeError(f"boolean mismatch at {pointer}")
    elif isinstance(expected, int) and isinstance(observed, int):
        if expected != observed:
            raise RuntimeError(f"integer mismatch at {pointer}")
    elif isinstance(expected, (int, float)) and isinstance(observed, (int, float)):
        if not math.isclose(expected, observed, rel_tol=REL_TOL, abs_tol=ABS_TOL):
            raise RuntimeError(f"numeric mismatch at {pointer}")
    elif isinstance(expected, dict) and isinstance(observed, dict):
        if expected.keys() != observed.keys():
            raise RuntimeError(f"keys differ at {pointer}")
        for key in expected:
            compare_record(expected[key], observed[key], f"{pointer}/{key}")
    elif isinstance(expected, list) and isinstance(observed, list):
        if len(expected) != len(observed):
            raise RuntimeError(f"length differs at {pointer}")
        for i, (a, b) in enumerate(zip(expected, observed)):
            compare_record(a, b, f"{pointer}/{i}")
    elif expected != observed:
        raise RuntimeError(f"value mismatch at {pointer}")


def calculate():
    exact = radial_reduction_identities()
    assert exact["symmetric_derivative_residual"] == exact["measure_residual"] == "0"
    constant_controls = []
    for length, angular_max in ((1., 32), (2., 16), (4., 8)):
        periodic = cylinder_stress(length, eta=0.)
        ap = cylinder_stress(length, eta=.5)
        reference = {"energy": periodic["energy"] - ap["energy"],
                     "rho": periodic["relative_energy_density"] - ap["relative_energy_density"],
                     "p_x": periodic["relative_radial_pressure"] - ap["relative_radial_pressure"],
                     "p_perp": periodic["relative_transverse_pressure"] - ap["relative_transverse_pressure"],
                     "null": periodic["radial_null_stress"] - ap["radial_null_stress"]}
        grids = []
        for points in (64, 128, 256):
            metric = StaticAxialMetric(length, np.ones(points), np.ones(points), np.ones(points))
            response = spin_shape_difference(metric, angular_max)
            mean = {"energy": response["energy_difference"],
                    "rho": float(np.mean(response["density_difference"])),
                    "p_x": float(np.mean(response["axial_pressure_difference"])),
                    "p_perp": float(np.mean(response["transverse_pressure_difference"])),
                    "null": float(np.mean(response["axial_null_difference"]))}
            assert np.ptp(response["axial_null_difference"]) < 1e-6
            grids.append({"points": points, "angular_max": angular_max,
                          "mean_response": mean,
                          "absolute_errors_from_Bessel_reference": {k: abs(mean[k] - reference[k]) for k in mean},
                          "maximum_spatial_null_variation": float(np.ptp(response["axial_null_difference"]))})
        for key in ("energy", "rho", "p_x", "p_perp", "null"):
            errors = [g["absolute_errors_from_Bessel_reference"][key] for g in grids]
            assert 3.5 < errors[0] / errors[1] < 4.5
            assert 3.5 < errors[1] / errors[2] < 4.5
        constant_controls.append({"length": length, "sphere_radius": 1., "Bessel_reference": reference,
                                  "grids": grids})

    base = smooth_metric(2., 64)
    base = replace(base, lapse=1 + .08 * np.cos(np.pi * base.x / 2),
                   radial_scale=1 + .06 * np.sin(np.pi * base.x / 2))
    response = spin_shape_difference(base, 4)
    variations = []
    for field, gradient in (("lapse", "gradient_lapse"), ("radial_scale", "gradient_radial_scale"),
                            ("sphere_radius", "gradient_sphere_radius")):
        for kind in ("smooth", "single_throat_node"):
            direction = (.4 * np.cos(np.pi * base.x / 2) + .2 * np.sin(np.pi * base.x)
                         if kind == "smooth" else np.eye(base.points)[base.points // 2])
            hf = float(np.dot(response[gradient], direction))
            differences = []
            for step in (.004, .002, .001):
                plus = replace(base, **{field: getattr(base, field) + step * direction})
                minus = replace(base, **{field: getattr(base, field) - step * direction})
                independent = (spin_shape_difference(plus, 4, False)["energy_difference"]
                               - spin_shape_difference(minus, 4, False)["energy_difference"]) / (2 * step)
                differences.append({"step": step, "centered_energy_derivative": independent,
                                    "absolute_error_from_polar_HF": abs(independent - hf)})
            assert differences[-1]["absolute_error_from_polar_HF"] < 2e-7
            variations.append({"field": field, "direction": kind, "polar_HF_derivative": hf,
                               "finite_differences": differences})

    families = []
    for radius_parameter in (2., 4.):
        grids = []
        for points in (64, 128, 256, 512):
            metric = smooth_metric(radius_parameter, points)
            np.testing.assert_allclose(metric.sphere_radius, areal_radius(metric.x, radius=radius_parameter), atol=1e-14)
            response = spin_shape_difference(metric, 8)
            summary = stress_summary(metric, response)
            assert summary["maximum_absolute_local_trace"] < 1e-11
            assert summary["maximum_absolute_local_weyl_identity"] < 1e-11
            assert abs(summary["global_lapse_homogeneity_residual"]) < 2e-7
            assert np.min(response["axial_null_difference"]) > 0
            grids.append(summary)
        for left, right in zip(grids, grids[1:]):
            ratio = left["rms_conservation_residual"] / right["rms_conservation_residual"]
            assert 3.4 < ratio < 4.6
        for key in ("energy_difference", "throat_axial_null_difference", "proper_axial_integral_null_difference"):
            changes = [abs(b[key] - a[key]) for a, b in zip(grids, grids[1:])]
            assert 3.4 < changes[0] / changes[1] < 4.6
            assert 3.4 < changes[1] / changes[2] < 4.6
        angular = []
        for maximum in (2, 4, 6, 8, 10):
            angular_metric = smooth_metric(radius_parameter, 256)
            angular_response = spin_shape_difference(angular_metric, maximum)
            angular.append({"angular_max": maximum,
                            "summary": stress_summary(angular_metric, angular_response),
                            "last_sector": angular_response["angular_sectors"][-1]})
        assert abs(angular[-1]["summary"]["throat_axial_null_difference"]
                   - angular[-2]["summary"]["throat_axial_null_difference"]) < 1e-7
        # Retain every finest-grid value alongside the independent geometric
        # Einstein residual; these are different quantities and units.
        samples = []
        required = eight_pi_G_stress(metric.x, radius=radius_parameter)
        for i in range(metric.points):
            samples.append({"x": float(metric.x[i]), "sphere_radius": float(metric.sphere_radius[i]),
                            "rho_difference": float(response["density_difference"][i]),
                            "p_x_difference": float(response["axial_pressure_difference"][i]),
                            "p_perp_difference": float(response["transverse_pressure_difference"][i]),
                            "null_difference": float(response["axial_null_difference"][i]),
                            "required_eight_pi_G_rho": float(required["eight_pi_G_rho"][i]),
                            "required_eight_pi_G_p_x": float(required["eight_pi_G_p_r"][i]),
                            "required_eight_pi_G_p_perp": float(required["eight_pi_G_p_t"][i]),
                            "required_eight_pi_G_null": float(required["eight_pi_G_rho_plus_p_r"][i])})
        families.append({"R": radius_parameter, "throat_radius": 1., "coordinate_period": 2 * radius_parameter,
                         "angular_max_for_radial_convergence": 8,
                         "radial_grids": grids, "angular_convergence_at_256_points": angular,
                         "finest_grid_angular_sectors": response["angular_sectors"],
                         "finest_grid_profile_sample_stride": 1, "profile_samples": samples,
                         "required_Einstein_null_comparison": {
                             "throat_eight_pi_G_null": float(required["eight_pi_G_rho_plus_p_r"][metric.points // 2]),
                             "proper_axial_integral_eight_pi_G_null": float(metric.spacing * np.sum(required["eight_pi_G_rho_plus_p_r"])),
                             "AP_minus_P_throat_null_change": -grids[-1]["throat_axial_null_difference"],
                             "interpretation": "AP lowers the null source relative to P, in the direction needed at the neck; Delta(P-AP) is not the absolute AP source, and G has not been fitted to equate magnitudes"},
                         "infinite_angular_or_spatial_error_bound_proved": False,
                         "tail_interpretation": "higher angular contributions reach subtractive SVD roundoff; their sign is not promoted into an exact tail result"})
    return {
        "schema": "nsc-shape-response-v1", "artifact_id": "NSC-4-SHAPE-RESPONSE",
        "classification": "same_smooth_varying_radius_full_angular_spin_structure_energy_and_local_stress_difference",
        "source_hashes": {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in (
            "src/recursive_horizons/nsc_shape_response.py", "scripts/check_nsc_shape_response.py",
            "src/recursive_horizons/nsc_covariant_measure.py", "src/recursive_horizons/nsc_smooth_geometry.py")},
        "authenticated_input_hashes": {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in (
            "results/nsc-4-covariant-measure.json", "results/nsc-4-smooth-geometry.json")},
        "reproduction_tolerance": {"float_atol": ABS_TOL, "float_rtol": REL_TOL,
                                   "reason": "subtractive SVD vacuum energies and local gradients; accommodates numerical BLAS variation",
                                   "integers_booleans_strings_keys": "exact", "excluded_fields": []},
        "conventions": {
            "signature": "+---", "units": "hbar=c=1; r_min=1 is the declared throat-length unit",
            "metric": "ds^2=N(x)^2dt^2-q(x)^2dx^2-r(x)^2dOmega_2^2",
            "state": "zero-temperature static ground state separately for each P/AP spin structure",
            "field": "one complex four-component massless Dirac field; both sphere signs included",
            "topology": "one smooth axial motif closed on S1; x is compact; different from unwrapped periodic array and compact transverse carrier y",
            "observable": "Delta=periodic minus antiperiodic on exactly the same local metric",
            "regulator": "spatial lattice and angular truncation removed by independent convergence checks; not a physical covariant finite Lambda for nonconstant lapse",
            "local_counterterms": "same local geometric terms and anomaly cancel in this spin-structure difference; absolute AP stress remains unspecified"},
        "radial_reduction": exact,
        "operator": {
            "continuum": "H=-i sigma2[c partial_x+c'/2]+v sigma1; c=N/q, v=N*kappa/r",
            "staggered_left": "A[e,i]=-sqrt(c_e*c_i)/h+(v_i+v_j)/4",
            "staggered_right": "A[e,j]=[sqrt(c_e*c_j)/h+(v_i+v_j)/4]*seam_phase; j=i+1 mod n",
            "edge_coefficient": "c_e=(c_i+c_j)/2", "seam_phases": [1, -1],
            "energy_counting": "E_eta=-4 sum_kappa>=1 kappa sum singular_values(A_kappa,eta)",
            "HF_identity": "d sum singular(A)=Re Tr[(U V_dagger)_polar_dagger dA]; vary every N_i,q_i,r_i before setting N=q=1",
            "lattice_resolution_guard": "max_i h*kappa*q_i/r_i < 2"},
        "stress_definitions": {
            "rho": "g_N/(4*pi*q*r^2*h)", "p_x": "-g_q/(4*pi*N*r^2*h)",
            "p_perp": "-g_r/(8*pi*N*q*r*h)",
            "trace": "Delta rho-Delta p_x-2 Delta p_perp=0",
            "local_Weyl_identity": "N_i g_N_i+q_i g_q_i+r_i g_r_i=0",
            "conservation": "p_x'+(N'/N)(rho+p_x)+2(r'/r)(p_x-p_perp)=0; centered derivative used as an independent continuum residual",
            "homogeneity": "sum_i N_i g_N_i=Delta E"},
        "constant_radius_Bessel_controls": constant_controls,
        "finite_variation_control": {
            "R": 2., "points": 64, "angular_max": 4,
            "background_lapse": "1+0.08*cos(pi*x/2)", "background_radial_scale": "1+0.06*sin(pi*x/2)",
            "background_radius": "sqrt(1+(sin(pi*x/4)/(pi/4))^2)",
            "smooth_direction": "0.4*cos(pi*x/2)+0.2*sin(pi*x)",
            "variations": variations},
        "smooth_families": families,
        "physical_connection": "changing axial spin structure P to AP lowers the computed throat null source on the actual varying-radius profile by Delta(rho+p_x); no local-density substitution from the cylinder is used",
        "smallest_absolute_stress_dependency": "derive the complete 4D covariant proper-time determinant and adopted compensator's lapse/radial/sphere variations in the chosen AP state on this same metric; the varying-radius local curvature terms no longer have the product cylinder's null cancellation, so their action-defined coefficients must be retained rather than set by fitting a neck residual",
        "nonclaims": {"absolute_AP_stress_negative_on_this_profile_established": False,
                      "metric_Einstein_residual_matched": False, "stationary_geometry_or_R_Omega_solved": False,
                      "unwrapped_periodic_array_stress_computed": False,
                      "physical_covariant_finite_cutoff_response_for_varying_lapse_computed": False,
                      "five_dimensional_carrier_or_particle_species_derived": False,
                      "new_to_world_Casimir_mechanism_discovered": False},
        "terminal": True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--check", action="store_true")
    group.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = calculate()
    if args.check:
        compare_record(json.loads(OUTPUT.read_text()), result)
        print("Smooth-metric spin-structure stress difference reproduced; every field checked.")
    elif args.output:
        with args.output.open("x") as stream:
            json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write("\n")
        print(args.output)
    else:
        print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
