#!/usr/bin/env python3
"""Reproduce the finite local metric-response identifiability calculation."""
import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from recursive_horizons.nsc_finite_terms import (  # noqa: E402
    BASIS, FIELDS, StaticAxialMetric, direct_curvature_identities, exact_identities,
    exact_throat_sensitivities, finite_response, general_control_metric,
    independent_energy, response_matrix, smooth_metric, ward_controls,
)

OUTPUT = ROOT / "results/nsc-8-finite-terms.json"
ABS_TOL, REL_TOL = 2e-6, 2e-8


def compare_record(expected, actual, path="$"):
    """Compare every scientific value, scope flag, expression and source hash."""
    if isinstance(expected, dict):
        if not isinstance(actual, dict) or expected.keys() != actual.keys():
            raise RuntimeError(f"record keys differ at {path}")
        for key in expected:
            compare_record(expected[key], actual[key], f"{path}.{key}")
    elif isinstance(expected, list):
        if not isinstance(actual, list) or len(expected) != len(actual):
            raise RuntimeError(f"record list differs at {path}")
        for i, (left, right) in enumerate(zip(expected, actual)):
            compare_record(left, right, f"{path}[{i}]")
    elif isinstance(expected, float):
        if isinstance(actual, bool) or not isinstance(actual, (float, int)) or not np.isclose(expected, actual, atol=ABS_TOL, rtol=REL_TOL):
            raise RuntimeError(f"numeric record differs at {path}: {expected!r} != {actual!r}")
    elif type(expected) is not type(actual) or expected != actual:
        raise RuntimeError(f"exact record differs at {path}")


def _finite_variations(metric, response):
    records = []
    angle = 2*np.pi*metric.x/metric.length
    for j, field in enumerate(FIELDS):
        # Two distinct smooth perturbations per retained metric component.
        for label, f in (("constant", np.ones(metric.points)),
                         ("mixed_harmonics", .4*np.cos(angle)+.2*np.sin(2*angle))):
            direction = getattr(metric, field)*f
            analytic = metric.spacing*response["gradients"][:, j]@direction
            differences = []
            for step in (.002, .001, .0005):
                plus = replace(metric, **{field: getattr(metric, field)+step*direction})
                minus = replace(metric, **{field: getattr(metric, field)-step*direction})
                observed = (independent_energy(plus)-independent_energy(minus))/(2*step)
                differences.append({"step": step, "centered_variation_by_basis": observed.tolist(),
                                    "absolute_error_by_basis": abs(observed-analytic).tolist()})
            np.testing.assert_allclose(observed, analytic, atol=3e-4, rtol=3e-6)
            records.append({"field": field, "direction": label,
                            "analytic_variation_by_basis": analytic.tolist(), "independent_controls": differences})
    return records


def _unit_and_weyl_controls(metric, response):
    factor = 1.7
    # Physical spatial dilation at fixed clock: energies scale as s^(3,1,-1,...).
    dilated = replace(metric, length=metric.length*factor, sphere_radius=metric.sphere_radius*factor)
    powers = np.array([3, 1, -1, -1, -1, -1])
    expected = response["energy"]*factor**powers
    observed = finite_response(dilated)["energy"]
    np.testing.assert_allclose(observed, expected, atol=1e-10)
    # Changing all units also changes M inversely: every energy then scales 1/s.
    units = finite_response(dilated, normalization_mass=1/factor)["energy"]
    np.testing.assert_allclose(units, response["energy"]/factor, atol=1e-10)
    angle = 2*np.pi*metric.x/metric.length
    omega = np.exp(.11*np.cos(angle)+.03*np.sin(2*angle))
    conformal = replace(metric, **{field: getattr(metric, field)*omega for field in FIELDS})
    changed = finite_response(conformal)
    np.testing.assert_allclose(changed["energy"][2], response["energy"][2], atol=1e-9)
    return {"spatial_dilation_factor": factor, "energy_scaling_powers_at_fixed_M": powers.tolist(),
            "spatial_dilation_energy_by_basis": observed.tolist(), "expected_spatial_dilation_energy_by_basis": expected.tolist(),
            "joint_unit_change_energy_by_basis": units.tolist(),
            "finite_local_weyl_factor": "exp(0.11*cos(2*pi*x/L)+0.03*sin(4*pi*x/L))",
            "C2_energy_before": float(response["energy"][2]), "C2_energy_after": float(changed["energy"][2]),
            "C2_maximum_pointwise_density_change": float(np.max(abs(changed["densities"][2]-response["densities"][2])))}


def _local_node_variations():
    """Direct local lattice-gradient control, separate from continuum probes."""
    metric = general_control_metric(64)
    response = finite_response(metric)
    direction = np.eye(metric.points)[metric.points//2]
    records = []
    for j, field in enumerate(FIELDS):
        analytic = metric.spacing*response["gradients"][:, j]@direction
        controls = []
        for step in (1e-5, 5e-6, 1e-6):
            plus = replace(metric, **{field: getattr(metric, field)+step*direction})
            minus = replace(metric, **{field: getattr(metric, field)-step*direction})
            observed = (independent_energy(plus)-independent_energy(minus))/(2*step)
            controls.append({"step": step, "independent_energy_variation_by_basis": observed.tolist(),
                             "absolute_error_by_basis": abs(observed-analytic).tolist()})
        np.testing.assert_allclose(observed, analytic, atol=5e-6, rtol=0)
        records.append({"field": field, "analytic_node_gradient_by_basis": analytic.tolist(),
                        "finite_variations": controls})
    return {"points": 64, "node": 32, "R": 2., "metric": "same analytic general N,q profile",
            "scope": "local derivative of discretized functional, not a nonsmooth continuum metric deformation",
            "variations": records}


def _samples(metric, response):
    return [{"x": float(metric.x[j]), "N": float(metric.lapse[j]),
             "q": float(metric.radial_scale[j]), "r": float(metric.sphere_radius[j]),
             "R": float(response["R"][j]),
             "Euler_Lagrange_density_by_basis_and_field": response["gradients"][:, :, j].tolist(),
             "rho_sensitivity_by_basis": response["rho"][:, j].tolist(),
             "p_x_sensitivity_by_basis": response["p_x"][:, j].tolist(),
             "p_perp_sensitivity_by_basis": response["p_perp"][:, j].tolist(),
             "null_sensitivity_by_basis": response["null"][:, j].tolist()}
            for j in range(metric.points)]


def calculate():
    exact = exact_identities()
    for key in ("weyl_contraction_residual", "euler_contraction_residual", "euler_primitive_residual", "pointwise_C2_local_weyl_residual"):
        assert exact[key] == "0"
    assert set(exact["rigid_weyl_density_residuals"]) == {"0"}
    direct = direct_curvature_identities()
    assert set(direct.values()) == {"0"}
    throat = exact_throat_sensitivities()
    cylinder = StaticAxialMetric(4., np.ones(64), np.ones(64), np.ones(64))
    cylinder_response = finite_response(cylinder)
    cylinder_matrix = response_matrix(cylinder, cylinder_response)
    assert cylinder_matrix["rank"] == 2
    cylinder_nulls = [[1, 1, 0, .25, 0, 0], [0, 0, 1, -1/3, 0, 0],
                      [0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 1]]
    cylinder_null_residuals = np.asarray(cylinder_matrix["matrix"])@np.asarray(cylinder_nulls).T
    assert np.max(abs(cylinder_null_residuals)) < 1e-11
    np.testing.assert_allclose(cylinder_response["energy"][2], 16*np.pi*4/3, atol=1e-12)
    np.testing.assert_allclose(cylinder_response["null"], 0., atol=1e-12)

    families = []
    for radius in (2., 4.):
        grids = []
        for points in (32, 64, 128):
            metric = smooth_metric(radius, points)
            response = finite_response(metric)
            matrix = response_matrix(metric, response)
            if points >= 64:
                assert matrix["rank"] == 4
            grids.append({"points": points, "energy_by_basis": response["energy"].tolist(),
                          "throat_null_sensitivity_by_basis": response["null"][:, points//2].tolist(),
                          "response_matrix": matrix, "ward_controls": ward_controls(metric, response)})
        previous = grids[-2]
        np.testing.assert_allclose(previous["energy_by_basis"], response["energy"], atol=2e-8)
        np.testing.assert_allclose(previous["response_matrix"]["matrix"], matrix["matrix"], atol=2e-7)
        coefficients = {sp.Symbol("a", positive=True): 1., sp.Symbol("k", positive=True): np.pi/(2*radius),
                        sp.Symbol("M", positive=True): 1.}
        expected_null = [float(sp.sympify(throat[name]["null"], locals={str(s): s for s in coefficients}).subs(coefficients)) for name in BASIS]
        np.testing.assert_allclose(response["null"][:, metric.points//2], expected_null, atol=1e-7)
        assert np.max(abs(np.array(grids[-1]["ward_controls"]["integrated_radial_gauge_response_by_basis"]))) < 1e-10
        families.append({"R": radius, "a": 1., "M": 1., "coordinate_period": metric.length,
                         "radial_convergence": grids, "exact_throat_null_sensitivities_evaluated": expected_null,
                         "sample_stride": 1, "samples": _samples(metric, response),
                         "proper_integral_null_sensitivity_by_basis": (metric.spacing*np.sum(response["null"]*metric.radial_scale, axis=1)).tolist()})

    general = general_control_metric(128)
    general_response = finite_response(general)
    general_matrix = response_matrix(general, general_response)
    assert general_matrix["rank"] == 4
    wards = ward_controls(general, general_response)
    assert max(wards["maximum_local_radial_diffeomorphism_residual_by_basis"]) < 1e-6
    assert max(wards["maximum_local_weyl_residual_by_basis"]) < 2e-6
    variation_controls = _finite_variations(general, general_response)

    source_paths = ("src/recursive_horizons/nsc_finite_terms.py", "scripts/check_nsc_finite_terms.py",
                    "tests/test_nsc_finite_terms.py", "docs/nsc-finite-terms.md",
                    "src/recursive_horizons/nsc_shape_response.py")
    input_paths = ("docs/nsc-covariant-measure.md", "docs/nsc-self-sourcing-progress.md",
                   "docs/nsc-smooth-geometry.md", "results/nsc-4-covariant-measure.json",
                   "results/nsc-4-shape-response.json", "results/nsc-4-smooth-geometry.json")
    return {
        "schema": "nsc-finite-terms-v1", "artifact_id": "NSC-8-FINITE-TERMS",
        "classification": "computed_local_covariant_metric_response_and_finite_coefficient_identifiability",
        "source_hashes": {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in source_paths},
        "authenticated_input_hashes": {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in input_paths},
        "reproduction_tolerance": {"float_atol": ABS_TOL, "float_rtol": REL_TOL,
                                   "reason": "local fourth-derivative roundoff; integrated sensitivities separately checked at 2e-7 or tighter",
                                   "integers_booleans_strings_keys": "exact", "excluded_fields": []},
        "conventions": {"signature": "+---", "units": "hbar=c=1; coordinates x,t carry length; N,q dimensionless; r has length; M has inverse length",
                        "metric": "ds^2=N(x)^2dt^2-q(x)^2dx^2-r(x)^2dOmega_2^2; x~x+L",
                        "Riemann": "R^a_bcd=partial_c Gamma^a_db-partial_d Gamma^a_cb+Gamma^a_ce Gamma^e_db-Gamma^a_de Gamma^e_cb",
                        "functional": "F_fin=sum_i c_i F_i; F_i=4*pi*int_0^L N*q*r^2 I_i dx; static S_fin=-int dt F_fin",
                        "curvature_sign": "R_L=-2/a^2 on constant cylinder; Wick-rotated positive metric has R_E=-R_L. Squared invariants agree; the M2R coefficient includes this sign convention.",
                        "basis": list(BASIS), "coefficients": ["c_0", "c_R", "c_C", "c_R2", "c_E", "c_box"],
                        "M": "normalization mass kept fixed through physical variations; M=1 chooses the declared length unit and is not an inserted fermion mass",
                        "domain": "torsionless 4D local metric sector through four derivatives on a closed smooth static axial cell; fixed topology and periodic variations",
                        "coefficient_meaning": "unresolved finite part of the same covariant operator/measure; no independent Einstein or counterterm action is added or fitted",
                        "state": "local terms are state independent; no absolute P/AP quantum expectation value is computed",
                        "numerics": "periodic real Fourier derivatives with zero Nyquist mode and trapezoidal integration; symbolic local jet gradients use exact discrete adjoints"},
        "exact_identities": exact, "independent_full_metric_Christoffel_contractions": direct,
        "variation_equations": {
            "Euler_Lagrange": "e_{if}=4*pi*(partial L_i/partial f - partial_x(partial L_i/partial f') + partial_x^2(partial L_i/partial f'')); L_i=N*q*r^2 I_i",
            "metric_residuals": "R_f(x)=R_f^(specified nonlocal/state/link part)(x)+sum_{i=0,R,C,R2} c_i e_{if}(x)=0 for f=N,q,r; include link equation from the same action",
            "stress": "rho_i=e_iN/(4*pi*q*r^2); p_x_i=-e_iq/(4*pi*N*r^2); p_perp_i=-e_ir/(8*pi*N*q*r)",
            "radial_coordinate_identity": "N'*e_N+r'*e_r-q*partial_x(e_q)=0; do not divide by r', which vanishes at the throat",
            "gauge_variation": "delta N=xi*N', delta q=xi*q'+q*xi', delta r=xi*r'; a row-direction identity, not an unidentifiable coefficient",
            "weyl_variation": "N*e_N+q*e_q+r*e_r=(4*density_M4,2*density_M2R,0,-12*density_boxR,0,0), where densities include 4*pi",
            "probe_equations": "b_A+sum_i J_Ai*c_i=0; J_Ai=int dx sum_f delta_A f*e_if; b_A must be predicted independently by the specified nonlocal/state/link functional",
            "implicit_stationary_sensitivity": "If an independently completed functional has a nonsingular physical Hessian H_AB at a stationary metric, du_star^A/dc_i=-(H^-1)^AB*J_Bi. H and u_star are not supplied by this calculation."},
        "cylinder_control": {"L": 4., "a": 1., "M": 1., "C2_energy": float(cylinder_response["energy"][2]),
                             "known_C2_energy": "16*pi*alpha*L/(3*a^2)", "response_matrix": cylinder_matrix,
                             "right_null_vectors_at_a_equals_M_equals_1": cylinder_nulls,
                             "maximum_right_null_residual": float(np.max(abs(cylinder_null_residuals))),
                             "interpretation": "rank two is an accidental constant-product probe degeneracy; it is not the general four-dimensional bulk coupling count"},
        "exact_throat_stress_sensitivities": throat,
        "smooth_profile": "r(x)=sqrt(a^2+(sin(k*x)/k)^2), k=pi/(2*R), x in [-R,R); N=q=1 only after variation",
        "smooth_families": families,
        "nonconstant_lapse_and_radial_control": {
            "R": 2., "points": 128, "N": "exp(0.08*cos(2*pi*x/L)+0.03*sin(4*pi*x/L))",
            "q": "exp(0.06*sin(2*pi*x/L)-0.02*cos(4*pi*x/L))",
            "energy_by_basis": general_response["energy"].tolist(), "response_matrix": general_matrix,
            "ward_controls": wards, "finite_metric_variations": variation_controls,
            "independent_local_node_variations": _local_node_variations(),
            "finite_unit_and_Weyl_controls": _unit_and_weyl_controls(general, general_response)},
        "identifiability": {
            "rank_on_each_declared_smooth_profile": 4, "bulk_coefficients_distinguishable": ["c_0", "c_R", "c_C", "c_R2"],
            "exact_bulk_right_null_vectors": [[0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 1]],
            "right_null_classification": "c_E is a fixed-topology Euler term and c_box a total derivative; both are boundary-only on this closed cell, not gauge couplings",
            "coarse_grid_caveat": "At R=4 with 32 points the raw discretized Euler column has a 3.41e-8 spurious singular value above the declared 1e-8 threshold. It falls below 1e-13 at 64 and 128 points; the exact continuum boundary identity fixes the physical nullspace.",
            "gauge_classification": "radial diffeomorphisms annihilate the functional for every coefficient. They relate metric residuals; they do not remove any of the four bulk columns.",
            "Weyl_classification": "C2 is locally Weyl invariant but has independent lapse/shape/null responses. Weyl invariance alone cannot fix c_C. Weyl is not declared a gauge symmetry of the incomplete full action.",
            "coefficients_fixed_by_this_calculation": [],
            "conditions_still_required": ["Specify the covariant regulated determinant and its normalization, finite subtraction/renormalization conditions and physical state on the same metric.",
                                           "Derive four independent finite bulk matching data from that same UV/operator measure, or supply four independent renormalization conditions. Rank four proves sensitivity, not their values.",
                                           "Specify compensator/link variations and any boundary/topology changes. c_E and c_box may then require boundary data despite vanishing here.",
                                           "Only after those inputs, solve all independent metric and link residuals and the state consistently; keep the Einstein coefficient counted exactly once."]},
        "conclusion": "The actual smooth-profile response distinguishes all four local bulk coefficients. Their unresolved values affect the independent absolute-stress equations, including the neck null equation; no stationary geometry follows from the already computed spin-structure difference or anomaly alone.",
        "nonclaims": {"complete_covariant_UV_measure_constructed": False, "absolute_quantum_stress_computed": False,
                      "stationary_metric_or_inheritance_scale_solved": False, "counterterm_fitted_to_neck": False,
                      "independent_gravity_action_added": False, "in_in_measure_constructed": False,
                      "global_infinite_tower_closed": False, "new_to_world_curvature_variation_discovered": False},
        "terminal": True,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--check", action="store_true")
    group.add_argument("--output", type=Path)
    args = parser.parse_args()
    # Refuse an existing output before the computation, and exclusively create
    # it afterward too, closing the ordinary check/create race.
    if args.output and args.output.exists():
        raise FileExistsError(args.output)
    result = calculate()
    if args.check:
        compare_record(json.loads(OUTPUT.read_text()), result)
        print("Finite local metric identifiability reproduced; all fields and hashes checked.")
    elif args.output:
        with args.output.open("x") as stream:
            json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write("\n")
        print(args.output)
    else:
        print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
