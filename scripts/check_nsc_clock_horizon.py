#!/usr/bin/env python3
"""Audit clock inheritance and observer geometry on the fixed PG benchmark."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
import sympy as sp

from check_nsc_scale_closure import compare

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results/nsc-5-clock-horizon.json"


def exact_zero(expression):
    residual = sp.simplify(expression)
    assert residual == 0, residual
    return str(residual)


def calculate():
    tau, rho = sp.symbols("tau rho", real=True)
    beta = sp.Function("beta")(rho)
    a_generic = 1 - beta**2
    metric = sp.Matrix([[a_generic, -beta], [-beta, -1]])
    inverse = metric.inv()
    coordinates = (tau, rho)
    # Angular acceleration components of the constant-angle horizon generator
    # vanish because g_tau_tau depends only on rho and g^{angular,rho}=0.
    connection = [sp.simplify(sum(inverse[i, j] * (
        2*sp.diff(metric[j, 0], tau) - sp.diff(metric[0, 0], coordinates[j]))/2
        for j in range(2))) for i in range(2)]
    connection_residuals = [
        exact_zero(connection[0] - beta*sp.diff(a_generic, rho)/2),
        exact_zero(connection[1] - a_generic*sp.diff(a_generic, rho)/2),
    ]
    inverse_residual = [[exact_zero(x) for x in row]
                        for row in (metric*inverse-sp.eye(2)).tolist()]
    assert metric.det() == -1
    ap = sp.symbols("A_prime_h", positive=True)
    horizon_connection = [sp.simplify(x.subs(sp.diff(beta, rho), -ap/2)
                                     .subs(beta, 1)) for x in connection]
    assert horizon_connection == [ap/2, 0]

    radius = sp.sqrt(1+rho**2)
    area = 4*sp.pi*radius**2
    expansions = [sp.simplify(sp.diff(area, rho)*(-beta+direction)/area)
                  for direction in (1, -1)]
    expansion_residuals = [exact_zero(theta - 2*rho*(-beta+direction)/(1+rho**2))
                          for theta, direction in zip(expansions, (1, -1))]
    null_residuals = [exact_zero((sp.Matrix([[1, -beta+direction]]) * metric
                                 * sp.Matrix([1, -beta+direction]))[0])
                      for direction in (1, -1)]
    product_residual = exact_zero(expansions[0]*expansions[1]
                                 + 4*a_generic*rho**2/(1+rho**2)**2)

    a = 1+3*rho+3*(1+rho**2)*(sp.atan(rho)-sp.pi/2)
    a_prime = sp.diff(a, rho)
    derivative_residual = exact_zero(a_prime-6*(1-rho*(sp.pi/2-sp.atan(rho))))
    a_numeric = sp.lambdify(rho, a, "numpy")
    ap_numeric = sp.lambdify(rho, a_prime, "numpy")
    horizon = float(brentq(a_numeric, 0., 8., xtol=5e-15))
    surface_gravity = float(ap_numeric(horizon)/2)
    assert abs(a_numeric(horizon)) < 1e-13 and surface_gravity > 0
    probes = []
    for point in (-8., -1., 0., 1., horizon, 8.):
        avalue = float(a_numeric(point))
        bvalue = float(np.sqrt(1-avalue))
        values = [float(2*point*(-bvalue+direction)/(1+point*point))
                  for direction in (1, -1)]
        probes.append({"rho": point, "A": avalue, "beta": bvalue,
                       "theta_plus": values[0], "theta_minus": values[1]})
    assert all(p["theta_plus"] > 0 and p["theta_minus"] > 0 for p in probes[:2])
    assert probes[2]["theta_plus"] == probes[2]["theta_minus"] == 0
    assert probes[3]["theta_plus"] < 0 and probes[3]["theta_minus"] < 0
    assert abs(probes[4]["theta_plus"]) < 1e-13 and probes[4]["theta_minus"] < 0
    assert probes[5]["theta_plus"] > 0 and probes[5]["theta_minus"] < 0

    kappa, lam = sp.symbols("kappa lambda", positive=True)
    affine = sp.exp(kappa*tau)
    inverse_affine = sp.log(lam)/kappa
    affine_residual = exact_zero(sp.diff(affine, tau, 2)-kappa*sp.diff(affine, tau))
    geodesic_residual = exact_zero(sp.diff(inverse_affine, lam, 2)
                                   + kappa*sp.diff(inverse_affine, lam)**2)
    assert sp.limit(affine, tau, -sp.oo) == 0
    affine_samples = [{"tau": t, "lambda_with_lambda_at_zero_one":
                       float(np.exp(surface_gravity*t))}
                      for t in (-100., -20., -4., 0., 4.)]

    # Independent exact asymptotic coefficients determine the interior rates.
    asymptotics = {"A_over_rho_squared": sp.limit(a/rho**2, rho, -sp.oo),
                   "A_prime_over_rho": sp.limit(a_prime/rho, rho, -sp.oo),
                   "areal_radius_squared_over_rho_squared":
                   sp.limit(radius**2/rho**2, rho, -sp.oo)}
    assert list(asymptotics.values()) == [-3*sp.pi, -6*sp.pi, 1]
    interior_probes = []
    for point in (0., -1., -4., -16., -64.):
        avalue, derivative = float(a_numeric(point)), float(ap_numeric(point))
        h_parallel = derivative/(2*np.sqrt(-avalue))
        h_sphere = -np.sqrt(-avalue)*point/(1+point*point)
        interior_probes.append({"rho": point, "H_parallel": float(h_parallel),
                                "H_sphere": float(h_sphere),
                                "shear_squared": float((h_parallel-h_sphere)**2/3)})
    assert interior_probes[0]["H_sphere"] == 0
    assert interior_probes[0]["shear_squared"] > 0
    assert interior_probes[-1]["shear_squared"] < interior_probes[-2]["shear_squared"]

    clock_examples = []
    for q, label, limit in ((sp.Rational(1, 2), "convergent", "2"),
                            (sp.Integer(1), "divergent", "+infinity"),
                            (sp.Integer(2), "divergent", "+infinity")):
        samples = [{"last_ancestor_index": n,
                    "total_ticks_exact": str(sum(q**m for m in range(n+1)))}
                   for n in (0, 1, 3, 7, 15)]
        n = sp.symbols("n", integer=True, nonnegative=True)
        total = n+1 if q == 1 else (1-q**(n+1))/(1-q)
        exact_limit = sp.limit(total, n, sp.oo)
        assert exact_limit == (sp.Integer(2) if q < 1 else sp.oo)
        clock_examples.append({"adjacent_clock_factor": str(q),
                               "each_local_interval_ticks": "1", "classification": label,
                               "infinite_sum": limit, "partial_sums": samples})

    return {
        "schema": "nsc-clock-horizon-v1", "artifact_id": "NSC-5-CLOCK-HORIZON",
        "classification": "exact_fixed_background_observer_map_and_clock_completeness_audit",
        "source_hashes": {name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                          for name in ("scripts/check_nsc_clock_horizon.py", "docs/nsc-clock-horizon.md",
                                       "scripts/check_nsc_dirac_tetrad.py")},
        "conventions": {"signature": "+---", "units": "c=L_throat=1",
                        "spacetime": "declared PG patch R_tau x R_rho x S²",
                        "metric": "d_tau²-(d_rho+beta*d_tau)²-r²*d_Omega²",
                        "radius": "sqrt(1+rho²)",
                        "A": "1+3rho+3(1+rho²)(atan(rho)-pi/2)", "beta": "sqrt(1-A)",
                        "null_normalization": "ell_plus/minus=partial_tau+(-beta+/-1)partial_rho"},
        "observer_map": {"sphere_area": "4pi(1+rho²)", "area_dimension": 2,
                         "spatial_slice_dimension": 3,
                         "null_residuals": null_residuals, "expansion_residuals": expansion_residuals,
                         "theta_plus": "2rho(1-beta)/(1+rho²)",
                         "theta_minus": "2rho(-1-beta)/(1+rho²)",
                         "product_identity_residual": product_residual,
                         "sign_regions": [
                             {"domain": "rho>rho_h", "plus": "+", "minus": "-", "type": "untrapped"},
                             {"domain": "rho=rho_h", "plus": "0", "minus": "-", "type": "outer_horizon"},
                             {"domain": "0<rho<rho_h", "plus": "-", "minus": "-", "type": "future_trapped"},
                             {"domain": "rho=0", "plus": "0", "minus": "0", "type": "minimal_areal_radius"},
                             {"domain": "rho<0", "plus": "+", "minus": "+", "type": "anti_trapped"}],
                         "analytic_sign_argument": "A_prime=6(1-rho*x)>0, x=pi/2-atan(rho); for rho>0 use atan(1/rho)<1/rho, while rho<=0 is immediate. A(0)<0 and A(+infinity)=1 imply a unique positive horizon. beta>0 follows from atan(y)>y/(1+y²).",
                         "A_prime_identity_residual": derivative_residual,
                         "geometry_probes": probes,
                         "horizon_rho": horizon,
                         "horizon_areal_radius": float(np.sqrt(1+horizon*horizon)),
                         "horizon_sphere_area": float(4*np.pi*(1+horizon*horizon)),
                         "parent_unresolved_point": "measurement limit; invariant area is unchanged",
                         "observer_cosmological_horizon": "requires the specified observer worldline and causal past/future; not identified with the parent horizon by coordinates"},
        "horizon_affine_parameter": {
            "radial_metric_determinant": "-1", "inverse_residual": inverse_residual,
            "connection_identity_residuals": connection_residuals,
            "Gamma_tau_tau_tau": "beta*A_prime/2", "Gamma_rho_tau_tau": "A*A_prime/2",
            "generator": "k=partial_tau at rho_h and fixed angles",
            "nonaffinity": "nabla_k k=(A_prime_h/2)k", "kappa_h": surface_gravity,
            "affine_parameter": "lambda=C*exp(kappa_h*tau)+lambda_0, C>0",
            "affine_equation_residual": affine_residual,
            "inverse_affine_geodesic_residual": geodesic_residual,
            "past_affine_endpoint": "lambda_0 at tau->-infinity",
            "finite_past_affine_interval_from_tau_zero_for_C_one": 1.,
            "normalized_samples": affine_samples,
            "declared_pg_patch_past_null_complete": False,
            "curvature_singularity_inferred": False,
            "extension_status": "an extension beyond this PG patch has not been constructed or tested here",
            "relation_to_dirac_flow": "coordinate-time characteristic completeness and norm-preserving evolution do not establish affine geodesic completeness"},
        "child_geometry": {
            "interior_proper_time": "dT=-d_rho/sqrt(-A) for A<0, with decreasing rho future-directed",
            "interior_metric": "dT²-(-A)dt²-r²d_Omega²; Kantowski-Sachs homogeneous slices R_t x S²",
            "directional_rates": {"H_parallel": "A_prime/(2sqrt(-A))", "H_sphere": "-sqrt(-A)rho/(1+rho²)"},
            "shear_squared": "(H_parallel-H_sphere)²/3",
            "asymptotic_coefficients": {name: str(value) for name, value in asymptotics.items()},
            "asymptotic_directional_H": "sqrt(3pi)", "asymptotic_Lambda": "9pi",
            "interior_probes": interior_probes,
            "exactly_isotropic_at_finite_rho": False,
            "interpretation": "asymptotically locally de Sitter; no claim of a global FLRW chart, compact spatial volume, or observer-independent horizon"},
        "clock_inheritance": {
            "indexing": "m=0 is local room; m=1,2,... are successive ancestors",
            "adjacent_map": "q_m=d(tau_hat_{m-1})/d(tau_hat_m)>0 for specified transported clocks",
            "local_assigned_ancestral_duration": "T_0=sum_{m=0}^infinity [product_{j=1}^m q_j] Delta(tau_hat_m)",
            "constant_interval_assumption": "q_j are constant over each mapped interval; otherwise integrate the pointwise composite clock map",
            "infinite_runtime_condition": "the positive series diverges; infinitely many rooms alone is insufficient",
            "exact_examples": clock_examples,
            "ruler_map": "d_ell_c=alpha*d_ell_p", "clock_map": "d_tau_c=beta_clock*d_tau_p",
            "dimensionless_parent_speed_in_child_units": "v=(alpha/beta_clock)*c_p/c_c",
            "pure_coordinate_rescaling": "v=1 when the null cone is unchanged",
            "crossing_duration_ratio_in_matched_units": "(L_p/c_p)/(L_c/c_c)=(L_p/L_c)/(c_p/c_c)",
            "larger_length_faster_speed_examples": [
                {"length_ratio": 2, "speed_ratio": 4, "duration_ratio_exact": "1/2"},
                {"length_ratio": 2, "speed_ratio": 2, "duration_ratio_exact": "1"},
                {"length_ratio": 4, "speed_ratio": 2, "duration_ratio_exact": "2"}],
            "local_proper_age": "elapsed proper time along a specified observer from a specified formation event remains meaningful",
            "clock_map_derived_from_operator": False,
            "global_age_status": "requires causal gluing and a clock convention; not ruled out by scale-dependent local clocks"},
        "entropy": {"protected_thesis_edited": False,
                    "working_direction": "retain original weak-gradient parent hypothesis",
                    "ordering_status": "unresolved until matter, gravitational/boundary and coarse-graining contributions define a generalized entropy"},
        "primary_sources": [
            {"url": "https://arxiv.org/abs/gr-qc/0611022", "use": "black-universe Kantowski-Sachs interior and de Sitter asymptotic prior art; does not derive a Dirac source"},
            {"url": "https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.96.251101", "use": "original regular phantom black-hole family; its matter source is not adopted here"}],
        "nonclaims": {"self_sourcing_stationary_solution": False,
                      "full_extension_constructed": False, "all_singularities_coordinate_artifacts": False,
                      "infinite_ancestral_proper_or_affine_runtime_proved": False,
                      "particle_identity_derived_from_external_localization": False,
                      "new_to_world_priority_established": False},
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
        print("Clock inheritance, null expansions and horizon affine-parameter audit reproduced; every field checked.")
    elif args.output:
        with args.output.open("x") as stream:
            json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write("\n")
        print(args.output)
    else:
        print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
