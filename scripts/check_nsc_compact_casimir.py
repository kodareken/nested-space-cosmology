#!/usr/bin/env python3
"""Reproduce only the new curved compact interaction and holonomy calculation."""
import argparse
from dataclasses import replace
import json
from pathlib import Path
import sys

import numpy as np
from scipy.integrate import quad
from scipy.special import k1, zeta

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from recursive_horizons.nsc_compact_casimir import (
    compact_response, holonomy_cutoff_energy, matched_local_terms,
)
from recursive_horizons.nsc_covariant_operator import euclidean_operator, smooth_metric
from recursive_horizons.nsc_finite_terms import periodic_derivative
from check_nsc_compact_boundary_action import authenticated_record
from check_nsc_vacuum_charge_matching import hashes

OUTPUT = ROOT / "results/development/compact-casimir.json"
SOURCES = ("scripts/check_nsc_compact_casimir.py",
           "src/recursive_horizons/nsc_compact_casimir.py",
           "docs/nsc-compact-casimir.md",
           "scripts/check_nsc_compact_boundary_action.py",
           "scripts/check_nsc_vacuum_charge_matching.py",
           "src/recursive_horizons/nsc_covariant_operator.py",
           "src/recursive_horizons/nsc_finite_terms.py")
INPUTS = ("results/development/compact-boundary-action.json",
          "results/nsc-9-covariant-source.json")


def native(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {key: native(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [native(item) for item in value]
    return value


def compare(expected, actual, path="$"):
    if isinstance(expected, dict):
        assert isinstance(actual, dict) and expected.keys() == actual.keys(), path
        for key in expected:
            compare(expected[key], actual[key], f"{path}/{key}")
    elif isinstance(expected, list):
        assert isinstance(actual, list) and len(expected) == len(actual), path
        for index, (a, b) in enumerate(zip(expected, actual)):
            compare(a, b, f"{path}/{index}")
    elif isinstance(expected, float):
        assert not isinstance(actual, bool) and isinstance(actual, (int, float)), path
        assert np.isclose(expected, actual, atol=3e-9, rtol=3e-8), (path, expected, actual)
    else:
        assert type(expected) is type(actual) and expected == actual, (path, expected, actual)


def summary(metric, result):
    i = metric.points // 2
    ward = (periodic_derivative(result["p_radial_4D"], metric.length)
            + periodic_derivative(metric.lapse, metric.length) / metric.lapse * result["radial_null_4D"]
            + 2 * periodic_derivative(metric.sphere_radius, metric.length) / metric.sphere_radius
            * (result["p_radial_4D"] - result["p_sphere_4D"]))
    scaling = (np.sum(result["metric_gradients"] * np.array([
        metric.lapse, metric.radial_scale, metric.sphere_radius]))
        + 2 * result["interval_derivative"])
    return {"points": metric.points, "energy": result["energy"],
            "interval_derivative": result["interval_derivative"],
            "neck": {name: float(result[name][i]) for name in
                     ("rho_4D", "p_radial_4D", "p_sphere_4D", "radial_null_4D")},
            "maximum_local_radial_conservation_residual": float(max(abs(ward))),
            "rigid_metric_interval_scaling_residual": float(scaling),
            "lapse_homogeneity_residual": result["lapse_homogeneity_residual"],
            "finite_angular_frequency_energy_tail_bound": result["finite_angular_frequency_energy_tail_bound"]}


def bessel_energy(metric, interval=2., angular_max=12, windings=24):
    assert np.all(metric.lapse == 1)
    total = 0.
    j = np.arange(1, windings + 1, dtype=float)
    for kappa in range(1, angular_max + 1):
        values = np.abs(np.linalg.eigvalsh(euclidean_operator(metric, 0., kappa)))
        total += 2 * 2 * kappa / np.pi * np.sum(
            values[:, None] * k1(2 * interval * values[:, None] * j[None, :]) / j[None, :])
    return float(total)


def calculate():
    for path in INPUTS:
        authenticated_record(path)
    ell, copies = 2., 2
    rank = 4 * copies
    flat = float(3 * rank * zeta(5) / (128 * np.pi ** 2 * ell ** 4))
    direct = quad(lambda p: -rank / (16 * np.pi ** 2) * p ** 3
                  * np.log(-np.expm1(-2 * ell * p)), 0, np.inf,
                  epsabs=1e-13, epsrel=1e-13)[0]
    twisted = quad(lambda p: -rank / (16 * np.pi ** 2) * p ** 3
                   * np.log1p(np.exp(-2 * ell * p)), 0, np.inf,
                   epsabs=1e-13, epsrel=1e-13)[0]
    assert abs(direct - flat) < 2e-14 and abs(twisted / flat + 15 / 16) < 2e-12
    metrics = {n: smooth_metric(n, general=True) for n in (16, 24, 32, 40)}
    responses = {n: compact_response(m, angular_max=12, frequency_points=64)
                 for n, m in metrics.items()}
    radial = [summary(metrics[n], responses[n]) for n in metrics]
    assert abs(radial[-1]["energy"] - radial[-2]["energy"]) < 1e-10
    assert abs(radial[-1]["neck"]["radial_null_4D"] - radial[-2]["neck"]["radial_null_4D"]) < 1e-10
    assert radial[-1]["maximum_local_radial_conservation_residual"] < 1e-10
    m, reference = metrics[32], responses[32]
    angular = []
    for count in (4, 8, 12, 16):
        r = reference if count == 12 else compact_response(m, angular_max=count, frequency_points=64)
        angular.append({"angular_max": count, **summary(m, r)})
    assert abs(angular[-1]["energy"] - angular[-2]["energy"]) < 1e-10
    frequency = []
    for count in (32, 48, 64):
        r = reference if count == 64 else compact_response(m, angular_max=12, frequency_points=count)
        frequency.append({"frequency_points": count, **summary(m, r)})
    assert abs(frequency[-1]["energy"] - frequency[-2]["energy"]) < 1e-10
    extended = compact_response(m, angular_max=12, frequency_points=80, extent_factor=1.25)
    assert abs(extended["energy"] - reference["energy"]) < 1e-10
    ultra = smooth_metric(32)
    ur = compact_response(ultra, angular_max=12, frequency_points=64)
    bessel = [{"windings": n, "energy": bessel_energy(ultra, windings=n)} for n in (16, 24)]
    assert abs(bessel[-1]["energy"] - ur["energy"]) < 1e-10
    local = matched_local_terms(ultra)
    general_local = matched_local_terms(m)

    # New profile derivatives, using the unchanged covariant operator.
    variation = []
    vm, vr = metrics[24], responses[24]
    angle = 2 * np.pi * vm.x / vm.length
    for field_index, field in enumerate(("lapse", "radial_scale", "sphere_radius")):
        direction = .2 * np.cos(angle) + .07 * np.sin(2 * angle + field_index)
        exact = float(np.dot(vr["metric_gradients"][field_index], getattr(vm, field) * direction))
        values = []
        for step in (2e-4, 1e-4):
            energies = [compact_response(replace(vm, **{field: getattr(vm, field) * np.exp(sign * step * direction)}),
                                        angular_max=12, frequency_points=64, gradients=False)["energy"]
                        for sign in (1, -1)]
            fd = (energies[0] - energies[1]) / (2 * step)
            values.append({"step": step, "finite_difference": fd, "error": abs(fd - exact)})
        assert values[-1]["error"] < 2e-9
        variation.append({"field": field, "matrix_derivative": exact, "finite_differences": values})
    step = 1e-5
    ep = compact_response(vm, interval=ell * np.exp(step), angular_max=12, frequency_points=64, gradients=False)["energy"]
    em = compact_response(vm, interval=ell * np.exp(-step), angular_max=12, frequency_points=64, gradients=False)["energy"]
    length_fd = (ep - em) / (2 * step)
    assert abs(length_fd - ell * vr["interval_derivative"]) < 2e-9

    domains = []
    for length, points, count in ((4., 32, 12), (6., 48, 16), (8., 64, 24), (12., 96, 32)):
        dm = smooth_metric(points, radius_parameter=length / 2)
        dr = ur if length == 4 else compact_response(dm, angular_max=count, frequency_points=64)
        domains.append({"axial_length": length, "points": points, "angular_max": count,
                        "energy": dr["energy"], "neck_null": float(dr["radial_null_4D"][points // 2])})
    assert domains[0]["neck_null"] < 0 < domains[-1]["neck_null"]
    pm = smooth_metric(33, eta=0., general=True)
    pr = compact_response(pm, angular_max=12, frequency_points=64)

    # Full five-dimensional free determinant, with its finite local ambiguity
    # cancelled between flat connections. This is not Gamma_int alone.
    hm = smooth_metric(64)
    holonomy = []
    for cutoff, angular_max, compact_max in ((2., 28, 12), (3., 40, 16), (4., 52, 22)):
        baseline = holonomy_cutoff_energy(hm, 0., cutoff=cutoff, angular_max=angular_max,
                                          compact_max=compact_max)["energy"]
        rows = []
        for phase in (0., .125, .25, .375, .5, .75, 1.):
            h = holonomy_cutoff_energy(hm, phase, cutoff=cutoff, angular_max=angular_max,
                                      compact_max=compact_max, derivatives=phase == .5)
            rows.append({"phase": phase, "difference_from_periodic": h.pop("energy") - baseline, **h})
        assert all(rows[i+1]["difference_from_periodic"] < rows[i]["difference_from_periodic"] for i in range(4))
        assert abs(rows[-1]["difference_from_periodic"]) < 2e-9
        assert abs(rows[4]["phase_derivative"]) < 2e-9 and rows[4]["phase_second_derivative"] > 0
        holonomy.append({"cutoff": cutoff, "points": 64, "angular_max": angular_max,
                         "compact_max": compact_max, "phases": rows})
    assert abs(holonomy[-1]["phases"][4]["difference_from_periodic"]
               - holonomy[-2]["phases"][4]["difference_from_periodic"]) < 1e-8
    # Independent changes at fixed cutoff: radial, compact and angular.
    refined = []
    for points, angular_max, compact_max in ((80, 40, 16), (64, 48, 16), (64, 40, 24)):
        mm = smooth_metric(points)
        hs = [holonomy_cutoff_energy(mm, p, cutoff=3., angular_max=angular_max, compact_max=compact_max)
              for p in (0., .5)]
        difference = hs[1]["energy"] - hs[0]["energy"]
        assert abs(difference - holonomy[1]["phases"][4]["difference_from_periodic"]) < 1e-8
        refined.append({"points": points, "angular_max": angular_max, "compact_max": compact_max,
                        "AP_minus_P_energy": difference})
    h = .003
    center = holonomy_cutoff_energy(hm, .5, cutoff=3., angular_max=40, compact_max=16)["energy"]
    neighbors = [holonomy_cutoff_energy(hm, .5 + sign * h, cutoff=3.,
                                        angular_max=40, compact_max=16)["energy"] for sign in (-1, 1)]
    phase_fd = (neighbors[0] - 2 * center + neighbors[1]) / h ** 2
    phase_hessian = holonomy[1]["phases"][4]["phase_second_derivative"]
    assert abs(phase_fd - phase_hessian) < 2e-4
    def remainder(result, terms):
        return {key: result[key] - terms["source"][key] for key in terms["source"]}
    return native({
        "schema": "NSC-COMPACT-CASIMIR-v1",
        "status": "finite curved interaction source and conditional one-loop axial holonomy saddle",
        "source_hashes": hashes(SOURCES), "input_hashes": hashes(INPUTS),
        "domain": {"transverse_interval": ell, "Dirac_copies": copies, "transverse_twist": "untwisted",
                   "zero_level_rank": 4, "massive_level_rank": 8, "KK_masses": "n*pi/ell",
                   "axial_length_reference": 4., "axial_phase_reference": .5, "neck_radius": 1.,
                   "curved_operator": "existing CovariantStaticMetric and paired-angular D_omega",
                   "state": "zero-temperature free vacuum, zero chemical potential; continuous Euclidean time frequency",
                   "scope_of_source": "effective 4D metric variation integrated over transverse interval",
                   "warp_cocycle_and_local_completion_included": False},
        "functional": {"interaction": "-1/2 Tr_two_copies log(1-exp(-2*ell*sqrt(D4^2)))",
                       "excluded_bulk_piece": "formally -ell/2 Tr_two_copies abs(D4), with its ultraviolet prescription; generally nonlocal in curved g4",
                       "derivative": "-ell*sign(lambda)/(exp(2*ell*abs(lambda))-1)",
                       "interval_force_sign": "outward for the untwisted nonzero spectrum",
                       "flat_volume_coefficient": flat, "independent_flat_integral": direct,
                       "twisted_flat_integral": twisted, "twisted_to_untwisted": twisted / flat,
                       "Q2": "3*zeta(5)/(8*ell^4)", "Q1": "zeta(3)/(4*ell^2)",
                       "higher_local_expansion_IR_singular": True},
        "ultrastatic_reference": {"summary": summary(ultra, ur), "full_source": ur,
                                  "matched_local_terms": local, "remaining_interaction": remainder(ur, local),
                                  "Bessel_frequency_integral": bessel},
        "general_metric_reference": {"summary": summary(m, reference), "full_source": reference,
                                     "matched_local_terms": general_local,
                                     "remaining_interaction": remainder(reference, general_local)},
        "radial_refinement": radial, "angular_refinement": angular,
        "frequency_refinement": frequency, "extended_frequency": summary(m, extended),
        "metric_variations": variation,
        "interval_variation": {"step_in_log_interval": step, "matrix_derivative": ell * vr["interval_derivative"],
                               "finite_difference": length_fd, "error": abs(length_fd - ell * vr["interval_derivative"])},
        "domain_family": {"observed": domains, "infinite_domain_bound_claimed": False,
                          "larger_cell_changes_global_geometry": True},
        "periodic_axial_control": {"points": 33, "energy": pr["energy"],
                                  "neck_null": float(pr["radial_null_4D"][pm.points // 2]),
                                  "state_selected_from_this_comparison": False},
        "holonomy": {"meaning": "effective spin plus flat axial U(1) phase modulo one",
                     "potential": "full free KK fermion modulus; local terms cancel in differences",
                     "hypotheses": "ultrastatic closed axial loop, positive r and q, free charged Dirac tower, unfixed holonomy",
                     "operator_partners": "-d_s^2+W^2+/-d_s W, W=kappa/r>0",
                     "positive_determinant": "Delta(omega^2+m_n^2)-2*cos(2*pi*alpha), Delta>2",
                     "minimum": "alpha=1/2 modulo one; conditional one-loop real potential",
                     "proof": "each -c log(Delta-2 cos(2 pi alpha)), c>0, decreases on (0,1/2)",
                     "cutoff_controls": holonomy, "separate_resolution_controls": refined,
                     "second_variation_control": {"phase_step": h, "matrix_hessian": phase_hessian,
                                                  "independent_energy_hessian": phase_fd},
                     "complete_state_or_recursive_return_path_derived": False},
        "scope": {"local_vacuum_coefficient_fixed": False, "determinant_phase_fixed": False,
                  "decompactified_bulk_metric_source_included": False,
                  "transmitting_interface_stress_fixed": False, "self_sourced_geometry_solved": False,
                  "radial_null_sign_universal": False, "cosmological_Q_computed": False,
                  "new_universal_Casimir_or_Wilson_line_mechanism_claimed": False},
        "comparison": {"fields": "all", "exact": "structure, types, strings and source/input hashes",
                       "float_atol": 3e-9, "float_rtol": 3e-8, "exceptions": []},
    })


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--check", action="store_true")
    group.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.output and args.output.exists():
        raise FileExistsError("refusing to overwrite recorded evidence")
    expected = json.loads(OUTPUT.read_text()) if args.check else None
    if args.check:
        compare(expected["source_hashes"], hashes(SOURCES), "$/source_hashes")
        compare(expected["input_hashes"], hashes(INPUTS), "$/input_hashes")
    result = calculate()
    if args.check:
        compare(expected, result)
        print("compact Casimir source: every field and dependency reproduced; prior generators were not run")
    elif args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("x") as stream:
            json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write("\n")
        print(args.output)
    else:
        print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
