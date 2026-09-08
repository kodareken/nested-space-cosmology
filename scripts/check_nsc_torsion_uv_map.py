#!/usr/bin/env python3
"""Match the five-dimensional proper-time bulk torsion coefficients.

Reuses the general-dimensional Bochner/a2 formulas and the universal heat
coefficient. Only two flat quadratic derivative polarizations are evaluated.
No complete curved a4, finite renormalization or physical state is inferred.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from recursive_horizons.nsc_covariant_identities import _matrices

OUTPUT = ROOT / "results/development/torsion-uv-map.json"
SOURCES = (
    "scripts/check_nsc_torsion_uv_map.py",
    "src/recursive_horizons/nsc_covariant_identities.py",
    "docs/nsc-torsion-uv-map.md",
)


def source_hashes():
    return {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in SOURCES}


def compare(expected, actual, path="$"):
    # All scientific outputs here are exact strings/integers, not floats.
    if json.dumps(expected, sort_keys=True) != json.dumps(actual, sort_keys=True):
        raise AssertionError(f"exact fields differ at {path}")


def derivative_polarization(indices):
    gamma4 = _matrices()[0]
    gamma5 = gamma4[0] * gamma4[1] * gamma4[2] * gamma4[3]
    clifford = tuple(sp.I * g for g in (*gamma4, gamma5))
    # T=tau(x0) dx_i wedge dx_j wedge dx_k. Set tau'=1;
    # each independent 3-form component has component norm squared 6*tau².
    connection_derivative = [sp.zeros(4) for _ in range(5)]
    for position, index in enumerate(indices):
        others = [value for value in indices if value != index]
        connection_derivative[index] = (sp.Rational(3, 2) * (-1) ** position
                                        * clifford[others[0]] * clifford[others[1]])
    potential_linear = sp.zeros(4)
    if 0 not in indices:
        potential_linear = -sp.Rational(3, 2) * clifford[0]
        for index in indices:
            potential_linear *= clifford[index]
    curvature_trace = 0
    for i in range(5):
        for j in range(5):
            omega = (connection_derivative[j] if i == 0 else sp.zeros(4))
            omega -= connection_derivative[i] if j == 0 else sp.zeros(4)
            curvature_trace += sp.trace(omega * omega)
    e_part = sp.simplify(sp.trace(potential_linear * potential_linear) / 2)
    omega_part = sp.simplify(curvature_trace / 12)
    return {"indices_zero_based": list(indices), "derivative_axis": 0,
            "half_trace_E_squared": str(e_part),
            "one_twelfth_trace_Omega_squared": str(omega_part),
            "quadratic_a4_coefficient": str(sp.simplify(e_part + omega_part))}


def calculate():
    longitudinal = derivative_polarization((0, 1, 2))
    transverse = derivative_polarization((1, 2, 3))
    if longitudinal["quadratic_a4_coefficient"] != "-3":
        raise AssertionError("longitudinal derivative normalization")
    if transverse["quadratic_a4_coefficient"] != "0":
        raise AssertionError("transverse derivative cancellation")
    cutoff = sp.Symbol("Lambda", positive=True)
    matching = sp.Symbol("nu_match", positive=True)
    a_r, a_t = sp.symbols("c_R c_T", real=True)
    c5 = (4 * sp.pi) ** sp.Rational(5, 2)
    einstein = (cutoff ** 3 - matching ** 3) / (9 * c5)
    stiffness = (cutoff ** 3 - matching ** 3) / c5
    if sp.simplify(stiffness / einstein) != 9:
        raise AssertionError("leading stiffness ratio")
    # Use differential-form norms for this derivative comparison.
    mass_coefficient = 6 * stiffness
    gradient_magnitude = 3 * (cutoff - matching) / c5
    derivative_scale = 2 * (cutoff ** 2 + cutoff * matching + matching ** 2)
    if sp.simplify(mass_coefficient / gradient_magnitude - derivative_scale) != 0:
        raise AssertionError("five-dimensional proper-time derivative scale")
    renormalized_ratio = (stiffness + a_t) / (einstein + a_r)
    mismatch = (a_t - 9 * a_r) / (einstein + a_r)
    if sp.simplify(renormalized_ratio - 9 - mismatch) != 0:
        raise AssertionError("finite matching dependence")
    beta_r = -cutoff * sp.diff(einstein, cutoff)
    beta_t = -cutoff * sp.diff(stiffness, cutoff)
    if sp.simplify(beta_t - 9 * beta_r) != 0:
        raise AssertionError("unfixed finite relative coefficient")
    return {
        "schema": "NSC-TORSION-UV-MAP-v2",
        "status": "development; bulk ultraviolet matching, not a completed quantum source",
        "sources": ["https://arxiv.org/html/1101.1424v3", "https://arxiv.org/html/hep-th/0306138v3"],
        "source_hashes": source_hashes(),
        "domain": "five-dimensional Riemannian bulk; local flat quadratic torsion test; no boundary terms",
        "regulator": "one-half integral from Lambda^-2 to nu_match^-2 of dt/t times Tr exp(-t D^2); 0<nu_match<Lambda",
        "matching_split": "remaining one-loop modulus uses proper times above nu_match^-2; the full determinant is not added again",
        "normalization_scale": "M is distinct from Lambda and nu_match; its full measure matching is not set here",
        "component_norm": "||K||_comp^2=sum_ABC K_ABC^2=6*|K|_form^2; physical torsion=2K",
        "heat_density_without_4pi_factor": {"a0": "4", "a2": "-R/3+3*||K||_comp^2",
                                             "flat_quadratic_integrated_a4": "-3*|delta K|_form^2"},
        "proper_time_weights": {"a0": "(Lambda^5-nu_match^5)/5", "a2": "(Lambda^3-nu_match^3)/3", "a4": "Lambda-nu_match"},
        "longitudinal_polarization": longitudinal,
        "transverse_polarization": transverse,
        "bulk_coefficients": {"common_denominator": "(4*pi)^(5/2)",
                              "volume": "4*(Lambda^5-nu_match^5)/(5*(4*pi)^(5/2))",
                              "Einstein_A5": "(Lambda^3-nu_match^3)/(9*(4*pi)^(5/2))",
                              "K_component_norm_squared": "9*A5",
                              "delta_K_form_norm_squared": "-3*(Lambda-nu_match)/(4*pi)^(5/2)",
                              "leading_stiffness_ratio": 9,
                              "derivative_coefficient_scale_squared": "2*(Lambda^2+Lambda*nu_match+nu_match^2)",
                              "leading_UV_derivative_scale_squared": "2*Lambda^2 as nu_match/Lambda -> 0",
                              "derivative_correction_relative_size": "q_E^2/[2*(Lambda^2+Lambda*nu_match+nu_match^2)]",
                              "first_KK_current_harmonic_unwarped_reference": "pi^2/[2*zeta*(1+u+u^2)], u=nu_match/Lambda",
                              "local_expansion_condition": "q_E, |K|, sqrt(||nabla K||), sqrt(||Riemann||) << nu_match; compact boundary terms separately required",
                              "leading_only_contact_prefactor_after_compact_reduction": "(4*pi)^(5/2)/[64*(Lambda^3-nu_match^3)*I3]"},
        "finite_matching": {"Einstein_coefficient": "A5+c_R", "torsion_coefficient": "9*A5+c_T",
                            "s_T": "(9*A5+c_T)/(A5+c_R)",
                            "s_T_minus_9": "(c_T-9*c_R)/(A5+c_R)",
                            "c_R": None, "c_T": None,
                            "leading_cutoff_flow_leaves_c_T_minus_9_c_R_unchanged": True},
        "scope": {"leading_bulk_ratio_matched": True,
                  "full_curved_5D_a4_computed": False,
                  "complete_finite_cutoff_kernel_computed": False,
                  "finite_coefficients_selected": False,
                  "Lorentzian_torsion_mass_or_health_proved": False,
                  "compact_boundary_action_matched": False,
                  "state_or_self_sourced_geometry_solved": False},
        "comparison": "all fields; exact symbolic strings, integers, types and authenticated source hashes",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--check", action="store_true")
    group.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.output and args.output.exists():
        raise FileExistsError("refusing to overwrite recorded evidence")
    expected = None
    if args.check:
        expected = json.loads(OUTPUT.read_text())
        compare(expected["source_hashes"], source_hashes(), "$/source_hashes")
    result = calculate()
    if args.check:
        compare(expected, result)
        print("five-dimensional torsion UV map: all exact fields reproduced")
    elif args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("x") as handle:
            json.dump(result, handle, indent=2, sort_keys=True, allow_nan=False)
            handle.write("\n")
        print(args.output)
    else:
        print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
