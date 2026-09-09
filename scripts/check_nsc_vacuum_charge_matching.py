#!/usr/bin/env python3
"""Check whether retained spectral coefficients admit the charged-throat seed.

This reuses heat and Reissner--Nordstrom--de Sitter formulas. It is a
compatibility test of the leading positive bulk contribution, not a claim
about all solutions of the completed NSC functional.
"""
import argparse
import hashlib
import json
from pathlib import Path

import sympy as sp
from scipy.integrate import quad

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results/development/vacuum-charge-matching.json"
SOURCES = ("scripts/check_nsc_vacuum_charge_matching.py", "docs/nsc-vacuum-charge-matching.md")
INPUTS = ("results/development/torsion-uv-map.json", "results/development/charged-sector.json",
          "results/development/compact-interaction.json")


def hashes(paths):
    return {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in paths}


def compare(expected, actual, path="$"):
    if isinstance(expected, dict):
        if not isinstance(actual, dict) or expected.keys() != actual.keys():
            raise AssertionError(f"keys differ at {path}")
        for key in expected:
            compare(expected[key], actual[key], f"{path}/{key}")
    elif isinstance(expected, list):
        if not isinstance(actual, list) or len(expected) != len(actual):
            raise AssertionError(f"list differs at {path}")
        for i, (a, b) in enumerate(zip(expected, actual)):
            compare(a, b, f"{path}/{i}")
    elif isinstance(expected, float):
        if isinstance(actual, bool) or not isinstance(actual, (int, float)):
            raise AssertionError(f"numeric type differs at {path}")
        if abs(expected - actual) > 3e-13 + 3e-13 * abs(expected):
            raise AssertionError(f"numeric value differs at {path}")
    elif type(expected) is not type(actual) or expected != actual:
        raise AssertionError(f"exact value differs at {path}")


def calculate():
    uv = json.loads((ROOT / INPUTS[0]).read_text())
    charged = json.loads((ROOT / INPUTS[1]).read_text())
    compact = json.loads((ROOT / INPUTS[2]).read_text())
    assert uv["heat_density_without_4pi_factor"]["a0"] == "4"
    assert charged["gauge_heat_coefficient"]["per_unit_charge_Dirac_a4_F_squared"] == "2/3"
    n, cutoff, f5, f3, f1, i5, i3, i1 = sp.symbols("N Lambda f5 f3 f1 I5 I3 I1", positive=True)
    q = sp.Symbol("q", integer=True, nonzero=True)
    c5 = (4 * sp.pi) ** sp.Rational(5, 2)
    volume = 4 * n * cutoff ** 5 * f5 * i5 / c5
    einstein = n * cutoff ** 3 * f3 * i3 / (3 * c5)
    gauge = 2 * n * cutoff * f1 * i1 / (3 * c5)
    newton = 1 / (16 * sp.pi * einstein)
    gauge_squared = 1 / (4 * gauge)
    vacuum_curvature = volume / (2 * einstein)
    charge_radius_squared = sp.pi * q ** 2 * newton / gauge_squared
    xi = sp.simplify(vacuum_curvature * charge_radius_squared)
    target = 3 * q ** 2 * f5 * f1 * i5 * i1 / (f3 ** 2 * i3 ** 2)
    assert sp.simplify(xi - target) == 0
    mellin_factor = sp.simplify(sp.gamma(sp.Rational(3, 2)) ** 2 /
                              (sp.gamma(sp.Rational(5, 2)) * sp.gamma(sp.Rational(1, 2))))
    assert mellin_factor == sp.Rational(1, 3)
    u = sp.Symbol("u", real=True)
    positive_polynomial = (1 - u) ** 4 * (4 * u ** 2 + 7 * u + 4)
    assert sp.expand(9 * (1-u**5) * (1-u) - 5 * (1-u**3)**2 - positive_polynomial) == 0
    radius, mass_length, rq2, lam = sp.symbols("r GM rQ2 lambda4", positive=True)
    lapse = 1 - 2 * mass_length / radius + rq2 / radius ** 2 - lam * radius ** 2 / 3
    assert sp.simplify(lapse + radius * sp.diff(lapse, radius)
                       - (1 - rq2 / radius ** 2 - lam * radius ** 2)) == 0
    z = sp.Symbol("z", real=True)
    assert sp.expand(sp.Rational(1, 4) - z * (1-z) - (z-sp.Rational(1, 2)) ** 2) == 0

    # Known conformal-curvature formula; only its one-dimensional integral
    # on the already specified Gaussian warp is evaluated here.
    a = sp.Rational(18, 1015)
    assert float(a) == -compact["domain"]["sigma_coefficient"]
    i5_dimensionless = sp.sqrt(sp.pi/(5*a)) * sp.erf(sp.sqrt(5*a))
    b = 3*a
    second_moment = sp.sqrt(sp.pi)*sp.erf(sp.sqrt(b))/(2*b**sp.Rational(3, 2))-sp.exp(-b)/b
    j_times_length = 4*a*a*second_moment
    boundary_times_length = -4*a*sp.exp(-3*a)
    bulk_factor = sp.simplify((12*j_times_length-8*boundary_times_length)/(12*i5_dimensionless))
    ghy_factor = sp.simplify(j_times_length/i5_dimensionless)
    af = float(a)
    import math
    direct_i5 = quad(lambda y: math.exp(-5*af*y*y), -1, 1, epsabs=1e-13, epsrel=1e-13)[0]
    direct_j = quad(lambda y: 4*af*af*y*y*math.exp(-3*af*y*y), -1, 1, epsabs=1e-13, epsrel=1e-13)[0]
    direct_b = -4*af*math.exp(-3*af)
    direct_bulk = (12*direct_j-8*direct_b)/(12*direct_i5)
    direct_ghy = direct_j/direct_i5
    bulk = float(sp.N(bulk_factor, 40))
    ghy = float(sp.N(ghy_factor, 40))
    assert abs(bulk-direct_bulk) < 2e-14 and abs(ghy-direct_ghy) < 2e-14
    maximum_bulk = 5 * bulk / 3
    return {
        "schema": "NSC-VACUUM-CHARGE-MATCHING-v2",
        "status": "leading-coefficient compatibility obstruction; full functional remains open",
        "source_hashes": hashes(SOURCES), "input_hashes": hashes(INPUTS),
        "assumptions": ["N identical unit-charge complex Dirac fields with a common positive scalar spectral profile",
                        "five-dimensional leading interior heat coefficients and a fixed Y-dependent common conformal warp",
                        "four-dimensional metric and gauge field independent of Y, A_Y=0; compact size held fixed in this coefficient projection",
                        "the displayed vacuum, Einstein and gauge channels do not include an independently evaluated scalar-link potential",
                        "standard local Lorentzian Einstein-Maxwell interpretation with positive A and C",
                        "nonzero integer magnetic flux q; rQ is a charge radius, not a particle Compton length",
                        "finite completion, bosonic/ghost loops, boundary action and higher-curvature terms are not silently included"],
        "coefficient_definitions": {"lorentzian_action": "integral sqrt(-g)*(A*R-C*F_mn*F^mn-V)",
                                    "GN": "1/(16*pi*A)", "gauge_coupling_squared": "1/(4*C)",
                                    "lambda4": "V/(2*A)", "rQ_squared": "q^2*C/(4*A)",
                                    "V_volume": "4*N*Lambda^5*f5*I5/(4*pi)^(5/2)",
                                    "A_bulk": "N*Lambda^3*f3*I3/[3*(4*pi)^(5/2)]",
                                    "C_bulk": "2*N*Lambda*f1*I1/[3*(4*pi)^(5/2)]"},
        "matching_invariant": {"Xi_volume": str(xi), "normalization_residual": "0",
                               "definition": "lambda4*rQ_squared=q^2*V*C/(8*A^2)",
                               "general_nonnegative_profile_lower_bound": "Xi_volume >= q^2",
                               "nonnegative_proper_time_weight_lower_bound": "Xi_volume >= 3*q^2",
                               "mellin_gamma_factor": str(mellin_factor),
                               "warp_moment_bound": "I5*I1/I3^2 >= 1",
                               "copy_count_and_overall_scale_cancel": True},
        "proper_time_window": {"domain": "0<u=nu_match/Lambda<1",
                               "Xi_volume": "(27/5)*q^2*(1-u^5)*(1-u)/(1-u^3)^2 * I5*I1/I3^2",
                               "positive_polynomial_identity": "9*(1-u^5)*(1-u)-5*(1-u^3)^2=(1-u)^4*(4*u^2+7*u+4)"},
        "RN_de_Sitter_control": {"lapse": str(lapse),
                                 "double_horizon_identity": "f+r*f_prime=1-rQ2/r^2-lambda4*r^2",
                                 "extremal_relation": "Xi=z*(1-z), z=lambda4*r_horizon^2",
                                 "bound": "Xi <= 1/4",
                                 "square_identity": "1/4-z*(1-z)=(z-1/2)^2"},
        "specified_warp_check": {"sigma": "-18*(Y/L_star)^2/1015",
                                  "A5_definition": "A_bulk/I3, including the same N Dirac copies",
                                  "scalar_curvature_formula": "R5=exp(-2*sigma)*(R4-8*sigma_second-12*sigma_prime^2)",
                                  "bulk_geometric_volume_contribution": "A5*(8*B-12*J), B=[exp(3*sigma)*sigma_prime]_ends, J=integral exp(3*sigma)*sigma_prime^2 dY",
                                  "bulk_ratio": "V_warp_bulk/V_volume=-(f3/f5)*K_bulk/zeta",
                                  "K_bulk": bulk, "adaptive_K_bulk": direct_bulk,
                                  "K_with_standard_GHY": ghy, "adaptive_K_with_standard_GHY": direct_ghy,
                                  "maximum_window_bulk_warp_correction_factor": maximum_bulk,
                                  "corrected_lower_bound_for_zeta_ge_1": 3*(1-maximum_bulk),
                                  "corrected_lower_bound_units": "Xi_retained >= this coefficient times q^2, retaining only the stated bulk terms",
                                  "standard_GHY_is_only_a_comparison": True,
                                  "NSC_boundary_action_derived": False},
        "required_completion": {"positive_vacuum_seed_condition": "V_full*C_full/A_full^2 <= 2/q^2",
                                "if_only_V_changes_from_positive_proper_time_baseline": "V_full/V_volume <= 1/(12*q^2) is necessary, not sufficient",
                                "no_finite_coefficient_assigned": True},
        "scope": {"full_NSC_excluded": False, "all_charged_solutions_excluded": False,
                  "controlled_direct_MMP_matching_from_raw_positive_bulk_supported": False,
                  "boson_ghost_cancellation_computed": False,
                  "complete_quantum_measure_or_vacuum_coefficient_fixed": False,
                  "new_to_world_priority_claim": False},
        "comparison": {"fields": "all", "exact": "identities, keys, types, source/input hashes",
                       "float_atol": 3e-13, "float_rtol": 3e-13, "exceptions": []},
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
        compare(expected["source_hashes"], hashes(SOURCES), "$/source_hashes")
        compare(expected["input_hashes"], hashes(INPUTS), "$/input_hashes")
    result = calculate()
    if args.check:
        compare(expected, result)
        print("vacuum/charge matching: all fields and source dependencies reproduced")
    elif args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("x") as handle:
            json.dump(result, handle, indent=2, sort_keys=True, allow_nan=False)
            handle.write("\n")
        print(args.output)
    else:
        print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
