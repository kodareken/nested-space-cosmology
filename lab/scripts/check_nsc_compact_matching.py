#!/usr/bin/env python3
"""Bind the computed compact determinant to canonical low-field matching."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
from scipy.linalg import eigvalsh
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))
from recursive_horizons.nsc_compact_matching import (
    flat_moments, local_coefficients, matched_weight, proper_ratio,
    spectral_moments, static_basis_coefficients, zero_limit,
)
from recursive_horizons.nsc_warped_source import WarpedCompactWeight
from check_nsc_compact_boundary_action import authenticated_record
from check_nsc_compact_casimir import compare, native
from check_nsc_vacuum_charge_matching import hashes

OUTPUT = ROOT/"results/development/compact-matching.json"
SOURCES = ("scripts/check_nsc_compact_matching.py",
           "src/recursive_horizons/nsc_compact_matching.py",
           "docs/nsc-compact-matching.md", "docs/nsc-covariant-source.md",
           "src/recursive_horizons/nsc_finite_terms.py",
           "scripts/check_nsc_compact_boundary_action.py",
           "scripts/check_nsc_compact_casimir.py",
           "scripts/check_nsc_vacuum_charge_matching.py")
INPUTS = ("results/development/warped-source.json",
          "results/nsc-9-covariant-source.json", "results/nsc-10-influence.json",
          "results/nsc-10-measure-normalization.json")


def calculate():
    for path in INPUTS:
        authenticated_record(path)
    exact = {}

    def zero(name, expression):
        assert sp.simplify(expression) == 0, name
        exact[name] = "0"

    q1, q2, q0, nu = sp.symbols("Q1 Q2 Q0 nu", positive=True)
    norm = (4*sp.pi)**2
    V, A, C = 2*(q2-nu**4/2)/norm, (q1-nu**2)/(6*norm), (q0-2*sp.log(nu))/(3*norm)
    zero("volume_matching_scale_cancellation", sp.diff(V+nu**4/norm, nu))
    zero("Einstein_matching_scale_cancellation", sp.diff(A+nu**2/(6*norm), nu))
    zero("gauge_matching_scale_cancellation", sp.diff(C, nu)+2/(3*norm*nu))
    y = sp.Symbol("y", positive=True)
    # Differential identity for E1(y/nu²), with positive spectrum/consistent IR treatment.
    zero("full_profile_matching_scale_cancellation",
         sp.diff(sp.expint(1, y/nu**2), nu)-2*sp.exp(-y/nu**2)/nu)
    riemann, ricci, scalar, box = sp.symbols("Riemann2 Ricci2 R_E boxR_E")
    weyl, euler = riemann-2*ricci+scalar**2/3, riemann-4*ricci+scalar**2
    zero("imported_Dirac_a4_basis_conversion",
         5*scalar**2-8*ricci-7*riemann-12*box-(-18*weyl+11*euler-12*box))

    flat = WarpedCompactWeight(path=0., modes=26, points=112)
    flat_numeric = spectral_moments(flat)
    f1, f2 = flat_moments(flat)
    assert abs(flat_numeric["Q1"]["weight"]-f1) < 2e-9
    assert abs(flat_numeric["Q2"]["weight"]-f2) < 2e-9

    convergence = []
    for modes in (18, 26, 34):
        weight = WarpedCompactWeight(modes=modes, points=112)
        moments = spectral_moments(weight)
        assert abs(moments["Q1_scaling_residual"]) < 2e-8
        assert abs(moments["Q2_scaling_residual"]) < 2e-8
        convergence.append({"compact_modes": modes, "moments": moments})
    moments = convergence[-1]["moments"]
    weight = WarpedCompactWeight(modes=34, points=112)
    for key in ("Q1", "Q2"):
        assert convergence[0]["moments"][key]["weight"] < convergence[1]["moments"][key]["weight"] < moments[key]["weight"]
        assert abs(moments[key]["weight"]-convergence[1]["moments"][key]["weight"]) < 2e-7
    refined = spectral_moments(WarpedCompactWeight(modes=34, points=144), extent=14., tolerance=1e-10)
    for key in ("Q1", "Q2"):
        assert max(abs(refined[key][f]-moments[key][f]) for f in refined[key]) < 2e-8

    # Proper-coordinate low spectrum; the old field/KK generators remain untouched.
    zero_data = zero_limit(weight, 1.)
    eigenvalues = eigvalsh(weight.k0, weight.mass)
    assert abs(eigenvalues[0]) < 1e-11
    low_exact = [(p*np.pi/zero_data["proper_length"])**2 for p in range(1, 4) for _ in (0, 1)]
    spectral_errors = eigenvalues[1:7]-low_exact
    assert max(abs(spectral_errors)) < 2e-5
    infrared = []
    for y in (1e-2, 1e-3, 1e-4):
        values = eigvalsh(weight.k0+np.sqrt(y)*weight.k1+y*weight.k2, weight.mass)
        point = matched_weight(weight, y, 1.)
        infrared.append({"D4_squared": y, "light_singular_ratio": float(values[0]/y),
                         "matched_weight": point,
                         "difference_from_zero_limit": point-zero_data["matched_weight_at_zero"]})
    assert abs(infrared[-1]["light_singular_ratio"]-zero_data["light_singular_slope"]) < 2e-6
    assert abs(infrared[-1]["difference_from_zero_limit"]) < 2e-4

    coefficients = [local_coefficients(weight, moments, nu) for nu in (.5, 1., 1.4)]
    base = coefficients[1]
    differences = []
    for other in (coefficients[0], coefficients[2]):
        n = other["matching_cutoff"]
        v = other["V_Dirac"]-base["V_Dirac"]+(n**4-1)/(4*np.pi)**2
        a = other["A_Dirac"]-base["A_Dirac"]+(n*n-1)/(6*(4*np.pi)**2)
        c = other["C_gauge_Dirac"]-base["C_gauge_Dirac"]+2*np.log(n)/(3*(4*np.pi)**2)
        assert max(abs(v), abs(a), abs(c)) < 1e-14
        differences.append({"matching_cutoff": n, "total_delta_V": v, "total_delta_A": a, "total_delta_C_gauge": c})

    # Match the positive/negative curvature conventions to the existing owner.
    mass, R_L, C2, R2, E4, box_L = sp.symbols("M R_L C2 R2 E4 box_L")
    v, a, cw, cr, ce, cb = sp.symbols("V A C_W C_R2 C_E C_box")
    coefficients_L = [v/mass**4, a/mass**2, cw, cr, ce, cb]
    energy_L = sum(c*b for c,b in zip(coefficients_L,(mass**4,mass**2*R_L,C2,R2,E4,box_L)))
    zero("static_metric_owner_sign_and_units", energy_L-(v-a*(-R_L)+cw*C2+cr*R2+ce*E4+cb*box_L))
    return native({
        "schema": "NSC-COMPACT-MATCHING-v1",
        "status": "finite Dirac coefficients matched to one retained canonical 4D field",
        "source_hashes": hashes(SOURCES), "input_hashes": hashes(INPUTS),
        "conventions": {
            "full_functional": "Gamma5=(1/2)Tr4 h_Lambda(D4²), two paired 5D Dirac copies",
            "light_functional": "Gamma_light,nu=(1/2)Tr4 E1(D4²/nu²), one canonical 4D Dirac",
            "complement": "Gamma_H,nu=Gamma5-Gamma_light,nu; exact spectral matching, not a sharp energy projection",
            "matching_cutoff_is_normalization_mass": False,
            "Mellin_definitions": "Q1=int H(y)dy, Q2=int y H(y)dy, Q0=H(0)",
            "coefficient_units": "V length^-4; A length^-2; curvature-squared and gauge coefficients dimensionless",
            "parameters": {"Lambda": 2., "ell": 2., "amplitude": 18/1015, "warp_path": 1.,
                           "basis_modes": 34, "basis_quadrature": 112, "momentum_extent": 12., "quadrature_tolerance": 1e-9},
            "validity": "local coefficients through 4D a4; external curvature, derivatives and gauge field must be small for truncation; compact warp treated by the full finite spectral weight",
        },
        "exact_matching_residuals": exact,
        "flat_moment_check": {"analytic_Q1": f1, "analytic_Q2": f2, "integrated": flat_numeric},
        "compact_convergence": convergence, "refined_moments": refined,
        "proper_coordinate_zero_limit": zero_data,
        "zero_tangent_spectrum": {"computed_first_six": eigenvalues[1:7], "exact_first_six": low_exact,
                                  "errors": spectral_errors},
        "infrared_matching": infrared,
        "matched_coefficients": coefficients,
        "matching_scale_cancellation": differences,
        "static_metric_owner": {"basis": ["M4","M2R_L","C2","R_L2","E4","box_L R_L"],
                                "reference_mass": 1., "coefficients": static_basis_coefficients(base),
                                "full_curved_remainder_recomputed": False},
        "state_completion": {
            "known_owner": "nsc_influence normalized light-field CTP functional",
            "required": "same-regulator causal continuation, light and heavy state/correlations, and completion beyond the free determinant",
            "equal_histories_or_state_normalization_fixes_missing_real_terms": False,
            "arbitrary_extra_vacuum_or_Einstein_coefficient_inserted": False,
        },
        "scope": {"new_universal_heat_coefficients_claimed": False,
                  "zero_Dirac_R2_coefficient_fixes_complete_R2_term": False,
                  "complete_G_Newton_or_gauge_coupling_predicted": False,
                  "complete_quantum_compensator_derived": False,
                  "joint_self_sourcing_or_cosmological_Q_solved": False,
                  "matching_cutoff_may_be_fitted_as_a_physical_parameter": False},
        "comparison": {"fields": "all", "float_atol": 3e-9, "float_rtol": 3e-8,
                       "exact": "keys, lengths, strings, non-float types and source/input hashes", "exceptions": []},
    })


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group()
    group.add_argument("--check",action="store_true")
    group.add_argument("--output",type=Path)
    args=parser.parse_args()
    if args.output and args.output.exists():
        raise FileExistsError("refusing to overwrite recorded evidence")
    expected=json.loads(OUTPUT.read_text()) if args.check else None
    if args.check:
        compare(expected["source_hashes"],hashes(SOURCES),"$/source_hashes")
        compare(expected["input_hashes"],hashes(INPUTS),"$/input_hashes")
    result=calculate()
    if args.check:
        compare(expected,result)
        print("compact matching: every field reproduced; prior generators were not run")
    elif args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        with args.output.open("x") as stream:
            json.dump(result,stream,indent=2,sort_keys=True,allow_nan=False)
            stream.write("\n")
        print(args.output)
    else:
        print(json.dumps(result,indent=2,allow_nan=False))


if __name__=="__main__":
    main()
