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
from math import pi, sqrt
from scipy.integrate import quad_vec
from scipy.special import erf, erfc, exp1
from recursive_horizons.nsc_warped_source import integrated_exp1
from recursive_horizons.nsc_warped_source import WarpedCompactWeight
from check_nsc_vacuum_charge_matching import hashes

OUTPUT = ROOT/"results/development/compact-matching.json"
SOURCES = ("scripts/check_nsc_compact_matching.py",
           "docs/nsc-compact-matching.md", "docs/nsc-covariant-source.md",
           "src/recursive_horizons/nsc_finite_terms.py",
           "scripts/check_nsc_vacuum_charge_matching.py")
INPUTS = ("results/development/warped-source.json",
          "results/nsc-9-covariant-source.json", "results/nsc-10-influence.json",
          "results/nsc-10-measure-normalization.json")




def proper_ratio(weight):
    """J=ell_proper/ell for the existing Gaussian warp."""
    a = weight.path*weight.amplitude
    return 1. if a == 0 else sqrt(pi)*erf(sqrt(a))/(2*sqrt(a))

def zero_limit(weight, matching_cutoff, maximum=24):
    """Finite limit of h_Lambda(y)-E1(y/nu²) as y decreases to zero.

    At lambda=0 the proper compact coordinate gives the exact nonzero
    singular values p*pi/ell_proper. The light singular value is
    epsilon0(y)=y/J²+O(y²), though its canonical Dirac mass stays zero.
    This limit is not a prescription to discard actual spacetime zero modes.
    """
    if not np.isfinite(matching_cutoff) or matching_cutoff <= 0:
        raise ValueError("positive matching cutoff required")
    if not isinstance(maximum, int) or maximum < 1:
        raise ValueError("at least one compact level required")
    ratio = proper_ratio(weight)
    length = weight.interval*ratio
    factor = (pi/(length*weight.cutoff))**2
    massive = 2*float(np.sum(exp1(factor*np.arange(1, maximum+1)**2)))
    value = 2*np.log(weight.cutoff*ratio/matching_cutoff)+massive
    # E1(c*p²)<=exp(-c*p²)/(c*N²), p>N. Integral comparison of the
    # decreasing Gaussian sum gives this bound, evaluated in floating point.
    tail = sqrt(pi)*erfc(sqrt(factor)*maximum)/(factor**1.5*maximum**2)
    return {"proper_ratio": ratio, "proper_length": length,
            "light_singular_slope": 1/ratio**2,
            "matched_weight_at_zero": float(value),
            "massive_zero_argument_weight": massive,
            "zero_limit_tail_bound": float(tail)}

def matched_weight(weight, squared, matching_cutoff):
    if not np.isfinite(squared) or squared < 0:
        raise ValueError("nonnegative squared spectral argument required")
    if squared == 0:
        return zero_limit(weight, matching_cutoff)["matched_weight_at_zero"]
    if not np.isfinite(matching_cutoff) or matching_cutoff <= 0:
        raise ValueError("positive matching cutoff required")
    return weight.response(squared, derivatives=False)["weight"]-float(exp1(squared/matching_cutoff**2))

def spectral_moments(weight, extent=12., tolerance=1e-9):
    """Q1[h]=int h(y)dy and Q2[h]=int y*h(y)dy, with first variations.

    Compact eigenvalues come from the existing owner; no angular or radial
    discretization is introduced. Omitted momentum bounds exclude retained
    Galerkin errors and are evaluated with ordinary floating-point arithmetic.
    """
    if extent < weight.cutoff or tolerance <= 0:
        raise ValueError("resolved positive integration extent/tolerance required")
    names = ("weight", "d_interval", "d_path", "d_cutoff")

    def integrand(k):
        response = weight.response(k*k)
        v = np.array([response[name] for name in names])
        return np.concatenate((2*k*v, 2*k**3*v))

    values, error = quad_vec(integrand, 0, extent, epsabs=tolerance,
                            epsrel=tolerance, norm="max")
    q1, q2 = dict(zip(names, map(float, values[:4]))), dict(zip(names, map(float, values[4:])))
    heat_bound = 1+weight.interval*weight.cutoff/sqrt(pi)
    tail = weight.cutoff**2*heat_bound*np.exp(-(extent/weight.cutoff)**2)
    return {"Q1": q1, "Q2": q2, "quadrature_error_estimate": float(error),
            "Q1_momentum_tail_bound": float(tail),
            "Q2_momentum_tail_bound": float(tail*(extent**2+weight.cutoff**2)),
            "Q1_scaling_residual": weight.interval*q1["d_interval"]-weight.cutoff*q1["d_cutoff"]+2*q1["weight"],
            "Q2_scaling_residual": weight.interval*q2["d_interval"]-weight.cutoff*q2["d_cutoff"]+4*q2["weight"]}

def local_coefficients(weight, moments, matching_cutoff):
    """Known 4D Dirac heat coefficients weighted by the calculated complement.

    Euclidean density: V-A*R_E+C_F*F²+C_W*C_E²+C_Euler*E4+C_box*box_E R_E.
    No R_E² term arises in this free Dirac a4 combination. An independent
    finite R² coefficient in the completed theory is not thereby fixed.
    """
    zero = zero_limit(weight, matching_cutoff)
    q1 = moments["Q1"]["weight"]-matching_cutoff**2
    q2 = moments["Q2"]["weight"]-matching_cutoff**4/2
    q0 = zero["matched_weight_at_zero"]
    norm = (4*pi)**2
    return {**zero, "matching_cutoff": matching_cutoff,
            "Q1_complement": q1, "Q2_complement": q2,
            "V_Dirac": 2*q2/norm, "A_Dirac": q1/(6*norm),
            "C_gauge_Dirac": q0/(3*norm),
            "C_Weyl_Dirac": -q0/(40*norm),
            "C_Euler_Dirac": 11*q0/(720*norm),
            "C_boxR_Dirac": -q0/(60*norm), "C_R2_Dirac": 0.,
            "d_nu_V": -4*matching_cutoff**3/norm,
            "d_nu_A": -matching_cutoff/(3*norm),
            "d_nu_C_gauge": -2/(3*norm*matching_cutoff)}

def flat_moments(weight):
    """Independent analytic momentum integral of the imported unwarped tower."""
    b2 = (np.arange(weight.modes+1)*pi/weight.interval)**2
    copies = np.array([1]+[2]*weight.modes)
    first = integrated_exp1(b2, weight.cutoff)
    second = (weight.cutoff**4*np.exp(-b2/weight.cutoff**2)-b2*first)/2
    return float(copies @ first), float(copies @ second)

def static_basis_coefficients(coefficients, basis_mass=1.):
    """Map to the existing finite_response energy convention, without rerunning it.

    Its basis is (M4,M2 R_L,C2,R_L2,E4,box_L R_L), with R_E=-R_L.
    Static squared invariants and box_E R_E=box_L R_L agree. This reference
    basis mass is neither the matching cutoff nor a determinant normalization.
    """
    if not np.isfinite(basis_mass) or basis_mass <= 0:
        raise ValueError("positive basis reference mass required")
    return [coefficients["V_Dirac"]/basis_mass**4,
            coefficients["A_Dirac"]/basis_mass**2,
            coefficients["C_Weyl_Dirac"], coefficients["C_R2_Dirac"],
            coefficients["C_Euler_Dirac"], coefficients["C_boxR_Dirac"]]

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

def authenticated_record(path, seen=None):
    seen = set() if seen is None else seen
    if path in seen:
        return json.loads((ROOT / path).read_text())
    seen.add(path)
    value = json.loads((ROOT / path).read_text())
    for field in ("source_hashes", "input_hashes"):
        declared = value.get(field, {})
        compare(declared, hashes(tuple(declared)), f"{path}/{field}")
    for dependency in value.get("input_hashes", {}):
        if dependency.endswith(".json"):
            authenticated_record(dependency, seen)
    return value

def compact_boundary_map(rank=8):
    """Mixed chiral D_E†D_E map; reflecting endpoints are not a transmitting wall."""
    if rank <= 0 or rank % 2:
        raise ValueError("positive even spinor rank required")
    return {"rank": rank, "a1": "0",
            "a2": "-rank/[12*(4*pi)^(5/2)] * (integral_M R + 2*integral_boundary K)",
            "a3_geometric": "rank/[128*(4*pi)^2] * integral_boundary (K^2-2*K_ab*K_ab)",
            "positive_operator": "D_E^dagger D_E with complementary adjoint domain",
            "bulk_boundary_ratio": 2,
            "full_source": "OPEN"}


def compact_casimir_protocol(interval):
    """Untwisted compact interaction and conditional axial holonomy saddle."""
    if not np.isfinite(interval) or interval <= 0:
        raise ValueError("positive compact interval required")
    from scipy.special import zeta
    return {"functional": "-1/2 Tr_two_copies log(1-exp(-2*ell*sqrt(D4^2)))",
            "Q2": float(3*zeta(5)/(8*interval**4)),
            "Q1": float(zeta(3)/(4*interval**2)),
            "holonomy_minimum": "alpha=1/2 modulo one",
            "hypotheses": "ultrastatic closed axial loop; positive r,q; free charged Dirac tower; unfixed holonomy",
            "complete_state_or_recursive_return_path_derived": False}


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
