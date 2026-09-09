#!/usr/bin/env python3
"""Match the finite-cutoff compact Dirac source to the physical Gaussian warp."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
from scipy.integrate import quad
from scipy.special import exp1
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from recursive_horizons.nsc_warped_source import (
    WarpedCompactWeight, omitted_tail_bounds, product_potential, source_potential,
)
from check_nsc_compact_boundary_action import authenticated_record
from check_nsc_compact_casimir import compare
from check_nsc_vacuum_charge_matching import hashes

OUTPUT = ROOT / "results/development/warped-source.json"
SOURCES = ("scripts/check_nsc_warped_source.py",
           "src/recursive_horizons/nsc_warped_source.py", "docs/nsc-warped-source.md",
           "scripts/check_nsc_compact_boundary_action.py",
           "scripts/check_nsc_compact_casimir.py",
           "scripts/check_nsc_vacuum_charge_matching.py",
           "docs/nsc-compact-mass-map.md")
INPUTS = ("results/development/horizon-source.json",
          "results/development/compact-boundary-action.json",
          "results/development/vacuum-charge-matching.json")
PHYSICAL_FIELDS = ("potential", "d_radius", "d_interval", "d_path", "d_cutoff")


def calculate():
    for path in INPUTS:
        authenticated_record(path)
    exact = {}

    def zero(name, expression):
        assert sp.simplify(expression) == 0, name
        exact[name] = "0"

    A4, C4, r, q, nu = sp.symbols("A4 C4 r q nu", positive=True)
    G = 1/(16*sp.pi*A4)
    zero("spherical_Einstein_term_counted_once", 2*G*(-8*sp.pi*A4)+1)
    zero("magnetic_source_normalization", 2*G*(2*sp.pi*C4*q*q/r**2)-sp.pi*G*q*q*4*C4/r**2)
    V = sp.Symbol("V")
    low = q*nu**2/(4*sp.pi)
    zero("matching_scale_cancels_between_low_and_rest", sp.diff(low+V-low, nu))

    # A real coefficient matrix still represents COMPLEX Dirac spinors.
    # Independent componentwise norm detects the weighted mixing term.
    f, g, fp, gp, lam = sp.symbols("f g fprime gprime lambda", real=True)
    sigma1, sigma3 = sp.Matrix([[0, 1], [1, 0]]), sp.diag(1, -1)
    Dpsi = lam*sigma1*sp.Matrix([f, sp.I*g])-sp.I*sigma3*sp.Matrix([fp, sp.I*gp])
    norm = sp.expand((Dpsi.conjugate().T*Dpsi)[0])
    zero("complex_spinor_quadratic_form", norm-((lam*g-fp)**2+(lam*f-gp)**2))
    trial = WarpedCompactWeight(modes=12)
    coefficients = np.zeros(25)
    coefficients[0] = coefficients[14] = 1.
    matrix_norm = float(coefficients @ (trial.k0+trial.k1+trial.k2) @ coefficients)
    no_mixing = float(coefficients @ (trial.k0+trial.k2) @ coefficients)

    def direct_norm(x):
        g = np.sqrt(2)*np.sin(2*np.pi*x)
        gp = np.sqrt(2)*np.pi*np.cos(2*np.pi*x)
        return np.exp((18/1015)*(2*x-1)**2)*(g*g+(1-gp)**2)

    direct, direct_error = quad(direct_norm, 0, 1, epsabs=1e-12)
    assert abs(direct-matrix_norm) < 1e-12
    assert abs(no_mixing-direct) > .01
    mixing_control = {"trial": "lambda=1, ell=2, f=1, g=sqrt(2)*sin(2*pi*x), chi=(f,i*g)",
                      "direct_complex_field_norm": direct, "matrix_norm": matrix_norm,
                      "quadrature_error_estimate": direct_error,
                      "norm_if_mixing_discarded": no_mixing,
                      "error_if_mixing_discarded": no_mixing-direct}

    # New compact evaluator recovers the imported tower; no old generator runs.
    flat_weight = WarpedCompactWeight(path=0., modes=12)
    spectral_controls = []
    for y in (.25, 1.3, 5.):
        b2 = (np.pi*np.arange(13)/2)**2
        multiplicity = np.array([1]+[2]*12)
        expected = float(multiplicity @ exp1((y+b2)/4))
        response = flat_weight.response(y)
        assert abs(response["weight"]-expected) < 2e-12
        assert abs(response["d_path"]-response["dual_heat_path_derivative"]) < 2e-12
        spectral_controls.append({"D4_squared": y, "expected_KK_weight": expected,
                                  "computed": response})

    convergence = []
    for modes in (12, 18, 26):
        response = source_potential(WarpedCompactWeight(modes=modes), angular_max=10)
        convergence.append({"compact_modes": modes, "response": response})
    base = convergence[-1]["response"]
    assert convergence[0]["response"]["potential"] < convergence[1]["response"]["potential"] < base["potential"]
    assert abs(base["scaling_residual"]) < 1e-8
    assert abs(base["path_variation_residual"]) < 1e-7
    assert abs(base["potential"]-convergence[-2]["response"]["potential"]) < 1e-6

    changes = {}
    variants = {
        "compact_quadrature_128": (dict(modes=26, points=128), {}),
        "angular_max_12": (dict(modes=26), dict(angular_max=12)),
        "momentum_extent_14": (dict(modes=26), dict(momentum_extent=14.)),
        "quadrature_tolerance_2e-10": (dict(modes=26), dict(quadrature_atol=2e-10)),
    }
    for name, (operator, quadrature) in variants.items():
        response = source_potential(WarpedCompactWeight(**operator), **quadrature)
        difference = {key: response[key]-base[key] for key in PHYSICAL_FIELDS}
        assert max(abs(x) for x in difference.values()) < 1e-7, name
        changes[name] = {"response": response, "difference_from_base": difference}

    # Independent finite variations of the same new potential, at fixed basis.
    controls = {"radius": 1., "interval": 2., "path": .8, "cutoff": 2.}
    point = source_potential(WarpedCompactWeight(modes=18, path=.8))
    finite_differences = {}
    for field in controls:
        step = 2e-4

        def evaluate(delta):
            args = dict(controls)
            args[field] += delta
            radius = args.pop("radius")
            return source_potential(WarpedCompactWeight(modes=18, **args), radius=radius)["potential"]

        derivative = (evaluate(step)-evaluate(-step))/(2*step)
        residual = derivative-point["d_"+field]
        assert abs(residual) < 4e-7, (field, residual)
        finite_differences[field] = {"step": step, "finite_difference": derivative,
                                     "functional_derivative": point["d_"+field], "residual": residual}

    flux_cases = {}
    for flux in (0, 1, 2):
        response = base if flux == 1 else source_potential(WarpedCompactWeight(modes=26), flux=flux)
        upper = product_potential(flux, 1., 2., 2., 14, 26)
        lower = product_potential(flux, 1., 2., 2*np.exp(-18/1015), 14, 26)
        assert lower < response["potential"] < upper
        flux_cases[str(flux)] = {"response": response, "flat_lower_cutoff_comparison": lower,
                                "flat_upper_comparison": upper,
                                "warp_change_from_product": response["potential"]-upper}

    rho = base["potential"]/(4*np.pi)
    sphere_pressure = -base["d_radius"]/(8*np.pi)
    return {
        "schema": "NSC-WARPED-SOURCE-v1",
        "status": "finite-cutoff two-copy Dirac modulus matched to the Gaussian compact warp",
        "source_hashes": hashes(SOURCES), "input_hashes": hashes(INPUTS),
        "domain": {
            "metric": "exp(2*s*sigma(Y))*(g4+dY^2); sigma=-18*(2Y/ell)^2/1015",
            "g4_for_potential": "Euclidean R2 x S2(r)",
            "Y_interval": "[-ell/2,ell/2], reflecting paired chiral domains",
            "quantum_fields": "two free massless complex 5D Dirac copies; unit U(1) charge; no link mass or torsion",
            "functional": "+(1/2) Tr integral_{Lambda^-2}^infinity dt exp(-t Ddagger D)/t",
            "fourD_reduction": "+(1/2) Tr4 h_Lambda(D4^2), with compact generalized eigenproblem",
            "potential_units": "inverse length squared per reference R2 area; hbar=c_light=1",
            "kernel_projector": "no L2 zero projector density on noncompact R2; separate prescription on a compact radial box",
            "parameters": {"cutoff": 2., "interval": 2., "radius": 1., "amplitude": 18/1015,
                           "warp_path": 1., "angular_max": 10, "compact_modes": 26,
                           "compact_points": 96, "momentum_extent": 12., "quadrature_atol": 2e-9},
        },
        "exact_matching_residuals": exact,
        "complex_spinor_mixing_control": mixing_control,
        "flat_spectral_controls": spectral_controls,
        "base_source": base,
        "compact_convergence": convergence,
        "independent_resolution_changes": changes,
        "first_variation_check_point": controls,
        "first_variations": finite_differences,
        "flux_cases": flux_cases,
        "tail_bounds": omitted_tail_bounds(1, 1., 2., 2., 10, 26, 12.),
        "bound_scope": "positive continuum-sector tail comparisons, evaluated in floating point; retained Galerkin eigenvalue errors are assessed separately and not bounded by those tails",
        "homogeneous_variational_source": {
            "rho": rho, "p_parallel": -rho, "p_sphere": sphere_pressure,
            "radial_null": 0., "sphere_null": rho+sphere_pressure,
            "scope": "Euclidean homogeneous projections integrated over compact Y; not the in-in stress on a varying Lorentzian throat",
        },
        "source_partition": {
            "low_canonical_LLL_potential": "abs(q)*nu^2/(4*pi)",
            "rest": "V5_Lambda - V_low_nu, before moving the same induced geometric functional to the left-hand side",
            "Einstein_sphere_potential": "-8*pi*A4; adding this to U while retaining the same induced Einstein term double-counts it",
            "full_U_from_this_potential_alone": False,
        },
        "scope": {
            "classical_KK_masses_changed_by_warp": False,
            "finite_cutoff_determinant_is_warp_invariant": False,
            "new_universal_heat_or_monopole_spectrum_claimed": False,
            "all_localized_a5_gauge_invariants_checked": False,
            "full_quantum_compensator_or_finite_coefficients_matched": False,
            "Lorentzian_state_or_transmission_phase_selected": False,
            "joint_self_sourced_geometry_solved": False,
            "cosmological_Q_or_dark_fraction_predicted": False,
        },
        "comparison": {"fields": "all", "float_atol": 3e-9, "float_rtol": 3e-8,
                       "exact": "keys, lengths, non-float types, strings and source/input hashes", "exceptions": []},
    }


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
        print("warped source: every field reproduced; prior generators were not run")
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
