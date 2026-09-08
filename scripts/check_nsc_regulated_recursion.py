#!/usr/bin/env python3
"""Reproduce finite regulated variations, frequency reduction and recursion."""
import argparse
from dataclasses import replace
import json
import hashlib
from pathlib import Path
import sys

import numpy as np
from scipy.integrate import quad
from scipy.linalg import expm
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from recursive_horizons.nsc_regulated import (  # noqa: E402
    ClosureSolution, OperatorConventions, RegulatedOperator,
    direct_chain, pg_principal_witness, recursive_response,
)
from check_nsc_scale_closure import compare  # noqa: E402

OUTPUT = ROOT / "results/nsc-3-regulated-recursion.json"


def calculate():
    conventions = OperatorConventions(cutoff=1.3, normalization=.7, throat_length=2., child_ratio=1.6)
    d = np.array([[.8, .3+.1j], [.3-.1j, -1.1]])
    op = RegulatedOperator(d, conventions)
    sigma = np.array([[1., .2j], [-.2j, -.4]])
    first = -(sigma @ d + d @ sigma) / 2
    second = (sigma @ sigma @ d + 2 * sigma @ d @ sigma + d @ sigma @ sigma) / 4
    exact_first, exact_second = op.variation(first, second)
    differences = []
    for step in (.002, .001, .0005):
        minus = op.scale_orbit(sigma, -step).action()
        plus = op.scale_orbit(sigma, step).action()
        differences.append({"step": step, "first": (plus-minus)/(2*step),
                            "second": (plus-2*op.action()+minus)/step**2})
    assert abs(differences[-1]["first"]-exact_first) < 1e-6
    assert abs(differences[-1]["second"]-exact_second) < 1e-6

    defect = op.scale_defect(sigma, .3)
    defect_integral, quadrature_error = quad(lambda t: op.scale_defect_derivative(sigma, t), 0, .3,
                                             epsabs=1e-12, epsrel=1e-12)
    assert abs(defect-defect_integral) < 1e-12
    p, phi, cutoff = .7, .4, 1.2
    relative = RegulatedOperator([[p, phi], [phi, -p]], OperatorConventions(cutoff=cutoff))
    generator = np.diag([1., -1.])
    rfirst = -(generator@relative.matrix+relative.matrix@generator)/2
    rsecond = relative.matrix.copy()
    rsecond[0, 1] = rsecond[1, 0] = 0
    relative_hessian = relative.variation(rfirst, rsecond)[1]
    expected_relative = 4*p*p/cutoff**2*np.exp(-(p*p+phi*phi)/cutoff**2)
    assert abs(relative_hessian-expected_relative) < 1e-12

    a = 1/conventions.cutoff**2
    spatial_squared = op.values**2
    integral, integration_error = quad(
        lambda t: float(np.exp(-t*spatial_squared).sum())/(4*np.sqrt(np.pi)*t**1.5),
        a, np.inf, epsabs=1e-12, epsrel=1e-12)
    closed = op.ultrastatic_action_per_time()
    assert abs(integral-closed) < 1e-12
    factor = 2.3
    rescaled = RegulatedOperator(factor*d, replace(conventions, cutoff=factor*conventions.cutoff,
                                  normalization=factor*conventions.normalization,
                                  throat_length=conventions.throat_length/factor))
    assert abs(rescaled.action()-op.action()) < 1e-12
    assert abs(rescaled.conventions.zeta-conventions.zeta) < 1e-12
    cutoff_step = 1e-5
    cutoff_fd = (RegulatedOperator(d, replace(conventions, cutoff=conventions.cutoff*np.exp(cutoff_step))).action()
                 - RegulatedOperator(d, replace(conventions, cutoff=conventions.cutoff*np.exp(-cutoff_step))).action())/(2*cutoff_step)
    assert abs(cutoff_fd-op.cutoff_derivative()) < 1e-9
    common_gradient, common_hessian = op.variation(-d, d)
    common_heat = float(np.exp(-(op.values / conventions.cutoff)**2).sum())
    assert abs(common_gradient-common_heat) < 1e-12 and common_gradient > 0

    local = np.array([[.8, .1], [.1, 1.2]])
    link = np.array([[.12, .03j], [.04, .1]])
    energy = .7+.2j
    chain_controls = []
    for omega in (1., 1.3, 2.):
        previous = None
        for depth in (2, 4, 8, 16, 32):
            reduced = recursive_response(local, link, energy, omega, depth)
            direct = direct_chain(local, link, energy, omega, depth)
            error = float(np.linalg.norm(reduced-direct))
            imaginary = (reduced-reduced.conj().T)/(2j)
            imaginary_min = float(np.linalg.eigvalsh(imaginary).min())
            assert error < 1e-11 and imaginary_min > 0
            row = {"omega": omega, "depth": depth, "direct_reduced_error": error,
                   "response_real": reduced.real.tolist(), "response_imag": reduced.imag.tolist(),
                   "minimum_imaginary_eigenvalue": imaginary_min,
                   "previous_depth_change": None if previous is None else float(np.linalg.norm(reduced-previous))}
            chain_controls.append(row)
            previous = reduced
    # Exact symbol proof at the throat, where beta²=3*pi/2>1.
    beta, kr, kt = sp.symbols('beta kr kt', positive=True)
    symbol = sp.Matrix([[-beta*kr+kr, kt], [kt, -beta*kr-kr]])
    determinant = sp.factor(symbol.det())
    null_symbol = sp.simplify(determinant.subs({kr: 1, kt: sp.sqrt(beta**2-1)}))
    assert null_symbol == 0
    closure = ClosureSolution("finite_matrix_recursive_response", max(x['direct_reduced_error'] for x in chain_controls),
                              ("covariant field measure and relative compensator", "stationary curved geometry and link",
                               "finite internal Dirac spectrum", "state-dependent retarded metric response"))
    return {
        "schema": "nsc-regulated-recursion-v1", "artifact_id": "NSC-3-REGULATED-RECURSION",
        "source_hashes": {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in (
            "src/recursive_horizons/nsc_regulated.py", "scripts/check_nsc_regulated_recursion.py")},
        "classification": "finite_regulated_variations_and_normalized_chain_closure_with_spatial_spacetime_obstruction",
        "normalization": {"action": "1/2 Tr E1(D²/Lambda²)+N[log(M/Lambda)+gamma_E/2]",
                          "measure_scope": "explicit finite orthonormal matrix basis; not continuum gravitational measure",
                          "cutoff": conventions.cutoff, "M": conventions.normalization,
                          "L_star": conventions.throat_length, "zeta": conventions.zeta,
                          "Omega": conventions.child_ratio, "zeta_equals_Omega_squared_assumed": False},
        "finite_variations": {"gradient": exact_first, "hessian": exact_second,
                              "centered_differences": differences, "relative_hessian": relative_hessian,
                              "relative_hessian_exact_formula": float(expected_relative),
                              "cutoff_gradient": op.cutoff_derivative(), "cutoff_finite_difference": cutoff_fd,
                              "unwrapped_negative_logdet_phase": op.negative_logdet_phase()},
        "finite_scale_cocycle": {"definition": "C(t)=Gamma_reg(D_t)-Gamma_reg(D_0)-t Tr Sigma",
                                 "derivative": "Tr Sigma[exp(-D_t²/Lambda²)-I]",
                                 "endpoint": defect, "independent_integral": defect_integral,
                                 "quadrature_error": quadrature_error,
                                 "physical_relative_compensator_derived": False},
        "frequency_reduction": {"hypothesis": "ultrastatic Euclidean product D_E²=-partial_tau²+H²",
                                "heat_trace_per_time": "Tr_sp exp(-t H²)/sqrt(4*pi*t)",
                                "action_per_time_formula": "sum[Lambda*exp(-e²/Lambda²)/(2sqrt(pi))-|e|*erfc(|e|/Lambda)/2]",
                                "closed_value": closed, "independent_integral": integral,
                                "quadrature_error": integration_error,
                                "valid_for_shifted_black_universe_without_further_derivation": False},
        "unit_change": {"factor": factor, "action_residual": rescaled.action()-op.action(),
                        "zeta_residual": rescaled.conventions.zeta-conventions.zeta},
        "finite_stationarity_obstruction": {
            "family": "D(t)=exp(-t)D with fixed cutoff, normalization and finite domain",
            "derivative_identity": "dGamma/dt=Tr exp(-D(t)^2/Lambda²)>0",
            "gradient": common_gradient, "independent_heat_trace": common_heat,
            "hessian": common_hessian,
            "isolated_stationary_common_scale_in_this_family": False,
            "scope": "finite prescribed determinant; excludes a separately derived continuum field measure and constrained metric variation",
        },
        "normalized_recursion": {"convention": "H_n=Omega^n H; B_n=Omega^n b; fixed dimensional E",
                                 "equation": "Gamma_n(x)=x I-H-(1/Omega)b Gamma_(n+1)(x/Omega)^-1 b_dagger",
                                 "local_matrix": local.tolist(), "link_real": link.real.tolist(),
                                 "link_imag": link.imag.tolist(), "energy": [energy.real, energy.imag],
                                 "controls": chain_controls, "tail_error_bound_proved": False,
                                 "link_derived_from_throat_in_this_control": False},
        "coordinate_time_symbol": {"determinant": str(determinant), "exact_null_residual": str(null_symbol),
                                   "numerical_witness": pg_principal_witness(),
                                   "consequence": "spatial H with lapse/shift is not elliptic in trapped region; elliptic spatial heat calculus cannot be imported as the full spacetime determinant"},
        "unresolved_equations": list(closure.unresolved_equations),
        "nonclaims": {"full_covariant_action_derived": False,
                      "physical_stationarity_solved": closure.physical_stationarity_solved,
                      "physical_gap_or_cosmological_scale_predicted": False,
                      "healthy_metric_or_child_transition_established": False,
                      "new_to_world_priority_established": False},
        "terminal": True,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--check', action='store_true')
    group.add_argument('--output', type=Path)
    args = parser.parse_args()
    result = calculate()
    if args.check:
        compare(json.loads(OUTPUT.read_text()), result)
        print('Regulated variations, frequency reduction and matrix recursion reproduced; all scientific fields compared.')
    elif args.output:
        with args.output.open('x') as stream:
            json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write('\n')
        print(args.output)
    else:
        print(json.dumps(result, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
