#!/usr/bin/env python3
"""Check ZETA1 units and shared-operator identities without importing its decisions.

This is a retrospective calculation on the declared factorized ZETA1 family.
It does not select a physical scale or evaluate a complete Lorentzian action.
Existing result files and generator sources are read-only inputs.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from math import pi
from pathlib import Path

import numpy as np
from scipy.linalg import eigh_tridiagonal, svdvals
from scipy.optimize import brentq
from scipy.special import exp1
import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results/nsc-2-zeta1-unit-closure-check.json"
DEFAULT_RESULT = OUTPUT
AUXILIARY_ROOT = ROOT / "provenance" / "nsc-2-zeta1-unit-closure-check"
ANCHOR = "ff2cf2722b966589b98a61accdbb6cee819a58c7"
FOLLOW_UP_COMMIT = "5f38712ca01ddd71e715fd265088925a73369aba"
# Declared so the public graph copies these generators without adding
# regenerated JSON as dependencies of this follow-up.
_AUTHENTICATED_SCRIPTS = (
    ROOT / "scripts/run_nsc_zeta1_regulated_determinant.py",
    ROOT / "scripts/run_nsc_zeta1_warped_y.py",
)
A = 3 * pi / 2
SOURCE_HASHES = {
    "scripts/run_nsc_zeta1_regulated_determinant.py":
        "01caf3d6061b1411a82ab030cdfe17fc8b0c5a9aa076d73ce8dbe589ff67276b",
    "scripts/run_nsc_zeta1_warped_y.py":
        "2206e30f746e3c6c0e04a65c9cead46122c6b97d9fd745b591f3243b43fab67b",
    "results/nsc-1-s-one-child-scale-correction.json":
        "19c50a14a9902b372a1868c147dd1506570b40a1ac3bd1f1200d03cfe5b8763b",
    "results/nsc-2-zeta1-recursion-map.json":
        "4b2c7c3cea91bc980a3de3750e89d32c1b3b105adbdaf7844be1ab70c9565774",
    "results/nsc-2-zeta1-anomaly-decomposition.json":
        "052750d49e969a1d82b683f6405a1568aff5f56352183eb9da948a6a2ab0b750",
}


def resolve_source(relative: str) -> Path:
    """Prefer frozen auxiliary copies during isolated reproduction.

    Historical v0.1.0 result JSON is not a graph dependency of this
    follow-up. Authenticated paths and hashes in the compact record remain
    the original laboratory paths.
    """
    auxiliary = AUXILIARY_ROOT / relative
    if auxiliary.is_file():
        return auxiliary
    return ROOT / relative


def authenticate() -> list[dict[str, str]]:
    inventory = []
    for relative, expected in SOURCE_HASHES.items():
        path = resolve_source(relative)
        if not path.is_file():
            raise RuntimeError(f"missing source: {relative}")
        observed = hashlib.sha256(path.read_bytes()).hexdigest()
        if observed != expected:
            raise RuntimeError(f"source drift: {relative}")
        inventory.append({"path": relative, "sha256": observed})
    return inventory


def exact_relations() -> dict:
    lam, zeta, a, omega, link, tail, parent = sp.symbols(
        "lambda zeta a Omega B Gamma_c K_p", positive=True
    )
    mu2 = 1 - a / zeta
    q_old = (lam + mu2) / zeta
    q_correct = lam / zeta + mu2
    # d E1(q) / dq = -exp(-q)/q.
    derivative = sp.factor(-sp.exp(-q_correct) * zeta * sp.diff(q_correct, zeta) / (2*q_correct))
    subtraction = sp.exp(-q_correct) / 2
    compensated = sp.simplify(derivative - subtraction)
    compensated_expected = -sp.exp(-q_correct) / (2*q_correct)
    residual = sp.simplify(compensated - compensated_expected)
    # Normalize a dimensional two-block inverse response by Lambda_parent.
    # Child block / Lambda_parent = Omega * Gamma_child(E/Lambda_child).
    block = sp.Matrix([[parent, link], [link, omega * tail]])
    exact_inverse = sp.factor(1 / block.inv()[0, 0])
    expected_inverse = parent - link**2 / (omega * tail)
    missing_factor = sp.factor(expected_inverse - (parent - link**2 / tail))
    # This finite complex-energy matrix also checks noncommuting block order.
    hp = sp.Matrix([[1, sp.Rational(1, 3)], [sp.Rational(1, 3), 2]])
    hc = sp.Matrix([[3, sp.I/4], [-sp.I/4, 4]])
    coupling = sp.Matrix([[sp.Rational(1, 5), sp.I/7], [sp.Rational(1, 6), 0]])
    energy = 2 + sp.I
    rp = energy * sp.eye(2) - hp
    rc = energy * sp.eye(2) - 2 * hc
    full = rp.row_join(-coupling).col_join((-coupling.conjugate().T).row_join(rc))
    schur = rp - coupling * rc.inv() * coupling.conjugate().T
    matrix_residual = (full.inv()[:2, :2] - schur.inv()).applyfunc(sp.simplify)
    # The strong particle/outside identity is verified directly, no stored labels.
    p, phi, e = sp.symbols("p Phi E", nonzero=True, real=True)
    d = sp.Matrix([[p, phi], [phi, -p]])
    visible = (e * sp.eye(2) - d).inv()[0, 0]
    gap_residual = sp.simplify(d**2 - (p**2 + phi**2)*sp.eye(2))
    visible_residual = sp.factor(1/visible - (e-p-phi**2/(e+p)))
    # Raw finite determinant invariance does not imply invariance of the
    # proper-time regulated determinant used in the later ZETA1 calculation.
    delta, cutoff = sp.symbols("delta Lambda", real=True, positive=True)
    radius = sp.sqrt(p**2*sp.cosh(delta)**2 + phi**2)
    squared_modes = [(-p*sp.sinh(delta) + sign*radius)**2/cutoff**2 for sign in (1,-1)]
    reg_hessian = 0
    for q in squared_modes:
        q0 = sp.simplify(q.subs(delta,0))
        q1 = sp.simplify(sp.diff(q,delta).subs(delta,0))
        q2 = sp.simplify(sp.diff(q,delta,2).subs(delta,0))
        reg_hessian += sp.exp(-q0)*((q0+1)*q1*q1/q0**2-q2/q0)/2
    reg_hessian = sp.simplify(reg_hessian)
    reg_expected = 4*p**2*sp.exp(-(p**2+phi**2)/cutoff**2)/cutoff**2
    assert sp.simplify(reg_hessian-reg_expected) == 0
    assert residual == 0 and sp.simplify(exact_inverse - expected_inverse) == 0
    assert matrix_residual == sp.zeros(2) and gap_residual == sp.zeros(2)
    assert visible_residual == 0
    return {
        "declared_units": {
            "lambda": "eigenvalue of the spatial squared operator times L_star^2",
            "zeta": "Lambda^2 * L_star^2",
            "mu_squared": "Phi^2 / Lambda^2 = 1 - (3*pi/2)/zeta",
            "gap_in_boundary_units": "Phi^2 * L_star^2 = zeta - 3*pi/2",
        },
        "old_heat_argument": str(sp.factor(q_old)),
        "unit_consistent_heat_argument": str(sp.factor(q_correct)),
        "argument_difference": str(sp.factor(q_correct - q_old)),
        "physical_mass_squared_implemented_by_old_argument": "Phi_declared^2 / zeta",
        "corrected_determinant_derivative": str(derivative),
        "inherited_cutoff_subtraction": str(subtraction),
        "corrected_subtracted_derivative": str(compensated_expected),
        "derivative_identity_residual": str(residual),
        "resolvent_scale_normalization": {
            "assumptions": "one common dimensional energy; first-order inverse response; b=B_dim/Lambda_parent",
            "dimensionless_child_block": "Omega * Gamma_c(x/Omega)",
            "parent_Schur_complement": "Gamma_p(x)=K_p(x)-(1/Omega)*b*Gamma_c(x/Omega)^(-1)*b_dagger",
            "symmetric_link_option": "b_sym=B_dim/sqrt(Lambda_parent*Lambda_child) absorbs 1/Omega",
            "unweighted_old_B_normalization_specified": False,
            "transformed_mode_argument": "E_child/Lambda_child=E_parent/Lambda_parent",
            "fixed_common_energy_argument": "E/Lambda_child=(E/Lambda_parent)/Omega",
            "scalar_residual": str(sp.simplify(exact_inverse - expected_inverse)),
            "difference_without_normalization": str(missing_factor),
            "noncommuting_finite_matrix_residual": str(matrix_residual),
            "normalization_selects_Omega": False,
        },
        "shared_Phi_identity": {
            "Dirac_squared": "(p^2+Phi^2)*I",
            "visible_inverse": "E-p-Phi^2/(E+p)",
            "spectral_square_residual": str(gap_residual),
            "Schur_residual": str(visible_residual),
            "electron_mass_determined": False,
        },
        "raw_versus_regulated_determinant": {
            "relative_operator": "D_delta=[[p*exp(-delta),Phi],[Phi,-p*exp(delta)]]",
            "raw_determinant": "-(p^2+Phi^2) independent of delta",
            "regulated_action": "one_half*sum E1(eigenvalue(D_delta)^2/Lambda^2)",
            "regulated_relative_Hessian_at_zero": str(reg_hessian),
            "Hessian_identity_residual": str(sp.simplify(reg_hessian-reg_expected)),
            "implication": "raw finite determinant invariance cannot justify discarding the regulated relative Hessian; a common regulator and anomaly must be specified",
        },
    }


def radial_spectrum(angular: int, intervals: int = 750) -> tuple:
    radius = 30.0
    h = radius / intervals
    rho = np.arange(-intervals+1, intervals, dtype=float) * h
    w = angular / np.sqrt(1 + rho*rho)
    dw = -angular*rho/(1 + rho*rho)**1.5
    joined, disconnected, margins = [], [], []
    for sign in (1, -1):
        diagonal = 2/h**2 + w*w + sign*dw
        off = np.full(len(diagonal)-1, -1/h**2)
        vals_j = eigh_tridiagonal(diagonal, off, eigvals_only=True)
        # A literal principal compression of the same joined matrix removes
        # the rho=0 node. No independent boundary spectrum is assumed.
        center = intervals-1
        vals_d = np.sort(np.concatenate([
            eigh_tridiagonal(diagonal[:center], off[:center-1], eigvals_only=True),
            eigh_tridiagonal(diagonal[center+1:], off[center+1:], eigvals_only=True),
        ]))
        margins.extend([np.min(vals_d-vals_j[:-1]), np.min(vals_j[1:]-vals_d)])
        joined.append(vals_j)
        disconnected.append(vals_d)
    return angular, np.concatenate(joined), np.concatenate(disconnected), float(min(margins))


def y_spectrum(intervals: int = 32) -> tuple[np.ndarray, np.ndarray]:
    h = 2/intervals
    nodes = np.linspace(-1, 1, intervals+1)
    edges = (nodes[:-1]+nodes[1:])/2
    incidence = (np.eye(intervals, intervals+1, k=1)-np.eye(intervals, intervals+1))/h
    weights = np.full(intervals+1, h)
    weights[[0, -1]] = h/2
    q = np.exp(9*edges**2/1015)[:, None] * incidence
    q *= (np.exp(9*nodes**2/1015)/np.sqrt(weights))[None, :]
    q *= np.sqrt(h)
    singular = np.sort(svdvals(q))
    return np.r_[0.0, singular**2], np.r_[1.0, np.full(intervals, 2.0)]


def relative_parts(zeta: float, radial: list, y: tuple, corrected: bool) -> np.ndarray:
    total = np.zeros(3)
    yvals, weights = y
    for angular, joined, disconnected, _ in radial:
        for radial_values, sign in ((joined, 1), (disconnected, -1)):
            lam = radial_values[:, None] + yvals[None, :]
            if corrected:
                quotient = 1 + (lam-A)/zeta
                gradient = (lam-A)/zeta
                coefficient = 1.0
            else:
                quotient = (lam+1)/zeta - A/zeta**2
                gradient = (lam+1)/zeta - 2*A/zeta**2
                coefficient = A/zeta**2
            if np.min(quotient) <= 0:
                raise RuntimeError("nonpositive heat argument")
            e = np.exp(-quotient)
            contributions = np.array([
                np.sum(weights*exp1(quotient))/2,
                np.sum(weights*e*gradient/quotient)/2,
                -coefficient*np.sum(weights*e/quotient)/2,
            ])
            total += sign * 4*angular * contributions
    return total


def calculate() -> dict:
    before = authenticate()
    algebra = exact_relations()
    with ThreadPoolExecutor(max_workers=4) as pool:
        radial = list(pool.map(radial_spectrum, range(1, 13)))
    y = y_spectrum()
    prior = json.loads(
        resolve_source("results/nsc-2-zeta1-anomaly-decomposition.json").read_text()
    )
    old_zeta = prior["prior_determinant_candidate"]["zeta"]
    old = relative_parts(old_zeta, radial, y, False)
    fixed = relative_parts(old_zeta, radial, y, True)
    historical_error = abs(old[2] - prior["prior_determinant_candidate"]["anomaly_compensated_child_link_derivative"])
    # Independent principal-compression gridding differs at rounding level
    # from separately generated historical half-grids. Record the discrepancy;
    # 1e-6 is a numerical replication check, not a scientific admission margin.
    assert historical_error < 1e-6
    bracket = (A+1e-5, 8.0)
    diagnostic_root = brentq(lambda z: relative_parts(z, radial, y, True)[1], *bracket, xtol=1e-11)
    root_parts = relative_parts(diagnostic_root, radial, y, True)
    h = 1e-4
    fd = (relative_parts(diagnostic_root*np.exp(h), radial, y, True)[0]
          - relative_parts(diagnostic_root*np.exp(-h), radial, y, True)[0])/(2*h)
    assert abs(fd-root_parts[1]) < 1e-5
    grid = np.geomspace(A+1e-5, 200.0, 257)
    with ThreadPoolExecutor(max_workers=4) as pool:
        rows = list(pool.map(lambda z: relative_parts(float(z), radial, y, True), grid))
    compensated = [float(row[2]) for row in rows]
    # This sign also follows without a scan from principal-submatrix
    # eigenvalue interlacing and strict decrease of exp(-q)/q for q>0.
    assert max(compensated) < 0
    interlacing_margin = min(row[3] for row in radial)
    assert interlacing_margin > -1e-8
    convergence = []
    for intervals in (375, 750, 1500):
        sectors = radial if intervals == 750 else [radial_spectrum(n, intervals) for n in range(1,13)]
        value = relative_parts(old_zeta, sectors, y, True)
        convergence.append({"radial_half_intervals": intervals, "radial_spacing": 30/intervals,
                            "subtracted_derivative": float(value[2])})
    assert before == authenticate()
    return {
        "artifact_id": "NSC-2-ZETA1-UNIT-CLOSURE-CHECK",
        "schema": "NSC-2-ZETA1-UNIT-CLOSURE-CHECK-v1",
        "source_commit": ANCHOR,
        "classification": "declared_scale_units_require_a_mass_term_correction_and_the_corrected_frozen_geometry_family_has_strictly_negative_subtracted_derivative",
        "authenticated_inputs": before,
        "exact_relations": algebra,
        "retained_numerical_family": {
            "radial_radius": 30, "radial_half_intervals": 750, "angular_max": 12,
            "y_intervals": 32, "y_interval": [-1,1], "angular_weight": "4*n as in inherited proxy",
            "warp": "sigma=-18*y^2/1015",
            "radial_plus_y_factorization_retained": True,
            "lapse_shift_spin_connection_and_nonseparable_operator_included": False,
        },
        "comparison_at_frozen_old_candidate": {
            "zeta": old_zeta,
            "historical_subtracted_derivative": float(old[2]),
            "historical_reproduction_absolute_error": historical_error,
            "corrected_determinant_derivative": float(fixed[1]),
            "corrected_subtracted_derivative": float(fixed[2]),
        },
        "corrected_determinant_diagnostic": {
            "zeta": diagnostic_root, "derivative_residual": float(root_parts[1]),
            "subtracted_derivative_at_root": float(root_parts[2]),
            "independent_finite_difference_derivative": float(fd),
            "physical_scale_selected": False,
        },
        "corrected_scan": {
            "count": len(grid), "zeta_min": float(grid[0]), "zeta_max": float(grid[-1]),
            "minimum_subtracted_derivative": min(compensated),
            "maximum_subtracted_derivative": max(compensated),
            "rows": [{"zeta": float(z), "subtracted_derivative": d} for z,d in zip(grid,compensated)],
        },
        "radial_convergence_at_old_candidate": convergence,
        "finite_family_sign_theorem": {
            "premise": "disconnected matrix is joined principal compression; geometry frozen; additive mass term; inherited cutoff subtraction",
            "eigenvalue_interlacing": "lambda_j<=nu_j<=lambda_(j+1)",
            "positive_decreasing_weight": "f(lambda)=exp[-(lambda/ zeta +mu^2)]/(lambda/zeta+mu^2)",
            "strict_trace_sign": "Tr f(H_joined)-Tr f(H_disconnected)>0 for each finite sector",
            "conclusion": "dGamma_subtracted/dlog(zeta)<0 for every zeta>3*pi/2 in this finite frozen-geometry family",
            "numerical_interlacing_margin": interlacing_margin,
            "continuum_or_general_recursive_theorem": False,
        },
        "next_required_equation": "derive the complete normalized curved K(E,zeta), B(E,zeta) and the anomaly-consistent variational derivative before selecting a physical recursive scale",
        "nonclaims": {
            "Nested_Space_hypothesis_proved_or_refuted": False,
            "complete_recursive_tail_solved": False,
            "physical_zeta_or_particle_mass_derived": False,
            "inherited_cutoff_subtraction_derived_as_full_local_anomaly": False,
            "new_to_world_physics_discovery_established": False,
        },
        "terminal": True,
    }


def compare(expected, observed, pointer="") -> None:
    if isinstance(expected, bool) or isinstance(observed, bool):
        if expected is not observed:
            raise RuntimeError(f"boolean mismatch at {pointer}")
    elif isinstance(expected, (int,float)) and isinstance(observed, (int,float)):
        if not np.isclose(expected,observed,rtol=1e-8,atol=1e-8):
            raise RuntimeError(f"numeric mismatch at {pointer}")
    elif isinstance(expected,dict) and isinstance(observed,dict):
        if expected.keys()!=observed.keys():
            raise RuntimeError(f"keys differ at {pointer}")
        for key in expected:
            compare(expected[key],observed[key],f"{pointer}/{key}")
    elif isinstance(expected,list) and isinstance(observed,list):
        if len(expected)!=len(observed):
            raise RuntimeError(f"length mismatch at {pointer}")
        for i,(a,b) in enumerate(zip(expected,observed)):
            compare(a,b,f"{pointer}/{i}")
    elif expected!=observed:
        raise RuntimeError(f"value mismatch at {pointer}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path)
    parser.add_argument("--check",action="store_true")
    args = parser.parse_args()
    if args.check and args.output:
        parser.error("choose --output or --check")
    record = calculate()
    payload = json.dumps(record, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()
    if args.check:
        compare(json.loads(DEFAULT_RESULT.read_bytes()), record)
        print("Scale closure calculation reproduced; every recorded numeric and exact field checked.")
    elif args.output:
        with args.output.open("xb") as stream:
            stream.write(payload)
        print(json.dumps({"output":str(args.output),"classification":record["classification"],
                          "comparison":record["comparison_at_frozen_old_candidate"],
                          "diagnostic":record["corrected_determinant_diagnostic"]},indent=2))
    else:
        if OUTPUT.exists():
            raise RuntimeError("ZETA1 unit-closure-check output already exists")
        OUTPUT.write_bytes(payload)
        print(json.dumps({"output":str(OUTPUT.relative_to(ROOT)),
                          "classification":record["classification"],
                          "follow_up_commit":FOLLOW_UP_COMMIT},indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
