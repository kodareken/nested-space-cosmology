#!/usr/bin/env python3
"""Reproduce necessary dark-fraction plateau conditions and analytic controls.

This is a requirement check, not an NSC cosmological solution. --output creates
its target exclusively; --check compares every field, including source hashes.
"""
from __future__ import annotations

import argparse
import cmath
import hashlib
import json
import math
from pathlib import Path

from scipy.integrate import quad
from scipy.optimize import brentq
import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results/nsc-5-plateau-conditions.json"
RTOL = 2e-10
ATOL = 3e-12
B0, D0, W0 = 0.05, 0.95, -14 / 19


def compare(expected, observed, pointer=""):
    """Check every value; provenance, formulas, scope and integers are exact."""
    if isinstance(expected, bool) or isinstance(observed, bool):
        if expected is not observed:
            raise RuntimeError(f"boolean mismatch at {pointer}")
    elif isinstance(expected, int) and isinstance(observed, int):
        if expected != observed:
            raise RuntimeError(f"integer mismatch at {pointer}")
    elif isinstance(expected, (int, float)) and isinstance(observed, (int, float)):
        if not math.isclose(expected, observed, rel_tol=RTOL, abs_tol=ATOL):
            raise RuntimeError(f"numeric mismatch at {pointer}")
    elif isinstance(expected, dict) and isinstance(observed, dict):
        if expected.keys() != observed.keys():
            raise RuntimeError(f"keys differ at {pointer}")
        for key in expected:
            compare(expected[key], observed[key], f"{pointer}/{key}")
    elif isinstance(expected, list) and isinstance(observed, list):
        if len(expected) != len(observed):
            raise RuntimeError(f"length differs at {pointer}")
        for i, (left, right) in enumerate(zip(expected, observed)):
            compare(left, right, f"{pointer}/{i}")
    elif type(expected) is not type(observed) or expected != observed:
        raise RuntimeError(f"value mismatch at {pointer}")


def exact_checks():
    b, d, h = sp.symbols("rho_b rho_d H", positive=True)
    w, transfer = sp.symbols("w_d Q", real=True)
    bp, dp = -3*b + transfer/h, -3*(1+w)*d - transfer/h
    ratio, fraction = d/b, d/(b+d)
    ratio_prime = sp.diff(ratio, b)*bp + sp.diff(ratio, d)*dp
    fraction_prime = sp.diff(fraction, b)*bp + sp.diff(fraction, d)*dp
    epsilon = transfer/(h*(b+d))
    deceleration = (1+3*w*fraction)/2
    residuals = {
        "ratio_from_two_continuity_equations": sp.simplify(
            ratio_prime + 3*w*ratio + (1+ratio)**2*epsilon),
        "fraction_from_two_continuity_equations": sp.simplify(
            fraction_prime + 3*w*fraction*(1-fraction) + epsilon),
        "fraction_from_deceleration": sp.simplify(
            fraction_prime - (1-fraction)*(1-2*deceleration) + epsilon),
    }
    n, k = sp.symbols("N k", nonnegative=True, positive=True)
    b0, d0 = sp.symbols("rho_b0 rho_d0", positive=True)
    w0 = sp.symbols("w0", real=True)
    density_b = b0*sp.exp(-3*n)
    density_d = d0*sp.exp(-3*n-3*w0*(1-sp.exp(-k*n))/k)
    wd = w0*sp.exp(-k*n)
    h2 = density_b+density_d
    q_h = -1-sp.diff(h2, n)/(2*h2)
    q_stress = (1+3*wd*density_d/h2)/2
    # Factor the dilution from the finite amplitude before taking limits;
    # mixed symbolic rates otherwise needlessly obscure the product limit.
    amplitude_limit = sp.limit(sp.simplify(h2*sp.exp(3*n)), n, sp.oo)
    tail = sp.symbols("tail", nonnegative=True)
    fraction_tail = d0*sp.exp(-3*w0*(1-tail)/k)/(b0+d0*sp.exp(-3*w0*(1-tail)/k))
    residuals.update({
        "exponential_history_baryon_conservation": sp.simplify(
            sp.diff(density_b, n)+3*density_b),
        "exponential_history_dark_conservation": sp.simplify(
            sp.diff(density_d, n)+3*(1+wd)*density_d),
        "exponential_history_Friedmann_vs_stress_deceleration": sp.simplify(q_h-q_stress),
        "observable_ratio_reconstruction": sp.simplify(
            h2*sp.exp(3*n)/b0-1-density_d/density_b),
        "exponential_history_ratio_limit": sp.simplify(
            sp.limit(density_d/density_b, n, sp.oo)-d0*sp.exp(-3*w0/k)/b0),
        "exponential_history_finite_amplitude_limit": sp.simplify(
            amplitude_limit-b0-d0*sp.exp(-3*w0/k)),
        "exponential_history_H_squared_limit": sp.limit(sp.exp(-3*n), n, sp.oo)*amplitude_limit,
        "exponential_history_deceleration_limit": sp.simplify(
            sp.limit((1+3*w0*tail*fraction_tail)/2, tail, 0)-sp.Rational(1, 2)),
    })
    r, wr, er = sp.symbols("r w_r epsilon_r", real=True)
    e = sp.symbols("epsilon", real=True)
    flow = -3*w*r-(1+r)**2*e
    derivative = sp.diff(flow, r)+sp.diff(flow, w)*wr+sp.diff(flow, e)*er
    residuals["conserved_baryon_fixed_point_derivative"] = sp.simplify(
        derivative.subs({e: 0, er: 0, w: 0})+3*r*wr)
    phi = (1+sp.sqrt(5))/2
    gamma = sp.I*phi
    residuals["spectral_fixed_point"] = sp.simplify(gamma-sp.I+1/gamma)
    residuals["spectral_contraction_derivative"] = sp.simplify(1/gamma**2+1/phi**2)
    assert all(value == 0 for value in residuals.values())
    return {
        "zero_residuals": {key: str(value) for key, value in residuals.items()},
        "equations": {
            "N": "ln(a); prime=d/dN; a=1 at the control snapshot",
            "energy_transfer_sign": "dot(rho_b)+3H*rho_b=Q; dot(rho_d)+3H*(rho_d+p_d)=-Q",
            "epsilon": "Q/[H*(rho_b+rho_d)]",
            "r_prime": "-3*w_d*r-(1+r)^2*epsilon",
            "f_prime": "-3*w_d*f*(1-f)-epsilon=(1-f)*(1-2*q_dec)-epsilon",
            "q_dec": "(1+3*w_d*f)/2=-1-d(ln H)/dN",
            "fixed_point_epsilon": "-3*w_d*r/(1+r)^2=(1-f)*(1-2*q_dec)",
            "autonomous_flow_derivative": str(derivative),
            "local_stability": "F'(r_star)<0 for scalar autonomous closure; full coupled physical Jacobian otherwise",
            "Q_zero_solution": "r(N)=r0*exp[-3*integral_0^N w_d(s) ds]",
            "finite_positive_ratio_condition": "integral_0^N w_d(s) ds converges to a finite real limit",
            "regular_Q_zero_fixed_point": "w_d(r_star)=0; F'(r_star)=-3*r_star*w_d'(r_star)",
            "regular_Q_zero_positive_root_stability": "w_d'(r_star)>0",
            "regular_Q_zero_late_deceleration": "q_dec -> 1/2; w_total -> 0",
            "flat_Q_zero_late_geometry": "H^2 proportional to a^(-3); a(t) proportional to t^(2/3); curvature -> 0",
            "continued_acceleration_at_fixed_fraction": "q_dec<0 requires epsilon_star=(1-f_star)*(1-2*q_dec)>0",
        },
        "regular_limit_assumption": "Limits are smooth asymptotic dynamical states; convergence of a function alone need not imply convergence of its derivative.",
    }


def history(n, k):
    """Analytic separately conserved fluid control in H0=1 units."""
    integrated_exponent = -3*W0*(-math.expm1(-k*n))/k
    ratio = 19*math.exp(integrated_exponent)
    fraction = ratio/(1+ratio)
    log_h2 = -3*n+math.log(B0+D0*math.exp(integrated_exponent))
    wd = W0*math.exp(-k*n)
    return {
        "N": float(n), "w_dark": wd, "ratio": ratio, "dark_fraction": fraction,
        "H_over_H0": math.exp(log_h2/2), "q_dec": (1+3*wd*fraction)/2,
        "log_r_over_r0": integrated_exponent,
    }


def exponential_controls():
    controls = []
    for k in (1, 4, 20):
        r_star = 19*math.exp(-3*W0/k)
        f_star = r_star/(1+r_star)
        samples, max_quad, max_q_error = [], 0., 0.
        for n in (0., .1, .25, .5, 1., 2., 4., 8.):
            row = history(n, k)
            integral, _ = quad(lambda s: W0*math.exp(-k*s), 0., n,
                               epsabs=2e-13, epsrel=2e-13)
            error = abs(row["log_r_over_r0"]+3*integral)
            max_quad = max(max_quad, error)
            delta = 1e-5/k
            # Complex-step differentiation of H avoids subtracting nearby
            # logarithms and is independent of the fluid-pressure expression.
            complex_n = complex(n, delta)
            complex_exponent = -3*W0*(1-cmath.exp(-k*complex_n))/k
            complex_log_h = (-3*complex_n+cmath.log(
                B0+D0*cmath.exp(complex_exponent)))/2
            q_difference = -1-complex_log_h.imag/delta
            max_q_error = max(max_q_error, abs(q_difference-row["q_dec"]))
            row["log_ratio_quadrature_error"] = error
            samples.append(row)
        future_integral, quadrature_error_estimate = quad(
            lambda s: W0*math.exp(-k*s), 0., math.inf,
            epsabs=2e-13, epsrel=2e-13)
        infinite_error = abs(future_integral-W0/k)
        end = brentq(lambda n: history(n, k)["q_dec"], 0., 4./k,
                     xtol=1e-14, rtol=1e-14)
        elapsed, _ = quad(lambda n: 1/history(n, k)["H_over_H0"], 0., end,
                          epsabs=2e-13, epsrel=2e-13)
        tails = []
        for decay_units in (0, 2, 4, 8, 12):
            n = decay_units/k
            log_tail = (-3*W0/k)*math.exp(-k*n)
            current_r = r_star*math.exp(-log_tail)
            exact_fraction_tail = (r_star*(-math.expm1(-log_tail))
                                   /((1+r_star)*(1+current_r)))
            # Logistic derivative f(1-f)<=1/4 gives a rigorous tail bound.
            bound = log_tail/4
            assert 0 <= exact_fraction_tail <= bound
            assert math.isclose(f_star-history(n, k)["dark_fraction"],
                                exact_fraction_tail, rel_tol=2e-8, abs_tol=4e-16)
            tails.append({"N": n, "remaining_log_ratio": log_tail,
                          "remaining_fraction": exact_fraction_tail,
                          "rigorous_fraction_tail_upper_bound": bound})
        assert max_quad < 2e-12 and infinite_error < 2e-12 and max_q_error < 5e-8
        first = samples[0]
        assert abs(first["q_dec"]+.55) < 1e-14
        assert first["H_over_H0"] == 1. and first["dark_fraction"] == .95
        controls.append({
            "k": k, "samples": samples,
            "r_star_exact": f"19*exp(42/(19*{k}))",
            "r_star": r_star, "f_star": f_star,
            "asymptotic_E_squared_a_cubed": B0*(1+r_star),
            "asymptotic_q_dec": .5,
            "asymptotic_H_over_H0": 0.,
            "acceleration_ends_N": end, "acceleration_ends_a": math.exp(end),
            "acceleration_ends_elapsed_H0_time": elapsed,
            "future_integral_w": future_integral,
            "future_integral_quadrature_error_estimate": quadrature_error_estimate,
            "maximum_finite_integral_identity_error": max_quad,
            "infinite_integral_identity_error": infinite_error,
            "maximum_independent_H_complex_step_q_error": max_q_error,
            "convergent_tails": tails,
        })
    return {
        "status": "illustrative future fluid histories; not derived NSC laws, fitted attractors, or perturbation-stability demonstrations",
        "shared_snapshot": {"rho_b0": B0, "rho_dark0": D0, "H0": 1.,
                            "f0": .95, "r0": 19, "w_dark0_exact": "-14/19",
                            "w_total0": -.7, "q_dec0": -.55},
        "w_history": "w_d(N)=(-14/19)*exp(-k*N), N>=0, Q=0",
        "analytic_density": "rho_d(N)=(19/20)*exp[-3N+42*(1-exp(-k*N))/(19*k)]",
        "analytic_ratio_limit": "r_star=19*exp[42/(19*k)]; f_star=r_star/(1+r_star)",
        "convergent_tail": "ln(r_star/r(N))=42*exp(-k*N)/(19*k); 0<=f_star-f(N)<=ln(r_star/r(N))/4",
        "meaning": "Identical present H, fraction and acceleration do not determine the future ratio or transition rate.",
        "histories": controls,
    }


def slow_decay_control():
    n_symbol = sp.symbols("N", nonnegative=True)
    h2 = sp.exp(-3*n_symbol)*(sp.Rational(1, 20)
                            +sp.Rational(19, 20)*(1+n_symbol)**3)
    hdot = sp.diff(h2, n_symbol)/2
    curvature_limits = {
        "H_squared": sp.limit(h2, n_symbol, sp.oo),
        "Ricci": sp.limit(6*(hdot+2*h2), n_symbol, sp.oo),
        "Kretschmann": sp.limit(12*((hdot+h2)**2+h2**2), n_symbol, sp.oo),
    }
    assert all(value == 0 for value in curvature_limits.values())
    assert sp.limit(-1-hdot/h2, n_symbol, sp.oo) == sp.Rational(1, 2)
    rows = []
    for n in (0., 1., 4., 16., 64.):
        integral, _ = quad(lambda s: -1/(1+s), 0., n,
                           epsabs=2e-13, epsrel=2e-13)
        exact = -math.log1p(n)
        assert abs(integral-exact) < 2e-12
        r = 19*(1+n)**3
        f = r/(1+r)
        rows.append({"N": n, "w_dark": -1/(1+n), "ratio": r,
                     "dark_fraction": f, "integral_w": integral,
                     "integral_identity_error": abs(integral-exact),
                     "q_dec": (1-3*f/(1+n))/2})
    return {
        "status": "independent analytic counterexample; not an NSC solution",
        "w_dark": "-1/(1+N)", "integral_w": "-ln(1+N) -> -infinity",
        "ratio": "19*(1+N)^3 -> infinity", "fraction_limit": 1.,
        "w_dark_limit": 0., "q_dec_limit": .5,
        "curvature_limits_exact": {key: str(value) for key, value in curvature_limits.items()},
        "meaning": "A pressure tending to zero and curvature flattening do not by themselves imply a finite ratio plateau.",
        "samples": rows,
    }


def lcdm_control():
    rows = []
    for n in (0., 1., 2., 4., 8.):
        matter = .3*math.exp(-3*n)
        h2 = matter+.7
        b = B0*math.exp(-3*n)
        d = .25*math.exp(-3*n)+.7
        ricci = 3*matter+12*.7
        kretschmann = 12*((.7-.5*matter)**2+h2**2)
        rows.append({"N": n, "H_over_H0": math.sqrt(h2),
                     "ratio": d/b, "dark_fraction": d/h2,
                     "q_dec": -1+1.5*matter/h2,
                     "Ricci_over_H0_squared": ricci,
                     "Kretschmann_over_H0_fourth": kretschmann})
    n = sp.symbols("N", real=True)
    h2 = sp.Rational(3, 10)*sp.exp(-3*n)+sp.Rational(7, 10)
    hdot = sp.diff(h2, n)/2
    ricci = 6*(hdot+2*h2)
    kretschmann = 12*((hdot+h2)**2+h2**2)
    assert sp.simplify(ricci-sp.Rational(9, 10)*sp.exp(-3*n)-sp.Rational(42, 5)) == 0
    assert sp.limit(ricci, n, sp.oo) == sp.Rational(42, 5)
    assert sp.limit(kretschmann, n, sp.oo) == sp.Rational(294, 25)
    return {
        "status": "rounded flat dust-plus-positive-Lambda control; no current parameter estimation",
        "snapshot": {"baryon_fraction": .05, "cold_dark_matter_fraction": .25,
                     "vacuum_fraction": .7},
        "ratio": "rho_dark/rho_b=5+14*exp(3N) -> infinity",
        "fraction_limit": 1., "H_limit_over_H0": math.sqrt(.7),
        "Ricci_limit_over_H0_squared_exact": "42/5",
        "Kretschmann_limit_over_H0_fourth_exact": "294/25",
        "curvature_convention": "signature (-,+,+,+), R=6*(dot(H)+2H^2), flat FLRW",
        "meaning": "A diverging density ratio has a diluting denominator and is not a curvature singularity; its fraction is bounded.",
        "samples": rows,
    }


def spectral_countermodel():
    phi = (1+math.sqrt(5))/2
    rows = []
    for depth in (0, 1, 2, 4, 8, 16, 32):
        imaginary = 1.
        for _ in range(depth):
            imaginary = 1+1/imaginary
        rows.append({"depth": depth, "Gamma_imaginary": imaginary,
                     "error_from_fixed_point": abs(imaginary-phi)})
    assert rows[-1]["error_from_fixed_point"] < 2e-13
    return {
        "status": "logical countermodel to recursion alone; not a combined self-sourced NSC solution",
        "parameters": {"Omega": 1, "b": 1, "K_at_i": "i"},
        "equation": "Gamma=i-1/Gamma",
        "physical_fixed_point": "Gamma_star=i*(1+sqrt(5))/2",
        "depth_iteration_derivative": "-1/phi^2",
        "depth_contraction_magnitude": 1/phi**2,
        "uncoupled_cosmic_history": "rho_b=exp(-3N), rho_dark=1, r=exp(3N), f -> 1",
        "cosmic_N_in_spectral_equation": False,
        "matter_density_derived_from_Gamma": False,
        "meaning": "A stable recursion fixed point in room depth imposes no cosmic-time ratio attractor without the metric stress and evolution equations.",
        "depth_samples": rows,
    }


def calculate():
    paths = ["scripts/check_nsc_plateau_conditions.py", "docs/nsc-plateau-conditions.md"]
    return {
        "artifact_id": "NSC-5-PLATEAU-CONDITIONS",
        "schema": "NSC-5-PLATEAU-CONDITIONS-v1",
        "scope": "Exact FLRW consistency requirements plus analytic future controls; full physical closure remains open.",
        "source_hashes": {path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest() for path in paths},
        "reproduction_tolerances": {"relative": RTOL, "absolute": ATOL,
                                    "non_numerical_fields": "exact; no excluded fields"},
        "assumptions": [
            "Homogeneous flat FLRW, constant G, positive densities, H>0, radiation negligible at late time.",
            "Visible means fixed-mass baryonic dust; dark is the combined effective unresolved stress, not a separately postulated substance.",
            "Total effective stress is conserved. Separately conserved baryons Q=0 are the baseline; a nonzero transfer must be derived.",
            "Asymptotic flattening means curvature invariants tend to zero, not a chosen coordinate chart or necessarily no expansion.",
            "The rounded 95% fraction is an illustrative current matching datum; no attractor at r=19 is asserted.",
        ],
        "exact": exact_checks(),
        "exponential_future_controls": exponential_controls(),
        "slow_decay_counterexample": slow_decay_control(),
        "LCDM_control": lcdm_control(),
        "spectral_countermodel": spectral_countermodel(),
        "observational_audit": {
            "required_input": "H(a), full state-defined metric-variation stress, and Q(a) from one common action and fixed parameters",
            "Q_zero_reconstructed_ratio": "r_H=E(a)^2*a^3/Omega_b0-1, E=H/H0",
            "Q_zero_reconstructed_fraction": "f_H=1-Omega_b0*a^(-3)/E(a)^2",
            "reconstructed_w_total": "-1-(2/3)*d(ln H)/dN",
            "reconstructed_w_dark": "w_total/f_H",
            "Q_zero_finite_plateau_condition": "E(a)^2*a^3 -> Omega_b0/(1-f_star); d(ln H)/dN -> -3/2",
            "independent_comparisons": [
                "Friedmann and Raychaudhuri residuals against the full stress",
                "Both continuity equations and the fraction derivative identity",
                "Present expansion and combined fraction, followed by predicted future shape and curvature",
                "Physical perturbations and full coupled attractor Jacobian for the same solution",
            ],
            "source_mapping": "No rho_dark, pressure, Q, or current observed density is identified with b, Phi, or a recursion eigenvalue by assertion.",
            "prediction_claim": False,
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--check", action="store_true")
    group.add_argument("--output", type=Path)
    args = parser.parse_args()
    record = calculate()
    if args.check:
        compare(json.loads(OUTPUT.read_text()), record)
        print("Plateau requirements and analytic controls reproduced; every field and source hash checked.")
    elif args.output:
        with args.output.open("x") as stream:
            json.dump(record, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write("\n")
        print(args.output)
    else:
        print(json.dumps({"artifact_id": record["artifact_id"],
                          "scope": record["scope"],
                          "future_controls": [{"k": row["k"], "f_star": row["f_star"],
                                               "acceleration_ends_N": row["acceleration_ends_N"]}
                                              for row in record["exponential_future_controls"]["histories"]]},
                         indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
