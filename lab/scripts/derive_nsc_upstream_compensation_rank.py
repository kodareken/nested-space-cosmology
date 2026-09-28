"""Exact algebra for one upstream response rank gate; no field propagation."""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results/development/nsc-upstream-compensation-rank.json"
INPUTS = (
    "results/development/nsc-retarded-radial-response.json",
    "results/development/nsc-mode-resolved-cauchy-state.json",
    "docs/nsc-upstream-metric-compensation.md",
    "docs/nsc-incoming-fourier-matching.md",
    "src/recursive_horizons/nsc_incoming_cauchy_jets.py",
    "src/recursive_horizons/nsc_lorentzian.py",
)


def guard_labels(m, ell, energy):
    if not (Q(3, 2) < m < Q(8, 5) and 2 < ell < Q(9, 4)
            and Q(1, 2) < energy < Q(3, 5)):
        raise ValueError("selected labels outside proved rational guards")


def exact_identities():
    import sympy as s

    m, ell, energy, a, r = s.symbols("m ell E a r", positive=True)
    ap, rp = s.symbols("a_prime r_prime", real=True)
    v = s.Matrix([[-m, 0, 0], [ell/r, 0, -ell/r**2],
                  [-energy/a, energy/a**2, 0]])
    inverse = s.Matrix([[-1/m, 0, 0], [-a/m, 0, a**2/energy],
                        [-r/m, -r**2/ell, 0]])
    derivative = s.Matrix([[0, 0, 0], [-ell*rp/r**2, 0, 2*ell*rp/r**3],
                           [energy*ap/a**2, -2*energy*ap/a**3, 0]])
    assert s.simplify(v*inverse - s.eye(3)) == s.zeros(3)
    assert s.simplify(v.det() + m*ell*energy/(a**2*r**2)) == 0
    assert s.simplify(v.diff(a)*ap + v.diff(r)*rp - derivative) == s.zeros(3)
    return {"inverse_identity": True, "determinant_identity": True,
            "derivative_identity": True}


def build_record():
    radial = json.loads((ROOT / INPUTS[0]).read_text())
    channel = json.loads((ROOT / INPUTS[1]).read_text())["channels"][14]
    assert radial["control"]["family"] == "14_1"
    assert channel["index"] == 14
    assert channel["compact_level"] == channel["angular_level"] == 1
    labels = {"mass": Q.from_float(channel["compact_mass"]),
              "angular": Q.from_float(channel["angular_eigenvalue"]),
              "energy": Q.from_float(radial["control"]["energies"][3])}
    guard_labels(labels["mass"], labels["angular"], labels["energy"])
    # The defining pi/2 and sqrt(5) labels satisfy the same guards.
    assert Q(3, 2) < Q(314159, 200000) < Q(11, 7) < Q(8, 5)
    assert 2**2 < 5 < Q(9, 4)**2
    center, halfwidth = Q(51, 50), Q(1, 200)
    left, right = center-halfwidth, center+halfwidth
    rho_outer = Q(103, 100)
    assert 1 < left < right < rho_outer
    a2_lower = 3*((1+rho_outer**2)*(Q(314159, 400000)-Q(3, 203))
                  - rho_outer)-1
    assert a2_lower > Q(16, 25)
    assert Q(3, 2)*Q(22, 7)-4 < 1  # a(1)^2 < 1
    assert 1+rho_outer**2 < Q(21, 10)
    assert Q(9, 4)**2 / 2 < Q(8, 5)**2  # ell/r < 8/5
    assert 2**3 > Q(5, 2)**2  # r^3 > 5/2 when r^2 >= 2
    inverse2 = Q(82, 45)+4+Q(441, 400)
    v2 = 2*Q(8, 5)**2+Q(3, 4)**2+Q(15, 16)**2+Q(9, 8)**2
    derivative_bounds = [Q(9, 4)/2, Q(3, 5)*Q(15, 8)/Q(4, 5)**2,
                         2*Q(3, 5)*Q(15, 8)/Q(4, 5)**3,
                         2*Q(9, 4)/Q(5, 2)]
    assert derivative_bounds == [Q(9, 8), Q(225, 128), Q(1125, 256), Q(9, 5)]
    vp2 = sum(x*x for x in derivative_bounds)
    generator2 = (2*Q(8, 5)**2+Q(3, 4)**2)/Q(4, 5)**2
    assert inverse2 < 7 and 7*Q(3, 8)**2 < 1
    assert v2 < 9 and vp2 < 36 and generator2 < 9
    variation = halfwidth*(2*3*3+6)
    sigma = Q(3, 8)-variation
    assert sigma == Q(51, 200) > 0
    return {
        "schema": "nsc-upstream-compensation-rank/v1",
        "status": "PASS: normalized one-energy diagonal response has full real rank",
        "proof_kind": "symbolic identities and exact rational inequalities; no propagation",
        "selected_channel": "14_1", "selected_input": {k: str(v) for k, v in labels.items()},
        "input_hashes": {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in INPUTS},
        "identities": exact_identities(),
        "support": {"center": str(center), "halfwidth": str(halfwidth),
                    "left": str(left), "right": str(right),
                    "bump": "any smooth common nonnegative nonzero bump supported in this interval"},
        "bounds": {"a_squared_strict_lower": str(a2_lower),
                   "inverse_F_squared_strict_upper": str(inverse2),
                   "V_F_squared_strict_upper": str(v2),
                   "V_prime_F_squared_strict_upper": str(vp2),
                   "G_norm_squared_strict_upper": str(generator2),
                   "V_prime_entry_absolute_upper": list(map(str, derivative_bounds)),
                   "center_sigma_min_strict_lower": "3/8",
                   "normalized_average_deviation_strict_upper": str(variation),
                   "normalized_response_sigma_min_strict_lower": str(sigma)},
        "actual_response": "sigma_min > mu_b * 51/200; mu_b = w_hat(0) integral(b/a0) > 0",
        "scope": {"clock": "fixed reference a0", "amplitudes_selected": False,
                  "bump_shape_selected": False, "field_solve_performed": False,
                  "same_metric_coefficients_required_across_all_channels": True,
                  "full_C0_matching": "OPEN", "full_action_stationarity": "OPEN"},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--write", action="store_true")
    action.add_argument("--check", action="store_true")
    args = parser.parse_args()
    record = build_record()
    if args.write:
        with OUTPUT.open("x") as stream:
            stream.write(json.dumps(record, indent=2, sort_keys=True)+"\n")
    else:
        assert json.loads(OUTPUT.read_text()) == record, "rank receipt differs from exact replay"
    print(record["status"])


if __name__ == "__main__":
    main()
