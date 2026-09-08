#!/usr/bin/env python3
"""Check a published torsion-flow projection before importing fixed points."""
import argparse
import hashlib
import json
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results/development/flow-compatibility.json"
SOURCES = ("scripts/check_nsc_flow_compatibility.py", "docs/nsc-flow-compatibility.md")


def hashes():
    return {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in SOURCES}


def calculate():
    # The source's independent contortion is antisymmetric in its last pair.
    b = {}
    for i in range(5):
        for j in range(5):
            b[i, j, j] = 0
            for k in range(j + 1, 5):
                b[i, j, k] = sp.Symbol(f"b_{i}{j}{k}", real=True)
                b[i, k, j] = -b[i, j, k]
    keys = list(b)
    skew = {key: (b[key] + b[key[1], key[2], key[0]] + b[key[2], key[0], key[1]]) / 3
            for key in keys}
    norm_b = sum(value ** 2 for value in b.values())
    swap = sum(b[i, j, k] * b[j, i, k] for i, j, k in keys)
    norm_skew = sum(value ** 2 for value in skew.values())
    residual = sp.expand(norm_skew - (norm_b - 2 * swap) / 3)
    orthogonal = sp.expand(sum(skew[key] * (b[key] - skew[key]) for key in keys))
    if residual != 0 or orthogonal != 0:
        raise AssertionError("connection projection convention")
    d = sp.Symbol("d", integer=True, positive=True)
    g1, g2, loop, fermion = sp.symbols("g1 g2 kappa_loop C_fermion", real=True)
    l1 = (d ** 2 - 7 * d - 12) / 4
    l2 = (d - 4) * (d + 1) / 4
    beta1 = -(d - 2) * g1 + loop * l1 * g1 - fermion
    beta2 = -(d - 2) * g2 + loop * l2 * g2 + 2 * fermion
    beta_u_on_subspace = sp.factor((2 * beta1 + beta2).subs(g2, -2 * g1))
    if sp.simplify(beta_u_on_subspace + 2 * (d + 2) * loop * g1) != 0:
        raise AssertionError("published flow projection")
    if sp.simplify(beta_u_on_subspace.subs(loop, 0)) != 0:
        raise AssertionError("fermion-loop subspace")
    c5 = sp.simplify(1 / ((4 * sp.pi) ** sp.Rational(5, 2) * sp.gamma(sp.Rational(5, 2))))
    if c5 != 1 / (24 * sp.pi ** 3):
        raise AssertionError("five-dimensional fermion coefficient")
    return {
        "schema": "NSC-PUBLISHED-FLOW-COMPATIBILITY-v2",
        "status": "exact source-compatibility test, not an NSC renormalization flow",
        "sources": ["https://arxiv.org/html/1506.02882v3", "https://arxiv.org/html/1101.1424v3"],
        "source_hashes": hashes(),
        "field_space": "source contortion beta_abc=-beta_acb; current candidate uses only K=beta_[abc]",
        "projection": {"norm_identity": "||K||^2=(I1-2*I2)/3", "norm_identity_residual": str(residual),
                       "orthogonality_residual": str(orthogonal),
                       "embedding_of_C_K_squared": "g1=C/3, g2=-2*C/3, g3=0",
                       "visible_subspace_condition": "u=2*g1+g2=0"},
        "imported_flow": {"dimensionless_gi": "gi/k^(d-2)",
                          "metric_lambda1": str(sp.factor(l1)), "metric_lambda2": str(sp.factor(l2)),
                          "fermion_offsets": "(-C_fermion,+2*C_fermion)",
                          "C_fermion_d5_one_Dirac": str(c5),
                          "beta_u_on_visible_subspace": str(beta_u_on_subspace),
                          "beta_u_d5_on_visible_subspace": str(sp.simplify(beta_u_on_subspace.subs(d, 5))),
                          "fermion_only_beta_u": "0"},
        "reuse_conditions": ["published generalized Palatini field space",
                             "connection has mass terms only and no regulating kernel",
                             "one-loop approximation; specified graviton gauge and cutoff",
                             "NSC proper-time regulator, boundary domain and constrained measure not matched"],
        "conclusions": {"quoted_Palatini_flow_preserves_selected_subspace": False,
                        "published_fermion_contribution_preserves_visible_subspace": True,
                        "this_is_a_derived_NSC_beta_function": False,
                        "finite_torsion_matching_datum_fixed": False,
                        "published_d4_fixed_point_is_NSC_d5_prediction": False},
        "next_action": "define the flow Hessian on metric, skew 3-form, physical fermions and links, with the actual measure/constraints before taking the trace",
        "comparison": "all exact fields and authenticated source hashes; no numerical evolution",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--check", action="store_true")
    group.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.output and args.output.exists():
        raise FileExistsError("refusing to overwrite evidence")
    expected = None
    if args.check:
        expected = json.loads(OUTPUT.read_text())
        if expected["source_hashes"] != hashes():
            raise AssertionError("source authentication failed")
    result = calculate()
    if args.check:
        if json.dumps(expected, sort_keys=True) != json.dumps(result, sort_keys=True):
            raise AssertionError("exact fields differ")
        print("published torsion flow: projection and all exact fields reproduced")
    elif args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open("x") as handle:
            json.dump(result, handle, indent=2, sort_keys=True)
            handle.write("\n")
        print(args.output)
    else:
        print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
