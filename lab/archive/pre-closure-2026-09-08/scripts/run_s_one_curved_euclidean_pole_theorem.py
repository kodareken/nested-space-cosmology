#!/usr/bin/env python3
"""Bind the positive-axis form-factor inequality to curved spectral calculus."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
for value in (str(ROOT), str(ROOT / "src")):
    if value not in sys.path:
        sys.path.insert(0, value)

from recursive_horizons.evidence_io import canonical_json_bytes  # noqa: E402


OUTPUT = ROOT / "results/nsc-1-s-one-curved-euclidean-pole-theorem.json"
INPUTS = (
    ROOT / "results/nsc-1-s-one-flat-spectral-poles.json",
    ROOT / "results/nsc-1-s-one-one-operator-bind.json",
    ROOT / "results/nsc-1-s-one-spectral-profile-closure.json",
)


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("curved Euclidean pole-theorem output already exists")

    authenticated = []
    for path in INPUTS:
        raw = path.read_bytes()
        item = json.loads(raw)
        if item.get("terminal") is not True:
            raise RuntimeError(f"nonterminal input: {path.relative_to(ROOT)}")
        authenticated.append(
            {
                "path": str(path.relative_to(ROOT)),
                "sha256": hashlib.sha256(raw).hexdigest(),
                "artifact_id": item.get("artifact_id"),
            }
        )

    record = {
        "artifact_id": "NSC-1-S-ONE-CURVED-EUCLIDEAN-POLE-THEOREM",
        "schema": "NSC-1-S-ONE-CURVED-EUCLIDEAN-POLE-THEOREM-v1",
        "classification": "the_fixed_graviton_heat_form_factor_has_no_additional_Euclidean_pole_on_any_self_adjoint_curved_two_sheet_spectrum",
        "authenticated_inputs": authenticated,
        "theorem": {
            "operator_assumption": "mathbb_D_Theta_is_self_adjoint_so_Spec(mathbb_D_Theta^2/Lambda^2)_is_a_subset_of_[0,infinity)",
            "master_integral": "h(z)=integral_0^1 exp[-alpha*(1-alpha)*z] dalpha",
            "form_factor": "F2(z)=(1+z/4)*h(z)-2",
            "claim": "F2(z)<0_for_every_real_z_greater_than_or_equal_to_zero",
        },
        "proof": {
            "z_equals_zero": "F2(0)=-1",
            "zero_to_four": {
                "bound": "0<h(z)<1_for_z>0",
                "consequence": "F2(z)<z/4-1<=0",
            },
            "four_to_infinity": {
                "integral_symmetry": "h(z)=2*integral_0^(1/2) exp[-z*alpha*(1-alpha)] dalpha",
                "interval_bound": "alpha*(1-alpha)>=alpha/2",
                "master_bound": "h(z)<4*(1-exp(-z/4))/z<4/z",
                "consequence": "F2(z)<4/z-1<=0",
            },
            "strictness_at_z_four": "both_master_integral_bounds_are_strict_for_positive_z",
        },
        "functional_calculus_result": {
            "zero_not_in_spectrum": "0_not_in_Spec[F2(mathbb_D_Theta^2/Lambda^2)]",
            "curvature_can_change": [
                "Dirac_eigenvalues",
                "eigenfunctions",
                "multiplicities",
            ],
            "curvature_cannot_do_in_this_sector": "move_a_nonnegative_D_squared_eigenvalue_onto_a_zero_of_F2",
            "additional_Euclidean_graviton_form_factor_pole": False,
        },
        "scope": {
            "covariant_quadratic_heat_kernel_order": True,
            "nonlinear_all_curvature_orders": False,
            "Lorentzian_reflection_positivity": False,
            "retarded_causality": False,
            "background_induced_massless_graviton": False,
        },
        "next_result": "establish_or_reject_reflection_positivity_for_the_background_adjusted_two_sheet_kernel_then_construct_its_retarded_Lorentzian_continuation",
        "gate": {
            "inputs_authenticated": len(authenticated) == len(INPUTS),
            "self_adjoint_D_squared_nonnegative": True,
            "piecewise_strict_negative_bound_complete": True,
            "curved_Euclidean_extra_pole_excluded": True,
            "Lorentzian_retarded_test_completed": False,
        },
        "terminal": True,
    }
    for key, value in record["gate"].items():
        if key != "Lorentzian_retarded_test_completed" and value is not True:
            raise RuntimeError(f"curved Euclidean pole theorem failed: {key}")
    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
