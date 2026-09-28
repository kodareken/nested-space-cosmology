#!/usr/bin/env python3
"""Test the bosonic room covariance against a Stieltjes positivity condition."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
for value in (str(ROOT), str(ROOT / "src")):
    if value not in sys.path:
        sys.path.insert(0, value)

from recursive_horizons.evidence_io import canonical_json_bytes  # noqa: E402


OUTPUT = ROOT / "results/nsc-1-s-one-reflection-positivity-obstruction.json"
INPUTS = (
    ROOT / "results/nsc-1-s-one-background-adjusted-graviton.json",
    ROOT / "results/nsc-1-s-one-spectral-room-bind.json",
)


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("reflection-positivity obstruction output already exists")

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

    z = sp.symbols("z")
    master = sum(
        (-z) ** order
        / sp.factorial(order)
        * sp.factorial(order) ** 2
        / sp.factorial(2 * order + 1)
        for order in range(12)
    )
    kernel = sp.series((1 + z / 4) * master - 1, z, 0, 9).removeO()
    covariance = sp.series(1 / kernel, z, 0, 7)
    weighted_covariance = sp.series(z / kernel, z, 0, 7)
    second_at_zero = sp.diff(weighted_covariance.removeO(), z, 2).subs(z, 0)

    record = {
        "artifact_id": "NSC-1-S-ONE-REFLECTION-POSITIVITY-OBSTRUCTION",
        "schema": "NSC-1-S-ONE-REFLECTION-POSITIVITY-OBSTRUCTION-v1",
        "classification": "the_standalone_bosonic_heat_trace_graviton_covariance_fails_a_necessary_positive_Stieltjes_condition",
        "authenticated_inputs": authenticated,
        "exact_series": {
            "K2": str(kernel),
            "C2_equals_1_over_K2": str(covariance),
            "g_equals_z_times_C2": str(weighted_covariance),
            "g_second_derivative_at_zero": str(second_at_zero),
        },
        "necessary_Stieltjes_condition": {
            "representation": "C(z)=a/z+b+integral_0^infinity d_rho(t)/(z+t), d_rho>=0",
            "weighted_function": "g(z)=z*C(z)",
            "derived_second_derivative": "g''(z)=-2*integral_0^infinity t*d_rho(t)/(z+t)^3<=0",
            "observed_sign": "g''(0)=228/175>0",
            "condition_passed": False,
        },
        "causal_owner": {
            "insufficient_object": "standalone_background_adjusted_bosonic_TT_heat_trace",
            "required_complete_object": "gauge_invariant_covariance_from_the_full_fermionic_plus_bosonic_two_sheet_spectral_action",
            "missing_terms": [
                "fermionic_determinant",
                "off_diagonal_parent_child_spectrum",
                "background_consistent_contact_and_subtraction_terms",
            ],
            "coefficient_retuning_is_the_fix": False,
        },
        "scope": {
            "bosonic_TT_two_point_test": True,
            "gauge_invariant_full_gravity_observable_tested": False,
            "contact_terms_fully_removed": False,
            "complete_one_action_reflection_positivity_rejected": False,
        },
        "next_result": "include_the_fermionic_determinant_and_off_diagonal_sheet_spectrum_then_repeat_the_Stieltjes_test_on_a_complete_gauge_invariant_two_point_observable",
        "gate": {
            "inputs_authenticated": len(authenticated) == len(INPUTS),
            "series_inversion_exact": sp.simplify(
                kernel * covariance.removeO() - 1
            ).series(z, 0, 7).removeO()
            == 0,
            "second_derivative_exact": second_at_zero == sp.Rational(228, 175),
            "standalone_bosonic_Stieltjes_condition_passed": False,
            "full_one_action_test_completed": False,
        },
        "terminal": True,
    }
    for key in ("inputs_authenticated", "series_inversion_exact", "second_derivative_exact"):
        if record["gate"][key] is not True:
            raise RuntimeError(f"reflection-positivity derivation failed: {key}")
    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
