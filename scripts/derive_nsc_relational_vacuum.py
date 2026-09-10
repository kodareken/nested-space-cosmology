#!/usr/bin/env python3
"""Apply the recursive zero-tadpole law and reopen the MMP vacuum gate."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from recursive_horizons.nsc_mmp_embedding import CommonCoefficients
from recursive_horizons.nsc_relational_vacuum import (
    RelationalVacuumNormalization,
    uniqueness_residual,
)


OUTPUT = ROOT / "results/development/relational-vacuum-normalization.json"
INPUTS = (
    "results/development/mmp-vacuum-owner.json",
    "results/development/mmp-embedding.json",
    "results/development/causal-common-functional.json",
)
SOURCES = (
    "src/recursive_horizons/nsc_relational_vacuum.py",
    "scripts/derive_nsc_relational_vacuum.py",
    "docs/nsc-relational-vacuum-normalization.md",
)


def hashes(paths):
    return {path: hashlib.sha256((ROOT / path).read_bytes()).hexdigest() for path in paths}


def compare(expected, actual, path="$"):
    if isinstance(expected, dict):
        if not isinstance(actual, dict) or expected.keys() != actual.keys():
            raise AssertionError(f"keys differ at {path}")
        for key in expected:
            compare(expected[key], actual[key], f"{path}/{key}")
    elif isinstance(expected, list):
        if not isinstance(actual, list) or len(expected) != len(actual):
            raise AssertionError(f"list differs at {path}")
        for index, (left, right) in enumerate(zip(expected, actual)):
            compare(left, right, f"{path}/{index}")
    elif isinstance(expected, float):
        if isinstance(actual, bool) or not isinstance(actual, (int, float)):
            raise AssertionError(f"numeric type differs at {path}")
        if abs(expected - actual) > 3e-13 + 3e-12 * abs(expected):
            raise AssertionError(f"numeric value differs at {path}")
    elif type(expected) is not type(actual) or expected != actual:
        raise AssertionError(f"exact value differs at {path}")


def calculate():
    owner, embedding, causal = (json.loads((ROOT / path).read_text()) for path in INPUTS)
    raw = owner["retained_partial_coefficients"]
    coefficients = CommonCoefficients(raw["A"], raw["C"], raw["V"])
    law = RelationalVacuumNormalization()
    normalized = law.normalize(coefficients)
    mmp = law.mmp_gate(coefficients)
    identities = uniqueness_residual()
    scale_defects = [law.scale_covariance_defect(raw["V"], 2.7, omega)
                     for omega in (0.5, 1.0, 3.0)]
    direct_sum = (
        law.counterterm(raw["V"], 2.0)
        + law.counterterm(raw["V"], 3.0)
        - law.counterterm(raw["V"], 5.0)
    )
    assert max(abs(value) for value in scale_defects) == 0
    assert abs(direct_sum) < 2e-16
    assert mmp["charged_seed_bound_passes"]
    assert mmp["asymptotically_flat_bulk_vacuum_passes"]
    assert all(mmp[name] for name in (
        "Einstein_coefficient_unchanged", "gauge_coefficient_unchanged",
        "charge_radius_unchanged", "throat_length_unchanged",
    ))
    return {
        "schema": "NSC-RELATIONAL-VACUUM-NORMALIZATION-v1",
        "status": (
            "recursive zero-tadpole normalization fixes the gravitating homogeneous "
            "V_full to zero while retaining gradient and boundary response"
        ),
        "source_hashes": hashes(SOURCES),
        "input_hashes": hashes(INPUTS),
        "law": {
            "projector": "P0 Gamma = homogeneous local zero-derivative metric tadpole",
            "normalized_action": "Gamma_rel=(I-P0) Gamma_one",
            "condition": "delta_g Gamma_rel=0 on every unlinked homogeneous room in inherited dimensionless units",
            "coefficient_action": "(V,A,C) maps to (0,A,C)",
            "CTP_term": "-V_full*(Vol_plus-Vol_minus), coefficient taken from the same complete action",
            "independent_coefficient_added": normalized.independent_coefficient_added,
            "applied_after_complete_same_action_assembly": True,
        },
        "classification": {
            "kind": "NSC recursive normalization law",
            "raw_spectral_a0_claimed_zero": False,
            "finite_cosmological_counterterm_fitted": False,
            "observational_value_used": False,
            "physical_statement": "homogeneous unlinked vacuum defines zero relational gravitational source",
            "consequence": "local a0 cannot select the inheritance scale; boundary, gradient and nonlocal terms must do so",
            "compact_modulus_rule": (
                "project only the local a0 density at fixed internal spectral data; "
                "retain the topology-dependent finite Casimir/link remainder"
            ),
        },
        "coefficient_result": normalized.to_dict(),
        "exact_identities": {
            **identities,
            "direct_sum_additivity_residual": direct_sum,
            "scale_covariance_defects": scale_defects,
            "equal_history_CTP_counterterm": law.ctp_counterterm(raw["V"], 4.2, 4.2),
        },
        "MMP_vacuum_gate": mmp,
        "gate": {
            "V_full_determined": True,
            "V_full": 0.0,
            "charged_seed_vacuum_condition_passes": True,
            "asymptotically_flat_MMP_vacuum_condition_passes": True,
            "A_and_C_preserved": True,
            "Casimir_link_and_curvature_terms_projected_out": False,
            "canonical_CTP_owner_preserved": causal["gate"]["finite_interface_executable"],
            "full_physical_MMP_embedding_complete": False,
        },
        "next_embedding_conditions": [
            "fix or declare the complete induced A_full and C_full normalization",
            "establish the physical closed magnetic return domain and covariance",
            "retain a massless lowest-Landau sector or compute the massive Casimir replacement",
        ],
        "independent_audit": {
            "a0_only_subtraction_preserves_A_and_C": True,
            "connected_minus_disconnected_would_cancel_all_local_Seeley_terms": True,
            "finite_topology_remainder_retained": True,
        },
        "untracked_normalization_or_Weyl_files_read_or_modified": False,
        "old_scientific_generators_rerun": False,
        "comparison": {
            "fields": "all", "exact": "keys, types, strings and source/input hashes",
            "float_atol": 3e-13, "float_rtol": 3e-12, "exceptions": [],
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--check", action="store_true")
    group.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = calculate()
    if args.check:
        compare(json.loads(OUTPUT.read_text()), result)
        print("relational vacuum normalization reproduced from authenticated inputs")
    elif args.output:
        if args.output.exists():
            raise FileExistsError("refusing to overwrite recorded evidence")
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
        print(args.output)
    else:
        print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
