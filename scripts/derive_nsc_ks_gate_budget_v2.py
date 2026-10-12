#!/usr/bin/env python3
"""Keep reusable baseline bounds separate from unresolved history errors."""
from hashlib import sha256
import argparse
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
from recursive_horizons.evidence_io import publish_exclusive_file

OUTPUT = ROOT / "results/development/nsc-ks-gate-budget-v2.json"
BASELINE = ROOT / "results/development/nsc-incoming-source-update-v5.json"
INVENTORY = ROOT / "results/development/nsc-ks-source-inventory.json"
LEGACY = ROOT / "results/development/nsc-ks-gate-budget-independent.json"
ALLOCATION = {
    "field_space_time": 5e-12,
    "changed_history_UV_tail": 5e-12,
    "baseline_low_subgap": 3e-12,
    "upstream": 2e-12,
    "energy_interpolation": 1e-12,
    "covered_regions": 1e-12,
    "phase_value": 1e-12,
    "between_node": 1e-12,
    "arithmetic": 1e-12,
}


def digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def assemble(baseline):
    """Import precisely the covered-domain result, never the full source."""
    budget = baseline["partial_error_budget"]
    if budget["constraint_order"] != ["N", "beta"] or budget["covered_domain_status"] != "PASS":
        raise ValueError("the declared V5 covered-domain certificate is required")
    if baseline["baseline"]["constraint_order"] != ["N", "beta"]:
        raise ValueError("baseline constraint order changed")
    covered = np.asarray(budget["covered_spectral_regions_action_error_upper"], float)
    if covered.shape != (2,) or not np.isfinite(covered).all() or np.any(covered < 0):
        raise ValueError("nonnegative finite covered bounds required")
    # These are bounds for the unchanged baseline. They say nothing about
    # the difference source, the newly evolved field, or omitted low/subgap.
    components = {
        name: {
            "allocation_per_component": allocation,
            "bound": None,
            "status": "OPEN: no verified bound for the current history",
        }
        for name, allocation in ALLOCATION.items()
    }
    components["covered_regions"].update({
        "bound": covered.tolist(),
        "status": ("REUSE: unchanged baseline covered-domain bound"
                   if np.all(covered <= ALLOCATION["covered_regions"])
                   else "REUSE: bound exceeds provisional allocation"),
        "allocation_excess": np.maximum(
            covered - ALLOCATION["covered_regions"], 0).tolist(),
        "scope": "only the unchanged V5 covered baseline spectral regions",
    })
    missing = [name for name, row in components.items() if row["bound"] is None]
    return {
        "schema": "NSC-KS-GATE-ERROR-BUDGET-v2",
        "status": "OPEN: retained baseline bound; history/source bounds missing",
        "constraint_order": ["N", "beta"],
        "components": components,
        "partial_known_error_sum": covered.tolist(),
        "full_error_sum": None,
        "missing_components": missing,
        "allocation_total": float(sum(ALLOCATION.values())),
        "residual_reserve": 1e-11,
        "physical_tolerance": 3e-11,
        "allocation_changed": False,
        "numerical_indicators_used_as_bounds": False,
        "physical_EXISTENCE_certificate": False,
        "physical_NONEXISTENCE_certificate": False,
        "named_gap": "field_source_and_continuous_residual_enclosures_missing",
    }


def compute():
    baseline = json.loads(BASELINE.read_text())
    legacy = json.loads(LEGACY.read_text())
    if legacy["allocation"] != ALLOCATION:
        raise ValueError("provisional allocation changed without a recorded decision")
    result = assemble(baseline)
    result["components"]["covered_regions"]["evidence"] = {
        "path": str(BASELINE.relative_to(ROOT)), "sha256": digest(BASELINE),
        "field": "partial_error_budget.covered_spectral_regions_action_error_upper",
    }
    result["source_hashes"] = {
        "scripts/derive_nsc_ks_gate_budget_v2.py": digest(__file__),
        "src/recursive_horizons/evidence_io.py": digest(
            ROOT / "src/recursive_horizons/evidence_io.py"),
    }
    result["input_hashes"] = {
        str(path.relative_to(ROOT)): digest(path) for path in (BASELINE, INVENTORY, LEGACY)
    }
    result["reuse_condition"] = (
        "same baseline action/source and spectral region allocation; no transfer "
        "of background bounds to changed-history matter or upstream propagation")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--record", action="store_true")
    modes.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = compute()
    if args.record:
        publish_exclusive_file(ROOT, str(OUTPUT.relative_to(ROOT)),
            (json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+"\n").encode())
    elif result != json.loads(OUTPUT.read_text()):
        raise ValueError("successor error-budget replay differs")
    print(json.dumps({name: result[name] for name in (
        "status", "partial_known_error_sum", "missing_components",
        "full_error_sum", "physical_EXISTENCE_certificate")}, indent=2))
