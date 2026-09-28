#!/usr/bin/env python3
"""Rank three saved histories under one corrected retained-source evaluator.

The ranking is of measured nodal residuals. It is not a continuum
certificate, a new optimizer step, or a transfer of one history's error
bounds to another history.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import platform
import sys

import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))
sys.path.insert(0, str(ROOT/"scripts"))
import derive_nsc_ks_gate_value as V
from recursive_horizons.evidence_io import publish_exclusive_file

OUTPUT = ROOT/"results/development/nsc-ks-gate-reanchor-v2.json"
LABELS = (
    "v2-current-dop-g1024-d48-metadata2",
    "v2-n64-best-dop-g1024-d48",
    "v2-trial-dop-g1024-d48",
)
CONTROLS = (
    "results/development/nsc-ks-family-comparison-v3-split-space-1024-2048.json",
    "results/development/nsc-ks-family-comparison-v3-split-cf4-dop.json",
    "results/development/nsc-ks-gate-budget-v4.json",
)


def digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def compute():
    records = [json.loads(V.output_paths(label)[1].read_text()) for label in LABELS]
    first = records[0]
    for row in records:
        if row["schema"] != "NSC-KS-GATE-VALUE-v2":
            raise ValueError("corrected bound evaluation required")
        if (row["completed_positive_families"] != 60
                or row["required_positive_families"] != 60
                or row["node_count"] != 129 or row["tangents"] != "zero"):
            raise ValueError("complete original-source 129-node values required")
        for key in ("solver", "source_identity", "merit_definition"):
            if row[key] != first[key]:
                raise ValueError("comparison changed source or numerical definition: " + key)
    measured = []
    hashes = {}
    for label, row in zip(LABELS, records):
        gradient = V.reconstruct_gradient(label, replay=False)
        maxima = np.max(np.abs(gradient), axis=0).tolist()
        if maxima != row["constraint_maxima"]:
            raise ValueError("reconstructed maxima differ")
        path = V.output_paths(label)[1]
        history = ROOT/row["history_path"]
        coefficients = json.loads(history.read_text())["history"]["coefficients"]
        measured.append({
            "label": label, "profile_identity": row["profile_identity"],
            "history_path": row["history_path"], "coefficients_per_function": len(coefficients[0]),
            "nodal_maxima": maxima, "merit": row["merit"],
            "value_record": {"path": str(path.relative_to(ROOT)), "sha256": digest(path)},
            "new_optimizer_step": False,
        })
        hashes[str(path.relative_to(ROOT))] = digest(path)
    best = min(measured, key=lambda r: r["merit"])
    for name in CONTROLS:
        hashes[name] = digest(ROOT/name)
    return {
        "schema": "NSC-KS-GATE-REANCHOR-v2",
        "status": "OPEN: comparable corrected nodal measurements; full error budget incomplete",
        "state_law": "C_Sigma[g]=U_g C_up U_g^dagger; C_up unchanged",
        "class": "delta r = chi(s)[s*w(z)+s^3*U(z)/6]",
        "interval": "S(1)+[0.12,0.18]", "node_count": 129,
        "source_identity": first["source_identity"], "solver": first["solver"],
        "merit_definition": first["merit_definition"], "rows": measured,
        "best_measured_label": best["label"], "best_measured_profile_identity": best["profile_identity"],
        "best_measured_merit": best["merit"], "physical_tolerance": 3e-11,
        "full_retained_source_covered": True,
        "comparison_is_of_certified_continuum_residuals": False,
        "history_specific_bounds_transferred": False,
        "physical_EXISTENCE_certificate": False, "physical_NONEXISTENCE_certificate": False,
        "certificate_numerics": None,
        "calibration_scope": "representative complete families on history 0b0e4ced; convergence indicators only",
        "runtime_versions": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
        "source_hashes": {str(Path(__file__).relative_to(ROOT)): digest(__file__),
                          "scripts/derive_nsc_ks_gate_value.py": digest(ROOT/"scripts/derive_nsc_ks_gate_value.py")},
        "input_hashes": hashes,
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    mode = p.add_mutually_exclusive_group(required=True)
    mode.add_argument("--record", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = p.parse_args()
    result = compute()
    if args.record:
        publish_exclusive_file(ROOT, str(OUTPUT.relative_to(ROOT)),
            (json.dumps(result, sort_keys=True, indent=2)+"\n").encode())
    elif result != json.loads(OUTPUT.read_text()):
        raise ValueError("reanchoring replay differs")
    print(json.dumps({key: result[key] for key in ("status", "best_measured_label", "best_measured_merit", "rows")}, indent=2))
