#!/usr/bin/env python3
"""Replayable coefficient proof with actual scientific dependencies only.

The v2/v3 development registers remain historical artifacts. Their recorded
working-tree HEAD was not a snapshot of the then-untracked source. This
successor recomputes the coefficient identities and pilot from their owners;
it neither imports those historical numbers nor rewrites their bytes.
"""
import argparse
from collections.abc import Mapping
from hashlib import sha256
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))
sys.path.insert(0, str(ROOT/"scripts"))
import derive_nsc_ks_current_uv_coefficient_control_v2 as P
from recursive_horizons.nsc_ks_finite_history_uv_coefficients import finite_history_uv_coefficient_report
from recursive_horizons.nsc_ks_evaluation_binding import implementation_hashes
from recursive_horizons.evidence_io import publish_exclusive_file

OUTPUT = ROOT/"results/development/nsc-ks-current-uv-coefficients-v4.json"
OWNERS = (
    "scripts/derive_nsc_ks_current_uv_coefficient_control_v4.py",
    "scripts/derive_nsc_ks_current_uv_coefficient_control_v2.py",
    "src/recursive_horizons/nsc_ks_finite_history_uv_coefficients.py",
)
INPUTS = (
    "results/development/nsc-ks-gate-history-lm-broyden.json",
    "results/development/nsc-ks-source-inventory.json",
    "results/development/nsc-ks-cutoff-bridge-control.json",
    "docs/nsc-source-cutoff-bridge.md",
    "docs/nsc-incoming-fixed-transfer.md",
    "requirements-validation.txt",
)


def plain(value):
    if isinstance(value, Mapping):
        return {key: plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [plain(item) for item in value]
    return value


def compute():
    report = plain(finite_history_uv_coefficient_report())
    pilot = P.current_history_pilot()
    pilot.pop("cpu_seconds", None)
    return {
        "schema": "NSC-KS-CURRENT-UV-COEFFICIENTS-v4",
        "status": "OPEN: leading paired vacuum coefficient verified; numerical tail remainder unbounded",
        "profile_identity": json.loads((ROOT/INPUTS[0]).read_text())["profile_identity"],
        "coefficient_proof": report,
        "pilot": pilot,
        "pilot_quadrature_error_enclosed": False,
        "numerical_C_M": None, "finite_cutoff_error_bound": None,
        "physical_EXISTENCE_certificate": False,
        "physical_NONEXISTENCE_certificate": False,
        "source_hashes": implementation_hashes(ROOT, OWNERS),
        "input_hashes": {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in INPUTS},
        "historical_note": (
            "Earlier v2/v3 development records are preserved. This proof is computed "
            "directly from the stated equations and data; those development records "
            "and their working-tree HEAD are not proof inputs."),
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    modes = p.add_mutually_exclusive_group(required=True)
    modes.add_argument("--record", action="store_true")
    modes.add_argument("--check", action="store_true")
    args = p.parse_args()
    if args.check and not OUTPUT.is_file():
        raise FileNotFoundError("UV coefficient record is missing; --check cannot create evidence")
    result = compute()
    if args.record:
        publish_exclusive_file(ROOT, str(OUTPUT.relative_to(ROOT)),
            (json.dumps(result, indent=2, sort_keys=True, allow_nan=False)+"\n").encode())
    elif result != json.loads(OUTPUT.read_text()):
        raise ValueError("UV coefficient replay differs")
    print(json.dumps({"status": result["status"],
        "jet_identities": len(result["coefficient_proof"]["identities"]["residuals"]),
        "sigma_identities": len(result["coefficient_proof"]["e_minus2_sigma"]["residuals"]),
        "angular_pairs": len(result["coefficient_proof"]["e_minus2_sigma"]["angular_pairs"]),
        "numerical_C_M": None, "physical_gate": "OPEN"}, indent=2))
