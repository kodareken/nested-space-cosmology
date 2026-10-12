#!/usr/bin/env python3
"""Compose the retained baseline and current-history interpolation enclosures."""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from recursive_horizons.evidence_io import publish_exclusive_file
from recursive_horizons.nsc_ks_ball_trajectory import restored_upper
from flint import ctx

BASE = ROOT/"results/development/nsc-ks-gate-budget-v2.json"
INTERPOLATION = ROOT/"results/development/nsc-ks-difference-interpolation-d48-v3.json"
OUTPUT = ROOT/"results/development/nsc-ks-gate-budget-v3.json"


def digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def outward(value):
    exact = Q(str(value))
    number = float(exact)
    if Q(number) < exact:
        number = math.nextafter(number, math.inf)
    return number


def compute():
    baseline = json.loads(BASE.read_text())
    interpolation = json.loads(INTERPOLATION.read_text())
    if baseline["schema"] != "NSC-KS-GATE-ERROR-BUDGET-v2":
        raise ValueError("versioned baseline budget required")
    if interpolation["schema"] != "NSC-KS-DIFFERENCE-ENERGY-INTERPOLATION-ACCURACY-v3":
        raise ValueError("difference-interpolant enclosure required")
    if (interpolation["reference_mode"] != "direct-original-energies"
            or interpolation["covered_positive_families"] != 60
            or interpolation["covered_signed_families"] != 120):
        raise ValueError("interpolation coverage or reference representation differs")
    for record in (baseline, interpolation):
        for name, expected in {**record["source_hashes"], **record["input_hashes"]}.items():
            if digest(ROOT/name) != expected:
                raise ValueError("budget dependency changed: " + name)
    bounds = []
    with ctx.workprec(160):
        for value in interpolation["total_N_beta_interpolation_error_upper"]:
            upper = restored_upper(value)
            number = float(upper)
            # Conversion to a report float must not weaken the proven upper.
            from flint import arb
            if not arb(number) >= upper:
                number = math.nextafter(number, math.inf)
            bounds.append(number)
    result = json.loads(json.dumps(baseline))
    result["schema"] = "NSC-KS-GATE-ERROR-BUDGET-v3"
    result["profile_identity"] = interpolation["profile_identity"]
    result["degree"] = interpolation["degree"]
    result["reference_mode"] = interpolation["reference_mode"]
    result["components"]["energy_interpolation"].update({
        "bound": bounds,
        "status": "ENCLOSED: ideal difference interpolation, not node or arithmetic error",
        "evidence": {"path": str(INTERPOLATION.relative_to(ROOT)),
                     "sha256": digest(INTERPOLATION),
                     "field": "total_N_beta_interpolation_error_upper"},
        "exact_upper_evidence": interpolation["total_N_beta_interpolation_error_upper"],
    })
    missing = [name for name, row in result["components"].items() if row["bound"] is None]
    result["missing_components"] = missing
    result["partial_known_error_sum"] = [
        math.nextafter(math.fsum(row["bound"][i] for row in result["components"].values()
                                if row["bound"] is not None), math.inf)
        for i in range(2)]
    result["full_error_sum"] = None
    result["status"] = "OPEN: baseline and ideal interpolation enclosed; seven components missing"
    result["source_hashes"] = {
        "scripts/derive_nsc_ks_gate_budget_v3.py": digest(__file__),
        "src/recursive_horizons/nsc_ks_ball_trajectory.py": digest(ROOT/"src/recursive_horizons/nsc_ks_ball_trajectory.py"),
        "requirements-validation.txt": digest(ROOT/"requirements-validation.txt"),
    }
    result["input_hashes"] = {str(path.relative_to(ROOT)): digest(path) for path in (BASE, INTERPOLATION)}
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    modes = p.add_mutually_exclusive_group(required=True)
    modes.add_argument("--record", action="store_true")
    modes.add_argument("--check", action="store_true")
    args = p.parse_args()
    result = compute()
    if args.record:
        publish_exclusive_file(ROOT, str(OUTPUT.relative_to(ROOT)),
            (json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+"\n").encode())
    elif result != json.loads(OUTPUT.read_text()):
        raise ValueError("component budget replay differs")
    print(json.dumps({key: result[key] for key in (
        "status", "partial_known_error_sum", "missing_components",
        "physical_EXISTENCE_certificate")}, indent=2))
