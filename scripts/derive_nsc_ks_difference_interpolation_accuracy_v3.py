#!/usr/bin/env python3
"""Ideal interpolation enclosure for U_g-U_ref, with direct reference values."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys
import time

from flint import arb, ctx

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
import derive_nsc_ks_energy_interpolation_accuracy_v2 as Full
from recursive_horizons.evidence_io import publish_exclusive_file
from recursive_horizons.nsc_ks_ball_trajectory import exact_upper, restored_upper


def difference_bound(full):
    """Transfer the stated full-operator lemma with a conservative factor four.

    ||d_E^n(U_g-U_ref)|| <= ||d_E^n U_g|| + ||d_E^n U_ref||.
    The common radius lower bound also covers U_ref, so the W remainder
    at most doubles. U_ref has zero envelope z derivative. Consequently
    the F and F_z interpolation errors each at most double. The owned
    quadratic stress-error polynomial has nonnegative coefficients and
    degree <=2 in these errors: multiplying its old bound by four is safe.
    The numerical reference solve error is separate and is not set to zero.
    """
    if full["schema"] != "NSC-KS-CURRENT-ENERGY-INTERPOLATION-ACCURACY-v2":
        raise ValueError("current-history full-operator interpolation lemma required")
    result = dict(full)
    result["schema"] = "NSC-KS-DIFFERENCE-ENERGY-INTERPOLATION-ACCURACY-v3"
    result["status"] = "COMPONENT: difference-interpolant enclosure; other gate errors OPEN"
    result["reference_mode"] = "direct-original-energies"
    result["transfer"] = {
        "W_error_factor": 2, "Wz_error_factor": 1,
        "F_error_factor_upper": 2, "Fz_error_factor_upper": 2,
        "quadratic_matter_error_factor_upper": 4,
        "reference_z_derivative": 0,
        "reference_numerical_error_bound": None,
        "measured_drift_subtracted": False,
    }
    rows, totals = [], [arb(0), arb(0)]
    with ctx.workprec(Full.BITS):
        for old in full["families"]:
            row = dict(old)
            row["full_operator_W_lemma"] = old["W_interpolation_remainder"]
            row["W_interpolation_remainder"] = exact_upper(
                2*restored_upper(old["W_interpolation_remainder"]))
            row["full_operator_matter_lemma"] = old["separate_signed_matter_error_uppers"]
            signed = {}
            for sign, previous in old["separate_signed_matter_error_uppers"].items():
                scaled = dict(previous)
                for i, name in enumerate(("N", "beta")):
                    scaled[name] = exact_upper(4*restored_upper(previous[name]))
                    totals[i] += restored_upper(scaled[name])
                signed[sign] = scaled
            row["separate_signed_matter_error_uppers"] = signed
            rows.append(row)
        result["families"] = rows
        result["total_N_beta_interpolation_error_upper"] = [exact_upper(value) for value in totals]
        result["within_provisional_allocation"] = [
            bool(value.upper() <= Full.ball(Full.Q(1, 10**12))) for value in totals]
    result["source_hashes"] = {
        **full["source_hashes"],
        str(Path(__file__).relative_to(ROOT)): sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--record", action="store_true")
    mode.add_argument("--check", action="store_true")
    parser.add_argument("--degree", type=int, default=48)
    parser.add_argument("--history", default=Full.HISTORY)
    args = parser.parse_args()
    output = ROOT/"results/development"/f"nsc-ks-difference-interpolation-d{args.degree}-v3.json"
    start = time.process_time()
    result = difference_bound(Full.compute(args.history, args.degree))
    if args.record:
        raw = (json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+"\n").encode()
        publish_exclusive_file(ROOT, str(output.relative_to(ROOT)), raw)
    elif result != json.loads(output.read_text()):
        raise ValueError("difference interpolation replay differs")
    print(json.dumps({
        "status": result["status"], "degree": args.degree,
        "N_beta_interpolation_upper": [float(restored_upper(v)) for v in result["total_N_beta_interpolation_error_upper"]],
        "within_allocation": result["within_provisional_allocation"],
        "CPU_seconds": time.process_time()-start,
    }, indent=2))
