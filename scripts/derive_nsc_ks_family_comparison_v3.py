#!/usr/bin/env python3
"""Compare calibrated families across an explicitly pinned metadata-only repair."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
from types import FunctionType

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))
sys.path.insert(0, str(ROOT/"scripts"))
import derive_nsc_ks_evaluator_accuracy_v2 as A
from derive_nsc_ks_value_metadata_successor import OWNER, unchanged_numerical_code
from recursive_horizons.evidence_io import publish_exclusive_file


def compare(left, right, producer_commit):
    pinned = subprocess.check_output(
        ["git", "-C", str(ROOT), "rev-parse", "--verify", producer_commit+"^{commit}"],
        text=True).strip()
    old_code = subprocess.check_output(["git", "-C", str(ROOT), "show", pinned+":"+OWNER])
    current_code = (ROOT/OWNER).read_bytes()
    unchanged_numerical_code(old_code, current_code)
    old_hash, current_hash = sha256(old_code).hexdigest(), sha256(current_code).hexdigest()
    provenance = {}

    def checked_load(name):
        row = A.load(name, verify_sources=False)
        checksums = dict(row["source_hashes"])
        for path, expected in checksums.items():
            actual = A.digest(ROOT/path)
            if path == OWNER and expected == old_hash:
                checksums[path] = current_hash
            elif actual != expected:
                raise ValueError("calibration dependency changed: " + path)
        provenance[name] = {"original_value_writer_sha256": row["source_hashes"][OWNER],
                            "metadata_only_equivalence_verified": True}
        # This normalization applies only to A.compare's code-equality test.
        # Original records and source hashes remain immutable and are linked
        # as inputs below. All numerical functions/imports and data matched.
        return {**row, "source_hashes": checksums}

    namespace = {**A.compare.__globals__, "load": checked_load}
    runner = FunctionType(A.compare.__code__, namespace, A.compare.__name__)
    result = runner(left, right)
    result["schema"] = "NSC-KS-EVALUATOR-COMPARISON-v3"
    result["profile_identity"] = A.load(left, verify_sources=False)["profile_identity"]
    result["writer_equivalence"] = {
        "pinned_producer_commit": pinned, "old_writer_sha256": old_hash,
        "current_writer_sha256": current_hash, "records": provenance,
        "numerical_arrays_modified": False,
    }
    paths = [Path(__file__), ROOT/OWNER, ROOT/"scripts/derive_nsc_ks_evaluator_accuracy_v2.py",
             ROOT/"scripts/derive_nsc_ks_value_metadata_successor.py"]
    result["source_hashes"] = {str(p.relative_to(ROOT)): A.digest(p) for p in paths}
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    modes = p.add_mutually_exclusive_group(required=True)
    modes.add_argument("--record", action="store_true")
    modes.add_argument("--check", action="store_true")
    p.add_argument("--left", required=True)
    p.add_argument("--right", required=True)
    p.add_argument("--name", required=True)
    p.add_argument("--producer-commit", required=True)
    args = p.parse_args()
    A.V.require_label(args.name)
    output = ROOT/"results/development"/("nsc-ks-family-comparison-v3-"+args.name+".json")
    result = compare(args.left, args.right, args.producer_commit)
    if args.record:
        publish_exclusive_file(ROOT, str(output.relative_to(ROOT)),
            (json.dumps(result, indent=2, sort_keys=True)+"\n").encode())
    elif result != json.loads(output.read_text()):
        raise ValueError("family comparison replay differs")
    print(json.dumps({key: result[key] for key in (
        "profile_identity", "families", "indicator_not_a_bound", "physical_EXISTENCE_certificate")}, indent=2))
