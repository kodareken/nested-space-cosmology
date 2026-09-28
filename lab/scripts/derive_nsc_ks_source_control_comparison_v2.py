#!/usr/bin/env python3
"""Replayable same-source refinement comparisons; no field evolution."""
from hashlib import sha256
import argparse
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))
import derive_nsc_ks_source_control_v2 as C
from recursive_horizons.evidence_io import publish_exclusive_file


def compute(left_name, right_name):
    result = C.compare(left_name, right_name)
    left, a = C.load(left_name)
    right, b = C.load(right_name)
    result["schema"] = "NSC-KS-ORIGINAL-SOURCE-COMPARISON-v2"
    result["profile_identity"] = left["profile_identity"]
    result["left_solver"], result["right_solver"] = left["solver"], right["solver"]
    paired = {}
    for family in left["families"]:
        group, angular = family["family"]
        prefix = f"{group}_{angular}/"
        delta = sum(a[prefix + sign + "/matter_change"] - b[prefix + sign + "/matter_change"]
                    for sign in ("positive", "negative"))
        paired[f"{group}_{angular}"] = np.max(np.abs(delta), axis=0).tolist()
    result["paired_matter_max_difference"] = paired
    # This triangle sum is still only a refinement indicator for the selected
    # energy rows; it is not an error bound on the original full source.
    result["selected_source_movement_indicator"] = np.sum(list(paired.values()), axis=0).tolist()
    result["constraint_order"] = ["N", "beta"]
    result["full_source_covered"] = False
    result["search_numerics_declared"] = False
    result["certificate_numerics_declared"] = False
    result["source_hashes"] = {
        str(Path(__file__).relative_to(ROOT)): C.digest(__file__),
        "scripts/derive_nsc_ks_source_control_v2.py": C.digest(ROOT / "scripts/derive_nsc_ks_source_control_v2.py"),
    }
    result["input_hashes"] = {
        str(C.output_paths(name)[0].relative_to(ROOT)): C.digest(C.output_paths(name)[0])
        for name in (left_name, right_name)
    }
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    mode = p.add_mutually_exclusive_group(required=True)
    mode.add_argument("--record", nargs=2, metavar=("LEFT", "RIGHT"))
    mode.add_argument("--check", action="store_true")
    p.add_argument("--name", required=True)
    args = p.parse_args()
    C.output_paths(args.name)  # validate the name
    path = ROOT / "results/development" / f"nsc-ks-source-comparison-v2-{args.name}.json"
    if args.record:
        result = compute(*args.record)
        raw = (json.dumps(result, sort_keys=True, indent=2, allow_nan=False)+"\n").encode()
        publish_exclusive_file(ROOT, str(path.relative_to(ROOT)), raw)
    else:
        saved = json.loads(path.read_text())
        result = compute(saved["left"], saved["right"])
        if result != saved:
            raise ValueError("source comparison replay differs")
    print(json.dumps({key: result[key] for key in (
        "left", "right", "paired_matter_max_difference",
        "selected_source_movement_indicator", "indicator_not_a_bound",
        "full_source_covered")}, indent=2))
