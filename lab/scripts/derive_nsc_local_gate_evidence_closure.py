#!/usr/bin/env python3
"""Build or replay the complete declared dependency graph of local evidence."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))
from recursive_horizons.evidence_io import publish_exclusive_file
from recursive_horizons.nsc_local_gate_evidence import build_closure, verify_closure

OWNERS = (
    "scripts/derive_nsc_local_gate_evidence_closure.py",
    "src/recursive_horizons/nsc_local_gate_evidence.py",
    "src/recursive_horizons/evidence_io.py",
    "requirements-validation.txt",
)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    modes = p.add_mutually_exclusive_group(required=True)
    modes.add_argument("--record", action="store_true")
    modes.add_argument("--check", action="store_true")
    p.add_argument("--name", required=True)
    p.add_argument("--root", action="append", default=[])
    p.add_argument("--pin", action="append", default=[])
    args = p.parse_args()
    if not re.fullmatch("[a-z0-9-]+", args.name):
        p.error("--name must contain lowercase letters, digits or hyphens")
    path = ROOT/"results/development"/("nsc-local-gate-evidence-"+args.name+".json")
    if args.record:
        if not args.root:
            p.error("--record requires at least one --root")
        manifest = build_closure(ROOT, [*args.root, *OWNERS], historical_commits=args.pin)
        manifest["scope"] = "declared dependency integrity; no physical result inferred"
        manifest["physical_claim_verified"] = False
        result = verify_closure(ROOT, manifest)
        raw = (json.dumps(manifest, indent=2, sort_keys=True)+"\n").encode()
        publish_exclusive_file(ROOT, str(path.relative_to(ROOT)), raw)
    else:
        if args.root or args.pin:
            p.error("--check uses the recorded roots and explicit historical pins")
        manifest = json.loads(path.read_text())
        result = verify_closure(ROOT, manifest)
    print(json.dumps({"integrity_status": result["integrity_status"],
        "nodes": len(result["verified_nodes"]),
        "historical_nodes": sum(n.get("git_commit") is not None for n in manifest["nodes"]),
        "manifest_sha256": sha256(path.read_bytes()).hexdigest(),
        "physical_claim_verified": False}, indent=2))


if __name__ == "__main__":
    main()
