#!/usr/bin/env python3
"""Compose Gate Budget v5 from authenticated PASS enclosures only.

Predecessor v4 is authenticated and left byte-for-byte unchanged. No
component is filled from an OPEN, ENCLOSED, REUSE or diagnostic record.
The current production record therefore remains OPEN with null bounds.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from recursive_horizons.nsc_ks_evaluation_binding import write_bytes_atomic
from recursive_horizons.nsc_ks_gate_budget_v5 import (
    SCHEMA,
    admit_budget,
    compose_open_budget,
    validate_open_budget,
)

PREDECESSOR_V5 = ROOT / "results/development/nsc-ks-gate-budget-v5.json"
OUTPUT = ROOT / "results/development/nsc-ks-gate-budget-v5-review1.json"
BUDGET_V4 = ROOT / "results/development/nsc-ks-gate-budget-v4.json"
ENCLOSURE = ROOT / "results/development/nsc-ks-continuous-constraint-enclosure-v2.json"
REANCHOR = ROOT / "results/development/nsc-ks-gate-reanchor-v2.json"
OWNERS = (
    "scripts/derive_nsc_ks_gate_budget_v5.py",
    "src/recursive_horizons/nsc_ks_gate_budget_v5.py",
    "tests/test_nsc_ks_gate_budget_v5.py",
    "docs/nsc-ks-gate-budget-v5.md",
)
INPUTS = (
    str(BUDGET_V4.relative_to(ROOT)),
    str(ENCLOSURE.relative_to(ROOT)),
    str(REANCHOR.relative_to(ROOT)),
    "scripts/derive_nsc_ks_gate_budget_v2.py",
    "scripts/derive_nsc_ks_gate_budget_v4.py",
    "src/recursive_horizons/nsc_local_gate_certificate_v2.py",
)


def digest(path):
    return sha256((ROOT / path).read_bytes()).hexdigest()


def calculate():
    predecessor = json.loads(BUDGET_V4.read_text())
    enclosure = json.loads(ENCLOSURE.read_text())
    reanchor = json.loads(REANCHOR.read_text())
    if enclosure["production"]["between_node"] is not None:
        raise ValueError("v5 may not take a null between-node production bound as PASS")
    if enclosure["production"]["arithmetic"] is not None:
        raise ValueError("v5 may not take a null arithmetic production bound as PASS")
    if enclosure["diagnostic"]["certificate_use"] is not False:
        raise ValueError("diagnostic enclosure cannot fill v5")
    result = compose_open_budget(
        predecessor=predecessor,
        pass_enclosures={},
        source_identity=reanchor["source_identity"],
        predecessor_sha256=digest(str(BUDGET_V4.relative_to(ROOT))),
    )
    result["enclosure_record"] = {
        "path": str(ENCLOSURE.relative_to(ROOT)),
        "sha256": digest(str(ENCLOSURE.relative_to(ROOT))),
        "between_node": enclosure["production"]["between_node"],
        "arithmetic": enclosure["production"]["arithmetic"],
        "certificate_use": enclosure["production"]["certificate_use"],
    }
    result["record_revision"] = 2
    result["supersedes_record"] = {
        "path": str(PREDECESSOR_V5.relative_to(ROOT)),
        "sha256": digest(str(PREDECESSOR_V5.relative_to(ROOT))),
    }
    result["source_hashes"] = {path: digest(path) for path in OWNERS}
    result["input_hashes"] = {path: digest(path) for path in INPUTS}
    validate_open_budget(result)
    return result


def check(record):
    validate_open_budget(record)
    for path, expected in {**record["source_hashes"], **record["input_hashes"]}.items():
        if digest(path) != expected:
            raise ValueError("gate-budget v5 dependency changed: " + path)
    if record["components"]["between_node"]["bound"] is not None:
        raise ValueError("v5 between_node must stay null without a PASS enclosure")
    if record["components"]["arithmetic"]["bound"] is not None:
        raise ValueError("v5 arithmetic must stay null without a PASS enclosure")
    if record["physical_EXISTENCE_certificate"] or record["physical_NONEXISTENCE_certificate"]:
        raise ValueError("OPEN budget may not be reported as EXISTENCE or NON_EXISTENCE")
    try:
        admit_budget(record)
    except ValueError:
        pass
    else:
        raise ValueError("current OPEN v5 budget must fail admission")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--record", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = parser.parse_args()
    computed = calculate()
    if args.check and not OUTPUT.is_file():
        raise FileNotFoundError(
            "gate-budget v5 record is missing; --check cannot create evidence")
    if args.record:
        payload = (
            json.dumps(computed, indent=2, sort_keys=True, allow_nan=False) + "\n"
        ).encode()
        if OUTPUT.is_file():
            if OUTPUT.read_bytes() != payload:
                raise ValueError("existing gate-budget v5 record differs")
        else:
            write_bytes_atomic(OUTPUT, payload, exclusive=True)
        value = computed
    else:
        value = json.loads(OUTPUT.read_text())
        if value != computed:
            raise ValueError("gate-budget v5 replay differs")
    check(value)
    print(json.dumps({
        "status": value["status"],
        "schema": value["schema"],
        "missing_components": value["missing_components"],
        "filled_from_pass": value["filled_from_pass"],
        "full_error_sum": value["full_error_sum"],
        "v4_bytes_preserved": value["v4_bytes_preserved"],
        "physical_EXISTENCE_certificate": value["physical_EXISTENCE_certificate"],
    }, indent=2))
