#!/usr/bin/env python3
"""Record/replay the 257-node continuous enclosure and remaining arithmetic.

Production between-node and arithmetic components stay null until the
owned directed inputs exist. The diagnostic block is a manufactured
Lipschitz/exact-rational control with certificate_use=false.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from recursive_horizons.nsc_ks_assembly_arithmetic import (
    SCHEMA as ARITHMETIC_SCHEMA,
    diagnostic_control_arithmetic,
    production_assembly_arithmetic,
    validate_assembly_arithmetic_record,
)
from recursive_horizons.nsc_ks_continuous_constraint_enclosure import (
    CANDIDATE_OWNER,
    CELL_COUNT,
    EVALUATOR_OWNERS,
    GRID_OWNER,
    INTERVAL_LABEL,
    MISSING_PRODUCTION_INPUTS,
    PROFILE_IDENTITY,
    SCHEMA,
    STATE_LAW,
    VERIFICATION_NODE_COUNT,
    EnclosureBinding,
    diagnostic_linear_enclosure,
    production_between_node_component,
    require_family,
    validate_continuous_enclosure_record,
)
from recursive_horizons.nsc_ks_evaluation_binding import write_bytes_atomic
from recursive_horizons.nsc_ks_profile_identity import profile_identity
from recursive_horizons.nsc_ks_source_envelope import physical_incoming_interval
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily
from recursive_horizons.nsc_local_gate_certificate_v2 import ERROR_COMPONENTS

PREDECESSOR = ROOT / "results/development/nsc-ks-continuous-constraint-enclosure.json"
OUTPUT = ROOT / "results/development/nsc-ks-continuous-constraint-enclosure-v2.json"
HISTORY = ROOT / "results/development/nsc-ks-gate-history-lm-broyden.json"
REANCHOR = ROOT / "results/development/nsc-ks-gate-reanchor-v2.json"
BUDGET_V4 = ROOT / "results/development/nsc-ks-gate-budget-v4.json"
PHASE = ROOT / "results/development/nsc-ks-current-phase-accuracy-v3-n129.json"
GEOMETRY = ROOT / "results/development/nsc-ks-between-node-geometry.json"
OWNERS = (
    "scripts/derive_nsc_ks_continuous_constraint_enclosure.py",
    "src/recursive_horizons/nsc_ks_continuous_constraint_enclosure.py",
    "src/recursive_horizons/nsc_ks_assembly_arithmetic.py",
    "src/recursive_horizons/nsc_ks_continuous_constraint_bound.py",
    "tests/test_nsc_ks_continuous_constraint_enclosure.py",
    "tests/test_nsc_ks_assembly_arithmetic.py",
    "docs/nsc-ks-continuous-constraint-enclosure.md",
)
INPUTS = (
    str(HISTORY.relative_to(ROOT)),
    str(REANCHOR.relative_to(ROOT)),
    str(BUDGET_V4.relative_to(ROOT)),
    str(PHASE.relative_to(ROOT)),
    str(GEOMETRY.relative_to(ROOT)),
    *EVALUATOR_OWNERS,
    "src/recursive_horizons/nsc_local_incoming_family.py",
    "src/recursive_horizons/nsc_ks_profile_identity.py",
)


def digest(path):
    return sha256((ROOT / path).read_bytes()).hexdigest()


def plain(value):
    if isinstance(value, dict):
        return {str(key): plain(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [plain(item) for item in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def load_family():
    history = json.loads(HISTORY.read_text())
    family = LocalIncomingFamily(np.asarray(history["history"]["coefficients"], float))
    identity = profile_identity(family, include_normal_window=True)
    if identity != PROFILE_IDENTITY or history["profile_identity"] != PROFILE_IDENTITY:
        raise ValueError("current history is not profile 0b0e4ced")
    require_family(family)
    return family, history


def calculate():
    family, history = load_family()
    reanchor = json.loads(REANCHOR.read_text())
    budget_v4 = json.loads(BUDGET_V4.read_text())
    phase = json.loads(PHASE.read_text())
    geometry = json.loads(GEOMETRY.read_text())
    if budget_v4["schema"] != "NSC-KS-GATE-ERROR-BUDGET-v4":
        raise ValueError("immutable v4 budget required")
    if budget_v4["profile_identity"] != PROFILE_IDENTITY:
        raise ValueError("v4 budget profile mismatch")
    if phase["history_identity"] != PROFILE_IDENTITY:
        raise ValueError("phase-value record profile mismatch")
    source_identity = reanchor["source_identity"]
    hashes = {path: digest(path) for path in (*OWNERS, *INPUTS)}
    binding = EnclosureBinding(
        PROFILE_IDENTITY, STATE_LAW, source_identity,
        physical_incoming_interval(),
        tuple(hashes.items()))
    between = production_between_node_component(binding=binding, family=family)
    arithmetic = production_assembly_arithmetic()
    diagnostic_grid = diagnostic_linear_enclosure(family, binding)
    diagnostic_arithmetic = dict(diagnostic_control_arithmetic())
    record = {
        "schema": SCHEMA,
        "record_revision": 2,
        "supersedes_record": {
            "path": str(PREDECESSOR.relative_to(ROOT)),
            "sha256": digest(str(PREDECESSOR.relative_to(ROOT))),
        },
        "accountable_author": "Douglas Ek",
        "status": (
            "OPEN: 257-node / 256-cell continuous enclosure and remaining "
            "assembly arithmetic defined; production bounds stay null until "
            "directed nodal, derivative and interval-array inputs exist"
        ),
        "profile_identity": PROFILE_IDENTITY,
        "state_law": STATE_LAW,
        "source_identity": source_identity,
        "interval": [float(v).hex() for v in family.interval],
        "interval_label": INTERVAL_LABEL,
        "verification_node_count": VERIFICATION_NODE_COUNT,
        "cell_count": CELL_COUNT,
        "grid_owner": GRID_OWNER,
        "candidate_owner": CANDIDATE_OWNER,
        "evaluator_owners": list(EVALUATOR_OWNERS),
        "nodes_hex": diagnostic_grid["nodes_hex"],
        "nodes_digest": diagnostic_grid["nodes_digest"],
        "binding": binding.as_record(),
        "production": {
            "between_node": between["bound"],
            "arithmetic": arithmetic["bound"],
            "between_node_status": between["status"],
            "arithmetic_status": arithmetic["status"],
            "certificate_use": False,
            "missing_inputs": list(between["missing_inputs"]) + list(
                arithmetic["missing_inputs"]),
        },
        "missing_scientific_inputs": list(MISSING_PRODUCTION_INPUTS),
        "error_components_not_filled": ["between_node", "arithmetic"],
        "nine_component_order": list(ERROR_COMPONENTS),
        "phase_value_arithmetic_excluded": True,
        "phase_record": {
            "path": str(PHASE.relative_to(ROOT)),
            "sha256": digest(str(PHASE.relative_to(ROOT))),
            "history_identity": phase["history_identity"],
            "node_count": phase["node_count"],
        },
        "geometry_between_node_is_not_full_residual": True,
        "geometry_record": {
            "path": str(GEOMETRY.relative_to(ROOT)),
            "sha256": digest(str(GEOMETRY.relative_to(ROOT))),
            "full_between_node_remainder": geometry["full_between_node_remainder"],
        },
        "predecessor_v4": {
            "path": str(BUDGET_V4.relative_to(ROOT)),
            "sha256": digest(str(BUDGET_V4.relative_to(ROOT))),
            "schema": budget_v4["schema"],
            "missing_components": budget_v4["missing_components"],
        },
        "history": {
            "path": str(HISTORY.relative_to(ROOT)),
            "sha256": digest(str(HISTORY.relative_to(ROOT))),
            "schema": history["schema"],
        },
        "diagnostic": {
            "continuous": diagnostic_grid,
            "arithmetic": diagnostic_arithmetic,
            "certificate_use": False,
            "scope": "synthetic/control outputs only",
        },
        "samples_are_not_proof": True,
        "dense_sampling_used_as_bound": False,
        "observed_drift_used_as_bound": False,
        "physical_EXISTENCE_certificate": False,
        "physical_NONEXISTENCE_certificate": False,
        "arithmetic_schema": ARITHMETIC_SCHEMA,
        "source_hashes": {path: digest(path) for path in OWNERS},
        "input_hashes": {path: digest(path) for path in INPUTS},
        "reproducer": (
            "python3 scripts/derive_nsc_ks_continuous_constraint_enclosure.py --check"
        ),
    }
    validate_continuous_enclosure_record(record)
    validate_assembly_arithmetic_record(arithmetic)
    validate_assembly_arithmetic_record(diagnostic_arithmetic)
    return plain(record)


def check(record):
    validate_continuous_enclosure_record(record)
    if record["production"]["between_node"] is not None:
        raise ValueError("recorder may not invent the production between-node bound")
    if record["production"]["arithmetic"] is not None:
        raise ValueError("recorder may not invent the production arithmetic bound")
    if record["diagnostic"]["certificate_use"] is not False:
        raise ValueError("diagnostic enclosure cannot be used as a certificate")
    if record["physical_EXISTENCE_certificate"] or record["physical_NONEXISTENCE_certificate"]:
        raise ValueError("OPEN calculation may not be reported as PASS or NON_EXISTENCE")
    for path, expected in {**record["source_hashes"], **record["input_hashes"]}.items():
        if digest(path) != expected:
            raise ValueError("continuous-enclosure dependency changed: " + path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--record", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = parser.parse_args()
    computed = calculate()
    if args.check and not OUTPUT.is_file():
        raise FileNotFoundError(
            "continuous enclosure record is missing; --check cannot create evidence")
    if args.record:
        payload = (
            json.dumps(computed, indent=2, sort_keys=True, allow_nan=False) + "\n"
        ).encode()
        if OUTPUT.is_file():
            if OUTPUT.read_bytes() != payload:
                raise ValueError("existing continuous enclosure record differs")
        else:
            write_bytes_atomic(OUTPUT, payload, exclusive=True)
        value = computed
    else:
        value = json.loads(OUTPUT.read_text())
        if value != computed:
            raise ValueError("continuous enclosure replay differs")
    check(value)
    print(json.dumps({
        "status": value["status"],
        "schema": value["schema"],
        "profile_identity": value["profile_identity"],
        "verification_node_count": value["verification_node_count"],
        "cell_count": value["cell_count"],
        "production": value["production"],
        "missing_scientific_inputs": value["missing_scientific_inputs"],
        "diagnostic_certificate_use": value["diagnostic"]["certificate_use"],
        "physical_EXISTENCE_certificate": value["physical_EXISTENCE_certificate"],
    }, indent=2))
