#!/usr/bin/env python3
"""Compose covered baseline, difference interpolation and nodal phase proofs.

Only scientific component records are computational inputs. The provisional
allocation is unchanged; the old incomplete search/bookkeeping record is
not a proof input for these three independently replayed components.
"""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
import json
import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))
sys.path.insert(0, str(ROOT/"scripts"))
import derive_nsc_ks_gate_budget_v2 as B
import derive_nsc_ks_current_phase_accuracy_v3 as P
from recursive_horizons.evidence_io import publish_exclusive_file
from recursive_horizons.nsc_local_gate_evidence import verify_closure

BASELINE = ROOT/"results/development/nsc-incoming-source-update-v5.json"
INTERPOLATION = ROOT/"results/development/nsc-ks-difference-interpolation-d48-v3.json"
PHASE = ROOT/"results/development/nsc-ks-current-phase-accuracy-v3-n129.json"
CLOSURE = ROOT/"results/development/nsc-local-gate-evidence-components-v3.json"
OUTPUT = ROOT/"results/development/nsc-ks-gate-budget-v4.json"


def digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def outward(value):
    exact = Q(value)
    rounded = float(exact)
    return math.nextafter(rounded, math.inf) if Q(rounded) < exact else rounded


def exact_packed(value):
    power = int(value["exponent"])
    return Q(int(value["mantissa"]))*Q(2)**power


def evidence(path, field):
    return {"path": str(path.relative_to(ROOT)), "sha256": digest(path), "field": field}


def compute():
    baseline, interpolation, phase = [json.loads(p.read_text()) for p in (BASELINE, INTERPOLATION, PHASE)]
    verify_closure(ROOT, json.loads(CLOSURE.read_text()))
    if (interpolation["schema"] != "NSC-KS-DIFFERENCE-ENERGY-INTERPOLATION-ACCURACY-v3"
            or interpolation["degree"] != 48
            or interpolation["reference_mode"] != "direct-original-energies"
            or interpolation["covered_positive_families"] != 60
            or interpolation["covered_signed_families"] != 120):
        raise ValueError("complete current difference-interpolation component required")
    if not P.nodal_component_pass(phase) or phase["node_count"] != 129:
        raise ValueError("all 129 current-history phase nodes must be enclosed")
    if phase["history_identity"] != interpolation["profile_identity"]:
        raise ValueError("component histories differ")
    for name, expected in phase["source_hashes"].items():
        if digest(ROOT/name) != expected:
            raise ValueError("phase source changed: " + name)
    directory, _ = P.paths("n129")
    P.authenticate_payloads(phase["payloads"], directory, tuple(range(129)))
    result = B.assemble(baseline)
    result["schema"] = "NSC-KS-GATE-ERROR-BUDGET-v4"
    result["profile_identity"] = interpolation["profile_identity"]
    result["degree"] = 48
    result["reference_mode"] = interpolation["reference_mode"]
    result["phase_nodes"] = 192
    result["phase_target_nodes_hex"] = phase["binding"]["target_node_hex"]
    result["components"]["covered_regions"]["evidence"] = evidence(
        BASELINE, "partial_error_budget.covered_spectral_regions_action_error_upper")
    interp = interpolation["total_N_beta_interpolation_error_upper"]
    phase_bounds = phase["quadrature_attempts"]["192"]["value_error_upper_N_beta"]
    for name, path, field, exact, scope in (
        ("energy_interpolation", INTERPOLATION, "total_N_beta_interpolation_error_upper", interp,
         "ideal difference interpolation; node solves and arithmetic excluded"),
        ("phase_value", PHASE, "quadrature_attempts.192.value_error_upper_N_beta", phase_bounds,
         "129 nodal phase values including phase contraction arithmetic; between-node excluded")):
        result["components"][name].update({
            "bound": [outward(exact_packed(v)) for v in exact],
            "status": "ENCLOSED: " + scope,
            "exact_upper_evidence": exact, "evidence": evidence(path, field),
        })
    components = result["components"]
    result["missing_components"] = [k for k, v in components.items() if v["bound"] is None]
    result["partial_known_error_sum"] = [outward(sum(
        (Q(row["bound"][i]) for row in components.values() if row["bound"] is not None), Q()))
        for i in range(2)]
    result["status"] = "OPEN: three error components enclosed; six remain missing"
    result["known_component_proofs_replayed"] = ["covered_regions", "energy_interpolation", "phase_value"]
    result["phase_continuous_remainder_included"] = False
    result["arithmetic_scope"] = "remaining geometry, field and assembly arithmetic; phase arithmetic is already in phase_value"
    owners = [Path(__file__), ROOT/"scripts/derive_nsc_ks_gate_budget_v2.py",
              ROOT/"scripts/derive_nsc_ks_current_phase_accuracy_v3.py"]
    result["source_hashes"] = {str(p.relative_to(ROOT)): digest(p) for p in owners}
    result["input_hashes"] = {str(p.relative_to(ROOT)): digest(p) for p in (BASELINE, INTERPOLATION, PHASE, CLOSURE)}
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
            (json.dumps(result, indent=2, sort_keys=True)+"\n").encode())
    elif result != json.loads(OUTPUT.read_text()):
        raise ValueError("component budget replay differs")
    print(json.dumps({key: result[key] for key in ("status", "partial_known_error_sum", "missing_components")}, indent=2))
