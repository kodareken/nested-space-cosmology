#!/usr/bin/env python3
"""Recover bound v2 fields after the missing payload identity writer defect.

This is an explicit successor conversion, not cache compatibility. The
original records stay immutable. Only save_family may differ from the pinned
producer, all other numerical/source dependencies must match, and every
stored operator is re-contracted against the unchanged complete source.
"""
import argparse
import ast
import copy
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
import derive_nsc_ks_gate_value as V
from recursive_horizons.nsc_ks_energy_propagator import KSEnergyPropagator
from recursive_horizons.nsc_ks_evaluation_binding import (
    PRODUCTION_MODE, authenticate_payload, binding_digest, compare_bindings,
    compare_family_cache, family_key_text, payload_descriptor)
from recursive_horizons.nsc_ks_source_envelope import (
    KSEnvelopeBinding, _sample_axial_profiles, usual_axial_support)

OWNER = "scripts/derive_nsc_ks_gate_value.py"


def unchanged_numerical_code(previous, current):
    """Only serialization may differ; imports/globals/all math stay identical."""
    def without_writer(source):
        tree = ast.parse(source)
        writers = [n for n in tree.body if isinstance(n, ast.FunctionDef)
                   and n.name == "save_family"]
        if len(writers) != 1:
            raise ValueError("one known family writer required")
        tree.body.remove(writers[0])
        return ast.dump(tree, include_attributes=False)
    if without_writer(previous) != without_writer(current):
        raise ValueError("numerical code changed outside the metadata writer")


def corrected_metadata_for_validation(record, arrays):
    """Recognize precisely the missing redundant identity, without relaxing IO."""
    metadata = dict(arrays["metadata"])
    if "evaluation_identity" in metadata:
        raise ValueError("conversion requires the known omitted identity field")
    if (metadata.get("binding_digest") != record.get("binding_digest")
            or record.get("evaluation_identity") != record.get("binding_digest")):
        raise ValueError("record and payload binding identity differ")
    metadata["evaluation_identity"] = metadata["binding_digest"]
    return {**arrays, "metadata": metadata}


def restore_operator(ctx, key, payload, metadata):
    with np.load(payload, allow_pickle=False) as handle:
        fields = [np.array(handle[name]) for name in V.KEYS]
    tagged = []
    for batch, channel in ctx["archive"].family_entries(key):
        normalized = V._channel_record("_", {"_": channel})
        tag = f'{batch.group}_{batch.angular_sign}_E{batch.energy_sign:+d}:{batch.panel_name}'
        tagged.append((tag, batch, normalized))
    group = V.H.group_signed_operator_families(tagged)[0]
    first = next(batch for _, batch, _ in group["applies"])
    metric, solver = ctx["family"].metric(), ctx["solver"]
    w, u = _sample_axial_profiles(metric.directions, ctx["grid"])
    binding = KSEnvelopeBinding(
        ctx["grid"], w, u, tuple(metric.amplitudes),
        tuple(d.inner_radius for d in metric.directions),
        tuple(d.outer_radius for d in metric.directions), usual_axial_support(),
        first.rho_up, 1.0, solver["rtol"], solver["atol"], solver["max_step"])
    return KSEnergyPropagator(
        V.R.interpolation_interval(ctx["archived"][key]), *fields,
        first.mass, first.angular, first.rho_up, binding, metadata["operator"]["diagnostics"])


def convert(old_label, new_label, producer_commit):
    if old_label == new_label:
        raise ValueError("a distinct successor label is required")
    old_directory, _ = V.output_paths(old_label)
    new_directory, new_output = V.output_paths(new_label)
    if new_output.exists():
        raise FileExistsError("successor value result already exists")
    pinned = subprocess.check_output(
        ["git", "-C", str(ROOT), "rev-parse", "--verify", producer_commit+"^{commit}"],
        text=True).strip()
    old_code = subprocess.check_output(["git", "-C", str(ROOT), "show", f"{pinned}:{OWNER}"])
    old = V.load_persisted_binding(old_directory)
    if sha256(old_code).hexdigest() != old["source_hashes"][OWNER]:
        raise ValueError("pinned producer differs from the original cache")
    unchanged_numerical_code(old_code, (ROOT/OWNER).read_bytes())
    for name, checksum in old["source_hashes"].items():
        if name != OWNER and V.digest(name) != checksum:
            raise ValueError("scientific dependency changed: " + name)
    settings = {name: float.fromhex(v) if isinstance(v, str) and v.startswith("0x") else v
                for name, v in old["solver"].items()}
    ctx = V.context(old["node_count"], settings, old["history_path"])
    expected = V.bind_context(ctx)
    updated = copy.deepcopy(old)
    for section in ("source_hashes", "implementation_hashes"):
        updated[section][OWNER] = V.digest(OWNER)
    compare_bindings(updated, expected)
    new_directory.mkdir(parents=True, exist_ok=True)
    V.persist_binding(new_directory, expected)
    rows = []
    for key in sorted(ctx["archived"]):
        original, payload = V.family_paths(old_directory, key)
        if not original.exists():
            if payload.exists():
                raise ValueError("orphan original payload")
            continue
        started = time.process_time()
        record = json.loads(original.read_text())
        arrays = authenticate_payload(record, payload, ROOT, target=ctx["target"])
        corrected = corrected_metadata_for_validation(record, arrays)
        compare_family_cache(record, old, key, mode=PRODUCTION_MODE, payload_arrays=corrected)
        operator = restore_operator(ctx, key, payload, corrected["metadata"])
        if operator.digest != record["operator_digest"]:
            raise ValueError("restored numerical operator identity differs")
        total, signed, _ = V.contract_family(ctx, key, operator)
        if (not np.array_equal(total, arrays["matter_change"])
                or signed != corrected["metadata"]["family_records"]):
            raise ValueError("fresh contraction differs from original arrays or inventory")
        successor, successor_payload = V.family_paths(new_directory, key)
        inheritance = {
            "reason": "missing redundant evaluation_identity in bound v2 payload",
            "original_record": payload_descriptor(original, ROOT),
            "original_payload": payload_descriptor(payload, ROOT),
            "original_binding": payload_descriptor(V.binding_file(old_directory), ROOT),
            "producer_commit": pinned,
            "conversion_owner": payload_descriptor(Path(__file__), ROOT),
            "numerical_arrays_changed": False,
            "source_contraction_replayed_exactly": True,
        }
        if successor.exists():
            saved = json.loads(successor.read_text())
            compare_family_cache(saved, expected, key, mode=PRODUCTION_MODE,
                payload_arrays=authenticate_payload(saved, successor_payload, ROOT, target=ctx["target"]))
            if saved.get("inherited_from") != inheritance:
                raise ValueError("successor inheritance changed")
        else:
            saved = V.save_family(ctx, new_directory, key, total, signed, operator,
                time.process_time()-started, inherited_from=inheritance)
            with np.load(successor_payload, allow_pickle=False) as fresh, np.load(payload, allow_pickle=False) as prior:
                for name in (*V.KEYS, "matter_change"):
                    if not np.array_equal(fresh[name], prior[name]):
                        raise ValueError("successor serialization changed numerical arrays")
        rows.append({"family": list(key), "original": V.digest(original), "successor": V.digest(successor)})
    inspected = V.inspect_existing_families(ctx, new_directory, mode=PRODUCTION_MODE)
    return {"status": "PASS: successor metadata and unchanged numerical arrays verified",
            "converted_families": len(rows), "remaining_families": len(inspected["pending"]),
            "evaluation_identity": binding_digest(expected), "new_operator_solves": 0,
            "physical_local_gate": "OPEN", "rows": rows}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--from-label", required=True)
    parser.add_argument("--to-label", required=True)
    parser.add_argument("--producer-commit", required=True)
    args = parser.parse_args()
    result = convert(args.from_label, args.to_label, args.producer_commit)
    print(json.dumps({k: v for k, v in result.items() if k != "rows"}, indent=2))
