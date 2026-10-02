#!/usr/bin/env python3
"""Freeze the numerical basis of existing nested-pair checkpoints; no evolution."""
from dataclasses import replace
import hashlib
import json
from pathlib import Path

import numpy as np
from recursive_horizons import nsc_nested_parent_child as model
from recursive_horizons import nsc_spherical_galerkin_coupling as galerkin

LAB = Path(__file__).resolve().parents[1]
OUT = LAB / "results/development/nsc-nested-parent-child-replay-basis-v1.json"
PAYLOAD = OUT.with_suffix(".npz")
INPUTS = [LAB / ("results/development/" + name + suffix)
          for name in ("nsc-nested-parent-child-v1", "nsc-nested-parent-child-confirmation-v2")
          for suffix in (".json", ".npz")]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def records():
    return [json.loads(path.read_text()) for path in INPUTS if path.suffix == ".json"]


def bindings():
    paths = INPUTS + [Path(__file__), LAB / "scripts/check_nsc_nested_parent_child_confirmation.py"]
    for record in records():
        paths += [LAB / name for name in record.get("source_bindings", record.get("source_bindings_before", {}))]
    return {p.relative_to(LAB).as_posix(): digest(p) for p in dict.fromkeys(paths)}


def frozen_pair(nf, weights, arrays, entry):
    """Use saved W and source/observer bytes; never reconstruct an SVD frame."""
    grid = galerkin.build_grid(nf, gauge="conformal")
    grid.fine = replace(grid.fine, occupations=np.array(weights))
    ref, src = arrays[f"nf{nf}_reference_columns"], arrays[f"nf{nf}_source_columns"]
    return model.NestedPair(grid, arrays[f"nf{nf}_W"],
        np.array(entry["coarse_indices"]), np.array(entry["child_indices"]),
        np.array(entry["parent_indices"]), ref[:nf], ref[nf:], ref,
        src[:nf], src[nf:], src, np.array(weights), entry["source"], entry["geometry"])


def validate(arrays, entries):
    checked, worst = 0, 0.
    for nf in (128, 256, 512):
        W, ref, src = (arrays[f"nf{nf}_{name}"] for name in ("W", "reference_columns", "source_columns"))
        if W.shape != (nf - 1, nf - 1) or np.max(abs(W.T @ W - np.eye(nf - 1))) > 1e-10:
            raise ValueError("replay basis is not the complete orthogonal canonical frame")
        for cols in (ref, src):
            if cols.shape != (2 * nf, 6) or np.max(abs(cols.conj().T @ cols - np.eye(6))) > 1e-10:
                raise ValueError("replay source/observer is not the preserved CAR frame")
        if not np.array_equal(src[:, 2:4], ref[:, 2:4]):
            raise ValueError("replay child source and observer differ")
    for record, payload in zip(records(), INPUTS[1::2]):
        if digest(payload) != record["payload_sha256"]:
            raise ValueError("resident scientific payload changed")
        with np.load(payload, allow_pickle=False) as stored:
            for name, case in record["results"].items():
                nf = case["nf"]
                pair = frozen_pair(nf, case["source"]["occupations"], arrays, entries[str(nf)])
                if not np.array_equal(stored[name + "_reference_columns"], pair.reference_columns):
                    raise ValueError("replay observer differs from the saved observer")
                initial = np.vstack((stored[name + "_phi0"][0], stored[name + "_phi1"][0]))
                if not np.array_equal(initial, pair.source_columns):
                    raise ValueError("replay preparation differs from the saved source")
                for i, row in enumerate(case["rows"]):
                    state = model.NestedState(*(stored[name + "_" + field][i] for field in model.STATE_NAMES))
                    metrics = model.metrics(pair, state)
                    for key in ("r_min", "r_max", "Q_min", "Q_max", "parent_proper_length",
                                "child_proper_length", "child_r_proper_mean", "parent_annulus_r_proper_mean"):
                        gap = abs(metrics[key] - row["metrics"][key])
                        worst = max(worst, gap)
                        if gap > 2e-9 * max(1., abs(row["metrics"][key])):
                            raise ValueError("replay basis changes saved full geometry: " + name + ":" + key)
                    if np.max(abs(np.array(metrics["clock_rates"]) - row["metrics"]["clock_rates"])) > 2e-9:
                        raise ValueError("replay basis changes saved normal clocks")
                    if i in (0, len(case["rows"]) - 1) and abs(model.energy(pair, state) - row["energy"]) > 2e-9:
                        raise ValueError("replay basis changes saved endpoint energy")
                    cols = np.vstack((state.phi0, state.phi1))
                    gram = cols.conj().T @ cols
                    eig = np.linalg.eigvalsh(np.sqrt(pair.weights)[:, None] * gram * np.sqrt(pair.weights)[None, :])
                    if eig.min() < -1e-10 or eig.max() > 1 + 1e-8:
                        raise ValueError("saved source lost CAR admissibility")
                    checked += 1
    return {"geometry_frames_checked": checked, "maximum_metric_gap": worst,
            "SVD_reconstructed_for_replay": False, "endpoint_energy_checked": True}


def run():
    if OUT.exists() or PAYLOAD.exists():
        raise FileExistsError("refusing to overwrite frozen replay-basis evidence")
    before, arrays, entries = bindings(), {}, {}
    for record in records():
        for name, expected in record.get("source_bindings", record.get("source_bindings_before", {})).items():
            if digest(LAB / name) != expected:
                raise ValueError("scientific source differs from its frozen producer: " + name)
    for nf in (128, 256, 512):
        pair = model.build_pair(nf)
        arrays[f"nf{nf}_W"] = pair.geometry_map
        arrays[f"nf{nf}_reference_columns"] = pair.reference_columns
        arrays[f"nf{nf}_source_columns"] = pair.source_columns
        entries[str(nf)] = {"coarse_indices": pair.geometry_coarse_indices.tolist(),
            "child_indices": pair.geometry_child_indices.tolist(),
            "parent_indices": pair.geometry_parent_indices.tolist(),
            "source": pair.source_metadata, "geometry": pair.geometry_metadata}
    report = validate(arrays, entries)
    if before != bindings():
        raise ValueError("frozen scientific inputs changed during basis capture")
    with PAYLOAD.open("xb") as stream:
        np.savez_compressed(stream, **arrays)
    record = {"schema": "NSC-NESTED-PARENT-CHILD-REPLAY-BASIS-v1",
        "source_bindings": before, "bases": entries, "validation": report,
        "payload_sha256": digest(PAYLOAD), "payload_bytes": PAYLOAD.stat().st_size,
        "trajectories_reexecuted": False, "original_payloads_changed": False}
    with OUT.open("x") as stream:
        json.dump(record, stream, indent=2, allow_nan=False)
        stream.write("\n")
    return record


def verify():
    record = json.loads(OUT.read_text())
    if digest(PAYLOAD) != record["payload_sha256"]:
        raise ValueError("frozen replay basis payload changed")
    for path, expected in record["source_bindings"].items():
        if digest(LAB / path) != expected:
            raise ValueError("replay basis source binding changed: " + path)
    with np.load(PAYLOAD, allow_pickle=False) as stored:
        report = validate(stored, record["bases"])
    return report


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    print(verify() if args.check else run())
