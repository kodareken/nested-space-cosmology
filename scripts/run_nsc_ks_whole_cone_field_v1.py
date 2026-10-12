#!/usr/bin/env python3
"""Production resumable whole-cone field campaign for representative family 14_1.

This driver binds the declared current history and the authenticated v4
profile payload, re-evolves family 14_1 from the unchanged retained
upstream source, and encloses streamed history-minus-reference cells.
It does not launch the expensive full scientific run unless
--run-family-14-1 is given explicitly. Diagnostic and bind modes cannot
emit a completed family record or a physical certificate.
"""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

import argparse
from dataclasses import replace
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from recursive_horizons.nsc_ks_current_field_campaign import (
    FREE_SPACE_FLOOR_BYTES, SCHEMA, STREAM_KIND_TRAJECTORY, V4_RELATIVE,
    bind_authenticated_v4_payload, bind_declared_history, bind_retained_source,
    production_field_config, production_stream, run_whole_cone_campaign,
    scientific_digest, scientific_record,
)
from recursive_horizons.nsc_ks_evaluation_binding import write_json_atomic
from recursive_horizons.nsc_ks_evaluation_binding import array_digest, implementation_hashes
from recursive_horizons.nsc_ks_source_envelope import usual_axial_support


DEFAULT_CHECKPOINT = ROOT / "results/development/nsc-ks-whole-cone-field-v1.checkpoint.json"
DEFAULT_RESULT = ROOT / "results/development/nsc-ks-whole-cone-field-v1.json"


def select_source_owner():
    import derive_nsc_ks_source_control_v2 as control
    return control.select_source


def bind_production(root=ROOT):
    payload, digest, v4 = bind_authenticated_v4_payload(root)
    family, identity, _history = bind_declared_history(root)
    selected = bind_retained_source(root, select_source_owner())
    config, grid = production_field_config(
        family, selected["source"], selected["batch"], digest,
        v4["settings"])
    target = family.collocation_nodes(config.settings["target_nodes"])
    config = bind_selected_inputs(config, selected, grid, target, root)
    return {
        "payload": payload,
        "payload_digest": digest,
        "v4_settings": dict(v4["settings"]),
        "family": family,
        "profile_identity": identity,
        "selected": selected,
        "config": config,
        "grid": grid,
        "target": target,
        "rho_up": float(selected["batch"].rho_up),
    }


def bind_selected_inputs(config, selected, grid, target, root):
    """Reject a changed preparation, covariance or implementation on resume."""
    import platform
    from importlib.metadata import version
    settings = dict(config.settings)
    settings["selected_input_binding"] = {
        "source": selected["source"].digest,
        "negative_source": selected["negative"].digest,
        "initial_columns": array_digest(selected["initial_columns"]),
        "archive_hashes": dict(selected["archive_hashes"]),
        "selections": selected["selections"],
        "channel": selected["channel"],
        "grid": array_digest(grid), "targets": array_digest(target),
    }
    settings["implementation_hashes"] = implementation_hashes(root, owners=(
        "scripts/run_nsc_ks_whole_cone_field_v1.py",
        "scripts/derive_nsc_ks_source_control_v2.py",
        "src/recursive_horizons/nsc_ks_current_field_campaign.py",
    ))
    settings["runtime_versions"] = {
        "python": platform.python_version(),
        **{name: version(name) for name in ("numpy", "scipy", "python-flint")},
    }
    return replace(config, settings=settings, settings_digest="")


def bind_report(bound):
    config = bound["config"]
    selected = bound["selected"]
    return {
        "schema": SCHEMA,
        "status": (
            "OPEN: production whole-cone bindings; no family 14_1 scientific "
            "run was launched; physical rho=1 source error unresolved"),
        "family": list(config.family),
        "profile_identity": bound["profile_identity"],
        "profile_payload_sha256": bound["payload_digest"],
        "profile_payload_recomputed": False,
        "source_via_select_source": True,
        "source_preparation_digest": selected["source"].digest,
        "source_columns": len(config.source_weights),
        "selections": selected["selections"],
        "mass": config.mass,
        "angular": config.angular,
        "rho_up": bound["rho_up"],
        "solver": {
            "integrator": "dop853",
            "rtol": config.settings["rtol"],
            "atol": config.settings["atol"],
            "max_step": config.settings["max_step"],
            "step_control": config.settings["step_control"],
            "tangents": config.settings["tangents"],
        },
        "config_digest": config.digest(),
        "axial_support": list(usual_axial_support()),
        "family_14_1_scientific_run": False,
        "all_history_cells": False,
        "physical_rho1_source_error": None,
        "physical_EXISTENCE_certificate": False,
        "physical_NONEXISTENCE_certificate": False,
        "certificate_use": False,
    }


def run_bound_campaign(bound, *, checkpoint, resume, max_new_cells, cpu_budget,
                       free_space_floor):
    selected = bound["selected"]
    batch = selected["batch"]
    stream = production_stream(
        selected["source"], selected["initial_columns"], bound["family"],
        bound["grid"], bound["target"], batch.mass, batch.angular, batch.rho_up)
    return run_whole_cone_campaign(
        bound["config"], bound["payload"], bound["family"], stream,
        checkpoint_path=checkpoint, resume=resume,
        max_new_cells=max_new_cells, cpu_budget=cpu_budget,
        free_space_floor=free_space_floor, stream_kind=STREAM_KIND_TRAJECTORY,
        source_via_select_source=True, continuous_source=selected,
        rho_up=bound["rho_up"])


def check_result(path):
    recorded = json.loads(Path(path).read_text())
    if recorded.get("schema") != SCHEMA:
        raise ValueError("unexpected whole-cone campaign schema")
    if recorded.get("certificate_use") or recorded.get("physical_EXISTENCE_certificate"):
        raise ValueError("campaign record cannot issue a physical certificate")
    if recorded.get("physical_rho1_source_error") is not None:
        raise ValueError("physical rho=1 source error is unresolved")
    if not recorded.get("status", "").startswith("OPEN"):
        raise ValueError("campaign record must remain OPEN")
    if scientific_digest(recorded) != recorded.get("scientific_digest"):
        raise ValueError("scientific digest changed")
    if "runtime" in scientific_record(recorded) or "forecast" in scientific_record(recorded):
        raise ValueError("runtime observations leaked into the scientific digest")
    if recorded.get("family_14_1_scientific_run"):
        coverage = recorded.get("coverage") or {}
        if not coverage.get("all_history_cells") or not coverage.get("stream_reached_rho1"):
            raise ValueError("family_14_1_scientific_run requires full cell coverage")
        if coverage.get("diagnostic") or coverage.get("resource_stopped"):
            raise ValueError("partial output cannot be a completed family record")
    return recorded


def display(value):
    print(json.dumps({
        "schema": value.get("schema"),
        "status": value.get("status"),
        "family": value.get("family"),
        "profile_identity": value.get("profile_identity"),
        "coverage": value.get("coverage"),
        "family_14_1_scientific_run": value.get("family_14_1_scientific_run"),
        "all_history_cells": value.get("all_history_cells"),
        "physical_rho1_source_error": value.get("physical_rho1_source_error"),
        "stop_reason": value.get("stop_reason"),
        "config_digest": value.get("config_digest"),
        "scientific_digest": value.get("scientific_digest"),
        "runtime": value.get("runtime"),
        "forecast": value.get("forecast"),
    }, indent=2, sort_keys=True))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--bind", action="store_true",
                       help="bind source, history and v4 payload; do not evolve")
    modes.add_argument("--diagnostic-max-new-cells", type=int, metavar="N",
                       help="enclose at most N new cells and keep a checkpoint")
    modes.add_argument("--run-family-14-1", action="store_true",
                       help="expensive full family 14_1 scientific run")
    modes.add_argument("--check", metavar="RESULT",
                       help="validate a campaign result record")
    parser.add_argument("--checkpoint", type=Path, default=DEFAULT_CHECKPOINT)
    parser.add_argument("--result", type=Path, default=DEFAULT_RESULT)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--cpu-budget", type=float, default=None)
    parser.add_argument("--write-result", action="store_true")
    args = parser.parse_args()
    if args.check:
        display(check_result(args.check))
        return
    if args.bind:
        report = bind_report(bind_production(ROOT))
        report["scientific_digest"] = scientific_digest(report)
        display(report)
        return
    if args.run_family_14_1 and args.diagnostic_max_new_cells:
        raise ValueError("diagnostic mode cannot launch the family scientific run")
    max_new_cells = args.diagnostic_max_new_cells
    if args.run_family_14_1:
        max_new_cells = None
    bound = bind_production(ROOT)
    mapping, _result = run_bound_campaign(
        bound, checkpoint=args.checkpoint, resume=args.resume,
        max_new_cells=max_new_cells, cpu_budget=args.cpu_budget,
        free_space_floor=FREE_SPACE_FLOOR_BYTES)
    if args.write_result:
        write_json_atomic(args.result, mapping, indent=2)
    display(mapping)


if __name__ == "__main__":
    main()
