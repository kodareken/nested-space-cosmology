#!/usr/bin/env python3
"""Read saved discovery stations and report the observer tidal field.

``--write PATH`` creates a summary once and
does not overwrite an existing file. ``--check PATH`` re-hashes every bound
input and producer. No station is evolved and no source file is written.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

from recursive_horizons.nsc_discovery_tidal import (
    measurement_report,
    jsonable,
    sha256_file,
)


def digest(report):
    """Short evidence block. The full station tables stay in the record."""
    comparisons = report["comparisons"]
    growth = comparisons["growth"]["worldline_x2"]["channels"]
    peak = comparisons["growth"]["peak_R0202"]["channels"]
    space = comparisons["space_nf256_versus_nf128_dt0.0005"]["3.0"]
    step = comparisons["timestep_nf256_dt0.001_versus_dt0.0005"]["3.0"]
    finer = comparisons.get("space_nf512_versus_nf256_dt0.0005")
    finer_tide = None if not finer else finer["3.0"]["R_0202"]["max_abs"]["relative_difference"]
    finer_w = None if not finer else finer["3.0"]["owned_W"]["max_abs"]["relative_difference"]
    return {
        "schema": report["schema"],
        "cpu_seconds": report["cpu_seconds"],
        "saved_state_hashes_unchanged": report["saved_state_hashes_unchanged"],
        "chi_substituted_for_curvature": report["chi_substituted_for_curvature"],
        "new_ode": report["new_ode"],
        "singularity_declared": report["singularity_declared"],
        "theory_declared_dead": report["theory_declared_dead"],
        "geometry_evolved": report["geometry_evolved"],
        "production_record_written": report["production_record_written"],
        "flat_minkowski_curvature_max_abs": report["exact"]["minkowski"]["flat_curvature_max_abs"],
        "flat_milne_curvature_max_abs": report["exact"]["milne"]["flat_curvature_max_abs"],
        "bertotti_robinson_R4": report["exact"]["bertotti_robinson"]["R4"],
        "bertotti_robinson_Ricci2": report["exact"]["bertotti_robinson"]["Ricci2"],
        "bertotti_robinson_K": report["exact"]["bertotti_robinson"]["K"],
        "worldline_R0202": growth["R_0202"]["values"],
        "worldline_R0101": growth["R_0101"]["values"],
        "worldline_K": growth["K"]["values"],
        "worldline_proper_time": comparisons["growth"]["worldline_x2"]["proper_time_x2"],
        "peak_R0202": peak["R_0202"]["values"],
        "peak_x": comparisons["growth"]["peak_R0202"]["x"],
        "clock_x2": comparisons["growth"]["clock_comparison"]["worldline_x2"],
        "space_T3_R0202_relative": space["R_0202"]["max_abs"]["relative_difference"],
        "space_T3_W_relative": space["owned_W"]["max_abs"]["relative_difference"],
        "space_T3_R4_relative": space["R4"]["max_abs"]["relative_difference"],
        "timestep_T3_R0202_relative": step["R_0202"]["relative_difference"],
        "nf512_status": report["confirmation_nf512"]["status"],
        "nf512_T3_R0202_relative": finer_tide,
        "nf512_T3_W_relative": finer_w,
        "frozen_geometry_jets_zero": comparisons["frozen"]["geometry_jets_zero"],
        "frozen_field_still_flows": comparisons["frozen"]["field_still_flows"],
        "finite_window": report["finite_window"]["statement"],
    }


def _bound_hashes(record):
    hashes = {}
    hashes.update(record.get("producer_hashes") or {})
    inputs = record.get("input_hashes") or {}
    hashes.update(inputs.get("before") or {})
    confirmation = record.get("confirmation_nf512") or {}
    hashes.update(confirmation.get("input_hashes_before") or {})
    for binding in (record.get("frozen_producer_bindings") or {}).values():
        hashes[binding["envelope"]] = binding["envelope_sha256"]
        hashes.update(binding["physics_hashes"])
    return hashes


def check(path):
    record = json.loads(Path(path).read_text(encoding="utf-8"))
    problems = []
    if record.get("schema") != "NSC-DISCOVERY-TIDAL-v1":
        problems.append("unexpected schema")
    if record.get("singularity_declared") is not False:
        problems.append("record declares a singularity")
    if record.get("theory_declared_dead") is not False:
        problems.append("record declares the theory dead")
    if record.get("chi_substituted_for_curvature") is not False:
        problems.append("record substitutes chi for curvature")
    if record.get("new_ode") is not False or record.get("geometry_evolved") is not False:
        problems.append("record claims a new evolution")
    if record.get("saved_state_hashes_unchanged") is not True:
        problems.append("record does not preserve saved-state hashes")
    if not record.get("producer_hashes") or not record.get("input_hashes", {}).get("before"):
        problems.append("record lacks source or input bindings")
    bindings = record.get("frozen_producer_bindings") or {}
    if not bindings:
        problems.append("record lacks frozen producer bindings")
    for name, digest_value in _bound_hashes(record).items():
        file_path = Path(name)
        if not file_path.is_file():
            problems.append("missing " + name)
            continue
        if sha256_file(file_path) != digest_value:
            problems.append("hash changed " + name)
    if problems:
        raise SystemExit("tidal check failed: " + "; ".join(problems))
    return {"ok": True, "checked_hashes": len(_bound_hashes(record)), "evolved": False, "bytes_written": 0}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--write", metavar="PATH", help="create a record once; never overwrite")
    modes.add_argument("--check", metavar="PATH", help="re-hash a previously written record")
    parser.add_argument("--nf512", action=argparse.BooleanOptionalAction, default=True,
                        help="include saved nf512 confirmation (default: selected)")
    parser.add_argument("--episode-dir", type=Path, default=None)
    parser.add_argument("--confirmation-dir", type=Path, default=None)
    args = parser.parse_args(argv)
    if args.check:
        print(json.dumps(check(args.check), indent=2, sort_keys=True))
        return 0
    destination = Path(args.write).expanduser().resolve()
    if destination.exists():
        raise SystemExit("refusing to overwrite " + str(destination))
    report = measurement_report(episode_dir=args.episode_dir, confirmation_dir=args.confirmation_dir,
                                include_nf512=args.nf512, extra_paths=(Path(__file__).resolve(),))
    payload = jsonable(report)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("x", encoding="utf-8") as stream:
        stream.write(json.dumps(payload, indent=2, allow_nan=False) + "\n")
    check(destination)
    print(json.dumps(digest(payload), indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
