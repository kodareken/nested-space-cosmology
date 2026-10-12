#!/usr/bin/env python3
"""Creation-only remaining timestep control from the exact stopped NF32 arm.

This separate35-CPU batch preserves the original30-CPU record and forecast.
It resumes the full autonomous leading RHS; it never prepares a new radius,
selects a source, or adjusts the prediction after observing the held source.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
import json
from pathlib import Path
import resource
import time

import numpy as np
from threadpoolctl import threadpool_limits

from recursive_horizons import nsc_discovery_homogeneous_turn as owner
from recursive_horizons import nsc_discovery_turn_fullfield as full
from recursive_horizons import provenance

SCHEMA = "NSC-DISCOVERY-HOMOGENEOUS-TURN-CONFIRMATION-v1"
PARENT_COMMIT = "8ef2da55af68be2eec10502f45d06638c16aafa8"
DEFAULT_SOURCE = owner.OUTPUT
DEFAULT_OUTPUT = full.leading.LAB/"results/development/nsc-discovery-homogeneous-turn-confirmation-v1"
CPU_LIMIT = 35.
OWN_PATHS = ("lab/scripts/confirm_nsc_discovery_homogeneous_turn.py",
             "lab/tests/test_nsc_discovery_turn_confirmation.py")
# Observed original artifacts, pinned by this independent consumer after the
# original run. These are not represented as a pre-run artifact attestation.
ORIGINAL_PINS = {
    "prediction.json":"9b1ca8bbdc1d0276e55fa13a23df46f0033e94489e340494028f3935022e04f4",
    "prediction.npz":"2c70d4560df8a818cfa52e934e5658ee47824d19e2e0996c0c1b423034fac8c9",
    "measurement.json":"a666d2f9bc863fc9e5696da94cd5e233bad8421595713c69598e0a880879d603",
    "measurement.npz":"f7561aca49463c1219edeb5526070b6f5d5fc963ec9aae5ff150bea5f59cdafd"}


def cpu_usage():
    children=resource.getrusage(resource.RUSAGE_CHILDREN)
    return time.process_time()+children.ru_utime+children.ru_stime


def authenticate_source(source,*,require_current=True):
    source = Path(source).resolve()
    for name,expected in ORIGINAL_PINS.items():
        if full.digest(source/name) != expected:
            raise ValueError("original observed artifact changed: "+name)
    locked, _forecast_arrays = owner._read(source,"prediction")
    measured, arrays = owner._read(source,"measurement")
    for record in (locked, measured):
        if record.get("producing_commit") != PARENT_COMMIT:
            raise ValueError("confirmation requires the original8ef2 producing commit")
    if (locked.get("mode") != "locked_prediction"
            or locked.get("locked_before_held_data") is not True
            or measured.get("mode") != "independent_fullfield_measurement"):
        raise ValueError("confirmation requires the original LOCKED forecast and independent measurement")
    if measured["prediction_json_sha256"] != full.digest(source/"prediction.json"):
        raise ValueError("original measurement/forecast binding changed")
    if measured["backend"]["binding"]["prediction_sha256"] != full.digest(source/"prediction.json"):
        raise ValueError("backend forecast lock binding changed")
    # Historical scientific records/documents authenticate against8ef2. Code
    # actually imported by this continuation must also retain its pinned bytes.
    for name, expected in measured["producers"].items():
        if require_current and name.startswith("lab/src/") and full.digest(full.ROOT/name) != expected:
            raise ValueError("executing original producer changed: "+name)
    coarse, stopped = measured["backend"]["events"]
    if (coarse["status"] != "FIRST_POSITIVE_TO_NEGATIVE_P_Q_EVENT" or coarse["step_cap"] != .001
            or stopped["status"] != "CPU_BUDGET_STOP" or stopped["step_cap"] != .0005
            or stopped["nf"] != full.NF or not stopped["armed_after_positive_p_Q"]):
        raise ValueError("requires the one unfinished fine-cap control")
    pins = {str(source/(name+ext)): full.digest(source/(name+ext))
            for name in ("prediction","measurement") for ext in (".json",".npz")}
    pins.update(measured["input_hashes"])
    return locked, measured, arrays, pins


def reconstruct_handoff(locked, measured, arrays):
    """Rebuild carriers only; preserve every saved state/source/frame byte."""
    weights = np.array(arrays["cap1_source_weights"],copy=True)
    if not np.array_equal(weights,locked["preparation"]["held_pair_occupations"]):
        raise ValueError("saved source weights differ from locked source")
    columns = tuple(np.array(arrays["initial_"+name],copy=True) for name in ("phi0","phi1"))
    pair = full.leading.nested.build_pair(full.NF,coarse_modes=1,child_details=2,
        source_layout="override",columns_override=columns,occupations=weights)
    W = np.array(arrays["cap1_W"],copy=True)
    if not np.array_equal(W,pair.geometry_map) or not np.array_equal(W,arrays["cap0_W"]):
        raise ValueError("saved canonical geometry frame mismatch")
    pair = replace(full.leading.backend.make_fft_pair(pair),geometry_map=W)
    state = full.leading.State(*(np.array(arrays["cap1_"+name],copy=True) for name in full.leading.FIELDS))
    old = measured["backend"]["events"][1]
    if not np.array_equal(arrays["cap1_normal_clocks"],old["normal_clocks"]):
        raise ValueError("saved clocks differ from the stopped ledger")
    if not old["time"] > 0 or full.momentum(pair,state) <= 0:
        raise ValueError("saved fine arm is not a pre-turn positive-momentum state")
    fine = full.leading.check_chart(pair,state)
    if state.Q.shape != (pair.grid.ng,) or fine.phi0.shape != (pair.grid.nq,6):
        raise ValueError("saved full carrier/state shape mismatch")
    return pair,state


def preview(source=DEFAULT_SOURCE):
    locked,measured,arrays,pins = authenticate_source(source)
    with threadpool_limits(limits=1),full.leading.backend.fft_thread_limit(1):
        _pair,_state = reconstruct_handoff(locked,measured,arrays)
    return {"schema":SCHEMA,"mode":"preflight","evolved":False,"bytes_written":0,
        "case":"remaining_timestep_control","source_producing_commit":PARENT_COMMIT,
        "resume_time":measured["backend"]["events"][1]["time"],
        "resume_steps":measured["backend"]["events"][1]["steps"],"nf":full.NF,"step_cap":.0005,
        "fresh_CPU_limit_seconds":CPU_LIMIT,"original_CPU_limit_seconds":measured["CPU_limit_seconds"],
        "original_aggregate_CPU_seconds":measured["aggregate_CPU_seconds"],
        "locked_predicted_radius_change":locked["predicted_radius_change"],"input_hashes":pins,
        "source_radius_eigenframe_reprepared":False,"physical_scope":"finite CAR leading EFT only"}


def continue_arm(pair,state,old,*,deadline,maximum_time=None,before_advance):
    """Existing full-field helper uses elapsed time; restore the original ledger."""
    remaining = full.MAX_TIME-old["time"] if maximum_time is None else maximum_time
    record,arrays = full._advance(pair,state,.0005,deadline=deadline,
                                  maximum_time=remaining,before_advance=before_advance)
    elapsed = record["time"]
    record["resume_elapsed_time"] = elapsed
    record["time"] += old["time"]
    record["diagnostics"]["time"] = record["time"]
    record["steps_this_batch"] = record["steps"]
    record["steps"] += old["steps"]
    clocks = arrays["normal_clocks"]+np.asarray(old["normal_clocks"])
    record["normal_clocks"] = clocks.tolist()
    arrays["normal_clocks"] = clocks
    record["integrated_coordinate_fieldwork"] += old["integrated_coordinate_fieldwork"]
    record["initial_energy"] = dict(old["initial_energy"])
    energy = record["diagnostics"]["energy"]
    record["field_energy_change"] = energy["field"]-old["initial_energy"]["field"]
    record["fieldwork_balance_residual"] = record["field_energy_change"]-record["integrated_coordinate_fieldwork"]
    record["total_energy_drift"] = energy["total"]-old["initial_energy"]["total"]
    record["prior_partial_CPU_seconds"] = old["CPU_seconds"]
    record["combined_fine_arm_CPU_seconds"] = old["CPU_seconds"]+record["CPU_seconds"]
    record["exact_saved_state_continuation"] = True
    return record,arrays


def producer_pins(measured):
    names = [name for name in measured["producers"] if name.startswith("lab/src/")]
    return {name:full.digest(full.ROOT/name) for name in names+list(OWN_PATHS)}


def run(source=DEFAULT_SOURCE,output=DEFAULT_OUTPUT,*,cpu_limit_seconds=CPU_LIMIT,producing_commit=None):
    if not 0 < cpu_limit_seconds <= CPU_LIMIT:
        raise ValueError("remaining timestep-control budget must be in(0,35]")
    if producing_commit is None:
        raise ValueError("new confirmation requires its explicit frozen producing commit")
    output = Path(output).resolve()
    if (output/"confirmation.json").exists() or (output/"confirmation.npz").exists():
        raise FileExistsError("confirmation output is immutable")
    start = cpu_usage()
    locked,measured,saved,pins = authenticate_source(source)
    producers = producer_pins(measured)
    for name,expected in producers.items():
        provenance.resolve_pinned_source_bytes(full.ROOT,name,expected,commit=producing_commit)
    def guard():
        for name,expected in pins.items():
            if full.digest(name) != expected:raise ValueError("confirmation input changed: "+name)
        for name,expected in producers.items():
            if full.digest(full.ROOT/name) != expected:raise ValueError("continuation producer changed: "+name)
    with threadpool_limits(limits=1),full.leading.backend.fft_thread_limit(1):
        pair,state = reconstruct_handoff(locked,measured,saved)
        old = measured["backend"]["events"][1]
        remaining_cpu=cpu_limit_seconds-(cpu_usage()-start)
        event,arrays = continue_arm(pair,state,old,deadline=time.process_time()+max(0.,remaining_cpu-.25),
                                    before_advance=guard)
    arrays.update({"resume_"+name:np.array(saved["cap1_"+name],copy=True)
                   for name in full.leading.FIELDS+("W","source_weights","normal_clocks")})
    comparison = None
    if event["status"] == "FIRST_POSITIVE_TO_NEGATIVE_P_Q_EVENT":
        coarse = measured["backend"]["events"][0]
        effect = locked["predicted_radius_change"]
        change = event["radius"]-locked["baseline"]["radius"]
        movement = event["radius"]-coarse["radius"]
        comparison = {"fine_radius":event["radius"],"fine_time":event["time"],
            "coarse_radius":coarse["radius"],"coarse_time":coarse["time"],
            "radius_timestep_movement":movement,"time_timestep_movement":event["time"]-coarse["time"],
            "abs_radius_movement_over_locked_predicted_effect":abs(movement/effect) if effect else None,
            "abs_radius_movement_over_coarse_radius":abs(movement/coarse["radius"]),
            "locked_predicted_change":effect,"measured_fine_change":change,
            "unchanged_prediction_error":event["radius"]-locked["predicted_held_radius"],
            "prediction_error_over_measured_effect":abs((event["radius"]-locked["predicted_held_radius"])/change) if change else None}
    guard()
    cpu = cpu_usage()-start
    if cpu > cpu_limit_seconds:raise RuntimeError("fresh confirmation CPU cap exceeded")
    record = {"schema":SCHEMA,"mode":"remaining_timestep_control","producing_commit":producing_commit,
        "producers":producers,"input_hashes":pins,"original_producing_commit":PARENT_COMMIT,
        "original_measurement_modified":False,"forecast_recalibrated":False,"event":event,
        "comparison":comparison,"resume_time":old["time"],"resume_steps":old["steps"],
        "resume_array_sha256":{name:full.array_digest(saved["cap1_"+name]) for name in
                               full.leading.FIELDS+("W","source_weights","normal_clocks")},
        "CPU_seconds":cpu,"fresh_CPU_limit_seconds":cpu_limit_seconds,
        "CPU_accounting":"process plus historical-authentication child CPU; single-thread trajectory",
        "original_CPU_limit_seconds":measured["CPU_limit_seconds"],
        "original_aggregate_CPU_seconds":measured["aggregate_CPU_seconds"],
        "original_backend_CPU_seconds":measured["backend"]["CPU_seconds"],
        "budget_policy":"separate declared35-CPU remaining timestep-control batch; original30-CPU record unchanged",
        "source_radius_eigenframe_reprepared":False,"original_initial_energy":old["initial_energy"],
        "prior_integrated_fieldwork":old["integrated_coordinate_fieldwork"],"prior_normal_clocks":old["normal_clocks"],
        "status":"COMPLETE" if comparison else event["status"],"no_further_refinement_declared":True,
        "physical_scope":"finite CAR leading EFT; strong-curvature validity and vacuum matching open",
        "holding_regeneration_or_bounce_claimed":False}
    return owner._write(output,"confirmation",record,arrays)


def check(output=DEFAULT_OUTPUT):
    record,arrays = owner._read(output,"confirmation")
    if record.get("schema") != SCHEMA:raise ValueError("not a homogeneous-turn confirmation")
    original = Path(next(name for name in record["input_hashes"] if name.endswith("measurement.json"))).parent
    _locked,_measured,saved,_pins = authenticate_source(original,require_current=False)
    for name,expected in record["resume_array_sha256"].items():
        if full.array_digest(arrays["resume_"+name]) != expected or not np.array_equal(arrays["resume_"+name],saved["cap1_"+name]):
            raise ValueError("confirmation exact handoff drift: "+name)
    return {"schema":SCHEMA,"mode":"read_only_check","ok":True,"evolved":False,"bytes_written":0,
            "status":record["status"],"comparison":record["comparison"],"CPU_seconds":record["CPU_seconds"]}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    modes=parser.add_mutually_exclusive_group()
    modes.add_argument("--run",action="store_true");modes.add_argument("--check",action="store_true")
    parser.add_argument("--source",type=Path,default=DEFAULT_SOURCE)
    parser.add_argument("--output",type=Path,default=DEFAULT_OUTPUT)
    parser.add_argument("--cpu-limit-seconds",type=float,default=CPU_LIMIT)
    parser.add_argument("--producing-commit")
    args=parser.parse_args()
    result=check(args.output) if args.check else run(args.source,args.output,
        cpu_limit_seconds=args.cpu_limit_seconds,producing_commit=args.producing_commit) if args.run else preview(args.source)
    print(json.dumps(owner._json(result),sort_keys=True,indent=2,allow_nan=False))


if __name__ == "__main__":main()
