#!/usr/bin/env python3
"""Create/check fixed-input64-cell SSPRK3 windows using the retained cache."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import time
import uuid

import numpy as np
from validation_source_successors import verify_saved_bindings, verify_saved_source

import record_nsc_atlas_ssprk3_refinement as previous
from analyze_nsc_atlas_spatial_refinement import diagnostics

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "nsc-atlas-spatial-refinement-v1"
RELATIVE_OUTPUT = "results/development/nsc-atlas-spatial-refinement-v1"
STATUS = "FIXED_INPUT_SPATIAL_COMPARISON"
SCOPE = ("64 versus 32 exterior/cap cells with the same 512-step SSPRK3 method, source snapshot, "
         "computational profile and cache; sampled spatial sensitivity, not a convergence order or error enclosure")
INPUT = previous.RELATIVE_OUTPUT+"/refinement.json"
INPUT_SHA256 = "ab95bdcd43b6df78e776f5282387e3f7398d8996a25d4e3de4e5a8226fc25e02"
MANIFEST_SHA256 = "f85a1491eeba5c77514936e48f69a4335d5a92edda2a0db2a0dc624a3701f041"
BASE_INPUT = "results/development/nsc-atlas-coupled-evolution-v1/coupled.json"
CASES = (("zero",64,512),("finite",64,512))
SETTINGS = dict(input=INPUT,cases=[dict(amplitude=a,cells=n,steps=k) for a,n,k in CASES],
    source_cache_nodes=1025,source_Hermite_knots=129,profile_training_rows=65,profile_withheld_rows=64,
    saved_time_nodes=129,exterior_quadrature_order=2,cap_quadrature_order=4,maximum_courant=.4,
    tube=previous.previous.TUBE,time_scheme=previous.previous.TIME_SCHEME,
    global_cpu_seconds=450.,global_evaluations=80000000,standalone_cpu_seconds=60.,standalone_evaluations=8000000,
    source_continuation=previous.SETTINGS["source_continuation"],clock_origin=previous.SETTINGS["clock_origin"],
    cache="reference the pinned parent's1025 exact rows; no reconstruction or quadrature",
    counts="measured shared-budget deltas; unsaved stage point totals are not independently reconstructed",
    cap_comparison="sum adjacent total physical dr integrals at common normalized radial cells; radii may differ")
SOURCES = tuple(dict.fromkeys((*previous.SOURCES,INPUT,previous.RELATIVE_OUTPUT+"/manifest.json",
    "scripts/analyze_nsc_atlas_spatial_refinement.py","tests/test_analyze_nsc_atlas_spatial_refinement.py",
    "scripts/record_nsc_atlas_spatial_refinement.py","tests/test_record_nsc_atlas_spatial_refinement.py")))
Budget = previous.Budget
plain,encoded,digest = previous.plain,previous.encoded,previous.digest
runtime_versions = previous.runtime_versions
_equal,_array,_count,_real = previous._equal,previous._array,previous._count,previous._real
subset,StoredRows = previous.subset,previous.StoredRows
physical_view,final_summary = previous.physical_view,previous.final_summary
fit_profile,background_stage = previous.fit_profile,previous.background_stage
evolve_background_ssprk3,quiet_departures = previous.evolve_background_ssprk3,previous.quiet_departures
_failure_safe = previous._failure_safe


def shared_budget(budget=None):
    if budget is None:
        budget = Budget(maximum_evaluations=80000000,cpu_seconds=450.)
    if (not isinstance(budget,Budget) or type(budget.maximum_evaluations) is not int
            or not 0 < budget.maximum_evaluations <= 80000000
            or isinstance(budget.cpu_seconds,bool) or not math.isfinite(budget.cpu_seconds)
            or not 0 < budget.cpu_seconds <= 450. or type(budget.evaluations) is not int or budget.evaluations < 0
            or isinstance(budget.started,bool) or not math.isfinite(budget.started)
            or not 0 <= budget.started <= time.process_time()):
        raise ValueError("invalid shared spatial-refinement allowance")
    budget.check()
    return budget


def source_bindings():
    return {name:digest((ROOT/name).read_bytes()) for name in SOURCES}


def load_inputs(budget):
    """Authenticate pinned bytes and86 bindings without recursive validation."""
    budget.check()
    directory = ROOT/previous.RELATIVE_OUTPUT
    if ({p.name for p in directory.iterdir()} != {"refinement.json","manifest.json"}
            or any(p.is_symlink() for p in directory.iterdir())):
        raise ValueError("inherited refinement inventory mismatch")
    raw = (directory/"refinement.json").read_bytes()
    mr = (directory/"manifest.json").read_bytes()
    manifest = json.loads(mr)
    bindings = manifest.get("source_bindings",{})
    if (digest(raw) != INPUT_SHA256 or digest(mr) != MANIFEST_SHA256
            or manifest.get("schema") != previous.SCHEMA+"-manifest"
            or manifest.get("outputs") != {"refinement.json":dict(sha256=digest(raw),bytes=len(raw))}
            or set(bindings) != set(previous.SOURCES) or len(bindings) != 86):
        raise ValueError("pinned refinement payload/manifest mismatch")
    base_raw = None
    for name,expected in bindings.items():
        path = (ROOT/name).resolve()
        if not path.is_relative_to(ROOT.resolve()):
            raise ValueError("bound input exits the repository")
        data = verify_saved_source(name, expected, ROOT)
        if name == BASE_INPUT:
            base_raw = data
    parent = json.loads(raw)
    if (parent.get("schema") != previous.SCHEMA or parent.get("settings") != previous.SETTINGS
            or parent.get("runtime") != runtime_versions() or base_raw is None):
        raise ValueError("inherited declaration/runtime mismatch")
    base = json.loads(base_raw)
    _equal(parent["geometry"],base["geometry"],"inherited geometry")
    if [(c["amplitude"],c["cells"],c["steps"]) for c in parent["cases"]] != [("zero",32,512),("finite",32,512)]:
        raise ValueError("require the two retained32-cell512-step cases")
    budget.check()
    identity = dict(path=INPUT,sha256=digest(raw),bytes=len(raw),manifest_sha256=digest(mr),
        authenticated_input_files=86,base_path=BASE_INPUT,base_sha256=bindings[BASE_INPUT],
        source_integrations=0,cache_reconstructions=0,trajectory_replays=0,scientific_validator_calls=0)
    return parent,base,identity


def cache_identity(parent):
    raw = encoded(parent["cache"])
    return dict(parent_payload_sha256=INPUT_SHA256,cache_sha256=digest(raw),cache_bytes=len(raw),
                nodes=1025,representation="exact parent JSON rows; no new interpolation")


def calculate(budget=None):
    budget = shared_budget(budget)
    initial = budget.evaluations
    phase,record = dict(name="hash-only input authentication"),{}
    try:
        parent,base,identity = load_inputs(budget)
        auth_charges = budget.evaluations-initial
        phase = dict(name="unchanged computational profile")
        before = budget.evaluations
        profile,summary = fit_profile(parent["profile_inputs"],budget)
        profile_charges = budget.evaluations-before
        _equal(summary,parent["profile"],"unchanged profile")
        record = dict(schema=SCHEMA,status=STATUS,scope=SCOPE,settings=SETTINGS,runtime=runtime_versions(),
            parent_input=identity,cache_identity=cache_identity(parent),geometry=parent["geometry"],
            profile_inputs=parent["profile_inputs"],profile=summary,cases=[],uniform_error_enclosure=False)
        cache = parent["cache"]
        x = _array(cache["x_nodes"],(1025,),"retained cache coordinates")
        for label,n,steps in CASES:
            phase = dict(name="new fixed-input spatial evolution",amplitude=label,cells=n,steps=steps)
            source,gauge = StoredRows(cache["source_rows"][label]),StoredRows(cache["gauge_rows"])
            before = budget.evaluations
            solution = evolve_background_ssprk3(quiet_departures(exterior_cells=n,cap_cells=n),x[::2],
                mid_nodes=x[1::2],eta=np.linspace(0,1,n+1),y_faces=np.linspace(0,1,n+1),
                source=source,gauge=gauge,alpha=base["initial_data"][label]["alpha"],tube=SETTINGS["tube"],
                budget=budget,profile=profile,save_every=4,maximum_courant=.4,
                exterior_quadrature_order=2,cap_quadrature_order=4)
            if solution["steps"] != steps or solution["time_scheme"] != SETTINGS["time_scheme"]:
                raise ValueError("spatial run returned a different time scheme")
            final = final_summary(solution["final_stage"])
            record["cases"].append(plain(dict(amplitude=label,cells=n,steps=steps,
                alpha=base["initial_data"][label]["alpha"],snapshots=solution["snapshots"],final=final,
                execution=dict(budget_start=before,budget_end=budget.evaluations,measured_evaluations=budget.evaluations-before,
                    stages=solution["stages"],source_cache_reads=source.calls,gauge_cache_reads=gauge.calls,
                    max_courant=solution["max_courant"],retries=solution["retries"],
                    final_stage_evaluations=8*(n+1)+2+final["cap_metric_evaluations"]+final["cap_rhs_evaluations"]))))
        phase = dict(name="saved spatial diagnostics")
        before = budget.evaluations
        record["diagnostics"] = diagnostics(record,parent,base,profile,budget)
        record["execution"] = dict(initial_evaluations=initial,evaluations=budget.evaluations,
            cpu_seconds=time.process_time()-budget.started,inherited_authentication_evaluations=auth_charges,
            profile_build_evaluations=profile_charges,diagnostic_evaluations=budget.evaluations-before,
            coupled_evolutions=len(CASES),time_steps=sum(k for _,_,k in CASES),stage_evaluations=sum(3*k+1 for _,_,k in CASES),
            source_integrations=0,cache_reconstructions=0,gauge_quadratures=0,old_trajectory_replays=0,
            profile_stage_extractions=0,stage_point_totals_are_measured=True)
        budget.check()
        return plain(record)
    except Exception as error:
        error.spatial_failure = plain(dict(phase=phase,exception=type(error).__name__,message=str(error),
            stage=getattr(error,"diagnostic",None),initial_evaluations=initial,evaluations=budget.evaluations,
            cpu_seconds=time.process_time()-budget.started,completed_cases=record.get("cases",[]),
            parent_input=record.get("parent_input")))
        raise


def validate_record(record,budget=None,*,require_creation=True):
    budget = shared_budget(budget)
    start = budget.evaluations
    try:
        if (record["schema"] != SCHEMA or record["status"] != STATUS or record["scope"] != SCOPE
                or record["uniform_error_enclosure"] is not False):
            raise ValueError("spatial declaration mismatch")
        _equal(record["settings"],SETTINGS,"settings")
        _equal(record["runtime"],runtime_versions(),"runtime")
        parent,base,identity = load_inputs(budget)
        auth_charges = budget.evaluations-start
        _equal(record["parent_input"],identity,"parent identity")
        _equal(record["cache_identity"],cache_identity(parent),"unchanged cache identity")
        _equal(record["geometry"],parent["geometry"],"geometry")
        _equal(record["profile_inputs"],parent["profile_inputs"],"profile inputs")
        before = budget.evaluations
        profile,summary = fit_profile(record["profile_inputs"],budget)
        profile_charges = budget.evaluations-before
        _equal(record["profile"],summary,"profile reconstruction")
        _equal(summary,parent["profile"],"unchanged profile")
        if len(record["cases"]) != len(CASES):
            raise ValueError("both ordered spatial cases required")
        cache = parent["cache"]
        x = _array(cache["x_nodes"],(1025,),"retained cache coordinates")
        initial = record["execution"]["initial_evaluations"]
        _count(initial,80000000,"initial shared charges")
        end = initial+auth_charges+profile_charges
        for case,(label,n,steps) in zip(record["cases"],CASES):
            alpha = base["initial_data"][label]["alpha"]
            if (case["amplitude"],case["cells"],case["steps"],case["alpha"]) != (label,n,steps,alpha):
                raise ValueError("case identity mismatch")
            if len(case["snapshots"]) != 129:
                raise ValueError("require129 saved cuts")
            eta = np.linspace(0,1,n+1)
            for j,(sample,at) in enumerate(zip(case["snapshots"],x[::8])):
                if sample["x"] != at or set(sample) != {"x","state","physical"}:
                    raise ValueError("saved row/time inventory mismatch")
                state = sample["state"]
                if set(state) != {"departure","wall","cap"}:
                    raise ValueError("raw state inventory mismatch")
                for group,names,shape in (("departure",("r","L","c"),(n-1,)),
                        ("wall",("r","L","T","s"),()),("cap",("Ip","Iq"),(n,))):
                    if set(state[group]) != set(names):
                        raise ValueError("raw state fields mismatch")
                    for key in names:
                        _array(state[group][key],shape,group+"."+key)
                expected = physical_view(state,at,subset(cache["source_rows"][label],8*j),
                    subset(cache["gauge_rows"],8*j),eta,alpha,profile,budget)
                _equal(sample["physical"],expected,"raw-to-physical reconstruction")
            _equal(case["snapshots"][0]["state"],quiet_departures(exterior_cells=n,cap_cells=n),"quiet departures")
            for key in ("T","s"):
                clocks = np.array([s["physical"]["wall"][key] for s in case["snapshots"]])
                if clocks[0] != 0 or np.any(np.diff(clocks) <= 0):
                    raise ValueError("clocks must increase from zero")
            before = budget.evaluations
            stage = background_stage(case["snapshots"][-1]["state"],x=x[-1],eta=eta,y_faces=eta,
                source=subset(cache["source_rows"][label],-1),gauge=subset(cache["gauge_rows"],-1),
                alpha=alpha,tube=SETTINGS["tube"],budget=budget,profile=profile,
                exterior_quadrature_order=2,cap_quadrature_order=4)
            final_charge = budget.evaluations-before
            _equal(case["final"],final_summary(stage),"final stage")
            execution = case["execution"]
            charge = execution["measured_evaluations"]
            _count(charge,80000000,"measured case points")
            if charge < (3*steps+1)*(8*(n+1)+2):
                raise ValueError("case charge below exterior/hybrid minimum")
            expected = dict(budget_start=end,budget_end=end+charge,measured_evaluations=charge,
                stages=3*steps+1,source_cache_reads=3*steps+1,gauge_cache_reads=3*steps+1,retries=0,
                final_stage_evaluations=final_charge)
            _equal({k:v for k,v in execution.items() if k != "max_courant"},expected,"case accounting")
            _real(execution["max_courant"],"Courant",high=.4)
            end += charge
        before = budget.evaluations
        _equal(record["diagnostics"],diagnostics(record,parent,base,profile,budget),"diagnostics")
        diagnostic_charges = budget.evaluations-before
        expected = dict(initial_evaluations=initial,evaluations=end+diagnostic_charges,
            inherited_authentication_evaluations=auth_charges,profile_build_evaluations=profile_charges,
            diagnostic_evaluations=diagnostic_charges,coupled_evolutions=len(CASES),time_steps=sum(k for _,_,k in CASES),
            stage_evaluations=sum(3*k+1 for _,_,k in CASES),source_integrations=0,cache_reconstructions=0,
            gauge_quadratures=0,old_trajectory_replays=0,profile_stage_extractions=0,stage_point_totals_are_measured=True)
        _equal({k:v for k,v in record["execution"].items() if k != "cpu_seconds"},expected,"execution accounting")
        _count(record["execution"]["evaluations"],80000000,"calculation points")
        _real(record["execution"]["cpu_seconds"],"calculation CPU",high=450.)
        budget.check()
        charges = budget.evaluations-start
        if require_creation:
            creation = record["creation_execution"]
            _equal({k:v for k,v in creation.items() if k != "cpu_seconds"},dict(
                calculation_evaluations=record["execution"]["evaluations"],validation_evaluations=charges,
                total_evaluations=record["execution"]["evaluations"]+charges),"creation accounting")
            _count(creation["total_evaluations"],80000000,"creation points")
            _real(creation["cpu_seconds"],"creation CPU",low=record["execution"]["cpu_seconds"],high=450.)
    except (KeyError,TypeError,IndexError,OverflowError) as error:
        raise ValueError("incomplete spatial-refinement record") from error
    return dict(authentication_evaluations=budget.evaluations-start,source_integrations=0,time_evolutions=0,
        cache_reconstructions=0,gauge_quadratures=0,current_final_radial_stages=len(CASES),
        unsaved_stage_point_totals_independently_recomputed=False)


def create(directory):
    started = time.process_time()
    budget = Budget(maximum_evaluations=80000000,cpu_seconds=450.,started=started)
    directory = Path(directory)
    if directory.exists():
        raise FileExistsError("creation-only record: choose a new directory or --check")
    bindings,runtime = source_bindings(),runtime_versions()
    try:
        record = calculate(budget)
        if bindings != source_bindings() or runtime != runtime_versions() or record["runtime"] != runtime:
            raise RuntimeError("spatial sources/runtime changed during calculation")
        auth = validate_record(record,budget,require_creation=False)
        if bindings != source_bindings() or runtime != runtime_versions():
            raise RuntimeError("spatial sources/runtime changed during authentication")
        record["creation_execution"] = dict(calculation_evaluations=record["execution"]["evaluations"],
            validation_evaluations=auth["authentication_evaluations"],total_evaluations=budget.evaluations,
            cpu_seconds=time.process_time()-started)
        raw = encoded(record)
        manifest = dict(schema=SCHEMA+"-manifest",source_bindings=bindings,
            outputs={"spatial.json":dict(sha256=digest(raw),bytes=len(raw))})
        budget.check()
    except Exception as error:
        receipt = dict(schema=SCHEMA+"-failed-attempt",status="FAILED_NOT_COMPLETED",settings=SETTINGS,
            source_bindings_before=bindings,runtime=runtime,
            failure=getattr(error,"spatial_failure",dict(exception=type(error).__name__,message=str(error),
                evaluations=budget.evaluations,cpu_seconds=time.process_time()-started)))
        path = ROOT/".build/spatial-refinement-failures"/(uuid.uuid4().hex+".json")
        path.parent.mkdir(parents=True,exist_ok=True)
        with path.open("xb") as stream:
            stream.write(encoded(_failure_safe(receipt)))
        error.failure_receipt = str(path)
        raise
    directory.mkdir(parents=True,exist_ok=False)
    (directory/"spatial.json").write_bytes(raw)
    (directory/"manifest.json").write_bytes(encoded(manifest))
    return dict(created=str(directory),execution=record["execution"],creation_execution=record["creation_execution"])


def check(directory):
    budget = Budget(maximum_evaluations=8000000,cpu_seconds=60.)
    directory = Path(directory)
    if ({p.name for p in directory.iterdir()} != {"spatial.json","manifest.json"}
            or any(p.is_symlink() for p in directory.iterdir())):
        raise ValueError("spatial output inventory mismatch")
    raw = (directory/"spatial.json").read_bytes()
    manifest = json.loads((directory/"manifest.json").read_bytes())
    bindings,runtime = source_bindings(),runtime_versions()
    if (manifest.get("schema") != SCHEMA+"-manifest" or not verify_saved_bindings(manifest.get("source_bindings"), ROOT, expected_paths=SOURCES)
            or manifest.get("outputs") != {"spatial.json":dict(sha256=digest(raw),bytes=len(raw))}):
        raise ValueError("spatial hash/source mismatch")
    result = validate_record(json.loads(raw),budget)
    if bindings != source_bindings() or runtime != runtime_versions():
        raise RuntimeError("spatial sources/runtime changed during authentication")
    budget.check()
    return dict(checked=str(directory),**result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--check",action="store_true")
    args = parser.parse_args()
    try:
        result = check(args.output) if args.check else create(args.output)
    except Exception as error:
        if hasattr(error,"failure_receipt"):
            print(json.dumps(dict(failed=str(error),diagnostic=error.failure_receipt)))
        raise
    print(json.dumps(result,sort_keys=True,allow_nan=False))


if __name__ == "__main__":
    main()
