#!/usr/bin/env python3
"""Source-bound terminal cap mass-functional derivative; no trajectory replay."""
from __future__ import annotations

import argparse
import json
import math
from numbers import Real
from pathlib import Path
import time
import uuid

import numpy as np

import diagnose_nsc_atlas_cap_consistency as previous
from analyze_nsc_atlas_cap_mass_direction import analyze
from validation_source_successors import verify_saved_bindings

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "nsc-atlas-cap-mass-direction-audit-v1"
RELATIVE_OUTPUT = "results/development/nsc-atlas-cap-mass-direction-v1.json"
CAP_INPUT = previous.RELATIVE_OUTPUT
CAP_SHA256 = "79d1e02917ba7c1a7d3326cd80612b0467f19b0b5b7b1e69b6e1e2f1f3ae9728"
LIMITS = dict(creation=dict(cpu_seconds=60., evaluated_points=4_000_000),
              check=dict(cpu_seconds=30., evaluated_points=2_000_000))
SCOPE = ("terminal cap mass-functional derivative along fixed implemented transport rates; "
         "backward step/order/branch sensitivities, not a trajectory or error enclosure")
OWN_SOURCES = (*previous.OWN_SOURCES,
    "scripts/analyze_nsc_atlas_cap_mass_direction.py",
    "tests/test_analyze_nsc_atlas_cap_mass_direction.py",
    "tests/test_nsc_atlas_cap_mass_direction_orchestration.py",
    "scripts/diagnose_nsc_atlas_cap_mass_direction.py",
    "tests/test_diagnose_nsc_atlas_cap_mass_direction.py")
Budget, plain, encoded, digest = previous.Budget, previous.plain, previous.encoded, previous.digest
runtime_versions = previous.runtime_versions
same_recomputed = previous.same_recomputed
EXPECTED_CALLS = dict(metric_calls=20, rhs_calls=4, profile_constructions=1,
    profile_fits=0, stage_calls=0, prefix_calls=0, hybrid_calls=0,
    source_solves=0, gauge_quadratures=0, evolution_calls=0)


def shared_budget(value=None, *, checking=False):
    limit = LIMITS["check" if checking else "creation"]
    if value is None:
        value = Budget(maximum_evaluations=limit["evaluated_points"], cpu_seconds=limit["cpu_seconds"])
    if (not isinstance(value,Budget) or type(value.maximum_evaluations) is not int
            or not 0 < value.maximum_evaluations <= limit["evaluated_points"]
            or type(value.evaluations) is not int or value.evaluations < 0):
        raise ValueError("invalid shared cap-mass-direction allowance")
    for name in ("started","cpu_seconds"):
        v = getattr(value,name)
        if (isinstance(v,(bool,np.bool_)) or not isinstance(v,Real) or not math.isfinite(v)
                or v < 0 or (name=="cpu_seconds" and not 0 < v <= limit["cpu_seconds"])):
            raise ValueError("invalid shared cap-mass-direction clock")
    if value.started > time.process_time():
        raise ValueError("shared clock cannot start in the future")
    value.check()
    return value


def _own_bindings(root):
    return {name:digest((Path(root)/name).read_bytes()) for name in OWN_SOURCES}


def authenticate(root,budget,bindings):
    """One pinned cap payload and one existing 104-file hash pass; no analysis.

    The prior audit's seven producer bindings reuse the combined 12-file
    inventory already read by the caller, rather than hashing them twice.
    """
    budget.check()
    path = Path(root)/CAP_INPUT
    if path.is_symlink():
        raise ValueError("pinned cap audit must be a regular retained file")
    raw = path.read_bytes()
    cap_digest = digest(raw)
    if cap_digest != CAP_SHA256:
        raise ValueError("pinned cap audit mismatch")
    cap = json.loads(raw)
    current,_,cache,identity = previous.authenticate(root,budget)
    expected = cap.get("audit_source_bindings")
    if (cap.get("schema") != previous.SCHEMA or cap.get("scope") != previous.SCOPE
            or cap.get("limits") != previous.LIMITS or cap.get("runtime") != runtime_versions()
            or cap.get("input") != identity):
        raise ValueError("cap audit declaration/input/source mismatch")
    verify_saved_bindings(expected,root,expected_paths=previous.OWN_SOURCES)
    budget.check()
    return current,cap,cache,dict(cap_path=CAP_INPUT,cap_sha256=cap_digest,cap_bytes=len(raw),
        cap_audit_source_bindings=expected,inherited=identity,cap_scientific_replays=0)


def _analysis(current,cap_audit,cache,budget):
    start = budget.evaluations
    result = plain(analyze(current,cap_audit,cache,budget))
    cost = budget.evaluations-start
    if type(result.get("evaluated_points")) is not int or result["evaluated_points"] != cost:
        raise ValueError("analysis point receipt disagrees with shared budget")
    execution = result.get("execution",{})
    if any(type(execution.get(k)) is not int or execution[k] != v for k,v in EXPECTED_CALLS.items()):
        raise ValueError("analysis call inventory mismatch")
    charges = execution.get("charges",{})
    if (set(charges) != {"input_and_profile","metrics","rhs","probe_analysis"}
            or any(type(v) is not int or v < 0 for v in charges.values()) or sum(charges.values()) != cost):
        raise ValueError("analysis charge inventory mismatch")
    encoded(result)  # rejects nonfinite values before any output is created
    budget.check()
    return result


def _execution(cost,initial,cpu):
    return dict(analysis_evaluations=cost,creation_validation_evaluations=cost,
        evaluated_points=2*cost,initial_evaluations=initial,cpu_seconds=cpu,
        hash_authentication_passes=1,source_integrations=0,gauge_quadratures=0,
        source_reconstructions=0,profile_fits=0,stage_calls=0,analysis_rhs_calls=4,time_steps=0,
        scientific_validator_calls=0)


def _failure_receipt(error,root,budget,initial,phase):
    receipt = Path(root)/".build/cap-mass-direction-failures"/(uuid.uuid4().hex+".json")
    receipt.parent.mkdir(parents=True,exist_ok=True)
    with receipt.open("xb") as stream:
        stream.write(encoded(dict(schema=SCHEMA+"-failed-attempt",phase=phase,
            reason=str(error),exception=type(error).__name__,evaluated_points=budget.evaluations-initial,
            cpu_seconds=time.process_time()-budget.started,retries=0,
            analysis_failure=getattr(error,"cap_mass_direction_failure",None))))
    error.failure_receipt = str(receipt)


def _validation(saved,current,cap_audit,cache,identity,bindings,budget,*,root=ROOT):
    if (saved.get("schema") != SCHEMA or saved.get("scope") != SCOPE
            or saved.get("limits") != LIMITS or saved.get("input") != identity
            or saved.get("runtime") != runtime_versions()):
        raise ValueError("cap-mass-direction declaration mismatch")
    if saved.get("audit_source_bindings") != bindings:
        verify_saved_bindings(saved.get("audit_source_bindings"),root,expected_paths=bindings)
    expected = _analysis(current,cap_audit,cache,budget)
    if not same_recomputed(saved.get("analysis"),plain(expected)):
        raise ValueError("cap-mass-direction saved reconstruction mismatch")
    receipt = saved.get("execution",{})
    initial,cpu = receipt.get("initial_evaluations"),receipt.get("cpu_seconds")
    cost = expected["evaluated_points"]
    if (type(initial) is not int or initial < 0 or initial+2*cost > LIMITS["creation"]["evaluated_points"]
            or isinstance(cpu,(bool,np.bool_)) or not isinstance(cpu,Real) or not math.isfinite(cpu)
            or not 0 <= cpu <= LIMITS["creation"]["cpu_seconds"]
            or set(receipt) != set(_execution(cost,initial,cpu))
            or any(type(v) is not int for k,v in receipt.items() if k != "cpu_seconds")
            or receipt != _execution(cost,initial,cpu)):
        raise ValueError("invalid cap-mass-direction execution receipt")
    return cost


def create(path, *, root=ROOT, budget=None):
    budget = shared_budget(budget)
    path = Path(path)
    if path.exists():
        raise FileExistsError("cap-mass-direction output is creation-only")
    initial = budget.evaluations
    phase = "producer bindings"
    try:
        bindings,runtime = _own_bindings(root),runtime_versions()
        phase = "input authentication"
        current,cap_audit,cache,identity = authenticate(root,budget,bindings)
        phase = "terminal mass-functional probes"
        result = _analysis(current,cap_audit,cache,budget)
        payload = dict(schema=SCHEMA,scope=SCOPE,limits=LIMITS,input=identity,
            audit_source_bindings=bindings,runtime=runtime,analysis=result,
            execution=_execution(result["evaluated_points"],initial,time.process_time()-budget.started))
        phase = "creation validation"
        _validation(payload,current,cap_audit,cache,identity,bindings,budget,root=root)
        if bindings != _own_bindings(root) or runtime != runtime_versions():
            raise ValueError("cap-mass-direction producer/runtime changed during calculation")
        if budget.evaluations-initial != payload["execution"]["evaluated_points"]:
            raise ValueError("unexpected charges outside the two declared analyses")
        payload["execution"]["cpu_seconds"] = time.process_time()-budget.started
        budget.check()
        raw = encoded(payload)
        budget.check()
    except Exception as error:
        _failure_receipt(error,root,budget,initial,phase)
        raise
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("xb") as stream:
        stream.write(raw)
    return payload


def check(path, *, root=ROOT, budget=None):
    budget = shared_budget(budget,checking=True)
    start = budget.evaluations
    phase = "check saved record and producer bindings"
    try:
        saved = json.loads(Path(path).read_bytes())
        bindings,runtime = _own_bindings(root),runtime_versions()
        phase = "check input authentication"
        current,cap_audit,cache,identity = authenticate(root,budget,bindings)
        phase = "check selected reconstruction"
        cost = _validation(saved,current,cap_audit,cache,identity,bindings,budget,root=root)
        if bindings != _own_bindings(root) or runtime != runtime_versions():
            raise ValueError("cap-mass-direction producer/runtime changed during authentication")
        if budget.evaluations-start != cost:
            raise ValueError("unexpected charges outside the declared analysis")
        budget.check()
    except Exception as error:
        _failure_receipt(error,root,budget,start,phase)
        raise
    return dict(checked=str(path),evaluated_points=cost,cpu_seconds=time.process_time()-budget.started,
        hash_authentication_passes=1,source_integrations=0,stage_calls=0,rhs_calls=4,time_steps=0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--check",action="store_true")
    args = parser.parse_args()
    try:
        result = check(args.output) if args.check else dict(created=str(args.output),execution=create(args.output)["execution"])
    except Exception as error:
        if hasattr(error,"failure_receipt"):
            print(json.dumps(dict(failed=str(error),diagnostic=error.failure_receipt)))
        raise
    print(json.dumps(result,sort_keys=True,allow_nan=False))


if __name__ == "__main__":
    main()
