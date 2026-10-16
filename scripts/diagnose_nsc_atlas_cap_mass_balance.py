#!/usr/bin/env python3
"""Record/check the fixed terminal ledger, authenticating old results without replay."""
import argparse
import json
from pathlib import Path
import time
import uuid

import analyze_nsc_atlas_cap_mass_balance as calculation
import diagnose_nsc_atlas_cap_mass_direction as previous

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "nsc-atlas-cap-mass-balance-audit-v1"
SCOPE = calculation.SCOPE
PRIOR = "results/summaries/nsc-atlas-cap-mass-direction-v1.json"
PRIOR_SHA256 = "9092a6f4d6dd2a58f7b3b31e0dc2d5c33b1a33c05755fe94501579bd808ebd7e"
LIMITS = dict(creation=dict(evaluated_points=12_000_000, cpu_seconds=90.),
              check=dict(evaluated_points=6_000_000, cpu_seconds=45.))
OWN_SOURCES = (*previous.OWN_SOURCES,
    "src/recursive_horizons/nsc_atlas_cap_mass_derivatives.py",
    "src/recursive_horizons/nsc_atlas_cap_mass_tangent.py",
    "scripts/analyze_nsc_atlas_cap_mass_balance.py",
    "scripts/diagnose_nsc_atlas_cap_mass_balance.py",
    "tests/test_nsc_atlas_cap_mass_derivatives.py",
    "tests/test_nsc_atlas_cap_mass_tangent.py",
    "tests/test_nsc_atlas_cap_mass_balance.py",
    "tests/control_budget.py")
Budget, encoded, digest, plain = previous.Budget, previous.encoded, previous.digest, previous.plain
runtime_versions, same_recomputed = previous.runtime_versions, previous.same_recomputed


def shared_budget(value=None, *, checking=False):
    limit = LIMITS["check" if checking else "creation"]
    if value is None:
        value = Budget(maximum_evaluations=limit["evaluated_points"], cpu_seconds=limit["cpu_seconds"])
    calculation.bounded_budget(value)
    if value.maximum_evaluations > limit["evaluated_points"] or value.cpu_seconds > limit["cpu_seconds"]:
        raise ValueError("cap mass-balance allowance exceeds the declared operation")
    return value


def bindings(root):
    return {name:digest((Path(root)/name).read_bytes()) for name in OWN_SOURCES}


def authenticate(root, budget, sources):
    current, cap, cache, inherited = previous.authenticate(root, budget, sources)
    path = Path(root)/PRIOR
    if path.is_symlink():
        raise ValueError("direction summary must be a regular retained file")
    raw = path.read_bytes()
    if digest(raw) != PRIOR_SHA256:
        raise ValueError("pinned direction summary mismatch")
    summary = json.loads(raw)
    budget.check()
    return (current, cap, cache, summary), dict(inherited=inherited,
        direction_summary=dict(path=PRIOR, sha256=PRIOR_SHA256, bytes=len(raw)), previous_analysis_replays=0)


def analysis(inputs, budget):
    start = budget.evaluations
    result = plain(calculation.analyze(*inputs, budget))
    cost = budget.evaluations-start
    if type(result.get("evaluated_points")) is not int or result["evaluated_points"] != cost:
        raise ValueError("ledger point receipt disagrees with shared allowance")
    execution = result.get("execution", {})
    if any(type(execution.get(k)) is not int or execution[k]!=v for k,v in calculation.EXPECTED_CALLS.items()):
        raise ValueError("ledger call inventory mismatch")
    charges = execution.get("charges", {})
    if (set(charges)!=set(calculation.CHARGES) or any(type(v) is not int or v<0 for v in charges.values())
            or sum(charges.values())!=cost):
        raise ValueError("ledger charge inventory mismatch")
    encoded(result)
    budget.check()
    return result


def receipt(cost, initial, cpu):
    return dict(analysis_evaluations=cost, creation_validation_evaluations=cost,
        evaluated_points=2*cost, initial_evaluations=initial, cpu_seconds=cpu,
        hash_authentication_passes=1, previous_analysis_replays=0, time_steps=0)


def validate(saved, inputs, identity, sources, budget):
    if (saved.get("schema")!=SCHEMA or saved.get("scope")!=SCOPE or saved.get("limits")!=LIMITS
            or saved.get("input")!=identity or saved.get("runtime")!=runtime_versions()
            or saved.get("audit_source_bindings")!=sources):
        raise ValueError("mass balance declaration, input or producer mismatch")
    expected = analysis(inputs, budget)
    if not same_recomputed(saved.get("analysis"), expected):
        raise ValueError("saved mass balance reconstruction mismatch")
    data = saved.get("execution", {})
    initial, cpu = data.get("initial_evaluations"), data.get("cpu_seconds")
    cost = expected["evaluated_points"]
    if (type(initial) is not int or initial<0 or initial+2*cost>LIMITS["creation"]["evaluated_points"]
            or type(cpu) not in (float,int) or not 0<=cpu<=LIMITS["creation"]["cpu_seconds"]
            or any(type(v) is not int for k,v in data.items() if k!="cpu_seconds")
            or data!=receipt(cost, initial, cpu)):
        raise ValueError("invalid mass balance execution receipt")
    return cost


def operation(path, *, checking=False, root=ROOT, budget=None):
    budget = shared_budget(budget, checking=checking)
    path = Path(path)
    if not checking and path.exists():
        raise FileExistsError("mass balance output is creation-only")
    start = budget.evaluations
    phase = "producer bindings"
    try:
        sources, runtime = bindings(root), runtime_versions()
        phase = "input authentication"
        inputs, identity = authenticate(root, budget, sources)
        phase = "check reconstruction" if checking else "mass derivative and ledger"
        if checking:
            payload = json.loads(path.read_bytes())
            cost = validate(payload, inputs, identity, sources, budget)
            expected_cost = cost
        else:
            result = analysis(inputs, budget)
            cost = result["evaluated_points"]
            payload = dict(schema=SCHEMA, scope=SCOPE, limits=LIMITS, input=identity,
                runtime=runtime, audit_source_bindings=sources, analysis=result,
                execution=receipt(cost, start, time.process_time()-budget.started))
            phase = "creation validation"
            validate(payload, inputs, identity, sources, budget)
            expected_cost = 2*cost
        if bindings(root)!=sources or runtime_versions()!=runtime:
            raise ValueError("mass balance producer/runtime changed during calculation")
        if budget.evaluations-start!=expected_cost:
            raise ValueError("unexpected charges outside declared analyses")
        budget.check()
        if not checking:
            payload["execution"]["cpu_seconds"] = time.process_time()-budget.started
            raw = encoded(payload)
            budget.check()
    except Exception as error:
        target = Path(root)/".build/cap-mass-balance-failures"/(uuid.uuid4().hex+".json")
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open("xb") as stream:
            stream.write(encoded(dict(schema=SCHEMA+"-failed-attempt", phase=phase,
                exception=type(error).__name__, reason=str(error), retries=0,
                evaluated_points=budget.evaluations-start, cpu_seconds=time.process_time()-budget.started,
                analysis_failure=getattr(error, "cap_mass_balance_failure", None))))
        error.failure_receipt = str(target)
        raise
    if checking:
        return dict(checked=str(path), evaluated_points=cost, cpu_seconds=time.process_time()-budget.started,
            hash_authentication_passes=1, metric_calls=4, rhs_calls=4, tangent_calls=4, time_steps=0)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("xb") as stream:
        stream.write(raw)
    return payload


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    try:
        result = operation(args.output, checking=args.check)
    except Exception as error:
        if hasattr(error, "failure_receipt"):
            print(json.dumps(dict(failed=str(error), diagnostic=error.failure_receipt)))
        raise
    print(json.dumps(result if args.check else dict(created=str(args.output), execution=result["execution"]), sort_keys=True))


if __name__ == "__main__":
    main()
