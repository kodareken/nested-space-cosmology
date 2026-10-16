#!/usr/bin/env python3
"""Six-stage local SSPRK3 ledger against pinned, unreplayed SSPRK2 results.

The Richardson endpoint is an asymptotic comparison, not a truth estimate.
Checking repeats only these six stages, with no source or trajectory replay.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import time

import numpy as np

import diagnose_nsc_atlas_background_step as previous
from recursive_horizons.nsc_atlas_background_ssprk3 import advance_background_ssprk3

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = "nsc-atlas-ssprk3-step-audit-v1"
STEP_INPUT = "results/development/nsc-atlas-background-step-audit-v1.json"
STEP_SHA256 = "a94a40a9a57f5b435b32d8d43005115650348f44d84a4ab67d49c3a181907aa4"
LIMITS = dict(cpu_seconds=30., evaluated_points=150000, stage_calls=6)
SOURCES = (*previous.SOURCES, STEP_INPUT,
    "src/recursive_horizons/nsc_atlas_background_ssprk3.py",
    "tests/test_nsc_atlas_background_ssprk3.py",
    "scripts/diagnose_nsc_atlas_ssprk3_step.py",
    "tests/test_diagnose_nsc_atlas_ssprk3_step.py")
CASES, CUT = (1, 4), 82
WEIGHTS = (1/6, 1/6, 2/3)
COMPONENTS = previous.COMPONENTS
METHODS = ("rk2_coarse", "rk2_two_halves", "rk2_richardson", "ssprk3")
producer = previous.producer
Budget, plain, encoded, digest = previous.Budget, previous.plain, previous.encoded, previous.digest
compact_state, frame, _cache = previous.compact_state, previous.frame, previous._cache


def authenticate(root=ROOT):
    """Authenticate original66 bindings and the pinned local record, without stages."""
    root = Path(root)
    record, base, identity = previous.previous.authenticate(root)
    raw = (root/STEP_INPUT).read_bytes()
    if digest(raw) != STEP_SHA256:
        raise ValueError("pinned SSPRK2 local-step payload mismatch")
    saved = json.loads(raw)
    if (saved.get("schema") != previous.SCHEMA or saved.get("input") != identity
            or saved.get("limits") != previous.LIMITS
            or saved.get("runtime") != producer.runtime_versions()
            or set(saved.get("audit_source_bindings", {})) != set(previous.SOURCES)):
        raise ValueError("pinned SSPRK2 declaration mismatch")
    for name, expected in saved["audit_source_bindings"].items():
        if digest((root/name).read_bytes()) != expected:
            raise ValueError("SSPRK2 source binding mismatch: "+name)
    settings = saved["settings"]
    if (settings["cases"] != list(CASES) or settings["cut"] != CUT
            or settings["cache_indices"] != [328, 329, 330]
            or len(saved["cases"]) != 2):
        raise ValueError("SSPRK2 local comparison coordinates differ")
    return record, base, saved, dict(original=identity, local_step=dict(
        path=STEP_INPUT, sha256=digest(raw), bytes=len(raw),
        source_bindings=saved["audit_source_bindings"], local_steps_replayed=0))


def weighted_step(initial, left, mid, right, *, source_rows, gauge_rows,
                  profile, alpha, budget, counter, label):
    """Advance three stages and close the five-part ledger at actual geometry."""
    before = budget.evaluations
    if counter["stage_calls"]+3 > LIMITS["stage_calls"]:
        raise ValueError("six-stage ceiling exhausted")
    indices = (left, right, mid)
    gauges = [_cache(gauge_rows, i) for i in indices]
    sources = [_cache(source_rows, i) for i in indices]
    times = [float(g["x"]) for g in gauges]
    if len(set(times)) != 3 or any(s["x"] != g["x"] for s, g in zip(sources, gauges)):
        raise ValueError("require exactly matched cached source/gauge coordinates")
    lookup_g = dict(zip(times, gauges)); lookup_s = dict(zip(times, sources))
    h = times[1]-times[0]
    ne, nc = len(initial["departure"]["c"])+1, len(initial["cap"]["Iq"])
    phase = "SSPRK3 advance"
    try:
        try:
            result = advance_background_ssprk3(initial, times[0], times[2], times[1],
                eta=np.linspace(0., 1., ne+1), y_faces=np.linspace(0., 1., nc+1),
                source=lambda x: lookup_s[x], gauge=lambda x: lookup_g[x],
                alpha=alpha, tube=producer.TUBE, budget=budget, profile=profile,
                maximum_courant=.4, exterior_quadrature_order=2, cap_quadrature_order=4)
        except Exception as error:
            diagnostic = getattr(error, "diagnostic", {})
            count = diagnostic.get("stage_calls", 0)
            if type(count) is int and 0 <= count <= 3:
                counter["stage_calls"] += count
            raise
        if (result["stage_calls"] != 3 or result["stage_times"] != times
                or len(result["stages"]) != 3 or len(result["stage_states"]) != 3):
            raise ValueError("SSPRK3 stage inventory mismatch")
        counter["stage_calls"] += 3
        counter["max_courant"] = max(counter["max_courant"], result["max_courant"])
        phase = "independent stage and accepted readouts"
        views = [compact_state(s, g, profile, alpha, budget)
                 for s, g in zip(result["stage_states"], gauges)]
        final = compact_state(result["accepted_state"], gauges[1], profile, alpha, budget)
        chains = []
        for stage, view in zip(result["stages"], views):
            producer._equal(stage["physical_wall"], view["wall"], "SSPRK3 independent scalar readout")
            chains.append(previous.previous.wall_chain(stage, budget))
        integrate = lambda values: h*sum(w*v for w, v in zip(WEIGHTS, values))
        chain = lambda key: integrate([c[key] for c in chains])
        wall_rate = lambda key: integrate([s["wall_rates"][key] for s in result["stages"]])
        w0, w1 = views[0]["wall"], final["wall"]
        Rrk, Trk = w0["r"]+wall_rate("r"), w0["T"]+wall_rate("T")
        chi_rk = views[0]["canonical_chi"]+chain("Jchi")
        RB = final["background_chi"]-views[0]["background_chi"]-chain("B_x")
        c_rk = frame(Rrk, chi_rk, alpha)
        c_radius = frame(w1["r"], chi_rk, alpha)
        c_background = frame(w1["r"], chi_rk+RB, alpha)
        parts = dict(exterior_flux=chain("Q_c"), cap_return_flux=chain("cap_c_return"),
            frame_remainder=c_rk-w0["c"]-chain("Jc"),
            reference_geometry_remainder=c_radius-c_rk,
            background_remainder=c_background-c_radius)
        increment = w1["c"]-w0["c"]
        geometry, stage_geometry = {}, {}
        for key, rate, trial in (("r", "r_x_total", Rrk), ("T", "T_x", Trk)):
            reference = float(gauges[1][key])-float(gauges[0][key])-integrate([float(g[rate]) for g in gauges])
            geometry[key] = dict(physical_rk=trial, accepted=w1[key],
                actual_defect=w1[key]-trial, reference_defect=reference, closure=w1[key]-trial-reference)
            # Y1 is Euler at the right key; Y2 has weights1/4,1/4 at the midpoint.
            rhs0, rhs1 = [float(s["wall_rates"][key]) for s in result["stages"][:2]]
            stage_geometry[key] = [
                views[1]["wall"][key]-w0[key]-h*rhs0,
                views[2]["wall"][key]-w0[key]-h*(rhs0+rhs1)/4]
        closures = dict(canonical=final["canonical_chi"]-chi_rk-RB,
            invariant=increment-sum(parts.values()),
            invariant_flux=chain("Jc")-chain("Q_c")-chain("cap_c_return"),
            departure=final["departure_chi"]-views[0]["departure_chi"]-chain("departure_chi_rate"))
        budget.check(54)
        row = dict(label=label, cache_indices=[left, mid, right], stage_times=times,
            x_left=times[0], x_mid=times[2], x_right=times[1], h=h, alpha=alpha,
            initial=views[0], right_predictor=views[1], midpoint_state=views[2], accepted=final,
            stages=[dict(chain=c, wall_rates=s["wall_rates"], courant=h*float(s["courant_rate"]))
                    for c, s in zip(chains, result["stages"])], weights=WEIGHTS,
            rk=dict(radius=Rrk, T=Trk, chi=chi_rk, c=c_rk), background_chain=RB,
            reference_geometry=geometry, stage_geometry_defects=stage_geometry,
            increment=increment, components=parts, closures=closures,
            evaluated_points=budget.evaluations-before)
        encoded(plain(row))  # Strict finite JSON before this result can be retained.
        return plain(row)
    except Exception as error:
        raise RuntimeError(f"SSPRK3 ledger {label}; phase={phase}; indices={left},{mid},{right}; "
            f"stage_calls={counter['stage_calls']}; cumulative_points={budget.evaluations}; "
            f"cpu_seconds={time.process_time()-budget.started:.9g}; {error}") from error


def endpoint_comparison(step, old, budget):
    """Algebra only; retained RK2 endpoints are never propagated again."""
    producer._equal(old["coarse"]["initial"], step["initial"], "common saved initial state")
    coarse, fine = old["coarse"]["accepted"]["wall"], old["halves"][1]["accepted"]["wall"]
    endpoints = dict(rk2_coarse=coarse, rk2_two_halves=fine,
        rk2_richardson={k:(4*fine[k]-coarse[k])/3 for k in coarse}, ssprk3=step["accepted"]["wall"])
    initial_c = step["initial"]["wall"]["c"]
    increments = {method:row["c"]-initial_c for method, row in endpoints.items()}
    components = dict(rk2_coarse=old["coarse"]["components"], rk2_two_halves=old["fine_components"],
        rk2_richardson={k:(4*old["fine_components"][k]-old["coarse"]["components"][k])/3 for k in COMPONENTS},
        ssprk3=step["components"])
    difference = {method:{k:endpoints["ssprk3"][k]-endpoints[method][k] for k in coarse}
                  for method in METHODS[:-1]}
    budget.check(65)
    return dict(endpoints=endpoints, increments=increments, components=components,
        ssprk3_minus_reference=difference,
        component_closures={m:increments[m]-sum(components[m].values()) for m in METHODS},
        richardson_scope="asymptotic (4*two_half_steps-coarse)/3 comparison; not truth or an error bound")


def _shared(budget):
    if budget is None:
        return Budget(maximum_evaluations=150000, cpu_seconds=30.)
    if (not isinstance(budget, Budget) or type(budget.maximum_evaluations) is not int
            or not 0 < budget.maximum_evaluations <= 150000 or type(budget.evaluations) is not int
            or budget.evaluations < 0 or isinstance(budget.cpu_seconds, (bool, np.bool_))
            or not math.isfinite(budget.cpu_seconds) or not 0 < budget.cpu_seconds <= 30
            or isinstance(budget.started, (bool, np.bool_)) or not math.isfinite(budget.started)
            or budget.started > time.process_time()):
        raise ValueError("invalid shared SSPRK3 diagnostic allowance")
    budget.check()
    return budget


def audit(root=ROOT, *, budget=None):
    budget = _shared(budget)
    initial_count = budget.evaluations
    counter = dict(stage_calls=0, max_courant=0.)
    phase = "input authentication"
    try:
        record, base, old, identity = authenticate(root)
        bindings = {name:digest((Path(root)/name).read_bytes()) for name in SOURCES}
        before = budget.evaluations
        profile, summary = producer.fit_profile(record["profile_inputs"], budget)
        producer._equal(summary, record["profile"], "retained computational profile")
        profile_points = budget.evaluations-before
        gauge = base["gauge"]["rows"]
        x = np.asarray(gauge["x"], dtype=float)
        if (x.shape != (513,) or np.any(np.diff(x) <= 0) or not np.all(np.isfinite(x))
                or not np.allclose(np.diff(x), np.diff(x)[0], rtol=1e-10, atol=1e-14)):
            raise ValueError("require original uniform513 cached coordinates")
        cases = []
        comparison_start = 0
        for index, case_id in enumerate(CASES):
            case, prior = record["cases"][case_id], old["cases"][index]
            label = "zero" if case_id == 1 else "finite"
            phase = f"case={case_id} amplitude={label}"
            if ((case["cells"], case["steps"], case["amplitude"]) != (32, 256, label)
                    or (prior["case"], prior["amplitude"], prior["alpha"]) != (case_id, label, case["alpha"])):
                raise ValueError("require matching saved32/256 amplitudes")
            sample = case["snapshots"][CUT]
            if sample["x"] != x[328]:
                raise ValueError("saved starting coordinate differs")
            row = weighted_step(sample["state"], 328, 329, 330,
                source_rows=base["initial_data"][label]["snapshot_rows"], gauge_rows=gauge,
                profile=profile, alpha=case["alpha"], budget=budget, counter=counter, label=label)
            producer._equal(row["initial"]["wall"], sample["physical"]["wall"], "retained starting wall")
            count = budget.evaluations
            comparison = endpoint_comparison(row, prior, budget)
            comparison_start += budget.evaluations-count
            cases.append(dict(case=case_id, amplitude=label, alpha=case["alpha"],
                              step=row, comparison=comparison))
        if counter["stage_calls"] != 6:
            raise ValueError("diagnostic needs exactly six coupled stages")
        paired = {}
        for method in METHODS:
            zero, finite = [c["comparison"] for c in cases]
            parts = {k:finite["components"][method][k]-zero["components"][method][k] for k in COMPONENTS}
            increment = finite["increments"][method]-zero["increments"][method]
            paired[method] = dict(endpoint={k:finite["endpoints"][method][k]-zero["endpoints"][method][k]
                for k in finite["endpoints"][method]}, increment=increment, components=parts,
                closure=increment-sum(parts.values()))
        paired_difference = {m:paired["ssprk3"]["increment"]-paired[m]["increment"] for m in METHODS[:-1]}
        budget.check(47)
        phase = "final input authentication"
        _, _, _, after = authenticate(root)
        if identity != after or bindings != {name:digest((Path(root)/name).read_bytes()) for name in SOURCES}:
            raise ValueError("SSPRK3 diagnostic input or producer changed during evaluation")
        budget.check()
        return plain(dict(schema=SCHEMA, input=identity, audit_source_bindings=bindings, limits=LIMITS,
            runtime=producer.runtime_versions(), settings=dict(cases=CASES, cut=CUT,
                cache_indices=[328, 329, 330], stage_order=[328, 330, 329], weights=WEIGHTS,
                time_scheme="SSPRK(3,3); one step per saved amplitude", maximum_courant=.4,
                exterior_quadrature_order=2, cap_quadrature_order=4),
            cases=cases, paired=paired, paired_ssprk3_minus_reference_increment=paired_difference,
            execution=dict(initial_evaluations=initial_count, evaluated_points=budget.evaluations,
                profile_fit_points=profile_points, local_step_points=sum(c["step"]["evaluated_points"] for c in cases),
                comparison_points=comparison_start+47, cpu_seconds=time.process_time()-budget.started,
                stage_calls=6, local_ssprk3_steps=2, old_ssprk2_steps_replayed=0,
                full_trajectory_replays=0, source_integrations=0, profile_stage_extractions=0,
                accepted_endpoint_extra_stage_calls=0, hash_authentication_passes=2,
                maximum_courant=counter["max_courant"], retries=0),
            scope="one local interval with unchanged physical inputs; weighted increment identities and sampled method sensitivity; Richardson is asymptotic, not truth; no global convergence or error enclosure"))
    except Exception as error:
        raise RuntimeError(f"SSPRK3 audit phase={phase}; stage_calls={counter['stage_calls']}; "
            f"cumulative_points={budget.evaluations}; cpu_seconds={time.process_time()-budget.started:.9g}; {error}") from error


def create(path, *, root=ROOT):
    path = Path(path)
    if path.exists():
        raise FileExistsError("SSPRK3 audit output is creation-only")
    result = audit(root)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("xb") as stream:
        stream.write(encoded(result))
    return result


def check(path, *, root=ROOT):
    saved = json.loads(Path(path).read_bytes())
    current = audit(root)
    for row in (saved, current):
        cpu = row["execution"].pop("cpu_seconds")
        if isinstance(cpu, bool) or not math.isfinite(cpu) or not 0 <= cpu <= 30:
            raise ValueError("invalid SSPRK3 diagnostic CPU")
    if encoded(saved) != encoded(current):
        raise ValueError("SSPRK3 local ledger recomputation mismatch")
    return dict(checked=str(path), stage_calls=6, local_ssprk3_steps=2, old_ssprk2_steps_replayed=0,
                full_trajectory_replays=0, source_integrations=0,
                evaluated_points=current["execution"]["evaluated_points"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.check:
        result = check(args.output)
    else:
        row = create(args.output)
        result = dict(created=str(args.output), execution=row["execution"],
                      paired_ssprk3_minus_reference_increment=row["paired_ssprk3_minus_reference_increment"])
    print(json.dumps(result, sort_keys=True, allow_nan=False))


if __name__ == "__main__":
    main()
