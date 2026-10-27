#!/usr/bin/env python3
"""Evaluate the fixed nested gravitational benchmark; compact data only."""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np

from recursive_horizons.nsc_atlas_gravitational_response import nested_body, response


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results/summaries/nsc-atlas-gravitational-response-v1.json"
SOURCES = ["src/recursive_horizons/nsc_atlas_gravitational_response.py",
           "scripts/check_nsc_atlas_gravitational_response.py",
           "tests/test_nsc_atlas_gravitational_response.py"]
TARGETS = (.4, 2.0)
CONTROLS = [(2e-8, .004, 1e-4), (5e-10, .002, 2e-5), (1e-11, .001, 4e-6)]
CURVE = np.linspace(.2, 2.4, 23)
MAX_CPU = 180.
MAX_RHS = 500_000


def json_value(value):
    if isinstance(value, complex):
        return {"real": float(value.real), "imag": float(value.imag)}
    if isinstance(value, dict):
        return {str(k): json_value(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_value(v) for v in value]
    if isinstance(value, np.generic):
        return json_value(value.item())
    return value


def task(args):
    omega, termination, level, check = args
    rtol, origin, horizon = CONTROLS[level]
    return response(nested_body(), omega, termination=termination, rtol=rtol,
                    origin_fraction=origin, horizon_fraction=horizon,
                    analytic=check == "analytic", method="RK45" if check == "rk45" else "DOP853")


def main():
    started = time.perf_counter()
    body = nested_body()
    shells = [body.inner.record(), body.outer.record()]
    for shell in shells:
        assert shell["Vsecond"] > 0 and shell["sigma"] > 0
        assert max(abs(shell[k]) for k in ("V", "Vprime", "israel_gap", "force_gap")) < 3e-13

    jobs = [(w, t, level, "series") for w in TARGETS for t in ("cap", "horizon")
            for level in range(3)]
    # These do not select parameters or tolerances: same fixed geometry and frequencies.
    checks = [(w, "cap", 2, "analytic") for w in TARGETS]
    checks += [(w, t, 2, "rk45") for w in TARGETS for t in ("cap", "horizon")]
    curve_jobs = [(float(w), t, 1, "series") for w in CURVE for t in ("cap", "horizon")]
    all_jobs = jobs+checks+curve_jobs
    with ProcessPoolExecutor(max_workers=4) as pool:
        values = list(pool.map(task, all_jobs))
    cpu = sum(v["cpu_seconds"] for v in values)
    rhs = sum(v["rhs_calls"] for v in values)
    if cpu > MAX_CPU or rhs > MAX_RHS:
        raise RuntimeError("aggregate finite calculation exceeded its preparation allowance")
    lookup = {key: value for key, value in zip(all_jobs, values)}
    targets = []
    for omega in TARGETS:
        pair = {}
        for termination in ("cap", "horizon"):
            levels = [lookup[(omega, termination, k, "series")] for k in range(3)]
            final = levels[-1]
            difference = abs(levels[1]["reflection"]-final["reflection"])
            coarse_difference = abs(levels[0]["reflection"]-levels[1]["reflection"])
            other = lookup[(omega, termination, 2, "rk45")]
            method_gap = abs(other["reflection"]-final["reflection"])
            analytic_gap = (abs(lookup[(omega, "cap", 2, "analytic")]["reflection"]
                                -final["reflection"]) if termination == "cap" else 0.)
            envelope = 2*max(difference, method_gap, analytic_gap, 1e-13)
            if not difference < coarse_difference or abs(final["energy_balance_gap"]) > 1e-8:
                raise RuntimeError("response error/flux preflight failed; diagnose before continuing")
            pair[termination] = dict(final=final, levels=levels,
                                     coarse_to_middle=coarse_difference,
                                     middle_to_fine=difference, other_integrator_gap=method_gap,
                                     analytic_cap_gap=analytic_gap,
                                     complex_reflection_error_indicator=envelope,
                                     absorbed_power_error_indicator=(2*abs(final["reflection"])*envelope
                                         +envelope**2+abs(final["energy_balance_gap"])))
        difference = abs(pair["cap"]["final"]["reflection"]
                         -pair["horizon"]["final"]["reflection"])
        if difference < 10*sum(pair[t]["complex_reflection_error_indicator"] for t in pair):
            raise RuntimeError("termination contrast is not numerically resolved")
        targets.append(dict(omega=omega, comparison=pair, complex_response_difference=difference))

    # IAU2015 nominal solar mass parameter; exact conversion convention, not a fitted mass.
    c = 299_792_458.
    nominal_GM_sun = 1.3271244e20
    seconds_per_solar_mass = nominal_GM_sun/c**3
    example_mass = 1_000_000.
    mass_seconds = example_mass*seconds_per_solar_mass
    first_time = body.causal_time()
    center_time = body.causal_time(center=True)
    record = dict(
        schema="NSC-ATLAS-GRAVITATIONAL-RESPONSE-v1", date="2026-10-27",
        theory="Einstein + regional four-form + canonical surface scalars; static finite candidate",
        domain="linear axial ell>=2; radial linear background stability; no polar/nonlinear stability claim",
        observational_validation=False,
        preparation=dict(outer_mass=1., outer_radius=2.5, intermediate_mass=.1,
                         inner_radius=.25, observer_radius=10., ell=2,
                         fixed_frequencies=list(TARGETS), matter_law="p=sigma",
                         controls=CONTROLS, comparison="same outer wall/parent/exterior; regular child vs ingoing termination"),
        regions=dict(core=asdict(body.core), parent=asdict(body.parent), exterior=asdict(body.exterior)),
        interfaces=shells,
        vacuum_products=dict(Lambda_parent_M2=3*body.parent.H2,
                             Lambda_child_intermediate_mass2=3*body.core.H2*body.parent.mass**2,
                             note="conditional geometric/state relations; radius and mass ratios supplied, no absolute vacuum selection"),
        targets=targets,
        curve=[dict(omega=float(w), cap=lookup[(float(w), "cap", 1, "series")]["reflection"],
                    horizon=lookup[(float(w), "horizon", 1, "series")]["reflection"],
                    absorbed=lookup[(float(w), "horizon", 1, "series")]["absorbed_power"])
               for w in CURVE],
        causal_times=dict(first_possible_difference_proper_over_M=first_time,
                          center_null_roundtrip_proper_over_M=center_time,
                          note="characteristic/front bounds, not fitted echo peaks"),
        physical_example=dict(nominal_solar_masses=example_mass,
                              nominal_GM_sun_m3_s2=nominal_GM_sun,
                              M_seconds=mass_seconds, M_meters=mass_seconds*c,
                              first_possible_difference_seconds=first_time*mass_seconds,
                              center_null_roundtrip_seconds=center_time*mass_seconds,
                              proper_frequencies_Hz=[w/(math.sqrt(body.exterior.f(body.observer))
                                                       *2*math.pi*mass_seconds) for w in TARGETS]),
        compute=dict(rhs_calls=rhs, scientific_cpu_seconds=cpu,
                     elapsed_seconds=time.perf_counter()-started, workers=4,
                     cpu_limit=MAX_CPU, rhs_limit=MAX_RHS,
                     cpu_note="sum of solver process CPU, excludes worker import/startup; no field grids"),
        source_sha256={p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCES},
    )
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(json_value(record), indent=2, allow_nan=False)+"\n", encoding="utf-8", newline="\n")
    print(json.dumps(dict(output=OUTPUT.relative_to(ROOT).as_posix(), cpu=cpu, rhs=rhs,
                          first_time=first_time,
                          absorbed_at_high_frequency=targets[-1]["comparison"]["horizon"]["final"]["absorbed_power"])))


if __name__ == "__main__":
    main()
