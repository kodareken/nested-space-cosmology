#!/usr/bin/env python3
"""Exact frozen-group comparison at the saved normal proper clock.

The sealed episode-v2 interpolated its frozen control between numerical
step times. This postprocessor evaluates the same available analytic frozen
group at the exact matched clock endpoint instead. It does not replace the
trajectory, observer, clock protocol, or a source and keeps the older
interpolated result as a measured temporal indicator.
"""
import json
from pathlib import Path

import numpy as np

import derive_nsc_spherical_conformal_episode as episode
import derive_nsc_spherical_conformal_continuation as continuation

LAB = Path(__file__).resolve().parents[1]
OUT = LAB / "results/development/nsc-spherical-conformal-clock-v2.json"
DRIVER = Path(__file__).resolve()


def run():
    core = json.loads(continuation.OUT.read_text())
    v1 = json.loads(episode.OUT.read_text())
    if not core["all_target_reached"]:
        raise ValueError("complete conformal successor required")
    bindings = {"episode_v2": episode.sha256(continuation.OUT), "payload_v2": episode.sha256(continuation.NPZ),
                "episode_v1": episode.sha256(episode.OUT), "analytic_frozen_producer": episode.sha256(episode.DRIVER),
                "driver": episode.sha256(DRIVER)}
    results = {}
    with np.load(continuation.NPZ, allow_pickle=False) as saved:
        for spec in episode.RUN_PLAN:
            original = episode.CauchyState(**{name: saved[f"nf{spec['nf']}_original_T0_{name}"] for name in episode.STATE_NAMES})
            final0, final1 = saved[spec["name"] + "_final_phi0"], saved[spec["name"] + "_final_phi1"]
            observer = np.vstack((original.phi0, original.phi1))[:, :2]
            clock = core["results"][spec["name"]]["final"]["proper_clock"]
            original_rate = v1["results"][spec["name"]]["initial"]["clock_rate"]
            frozen_time = clock / original_rate
            if frozen_time > core["results"][spec["name"]]["time"] + 1e-12:
                raise ValueError("common clock would require interpolating the coupled endpoint")
            frozen0, frozen1 = episode.frozen_columns(original, episode.galerkin.PERIOD, episode.galerkin.KAPPA, frozen_time)
            coupled = episode.occupation(observer, final0, final1, episode.galerkin.OCCUPATIONS)
            frozen = episode.occupation(observer, frozen0, frozen1, episode.galerkin.OCCUPATIONS)
            difference = coupled - frozen
            legacy = np.asarray(core["results"][spec["name"]]["proper_clock_comparison"]["difference"])
            results[spec["name"]] = {"proper_clock": clock, "original_clock_rate": original_rate,
                                    "frozen_coordinate_time": frozen_time, "coupled_coordinate_time": core["results"][spec["name"]]["time"],
                                    "coupled_occupation": coupled, "frozen_occupation": frozen, "difference": difference,
                                    "interpolated_control_difference": legacy, "interpolated_minus_exact": legacy - difference}
    primary = results["nf512_dt_0.0005"]["difference"]
    time_gap = abs(primary - results["nf512_dt_0.0010"]["difference"])
    space_gap = abs(primary - results["nf256_dt_0.0005"]["difference"])
    record = {"schema": "NSC-SPHERICAL-CONFORMAL-CLOCK-v2", "gauge": "conformal", "source_bindings": bindings,
              "protocol": "saved coupled endpoint clock tau(x=1), exact analytic frozen group at coordinate time tau/N0(x=1), original T0 mode observer",
              "trajectory_evolved": False, "frozen_time_interpolation_used": False, "coupled_time_interpolation_used": False,
              "results": results, "time_indicator": time_gap, "space_indicator": space_gap,
              "time_fraction_of_effect": time_gap / abs(primary), "space_fraction_of_effect": space_gap / abs(primary),
              "one_percent_indicators": bool(np.max(time_gap / abs(primary)) < .01 and np.max(space_gap / abs(primary)) < .01),
              "clock_quadrature_certified": False, "observable_error_certified": False, "renewal": False}
    record["source_unchanged"] = episode.sha256(continuation.OUT) == bindings["episode_v2"] and episode.sha256(continuation.NPZ) == bindings["payload_v2"]
    OUT.write_text(json.dumps(episode.jsonable(record), indent=2, allow_nan=False) + "\n")
    print("saved", OUT, "difference", primary, "timefractions", record["time_fraction_of_effect"], "spacefractions", record["space_fraction_of_effect"], flush=True)
    return record


if __name__ == "__main__":
    run()
