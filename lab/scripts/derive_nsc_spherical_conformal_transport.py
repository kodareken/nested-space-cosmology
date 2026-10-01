#!/usr/bin/env python3
"""Surface-resolved turnover assessment of the completed conformal successor.

Entry and exit surfaces are assessed separately from net flux. Opposing
nonzero entry/exit is not vetoed because the net content-flux contribution
cancels. The finite measurement has empirical space/time/frame indicators;
it is not an observable-error certificate or a continuum transport proof.
"""
from __future__ import annotations
import json
from pathlib import Path
import sys
import numpy as np

import derive_nsc_spherical_conformal_episode as episode
import derive_nsc_spherical_conformal_continuation as continuation

LAB = Path(__file__).resolve().parents[1]
OUT = LAB / "results/development/nsc-spherical-conformal-transport-v2.json"
DRIVER = Path(__file__).resolve()
CHANNELS = ("normal_left_outward_flux", "normal_right_outward_flux", "normal_boundary_flux")
PRIMARY = "nf512_dt_0.0005"
TIME = "nf512_dt_0.0010"
SPACE = "nf256_dt_0.0005"
ONE_PERCENT = .01


def integrate(times, values):
    values, times = np.asarray(values, dtype=float), np.asarray(times, dtype=float)
    dt = np.diff(times)
    signed = float(np.sum(.5 * dt * (values[:-1] + values[1:])))
    absolute = float(np.sum(.5 * dt * (abs(values[:-1]) + abs(values[1:]))))
    uniform = len(dt) > 0 and max(abs(dt - dt[0])) < 1e-12
    if uniform:
        signed_indicator = episode.integrate(values, float(dt[0]))
        absolute_indicator = episode.integrate(abs(values), float(dt[0]))
    else:
        signed_indicator = {"indicator": None, "simpson": None}
        absolute_indicator = {"indicator": None, "simpson": None}
    return {"signed_integral": signed, "absolute_integral": absolute,
            "signed_simpson": signed_indicator["simpson"], "absolute_simpson": absolute_indicator["simpson"],
            "signed_frame_indicator": signed_indicator["indicator"], "absolute_frame_indicator": absolute_indicator["indicator"],
            "peak_absolute": float(np.max(abs(values))), "endpoint": float(values[-1])}


def assess_surfaces(physical_series):
    assessments = []
    for window in range(4):
        for channel in CHANNELS:
            per_case = {}
            for name, rows in physical_series.items():
                times = [row["time"] for row in rows]
                values = [row["windows"][window][channel] for row in rows]
                per_case[name] = integrate(times, values)
            primary = per_case[PRIMARY]
            effect = primary["absolute_integral"]
            time_indicator = abs(effect - per_case[TIME]["absolute_integral"])
            space_indicator = abs(effect - per_case[SPACE]["absolute_integral"])
            frame_indicator = primary["absolute_frame_indicator"]
            resolved = effect > 0 and frame_indicator is not None and max(time_indicator, space_indicator, frame_indicator) < ONE_PERCENT * effect
            assessments.append({"interval": [2 * window, 2 * window + 2], "channel": channel, "per_case": per_case,
                                "absolute_integral_effect": effect, "time_indicator": time_indicator,
                                "space_indicator": space_indicator, "frame_quadrature_indicator": frame_indicator,
                                "time_fraction_of_effect": None if effect == 0 else time_indicator / effect,
                                "space_fraction_of_effect": None if effect == 0 else space_indicator / effect,
                                "frame_fraction_of_effect": None if effect == 0 or frame_indicator is None else frame_indicator / effect,
                                "resolved_one_percent_indicators": bool(resolved), "certified": False})
    return assessments


def run():
    episode.refuse_existing_outputs(OUT)
    raw = json.loads(continuation.OUT.read_text())
    if raw["gauge"] != "conformal" or not raw["all_target_reached"]:
        raise ValueError("completed same-source conformal successor required")
    before = {"episode_v2": episode.sha256(continuation.OUT), "payload_v2": episode.sha256(continuation.NPZ),
              "producer_v2": episode.sha256(continuation.DRIVER), "surface_producer": episode.sha256(DRIVER)}
    assessments = assess_surfaces(raw["physical_series"])
    actual_surfaces = [row for row in assessments if row["channel"] != "normal_boundary_flux"]
    resolved = any(row["resolved_one_percent_indicators"] for row in actual_surfaces)
    localization = {name: {"shell_leader_unchanged": value["shell_leader_unchanged"],
                           "probability_leader_unchanged": value["probability_leader_unchanged"],
                           "shell_leader_share_min": value["normal_shell_leader_share_min"],
                           "probability_leader_share_min": value["probability_leader_share_min"]}
                    for name, value in raw["physical_results"].items()}
    maintained = all(value["shell_leader_unchanged"] and value["probability_leader_unchanged"] for value in localization.values())
    record = {"schema": "NSC-SPHERICAL-CONFORMAL-TRANSPORT-v2", "gauge": "conformal", "start": continuation.START, "end": continuation.TARGET,
              "source_bindings": before, "surface_assessments": assessments, "localization": localization,
              "resolved_regional_surface_exchange": resolved, "maintained_localization": maintained,
              "net_zero_is_a_veto": False, "flow_effect": "time integral of absolute normal energy flux at each natural width-2 surface; signed/net quantities retained",
              "refinement_domain": "nf256/nf512, dtcaps .001/.0005, .005 saved physical frames, same continuous original observer and Cauchy trajectory",
              "criterion": "each surface uses its own absolute-flow effect and space/time/Simpson-trapezoid indicators below one percent",
              "continuum_error_certified": False, "observable_error_certified": False, "renewal": False,
              "verdict": "MEASURED_MAINTAINED_LOCALIZATION_WITH_RESOLVED_SURFACE_FLOW" if maintained and resolved else "SURFACE_EXCHANGE_UNRESOLVED",
              "producer_net_assessment_scope": "the sealed episode-v2 boundary_assessment measures net content-flux contribution; this record adds each surface without rewriting it"}
    record["source_unchanged"] = episode.sha256(continuation.OUT) == before["episode_v2"] and episode.sha256(continuation.NPZ) == before["payload_v2"]
    OUT.write_text(json.dumps(episode.jsonable(record), indent=2, allow_nan=False) + "\n")
    print("saved", OUT, record["verdict"], flush=True)
    for row in actual_surfaces:
        print(row["interval"], row["channel"], "effect", row["absolute_integral_effect"], "fractions", row["time_fraction_of_effect"], row["space_fraction_of_effect"], row["frame_fraction_of_effect"], "resolved", row["resolved_one_percent_indicators"], flush=True)
    return record


def verify_saved(*, source_ref=None):
    before = episode.sha256(OUT)
    record = json.loads(OUT.read_text())
    if record.get("schema") != "NSC-SPHERICAL-CONFORMAL-TRANSPORT-v2":
        raise ValueError("transport schema differs")
    episode.check_recorded_sources(record.get("source_bindings"),
        {"episode_v2": continuation.OUT, "payload_v2": continuation.NPZ,
         "producer_v2": continuation.DRIVER, "surface_producer": DRIVER}, source_ref=source_ref)
    raw = json.loads(continuation.OUT.read_text())
    episode.check_saved_values(assess_surfaces(raw["physical_series"]), record["surface_assessments"], "surface assessment")
    if record.get("net_zero_is_a_veto") is not False or record.get("observable_error_certified") is not False:
        raise ValueError("transport claim domain differs")
    if episode.sha256(OUT) != before:
        raise RuntimeError("read-only check changed transport bytes")
    return {"status": record["verdict"], "wrote": False, "source_ref": source_ref}


if __name__ == "__main__":
    if sys.argv[1:] == ["--check"]:
        print(json.dumps(verify_saved(source_ref=episode.SEALED_SOURCE_REF)))
    elif len(sys.argv) == 1:
        run()
    else:
        raise SystemExit("Use --check for sealed evidence, or no arguments to create new outputs")
