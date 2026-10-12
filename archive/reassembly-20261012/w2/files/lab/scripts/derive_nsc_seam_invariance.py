"""Seam-profile invariance across the three charts.

Round-32 pinned the seam-8 fee as chart-invariant. This driver takes the
remaining step: the whole five-seam profile at mu = 1.002, 1.010, 1.040
with the identical tangent probe, and records the pairwise maximum
deviation. If the curves coincide, the outside scaffold is complete and
the programme's structure questions become purely interior ones.

    python scripts/lab.py scripts/derive_nsc_seam_invariance.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))
if str(LAB / "scripts") not in sys.path:
    sys.path.insert(0, str(LAB / "scripts"))

import numpy as np

import derive_nsc_rn_parent as driver
from recursive_horizons import nsc_influence_chain as chain
from recursive_horizons import nsc_rn_pg as pg
from recursive_horizons import nsc_rn_response as response
from recursive_horizons import nsc_rn_source as source

OUTPUT = LAB / "results" / "development" / "nsc-seam-invariance-v1.json"
SCHEMA = "NSC-SEAM-INVARIANCE-v1"
SEAMS = (4.0, 8.0, 16.0, 24.0, 30.0)
OUTER_START = 24.0
MASSES = (1.002, 1.010, 1.040)


def _window(radius, left, right, sharpness=8.0):
    return 0.5 * (1.0 + np.tanh(sharpness * (radius - left))) \
        * 0.5 * (1.0 + np.tanh(sharpness * (right - radius)))


def probe(mass):
    geometry = driver.case_geometry(mass)
    grid = pg.build_grid(
        geometry["mass"], r_m=geometry["r_m"], points=geometry["route"]["points"],
        r_out=geometry["route"]["outer"],
    )
    packet = source.prepare_radial_packet(
        mass_over_rm=mass, energy_target=1.0e-3, rank=driver.PILOT_RANK,
        model=grid["model"], radii=grid["radius"], weights=grid["weights"],
        derivative=grid["derivative"],
    )
    state = pg.make_state(
        grid, np.asarray(packet["phi"]),
        inner_mass=geometry["mass"],
        killing_frequency=float(packet["radial_moments"]["requested_omega"]),
        occupations=np.asarray(packet["occupations"], dtype=float),
    )
    radius = grid["radius"]
    weights = grid["weights"]
    direction = {
        "dphi": np.asarray(packet["phi"]) * _window(
            radius, OUTER_START, float(geometry["route"]["outer"]))[None, :, None],
        "d_occupations": np.zeros_like(state["occupations"]),
        "d_inner_mass": 0.0,
    }
    tangent = response.stage_tangent(state, grid, direction)
    profile = {
        float(seam): chain.spatial_fee(radius, weights, tangent["d_lapse"], seam)[
            "l1_mass_response_fraction_inside"]
        for seam in SEAMS
    }
    return profile


def main():
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    profiles = {}
    for mass in MASSES:
        profiles["mu-{:.3f}".format(mass)] = probe(mass)
    deviations = []
    keys = list(profiles)
    for index in range(len(keys)):
        for other in range(index + 1, len(keys)):
            a = profiles[keys[index]]
            b = profiles[keys[other]]
            deviations.append({
                "pair": [keys[index], keys[other]],
                "max_abs_deviation": max(abs(a[seam] - b[seam]) for seam in SEAMS),
                "per_seam": {
                    str(seam): abs(a[seam] - b[seam]) for seam in SEAMS
                },
            })
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "profiles": profiles,
        "pairwise": deviations,
        "max_pairwise_deviation": max(item["max_abs_deviation"] for item in deviations),
        "verdict": (
            "PROFILES_COINCIDE" if max(item["max_abs_deviation"] for item in deviations) < 0.01
            else "PROFILES_DIVERGE"
        ),
        "reading": (
            "identical probe, three charts; the inverted mass direction and "
            "the identical outside scaffold both predict coincidence"
        ),
        "notes": [
            "the round-11 recorded profile for mu=1.040 is recomputed here "
            "with the same code path, not imported",
        ],
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_seam_invariance.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_seam_invariance.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    for key, profile in profiles.items():
        parts = "  ".join(
            "f({:g})={:.4f}".format(seam, profile[seam]) for seam in SEAMS
        )
        print("  {:<10} {}".format(key, parts))
    print("  max pairwise deviation: {:.4f}   verdict: {}".format(
        record["max_pairwise_deviation"], record["verdict"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())