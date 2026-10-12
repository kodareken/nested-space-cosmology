"""Seam-8 lattice: does the recorded fee anchor at 1/Omega where the wall locks?

Round-29 joined the lapse fee f=0.2610 (mu=1.040, seam 8) to the
Schroeder Omega=4 anchor (+15.2% over-exclosure) and pre-registered the
convergence hypothesis: if the fee's anchor is the wall's converter lock
(recorded only at mu=1.040), then the light-wall cases should sit
FARTHER from 1/Omega = 0.25 than the locked case. This driver computes
the same seam-8 tangent fee on all three neutral charts and measures the
distance of each to 0.25. The hypothesis is declared before the run.

    python scripts/lab.py scripts/derive_nsc_fee_convergence.py
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

OUTPUT = LAB / "results" / "development" / "nsc-fee-convergence-v1.json"
SCHEMA = "NSC-FEE-CONVERGENCE-v1"
SEAM = 8.0
OUTER_START = 24.0
MASSES = (1.002, 1.010, 1.040)
OMEGA_INVERSE = 0.25


def _window(radius, left, right, sharpness=8.0):
    return 0.5 * (1.0 + np.tanh(sharpness * (radius - left))) \
        * 0.5 * (1.0 + np.tanh(sharpness * (right - radius)))


def main():
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    rows = {}
    for mass in MASSES:
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
        fee = chain.spatial_fee(
            radius, weights, tangent["d_lapse"], SEAM,
        )["l1_mass_response_fraction_inside"]
        rows["mu-{:.3f}".format(mass)] = {
            "mass_over_rm": mass,
            "seam": SEAM,
            "fee": fee,
            "distance_to_one_over_omega": abs(fee - OMEGA_INVERSE),
        }
    nearest = min(rows, key=lambda key: rows[key]["distance_to_one_over_omega"])
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "one_over_omega": OMEGA_INVERSE,
        "hypothesis_pre_registered": (
            "if the fee's anchor is the wall's converter lock, the locked "
            "chart (mu=1.040) is closest to 1/Omega among the three"
        ),
        "rows": rows,
        "nearest_case": nearest,
        "verdict": (
            "LOCK_IS_ANCHOR" if nearest == "mu-1.040" else "LOCK_NOT_ANCHOR"
        ),
        "notes": [
            "same tangent probe as the round-10/11 seam record, seam 8",
            "1/Omega is exact; the round-1 value 0.2517 was a rounding of "
            "the recorded Omega=4 resolvent",
        ],
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_fee_convergence.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_fee_convergence.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    for key, row in rows.items():
        print("  {:<10} f(seam 8) = {:.4f}   distance to 0.25 = {:.4f}".format(
            key, row["fee"], row["distance_to_one_over_omega"]))
    print("  nearest: {}   verdict: {}".format(nearest, record["verdict"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())