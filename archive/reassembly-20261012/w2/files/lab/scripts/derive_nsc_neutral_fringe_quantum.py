"""The neutral fringe staircase: does the charged quantum recur?

Rounds 44-47 resolved the charged fringe tilt into nodal quanta with a
count ladder. The symmetry partner is the neutral chart: a fine-mu scan
at seams 4 and 8 (12.5% and 25% of the 32 r_m cap), fixed 1.04 probe,
and the measured step quantum. Quantization is compared chart-free as
quantum * node count = the seam node's magnitude in units of the mean;
the charged value was 0.00131 * 876 = 1.148.

    python scripts/lab.py scripts/derive_nsc_neutral_fringe_quantum.py
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

OUTPUT = LAB / "results" / "development" / "nsc-neutral-fringe-quantum-v1.json"
SCHEMA = "NSC-NEUTRAL-FRINGE-QUANTUM-v1"
MASSES = (1.015, 1.020, 1.025, 1.030, 1.035, 1.040, 1.045)
SEAMS = (4.0, 8.0)
POINTS = 3921


def _window(radius, left, right, sharpness=8.0):
    return 0.5 * (1.0 + np.tanh(sharpness * (radius - left))) \
        * 0.5 * (1.0 + np.tanh(sharpness * (right - radius)))


def main():
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    geometry = driver.case_geometry(1.040)
    template_grid = pg.build_grid(
        geometry["mass"], r_m=geometry["r_m"], points=POINTS, r_out=32.0,
    )
    packet = source.prepare_radial_packet(
        mass_over_rm=1.040, energy_target=1.0e-3, rank=driver.PILOT_RANK,
        model=template_grid["model"], radii=template_grid["radius"],
        weights=template_grid["weights"], derivative=template_grid["derivative"],
    )
    series = {seam: [] for seam in SEAMS}
    for mass in MASSES:
        grid = pg.build_grid(mass, r_m=1.0, points=POINTS, r_out=32.0)
        state = pg.make_state(
            grid, np.asarray(packet["phi"]), inner_mass=mass,
            killing_frequency=float(packet["radial_moments"]["requested_omega"]),
            occupations=np.asarray(packet["occupations"], dtype=float),
        )
        radius = grid["radius"]
        direction = {
            "dphi": np.asarray(packet["phi"]) * _window(radius, 24.0, 32.0)[None, :, None],
            "d_occupations": np.zeros_like(state["occupations"]),
            "d_inner_mass": 0.0,
        }
        tangent = response.stage_tangent(state, grid, direction)
        for seam in SEAMS:
            series[seam].append(chain.spatial_fee(
                radius, grid["weights"], tangent["d_lapse"], seam,
            )["l1_mass_response_fraction_inside"])
    seams = {}
    for seam in SEAMS:
        values = series[seam]
        steps = [values[i + 1] - values[i] for i in range(len(values) - 1)]
        active = [step for step in steps if abs(step) > 1e-5]
        quantum = max(abs(step) for step in active) if active else None
        equal = all(abs(abs(step) - quantum) <= 0.25 * quantum for step in active) \
            if active else True
        seams[str(seam)] = {
            "series": values,
            "steps": steps,
            "active_steps": active,
            "quantum": quantum,
            "quantum_times_points": None if quantum is None else quantum * POINTS,
        }
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "seams": seams,
        "charged_reference": {
            "quantum": 0.00131,
            "points": 876,
            "quantum_times_points": 0.00131 * 876,
        },
        "reading": (
            "chart-free quantization compares quantum * points: the seam "
            "node's magnitude in units of the mean node"
        ),
        "notes": [
            "fixed 1.04 probe across the fine-mu grids; seam 4 is the "
            "12.5%-cap partner of the charged scan",
        ],
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_neutral_fringe_quantum.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_neutral_fringe_quantum.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    for seam, entry in seams.items():
        print("  seam {:<3} active {!s} quantum {!s} q*N {!s}".format(
            seam, [round(v, 6) for v in entry["active_steps"]],
            None if entry["quantum"] is None else "{:.2e}".format(entry["quantum"]),
            None if entry["quantum_times_points"] is None
            else "{:.3f}".format(entry["quantum_times_points"])))
    print("  charged q*N = 1.148")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())