"""Decomposing p_Q into its lapse, radius and shift ingredients.

Round-63 recorded the blocker: the two-ends v2 record stores only
reduced qualities, so p_Q's provenance could not be derived. The
payloads, however, hold everything the reduction needs and the cone-map
owner rebuilds the metric from them. This driver reads the last-station
payload of each case, rebuilds grid/lapse/shift, evaluates the room
qualities and decomposes the recorded Q ratio's gap power into the
three factors of Q = N^2 spacing/(r (N+|beta|)). The measured p_Q =
0.1173 is the reconciliation target; nothing is fitted.

    python scripts/lab.py scripts/derive_nsc_q_decomposition.py
"""
from __future__ import annotations

import glob
import json
import math
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

import numpy as np

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402
from recursive_horizons import nsc_rn_cone_map as cone  # noqa: E402
from recursive_horizons import nsc_rn_two_ends as ends  # noqa: E402
from recursive_horizons import nsc_rn_source as source  # noqa: E402

MARCH = LAB / "results" / "development" / "nsc-rn-inward-march-v1"
COVARIANCE = LAB / "results" / "development" / "nsc-two-ends-covariance-v1.json"
OUTPUT = LAB / "results" / "development" / "nsc-q-decomposition-v1.json"
SCHEMA = "NSC-Q-DECOMPOSITION-v1"
CASES = ("mu-1.002-eps-0.001", "mu-1.010-eps-0.001", "mu-1.040-eps-0.001")


def main():
    if not MARCH.is_dir() or not COVARIANCE.is_file():
        raise SystemExit("missing march campaign or covariance record")
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    covariance = json.loads(COVARIANCE.read_text(encoding="utf-8"))
    measured_p_q = covariance["q_gap_power_post_hoc"]
    rows = []
    for name in CASES:
        stations = sorted(
            glob.glob(str(MARCH / name / "station-*.npz")),
            key=lambda path: float(path.rsplit("station-", 1)[1][:-4]),
        )
        payload = np.load(stations[-1])
        grid = cone.station_grid(payload)
        metric = cone.coupled_metric(payload, grid)
        geometry = metric["geometry"]
        radius = np.asarray(grid["radius"], dtype=float)
        lapse = np.asarray(geometry["lapse"], dtype=float)
        shift = np.asarray(geometry["shift"], dtype=float)
        spacing = float(grid["spacing"])
        quality = ends.local_room_scales(radius, lapse, shift, spacing)["room_quality"]
        inner = int(0)
        outer = int(radius.shape[0] - 1)
        rows.append({
            "case": name,
            "q_inner": float(quality[inner]),
            "q_outer": float(quality[outer]),
            "q_ratio": float(quality[inner] / quality[outer]),
            "n_factor_squared": float((lapse[inner] / lapse[outer]) ** 2),
            "r_factor": float(radius[outer] / radius[inner]),
            "beta_factor": float(
                (lapse[outer] + abs(shift[outer])) / (lapse[inner] + abs(shift[inner]))),
        })
    light = rows[0]
    heavy = rows[2]
    ln_gap = math.log(1.0)
    model_light = source.action_reference(1.002)
    model_heavy = source.action_reference(1.040)
    gap_light = float(model_light.horizons["r_plus"]) - float(model_light.horizons["r_minus"])
    gap_heavy = float(model_heavy.horizons["r_plus"]) - float(model_heavy.horizons["r_minus"])
    ln_gap_ratio = math.log(gap_heavy / gap_light)
    def power(ratio_value):
        return math.log(ratio_value) / ln_gap_ratio
    p_n = power(heavy["n_factor_squared"] / light["n_factor_squared"]) / 2.0
    p_r = power(heavy["r_factor"] / light["r_factor"])
    p_beta = power(heavy["beta_factor"] / light["beta_factor"])
    components = {
        "p_N_single": p_n,
        "p_r": p_r,
        "p_beta_factor": p_beta,
        "p_sum": power(heavy["q_ratio"] / light["q_ratio"]),
    }
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "rows": rows,
        "component_powers": components,
        "measured_p_Q_reference": measured_p_q,
        "reconciliation": {
            "sum_minus_measured": components["p_sum"] - measured_p_q,
        },
        "notes": [
            "powers are two-point logs between the light and heavy cases",
            "the sum power is the full ratio's power: it must reconcile "
            "with p_Q measured in round 52",
        ],
        "input_sha256": {
            "lab/results/development/nsc-two-ends-covariance-v1.json":
                chain.path_sha256(COVARIANCE),
        },
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_q_decomposition.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_q_decomposition.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    for row in rows:
        print("  {:<20} Q {:.4f} N2 {:.6f} r {:.5f} beta {:.6f}".format(
            row["case"], row["q_ratio"], row["n_factor_squared"], row["r_factor"],
            row["beta_factor"]))
    print("  powers: p_N {:.4f} p_r {:.4f} p_beta {:.4f}  sum {:.4f} vs measured {:.4f}".format(
        p_n, p_r, p_beta, components["p_sum"], measured_p_q))
    print("  reconciliation: {:+.4f}".format(record["reconciliation"]["sum_minus_measured"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())