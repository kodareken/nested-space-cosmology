"""Lapse-channel interior awareness versus the infall clock.

At the round-35 target seam (28.17 r_m) the lapse channel leaves 4.9%
of its response outside while the vacuum clock leaves ~11.9% of its
budget outside. This driver makes the comparison systematic: at every
measured seam it records f(s) (L1 lapse-response fraction inside) and
T_below(s) (infall-time fraction inside), then the bias index

    beta(s) = (1 - f(s)) / T_above(s)

the lapse weight left outside per unit of clock weight left outside.
beta < 1 means the lapse channel is interior-aware (interior-biased)
relative to the clock. The seam where beta = 1 is the fairness point.

    python scripts/lab.py scripts/derive_nsc_time_awareness.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

import numpy as np

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402
from recursive_horizons import nsc_rn_source as source  # noqa: E402

INVARIANCE = LAB / "results" / "development" / "nsc-seam-invariance-v1.json"
BAND_LOCATION = LAB / "results" / "development" / "nsc-cosmic-band-location-v1.json"
OUTPUT = LAB / "results" / "development" / "nsc-time-awareness-v1.json"
SCHEMA = "NSC-TIME-AWARENESS-v1"
SEAMS = (4.0, 8.0, 16.0, 24.0, 30.0)
CAP = 32.0
MASS = 1.040


def mean_profile(profiles, seams, extra=None):
    profile = {
        float(seam_key): sum(float(values[seam_key]) for values in profiles.values())
        / len(profiles)
        for seam_key in profiles[next(iter(profiles))]
    }
    if extra is not None:
        import math
        pairs = sorted(profile.items())
        extra_value = None
        for (s_left, f_left), (s_right, f_right) in zip(pairs, pairs[1:]):
            if s_left <= extra <= s_right:
                t = (extra - s_left) / (s_right - s_left)
                ln_wanted = math.log(f_left) + t * (math.log(f_right) - math.log(f_left))
                extra_value = math.exp(ln_wanted)
                break
        if extra_value is not None:
            profile[extra] = extra_value
    return profile


def main():
    if not INVARIANCE.is_file() or not BAND_LOCATION.is_file():
        raise SystemExit("missing invariance or band-location record")
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    invariance = json.loads(INVARIANCE.read_text(encoding="utf-8"))
    band = json.loads(BAND_LOCATION.read_text(encoding="utf-8"))
    target_seam = band["mapping"]["seam_of_target_0_049"]
    profile = mean_profile(invariance["profiles"], SEAMS, extra=target_seam)
    model = source.action_reference(MASS)
    r_plus_h = float(model.horizons["r_plus"])
    fine = np.linspace(r_plus_h, CAP, 3200)
    fval = 1.0 - 2.0 * float(model.mass) / fine + float(model.magnetic_r2) / fine ** 2
    weight = 1.0 / (1.0 + np.sqrt(np.maximum(0.0, fval)))
    total = float(np.trapezoid(weight, fine))
    rows = []
    for seam in sorted(profile):
        inside = fine < seam
        t_below = float(np.trapezoid(weight[inside], fine[inside])) / total
        f_seam = profile[seam]
        t_above = 1.0 - t_below
        beta = (1.0 - f_seam) / t_above if t_above > 1e-15 else None
        rows.append({
            "seam": seam,
            "f_lapse_inside": f_seam,
            "T_below_clock": t_below,
            "T_above_clock": t_above,
            "beta": beta,
        })
    fairness = None
    for left, right in zip(rows, rows[1:]):
        if (left["beta"] - 1.0) * (right["beta"] - 1.0) < 0.0:
            t = (1.0 - left["beta"]) / (right["beta"] - left["beta"])
            fairness = left["seam"] + t * (right["seam"] - left["seam"])
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "rows": rows,
        "fairness_seam": fairness,
        "target_seam": target_seam,
        "reading": (
            "beta = lapse weight outside per clock weight outside; beta<1 "
            "is interior awareness of the lapse channel relative to the clock"
        ),
        "notes": [
            "profile f(s) is the round-33 mean across the three charts; the "
            "clock is the vacuum infall quadrature of mu=1.040",
        ],
        "input_sha256": {
            "lab/results/development/nsc-seam-invariance-v1.json":
                chain.path_sha256(INVARIANCE),
            "lab/results/development/nsc-cosmic-band-location-v1.json":
                chain.path_sha256(BAND_LOCATION),
        },
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_time_awareness.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_time_awareness.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    for row in rows:
        print("  seam {:<8} f_inside {:.4f} T_below {:.4f} T_above {:.4f} beta {:.3f}".format(
            row["seam"], row["f_lapse_inside"], row["T_below_clock"],
            row["T_above_clock"], row["beta"],
        ))
    print("  fairness seam (beta=1): {!s}".format(
        "none" if fairness is None else "{:.3f}".format(fairness)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())