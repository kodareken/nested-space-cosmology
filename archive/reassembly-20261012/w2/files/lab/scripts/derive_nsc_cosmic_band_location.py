"""Where does the cosmic local band live on this chart?

Every exterior channel of the founding riddle was refuted, and the
interior is the sink side of an exact one-way funnel. The remaining
positive object is a location question: at which seam does the measured
lapse-profile leave exactly the pre-registered local fraction
[0.03, 0.10] (BBN-anchored target 0.049)? This driver maps the band
onto the seam curve from the round-33 invariance record (log-linear
between the measured points) for all three charts. The mapping is pure
algebra on recorded numbers; the cosmic reading is labeled hypothesis.

    python scripts/lab.py scripts/derive_nsc_cosmic_band_location.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402

INVARIANCE = LAB / "results" / "development" / "nsc-seam-invariance-v1.json"
OUTPUT = LAB / "results" / "development" / "nsc-cosmic-band-location-v1.json"
SCHEMA = "NSC-COSMIC-BAND-LOCATION-v1"
SEAMS = (4.0, 8.0, 16.0, 24.0, 30.0)
CAP = 32.0
TARGET = 0.049
BAND = (0.03, 0.10)


def seam_for_fraction(profile, wanted_local_fraction):
    """Log-linear interpolation of local = 1 - f to the seam where the
    profile holds exactly the wanted local fraction."""
    import math
    wanted_f = 1.0 - wanted_local_fraction
    pairs = sorted(profile.items())
    for (s_left, f_left), (s_right, f_right) in zip(pairs, pairs[1:]):
        if (f_left - wanted_f) * (f_right - wanted_f) <= 0.0:
            if abs(f_right - f_left) < 1e-15:
                return 0.5 * (s_left + s_right)
            ln_left = math.log(f_left)
            ln_right = math.log(f_right)
            ln_wanted = math.log(wanted_f)
            t = (ln_wanted - ln_left) / (ln_right - ln_left) if ln_right != ln_left else 0.5
            return s_left + t * (s_right - s_left)
    raise ValueError("fraction not bracketed by the profile")


def main():
    if not INVARIANCE.is_file():
        raise SystemExit("missing seam invariance record: " + str(INVARIANCE))
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    record = json.loads(INVARIANCE.read_text(encoding="utf-8"))
    profiles = {}
    for key, profile in record["profiles"].items():
        profiles[key] = {float(seam_key): float(value) for seam_key, value in profile.items()}
    mean_profile = {
        seam: sum(profiles[key][seam] for key in profiles) / len(profiles) for seam in SEAMS
    }
    band_lo_seam = seam_for_fraction(mean_profile, BAND[0])
    target_seam = seam_for_fraction(mean_profile, TARGET)
    band_hi_seam = seam_for_fraction(mean_profile, BAND[1])
    per_case = {}
    for key, profile in profiles.items():
        per_case[key] = {
            "local_0_03_seam": seam_for_fraction(profile, BAND[0]),
            "local_0_049_seam": seam_for_fraction(profile, TARGET),
            "local_0_10_seam": seam_for_fraction(profile, BAND[1]),
        }
    body = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "mapping": {
            "local_fraction_band": list(BAND),
            "target_local_fraction": TARGET,
            "seam_of_band_low_0_03": band_lo_seam,
            "seam_of_target_0_049": target_seam,
            "seam_of_band_high_0_10": band_hi_seam,
            "band_seam_interval": [band_lo_seam, band_hi_seam],
            "band_cap_fractions": [band_lo_seam / CAP, band_hi_seam / CAP],
            "distance_from_cap": [CAP - band_hi_seam, CAP - band_lo_seam],
        },
        "per_case": per_case,
        "status": "observation-level mapping; the cosmic location reading is a hypothesis",
        "notes": [
            "the band [0.03, 0.10] is the pre-registered BBN-adjacent "
            "acceptance of the influence question, unchanged since round 1",
            "interpolation is log-linear in f between the measured seams; "
            "no further model enters",
        ],
        "input_sha256": {
            "lab/results/development/nsc-seam-invariance-v1.json":
                chain.path_sha256(INVARIANCE),
        },
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_cosmic_band_location.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_cosmic_band_location.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(body, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    print("  band seam interval: [{:.3f}, {:.3f}] r_m   target 0.049 at {:.3f} r_m".format(
        band_lo_seam, band_hi_seam, target_seam))
    print("  cap fractions: [{:.2%}, {:.2%}]   distance from cap [{:.2f}, {:.2f}] r_m".format(
        band_lo_seam / CAP, band_hi_seam / CAP, CAP - band_hi_seam, CAP - band_lo_seam))
    for key, entry in per_case.items():
        print("  {:<10} 0.03@ {:.3f}  0.049@ {:.3f}  0.10@ {:.3f}".format(
            key, entry["local_0_03_seam"], entry["local_0_049_seam"],
            entry["local_0_10_seam"],
        ))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())