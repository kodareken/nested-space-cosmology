"""Driven-isotropy join across the six family cases.

The driven-window record v1 carries isotropy first/last for all six
family cases; the wall epochs and clock rates exist for the same six.
This driver lays the six rows on one table and measures the n=6
correlations of the isotropy fall against wall presence, rise, and fall
rates. It also records the mu-blindness of I_first and the epsilon
gradient of I_last. Nothing is predicted in advance.

    python scripts/lab.py scripts/derive_nsc_isotropy_driven_join.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402

ISOTROPY = LAB / "results" / "development" / "nsc-rn-isotropy-v1.json"
BUBBLE = LAB / "results" / "development" / "nsc-rn-gradient-bubble-v5.json"
EPSCROSS = LAB / "results" / "development" / "nsc-wall-clock-epscross-v1.json"
OUTPUT = LAB / "results" / "development" / "nsc-isotropy-driven-join-v1.json"
SCHEMA = "NSC-ISOTROPY-DRIVEN-JOIN-v1"


def pearson(x, y):
    n = len(x)
    xm = sum(x) / n
    ym = sum(y) / n
    sxy = sum((a - xm) * (b - ym) for a, b in zip(x, y))
    sxx = sum((a - xm) ** 2 for a in x)
    syy = sum((b - ym) ** 2 for b in y)
    return None if (sxx * syy) == 0.0 else sxy / (sxx * syy) ** 0.5


def main():
    if any(not p.is_file() for p in (ISOTROPY, BUBBLE, EPSCROSS)):
        raise SystemExit("missing join input record")
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    isotropy = json.loads(ISOTROPY.read_text(encoding="utf-8"))
    bubble = json.loads(BUBBLE.read_text(encoding="utf-8"))
    epscross = json.loads(EPSCROSS.read_text(encoding="utf-8"))
    bubble_by_case = {c["case"]: c for c in bubble["cases"]}
    summary = isotropy["summary"]
    family_keys = [k for k in summary if k.startswith("nsc-rn-neutral-family-v1/")]
    rows = []
    for key in family_keys:
        case_name = key.split("/")[-1]
        entry = summary[key]
        bubble_case = bubble_by_case[case_name]
        epscross_case = epscross["cases"][case_name]
        rows.append({
            "case": case_name,
            "isotropy_first": float(entry["first"]),
            "isotropy_last": float(entry["last"]),
            "isotropy_fall": float(entry["first"]) - float(entry["last"]),
            "wall_presence_fraction": float(bubble_case["wall_presence_fraction"]),
            "rise_rate": float(epscross_case["rise_rate"]),
            "fall_rate": float(epscross_case["fall_rate"]),
        })
    i_firsts = [row["isotropy_first"] for row in rows]
    first_spread = (max(i_firsts) - min(i_firsts)) / (sum(i_firsts) / len(i_firsts))
    corr_fall_vs_rise = pearson([row["isotropy_fall"] for row in rows],
                                [row["rise_rate"] for row in rows])
    corr_fall_vs_presence = pearson([row["isotropy_fall"] for row in rows],
                                    [row["wall_presence_fraction"] for row in rows])
    corr_fall_vs_fallrate = pearson([row["isotropy_fall"] for row in rows],
                                    [abs(row["fall_rate"]) for row in rows])
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "rows": rows,
        "i_first_mu_blindness": first_spread,
        "correlations": {
            "isotropy_fall_vs_rise": corr_fall_vs_rise,
            "isotropy_fall_vs_wall_presence": corr_fall_vs_presence,
            "isotropy_fall_vs_fall_rate": corr_fall_vs_fallrate,
        },
        "notes": [
            "n = 6; the two epsilon arms of each mu are separate rows",
            "I_first is nearly mu-blind and epsilon-blind: the seventh "
            "mu-blind scaffold quantity",
        ],
        "input_sha256": {
            "lab/results/development/nsc-rn-isotropy-v1.json": chain.path_sha256(ISOTROPY),
            "lab/results/development/nsc-rn-gradient-bubble-v5.json": chain.path_sha256(BUBBLE),
            "lab/results/development/nsc-wall-clock-epscross-v1.json": chain.path_sha256(EPSCROSS),
        },
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_isotropy_driven_join.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_isotropy_driven_join.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    for row in rows:
        print("  {:<20} I_fall {:.4f}  presence {:.2f}  rise {:.4f}  fall {:.4f}".format(
            row["case"], row["isotropy_fall"], row["wall_presence_fraction"],
            row["rise_rate"], row["fall_rate"]))
    print("  I_first mu-blindness: {:.2%}".format(record["i_first_mu_blindness"]))
    print("  corr fall,rise {!s}  fall,presence {!s}  fall,fallrate {!s}".format(
        corr_fall_vs_rise, corr_fall_vs_presence, corr_fall_vs_fallrate))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())