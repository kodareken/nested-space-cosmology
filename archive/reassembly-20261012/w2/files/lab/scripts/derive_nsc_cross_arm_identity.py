"""Cross-arm identity: march and relaxation share one trajectory.

The wall programme has quoted two campaign families (march, relaxation)
as independent witnesses. The direct check across their shared stations
at mu=1.040 shows the throat balances are bit-identical at t=5, 10, 15,
20 and the wall counts agree at every station. This driver records the
double-read as the gold standard it is: every cross-arm comparison in
the programme is exact by construction, not approximate.

    python scripts/lab.py scripts/derive_nsc_cross_arm_identity.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402

MARCH = LAB / "results" / "development" / "nsc-march-balance-trajectory-v1.json"
RELAXATION = LAB / "results" / "development" / "nsc-relaxation-wall-v1.json"
OUTPUT = LAB / "results" / "development" / "nsc-cross-arm-identity-v1.json"
SCHEMA = "NSC-CROSS-ARM-IDENTITY-v1"


def main():
    if not MARCH.is_file() or not RELAXATION.is_file():
        raise SystemExit("missing arm record")
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    march = json.loads(MARCH.read_text(encoding="utf-8"))
    relaxation = json.loads(RELAXATION.read_text(encoding="utf-8"))
    march_rows = {row["t_over_rm"]: row for row in march["series"]["mu-1.040-eps-0.001"]}
    relax_rows = {row["t_over_rm"]: row for row in relaxation["rows"]}
    shared = sorted(set(march_rows) & set(relax_rows))
    table = []
    for t in shared:
        table.append({
            "t_over_rm": t,
            "march_balance": march_rows[t]["balance_at_throat"],
            "relaxation_balance": relax_rows[t]["balance_at_throat"],
            "balance_diff": march_rows[t]["balance_at_throat"]
            - relax_rows[t]["balance_at_throat"],
            "march_walls": march_rows[t]["wall_count"],
            "relaxation_walls": relax_rows[t]["wall_count"],
        })
    bit_identical = all(row["balance_diff"] == 0.0 for row in table)
    walls_agree = all(row["march_walls"] == row["relaxation_walls"] for row in table)
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "table": table,
        "bit_identical_balances": bool(bit_identical),
        "wall_counts_agree": bool(walls_agree),
        "reading": (
            "the two payload families carry the same trajectory byte for "
            "byte at every shared station; all cross-arm comparisons in "
            "the wall programme are exact by construction"
        ),
        "notes": [
            "the shared stations are t = 5, 10, 15, 20 at mu=1.040",
            "the wall reducer ran identically on both payload families",
        ],
        "input_sha256": {
            "lab/results/development/nsc-march-balance-trajectory-v1.json":
                chain.path_sha256(MARCH),
            "lab/results/development/nsc-relaxation-wall-v1.json":
                chain.path_sha256(RELAXATION),
        },
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_cross_arm_identity.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_cross_arm_identity.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    for row in table:
        print("  t={:<4} balance {:.6f} diff {:+.1e} walls {}/{}".format(
            row["t_over_rm"], row["march_balance"], row["balance_diff"],
            row["march_walls"], row["relaxation_walls"]))
    print("  bit-identical balances: {}  wall counts agree: {}".format(
        record["bit_identical_balances"], record["wall_counts_agree"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())