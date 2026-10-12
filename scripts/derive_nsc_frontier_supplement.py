"""Frontier supplement: t-extensions are operational, the 4th gap stays gated.

The round-63/78 frontier records listed station extensions as
production. Round-88 re-examined the registry: the inward-march runner
marches any station list on the three legal masses, and a t30 campaign
for the light case was launched and is integrating. This supplement
records the precise new shape of the frontier: the exponent-provenance
cell remains blocked only by the packet preparer's mass pin; every
t-extension on the three recorded masses is now operational machinery,
with the t30 campaign as its first new instance.

    python scripts/lab.py scripts/derive_nsc_frontier_supplement.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402

OUTPUT = LAB / "results" / "development" / "nsc-frontier-supplement-v1.json"
SCHEMA = "NSC-FRONTIER-SUPPLEMENT-v1"


def main():
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    prepare = (LAB / "results" / "development" / "nsc-rn-inward-march-t30-v1"
               / "mu-1.002-eps-0.001" / "prepare.json")
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "supplements": ["nsc-frontier-amendment-v1.json"],
        "cells_now": {
            "exponent_provenance_fourth_gap": (
                "still blocked: source.prepare_radial_packet accepts only "
                "mass 1.002/1.01/1.04; no fourth gap point possible from "
                "recorded bytes or the march runner"
            ),
            "t_extensions_on_recorded_masses": (
                "now operational: derive_nsc_rn_inward_march.py marches "
                "any station list on the legal masses; the t30 campaign "
                "for mu=1.002 was launched in round 88 and is integrating"
            ),
            "paper_level_coupling": "unchanged, out of lab domain",
        },
        "launch_trail": {
            "campaign": "nsc-rn-inward-march-t30-v1",
            "prepare_exists": bool(prepare.is_file()),
            "stations": [25.0, 30.0],
        },
        "notes": [
            "the freeze-gate refusal of round 71 applies to the campaign "
            "CLI; the march runner is a separate operating corridor",
        ],
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_frontier_supplement.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_frontier_supplement.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    print("  prepare exists:", record["launch_trail"]["prepare_exists"])
    print("  4th gap:", "blocked (preparer pin)")
    print("  t-extensions:", "operational (march runner)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())