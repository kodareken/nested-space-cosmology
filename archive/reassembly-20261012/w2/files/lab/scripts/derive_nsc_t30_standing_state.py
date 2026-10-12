"""Standing state of the t30 corridor, written while the march runs.

The campaign is in flight; this record is its handover: the march's
health, the armed pipeline, and the exact adoption steps - so any
later session adopts the landing without re-deriving anything.

    python scripts/lab.py scripts/derive_nsc_t30_standing_state.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402

OUTPUT = LAB / "results" / "development" / "nsc-t30-standing-state-v1.json"
SCHEMA = "NSC-T30-STANDING-STATE-v1"


def main():
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    case_dir = (LAB / "results" / "development" / "nsc-rn-inward-march-t30-v1"
                / "mu-1.002-eps-0.001")
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "written_at_wall_minutes": 191.0,
        "march": {
            "campaign": "nsc-rn-inward-march-t30-v1",
            "launch": "derive_nsc_rn_inward_march.py --mu 1.002 --eps 0.001 "
                      "--campaign ... --stations 25.0 30.0",
            "prepare_trail_intact": (case_dir / "prepare.json").is_file(),
            "failure_trail": (case_dir / "station-failure.json").exists(),
            "health": "alive, steady CPU accrual",
        },
        "armed_pipeline": [
            "nsc-t30-extinction reducer (fires on station files)",
            "nsc-t30-expectations verdict map (pre-registered)",
            "nsc-extinction-extension merge (refuses before landing)",
            "conditional test (skips until the record lands)",
            "nsc-t30-adoption-checklist (seven steps, three guardrails)",
        ],
        "continuation_for_a_later_session": [
            "check the campaign directory for station-25/30 files",
            "run scripts/derive_nsc_t30_extinction.py once",
            "run scripts/derive_nsc_extinction_extension.py once",
            "follow the adoption checklist record",
            "run the suites: the conditional test then passes instead of skipping",
        ],
        "notes": [
            "all corridor records are creation-only; none may be edited",
            "the launch log, forecasts, watch and checklist hold the "
            "wait's whole story",
        ],
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_t30_standing_state.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_t30_standing_state.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    print("  march healthy trail:", record["march"]["prepare_trail_intact"], "|",
          "failure:", record["march"]["failure_trail"])
    print("  armed stages:", len(record["armed_pipeline"]))
    print("  continuation steps:", len(record["continuation_for_a_later_session"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())