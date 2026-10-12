"""Pre-registered expectations for the t30 extinction verdict.

Written before the march lands, this record maps each of the three
possible verdicts to the ledger changes it will mandate, so the
adoption step that follows the landing is mechanical and its framing
was fixed in advance (no post-hoc reading).

    python scripts/lab.py scripts/derive_nsc_t30_expectations.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402

OUTPUT = LAB / "results" / "development" / "nsc-t30-expectations-v1.json"
SCHEMA = "NSC-T30-EXPECTATIONS-v1"


def main():
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    landing = {
        "DIES_IN_20_25": (
            "the light wall's first death window is pinned between the march "
            "v1 coverage and t=25: the extinction ledger gains its light-side "
            "endpoint, the lifespan's mu-ordering gains its light pole, and "
            "the two-domain comparison becomes complete"
        ),
        "DIES_IN_25_30": (
            "the light wall survives past t=25 but dies before t=30: lifespan "
            "exceeds all march-v1 coverage and the death window tightens to "
            "a five-unit band"
        ),
        "WALL_SURVIVES_TO_30": (
            "the light wall lives at least to t=30: the light-side lifespan "
            "extends beyond every recorded era, the deep-domain amplitude "
            "law's quiet breathing gains its longest recorded member, and a "
            "longer campaign becomes the named next question"
        ),
    }
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "pre_registered_before_landing": True,
        "campaign": "nsc-rn-inward-march-t30-v1",
        "verdict_map": landing,
        "allowed_verdicts": sorted(landing),
        "adoption_steps_when_landed": [
            "run scripts/derive_nsc_t30_extinction.py",
            "read the verdict from nsc-t30-extinction-v1.json",
            "summarise per this pre-registration; no new framing",
        ],
        "notes": [
            "written at round 98 while the march integrating continues",
        ],
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_t30_expectations.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_t30_expectations.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    print("  registered verdict map for:", record["allowed_verdicts"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())