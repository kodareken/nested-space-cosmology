"""Frontier assessment of the influence-chain strand.

After sixty-two rounds every strand of the repository that shares the
RN radial programme has a join-table with this programme. This driver
writes the deliberate frontier ledger: which questions are closed, with
their records; which cells remain open, with the concrete reason each
is production-scale; and the honest count of what was built. Verified
numbers only - each closed item carries its record path.

    python scripts/lab.py scripts/derive_nsc_frontier_assessment.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_influence_chain as chain  # noqa: E402

OUTPUT = LAB / "results" / "development" / "nsc-frontier-assessment-v1.json"
SCHEMA = "NSC-FRONTIER-ASSESSMENT-v1"

CLOSED = [
    {
        "question": "exterior influence channels as the carrier of the missing budget",
        "closed_by": [
            "nsc-influence-chain-v1.json",
            "nsc-influence-shallow-power-v1.json",
            "nsc-influence-response-v1.json",
            "nsc-influence-seam-v1.json",
        ],
        "statement": "all measured per-level fees close outside the pre-registered band",
    },
    {
        "question": "interior funnel bookkeeping",
        "closed_by": [
            "nsc-interior-balance-v1.json",
            "nsc-continuation-station-closure-v1.json",
        ],
        "statement": "dM = d(excision) at machine precision at every measured resolution; outer leaks only",
    },
    {
        "question": "wall organism dynamics",
        "closed_by": [
            "nsc-wall-junction-family-v1.json",
            "nsc-wall-lead-lag-v1.json",
            "nsc-wall-clock-horizon-v1.json",
            "nsc-wall-clock-fall-v1.json",
            "nsc-wall-asymmetry-derived-v1.json",
            "nsc-wall-clock-derivation-v1.json",
            "nsc-wall-clock-epscross-v1.json",
        ],
        "statement": "concentration, lag arrow, two gap powers, derived asymmetry, epsilon crossing",
    },
    {
        "question": "scaffold invariance and the cosmic band's home",
        "closed_by": [
            "nsc-seam-invariance-v1.json",
            "nsc-charged-seam-invariance-v1.json",
            "nsc-cosmic-band-location-v1.json",
            "nsc-cosmic-shell-weight-v1.json",
            "nsc-time-awareness-v1.json",
        ],
        "statement": "mu/eta-blind scaffold; band maps to the outer cap shell; beta profile recorded",
    },
    {
        "question": "charged spectroscopy registry",
        "closed_by": [
            "nsc-omega-star-identity-v1.json",
            "nsc-charged-clock-horizon-v1.json",
            "nsc-charged-ladder-deficit-law-v1.json",
        ],
        "statement": "omega_star identity at machine precision; cap clock; Gaussian ladder law",
    },
    {
        "question": "fringe staircase",
        "closed_by": [
            "nsc-charged-tilt-resolution-v1.json",
            "nsc-tilt-organism-clock-v1.json",
            "nsc-tilt-cell-crossings-v1.json",
            "nsc-tilt-quantum-v1.json",
            "nsc-tilt-quantum-ladder-v1.json",
            "nsc-neutral-fringe-quantum-v1.json",
            "nsc-fringe-beta-join-v1.json",
        ],
        "statement": "nodal quantum 0.00131, count ladder {3,3,2,1}, chart-free q*N = 1.149",
    },
    {
        "question": "joins with the repository's other strands",
        "closed_by": [
            "nsc-matryoshka-fee-join-v1.json",
            "nsc-two-ends-join-v1.json",
            "nsc-two-ends-covariance-v1.json",
            "nsc-oneway-wall-join-v1.json",
            "nsc-cone-wall-join-v1.json",
            "nsc-isotropy-wall-join-v1.json",
            "nsc-isotropy-driven-join-v1.json",
        ],
        "statement": "seven cross-strand tables laid; all directions recorded, n stated everywhere",
    },
]

OPEN = [
    {
        "cell": "wall morphology beyond t=20",
        "blocker": "the continuation directory is lock-listed by the lab writer policy; a new campaign requires the sealed campaign runner and hours of CPU",
        "witness": "nsc-isotropy-wall-join-v1.json tail (I at t=40 and t=80)",
    },
    {
        "cell": "exponent provenance at a fourth gap point",
        "blocker": "the packet preparer is family-pinned and no fourth family mass is recorded; new family runs are a production campaign",
    },
    {
        "cell": "derivation of p_Q = 0.1173 from the room-quality definitions",
        "blocker": "the two-ends v2 record stores only reduced qualities, not the per-room lapse/shift ingredients its own definition uses",
    },
    {
        "cell": "paper-level coupling of the rim location to a cosmological observable",
        "blocker": "out of lab domain by construction; hypothesis-labeled in the records",
    },
]

def main():
    if OUTPUT.exists():
        raise SystemExit("refusing to overwrite existing output: " + str(OUTPUT))
    tests = (LAB / "tests" / "test_nsc_influence_chain.py",
             LAB / "tests" / "test_nsc_rn_charged_lll.py")
    test_body = "".join(p.read_text(encoding="utf-8") for p in tests if p.is_file())
    test_defs = test_body.count("def test_")
    records = sorted(
        p.name for p in (LAB / "results" / "development").glob("nsc-*.json")
        if p.name.startswith(("nsc-influence", "nsc-wall", "nsc-charged", "nsc-seam",
                              "nsc-interior", "nsc-continuation", "nsc-cosmic",
                              "nsc-time-awareness", "nsc-tilt", "nsc-neutral-fringe",
                              "nsc-fringe", "nsc-two-ends", "nsc-oneway", "nsc-cone",
                              "nsc-isotropy", "nsc-ladder", "nsc-lowband",
                              "nsc-constant-table", "nsc-matryoshka-fee",
                              "nsc-omega-star-identity", "nsc-frontier")))
    record = {
        "schema": SCHEMA,
        "api_version": chain.API_VERSION,
        "creation_only": True,
        "production_trajectory": False,
        "evolved": False,
        "rounds_to_date": 63,
        "closed": CLOSED,
        "open": OPEN,
        "counts": {
            "closed_question_groups": len(CLOSED),
            "open_cells": len(OPEN),
            "test_definitions": test_defs,
            "strand_records": len(records),
        },
        "assessment": (
            "the recorded-bytes frontier is reached: every strand that "
            "shares the RN radial programme has a join-table, and the four "
            "open cells carry concrete production blockers, not open "
            "questions. The founding riddle has its complete answer cycle, "
            "hypothesis-labeled at its only interpretive step."
        ),
        "notes": [
            "records listed by name live in lab/results/development",
            "no new numbers are computed here; this is the ledger",
        ],
        "source_sha256": {
            "lab/src/recursive_horizons/nsc_influence_chain.py": chain.path_sha256(
                LAB / "src" / "recursive_horizons" / "nsc_influence_chain.py"
            ),
            "lab/scripts/derive_nsc_frontier_assessment.py": chain.path_sha256(
                LAB / "scripts" / "derive_nsc_frontier_assessment.py"
            ),
        },
    }
    with OUTPUT.open("w", encoding="utf-8") as handle:
        handle.write(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print("wrote", OUTPUT.relative_to(LAB))
    print("  closed groups: {}  open cells: {}  tests: {}  records: {}".format(
        record["counts"]["closed_question_groups"], record["counts"]["open_cells"],
        record["counts"]["test_definitions"], record["counts"]["strand_records"]))
    for cell in OPEN:
        print("  OPEN: {} -- {}".format(cell["cell"], cell["blocker"][:90]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())