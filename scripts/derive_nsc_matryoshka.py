"""Write the Matryoshka comparison record.

Enumerative and recursion-structure comparison between the cosmohedron
Matryoshkas (arXiv:2603.03425) and the NSC chain nesting. Creation only;
the record is refused if it already exists. No external compute.

Run from the repository root:

    python scripts/lab.py scripts/derive_nsc_matryoshka.py
"""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

LAB = Path(__file__).resolve().parents[1]
if str(LAB / "src") not in sys.path:
    sys.path.insert(0, str(LAB / "src"))

from recursive_horizons import nsc_matryoshka as matryoshka  # noqa: E402

RECORD = LAB / "results/development/nsc-matryoshka-comparison-v1.json"


def _producing_commit() -> str | None:
    try:
        out = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=LAB, capture_output=True, text=True, check=True,
        )
        return out.stdout.strip()
    except Exception:
        return None


def _sha256_of_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def main() -> int:
    if RECORD.exists():
        print(f"refusing to overwrite sealed record: {RECORD}", file=sys.stderr)
        return 1

    body = matryoshka.comparison()
    record = {
        "schema": "NSC-MATRYOSHKA-v1",
        "api_version": matryoshka.API_VERSION,
        "creation_only": True,
        "producing_commit": _producing_commit(),
        "source_paper": "arXiv:2603.03425 (Ardila-Mantilla, Arkani-Hamed, Figueiredo, Vazao)",
        "comparison": body,
    }

    RECORD.write_text(json.dumps(record, indent=1, sort_keys=True) + "\n", encoding="utf-8")
    print(f"wrote {RECORD.relative_to(LAB)}")
    print(f"  matryoshka prefix matches paper  : {body['matryoshka_prefix_matches_paper']}")
    print(f"  schroder prefix matches reference: {body['schroder_prefix_matches_reference']}")
    print(f"  matryoshka growth ratios         : {body['matryoshka_growth_ratios']}")
    print(f"  chain growth ratios              : {body['chain_count_growth_ratios']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())