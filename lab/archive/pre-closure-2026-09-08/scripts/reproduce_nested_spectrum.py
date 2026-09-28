#!/usr/bin/env python3
"""Regenerate the bounded NGS-0 scale, recurrence, and probability controls."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

from recursive_horizons.nested_spectrum import (  # noqa: E402
    eventual_occurrence_probability,
    finite_band_selection,
    required_nonclaims,
    two_level_recurrence_ratio,
)


def record() -> dict[str, object]:
    return {
        "schema_version": 1,
        "project_version": "0.11.0",
        "model_id": "NGS-0",
        "artifact": "nested_gradient_spectral_recurrence_probability_controls",
        "classification": "conditional_mathematical_controls_not_nested_cosmology",
        "finite_band_selection": finite_band_selection(alpha=2.0, beta=1.0),
        "recurrence_controls": {
            "golden_fixture": two_level_recurrence_ratio(a=1.0, b=1.0),
            "non_golden_fixture": two_level_recurrence_ratio(a=2.0, b=1.0),
        },
        "eternal_opportunity_control": eventual_occurrence_probability(
            probability_per_trial=1.0e-6,
            trial_count=1_000_000,
        ),
        "nonclaims": required_nonclaims(),
        "scope": (
            "NGS-0 executes a spatial finite-band selector, a conditional two-level "
            "recurrence eigenvalue, and a finite independent-trial probability. It "
            "does not derive a covariant FGC action, nested spacetime, dark sector, "
            "universal hertz, cosmological golden ratio, inevitable observers, or plan."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=REPOSITORY / "results" / "nested-gradient-spectrum.json",
    )
    output = parser.parse_args().output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    payload = json.dumps(record(), indent=2, sort_keys=True, allow_nan=False) + "\n"
    temporary = output.with_suffix(output.suffix + ".tmp")
    temporary.write_text(payload, encoding="utf-8")
    temporary.replace(output)
    print(f"wrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
