#!/usr/bin/env python3
"""Regenerate the GMF-1B-PF1 focusing/thermal-closure preflight artifact."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

from recursive_horizons.gmf1b_preflight import (  # noqa: E402
    EC1BounceSpec,
    canonical_scalar_focusing_gate,
    ec1_bounce_thermal_gate,
    thermal_ec_state,
)


def record() -> dict[str, object]:
    canonical = canonical_scalar_focusing_gate(
        gravitational_constant=1.0,
        initial_expansion=-1.0,
        shear_squared=0.25,
        k_dot_phi=0.5,
    )
    ec_gate = ec1_bounce_thermal_gate(EC1BounceSpec(A=10.0, B=1.0))
    return {
        "schema_version": 1,
        "project_version": "0.11.0",
        "model_id": "GMF-1B-PF1",
        "artifact": "gmf1b_focusing_and_thermal_spin_fluid_preflight",
        "classification": "conditional_material_law_and_null_focusing_rejection_gates",
        "canonical_scalar_null_focusing": canonical,
        "thermal_ec_closure_regimes": {
            "causal_example": thermal_ec_state(0.1),
            "degenerate_example": thermal_ec_state(2.0 / 9.0),
            "gradient_unstable_example": thermal_ec_state(0.4),
            "singular_example": thermal_ec_state(2.0 / 3.0),
            "superluminal_example": thermal_ec_state(0.8),
        },
        "ec1_bounce_gate": ec_gate,
        "scope": "Rejects only the naive thermal averaged perfect-fluid closure as GMF-1B material law and the stated canonical-scalar null-defocusing route; not full Einstein-Cartan-Dirac, EC-1 background algebra, noncanonical/modified gravity, or topology-changing physics.",
        "global_solution_constructed": False,
        "general_recursive_horizons_no_go_proven": False,
        "canonical_scalar_regular_bridge_supported_under_stated_assumptions": False,
        "full_einstein_cartan_dirac_rejected": False,
        "ec1_background_algebra_rejected": False,
        "external_Q_derived": False,
        "late_dark_energy_derived": False,
        "variable_c_derived": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=REPOSITORY / "results" / "gmf-1b-preflight.json")
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
