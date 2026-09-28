#!/usr/bin/env python3
"""Regenerate the bounded EC-1 homogeneous spin-fluid artifact."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

from recursive_horizons.einstein_cartan import (  # noqa: E402
    EinsteinCartanSpec,
    boundary_areal_radius,
    boundary_misner_sharp_mass,
    boundary_work_residual,
    causal_sequence,
    continuity_residual,
    curvature_invariants,
    fluid_state,
    friedmann_F,
    half_cycle_convergence,
    scale_factor_acceleration,
    static_parent_mass_drift,
    turning_points,
)


def record() -> dict[str, object]:
    spec = EinsteinCartanSpec(A=10.0, B=1.0)
    roots = turning_points(spec)
    a_min = float(roots["a_min"])
    a_max = float(roots["a_max"])
    probe = (a_min + a_max) / 2.0
    sequence = causal_sequence(spec)
    return {
        "schema_version": 1,
        "project_version": "0.11.0",
        "model_id": "EC-1",
        "artifact": "einstein_cartan_homogeneous_spin_fluid_cap",
        "classification": "reduced_action_motivated_homogeneous_effective_fluid_benchmark",
        "scope": "Closed-FLRW random-spin ultrarelativistic reduction only; not full Dirac/tetrad evolution, a resolved vacuum boundary, or GMF-1B.",
        "parameters": {
            "A": spec.A,
            "B": spec.B,
            "chi_boundary": spec.chi_boundary,
            "gravitational_constant": spec.gravitational_constant,
            "kappa": spec.kappa,
        },
        "turning_points": roots,
        "bounce": {
            "friedmann_F": friedmann_F(spec, a_min),
            "acceleration": scale_factor_acceleration(spec, a_min),
            "fluid": fluid_state(spec, a_min),
            "curvature_invariants": curvature_invariants(spec, a_min),
            "continuity_residual_expansion": continuity_residual(spec, a_min, "expansion"),
        },
        "interior_probe": {
            "scale_factor": probe,
            "friedmann_F": friedmann_F(spec, probe),
            "fluid": fluid_state(spec, probe),
            "curvature_invariants": curvature_invariants(spec, probe),
            "continuity_residual_collapse": continuity_residual(spec, probe, "collapse"),
            "continuity_residual_expansion": continuity_residual(spec, probe, "expansion"),
        },
        "null_expansion_sequence": sequence,
        "proper_time": half_cycle_convergence(spec),
        "boundary": {
            "areal_radius_at_bounce": boundary_areal_radius(spec, a_min),
            "areal_radius_at_maximum": boundary_areal_radius(spec, a_max),
            "misner_sharp_mass_at_bounce": boundary_misner_sharp_mass(spec, a_min),
            "misner_sharp_mass_at_maximum": boundary_misner_sharp_mass(spec, a_max),
            "work_residual_collapse": boundary_work_residual(spec, probe, "collapse"),
            "work_residual_expansion": boundary_work_residual(spec, probe, "expansion"),
            "static_parent_mass_drift": static_parent_mass_drift(spec),
            "fixed_parent_crossing_comparator_only": True,
        },
        "components": {
            "radiation": "w=1/3, rho proportional to a^-4",
            "negative_torsion": "w=1, rho proportional to -a^-6",
            "q_internal": 0.0,
            "late_w_minus_one_component": False,
        },
        "global_solution_constructed": False,
        "smooth_static_vacuum_boundary_proven": False,
        "child_topology_proven": False,
        "derived_external_Q": False,
        "late_dark_energy_derived": False,
        "variable_c_derived": False,
        "radial_stability_proven": False,
        "microscopic_dirac_stability_or_eft_control_proven": False,
        "homogeneous_spin_fluid_closure_assumed": True,
        "single_interior_metric_causal_sequence": sequence[
            "single_interior_metric_causal_sequence"
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=REPOSITORY / "results" / "einstein-cartan-collapse.json",
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
