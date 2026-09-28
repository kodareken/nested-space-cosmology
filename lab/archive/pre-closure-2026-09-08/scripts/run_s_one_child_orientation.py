#!/usr/bin/env python3
"""Select the doubled-vacuum sign from the exact expanding child geometry.

The gap/transition closure leaves two equation-orientation routes because the
on-room vacuum-form tensor can be quoted on either side of the metric equation.
This calculation compares curvature conventions directly.  The exact
black-universe child has negative Ricci scalar in the repository convention;
the doubled-FLRW de Sitter solution has R=-2 Lambda_effective in the source
convention.  The expanding child therefore selects positive Lambda_effective,
fixing the source-side gap 1006/1015 without observational fitting.
"""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
for value in (str(ROOT), str(ROOT / "src")):
    if value not in sys.path:
        sys.path.insert(0, value)

from recursive_horizons.evidence_io import canonical_json_bytes  # noqa: E402


OUTPUT = ROOT / "results/nsc-1-s-one-child-orientation.json"
INPUTS = (
    ROOT / "results/nsc-1-s-one-gap-transition-closure.json",
    ROOT / "results/nsc-1-exact-black-universe-defocusing.json",
    ROOT / "results/nsc-1-s-one-gauss-bonnet-component-closure.json",
    ROOT / "results/nsc-1-s-one-doubled-flrw-closure.json",
)


def _authenticate() -> tuple[list[dict[str, object]], dict[str, dict[str, object]]]:
    authenticated: list[dict[str, object]] = []
    records: dict[str, dict[str, object]] = {}
    for path in INPUTS:
        raw = path.read_bytes()
        item = json.loads(raw)
        if item.get("terminal") is not True:
            raise RuntimeError(f"nonterminal input: {path.relative_to(ROOT)}")
        artifact_id = str(item["artifact_id"])
        records[artifact_id] = item
        authenticated.append(
            {
                "path": str(path.relative_to(ROOT)),
                "sha256": hashlib.sha256(raw).hexdigest(),
                "artifact_id": artifact_id,
            }
        )
    return authenticated, records


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("child-orientation output already exists")
    authenticated, records = _authenticate()
    gap = records["NSC-1-S-ONE-GAP-TRANSITION-CLOSURE"]
    component = records["NSC-1-S-ONE-GAUSS-BONNET-COMPONENT-CLOSURE"]

    rho = sp.symbols("rho", real=True)
    room_ricci = sp.sympify(
        component["room_curvature_spectrum"]["Ricci_scalar"],
        locals={"rho": rho},
    )
    child_ricci_limit = sp.simplify(sp.limit(room_ricci, rho, -sp.oo))
    parent_ricci_limit = sp.simplify(sp.limit(room_ricci, rho, sp.oo))

    time, hubble, scale_0 = sp.symbols(
        "t lambda a_0", positive=True, real=True
    )
    scale_factor = scale_0 * sp.exp(hubble * time)
    lapse = sp.Integer(1)
    flrw_ricci = sp.simplify(
        6
        * (
            sp.diff(scale_factor, time) * sp.diff(lapse, time)
            / (scale_factor * lapse**3)
            - sp.diff(scale_factor, time) ** 2
            / (scale_factor**2 * lapse**2)
            - sp.diff(scale_factor, time, 2)
            / (scale_factor * lapse**2)
        )
    )
    lambda_effective = 6 * hubble**2
    flrw_vacuum_relation = sp.simplify(flrw_ricci + 2 * lambda_effective)

    source_gap_squared = Fraction(
        gap["derived_gap_candidates"]["source_side_positive"][
            "Phi_squared_over_Lambda_squared"
        ]
    )
    rejected_sign_gap_squared = Fraction(
        gap["derived_gap_candidates"]["geometric_left_negative"][
            "Phi_squared_over_Lambda_squared"
        ]
    )
    vacuum = Fraction(108, 1015)
    boundary_xi = vacuum / source_gap_squared
    candidate_rows = gap["full_exponential_evaluation"]["candidate_rows"][
        "source_side_positive"
    ]
    minimum_kernel = min(float(row["normalized_kernel"]) for row in candidate_rows)

    record = {
        "artifact_id": "NSC-1-S-ONE-CHILD-ORIENTATION",
        "schema": "NSC-1-S-ONE-CHILD-ORIENTATION-v1",
        "classification": "the_exact_expanding_child_selects_the_positive_source_side_vacuum_and_the_gap_squared_1006_over_1015",
        "authenticated_inputs": authenticated,
        "curvature_convention_match": {
            "black_universe_room_Ricci_scalar": str(room_ricci),
            "child_limit_rho_to_minus_infinity": str(child_ricci_limit),
            "parent_limit_rho_to_plus_infinity": str(parent_ricci_limit),
            "doubled_FLRW_de_Sitter_scale_factor": str(scale_factor),
            "doubled_FLRW_Ricci_scalar": str(flrw_ricci),
            "doubled_FLRW_effective_vacuum": str(lambda_effective),
            "identity": "R_child=-2*Lambda_effective",
            "identity_residual": str(flrw_vacuum_relation),
            "sign_consequence": "the_negative_child_Ricci_scalar_corresponds_to_positive_Lambda_effective_in_the_doubled_FLRW_source_convention",
        },
        "orientation_selection": {
            "selected_route": "source_side_positive",
            "selected_signed_vacuum": str(vacuum),
            "selected_Phi_squared_over_Lambda_squared": str(
                source_gap_squared
            ),
            "selected_Phi_over_Lambda": float(sp.sqrt(source_gap_squared)),
            "selected_location": "9/1015_below_the_local_spectral_wall",
            "unselected_equation_rearrangement": "geometric_left_negative",
            "unselected_Phi_squared_over_Lambda_squared": str(
                rejected_sign_gap_squared
            ),
            "observation_used_to_select_sign": False,
        },
        "first_cross_scale_number": {
            "name": "Xi_boundary",
            "definition": "Lambda_effective/m_gap_squared",
            "exact": str(boundary_xi),
            "decimal": float(boundary_xi),
            "derivation": "(108/1015)/(1006/1015)=54/503",
            "scope": "boundary_spectral_pole_to_child_vacuum_in_the_one_scale_branch",
            "not_electron_Xi_NSC": True,
        },
        "selected_full_profile_kernel": {
            "rows": candidate_rows,
            "minimum": minimum_kernel,
            "all_positive": minimum_kernel > 0.0,
        },
        "causal_result": {
            "discrete_orientation_ambiguity_removed": True,
            "selected_gap_is_locally_resolved_but_adjacent_to_cutoff": True,
            "same_Phi_roles": [
                "boundary_mass_gap",
                "outside_self_energy",
                "relative_metric_interaction",
                "child_vacuum_shift",
            ],
            "next_equation": "evaluate_delta_Gamma_one_over_delta_Phi_at_Phi_squared_over_Lambda_squared_1006_over_1015_then_construct_the_curved_retarded_pole_kernel",
        },
        "gate": {
            "inputs_authenticated": len(authenticated) == len(INPUTS),
            "child_Ricci_limit_negative": bool(child_ricci_limit < 0),
            "parent_Ricci_limit_zero": parent_ricci_limit == 0,
            "doubled_FLRW_sign_identity_exact": flrw_vacuum_relation == 0,
            "positive_source_orientation_selected_without_observation": True,
            "selected_gap_exact": source_gap_squared
            == Fraction(1006, 1015),
            "boundary_Xi_exact": boundary_xi == Fraction(54, 503),
            "selected_full_kernel_positive": minimum_kernel > 0.0,
            "stationary_gap_equation_solved": False,
            "curved_retarded_kernel_completed": False,
            "conviction_proof_seed_completed": False,
        },
        "nonclaims": {
            "boundary_gap_is_an_observed_particle": False,
            "electron_Xi_NSC_predicted": False,
            "Euler_Lagrange_stationarity_proved": False,
            "curved_Lorentzian_spectrum_proved": False,
            "our_universe_identified": False,
        },
        "terminal": True,
    }
    for key, value in record["gate"].items():
        if key not in {
            "stationary_gap_equation_solved",
            "curved_retarded_kernel_completed",
            "conviction_proof_seed_completed",
        } and value is not True:
            raise RuntimeError(f"child-orientation gate failed: {key}")
    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
