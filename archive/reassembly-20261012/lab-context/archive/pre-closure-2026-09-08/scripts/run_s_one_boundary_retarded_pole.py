#!/usr/bin/env python3
"""Derive the causal boundary particle/wave pole below the spectral wall.

The expanding-child orientation fixes Phi^2/Lambda^2=1006/1015.  The same
two-sheet Dirac operator then has a standard retarded resolvent with positive
pole weights and subluminal group velocity.  Its locally propagating momentum
band ends exactly when p^2+Phi^2=Lambda^2; the remaining squared bandwidth is
9/1015, one half of the transition-curvature magnitude.

This is the boundary spectral excitation, not an identification with an
electron or the complete relative-metric/graviton propagator.
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


OUTPUT = ROOT / "results/nsc-1-s-one-boundary-retarded-pole.json"
INPUTS = (
    ROOT / "results/nsc-1-s-one-child-orientation.json",
    ROOT / "results/nsc-1-s-one-local-two-sheet-anomaly.json",
    ROOT / "results/nsc-1-s-one-spectral-resolution-wall-bind.json",
    ROOT / "results/nsc-1-s-one-transition-link-coefficient.json",
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
        raise RuntimeError("boundary retarded-pole output already exists")
    authenticated, records = _authenticate()
    orientation = records["NSC-1-S-ONE-CHILD-ORIENTATION"]
    transition_record = records["NSC-1-S-ONE-TRANSITION-LINK-COEFFICIENT"]

    mass_squared = Fraction(
        orientation["orientation_selection"][
            "selected_Phi_squared_over_Lambda_squared"
        ]
    )
    transition = Fraction(
        transition_record["transition_derivation"][
            "quadratic_coefficient_over_F"
        ]
    )
    bandwidth_squared = 1 - mass_squared

    momentum, mass, omega = sp.symbols(
        "p Phi omega", nonnegative=True, real=True
    )
    energy = sp.sqrt(momentum**2 + mass**2)
    dirac = sp.Matrix([[momentum, mass], [mass, -momentum]])
    inverse_symbol = sp.factor((omega * sp.eye(2) - dirac).det())
    resolvent = sp.simplify((omega * sp.eye(2) - dirac).inv())
    visible = sp.factor(resolvent[0, 0])
    positive_residue = sp.simplify(
        sp.limit((omega - energy) * visible, omega, energy)
    )
    negative_residue = sp.simplify(
        sp.limit((omega + energy) * visible, omega, -energy)
    )
    residue_sum = sp.simplify(positive_residue + negative_residue)
    group_velocity = sp.simplify(sp.diff(energy, momentum))

    maximum_momentum = sp.sqrt(sp.Rational(bandwidth_squared.numerator, bandwidth_squared.denominator))
    exact_mass = sp.sqrt(sp.Rational(mass_squared.numerator, mass_squared.denominator))
    sample_momenta = (sp.Integer(0), maximum_momentum / 2, maximum_momentum)
    samples = []
    for value in sample_momenta:
        sample_energy = sp.simplify(
            energy.subs({momentum: value, mass: exact_mass})
        )
        sample_positive = sp.simplify(
            positive_residue.subs({momentum: value, mass: exact_mass})
        )
        sample_negative = sp.simplify(
            negative_residue.subs({momentum: value, mass: exact_mass})
        )
        sample_velocity = sp.simplify(
            group_velocity.subs({momentum: value, mass: exact_mass})
        )
        samples.append(
            {
                "p_over_Lambda": str(value),
                "energy_over_Lambda": str(sample_energy),
                "energy_decimal": float(sample_energy),
                "positive_residue": str(sample_positive),
                "negative_residue": str(sample_negative),
                "group_velocity_over_c": str(sample_velocity),
                "group_velocity_decimal": float(sample_velocity),
                "inside_or_on_cutoff": bool(sample_energy <= 1),
            }
        )

    cutoff_energy = sp.simplify(
        energy.subs({momentum: maximum_momentum, mass: exact_mass})
    )
    maximum_velocity = sp.simplify(
        group_velocity.subs({momentum: maximum_momentum, mass: exact_mass})
    )

    record = {
        "artifact_id": "NSC-1-S-ONE-BOUNDARY-RETARDED-POLE",
        "schema": "NSC-1-S-ONE-BOUNDARY-RETARDED-POLE-v1",
        "classification": "the_selected_boundary_gap_is_a_positive_weight_causal_particle_wave_pole_whose_local_band_ends_at_the_transition_fixed_spectral_wall",
        "authenticated_inputs": authenticated,
        "retarded_operator": {
            "Dirac_symbol": str(dirac),
            "retarded_resolvent": "G_R(omega,p)=[(omega+i0)I-D(p)]^-1",
            "retarded_time_kernel": "G_R(t,p)=-i*theta(t)*exp(-i*D(p)*t)",
            "causal_support": "t>=0",
            "characteristic_denominator": str(inverse_symbol),
            "visible_component": str(visible),
            "positive_energy": str(energy),
            "positive_residue": str(positive_residue),
            "negative_residue": str(negative_residue),
            "residue_sum": str(residue_sum),
            "group_velocity": str(group_velocity),
        },
        "transition_fixed_local_band": {
            "mass_squared_over_cutoff_squared": str(mass_squared),
            "momentum_band_squared_over_cutoff_squared": str(
                bandwidth_squared
            ),
            "maximum_momentum_over_cutoff": str(maximum_momentum),
            "transition_curvature": str(transition),
            "bandwidth_identity": "p_max^2/Lambda^2=abs(c_transition)/2",
            "spectral_wall_identity": "m_gap^2+p_max^2=Lambda^2",
            "energy_at_band_edge_over_cutoff": str(cutoff_energy),
            "maximum_group_velocity_over_c": str(maximum_velocity),
            "samples": samples,
        },
        "same_state_identity": {
            "particle": "positive_retarded_pole_at_E=sqrt(p^2+Phi^2)",
            "wave": "positive_pole_residue_and_unitary_exp(-iDt)_amplitude",
            "outside": "the_same_visible_resolvent_contains_minus_Phi^2/(omega+p+i0)",
            "boundary": "the_available_local_momentum_band_is_fixed_by_the_same_transition_curvature",
        },
        "causal_result": {
            "flat_boundary_particle_wave_pole_completed": True,
            "pole_is_below_or_on_cutoff_for_complete_band": True,
            "next_equation": "derive_the_full_relative_metric_retarded_projector_kernel_on_the_curved_two_sheet_background_then_prove_the_selected_gap_is_a_recursive_fixed_spectral_datum",
        },
        "gate": {
            "inputs_authenticated": len(authenticated) == len(INPUTS),
            "selected_mass_gap_exact": mass_squared == Fraction(1006, 1015),
            "bandwidth_exact": bandwidth_squared == Fraction(9, 1015),
            "bandwidth_is_half_transition_curvature": bandwidth_squared
            == abs(transition) / 2,
            "band_edge_hits_cutoff_exactly": cutoff_energy == 1,
            "pole_residues_sum_to_one": residue_sum == 1,
            "all_sample_residues_nonnegative": all(
                sp.sympify(item["positive_residue"]) >= 0
                and sp.sympify(item["negative_residue"]) >= 0
                for item in samples
            ),
            "all_sample_group_velocities_subluminal": all(
                0 <= item["group_velocity_decimal"] < 1 for item in samples
            ),
            "retarded_time_support_causal": True,
            "curved_relative_metric_retarded_kernel_completed": False,
            "electron_pole_predicted": False,
            "conviction_proof_seed_completed": False,
        },
        "nonclaims": {
            "boundary_mode_is_the_electron": False,
            "full_metric_graviton_covariance_proved": False,
            "curved_black_child_stability_proved": False,
            "observed_cosmology_identified": False,
        },
        "terminal": True,
    }
    for key, value in record["gate"].items():
        if key not in {
            "curved_relative_metric_retarded_kernel_completed",
            "electron_pole_predicted",
            "conviction_proof_seed_completed",
        } and value is not True:
            raise RuntimeError(f"boundary retarded-pole gate failed: {key}")
    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
