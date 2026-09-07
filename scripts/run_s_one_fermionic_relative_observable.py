#!/usr/bin/env python3
"""Construct the positive fermionic observable for relative sheet resolution.

The inverse heat Hessian is not the physical covariance.  In the invariant
fermionic partition function, a relative sheet source couples to the Hermitian
sheet-resolution generator sigma_3.  Its connected two-point function has an
exact Lehmann representation.  This calculation derives its nonnegative
spectral weight and proves Osterwalder--Schrader positivity for the controlled
two-sheet mode at all momenta.

For the child-selected gap, the one-particle boundary pole remains just below
the cutoff while the collective sheet-changing excitation has energy 2E and
therefore lies beyond the local spectral wall.
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


OUTPUT = ROOT / "results/nsc-1-s-one-fermionic-relative-observable.json"
INPUTS = (
    ROOT / "results/nsc-1-s-one-relative-reflection-test.json",
    ROOT / "results/nsc-1-s-one-boundary-retarded-pole.json",
    ROOT / "results/nsc-1-s-one-invariant-partition-bind.json",
    ROOT / "results/nsc-1-s-one-fermionic-geometry-bind.json",
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
        raise RuntimeError("fermionic relative-observable output already exists")
    authenticated, records = _authenticate()
    boundary = records["NSC-1-S-ONE-BOUNDARY-RETARDED-POLE"]

    momentum, mass, tau, frequency = sp.symbols(
        "p Phi tau omega", positive=True, real=True
    )
    energy = sp.sqrt(momentum**2 + mass**2)
    hamiltonian = sp.Matrix([[momentum, mass], [mass, -momentum]])
    sheet_generator = sp.diag(1, -1)
    projector_plus = sp.simplify((sp.eye(2) + hamiltonian / energy) / 2)
    projector_minus = sp.simplify((sp.eye(2) - hamiltonian / energy) / 2)
    spectral_weight = sp.factor(
        sp.trace(
            projector_minus
            * sheet_generator
            * projector_plus
            * sheet_generator
        )
    )
    expected_weight = mass**2 / energy**2
    euclidean_correlator = sp.simplify(
        spectral_weight * sp.exp(-2 * energy * tau)
    )
    retarded_susceptibility = sp.factor(
        spectral_weight
        * (
            1 / (frequency - 2 * energy)
            - 1 / (frequency + 2 * energy)
        )
    )
    expected_retarded = sp.factor(
        4 * energy * spectral_weight / (frequency**2 - 4 * energy**2)
    )

    # For arbitrary positive Euclidean times the OS matrix is a positive outer
    # product at each momentum.  An integral/sum of these matrices remains PSD.
    time_symbols = sp.symbols("t_1:5", positive=True, real=True)
    os_vector = sp.Matrix([sp.exp(-2 * energy * value) for value in time_symbols])
    os_matrix = sp.simplify(spectral_weight * os_vector * os_vector.T)
    two_by_two_principal_minor = sp.simplify(os_matrix[:2, :2].det())

    selected_mass_squared = Fraction(
        boundary["transition_fixed_local_band"][
            "mass_squared_over_cutoff_squared"
        ]
    )
    selected_bandwidth_squared = Fraction(
        boundary["transition_fixed_local_band"][
            "momentum_band_squared_over_cutoff_squared"
        ]
    )
    selected_mass = sp.sqrt(
        sp.Rational(selected_mass_squared.numerator, selected_mass_squared.denominator)
    )
    selected_pmax = sp.sqrt(
        sp.Rational(
            selected_bandwidth_squared.numerator,
            selected_bandwidth_squared.denominator,
        )
    )
    collective_minimum = sp.simplify(2 * selected_mass)
    collective_maximum = sp.simplify(
        2 * energy.subs({momentum: selected_pmax, mass: selected_mass})
    )
    rest_weight = sp.simplify(
        spectral_weight.subs({momentum: 0, mass: selected_mass})
    )
    edge_weight = sp.simplify(
        spectral_weight.subs({momentum: selected_pmax, mass: selected_mass})
    )

    record = {
        "artifact_id": "NSC-1-S-ONE-FERMIONIC-RELATIVE-OBSERVABLE",
        "schema": "NSC-1-S-ONE-FERMIONIC-RELATIVE-OBSERVABLE-v1",
        "classification": "the_anomaly_consistent_fermionic_sheet_resolution_observable_has_a_positive_spectral_measure_and_places_the_collective_transition_above_the_local_cutoff",
        "authenticated_inputs": authenticated,
        "exact_fermionic_observable": {
            "Hamiltonian": str(hamiltonian),
            "relative_sheet_generator": str(sheet_generator),
            "positive_projector": str(projector_plus),
            "negative_projector": str(projector_minus),
            "Lehmann_weight": str(spectral_weight),
            "expected_weight": str(expected_weight),
            "Euclidean_connected_correlator": str(euclidean_correlator),
            "retarded_susceptibility_without_i0_notation": str(
                retarded_susceptibility
            ),
            "retarded_closed_form": str(expected_retarded),
            "spectral_measure": "rho_minus(Omega,p)=weight(p)*delta(Omega-2*E(p)); weight>=0",
        },
        "reflection_positivity": {
            "OS_matrix": str(os_matrix),
            "factorization": "M_ij(p)=weight(p)*v_i(p)*v_j(p)",
            "rank_per_momentum": 1,
            "two_by_two_principal_minor": str(two_by_two_principal_minor),
            "proof": "every_momentum_matrix_is_a_nonnegative_outer_product_and_positive_sums_or_integrals_preserve_positive_semidefiniteness",
            "heat_inverse_reused_as_physical_covariance": False,
        },
        "selected_transition_spectrum": {
            "Phi_squared_over_Lambda_squared": str(selected_mass_squared),
            "single_particle_energy_band_over_Lambda": [
                str(selected_mass),
                "1",
            ],
            "collective_sheet_excitation_band_over_Lambda": [
                str(collective_minimum),
                str(collective_maximum),
            ],
            "collective_minimum_decimal": float(collective_minimum),
            "collective_maximum_decimal": float(collective_maximum),
            "collective_mode_above_local_cutoff": bool(collective_minimum > 1),
            "rest_spectral_weight": str(rest_weight),
            "band_edge_spectral_weight": str(edge_weight),
        },
        "causal_result": {
            "measured_heat_inverse_nonpass_preserved": True,
            "physical_fermionic_relative_observable_positive": True,
            "single_boundary_pole_locally_resolved": True,
            "collective_parent_child_mode_locally_unresolved": True,
            "meaning": "the_local_room_can_resolve_the_boundary_particle_pole_but_the_correlated_sheet_change_belongs_to_the_outside_transition_sector",
            "next_equation": "couple_this_positive_sheet_generator_correlator_to_the_gauge_invariant_relative_metric_projectors_on_the_curved_parent_child_background",
        },
        "gate": {
            "inputs_authenticated": len(authenticated) == len(INPUTS),
            "Lehmann_weight_exact": sp.simplify(
                spectral_weight - expected_weight
            )
            == 0,
            "Lehmann_weight_nonnegative": bool(spectral_weight >= 0),
            "retarded_form_exact": sp.simplify(
                retarded_susceptibility - expected_retarded
            )
            == 0,
            "OS_outer_product_exact": two_by_two_principal_minor == 0,
            "OS_reflection_positive_all_momenta": True,
            "selected_collective_mode_above_cutoff": bool(
                collective_minimum > 1
            ),
            "selected_collective_band_edge_exact": collective_maximum == 2,
            "curved_metric_projector_coupling_completed": False,
            "electron_spectrum_predicted": False,
            "conviction_proof_seed_completed": False,
        },
        "nonclaims": {
            "sheet_generator_is_the_complete_graviton": False,
            "boundary_pole_is_the_electron": False,
            "curved_parent_child_linear_stability_proved": False,
            "observed_dark_response_predicted": False,
        },
        "terminal": True,
    }
    for key, value in record["gate"].items():
        if key not in {
            "curved_metric_projector_coupling_completed",
            "electron_spectrum_predicted",
            "conviction_proof_seed_completed",
        } and value is not True:
            raise RuntimeError(f"fermionic relative observable failed: {key}")
    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
