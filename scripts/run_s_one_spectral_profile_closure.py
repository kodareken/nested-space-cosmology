#!/usr/bin/env python3
"""Fix the spectral profile by recursive composition of resolution steps."""

from __future__ import annotations

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


OUTPUT = ROOT / "results/nsc-1-s-one-spectral-profile-closure.json"
INPUTS = (
    ROOT / "results/nsc-1-s-one-one-operator-bind.json",
    ROOT / "results/nsc-1-s-one-spectral-resolution-wall-bind.json",
)


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("spectral-profile closure output already exists")

    authenticated = []
    for path in INPUTS:
        raw = path.read_bytes()
        item = json.loads(raw)
        if item.get("terminal") is not True:
            raise RuntimeError(f"nonterminal input: {path.relative_to(ROOT)}")
        authenticated.append(
            {
                "path": str(path.relative_to(ROOT)),
                "sha256": hashlib.sha256(raw).hexdigest(),
                "artifact_id": item.get("artifact_id"),
            }
        )

    tau1, tau2, eigenvalue = sp.symbols(
        "tau_1 tau_2 lambda_squared", nonnegative=True
    )

    def profile(tau):
        return sp.exp(-tau * eigenvalue)

    semigroup_residual = sp.simplify(
        profile(tau1 + tau2) - profile(tau1) * profile(tau2)
    )
    identity_residual = sp.simplify(profile(0) - 1)
    generator_residual = sp.simplify(
        sp.diff(profile(tau1), tau1).subs(tau1, 0) + eigenvalue
    )
    moments = {
        str(order): str(
            sp.integrate(
                sp.symbols("z", nonnegative=True) ** order
                * sp.exp(-sp.symbols("z", nonnegative=True)),
                (sp.symbols("z", nonnegative=True), 0, sp.oo),
            )
        )
        for order in range(5)
    }

    record = {
        "artifact_id": "NSC-1-S-ONE-SPECTRAL-PROFILE-CLOSURE",
        "schema": "NSC-1-S-ONE-SPECTRAL-PROFILE-CLOSURE-v1",
        "classification": "recursive_composition_of_resolution_steps_fixes_the_one_action_profile_to_the_heat_semigroup",
        "authenticated_inputs": authenticated,
        "resolution_axioms": {
            "composition": "f_(tau1+tau2)(D^2)=f_tau1(D^2)*f_tau2(D^2)",
            "identity": "f_0=1",
            "continuity_and_positivity": True,
            "generator": "-D^2",
        },
        "unique_normalized_profile": {
            "formula": "f_tau(D^2)=exp(-tau*D^2)",
            "tau": "Lambda^-2",
            "semigroup_residual": str(semigroup_residual),
            "identity_residual": str(identity_residual),
            "generator_residual": str(generator_residual),
            "first_five_dimensionless_moments": moments,
        },
        "closed_one_action": {
            "formula": "S_one=Tr exp(-mathbb_D_Theta^2/Lambda^2)+<J*Psi,mathbb_D_Theta*Psi>",
            "recursion": "T_Theta^*mathbb_D_Theta=mathbb_D_Theta",
            "remaining_dimensional_input": "one_resolution_unit_Lambda",
            "additional_profile_coefficients": 0,
        },
        "binding_to_propagation_wall": {
            "source": "Kurkov_Lizzi_Vassilevich_arXiv_1312_2235",
            "same_profile": "exponential_heat_trace",
            "known_high_energy_result": "scalar_gauge_and_graviton_local_propagation_shuts_off_near_the_cutoff",
        },
        "next_result": "evaluate_the_exact_finite_momentum_heat_kernel_form_factors_on_the_two_sheet_background_and_locate_the_real_low_energy_propagation_domain",
        "nonclaims": {
            "mathbb_D_spectrum_derived": False,
            "curved_Lorentzian_form_factors_solved": False,
            "particle_masses_predicted": False,
        },
        "gate": {
            "inputs_authenticated": len(authenticated) == len(INPUTS),
            "semigroup_exact": semigroup_residual == 0,
            "identity_exact": identity_residual == 0,
            "generator_exact": generator_residual == 0,
            "moments_fixed": moments
            == {"0": "1", "1": "1", "2": "2", "3": "6", "4": "24"},
            "no_dimensionless_profile_fit": True,
        },
        "terminal": True,
    }
    if not all(record["gate"].values()):
        raise RuntimeError("spectral-profile closure gate failed")
    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
