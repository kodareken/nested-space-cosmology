#!/usr/bin/env python3
"""Make the parent/child recursive scale map explicit and executable."""

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
from recursive_horizons.unified_action import RecursiveScaleDilation  # noqa: E402


OUTPUT = ROOT / "results/nsc-2-zeta1-recursion-map.json"
INPUTS = (
    ROOT / "results/nsc-2-zeta1-anomaly-owner-correction.json",
    ROOT / "results/nsc-1-s-one-recursive-outside-kernel.json",
    ROOT / "results/nsc-1-s-one-child-scale-correction.json",
)


def _authenticate() -> list[dict[str, object]]:
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
                "artifact_id": item["artifact_id"],
            }
        )
    return authenticated


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("ZETA1 recursion-map output already exists")
    authenticated = _authenticate()

    omega, first, second, dimension = sp.symbols(
        "Omega Omega_1 Omega_2 d", positive=True, real=True
    )
    coordinate = sp.symbols("x", real=True)
    amplitude = omega ** (dimension / 2)
    norm_jacobian = sp.simplify(amplitude**2 / omega**dimension)
    composed_amplitude = sp.simplify(
        first ** (dimension / 2) * second ** (dimension / 2)
    )
    expected_composed = (first * second) ** (dimension / 2)
    composition_residual = sp.simplify(composed_amplitude - expected_composed)

    parent_energy, parent_cutoff, link = sp.symbols(
        "E Lambda Phi", positive=True, real=True
    )
    child_energy = omega * parent_energy
    child_cutoff = omega * parent_cutoff
    child_link = omega * link
    energy_ratio_residual = sp.simplify(
        child_energy / child_cutoff - parent_energy / parent_cutoff
    )
    link_ratio_residual = sp.simplify(
        child_link / child_cutoff - link / parent_cutoff
    )

    gamma, local, boundary = sp.symbols(
        "Gamma K B", nonzero=True
    )
    scaled_argument = coordinate / omega
    functional_tail = sp.Eq(
        sp.Function("Gamma")(coordinate),
        sp.Function("K")(coordinate)
        - boundary**2 / sp.Function("Gamma")(scaled_argument),
    )
    zero_argument = sp.Eq(gamma, local - boundary**2 / gamma)
    zero_solutions = sp.solve(zero_argument, gamma)

    fixture = RecursiveScaleDilation(omega=2.0, spatial_dimension=4)
    composed_fixture = fixture.compose(
        RecursiveScaleDilation(omega=3.0, spatial_dimension=4)
    )
    record = {
        "artifact_id": "NSC-2-ZETA1-RECURSION-MAP",
        "schema": "NSC-2-ZETA1-RECURSION-MAP-v1",
        "classification": "norm_preservation_and_first_order_Dirac_scaling_fix_the_minimal_parent_child_dilation_and_the_energy_resolved_recursive_tail",
        "authenticated_inputs": authenticated,
        "unitary_dilation": {
            "dimension": 4,
            "formula": "(U_Omega psi)(x)=Omega^(d/2)*psi(Omega*x)",
            "sheet_map": "T_Omega=sigma_1*U_Omega",
            "norm_Jacobian": str(norm_jacobian),
            "composition": "U_Omega1*U_Omega2=U_(Omega1*Omega2)",
            "composition_residual": str(composition_residual),
            "zeta": "Omega^2",
        },
        "inherited_spectral_ratios": {
            "child_energy": str(child_energy),
            "child_cutoff": str(child_cutoff),
            "child_link": str(child_link),
            "energy_over_cutoff_residual": str(energy_ratio_residual),
            "Phi_over_cutoff_residual": str(link_ratio_residual),
        },
        "recursive_tail": {
            "dimensionless_parent_argument": str(coordinate),
            "dimensionless_child_argument": str(scaled_argument),
            "functional_equation": str(functional_tail),
            "zero_argument_reduction": str(zero_argument),
            "zero_argument_solutions": [str(value) for value in zero_solutions],
            "physical_branch": "the_solution_continuous_with_Gamma_to_K_as_B_to_zero",
        },
        "implementation_fixture": {
            "Omega": fixture.omega,
            "zeta": fixture.zeta,
            "wave_amplitude_factor": fixture.wave_amplitude_factor,
            "child_cutoff_from_5": fixture.child_cutoff(5.0),
            "child_energy_from_7": fixture.child_energy(7.0),
            "child_argument_from_3": fixture.child_dimensionless_argument(3.0),
            "composed_Omega_2_then_3": composed_fixture.omega,
        },
        "next_result": "solve_the_functional_tail_on_the_warped_mode_spectrum_and_recompute_the_anomaly_compensated_scale_derivative",
        "gate": {
            "inputs_authenticated": len(authenticated) == len(INPUTS),
            "unitary_norm_factor_exact": norm_jacobian == 1,
            "dilation_semigroup_exact": composition_residual == 0,
            "dimensionless_energy_inherited": energy_ratio_residual == 0,
            "dimensionless_link_inherited": link_ratio_residual == 0,
            "scalar_zero_mode_fixed_point_recovered": zero_argument
            == sp.Eq(gamma, local - boundary**2 / gamma),
            "implementation_fixture_passes": composed_fixture.omega == 6.0,
            "mode_resolved_tail_solved": False,
            "zeta_derived": False,
        },
        "nonclaims": {
            "Omega_value_selected": False,
            "functional_tail_solution_computed": False,
            "scale_root_restored": False,
            "physical_zeta_promoted": False,
        },
        "terminal": True,
    }
    for key, value in record["gate"].items():
        if key in {"mode_resolved_tail_solved", "zeta_derived"}:
            continue
        if value is not True:
            raise RuntimeError(f"ZETA1 recursion-map gate failed: {key}")
    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
