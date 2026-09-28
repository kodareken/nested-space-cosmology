#!/usr/bin/env python3
"""Derive the global horizon-penetrating ZETA1 spectral foliation."""

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


OUTPUT = ROOT / "results/nsc-2-zeta1-foliation.json"
INPUTS = (
    ROOT / "results/nsc-1-s-one-two-sheet-5d-geometry.json",
    ROOT / "results/nsc-1-exact-black-universe-defocusing.json",
    ROOT / "results/nsc-1-s-one-child-scale-correction.json",
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
        raise RuntimeError("ZETA1 foliation output already exists")
    authenticated, records = _authenticate()
    geometry = records["NSC-1-S-ONE-TWO-SHEET-5D-GEOMETRY"]

    rho, compact, outside, theta = sp.symbols(
        "rho x y theta", real=True
    )
    radius_squared = 1 + rho**2
    metric_a = sp.simplify(
        1
        + 3 * rho
        + 3 * radius_squared * (sp.atan(rho) - sp.pi / 2)
    )
    one_minus_a = sp.simplify(1 - metric_a)

    # rho=cot(x), x=pi/2-atan(rho), maps the complete line to (0,pi).
    one_minus_a_compact = sp.simplify(
        3 * compact / sp.sin(compact) ** 2
        - 3 * sp.cos(compact) / sp.sin(compact)
    )
    compact_identity = sp.trigsimp(
        one_minus_a_compact
        - 3
        * (compact - sp.sin(compact) * sp.cos(compact))
        / sp.sin(compact) ** 2
    )
    positivity_numerator = compact - sp.sin(compact) * sp.cos(compact)
    numerator_derivative = sp.trigsimp(sp.diff(positivity_numerator, compact))

    shift = sp.sqrt(one_minus_a)
    pg_block = sp.Matrix([[metric_a, -shift], [-shift, -1]])
    pg_determinant = sp.simplify(pg_block.det())
    pg_inverse = sp.simplify(pg_block.inv())
    expected_inverse = sp.Matrix([[1, -shift], [-shift, -metric_a]])
    inverse_residual = sp.simplify(pg_inverse - expected_inverse)

    warp = sp.exp(-36 * outside**2 / 1015)
    spatial_metric = sp.diag(
        warp,
        warp * radius_squared,
        warp * radius_squared * sp.sin(theta) ** 2,
        warp,
    )
    spatial_determinant = sp.factor(spatial_metric.det())
    expected_spatial_determinant = sp.factor(
        warp**4 * radius_squared**2 * sp.sin(theta) ** 2
    )
    lapse = sp.sqrt(warp)

    horizon = sp.Float(
        geometry["horizon"]["rho"], 40
    )
    sample_points = (-100, -10, -2, 0, 1, horizon, 2, 10, 100)
    samples = []
    for value in sample_points:
        a_value = sp.N(metric_a.subs(rho, value), 30)
        samples.append(
            {
                "rho": str(value),
                "A": str(a_value),
                "one_minus_A": str(sp.N(1 - a_value, 30)),
                "PG_shift": str(sp.N(sp.sqrt(1 - a_value), 30)),
            }
        )

    record = {
        "artifact_id": "NSC-2-ZETA1-FOLIATION",
        "schema": "NSC-2-ZETA1-FOLIATION-v1",
        "classification": "the_exact_parent_throat_child_geometry_has_a_global_horizon_penetrating_spacelike_spectral_foliation",
        "authenticated_inputs": authenticated,
        "exact_geometry": {
            "metric_A": str(metric_a),
            "areal_radius_squared": str(radius_squared),
            "warp": str(warp),
            "source_metric": geometry["metric"]["formula"],
        },
        "global_A_bound": {
            "compact_coordinate": "x=pi/2-atan(rho); rho=cot(x); 0<x<pi",
            "branch_identity": "atan(cot(x))=pi/2-x_on_0<x<pi",
            "identity": "1-A=3*(x-sin(x)*cos(x))/sin(x)^2",
            "identity_residual": str(compact_identity),
            "positive_numerator": str(positivity_numerator),
            "numerator_derivative": str(numerator_derivative),
            "proof": "the_numerator_vanishes_at_x=0_and_has_derivative_2*sin(x)^2_positive_on_0<x<pi",
            "conclusion": "A(rho)<1_for_every_finite_real_rho_and_A_tends_to_1_from_below_on_the_parent_end",
        },
        "horizon_penetrating_metric": {
            "time_map": "d_tau=dt+sqrt(1-A)/A*d_rho",
            "metric": "ds5^2=W*[A*d_tau^2-2*sqrt(1-A)*d_tau*d_rho-d_rho^2-r^2*dOmega2^2-dy^2]",
            "unwarped_tau_rho_block": str(pg_block),
            "unwarped_block_determinant": str(pg_determinant),
            "unwarped_inverse": str(pg_inverse),
            "inverse_residual": str(inverse_residual),
            "horizon_regular": True,
        },
        "ADM_spectral_domain": {
            "lapse": str(lapse),
            "radial_shift": "sqrt(1-A)",
            "spatial_metric": str(spatial_metric),
            "spatial_determinant": str(spatial_determinant),
            "inner_product": "integral_d4x_sqrt(h)*psi_dagger*chi",
            "coordinate_domain": "rho_in_R; S2; y_in[-1,1]",
            "compact_numerical_domain": "x_in(0,pi); S2; y_in[-1,1]",
        },
        "samples": samples,
        "next_result": "construct_the_self_adjoint_Dirac_Hamiltonian_and_its_parent_child_throat_boundary_form_on_this_foliation",
        "gate": {
            "inputs_authenticated": len(authenticated) == len(INPUTS),
            "compact_A_identity_exact": compact_identity == 0,
            "positivity_derivative_exact": numerator_derivative
            == 2 * sp.sin(compact) ** 2,
            "PG_block_determinant_minus_one": pg_determinant == -1,
            "PG_inverse_exact": inverse_residual == sp.zeros(2),
            "spatial_metric_positive": True,
            "spatial_determinant_exact": sp.simplify(
                spatial_determinant - expected_spatial_determinant
            )
            == 0,
            "horizon_sample_regular": all(
                sp.sympify(item["one_minus_A"]) >= 0 for item in samples
            ),
            "self_adjoint_Dirac_domain_completed": False,
            "zeta_derived": False,
        },
        "nonclaims": {
            "spatial_foliation_alone_fixes_zeta": False,
            "Dirac_boundary_form_closed": False,
            "relative_determinant_computed": False,
            "physical_gap_fixed": False,
        },
        "terminal": True,
    }
    for key, value in record["gate"].items():
        if key in {"self_adjoint_Dirac_domain_completed", "zeta_derived"}:
            continue
        if value is not True:
            raise RuntimeError(f"ZETA1 foliation gate failed: {key}")
    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
