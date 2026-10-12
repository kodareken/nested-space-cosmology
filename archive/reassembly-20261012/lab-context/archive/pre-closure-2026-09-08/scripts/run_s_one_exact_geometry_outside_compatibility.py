#!/usr/bin/env python3
"""Test whether the old exact phantom geometry can be relabeled as outside stress."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import mpmath as mp
import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
for value in (str(ROOT), str(ROOT / "src")):
    if value not in sys.path:
        sys.path.insert(0, value)

from recursive_horizons.evidence_io import canonical_json_bytes  # noqa: E402


OUTPUT = ROOT / "results/nsc-1-s-one-exact-geometry-outside-compatibility.json"
SOURCE = ROOT / "results/nsc-1-exact-black-universe-defocusing.json"


def run() -> dict[str, object]:
    if OUTPUT.exists():
        raise RuntimeError("exact-geometry outside-compatibility output already exists")

    source_raw = SOURCE.read_bytes()
    source = json.loads(source_raw)
    if source.get("terminal") is not True:
        raise RuntimeError("exact black-universe source is not terminal")

    u = sp.symbols("u", real=True)
    y = sp.Function("y")(u)
    radius = sp.sqrt(u**2 + 1)
    metric_b = (
        -3 * sp.pi / 2
        + 1 / (1 + u**2)
        + 3 * (u / (1 + u**2) + sp.atan(u))
    )
    metric_a = sp.simplify(metric_b * radius**2)

    # Bronnikov-Donskoy Eqs. (21)-(25), with y=phi'^2 for a normal scalar.
    outside_f = sp.simplify(2 * sp.diff(radius, u, 2) / radius + y)
    outside_p = sp.simplify(
        (
            metric_a * sp.diff(radius**2, u, 2)
            - radius**2 * sp.diff(metric_a, u, 2)
            - 2
            - sp.Rational(3, 2) * metric_a * outside_f
        )
        / 2
    )
    conservation = sp.simplify(
        sp.diff(outside_p * radius**4, u)
        - sp.Rational(1, 2)
        * radius**6
        * outside_f
        * sp.diff(metric_b, u)
    )
    y_prime_coefficient = sp.simplify(
        sp.diff(conservation, sp.diff(y, u))
    )
    y_coefficient = sp.simplify(sp.diff(conservation, y))
    independent = sp.simplify(
        conservation
        - y_prime_coefficient * sp.diff(y, u)
        - y_coefficient * y
    )

    exact_regular_y = -2 / (1 + u**2) ** 2
    exact_residual = sp.simplify(
        conservation.subs(
            {
                y: exact_regular_y,
                sp.diff(y, u): sp.diff(exact_regular_y, u),
            }
        )
    )
    exact_f = sp.simplify(outside_f.subs(y, exact_regular_y))
    exact_p = sp.simplify(outside_p.subs(y, exact_regular_y))

    mp.mp.dps = 60

    def a_numeric(value: mp.mpf) -> mp.mpf:
        return (
            -3 * mp.pi / 2
            + 1 / (1 + value**2)
            + 3 * (value / (1 + value**2) + mp.atan(value))
        ) * (1 + value**2)

    horizon = mp.findroot(a_numeric, (mp.mpf(1), mp.mpf(2)))
    coefficient_functions = sp.lambdify(
        u,
        (y_prime_coefficient, y_coefficient, independent),
        "mpmath",
    )
    derivative_coefficient, algebraic_coefficient, algebraic_source = (
        coefficient_functions(horizon)
    )
    horizon_regular_y = -algebraic_source / algebraic_coefficient
    expected_horizon_y = -2 / (1 + horizon**2) ** 2

    record = {
        "artifact_id": "NSC-1-S-ONE-EXACT-GEOMETRY-OUTSIDE-COMPATIBILITY",
        "schema": "NSC-1-S-ONE-EXACT-GEOMETRY-OUTSIDE-COMPATIBILITY-v1",
        "classification": "the_old_exact_phantom_geometry_cannot_be_turned_into_the_positive_local_matter_solution_by_relabeling_its_source_as_outside_geometry",
        "authenticated_input": {
            "path": str(SOURCE.relative_to(ROOT)),
            "sha256": hashlib.sha256(source_raw).hexdigest(),
            "artifact_id": source.get("artifact_id"),
        },
        "imported_equations": {
            "source": "Bronnikov_and_Donskoy_arXiv_0910.4930",
            "equation_21": "E^nu_mu=diag(-P-A*f,-P,P+A*f/2,P+A*f/2)",
            "equation_22": "(P*r^4)'=r^6*f*(A/r^2)'/2",
            "equation_24": "2*r''/r=-phi'^2+f",
            "equation_25": "A*(r^2)''-r^2*A''-2=2*P+3*A*f/2",
            "normal_scalar_variable": "y(u)=phi'(u)^2 must be nonnegative",
        },
        "compatibility_equation": {
            "radius": str(radius),
            "metric_A": str(metric_a),
            "outside_f_from_equation_24": str(outside_f),
            "outside_P_from_equation_25": str(outside_p),
            "conservation_ODE": str(conservation),
            "coefficient_of_y_prime": str(y_prime_coefficient),
            "coefficient_of_y": str(y_coefficient),
            "independent_source": str(independent),
        },
        "horizon_regularity": {
            "horizon_u": mp.nstr(horizon, 40),
            "coefficient_of_y_prime_at_horizon": mp.nstr(
                derivative_coefficient, 12
            ),
            "coefficient_of_y_at_horizon": mp.nstr(
                algebraic_coefficient, 30
            ),
            "regularity_condition": "coefficient_y*y_horizon+source=0 because A=0 removes the y_prime term",
            "forced_y_at_horizon": mp.nstr(horizon_regular_y, 30),
            "exact_solution_y_at_horizon": mp.nstr(expected_horizon_y, 30),
            "normal_scalar_condition_satisfied": horizon_regular_y >= 0,
        },
        "exact_solution": {
            "y_phi_prime_squared": str(exact_regular_y),
            "conservation_residual": str(exact_residual),
            "outside_f": str(exact_f),
            "outside_P": str(exact_p),
            "meaning": "the horizon-regular continuation returns the original negative-kinetic scalar and switches off the outside Weyl profile",
        },
        "causal_result": {
            "wrong_move": "retain_the_old_exact_metric_and_only_rename_its_phantom_source",
            "measured_obstruction": "the_horizon_regularity_condition_forces_phi_prime_squared_negative",
            "required_move": "use_the_distinct_published_positive_scalar_brane_black_universe_geometry_where_E_mu_nu_is_nonzero_then_close_its_bulk",
            "theory_scope": "this_rejects_one_transplant_not_the_nested_outside_geometry_mechanism",
        },
        "next_result": "bind_the_smooth_positive_scalar_RS2_black_universe_background_then_solve_or_import_a_regular_bulk_extension_whose_projected_Weyl_tensor_generates_its_E_mu_nu",
        "gate": {
            "source_terminal_and_hashed": source.get("terminal") is True,
            "exact_conservation_solution": exact_residual == 0,
            "outside_profile_switches_off": exact_f == 0 and exact_p == 0,
            "horizon_regular_value_matches_exact_solution": abs(
                horizon_regular_y - expected_horizon_y
            )
            < mp.mpf("1e-50"),
            "normal_scalar_transplant_passes": horizon_regular_y >= 0,
        },
        "terminal": True,
    }
    required = (
        "source_terminal_and_hashed",
        "exact_conservation_solution",
        "outside_profile_switches_off",
        "horizon_regular_value_matches_exact_solution",
    )
    if not all(record["gate"][name] for name in required):
        raise RuntimeError("exact-geometry compatibility derivation failed")
    if record["gate"]["normal_scalar_transplant_passes"]:
        raise RuntimeError("unexpected positive normal-scalar transplant")
    OUTPUT.write_bytes(canonical_json_bytes(record))
    return record


def main() -> int:
    record = run()
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
