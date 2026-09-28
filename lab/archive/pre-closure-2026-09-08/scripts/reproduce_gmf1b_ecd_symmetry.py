#!/usr/bin/env python3
"""Regenerate the bounded GMF-1B-ECD-SYM1 source-certificate artifact."""

from __future__ import annotations

import argparse
import json
from math import pi, sqrt
from pathlib import Path
import sys

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

from recursive_horizons.gmf1b_ecd_symmetry import (  # noqa: E402
    audit_chiral_annular_packet,
    double_null_metric_null_cone,
    finster_equal_profile_axial_control,
    geometric_normalization,
    minimal_ecd_contact_scalar,
    required_nonclaims,
    ventrella_chiral_axial_bilinears,
    ventrella_gamma_matrix_certificate,
)


def record() -> dict[str, object]:
    kappa, F, G, radial_metric = 1.0, 0.75, 0.5, 1.0
    fixture_radius = 1.0 / (2.0 * sqrt(pi))
    axial = ventrella_chiral_axial_bilinears(F, G, fixture_radius, radial_metric)
    samples = [
        audit_chiral_annular_packet(0.5 + 2.0 * index / 256.0, 1.0, 2.0, F, G, radial_metric, kappa)
        for index in range(257)
    ]
    max_abs_axial_squared = max(abs(float(item["axial"]["axial_squared"])) for item in samples)
    return {
        "schema_version": 1,
        "project_version": "0.11.0",
        "model_id": "GMF-1B-ECD-SYM1",
        "artifact": "minimal_ecd_spherical_spinor_pair_axial_current_and_taper_preflight",
        "classification": "local_action_locked_axial_current_and_metric_null_cone_gate",
        "scope": "Local source certificate only; not an Einstein-Cartan evolution, bounce, constraint solution, horizon solution, global matching construction, or child-universe claim.",
        "conventions": {
            "metric_signature": "-+++",
            "orientation_0123": 1,
            "torsion_treatment": "minimal_ecks_cartan_equation_eliminated_algebraically",
            "gravitational_coupling_kappa": kappa,
            "contact_action_coefficient_times_kappa": 3.0 / 16.0,
            "hehl_datta_cubic_coefficient_times_kappa": 3.0 / 8.0,
        },
        "unit_normalization_algebraic_fixture": {
            "normalization": 1.0,
            "F": F,
            "G": G,
            "A_hat_0": 13.0 / 8.0,
            "A_hat_r": 5.0 / 8.0,
            "axial_squared": -9.0 / 4.0,
            "contact_scalar": -27.0 / 64.0,
        },
        "geometric_vc_fixture": {
            "areal_radius": fixture_radius,
            "radial_metric": radial_metric,
            "normalization": geometric_normalization(fixture_radius, radial_metric),
            "normalization_is_unit_within_binary64": abs(geometric_normalization(fixture_radius, radial_metric) - 1.0) <= 1.0e-14,
        },
        "finster_full_dirac_singlet_control": finster_equal_profile_axial_control(
            0.75 + 0.2j,
            0.5 - 0.3j,
            normalization=1.25,
            theta=0.7,
            phi=-0.4,
        ),
        "ventrella_left_chiral_pair": {
            "axial_current_survives_spherical_pair_sum": axial["axial_current_nonzero"],
            "algebraic_torsion_source_survives_spherical_pair_sum": axial[
                "algebraic_torsion_source_nonzero"
            ],
            "axial_contact_invariant_nonzero_when_both_amplitudes_populated": axial[
                "axial_contact_invariant_nonzero"
            ],
            "midpoint_fixture": {
                "F": F,
                "G": G,
                "A_hat_0": axial["A_hat_0"],
                "A_hat_r": axial["A_hat_r"],
                "A_hat_theta": axial["A_hat_theta"],
                "A_hat_phi": axial["A_hat_phi"],
                "axial_squared": axial["axial_squared"],
                "contact_scalar": minimal_ecd_contact_scalar(kappa, axial["axial_squared"]),
            },
            "gamma_matrix_certificate": ventrella_gamma_matrix_certificate(F, G, fixture_radius, radial_metric, 0.7, -0.4),
            "one_amplitude_zero_controls": {
                "F_zero": ventrella_chiral_axial_bilinears(0.0, G, fixture_radius, radial_metric),
                "G_zero": ventrella_chiral_axial_bilinears(F, 0.0, fixture_radius, radial_metric),
            },
        },
        "principal_cone_gate": {
            "outgoing_null": double_null_metric_null_cone(1.0, 1.0, 0.0),
            "ingoing_null": double_null_metric_null_cone(1.0, 0.0, 1.0),
            "non_null_control": double_null_metric_null_cone(1.0, 1.0, 1.0),
            "thermal_perfect_fluid_dp_drho_used_as_ecd_characteristic_speed": False,
        },
        "annular_packet_gate": {
            "radius_inner": 1.0,
            "radius_outer": 2.0,
            "sample_count": 257,
            "c_infinity_compact_support": True,
            "centre_spinor_zero_buffer": True,
            "outer_spinor_zero_buffer": True,
            "sampled_axial_bilinears_and_contact_scalar_finite": True,
            "nonzero_interior_axial_contact_invariant_present": max_abs_axial_squared > 0.0,
            "spinor_boundary_axial_current_zero": all(
                audit_chiral_annular_packet(radius, 1.0, 2.0, F, G, radial_metric, kappa)[
                    "axial_current_vanishes_at_this_outside_support_point"
                ]
                for radius in (1.0, 2.0)
            ),
            "maximum_abs_axial_squared_on_samples": max_abs_axial_squared,
        },
        **required_nonclaims(),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=REPOSITORY / "results" / "gmf-1b-ecd-symmetry.json")
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
