from __future__ import annotations

from dataclasses import replace
from pathlib import Path
import sys
import unittest

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.covariant_principal_health import (  # noqa: E402
    CovariantPrincipalBackground,
    background_from_spherical_state,
    esf_background,
)
from recursive_horizons.fgc.evolution.initial_state_bridge import (  # noqa: E402
    solve_initial_second_jet,
)
from recursive_horizons.fgc.evolution.multidirectional_health import (  # noqa: E402
    WeakCouplingThresholds,
    canonical_health_monitor_values,
    coefficient_deformation_bounds,
    esf_reference_cluster_certificate,
    weak_coupling_health_certificate,
)
from recursive_horizons.fgc.initial_data_family import (  # noqa: E402
    InitialDataParameters,
    solve_initial_data,
)


class FGCMultidirectionalHealthTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        parameters = InitialDataParameters(
            chi_amplitude=2.0,
            chi_half_width=2.0,
            phi_amplitude=1.0 / 131072.0,
        )
        solution = solve_initial_data(parameters, step_count=1024)
        index = min(
            range(2, len(solution.points) - 2),
            key=lambda item: abs(solution.points[item].radius - 12.0),
        )
        jet = solve_initial_second_jet(solution, index)
        cls.background = background_from_spherical_state(jet.accelerated_state)

    def test_flat_reference_has_six_quantified_clusters(self) -> None:
        certificate = esf_reference_cluster_certificate()
        self.assertEqual(len(certificate["clusters"]), 6)
        self.assertTrue(
            certificate["six_contours_partition_all_24_characteristics"]
        )
        for cluster in certificate["clusters"].values():
            self.assertEqual(cluster["eigenvalue_count"], 4)
            self.assertGreater(
                cluster["certified_lower_bound"],
                cluster["locked_resolvent_singular_lower_bound"],
            )

    def test_real_id1_point_passes_all_covector_sufficient_envelope(self) -> None:
        certificate = weak_coupling_health_certificate(self.background)
        self.assertTrue(
            certificate[
                "quantitative_all_covector_weak_coupling_health_envelope_passed"
            ]
        )
        self.assertTrue(all(certificate["predicates"].values()))
        self.assertGreater(certificate["physical_energy_coercivity_lower"], 0.0)
        monitor = canonical_health_monitor_values(certificate)
        self.assertGreater(monitor["physical_energy_coercivity_lower"], 0.0)
        self.assertGreater(
            monitor["minimum_cluster_resolvent_singular_margin"], 0.0
        )

    def test_norm_bound_dominates_direct_directional_deformation(self) -> None:
        bounds = coefficient_deformation_bounds(self.background)
        from recursive_horizons.fgc.evolution.covariant_principal_health import (
            normalized_first_order_matrix_from_coefficients,
            principal_coefficient_tensors,
        )

        candidate = principal_coefficient_tensors(self.background)
        reference = principal_coefficient_tensors(esf_background())
        for angle in np.linspace(0.0, 2.0 * np.pi, 33, endpoint=False):
            direction = (np.cos(angle), np.sin(angle), 0.0)
            candidate_matrix, _ = normalized_first_order_matrix_from_coefficients(
                candidate, direction
            )
            reference_matrix, _ = normalized_first_order_matrix_from_coefficients(
                reference, direction
            )
            direct = float(
                np.linalg.norm(candidate_matrix - reference_matrix, ord=2)
            )
            self.assertLessEqual(direct, bounds.companion_operator_2)

    def test_envelope_fails_closed_when_deformation_is_too_large(self) -> None:
        strong = CovariantPrincipalBackground(
            effective_planck_squared=4.0,
            effective_planck_prime=0.5,
            gb_coupling_prime=0.0,
            riemann_lower=np.zeros((4, 4, 4, 4)),
            hessian_gb_lower=np.zeros((4, 4)),
        )
        certificate = weak_coupling_health_certificate(strong)
        self.assertFalse(
            certificate[
                "quantitative_all_covector_weak_coupling_health_envelope_passed"
            ]
        )
        self.assertFalse(
            certificate["predicates"][
                "all_direction_companion_deformation_inside_ball"
            ]
        )
        with self.assertRaisesRegex(ValueError, "failed certificate"):
            canonical_health_monitor_values(certificate)

        limits = replace(
            WeakCouplingThresholds(),
            companion_deformation_maximum=1.0e-12,
        )
        self.assertFalse(
            weak_coupling_health_certificate(
                self.background, thresholds=limits
            )["quantitative_all_covector_weak_coupling_health_envelope_passed"]
        )


if __name__ == "__main__":
    unittest.main()
