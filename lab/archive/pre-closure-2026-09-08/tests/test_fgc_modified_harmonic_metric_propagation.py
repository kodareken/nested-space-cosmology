from __future__ import annotations

from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest
from unittest.mock import patch


sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from scripts.reproduce_fgc_hyp1_fo1_rc1 import load_config  # noqa: E402
from recursive_horizons.fgc.modified_harmonic_constraints import (  # noqa: E402
    activated_compatible_state,
)
from recursive_horizons.fgc.modified_harmonic_metric_propagation import (  # noqa: E402
    SphericalThirdJetState,
    metric_derived_gauge_jet,
    metric_derived_gauge_propagation_identity,
)
from recursive_horizons.fgc.reference_connection import (  # noqa: E402
    flat_spherical_annulus_reference,
)
from recursive_horizons.fgc.reference_connection_second import (  # noqa: E402
    flat_spherical_annulus_connection_second_derivative,
)
from recursive_horizons.fgc.spherical_reduction import (  # noqa: E402
    BASE_FIELD_ORDER,
    state_from_generalized_adm_pg_fixture,
)


THIRD_ORDER = ("dttt", "dttr", "dtrr", "drrr")


def third_table(*, zero: bool = False):
    return {
        field: {
            "dttt": Q(0) if zero else Q(index, 17),
            "dttr": Q(0) if zero else Q(-index, 19),
            "dtrr": Q(0) if zero else Q(index, 23),
            "drrr": Q(0) if zero else Q(-index, 29),
        }
        for index, field in enumerate(BASE_FIELD_ORDER, 1)
    }


class MetricDerivedGaugePropagationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        predecessor = load_config()
        cls.reference = flat_spherical_annulus_reference(
            radial_domain_minimum=Q(1, 2)
        )
        cls.flat_state = state_from_generalized_adm_pg_fixture(
            predecessor.flat_fixture
        )
        cls.flat_third = SphericalThirdJetState.from_spherical_state(
            cls.flat_state, third_table(zero=True)
        )
        activated = activated_compatible_state(predecessor.flat_fixture)
        cls.activated_state = activated["state"]
        cls.activated_third = SphericalThirdJetState.from_spherical_state(
            cls.activated_state, third_table()
        )
        cls.flat = metric_derived_gauge_propagation_identity(
            cls.flat_third,
            reference=cls.reference,
            coordinate_radius=Q(4),
            tilde_normal_factor=Q(4),
            hat_normal_factor=Q(9),
        )
        cls.activated = metric_derived_gauge_propagation_identity(
            cls.activated_third,
            reference=cls.reference,
            coordinate_radius=Q(4),
            tilde_normal_factor=Q(4),
            hat_normal_factor=Q(9),
        )

    def test_flat_reference_control_has_zero_metric_derived_gauge_two_jet(self) -> None:
        gauge = self.flat["metric_derived_gauge"]
        self.assertEqual(gauge.gauge_jet.t, type(gauge.gauge_jet.t)(0))
        self.assertEqual(gauge.gauge_jet.r, type(gauge.gauge_jet.r)(0))
        self.assertEqual(
            self.flat["metric_derived_extension_operator"]["coordinate_divergence"],
            (Q(0),) * 4,
        )
        self.assertEqual(self.flat["full_modified_metric_divergence_up"], (Q(0),) * 4)
        self.assertTrue(
            self.flat["local_metric_derived_gauge_subsidiary_identity_derived"]
        )

    def test_activated_third_jet_closes_noether_and_metric_gauge_routes(self) -> None:
        result = self.activated
        gauge = result["metric_derived_gauge"]
        operator = result["metric_derived_extension_operator"]
        noether = result["unredefined_noether"]
        self.assertTrue(gauge.physical_second.mixed_partials_agree_exactly)
        self.assertTrue(gauge.mixed_partials_agree_exactly)
        self.assertTrue(
            gauge.reference_second_connection_contribution_zero_by_spherical_contraction
        )
        self.assertEqual(
            gauge.reference_second_connection_constraint_contribution,
            (((Q(0),) * 4,) * 2,) * 2,
        )
        self.assertNotEqual(gauge.gauge_jet.t.dtt, 0)
        self.assertNotEqual(gauge.gauge_jet.r.drr, 0)
        self.assertTrue(operator["two_routes_agree_exactly"])
        self.assertNotEqual(operator["coordinate_divergence"][:2], (Q(0), Q(0)))
        self.assertEqual(noether["divergence_plus_scalar_source_down"], (Q(0),) * 4)
        self.assertEqual(noether["divergence_plus_scalar_source_up"], (Q(0),) * 4)
        self.assertNotEqual(noether["divergence_up"][1], 0)
        self.assertEqual(result["subsidiary_identity_residual"], (Q(0),) * 4)
        self.assertTrue(result["closed_linear_second_order_operator_in_metric_derived_C"])
        self.assertTrue(all(value is False for value in result["nonclaims"].values()))

    def test_third_metric_and_second_reference_derivatives_are_live_inputs(self) -> None:
        changed_table = third_table()
        changed_table["h_tt"] = dict(changed_table["h_tt"])
        changed_table["h_tt"]["dttt"] += 1
        changed_state = SphericalThirdJetState.from_spherical_state(
            self.activated_state, changed_table
        )
        changed = metric_derived_gauge_jet(
            changed_state,
            reference=self.reference,
            coordinate_radius=Q(4),
            tilde_normal_factor=Q(4),
        )
        self.assertNotEqual(
            changed.gauge_jet.t.dtt,
            self.activated["metric_derived_gauge"].gauge_jet.t.dtt,
        )

        reference_second = flat_spherical_annulus_connection_second_derivative(
            self.reference, coordinate_radius=Q(4)
        )
        corrupted = [
            [
                [
                    [list(lower) for lower in upper]
                    for upper in second_direction
                ]
                for second_direction in first_direction
            ]
            for first_direction in reference_second
        ]
        # The real nonzero dd-bar-Gamma entries contract away in the spherical
        # C^t,C^r sector.  Inject a forbidden upper-r/angular-trace entry to
        # prove that the complete tensor is nevertheless consumed rather than
        # silently omitted by the implementation.
        corrupted[1][1][1][2][2] = Q(1)
        corrupted_data = tuple(
            tuple(
                tuple(
                    tuple(tuple(lower) for lower in upper)
                    for upper in second_direction
                )
                for second_direction in first_direction
            )
            for first_direction in corrupted
        )
        with patch(
            "recursive_horizons.fgc.modified_harmonic_metric_propagation.flat_spherical_annulus_connection_second_derivative",
            return_value=corrupted_data,
        ):
            corrupted_gauge = metric_derived_gauge_jet(
                self.activated_third,
                reference=self.reference,
                coordinate_radius=Q(4),
                tilde_normal_factor=Q(4),
            )
        self.assertNotEqual(
            corrupted_gauge.gauge_jet.r.drr,
            self.activated["metric_derived_gauge"].gauge_jet.r.drr,
        )
        self.assertFalse(
            corrupted_gauge.reference_second_connection_contribution_zero_by_spherical_contraction
        )

    def test_input_and_theorem_boundaries_fail_closed(self) -> None:
        missing = third_table()
        del missing["chi"]
        with self.assertRaisesRegex(ValueError, "exactly the six"):
            SphericalThirdJetState.from_spherical_state(self.activated_state, missing)
        malformed = third_table()
        del malformed["phi"]["drrr"]
        with self.assertRaisesRegex(ValueError, "exactly the four"):
            SphericalThirdJetState.from_spherical_state(self.activated_state, malformed)
        inexact = third_table()
        inexact["phi"]["dttt"] = 0.5
        with self.assertRaisesRegex(TypeError, "exact rational"):
            SphericalThirdJetState.from_spherical_state(self.activated_state, inexact)
        with self.assertRaisesRegex(ValueError, "1 < tilde"):
            metric_derived_gauge_propagation_identity(
                self.flat_third,
                reference=self.reference,
                coordinate_radius=Q(4),
                tilde_normal_factor=Q(9),
                hat_normal_factor=Q(4),
            )


if __name__ == "__main__":
    unittest.main()
