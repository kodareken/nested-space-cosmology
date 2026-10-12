from __future__ import annotations

from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from recursive_horizons.fgc.modified_harmonic_reduction_subsidiary import (  # noqa: E402
    ReductionDifferentialFieldJet,
    ReductionDifferentialState,
    reduction_subsidiary_identity,
)
from recursive_horizons.fgc.spherical_reduction import BASE_FIELD_ORDER  # noqa: E402


class ReductionSubsidiaryTests(unittest.TestCase):
    def test_complete_off_shell_identity_for_six_independent_fields(self) -> None:
        state = ReductionDifferentialState(
            **{
                field: ReductionDifferentialFieldJet(
                    index,
                    -index,
                    2 * index,
                    -3 * index,
                    5 * index,
                    -7 * index,
                    11 * index,
                )
                for index, field in enumerate(BASE_FIELD_ORDER, 1)
            }
        )
        out = reduction_subsidiary_identity(state)
        self.assertEqual(
            out["off_shell_identity_residual_d_t_C_plus_d_r_D_minus_K"],
            (Q(0),) * 6,
        )
        self.assertTrue(out["off_shell_identity_exact"])
        self.assertFalse(out["on_differentiable_kinematic_shell"])
        self.assertIsNone(out["on_shell_propagation_residual_d_t_C"])
        self.assertTrue(out["complete_kinematic_reduction_subsidiary_system_derived"])

    def test_on_shell_constraints_are_time_constant_and_have_no_incoming_mode(self) -> None:
        # C itself need not vanish for the identity to show it is constant;
        # setting its initial value to zero supplies the usual propagation
        # consequence without hiding that initial-data premise.
        state = ReductionDifferentialState(
            **{
                field: ReductionDifferentialFieldJet(
                    u_t=index,
                    p=index,
                    q=2 * index,
                    u_r=-index,
                    q_t=3 * index,
                    p_r=3 * index,
                    u_tr=3 * index,
                )
                for index, field in enumerate(BASE_FIELD_ORDER, 1)
            }
        )
        out = reduction_subsidiary_identity(state)
        self.assertTrue(out["on_differentiable_kinematic_shell"])
        self.assertEqual(out["on_shell_propagation_residual_d_t_C"], (Q(0),) * 6)
        self.assertTrue(out["on_shell_radial_constraint_is_time_constant"])
        self.assertEqual(out["subsidiary_characteristic_speeds"], (Q(0),) * 6)
        self.assertEqual(out["independent_spatial_curl_reduction_constraints"], 0)
        self.assertEqual(out["incoming_reduction_constraint_fields_at_inner_boundary"], ())
        self.assertEqual(out["incoming_reduction_constraint_fields_at_outer_boundary"], ())
        self.assertTrue(all(value is False for value in out["nonclaims"].values()))

    def test_exact_input_contract_fails_closed(self) -> None:
        with self.assertRaisesRegex(TypeError, "exact rational"):
            ReductionDifferentialFieldJet(0.5, 0, 0, 0, 0, 0, 0)
        with self.assertRaisesRegex(ValueError, "seven"):
            ReductionDifferentialFieldJet.from_mapping({"u_t": 0})
        with self.assertRaisesRegex(ValueError, "exactly six"):
            ReductionDifferentialState.from_mapping({})
        with self.assertRaisesRegex(TypeError, "ReductionDifferentialState"):
            reduction_subsidiary_identity(object())  # type: ignore[arg-type]


if __name__ == "__main__":
    unittest.main()
