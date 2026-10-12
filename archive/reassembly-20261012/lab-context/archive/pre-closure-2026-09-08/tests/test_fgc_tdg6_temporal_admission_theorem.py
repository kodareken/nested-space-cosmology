from __future__ import annotations

import sys
import unittest
from fractions import Fraction
from pathlib import Path


REPOSITORY = Path(__file__).resolve().parents[1]
SOURCE = REPOSITORY / "src"
if str(SOURCE) not in sys.path:
    sys.path.insert(0, str(SOURCE))

from recursive_horizons.fgc.evolution.tdg6_temporal_admission_theorem import (  # noqa: E402
    IndependentMagnitudeInterval,
    TDG6_BINDER_CHANNELS,
    channel_debit_certificate,
    direct_restriction_matrix,
    independent_channel_decision,
    interval_order_certificate,
    quarter_restriction_certificate,
    reduce_complete_channel_ledger,
    restrict_cubic_directly,
    runtime_extension_certificate,
    tdg6_independent_binder_certificate,
)


class TDG6TemporalAdmissionTheoremTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.tdg5_source = (
            REPOSITORY
            / "src/recursive_horizons/fgc/evolution/tdg5_stage_complete_refinement_runtime.py"
        ).read_text(encoding="utf-8")
        cls.proto7_source = (
            REPOSITORY / "scripts/run_fgc_gr0_calibration_v7.py"
        ).read_text(encoding="utf-8")
        cls.proto13_source = (
            REPOSITORY / "scripts/run_fgc_gr0_calibration_v13.py"
        ).read_text(encoding="utf-8")

    def test_direct_quarter_restrictions_are_exact_and_compositional(self) -> None:
        result = quarter_restriction_certificate()
        self.assertTrue(result["four_quarter_maps_are_exact_and_injective"])
        self.assertEqual(result["quarter_determinants"], ["1/4096"] * 4)
        self.assertEqual(result["basis_control_count"], 16)
        coefficients = (Fraction(7, 3), Fraction(-2, 5), Fraction(11, 7), Fraction(4, 9))
        for index in range(4):
            restricted = restrict_cubic_directly(
                coefficients, divisor=4, interval_index=index
            )
            for local in (Fraction(0), Fraction(1, 3), Fraction(1)):
                original_point = (Fraction(index) + local) / 4
                original = sum(
                    coefficient * original_point**power
                    for power, coefficient in enumerate(coefficients)
                )
                mapped = sum(
                    coefficient * local**power
                    for power, coefficient in enumerate(restricted)
                )
                self.assertEqual(mapped, original)

    def test_invalid_restriction_partitions_fail_closed(self) -> None:
        with self.assertRaises(ValueError):
            direct_restriction_matrix(divisor=0, interval_index=0)
        with self.assertRaises(ValueError):
            direct_restriction_matrix(divisor=4, interval_index=4)
        with self.assertRaises(ValueError):
            direct_restriction_matrix(divisor=4, interval_index=-1)

    def test_interval_order_implication_and_boundary_are_exact(self) -> None:
        result = interval_order_certificate()
        self.assertTrue(result["controls_passed"])
        self.assertTrue(result["exact_boundary_passes"])
        self.assertTrue(result["one_rational_unit_below_boundary_fails"])
        self.assertEqual(result["squared_multiplier_rederived_as_2_to_the_power_2p"], 8)
        self.assertTrue(result["universal_nonnegative_interval_implication"])

    def test_enclosure_and_nonconvergent_classes_keep_distinct_meanings(self) -> None:
        enclosure = independent_channel_decision(
            IndependentMagnitudeInterval(Fraction(0), Fraction(1, 1024)),
            IndependentMagnitudeInterval(Fraction(0), Fraction(1, 2048)),
        )
        self.assertTrue(enclosure.passed)
        self.assertFalse(enclosure.order_resolved)
        self.assertEqual(enclosure.debit, Fraction(1, 2048))
        failure = independent_channel_decision(
            IndependentMagnitudeInterval(Fraction(1, 10**12), Fraction(1, 10**12)),
            IndependentMagnitudeInterval(Fraction(9, 10**13), Fraction(9, 10**13)),
        )
        self.assertFalse(failure.passed)
        self.assertTrue(failure.retryable)
        self.assertEqual(failure.classification, "resolved_order_failure")

    def test_complete_channel_reduction_is_ordered_all_of(self) -> None:
        result = channel_debit_certificate()
        self.assertTrue(result["controls_passed"])
        self.assertTrue(result["one_failed_channel_vetoes_complete_admission"])
        decision = independent_channel_decision(
            IndependentMagnitudeInterval(Fraction(8), Fraction(8)),
            IndependentMagnitudeInterval(Fraction(1), Fraction(1)),
        )
        complete = {channel: decision for channel in TDG6_BINDER_CHANNELS}
        self.assertTrue(reduce_complete_channel_ledger(complete)["all_channels_pass"])
        with self.assertRaises(ValueError):
            reduce_complete_channel_ledger(dict(list(complete.items())[:-1]))
        with self.assertRaises(ValueError):
            reduce_complete_channel_ledger(dict(reversed(tuple(complete.items()))))
        with self.assertRaises(ValueError):
            IndependentMagnitudeInterval(Fraction(-1), Fraction(0))

    def test_immutable_runtime_exposes_every_required_extension_surface(self) -> None:
        result = runtime_extension_certificate(
            tdg5_runtime_source=self.tdg5_source,
            proto7_runner_source=self.proto7_source,
            proto13_runner_source=self.proto13_source,
        )
        self.assertTrue(result["production_compositor_extension_is_source_shape_feasible"])
        self.assertTrue(result["TDG5_shadow_attempt_deliberately_has_no_tracer_preaccept"])
        self.assertTrue(result["PROTO7_exposes_tracer_preview_commit_and_rollback_inputs"])
        self.assertTrue(result["PROTO13_checkpoint_exposes_state_monitor_causal_and_tracer_extension_surface"])
        self.assertFalse(result["production_compositor_implemented"])

    def test_runtime_shape_mutations_fail_the_feasibility_audit(self) -> None:
        no_tracer_preview = runtime_extension_certificate(
            tdg5_runtime_source=self.tdg5_source,
            proto7_runner_source=self.proto7_source.replace(
                "preview_advance", "preview_removed", 1
            ),
            proto13_runner_source=self.proto13_source,
        )
        self.assertFalse(no_tracer_preview["production_compositor_extension_is_source_shape_feasible"])
        wrong_commit = runtime_extension_certificate(
            tdg5_runtime_source=self.tdg5_source.replace(
                "prepared.fine_monitor_state", "prepared.original_monitor_state", 1
            ),
            proto7_runner_source=self.proto7_source,
            proto13_runner_source=self.proto13_source,
        )
        self.assertFalse(wrong_commit["production_compositor_extension_is_source_shape_feasible"])
        no_checkpoint_writer = runtime_extension_certificate(
            tdg5_runtime_source=self.tdg5_source,
            proto7_runner_source=self.proto7_source,
            proto13_runner_source=self.proto13_source.replace(
                "_write_checkpoint", "_write_removed"
            ),
        )
        self.assertFalse(no_checkpoint_writer["production_compositor_extension_is_source_shape_feasible"])

    def test_binder_authorizes_implementation_but_no_trajectory(self) -> None:
        result = tdg6_independent_binder_certificate(
            tdg5_runtime_source=self.tdg5_source,
            proto7_runner_source=self.proto7_source,
            proto13_runner_source=self.proto13_source,
        )
        self.assertTrue(result["all_independent_binder_controls_pass"])
        self.assertTrue(result["TDG6_independent_binder_completed"])
        self.assertTrue(result["production_compositor_implementation_authorized"])
        self.assertFalse(result["production_compositor_implemented"])
        self.assertFalse(result["production_state_advanced"])
        self.assertFalse(result["FGCQR_trajectory_read"])

    def test_theorem_module_does_not_import_freeze_design(self) -> None:
        source = (
            REPOSITORY
            / "src/recursive_horizons/fgc/evolution/tdg6_temporal_admission_theorem.py"
        ).read_text(encoding="utf-8")
        self.assertNotIn("tdg6_temporal_admission_design", source)


if __name__ == "__main__":
    unittest.main()
