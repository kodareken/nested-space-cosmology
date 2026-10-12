from __future__ import annotations

from copy import deepcopy
import json
import unittest

from scripts import reproduce_fgc_tdg6_pref21 as pref21


class FGCTDG6PREF21ReproductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.record = pref21.verify_canonical()

    def test_canonical_independent_binder_reproduces(self) -> None:
        self.assertEqual(self.record["artifact_id"], pref21.ARTIFACT_ID)
        self.assertEqual(
            self.record["gate_status"],
            "PASS_PRODUCTION_COMPOSITOR_IMPLEMENTATION_AUTHORIZED",
        )
        binder = self.record["artifact_payload"]["independent_binder"]
        self.assertTrue(binder["all_independent_binder_controls_pass"])
        self.assertTrue(binder["TDG6_independent_binder_completed"])
        self.assertTrue(binder["production_compositor_implementation_authorized"])
        self.assertFalse(binder["production_compositor_implemented"])

    def test_quarter_interval_and_channel_results_are_exact(self) -> None:
        binder = self.record["artifact_payload"]["independent_binder"]
        quarter = binder["independent_quarter_restriction"]
        interval = binder["independent_interval_order_theorem"]
        channel = binder["independent_channel_debit_reduction"]
        self.assertEqual(quarter["quarter_determinants"], ["1/4096"] * 4)
        self.assertTrue(quarter["direct_quarter_maps_equal_two_composed_half_maps"])
        self.assertEqual(interval["minimum_observed_order"], "3/2")
        self.assertEqual(interval["squared_multiplier_rederived_as_2_to_the_power_2p"], 8)
        self.assertTrue(interval["exact_boundary_passes"])
        self.assertTrue(interval["one_rational_unit_below_boundary_fails"])
        self.assertEqual(channel["channel_count"], 18)
        self.assertTrue(channel["one_failed_channel_vetoes_complete_admission"])
        self.assertTrue(channel["accumulation_uses_no_cancellation"])

    def test_runtime_feasibility_keeps_missing_production_work_visible(self) -> None:
        runtime = self.record["artifact_payload"]["independent_binder"][
            "immutable_runtime_extension_audit"
        ]
        self.assertTrue(runtime["production_compositor_extension_is_source_shape_feasible"])
        self.assertTrue(runtime["TDG5_shadow_attempt_deliberately_has_no_tracer_preaccept"])
        self.assertTrue(runtime["TDG6_debit_and_temporal_retry_fields_are_not_yet_in_PROTO13"])
        self.assertFalse(runtime["production_compositor_implemented"])
        self.assertFalse(runtime["production_trajectory_authorized"])

    def test_authorization_stops_before_trajectory_or_physics(self) -> None:
        payload = self.record["artifact_payload"]
        authorization = payload["implementation_authorization"]
        boundary = payload["claim_boundary"]
        self.assertTrue(boundary["TDG6_independent_binder_completed"])
        self.assertTrue(boundary["production_compositor_implementation_authorized"])
        self.assertFalse(boundary["production_compositor_implemented"])
        self.assertFalse(authorization["production_trajectory_authorized"])
        self.assertTrue(authorization["historical_PROTO13_stop_preserved"])
        self.assertFalse(authorization["historical_CAL10_result_reclassified"])
        for name in (
            "PROTO14_frozen",
            "fresh_GR0_dynamic_calibration_completed",
            "GR0_case_eligible",
            "classical_spherical_diagnostic_authorized",
            "SGBL_execution_authorized",
            "FGCQR_holdout_execution_authorized",
            "physical_transition_claim_authorized",
        ):
            self.assertFalse(boundary[name], name)

    def test_every_embedded_mutation_control_passes(self) -> None:
        mutations = self.record["artifact_payload"]["mutation_controls"]
        self.assertGreaterEqual(len(mutations), 18)
        self.assertEqual(set(mutations.values()), {True})

    def test_config_mutations_fail_closed(self) -> None:
        config = pref21.load_config()
        attacks = []

        value = deepcopy(config)
        value["immutable_lineage"]["freeze_result_sha256"] = "0" * 64
        attacks.append(value)

        value = deepcopy(config)
        value["quarter_restriction"]["quarter_restriction_determinant"] = "1/64"
        attacks.append(value)

        value = deepcopy(config)
        value["interval_order"]["squared_multiplier"] = 7
        attacks.append(value)

        value = deepcopy(config)
        value["channel_debit"]["complete_admission_is_all_of"] = False
        attacks.append(value)

        value = deepcopy(config)
        value["runtime_feasibility"]["production_compositor_implemented"] = True
        attacks.append(value)

        value = deepcopy(config)
        value["implementation_authorization"]["production_trajectory_authorized"] = True
        attacks.append(value)

        value = deepcopy(config)
        value["claims"]["GR0_case_eligible"] = True
        attacks.append(value)

        value = deepcopy(config)
        value["claims"]["FGCQR_holdout_execution_authorized"] = True
        attacks.append(value)

        for attacked in attacks:
            with self.assertRaises(ValueError):
                pref21.validate_config_data(attacked)

    def test_tracked_result_is_canonical(self) -> None:
        raw = pref21.DEFAULT_OUTPUT.read_bytes()
        stored = json.loads(raw)
        self.assertEqual(raw, pref21._canonical_bytes(stored))


if __name__ == "__main__":
    unittest.main()
