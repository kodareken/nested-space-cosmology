from __future__ import annotations

from copy import deepcopy
from fractions import Fraction
import json
from pathlib import Path
import sys
import tomllib
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.cal7_center_roundoff_diagnosis import (  # noqa: E402
    diagnose_campaign_pair,
    initial_source_floor_record,
    minkowski_reference_profiles,
)
from recursive_horizons.fgc.evolution.gr0_calibration import (  # noqa: E402
    construct_gr0_grid_initial_data,
)
from recursive_horizons.fgc.initial_data_preflight import (  # noqa: E402
    PulseParameters,
)


Q = Fraction


class FGCCAL7CenterRoundoffDiagnosisTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        with (REPOSITORY / "configs/fgc/fgc-1-cal7-run1.toml").open("rb") as handle:
            cls.plan = tomllib.load(handle)
        cls.events = [
            json.loads(line)
            for line in (
                REPOSITORY
                / "runs/fgc-2-sf1/proto10/calibration/events.jsonl"
            ).read_text(encoding="utf-8").splitlines()
        ]
        cls.campaign = json.loads(
            (
                REPOSITORY
                / "runs/fgc-2-sf1/proto10/calibration/campaign-result.json"
            ).read_text(encoding="utf-8")
        )

    def test_reference_profiles_are_exact_minkowski_adm_data(self) -> None:
        import numpy as np

        radii = np.linspace(0.0, 1.0, 9)
        values, derivative = minkowski_reference_profiles(radii)
        self.assertTrue(np.array_equal(values[:, 0], np.ones(9)))
        self.assertTrue(np.array_equal(values[:, 2], np.ones(9)))
        self.assertTrue(np.array_equal(values[:, 3], radii))
        self.assertEqual(int(np.count_nonzero(values[:, [1, 4, 5]])), 0)
        self.assertTrue(np.array_equal(derivative[:, 3], np.ones(9)))
        self.assertEqual(int(np.count_nonzero(derivative[:, [0, 1, 2, 4, 5]])), 0)

    def test_campaign_pair_is_a_repeated_source_only_obstruction(self) -> None:
        diagnosis = diagnose_campaign_pair(self.events, self.campaign)
        self.assertEqual(diagnosis["event_count"], 264)
        self.assertEqual(diagnosis["source_only_rejection_count"], 262)
        self.assertTrue(diagnosis["all_rejections_exactly_rolled_back"])
        self.assertTrue(diagnosis["no_non_source_failure_observed"])
        self.assertLess(diagnosis["terminal_time_absolute_difference"], 2.0e-13)
        self.assertLess(diagnosis["terminal_residual_relative_difference"], 4.0e-11)
        for amplitude in ("5/2", "3"):
            self.assertEqual(
                diagnosis["terminal"][amplitude]["accepted_step_index"],
                15,
            )

    def test_relabelled_non_source_failure_fails_closed(self) -> None:
        mutated = deepcopy(self.events)
        rejected = next(
            item
            for item in mutated
            if item["event_type"] == "rejected_unaccepted_source_only_proposal"
        )
        rejected["complete_failure_set"].append("kinetic_condition_limit")
        with self.assertRaisesRegex(ValueError, "rejection semantics"):
            diagnose_campaign_pair(mutated, self.campaign)

    def test_well_balanced_control_removes_the_primary_h_minus_two_floor(self) -> None:
        physical = self.plan["physical_inputs"]
        thresholds = self.plan["universal_thresholds"]
        method = self.plan["primary_method"]
        parameters = PulseParameters(
            chi_amplitude=float(Q("5/2")),
            center=float(Q(physical["chi_center"])),
            half_width=float(Q(physical["chi_half_width"])),
            phi_amplitude=float(Q(physical["phi_seed_amplitude"])),
            planck_mass=float(Q(physical["planck_mass"])),
            scalar_mass=float(Q(physical["scalar_mass"])),
            quartic_coupling=float(Q(physical["quartic_coupling"])),
        )
        initial = construct_gr0_grid_initial_data(
            parameters,
            point_count=8193,
            outer_radius=float(Q(physical["outer_radius"])),
            constraint_method=method["constraint_solve_method"],
            diagnostic_spatial_order=method["spatial_order"],
        )
        common = dict(
            initial=initial,
            amplitude="5/2",
            method="RK4",
            spatial_order=4,
            residual_tolerance=float(Q(thresholds["source_residual_infinity_max"])),
            condition_number_maximum=float(Q(thresholds["kinetic_condition_number_max"])),
        )
        raw = initial_source_floor_record(
            **common,
            representation="PROTO10_raw",
        )
        balanced = initial_source_floor_record(
            **common,
            representation="well_balanced_control",
        )
        self.assertEqual(raw.maximum_grid_index, 1)
        self.assertEqual(raw.maximum_field, "alpha")
        self.assertAlmostEqual(raw.roundoff_coefficient, 83.0 / 96.0, places=12)
        self.assertAlmostEqual(raw.residual_infinity / balanced.residual_infinity, 10.375, places=11)
        self.assertLess(balanced.residual_infinity, 1.0e-12)


if __name__ == "__main__":
    unittest.main()
