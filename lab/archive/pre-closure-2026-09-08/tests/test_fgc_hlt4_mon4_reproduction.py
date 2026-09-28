from __future__ import annotations

from pathlib import Path
import sys
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts.reproduce_fgc_hlt4_mon4 import (  # noqa: E402
    DEFAULT_CONFIG,
    DEFAULT_OUTPUT,
    _canonical,
    _validate_run_plan,
    load_canonical_result,
    load_config,
    record,
    verify_canonical,
)


class FGCHLT4MON4ReproductionTests(unittest.TestCase):
    @staticmethod
    def _postlaunch_record() -> dict[str, object]:
        observed = load_canonical_result(DEFAULT_OUTPUT)
        return record(
            DEFAULT_CONFIG,
            historical_namespace_evidence=observed["artifact_payload"][
                "namespace_precondition"
            ],
        )

    def test_canonical_result_reproduces_and_stays_pretrajectory(self) -> None:
        observed = load_canonical_result(DEFAULT_OUTPUT)
        verify_canonical(DEFAULT_OUTPUT, DEFAULT_CONFIG)
        self.assertEqual(observed, self._postlaunch_record())
        gates = observed["gate_status"]
        self.assertTrue(gates["PROTO6_successor_runtime_compositor_implemented"])
        self.assertTrue(gates["PROTO6_fresh_GR0_dynamic_calibration_authorized"])
        self.assertFalse(gates["FGCQR_holdout_execution_authorized"])
        payload = observed["artifact_payload"]
        self.assertTrue(payload["runtime_controls"]["all_controls_passed"])
        self.assertTrue(
            payload["accepted_source_operator_equivalence"]["all_controls_passed"]
        )
        self.assertEqual(len(payload["frozen_run_inputs"]), 12)
        self.assertEqual(len(payload["accepted_state_source_prechecks"]), 12)
        self.assertTrue(
            all(
                item["accepted_state_source_gate_passed"]
                and not item["trajectory_advanced"]
                for item in payload["accepted_state_source_prechecks"]
            )
        )
        self.assertEqual(len(payload["initial_common_events"]), 4)
        self.assertTrue(all(value is False for value in observed["nonclaims"].values()))

    def test_authorization_mutations_fail_closed(self) -> None:
        source = DEFAULT_CONFIG.read_text(encoding="utf-8")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            weakened = root / "weakened.toml"
            weakened.write_text(
                source.replace(
                    "only_source_only_internal_RK_stage_failure_is_retryable = true",
                    "only_source_only_internal_RK_stage_failure_is_retryable = false",
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "runtime composition"):
                load_config(weakened)

            promoted = root / "promoted.toml"
            promoted.write_text(
                source.replace(
                    "FGCQR_holdout_execution_authorized = false",
                    "FGCQR_holdout_execution_authorized = true",
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "claims"):
                load_config(promoted)

    def test_run_plan_cannot_change_physical_inputs_or_retry_budget(self) -> None:
        source_path = REPOSITORY / "configs/fgc/fgc-1-cal3-run1.toml"
        source = source_path.read_text(encoding="utf-8")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            physical = root / "physical.toml"
            physical.write_text(
                source.replace('chi_center = "12"', 'chi_center = "13"'),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "changes more"):
                _validate_run_plan(physical)

            retries = root / "retries.toml"
            retries.write_text(
                source.replace(
                    "maximum_source_retries_per_step = 32",
                    "maximum_source_retries_per_step = 31",
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "numerical ownership"):
                _validate_run_plan(retries)

    def test_noncanonical_and_duplicate_results_fail(self) -> None:
        payload = self._postlaunch_record()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            noncanonical = root / "noncanonical.json"
            noncanonical.write_text(_canonical(payload).rstrip(), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "canonical sorted"):
                load_canonical_result(noncanonical)
            duplicate = root / "duplicate.json"
            duplicate.write_text('{"a":1,"a":2}\n', encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "duplicate JSON key"):
                load_canonical_result(duplicate)


if __name__ == "__main__":
    unittest.main()
