from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from recursive_horizons.fgc.evolution import tdg9_ar1_pref1_binder as binder


ROOT = Path(__file__).resolve().parents[1]
CONFIG = (ROOT / binder.CONFIG_PATH).read_bytes()
RESULT = (ROOT / binder.RESULT_PATH).read_bytes()


class TDG9AR1PREF1BinderTests(unittest.TestCase):
    def test_compact_result_is_store_blind_and_valid(self) -> None:
        forbidden = AssertionError("compact verification opened a live input")
        with (
            patch.object(binder, "_raw_snapshot", side_effect=forbidden),
            patch.object(binder, "_authenticate_pref2", side_effect=forbidden),
            patch.object(binder, "bind_raw_diagnostic", side_effect=forbidden),
        ):
            value = binder.validate_compact_result(CONFIG, RESULT)
        payload = value["artifact_payload"]
        self.assertEqual(payload["independent_replay"]["evaluator_calls"], 108)
        self.assertEqual(payload["conclusion"]["failed_channel_count"], 10)

    def test_compact_reproducer_has_no_live_flag_or_raw_requirement(self) -> None:
        source = (
            ROOT / "scripts/reproduce_fgc_tdg9_ar1_pref1.py"
        ).read_text(encoding="utf-8")
        self.assertIn("if arguments.live:", source)
        self.assertIn("else:\n        raw = output.read_bytes()", source)
        check_source = (
            ROOT / "scripts/repo_checks/temporal_diagnostics.py"
        ).read_text(encoding="utf-8")
        compact_start = check_source.index("def _check_tdg9_ar1_pref1_result")
        compact_end = check_source.index(
            "\ndef _check_tdg9_loc1_frz1_result", compact_start
        )
        compact = check_source[compact_start:compact_end]
        self.assertNotIn("--live", compact)
        self.assertNotIn(binder.RAW_NAMESPACE, compact)

    def test_compact_result_rejects_one_bit_classification_and_claim_promotion(self) -> None:
        value = json.loads(RESULT)
        mutations = (
            lambda item: item["artifact_payload"]["independent_replay"].__setitem__("classification", "arithmetic_discriminator_inconclusive"),
            lambda item: item["artifact_payload"]["claims"].__setitem__("candidate_execution_authorized", True),
            lambda item: item["artifact_payload"]["independent_replay"]["retries"][0]["channels"][0].__setitem__("row_stream_sha256", "0" * 64),
        )
        for mutation in mutations:
            changed = deepcopy(value)
            mutation(changed)
            with self.assertRaises(binder.TDG9AR1PREF1Error):
                binder.validate_compact_result(CONFIG, binder.canonical_result(changed))

    def test_aggregate_precedence_is_failure_then_resource_then_pass(self) -> None:
        def retries(*classes: str) -> list[dict[str, object]]:
            return [{"channels": [{"classification": item} for item in classes]}]
        self.assertEqual(binder._reduce(retries("resource_inconclusive", "sufficient_contraction_failure")), "legacy_enclosure_not_sole_owner_on_frozen_samples")
        self.assertEqual(binder._reduce(retries("sufficient_contraction_pass", "resource_inconclusive")), "arithmetic_discriminator_inconclusive")
        self.assertEqual(binder._reduce(retries("exact_zero", "sufficient_contraction_pass")), "legacy_enclosure_owned_nonadmission_on_all_frozen_samples")

    def test_owned_array_and_width_contract_reject_mutations(self) -> None:
        good = np.zeros((3, 2044, 6), dtype="<f8", order="C")
        self.assertIs(binder._owned(good, "good"), good)
        for bad in (
            good[:, :-1, :],
            np.zeros((3, 2044, 6), dtype=">f8"),
            np.asfortranarray(good),
            np.full((3, 2044, 6), np.nan),
        ):
            with self.assertRaises(binder.TDG9AR1PREF1Error):
                binder._owned(bad, "bad")
        class Proposal:
            initial_state = good
            candidate_state = good
            initial_time = 1.0
            final_time = 2.0
        with (
            patch.object(binder.tdg5, "_proposal_endpoint_records", return_value=(good, good)),
            patch.object(binder.tdg5, "_owned_state", side_effect=lambda value: value),
            patch.object(binder.tdg5, "_owned_rhs", side_effect=lambda value: value),
        ):
            self.assertEqual(binder._segment(Proposal()).width, 1.0)
            for final in (1.0, 0.0, float("nan"), np.float64(2.0)):
                Proposal.final_time = final
                with self.assertRaises(binder.TDG9AR1PREF1Error):
                    binder._segment(Proposal())

    def test_nofollow_reader_rejects_symlink_and_noncanonical_raw(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "real").write_text("x", encoding="ascii")
            (root / "link").symlink_to(root / "real")
            with self.assertRaises(binder.TDG9AR1PREF1Error):
                binder._read_leaf(root, "link")
            (root / "directory").mkdir()
            (root / "directory" / "inner").symlink_to(root)
            with self.assertRaises(binder.TDG9AR1PREF1Error):
                binder._read_leaf(root, "directory/inner/real")
        value = json.loads(RESULT)
        noncanonical = json.dumps(value).encode("ascii")
        with self.assertRaises(binder.TDG9AR1PREF1Error):
            binder.validate_compact_result(CONFIG, noncanonical)

    def test_binder_imports_neither_runner_nor_persisted_retry_wrapper(self) -> None:
        source = (ROOT / "src/recursive_horizons/fgc/evolution/tdg9_ar1_pref1_binder.py").read_text(encoding="utf-8")
        self.assertNotIn("run_fgc_tdg9_ar1", source)
        self.assertNotIn("tdg8_persisted_retry_replay", source)
        probe = subprocess.run(
            [
                sys.executable,
                "-I",
                "-B",
                "-c",
                (
                    "import sys; sys.path.insert(0, 'src'); "
                    "import recursive_horizons.fgc.evolution.tdg9_ar1_pref1_binder; "
                    "bad=[name for name in sys.modules if "
                    "'tdg8_persisted_retry_replay' in name or "
                    "'run_fgc_tdg9_ar1' in name]; "
                    "assert not bad, bad"
                ),
            ],
            cwd=ROOT,
            capture_output=True,
            text=True,
        )
        self.assertEqual(probe.returncode, 0, probe.stdout + probe.stderr)


if __name__ == "__main__":
    unittest.main()
