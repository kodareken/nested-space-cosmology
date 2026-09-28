from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts import reproduce_fgc_cal11_pref22 as pref22  # noqa: E402
from recursive_horizons.fgc.evolution.cal11_proto14_campaign_diagnosis import (  # noqa: E402
    _checkpoint_payload,
    load_canonical_jsonl,
    source_pinned_static_replay,
)


class FGCCAL11PREF22ReproductionTests(unittest.TestCase):
    def test_canonical_reproduction_reconstructs_raw_bundle(self) -> None:
        completed = subprocess.run(
            [sys.executable, "scripts/reproduce_fgc_cal11_pref22.py", "--check"],
            cwd=REPOSITORY, check=True, capture_output=True, text=True,
        )
        self.assertIn('"artifact_id": "FGC-1-CAL11-PREF22"', completed.stdout)
        self.assertIn(
            '"raw_bundle_loaded_for_this_verification": true', completed.stdout
        )

    def test_tracked_result_is_canonical_and_closed(self) -> None:
        raw = pref22.OUTPUT.read_bytes()
        record = json.loads(raw)
        self.assertEqual(raw, pref22.canonical(record))
        self.assertEqual(
            record["classification"],
            "completed_PROTO14_invalid_runtime_terminal_result",
        )
        self.assertFalse(record["gate_status"]["GR0_case_eligible"])
        self.assertFalse(
            record["gate_status"]["FGCQR_holdout_execution_authorized"]
        )
        replay = record["artifact_payload"]["terminal_checkpoint_replay"]
        self.assertEqual(replay["completed_common_event_index"], 23)
        self.assertTrue(replay["all_members_zero_progress"])

    def test_complete_absence_is_portable_but_partial_bundle_fails(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            paths = {name: root / name for name in ("a", "b", "c", "d")}
            self.assertFalse(pref22._raw_bundle_is_complete_or_absent(paths))
            paths["a"].write_bytes(b"present")
            with self.assertRaisesRegex(ValueError, "partial"):
                pref22._raw_bundle_is_complete_or_absent(paths)
            for name in ("b", "c", "d"):
                paths[name].write_bytes(b"present")
            self.assertTrue(pref22._raw_bundle_is_complete_or_absent(paths))

    def test_portable_absence_reconstructs_identical_compact_record(self) -> None:
        config = pref22.load_config()
        present = pref22.build(config)
        with mock.patch.object(
            pref22, "_raw_bundle_is_complete_or_absent", return_value=False
        ):
            absent = pref22.build(config)
        self.assertEqual(absent, present)

    def test_claim_promotion_fails_closed(self) -> None:
        mutated = deepcopy(pref22.load_config())
        mutated["claims"]["FGCQR_holdout_execution_authorized"] = True
        with self.assertRaisesRegex(ValueError, "configuration contract"):
            pref22.validate_config_data(mutated)

    def test_terminal_identity_mutation_fails_closed(self) -> None:
        mutated = deepcopy(pref22.load_config())
        mutated["terminal_identity"]["error_message"] = "rounded into a pass"
        with self.assertRaisesRegex(ValueError, "configuration contract"):
            pref22.validate_config_data(mutated)

    def test_authorization_hash_or_path_mutation_fails_closed(self) -> None:
        mutated = deepcopy(pref22.load_config())
        mutated["authorization_lineage"]["tdg6_runtime_module"] = (
            "src/recursive_horizons/fgc/evolution/proto14_runtime.py"
        )
        with self.assertRaisesRegex(ValueError, "configuration contract"):
            pref22.validate_config_data(mutated)

    def test_source_pinned_guard_mutation_fails_closed(self) -> None:
        config = pref22.load_config()
        source = pref22._git_blob(
            config["immutable_campaign"]["authorization_commit"],
            config["source_pinned_replay"]["runtime_source"],
        )
        with self.assertRaisesRegex(ValueError, "source hash"):
            source_pinned_static_replay(
                source_payload=source.replace(b"bitwise uniform", b"bitwise changed"),
                expected=config["source_pinned_replay"],
            )

    def test_noncanonical_event_log_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            event_log = Path(directory) / "events.jsonl"
            event_log.write_bytes(b'{"b": 2, "a": 1}\n')
            with self.assertRaisesRegex(ValueError, "canonical JSONL"):
                load_canonical_jsonl(event_log)

    def test_missing_checkpoint_array_fails_closed(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            checkpoint = Path(directory) / "checkpoint.npz"
            np.savez(
                checkpoint,
                metadata_utf8=np.frombuffer(b"{}", dtype=np.uint8),
                event_log_utf8=np.frombuffer(b"", dtype=np.uint8),
            )
            with self.assertRaisesRegex(ValueError, "array schema"):
                _checkpoint_payload(checkpoint)

    def test_binder_exposes_no_campaign_or_resume_control(self) -> None:
        source = (REPOSITORY / "scripts/reproduce_fgc_cal11_pref22.py").read_text()
        diagnosis = (
            REPOSITORY
            / "src/recursive_horizons/fgc/evolution/cal11_proto14_campaign_diagnosis.py"
        ).read_text()
        self.assertNotIn("run_campaign", source)
        self.assertNotIn('add_argument("--resume"', source)
        self.assertIn('"output_namespace_creation": False', diagnosis)


if __name__ == "__main__":
    unittest.main()
