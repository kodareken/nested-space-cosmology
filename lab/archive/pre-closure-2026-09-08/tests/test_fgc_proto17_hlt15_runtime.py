from __future__ import annotations

import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src")]
import recursive_horizons.fgc.evolution.proto17_hlt15_runtime as runtime_module
from recursive_horizons.fgc.evolution.proto17_hlt15_runtime import (
    Proto17HLT15Error,
    Proto17HLT15Runtime,
)
from recursive_horizons.fgc.evolution.proto18_auth1_inputs import build_calibration_genesis
from recursive_horizons.fgc.evolution.proto17_pure_construction import digest


class HLT15RuntimeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        record = json.loads((ROOT / "results/fgc-1-pro18-pref26.json").read_text())
        hashes = {name: "b" * 64 for name in (
            "protocol_config_sha256", "protocol_freeze_result_sha256", "hlt13_result_sha256",
            "run_plan_sha256", "runtime_module_sha256", "adapter_module_sha256", "runner_sha256",
        )}
        cls.spec = build_calibration_genesis(
            record["artifact_payload"]["AUTH1_input_evidence"], campaign_id="test-gen1",
            namespace="runs/fgc-2-sf1/proto17/calibration", source_hashes=hashes,
            numerical_environment={"test": True},
        ).genesis_spec

    @classmethod
    def binding(cls):
        return {
            "artifact_id": "FGC-1-PRO18-AUTH1",
            "authorization_commit": "a" * 40,
            "authorization_result_path": "results/fgc-1-pro18-auth1.json",
            "authorization_result_sha256": "c" * 64,
            "genesis_spec_sha256": digest(cls.spec),
            "source_closure_sha256": "d" * 64,
            "numerical_environment_sha256": "e" * 64,
            "pref26_auth1_input_evidence_sha256": "f" * 64,
        }

    def test_atomic_generation_zero_reopens_binds_authority_and_refuses_resume(self):
        with TemporaryDirectory() as temporary:
            repository = Path(temporary)
            runtime = Proto17HLT15Runtime(
                repository, self.spec, authority_binding=self.binding(),
            )
            observed = runtime.materialize_generation_zero()
            self.assertEqual(observed["checkpoint"], runtime.genesis)
            self.assertEqual(runtime.reopen_verify(), observed)
            self.assertEqual(observed["receipt"]["external_authority"], self.binding())
            self.assertTrue(observed["receipt"]["claims"]["generation_zero_materialized"])
            self.assertFalse(observed["receipt"]["claims"]["state_advanced"])
            with self.assertRaises(Proto17HLT15Error):
                runtime.materialize_generation_zero()

    def test_pre_adoption_fault_cleans_staging_without_creating_target(self):
        with TemporaryDirectory() as temporary:
            repository = Path(temporary)
            runtime = Proto17HLT15Runtime(
                repository, self.spec, authority_binding=self.binding(),
            )
            with patch.object(
                runtime_module, "_rename_exclusive",
                side_effect=Proto17HLT15Error("injected adoption fault"),
            ):
                with self.assertRaisesRegex(Proto17HLT15Error, "injected"):
                    runtime.materialize_generation_zero()
            target = repository / self.spec["namespace"]
            self.assertFalse(target.exists())
            self.assertEqual(list(target.parent.glob(".calibration.hlt15-stage-*")), [])

    def test_post_adoption_verification_fault_preserves_forensic_store(self):
        with TemporaryDirectory() as temporary:
            repository = Path(temporary)
            runtime = Proto17HLT15Runtime(
                repository, self.spec, authority_binding=self.binding(),
            )
            with patch.object(
                runtime, "reopen_verify",
                side_effect=Proto17HLT15Error("injected reopen fault"),
            ):
                with self.assertRaisesRegex(Proto17HLT15Error, "injected"):
                    runtime.materialize_generation_zero()
            target = repository / self.spec["namespace"]
            self.assertTrue(target.is_dir())
            self.assertTrue((target / "receipts/generation-zero.json").is_file())

    def test_authority_mutation_and_advance_surfaces_fail_closed(self):
        with TemporaryDirectory() as temporary:
            bad = self.binding()
            bad["genesis_spec_sha256"] = "0" * 64
            with self.assertRaises(Proto17HLT15Error):
                Proto17HLT15Runtime(Path(temporary), self.spec, authority_binding=bad)
        runtime = Proto17HLT15Runtime(
            ROOT, self.spec, authority_binding=self.binding(),
        )
        for name in ("advance", "commit_common_event", "recover", "resume", "replay_six_semantic_accepted_records"):
            self.assertFalse(hasattr(runtime, name), name)


if __name__ == "__main__":
    unittest.main()
