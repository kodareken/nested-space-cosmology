from __future__ import annotations

from copy import deepcopy
from hashlib import sha256
from pathlib import Path
import tempfile
import unittest

from recursive_horizons.fgc.evolution.hlt16_mon16_artifact import (
    HLT16Mon16ArtifactError, PENDING_CLASSIFICATION, PENDING_SHA256,
    SEALED_CLASSIFICATION, build_mon16, validate_mon16_config,
)


def digest(value: bytes) -> str:
    return sha256(value).hexdigest()


class HLT16Mon16ArtifactTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.files = {
            "configs/fgc/fgc-1-pro19-frz1.toml": b"pro19\n",
            "results/fgc-1-pro18-pref27.json": b"pref27\n",
            "src/runtime.py": b"source\n",
            "tests/runtime_test.py": b"test\n",
            "src/status.py": b"status\n",
            "scripts/runner.py": b"runner\n",
            "configs/fgc/manifest.toml": b"manifest\n",
        }
        for relative, data in self.files.items():
            path = self.root / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(data)
        self.config = {
            "schema_version": 1,
            "artifact_id": "FGC-1-HLT16-MON16",
            "target_protocol": "FGC-2-SF1-PROTO18",
            "predecessors": [
                {"artifact_id": "FGC-1-PRO19-FRZ1", "path": "configs/fgc/fgc-1-pro19-frz1.toml",
                 "sha256": digest(self.files["configs/fgc/fgc-1-pro19-frz1.toml"])},
                {"artifact_id": "FGC-1-PRO18-PREF27", "path": "results/fgc-1-pro18-pref27.json",
                 "sha256": digest(self.files["results/fgc-1-pro18-pref27.json"])},
            ],
            "auth1_environment": {"python_version": "3.14.3", "numpy_version": "2.5.1", "system": "Darwin", "machine": "arm64"},
            "inventory": [
                {"role": "source", "path": "src/runtime.py", "sha256": PENDING_SHA256},
                {"role": "test", "path": "tests/runtime_test.py", "sha256": PENDING_SHA256},
                {"role": "status", "path": "src/status.py", "sha256": PENDING_SHA256},
                {"role": "runner", "path": "scripts/runner.py", "sha256": PENDING_SHA256},
                {"role": "manifest", "path": "configs/fgc/manifest.toml", "sha256": PENDING_SHA256},
            ],
            "claims": {
                "durable_runtime_implemented": True, "durable_runtime_qualified": True,
                "one_event_preflight_authorized": False, "state_advance_authorized": False,
                "GR0_calibration_completed": False, "candidate_execution_authorized": False,
                "physical_transition_claim_authorized": False,
            },
        }

    def tearDown(self) -> None:
        self.temp.cleanup()

    def test_pending_inventory_records_observed_hashes_without_authorization(self) -> None:
        result = build_mon16(self.root, self.config)
        self.assertEqual(result["classification"], PENDING_CLASSIFICATION)
        payload = result["artifact_payload"]
        self.assertEqual(payload["coordinator_refresh_required_paths"], sorted(item["path"] for item in self.config["inventory"]))
        self.assertTrue(all(item["binding_status"] == "pending_coordinator_refresh" for item in payload["inventory"]))
        self.assertFalse(payload["claims"]["one_event_preflight_authorized"])
        self.assertFalse(payload["explicit_nonexecution"]["live_run_store_opened"])
        self.assertEqual(
            next(item for item in payload["inventory"] if item["path"] == "src/runtime.py")["observed_sha256"],
            digest(self.files["src/runtime.py"]),
        )

    def test_sealed_refresh_requires_exact_hash_and_clears_only_that_entry(self) -> None:
        config = deepcopy(self.config)
        config["inventory"][0]["sha256"] = digest(self.files["src/runtime.py"])
        result = build_mon16(self.root, config)
        item = next(item for item in result["artifact_payload"]["inventory"] if item["path"] == "src/runtime.py")
        self.assertEqual(item["binding_status"], "sealed")
        self.assertNotIn("src/runtime.py", result["artifact_payload"]["coordinator_refresh_required_paths"])
        config["inventory"][0]["sha256"] = "0" * 64
        with self.assertRaises(HLT16Mon16ArtifactError):
            build_mon16(self.root, config)

    def test_fully_sealed_inventory_has_distinct_sealed_classification(self) -> None:
        config = deepcopy(self.config)
        for item in config["inventory"]:
            item["sha256"] = digest(self.files[item["path"]])
        result = build_mon16(self.root, config)
        self.assertEqual(result["classification"], SEALED_CLASSIFICATION)
        self.assertEqual(result["artifact_payload"]["coordinator_refresh_required_paths"], [])
        self.assertTrue(all(
            item["binding_status"] == "sealed"
            for item in result["artifact_payload"]["inventory"]
        ))

    def test_config_rejects_runtime_payload_paths_and_claim_promotion(self) -> None:
        config = deepcopy(self.config)
        config["inventory"][0]["path"] = "runs/fgc-2-sf1/array.npz"
        with self.assertRaises(HLT16Mon16ArtifactError):
            validate_mon16_config(config)
        config = deepcopy(self.config)
        config["claims"]["state_advance_authorized"] = True
        with self.assertRaises(HLT16Mon16ArtifactError):
            validate_mon16_config(config)

    def test_predecessor_must_be_sealed_and_inventory_must_cover_roles(self) -> None:
        config = deepcopy(self.config)
        config["predecessors"][0]["sha256"] = PENDING_SHA256
        with self.assertRaises(HLT16Mon16ArtifactError):
            validate_mon16_config(config)

    def test_no_follow_reader_rejects_inner_and_leaf_symlinks(self) -> None:
        config = deepcopy(self.config)
        import os
        safe = self.root / "safe"
        safe.mkdir()
        (safe / "source.py").write_bytes(b"safe\n")
        (self.root / "src" / "runtime.py").unlink()
        os.symlink(safe / "source.py", self.root / "src" / "runtime.py")
        with self.assertRaises(HLT16Mon16ArtifactError):
            build_mon16(self.root, config)
        (self.root / "src" / "runtime.py").unlink()
        os.symlink(safe, self.root / "src" / "linked")
        config["inventory"][0]["path"] = "src/linked/source.py"
        with self.assertRaises(HLT16Mon16ArtifactError):
            build_mon16(self.root, config)
        config = deepcopy(self.config)
        config["inventory"] = config["inventory"][:-1]
        with self.assertRaises(HLT16Mon16ArtifactError):
            validate_mon16_config(config)


if __name__ == "__main__":
    unittest.main()
