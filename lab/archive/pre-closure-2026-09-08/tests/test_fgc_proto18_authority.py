from __future__ import annotations

from copy import deepcopy
from pathlib import Path
import tomllib
import unittest

from recursive_horizons.fgc.evolution.proto17_pure_construction import GENESIS_INPUT_FIELDS
from recursive_horizons.fgc.evolution.proto18_authority import (
    Proto18AuthorityError,
    build_authority,
    verify_authority_payload,
)


ROOT = Path(__file__).resolve().parents[1]
CONFIG = tomllib.loads((ROOT / "configs/fgc/fgc-1-pro18-auth1.toml").read_text("utf-8"))


class Proto18AuthorityTests(unittest.TestCase):
    def setUp(self) -> None:
        self.config = deepcopy(CONFIG)
        self.result = build_authority(ROOT, self.config)

    def test_constructs_complete_calibration_only_genesis(self) -> None:
        payload = self.result["artifact_payload"]
        self.assertEqual(payload["genesis_spec"]["namespace"], "runs/fgc-2-sf1/proto17/calibration")
        self.assertEqual(payload["claims"]["HLT15_authorized"], False)
        self.assertEqual(payload["claims"]["future_output_roots_observed"], False)
        self.assertEqual(payload["claims"]["candidate_execution_authorized"], False)
        self.assertNotIn("authorization_commit", str(self.result))
        self.assertEqual(
            payload["authority_inputs"]["raw_source_config"]["path"],
            "configs/fgc/fgc-1-pro18-frz1.toml",
        )
        self.assertIn(
            "src/recursive_horizons/fgc/evolution/proto18_authority.py",
            {item["path"] for item in payload["source_closure"]},
        )

    def test_every_genesis_input_mutation_fails_closed(self) -> None:
        for field in GENESIS_INPUT_FIELDS:
            with self.subTest(field=field):
                mutated = deepcopy(self.result)
                value = mutated["artifact_payload"]["genesis_spec"][field]
                if isinstance(value, bool):
                    changed = not value
                elif isinstance(value, int):
                    changed = value + 1
                elif isinstance(value, str):
                    changed = "f" * 64 if len(value) == 64 else f"MUTATED-{value}"
                elif isinstance(value, list):
                    changed = list(reversed(value)) if len(value) > 1 else ["MUTATED"]
                elif isinstance(value, dict):
                    changed = deepcopy(value); changed["__mutation__"] = True
                else:  # pragma: no cover - GENESIS_INPUT_FIELDS fixes the vocabulary
                    self.fail(f"unsupported genesis field type for {field}")
                mutated["artifact_payload"]["genesis_spec"][field] = changed
                with self.assertRaisesRegex(Proto18AuthorityError, "AUTH1_GENESIS_INPUT_DRIFT"):
                    verify_authority_payload(ROOT, self.config, mutated)

    def test_holdout_namespace_is_forbidden_without_observing_it(self) -> None:
        self.config["genesis_bindings"]["namespace"] = "runs/fgc-2-sf1/proto17/holdout"
        with self.assertRaisesRegex(Proto18AuthorityError, "AUTH1_HOLDOUT_NAMESPACE_FORBIDDEN"):
            build_authority(ROOT, self.config)

    def test_missing_closure_source_fails(self) -> None:
        self.config["source_closure"] = self.config["source_closure"][:-1]
        with self.assertRaisesRegex(Proto18AuthorityError, "AUTH1_SOURCE_CLOSURE_DRIFT"):
            build_authority(ROOT, self.config)

    def test_source_hash_environment_and_claim_drift_fail(self) -> None:
        for mutation, stop in (
            (lambda c: c["source_closure"].__setitem__(0, {**c["source_closure"][0], "sha256": "0" * 64}), "AUTH1_SOURCE_CLOSURE_DRIFT"),
            (lambda c: c["numerical_environment"]["contract"].__setitem__("machine", "x86_64"), "AUTH1_ENVIRONMENT_DRIFT"),
            (lambda c: c["runtime_bindings"]["raw_source_config"].__setitem__("sha256", "0" * 64), "AUTH1_HASH_DRIFT"),
            (lambda c: c["claims"].__setitem__("HLT15_authorized", True), "AUTH1_CLAIM_DRIFT"),
        ):
            with self.subTest(stop=stop):
                config = deepcopy(CONFIG)
                mutation(config)
                with self.assertRaisesRegex(Proto18AuthorityError, stop):
                    build_authority(ROOT, config)
