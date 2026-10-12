from __future__ import annotations

from copy import deepcopy
from dataclasses import replace
from pathlib import Path
import sys
import tomllib
import unittest
from unittest.mock import patch


ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]

from recursive_horizons.fgc.evolution import proto18_pref26_binder as binder
from recursive_horizons.fgc.evolution.proto18_historical_replay import (
    Proto18HistoricalReplayError,
    replay_proto12_history,
)
from recursive_horizons.fgc.evolution.proto18_pref26_binder import (
    Proto18Pref26BinderError,
    bind_pref26,
)
from recursive_horizons.fgc.evolution.proto18_production_inputs import (
    VerifiedProductionCorpus,
    verify_production_corpus,
)
from scripts import reproduce_fgc_pro18_pref26 as reproduction


def toml(path: Path) -> dict:
    with path.open("rb") as handle:
        return tomllib.load(handle)


PREF26 = ROOT / "configs/fgc/fgc-1-pro18-pref26.toml"
PRO18 = ROOT / "configs/fgc/fgc-1-pro18-frz1.toml"


class Pref26BinderTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.pref26 = toml(PREF26)
        cls.pro18 = toml(PRO18)
        # This is read-only production evidence.  It never addresses a
        # PROTO17 output root, and is reused by the mutation-only tests.
        cls.corpus = verify_production_corpus(ROOT, cls.pro18)

    def test_real_read_only_happy_path_completes_only_evidence_duties(self) -> None:
        evidence = bind_pref26(ROOT, self.pref26, self.pro18)
        self.assertEqual(
            tuple(evidence.source_archive_manifest["sources"][index]["source_id"] for index in range(2)),
            ("PROTO12", "RSP2"),
        )
        self.assertTrue(evidence.duty_status["raw_bundle_byte_and_hash_recomputation"]["completed"])
        self.assertTrue(evidence.duty_status["semantic_historical_journal_payload_replay"]["completed"])
        self.assertFalse(evidence.duty_status["external_Git_authority_blob_and_live_import_validation"]["completed"])
        self.assertFalse(evidence.duty_status["namespace_reuse_and_foreign_store_rejection"]["completed"])
        self.assertTrue(self.pref26["scope"]["PREF26_read_only_source_execution_authorized"])
        self.assertFalse(self.pref26["scope"]["future_PROTO17_output_root_observation"])
        self.assertFalse(self.pref26["scope"]["future_namespace_created"])

    def test_exact_contract_and_typed_stop_mutations_fail_closed(self) -> None:
        for mutate, stop in (
            (lambda value: value["typed_stops"].__setitem__("source_byte_or_hash_drift", "WRONG"), "PREF26_INVALID_OR_NONCONVERGED_BINDER"),
            (lambda value: value["resource_bounds"].__setitem__("maximum_total_source_bytes", 1), "PREF26_SOURCE_RESOURCE_BOUND_EXCEEDED"),
            (lambda value: value["scope"].__setitem__("PREF26_read_only_source_execution_authorized", False), "PREF26_AUTHORITY_OR_LAUNCH_PROMOTION_FORBIDDEN"),
        ):
            altered = deepcopy(self.pref26)
            mutate(altered)
            with self.assertRaises(Proto18Pref26BinderError) as caught:
                bind_pref26(ROOT, altered, self.pro18)
            self.assertEqual(caught.exception.stop_id, stop)
        altered = deepcopy(self.pref26)
        altered["future_chain"]["AUTH1"] = "OTHER"
        with self.assertRaises(Proto18Pref26BinderError) as caught:
            reproduction._validate_config(altered)
        self.assertEqual(caught.exception.stop_id, "PREF26_AUTHORITY_OR_LAUNCH_PROMOTION_FORBIDDEN")

    def test_source_hash_and_resource_drift_map_to_frozen_stops(self) -> None:
        original = self.corpus.archives["PROTO12"]
        hash_drift = replace(
            original,
            archive_snapshot=replace(original.archive_snapshot, sha256="0" * 64),
        )
        altered = VerifiedProductionCorpus(
            archives={**self.corpus.archives, "PROTO12": hash_drift}, members=self.corpus.members,
        )
        with patch.object(binder, "verify_production_corpus", return_value=altered):
            with self.assertRaises(Proto18Pref26BinderError) as caught:
                bind_pref26(ROOT, self.pref26, self.pro18)
        self.assertEqual(caught.exception.stop_id, "PREF26_SOURCE_BYTE_OR_HASH_DRIFT")

        resource_drift = replace(original, total_uncompressed_bytes=original.total_uncompressed_bytes + 1)
        altered = VerifiedProductionCorpus(
            archives={**self.corpus.archives, "PROTO12": resource_drift}, members=self.corpus.members,
        )
        with patch.object(binder, "verify_production_corpus", return_value=altered):
            with self.assertRaises(Proto18Pref26BinderError) as caught:
                bind_pref26(ROOT, self.pref26, self.pro18)
        self.assertEqual(caught.exception.stop_id, "PREF26_SOURCE_RESOURCE_BOUND_EXCEEDED")

    def test_historical_final_hash_mismatch_is_rejected(self) -> None:
        archive = self.corpus.archives["PROTO12"]
        selected = {
            key: member.physical_state_sha256
            for key, member in self.corpus.members.items() if key != "SSPRK3-16385"
        }
        selected["RK4-2049"] = "0" * 64
        with self.assertRaises(Proto18HistoricalReplayError):
            replay_proto12_history(
                archive.event_snapshot.payload, archive.embedded_event_log_payload,
                archive.metadata, selected,
            )


if __name__ == "__main__":
    unittest.main()
