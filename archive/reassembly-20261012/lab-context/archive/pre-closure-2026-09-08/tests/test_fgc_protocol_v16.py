from copy import deepcopy
from pathlib import Path
import sys
import tomllib
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]
from recursive_horizons.fgc.evolution.protocol_v16 import (
    FALSE_CLAIMS,
    GENESIS_SPEC_FIELDS,
    MEMBER_DESCRIPTOR_FIELDS,
    RECEIPT_PAYLOAD_FIELDS,
    TRUE_CLAIMS,
    validate_sf1_protocol_v16,
)


class ProtocolV16Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with (ROOT / "configs/fgc/fgc-2-sf1-protocol-v16.toml").open("rb") as f:
            cls.p = tomllib.load(f)

    def test_contract(self):
        validate_sf1_protocol_v16(self.p)
        self.assertEqual(
            self.p["event_commit"]["journal_record_kind"], "COMMON_EVENT_COMMIT"
        )
        self.assertEqual(
            self.p["event_commit"]["receipt_payload_required_fields"],
            RECEIPT_PAYLOAD_FIELDS,
        )
        self.assertTrue(
            self.p["event_commit"][
                "self_successor_cursor_set_and_checkpoint_hashes_forbidden_in_receipt_payload"
            ]
        )
        self.assertTrue(
            self.p["event_commit"][
                "new_checkpoint_binds_derived_successor_cursor_set_after_receipt_hash_exists"
            ]
        )
        self.assertTrue(
            self.p["event_commit"]["receipt_sha256_is_outer_journal_record_sha256"]
        )
        self.assertEqual(
            self.p["trusted_genesis"]["genesis_spec_required_fields"],
            GENESIS_SPEC_FIELDS,
        )
        self.assertEqual(
            self.p["trusted_genesis"]["member_descriptor_required_fields"],
            MEMBER_DESCRIPTOR_FIELDS,
        )
        self.assertEqual(self.p["claims"].keys(), TRUE_CLAIMS | FALSE_CLAIMS)

    def test_mutations_fail(self):
        for fn in (
            lambda x: x["event_commit"].__setitem__("journal_record_kind", "X"),
            lambda x: x["event_commit"].__setitem__(
                "receipt_payload_required_fields", ["common_event_receipt_sha256"]
            ),
            lambda x: x["event_commit"].__setitem__(
                "frozen_common_event_cadence", "1/8"
            ),
            lambda x: x["event_commit"].__setitem__(
                "completed_common_event_index_is_previous_plus_one", False
            ),
            lambda x: x["trusted_genesis"].__setitem__("authority_artifact_id", "X"),
            lambda x: x["trusted_genesis"].__setitem__(
                "launch_authority_tuple_required_fields", []
            ),
            lambda x: x["trusted_genesis"].__setitem__(
                "genesis_spec_required_fields",
                ["genesis_spec_sha256"],
            ),
            lambda x: x["trusted_genesis"].__setitem__(
                "member_descriptor_required_fields",
                [],
            ),
            lambda x: x["claims"].__setitem__(
                "FGCQR_holdout_execution_authorized", True
            ),
        ):
            v = deepcopy(self.p)
            fn(v)
            with self.assertRaises((TypeError, ValueError)):
                validate_sf1_protocol_v16(v)
