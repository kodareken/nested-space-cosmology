from copy import deepcopy
from pathlib import Path
import sys
import tomllib
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src")]

from recursive_horizons.fgc.evolution.protocol_v17 import (
    CADENCE_RATIONAL,
    GENESIS_INPUT_FIELDS,
    MEMBER_DESCRIPTOR_FIELDS,
    MEMBER_KEYS,
    RECEIPT_FIELDS,
    SF1_PROTOCOL_V17_ARTIFACT_ID,
    STATE_OBJECT_FIELDS,
    ARRAY_DESCRIPTOR_FIELDS,
    LEDGER_FIELDS,
    validate_proto17_schema,
)


def schema_fixture():
    return {
        "schema_version": 1,
        "artifact_id": SF1_PROTOCOL_V17_ARTIFACT_ID,
        "protocol_version": 17,
        "frozen": True,
        "member_keys": list(MEMBER_KEYS),
        "cadence": CADENCE_RATIONAL,
        "genesis_input_fields": list(GENESIS_INPUT_FIELDS),
        "member_descriptor_fields": list(MEMBER_DESCRIPTOR_FIELDS),
        "receipt_fields": list(RECEIPT_FIELDS),
        "claims": {
            "PROTO17_pure_reference_oracle_implemented": True,
            "PROTO17_runtime_implemented": False,
            "PROTO17_namespace_authorized": False,
            "PROTO17_trajectory_read": False,
        },
    }


class ProtocolV17Tests(unittest.TestCase):
    def test_schema(self):
        self.assertEqual(validate_proto17_schema(schema_fixture())["protocol_version"], 17)

    def test_live_frozen_toml_vocabulary_is_byte_identical(self):
        with (ROOT / "configs/fgc/fgc-2-sf1-protocol-v17.toml").open("rb") as handle:
            frozen = tomllib.load(handle)["trusted_genesis"]
        self.assertEqual(tuple(frozen["genesis_spec_required_fields"]), GENESIS_INPUT_FIELDS)
        self.assertEqual(tuple(frozen["member_descriptor_required_fields"]), MEMBER_DESCRIPTOR_FIELDS)
        self.assertEqual(tuple(frozen["state_object_required_fields"]), STATE_OBJECT_FIELDS)
        self.assertEqual(tuple(frozen["physical_array_required_fields"]), ARRAY_DESCRIPTOR_FIELDS)
        self.assertEqual(tuple(frozen["tdg6_ledger_required_fields"]), LEDGER_FIELDS)

    def test_schema_mutations_fail(self):
        for mutate in (
            lambda x: x.__setitem__("cadence", "1/8"),
            lambda x: x.__setitem__("member_keys", list(reversed(MEMBER_KEYS))),
            lambda x: x.__setitem__("genesis_input_fields", []),
            lambda x: x["claims"].__setitem__("PROTO17_runtime_implemented", True),
        ):
            fixture = deepcopy(schema_fixture())
            mutate(fixture)
            with self.assertRaises(ValueError):
                validate_proto17_schema(fixture)
