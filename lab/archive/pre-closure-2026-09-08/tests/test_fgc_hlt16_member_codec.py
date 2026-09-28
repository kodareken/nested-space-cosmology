from __future__ import annotations

import os
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import unittest

ROOT=Path(__file__).resolve().parents[1]; sys.path[:0]=[str(ROOT/"src")]
from recursive_horizons.fgc.evolution.hlt16_member_codec import (
    HLT16MemberCodecError,
    RUNTIME_PROTOCOL,
    _validated_proto15_cursor,
    encode_member,
    restore_member,
)
from recursive_horizons.fgc.evolution.hlt16_state_store import (
    HLT16StateStore,
    authorize_gen0_import,
)
from tests.fgc_gen0_fixture import reconstruct_sealed_gen0_members
from recursive_horizons.fgc.evolution import proto15_runtime as p15


class HLT16MemberCodecPureTests(unittest.TestCase):
    def test_persisted_proto18_cursor_hash_is_verified_not_recomputed(self):
        ledger = p15.TDG6TemporalLedger.zero(initial_time=23 / 16)
        cursor = p15._cursor({
            "protocol_artifact_id": RUNTIME_PROTOCOL,
            "campaign_id": "HLT16-TEST",
            "member_key": "RK4-2049",
            "method": "RK4",
            "point_count": 2049,
            "committed_common_event_index": 23,
            "accepted_boundary_time": p15._time_identity(23 / 16),
            "accepted_state_sha256": "a" * 64,
            "previous_step_index": 1,
            "previous_transaction_serial": 1,
            "TDG6_ledger_sha256": p15._hash(p15._ledger_mapping(ledger)),
            "cursor_generation": 0,
            "event_target_time": p15._time_identity(3 / 2),
            "attempt_serial": 0,
            "attempt_id": "HLT16-TEST:RK4-2049:23:0",
            "mode": "FRESH_READY",
            "retry_successor_payload_or_none": None,
            "journal_tip_sha256": "0" * 64,
            "cursor_chain_parent_sha256": "b" * 64,
        })
        self.assertEqual(_validated_proto15_cursor(cursor.payload).sha256, cursor.sha256)
        changed = dict(cursor.payload)
        changed["cursor_chain_sha256"] = "f" * 64
        with self.assertRaises(HLT16MemberCodecError):
            _validated_proto15_cursor(changed)


@unittest.skipUnless(os.environ.get("HLT16_REAL_GEN0") == "1", "real sealed GEN0 codec audit is opt-in")
class HLT16MemberCodecTests(unittest.TestCase):
    def test_all_six_sealed_members_round_trip_bitwise(self):
        inputs=reconstruct_sealed_gen0_members(ROOT)
        for key,member in inputs.members.items():
            cursor=inputs.evidence.checkpoint["cursors"][key]
            snapshot=encode_member(member,cursor,protocol_artifact_id=RUNTIME_PROTOCOL,campaign_id=cursor["campaign_id"],generation=0)
            before={name:value.tobytes() for name,value in snapshot.arrays.items()}
            static=(member.tracers.labels.tobytes(),member.tracers.cutoff,member.tracers.outer_radius)
            restore_member(member,snapshot)
            again=encode_member(member,cursor,protocol_artifact_id=RUNTIME_PROTOCOL,campaign_id=cursor["campaign_id"],generation=0)
            self.assertEqual(before,{name:value.tobytes() for name,value in again.arrays.items()})
            self.assertEqual(static,(member.tracers.labels.tobytes(),member.tracers.cutoff,member.tracers.outer_radius))

    def test_cursor_mutation_fails(self):
        inputs=reconstruct_sealed_gen0_members(ROOT); key=next(iter(inputs.members)); cursor=dict(inputs.evidence.checkpoint["cursors"][key]); cursor["member_key"]="FGC-QR"
        with self.assertRaises(HLT16MemberCodecError): encode_member(inputs.members[key],cursor,protocol_artifact_id=RUNTIME_PROTOCOL,campaign_id=cursor["campaign_id"],generation=0)

    def test_protocol_and_template_mutations_fail(self):
        inputs=reconstruct_sealed_gen0_members(ROOT); key=next(iter(inputs.members)); member=inputs.members[key]
        cursor=inputs.evidence.checkpoint["cursors"][key]
        with self.assertRaises(HLT16MemberCodecError):
            encode_member(member,cursor,protocol_artifact_id=cursor["protocol_artifact_id"],campaign_id=cursor["campaign_id"],generation=0)
        snapshot=encode_member(member,cursor,protocol_artifact_id=RUNTIME_PROTOCOL,campaign_id=cursor["campaign_id"],generation=0)
        changed=dict(snapshot.metadata); changed["template_sha256"]="f"*64
        with self.assertRaises(HLT16MemberCodecError):
            restore_member(member,type(snapshot)(snapshot.arrays,changed))

    def test_all_six_codec_snapshots_persist_load_and_restore(self):
        inputs=reconstruct_sealed_gen0_members(ROOT)
        with TemporaryDirectory() as directory:
            store=HLT16StateStore(Path(directory))
            for key,member in inputs.members.items():
                cursor=inputs.evidence.checkpoint["cursors"][key]
                snapshot=encode_member(member,cursor,protocol_artifact_id=RUNTIME_PROTOCOL,campaign_id=cursor["campaign_id"],generation=0)
                persisted=store.persist_snapshot(snapshot,gen0_authority=authorize_gen0_import(snapshot))
                loaded=store.load_state(persisted.descriptor_sha256)
                self.assertEqual({k:v.tobytes() for k,v in loaded.arrays.items()},{k:v.tobytes() for k,v in snapshot.arrays.items()})
                restore_member(member,type(snapshot)(loaded.arrays,loaded.descriptor["metadata"]))


if __name__ == "__main__": unittest.main()
