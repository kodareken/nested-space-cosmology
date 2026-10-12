from __future__ import annotations

from copy import deepcopy
from pathlib import Path
import sys
import tomllib
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src")]

from recursive_horizons.fgc.evolution.proto19_progression_freeze import Proto19FreezeError, build_freeze

CONFIG = tomllib.loads((ROOT / "configs/fgc/fgc-1-pro19-frz1.toml").read_text("utf-8"))


class Proto19FreezeTests(unittest.TestCase):
    def valid(self): return deepcopy(CONFIG)
    def reject(self, config):
        with self.assertRaises(Proto19FreezeError): build_freeze(ROOT, config)

    def test_exact_contract_constructs_without_runtime_or_array_access(self):
        result = build_freeze(ROOT, self.valid())
        self.assertEqual(result["artifact_id"], "FGC-1-PRO19-FRZ1")
        self.assertFalse(result["artifact_payload"]["claims"]["candidate_execution_authorized"])

    def test_identity_foundation_and_protocol_self_binding_mutations_fail(self):
        for mutate in (lambda c:c.__setitem__("artifact_id","OTHER"), lambda c:c.__setitem__("protocol_artifact_id","OTHER"), lambda c:c.__setitem__("sealed_foundation_commit","0"*40), lambda c:c["bindings"][1].__setitem__("sha256","0"*64)):
            with self.subTest(mutate=mutate): c=self.valid(); mutate(c); self.reject(c)

    def test_pref27_checkpoint_receipt_tree_cursor_and_terminal_facts_fail_one_bit(self):
        for field,value in (("checkpoint_sha256","0"*64),("receipt_sha256","0"*64),("store_tree_sha256","0"*64),("generation",1),("terminal_lock",True),("cursor_mode","RETRY_PENDING")):
            with self.subTest(field=field): c=self.valid(); c["sealed_genesis"][field]=value; self.reject(c)

    def test_every_zero_ledger_and_retry_premise_mutation_fails(self):
        for field in ("zero_accepted_macro_steps","zero_source_retries","zero_CFL_retries","zero_temporal_retries","zero_debits"):
            with self.subTest(field=field): c=self.valid(); c["sealed_genesis"][field]=False; self.reject(c)
        c=self.valid(); c["sealed_genesis"]["journal_tip_sha256"]="1"+"0"*63; self.reject(c)

    def test_time_event_member_and_descriptor_set_mutations_fail(self):
        mutations=(lambda c:c["first_edge"].__setitem__("start_time","1"),lambda c:c["first_edge"].__setitem__("target_time","2"),lambda c:c["first_edge"].__setitem__("common_event_index",24),lambda c:c["first_edge"].__setitem__("event_count",2),lambda c:c["first_edge"].__setitem__("members",list(reversed(c["first_edge"]["members"]))),lambda c:c["first_edge"]["state_hashes"].__setitem__(0,"0"*64))
        for mutate in mutations:
            with self.subTest(mutate=mutate): c=self.valid(); mutate(c); self.reject(c)

    def test_selector_source_and_environment_compact_evidence_mutations_fail(self):
        for field,value in (("pref26_evidence_sha256","0"*64),("proto12_container_sha256","0"*64),("rsp2_container_sha256","0"*64),("python_version","0"),("numpy_version","0"),("system","Other"),("machine","x86_64")):
            with self.subTest(field=field): c=self.valid(); c["compact_input_evidence"][field]=value; self.reject(c)

    def test_every_scope_and_claim_prohibition_rejects_promotion(self):
        for field in ("arrays_inspected","namespace_created","state_advanced","candidate_execution","holdout_opened","trajectory_read"):
            with self.subTest(scope=field): c=self.valid(); c["scope"][field]=True; self.reject(c)
        for field in ("PROTO18_runtime_implemented","PROTO18_preflight_authorized","PROTO18_state_advance_authorized","fresh_GR0_dynamic_calibration_authorized","GR0_case_eligible","SGBL_execution_authorized","FGCQR_holdout_execution_authorized","DEF1_execution_authorized","ROB1_authorized","candidate_execution_authorized","physical_transition_claim_authorized"):
            with self.subTest(claim=field): c=self.valid(); c["claims"][field]=True; self.reject(c)

    def test_candidate_branch_holdout_and_run_field_mutations_fail(self):
        for field,value in (("branch","FGC-QR"),("amplitude","5"),("candidate_forbidden",False),("raw_import_permitted",True),("run_permitted",True)):
            with self.subTest(field=field): c=self.valid(); c["first_edge"][field]=value; self.reject(c)
