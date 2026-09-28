from __future__ import annotations
from copy import deepcopy
from pathlib import Path
import ast
import json
import os
import sys
from tempfile import TemporaryDirectory
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src")]

AUTHORIZATION = "6be5a1b88ef5e1448599beb111b96b5ddb8a79c8"

from recursive_horizons.fgc.evolution.proto19_progression_contract import (
    MemberPlan,
    Proto19ProgressionContractError,
    StopPolicy,
    TimeIdentity,
    _read_nofollow,
    construct_first_event,
)
class TestContract(unittest.TestCase):
 def test_exact_plan(self):
  p=construct_first_event(ROOT, authorization_commit=AUTHORIZATION); self.assertEqual(p.members[0].member_key,'RK4-2049'); self.assertEqual(len(p.members),6); self.assertEqual(p.start_time.rational,'23/16')
 def test_mutations(self):
  p=construct_first_event(ROOT, authorization_commit=AUTHORIZATION)
  for change in (lambda:TimeIdentity('1','0x1.8000000000000p+0'),lambda:MemberPlan('FGC-QR-1','0'*64,'FGC-QR',1),lambda:StopPolicy(1,32,32,('x',))):
   with self.assertRaises(Proto19ProgressionContractError): change()
  q=deepcopy(p.mapping()); q['branch']='FGC-QR'; self.assertNotEqual(q,p.mapping())
  with self.assertRaises(Proto19ProgressionContractError):
   construct_first_event(ROOT, authorization_commit="0")
 def test_ast_has_no_runtime_or_payload_import(self):
  tree=ast.parse((ROOT/'src/recursive_horizons/fgc/evolution/proto19_progression_contract.py').read_text())
  text=ast.unparse(tree); self.assertNotIn('numerical_engine',text); self.assertNotIn('EvolutionState',text); self.assertNotIn('np.load',text)
 def test_captured_evidence_is_exact_and_closed(self):
  paths=("results/fgc-1-pro19-frz1.json","results/fgc-1-pro18-auth1.json")
  captured={path:(ROOT/path).read_bytes() for path in paths}
  direct=construct_first_event(ROOT,authorization_commit=AUTHORIZATION)
  supplied=construct_first_event(ROOT,authorization_commit=AUTHORIZATION,evidence_bytes=captured)
  self.assertEqual(direct,supplied)
  with self.assertRaises(Proto19ProgressionContractError):
   construct_first_event(ROOT,authorization_commit=AUTHORIZATION,evidence_bytes={paths[0]:captured[paths[0]]})
  changed=json.loads(captured[paths[0]]); changed["artifact_id"]="FGC-QR"
  with self.assertRaises(Proto19ProgressionContractError):
   construct_first_event(ROOT,authorization_commit=AUTHORIZATION,evidence_bytes={**captured,paths[0]:json.dumps(changed).encode()})
 def test_nofollow_reader_rejects_leaf_and_inner_symlinks(self):
  with TemporaryDirectory() as directory:
   root=Path(directory); (root/'real').mkdir(); (root/'real'/'value').write_bytes(b'{}')
   os.symlink(root/'real'/'value',root/'leaf')
   os.symlink(root/'real',root/'inner')
   for path in ('leaf','inner/value'):
    with self.assertRaises(Proto19ProgressionContractError): _read_nofollow(root,path,'evidence')
