from __future__ import annotations
from pathlib import Path
import sys, tempfile, unittest
sys.path[:0]=[str(Path(__file__).resolve().parents[1]),str(Path(__file__).resolve().parents[1]/"src")]
from scripts.reproduce_fgc_hyp1_con2_mprop1 import DEFAULT_CONFIG, DEFAULT_OUTPUT, _canonical_json_text, load_canonical_result, load_config, record
class CON2MPROP1Tests(unittest.TestCase):
 @classmethod
 def setUpClass(cls): cls.payload=record()
 def test_record(self):
  self.assertEqual(self.payload["artifact_id"],"FGC-1-HYP1-CON2-MPROP1"); self.assertTrue(all(self.payload["gate_status"].values())); self.assertTrue(all(v is False for v in self.payload["nonclaims"].values()))
  for control in self.payload["controls"].values(): self.assertTrue(control["metric_derived_gauge_subsidiary"]["local_metric_derived_gauge_subsidiary_identity_derived"]); self.assertTrue(control["reduction_subsidiary"]["on_shell_radial_constraint_is_time_constant"])
 def test_result_canonical(self):
  self.assertEqual(_canonical_json_text(self.payload),DEFAULT_OUTPUT.read_text(encoding="utf-8")); self.assertEqual(load_canonical_result()["artifact_id"],self.payload["artifact_id"])
 def test_mutations_fail(self):
  source=DEFAULT_CONFIG.read_text()
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/"bad.toml"; p.write_text(source+"\nunknown=true\n")
   with self.assertRaisesRegex(ValueError,"keys differ"): load_config(p)
   p.write_text(source.replace("evolution_authorized = false","evolution_authorized = true"))
   with self.assertRaisesRegex(ValueError,"open gate"): load_config(p)
   p.write_text(source.replace('coordinate_radius = "4"','coordinate_radius = "04"'))
   with self.assertRaisesRegex(ValueError,"canonical rational"): load_config(p)
   p.write_text(source.replace('fo1_config = "configs/fgc/fgc-1-hyp1-fo1-rc1.toml"','fo1_config = "../bad.toml"'))
   with self.assertRaisesRegex(ValueError,"stay inside|canonical"): load_config(p)
   p.write_text(source.replace('fo1_result = "results/fgc-1-hyp1-fo1-rc1.json"','fo1_result = "results/fgc-1-hyp1-dom1-qift1.json"'))
   with self.assertRaisesRegex(ValueError,"provenance differs"): load_config(p)
   duplicate=Path(d)/"duplicate.json"; duplicate.write_text(_canonical_json_text(self.payload).replace('{\n  "artifact_id"','{\n  "artifact_id": "x",\n  "artifact_id"',1))
   with self.assertRaisesRegex(ValueError,"unique-key"): load_canonical_result(duplicate)
if __name__=="__main__": unittest.main()
