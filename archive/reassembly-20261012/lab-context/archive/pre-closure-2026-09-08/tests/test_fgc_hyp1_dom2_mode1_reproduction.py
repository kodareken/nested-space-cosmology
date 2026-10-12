from __future__ import annotations
from pathlib import Path
import sys, tempfile, unittest
sys.path[:0]=[str(Path(__file__).resolve().parents[1]),str(Path(__file__).resolve().parents[1]/"src")]
from scripts.reproduce_fgc_hyp1_dom2_mode1 import DEFAULT_CONFIG, DEFAULT_OUTPUT, _canonical, load_canonical_result, load_config, record

class MODE1ReproductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls): cls.payload=record()
    def test_record(self):
        self.assertEqual(self.payload["artifact_id"],"FGC-1-HYP1-DOM2-MODE1")
        self.assertTrue(all(item["ref1_branch_equals_mhg1"] for item in self.payload["controls"].values()))
        self.assertTrue(self.payload["controls"]["comp1_exact_root"]["auxiliary_polynomial_mode_bases"]["tilde"]["all_residues_zero"])
        root=self.payload["exact_comp1_acceleration_root"]
        self.assertEqual(len(root["root_acceleration"]),6)
        self.assertTrue(root["solved_full_residual_zero"])
        self.assertEqual(root["solved_full_residual_vector"],["0"]*6)
        self.assertTrue(root["root_strictly_inside_qift_acceleration_box"])
        self.assertTrue(all(v is False for v in self.payload["nonclaims"].values()))
    @unittest.skipUnless(DEFAULT_OUTPUT.is_file(),"result generation belongs to integration")
    def test_canonical_result(self): self.assertEqual(self.payload,load_canonical_result())
    def test_mutations_fail_closed(self):
        source=DEFAULT_CONFIG.read_text()
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)
            for name,text,needle in (("unknown.toml",source+"\nunknown=true\n","keys differ"),("path.toml",source.replace('first_order_config = "configs/fgc/fgc-1-hyp1-fo1-rc1.toml"','first_order_config = "../bad.toml"'),"repository"),("claim.toml",source.replace("evolution_authorized = false","evolution_authorized = true"),"proof/open-gate"),("poly.toml",source.replace("exact_rref_nullspace_at_exact_auxiliary_polynomial_roots","false"),"formulation")):
                f=p/name; f.write_text(text)
                with self.assertRaises(ValueError): load_config(f)
            bad=p/"bad.json"; bad.write_text(_canonical(self.payload).replace('{\n  "artifact_id"','{\n  "artifact_id": "FGC-1-HYP1-DOM2-MODE1",\n  "artifact_id"',1))
            with self.assertRaises(ValueError): load_canonical_result(bad)
            non=p/"non.json"; non.write_text(_canonical(self.payload).rstrip())
            with self.assertRaises(ValueError): load_canonical_result(non)
if __name__=="__main__": unittest.main()
