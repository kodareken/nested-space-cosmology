from __future__ import annotations
from fractions import Fraction
from pathlib import Path
import sys, tempfile, unittest
sys.path[:0]=[str(Path(__file__).resolve().parents[1]),str(Path(__file__).resolve().parents[1]/"src")]
from scripts.reproduce_fgc_hyp1_dom3_uhyp1 import DEFAULT_CONFIG, DEFAULT_OUTPUT, _canonical_json_text, load_canonical_result, load_config, record

class UHYP1ReproductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls): cls.payload=load_canonical_result() if DEFAULT_OUTPUT.is_file() else record()
    def test_real_certificate_and_scope(self):
        out=self.payload["uniform_radial_certificate"]
        self.assertTrue(out["eigenframe"]["all_enclosed_frames_invertible"])
        self.assertLess(Fraction(out["eigenframe"]["neumann_inverse"]["rho_infinity"]), 1)
        self.assertEqual(len(self.payload["full_uniform_radial_certificate_sha256"]),64)
        self.assertGreater(self.payload["full_uniform_radial_certificate_canonical_json_bytes"],1_000_000)
        self.assertTrue(all(v is False for v in self.payload["nonclaims"].values()))
    @unittest.skipUnless(DEFAULT_OUTPUT.is_file(),"generation is owned by integration")
    def test_canonical_regeneration_equality(self): self.assertEqual(self.payload,load_canonical_result())
    def test_mutations_fail_closed(self):
        source=DEFAULT_CONFIG.read_text(encoding="utf-8")
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            p=root/"bad.toml"; p.write_text(source.replace('parameter_half_width = "1/','parameter_half_width = "2/'),encoding="utf-8")
            with self.assertRaisesRegex(ValueError,"canonical rational|widths"): load_config(p)
            p.write_text(source.replace("evolution_authorized = false","evolution_authorized = true"),encoding="utf-8")
            with self.assertRaisesRegex(ValueError,"open gate"): load_config(p)
            p=root/"bad.json"; p.write_text(_canonical_json_text(self.payload).rstrip(),encoding="utf-8")
            with self.assertRaisesRegex(ValueError,"canonical sorted"): load_canonical_result(p)
if __name__=="__main__": unittest.main()
