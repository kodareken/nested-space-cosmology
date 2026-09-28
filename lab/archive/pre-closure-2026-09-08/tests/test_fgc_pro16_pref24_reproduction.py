from pathlib import Path
import sys, unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]
from scripts import reproduce_fgc_pro16_pref24 as pref


class PREF24Tests(unittest.TestCase):
    def test_reproduces(self):
        pref.verify()
        r = pref.record()
        self.assertTrue(r["gate_status"]["PROTO17_successor_freeze_design_authorized"])
        self.assertFalse(
            r["gate_status"][
                "HLT14_implementation_and_synthetic_qualification_authorized"
            ]
        )
        self.assertFalse(r["gate_status"]["PROTO16_pretrajectory_runtime_authorized"])
