from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]
from scripts import reproduce_fgc_pro16_frz1 as pro16


class PRO16FreezeTests(unittest.TestCase):
    def test_canonical_compact_freeze(self):
        pro16.verify()
        r = pro16.record()
        self.assertTrue(r["gate_status"]["PROTO16_frozen"])
        self.assertFalse(r["gate_status"]["PROTO16_pretrajectory_runtime_authorized"])
        self.assertEqual(
            r["artifact_payload"]["trusted_genesis_contract"]["authority_artifact_id"],
            "FGC-1-HLT14-MON14",
        )
        self.assertTrue(
            all(
                not (ROOT / x["path"]).exists()
                for x in r["artifact_payload"]["namespace_precondition"]["records"]
            )
        )

    def test_mutations_are_bound(self):
        self.assertTrue(
            all(pro16.record()["artifact_payload"]["mutation_controls"].values())
        )
