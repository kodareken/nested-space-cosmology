from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import sys
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

from recursive_horizons.fgc.eft_run_envelope import (  # noqa: E402
    EFT1_OPEN1_ARTIFACT_ID,
    retained_eft_open_run_audit,
)


def _load(name: str):
    return json.loads((REPOSITORY / "results" / name).read_text(encoding="utf-8"))


class FGCEFT1OpenRunEnvelopeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.eft0 = _load("fgc-1-eft0-led1.json")
        cls.uhyp1 = _load("fgc-1-hyp1-dom3-uhyp1.json")
        cls.bnd1 = _load("fgc-1-hyp1-bnd1-md1.json")
        cls.con3 = _load("fgc-1-hyp1-con3-cau1.json")

    def audit(self, **changes):
        values = {
            "eft0": self.eft0,
            "uhyp1": self.uhyp1,
            "bnd1": self.bnd1,
            "con3": self.con3,
        }
        values.update(changes)
        return retained_eft_open_run_audit(**values)

    def test_current_evidence_is_preserved_but_run_fails_closed(self) -> None:
        output = self.audit()
        self.assertEqual(output["artifact_id"], EFT1_OPEN1_ARTIFACT_ID)
        self.assertEqual(output["passed_predicate_count"], 3)
        self.assertEqual(output["required_predicate_count"], 12)
        self.assertFalse(output["retained_EFT_open_run_envelope_passed"])
        self.assertTrue(all(output["partial_positive_evidence"].values()))
        self.assertTrue(all(value is False for value in output["stop_mask"].values()))
        self.assertTrue(all(value is False for value in output["nonclaims"].values()))
        self.assertIn(
            "invariant_proper_frequency_and_wavenumber_bound",
            output["missing_predicate_ids"],
        )
        self.assertIn(
            "constraint_complete_local_existence_or_IBVP",
            output["missing_predicate_ids"],
        )

    def test_predecessor_promotions_and_missing_positive_evidence_fail_closed(self) -> None:
        promoted = deepcopy(self.uhyp1)
        promoted["nonclaims"]["evolution_authorized"] = True
        with self.assertRaisesRegex(ValueError, "promoted nonclaim"):
            self.audit(uhyp1=promoted)

        missing = deepcopy(self.bnd1)
        missing["gate_status"][
            "uniform_frozen_radial_main_system_boundary_dissipation_passed"
        ] = False
        with self.assertRaisesRegex(ValueError, "must be true"):
            self.audit(bnd1=missing)

        circular = deepcopy(self.eft0)
        circular["conditional_local_ledger_certificate"]["cutoff"][
            "inferred_from_fixture"
        ] = True
        with self.assertRaisesRegex(ValueError, "must not be inferred"):
            self.audit(eft0=circular)


if __name__ == "__main__":
    unittest.main()
