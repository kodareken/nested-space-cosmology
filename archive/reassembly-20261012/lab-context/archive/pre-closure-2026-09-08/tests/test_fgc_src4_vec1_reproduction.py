from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest
from unittest import mock


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts import reproduce_fgc_src4_vec1 as reproduce  # noqa: E402


class FGCSRC4VEC1ReproductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.record = reproduce.reproduce()
        cls.stored = json.loads(
            (REPOSITORY / "results/fgc-1-src4-vec1.json").read_text(
                encoding="utf-8"
            )
        )

    def test_stored_result_reproduces_exactly(self) -> None:
        self.assertEqual(self.record, self.stored)
        self.assertTrue(
            self.record["gate_status"][
                "SRC4_runtime_source_instrument_authorized_for_a_separately_frozen_GR0_protocol"
            ]
        )
        self.assertFalse(self.record["gate_status"]["PROTO12_frozen"])

    def test_exact_differential_and_full_grid_controls_pass(self) -> None:
        payload = self.record["artifact_payload"]
        self.assertTrue(
            payload["exact_reference_control"]["complete_residual_bitwise_zero"]
        )
        self.assertTrue(
            payload["PREF11_exact_dyadic_control"]["strict_raw_gate_passed"]
        )
        self.assertTrue(
            payload["CAP1_full_grid_control"][
                "both_strict_raw_gates_passed_in_frozen_CAP1_evidence"
            ]
        )
        self.assertFalse(
            payload["tensor_contraction_differential"][
                "point_and_seed_axes_contracted"
            ]
        )

    def test_hash_mismatched_optional_raw_fixture_fails_closed(self) -> None:
        original = reproduce._sha

        def changed(path: Path) -> str:
            if path.name == "affine-wall-fixture.npz":
                return "0" * 64
            return original(path)

        with mock.patch.object(reproduce, "_sha", side_effect=changed):
            with self.assertRaisesRegex(ValueError, "CAP1 raw fixture hash"):
                reproduce.reproduce()

    def test_optional_raw_absence_preserves_the_canonical_record(self) -> None:
        with mock.patch.object(reproduce, "_optional_raw_exists", return_value=False):
            without_raw = reproduce.reproduce()
        self.assertEqual(without_raw, self.record)


if __name__ == "__main__":
    unittest.main()
