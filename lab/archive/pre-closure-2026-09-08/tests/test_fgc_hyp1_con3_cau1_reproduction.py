from __future__ import annotations

from pathlib import Path
import sys
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts.reproduce_fgc_hyp1_con3_cau1 import (  # noqa: E402
    DEFAULT_CONFIG,
    DEFAULT_OUTPUT,
    _canonical,
    load_canonical_result,
    load_config,
    record,
)


class CON3CAU1ReproductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = load_canonical_result(DEFAULT_OUTPUT)

    def test_conditional_theorem_exact_point_and_scope(self) -> None:
        certificate = self.payload["conditional_gauge_Cauchy_certificate"]
        self.assertTrue(certificate["subsidiary_operator"]["normally_hyperbolic_at_COMP1"])
        self.assertTrue(
            certificate["conditional_theorem"][
                "conditional_boundary_free_Cauchy_uniqueness_statement_derived"
            ]
        )
        self.assertTrue(certificate["compatible_point_witness"]["physical_H_zero"])
        self.assertTrue(certificate["compatible_point_witness"]["physical_M_zero"])
        self.assertTrue(all(value is False for value in self.payload["nonclaims"].values()))

    def test_canonical_reproduction_equality(self) -> None:
        self.assertEqual(self.payload, record(DEFAULT_CONFIG))

    def test_open_gate_and_canonical_mutations_fail_closed(self) -> None:
        source = DEFAULT_CONFIG.read_text(encoding="utf-8")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bad_config = root / "bad.toml"
            bad_config.write_text(
                source.replace(
                    "constraint_preserving_ACT1_IBVP_proven = false",
                    "constraint_preserving_ACT1_IBVP_proven = true",
                ),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "open gate"):
                load_config(bad_config)
            bad_result = root / "bad.json"
            bad_result.write_text(_canonical(self.payload).rstrip(), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "canonical sorted"):
                load_canonical_result(bad_result)


if __name__ == "__main__":
    unittest.main()
