from __future__ import annotations

from pathlib import Path
import sys
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts.reproduce_fgc_def0_obs1 import (  # noqa: E402
    DEFAULT_CONFIG,
    DEFAULT_OUTPUT,
    _canonical,
    load_canonical_result,
    load_config,
    record,
)


class DEF0OBS1ReproductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = load_canonical_result(DEFAULT_OUTPUT)

    def test_exact_controls_and_nonclaims(self) -> None:
        controls = self.payload["controls"]
        self.assertEqual(controls["minkowski_spherical"]["raychaudhuri"]["complete_rhs"], "-1/8")
        self.assertEqual(controls["contracting_flat_de_sitter"]["null_frame"]["round_sphere_classification"], "trapped")
        self.assertEqual(controls["synthetic_negative_Rkk"]["raychaudhuri"]["complete_rhs"], "1/2")
        self.assertTrue(controls["synthetic_negative_Rkk"]["raychaudhuri"]["locally_defocusing_at_point"])
        self.assertTrue(all(value is False for value in self.payload["nonclaims"].values()))

    def test_canonical_reproduction_equality(self) -> None:
        self.assertEqual(self.payload, record(DEFAULT_CONFIG))

    def test_open_gate_and_canonical_mutations_fail_closed(self) -> None:
        source = DEFAULT_CONFIG.read_text(encoding="utf-8")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bad_config = root / "bad.toml"
            bad_config.write_text(source.replace("physical_wall_derived = false", "physical_wall_derived = true"), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "open gate"):
                load_config(bad_config)
            bad_result = root / "bad.json"
            bad_result.write_text(_canonical(self.payload).rstrip(), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "canonical sorted"):
                load_canonical_result(bad_result)


if __name__ == "__main__":
    unittest.main()
