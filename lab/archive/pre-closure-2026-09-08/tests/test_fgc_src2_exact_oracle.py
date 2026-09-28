from __future__ import annotations

import json
from pathlib import Path
import sys
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from recursive_horizons.fgc.evolution.src2_exact_oracle import (  # noqa: E402
    diagnose_exact_oracle,
    exact_dyadic,
)


FIXTURE_PATH = REPOSITORY / "configs/fgc/fgc-1-src2-pref11-point0.json"
COMPACT_POINT0 = json.loads(FIXTURE_PATH.read_text(encoding="utf-8"))


def _row(name: str) -> tuple[float, ...]:
    source = (
        COMPACT_POINT0["baseline_acceleration"]
        if name == "baseline_acceleration"
        else COMPACT_POINT0["lower_jet"][name]
    )
    return tuple(float.fromhex(item) for item in source)


class FGCSRC2ExactOracleTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.record = diagnose_exact_oracle(
            _row("u"), _row("p"), _row("q"), _row("p_r"), _row("q_r"),
            float.fromhex(COMPACT_POINT0["coordinate_radius"]),
            _row("baseline_acceleration"),
        )

    def test_binary64_inputs_become_exact_dyadics(self) -> None:
        value = float.fromhex("0x1.5af0e4cec54d8p-57")
        self.assertEqual(exact_dyadic(value).denominator & (exact_dyadic(value).denominator - 1), 0)

    def test_exact_affine_root_really_zeros_the_equations(self) -> None:
        self.assertTrue(self.record["exact_affine_nonsingular"])
        self.assertTrue(self.record["exact_root_residual_is_zero"])
        self.assertEqual(self.record["rounded_exact_root_binary64_hex"], [
            "-0x1.754c79ad5b587p-46", "0x1.4314b747bd5b8p-42", "0x1.d781691ca502ep-43",
            "-0x1.bd9644f80d8fap-49", "0x0.0p+0", "0x0.0p+0",
        ])
        self.assertLess(self.record["rounded_exact_root_residual_infinity"], 5.1e-28)

    def test_both_float_paths_mismeasure_the_rounded_exact_root(self) -> None:
        values = self.record["floating_residual_infinities"]
        self.assertGreater(values["rounded_exact_root_direct"], 7.0e-12)
        self.assertGreater(values["rounded_exact_root_generic"], 6.0e-12)
        self.assertLess(values["baseline_generic"], 1.0e-12)
        self.assertGreater(values["baseline_direct"], 1.0e-12)

    def test_one_bit_fixture_mutation_changes_the_exact_diagnosis(self) -> None:
        mutated = list(_row("baseline_acceleration"))
        mutated[0] = float.fromhex("-0x1.3b4d60d8187f3p-47")
        changed = diagnose_exact_oracle(
            _row("u"), _row("p"), _row("q"), _row("p_r"), _row("q_r"),
            float.fromhex(COMPACT_POINT0["coordinate_radius"]), tuple(mutated),
        )
        self.assertNotEqual(
            changed["baseline_exact_residual"], self.record["baseline_exact_residual"]
        )


if __name__ == "__main__":
    unittest.main()
