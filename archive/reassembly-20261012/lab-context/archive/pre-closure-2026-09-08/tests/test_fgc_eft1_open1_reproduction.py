from __future__ import annotations

from pathlib import Path
import sys
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(REPOSITORY), str(REPOSITORY / "src")]

from scripts.reproduce_fgc_eft1_open1 import (  # noqa: E402
    DEFAULT_CONFIG,
    DEFAULT_OUTPUT,
    _canonical,
    load_canonical_result,
    load_config,
    record,
)


class EFT1OPEN1ReproductionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.payload = load_canonical_result(DEFAULT_OUTPUT)

    def test_exact_blocker_ledger_and_scope(self) -> None:
        audit = self.payload["open_run_authorization_audit"]
        self.assertEqual(audit["passed_predicate_count"], 3)
        self.assertEqual(audit["required_predicate_count"], 12)
        self.assertFalse(audit["retained_EFT_open_run_envelope_passed"])
        self.assertFalse(self.payload["gate_status"]["evolution_authorized"])
        self.assertTrue(all(value is False for value in self.payload["nonclaims"].values()))

    def test_canonical_reproduction_equality(self) -> None:
        self.assertEqual(self.payload, record(DEFAULT_CONFIG))

    def test_claim_and_canonical_mutations_fail_closed(self) -> None:
        source = DEFAULT_CONFIG.read_text(encoding="utf-8")
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bad_config = root / "bad.toml"
            bad_config.write_text(
                source.replace("evolution_authorized = false", "evolution_authorized = true"),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "fail-closed"):
                load_config(bad_config)

            bad_result = root / "bad.json"
            bad_result.write_text(_canonical(self.payload).rstrip(), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "canonical sorted"):
                load_canonical_result(bad_result)


if __name__ == "__main__":
    unittest.main()
