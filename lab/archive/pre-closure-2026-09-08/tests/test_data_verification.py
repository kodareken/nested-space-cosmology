"""Focused checks for the collected-data verification contract."""

from __future__ import annotations

import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("verify_data", REPOSITORY / "scripts" / "verify_data.py")
assert SPEC is not None and SPEC.loader is not None
VERIFY = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(VERIFY)


class DataVerificationTests(unittest.TestCase):
    def test_tracked_artifacts_verify_with_semantic_contracts(self) -> None:
        passes, failures, verified, skipped = VERIFY.verify_manifest(REPOSITORY / "collected-data")
        self.assertFalse(failures, failures)
        optional_planck = ("smica_2048.fits", "mask_common.fits")
        optional_present = sum(
            (REPOSITORY / "collected-data" / filename).is_file()
            for filename in optional_planck
        )
        self.assertEqual((verified, skipped), (9 + optional_present, 2 - optional_present))
        self.assertTrue(any("DESI DR2 BAO 13-element mean vector" in item for item in passes))
        self.assertTrue(any("Fermi GBM trigger-entry FITS" in item for item in passes))
        self.assertFalse(any("expected FITS signature" in item for item in passes))

    def test_optional_files_are_skipped_unless_require_is_set(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest = root / "manifest.csv"
            with manifest.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=["local_filename", "bytes", "sha256", "source_url", "required_in_clone"])
                writer.writeheader()
                writer.writerow({"local_filename": "optional.txt", "bytes": "1", "sha256": "0" * 64, "source_url": "https://example.invalid/optional.txt", "required_in_clone": "no"})
            _, failures, _, skipped = VERIFY.verify_manifest(root, manifest=manifest)
            self.assertFalse(failures)
            self.assertEqual(skipped, 1)
            _, failures, _, _ = VERIFY.verify_manifest(root, require_all=True, manifest=manifest)
            self.assertEqual(failures, ["optional.txt: required file is absent"])

    def test_desi_covariance_wrong_shape_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "desi_gaussian_bao_ALL_GCcomb_cov.txt"
            path.write_text("1 0\n0 1\n", encoding="utf-8")
            label, failures = VERIFY.verify_semantics(path)
            self.assertEqual(label, "DESI DR2 BAO 13 by 13 covariance")
            self.assertTrue(failures)

    def test_fermi_signature_catches_wrong_primary_identity(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "glg_tcat_all_bn170817529_v03.fit"
            path.write_bytes(b"not a FITS file")
            label, failures = VERIFY.verify_fits(path)
            self.assertEqual(label, "Fermi GBM trigger-entry FITS")
            self.assertTrue(failures)


if __name__ == "__main__":
    unittest.main()
