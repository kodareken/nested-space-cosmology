from __future__ import annotations

import json
import math
import sys
import tempfile
import unittest
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

from recursive_horizons.bao import (  # noqa: E402
    C_KM_S,
    BAODataset,
    BAOPoint,
    analytic_amplitude_profile,
    default_data_directory,
    e_squared,
    inverse_e_integral,
    load_bao_dataset,
    reproduce_desi_dr2_bao,
)


class BAODataTests(unittest.TestCase):
    def test_public_vector_covariance_and_profile_are_valid(self) -> None:
        directory = default_data_directory()
        dataset = load_bao_dataset(
            directory / "desi_gaussian_bao_ALL_GCcomb_mean.txt",
            directory / "desi_gaussian_bao_ALL_GCcomb_cov.txt",
        )
        self.assertEqual(len(dataset.points), 13)
        self.assertEqual(len(dataset.covariance), 13)
        self.assertEqual(dataset.points[0].observable, "DV_over_rs")
        self.assertEqual(dataset.points[-1].observable, "DM_over_rs")
        profile = analytic_amplitude_profile(dataset, 0.2976935, -0.9118992)
        self.assertGreater(profile.amplitude_c_over_h0rd, 0.0)
        self.assertAlmostEqual(profile.h0_rd_km_s * profile.amplitude_c_over_h0rd, C_KM_S)

    def test_invalid_labels_and_indefinite_covariance_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            mean = root / "mean.txt"
            covariance = root / "cov.txt"
            mean.write_text("0.3 1.0 not_an_observable\n" * 13, encoding="utf-8")
            covariance.write_text("\n".join(" ".join("1" if i == j else "0" for j in range(13)) for i in range(13)), encoding="utf-8")
            with self.assertRaises(ValueError):
                load_bao_dataset(mean, covariance)
            mean.write_text("\n".join(f"{0.1 + i / 10:.1f} 1.0 DM_over_rs" for i in range(13)), encoding="utf-8")
            covariance.write_text("\n".join(" ".join("1" for _ in range(13)) for _ in range(13)), encoding="utf-8")
            with self.assertRaises(ValueError):
                load_bao_dataset(mean, covariance)


class BAOBackgroundTests(unittest.TestCase):
    def test_flat_constant_w_background_and_simpson_refinement(self) -> None:
        self.assertAlmostEqual(e_squared(0.0, 0.3, -1.0), 1.0)
        self.assertAlmostEqual(e_squared(1.0, 0.3, -1.0), 3.1)
        coarse = inverse_e_integral(2.33, 0.3, -0.9, 128)
        fine = inverse_e_integral(2.33, 0.3, -0.9, 256)
        self.assertLess(abs(fine - coarse), 1.0e-9)
        with self.assertRaises(ValueError):
            e_squared(-0.1, 0.3, -1.0)
        with self.assertRaises(ValueError):
            e_squared(1.0, True, -1.0)  # type: ignore[arg-type]
        with self.assertRaises(ValueError):
            e_squared(1.0, 0.3, True)  # type: ignore[arg-type]
        with self.assertRaises(ValueError):
            e_squared(1.0, 0.3, "-1")  # type: ignore[arg-type]
        with self.assertRaises(ValueError):
            inverse_e_integral(1.0, 0.3, -1.0, 127)


class BAOArtifactTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.record = reproduce_desi_dr2_bao()

    def test_profile_result_has_expected_robust_ranges(self) -> None:
        best = self.record["best_fit"]
        self.assertTrue(-0.93 < best["w"] < -0.90)
        self.assertTrue(0.29 < best["omega_m"] < 0.31)
        self.assertTrue(9900.0 < best["H0rd_km_s"] < 10050.0)
        self.assertTrue(8.5 < best["chi2"] < 9.6)
        candidates = self.record["fixed_w_candidates"]
        self.assertTrue(1.0 < candidates["-1.0"]["delta_chi2_from_best"] < 1.5)
        self.assertTrue(10.0 < candidates["-0.6666666666666666"]["delta_chi2_from_best"] < 12.5)
        self.assertTrue(65.0 < candidates["-0.3333333333333333"]["delta_chi2_from_best"] < 75.0)

    def test_profile_is_converged_and_json_safe(self) -> None:
        refinement = self.record["convergence_refinement"]
        self.assertLess(refinement["abs_best_w_difference"], 2.0e-6)
        self.assertLess(refinement["abs_best_chi2_difference"], 1.0e-6)
        self.assertEqual(self.record["dof"], 10)
        self.assertEqual(self.record["project_version"], "0.11.0")
        self.assertIn("declared_search", self.record)
        self.assertNotIn("preregistered_search", self.record)
        self.assertEqual(len(self.record["data_points"]), 13)
        self.assertIn("not DESI+CMB/SN", self.record["scope"])
        json.dumps(self.record, allow_nan=False)


if __name__ == "__main__":
    unittest.main()
