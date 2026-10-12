from __future__ import annotations

from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

sys.path[:0] = [str(Path(__file__).resolve().parents[1]), str(Path(__file__).resolve().parents[1] / "src")]

from scripts.reproduce_fgc_hyp1_fo1_rc1 import load_config
from recursive_horizons.fgc.modified_harmonic_constraints import activated_compatible_state
from recursive_horizons.fgc.modified_harmonic_uniform_domain import uniform_radial_hyperbolicity_certificate


class UHYP1Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cfg = load_config()
        cls.datum = activated_compatible_state(cfg.flat_fixture)

    def test_nonzero_compact_radial_frame(self) -> None:
        out = uniform_radial_hyperbolicity_certificate(self.datum["state"], reference=self.datum["reference"])
        self.assertTrue(out["domain"]["all_30_parameter_axes_nonzero"])
        self.assertTrue(out["domain"]["all_6_acceleration_axes_nonzero"])
        self.assertTrue(out["eigenframe"]["all_enclosed_frames_invertible"])
        self.assertLess(out["eigenframe"]["neumann_inverse"]["rho_infinity"], 1)
        self.assertEqual(len(out["eigenframe"]["interval_columns"]), 12)
        self.assertEqual(out["eigenframe"]["neumann_inverse"]["dimension"], 12)
        for mode in out["modes"]["hat"]:
            self.assertLess(mode["linear_pivot_krawczyk"]["rho_infinity_upper_bound"], 1)
            self.assertTrue(mode["linear_pivot_krawczyk"]["krawczyk_image_strictly_inside_displacement_box"])
        for mode, physical in zip(out["modes"]["regulator"], out["cone_quadratics"]["root_brackets"]["physical"], strict=True):
            self.assertLess(mode["nonlinear_full_five_row_krawczyk"]["rho_infinity_upper_bound"], 1)
            self.assertTrue(mode["nonlinear_full_five_row_krawczyk"]["krawczyk_image_strictly_inside_displacement_box"])
            self.assertTrue(mode["speed_box"].upper < physical["root_interval"].lower or physical["root_interval"].upper < mode["speed_box"].lower)
        routes = out["modes"]["correlated_eigenmode_routes"]
        self.assertEqual(set(routes), {"tilde", "hat", "physical_chi", "regulator"})
        self.assertTrue(all("interval-residual containment" not in route or "never" in route for route in routes.values()))
        sym = out["radial_symmetrizer"]
        self.assertGreater(sym["euclidean_coercivity_lower_bound"], 0)
        self.assertGreater(sym["euclidean_coercivity_upper_bound"], sym["euclidean_coercivity_lower_bound"])
        self.assertTrue(all(value is False for value in out["nonclaims"].values()))

    def test_widths_fail_closed(self) -> None:
        with self.assertRaises(ValueError):
            uniform_radial_hyperbolicity_certificate(self.datum["state"], reference=self.datum["reference"], parameter_half_width=Q(0))
        with self.assertRaises(ValueError):
            uniform_radial_hyperbolicity_certificate(self.datum["state"], reference=self.datum["reference"], root_bracket_half_width=Q(0))


if __name__ == "__main__":
    unittest.main()
