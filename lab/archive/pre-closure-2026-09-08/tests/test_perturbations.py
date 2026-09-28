from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path

REPOSITORY = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY / "src"))

from recursive_horizons.perturbations import (  # noqa: E402
    closed_desitter_scale_factor,
    core_transfer_matrix,
    scale_factor_potential,
    spectator_frequency_squared,
    tensor_frequency_squared,
    transfer_convergence,
)


class ClosedDeSitterFormulaTests(unittest.TestCase):
    def test_scale_factor_and_potential(self) -> None:
        eta = 0.4
        self.assertAlmostEqual(
            closed_desitter_scale_factor(eta, 3.0), 3.0 / math.cos(eta), places=14
        )
        self.assertAlmostEqual(
            scale_factor_potential(eta), 2.0 / math.cos(eta) ** 2 - 1.0, places=14
        )

    def test_scalar_frequency_formula(self) -> None:
        eta = 0.25
        mass_radius = 1.7
        expected = (4 + 1) ** 2 + (mass_radius**2 - 2.0) / math.cos(eta) ** 2
        self.assertAlmostEqual(
            spectator_frequency_squared(4, eta, mass_radius), expected, places=13
        )

    def test_tensor_frequency_formula(self) -> None:
        eta = 0.25
        self.assertAlmostEqual(
            tensor_frequency_squared(5, eta), 5**2 - 2.0 / math.cos(eta) ** 2, places=13
        )


class TransferTests(unittest.TestCase):
    def test_zero_width_transfer_is_identity(self) -> None:
        matrix = core_transfer_matrix(
            sector="tensor", n=3, eta_core=0.0, steps=1
        )
        self.assertEqual((matrix.m11, matrix.m12, matrix.m21, matrix.m22), (1.0, 0.0, 0.0, 1.0))

    def test_scalar_constant_frequency_has_known_rotation(self) -> None:
        # mL=sqrt(2) cancels the sec^2 term, leaving omega=n+1 exactly.
        n = 2
        eta_core = 0.7
        omega = n + 1
        elapsed = 2.0 * eta_core
        matrix = core_transfer_matrix(
            sector="scalar",
            n=n,
            eta_core=eta_core,
            mass_radius=math.sqrt(2.0),
            steps=4_096,
        )
        self.assertAlmostEqual(matrix.m11, math.cos(omega * elapsed), places=12)
        self.assertAlmostEqual(matrix.m12, math.sin(omega * elapsed) / omega, places=12)
        self.assertAlmostEqual(matrix.m21, -omega * math.sin(omega * elapsed), places=12)
        self.assertAlmostEqual(matrix.m22, math.cos(omega * elapsed), places=12)

    def test_wronskian_is_preserved(self) -> None:
        for sector, n in (("scalar", 1), ("tensor", 3)):
            with self.subTest(sector=sector):
                matrix = core_transfer_matrix(
                    sector=sector, n=n, eta_core=0.8, steps=4_096
                )
                self.assertAlmostEqual(matrix.determinant, 1.0, places=10)

    def test_refinement_converges(self) -> None:
        records = transfer_convergence(
            sector="tensor", n=4, eta_core=0.7, initial_steps=128, refinements=4
        )
        self.assertEqual([record.steps for record in records], [128, 256, 512, 1024])
        self.assertLess(records[-1].determinant_error, 1e-10)
        self.assertIsNotNone(records[1].difference_from_previous)
        self.assertLess(
            float(records[-1].difference_from_previous),
            float(records[1].difference_from_previous),
        )


class ValidationTests(unittest.TestCase):
    def test_invalid_harmonics_and_domains_are_rejected(self) -> None:
        with self.assertRaises(ValueError):
            spectator_frequency_squared(-1, 0.0)
        with self.assertRaises(ValueError):
            tensor_frequency_squared(2, 0.0)
        with self.assertRaises(ValueError):
            closed_desitter_scale_factor(math.pi / 2.0)
        with self.assertRaises(ValueError):
            core_transfer_matrix(sector="vector", n=1, eta_core=0.2)
        with self.assertRaises(ValueError):
            core_transfer_matrix(sector="tensor", n=3, eta_core=0.2, steps=0)
        with self.assertRaises(ValueError):
            transfer_convergence(
                sector="scalar", n=0, eta_core=0.2, refinements=0
            )


if __name__ == "__main__":
    unittest.main()
