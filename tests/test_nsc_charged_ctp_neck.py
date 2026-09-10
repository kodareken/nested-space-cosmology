import unittest

from recursive_horizons.nsc_charged_ctp_neck import (
    ChargedCTPNeckConfig,
    _magnetic_levels,
)


class ChargedCTPNeckTests(unittest.TestCase):
    def test_q4_magnetic_spectrum_and_configuration(self):
        masses, degeneracies = _magnetic_levels(4, 3)
        self.assertEqual(degeneracies.tolist(), [12.0, 16.0, 20.0])
        self.assertAlmostEqual(masses[0]**2, 5.0)
        self.assertAlmostEqual(masses[1]**2, 12.0)
        self.assertAlmostEqual(masses[2]**2, 21.0)
        with self.assertRaises(ValueError):
            ChargedCTPNeckConfig(4, 1, 1, 1, evolution_steps=20)


if __name__ == "__main__":
    unittest.main()

