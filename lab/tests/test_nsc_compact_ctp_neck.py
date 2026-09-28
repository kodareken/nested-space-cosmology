import unittest
from math import pi

import numpy as np

from recursive_horizons.nsc_angular_stress import (
    adiabatic_coefficients,
    profile_jets,
)
from recursive_horizons.nsc_compact_ctp_neck import (
    _su2_step,
    adiabatic_bloch,
)


class CompactCTPNeckTests(unittest.TestCase):
    def test_three_axis_step_is_unitary(self):
        upper = np.array([0.7+0.2j, -0.1+0.3j])
        lower = np.array([0.2-0.4j, 0.6+0.1j])
        before = abs(upper)**2+abs(lower)**2
        changed = _su2_step(upper, lower, 0.4, -0.3, 0.8, 0.17)
        np.testing.assert_allclose(
            abs(changed[0])**2+abs(changed[1])**2, before,
            atol=2e-15, rtol=0,
        )

    def test_massless_limit_reproduces_existing_adiabatic_subtraction(self):
        frequency = 1.0
        jets = profile_jets(pi/2)
        momentum = -frequency/np.sqrt(jets[0])
        bloch = adiabatic_bloch(pi/2, 0.0, 1.0, frequency)
        leading = np.array([0.0, 1.0, momentum])
        leading /= np.linalg.norm(leading)
        correction = bloch+leading
        energy_projection = np.dot(
            np.array([0.0, 1.0, momentum]), correction
        )
        pressure_projection = momentum*correction[2]
        e2, p2, e4, p4 = adiabatic_coefficients(
            np.array([[momentum]]), 1.0, jets
        )
        self.assertAlmostEqual(
            energy_projection, float((e2+e4)[0, 0]), places=13
        )
        self.assertAlmostEqual(
            pressure_projection, float((p2+p4)[0, 0]), places=13
        )


if __name__ == "__main__":
    unittest.main()
