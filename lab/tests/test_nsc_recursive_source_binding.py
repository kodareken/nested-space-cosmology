import unittest
from math import pi, sqrt

from recursive_horizons.nsc_recursive_source_binding import (
    bind_recursive_lll_to_neck,
    inherited_state_constants,
)


class RecursiveSourceBindingTests(unittest.TestCase):
    def test_inherited_state_power_and_neck_sign(self):
        omega = 3.973074368754331
        kappa = 0.23832579963401956
        outgoing, incoming = inherited_state_constants(4, kappa, omega)
        source = bind_recursive_lll_to_neck(
            magnetic_flux=4,
            omega=omega,
            surface_gravity=kappa,
            metric_A=1 - 3*pi/2,
            metric_A_prime=6,
            metric_A_second=-3*pi,
            beta=sqrt(3*pi/2),
            radius=1,
        )
        self.assertAlmostEqual(incoming/outgoing, omega**2, places=12)
        self.assertAlmostEqual(
            source.parent_killing_power, outgoing-incoming, places=12
        )
        self.assertGreater(source.pg_source["fourD_null_plus"], 0)
        self.assertGreater(source.pg_source["fourD_null_minus"], 0)
        self.assertFalse(source.source_null_signs_match)


if __name__ == "__main__":
    unittest.main()

