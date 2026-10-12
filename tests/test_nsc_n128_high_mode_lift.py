import unittest

import numpy as np

from recursive_horizons.nsc_evolved_incoming_constraints import compatible_history_slots
from recursive_horizons.nsc_ks_profile_identity import profile_identity

import derive_nsc_ks_coupled_newton_n64_truncated_iterate as I64
from derive_nsc_ks_n128_high_modes import (
    LIFTED_IDENTITY, LocalIncomingFamily128, family_order, next_metric,
)


class N128HighModeLiftTests(unittest.TestCase):
    def test_zero_pad_keeps_the_n64_history_and_adds_matching_high_slots(self):
        coeff = np.zeros((2, 64))
        coeff[0, :3] = (0.01, -0.002, 0.0004)
        coeff[1, :2] = (-0.02, 0.001)
        low = I64.LocalIncomingFamily64(coeff)
        padded = np.zeros((2, 128))
        padded[:, :64] = coeff
        high = LocalIncomingFamily128(padded)
        z = low.collocation_nodes(16)
        slots_low, tang_low = compatible_history_slots(low.metric(), z, 128)
        slots_high, tang_high = compatible_history_slots(high.metric(), z, 256)
        np.testing.assert_array_equal(slots_low, slots_high)
        np.testing.assert_array_equal(tang_low[:64], tang_high[:64])
        np.testing.assert_array_equal(tang_low[64:], tang_high[128:192])
        metric = next_metric(high)
        _slots, tang_metric = compatible_history_slots(metric, z, 129)
        np.testing.assert_array_equal(tang_metric[1:65], tang_high[64:128])
        np.testing.assert_array_equal(tang_metric[65:], tang_high[192:])
        self.assertEqual(high.radius_lower_bound(), low.radius_lower_bound())
        self.assertGreater(high.radius_lower_bound(), 0.0)

    def test_family_order_interleaves_reused_and_new_blocks(self):
        low = np.zeros((128, 2, 2))
        extra = np.zeros((128, 2, 2))
        low[0, 0, 0] = 1.0
        low[64, 1, 1] = 2.0
        extra[0, 0, 1] = 3.0
        extra[64, 1, 0] = 4.0
        full = family_order(low, extra)
        self.assertEqual(full.shape, (256, 2, 2))
        self.assertEqual(full[0, 0, 0], 1.0)
        self.assertEqual(full[64, 0, 1], 3.0)
        self.assertEqual(full[128, 1, 1], 2.0)
        self.assertEqual(full[192, 1, 0], 4.0)

    def test_base_allowlist_still_rejects_128_and_the_lift_identity_is_stable(self):
        with self.assertRaises(ValueError):
            I64.LocalIncomingFamily(np.zeros((2, 128)))
        with self.assertRaises(ValueError):
            LocalIncomingFamily128(np.zeros((2, 64)))
        padded = np.zeros((2, 128))
        self.assertNotEqual(
            profile_identity(LocalIncomingFamily128(padded), include_normal_window=True),
            LIFTED_IDENTITY)


if __name__ == '__main__':
    unittest.main()
