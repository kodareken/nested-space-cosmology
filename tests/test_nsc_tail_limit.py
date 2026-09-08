import unittest
import mpmath as mp
from recursive_horizons.nsc_tail_limit import finite_response,weyl_disk,motif_hoppings


class TailLimitTests(unittest.TestCase):
    def test_passive_sign_is_necessary_for_disk_bound(self):
        with mp.workdps(60):
            z=mp.mpc('.3','.4');center,radius,_=weyl_disk([],z)
            passive=finite_response([],z,-mp.j)
            active=finite_response([],z,mp.j)
            self.assertLess(abs(passive-center),radius)
            self.assertGreater(abs(active-center),radius)

    def test_weyl_disk_matches_free_two_site_end_potentials(self):
        with mp.workdps(60):
            z=mp.j;center,radius,_=weyl_disk([mp.mpf(1)],z)
            self.assertLess(abs(center+mp.mpc(0,.75)),mp.mpf('1e-50'))
            self.assertEqual(radius,mp.mpf('.25'))
            for terminal in (-2,0,2):
                self.assertLess(abs(abs(finite_response([mp.mpf(1)],z,terminal)-center)-radius),mp.mpf('1e-50'))

    def test_smooth_chain_zero_energy_threshold_has_both_regimes(self):
        motif,base=motif_hoppings()
        q=motif['zero_energy_transfer']
        self.assertTrue(all(base>0))
        self.assertLess(4*q,1)
        self.assertGreater(50*q,1)


if __name__=='__main__':unittest.main()
