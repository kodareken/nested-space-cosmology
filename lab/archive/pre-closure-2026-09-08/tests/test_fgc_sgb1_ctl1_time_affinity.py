"""Formal highest-time cancellations; no source solve or trajectory.

At fixed lower jets the acceleration part of Riemann has one dt in each
antisymmetric pair. In a normal orthonormal frame it is the arbitrary
symmetric electric block R_0i0j=A_ij. Six independent formal atoms below
cover this whole block, not a list of numerical samples. Tensor covariance
and a timelike dt extend the normal-frame identities to any ADM chart.
"""

from __future__ import annotations

from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "src")]

from recursive_horizons.fgc.spherical_reduction import (  # noqa: E402
    _gb,
    _inv,
    _p,
    _ricci,
    _scalar,
)
from tests.test_fgc_sgb1_ctl1_constraints import _Polynomial as P  # noqa: E402


def _metric():
    return [
        [Q(-1 if a == 0 else 1) if a == b else Q(0) for b in range(4)]
        for a in range(4)
    ]


def _electric_curvature():
    tensor = [[[[P() for _ in range(4)] for _ in range(4)]
               for _ in range(4)] for _ in range(4)]
    for i, j in product(range(1, 4), repeat=2):
        atom = P.atom(f"A{min(i, j)}{max(i, j)}")
        for a, b, first in ((0, i, 1), (i, 0, -1)):
            for c, d, second in ((0, j, 1), (j, 0, -1)):
                tensor[a][b][c][d] = first * second * atom
    return tensor


def _dual_and_gb(riemann, metric):
    inverse = _inv(metric)
    ricci = _ricci(riemann, inverse)
    scalar = _scalar(ricci, inverse)
    return (
        _p(riemann, ricci, scalar, metric),
        _gb(riemann, ricci, scalar, inverse),
        inverse,
    )


def _pullback(riemann, coframe):
    # Sparse tensor pullback: retain the free polynomial atoms and avoid
    # spending time multiplying the many exact zero coframe entries.
    columns = [tuple((a, coframe[a][b]) for a in range(4) if coframe[a][b])
               for b in range(4)]
    transformed = [[[[P() for _ in range(4)] for _ in range(4)]
                    for _ in range(4)] for _ in range(4)]
    for a, b, c, d in product(range(4), repeat=4):
        transformed[a][b][c][d] = sum(
            (x * y * z * w * riemann[i][j][k][m]
             for (i, x), (j, y), (k, z), (m, w)
             in product(columns[a], columns[b], columns[c], columns[d])),
            P(),
        )
    return transformed


class SGBLHighestTimeAffinityTests(unittest.TestCase):
    def test_six_atom_electric_curvature_has_zero_gb_quadratic_form(self):
        riemann = _electric_curvature()
        atoms = {
            atom
            for a, b, c, d in product(range(4), repeat=4)
            for monomial, _ in riemann[a][b][c][d].formal.terms
            for atom in monomial
        }
        self.assertEqual(atoms, {"A11", "A12", "A13", "A22", "A23", "A33"})
        dual, gb, _ = _dual_and_gb(riemann, _metric())
        self.assertEqual(gb.formal.terms, ())
        for a, b in product(range(4), repeat=2):
            self.assertEqual(dual[a][0][b][0], 0)
        # P itself is not zero; only its twice-normal contraction vanishes.
        self.assertEqual(dual[1][2][1][2], -P.atom("A33"))
        scalar_acceleration = P.atom("phi_tt")
        for a, b in product(range(4), repeat=2):
            self.assertEqual(dual[a][0][b][0] * scalar_acceleration, 0)

    def test_nonzero_shift_requires_raised_time_covector_not_coordinate_zero(self):
        # e^0=2dt, e^1=(3/2)(dr+dt/3), e^2=4dx^2, e^3=4dx^3.
        # This gives lapse2, shift1/3, lambda3/2 and det(h)=-9.
        coframe = (
            (Q(2), Q(0), Q(0), Q(0)),
            (Q(1, 2), Q(3, 2), Q(0), Q(0)),
            (Q(0), Q(0), Q(4), Q(0)),
            (Q(0), Q(0), Q(0), Q(4)),
        )
        eta = _metric()
        metric = [
            [sum(eta[i][j] * coframe[i][a] * coframe[j][b]
                 for i, j in product(range(4), repeat=2)) for b in range(4)]
            for a in range(4)
        ]
        riemann = _pullback(_electric_curvature(), coframe)
        dual, gb, inverse = _dual_and_gb(riemann, metric)
        self.assertEqual(gb.formal.terms, ())
        self.assertEqual(inverse[0][0], -Q(1, 4))
        self.assertEqual(inverse[1][0], Q(1, 12))
        self.assertEqual(dual[2][0][2][0], -4 * P.atom("A33"))
        for a, b in product(range(4), repeat=2):
            raised_electric = sum(
                dual[a][c][b][d] * inverse[c][0] * inverse[d][0]
                for c, d in product(range(4), repeat=2)
            )
            self.assertEqual(raised_electric, 0)
            self.assertEqual(raised_electric * P.atom("phi_tt"), 0)


if __name__ == "__main__":
    unittest.main()
