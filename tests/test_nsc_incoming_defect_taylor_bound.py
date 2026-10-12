"""Normalized Taylor remainder, original recurrence and enclosure controls."""
import json
from pathlib import Path

import mpmath as mp
import numpy as np
import pytest

from recursive_horizons.nsc_incoming_defect_taylor_bound import (
    defect_coefficient_jets, centered_coefficient_enclosure, intersect_complex,
    CenteredDefectGeometry,
)
from recursive_horizons.nsc_incoming_vacuum_tail_bound import (
    _precision, _lo, _hi, _range, _geometry_jets, _defect_coefficients,
)
from recursive_horizons.nsc_pg_high_energy import riccati_coefficients

ROOT = Path(__file__).resolve().parents[1]


def _contains(outer, inner):
    return (_lo(outer.real) <= _lo(inner.real) <= _hi(inner.real) <= _hi(outer.real)
            and _lo(outer.imag) <= _lo(inner.imag) <= _hi(inner.imag) <= _hi(outer.imag))


def test_normalized_complex_Taylor_remainder_keeps_factorials():
    # f(x)=(1+2i)*x^4, center0, fourth normalized coefficient=(1+2i).
    # An erroneous second factorial division would fail at either endpoint.
    with _precision(40):
        h = _range(-mp.mpf('0.25'), mp.mpf('0.25'))
        c = mp.iv.mpc(1, 2)
        whole = [c*h**4, 4*c*h**3, 6*c*h*h, 4*c*h, c]
        center = [mp.iv.mpc(0)]*4+[c]
        enclosure = centered_coefficient_enclosure(whole, center, h, 4)
        assert _contains(enclosure, c*mp.iv.mpf('0.25')**4)
        assert _contains(enclosure, mp.iv.mpc(0))
        with pytest.raises(ArithmeticError, match='empty'):
            intersect_complex(mp.iv.mpc(1, 0), mp.iv.mpc(2, 0))


def test_zeroth_defect_matches_previous_interval_owner():
    with _precision(40):
        rho = _range(mp.mpf('1.2'), mp.mpf('1.21'))
        old = _geometry_jets(rho, 16)
        deep = _geometry_jets(rho, 20)
        mass, angular = mp.iv.pi/2, mp.iv.sqrt(5)
        expected = _defect_coefficients(*old, mass, angular)
        actual = defect_coefficient_jets(*deep, mass, angular, depth=4)
        for previous, current in zip(expected, actual):
            # Parenthesization can tighten interval products; their overlap
            # must contain the same point value and not change its meaning.
            intersect_complex(previous, current[0])
        point_geometry = _geometry_jets(mp.iv.mpf('1.205'), 20)
        point = defect_coefficient_jets(*point_geometry, mass, angular, depth=4)
        for interval, jet in zip(actual, point):
            assert _contains(interval[0], jet[0])
        with pytest.raises(ValueError, match='geometry jets'):
            defect_coefficient_jets(*old, mass, angular, depth=4)


def test_centered_enclosure_contains_independent_owned_point_coefficients():
    with _precision(40):
        rho = _range(mp.mpf('1.29'), mp.mpf('1.31')); center = mp.iv.mpf('1.3')
        m, ell = np.pi/2, np.sqrt(5)
        interval = defect_coefficient_jets(*_geometry_jets(rho, 20), mp.iv.mpf(m), mp.iv.mpf(ell), depth=4)
        point = defect_coefficient_jets(*_geometry_jets(center, 20), mp.iv.mpf(m), mp.iv.mpf(ell), depth=4)
        enclosed = [centered_coefficient_enclosure(a, b, rho-center, 4) for a, b in zip(interval, point)]
        for x in (1.29, 1.295, 1.3, 1.305, 1.31):
            c, dc = riccati_coefficients([x], m, ell, 16); c, dc = c[0], dc[0]
            A = 1+3*x+3*(1+x*x)*(np.arctan(x)-np.pi/2)
            Ap = 6+6*x*(np.arctan(x)-np.pi/2); U = ell/np.sqrt(1+x*x)-1j*m
            for p, enclosure in zip(range(16, 33), enclosed):
                value = -A*U*sum(c[j-1]*c[p-j-1] for j in range(1, 17) if 1 <= p-j <= 16)
                if p == 16: value -= 1j*(A*dc[-1]+Ap*c[-1]/2)
                assert _lo(enclosure.real) <= value.real <= _hi(enclosure.real)
                assert _lo(enclosure.imag) <= value.imag <= _hi(enclosure.imag)


def test_centered_cells_cover_same_horizon_and_have_valid_centers():
    config = json.loads((ROOT/'results/development/nsc-compact-matched-restart.json').read_text())['scattering_provenance']['config']
    owner = CenteredDefectGeometry(config, intervals=16, depth=4)
    with _precision(40):
        assert _hi(owner.cells[0][0]) >= _hi(owner.horizon)
        assert _lo(owner.cells[-1][0]) <= 1
        for i, (rho, center, _, _) in enumerate(owner.cells):
            assert _lo(rho) <= _lo(center) <= _hi(center) <= _hi(rho)
            if i: assert _hi(rho) >= _lo(owner.cells[i-1][0])
    with pytest.raises(ValueError, match='Taylor depth'):
        CenteredDefectGeometry(config, depth=16)
