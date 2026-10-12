"""Required-jet truncation, normalized remainder and retained-domain checks."""
import json
from pathlib import Path

import mpmath as mp
import pytest

from recursive_horizons.nsc_incoming_centered_order24 import trimmed_defect_jets, CenteredOrder24Geometry
from recursive_horizons.nsc_incoming_defect_taylor_bound import (
    defect_coefficient_jets, centered_coefficient_enclosure, intersect_complex,
)
from recursive_horizons.nsc_incoming_vacuum_tail_bound import (
    _precision, _lo, _hi, _range, _geometry_jets, _defect_coefficients,
)

ROOT = Path(__file__).resolve().parents[1]


def contains(a, b):
    return all(_lo(x) <= _lo(y) <= _hi(y) <= _hi(x)
               for x, y in ((a.real, b.real), (a.imag, b.imag)))


def test_trimmed16_all_required_jets_equal_existing_untrimmed_owner():
    with _precision(40):
        geometry = _geometry_jets(_range(mp.mpf('1.2'), mp.mpf('1.21')), 20)
        args = (*geometry, mp.iv.pi/2, mp.iv.sqrt(5))
        old = defect_coefficient_jets(*args, depth=4)
        new = trimmed_defect_jets(*args, order=16)
        assert [[v._mpci_ for v in row] for row in new] == [[v._mpci_ for v in row] for row in old]


def test_trimmed24_zeroth_matches_generic_and_minimum_depth_is_sufficient():
    with _precision(40):
        geometry = _geometry_jets(_range(mp.mpf('1.25'), mp.mpf('1.251')), 28)
        args = (*geometry, mp.iv.pi/2, mp.iv.sqrt(5))
        new = trimmed_defect_jets(*args)
        generic = _defect_coefficients(*args, order=24)
        for current, old in zip(new, generic):
            assert current[0]._mpci_ == old._mpci_
        minimum = trimmed_defect_jets(geometry[0][:29], geometry[1][:29], *args[2:])
        assert [[v._mpci_ for v in row] for row in new] == [[v._mpci_ for v in row] for row in minimum]
        with pytest.raises(ValueError, match='geometry'):
            trimmed_defect_jets(geometry[0][:28], geometry[1][:28], *args[2:])


def test_complex_remainder_factorial_and_intersection():
    with _precision(40):
        h = _range(mp.mpf('-0.2'), mp.mpf('0.3')); c = mp.iv.mpc(2, -3)
        whole = [c*h**4, 4*c*h**3, 6*c*h*h, 4*c*h, c]
        point = [mp.iv.mpc(0)]*4+[c]
        result = centered_coefficient_enclosure(whole, point, h, 4)
        assert contains(whole[0], result)
        for x in ('-0.2', '0', '0.3'):
            assert contains(result, c*mp.iv.mpf(mp.mpf(x))**4)
        with pytest.raises(ArithmeticError, match='empty'):
            intersect_complex(mp.iv.mpc(1), mp.iv.mpc(2))


def test_same128_cells_cover_collar_and_centers():
    config = json.loads((ROOT/'results/development/nsc-compact-matched-restart.json').read_text())['scattering_provenance']['config']
    owner = CenteredOrder24Geometry(config)
    assert len(owner.cells) == 128
    with _precision(40):
        assert _hi(owner.cells[0][0]) >= _hi(owner.horizon)
        assert _lo(owner.cells[-1][0]) <= 1
        for i, (rho, center) in enumerate(owner.cells):
            assert _lo(rho) <= _lo(center) <= _hi(center) <= _hi(rho)
            if i: assert _hi(rho) >= _lo(owner.cells[i-1][0])
    with pytest.raises(ValueError, match='group14'):
        owner.bound_channel({'index': 12})


def test_no_unrequested_approximation_or_remainder_depth():
    with pytest.raises(ValueError, match='order16'):
        trimmed_defect_jets([], [], 1, 1, order=32)
    with pytest.raises(ValueError, match='depth4'):
        trimmed_defect_jets([], [], 1, 1, depth=8)
