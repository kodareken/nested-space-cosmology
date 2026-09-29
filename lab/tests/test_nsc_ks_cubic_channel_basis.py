"""Shared cubic channels: manufactured integrals, group14 jets, one stored pilot."""
import json
import time
from pathlib import Path

import numpy as np
import pytest
from flint import arb, arb_series, ctx

from recursive_horizons import nsc_ks_ball_geometry as G
from recursive_horizons.nsc_ks_ball_trajectory import restored_upper
from recursive_horizons.nsc_ks_cubic_channel_basis import (
    CubicChannelBudgetExceeded,
    assemble_cubic_channels,
    enclose_cubic_channel_basis,
    integrate_shared_cell,
    reconstructed_massive_forcing,
    split_massless_channels,
    unit_cubic_geometry,
)
from recursive_horizons.nsc_ks_cubic_uv_enclosure import CubicGeometry
from recursive_horizons.nsc_ks_current_uv_transport import ARCHIVE_RHO_UP
from recursive_horizons.nsc_ks_massive_cubic_uv import (
    MassiveCubicGeometry,
    enclose_massive_characteristic,
    original_group14_labels,
    require_original_history,
)
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily


ROOT = Path(__file__).resolve().parents[1]
RECORD = ROOT / 'results/development/nsc-ks-massive-cubic-uv-v1.json'
COEFFICIENTS = ('q_z', 'q_zz', 'J3', 'N_bracket3', 'surface_mass_correction')


def _series(coefficients, order):
    return arb_series([arb(value) for value in coefficients], prec=order + 1)


def _overlaps_stored(ball, item):
    with ctx.workprec(256):
        lower = restored_upper(item['lower'])
        upper = restored_upper(item['upper'])
    if ball.upper() < lower or ball.lower() > upper:
        return False
    return True


def _display(ball):
    return {'mid': float(ball.mid()), 'rad': float(ball.rad()), 'ball': str(ball)}


def test_manufactured_scalar_and_series_cells_match_the_integrals():
    """Constant and polynomial forcings, both characteristic signs.

    The Taylor cell is exact for these polynomials. The uniform cell only has
    to contain that value: its tube is wider than the nested integral. A frozen
    end-of-cell Qzz is not the integral.
    """
    with ctx.workprec(160):
        order = 4
        expected = {
            1: (arb(-1), arb(-2), arb(-4), arb(0), arb(0)),
            -1: (arb(-1), arb(-2), arb(-10), arb(6), arb(-2)),
        }
        frozen_wrong_j4 = {1: arb(3), -1: arb(9)}
        for sign, exact in expected.items():
            scalar = split_massless_channels(
                arb(1), arb(2), arb(7), arb(3), arb(1), sign, 0)
            assert scalar['B'] == (3 if sign == 1 else -3)
            assert scalar['A'] == (4 if sign == 1 else 10)
            point = split_massless_channels(
                _series([1], order), _series([2], order), _series([7], order),
                _series([3], order), _series([1], order), sign, order)
            tight = integrate_shared_cell(
                (0, 0, 0, 0, 0), point, point, arb(1), arb(0), order)
            wide = integrate_shared_cell(
                (0, 0, 0, 0, 0), scalar, None, arb(1), arb(0), 0)
            for got, want in zip(tight, exact):
                assert got.contains(want)
                assert got.rad() < arb('1e-30')
            for got, want in zip(wide, exact):
                assert got.contains(want)
            for left, right in zip(tight, wide):
                assert right.contains(left)
            assert not tight[3].contains(frozen_wrong_j4[sign])
        variable = split_massless_channels(
            _series([1], order), _series([0], order), _series([0], order),
            _series([0], order), _series([1, 1], order), 1, order)
        stepped = integrate_shared_cell(
            (0, 0, 0, 0, 0), variable, variable, arb(1), arb(0), order)
        assert stepped[0].contains(-1) and stepped[0].rad() < arb('1e-30')
        assert stepped[4].contains(-arb(1) / 3)
        assert stepped[4].rad() < arb('1e-30')
        for index in (1, 2, 3):
            assert stepped[index].contains(0)
            assert stepped[index].rad() < arb('1e-30')


def test_assembly_squares_mass_and_is_even_in_angular():
    with ctx.workprec(160):
        basis = {
            'Qz': arb(2), 'Qzz': arb(3), 'J2': arb(5), 'J4': arb(7), 'JM': arb(11),
            'N2': arb(13), 'N4': arb(17), 'NM': arb(19), 'surface_a': arb(23),
            'sign': -1, 'bits': 160, 'massless_A': True, 'C_M': None,
            'higher_uv_remainder_C_M': None, 'entire_incoming_interval': False,
        }
        got = assemble_cubic_channels(basis, 3, 2)
        assert got['q_z'].contains(18) and got['q_zz'].contains(27)
        assert got['J3'].contains(1008) and got['N_bracket3'].contains(2178)
        # -s a m^2 qz = -(-1)*23*4*18
        assert got['surface_mass_correction'].contains(1656)
        for key in ('q_z', 'q_zz', 'J3', 'N_bracket3', 'surface_mass_correction'):
            assert got[key].rad() < arb('1e-30')
        flipped = assemble_cubic_channels(basis, -3, -2)
        for key in ('q_z', 'q_zz', 'J3', 'N_bracket3', 'surface_mass_correction'):
            assert flipped[key].overlaps(got[key])
        assert got['C_M'] is None and got['entire_incoming_interval'] is False
        assert got['physical_local_gate'] == 'OPEN'
        refused = dict(basis)
        refused['C_M'] = arb(1)
        with pytest.raises(ValueError, match='C_M'):
            assemble_cubic_channels(refused, 3, 2)
        refused = dict(basis)
        refused['entire_incoming_interval'] = True
        with pytest.raises(ValueError, match='incoming interval'):
            assemble_cubic_channels(refused, 3, 2)
        refused = dict(basis)
        refused['massless_A'] = False
        with pytest.raises(ValueError, match='massless'):
            assemble_cubic_channels(refused, 3, 2)


def test_massive_geometry_nonunit_angular_and_expired_budget_emit_nothing():
    family, _ = require_original_history()
    before = (ctx.prec, ctx.cap)
    model = unit_cubic_geometry(family)
    incoming = arb(family.center)
    with pytest.raises(TypeError, match='double-counts'):
        enclose_cubic_channel_basis(
            MassiveCubicGeometry(family, 1, arb.pi() / 2, bits=160),
            incoming, 1, ARCHIVE_RHO_UP, cells=1, order=1)
    with pytest.raises(ValueError, match='unit angular'):
        enclose_cubic_channel_basis(
            CubicGeometry(family, arb(5).sqrt(), bits=160),
            incoming, 1, ARCHIVE_RHO_UP, cells=1, order=1)
    with pytest.raises(CubicChannelBudgetExceeded, match='no enclosure returned'):
        enclose_cubic_channel_basis(
            model, incoming, 1, ARCHIVE_RHO_UP, cells=1024, order=8,
            cpu_deadline=time.process_time() - 1)
    assert (ctx.prec, ctx.cap) == before


def test_subdivision_failure_emits_no_enclosure_and_restores_context():
    family, _ = require_original_history()
    model = unit_cubic_geometry(family)
    before = (ctx.prec, ctx.cap)
    original = CubicGeometry.forcing_series

    def fail(self, rho, z_in, sign, order):
        raise G.SubdivisionNeeded('positive bivariate radius not enclosed')

    CubicGeometry.forcing_series = fail
    try:
        with pytest.raises(G.SubdivisionNeeded, match='no enclosure emitted'):
            enclose_cubic_channel_basis(
                model, arb(family.center), -1, ARCHIVE_RHO_UP,
                cells=1, max_depth=0, order=1)
    finally:
        CubicGeometry.forcing_series = original
    assert (ctx.prec, ctx.cap) == before


def test_domain_rules_z_box_and_zero_history_are_not_a_source():
    family, identity = require_original_history()
    model = unit_cubic_geometry(family)
    incoming = arb(family.center)
    before = (ctx.prec, ctx.cap)
    with pytest.raises(ValueError):
        enclose_cubic_channel_basis(
            model, incoming, 0, ARCHIVE_RHO_UP, cells=2, order=1)
    with pytest.raises(ValueError, match='exact upstream'):
        enclose_cubic_channel_basis(
            model, incoming, 1, arb('1.02 +/- 0.001'), cells=2, order=1)
    with pytest.raises(ValueError, match='flat'):
        enclose_cubic_channel_basis(model, incoming, 1, 1.005, cells=2, order=1)
    with pytest.raises(ValueError, match='cell count'):
        enclose_cubic_channel_basis(
            model, incoming, 1, ARCHIVE_RHO_UP, cells=True, order=1)
    with pytest.raises(ValueError, match='Taylor order'):
        enclose_cubic_channel_basis(
            model, incoming, 1, ARCHIVE_RHO_UP, cells=2, order=14)
    zero = unit_cubic_geometry(LocalIncomingFamily(np.zeros((2, 8))))
    center = arb(zero.family.center)
    box = center.union(center + arb(10)**-6)
    for sign in (1, -1):
        basis = enclose_cubic_channel_basis(
            zero, box if sign == 1 else center, sign, ARCHIVE_RHO_UP,
            cells=2, order=1)
        assert basis['profile_identity']
        assert basis['entire_incoming_interval'] is False
        assert basis['z_box_is_full_incoming_interval'] is False
        assert basis['C_M'] is None
        assert basis['higher_uv_remainder_bound'] is None
        assert basis['higher_uv_remainder_C_M'] is None
        assert basis['uniform_C4_on_I'] is None
        assert basis['complete_UV_tail'] is None
        assert basis['physical_local_gate'] == 'OPEN'
        assert basis['unit_angular_is_source_species'] is False
        assert basis['massless_A'] is True
        assert basis['accepted_cells'] >= 2
        if sign == 1:
            assert basis['z_domain'].rad() > 0
        assembled = assemble_cubic_channels(basis, arb(5).sqrt(), arb.pi() / 2)
        for key in ('Qz', 'Qzz', 'J2', 'J4', 'JM', 'N2', 'N4', 'NM'):
            assert basis[key].contains(0)
        for key in COEFFICIENTS:
            assert assembled[key].contains(0)
            assert assembled[key].is_finite()
        assert assembled['entire_incoming_interval'] is False
        assert assembled['C_M'] is None
    assert identity.startswith('0b0e4ced')
    assert (ctx.prec, ctx.cap) == before


def test_group14_forcing_jets_match_and_massive_A_double_counts():
    family, _ = require_original_history()
    before = (ctx.prec, ctx.cap)
    with ctx.workprec(160):
        model = unit_cubic_geometry(family)
        angular, mass = arb(5).sqrt(), arb.pi() / 2
        massive = MassiveCubicGeometry(family, angular, mass, bits=160)
        negative = MassiveCubicGeometry(family, -angular, mass, bits=160)
        rho, incoming = arb('1.018'), arb(family.center)
        couplings = []
        for sign in (1, -1):
            for order in (0, 4):
                got = reconstructed_massive_forcing(
                    model, rho, incoming, sign, order, angular, mass)
                if order == 0:
                    raw = massive.forcing(rho, incoming, sign)
                    even = negative.forcing(rho, incoming, sign)
                    for name in vars(raw):
                        assert getattr(got, name).overlaps(getattr(raw, name))
                        assert getattr(raw, name).overlaps(getattr(even, name))
                    couplings.append(raw.coupling)
                else:
                    raw = massive.forcing_series(rho, incoming, sign, order)
                    even = negative.forcing_series(rho, incoming, sign, order)
                    for name in raw:
                        for left, right in zip(
                                G._coeffs(got[name], order), G._coeffs(raw[name], order)):
                            assert (left - right).contains(0)
                        for left, right in zip(
                                G._coeffs(raw[name], order), G._coeffs(even[name], order)):
                            assert (left - right).contains(0)
            axial = G.background_series(rho, 0).a[0]
            base = CubicGeometry.forcing(model, rho, incoming, sign)
            extra = mass**2 * axial**2 * base.qz
            bad = split_massless_channels(
                base.qz, base.qzz, base.current + extra, base.coupling, axial, sign, 0)
            good = split_massless_channels(
                base.qz, base.qzz, base.current, base.coupling, axial, sign, 0)
            assert (bad['A'] - good['A'] - extra).contains(0)
            assert not extra.contains(0)
            correct = reconstructed_massive_forcing(
                model, rho, incoming, sign, 0, angular, mass)
            contaminated = correct.current + angular**2 * extra
            raw = massive.forcing(rho, incoming, sign)
            assert (contaminated - raw.current - angular**2 * extra).contains(0)
            assert not (contaminated - raw.current).contains(0)
        assert not couplings[0].overlaps(couplings[1])
    assert (ctx.prec, ctx.cap) == before


def test_short_path_both_signs_overlap_the_direct_massive_integrator():
    family, identity = require_original_history()
    with ctx.workprec(160):
        model = unit_cubic_geometry(family)
        angular, mass = arb(5).sqrt(), arb.pi() / 2
        massive = MassiveCubicGeometry(family, angular, mass, bits=160)
        incoming = arb(family.center)
        brackets = []
        for sign in (1, -1):
            basis = enclose_cubic_channel_basis(
                model, incoming, sign, ARCHIVE_RHO_UP, cells=8, order=3)
            assert basis['profile_identity'] == identity
            assert basis['massless_A'] is True
            assert basis['C_M'] is None
            assert basis['entire_incoming_interval'] is False
            got = assemble_cubic_channels(basis, angular, mass)
            direct = enclose_massive_characteristic(
                massive, incoming, sign, ARCHIVE_RHO_UP, cells=8, order=3)
            for key in COEFFICIENTS:
                assert got[key].overlaps(direct[key]), (
                    key, sign, str(got[key]), str(direct[key]))
            assert got['C_M'] is None and direct['C_M'] is None
            assert got['entire_incoming_interval'] is False
            assert direct['entire_incoming_interval'] is False
            brackets.append(got['N_bracket3'])
        assert len(brackets) == 2


def test_group14_center_pilot_overlaps_stored_enclosure():
    """One full path, both signs, 1024 cells, order 8. Overlap is not the proof."""
    family, identity = require_original_history()
    labels = original_group14_labels()
    record = json.loads(RECORD.read_text())
    assert record['higher_uv_remainder_C_M'] is None
    assert record['coverage']['entire_incoming_interval'] is False
    assert record['coverage']['all_source_families'] is False
    assert len(record['identities']['residuals']) == 74
    assert set(record['identities']['residuals'].values()) == {'0'}
    stored = {row['sign']: row for row in record['cases']}
    assert set(stored) == {1, -1}
    model = unit_cubic_geometry(family, bits=160)
    incoming = arb(family.center)
    deadline = time.process_time() + 90
    started = time.process_time()
    rows = {}
    for sign in (1, -1):
        basis = enclose_cubic_channel_basis(
            model, incoming, sign, ARCHIVE_RHO_UP, cells=1024, max_depth=16,
            order=8, cpu_deadline=deadline)
        assembled = assemble_cubic_channels(
            basis, labels['angular_ball'], labels['mass_ball'])
        flipped = assemble_cubic_channels(
            basis, -labels['angular_ball'], labels['mass_ball'])
        rows[sign] = (basis, assembled, flipped)
    cpu = time.process_time() - started
    summary = {
        'pilot_cpu_seconds': cpu,
        'cells': 1024,
        'order': 8,
        'rho_up': float(ARCHIVE_RHO_UP),
        'profile_identity': identity,
        'overlap_is_factorization_proof': False,
        'entire_incoming_interval': False,
        'C_M': None,
        'signs': {},
    }
    for sign, (basis, assembled, flipped) in rows.items():
        assert basis['accepted_cells'] == 1024
        assert basis['profile_identity'] == identity
        assert basis['C_M'] is None
        assert basis['higher_uv_remainder_bound'] is None
        assert basis['entire_incoming_interval'] is False
        assert basis['z_box_is_full_incoming_interval'] is False
        assert assembled['C_M'] is None
        assert assembled['entire_incoming_interval'] is False
        assert assembled['physical_local_gate'] == 'OPEN'
        case = stored[sign]
        assert case['cells'] == 1024 and case['order'] == 8
        widths = {}
        for key in COEFFICIENTS:
            ball = assembled[key]
            assert ball.is_finite()
            assert _overlaps_stored(ball, case['values'][key]), (sign, key, str(ball))
            assert flipped[key].overlaps(ball)
            with ctx.workprec(256):
                lower = restored_upper(case['values'][key]['lower'])
                upper = restored_upper(case['values'][key]['upper'])
                stored_rad = (upper - lower) / 2
            widths[key] = _display(ball)
            widths[key]['stored_rad'] = float(stored_rad)
            widths[key]['rad_over_stored'] = float(ball.rad() / stored_rad)
            assert ball.rad() < arb('1e-8')
            assert widths[key]['rad_over_stored'] < 8
        summary['signs'][str(sign)] = {
            'accepted_cells': basis['accepted_cells'],
            'cpu_seconds': basis['cpu_seconds'],
            'widths': widths,
        }
    print(json.dumps(summary, indent=2, sort_keys=True))
    assert cpu <= 90
    assert not rows[1][1]['N_bracket3'].overlaps(rows[-1][1]['N_bracket3'])


@pytest.mark.parametrize('missing', [None, False, True])
def test_missing_channel_data_are_not_zero(missing):
    with pytest.raises(ValueError, match='explicit'):
        split_massless_channels(missing,0,0,0,1,1)
    channels=split_massless_channels(0,0,0,0,1,1)
    with pytest.raises(ValueError, match='explicit'):
        integrate_shared_cell((missing,arb(0),arb(0),arb(0),arb(0)),channels,None,2,1,0)
    basis={k:arb(0) for k in ('Qz','Qzz','J2','J4','JM','N2','N4','NM')}
    basis.update(surface_a=arb(1),sign=1,bits=160,massless_A=True,C_M=None,
                 higher_uv_remainder_C_M=None,entire_incoming_interval=False)
    for key in ('J2','J4','JM','Qzz','surface_a'):
        bad=dict(basis);bad[key]=missing
        with pytest.raises(ValueError, match='explicit'):assemble_cubic_channels(bad,1,1)
