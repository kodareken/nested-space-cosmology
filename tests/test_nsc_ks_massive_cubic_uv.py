"""Exact massive cubic identities and the unchanged directed integrator."""
import numpy as np
import pytest
from flint import arb, ctx

from recursive_horizons.nsc_ks_cubic_uv_current import _require_massless_angular
from recursive_horizons.nsc_ks_cubic_uv_enclosure import (
    CubicGeometry, enclose_characteristic,
)
from recursive_horizons.nsc_ks_current_uv_transport import ARCHIVE_RHO_UP
from recursive_horizons.nsc_ks_massive_cubic_uv import (
    MASSLESS_J3, MassiveCubicGeometry, enclose_massive_characteristic,
    massive_cubic_identities, original_group14_labels, require_original_history,
)
from recursive_horizons.nsc_local_incoming_family import LocalIncomingFamily


REQUIRED = (
    'minor_A2_from_L0', 'minor_A3_from_L0', 'L_q', 'L_h_parent',
    'L_h_massless_limit', 'L_h_nonzero_unit_jet', 'L_h3_from_V',
    'continuity_n4', 'L_h3z_minus_s_n4', 'Re_lC_independent_of_h',
    'L_qz', 'L_qzz', 'L_J3', 'L_J3_massless_reduction',
    'mass_correction_not_zero', 'J3_pauli', 'J3_h_coefficient',
    'J3_hz_coefficient', 'raw_h_coefficient_is_q_z', 'raw_hz_coefficient_is_q',
    'N_minus_s_J3_over_a', 'Sigma_difference', 'Sigma_common_k_cancels',
    'J_factor', 'N_factor', 'Qz_factor', 'Qzz_factor', 'even_in_ell',
)


def test_generic_residuals_keep_h_until_it_cancels_and_factor_three_powers():
    proven = massive_cubic_identities()
    residuals = proven['residuals']
    assert residuals and set(residuals.values()) == {'0'}
    for sign in (1, -1):
        for name in REQUIRED:
            assert residuals['s%+d/%s' % (sign, name)] == '0'
    assert proven['massless_Ts_J3'] == MASSLESS_J3
    assert proven['Ts_J3_massive'].endswith('+ m^2 (s q_zz - 4 b b_z)')
    assert ' - s a m^2 delta(q_z)' in proven['Sigma_delta_N_bracket']
    assert proven['massless_h_equation_imposed'] is False
    assert proven['h_cancelled_from_J3_and_Re_lC'] is True
    assert proven['common_k_cancels_on_Sigma'] is True
    assert proven['L_h_parent'] == '-m ell/4 (s a^2 r_rho + r_z)/r^2'
    factor = proven['factorization']
    assert factor['J3'] == 'ell^2 J2 + ell^4 J4 + m^2 ell^2 JM'
    assert factor['N_bracket'] == 'ell^2 N2 + ell^4 N4 + m^2 ell^2 NM'
    assert factor['JM_prime'] == '-4 b0 b0_z + s Qzz'
    assert factor['J4_prime'] == '-16 b0^3 b0_z/a^2 + (4 s b0^2/a^2) Qzz'
    assert proven['action_N'] == '-mu/(2 pi) N_bracket'
    assert proven['action_beta'] == '+mu/(2 pi) J3'
    assert proven['group14_label_action_signs'] == {
        's+1/N': 2.0, 's+1/beta': -2.0, 's-1/N': 2.0, 's-1/beta': 2.0}
    for name in ('higher_uv_remainder_C_M', 'uniform_C4_on_I', 'complete_UV_tail'):
        assert proven[name] is None
    assert proven['physical_local_gate'] == 'OPEN'
    assert proven['one_point_or_formal_coverage_is_full_source_tail'] is False
    assert proven['factorization_replaces_a_finished_family_integral'] is False


def test_massless_guard_still_rejects_the_original_massive_channel():
    with pytest.raises(ValueError, match='massive'):
        _require_massless_angular(np.pi / 2, np.sqrt(5))


def test_zero_mass_forcing_matches_the_historical_owner():
    family, _ = require_original_history()
    with ctx.workprec(160):
        angular = arb(5).sqrt()
        owner = CubicGeometry(family, angular, bits=160)
        wrapped = MassiveCubicGeometry(family, angular, arb(0), bits=160)
        rho, z = arb('1.018'), arb(family.center)
        for sign in (1, -1):
            base = owner.forcing(rho, z, sign)
            same = wrapped.forcing(rho, z, sign)
            for name in vars(base):
                assert getattr(same, name).overlaps(getattr(base, name))
            jets = wrapped.forcing_series(rho, z, sign, 3)
            owner_jets = owner.forcing_series(rho, z, sign, 3)
            for name in owner_jets:
                for left, right in zip(jets[name].coeffs(), owner_jets[name].coeffs()):
                    assert (left - right).contains(0)


def test_massive_additions_signs_and_angular_evenness():
    family, _ = require_original_history()
    with ctx.workprec(160):
        angular = arb(5).sqrt()
        mass = arb.pi() / 2
        owner = CubicGeometry(family, angular, bits=160)
        positive = MassiveCubicGeometry(family, angular, mass, bits=160)
        negative = MassiveCubicGeometry(family, -angular, mass, bits=160)
        rho, z = arb('1.018'), arb(family.center)
        axial = owner.jets(rho, z)['a']
        for sign in (1, -1):
            base = owner.forcing(rho, z, sign)
            moved = positive.forcing(rho, z, sign)
            assert (moved.qz - base.qz).contains(0)
            assert (moved.qzz - base.qzz).contains(0)
            assert (moved.coupling - base.coupling - sign * mass**2).contains(0)
            assert (moved.current - base.current - mass**2 * axial**2 * base.qz).contains(0)
            other = negative.forcing(rho, z, sign)
            for name in vars(moved):
                assert getattr(moved, name).overlaps(getattr(other, name))
        with pytest.raises(ValueError):
            MassiveCubicGeometry(family, angular, arb(-1), bits=160)
        with pytest.raises(ValueError, match='nonnegative'):
            MassiveCubicGeometry(family, angular, arb(0, '.1'), bits=160)
        with pytest.raises(TypeError):
            enclose_massive_characteristic(owner, z, 1, ARCHIVE_RHO_UP, cells=1, order=0)


def test_surface_correction_and_null_tail_flags_on_a_short_path():
    family, _ = require_original_history()
    with ctx.workprec(160):
        model = MassiveCubicGeometry(
            family, arb(5).sqrt(), arb.pi() / 2, bits=160)
        z = arb(family.center)
        raw = enclose_characteristic(model, z, 1, ARCHIVE_RHO_UP, cells=8, order=3)
        done = enclose_massive_characteristic(
            model, z, 1, ARCHIVE_RHO_UP, cells=8, order=3)
        assert (done['N_bracket3'] - raw['N_bracket3']).overlaps(
            done['surface_mass_correction'])
        assert done['surface_mass_correction'].overlaps(
            -done['mass']**2 * model.jets(arb(1), z)['a'] * done['q_z'])
        assert done['C_M'] is None
        assert done['complete_UV_tail'] is None
        assert done['higher_uv_remainder_bound'] is None
        assert done['entire_incoming_interval'] is False
        assert done['physical_local_gate'] == 'OPEN'
        zero = enclose_massive_characteristic(
            MassiveCubicGeometry(LocalIncomingFamily(np.zeros((2, 8))),
                                 np.sqrt(5), np.pi / 2),
            arb('1.2'), -1, ARCHIVE_RHO_UP, cells=4, order=2)
        for key in ('q_z', 'q_zz', 'J3', 'N_bracket3'):
            assert zero[key].contains(0)
        assert zero['complete_UV_tail'] is None
        assert zero['physical_local_gate'] == 'OPEN'


def test_group14_metadata_is_mu_12_cutoff_160_positive_pi_over_2():
    labels = original_group14_labels()
    assert labels['group'] == 14
    assert labels['positive_mass'] is True
    assert labels['multiplicity_per_signed_family'] == 12
    assert labels['degeneracy'] == 12
    assert labels['cutoff'] == 160
    assert labels['quadrature_measure_by_sign'] == {1: 160.0, -1: 160.0}
    with ctx.workprec(160):
        assert labels['mass_ball'].contains(arb.pi() / 2)
        assert labels['angular_ball'].contains(arb(5).sqrt())
    _, identity = require_original_history()
    assert identity.startswith('0b0e4ced')
