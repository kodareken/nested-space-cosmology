"""Finite primitive, source-law, action-unit and endpoint ownership controls."""
from copy import deepcopy
from fractions import Fraction
import json
from pathlib import Path

import mpmath as mp
import numpy as np
import pytest

from recursive_horizons.nsc_incoming_middle_bound import (
    authenticated_middle_inputs, finite_power_primitive, validate_middle_interval,
    vacuum_middle_bound, thermal_middle_difference_bound, constraint_action_bound,
    STATIONARITY_TOLERANCE,
)
from recursive_horizons.nsc_incoming_vacuum_tail_bound import _precision, _lo, _hi
from recursive_horizons.nsc_incoming_cauchy_jets import incoming_cauchy_jets
from recursive_horizons.nsc_incoming_joint_constraints import source_action_gradient


ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope='module')
def owned(): return authenticated_middle_inputs(ROOT)


def test_finite_primitives_enclose_exact_rational_integrals():
    with _precision(40):
        for left, right in ((16, 160), (40, 160), (40, 320)):
            for power in (16, 24, 32):
                for moment in (0, 1):
                    p = power-moment-1
                    exact = (Fraction(1, left**p)-Fraction(1, right**p))/p
                    value = mp.mpf(exact.numerator)/exact.denominator
                    bound = finite_power_primitive(power, moment, left, right)
                    assert 0 < _lo(bound) <= value <= _hi(bound)
    with pytest.raises(ValueError, match='finite interval'):
        finite_power_primitive(16, 0, 160., 40.)


def test_authenticated_selected_interval_and_group_rejections(owned):
    assert len(owned['middle_intervals']) == 32
    assert owned['endpoint_measure_residual'] <= 3e-11
    for row in owned['middle_intervals']:
        group = row['group']
        assert row['left'] == (16. if group in (13, 14) else 40.)
        assert row['right'] == (320. if group in (10, 11, 12, 31, 32) else 160.)
    channel = owned['channels'][22]
    metadata = {'channel': channel, 'angular_sign': 1, 'left': 40., 'right': 160.}
    assert validate_middle_interval(channel, 1, metadata, 160.) == (40., 160.)
    for field, changed in (('left', 16.), ('right', 320.), ('angular_sign', -1)):
        wrong = deepcopy(metadata); wrong[field] = changed
        with pytest.raises(ValueError, match='endpoints'):
            validate_middle_interval(channel, 1, wrong, 160.)
    with pytest.raises(ValueError, match='endpoint'):
        vacuum_middle_bound(channel, owned['radial_groups'][21], 40., 320.)


def test_action_bound_has_existing_source_action_normalization():
    stress = [1e-8, 2e-8, 3e-8, 4e-8]
    expected = abs(source_action_gradient(incoming_cauchy_jets(), stress)[:2])
    upper = np.array(constraint_action_bound(stress))
    assert np.all(upper >= expected)
    assert np.max(abs(upper-expected)/expected) < 2e-15
    assert constraint_action_bound([1., 0., 0., 0.])[1] == 0.
    assert STATIONARITY_TOLERANCE == 3e-11


def test_thermal_trace_triangle_bound_matches_independent_exponential_integral(owned):
    channel = owned['channels'][14]; config = owned['config']
    upper = thermal_middle_difference_bound(channel, config, 16., 160.)
    # Independent quadrature of the already-owned scalar bound; no mode or
    # radial recurrence. This checks the two-departure factor and units.
    with mp.workdps(90):
        a, r = mp.sqrt(3*mp.pi/2-4), mp.sqrt(2)
        m, ell = mp.mpf(channel['compact_mass']), mp.mpf(channel['angular_eigenvalue'])
        kappa, omega = mp.mpf(config['surface_gravity']), mp.mpf(config['omega'])
        factor = channel['copy_count']*channel['degeneracy']/(4*mp.pi**2*r*r*a)
        b = lambda E: 2*mp.exp(-mp.pi*E/kappa)+mp.exp(-2*mp.pi*E/(omega*kappa))
        current = 2*factor*mp.quad(lambda E: E/a*b(E), [16, 17, 20, 40, 160])
        density = 2*factor*mp.quad(lambda E: (E/a+mp.sqrt(m*m+(ell/r)**2))*b(E), [16, 17, 20, 40, 160])
        assert current <= upper[2] and density <= upper[0]
        assert abs(mp.mpf(upper[2])/current-1) < mp.mpf('2e-14')
        assert abs(mp.mpf(upper[0])/density-1) < mp.mpf('2e-14')


def test_record_is_valid_enclosure_without_nonexistence_claim(owned):
    record = json.loads((ROOT/'results/development/nsc-incoming-middle-bound.json').read_text())
    assert record['valid_enclosure']
    assert not record['physical_NON_EXISTENCE_claimed']
    assert record['group_count'] == 32
    assert record['aggregate']['stationarity_tolerance'] == 3e-11
    for group in record['groups']:
        index = group['group']; interval = group['middle_interval']
        actual = vacuum_middle_bound(owned['channels'][index], owned['radial_groups'][index-1],
                                     interval['left'], interval['right'])
        assert actual == group['vacuum_stress_error_upper']
    assert not record['scope']['quadrature_error_covered']
    assert not record['scope']['low_subgap_modal_accuracy_covered']
    assert not record['scope']['radial_coefficient_preparation_rerun']
