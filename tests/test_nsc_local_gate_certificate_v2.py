"""Certificate-v2 decisions cannot promote missing bounds or optimizer stalls."""
from fractions import Fraction

import pytest

from recursive_horizons.nsc_local_gate_certificate_v2 import (
    ERROR_COMPONENTS, REQUIRED_MUTATIONS, TOLERANCE, build_certificate,
    existence_decision, nonexistence_decision)


def components(value):
    return {name: [value, value] for name in ERROR_COMPONENTS}


def mutations(value=True):
    return {name: value for name in REQUIRED_MUTATIONS}


def test_existence_requires_every_component_and_negative_control():
    result = existence_decision([Fraction(1, 10**12)] * 2,
                                components(Fraction(1, 10**13)), mutations())
    assert result['verdict'] == 'EXISTENCE'
    missing = components(Fraction(1, 10**13))
    missing['upstream'] = None
    assert existence_decision([0, 0], missing, mutations())['verdict'] == 'OPEN'
    failed = mutations()
    failed['drop-coherence'] = False
    assert existence_decision([0, 0], components(0), failed)['verdict'] == 'OPEN'


def test_residual_plus_error_is_componentwise_and_exact():
    each = TOLERANCE / (2 * len(ERROR_COMPONENTS))
    result = existence_decision([TOLERANCE / 2] * 2, components(each), mutations())
    assert result['verdict'] == 'EXISTENCE'
    too_large = existence_decision([TOLERANCE] * 2, components(each), mutations())
    assert too_large['verdict'] == 'OPEN'


def test_nonexistence_requires_strict_class_wide_relation_not_stall():
    obstruction = {
        'necessary_relation': True,
        'entire_declared_class': True,
        'separation_lower': Fraction(5, 10**11),
        'error_upper': Fraction(1, 10**11),
        'optimizer_failure_used': False,
    }
    assert nonexistence_decision(obstruction)['verdict'] == 'NON_EXISTENCE'
    obstruction['optimizer_failure_used'] = True
    assert nonexistence_decision(obstruction)['verdict'] == 'OPEN'
    obstruction['optimizer_failure_used'] = False
    obstruction['error_upper'] = obstruction['separation_lower']
    assert nonexistence_decision(obstruction)['verdict'] == 'OPEN'


def test_certificate_serializes_exact_fractions_and_open_state():
    certificate = build_certificate(
        mode='existence', profile_identity='p', source_identity='s', interval=[0, 1],
        state_law='C=U C_up U^dagger', history_class='w,U', numerical_settings={},
        coverage={}, provenance={}, residual_upper=[0, 0],
        components={**components(0), 'arithmetic': None}, mutations=mutations())
    assert certificate['schema'].endswith('-v2')
    assert certificate['verdict'] == 'OPEN'
    assert certificate['open_is_not_done']
    assert certificate['failed_optimization_is_nonexistence'] is False


def test_wrong_error_schema_rejected():
    with pytest.raises(ValueError, match='nine-component'):
        existence_decision([0, 0], {'field_space_time': [0, 0]}, mutations())
