"""Mutation controls on a captured-cell record, without expensive enclosure."""
from copy import deepcopy
import importlib
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
P = importlib.import_module('validate_nsc_ks_field_ramp_pilot')


def record():
    pair = [{'mantissa': '1', 'exponent': -50}, {'mantissa': '1', 'exponent': -40}]
    result = {'schema': 'NSC-KS-FIELD-RAMP-CELL-PILOT-v1', 'all_history_cells': False,
              'full_family_source_coverage': False, 'physical_rho1_source_error': None,
              'physical_local_gate': 'OPEN',
              'cell': {'cell_index': 150, 'total': pair, 'integral': deepcopy(pair)}}
    result['scientific_digest'] = P.proof_digest(result)
    return result


def test_altered_bound_fails_integrity():
    value = record()
    P.validate_pilot_shape(value)
    value['cell']['total'][0]['exponent'] -= 20
    with pytest.raises(ValueError, match='digest'):
        P.validate_pilot_shape(value)


def test_reintroducing_cell_width_is_rejected_even_with_updated_digest():
    value = record()
    value['cell']['integral'][0]['exponent'] -= 13
    value['scientific_digest'] = P.proof_digest(value)
    with pytest.raises(ValueError, match='cell width'):
        P.validate_pilot_shape(value)


@pytest.mark.parametrize('name,value', [('all_history_cells', True),
    ('full_family_source_coverage', True), ('physical_rho1_source_error', 0),
    ('physical_local_gate', 'EXISTENCE')])
def test_single_cell_cannot_claim_closure(name, value):
    changed = record()
    changed[name] = value
    changed['scientific_digest'] = P.proof_digest(changed)
    with pytest.raises(ValueError, match='physical gate'):
        P.validate_pilot_shape(changed)
