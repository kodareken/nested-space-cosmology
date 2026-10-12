import importlib.util
import json
from copy import deepcopy
from pathlib import Path
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
spec = importlib.util.spec_from_file_location('source_update_v5', ROOT/'scripts/derive_nsc_incoming_source_update_v5.py')
script = importlib.util.module_from_spec(spec); spec.loader.exec_module(script)


def inputs():
    return [json.loads((ROOT/p).read_text()) for p in script.INPUTS[:2]]


def test_reject_overlap_wrong_cell_or_missing_thermal_remainder():
    previous, window = inputs()
    assert script.validate_window(previous, window)['32']['joined'] == [24, 'infinity']
    wrong = deepcopy(previous); wrong['coverage']['32']['joined'][0] = 24
    with pytest.raises(ValueError, match='overlap'):
        script.validate_window(wrong, window)
    wrong = deepcopy(window); wrong['energy_interval'] = [16, 24]
    with pytest.raises(ValueError, match='owned'):
        script.validate_window(previous, wrong)
    wrong = deepcopy(window); wrong['positive_thermal_stress_error_upper'][2] = 0
    with pytest.raises(ValueError, match='thermal'):
        script.validate_window(previous, wrong)


def test_authenticated_replay_keeps_other_source_parts_and_open_gates():
    got = json.loads(json.dumps(script.make_record()))
    assert got == json.loads(script.OUTPUT.read_text())
    previous, _ = inputs()
    for key, value in previous['matter']['parts'].items():
        assert got['matter']['parts'][key] == value
    assert got['partial_error_budget']['covered_domain_status'] == 'PASS'
    assert got['partial_error_budget']['full_source_error_bound'] is None
    assert got['baseline']['physical_constraint_status'] == 'OPEN'
    assert not got['scope']['scientific_producers_run']
    assert not got['scope']['metric_evolution']
