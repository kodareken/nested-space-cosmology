import importlib.util
import json
from copy import deepcopy
from pathlib import Path
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
spec = importlib.util.spec_from_file_location('source_update_v4', ROOT/'scripts/derive_nsc_incoming_source_update_v4.py')
script = importlib.util.module_from_spec(spec); spec.loader.exec_module(script)


def inputs():
    return [json.loads((ROOT/p).read_text()) for p in script.INPUTS[:2]]


def test_tail_union_rejects_gap_duplicate_and_changed_source():
    previous, tail = inputs()
    assert script.join_coverage(previous, tail)[13]['joined'] == [16., 'infinity']
    bad = deepcopy(tail); bad['groups'][0]['lower_energy'] -= 1
    with pytest.raises(ValueError, match='meet'):
        script.join_coverage(previous, bad)
    bad = deepcopy(tail); bad['groups'][0] = bad['groups'][1]
    with pytest.raises(ValueError, match='distinct'):
        script.join_coverage(previous, bad)
    bad = deepcopy(tail); bad['aggregate']['archived_total_order13_source'][0] = 0
    with pytest.raises(ValueError, match='exact source'):
        script.join_coverage(previous, bad)


def test_authentic_replay_preserves_sources_and_unresolved_gates():
    got = json.loads(json.dumps(script.make_record()))
    assert got == json.loads(script.OUTPUT.read_text())
    previous, _ = inputs()
    assert got['matter']['stress_approximant'] == previous['matter']['stress_approximant']
    assert got['matter']['parts'] == previous['matter']['parts']
    assert 'tail numerical integration, distinct from physical mode error' not in got['matter']['unresolved']
    assert got['baseline'] == previous['baseline']
    assert got['partial_error_budget']['covered_domain_status'] == 'PASS'
    assert got['partial_error_budget']['full_source_error_bound'] is None
    assert not got['scope']['thermal_declared_zero']
    assert not got['scope']['metric_evolution']
