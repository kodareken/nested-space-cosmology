import importlib.util
import json
from copy import deepcopy
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('source_update_v3', ROOT/'scripts/derive_nsc_incoming_source_update_v3.py')
script = importlib.util.module_from_spec(spec); spec.loader.exec_module(script)


def inputs():
    records = {n: json.loads((ROOT/f'results/development/nsc-incoming-{n}.json').read_text()) for n in script.NAMES}
    return records['high-band-quadrature-batch']['groups'], script.physical_pieces(records)


def test_direct_unions_reject_missing_coverage_old_order_and_double_count():
    direct, physical = inputs()
    script.validate_coverage(direct, physical)
    wrong = deepcopy(physical); wrong[12].append(wrong[12][0])
    with pytest.raises(ValueError, match='overlap'):
        script.validate_coverage(direct, wrong)
    wrong = deepcopy(direct); wrong['13']['interval'] = [32, 160]
    with pytest.raises(ValueError, match='interval'):
        script.validate_coverage(wrong, physical)
    wrong = deepcopy(direct); wrong['1']['source_order'] = 16
    with pytest.raises(ValueError, match='order24'):
        script.validate_coverage(wrong, physical)
    wrong = deepcopy(physical); del wrong[32]
    with pytest.raises(ValueError, match='exactly32'):
        script.validate_coverage(direct, wrong)


def test_authenticated_composition_replays_and_preserves_open_gates():
    got = json.loads(json.dumps(script.make_record()))
    assert got == json.loads(script.OUTPUT.read_text())
    budget = got['partial_error_budget']
    assert budget['finite_high_status'] == 'PASS'
    assert budget['full_source_error_bound'] is None
    assert budget['tail_numerical_integration_error_bound'] is None
    assert got['baseline']['physical_constraint_status'] == 'OPEN'
    assert not got['scope']['prior_high_corrections_added']
    assert not got['scope']['thermal_declared_zero']
    assert not got['scope']['metric_evolution']
