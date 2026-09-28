import importlib.util
import json
from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('source_update_v2',ROOT/'scripts/derive_nsc_incoming_source_update_v2.py')
script=importlib.util.module_from_spec(spec);spec.loader.exec_module(script)


def test_adjacent_windows_count_once_and_reject_modified_interval():
    a=json.loads((ROOT/'results/development/nsc-incoming-window-reuse.json').read_text())
    b=json.loads((ROOT/'results/development/nsc-incoming-group32-source.json').read_text())
    windows={**a['windows'],**b['windows']};coverage=script.validate_windows(windows)
    assert sorted(coverage[32])==[[32.,40.],[40.,320.]]
    wrong=deepcopy(windows);wrong['group32_middle']['interval']=[32.,320.]
    with pytest.raises(ValueError,match='window'):
        script.validate_windows(wrong)
    wrong=deepcopy(windows);wrong['union']=wrong['group32_middle']
    with pytest.raises(ValueError,match='four explicit'):
        script.validate_windows(wrong)


def test_source_update_replays_and_does_not_promote_component_pass():
    got=json.loads(json.dumps(script.make_record()))
    want=json.loads((ROOT/'results/development/nsc-incoming-source-update-v2.json').read_text())
    assert got==want
    assert got['partial_error_budget']['full_source_error_bound'] is None
    assert not got['scope']['metric_evolution']
    assert got['baseline']['physical_constraint_status']=='OPEN'
    previous=json.loads((ROOT/'results/development/nsc-incoming-source-update.json').read_text())
    np.testing.assert_allclose(np.asarray(got['baseline']['action_gradient_approximant'])-previous['baseline']['action_gradient_approximant'],got['additional_action_gradient'][:2],rtol=0,atol=3e-13)
