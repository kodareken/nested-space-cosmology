import importlib.util
import json
from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest

from recursive_horizons.nsc_incoming_source_update import refined_matter_source
from recursive_horizons.nsc_incoming_joint_constraints import source_action_gradient
from recursive_horizons.nsc_incoming_cauchy_jets import incoming_cauchy_jets

ROOT=Path(__file__).resolve().parents[1]
NAMES=('nsc-incoming-spectral-source','nsc-incoming-subgap-source','nsc-incoming-source-tail',
       'nsc-incoming-middle-order-correction','nsc-incoming-retained-order-correction','nsc-incoming-low-high-source')


def records():return[json.loads((ROOT/f'results/development/{n}.json').read_text()) for n in NAMES]


def test_each_interval_is_composed_once_in_same_action_units():
    inputs=records();result=refined_matter_source(*inputs)
    increment=np.asarray(inputs[3]['refined48_node_vacuum_correction'])+inputs[4]['vacuum_order24_minus16']+np.asarray(inputs[5]['new_vacuum_minus_archived_LOW32'])
    np.testing.assert_array_equal(increment,result['total_stress_increment'])
    expected=np.asarray(inputs[3]['source_action_gradient_correction'])+inputs[4]['action_gradient_correction']+np.asarray(inputs[5]['additive_action_gradient_correction'])
    np.testing.assert_allclose(result['action_gradient_increment'],expected,rtol=0,atol=3e-23)
    assert not result['scope']['local_light_geometry_included']
    assert result['full_source_error_bound'] is None


def test_wrong_interval_cannot_be_applied_as_the_owned_update():
    inputs=records();wrong=deepcopy(inputs)
    wrong[-1]['interval']=[24.,32.]
    with pytest.raises(ValueError,match='group22 low interval'):
        refined_matter_source(*wrong)
    wrong=deepcopy(inputs);wrong[3]['new_order']=16
    with pytest.raises(ValueError,match='order24-minus16'):
        refined_matter_source(*wrong)


def test_full_source_and_constraint_gate_stay_open_after_component_passes():
    spec=importlib.util.spec_from_file_location('source_update_record',ROOT/'scripts/derive_nsc_incoming_source_update.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    current=json.loads(json.dumps(module.make_record(),default=module.encode))
    stored=json.loads((ROOT/'results/development/nsc-incoming-source-update.json').read_text())
    assert current==stored
    assert current['partial_error_budget']['full_source_error_bound'] is None
    assert not current['partial_error_budget']['other_low_subgap_errors_bounded']
    assert current['scope']['local_light_geometry_count']==1
    assert current['baseline']['physical_constraint_status']=='OPEN'
    assert current['baseline']['certified_constraint_residual'] is None
    assert not current['scope']['metric_evolution']
