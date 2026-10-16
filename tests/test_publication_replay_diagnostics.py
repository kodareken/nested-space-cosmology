"""Failure descriptions preserve rejection and expose only selected operands."""
import json
from types import SimpleNamespace

import numpy as np
import pytest

from check_publication import replay_failure_message
from record_nsc_atlas_nonlinear_rectangle import _same_recomputed


def test_failed_equal_reports_rejected_leaf_using_original_comparator():
    namespace = dict(previous=SimpleNamespace(_same_recomputed=_same_recomputed), plain=lambda x:x)
    exec(compile('def _equal(actual, expected, label):\n raise ValueError(label+" mismatch")\n',
                 'record_nsc_atlas_coupled_evolution.py', 'exec'), namespace)
    with pytest.raises(ValueError) as caught:
        namespace['_equal']({'small':1.+1e-14,'nested':[1.,2.]},
                            {'small':1.,'nested':[1.,2.01]},'diagnostics')
    original = str(caught.value)
    message = replay_failure_message(caught.value)
    assert str(caught.value) == original
    rows = json.loads(message.split(': ',1)[1])
    assert len(rows)==1 and rows[0]['path']=='diagnostics.nested[1]'
    assert float(rows[0]['absolute_difference']) == pytest.approx(.01)


def test_forcing_failure_separates_curve_and_internal_grid():
    namespace = dict(WALL_COORDINATE_REPLAY_ATOL=2e-13,np=np)
    source = '''def from_records():
 t=.5
 wall=1.
 recorded_wall=1.+1e-13
 row={"U":np.array([0.,.5+1e-12,recorded_wall])}
 expected=np.array([0.,recorded_wall/2,recorded_wall])
 raise ValueError("recorded forcing coordinates differ from the original wall")
'''
    exec(compile(source,'nsc_atlas_forcing_cache.py','exec'),namespace)
    with pytest.raises(ValueError) as caught: namespace['from_records']()
    row = json.loads(replay_failure_message(caught.value).split(': ',1)[1])[0]
    assert float(row['wall_difference']) < float(row['wall_allowance'])
    assert row['internal_grid_index']==1
    assert float(row['internal_grid_difference']) > float(row['internal_grid_allowance'])


def test_unrelated_frame_locals_are_not_reported():
    private_local = 'unrelated content'
    try: raise ValueError('ordinary failure')
    except ValueError as error: assert replay_failure_message(error)=='ordinary failure'
