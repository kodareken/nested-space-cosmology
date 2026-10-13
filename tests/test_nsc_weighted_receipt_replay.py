"""Finite-recipe replay bounds preserve exact-zero channels and scope."""
from copy import deepcopy
import json
from pathlib import Path
import numpy as np
import pytest
from recursive_horizons.nsc_weighted_receipt_replay import (
    weighted_bound, require_values, require_record)


def test_weighted_bound_tracks_uncancelled_operands_and_rejects_larger_drift():
    weights = np.array([.5,.5])
    kernels = np.array([[1e12,1.,0.,1e-100],[-1e12,1.,0.,-1e-100]])
    bound = weighted_bound(weights,kernels,2.)
    assert bound[0] > 1e-4
    assert bound[3] > 0 and bound[3] < 1e-112
    assert bound[2] == 0
    original = 2*np.einsum('n,nv->v',weights,kernels)
    changed = original.copy();changed[0] += 8*bound[0]
    with pytest.raises(ValueError,match='arithmetic envelope'):
        require_values(original,changed,bound)
    changed = original.copy();changed[2] = 1e-300
    with pytest.raises(ValueError,match='arithmetic envelope'):
        require_values(original,changed,bound)


def test_reduction_guard_never_substitutes_status_or_provenance():
    saved={'source':[1.,0.,0.,0.],'status':'OPEN','sha256':'a'*64,'tolerance':3e-11}
    bounds={('source',):np.array([1e-14,0.,0.,0.])}
    current=deepcopy(saved);current['source'][0] += 1e-14
    assert require_record(saved,current,bounds) is saved
    current['status']='PASS'
    with pytest.raises(ValueError,match='outside declared'):
        require_record(saved,current,bounds)


def test_receipt_guards_reject_numeric_text_complex_and_bool_aliases():
    with pytest.raises(ValueError,match='numeric type'):
        require_record({'value':1.0},{'value':'1.0'},{('value',):1e-15})
    with pytest.raises(ValueError,match='strict finite real'):
        weighted_bound([1.],np.ones((1,4),complex),1.)
    with pytest.raises(ValueError,match='outside declared'):
        require_record({'value':1.,'scope':True},{'value':1.,'scope':1},{('value',):1e-15})


def test_coherent_matter_receipt_rejects_an_out_of_envelope_nodal_change():
    from recursive_horizons.nsc_frozen_sources import frozen_source_context
    from recursive_horizons.nsc_ks_source_envelope import comparison_control as C
    from recursive_horizons.nsc_ks_matter_receipt_replay import matter_bounds, require_envelope_record
    with frozen_source_context('329ef809e466a86ed4a290ee32c19485f85add9d'):
        saved=C.check()
        with np.load(C.ROOT/C.ARTIFACT,allow_pickle=False) as payload:
            arrays={key:payload[key] for key in payload.files}
        meta=json.loads(arrays['metadata_json'].tobytes())
        current=C.record(arrays,meta,saved['payload'])
        bound,_=matter_bounds(arrays,'n128',meta['parameters'])
        row,column=np.unravel_index(np.argmax(bound),bound.shape)
        corrupted=deepcopy(current)
        corrupted['measurements']['per_case']['n128']['matter'][row][column] += float(8*bound[row,column])
        with pytest.raises(ValueError,match='arithmetic envelope'):
            require_envelope_record(saved,corrupted,arrays,meta)
        corrupted=deepcopy(current);corrupted['scope']['physical_local_gate']='PASS'
        with pytest.raises(ValueError,match='outside declared'):
            require_envelope_record(saved,corrupted,arrays,meta)
        bindings=C.successor_source_hashes()
        assert 'src/recursive_horizons/nsc_ks_source_envelope.py' in bindings


def test_germ_receipt_rejects_a_coefficient_change_beyond_its_reference_sum():
    from recursive_horizons.nsc_weighted_receipt_replay import germ_receipt_bounds
    root=Path(__file__).resolve().parents[1]
    saved=json.loads((root/'results/development/nsc-incoming-constraint-germ.json').read_text())
    payload=json.loads((root/saved['payload']['path']).read_text())
    joint=json.loads((root/'results/development/nsc-incoming-joint-constraints.json').read_text())
    baseline=json.loads((root/'results/development/nsc-incoming-local-constraints.json').read_text())['baseline']
    bounds=germ_receipt_bounds(saved,payload,joint,baseline)
    changed=deepcopy(saved);changed['numerical_coefficients']['c_u'] += float(8*bounds['numerical_coefficients','c_u'])
    with pytest.raises(ValueError,match='arithmetic envelope'):
        require_record(saved,changed,bounds)
    current=deepcopy(saved);current['sha256']='b'*64
    with pytest.raises(ValueError,match='outside declared'):
        require_record(saved,current,bounds)
