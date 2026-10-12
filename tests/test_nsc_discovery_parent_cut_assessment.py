"""Assessment refusals and original-column matched-clock identities; no evolution."""
from pathlib import Path
import numpy as np
import pytest
import assess_nsc_discovery_parent_cut_response as assess


def test_column_sum_retains_original_weights_and_proper_clock_correction():
    c=assess.cut;r=c.response;le=c.leading;p=c.parent
    grid=r.install_native_fft(p.galerkin.build_grid(16,quadrature=64,gauge='conformal'))
    x=np.arange(16)*grid.length/16
    phi0=np.column_stack((np.exp(1j*np.pi*x/8),np.exp(3j*np.pi*x/8)))/np.sqrt(32)
    phi1=1j*phi0;source=np.vstack((phi0,phi1))
    pair=r.make_holder(grid,weights=np.array([.2,.1]),reference_columns=source,source_columns=source)
    state=le.encode(pair,le.State(np.ones(grid.ng),np.ones(grid.ng),np.full(grid.ng,.01),np.zeros(grid.ng),phi0,phi1))
    zero=[np.zeros_like(getattr(state,n)) for n in le.FIELDS]
    tangent=r.ParentTangent(*zero,np.array([0.,.1]))
    content,delta=assess.columns(pair,state,tangent,.3)
    total,dt=r.child_regional_content(pair,state,tangent)
    _,slope=r.child_regional_content(pair,state,r.velocity_tangent(pair,state))
    expected=r.centre_clock_variation(dt,slope,.3,r.centre_clock(pair,state))
    assert content.sum()==pytest.approx(total,abs=1e-14)
    assert delta.sum()==pytest.approx(expected,abs=1e-14)
    _,unshifted=assess.columns(pair,state,tangent,0.)
    assert unshifted[0]==0.  # c0 fixed, no field tangent, no clock displacement.
    assert unshifted[1]==pytest.approx(content[1])


def test_authentication_rejects_source_digest_and_input_drift(tmp_path,monkeypatch):
    path=tmp_path/'baseline';path.write_text('original')
    record={'producing_commit':assess.COMMIT,'producers':{'owner':'accepted'},'inputs':{str(path):assess.digest(path)}}
    calls=[]
    def verify(root,path,sha,commit):
        calls.append((path,sha,commit))
        if sha!='accepted':raise ValueError('source digest failed')
    monkeypatch.setattr(assess.cut.provenance,'resolve_pinned_source_bytes',verify)
    assess.authenticate(record);assert calls==[('owner','accepted',assess.COMMIT)]
    record['producers']['owner']='corrupt'
    with pytest.raises(ValueError,match='source digest'):assess.authenticate(record)
    record['producers']['owner']='accepted';path.write_text('changed')
    with pytest.raises(ValueError,match='input hash'):assess.authenticate(record)


def test_effect_error_uses_effect_scale_and_zero_is_undefined():
    row=assess.effect(.01,.009)
    assert row['relative_to_predicted_effect']==pytest.approx(.1)
    assert row['relative_to_measured_effect']==pytest.approx(1/9)
    assert assess.effect(0.,0.)['relative_to_predicted_effect'] is None


def test_creation_only_readonly_writer_and_guards(tmp_path,monkeypatch):
    path=tmp_path/'record.json';sealed=tmp_path/'inputs'
    assess.write_record(path,{'value':3},sealed)
    assert path.stat().st_mode&0o222==0
    with pytest.raises(FileExistsError):assess.write_record(path,{},sealed)
    with pytest.raises(ValueError,match='sealed'):assess.write_record(sealed/'new.json',{},sealed)
    monkeypatch.setattr(assess,'LIMIT',2)
    with pytest.raises(ValueError,match='64 MiB'):assess.write_record(tmp_path/'big.json',{'x':4},sealed)
    assert not (tmp_path/'big.json').exists()
