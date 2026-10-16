"""Independent saved-array controls; no source, cache, metric or evolution replay."""
import copy
import json
import time

import numpy as np
import pytest

import diagnose_nsc_atlas_constraint_transport as tool


@pytest.fixture(scope="module")
def budget():
    shared = tool.Budget(maximum_evaluations=300000,cpu_seconds=30.)
    yield shared
    shared.check()
    print("constraint transport controls",shared.evaluations,time.process_time()-shared.started)


def example(budget,nt=9):
    x,eta = np.linspace(0,1,nt),np.linspace(0,1,5)
    W = 1+4*x
    src = 2+x*x
    Es = 2-.5*x
    r = src[:,None]+.2*W[:,None]*eta[None,1:-1]
    extra = .03*x[1:-1,None]*eta[None,1:-1]
    Ex = (src[1:-1,None]*Es[1:-1,None]+extra)/r[1:-1]
    # An arbitrary incoming wall-origin value must remain visible.
    EU = np.broadcast_to(eta[None,1:-1]**2,Ex.shape).copy()
    return tool.decompose_grid(x,eta,W,r,Ex,EU,src,Es,fronts=[.31,.73],W0=1.,budget=budget)


def test_nonzero_inherited_transport_and_independent_residual(budget):
    row = example(budget)
    x,eta = row["interior_x"],row["interior_eta"]
    expected = .03*x[:,None]*eta[None,:]
    np.testing.assert_allclose(row["arrays"]["transport"],expected,atol=9e-16)
    np.testing.assert_allclose(row["arrays"]["residual"],expected/row["radius"],atol=4e-16)
    inherited = (2+x*x)[:,None]*(2-.5*x)[:,None]/row["radius"]
    np.testing.assert_allclose(row["arrays"]["inherited"],inherited,atol=0)
    assert row["summaries"]["inherited"]["all_interior"]["minimum"]>0
    assert row["incoming_constraint"]["wall_expected"]=="UNKNOWN"
    assert row["summaries"]["E_Ubar"]["wall_origin"]["maximum_absolute"]>0


def test_true_inverse_radius_transport_zero_remainder_with_nonzero_source(budget):
    x,eta = np.linspace(0,1,5),np.linspace(0,1,5)
    W=1+4*x;src=2+x;Es=np.full(5,-.4)
    r=src[:,None]+W[:,None]*eta[None,1:-1]
    Ex=src[1:-1,None]*Es[1:-1,None]/r[1:-1]
    row=tool.decompose_grid(x,eta,W,r,Ex,np.zeros_like(Ex),src,Es,
        fronts=[.31,.73],W0=1.,budget=budget)
    np.testing.assert_array_equal(row["arrays"]["residual"],np.zeros_like(Ex))
    np.testing.assert_allclose(row["arrays"]["transport"],0,atol=3e-16)
    assert np.max(abs(row["arrays"]["E_x"]))>.2


def test_source_jet_expression_has_correct_scalar_factor_and_curvature(budget):
    x=np.linspace(.1,.7,7);r=2+x*x;rx=2*x;k=.3;alpha=.2
    source=dict(r=r,r_x=rx,r_xx=np.full(7,2.),log_C_x=np.full(7,.25),
        c=k*r*r,c_x=2*k*r*rx)
    expected=2-.25*rx+alpha*r*k*k*rx*rx/(2*(1+alpha*k*k*r*r/6))
    np.testing.assert_allclose(tool.source_constraint(source,alpha,budget),expected,rtol=1e-15)
    np.testing.assert_array_equal(tool.source_constraint(source,0.,budget),2-.25*rx)
    source["r_xx"]+=.125
    np.testing.assert_allclose(tool.source_constraint(source,alpha,budget)-expected,.125,atol=3e-16)


def test_three_node_origin_partition_and_separator_equality(budget):
    row=example(budget,5)
    # First interior time has W=2: upper endpoint of first stencil=1,
    # lower endpoint of last stencil=1. Only the former is bottom-origin.
    masks=row["masks"]
    np.testing.assert_array_equal(masks["bottom_origin"][0],[True,False,False])
    np.testing.assert_array_equal(masks["wall_origin"][0],[False,False,False])
    np.testing.assert_array_equal(masks["separator"][0],[False,True,True])
    np.testing.assert_array_equal(masks["wall_origin"][1],[False,False,True])
    np.testing.assert_array_equal(masks["bottom_origin"].astype(int)+masks["wall_origin"]+masks["separator"],1)


def matched_pair(budget):
    fine=example(budget)
    x=fine["x"][::2];eta=fine["eta"];W=fine["W"][::2].copy()
    source=fine["source_r"][::2];Es=fine["source_E"][::2]
    radius=source[:,None]+.2*W[:,None]*eta[None,1:-1]
    delta=.01*x[1:-1,None]*eta[None,1:-1]
    Ex=fine["arrays"]["E_x"][1::2]+delta
    EU=fine["arrays"]["E_Ubar"][1::2].copy()
    coarse=tool.decompose_grid(x,eta,W,radius,Ex,EU,source,Es,fronts=[.31,.73],W0=1.,budget=budget)
    return fine,coarse,delta


def test_coarsening_comparison_holds_radius_source_and_uses_union_masks(budget):
    fine,coarse,delta=matched_pair(budget)
    row=tool.compare_sampling(fine,coarse,budget)
    np.testing.assert_array_equal(row["fine_interior_indices"],[1,3,5])
    np.testing.assert_allclose(row["difference_coarse_minus_fine"]["E_x"],delta,atol=2e-16)
    np.testing.assert_allclose(row["difference_coarse_minus_fine"]["transport"],coarse["radius"]*delta,atol=1e-15)
    np.testing.assert_array_equal(row["difference_coarse_minus_fine"]["inherited"],0)
    assert not fine["masks"]["front_near"][3,0]
    assert row["union_masks"]["front_near"][1,0]
    assert row["E_Ubar_same_spatial_stencil_exact"]
    location=row["summaries"]["E_x"]["all_interior"]["location"]
    assert location["x"]==.75 and location["eta"]==.75


@pytest.mark.parametrize("key",("radius","W","interior_Ubar","E_Ubar"))
def test_matched_stencil_inputs_cannot_change_silently(budget,key):
    fine,coarse,_=matched_pair(budget)
    if key=="E_Ubar": coarse["arrays"][key][0,0]+=.01
    elif key=="W": coarse[key][1]+=.01
    else: coarse[key][0,0]+=.01
    with pytest.raises(ValueError):tool.compare_sampling(fine,coarse,budget)


def test_grid_shape_radius_and_original_width_guards(budget):
    x=np.linspace(0,1,5);eta=np.linspace(0,1,5);r=np.full((5,3),2.)
    common=dict(fronts=[.31,.73],W0=1.,budget=budget)
    with pytest.raises(ValueError,match="expanding"):
        tool.decompose_grid(x,eta,np.ones(5),r,np.zeros((3,3)),np.zeros((3,3)),np.ones(5),np.zeros(5),**common)
    with pytest.raises(ValueError,match="expanding"):
        tool.decompose_grid(x,eta,2+x,r,np.zeros((3,3)),np.zeros((3,3)),np.ones(5),np.zeros(5),**common)
    with pytest.raises(ValueError,match="shape"):
        tool.decompose_grid(x,eta,1+x,r,np.zeros((4,3)),np.zeros((3,3)),np.ones(5),np.zeros(5),**common)
    r[2,1]=-1
    with pytest.raises(ValueError,match="positive"):
        tool.decompose_grid(x,eta,1+x,r,np.zeros((3,3)),np.zeros((3,3)),np.ones(5),np.zeros(5),**common)


def test_empty_origin_summary_has_no_fabricated_location(budget):
    values=np.ones((2,3));mask=np.zeros((2,3),dtype=bool)
    row=tool.summary(values,mask,np.array([1.,2.]),np.array([.2,.4,.6]),values,budget)
    assert row==dict(count=0,minimum=None,maximum=None,maximum_absolute=None,rms=None,location=None)


def synthetic_record():
    x=np.linspace(.4,.5,1025);saved=x[::8];W=1+4*(x-.4)
    r=2+x*x;src=dict(r=r,r_x=2*x,r_xx=np.full(1025,2.),log_C_x=np.zeros(1025),
        c=np.zeros(1025),c_x=np.zeros(1025),source_constraint=np.full(1025,2.),x=x)
    record=dict(geometry=dict(support_x=[.431,.473]),cache=dict(x_nodes=x,gauge_rows=dict(x=x,W=W),
        source_rows={k:copy.deepcopy(src) for k in ("zero","finite")}),cases=[],diagnostics=dict(histories=[]))
    for label,alpha in (("zero",0.),("finite",.01)):
        snapshots=[dict(x=float(at),physical=dict(increment=dict(r=(.1*W[8*i]*np.linspace(0,1,33)[1:-1]).tolist()))) for i,at in enumerate(saved)]
        record["cases"].append(dict(amplitude=label,alpha=alpha,cells=32,steps=512,snapshots=snapshots))
        record["diagnostics"]["histories"].append({"departure"+str(nt):dict(source_constraint=np.full(nt,2.),
            E_x=np.zeros((nt-2,31)),E_Ubar=np.zeros((nt-2,31)),settings=dict(time_nodes=nt,eta_nodes=33,analytic_vacuum_jets=True)) for nt in (129,65)})
    return record,{k:saved.copy() for k in ("zero","finite")}


def test_full_saved_indices_and_ancestor_knots_without_large_arithmetic_replay(budget,monkeypatch):
    record,knots=synthetic_record();calls=[]
    def decompose(x,eta,W,r,Ex,EU,rs,Es,**kwargs):
        assert kwargs["budget"] is budget
        calls.append((x.copy(),r.copy(),rs.copy(),Es.copy()))
        return dict(time_nodes=len(x))
    monkeypatch.setattr(tool,"decompose_grid",decompose)
    monkeypatch.setattr(tool,"compare_sampling",lambda f,c,b:dict(sizes=[f["time_nodes"],c["time_nodes"]]))
    result=tool.analyze(record,knots,budget)
    assert [len(c[0]) for c in calls]==[129,65,129,65]
    np.testing.assert_array_equal(calls[1][0],calls[0][0][::2])
    np.testing.assert_array_equal(calls[1][1],calls[0][1][::2])
    np.testing.assert_array_equal(calls[0][1][0],calls[0][2][0]+.1*np.linspace(0,1,33)[1:-1])
    assert result[1]["source"]["all_saved_cuts_are_original_Hermite_knots"]
    record["diagnostics"]["histories"][0]["departure65"]["source_constraint"][1]+=.01
    with pytest.raises(ValueError,match="same source residual"):
        tool.analyze(record,knots,budget)
    knots["zero"][2]=np.nextafter(knots["zero"][2],1.)
    with pytest.raises(ValueError,match="Hermite knots"):
        tool.analyze(record,knots,budget)


def test_pinned_hashes_all86_bindings_and_no_scientific_validator(tmp_path,budget,monkeypatch):
    names=tuple(["bound"+str(i)+".txt" for i in range(85)]+[tool.BASE_INPUT])
    monkeypatch.setattr(tool.producer,"SOURCES",names)
    contents={name:name.encode() for name in names}
    contents[tool.BASE_INPUT]=tool.encoded(dict(initial_data={k:dict(exact_rows=dict(x=np.linspace(.4,.5,257).tolist())) for k in ("zero","finite")}))
    for name,data in contents.items():
        path=tmp_path/name;path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(data)
    raw=tool.encoded(dict(schema=tool.producer.SCHEMA,settings=tool.producer.SETTINGS,runtime=tool.runtime_versions()))
    manifest=dict(schema=tool.producer.SCHEMA+"-manifest",source_bindings={k:tool.digest(v) for k,v in contents.items()},
        outputs={"refinement.json":dict(sha256=tool.digest(raw),bytes=len(raw))})
    mr=tool.encoded(manifest)
    monkeypatch.setattr(tool,"INPUT_SHA256",tool.digest(raw));monkeypatch.setattr(tool,"MANIFEST_SHA256",tool.digest(mr))
    target=tmp_path/tool.INPUT;target.parent.mkdir(parents=True)
    target.write_bytes(raw);target.with_name("manifest.json").write_bytes(mr)
    def forbidden(*args,**kwargs):raise AssertionError("scientific replay prohibited")
    monkeypatch.setattr(tool.producer,"validate_record",forbidden)
    assert tool.authenticate(tmp_path,budget)[2]["authenticated_input_files"]==86
    (tmp_path/names[0]).write_bytes(b"changed")
    with pytest.raises(ValueError,match=r"^saved data/source binding changed: bound0\.txt$"):
        tool.authenticate(tmp_path,budget)
    manifest["source_bindings"][names[0]]=tool.digest(b"changed")
    target.with_name("manifest.json").write_bytes(tool.encoded(manifest))
    with pytest.raises(ValueError,match="pinned"):tool.authenticate(tmp_path,budget)


def test_audit_uses_one_budget_and_only_saved_analysis(tmp_path,budget,monkeypatch):
    for name in tool.OWN_SOURCES:
        path=tmp_path/name;path.parent.mkdir(parents=True,exist_ok=True);path.write_text(name)
    calls=[];identity=dict(pinned=True)
    def authenticate(root,b):calls.append(b);return {},{},identity
    monkeypatch.setattr(tool,"authenticate",authenticate)
    def analyze(record,knots,b):
        assert b is budget;b.check(7);return [dict(manufactured=True)]
    monkeypatch.setattr(tool,"analyze",analyze)
    before,clock=budget.evaluations,budget.started
    row=tool.audit(tmp_path,budget=budget)
    assert calls==[budget,budget] and budget.started==clock and budget.evaluations-before==7
    assert row["execution"]["evaluated_points"]==7
    assert row["execution"]["stage_calls"]==row["execution"]["source_evaluations"]==0


def test_creation_check_tampering_and_failure_receipt(tmp_path,budget,monkeypatch):
    original=tool.shared_budget
    monkeypatch.setattr(tool,"shared_budget",lambda value=None:budget if value is None else original(value))
    def audit(root,*,budget):
        budget.check(3)
        return dict(schema=tool.SCHEMA,value=[1.,-2.],execution=dict(evaluated_points=3,cpu_seconds=time.process_time()-budget.started))
    monkeypatch.setattr(tool,"audit",audit)
    output=tmp_path/"audit.json"
    tool.create(output,root=tmp_path)
    assert tool.check(output,root=tmp_path)["stage_calls"]==0
    with pytest.raises(FileExistsError):tool.create(output,root=tmp_path)
    row=json.loads(output.read_text());row["value"][1]=-3.;output.write_bytes(tool.encoded(row))
    with pytest.raises(ValueError,match="recomputation"):tool.check(output,root=tmp_path)
    def fail(*args,**kwargs):raise ValueError("manufactured first failure")
    monkeypatch.setattr(tool,"audit",fail)
    rejected=tmp_path/"failed.json"
    with pytest.raises(ValueError) as error:tool.create(rejected,root=tmp_path)
    assert not rejected.exists()
    assert json.loads(open(error.value.failure_receipt).read())["status"]=="FAILED_NOT_COMPLETED"


@pytest.mark.parametrize("field,value",(("cpu_seconds",float("nan")),("maximum_evaluations",1000001),
    ("started",True),("started",-1.),("evaluations",-1)))
def test_invalid_allowance_does_not_reset_or_expand(field,value):
    bad=tool.Budget(maximum_evaluations=20,cpu_seconds=30.)
    setattr(bad,field,value)
    with pytest.raises(ValueError):tool.shared_budget(bad)


def test_exhaustion_and_nonfinite_values_stop_immediately(budget):
    budget.check(2)  # include the deliberately exhausted child request in the aggregate allowance
    child=tool.Budget(maximum_evaluations=1,cpu_seconds=30.)
    with pytest.raises(RuntimeError,match="budget exhausted"):
        tool.array([1.,2.],(2,),"deliberate exhaustion",child)
    with pytest.raises(ValueError,match="nonfinite"):
        tool.array([float("nan")],(1,),"bad",budget)
    with pytest.raises(ValueError,match="shape"):
        tool.array([True,1.],(2,),"bad",budget)
