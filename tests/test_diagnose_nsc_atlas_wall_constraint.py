"""Bounded independent wall-jet and causal-transport controls; no evolution."""
import json
import time

import numpy as np
import pytest

import diagnose_nsc_atlas_wall_constraint as tool


@pytest.fixture(scope="module")
def budget():
    shared = tool.Budget(maximum_evaluations=200000,cpu_seconds=30.)
    yield shared
    shared.check()
    print("wall constraint controls",shared.evaluations,time.process_time()-shared.started)


@pytest.mark.parametrize("points,stride",((3,1),(3,2),(5,1),(5,2)))
def test_backward_polynomial_oracles(points,stride,budget):
    eta=np.linspace(0,1,33)
    # Constants and low-degree polynomial derivatives are independent exact oracles.
    first_degree=2 if points==3 else 4
    second_degree=3 if points==3 else 4
    values=np.vstack((eta**first_degree,eta**second_degree,np.ones(33)))
    row=tool.backward_jets(values,points=points,stride=stride,budget=budget)
    np.testing.assert_allclose(row["first"][0],first_degree,atol=4e-13)
    np.testing.assert_allclose(row["second"][1],second_degree*(second_degree-1),atol=2e-12)
    np.testing.assert_allclose([row["first"][2],row["second"][2]],0,atol=1e-12)
    assert row["first_indices"][0]==row["second_indices"][0]==32
    assert np.min(row["second_indices"])==32-stride*(len(row["second_indices"])-1)


def test_full_wall_constraint_has_physical_scale_and_scalar_factor(budget):
    eta=np.linspace(0,1,33);W=np.array([1.,2.,3.]);alpha=.3
    Z=dict(r=np.broadcast_to(.2*eta**3,(3,33)),L=np.broadcast_to(.4*eta**2,(3,33)),
           c=np.broadcast_to(.5*eta**4,(3,33)))
    vac=dict(r_Ubar=np.array([.1,.2,.3]),r_Ubar_Ubar=np.array([.01,.04,.09]),
             log_C_Ubar=np.array([.1,.2,.3]))
    wall=dict(r=np.array([2.,2.5,3.]),L=np.zeros(3),c=np.array([.7,.8,.9]))
    got=tool.wall_estimate(Z,wall,vac,W,alpha,points=5,stride=2,budget=budget)
    rU=vac["r_Ubar"]+.6/W;rUU=vac["r_Ubar_Ubar"]+1.2/W**2
    LU=vac["log_C_Ubar"]+.8/W;cU=2/W
    u=wall["c"]/wall["r"];uU=(cU-u*rU)/wall["r"]
    expected=rUU-LU*rU+alpha*wall["r"]*uU*uU/(2*(1+alpha*u*u/6))
    np.testing.assert_allclose(got["E_Ubar"],expected,rtol=2e-12,atol=3e-13)
    np.testing.assert_allclose(got["B"],wall["r"]*expected,rtol=2e-12,atol=1e-12)
    np.testing.assert_allclose(sum(got["terms"].values()),got["E_Ubar"],atol=0)
    assert np.min(got["terms"]["scalar"])>0


def test_analytic_vacuum_zero_departure_recovers_constraint(budget):
    x=np.array([.4,.45,.49]);W=np.array([.1,.15,.2])
    vacuum=tool.vacuum_slice(W,x,budget=budget)
    wall=dict(r=vacuum["r"],L=vacuum["L"],c=np.zeros(3))
    Z={k:np.zeros((3,33)) for k in ("r","L","c")}
    for _,points,stride in tool.VARIANTS:
        got=tool.wall_estimate(Z,wall,vacuum,W,.01,points=points,stride=stride,budget=budget)
        np.testing.assert_array_equal(got["E_Ubar"],np.zeros(3))


def test_causal_linear_birth_has_exact_nodes_and_no_future_bracket(budget):
    W=np.array([1.,2.,3.,4.,5.]);U=np.array([[1.5,2.,.9],[2.5,3.,6.]])
    obs=np.array([1,2]);mask=np.ones_like(U,dtype=bool)
    fine=tool.causal_plan(W,U,obs,mask,sampling=1,budget=budget)
    coarse=tool.causal_plan(W,U,obs,mask,sampling=2,budget=budget)
    np.testing.assert_array_equal(fine["supported"],[[True,True,False],[True,True,False]])
    np.testing.assert_array_equal(coarse["supported"],[[False,False,False],[True,True,False]])
    assert fine["lower_index"][0,1]==fine["upper_index"][0,1]==1
    assert coarse["lower_index"][1,1]==coarse["upper_index"][1,1]==2
    B=3*W-2
    pred=tool.transport_values(B,fine,budget)
    np.testing.assert_allclose(pred[fine["supported"]],(3*U-2)[fine["supported"]],atol=0)
    values=tool.masked(pred,coarse["supported"])
    assert values[0,0] is None and values[1,0]==5.5
    with pytest.raises(ValueError,match="causal"):
        tool.causal_plan(W,U,np.array([1.5,2.]),mask,sampling=1,budget=budget)


def test_spatially_disjoint_stencils_end_at_centre22():
    wall=np.ones((3,31),dtype=bool);common=wall.copy();common[1,0]=False
    masks=tool.subset_masks(wall,common)
    np.testing.assert_array_equal(masks["disjoint_wall_stencils"][0],np.arange(1,32)<=22)
    assert not masks["disjoint_wall_stencils"][1,0]
    np.testing.assert_array_equal(masks["disjoint_wall_stencils"]|masks["overlapping_wall_stencils"],common)
    assert not np.any(masks["disjoint_wall_stencils"]&masks["overlapping_wall_stencils"])


def test_nonzero_transport_and_comparison_axes_are_kept_separate(budget):
    W=np.array([1.,2.,3.,4.,5.]);x=np.linspace(0,1,5);eta=np.linspace(0,1,33)[1:-1]
    U=W[1:-1,None]*eta[None,:];r=2+U/10
    trueB=2*U-1;remainder=.03*U
    grid=dict(interior_x=x[1:-1],interior_eta=eta,interior_Ubar=U,radius=r,W=W,
        arrays=dict(E_Ubar=(trueB+remainder)/r),
        masks=dict(wall_origin=W[1:-1,None]*np.linspace(0,1,33)[None,:-2]>W[0]))
    estimates={name:dict(B=2*W-1+i*.125) for i,(name,_,_) in enumerate(tool.VARIANTS)}
    row=tool.case_transport(estimates,grid,W,budget)
    supported=row["common_supported"]
    baseline=tool.VARIANTS[0][0]
    got=np.asarray(row["transported"][baseline]["129"]["residual"],dtype=float)
    np.testing.assert_allclose(got[supported],remainder[supported],atol=9e-16)
    for name,_,_ in tool.VARIANTS:
        change=np.asarray(row["sampling_comparisons"][name]["B_65_minus_129"],dtype=float)
        np.testing.assert_allclose(change[supported],0,atol=9e-16)
    for i,(name,_,_) in enumerate(tool.VARIANTS[1:],1):
        for sampling in ("129","65"):
            change=np.asarray(row["spatial_comparisons"][sampling][name]["B_minus_baseline"],dtype=float)
            np.testing.assert_allclose(change[supported],i*.125,atol=9e-16)
    assert row["coverage"]["65"]["unsupported_count"]>0
    assert row["coverage"]["129"]["unsupported_count"]==0
    assert row["coverage"]["65"]["first_unsupported"]["time_grid_index"]==1
    assert row["transported"][baseline]["129"]["summaries"]["residual"]["all_wall_origin"]["count"]==np.count_nonzero(supported)
    domain=row["comparison_mask_coverage"]["disjoint_wall_stencils"]
    assert domain["count"]>0 and domain["eta_interval"][1]<=22/32
    assert domain["x_interval"][0]>=x[1] and domain["x_interval"][1]<=x[-2]


def synthetic_history(budget):
    cachex=np.linspace(.4,.5,1025);x=cachex[::8];W=.1+x-.4
    vacuum=tool.vacuum_slice(W,x,budget=budget)
    gauge=dict(W=np.interp(cachex,x,W))
    for k in ("r","L","r_Ubar","log_C_Ubar"):
        gauge[k]=np.interp(cachex,x,vacuum[k])
    record=dict(cache=dict(x_nodes=cachex,gauge_rows=gauge,source_rows={}),cases=[])
    retained=dict(cases=[])
    eta=np.linspace(0,1,33)
    for label,alpha in (("zero",0.),("finite",.001)):
        source=dict(x=cachex,r=cachex+1e-5*cachex**2,log_C=.01*cachex,c=cachex**2)
        record["cache"]["source_rows"][label]=source
        offsets=dict(r=source["r"][::8]-x,L=source["log_C"][::8],c=source["c"][::8])
        snapshots=[]
        for i,t in enumerate(x):
            middle={k:float(i/128)*eta[1:-1]**2 for k in ("r","L","c")}
            wall={k:(vacuum[k][i] if k!="c" else 0)+offsets[k][i]+i/128 for k in middle}
            snapshots.append(dict(x=t,state=dict(departure=middle),physical=dict(wall=wall,increment=dict(r=np.ones(31)))))
        record["cases"].append(dict(amplitude=label,alpha=alpha,cells=32,steps=512,snapshots=snapshots))
        radius=source["r"][::8][1:-1,None]+np.ones((127,31))
        retained["cases"].append(dict(amplitude=label,alpha=alpha,grids={"129":dict(x=x,radius=radius)}))
    return record,retained


def test_reconstruct_wall_departure_from_source_offset_and_one_shared_vacuum(budget,monkeypatch):
    record,retained=synthetic_history(budget)
    original=tool.vacuum_slice;calls=[]
    def counted(U,x,*,budget):calls.append(len(U));return original(U,x,budget=budget)
    monkeypatch.setattr(tool,"vacuum_slice",counted)
    seen=[]
    def estimate(Z,wall,vacuum,W,alpha,*,points,stride,budget):
        seen.append((alpha,points,stride))
        expected=np.arange(129)/128
        for k in ("r","L","c"):
            np.testing.assert_allclose(Z[k][:,-1],expected,atol=3e-16)
            np.testing.assert_array_equal(Z[k][:,0],np.zeros(129))
        return dict(B=np.zeros(129))
    monkeypatch.setattr(tool,"wall_estimate",estimate)
    monkeypatch.setattr(tool,"case_transport",lambda *args:dict(manufactured=True))
    result=tool.analyze(record,retained,budget)
    assert calls==[129] and len(seen)==8 and len(result["cases"])==2
    record["cases"][1]["snapshots"][2]["x"]+=.001
    with pytest.raises(ValueError,match="saved wall cuts"):
        tool.analyze(record,retained,budget)


def test_pinned_transport_result_and_producer_changes_reject_without_replay(tmp_path,budget,monkeypatch):
    identity=dict(frozen=True);record=dict(fake=True);seen=[]
    def auth(root,b):seen.append(b);return record,None,identity
    monkeypatch.setattr(tool.prior,"authenticate",auth)
    bindings={}
    for name in tool.prior.OWN_SOURCES:
        path=tmp_path/name;path.parent.mkdir(parents=True,exist_ok=True);path.write_text(name)
        bindings[name]=tool.digest(path.read_bytes())
    row=dict(schema=tool.prior.SCHEMA,input=identity,runtime=tool.runtime_versions(),audit_source_bindings=bindings)
    raw=tool.encoded(row);monkeypatch.setattr(tool,"TRANSPORT_SHA256",tool.digest(raw))
    target=tmp_path/tool.TRANSPORT_INPUT;target.parent.mkdir(parents=True);target.write_bytes(raw)
    def forbidden(*args,**kwargs):raise AssertionError("prior diagnostic replay forbidden")
    monkeypatch.setattr(tool.prior,"audit",forbidden)
    assert tool.authenticate(tmp_path,budget)[0] is record and seen==[budget]
    name=tool.prior.OWN_SOURCES[0];(tmp_path/name).write_text("changed")
    with pytest.raises(ValueError,match="producer changed"):tool.authenticate(tmp_path,budget)
    row["audit_source_bindings"][name]=tool.digest(b"changed");target.write_bytes(tool.encoded(row))
    with pytest.raises(ValueError,match="pinned"):tool.authenticate(tmp_path,budget)


def test_shared_audit_budget_and_source_drift_after_analysis(tmp_path,budget,monkeypatch):
    for name in tool.OWN_SOURCES:
        path=tmp_path/name;path.parent.mkdir(parents=True,exist_ok=True);path.write_text(name)
    calls=[]
    def auth(root,b):calls.append(b);return {},{},dict(fixed=True)
    monkeypatch.setattr(tool,"authenticate",auth)
    def analyze(a,b,c):
        assert c is budget;c.check(9);return dict(manufactured=True)
    monkeypatch.setattr(tool,"analyze",analyze)
    clock=budget.started;before=budget.evaluations
    row=tool.audit(tmp_path,budget=budget)
    assert calls==[budget,budget] and budget.started==clock and budget.evaluations-before==9
    assert row["execution"]["evaluated_points"]==9 and row["execution"]["time_steps"]==0
    def drift(a,b,c):
        (tmp_path/tool.OWN_SOURCES[0]).write_text("changed");c.check(1);return {}
    monkeypatch.setattr(tool,"analyze",drift)
    with pytest.raises(ValueError,match="changed during"):tool.audit(tmp_path,budget=budget)


def test_creation_roundtrip_tamper_first_failure_and_no_overwrite(tmp_path,budget,monkeypatch):
    original=tool.shared_budget
    monkeypatch.setattr(tool,"shared_budget",lambda value=None:budget if value is None else original(value))
    def audit(root,*,budget):
        budget.check(3)
        return dict(schema=tool.SCHEMA,value=[1.,None,-2.],execution=dict(evaluated_points=3,cpu_seconds=time.process_time()-budget.started))
    monkeypatch.setattr(tool,"audit",audit)
    path=tmp_path/"wall.json";tool.create(path,root=tmp_path)
    assert tool.check(path,root=tmp_path)["stage_calls"]==0
    with pytest.raises(FileExistsError):tool.create(path,root=tmp_path)
    row=json.loads(path.read_text());row["value"][2]=2.;path.write_bytes(tool.encoded(row))
    with pytest.raises(ValueError,match="recomputation"):tool.check(path,root=tmp_path)
    def fail(*args,**kwargs):raise ValueError("manufactured failure")
    monkeypatch.setattr(tool,"audit",fail)
    rejected=tmp_path/"failed.json"
    with pytest.raises(ValueError) as error:tool.create(rejected,root=tmp_path)
    assert not rejected.exists()
    assert json.loads(open(error.value.failure_receipt).read())["status"]=="FAILED_NOT_COMPLETED"


@pytest.mark.parametrize("field,value",(("cpu_seconds",float("nan")),("maximum_evaluations",2000001),
    ("started",True),("started",-1.),("evaluations",-1)))
def test_invalid_shared_budget_cannot_disable_allowance(field,value):
    bad=tool.Budget(maximum_evaluations=10,cpu_seconds=30.)
    setattr(bad,field,value)
    with pytest.raises(ValueError):tool.shared_budget(bad)


def test_exhaustion_and_invalid_stencils_stop(budget):
    budget.check(66)
    child=tool.Budget(maximum_evaluations=1,cpu_seconds=30.)
    with pytest.raises(RuntimeError,match="budget exhausted"):
        tool.backward_jets(np.zeros((2,33)),points=3,stride=1,budget=child)
    with pytest.raises(ValueError):tool.backward_jets(np.zeros((2,33)),points=3,stride=True,budget=budget)
    with pytest.raises(ValueError):tool.backward_jets(np.zeros((2,4)),points=5,stride=2,budget=budget)
    with pytest.raises(ValueError):tool.backward_jets(np.full((2,33),float("nan")),points=5,stride=2,budget=budget)
