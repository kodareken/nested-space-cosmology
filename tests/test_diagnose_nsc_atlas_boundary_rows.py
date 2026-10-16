"""Selected-row algebra/provenance controls, with no retained physical replay."""
import copy
import json
import time

import numpy as np
import pytest

import diagnose_nsc_atlas_boundary_rows as tool


@pytest.fixture(scope="module")
def budget():
    shared=tool.Budget(maximum_evaluations=100000,cpu_seconds=30.)
    yield shared
    shared.check()
    print("boundary-row controls",shared.evaluations,"points",time.process_time()-shared.started,"process CPU seconds")


def test_polynomial_gradients_and_signed_transition_identity(budget):
    eta=np.linspace(0,1,65)
    z=np.stack((0*eta+.25,2*eta-.3,3*eta**2-eta))
    before=z.copy(); row=tool.gradient_choices(z,budget)
    np.testing.assert_allclose(row["quadratic"],np.stack((0*eta,0*eta+2,6*eta-1)),atol=5e-14)
    np.testing.assert_allclose(row["identity_defect"],0,atol=5e-14)
    np.testing.assert_array_equal(z,before)
    z=np.stack((eta**3,2*eta**3,-eta**3)); row=tool.gradient_choices(z,budget)
    h=1/64; third=np.array([6.,12.,-6.])*h**3
    np.testing.assert_array_equal(row["third_difference"],third)
    np.testing.assert_array_equal(row["quadratic_minus_alternative"],np.stack((-third/(2*h),third/(2*h)),axis=1))
    assert np.any(row["quadratic_minus_first"]!=0)


def test_penultimate_impulse_wall_curvature_and_stride_omission(budget):
    h=1/64; width=.3; amplitude=2**-20
    z=np.zeros((3,65)); z[:,63]=amplitude
    for points,coefficient in ((3,-5.),(5,-26/3)):
        native=tool.spatial.backward_jets(z,points=points,stride=1,budget=budget)
        wide=tool.spatial.backward_jets(z,points=points,stride=2,budget=budget)
        np.testing.assert_allclose(native["second"]/width**2,coefficient*amplitude/(h*width)**2,rtol=2e-15)
        np.testing.assert_array_equal(wide["second"],0.)
        assert 63 in native["second_indices"] and 63 not in wide["second_indices"]


def test_source_offsets_five_term_scalar_factor_and_vacuum_cancellation(budget):
    U=np.array([.1,.2,.3]); W=.4; a=.03; b=-.04; alpha=.2
    # Exponential vacuum r0 makes its nonzero curvature/lapse cancel exactly.
    r0=np.exp(U); L0=U; offset_r=.2; offset_c=.07
    r=r0+offset_r+a*U**2; c=offset_c+.06*U**2
    vacuum=dict(r_Ubar=r0,r_Ubar_Ubar=r0,log_C_Ubar=np.ones(3))
    first=dict(r=2*a*U*W,L=2*b*U*W,c=.12*U*W)
    got=tool.constraint_terms(dict(r=r,c=c),vacuum,first,np.full(3,2*a*W**2),W,alpha,budget)
    rU=r0+2*a*U; LU=1+2*b*U; u=c/r; uU=(.12*U-u*rU)/r
    expected=r0+2*a-LU*rU+alpha*r*uU*uU/(2*(1+alpha*u*u/6))
    np.testing.assert_allclose(got["E_Ubar"],expected,atol=3e-16)
    np.testing.assert_allclose(got["terms"]["background_lapse"],-2*a*U,atol=1e-17)
    np.testing.assert_allclose(got["terms"]["background_radius"],-r0*2*b*U,atol=1e-17)
    np.testing.assert_allclose(got["terms"]["departure_product"],-4*a*b*U*U,atol=1e-18)
    np.testing.assert_allclose(got["algebraic_identity_defect"],0,atol=3e-16)
    zero=tool.constraint_terms(dict(r=r0,c=np.zeros(3)),vacuum,
        {k:np.zeros(3) for k in tool.FIELDS},np.zeros(3),W,0.,budget)
    np.testing.assert_array_equal(zero["E_Ubar"],0.)
    np.testing.assert_array_equal(zero["full_E_Ubar"],0.)


def manufactured_inputs(monkeypatch,budget):
    """Only selected snapshots carry data; all other rows are unused sentinels."""
    xc=np.linspace(.4,.48,1025); x=xc[::8]; eta=np.linspace(0,1,65); W=.2+.2*(xc-.4)
    cache=dict(x_nodes=xc,gauge_rows=dict(x=xc,W=W,W_x=np.full(1025,.2),W_xx=np.zeros(1025)),source_rows={})
    records=[]; vacuum_calls=[]
    def vacuum(U,xx,*,budget):
        U=np.asarray(U);xx=np.asarray(xx);budget.check(U.size);vacuum_calls.append(U.size)
        return dict(r=xx+U,L=np.zeros_like(U),r_Ubar=np.ones_like(U),
                    r_Ubar_Ubar=np.zeros_like(U),log_C_Ubar=np.zeros_like(U))
    monkeypatch.setattr(tool.spatial,"vacuum_slice",vacuum)
    for label in ("zero","finite"):
        offset=.01 if label=="finite" else 0.
        cache["source_rows"][label]=dict(x=xc,r=xc+offset,log_C=np.full(1025,offset/2),c=np.full(1025,.02))
    for method in ("quadratic","first_order"):
        record=dict(cases=[],diagnostics=dict(constraint_transport=[],wall_curvature=dict(cases=[])))
        for ci,label in enumerate(("zero","finite")):
            alpha=0. if ci==0 else .001; offset=.01 if ci else 0.
            multiplier=(1 if method=="quadratic" else 2)*(ci+1)
            coefficients={k:multiplier*v for k,v in zip(tool.FIELDS,(2**-12,-2**-13,2**-11))}
            snapshots=[dict(x=float(t)) for t in x]
            EU=np.zeros((127,63)); wallarrays={k:np.zeros(129) for k in ("curvature","lapse_product","scalar","E_Ubar","weighted")}
            for cut in tool.CUTS:
                t=x[cut];w=W[8*cut];U=w*eta
                Z={k:a*eta**2 for k,a in coefficients.items()}
                physical=dict(r=t+U+offset+Z["r"],L=offset/2+Z["L"],c=.02+Z["c"])
                snapshots[cut]=dict(x=float(t),state=dict(departure={k:v[1:-1] for k,v in Z.items()}),
                    physical=dict(wall={k:v[-1] for k,v in physical.items()},
                        increment=dict(r=physical["r"][1:-1]-(t+offset),
                                       L=physical["L"][1:-1]-offset/2,c=physical["c"][1:-1]-.02)))
                rU=1+2*coefficients["r"]*eta/w; LU=2*coefficients["L"]*eta/w
                cU=2*coefficients["c"]*eta/w; r=physical["r"]; u=physical["c"]/r; uU=(cU-u*rU)/r
                curvature=np.full(65,2*coefficients["r"]/w**2)
                lapse=-LU*rU; scalar=alpha*r*uU*uU/(2*(1+alpha*u*u/6)); E=curvature+lapse+scalar
                if cut<128: EU[cut-1]=E[1:-1]
                for k,v in dict(curvature=curvature,lapse_product=lapse,scalar=scalar,E_Ubar=E,weighted=r*E).items():
                    wallarrays[k][cut]=v[-1]
                budget.check(65*8)
            record["cases"].append(dict(amplitude=label,cells=64,steps=512,alpha=alpha,snapshots=snapshots))
            grid=dict(eta=eta[1:-1],arrays=dict(E_Ubar=EU))
            record["diagnostics"]["constraint_transport"].append(dict(grids={"129" if method=="quadratic" else "native64_129":grid}))
            estimates={f"p{p}_s{s}":dict(arrays=copy.deepcopy(wallarrays)) for p in (3,5) for s in (1,2)}
            record["diagnostics"]["wall_curvature"]["cases"].append(dict(estimates=estimates if method=="quadratic" else dict(fine64=estimates)))
        records.append(record)
    return *records,cache,vacuum_calls


def test_fixed_three_cut_orchestration_same_data_and_same_stencil(monkeypatch,budget):
    current,old,cache,calls=manufactured_inputs(monkeypatch,budget)
    before=copy.deepcopy(current); start=budget.evaluations
    got=tool.analyze(current,old,cache,budget)
    assert calls==[54] and got["analytic_vacuum_points"]==54
    assert got["evaluated_points"]==budget.evaluations-start
    assert got["cuts"]==[63,68,128] and got["patch_indices"]["wall"]==list(range(56,65))
    assert tool.encoded(current)==tool.encoded(before)
    assert [(c["method"],c["amplitude"]) for c in got["cases"]]==[
        (m,a) for m in ("first_order","quadratic") for a in ("zero","finite")]
    for case in got["cases"]:
        for patch in case["patches"].values():
            assert np.shape(patch["same_data_gradients"]["quadratic"])==(3,3,9)
            assert patch["retained_E_Ubar"][-1] is None
            np.testing.assert_allclose(patch["retained_difference"][:2],0,atol=1e-12)
        np.testing.assert_allclose(case["closure_identity"]["identity_defect"],0,atol=1e-15)
    # The two evolved arrays differ; both are evaluated with the same centered diagnostic.
    for ci,difference in enumerate(got["same_stencil_method_difference"]):
        expected=(np.asarray(got["cases"][ci+2]["patches"]["wall"]["centered"]["E_Ubar"])
                  -np.asarray(got["cases"][ci]["patches"]["wall"]["centered"]["E_Ubar"]))
        np.testing.assert_array_equal(difference["patches"]["wall"]["E_Ubar"],expected)
        assert np.max(abs(expected))>0


def test_local_invalid_data_and_shared_exhaustion(budget):
    with pytest.raises(ValueError): tool.gradient_choices(np.zeros((65,3)),budget)
    with pytest.raises(ValueError): tool.gradient_choices(np.full((3,65),np.nan),budget)
    bad=tool.Budget(started=time.process_time()+100.,cpu_seconds=30.)
    with pytest.raises(ValueError,match="future"): tool.shared_budget(bad)
    class Tiny(tool.Budget):
        def check(self,count=0):
            budget.check(count); super().check(count)
    tiny=Tiny(maximum_evaluations=2,cpu_seconds=30.,started=budget.started)
    with pytest.raises(RuntimeError,match="exhausted"):
        tool.gradient_choices(np.zeros((3,65)),tiny)
    assert tiny.evaluations==195


def auth_fixture(tmp_path,monkeypatch):
    own={name:b"test producer" for name in tool.OWN_SOURCES}
    identity=dict(test="cache")
    baseline=dict(geometry={},profile_inputs={},cache_identity=identity)
    ancestor=dict(cache={})
    files={tool.producer.INPUT:tool.encoded(baseline),tool.producer.CACHE_INPUT:tool.encoded(ancestor)}
    files.update({f"fixture/{i}.txt":b"frozen" for i in range(102)})
    for name,data in {**files,**own}.items():
        path=tmp_path/name; path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(data)
    current=dict(baseline,schema=tool.producer.SCHEMA,settings=tool.producer.SETTINGS,runtime=tool.runtime_versions())
    raw=tool.encoded(current)
    manifest=dict(schema=tool.producer.SCHEMA+"-manifest",outputs={"operator.json":dict(sha256=tool.digest(raw),bytes=len(raw))},
                  source_bindings={name:tool.digest(data) for name,data in files.items()})
    mr=tool.encoded(manifest);path=tmp_path/tool.INPUT;path.parent.mkdir(parents=True,exist_ok=True)
    path.write_bytes(raw);path.with_name("manifest.json").write_bytes(mr)
    monkeypatch.setattr(tool,"INPUT_SHA256",tool.digest(raw));monkeypatch.setattr(tool,"MANIFEST_SHA256",tool.digest(mr))
    monkeypatch.setattr(tool.producer,"SOURCES",tuple(files))
    monkeypatch.setattr(tool.producer,"cache_identity",lambda r:identity)
    return files


def test_hash_only_authentication_and_binding_mutation(tmp_path,monkeypatch,budget):
    files=auth_fixture(tmp_path,monkeypatch)
    def forbidden(*args,**kwargs): raise AssertionError("scientific replay")
    monkeypatch.setattr(tool.producer,"calculate",forbidden)
    monkeypatch.setattr(tool.producer,"load_inputs",forbidden)
    _,_,_,identity=tool.authenticate(tmp_path,budget)
    assert identity["authenticated_input_files"]==104 and not identity["recursive_validation"]
    (tmp_path/"fixture/0.txt").write_bytes(b"changed")
    with pytest.raises(ValueError,match="binding changed"): tool.authenticate(tmp_path,budget)


def test_creation_check_used_budget_and_output_mutation(tmp_path,monkeypatch,budget):
    for name in tool.OWN_SOURCES:
        p=tmp_path/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_text("owned")
    identity=dict(test="authenticated"); inputs=({}, {}, {})
    monkeypatch.setattr(tool,"authenticate",lambda root,b: (*inputs,identity))
    def analyze(a,b,c,shared):
        assert shared is budget
        shared.check(7);return dict(selected=[63,68,128],evaluated_points=7)
    monkeypatch.setattr(tool,"analyze",analyze)
    path=tmp_path/"result.json";initial=budget.evaluations
    created=tool.create(path,root=tmp_path,budget=budget)
    assert created["execution"]["initial_evaluations"]==initial
    assert budget.evaluations-initial==14 and created["execution"]["evaluated_points"]==14
    assert tool.check(path,root=tmp_path,budget=budget)["evaluated_points"]==7
    with pytest.raises(FileExistsError): tool.create(path,root=tmp_path,budget=budget)
    changed=json.loads(path.read_text());changed["analysis"]["selected"][0]=62
    path.write_bytes(tool.encoded(changed))
    with pytest.raises(ValueError,match="arithmetic mismatch"): tool.check(path,root=tmp_path,budget=budget)
    path.write_bytes(tool.encoded(created))
    (tmp_path/tool.OWN_SOURCES[0]).write_text("changed")
    with pytest.raises(ValueError,match="declaration mismatch"): tool.check(path,root=tmp_path,budget=budget)


def test_failure_receipt_stops_without_retry(tmp_path,monkeypatch,budget):
    calls=[]
    def fail(root,shared): calls.append(1);raise ValueError("injected bad input")
    monkeypatch.setattr(tool,"authenticate",fail)
    with pytest.raises(ValueError,match="injected") as error:
        tool.create(tmp_path/"result.json",root=tmp_path,budget=budget)
    assert calls==[1] and not (tmp_path/"result.json").exists()
    receipt=json.loads(open(error.value.failure_receipt).read())
    assert receipt["phase"]=="input authentication" and receipt["retries"]==0
