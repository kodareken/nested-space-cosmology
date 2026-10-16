"""Fixed-case dispatch and authenticated ledger lifecycle; synthetic inputs only."""
import copy
import json
from pathlib import Path

import numpy as np
import pytest

import analyze_nsc_atlas_cap_mass_balance as calc
import diagnose_nsc_atlas_cap_mass_balance as record
from test_nsc_atlas_cap_mass_derivatives import budget, report_allowance
from test_nsc_atlas_cap_mass_direction_orchestration import inputs, manufactured_metric, manufactured_fields
from test_diagnose_nsc_atlas_cap_mass_direction import pinned_fixture


def install(monkeypatch):
    seen = dict(metric=[], rhs=[], tangent=[])
    def metric(*a, **kw):
        seen["metric"].append((a[1],kw["T"],kw["alpha"],kw["quadrature_order"]))
        return manufactured_metric(*a, **kw)
    def rhs(m, *, T_x, radius_x, incoming_w_plus, budget):
        seen["rhs"].append((T_x,radius_x,incoming_w_plus))
        budget.check(128)
        return dict(rates=dict(Ip=np.full(64,.001/64), Iq=np.full(64,-.002/64)),
                    p_flux=np.zeros(65), q_flux=np.zeros(65))
    def tangent(m, direction, *, budget):
        seen["tangent"].append(copy.deepcopy(direction))
        R,T = m["radius"],m["T"]
        P,Q = (np.sum(m["state"][k]) for k in ("Ip","Iq"))
        Nr = 2*R*T+m["quadrature_order"]*1e-5
        Nt = R*R+3*T*T
        D = Nr*direction["R"]+Nt*direction["T"]+Q*.001-P*.002
        budget.check(260)
        return dict(forward=D,reverse=D,N_R=Nr,N_T=Nt,Ip_weights=np.full(64,Q),Iq_weights=np.full(64,P),
                    contraction_roundoff_allowance=1e-12)
    def ledger(m, rhs, tan, *, R_x, T_x, budget):
        radius,clock = tan["N_R"]*R_x,tan["N_T"]*T_x
        interior = tan["forward"]-radius-clock
        budget.check(16)
        return dict(radius=radius,clock=clock,centre=0.,wall=0.,interior=interior,reaction=0.,
            total=tan["forward"],telescoping_defect=0.,interior_split_defect=0.,
            maximum_face_split_defect=0.,interior_central=interior-.001,interior_upwind=.001,
            incoming_mass_rate_derivative=.02)
    for name,fn in dict(background_cap_metric=metric,background_cap_rhs=rhs,mass_tangent=tangent,
                        mass_balance_ledger=ledger,background_cap_profile=manufactured_fields).items():
        monkeypatch.setattr(calc,name,fn)
    return seen


def fixtures(budget):
    current,cap,cache = inputs(budget)
    prior = dict(data=dict(schema="nsc-atlas-cap-mass-direction-analysis-v1", cases=[dict(
        amplitude=c["amplitude"], alpha=c["alpha"], orders={str(q):dict(widths=[dict(width=1e-5,D_n=.3)])
        for q in (4,8)}) for c in current["cases"]]))
    return current,cap,cache,prior


def test_four_frozen_ledgers_use_saved_rates_and_keep_inputs(monkeypatch,budget):
    data = fixtures(budget)
    before = copy.deepcopy(data)
    seen = install(monkeypatch)
    start = budget.evaluations
    result = calc.analyze(*data,budget)
    assert result["evaluated_points"] == budget.evaluations-start == sum(result["execution"]["charges"].values())
    for k,v in calc.EXPECTED_CALLS.items(): assert result["execution"][k] == v
    assert [len(v) for v in seen.values()] == [4,4,4]
    for m,r,d in zip(seen["metric"],seen["rhs"],seen["tangent"]):
        assert m[:2] == (.49,.08)
        assert r == (1.2,(.3+m[3]*.001)/2,.25+m[2])
        assert d["R"] == r[1] and d["T"] == r[0]
    assert data == before
    for c in result["cases"]:
        for row in c["orders"].values():
            assert row["B_mass_without_upwind"] == pytest.approx(row["B_mass"]-.001)
            assert row["actual_B_mass"] == c["alpha"]*row["B_mass"]
            for radial in row["radial"]: assert abs(radial["identity_defect"]) < 1e-14


def test_identity_failure_stops_at_first_case(monkeypatch,budget):
    data = fixtures(budget); seen = install(monkeypatch)
    original = calc.mass_balance_ledger
    def fail(*args,**kwargs):
        row = original(*args,**kwargs); row["total"] += .1
        return row
    monkeypatch.setattr(calc,"mass_balance_ledger",fail)
    with pytest.raises(ValueError,match="identity failed") as caught:
        calc.analyze(*data,budget)
    assert [len(v) for v in seen.values()] == [1,1,1]
    assert caught.value.cap_mass_balance_failure["phase"]["name"] == "departure flux ledger"
    evidence = caught.value.cap_mass_balance_failure["identity_failure"]
    assert abs(evidence["defects"]["derivative_ledger"]) == pytest.approx(.1)
    assert evidence["rounding_allowance"] < 1e-12


@pytest.mark.parametrize("difference,accepted", [(1e-14,True),(.01,False)])
def test_saved_reconstruction_tolerance_is_bounded(difference,accepted):
    if accepted: calc.same_values([1.],[1.+difference],"fixture")
    else:
        with pytest.raises(ValueError,match="reconstruction differs"):
            calc.same_values([1.],[1.+difference],"fixture")


def owned(root):
    for name in record.OWN_SOURCES:
        path = root/name; path.parent.mkdir(parents=True,exist_ok=True)
        path.write_bytes(b"synthetic frozen producer")


def synthetic_record(monkeypatch,root,budget):
    owned(root)
    calls = []
    def auth(root,shared,sources):
        assert shared is budget
        calls.append("auth")
        return (1,2,3,4),dict(fixture=True)
    def analyze(*args):
        assert args == (1,2,3,4,budget)
        calls.append("analysis"); budget.check(7)
        return dict(mass=.2,evaluated_points=7,execution=dict(calc.EXPECTED_CALLS,
            charges=dict(input_and_profile=1,metrics=1,rhs=1,tangents=3,ledger_analysis=1)))
    monkeypatch.setattr(record,"authenticate",auth)
    monkeypatch.setattr(calc,"analyze",analyze)
    return calls


def test_create_check_and_no_overwrite(monkeypatch,tmp_path,budget):
    calls = synthetic_record(monkeypatch,tmp_path,budget)
    path = tmp_path/"result.json"
    saved = record.operation(path,root=tmp_path,budget=budget)
    assert calls == ["auth","analysis","analysis"]
    assert saved["execution"]["evaluated_points"] == 14
    assert record.operation(path,checking=True,root=tmp_path,budget=budget)["evaluated_points"] == 7
    raw = path.read_bytes()
    with pytest.raises(FileExistsError): record.operation(path,root=tmp_path,budget=budget)
    assert raw == path.read_bytes()


@pytest.mark.parametrize("kind",["mass","count","producer"])
def test_record_drift_rejected_with_receipt(monkeypatch,tmp_path,budget,kind):
    calls = synthetic_record(monkeypatch,tmp_path,budget)
    path = tmp_path/"result.json"
    saved = record.operation(path,root=tmp_path,budget=budget)
    if kind=="mass": saved["analysis"]["mass"] += .01
    elif kind=="count": saved["execution"]["hash_authentication_passes"] = True
    else: (tmp_path/record.OWN_SOURCES[-1]).write_bytes(b"producer drift")
    raw = record.encoded(saved); path.write_bytes(raw)
    with pytest.raises(ValueError) as caught:
        record.operation(path,checking=True,root=tmp_path,budget=budget)
    failure = json.loads(Path(caught.value.failure_receipt).read_bytes())
    assert failure["phase"] == "check reconstruction" and failure["retries"] == 0
    assert path.read_bytes() == raw
    assert calls.count("analysis") == (2 if kind=="producer" else 3)


def test_hash_only_prior_authentication_and_summary_pin(monkeypatch,tmp_path,budget):
    _,cap = pinned_fixture(tmp_path,monkeypatch,budget)
    owned(tmp_path)
    path = tmp_path/record.PRIOR; path.parent.mkdir(parents=True,exist_ok=True)
    raw = record.encoded(dict(data=dict(fixture=True))); path.write_bytes(raw)
    monkeypatch.setattr(record,"PRIOR_SHA256",record.digest(raw))
    def forbidden(*a,**kw): raise AssertionError("previous scientific replay")
    for name in ("analyze","check","_validation"):
        monkeypatch.setattr(record.previous,name,forbidden)
    start = budget.evaluations
    data,identity = record.authenticate(tmp_path,budget,record.bindings(tmp_path))
    assert data[1] == cap and data[3]["data"]["fixture"]
    assert budget.evaluations == start and identity["previous_analysis_replays"] == 0
    path.write_bytes(record.encoded(dict(data=dict(fixture=False))))
    with pytest.raises(ValueError,match="summary mismatch"):
        record.authenticate(tmp_path,budget,record.bindings(tmp_path))
