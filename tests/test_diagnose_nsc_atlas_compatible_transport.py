"""Manufactured transport oracles and input authentication; no physical run."""
import copy
import json
import time

import numpy as np
import pytest

import diagnose_nsc_atlas_compatible_transport as tool


@pytest.fixture(scope="module")
def budget():
    value = tool.Budget(maximum_evaluations=20000, cpu_seconds=30.)
    yield value
    value.check()
    print("transport audit controls", value.evaluations, time.process_time()-value.started)


def cap_case():
    R, n = .45, 8
    qa = np.arctanh(R)
    samples = []
    for i, T in enumerate(np.linspace(0., .2, 9)):
        r = R*np.linspace(0., 1., n+1)
        primitive = np.maximum(0., T+np.arctanh(r)-qa)
        samples.append(dict(x=.4+i*.01, physical=dict(wall=dict(T=T, r=R, c=T)),
                            state=dict(cap=dict(Iq=np.diff(primitive).tolist()))))
    return dict(alpha=0., amplitude="zero", cells=n, steps=16, snapshots=samples)


def test_exact_incoming_cap_and_quiet_extension(budget):
    row = tool.cap_pullback(cap_case(), budget)
    assert row["reflected_branch_upper_label"] < row["q_wall"][0]
    assert row["maxima"]["linear"]["absolute"] < 2e-16
    assert row["maxima"]["pchip"]["absolute"] < 2e-16
    assert row["rows"][0]["quiet_points"] == 8


def test_cap_oracle_detects_local_smearing_despite_exact_wall(budget):
    case = cap_case()
    case["snapshots"][4]["state"]["cap"]["Iq"][2] += .03
    case["snapshots"][4]["state"]["cap"]["Iq"][3] -= .03
    row = tool.cap_pullback(case, budget)
    assert row["maxima"]["pchip"]["absolute"] == pytest.approx(.03)
    assert row["maxima"]["pchip"]["face"] == 3
    assert row["maxima"]["pchip"]["snapshot"] == 4


@pytest.mark.parametrize("change", ("alpha", "reflection", "wall", "arrival"))
def test_cap_domain_and_matching_controls(change, budget):
    case = cap_case()
    if change == "alpha":
        case["alpha"] = .1
    elif change == "reflection":
        case["snapshots"][-1]["physical"]["wall"]["T"] = 1.
    elif change == "wall":
        case["snapshots"][3]["physical"]["wall"]["c"] += .01
    else:
        case["snapshots"][2]["physical"]["wall"]["T"] = 0.
    with pytest.raises(ValueError):
        tool.cap_pullback(case, budget)


def manufactured(n, budget):
    # Closed physical prefixes: r,C constant and canonical f linear in U.
    U, x, r, C, alpha = np.linspace(0., .4, n), .45, .65, 1.3, .3
    v = tool.vacuum_slice(U, x, budget=budget)
    f0, b, R0, S0, L0 = .3, 1.2, .8, -.4, .15
    f = f0+b*U
    scale = np.sqrt(alpha/6)
    u, rootD = np.sinh(scale*f)/scale, np.cosh(scale*f)
    X = dict(r=np.full(n, r), L=np.full(n, np.log(C)), c=r*u)
    Z = {key: X[key]-(v[key] if key != "c" else 0.) for key in tool.FIELDS}
    exact_jets = dict(r_Ubar=-v["r_Ubar"], log_C_Ubar=-v["log_C_Ubar"], c_Ubar=r*rootD*b)
    source = dict(r=r, log_C=np.log(C), c=X["c"][0], r_x=R0,
                  log_C_x=L0, c_x=u[0]*R0+r*rootD[0]*S0)
    R = R0-C*U/(4*r)
    S = S0-b*(R0*U-C*U**2/(8*r))/r
    exact = dict(r_x=R, f_x=S, c_x=u*R+r*rootD*S,
        log_C_x=L0+C*U/(2*r*r)-alpha*b*(S0*U-b*R0*U**2/(2*r)+b*C*U**3/(24*r*r)))
    return U, X, Z, exact_jets, source, exact


def test_closed_prefix_oracle_separates_jets_from_quadrature(budget):
    errors = {(method, order): [] for method in ("exact", "first", "second") for order in (2, 4)}
    for n in (9, 17, 33):
        U, X, Z, exact_jets, source, exact = manufactured(n, budget)
        for method in ("exact", "first", "second"):
            jets = exact_jets if method == "exact" else tool.departure_jets(Z, U, method)
            for order in (2, 4):
                row = tool.departure_slice(U, X, jets, source, x=.45, alpha=.3,
                    tube=tool.producer.TUBE, budget=budget, quadrature_order=order)
                errors[method, order].append(max(float(np.max(abs(row[key]-value)))
                                                 for key, value in exact.items()))
    print("crossed manufactured errors", errors)
    for order in (2, 4):
        first, second = np.array(errors["first", order]), np.array(errors["second", order])
        assert np.all(first[:-1] > first[1:])
        if order == 4:
            # The higher-order prefix isolates the first-order jet error.
            # With trapezoids the two leading errors partially cancel.
            assert np.all(first[:-1]/first[1:] > 1.8)
            assert np.all(first[:-1]/first[1:] < 2.3)
        assert np.all(second[:-1]/second[1:] > 3.4)
        assert second[-1] < first[-1]/5
    for order, minimum in ((2, 3.7), (4, 12.)):
        exact = np.array(errors["exact", order])
        assert np.all(exact[:-1]/exact[1:] > minimum)


def test_derivative_methods_and_endpoint_sign():
    U = np.linspace(0., .4, 9)
    Z = dict(r=U**2, L=3*U**2, c=-U**2)
    first, second = (tool.departure_jets(Z, U, k) for k in ("first", "second"))
    np.testing.assert_allclose(second["r_Ubar"], 2*U, atol=3e-15)
    assert first["r_Ubar"][-1] == pytest.approx(2*U[-1]-(U[1]-U[0]))
    Wx = 1.7
    d1, d2 = .3-Wx*first["r_Ubar"][-1], .3-Wx*second["r_Ubar"][-1]
    assert d2-d1 == pytest.approx(-Wx*(second["r_Ubar"][-1]-first["r_Ubar"][-1]))
    with pytest.raises(ValueError):
        tool.departure_jets(Z, U[::-1], "second")


def test_crossed_stage_freezes_cap_and_radial_longitudinal_jet(budget):
    n, x = 8, .44
    eta, W, Wx = np.linspace(0., 1., n+1), .2, 1.2
    v = tool.vacuum_slice(eta*W, x, budget=budget)
    gauge = dict(x=x, W=W, W_x=Wx, eta_plus=-.5*np.log(Wx), eta_plus_x=.1,
                 r=v["r"][-1], L=v["L"][-1], T=.01, s=.01)
    source = dict(x=x, r=x, r_x=1., log_C=0., log_C_x=0., c=0., c_x=0.)
    state = dict(departure={key: (scale*eta[1:-1]*(1-eta[1:-1])).tolist()
                           for key, scale in (("r", .001), ("L", .002), ("c", 0.))},
                 wall=dict(r=0., L=0., T=0., s=0.), cap=dict(Ip=[0.]*n, Iq=[0.]*n))
    physical = dict(increment={key: (v[key][1:-1]-source["log_C" if key == "L" else key]
                         +np.array(state["departure"][key])).tolist() if key != "c"
                         else state["departure"][key] for key in tool.FIELDS},
                    wall={**{k: gauge[k] for k in ("r", "L", "T", "s")}, "c": 0.})
    case = dict(alpha=0., amplitude="finite", cells=n,
                snapshots=[dict(x=x, state=state, physical=physical)])
    parent = dict(gauge=dict(rows={k: [v] for k, v in gauge.items()}),
                  initial_data=dict(finite=dict(snapshot_rows={k: [v] for k, v in source.items()})))
    result = tool.crossed_stage(case, parent, 0, budget)
    assert result["cap_radial_reconstructions"] == 1
    for order in (2, 4):
        a, b = result["results"][f"first_{order}"], result["results"][f"second_{order}"]
        assert a["longitudinal"]["r"] == b["longitudinal"]["r"]
        assert a["wall_rates"]["r"] == b["wall_rates"]["r"]
        assert abs(result["radial_jet_change_identity"][str(order)]) < 2e-15
        for key in tool.FIELDS:
            assert b["required_transverse"][key]-b["transverse"][key] == pytest.approx(
                b["wall_kinematic_defect"][key]/Wx, abs=2e-15)


def binding_fixture(tmp_path, monkeypatch):
    parent_name = tool.producer.previous.RELATIVE_OUTPUT+"/coupled.json"
    names = (parent_name, *[f"bound-{i}.txt" for i in range(50)])
    bindings = {}
    for name in names:
        path = tmp_path/name
        path.parent.mkdir(parents=True, exist_ok=True)
        raw = b"{}" if name == parent_name else name.encode()
        path.write_bytes(raw)
        bindings[name] = tool.digest(raw)
    record = dict(schema=tool.producer.SCHEMA, settings=tool.producer.SETTINGS,
                  parent_input=dict(path=parent_name, sha256=bindings[parent_name]))
    raw = tool.encoded(record)
    path = tmp_path/tool.INPUT
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(raw)
    manifest = dict(schema=tool.producer.SCHEMA+"-manifest", source_bindings=bindings,
                    outputs={"compatible.json": dict(sha256=tool.digest(raw), bytes=len(raw))})
    path.with_name("manifest.json").write_bytes(tool.encoded(manifest))
    monkeypatch.setattr(tool, "INPUT_SHA256", tool.digest(raw))
    monkeypatch.setattr(tool, "INPUT_MANIFEST_SHA256", tool.digest(tool.encoded(manifest)))
    monkeypatch.setattr(tool.producer, "SOURCES", names)
    return path, manifest


def test_hash_authentication_never_calls_scientific_validator(tmp_path, monkeypatch):
    binding_fixture(tmp_path, monkeypatch)
    def forbidden(*args, **kwargs):
        raise AssertionError("no physical or recursive replay")
    for name in ("validate_record", "evolve_compatible", "load_parent"):
        monkeypatch.setattr(tool.producer, name, forbidden)
    _, _, identity = tool.authenticate(tmp_path)
    assert identity["authenticated_source_files"] == 51
    assert identity["recursive_validation"] is False


@pytest.mark.parametrize("target", ("payload", "rehashed_payload", "binding", "rehashed_binding", "manifest"))
def test_modified_inputs_rejected(tmp_path, monkeypatch, target):
    path, manifest = binding_fixture(tmp_path, monkeypatch)
    if target in ("payload", "rehashed_payload"):
        raw = path.read_bytes()+b" "
        path.write_bytes(raw)
        if target == "rehashed_payload":
            manifest["outputs"]["compatible.json"] = dict(sha256=tool.digest(raw), bytes=len(raw))
    elif target in ("binding", "rehashed_binding"):
        (tmp_path/"bound-0.txt").write_text("modified")
        if target == "rehashed_binding":
            manifest["source_bindings"]["bound-0.txt"] = tool.digest((tmp_path/"bound-0.txt").read_bytes())
    else:
        manifest["schema"] = "wrong"
    path.with_name("manifest.json").write_bytes(tool.encoded(manifest))
    with pytest.raises(ValueError):
        tool.authenticate(tmp_path)


def test_audit_orchestration_shared_budget_no_replay(tmp_path, monkeypatch):
    cases = [dict(marker=i) for i in range(5)]
    calls = []
    monkeypatch.setattr(tool, "authenticate", lambda root: (dict(cases=cases), {}, dict(fixed=True)))
    for name in tool.OWN_SOURCES:
        p = tmp_path/name
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text("fixture")
    def cap(case, budget):
        calls.append(("cap", case["marker"]))
        budget.check(10)
        return {}
    def stage(case, parent, index, budget):
        calls.append(("stage", index))
        budget.check(20)
        return {}
    monkeypatch.setattr(tool, "cap_pullback", cap)
    monkeypatch.setattr(tool, "crossed_stage", stage)
    budget = tool.Budget(maximum_evaluations=1000, cpu_seconds=30.)
    result = tool.audit(tmp_path, budget=budget)
    assert result["execution"]["evaluated_points"] == 100
    assert result["execution"]["source_integrations"] == result["execution"]["time_steps"] == 0
    assert calls == [("cap", 0), ("cap", 1)]+[("stage", i) for i in tool.CUTS]
    with pytest.raises(RuntimeError, match="budget exhausted"):
        tool.audit(tmp_path, budget=tool.Budget(maximum_evaluations=99, cpu_seconds=30.))


def test_creation_only_and_stdout_result(tmp_path, monkeypatch):
    monkeypatch.setattr(tool, "audit", lambda root: dict(synthetic=True))
    output = tmp_path/"audit.json"
    assert tool.create(output) == dict(synthetic=True)
    assert json.loads(output.read_bytes()) == dict(synthetic=True)
    with pytest.raises(FileExistsError):
        tool.create(output)
    assert tool.create() == dict(synthetic=True)


@pytest.mark.parametrize("field,value", (("cpu_seconds", float("nan")), ("cpu_seconds", 31.),
    ("cpu_seconds", True), ("maximum_evaluations", float("nan")), ("maximum_evaluations", True),
    ("maximum_evaluations", 0), ("evaluations", -1), ("started", float("nan"))))
def test_caller_cannot_disable_budget_with_nonfinite_or_invalid_fields(field, value):
    budget = tool.Budget(maximum_evaluations=1000, cpu_seconds=30.)
    setattr(budget, field, value)
    with pytest.raises(ValueError, match="allowance"):
        tool.audit(budget=budget)
