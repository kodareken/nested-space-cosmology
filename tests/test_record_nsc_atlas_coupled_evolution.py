"""Synthetic five-case evidence/authentication controls; no physical solve."""
import copy
import json
import time

import numpy as np
import pytest
from scipy.special import lambertw

import record_nsc_atlas_coupled_evolution as tool
from recursive_horizons.nsc_atlas_characteristic_initialization import source_jet


GEOMETRY = dict(U0=-.32130105902052203, Ua=-.2992702973503126,
    Ub=-.2747468030926288, Va=-2.565695357068247, Vb=-2.262011613052105,
    order=64, support_V=[-2.489774421064211, -2.3379325490561405],
    support_x=[.42269179059454887, .4631925883576189], x_corner=.4, x_last=.4815042082576364)


def test_recomputed_numerics_allow_roundoff_with_exact_structure_and_counts():
    tool._equal({"metric": np.array([1., 0.]), "count": 2, "flag": False},
                {"metric": [1.+1e-13, 1e-14], "count": 2, "flag": False}, "recomputed metric")
    for actual, expected in (({"count": 2.}, {"count": 2}),
                             ({"flag": 0}, {"flag": False}),
                             ([0.], [0., 0.]),
                             ([1.+1e-8], [1.]),
                             ([float("nan")], [0.]),
                             ([float("inf")], [0.])):
        with pytest.raises(ValueError, match="recomputed"):
            tool._equal(actual, expected, "recomputed")


def test_snapshot_replay_conditions_residuals_but_not_copied_coordinates():
    from recursive_horizons.nsc_atlas_validation import source_snapshot_roundoff_tolerances
    x = np.linspace(.4, GEOMETRY["x_last"], 129)
    snapshot = object.__new__(tool.SourceSnapshot)
    snapshot._geometry, snapshot._alpha = GEOMETRY, 1e-4
    snapshot._bounds = dict(r_min=.2, v_max=2., F_min=.5)
    h, hx = 30*(x-.4)**3, 90*(x-.4)**2
    preliminary = snapshot._frame(x, h, hx, np.zeros_like(x))
    row = snapshot._frame(x, h, hx, preliminary["source_ode_h_xx"])
    knots = dict(x=x, h=h, h_x=hx)
    tolerances = source_snapshot_roundoff_tolerances(row, GEOMETRY, 1e-4, source_knots=knots)
    changed = copy.deepcopy(row)
    index = int(np.argmax(tolerances["source_h_xx_defect"]))
    changed["source_h_xx_defect"][index] += tolerances["source_h_xx_defect"][index]/2
    tool._source_snapshot_equal(changed, row, GEOMETRY, 1e-4, source_knots=knots)
    changed["source_h_xx_defect"][index] += 1e-4
    with pytest.raises(ValueError, match="source_h_xx_defect.*max_abs_difference=.*allowance="):
        tool._source_snapshot_equal(changed, row, GEOMETRY, 1e-4, source_knots=knots)
    changed = copy.deepcopy(row); changed["x"][2] = np.nextafter(x[2], np.inf)
    with pytest.raises(ValueError, match="exact coordinates"):
        tool._source_snapshot_equal(changed, row, GEOMETRY, 1e-4, source_knots=knots)
    changed = copy.deepcopy(row); changed["c"][0] = 1e-16
    with pytest.raises(ValueError, match="mismatch: c;"):
        tool._source_snapshot_equal(changed, row, GEOMETRY, 1e-4, source_knots=knots)
    changed = copy.deepcopy(row); changed["C"] = list(changed["C"]); changed["C"][0] = True
    with pytest.raises(ValueError, match="shape/type/nonfinite"):
        tool._source_snapshot_equal(changed, row, GEOMETRY, 1e-4, source_knots=knots)


@pytest.fixture(autouse=True)
def bounded_control():
    started = time.process_time()
    yield
    assert time.process_time()-started < 30.


def synthetic_gauge(x):
    # Algebraically consistent pointwise reference rows, not a quadrature
    # claim. The synthetic radius sequence/clock integrals are test data.
    R = np.linspace(.45, .55, len(x))
    bm, bp = .2/R**2+.75*R, .2/R**2-1.25*R
    v = np.sqrt(bm*bm-1+R*R)
    dv = (25*R/8-3/(10*R**2)-4/(25*R**5))/(2*v)
    dbp = -.4/R**3-1.25
    lam = -x/(.8*(.8-x))
    a = -1/(v*(v-bp))
    Rx = 1.6*lam/a
    Rxx = Rx*((-1/(.8-x)**2)/lam-(-dv/v-(dv-dbp)/(v-bp))*Rx)
    V = (1-x/.8)*np.exp(x/.8)/GEOMETRY["U0"]
    U = (1-R/.8)*np.exp(R/.8)/V
    U[0], U[-1] = GEOMETRY["Ua"], GEOMETRY["Ub"]
    CJ = -3.2/GEOMETRY["U0"]
    Z = -R*Rx/(.8*(.8-R))-lam
    Zx = -Rx*Rx/(.8-R)**2-R*Rxx/(.8*(.8-R))+1/(.8-x)**2
    W, Wx, Wxx = CJ*(U-GEOMETRY["U0"]), CJ*U*Z, CJ*U*(Z*Z+Zx)
    C = x/R*np.exp((x-R)/.8)
    sx = Rx/v
    rows = dict(x=x, r=R, U=U, V=V, W=W, W_x=Wx, W_xx=Wxx,
        r_x_total=Rx, r_xx_total=Rxx, eta_plus=-.5*np.log(Wx), eta_plus_x=-.5*Wxx/Wx,
        L=np.log(x/R)+(x-R)/.8, C=C, s_x=sx, T_x=bm*sx/(1-R*R),
        T=3*(x-x[0]), s=2*(x-x[0]), inverse_residual=np.zeros(len(x)),
        normalization_defect=C*Wx/(sx*sx)-1)
    rU, rx = C*(.8-x)/(4*x), (1-W/3.2)*C
    rows.update(r_Ubar=rU, r_x=rx, log_C_Ubar=-(1/R+1/.8)*rU,
        log_C_x=1/x+1/.8-(1/R+1/.8)*rx, L_x_total=1/x+1/.8-(1/R+1/.8)*Rx,
        T_s=bm/(1-R*R), x_s=1/sx, Ubar_s=Wx/sx, beta_minus=bm, beta_plus=bp,
        r_s=v, c=np.zeros(len(x)), c_Ubar=np.zeros(len(x)), c_x=np.zeros(len(x)),
        inverse_iterations=np.zeros(len(x), dtype=int))
    summary = dict(query_points=513, cache_hits=0, inverse_calls=511, quadrature_order=64,
        source_integrations=0, regional_marches=0, comparison_is_rigorous_enclosure=False,
        quadrature_nodes=64*512, evaluations=3+513+64*512, maximum_inverse_residual=0.,
        maximum_inverse_iterations=0, endpoint_integral_defect=0.,
        schema="nsc-atlas-reference-wall-gauge-v1", domain=dict(x_a=.4, x_b=x[-1],
            r_a=.45, r_b=.55, W_a=W[0], W_b=W[-1]))
    return dict(rows=rows, summary=summary)


def synthetic_source(x, alpha, budget):
    xexact = x[::2]
    t = xexact-.4
    h, hx = 20*t**3, 60*t*t
    p, px, _ = source_jet(xexact, GEOMETRY)
    quiet = (xexact <= GEOMETRY["support_x"][0]) | (xexact >= GEOMETRY["support_x"][1])
    g, gx = np.where(quiet, 0., xexact*p), np.where(quiet, 0., p+xexact*px)
    rJ, rJx = xexact+alpha*h, 1+alpha*hx
    F = 1-alpha*g*g/(6*rJ*rJ)
    r = rJ*np.sqrt(F)
    rx = (rJ*rJx-alpha*g*gx/6)/r
    pj, pjx = g/rJ, (gx-g*rJx/rJ)/rJ
    a = (.8/xexact-1)/4
    exact = dict(x=xexact, h=h, h_x=hx, rJ=rJ, rJ_x=rJx, r=r, C=F,
        log_C=np.log1p(-alpha*g*g/(6*rJ*rJ)), c=g, c_x=gx, r_x=rx,
        log_C_x=-alpha*pj*pjx/(3*F), r_Ubar=a, c_Ubar=-t*t,
        log_C_Ubar=-(1/xexact+1/.8)*a)
    check = tool.SourceSnapshot(GEOMETRY, alpha=alpha, x_nodes=xexact, source_rows=exact, budget=budget)
    base = check._frame(xexact, h, hx, np.zeros(257))
    frame = check._frame(xexact, h, hx, base["source_ode_h_xx"])
    for key in ("r_xx", "source_constraint", "source_constraint_jordan"):
        exact[key] = frame[key]
    exact["m"] = r*(1+4*a*rx/F)/2
    snapshot = tool.SourceSnapshot(GEOMETRY, alpha=alpha, x_nodes=x[::4],
        source_rows=tool.subset(exact, slice(None, None, 2)), budget=budget)
    cached = snapshot.sample(x)
    W0 = -3.2/GEOMETRY["U0"]*(GEOMETRY["Ua"]-GEOMETRY["U0"])
    ub = np.linspace(0., W0, 33)
    vac = tool.previous.schwarzschild(ub, [.4])
    br, bC, bL = (vac[key][:, 0] for key in ("r", "C", "L"))
    bottom = {key: np.zeros(33) for key in tool.previous.BOTTOM_KEYS}
    bottom.update(Ubar=ub, U=GEOMETRY["U0"]+ub/(-3.2/GEOMETRY["U0"]), x=np.full(33, .4),
        rJ=br, rJ_x=(1-ub/3.2)*bC, r=br, C=bC, log_C=bL,
        r_Ubar=bC/4, r_x=(1-ub/3.2)*bC,
        log_C_Ubar=-(1/br+1/.8)*bC/4,
        log_C_x=1/.4+1/.8-(1/br+1/.8)*(1-ub/3.2)*bC,
        m=np.full(33, .4))
    return dict(alpha=alpha, domain=snapshot.domain, exact_rows=exact, bottom=bottom,
        snapshot_rows=cached, snapshot_summary=snapshot.summary(), source_nfev=100,
        source_max_step=(x[-1]-x[0])/32, source_execution=dict(evaluations=393, cpu_seconds=.001))


@pytest.fixture(scope="module")
def record():
    budget = tool.Budget(maximum_evaluations=8000000, cpu_seconds=30.)
    x = np.linspace(.4, GEOMETRY["x_last"], 513)
    bounds = tool.source_amplitude_bound(GEOMETRY)
    initial = {label: synthetic_source(x, alpha, budget)
        for label, alpha in (("zero", 0.), ("finite", bounds["test_alpha"]))}
    gauge = synthetic_gauge(x)
    result = dict(schema=tool.SCHEMA, status=tool.STATUS, scope=tool.SCOPE, settings=tool.SETTINGS,
        runtime=tool.runtime_versions(), geometry=GEOMETRY, source_bounds=bounds, initial_data=initial,
        gauge=gauge, cases=[], refinement_is_rigorous_error_enclosure=False)
    construction = 963+gauge["summary"]["evaluations"]+sum(
        data["source_execution"]["evaluations"]+data["snapshot_summary"]["evaluations"] for data in initial.values())
    end = construction
    for label, n, steps in tool.CASES:
        data, samples = initial[label], []
        eta = np.linspace(0., 1., n+1)
        source_values = dict(r=data["snapshot_rows"]["r"], L=data["snapshot_rows"]["log_C"], c=data["snapshot_rows"]["c"])
        for i in range(0, 513, 4):
            U = eta*gauge["rows"]["W"][i]
            product = (1-U/3.2)*(1-x[i]/.8)*np.exp(x[i]/.8)
            r = .8*(1+lambertw(-product/np.e).real)
            L = np.log(x[i]/r)+(x[i]-r)/.8
            fields = dict(r=r, L=L, c=np.zeros(n+1))
            for key in fields:
                fields[key] += (1-eta)*(source_values[key][i]-fields[key][0])
            # Controlled synthetic resolution effects keep comparisons nonzero.
            fields["r"] += (x[i]-.4)*eta/n*(.01+(data["alpha"]*.2))
            state = dict(increment={key: fields[key][1:-1]-source_values[key][i] for key in fields},
                wall={**{key: float(fields[key][-1]) for key in fields},
                      "T": float(gauge["rows"]["T"][i]), "s": float(gauge["rows"]["s"][i])},
                cap={key: np.zeros(n) for key in ("f", "Pi", "Q")})
            samples.append(dict(x=float(x[i]), state=state))
        samples[0]["state"] = tool.quiet_initial_state(tool.subset(data["bottom"], slice(None, None, 32//n)),
            tool.subset(data["snapshot_rows"], 0), cap_cells=n)
        stage = tool.coupled_stage(samples[-1]["state"], eta=eta, y_faces=eta,
            source=tool.subset(data["snapshot_rows"], -1), gauge=tool.subset(gauge["rows"], -1),
            alpha=data["alpha"], tube=tool.TUBE, budget=budget, diagnostics=True)
        charge = (2*steps+1)*(32*n+7)+3*n
        result["cases"].append(dict(amplitude=label, cells=n, steps=steps, alpha=data["alpha"],
            snapshots=samples, final=tool.final_summary(stage), execution=dict(budget_start=end,
                budget_end=end+charge, stages=2*steps+1, source_cache_reads=2*steps+2,
                gauge_cache_reads=2*steps+1, max_courant=.1, retries=0)))
        end += charge
    start = budget.evaluations
    result["diagnostics"] = tool.diagnostics(result, budget)
    charges = budget.evaluations-start
    result["execution"] = dict(cpu_seconds=.02, evaluations=end+charges,
        construction_evaluations=construction, diagnostic_evaluations=charges,
        source_integrations=2, source_segments=6, gauge_constructions=1, gauge_queries=513,
        snapshot_constructions=2, snapshot_queries=1026, coupled_evolutions=5,
        time_steps=sum(k for _, _, k in tool.CASES), stage_evaluations=sum(2*k+1 for _, _, k in tool.CASES),
        historical_record_replays=0, old_rectangle_solves=0)
    result = tool.plain(result)
    auth = tool.validate_record(result, require_creation=False)
    result["creation_execution"] = dict(calculation_evaluations=result["execution"]["evaluations"],
        validation_evaluations=auth["authentication_evaluations"],
        total_evaluations=result["execution"]["evaluations"]+auth["authentication_evaluations"], cpu_seconds=.2)
    return result


def test_synthetic_record_reconstructs_without_physical_calls(record, monkeypatch):
    def forbidden(*a, **kw):
        pytest.fail("physical solve/quadrature in authentication")
    monkeypatch.setattr(tool, "ExactSource", forbidden)
    monkeypatch.setattr(tool, "ReferenceWallGauge", forbidden)
    monkeypatch.setattr(tool, "evolve_coupled", forbidden)
    result = tool.validate_record(record)
    assert result["source_integrations"] == result["time_evolutions"] == result["reference_quadrature"] == 0
    assert result["final_cap_radial_reconstructions"] == 5


@pytest.mark.parametrize("change", ["settings", "source", "curvature", "gauge", "clock", "initial", "snapshot", "final", "counts", "creation", "comparison", "overlap", "finite"])
def test_tampered_record_rejected(record, change):
    value = copy.deepcopy(record)
    if change == "settings": value["settings"]["source_snapshot_knots"] = 257
    elif change == "source": value["initial_data"]["finite"]["exact_rows"]["c"][25] += .01
    elif change == "curvature": value["initial_data"]["finite"]["snapshot_rows"]["r_xx"][99] += .1
    elif change == "gauge": value["gauge"]["rows"]["W_x"][99] *= 1.01
    elif change == "clock": value["cases"][0]["snapshots"][2]["state"]["wall"]["T"] = 0.
    elif change == "initial": value["cases"][0]["snapshots"][0]["state"]["cap"]["f"][0] = .01
    elif change == "snapshot": value["cases"][0]["snapshots"][50]["state"]["increment"]["r"][2] += .001
    elif change == "final": value["cases"][2]["final"]["cap_faces"]["A_faces"][1] += .001
    elif change == "counts": value["cases"][1]["execution"]["budget_start"] += 1
    elif change == "creation": value["creation_execution"]["total_evaluations"] -= 1
    elif change == "comparison": value["diagnostics"]["paired_spatial"]["wall_arrays"]["r"][20] += .001
    elif change == "overlap": value["diagnostics"]["original_rectangle_overlap"]["cases"][0]["difference_linear"]["r"][1][1] += .001
    elif change == "finite": value["cases"][0]["snapshots"][20]["state"]["cap"]["Pi"][2] = float("nan")
    with pytest.raises(ValueError):
        tool.validate_record(value)


def test_weighted_cap_restriction_is_not_plain_Pi_average():
    table = dict(exterior={"r": np.zeros((2, 5))}, wall={"r": np.ones(2)},
        cap={key: np.broadcast_to(np.arange(4), (2, 4)) for key in ("f", "Pi", "Q")})
    restricted = tool._restricted(table)
    np.testing.assert_allclose(restricted["cap"]["f"], [[.5, 2.5]]*2)
    assert restricted["cap"]["Pi"][0, 0] == pytest.approx(7/8)
    assert restricted["cap"]["Q"][0, 0] == pytest.approx(3/4)


def test_cache_has_exact_keys_and_no_fallback():
    rows = tool.StoredRows(dict(x=[.4, .5], r=[1., 2.]))
    assert rows(.4)["r"] == 1.
    with pytest.raises(ValueError, match="absent"):
        rows(np.nextafter(.4, 1.))
    assert rows.calls == 1


def fake_calculation(record):
    def calculate(budget):
        budget.check(record["execution"]["evaluations"])
        result = copy.deepcopy(record)
        result.pop("creation_execution", None)
        result["execution"]["cpu_seconds"] = 0.
        return result
    return calculate


def test_creation_only_roundtrip_shared_budget_and_hash_checks(record, tmp_path, monkeypatch):
    monkeypatch.setattr(tool, "calculate", fake_calculation(record))
    directory = tmp_path/"evidence"
    created = tool.create(directory)
    assert created["creation_execution"]["total_evaluations"] == record["creation_execution"]["total_evaluations"]
    assert tool.check(directory)["time_evolutions"] == 0
    with pytest.raises(FileExistsError):
        tool.create(directory)
    raw = (directory/"coupled.json").read_bytes()
    (directory/"coupled.json").write_bytes(raw+b" ")
    with pytest.raises(ValueError, match="hash"):
        tool.check(directory)


def test_source_drift_fails_before_writing_result(record, tmp_path, monkeypatch):
    calls = [0]
    def bindings():
        calls[0] += 1
        return {"source": str(calls[0])}
    monkeypatch.setattr(tool, "source_bindings", bindings)
    monkeypatch.setattr(tool, "calculate", fake_calculation(record))
    monkeypatch.setattr(tool, "ROOT", tmp_path)
    with pytest.raises(RuntimeError, match="changed") as error:
        tool.create(tmp_path/"no-result")
    assert not (tmp_path/"no-result").exists()
    receipt = json.loads(open(error.value.failure_receipt).read())
    assert receipt["status"] == "FAILED_NOT_COMPLETED"


def test_calculation_and_authentication_share_one_ceiling(record, tmp_path, monkeypatch):
    seen = []
    def calculate(budget):
        seen.append(budget)
        budget.check(7999999)
        return copy.deepcopy(record)
    def validate(value, budget=None, **kwargs):
        assert budget is seen[0]
        budget.check(2)
    monkeypatch.setattr(tool, "calculate", calculate)
    monkeypatch.setattr(tool, "validate_record", validate)
    monkeypatch.setattr(tool, "ROOT", tmp_path)
    monkeypatch.setattr(tool, "source_bindings", lambda: {"fixed": "hash"})
    with pytest.raises(RuntimeError, match="budget exhausted"):
        tool.create(tmp_path/"no-result")
    assert seen[0].evaluations == 8000001
    assert not (tmp_path/"no-result").exists()


def test_global_setup_reaches_bounded_source_child_without_physical_work(monkeypatch):
    budget = tool.Budget(maximum_evaluations=8000000, cpu_seconds=180.)
    def geometry(order, owner):
        assert owner is budget and order == 64
        owner.check(960)
        return copy.deepcopy(GEOMETRY)
    def source(*args, **kwargs):
        child = kwargs["budget"]
        assert child.parent is budget and child.maximum_evaluations == 20000 and child.cpu_seconds == 30.
        raise RuntimeError("intentional stop before source construction")
    monkeypatch.setattr(tool, "source_geometry", geometry)
    monkeypatch.setattr(tool, "ExactSource", source)
    with pytest.raises(RuntimeError, match="intentional stop"):
        tool.calculate(budget)
    assert budget.evaluations == 963
