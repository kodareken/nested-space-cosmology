"""Synthetic record/authentication controls: no source or regional solve."""
import copy
import json
import time

import numpy as np
import pytest

import record_nsc_atlas_nonlinear_rectangle as tool


@pytest.fixture(autouse=True)
def finite_test_cpu():
    started = time.process_time()
    yield
    assert time.process_time()-started < 30


def synthetic():
    # The original retained source labels; no reference quadrature is run.
    geometry = dict(U0=-.32130105902052203, Ua=-.2992702973503126,
        Ub=-.2747468030926288, Va=-2.565695357068247, Vb=-2.262011613052105,
        order=64, support_V=[-2.489774421064211, -2.3379325490561405],
        support_x=[.42269179059454887, .4631925883576189], x_corner=.4, x_last=.4815042082576364)
    bounds = tool.source_amplitude_bound(geometry)
    CJ0 = -3.2/geometry["U0"]
    Ubar = np.linspace(0., CJ0*(geometry["Ua"]-geometry["U0"]), 33)
    x = np.linspace(.4, geometry["x_last"], 33)
    domain = {key: geometry[key] for key in ("U0", "Ua", "Ub", "Va", "Vb")}
    domain.update(x_a=.4, x_b=x[-1], Ubar_a=0., Ubar_b=Ubar[-1], CJ0=CJ0)
    p, px, pxx = tool.source_jet(x, geometry)
    quiet = (x <= geometry["support_x"][0]) | (x >= geometry["support_x"][1])
    g, gx, gxx = (np.where(quiet, 0., value) for value in (x*p, p+x*px, 2*px+x*pxx))
    initial = {}
    for label, alpha in (("zero", 0.), ("finite", bounds["alpha_star"]/2)):
        # Pointwise frame-consistent synthetic boundary rows. These are not
        # claimed to be an integrated solution; calculate() is replaced below.
        F = 1-alpha*g*g/(6*x*x)
        r = x*np.sqrt(F)
        rx = (x-alpha*g*gx/6)/r
        a = (.8/x-1)/4
        Q = (g*gxx-2*gx*gx)/(6*x)+g*gx/(3*x*x)
        rxx = (1+x*alpha*Q-alpha*(gx*gx+g*gxx)/6-rx*rx)/r
        pJ, pJx = g/x, (gx-g/x)/x
        Lx = -alpha*pJ*pJx/(3*F)
        source = {key: np.zeros(33) for key in tool.SOURCE_KEYS}
        source.update(x=x, rJ=x, rJ_x=np.ones(33), r=r, C=F,
            log_C=np.log1p(-alpha*g*g/(6*x*x)), c=g, r_Ubar=a, r_x=rx,
            c_Ubar=-.1*(x-.4), c_x=gx, log_C_Ubar=-(1/r+1/.8)*a, log_C_x=Lx,
            m=r*(1+4*a*rx/F)/2, r_xx=rxx,
            source_constraint=rxx-Lx*rx+alpha*r*F*((gx-g/r*rx)/r)**2/2)
        exact = tool.schwarzschild(Ubar, np.array([.4]))
        br, bL, bC = (exact[key][:, 0] for key in ("r", "L", "C"))
        bU, bx = bC*(.8-.4)/(4*.4), (1-Ubar/3.2)*bC
        bottom = {key: np.zeros(33) for key in tool.BOTTOM_KEYS}
        bottom.update(Ubar=Ubar, U=geometry["U0"]+Ubar/CJ0, x=np.full(33, .4),
            rJ=br, rJ_x=bx, r=br, C=bC, log_C=bL, r_Ubar=bU, r_x=bx,
            log_C_Ubar=-(1/br+1/.8)*bU,
            log_C_x=1/.4+1/.8-(1/br+1/.8)*bx,
            m=br*(1+4*bU*bx/bC)/2)
        initial[label] = dict(alpha=alpha, domain=domain, source=source, bottom=bottom,
            max_step=bounds["length"]/32, nfev=100,
            execution=dict(evaluations=169+(3 if label == "zero" else 0), cpu_seconds=.01))
    grids = []
    charge = 960+sum(item["execution"]["evaluations"] for item in initial.values())
    for label, cells in tool.CASES:
        step, data = 32//cells, initial[label]
        ua, xa = Ubar[::step], x[::step]
        source, bottom = tool._boundary_states(data["source"], step), tool._boundary_states(data["bottom"], step)
        exact = tool.schwarzschild(ua, xa)
        states = {"r": exact["r"].copy(), "L": exact["L"].copy(),
                  "c": np.broadcast_to(source["c"], exact["r"].shape).copy()}
        for key in ("r", "L"):
            states[key] += (1-ua[:, None]/ua[-1])*(source[key]-states[key][0])[None, :]
        for key in states:
            states[key][:, 0] = bottom[key]
            states[key][0, :] = source[key]
        cell_data = {key: np.full((cells, cells), value, dtype=int) for key, value in
                     (("updates", 1), ("backtrack_trials", 1), ("rhs_evaluations", 3), ("jacobian_evaluations", 2))}
        cell_data.update({key: np.full((cells, cells), value) for key, value in
                          (("residual_norm", .1), ("correction_norm", .1), ("condition_number", 1.1))})
        cell_data["seed_kind"] = np.full((cells, cells), "affine_predictor")
        summary = dict(cell_count=cells*cells, budget_start=charge, budget_end=charge+5*cells*cells,
            seed_counts={"affine_predictor": cells*cells},
            backtrack_trials=cells*cells, rhs_evaluations=3*cells*cells, jacobian_evaluations=2*cells*cells,
            max_updates=1, max_residual_norm=.1, max_correction_norm=.1, max_condition_number=1.1)
        charge = summary["budget_end"]
        grids.append(dict(Ubar=ua, x=xa, states=states, cells=cell_data, summary=summary,
            source=source, bottom=bottom, amplitude=label, cell_count_per_axis=cells,
            settings=dict(alpha=data["alpha"], tube=tool.TUBE, cell_options=tool.CELL_OPTIONS,
                axis_order=["Ubar", "x"], corner_tolerance=5e-13, corner_choice="source",
                corner_difference={key: float(source[key][0]-bottom[key][0]) for key in source})))
    record = tool.plain(dict(schema=tool.SCHEMA, status=tool.STATUS, scope=tool.SCOPE,
        settings=tool.SETTINGS, runtime=tool.runtime_versions(), geometry=geometry,
        source_bounds=bounds, initial_data=initial, grids=grids,
        refinement_is_rigorous_error_enclosure=False,
        execution=dict(cpu_seconds=.1, evaluations=charge, geometry_evaluations=960,
            source_integrations=2, source_segments=6, rectangle_marches=4, successful_cells=2560,
            leading_field_solves=0, wall_solves=0, cap_solves=0)))
    record["diagnostics"] = tool.diagnostics(record)
    return record


@pytest.fixture(scope="module")
def base():
    return synthetic()


@pytest.fixture
def isolated(tmp_path, monkeypatch, base):
    root = tmp_path/"repository"
    for name in tool.SOURCES:
        path = root/name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(("fixture "+name).encode())
    monkeypatch.setattr(tool, "ROOT", root)
    monkeypatch.setattr(tool, "calculate", lambda: copy.deepcopy(base))
    return root/"result"


def test_four_grid_authentication_accepts_recorded_roundoff_corner(base):
    assert base["grids"][0]["settings"]["corner_difference"]["r"] != 0
    assert tool.validate_record(base)
    assert base["diagnostics"]["grids"][0]["boundary_summary"]["source"]["postpulse"]["c_Ubar"]["count"] > 0
    assert len(base["diagnostics"]["quiet_field_controls"]) == 4


def test_creation_only_and_check_never_solve(isolated, monkeypatch):
    tool.create(isolated)
    before = {p.name: p.read_bytes() for p in isolated.iterdir()}
    def forbidden(*args, **kwargs):
        pytest.fail("physical solve during authentication")
    monkeypatch.setattr(tool, "calculate", forbidden)
    monkeypatch.setattr(tool, "ExactSource", forbidden)
    monkeypatch.setattr(tool, "source_geometry", forbidden)
    import recursive_horizons.nsc_atlas_nonlinear_rectangle as rectangle
    monkeypatch.setattr(rectangle, "solve_rectangle", forbidden)
    monkeypatch.setattr(rectangle, "solve_characteristic_cell", forbidden)
    assert tool.check(isolated)["cell_solves"] == 0
    with pytest.raises(FileExistsError, match="creation-only"):
        tool.create(isolated)
    assert before == {p.name: p.read_bytes() for p in isolated.iterdir()}


@pytest.mark.parametrize("change", ["source", "bytes", "inventory", "rehash-state", "rehash-summary",
                                  "rehash-settings", "rehash-boundary", "rehash-counts"])
def test_tampering_rejected_after_rehash_too(isolated, change):
    tool.create(isolated)
    payload = isolated/"rectangle.json"
    if change == "source":
        (tool.ROOT/tool.SOURCES[0]).write_bytes(b"changed")
    elif change == "bytes":
        payload.write_bytes(payload.read_bytes()+b"\n")
    elif change == "inventory":
        (isolated/"extra").write_bytes(b"unexpected")
    else:
        record = json.loads(payload.read_bytes())
        if change == "rehash-state":
            record["grids"][2]["states"]["c"][6][8] += .001
        elif change == "rehash-summary":
            record["grids"][1]["summary"]["budget_end"] -= 1
        elif change == "rehash-settings":
            record["grids"][1]["settings"]["cell_options"]["rtol"] = 1e-10
        elif change == "rehash-boundary":
            record["initial_data"]["finite"]["source"]["log_C_x"][16] += .01
        else:
            record["execution"]["evaluations"] = 60001
        raw = tool.encoded(record)
        payload.write_bytes(raw)
        manifest_path = isolated/"manifest.json"
        manifest = json.loads(manifest_path.read_bytes())
        manifest["outputs"]["rectangle.json"] = dict(sha256=tool.digest(raw), bytes=len(raw))
        manifest_path.write_bytes(tool.encoded(manifest))
    with pytest.raises(ValueError):
        tool.check(isolated)


def test_failed_attempt_preserves_diagnosis_without_completed_output(isolated, monkeypatch):
    def failed():
        error = RuntimeError("bounded cell stopped")
        error.rectangle_failure = dict(phase=dict(name="march", amplitude="finite", cells=32),
            diagnostics=dict(indices=[4, 7], reason="newton_update_limit", condition_number=float("inf")),
            execution=dict(evaluations=4500, cpu_seconds=2.))
        raise error
    monkeypatch.setattr(tool, "calculate", failed)
    with pytest.raises(RuntimeError) as stopped:
        tool.create(isolated)
    assert not isolated.exists()
    receipt = json.loads(tool.Path(stopped.value.failure_receipt).read_bytes())
    assert receipt["status"] == "FAILED_NOT_COMPLETED"
    assert receipt["failure"]["diagnostics"]["indices"] == [4, 7]
    assert receipt["failure"]["diagnostics"]["condition_number"] == {"nonfinite": "inf"}
    assert receipt["source_bindings_before"] == receipt["source_bindings_after"]


def test_source_change_during_calculation_leaves_only_failure_receipt(isolated, monkeypatch, base):
    def changed():
        (tool.ROOT/tool.SOURCES[0]).write_bytes(b"changed mid-calculation")
        return copy.deepcopy(base)
    monkeypatch.setattr(tool, "calculate", changed)
    with pytest.raises(RuntimeError, match="inputs or runtime changed"):
        tool.create(isolated)
    assert not isolated.exists()
    receipt = json.loads(next((tool.ROOT/".build/nonlinear-rectangle-failures").glob("*.json")).read_bytes())
    assert receipt["source_bindings_before"] != receipt["source_bindings_after"]


def test_child_budget_debits_parent_and_both_limits_remain_live(monkeypatch):
    now = [100.]
    monkeypatch.setattr(tool.time, "process_time", lambda: now[0])
    parent = tool.Budget(maximum_evaluations=60000, cpu_seconds=120., started=100.)
    child = tool.SourceBudget(parent)
    child.started = now[0]
    child.check(13)
    assert child.evaluations == parent.evaluations == 13
    now[0] = 131.
    with pytest.raises(RuntimeError, match="child budget"):
        child.check()
    now[0] = 132.
    second = tool.SourceBudget(parent)
    second.started = now[0]
    second.check(19999)
    with pytest.raises(RuntimeError, match="child budget"):
        second.check(2)
    assert parent.evaluations == 20014  # The failing charge is not erased.
    second.started = now[0]
    second.evaluations = 0
    parent.evaluations = 60000
    with pytest.raises(RuntimeError, match="budget exhausted"):
        second.check(1)


def test_publication_hook_requires_named_record_and_no_fallback(monkeypatch):
    import check_publication as publication
    visited = []
    monkeypatch.setattr(tool, "check", lambda path: visited.append(path))
    errors = []
    publication.check_nonlinear_rectangle_evidence(errors)
    assert visited == [publication.ROOT/tool.RELATIVE_OUTPUT] and errors == []
    def missing(path):
        raise FileNotFoundError("missing record")
    monkeypatch.setattr(tool, "check", missing)
    publication.check_nonlinear_rectangle_evidence(errors)
    assert len(errors) == 1 and "nonlinear exterior rectangle" in errors[0]
