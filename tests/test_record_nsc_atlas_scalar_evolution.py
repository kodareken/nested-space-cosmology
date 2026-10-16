"""Saved-input restoration and discrete readout authentication, with no solve."""
from copy import deepcopy
import json
from types import SimpleNamespace

import numpy as np
import pytest

import record_nsc_atlas_scalar_evolution as tool


def restored_fixture():
    from recursive_horizons.nsc_atlas_characteristic_initialization import source_geometry
    geometry = source_geometry()
    d = {key: geometry[key] for key in tool.cache_record.DOMAIN_KEYS}
    du, dv = d["Ub"]-d["Ua"], d["Vb"]-d["Va"]
    def wall(V, derivative=0):
        values = np.asarray(V)
        return np.full_like(values, du/dv) if derivative else d["Ua"]+du*(values-d["Va"])/dv
    start = d["Ua"]+du/4
    lift = SimpleNamespace(B=lambda U: np.maximum(np.asarray(U)-start, 0.)**4,
        b_at_U=lambda U: 4*np.maximum(np.asarray(U)-start, 0.)**3,
        knot_times=tool.clock_record.predecessor.FINE_TIMES)
    class Forcing:
        def __init__(self, count):
            self.count = count
        def eval(self, U, t):
            U = np.asarray(U)
            W = float(wall(d["Va"]+dv*t))
            eta = (U-d["U0"])/(W-d["U0"])
            h = max(0., t-.25)**4
            F = h*(1-eta)*(1+1/self.count)
            return dict(U=U, F=F, Fsource=float(F[0]), G=h*(1+eta**2),
                        g1=h, g1_V=4*max(0., t-.25)**3/dv)
    maxes = [{key: value for key in ("temporal_F_defect", "combined_F_defect",
        "temporal_G_defect", "combined_G_defect", "g1_defect")} for value in (.1, .01)]
    cache = {"domain": d, "diagnostics": {"maximum_absolute_defects": maxes}, "solutions": [
        {"validation_samples": [{"interpolated": {"F": [value, 0.], "G": [value, value], "g1": value}}]}
        for value in (.1, .01)]}
    return SimpleNamespace(domain=d, geometry={"wall": wall}, lift=lift,
        forcing={count: Forcing(count) for count in (17, 33)}, cache=cache,
        receipt={"files": {"cache": "a"*64, "clock": "b"*64}, "clock_time_panels": 32,
            "clock_sampled_b_target_met": True, "restored_forcing_slices": 65, "restored_boundary_jets": 81,
            "new_leading_field_solves": 0, "new_wall_solves": 0},
        geometry_defects=dict(wall_coordinate=0., wall_derivative=0., arrival_coordinate=0.),
        geometry_tolerances=dict(wall_coordinate=2e-13, wall_derivative=2e-11, arrival_coordinate=2e-13))


def synthetic(restored):
    from recursive_horizons.nsc_atlas_scalar_evolution import SCHEMA, SCOPE, reference_kernel
    d = restored.domain
    solutions = []
    for (count, nodes), ode in zip(tool.CASES, tool.ODE_CASES):
        grid = np.linspace(d["U0"], d["Ub"], nodes)
        rows, interiors = [], []
        for t in tool.TIMES:
            V = min(d["Vb"], max(d["Va"], d["Va"]+(d["Vb"]-d["Va"])*t))
            wall = min(d["Ub"], max(d["Ua"], float(restored.geometry["wall"](V))))
            wall_V = float(restored.geometry["wall"](V, 1))
            U = np.minimum(grid, wall)
            eta = (U-d["U0"])/(wall-d["U0"])
            z = max(0., t-.25)**4*(2*eta-eta*eta)*(1+1/nodes+1/count)
            supplied = restored.forcing[count].eval(U, t)
            B = restored.lift.B(U)
            y = z+supplied["g1"]+B
            integrand = reference_kernel(U, V)*y+supplied["G"]
            zV = float(supplied["F"][-1]-supplied["Fsource"]-np.sum((integrand[:-1]+integrand[1:])*np.diff(U)/2))
            b = float(restored.lift.b_at_U(wall))
            rows.append(dict(t=t, V=V, U_wall=wall, U_wall_V=wall_V, g1=supplied["g1"],
                g1_V=supplied["g1_V"], B_wall=float(B[-1]), b_wall=b, z_wall=float(z[-1]),
                y_wall=float(y[-1]), z_wall_V=zV, y_wall_V=supplied["g1_V"]+zV+wall_V*b,
                source_defect=0., source_rate_defect=0., wall_flux=float(supplied["F"][-1]), inactive_state_spread=0.))
            interiors.append(dict(U=U.tolist(), V=np.full(nodes, V).tolist(), z=z.tolist(), y=y.tolist()))
        cuts = sorted(set([*np.linspace(0, 1, count), *restored.lift.knot_times,
            *((grid[(grid > d["Ua"]) & (grid < d["Ub"]) ]-d["Ua"])/(d["Ub"]-d["Ua"]))]))
        solutions.append(dict(schema=SCHEMA, scope=SCOPE, domain=d, settings=dict(grid_nodes=nodes, **ode),
            cache_times=count, wall_samples=rows, interior_samples=interiors, rhs_calls=100,
            readout_rhs_calls=17, kernel_values=nodes*117, segments=len(cuts)-1, split_times=cuts,
            maximum_split_merge_width=0., maximum_absolute_source_defect=0.,
            maximum_absolute_source_rate_defect=0., maximum_absolute_wall_flux=max(abs(row["wall_flux"]) for row in rows)))
    return dict(schema=tool.SCHEMA, status=tool.STATUS, scope=tool.SCOPE,
        runtime={"python": "3.14.8", "numpy": "2.5.1", "scipy": "1.17.1"}, domain=d,
        source_scale=1., geometry_knots=257, clock_time_panels=32, inputs=restored.receipt,
        geometry_restoration_defects=restored.geometry_defects, solutions=solutions,
        diagnostics=tool.comparison_diagnostics(solutions), sampled_input_indicators=tool.input_indicators(restored.cache),
        refinement_is_rigorous_error_enclosure=False, execution=dict(cpu_seconds=1., cpu_seconds_max=120.,
            scalar_kernel_values=sum(row["kernel_values"] for row in solutions), scalar_kernel_values_max=30_000_000,
            rhs_calls=400, readout_rhs_calls=68, correction_solves=4, leading_field_solves=0,
            wall_solves=0, new_forcing_slices=0, new_boundary_jets=0))


@pytest.fixture
def record(monkeypatch):
    restored = restored_fixture()
    monkeypatch.setattr(tool, "restore_inputs", lambda: restored)
    return synthetic(restored)


@pytest.fixture
def isolated(tmp_path, monkeypatch, record):
    root = tmp_path/"repository"
    for name in tool.SOURCES:
        path = root/name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(("fixture: "+name).encode())
    monkeypatch.setattr(tool, "ROOT", root)
    monkeypatch.setattr(tool, "runtime_versions", lambda: dict(record["runtime"]))
    monkeypatch.setattr(tool, "compare", lambda: deepcopy(record))
    return root/"result"


def test_axes_and_discrete_readout_identities_validate_without_evolving(record):
    assert tool.validate_record(record)
    axes = record["diagnostics"]["axes"]
    assert [axis["axis"] for axis in axes] == ["cache_time", "correction_grid", "ODE_tolerance"]
    assert all(axis["maximum_absolute_differences"]["y_wall"] > 0 for axis in axes[:2])
    assert axes[2]["maximum_absolute_differences"]["y_wall"] == 0
    assert record["diagnostics"]["ODE_tolerance_refined"] is True
    assert record["execution"]["readout_rhs_calls"] == 68
    assert record["sampled_input_indicators"]["uniform_bound"] is False


@pytest.mark.parametrize("path,value", [
    (("solutions", 0, "cache_times"), 33),
    (("solutions", 2, "settings", "grid_nodes"), 129),
    (("solutions", 3, "settings", "tolerance"), 2e-6),
    (("solutions", 1, "settings", "atol"), 1e-9),
    (("solutions", 1, "wall_samples", 10, "y_wall_V"), 99.),
    (("solutions", 1, "interior_samples", 10, "y", 50), 99.),
    (("solutions", 2, "split_times", 4), .13),
    (("solutions", 0, "rhs_calls"), 6001),
    (("execution", "leading_field_solves"), 1),
    (("execution", "scalar_kernel_values"), 30_000_001),
    (("execution", "cpu_seconds"), 121.),
    (("inputs", "clock_sampled_b_target_met"), False),
    (("sampled_input_indicators", "uniform_bound"), True),
])
def test_changed_axes_readout_algebra_and_budget_are_rejected(record, path, value):
    changed = deepcopy(record)
    target = changed
    for key in path[:-1]:
        target = target[key]
    target[path[-1]] = value
    with pytest.raises(ValueError):
        tool.validate_record(changed)


def test_creation_only_check_never_propagates(isolated, monkeypatch):
    from recursive_horizons import nsc_atlas_scalar_evolution as solver
    tool.create(isolated)
    before = {path.name: path.read_bytes() for path in isolated.iterdir()}
    monkeypatch.setattr(tool, "compare", lambda: pytest.fail("unexpected comparison"))
    monkeypatch.setattr(solver, "solve_scalar_correction", lambda *a, **k: pytest.fail("unexpected propagation"))
    assert tool.check(isolated)["propagated_scalar_solves"] == 0
    with pytest.raises(FileExistsError):
        tool.create(isolated)
    assert before == {path.name: path.read_bytes() for path in isolated.iterdir()}
    (tool.ROOT/"src/recursive_horizons/nsc_atlas_scalar_evolution.py").write_bytes(b"changed solver")
    with pytest.raises(ValueError, match="source binding"):
        tool.check(isolated)


def test_rehashed_retained_state_tamper_fails(isolated):
    tool.create(isolated)
    path = isolated/"evolution.json"
    record = json.loads(path.read_bytes())
    record["solutions"][1]["interior_samples"][10]["z"][50] += 1.
    raw = tool.encoded(record)
    path.write_bytes(raw)
    manifest_path = isolated/"manifest.json"
    manifest = json.loads(manifest_path.read_bytes())
    manifest["outputs"]["evolution.json"] = {"sha256": tool.digest(raw), "bytes": len(raw)}
    manifest_path.write_bytes(tool.encoded(manifest))
    with pytest.raises(ValueError):
        tool.check(isolated)


def test_source_drift_during_compare_does_not_write(isolated, monkeypatch, record):
    def drift():
        (tool.ROOT/"src/recursive_horizons/nsc_atlas_scalar_evolution.py").write_bytes(b"drift")
        return deepcopy(record)
    monkeypatch.setattr(tool, "compare", drift)
    with pytest.raises(RuntimeError, match="inputs or runtime changed"):
        tool.create(isolated)
    assert not isolated.exists()


def test_clock_gate_compares_shared_coordinates_in_different_domain_shapes(tmp_path, monkeypatch):
    d = dict(U0=-.32, Ua=-.30, Ub=-.27, Va=-2.5, Vb=-2.25)
    cache = {"domain": d}
    clock = {"domain": {**d, "qa": 1., "qb": 1.2, "Ta": .1, "Tb": .2, "qmin": 0., "xa": .45, "xb": .55},
             "solutions": [{"wall_time_panels": 16}, {"wall_time_panels": 32}],
             "diagnostics": {"sampled_b_target_met": [False, True]}}
    monkeypatch.setattr(tool, "ROOT", tmp_path)
    monkeypatch.setattr(tool.cache_record, "check", lambda *_: {"physical_solves": 0})
    monkeypatch.setattr(tool.clock_record, "check", lambda *_: {"physical_solves": 0})
    for name, value, filename in ((tool.CACHE, cache, "cache.json"), (tool.CLOCK, clock, "clock.json")):
        directory = tmp_path/name
        directory.mkdir(parents=True)
        (directory/filename).write_bytes(tool.encoded(value))
        (directory/"manifest.json").write_bytes(b"{}")
    assert tool.input_snapshot()[2]["clock_sampled_b_target_met"] is True
    clock["diagnostics"]["sampled_b_target_met"][1] = False
    path = tmp_path/tool.CLOCK/"clock.json"
    path.write_bytes(tool.encoded(clock))
    with pytest.raises(ValueError, match="clock32"):
        tool.input_snapshot()
    clock["diagnostics"]["sampled_b_target_met"][1] = True
    clock["domain"]["Ua"] += .001
    path.write_bytes(tool.encoded(clock))
    with pytest.raises(ValueError, match="domains"):
        tool.input_snapshot()


def test_restore_uses_full_geometry_and_shared_budget_without_physical_fallback(monkeypatch):
    from recursive_horizons import nsc_atlas_prepared_response as prepared
    from recursive_horizons import nsc_atlas_forcing_cache as cache_module
    from recursive_horizons import nsc_atlas_wall_evolution as wall_module
    from recursive_horizons.nsc_atlas_characteristic_initialization import source_geometry
    source = source_geometry()
    d = {key: source[key] for key in tool.cache_record.DOMAIN_KEYS}
    geometry = prepared._geometry(257)
    direct = []
    eta = np.linspace(0, 1, 129)
    for t in np.linspace(0, 1, 33):
        V = d["Va"]+(d["Vb"]-d["Va"])*t
        wall = float(geometry["wall"](V))
        U = d["U0"]+eta*(wall-d["U0"])
        U[0], U[-1] = d["U0"], wall
        direct.append(dict(t=float(t), V=V, U_wall=wall, eta=eta.tolist(), U=U.tolist(),
            F=np.zeros(129).tolist(), G=np.zeros(129).tolist(), g1=0., g1_V=0.,
            source_p_defect=0., source_p_V_defect=0., radial_nodes=1, field_points=1, field_batches=1,
            slice_costs=dict(kernel_values=0, metric_values=0, source_nodes=0)))
    boundary = []
    for t in tool.clock_record.predecessor.DIRECT_TIMES:
        V = d["Va"]+(d["Vb"]-d["Va"])*t
        boundary.append(dict(t=t, V=V, U=float(geometry["wall"](V)), U_V=float(geometry["wall"](V, 1)),
            q_plus=float(geometry["arrival"](V)), b=0., b_U=0., b_V=0.))
    cache = dict(domain=d, source_geometry=source, direct_samples=direct)
    clock = dict(domain={**d, "qa": 1., "qb": 2.}, solutions=[{}, {"boundary_samples": boundary}])
    monkeypatch.setattr(tool, "input_snapshot", lambda: (cache, clock, {"synthetic": True}))
    monkeypatch.setattr(prepared, "solve_prepared_field", lambda *a, **k: pytest.fail("new leading solve"))
    monkeypatch.setattr(wall_module, "solve_wall_evolution", lambda *a, **k: pytest.fail("new wall solve"))
    monkeypatch.setattr(cache_module, "ScalarForcing", lambda *a, **k: pytest.fail("new forcing slice"))
    restored = tool.restore_inputs()
    assert restored.lift.budget is restored.forcing[17].budget is restored.forcing[33].budget
    V = d["Va"]+.373*(d["Vb"]-d["Va"])
    assert restored.geometry["wall"](V, 1) > 0
    assert restored.lift.B_at_wall(V) == 0.
    assert restored.forcing[33].eval(np.array([d["U0"]]), .373)["g1"] == 0.
    assert restored.shared.kernel_values == 0
    cache["direct_samples"][5]["U_wall"] += 1e-6
    with pytest.raises(ValueError, match="reference geometry"):
        tool.restore_inputs()


def test_geometry_replay_allowance_tracks_hermite_derivative_conditioning():
    from scipy.interpolate import CubicHermiteSpline
    x = np.linspace(-2.5, -2.2, 257)
    values = -.3+.1*(x-x[0])
    slopes = np.full_like(x, .1)
    original = CubicHermiteSpline(x, values, slopes)
    changed = values.copy()
    changed[128] += 8*np.finfo(float).eps
    replay = CubicHermiteSpline(x, changed, slopes)
    query = (x[127]+x[128])/2
    defect = abs(float(replay(query, 1)-original(query, 1)))
    allowance = tool.reference_derivative_roundoff(replay)
    assert 2e-13 < defect < allowance < 1e-8
    assert abs(float(replay(query)-original(query))) < 2e-13
    limits = tool._geometry_tolerances({"wall": replay})
    saved = dict(wall_coordinate=0., wall_derivative=0., arrival_coordinate=0.)
    current = {**saved, "wall_derivative": defect}
    tool._check_geometry_receipt(saved, current, limits)
    with pytest.raises(ValueError, match="geometry arithmetic"):
        tool._check_geometry_receipt(saved, {**current, "wall_derivative": 1e-7}, limits)
    with pytest.raises(ValueError, match="geometry arithmetic"):
        tool._check_geometry_receipt(saved, {**current, "wall_coordinate": 1e-10}, limits)


def test_wall_derivative_replay_allowance_is_propagated_only_to_its_readouts(monkeypatch):
    restored = restored_fixture()
    saved = synthetic(restored)
    wall = restored.geometry["wall"]
    derivative_shift = restored.geometry_tolerances["wall_derivative"]/4
    restored.geometry["wall"] = lambda V, derivative=0: wall(V, derivative)+(derivative_shift if derivative else 0.)
    restored.geometry_defects["wall_derivative"] = derivative_shift
    monkeypatch.setattr(tool, "restore_inputs", lambda: restored)
    assert tool.validate_record(saved)
    damaged = deepcopy(saved)
    damaged["solutions"][0]["wall_samples"][10]["U_wall_V"] += 1e-7
    with pytest.raises(ValueError, match="wall-rate"):
        tool.validate_record(damaged)
