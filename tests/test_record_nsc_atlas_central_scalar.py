"""Saved reconstruction and provenance controls, without incoming integration."""
from copy import deepcopy
import json
import math

import pytest

import record_nsc_atlas_central_scalar as tool


def test_arrival_tangent_replay_keeps_retained_value_and_rejects_corruption():
    import numpy as np
    from scipy.interpolate import CubicHermiteSpline
    x = np.linspace(-2.5, -2.2, 257)
    q = .5+2*(x-x[0])
    original = CubicHermiteSpline(x, q, np.full_like(x, 2.))
    changed = q.copy()
    changed[128] += 8*np.finfo(float).eps
    replay = CubicHermiteSpline(x, changed, np.full_like(x, 2.))
    query = (x[127]+x[128])/2
    retained, current = float(original(query, 1)), float(replay(query, 1))
    allowance = tool.exterior_record.reference_derivative_roundoff(replay)
    assert retained != current
    assert tool._check_arrival_derivative(retained, current, allowance) == retained
    for invalid in (retained+1e-7, float("nan"), -retained):
        with pytest.raises(ValueError, match="arrival"):
            tool._check_arrival_derivative(invalid, current, allowance)


def synthetic(inputs):
    solutions = []
    for order in tool.RADIAL_ORDERS:
        samples = []
        for i, t in enumerate(tool.TIMES):
            c, w = inputs["caps"][i], inputs["walls"][i]
            q, x, d = c["q0"], c["x"], inputs["domain"]
            qV = inputs["arrival_spline_derivatives"][i]
            profile = dict(f=w["w"], f_prime=w["j"]/(qV*w["Vdot"]), f_second=0. if i == 0 else .1)
            radii, optical = tool._quadrature_shape(q, x, inputs["profile_edges"], order)
            integral = dict(q=q, x=x, u_wall=q-2*math.atanh(x), reflected_profile_defect=0.,
                central_profile_defect=0., analytic_regular_center_K1=0., first_arrival_margin=q-d["qa"],
                last_arrival_margin=d["qb"]-q, no_return_margin=d["qa"]-q+2*math.atanh(x),
                ray_nodes=len(radii), optical_radial_nodes=optical, center_samples=[])
            for kind, value in (("forcing", t*(1+1/order)), ("optical_time", t/2*(1+1/order)),
                                 ("lambda_plus", c["ray_integral"]+.001*t/order)):
                value = 0. if i == 0 else value
                integral.update({kind+"_integral": value, kind+"_positive": max(0., value),
                    kind+"_negative": min(0., value), kind+"_sample_min": value/x,
                    kind+"_sample_max": value/x})
            for r in radii[:2]:
                K = float(r*r*(0. if i == 0 else t))
                integral["center_samples"].append(dict(r=float(r), T=q-math.atanh(r),
                    W1=4*K/(1-r*r), K1=K, K1_over_r_squared=K/(r*r)))
            controls = tool.reconstruct_rows(inputs, i, profile, integral, qV)
            samples.append(dict(t=t, profile=profile, arrival_spline_derivative=qV,
                leading_compatibility=tool.leading_compatibility(inputs, i, profile, qV), integrals=integral,
                reconstruction=controls[3], exterior_controls=controls,
                saved_cap_arrival_difference=controls[3]["Delta_C_fixed_s"]-c["Delta_C_fixed_s"]))
        solutions.append(dict(radial_order=order, ray_order=16, samples=samples))
    return dict(schema=tool.SCHEMA, status=tool.STATUS, scope=tool.SCOPE,
        runtime={"python": "3.14.8", "numpy": "2.5.1", "scipy": "1.17.1"}, settings=dict(tool.SETTINGS),
        source_scale=1., domain=inputs["domain"], inputs=inputs["receipt"], solutions=solutions,
        diagnostics=tool.diagnostics(solutions), refinement_is_rigorous_error_enclosure=False,
        profile_evaluation=dict(counts=dict(profile_requested=100, wall_points=80, second_points=10, cache_hits=20),
            maximum_inverse_bracket_width=0.), execution=dict(cpu_seconds=1., cpu_seconds_max=120.,
            scalar_kernel_values=10000, scalar_kernel_values_max=30_000_000, metric_values=1000,
            source_quadrature_nodes=0, field_rhs_calls=10, leading_field_solves=1, wall_solves=0,
            cap_clock_solves=0, scalar_propagations=0, incoming_ray_evaluations=10,
            incoming_ray_nodes=sum(row["integrals"]["ray_nodes"] for s in solutions for row in s["samples"]),
            optical_radial_nodes=sum(row["integrals"]["optical_radial_nodes"] for s in solutions for row in s["samples"])))


@pytest.fixture(scope="module")
def authenticated_inputs():
    # These immutable predecessor records are already in the producer's
    # source inventory. Authentication performs no physical integration.
    return tool.input_snapshot()


@pytest.fixture
def record(monkeypatch, authenticated_inputs):
    inputs = deepcopy(authenticated_inputs)
    monkeypatch.setattr(tool, "input_snapshot", lambda: inputs)
    return synthetic(inputs)


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


def test_common_saved_inputs_and_central_algebra_have_separate_axes(record):
    assert tool.validate_record(record)
    assert record["inputs"]["cap_saved_wall_time_panels"] == 8
    assert record["inputs"]["current_wall_time_panels"] == 32
    assert record["diagnostics"]["maximum_absolute_quiet_recovery_defect"] <= 1e-24
    assert [row["axis"] for row in record["diagnostics"]["exterior_axes_at_radial32"]] == [
        "cache_time", "correction_grid", "ODE_tolerance"]
    for s in record["solutions"]:
        for row in s["samples"]:
            r = row["reconstruction"]
            assert r["p1_center"] == pytest.approx(sum(r[key] for key in ("exterior_trace_term",
                "arrival_rate_term", "arrival_shift_term", "outgoing_wall_term", "incoming_forcing_term")))


@pytest.mark.parametrize("path,value", [
    (("settings", "field_nodes"), 257),
    (("solutions", 0, "radial_order"), 32),
    (("solutions", 0, "samples", 2, "profile", "f"), 100.),
    (("solutions", 1, "samples", 2, "arrival_spline_derivative"), 1.),
    (("solutions", 1, "samples", 2, "integrals", "forcing_integral"), 1.),
    (("solutions", 1, "samples", 2, "integrals", "ray_nodes"), 1),
    (("solutions", 1, "samples", 2, "integrals", "center_samples", 0, "K1"), 1.),
    (("solutions", 1, "samples", 2, "reconstruction", "Delta_C_fixed_V"), 1.),
    (("solutions", 1, "samples", 2, "exterior_controls", 0, "p1_center"), 1.),
    (("execution", "scalar_kernel_values"), 30_000_001),
    (("execution", "cpu_seconds"), 121.),
    (("execution", "wall_solves"), 1),
    (("execution", "leading_field_solves"), True),
    (("refinement_is_rigorous_error_enclosure",), True),
])
def test_changed_settings_input_terms_and_budgets_fail(record, path, value):
    changed = deepcopy(record)
    target = changed
    for key in path[:-1]:
        target = target[key]
    target[path[-1]] = value
    with pytest.raises(ValueError):
        tool.validate_record(changed)


def test_creation_only_check_never_integrates(isolated, monkeypatch):
    from recursive_horizons import nsc_atlas_central_scalar as central
    from recursive_horizons import nsc_atlas_prepared_response as prepared
    tool.create(isolated)
    before = {path.name: path.read_bytes() for path in isolated.iterdir()}
    monkeypatch.setattr(tool, "compare", lambda: pytest.fail("unexpected comparison"))
    monkeypatch.setattr(central, "incoming_response_integrals", lambda *a, **k: pytest.fail("unexpected ray integral"))
    monkeypatch.setattr(prepared, "solve_prepared_field", lambda *a, **k: pytest.fail("unexpected leading solve"))
    assert tool.check(isolated)["physical_solves"] == 0
    with pytest.raises(FileExistsError):
        tool.create(isolated)
    assert before == {path.name: path.read_bytes() for path in isolated.iterdir()}
    (tool.ROOT/"src/recursive_horizons/nsc_atlas_central_scalar.py").write_bytes(b"source drift")
    with pytest.raises(ValueError, match="source binding"):
        tool.check(isolated)


def test_rehashed_reconstruction_tamper_fails(isolated):
    tool.create(isolated)
    path = isolated/"central.json"
    record = json.loads(path.read_bytes())
    record["solutions"][1]["samples"][2]["reconstruction"]["p1_center"] += 1.
    raw = tool.encoded(record)
    path.write_bytes(raw)
    manifest_path = isolated/"manifest.json"
    manifest = json.loads(manifest_path.read_bytes())
    manifest["outputs"]["central.json"] = {"sha256": tool.digest(raw), "bytes": len(raw)}
    manifest_path.write_bytes(tool.encoded(manifest))
    with pytest.raises(ValueError):
        tool.check(isolated)


def test_source_drift_during_compare_leaves_no_record(isolated, monkeypatch, record):
    def drift():
        (tool.ROOT/"docs/nsc-atlas-central-scalar-response.md").write_bytes(b"drift")
        return deepcopy(record)
    monkeypatch.setattr(tool, "compare", drift)
    with pytest.raises(RuntimeError, match="inputs or runtime changed"):
        tool.create(isolated)
    assert not isolated.exists()


def test_quiet_recovery_is_not_replaced_by_rehashed_algebra(record):
    changed = deepcopy(record)
    # A finite but invented incident quiet-corner contribution cannot become
    # acceptable merely by maintaining signed-integral bookkeeping.
    row = changed["solutions"][0]["samples"][0]
    row["profile"]["f_second"] = 1.
    with pytest.raises(ValueError, match="quiet|reconstruction"):
        tool.validate_record(changed)
