"""PG-normal observer clock. No field march and no result record."""
import hashlib
import math
import os
from pathlib import Path
import time

for _name in (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(_name, "1")

import numpy as np
import pytest

from recursive_horizons import nsc_rn_observer_clock as clock
from recursive_horizons import nsc_rn_pg as pg
from recursive_horizons import nsc_rn_source as source


ROOT = Path(__file__).resolve().parents[2]
CPU_LIMIT_SECONDS = 60.0
EVIDENCE = []
FROZEN_SHA256 = {
    "lab/src/recursive_horizons/nsc_rn_pg.py":
        "96660c727a90024dcfaf2cdd802879d6d55e5739f03ffeb29f2660ed49ffbda0",
    "lab/src/recursive_horizons/nsc_rn_source.py":
        "dd3e63cdc22042d7208ceca5635dbc4983b16aeebdb0de7ce2500eae4928a1c4",
    "lab/src/recursive_horizons/nsc_rn_observables.py":
        "cd2e3b58129acd92a513f7e03b4db283ff15011293c2ed0ae12091764e6a3307",
    "lab/src/recursive_horizons/nsc_rn_reference.py":
        "4e61faaa18ba7f82461a6e2ac7e19590786a24d7cce79cf30241079a22bbbecc",
    "lab/scripts/derive_nsc_rn_parent.py":
        "8207f06e83a2dc8118a3dc366706a07c2f0dbc3d987406afbd7cc92272211cdf",
    "lab/docs/nsc-rn-parent-episode.md":
        "afc5b7729c31967308b8e7a83f6d72e42443979c96aba1e97af95baec0b10283",
}


def _record(name, **values):
    EVIDENCE.append({"check": name, **values})
    print(f"EVIDENCE {name} {values}")


@pytest.fixture(scope="session", autouse=True)
def _cpu_budget():
    started = time.process_time()
    yield
    used = time.process_time() - started
    print(f"EVIDENCE cpu_seconds {used}")
    assert used < CPU_LIMIT_SECONDS


def _frozen_digests():
    """Digest each pinned file with line endings normalised to LF.

    The checkout is used from macOS and from Windows. A bare byte hash of the
    same content differs between them, and the pins were taken normalised, so
    normalise on read rather than pinning one platform's line endings.
    """
    found = {}
    for relative, digest in FROZEN_SHA256.items():
        data = (ROOT / relative).read_bytes().replace(b"\r\n", b"\n")
        found[relative] = hashlib.sha256(data).hexdigest()
        assert found[relative] == digest, relative
    return found


def test_frozen_producer_files_are_untouched():
    digests = _frozen_digests()
    _record("frozen_producers", files=len(digests))


def test_flat_metric_normal_and_static_clocks():
    radius = np.linspace(1.0, 5.0, 9)
    lapse = np.ones_like(radius)
    shift = np.zeros_like(radius)
    domain = (float(radius[0]), float(radius[-1]))
    observer = clock.make_normal_observer(3.0, domain, t=0.25, tau=0.0, label="flat")
    sample = clock.sample_metric(radius, lapse, shift, observer["r"])
    rate = clock.normal_observer_rate(sample, delta_r=0.0)
    assert rate["dr_dt"] == 0.0
    assert rate["dtau_dt"] == 1.0
    assert rate["normalization"] == pytest.approx(1.0)
    moved = observer
    for _ in range(4):
        moved = clock.rk4_normal_step(moved, 0.1, evaluator=_flat_evaluator)
    assert moved["status"] == "active"
    assert moved["r"] == pytest.approx(3.0)
    assert moved["t"] == pytest.approx(0.65)
    assert moved["tau"] == pytest.approx(0.4)
    static = clock.make_static_observer(3.0, domain, t=0.0, tau=0.0)
    static = clock.rk4_static_step(static, 0.3, evaluator=_flat_evaluator)
    assert static["r"] == pytest.approx(3.0)
    assert static["tau"] == pytest.approx(0.3)
    assert static["static_stage_rates"] == [1.0, 1.0, 1.0, 1.0]
    _record("flat", tau=moved["tau"], static_tau=static["tau"], normalization=rate["normalization"])


def _flat_evaluator(time_value, radius):
    return {"N": 1.0, "beta": 0.0, "N_r": 0.0, "beta_r": 0.0, "delta_N": 0.0, "delta_beta": 0.0}


def test_static_rn_proper_and_geodesic_normal_relation():
    grid = pg.build_grid(1.01, r_m=1.0, points=33, r_out=32.0, layout="dense")
    model = grid["model"]
    positions = clock.initial_observers(grid, t0=0.5, tau0=0.0)
    parent = positions["parent_packet_centre"]
    assert parent["label"] == "parent_packet_centre"
    assert parent["t"] == pytest.approx(0.5)
    assert parent["r"] == pytest.approx(12.0 * grid["r_m"])
    assert parent["status"] == "active"
    assert grid["parent_interval"][0] < parent["r"] < grid["parent_interval"][1]
    centre = parent["r"]
    beta = float(model.beta_rn(centre))
    beta_r = float(model.beta_derivative(centre))
    metric_tt = 1.0 - beta ** 2
    assert metric_tt == pytest.approx(float(model.f(centre)))
    assert metric_tt > 0.0
    sample = {
        "status": "active",
        "r": centre,
        "N": 1.0,
        "beta": beta,
        "N_r": 0.0,
        "beta_r": beta_r,
        "delta_N": 0.0,
        "delta_beta": 0.0,
    }
    normal = clock.normal_observer_rate(sample)
    static = clock.static_areal_rate(sample)
    nodal_lapse = np.ones_like(grid["radius"])
    nodal_shift = np.asarray(model.beta_rn(grid["radius"]), dtype=float)
    node_index = int(np.argmin(np.abs(grid["radius"] - centre)))
    node_radius = float(grid["radius"][node_index])
    nodal_sample = clock.sample_metric(grid["radius"], nodal_lapse, nodal_shift, node_radius)
    assert nodal_sample["beta"] == pytest.approx(float(nodal_shift[node_index]), abs=1.0e-12)
    assert nodal_sample["N"] == pytest.approx(1.0)
    assert normal["dtau_dt"] == pytest.approx(1.0)
    assert normal["dr_dt"] == pytest.approx(-beta)
    assert normal["normalization"] == pytest.approx(1.0)
    assert static["status"] == "active"
    assert static["rate"] == pytest.approx(math.sqrt(metric_tt))
    assert static["rate"] ** 2 + beta ** 2 == pytest.approx(normal["dtau_dt"] ** 2)
    killing_energy = (1.0 - beta ** 2) * 1.0 + (-beta) * (-beta)
    assert killing_energy == pytest.approx(1.0)

    def _vacuum(time_value, radius):
        return {
            "N": 1.0,
            "beta": float(model.beta_rn(radius)),
            "N_r": 0.0,
            "beta_r": float(model.beta_derivative(radius)),
            "delta_N": 0.0,
            "delta_beta": 0.0,
        }

    coarse = clock.integrate_normal(parent, 0.05, 8, evaluator=_vacuum)
    fine = clock.integrate_normal(parent, 0.005, 80, evaluator=_vacuum)
    assert coarse["stopped"] == "active"
    assert fine["stopped"] == "active"
    assert coarse["observer"]["tau"] - 0.0 == pytest.approx(0.4)
    assert fine["observer"]["tau"] == pytest.approx(coarse["observer"]["tau"])
    assert coarse["observer"]["r"] < parent["r"]
    assert abs(coarse["observer"]["r"] - fine["observer"]["r"]) < 1.0e-8
    excision = float(grid["radius"][0])
    excision_shift = float(model.beta_rn(excision))
    excision_gtt = 1.0 - excision_shift ** 2
    assert excision_gtt < 0.0
    excision_sample = clock.sample_metric(grid["radius"], nodal_lapse, nodal_shift, excision)
    assert clock.static_areal_rate(excision_sample)["status"] == "static_not_timelike"
    _record(
        "static_rn",
        centre=centre,
        beta=beta,
        static_rate=static["rate"],
        killing_energy=killing_energy,
        radius_error=abs(coarse["observer"]["r"] - fine["observer"]["r"]),
        excision_g_tt=excision_gtt,
        coarse_r=coarse["observer"]["r"],
    )


def test_historical_scalars_match_producer_and_do_not_move_the_field():
    grid = pg.build_grid(1.01, r_m=1.0, points=17, r_out=32.0, layout="dense")
    state = pg.make_state(grid, killing_frequency=1.0, inner_mass=grid["mass"])
    first = pg.stage_rates(state, grid)
    labeled = clock.historical_producer_rates(first["lapse"], first["shift"])
    assert labeled["normal_rate"] == first["normal_rate"]
    assert labeled["killing_rate"] == first["killing_rate"]
    assert labeled["normal_rate"] == float(first["lapse"][0])
    outer_gtt = float(first["lapse"][-1] ** 2 - first["shift"][-1] ** 2)
    assert outer_gtt > 0.0
    assert labeled["killing_rate"] == math.sqrt(outer_gtt)
    assert labeled["normal_is_pg_normal_proper_time"] is False
    assert labeled["normal_is_child_clock"] is False
    assert labeled["affects_phi_or_metric"] is False
    excision_gtt = float(first["lapse"][0] ** 2 - first["shift"][0] ** 2)
    assert excision_gtt < 0.0
    outer_sample = clock.sample_metric(grid["radius"], first["lapse"], first["shift"], grid["radius"][-1])
    assert clock.static_areal_rate(outer_sample)["rate"] == pytest.approx(first["killing_rate"])
    clocks = {"t": 1.25, "killing": 4.5, "normal": -2.0}
    renamed = clock.label_saved_clocks(clocks)
    assert renamed["t"] == clocks["t"]
    assert renamed["killing"] == clocks["killing"]
    assert renamed["normal"] == clocks["normal"]
    assert renamed["normal_is_pg_normal_proper_time"] is False
    state["clocks"] = {"t": 30.0, "killing": -8.0, "normal": 12.0}
    second = pg.stage_rates(state, grid)
    assert np.max(np.abs(first["phi_rate"] - second["phi_rate"])) == 0.0
    assert np.max(np.abs(first["lapse"] - second["lapse"])) == 0.0
    assert np.max(np.abs(first["shift"] - second["shift"])) == 0.0
    manufactured_lapse = np.linspace(1.4, 1.0, grid["points"])
    manufactured_shift = np.full(grid["points"], 0.1)
    manufactured = clock.historical_producer_rates(manufactured_lapse, manufactured_shift)
    assert manufactured["normal_rate"] == pytest.approx(1.4)
    assert manufactured["normal_rate"] != pytest.approx(float(manufactured_lapse[grid["points"] // 2]))
    clamped = clock.historical_producer_rates(np.array([1.0, 1.0]), np.array([2.0, 2.0]))
    assert clamped["killing_rate"] == 0.0
    assert clamped["killing_branch"] == "producer_clamped_to_zero"
    _record(
        "historical_scalars",
        normal_rate=labeled["normal_rate"],
        killing_rate=labeled["killing_rate"],
        excision_g_tt=excision_gtt,
        outer_g_tt=outer_gtt,
        phi_rate_changed=False,
    )


def test_manufactured_varying_metric_trajectory():
    # β = 0.15 + 0.05 t, N = 0.8 + 0.1 t + 0.02 r. RK4 is exact for this polynomial.
    r0 = 4.0
    domain = (0.5, 8.0)

    def _evaluator(time_value, radius):
        time_value = float(time_value)
        radius = float(radius)
        return {
            "N": 0.8 + 0.1 * time_value + 0.02 * radius,
            "beta": 0.15 + 0.05 * time_value,
            "N_r": 0.02,
            "beta_r": 0.0,
            "delta_N": 0.0,
            "delta_beta": 0.0,
        }

    def _exact(time_value):
        radius = r0 - 0.15 * time_value - 0.025 * time_value ** 2
        proper = (0.8 + 0.02 * r0) * time_value + 0.0485 * time_value ** 2 - (0.0005 / 3.0) * time_value ** 3
        return radius, proper

    observer = clock.make_normal_observer(r0, domain, t=0.0, tau=0.0)
    result = clock.integrate_normal(observer, 0.2, 3, evaluator=_evaluator)
    final = result["observer"]
    radius, proper = _exact(final["t"])
    radius_error = abs(final["r"] - radius)
    proper_error = abs(final["tau"] - proper)
    assert final["status"] == "active"
    assert radius_error < 1.0e-12
    assert proper_error < 1.0e-12
    _record(
        "manufactured_trajectory",
        t=final["t"],
        radius_error=radius_error,
        proper_error=proper_error,
        r=final["r"],
        tau=final["tau"],
    )


def test_domain_exit_preserves_last_valid_state():
    domain = (1.0, 5.0)
    observer = clock.make_normal_observer(1.4, domain, t=0.0, tau=0.2)

    def _inward(time_value, radius):
        return {"N": 1.0, "beta": 3.0, "N_r": 0.0, "beta_r": 0.0}

    refused = clock.rk4_normal_step(observer, 0.2, evaluator=_inward)
    assert refused["status"] == "domain_exit"
    assert refused["r"] == observer["r"]
    assert refused["t"] == observer["t"]
    assert refused["tau"] == observer["tau"]
    assert refused["last_valid"]["r"] == observer["r"]
    assert refused["exit"]["attempted_r"] < domain[0]
    assert refused["clock_advanced"] is False
    held = clock.rk4_normal_step(refused, 0.2, evaluator=_inward)
    assert held["r"] == refused["r"]
    assert held["tau"] == refused["tau"]
    assert held["t"] == refused["t"]
    calls = {"count": 0}

    def _counted(time_value, radius):
        calls["count"] += 1
        return _inward(time_value, radius)

    clock.rk4_normal_step(refused, 0.05, evaluator=_counted)
    assert calls["count"] == 0
    advanced = clock.rk4_normal_step(observer, 0.05, evaluator=_inward)
    assert advanced["status"] == "active"
    assert advanced["r"] == pytest.approx(1.25)
    assert advanced["tau"] == pytest.approx(0.25)

    def _negative_lapse(time_value, radius):
        return {"N": -0.2, "beta": 0.0, "N_r": 0.0, "beta_r": 0.0}

    refused_lapse = clock.rk4_normal_step(observer, 0.05, evaluator=_negative_lapse)
    assert refused_lapse["status"] == "lapse_not_positive"
    assert refused_lapse["r"] == observer["r"]
    assert refused_lapse["tau"] == observer["tau"]
    assert refused_lapse["t"] == observer["t"]
    _record(
        "domain_exit",
        kept_r=refused["r"],
        attempted_r=refused["exit"]["attempted_r"],
        advanced_r=advanced["r"],
    )


def test_child_outer_edge_exits_and_parent_stays():
    grid = pg.build_grid(1.01, r_m=1.0, points=33, r_out=32.0, layout="dense")
    model = grid["model"]
    positions = clock.initial_observers(grid, t0=0.0)
    child = positions["child_outer_edge"]
    parent = positions["parent_packet_centre"]
    assert child["r"] == pytest.approx(grid["child_interval"][1])
    assert child["status"] == "active"
    assert child["r"] > float(grid["radius"][0])

    def _vacuum(time_value, radius):
        return {
            "N": 1.0,
            "beta": float(model.beta_rn(radius)),
            "N_r": 0.0,
            "beta_r": float(model.beta_derivative(radius)),
        }

    exited = clock.integrate_normal(child, 0.05, 20, evaluator=_vacuum)
    assert exited["stopped"] == "domain_exit"
    assert exited["accepted_steps"] >= 1
    assert exited["observer"]["r"] >= float(grid["radius"][0])
    assert exited["observer"]["tau"] == pytest.approx(exited["accepted_steps"] * 0.05)
    again = clock.rk4_normal_step(exited["observer"], 0.05, evaluator=_vacuum)
    assert again["r"] == exited["observer"]["r"]
    assert again["t"] == exited["observer"]["t"]
    assert again["tau"] == exited["observer"]["tau"]
    stayed = clock.integrate_normal(parent, 0.05, 6, evaluator=_vacuum)
    assert stayed["stopped"] == "active"
    assert stayed["observer"]["r"] > grid["child_interval"][1]
    _record(
        "child_exit",
        accepted_steps=exited["accepted_steps"],
        exit_r=exited["observer"]["r"],
        exit_tau=exited["observer"]["tau"],
        parent_r=stayed["observer"]["r"],
    )


def test_static_spacelike_refusal_preserves_last_valid_state():
    domain = (1.0, 20.0)
    static = clock.make_static_observer(12.0, domain, t=0.0, tau=0.0)

    def _timelike(time_value, radius):
        return {"N": 1.0, "beta": 0.2, "N_r": 0.0, "beta_r": 0.0}

    valid = clock.rk4_static_step(static, 0.1, evaluator=_timelike)
    assert valid["status"] == "active"
    assert valid["tau"] == pytest.approx(0.1 * math.sqrt(1.0 - 0.04))
    kept_tau = valid["tau"]

    def _spacelike(time_value, radius):
        return {"N": 1.0, "beta": 2.0, "N_r": 0.0, "beta_r": 0.0}

    refused = clock.rk4_static_step(valid, 0.1, evaluator=_spacelike)
    assert refused["status"] == "static_not_timelike"
    assert refused["tau"] == kept_tau
    assert refused["t"] == valid["t"]
    assert refused["r"] == valid["r"]
    assert refused["last_valid"]["tau"] == kept_tau
    assert refused["clock_advanced"] is False
    sample = {
        "status": "active", "N": 1.0, "beta": 1.0, "r": 12.0,
    }
    null_curve = clock.static_areal_rate(sample)
    assert null_curve["status"] == "static_not_timelike"
    assert null_curve["rate"] is None
    assert null_curve["clamped_to_zero"] is False
    assert null_curve["g_tt"] == pytest.approx(0.0)
    _record("static_refusal", kept_tau=kept_tau, null_g_tt=null_curve["g_tt"])


def test_interpolation_ties_edges_and_tangent_at_two_amplitudes():
    radius = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
    values = np.array([0.0, 1.0, 3.0, 6.0, 10.0])
    lapse = np.ones_like(radius)
    left = clock.sample_metric(radius, lapse, values, 0.0)
    right = clock.sample_metric(radius, lapse, values, 4.0)
    node = clock.sample_metric(radius, lapse, values, 2.0)
    opening = clock.sample_metric(radius, lapse, values, 2.5)
    outside = clock.sample_metric(radius, lapse, values, 4.1)
    assert left["edge"] == "left"
    assert left["beta_r"] == pytest.approx(1.0)
    assert left["tie"] is False
    assert right["edge"] == "right"
    assert right["slope_side"] == "left_secant_of_last_interval"
    assert right["beta_r"] == pytest.approx(4.0)
    assert node["tie"] is True
    assert node["beta"] == pytest.approx(3.0)
    assert node["beta_r"] == pytest.approx(3.0)
    assert node["beta_r_left"] == pytest.approx(2.0)
    assert node["beta_r_right"] == pytest.approx(3.0)
    straddling = (
        clock.sample_metric(radius, lapse, values, 2.2)["beta"]
        - clock.sample_metric(radius, lapse, values, 1.8)["beta"]
    ) / 0.4
    assert straddling == pytest.approx(0.5 * (2.0 + 3.0))
    assert straddling != pytest.approx(node["beta_r"])
    one_sided = (clock.sample_metric(radius, lapse, values, 2.2)["beta"] - node["beta"]) / 0.2
    assert one_sided == pytest.approx(node["beta_r"])
    assert opening["tie"] is False
    assert opening["beta_r"] == pytest.approx(3.0)
    assert outside["status"] == "domain_exit"
    assert outside["extrapolated"] is False
    assert outside["N"] is None

    nodes = np.linspace(2.0, 10.0, 17)
    base_lapse = np.ones_like(nodes)
    base_shift = 0.1 * nodes
    delta_lapse = 0.05 * nodes
    delta_shift = nodes.copy()
    domain = (float(nodes[0]), float(nodes[-1]))
    origin = 4.0
    step = 0.05
    count = 4

    def _frames(amplitude):
        lapse = base_lapse + amplitude * delta_lapse
        shift = base_shift + amplitude * delta_shift

        def _evaluator(time_value, radius_value, _lapse=lapse, _shift=shift):
            return clock.sample_metric(nodes, _lapse, _shift, radius_value)

        return _evaluator

    base_rate = clock.normal_observer_rate(clock.sample_metric(
        nodes, base_lapse, base_shift, origin,
        delta_lapse=delta_lapse, delta_shift=delta_shift,
    ))
    assert base_rate["dr_dt"] == pytest.approx(-0.1 * origin)
    assert base_rate["dtau_dt"] == pytest.approx(1.0)
    assert base_rate["d_delta_r_dt"] == pytest.approx(-origin)
    assert base_rate["d_delta_tau_dt"] == pytest.approx(0.05 * origin)
    tangent = clock.make_normal_observer(origin, domain)
    frame = {
        "radius": nodes,
        "lapse": base_lapse,
        "shift": base_shift,
        "delta_lapse": delta_lapse,
        "delta_shift": delta_shift,
    }
    for _ in range(count):
        tangent = clock.rk4_normal_step(tangent, step, stages=(frame, frame, frame, frame))
    radius_errors = {}
    proper_errors = {}
    for amplitude in (1.0e-2, 1.0e-3):
        ahead = clock.make_normal_observer(origin, domain)
        behind = clock.make_normal_observer(origin, domain)
        for _ in range(count):
            ahead = clock.rk4_normal_step(ahead, step, evaluator=_frames(amplitude))
            behind = clock.rk4_normal_step(behind, step, evaluator=_frames(-amplitude))
        radius_errors[amplitude] = abs((ahead["r"] - behind["r"]) / (2.0 * amplitude) - tangent["delta_r"])
        proper_errors[amplitude] = abs((ahead["tau"] - behind["tau"]) / (2.0 * amplitude) - tangent["delta_tau"])
        assert ahead["status"] == "active"
    assert radius_errors[1.0e-2] > 10.0 * radius_errors[1.0e-3]
    assert proper_errors[1.0e-2] > 10.0 * proper_errors[1.0e-3]
    assert radius_errors[1.0e-3] < 1.0e-6
    assert proper_errors[1.0e-3] < 1.0e-6
    _record(
        "tangent",
        radius_error_1em2=radius_errors[1.0e-2],
        radius_error_1em3=radius_errors[1.0e-3],
        radius_error_ratio=radius_errors[1.0e-2] / radius_errors[1.0e-3],
        proper_error_1em2=proper_errors[1.0e-2],
        proper_error_1em3=proper_errors[1.0e-3],
        proper_error_ratio=proper_errors[1.0e-2] / proper_errors[1.0e-3],
        delta_r=tangent["delta_r"],
        delta_tau=tangent["delta_tau"],
        tie_right_slope=node["beta_r"],
        tie_left_slope=node["beta_r_left"],
        straddling_slope=straddling,
    )


def test_saved_history_stride_is_a_sampling_indicator_only():
    nodes = np.linspace(1.0, 8.0, 15)
    times = np.arange(0.0, 2.0 + 1.0e-12, 0.5)
    samples = []
    for time_value in times:
        beta = (0.1 + 0.04 * time_value ** 2) * np.ones_like(nodes)
        samples.append({
            "t": float(time_value),
            "radius": nodes,
            "lapse": np.ones_like(nodes),
            "shift": beta,
        })
    reconstructed = clock.reconstruct_saved_history(samples, r0=4.0, tau0=0.0, strides=(0.5, 1.0), substeps=4)
    assert reconstructed["exact_old_run_proper_clock"] is False
    assert reconstructed["claim"] == "sampling_indicator_only"
    fine = reconstructed["by_stride"][0.5]["observer"]
    coarse = reconstructed["by_stride"][1.0]["observer"]
    exact_r = 4.0 - (0.1 * 2.0 + 0.04 * (2.0 ** 3) / 3.0)
    fine_error = abs(fine["r"] - exact_r)
    coarse_error = abs(coarse["r"] - exact_r)
    assert fine["status"] == "active" and coarse["status"] == "active"
    assert abs(reconstructed["sampling_difference"]["delta_r"]) > 1.0e-4
    assert fine_error < coarse_error
    scope = clock.saved_history_scope()
    assert scope["exact_old_run_proper_clock_recoverable"] is False
    assert scope["evolved_lapse_and_shift_stored"] is False
    assert scope["rk4_stage_history_stored"] is False
    _record(
        "saved_history",
        fine_radius_error=fine_error,
        coarse_radius_error=coarse_error,
        stride_delta_r=reconstructed["sampling_difference"]["delta_r"],
        fine_tau=fine["tau"],
        exact_old_run=False,
    )


def test_stage_order_accepts_coupled_metric_frames():
    nodes = np.linspace(2.0, 6.0, 9)
    domain = (float(nodes[0]), float(nodes[-1]))
    observer = clock.make_normal_observer(4.0, domain)
    frames = []
    for index, beta in enumerate((0.2, 0.3, 0.3, 0.4)):
        frames.append({
            "radius": nodes,
            "lapse": np.full(nodes.shape, 1.0 + 0.01 * index),
            "shift": np.full(nodes.shape, beta),
        })
    # The two midpoint frames share β = 0.3 but the lapse stamps differ.
    # A stage-blind average would not reproduce the RK4 weights on the lapse.
    bare_frames = [{"lapse": frame["lapse"], "shift": frame["shift"]} for frame in frames]
    updated = clock.rk4_normal_step(observer, 0.2, stages=bare_frames, radius=nodes)
    lapse_rates = [1.0, 1.01, 1.02, 1.03]
    expected_tau = 0.2 / 6.0 * (lapse_rates[0] + 2.0 * lapse_rates[1] + 2.0 * lapse_rates[2] + lapse_rates[3])
    expected_r = 4.0 - 0.2 / 6.0 * (0.2 + 2.0 * 0.3 + 2.0 * 0.3 + 0.4)
    assert updated["tau"] == pytest.approx(expected_tau)
    assert updated["r"] == pytest.approx(expected_r)
    assert updated["writes_producer_clocks"] is False
    _record("rk4_stages", tau_error=abs(updated["tau"] - expected_tau), r_error=abs(updated["r"] - expected_r))
