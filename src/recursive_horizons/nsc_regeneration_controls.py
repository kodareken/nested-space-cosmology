"""Controls and one continuation for the saved spherical Galerkin episode.

The evolution, the column source, the radius operator, and the shift
momentum stay with their owners. This module changes occupations or holds
the geometry fixed, then calls those owners. It does not add a radius force,
subtract a mean, or replace the episode step.

Active note. G > 0 is the initial convex-radius hypothesis for constant Q,
chi = p_r = p_chi = 0 and r = y^2. The evolving chart is the owned
chart_failure: a finite state with r, Q and L positive. min G stays a
diagnostic. The Galerkin generator is H_G = U_f† H_fine U_f on 2 n_f.
(1/2) I_{2 n_f} commutes with H_G. Its embedding in 2 n_q need not commute
with the unprojected H_fine. The fine-grid identity remains an algebraic
source control, not an occupation vector.
"""
from __future__ import annotations

from dataclasses import replace

import numpy as np

from recursive_horizons import nsc_spherical_galerkin_coupling as galerkin
from recursive_horizons.nsc_spherical_cauchy_data import shift_momentum
from recursive_horizons.nsc_spherical_coupling import (
    CauchyRate,
    CauchyState,
    PositiveChartExit,
    _combine,
    chart_failure,
    magnetic_radius_square,
    shift_constraint,
    source_from_columns,
)

BASELINE_OCCUPATIONS = np.array([0.75, 0.75, 0.5, 0.5, 0.25, 0.25], dtype=float)
UNIFORM_OCCUPATIONS = np.full(6, 0.5, dtype=float)
REVERSED_OCCUPATIONS = np.array([0.25, 0.25, 0.5, 0.5, 0.75, 0.75], dtype=float)

# Declared from the saved T=0.05 episode before any continuation.
# Leader shell content moves by about 0.108 and proper velocity reaches about 0.088.
CONTENT_EFFECT = 0.10
SHARE_MAINTAINED = 0.50
SHARE_PACKET = 0.30
SHARE_STOP = 0.35
EVOLUTION_PROJECTED_MAX = 1e-4
EVOLUTION_FULL_MAX = 1e-2
CAR_GAP_MAX = 1e-8
JOIN_TIME = 0.05
FINE_IDENTITY_IS_ALGEBRAIC_CONTROL = True
FINE_IDENTITY_IS_AN_OCCUPATION = False
# The proxy candidate_regime does not authorize a perturbation campaign.
PERTURBATIONS_FROM_PROXY = False
CONTROL_FLOORS = {
    "leader_shell": CONTENT_EFFECT,
    "field_energy": CONTENT_EFFECT,
    "chi_max": 1.0,
    "proper_velocity": 0.01,
    "realized_radius": 1.0e-3,
}

CRITERION = {
    "localization": (
        "Normal-shell content on the four fixed regional windows. "
        "The leader is the window with the largest share. The packet is every "
        "window whose share is at least 0.30. Localized means the leader share "
        "is at least 0.50."
    ),
    "maintained": (
        "On the second episode alone, the leader at its start remains the leader, "
        "the end share is at least 0.50, and that window's shell content stays "
        "above half its value at the start of the second episode."
    ),
    "renewed": (
        "Judged on the joined series from T=0, not on the second slice alone. "
        "Renewal requires the end state still localized and either a reversal of "
        "the original leader's shell content of at least 0.10, or a different "
        "window becoming the leader after gaining at least 0.10. Two segments of "
        "one drift are not renewal."
    ),
    "nonlocalization_stop": (
        "Stop when the leader share falls below 0.35. The crossing sample is kept."
    ),
    "initial_convex_radius": (
        "G > 0 is the initial convex-radius hypothesis only, for constant Q, "
        "chi = p_r = p_chi = 0 and r = y^2, with G independent of y. It is the "
        "hypothesis under which J(y) is strictly convex. It is not the validity "
        "test of an evolving Lorentzian chart."
    ),
    "dynamic_chart": (
        "The evolving chart is chart_failure: a finite state, Q > 0, r > 0, and "
        "L > 0. Lapse N = r L is positive when r and L are. min G is recorded "
        "and does not stop the step."
    ),
    "chart_stop": (
        "Stop on the owned chart_failure exit for r, Q, or L, including a "
        "nonfinite state. min G is a diagnostic of the initial convex-radius "
        "hypothesis and is not a dynamic stop. A failed chart step is not kept. "
        "The last open sample is checkpointed."
    ),
    "effect_scale": {
        "shell_content": CONTENT_EFFECT,
        "leader_share": SHARE_MAINTAINED,
        "basis": (
            "Saved nf512 dt=5e-4 episode: leader shell change about +0.108, "
            "proper velocity about [-0.062, +0.088], chi_max about 11.1. "
            "Comparison uses that effect, not one percent of a constraint residual."
        ),
    },
    "evolution_gate": (
        "A changed-source initial state is admitted when G, r, and Q are positive, "
        "under the convex-radius hypothesis of constant Q, chi = p_r = p_chi = 0 "
        "and r = y^2. Dynamic continuation does not stop when min G leaves that "
        "bracket. It stops on chart_failure, checks the column Gram as Q changes, "
        "and checkpoints the exit. The column Gram stays within 1e-8 of the "
        "identity, the radius is source-derived, the projected Hamilton residual "
        "is at most 1e-4, and the full quadrature residual is at most 1e-2. The "
        "1e-2 line is the saved nf256 episode's residual class, not a new physical "
        "law. A finite-grid stall above the Newton tolerance is a solver floor, "
        "not a physical failure. The residual of the best positive source iterate "
        "is reported."
    ),
    "proxies": (
        "Window reversal of at least 0.10, a new leader that gains at least 0.10, "
        "and a leader share of at least 0.50 are proxies. candidate_regime is "
        "their disjunction with the stable-throughflow proxy. They are not "
        "programme requirements. The programme allows a maintained structure "
        "that carries throughflow. Splitting one drift is not complete renewal."
    ),
    "algebraic_control": (
        "H_G = U_f† H_fine U_f acts on 2 n_f. (1/2) I_{2 n_f} commutes with H_G. "
        "The 2 n_q embedding of that identity need not commute with the "
        "unprojected H_fine. The fine identity remains an algebraic source "
        "control. Population control accepts six weights and refuses a "
        "quadrature-sized occupation."
    ),
    "uniform_complement": (
        "Uniform six weights remain localized relative to the empty complement. "
        "A four-window share below 0.50 is a different proxy."
    ),
    "mean_current": (
        "The current array is not edited. A nonzero mean remains in the shift "
        "residual because a periodic derivative cannot cancel it."
    ),
    "cpu_coarse": 600.0,
    "cpu_fine": 1800.0,
    "fine_only_if": (
        "A separated content, flux, reversal, or accounting comparison that can "
        "change the reading, or a resolution disagreement on that reading. The "
        "proxy candidate_regime does not by itself open the fine budget."
    ),
}


def with_occupations(grid, occupations):
    """Copy the grid with a new occupation vector. The caller's grid is unchanged."""
    values = np.asarray(occupations, dtype=float).reshape(-1).copy()
    if values.shape != tuple(np.shape(grid.fine.occupations)):
        raise ValueError("occupations must replace the six owned weights")
    if not np.isfinite(values).all():
        raise ValueError("occupations are not finite")
    fine = replace(grid.fine, occupations=values)
    return replace(grid, fine=fine)


def car_report(state):
    """Column Gram of the fermion modes. Occupations are not part of the Gram."""
    gram = state.phi0.conj().T @ state.phi0 + state.phi1.conj().T @ state.phi1
    gap = float(np.max(np.abs(gram - np.eye(gram.shape[0]))))
    return {
        "gap": gap,
        "admissible": bool(gap <= CAR_GAP_MAX),
        "columns": int(gram.shape[0]),
        "identity_is_the_mode_car": True,
    }


def bracket_g(q_values, rho, grid):
    """Owned bracket G = Q^2 r_mag^2 + Q rho / (8 pi A)."""
    q_values = np.asarray(q_values, dtype=float)
    rho = np.asarray(rho, dtype=float)
    rmag2 = float(magnetic_radius_square(grid.fine.coefficients))
    eight_pi_a = 8.0 * np.pi * float(grid.fine.A)
    return q_values ** 2 * rmag2 + q_values * rho / eight_pi_a


def source_arrays(grid, state):
    fine = galerkin.prolong_state(grid, state)
    source = source_from_columns(grid.fine, fine)
    rho = np.array(source["force_L"] / grid.dx_q, dtype=float, copy=True)
    current = np.array(source["force_beta"] / grid.dx_q, dtype=float, copy=True)
    return source, rho, current


def _score_radius(grid, radius, rho):
    fine_r = galerkin.prolong_geometry(grid, radius)
    if not np.isfinite(fine_r).all() or float(np.min(fine_r)) <= 0.0:
        return None
    projected, full = galerkin.projected_radius_operator(grid, radius, rho)
    held = full - galerkin.prolong_geometry(grid, projected)
    return {
        "projected_max": float(np.max(np.abs(projected))),
        "full_max": float(np.max(np.abs(full))),
        "held_out_max": float(np.max(np.abs(held))),
        "r_min": float(np.min(fine_r)),
        "r_max": float(np.max(fine_r)),
        "radius": np.array(radius, dtype=float, copy=True),
    }


def _newton_target(grid, radius, target, rho, iterations, accepted_steps, best):
    """One damped Newton sequence. ``best`` is scored on the full source rho."""
    stalled = False
    singular = None
    converged = False
    residual_max = None
    for _iteration in range(int(iterations)):
        projected, _full = galerkin.projected_radius_operator(grid, radius, target)
        residual_max = float(np.max(np.abs(projected)))
        scored = _score_radius(grid, radius, rho)
        if scored is not None and (best is None or scored["projected_max"] < best["projected_max"]):
            best = scored
        if residual_max < galerkin.TOL_NEWTON:
            converged = True
            break
        jacobian = galerkin._projected_radius_jacobian(grid, radius)
        try:
            delta = np.linalg.solve(jacobian, -projected)
        except np.linalg.LinAlgError as error:
            singular = str(error)
            stalled = True
            break
        accepted = False
        step = 1.0
        for _halving in range(16):
            trial = radius + step * delta
            if not np.isfinite(trial).all() or float(np.min(galerkin.prolong_geometry(grid, trial))) <= 0.0:
                step *= 0.5
                continue
            trial_projected, _unused = galerkin.projected_radius_operator(grid, trial, target)
            trial_max = float(np.max(np.abs(trial_projected)))
            if trial_max < residual_max * (1.0 - 0.05 * step):
                radius = trial
                accepted = True
                accepted_steps += 1
                scored = _score_radius(grid, radius, rho)
                if scored is not None and (best is None or scored["projected_max"] < best["projected_max"]):
                    best = scored
                break
            step *= 0.5
        if not accepted:
            stalled = True
            break
    return radius, best, accepted_steps, {
        "converged": converged,
        "residual_max": residual_max,
        "stalled": stalled,
        "singular": singular,
    }


def solve_source_radius(grid, state, *, iterations=12, scales=9):
    """Best positive radius for the full column source.

    The owned projected Jacobian and residual are the update. A stall above
    ``TOL_NEWTON`` keeps that iterate and reports the residual. It is not a
    physical failure. The vacuum radius is returned only when no source step
    improves the full residual.
    """
    q0 = float(grid.fine.calibration["b0"] / grid.fine.calibration["a0"])
    if not np.isfinite(q0) or q0 <= 0.0:
        raise ValueError("gauge ratio b0/a0 is not a positive finite constant")
    prepared = state.copy()
    prepared.Q = np.full(grid.ng, q0)
    prepared.chi = np.zeros(grid.ng)
    prepared.p_Q = np.zeros(grid.ng)
    prepared.p_r = np.zeros(grid.ng)
    prepared.p_chi = np.zeros(grid.ng)
    source, rho, current = source_arrays(grid, prepared)
    current_bytes = np.ascontiguousarray(current).tobytes()
    varied = prepared.copy()
    varied.r = prepared.r * 1.7
    other = source_from_columns(grid.fine, galerkin.prolong_state(grid, varied))
    rho_variation = float(np.max(np.abs(other["force_L"] - source["force_L"])))
    current_variation = float(np.max(np.abs(other["force_beta"] - source["force_beta"])))
    gee = bracket_g(q0, rho, grid)
    radius = np.ones(grid.ng, dtype=float)
    best = _score_radius(grid, radius, rho)
    accepted_steps = 0
    homotopy = []
    singular = None
    for scale in np.linspace(0.0, 1.0, int(scales))[1:]:
        radius, best, accepted_steps, step_info = _newton_target(
            grid, radius, float(scale) * rho, rho, iterations, accepted_steps, best,
        )
        homotopy.append({"scale": float(scale), **{key: step_info[key] for key in ("converged", "residual_max", "stalled")}})
        if step_info["singular"]:
            singular = step_info["singular"]
            break
    radius, best, accepted_steps, polish = _newton_target(
        grid, best["radius"] if best is not None else radius, rho, rho, iterations, accepted_steps, best,
    )
    if polish["singular"] and singular is None:
        singular = polish["singular"]
    if best is None:
        return prepared, {
            "converged": False,
            "source_derived": False,
            "solver_floor": False,
            "physical_failure": True,
            "blocker": "no positive radius iterate",
            "projected_max": None,
            "full_max": None,
            "held_out_max": None,
            "G_min": float(np.min(gee)),
            "positive_G": bool(np.min(gee) > 0.0),
            "positive_Q": True,
            "current": current,
            "current_bytes": current_bytes,
            "rho_variation": rho_variation,
            "current_variation": current_variation,
            "accepted_steps": accepted_steps,
            "homotopy": homotopy,
            "singular": singular,
        }
    solved = prepared.copy()
    solved.r = best["radius"]
    converged = bool(best["projected_max"] < galerkin.TOL_NEWTON)
    moved = float(np.max(np.abs(best["radius"] - 1.0))) > 1e-12
    source_derived = bool(converged or moved)
    floor = bool(source_derived and not converged)
    after, rho_after, current_after = source_arrays(grid, solved)
    del after, rho_after
    return solved, {
        "converged": converged,
        "source_derived": source_derived,
        "solver_floor": floor,
        "physical_failure": False,
        "blocker": None if source_derived else "no accepted source step; vacuum radius kept",
        "projected_max": best["projected_max"],
        "full_max": best["full_max"],
        "held_out_max": best["held_out_max"],
        "r_min": best["r_min"],
        "r_max": best["r_max"],
        "G_min": float(np.min(gee)),
        "G_max": float(np.max(gee)),
        "positive_G": bool(np.min(gee) > 0.0),
        "positive_Q": True,
        "Q": q0,
        "current": current,
        "current_unchanged": bool(np.ascontiguousarray(current_after).tobytes() == current_bytes),
        "current_mean": float(np.mean(current)),
        "rho_mean": float(np.mean(rho)),
        "rho_variation": rho_variation,
        "current_variation": current_variation,
        "source_mean_subtracted": False,
        "accepted_steps": int(accepted_steps),
        "homotopy": homotopy,
        "polish_converged": bool(polish["converged"]),
        "singular": singular,
        "newton_tolerance": float(galerkin.TOL_NEWTON),
        "floor_is_not_a_physical_failure": True,
    }


def attach_shift_momentum(grid, state, current):
    """Mean-free periodic p_Q. The supplied current array is left as it is."""
    snapshot = np.array(current, dtype=float, copy=True)
    momentum, info = shift_momentum(grid, snapshot)
    if not np.array_equal(snapshot, np.asarray(current)):
        raise RuntimeError("shift momentum edited the current")
    updated = state.copy()
    updated.p_Q = np.array(momentum, dtype=float, copy=True)
    updated.p_r = np.zeros(grid.ng)
    updated.p_chi = np.zeros(grid.ng)
    updated.chi = np.zeros(grid.ng)
    q0 = float(grid.fine.calibration["b0"] / grid.fine.calibration["a0"])
    updated.Q = np.full(grid.ng, q0)
    info = dict(info)
    info["source_mean_subtracted"] = False
    info["current_array_unchanged"] = True
    return updated, info


def prepare_population(grid, columns_state, occupations):
    """Same mode columns, new occupations, source recomputed, radius and p_Q solved."""
    before_phi0 = np.ascontiguousarray(columns_state.phi0).tobytes()
    before_phi1 = np.ascontiguousarray(columns_state.phi1).tobytes()
    owned = with_occupations(grid, occupations)
    state = galerkin.blank_state(owned, columns_state.phi0, columns_state.phi1)
    solved, radius_info = solve_source_radius(owned, state)
    solved, momentum_info = attach_shift_momentum(owned, solved, radius_info["current"])
    car = car_report(solved)
    columns_unchanged = bool(
        np.ascontiguousarray(solved.phi0).tobytes() == before_phi0
        and np.ascontiguousarray(solved.phi1).tobytes() == before_phi1
        and np.ascontiguousarray(columns_state.phi0).tobytes() == before_phi0
        and np.ascontiguousarray(columns_state.phi1).tobytes() == before_phi1
    )
    r_min = radius_info.get("r_min")
    chart = bool(
        radius_info.get("positive_G")
        and r_min is not None
        and r_min > 0.0
        and float(np.min(solved.Q)) > 0.0
        and car["admissible"]
        and not radius_info["physical_failure"]
    )
    projected = radius_info.get("projected_max")
    full_residual = radius_info.get("full_max")
    evolve = bool(
        chart
        and radius_info.get("source_derived")
        and projected is not None
        and full_residual is not None
        and projected <= EVOLUTION_PROJECTED_MAX
        and full_residual <= EVOLUTION_FULL_MAX
    )
    fine = galerkin.prolong_state(owned, solved)
    shift_residual = shift_constraint(owned.fine, fine) + np.asarray(radius_info["current"], dtype=float)
    shift_mean_gap = abs(float(np.mean(shift_residual)) - float(np.mean(radius_info["current"])))
    report = {
        "occupations": [float(value) for value in owned.fine.occupations],
        "columns_unchanged": columns_unchanged,
        "caller_columns_unchanged": True,
        "car": car,
        "positive_G": bool(radius_info.get("positive_G")),
        "positive_r": bool(r_min is not None and r_min > 0.0),
        "positive_Q": bool(float(np.min(solved.Q)) > 0.0),
        "G_min": radius_info.get("G_min"),
        "G_max": radius_info.get("G_max"),
        "r_min": r_min,
        "r_max": radius_info.get("r_max"),
        "Q": radius_info.get("Q"),
        "projected_hamilton_max": projected,
        "full_hamilton_max": radius_info.get("full_max"),
        "held_out_hamilton_max": radius_info.get("held_out_max"),
        "converged_to_newton_tolerance": bool(radius_info.get("converged")),
        "solver_floor": bool(radius_info.get("solver_floor")),
        "physical_failure": bool(radius_info.get("physical_failure")),
        "source_derived": bool(radius_info["source_derived"]),
        "source_mean_subtracted": False,
        "current_mean": radius_info.get("current_mean"),
        "current_unchanged": bool(radius_info.get("current_unchanged")),
        "shift_mean_gap": shift_mean_gap,
        "mean_retained_in_shift_residual": True,
        "rho_mean": radius_info.get("rho_mean"),
        "rho_independent_of_r_max": radius_info.get("rho_variation"),
        "current_independent_of_r_max": radius_info.get("current_variation"),
        "momentum": {
            "current_mean": momentum_info.get("current_mean"),
            "p_Q_max": momentum_info.get("p_Q_max"),
            "p_Q_mean": momentum_info.get("p_Q_mean"),
            "antiderivative_accepted": momentum_info.get("antiderivative_accepted"),
            "source_mean_subtracted": False,
            "derivative_gap": momentum_info.get("derivative_gap"),
        },
        "chart_admissible": chart,
        "evolve": evolve,
        "evolution_gate": {"projected": EVOLUTION_PROJECTED_MAX, "full": EVOLUTION_FULL_MAX},
        "blocker": radius_info.get("blocker"),
        "accepted_steps": radius_info.get("accepted_steps"),
    }
    if not columns_unchanged:
        report["evolve"] = False
        report["blocker"] = "mode columns changed during preparation"
    return owned, solved, report


def zero_geometry_rate(state, rate):
    """Same column rates. Geometry velocities and the applied fieldwork are zero."""
    zero = np.zeros_like(state.Q)
    return CauchyRate(
        zero, zero, zero, zero, zero, zero,
        rate.phi0, rate.phi1,
        rate.force_L, rate.force_Q, rate.force_beta,
        0.0,
    )


def frozen_rk4_step(grid, state, dt):
    """RK4 of the owned column rates at fixed geometry. No extra force."""
    def stage(current):
        full = galerkin.rates(grid, current, include_matter_force=True)
        return zero_geometry_rate(current, full), full

    k1, full1 = stage(state)
    k2, _full2 = stage(_combine(state, k1, 0.5 * dt))
    k3, _full3 = stage(_combine(state, k2, 0.5 * dt))
    k4, _full4 = stage(_combine(state, k3, dt))

    def avg(one, two, three, four):
        return (one + 2.0 * two + 2.0 * three + four) / 6.0

    combined = CauchyRate(
        avg(k1.Q, k2.Q, k3.Q, k4.Q),
        avg(k1.r, k2.r, k3.r, k4.r),
        avg(k1.chi, k2.chi, k3.chi, k4.chi),
        avg(k1.p_Q, k2.p_Q, k3.p_Q, k4.p_Q),
        avg(k1.p_r, k2.p_r, k3.p_r, k4.p_r),
        avg(k1.p_chi, k2.p_chi, k3.p_chi, k4.p_chi),
        avg(k1.phi0, k2.phi0, k3.phi0, k4.phi0),
        avg(k1.phi1, k2.phi1, k3.phi1, k4.phi1),
        k1.force_L, k1.force_Q, k1.force_beta, 0.0,
    )
    updated = _combine(state, combined, dt)
    return updated, float(full1.fieldwork_power)


def geometry_bytes(state):
    return tuple(
        np.ascontiguousarray(array).tobytes()
        for array in (state.Q, state.r, state.chi, state.p_Q, state.p_r, state.p_chi)
    )


def state_sha256(state):
    import hashlib
    digest = hashlib.sha256()
    for array in (state.Q, state.r, state.chi, state.p_Q, state.p_r, state.p_chi, state.phi0, state.phi1):
        digest.update(np.ascontiguousarray(array).tobytes())
    return digest.hexdigest()


def reversal_amplitude(values):
    """Largest adverse move against the net trend. A monotone series scores zero."""
    values = np.asarray(values, dtype=float)
    if values.size < 2:
        return 0.0
    trend = float(values[-1] - values[0])
    if trend >= 0.0:
        peak = np.maximum.accumulate(values)
        return float(np.max(peak - values))
    trough = np.minimum.accumulate(values)
    return float(np.max(values - trough))


def assess_windows(content, flux, content_effect=CONTENT_EFFECT, share_floor=SHARE_MAINTAINED, packet_share=SHARE_PACKET):
    """Localization, maintenance, and renewal from shell content and proper flux."""
    content = np.asarray(content, dtype=float)
    flux = np.asarray(flux, dtype=float)
    if content.ndim != 2 or content.shape[1] != 4 or flux.shape != content.shape:
        raise ValueError("content and flux must be samples by four windows")
    totals = np.sum(content, axis=1)
    if np.any(totals == 0.0):
        raise ValueError("shell content vanished")
    shares = content / totals[:, None]
    leaders = np.argmax(shares, axis=1)
    start_leader = int(leaders[0])
    end_leader = int(leaders[-1])
    leader_content = content[:, start_leader]
    reversal = reversal_amplitude(leader_content)
    end_share = float(shares[-1, end_leader])
    localized_end = bool(end_share >= share_floor)
    same_leader = bool(np.all(leaders == start_leader))
    maintained = bool(
        same_leader
        and localized_end
        and float(leader_content[-1]) > 0.5 * float(leader_content[0])
    )
    rebuilt = bool(
        reversal >= content_effect
        and localized_end
        and float(shares[-1, start_leader]) >= share_floor
    )
    gained = content[-1] - content[0]
    new_leader = bool(
        end_leader != start_leader
        and localized_end
        and float(gained[end_leader]) >= content_effect
    )
    packet = shares >= packet_share
    packet_flux = np.sum(np.where(packet, flux, 0.0), axis=1)
    reservoir_flux = np.sum(np.where(packet, 0.0, flux), axis=1)
    renewed = bool(rebuilt or new_leader)
    content_change = float(leader_content[-1] - leader_content[0])
    stable_throughflow = bool(
        localized_end
        and same_leader
        and abs(content_change) < content_effect
        and abs(float(packet_flux[-1])) >= content_effect
    )
    packet_end = float(packet_flux[-1])
    reservoir_end = float(reservoir_flux[-1])
    flux_sum = float(packet_end + reservoir_end)
    flux_scale = max(1.0, abs(packet_end), abs(reservoir_end))
    reports = {
        "content": {
            "maintained": maintained,
            "localized_end": localized_end,
            "same_leader": same_leader,
            "start_leader": start_leader,
            "end_leader": end_leader,
            "end_share": end_share,
            "leader_content_change": content_change,
        },
        "flux": {
            "packet_flux_end": packet_end,
            "reservoir_flux_end": reservoir_end,
            "flux_sum_end": flux_sum,
            "balanced": bool(abs(flux_sum) <= 1e-9 * flux_scale),
            "throughflow_magnitude": abs(packet_end),
        },
        "reversal": {
            "reversal": reversal,
            "one_drift": bool(not renewed and reversal < content_effect),
            "rebuilt": rebuilt,
            "new_leader": new_leader,
            "renewed": renewed,
            "splitting_one_drift_is_complete_renewal": False,
        },
        "accounting": {
            "reported_with_windows": False,
            "role": "Energy exchange stays with the caller balance.",
        },
    }
    return {
        "start_leader": start_leader,
        "end_leader": end_leader,
        "same_leader": same_leader,
        "end_share": end_share,
        "start_share": float(shares[0, start_leader]),
        "leader_content_start": float(leader_content[0]),
        "leader_content_end": float(leader_content[-1]),
        "leader_content_change": content_change,
        "reversal": reversal,
        "localized_end": localized_end,
        "maintained": maintained,
        "rebuilt": rebuilt,
        "new_leader": new_leader,
        "renewed": renewed,
        "one_drift": bool(not renewed and reversal < content_effect),
        "stable_throughflow": stable_throughflow,
        "candidate_regime": bool(renewed or stable_throughflow),
        "candidate_regime_is_programme_requirement": False,
        "packet_flux_end": packet_end,
        "reservoir_flux_end": reservoir_end,
        "flux_sum_end": flux_sum,
        "content_effect": float(content_effect),
        "splitting_one_drift_is_renewal": False,
        "splitting_one_drift_is_complete_renewal": False,
        "maintained_with_throughflow": bool(maintained and abs(packet_end) >= float(content_effect)),
        "reports": reports,
    }


def separate_ledgers(content, flux, accounting=None, content_effect=CONTENT_EFFECT, share_floor=SHARE_MAINTAINED, packet_share=SHARE_PACKET):
    """Content, flux, reversal, and accounting as four readings.

    ``candidate_regime`` stays available on the window proxy and is not a
    programme requirement. A maintained structure may carry throughflow
    while its shell content is still drifting.
    """
    windows = assess_windows(
        content, flux,
        content_effect=content_effect, share_floor=share_floor, packet_share=packet_share,
    )
    return {
        "content": windows["reports"]["content"],
        "flux": windows["reports"]["flux"],
        "reversal": windows["reports"]["reversal"],
        "accounting": dict(accounting) if accounting is not None else dict(windows["reports"]["accounting"]),
        "proxies": {
            "candidate_regime": windows["candidate_regime"],
            "stable_throughflow": windows["stable_throughflow"],
            "renewed": windows["renewed"],
            "reversal_floor": float(content_effect),
            "share_floor": float(share_floor),
            "programme_requirements": False,
        },
        "maintained_with_throughflow": windows["maintained_with_throughflow"],
        "splitting_one_drift_is_complete_renewal": False,
    }


def dynamic_stop_reason(failure, gram, *, delocalized=False):
    """Dynamic exit from chart_failure or the column Gram.

    G is not an argument. A negative bracket does not stop the step.
    """
    if failure:
        return str(failure), True
    if gram is not None and not bool(gram.get("admissible")):
        return "GRAM_LEFT_IDENTITY", False
    if delocalized:
        return "NONLOCALIZATION", False
    return None, False


def uniform_complement_localization(occupations=None):
    """Six positive weights against an empty quadrature complement."""
    values = np.asarray(UNIFORM_OCCUPATIONS if occupations is None else occupations, dtype=float)
    return {
        "weights": [float(value) for value in values],
        "support": int(values.size),
        "complement_occupation": 0.0,
        "localized_relative_to_empty_complement": bool(values.size == 6 and np.all(values > 0.0)),
        "window_share_proxy_is_a_different_question": True,
    }


def galerkin_half_identity_report(grid, state=None):
    """H_G on 2 n_f against the unprojected fine embedding.

    (1/2) I_{2 n_f} commutes with H_G = U_f† H_fine U_f. The same identity
    pushed to 2 n_q is a different operator and need not commute with H_fine.
    """
    from recursive_horizons.nsc_conformal_adm_source import direct_hamiltonian
    from recursive_horizons.nsc_covariant_operator import CovariantStaticMetric
    from recursive_horizons.nsc_spherical_coupling import apply_dirac

    if state is None:
        phi0, phi1, _preparation, problems = galerkin.load_physical_columns(grid.nf)
        if problems:
            raise ValueError("physical columns are unresolved: " + ", ".join(problems))
        state = galerkin.blank_state(grid, phi0, phi1)
    fine = galerkin.prolong_state(grid, state)
    system = grid.fine
    metric = CovariantStaticMetric(
        system.length,
        system.length_density * fine.r,
        fine.Q * fine.r,
        fine.r,
        eta=0.5,
    )
    hamiltonian = direct_hamiltonian(metric, system.shift, system.kappa)
    points = int(system.points)
    fermions = int(grid.nf)
    isometry = grid.U_f
    embedding = np.zeros((2 * points, 2 * fermions), dtype=complex)
    embedding[:points, :fermions] = isometry
    embedding[points:, fermions:] = isometry
    projected = embedding.conj().T @ hamiltonian @ embedding
    half = 0.5 * np.eye(2 * fermions, dtype=complex)
    galerkin_commutator = projected @ half - half @ projected
    embedded = embedding @ half @ embedding.conj().T
    fine_commutator = hamiltonian @ embedded - embedded @ hamiltonian
    eye = np.eye(fermions, dtype=complex)
    zero = np.zeros((fermions, fermions), dtype=complex)
    band0 = np.concatenate((eye, zero), axis=1)
    band1 = np.concatenate((zero, eye), axis=1)
    image0, image1 = apply_dirac(
        isometry @ band0, isometry @ band1,
        system.length_density, fine.Q, system.shift, system.kappa, system.momentum,
    )
    pulled = np.vstack((isometry.conj().T @ image0, isometry.conj().T @ image1))
    galerkin_norm = float(np.linalg.norm(galerkin_commutator, ord="fro"))
    embedding_norm = float(np.linalg.norm(fine_commutator, ord="fro"))
    return {
        "nf": fermions,
        "nq": points,
        "dimension_H_G": int(projected.shape[0]),
        "dimension_H_fine": int(hamiltonian.shape[0]),
        "pullback_gap": float(np.max(np.abs(pulled - projected))),
        "half_identity_commutator_H_G": galerkin_norm,
        "half_identity_commutes_with_H_G": bool(galerkin_norm <= 1e-8),
        "embedding_commutator_H_fine": embedding_norm,
        "embedding_commutes_with_unprojected_H_fine": bool(embedding_norm <= 1e-8),
        "fine_identity_is_algebraic_control": FINE_IDENTITY_IS_ALGEBRAIC_CONTROL,
        "fine_identity_is_an_occupation": FINE_IDENTITY_IS_AN_OCCUPATION,
    }


def continuation_check(grid, state, dt, *, frozen=False):
    """One owned RK4 step with Gram and chart_failure. G is diagnostic only."""
    before_q = np.array(state.Q, dtype=float, copy=True)
    try:
        if frozen:
            candidate, refused = frozen_rk4_step(grid, state, dt)
        else:
            candidate = galerkin.rk4_step(grid, state, dt)
            refused = 0.0
    except PositiveChartExit as exit_chart:
        return {
            "accepted": False,
            "chart_stop": True,
            "stop_reason": str(exit_chart.reason),
            "state": state,
            "Q_change": 0.0,
            "gram": car_report(state),
            "positive_G": None,
            "G_min": None,
            "kept_failed_step": False,
        }
    fine = galerkin.prolong_state(grid, candidate)
    failure = chart_failure(grid.fine, fine)
    gram = car_report(candidate)
    q_change = float(np.max(np.abs(np.asarray(candidate.Q, dtype=float) - before_q)))
    if failure:
        reason, chart_stop = dynamic_stop_reason(failure, gram)
        return {
            "accepted": False,
            "chart_stop": chart_stop,
            "stop_reason": reason,
            "state": state,
            "Q_change": q_change,
            "gram": gram,
            "positive_G": None,
            "G_min": None,
            "kept_failed_step": False,
        }
    _source, rho, _current = source_arrays(grid, candidate)
    gee = bracket_g(fine.Q, rho, grid)
    reason, chart_stop = dynamic_stop_reason(None, gram)
    accepted = reason is None
    return {
        "accepted": accepted,
        "chart_stop": chart_stop,
        "stop_reason": reason,
        "state": candidate if accepted or reason == "GRAM_LEFT_IDENTITY" else state,
        "Q_change": q_change,
        "gram": gram,
        "positive_G": bool(float(np.min(gee)) > 0.0),
        "G_min": float(np.min(gee)),
        "positive_r": bool(float(np.min(fine.r)) > 0.0),
        "positive_Q": bool(float(np.min(fine.Q)) > 0.0),
        "positive_L": bool(float(np.min(grid.fine.length_density)) > 0.0),
        "kept_failed_step": False,
    }


def chart_from_fields(grid, fields):
    """Bracket diagnostic plus the owned positive chart.

    ``positive_G`` records the convex-radius bracket. It is not a dynamic stop.
    """
    gee = bracket_g(fields["Q"], fields["rho"], grid)
    lapse = np.asarray(fields["r"], dtype=float) * np.asarray(grid.fine.length_density, dtype=float)
    return {
        "G_min": float(np.min(gee)),
        "G_max": float(np.max(gee)),
        "positive_G": bool(np.min(gee) > 0.0),
        "positive_r": bool(np.min(fields["r"]) > 0.0),
        "positive_Q": bool(np.min(fields["Q"]) > 0.0),
        "positive_L": bool(np.min(grid.fine.length_density) > 0.0),
        "positive_lapse": bool(np.min(lapse) > 0.0),
        "g_is_dynamic_stop": False,
    }


def shell_views(sample):
    content = np.asarray(sample["window_normal"], dtype=float)
    flux = np.asarray(sample["window_proper_flux"], dtype=float)
    total = float(np.sum(content))
    shares = content / total
    leader = int(np.argmax(shares))
    packet = shares >= SHARE_PACKET
    return {
        "content": content,
        "flux": flux,
        "shares": shares,
        "leader": leader,
        "leader_share": float(shares[leader]),
        "packet_flux": float(np.sum(flux[packet])),
        "reservoir_flux": float(np.sum(flux[~packet])),
    }


def mode_energies(sample, leader):
    totals = np.zeros(6, dtype=float)
    leader_energy = np.zeros(6, dtype=float)
    occupations = np.zeros(6, dtype=float)
    for index, row in enumerate(sample["windows"]):
        for column_index, column in enumerate(row["mode_columns"]):
            totals[column_index] += float(column["hamiltonian_energy"])
            occupations[column_index] = float(column["occupation"])
            if index == leader:
                leader_energy[column_index] = float(column["hamiltonian_energy"])
    return occupations, totals, leader_energy


def proper_clock_rates(grid, fields, windows, leader):
    lapse = np.asarray(fields["r"], dtype=float) * np.asarray(grid.fine.length_density, dtype=float)
    window = np.asarray(windows[leader], dtype=float)
    mass = float(np.sum(window))
    leader_rate = float(np.sum(window * lapse) / mass)
    return float(np.mean(lapse)), leader_rate


def load_episode_final(npz_path, run_name):
    """Bitwise Cauchy state from a saved episode final, without a new solve."""
    with np.load(npz_path, allow_pickle=False) as data:
        state = CauchyState(
            Q=np.array(data[run_name + "_final_Q"], dtype=float, copy=True),
            r=np.array(data[run_name + "_final_r"], dtype=float, copy=True),
            chi=np.array(data[run_name + "_final_chi"], dtype=float, copy=True),
            p_Q=np.array(data[run_name + "_final_p_Q"], dtype=float, copy=True),
            p_r=np.array(data[run_name + "_final_p_r"], dtype=float, copy=True),
            p_chi=np.array(data[run_name + "_final_p_chi"], dtype=float, copy=True),
            phi0=np.array(data[run_name + "_final_phi0"], copy=True),
            phi1=np.array(data[run_name + "_final_phi1"], copy=True),
        )
    return state


def imbalance_perturbation(occupations, fraction):
    """Scale the deviation from one half. fraction is +0.05 or -0.05."""
    base = np.asarray(occupations, dtype=float)
    return 0.5 + (1.0 + float(fraction)) * (base - 0.5)


def join_series(first_time, first_values, second_time, second_values, join=JOIN_TIME):
    first_time = np.asarray(first_time, dtype=float)
    second_time = np.asarray(second_time, dtype=float)
    first_values = np.asarray(first_values, dtype=float)
    second_values = np.asarray(second_values, dtype=float)
    keep = first_time < float(join) - 1e-12
    return np.concatenate([first_time[keep], second_time]), np.concatenate(
        [first_values[keep], second_values], axis=0
    )
