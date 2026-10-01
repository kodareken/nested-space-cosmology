#!/usr/bin/env python3
"""Probe: prolong nf128 radius as the full-source nf256 Newton seed."""
import os
import time

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"

import numpy as np

from recursive_horizons import nsc_spherical_galerkin_coupling as galerkin

print("longdouble_is_float64", np.dtype(np.longdouble) == np.dtype(np.float64), np.finfo(np.longdouble).eps)


def embed_radius(radius, ng_coarse, ng_fine, length):
    modes_c = galerkin.geometry_modes(ng_coarse)
    modes_f = galerkin.geometry_modes(ng_fine)
    coarse = np.arange(ng_coarse, dtype=float) * (length / ng_coarse)
    fine = np.arange(ng_fine, dtype=float) * (length / ng_fine)
    analysis = np.exp(-2j * np.pi * coarse[:, None] * modes_c[None, :] / length) / ng_coarse
    coefficients = analysis.T @ np.asarray(radius, dtype=float)
    coeff_f = np.zeros(ng_fine, dtype=complex)
    lookup = {float(mode): index for index, mode in enumerate(modes_f)}
    for coefficient, mode in zip(coefficients, modes_c):
        coeff_f[lookup[float(mode)]] = coefficient
    synthesis = np.exp(2j * np.pi * fine[:, None] * modes_f[None, :] / length)
    nodal = synthesis @ coeff_f
    return np.real(nodal), float(np.max(np.abs(np.imag(nodal)))), float(np.max(np.abs(coeff_f[np.abs(modes_f) > ng_coarse // 2])))


def newton_full(grid, radius, rho, rounds=8):
    history = []
    for iteration in range(rounds):
        projected, full = galerkin.projected_radius_operator(grid, radius, rho)
        residual_max = float(np.max(np.abs(projected)))
        full_max = float(np.max(np.abs(full)))
        history.append({"iteration": iteration, "projected": residual_max, "full_operator": full_max})
        print("newton", history[-1], flush=True)
        if residual_max < galerkin.TOL_NEWTON:
            return radius, history, True
        jacobian = galerkin._projected_radius_jacobian(grid, radius)
        delta = np.linalg.solve(jacobian, -projected)
        linear = jacobian @ delta + projected
        print(
            "linear_residual",
            float(np.max(np.abs(linear))),
            "delta_max",
            float(np.max(np.abs(delta))),
            "cond_est",
            float(np.linalg.cond(jacobian)),
            flush=True,
        )
        accepted = False
        step = 1.0
        for _halving in range(16):
            trial = radius + step * delta
            if not np.isfinite(trial).all() or np.min(galerkin.prolong_geometry(grid, trial)) <= 0:
                step *= 0.5
                continue
            trial_projected, _ = galerkin.projected_radius_operator(grid, trial, rho)
            trial_max = float(np.max(np.abs(trial_projected)))
            if trial_max < residual_max * (1 - 0.05 * step):
                radius = trial
                accepted = True
                print("accepted", step, trial_max, flush=True)
                break
            step *= 0.5
        if not accepted:
            print("stalled", residual_max, flush=True)
            return radius, history, False
    projected, _full = galerkin.projected_radius_operator(grid, radius, rho)
    residual_max = float(np.max(np.abs(projected)))
    history.append({"iteration": rounds, "projected": residual_max})
    return radius, history, residual_max < galerkin.TOL_NEWTON


def main():
    wall = time.perf_counter()
    cpu = time.process_time()
    loaded = galerkin.load_physical_columns(128)
    print("prep128", loaded[3], flush=True)
    grid128 = galerkin.build_grid(128, quadrature=512)
    state128 = galerkin.blank_state(grid128, loaded[0], loaded[1])
    solved128, info128 = galerkin.solve_initial_radius(grid128, state128)
    print(
        "nf128",
        info128["converged"],
        info128["residual_max"],
        info128["r_min_coarse"],
        info128["r_max_coarse"],
        "cpu",
        time.process_time() - cpu,
        flush=True,
    )
    loaded256 = galerkin.load_physical_columns(256)
    print("prep256", loaded256[3], flush=True)
    grid = galerkin.build_grid(256, quadrature=1024)
    state = galerkin.blank_state(grid, loaded256[0], loaded256[1])
    fine_state = galerkin.prolong_state(grid, state)
    source = galerkin.source_from_columns(grid.fine, fine_state)
    rho = source["force_L"] / grid.dx_q
    print("rho_max", float(np.max(np.abs(rho))), "nq", grid.nq, "ng", grid.ng, flush=True)
    radius, imag, high = embed_radius(solved128.r, 127, 255, grid.length)
    print("embed_imag", imag, "high_mode_coeff", high, "r_min", float(np.min(radius)), "r_max", float(np.max(radius)), flush=True)
    seed_projected, seed_full = galerkin.projected_radius_operator(grid, radius, rho)
    print(
        "seed_projected",
        float(np.max(np.abs(seed_projected))),
        "seed_full_operator",
        float(np.max(np.abs(seed_full))),
        flush=True,
    )
    state.r = radius
    diagnostics = galerkin.constraint_diagnostics(grid, state)
    print(
        "seed_full_hamilton",
        diagnostics["full_hamilton_max"],
        "seed_projected_hamilton",
        diagnostics["projected_hamilton_max"],
        "positive",
        diagnostics["positive_r"],
        diagnostics["positive_Q"],
        "r_min",
        diagnostics["r_min_quadrature"],
        flush=True,
    )
    radius, history, converged = newton_full(grid, radius, rho)
    state.r = radius
    diagnostics = galerkin.constraint_diagnostics(grid, state)
    print(
        "after_converged_flag",
        converged,
        "full_hamilton",
        diagnostics["full_hamilton_max"],
        "projected_hamilton",
        diagnostics["projected_hamilton_max"],
        "full_momentum",
        diagnostics["full_momentum_max"],
        "held_h",
        diagnostics["held_out_hamilton_max"],
        "positive",
        diagnostics["positive_r"],
        diagnostics["r_min"],
        diagnostics["r_min_quadrature"],
        flush=True,
    )
    print("wall", time.perf_counter() - wall, "cpu", time.process_time() - cpu, flush=True)


if __name__ == "__main__":
    main()
