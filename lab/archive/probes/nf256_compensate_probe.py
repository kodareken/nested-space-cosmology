#!/usr/bin/env python3
"""Diagnose cancellation and compensated residual at the nf256 full-source seed."""
import math
import os

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"

import numpy as np

from recursive_horizons import nsc_spherical_galerkin_coupling as galerkin
from recursive_horizons.nsc_spherical_feedback_action import feedback_F, feedback_V, feedback_Z

print("longdouble_eq_float64", np.dtype(np.longdouble) == np.dtype(np.float64))


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
    return np.real(synthesis @ coeff_f)


def fsum_matvec(matrix, vector):
    products = np.asarray(matrix, dtype=np.float64) * np.asarray(vector, dtype=np.float64)
    out = np.empty(matrix.shape[0], dtype=np.float64)
    for index in range(matrix.shape[0]):
        out[index] = math.fsum(products[index])
    return out


def split_residual(system, radius_fine, rho):
    q_value = system.calibration["b0"] / system.calibration["a0"]
    q = np.full(system.points, q_value)
    derivative = fsum_matvec(system.derivative, radius_fine)
    derivative_blas = system.derivative @ radius_fine
    force = feedback_F(radius_fine, np.zeros(system.points), system.A, system.C_W)
    f_over_q = fsum_matvec(system.derivative, force) / q
    potential = feedback_V(radius_fine, np.zeros(system.points), system.A, system.C_W, system.C_F, system.flux)
    z_value = float(feedback_Z(system.A))
    term_d = z_value * derivative ** 2 / q
    term_v = -q * potential
    term_f = -2 * fsum_matvec(system.derivative, f_over_q)
    compensated = term_d + term_v + term_f + rho
    blas = galerkin._radius_residual(system, radius_fine, rho)
    # Same formula with BLAS matvecs, to prove the split matches the module.
    f_over_q_blas = (system.derivative @ force) / q
    mine = (
        z_value * derivative_blas ** 2 / q
        - q * potential
        - 2 * (system.derivative @ f_over_q_blas)
        + rho
    )
    return {
        "term_d_max": float(np.max(np.abs(term_d))),
        "term_v_max": float(np.max(np.abs(term_v))),
        "term_f_max": float(np.max(np.abs(term_f))),
        "rho_max": float(np.max(np.abs(rho))),
        "compensated_max": float(np.max(np.abs(compensated))),
        "blas_max": float(np.max(np.abs(blas))),
        "formula_match_max": float(np.max(np.abs(mine - blas))),
        "compensated_minus_blas_max": float(np.max(np.abs(compensated - blas))),
        "derivative_minus_blas_max": float(np.max(np.abs(derivative - derivative_blas))),
        "compensated": compensated,
        "blas": blas,
    }


def project(grid, fine_residual):
    blas = galerkin.pull_geometry(grid, fine_residual)
    compensated = grid.weight * fsum_matvec(grid.A_g.T, fine_residual)
    return blas, compensated


def main():
    loaded = galerkin.load_physical_columns(128)
    grid128 = galerkin.build_grid(128, quadrature=512)
    solved128, info128 = galerkin.solve_initial_radius(grid128, galerkin.blank_state(grid128, loaded[0], loaded[1]))
    print("nf128", info128["converged"], info128["residual_max"])
    loaded256 = galerkin.load_physical_columns(256)
    grid = galerkin.build_grid(256, quadrature=1024)
    state = galerkin.blank_state(grid, loaded256[0], loaded256[1])
    rho = galerkin.source_from_columns(grid.fine, galerkin.prolong_state(grid, state))["force_L"] / grid.dx_q
    radius = embed_radius(solved128.r, 127, 255, grid.length)
    best = None
    for iteration in range(6):
        projected, _full = galerkin.projected_radius_operator(grid, radius, rho)
        residual_max = float(np.max(np.abs(projected)))
        print("iter", iteration, residual_max)
        if best is None or residual_max < best[0]:
            best = (residual_max, radius.copy())
        if residual_max < galerkin.TOL_NEWTON:
            break
        jacobian = galerkin._projected_radius_jacobian(grid, radius)
        delta = np.linalg.solve(jacobian, -projected)
        print("delta", float(np.max(np.abs(delta))))
        trial = radius + delta
        if np.min(galerkin.prolong_geometry(grid, trial)) <= 0:
            print("left chart")
            break
        trial_projected, _ = galerkin.projected_radius_operator(grid, trial, rho)
        trial_max = float(np.max(np.abs(trial_projected)))
        print("trial", trial_max)
        if trial_max < residual_max:
            radius = trial
        else:
            print("no decrease")
            break
        if float(np.max(np.abs(delta))) < 1e-12:
            break
    radius = best[1]
    print("best_blas_projected", best[0])
    fine = galerkin.prolong_geometry(grid, radius)
    report = split_residual(grid.fine, fine, rho)
    blas_p, comp_p = project(grid, report["blas"])
    _blas_from_comp, comp_of_comp = project(grid, report["compensated"])
    print("formula_match", report["formula_match_max"])
    print("term_d", report["term_d_max"], "term_v", report["term_v_max"], "term_f", report["term_f_max"], "rho", report["rho_max"])
    print("blas_full", report["blas_max"], "comp_full", report["compensated_max"])
    print("comp_minus_blas_full", report["compensated_minus_blas_max"], "derivative_gap", report["derivative_minus_blas_max"])
    print("blas_projected", float(np.max(np.abs(blas_p))), "comp_projected_of_blas_fine", float(np.max(np.abs(comp_p))))
    print("comp_projected", float(np.max(np.abs(comp_of_comp))))
    # One iterative refinement using the compensated projected residual and the same Jacobian.
    jacobian = galerkin._projected_radius_jacobian(grid, radius)
    comp_residual = comp_of_comp
    delta = np.linalg.solve(jacobian, -comp_residual)
    # Compensated linear residual of that step.
    linear = fsum_matvec(jacobian, delta) + comp_residual
    correction = np.linalg.solve(jacobian, -linear)
    print("polish_delta", float(np.max(np.abs(delta))), "linear_comp", float(np.max(np.abs(linear))), "correction", float(np.max(np.abs(correction))))
    polished = radius + delta + correction
    fine2 = galerkin.prolong_geometry(grid, polished)
    report2 = split_residual(grid.fine, fine2, rho)
    print("polished_blas_full", report2["blas_max"], "polished_comp_full", report2["compensated_max"])
    blas2, _ = project(grid, report2["blas"])
    _b, comp2 = project(grid, report2["compensated"])
    print("polished_blas_projected", float(np.max(np.abs(blas2))), "polished_comp_projected", float(np.max(np.abs(comp2))))
    state.r = radius
    diagnostics = galerkin.constraint_diagnostics(grid, state)
    print("diagnostics", {key: diagnostics[key] for key in (
        "projected_hamilton_max", "projected_momentum_max", "full_hamilton_max", "full_momentum_max",
        "held_out_hamilton_max", "held_out_momentum_max", "positive_r", "positive_Q", "r_min_quadrature",
    )})
    print("unresolved", diagnostics["hamilton_modes"]["unresolved_fraction"], "modes_max", diagnostics["hamilton_modes"]["unresolved_max_mode"])


if __name__ == "__main__":
    main()
