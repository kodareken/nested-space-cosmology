#!/usr/bin/env python3
"""Actual metric curvature of the sealed conformal v2 Cauchy trajectory.

Qddot is the exact directional Jacobian of the actual projected bilinear
Qrate, evaluated along the actual projected state rate. Both projections and
Ldot=Qdot are retained. Rh is built from the metric jets; chi+2 is compared
afterward and never substituted into that calculation.
"""
from __future__ import annotations
import json
import os
from pathlib import Path
import time

for _name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_name] = "1"

import numpy as np
import derive_nsc_spherical_conformal_episode as episode
import derive_nsc_spherical_conformal_continuation as continuation
from recursive_horizons import nsc_spherical_galerkin_coupling as galerkin
from recursive_horizons import nsc_spherical_episode_assessment as metric
from recursive_horizons.nsc_spherical_feedback_action import partial_F
from recursive_horizons.nsc_spherical_coupling import CauchyState

LAB = Path(__file__).resolve().parents[1]
OUT = LAB / "results/development/nsc-spherical-conformal-curvature-v1.json"
NPZ = LAB / "results/development/nsc-spherical-conformal-curvature-v1.npz"
DRIVER = Path(__file__).resolve()


def actual_q_second_rate(grid, state, rate=None, bundle=None):
    if grid.gauge != "conformal":
        raise ValueError("the bilinear Qrate identity requires conformal gauge")
    if rate is None or bundle is None:
        rate, bundle = galerkin.compose_fine_hamiltonian(grid, state)
    fine = bundle["fine_state"]
    _fr, f_chi = partial_F(fine.r, grid.fine.A, grid.fine.C_W)
    f_chi = float(f_chi)
    q_dot = galerkin.prolong_geometry(grid, rate.Q)
    p_chi_dot = galerkin.prolong_geometry(grid, rate.p_chi)
    directional = (q_dot * fine.p_chi + fine.Q * p_chi_dot) / (2 * f_chi)
    q_ddot = galerkin.prolong_geometry(grid, galerkin.pull_geometry(grid, directional))
    return q_dot, q_ddot, directional


def christoffel_crosscheck(Q, Q_t, Q_tt, Q_x, Q_tx, Q_xx, nodes):
    errors = []
    for node in nodes:
        q, qt, qtt, qx, qtx, qxx = (float(values[node]) for values in (Q, Q_t, Q_tt, Q_x, Q_tx, Q_xx))
        signature = np.array((1., -1.))
        g = np.diag(signature * q ** 2)
        dg = np.zeros((2, 2, 2))
        ddg = np.zeros((2, 2, 2, 2))
        for i in range(2):
            dg[i, i, 0] = signature[i] * 2 * q * qt
            dg[i, i, 1] = signature[i] * 2 * q * qx
            ddg[0, i, i, 0] = signature[i] * (2 * qt ** 2 + 2 * q * qtt)
            ddg[0, i, i, 1] = ddg[1, i, i, 0] = signature[i] * (2 * qt * qx + 2 * q * qtx)
            ddg[1, i, i, 1] = signature[i] * (2 * qx ** 2 + 2 * q * qxx)
        actual = metric.ricci_scalar_from_christoffel(g, dg, ddg)
        explicit = 2 / q ** 2 * (qxx / q - qx ** 2 / q ** 2 - qtt / q + qt ** 2 / q ** 2)
        errors.append(abs(actual - explicit))
    return max(errors)


def evaluate(grid, state):
    rate, bundle = galerkin.compose_fine_hamiltonian(grid, state)
    fine = bundle["fine_state"]
    qt, qtt, before_final_projection = actual_q_second_rate(grid, state, rate, bundle)
    zero = np.zeros(grid.nq)
    curvature = metric.direct_rh_grid(fine.Q, fine.Q, zero, qt, qtt, grid.length, L_dot=qt, beta_dot=zero)["R_h"]
    actual_weyl = metric.weyl_actual(curvature, fine.r)
    proxy = metric.weyl_proxy_from_chi(fine.chi, fine.r)
    qx = grid.derivative @ fine.Q
    qxx = grid.derivative @ qx
    qtx = grid.derivative @ qt
    formula = 2 / fine.Q ** 2 * (qxx / fine.Q - qx ** 2 / fine.Q ** 2 - qtt / fine.Q + qt ** 2 / fine.Q ** 2)
    nodes = (0, grid.nq // 8, grid.nq // 4, 3 * grid.nq // 8, grid.nq // 2)
    report = {"R_h_min": float(np.min(curvature)), "R_h_max": float(np.max(curvature)),
              "R_h_l2": float(np.sqrt(grid.dx_q * np.sum(curvature ** 2))),
              "weyl_C2_max": float(np.max(actual_weyl)), "weyl_C2_l2": float(np.sqrt(grid.dx_q * np.sum(actual_weyl ** 2))),
              "chi_proxy_max": float(np.max(proxy)), "R_h_minus_chi2_max": float(np.max(abs(curvature - fine.chi - 2))),
              "actual_minus_proxy_C2_max": float(np.max(abs(actual_weyl - proxy))),
              "direct_formula_gap": float(np.max(abs(curvature - formula))),
              "christoffel_contraction_gap": christoffel_crosscheck(fine.Q, qt, qtt, qx, qtx, qxx, nodes),
              "Qddot_projection_gap_max": float(np.max(abs(qtt - before_final_projection))),
              "r_min": float(np.min(fine.r)), "Q_min": float(np.min(fine.Q)),
              "positive_chart": bool(np.min(fine.r) > 0 and np.min(fine.Q) > 0)}
    profiles = {"R_h": curvature, "weyl_C2": actual_weyl, "Q_dot": qt, "Q_ddot": qtt, "r": fine.r, "chi": fine.chi}
    return report, profiles


def run():
    start_cpu = time.process_time()
    bindings = {"episode_v2": episode.sha256(continuation.OUT), "payload_v2": episode.sha256(continuation.NPZ),
                "galerkin": episode.sha256(galerkin._MODULE_PATH), "metric_helper": episode.sha256(Path(metric.__file__)),
                "action": episode.sha256(galerkin._ACTION_PATH), "driver": episode.sha256(DRIVER)}
    record = {"schema": "NSC-SPHERICAL-CONFORMAL-CURVATURE-v1", "gauge": "conformal", "source_bindings": bindings,
              "geometry_evolved": False, "Cauchy_state_changed": False,
              "Qddot": "A_g pull_g[(Qdot_lift*pchi_fine+Q_fine*pchi_dot_lift)/(2Fchi)], exact actual-ODE Qrate Jacobian",
              "lapse_rate": "Ldot=Qdot_lift", "shift_rate": "betadot=0", "auxiliary_substituted": False,
              "R_h": "metric jets of h=Q^2 diag(1,-1), with actual projected Qdot/Qddot",
              "weyl_C2": "(R_h-2)^2/(3r^4), owned conformal spherical invariant identity",
              "continuum_curvature_certified": False, "series": {}, "results": {}}
    arrays = {}
    with np.load(continuation.NPZ, allow_pickle=False) as stored:
        for nf in (256, 512):
            grid = galerkin.build_grid(nf, quadrature=4 * nf, gauge="conformal")
            for spec in (row for row in episode.RUN_PLAN if row["nf"] == nf):
                times = stored[spec["name"] + "_frame_times"]
                states = {name: stored[spec["name"] + "_frames_" + name] for name in episode.STATE_NAMES}
                rows = []
                for index, t in enumerate(times):
                    state = CauchyState(**{name: values[index] for name, values in states.items()})
                    row, profiles = evaluate(grid, state)
                    row["time"] = float(t)
                    rows.append(row)
                    if index in (0, len(times) - 1):
                        stage = "initial" if index == 0 else "final"
                        for name, value in profiles.items():
                            arrays[spec["name"] + "_" + stage + "_" + name] = value
                record["series"][spec["name"]] = rows
                record["results"][spec["name"]] = {"initial": rows[0], "final": rows[-1],
                                                    "maximum_auxiliary_gap": max(row["R_h_minus_chi2_max"] for row in rows),
                                                    "maximum_christoffel_gap": max(row["christoffel_contraction_gap"] for row in rows)}
                print("curvature", spec["name"], rows[-1], flush=True)
    record["cpu_seconds"] = time.process_time() - start_cpu
    record["source_unchanged"] = bindings["episode_v2"] == episode.sha256(continuation.OUT) and bindings["payload_v2"] == episode.sha256(continuation.NPZ)
    with NPZ.open("wb") as handle:
        np.savez_compressed(handle, **arrays)
    record["payload_sha256"] = episode.sha256(NPZ)
    OUT.write_text(json.dumps(episode.jsonable(record), indent=2, allow_nan=False) + "\n")
    record["payload_bytes"] = OUT.stat().st_size + NPZ.stat().st_size
    OUT.write_text(json.dumps(episode.jsonable(record), indent=2, allow_nan=False) + "\n")
    print("saved", OUT, "CPU", record["cpu_seconds"], "bytes", record["payload_bytes"], flush=True)
    return record


if __name__ == "__main__":
    run()
