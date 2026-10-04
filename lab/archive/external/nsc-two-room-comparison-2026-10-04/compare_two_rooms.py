#!/usr/bin/env python3
"""Read-only NSC benchmark calculation: ray transfer and child Dirac modes.

Reuses the fixed Bronnikov--Fabris geometry recorded in the NSC repository at
7656eec34141ceca22c073bca997acefb317be05. This is NOT a rerun of the repository's
self-consistent state--geometry campaigns and does not alter the repository.

Run: python compare_two_rooms.py --output result.json
Dependencies: numpy, scipy, mpmath. Existing output files are never replaced.
Units: c = hbar = L_throat = 1. The geometry also sets the mass parameter m=1.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
from pathlib import Path

import mpmath as mp
import numpy as np
import scipy
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq

COMMIT = "7656eec34141ceca22c073bca997acefb317be05"
REPO = "https://github.com/kodareken/nested-space-cosmology"
REFERENCE_PATHS = (
    "docs/nsc-clock-horizon.md", "docs/nsc-dirac-tetrad.md",
    "results/nsc-5-clock-horizon.json",
)
SIGMA1 = np.array([[0, 1], [1, 0]], dtype=complex)
SIGMA3 = np.diag([1., -1.]).astype(complex)


def geometry(rho: float) -> tuple[float, float, float, float]:
    """Return A, A', beta, r. atan2 avoids subtracting pi/2 for rho>0."""
    x = math.atan2(1.0, float(rho))
    radius2 = 1.0 + rho*rho
    A = 1.0 + 3.0*rho - 3.0*radius2*x
    Ap = 6.0*(1.0-rho*x)
    if 1.0-A <= 0:
        raise ValueError("The retained PG patch requires beta>0.")
    return A, Ap, math.sqrt(1.0-A), math.sqrt(radius2)


def local_ray_frequency(rho: float, energy: float = 1.0) -> float:
    """Positive frequency measured by u=partial_tau-beta partial_rho.

    A future inward radial null ray has conserved Killing frequency E,
    p_rho=-E/(1+beta), and local frequency E+beta*p_rho = E/(1+beta).
    E=1 is a normalization of the ray; ratios are independent of its value.
    The ray approximation can be rescaled to large physical frequencies.
    """
    if energy <= 0:
        raise ValueError("energy must be positive")
    return energy/(1.0+geometry(rho)[2])


def child_data(rho: float) -> dict[str, float]:
    """Geometric quantities in the expanding KS region rho<=0."""
    if rho > 0:
        raise ValueError("Use the expanding child branch rho<=0.")
    A, Ap, beta, radius = geometry(rho)
    if A >= 0:
        raise ValueError("The KS chart requires A<0.")
    a = math.sqrt(-A)
    hp = Ap/(2.0*a)
    ht = -a*rho/(radius*radius)
    proper_time, quadrature_error = quad(
        lambda x: 1.0/math.sqrt(-geometry(x)[0]), rho, 0.0,
        epsabs=2e-13, epsrel=2e-13,
    )
    return {
        "rho": float(rho), "proper_time_from_neck": proper_time,
        "proper_time_quadrature_indicator": quadrature_error,
        "a_parallel": a, "sphere_radius": radius,
        "H_parallel": hp, "H_sphere": ht,
        "shear_squared": (hp-ht)**2/3.0,
        "angular_frequency_k0_kappa1": 1.0/radius,
        "frequency_k1_kappa1": math.hypot(1.0/a, 1.0/radius),
        "frequency_ratio_k1_to_k0_kappa1": math.hypot(radius/a, 1.0),
        "normalized_longitudinal_coefficient_r_over_a": radius/a,
    }


def child_hamiltonian(rho: float, k: float, kappa: int) -> np.ndarray:
    """Exact instantaneous reduced free massless Dirac Hamiltonian in KS time.

    In ds^2=dT^2-a(T)^2 dz^2-r(T)^2 dOmega^2, rescale the spinor by
    r*sqrt(a) to remove the homogeneous spin-connection dilution term.
    In a constant Pauli basis, each angular/longitudinal block is
    H=(k/a)*sigma3+(kappa/r)*sigma1. k is continuous; no radial box is added.
    kappa is a nonzero integer for the untwisted unit-sphere Dirac operator.
    """
    if kappa == 0 or int(kappa) != kappa:
        raise ValueError("kappa must be a nonzero integer")
    A, _, _, radius = geometry(rho)
    if A >= 0:
        raise ValueError("Child Hamiltonian requires A<0.")
    return (k/math.sqrt(-A))*SIGMA3+(kappa/radius)*SIGMA1


def integrate_child_mode(k: float, kappa: int, rho_end: float,
                         rtol: float, atol: float, max_step: float) -> dict:
    """Evolve from the positive instantaneous neck eigenvector.

    The final negative-instantaneous-frequency projection is a mode-basis
    diagnostic, NOT an identified particle-production probability.
    """
    _, eigenvectors = np.linalg.eigh(child_hamiltonian(0.0, k, kappa))
    initial = eigenvectors[:, 1].astype(complex)

    def rhs(rho: float, psi: np.ndarray) -> np.ndarray:
        # dT/drho=-1/a and dpsi/dT=-iH psi, hence dpsi/drho=+iH psi/a.
        a = math.sqrt(-geometry(rho)[0])
        return 1j*(child_hamiltonian(rho, k, kappa) @ psi)/a

    solution = solve_ivp(
        rhs, (0.0, rho_end), initial, method="DOP853", rtol=rtol, atol=atol,
        dense_output=True, max_step=max_step,
    )
    if not solution.success:
        raise RuntimeError(solution.message)
    endpoint = solution.y[:, -1]
    eigenvalues, eigenvectors = np.linalg.eigh(child_hamiltonian(rho_end, k, kappa))
    sample = solution.sol(np.linspace(0.0, rho_end, 1001))
    norm_defect = np.max(np.abs(np.sum(np.abs(sample)**2, axis=0)-1.0))
    return {
        "k": k, "kappa": kappa, "rho_start": 0.0, "rho_end": rho_end,
        "initial_state": "positive instantaneous eigenstate at the minimum sphere",
        "negative_instantaneous_projection": float(abs(np.vdot(eigenvectors[:, 0], endpoint))**2),
        "positive_instantaneous_projection": float(abs(np.vdot(eigenvectors[:, 1], endpoint))**2),
        "final_instantaneous_eigenvalues": eigenvalues.tolist(),
        "sampled_norm_defect": float(norm_defect),
        "rhs_evaluations": solution.nfev,
        "spinor_real": endpoint.real.tolist(), "spinor_imag": endpoint.imag.tolist(),
        "rtol": rtol, "atol": atol, "max_step": max_step,
    }


def integrate_ray(rtol: float, atol: float, max_step: float) -> dict:
    rho_start, rho_end, E = 8.0, -64.0, 1.0
    p0 = -local_ray_frequency(rho_start, E)

    def rhs(tau: float, state: np.ndarray) -> list[float]:
        rho, p = state
        _, Ap, beta, _ = geometry(float(rho))
        beta_prime = -Ap/(2.0*beta)
        return [-1.0-beta, beta_prime*p]

    def at_end(tau: float, state: np.ndarray) -> float:
        return float(state[0]-rho_end)

    at_end.terminal = True
    at_end.direction = -1
    solution = solve_ivp(
        rhs, (0., 20.), [rho_start, p0], events=at_end, dense_output=True,
        method="DOP853", rtol=rtol, atol=atol, max_step=max_step,
    )
    if not solution.success or len(solution.t_events[0]) != 1:
        raise RuntimeError("Ray failed to reach the child endpoint.")
    tau_end = float(solution.t_events[0][0])
    sample = solution.sol(np.linspace(0., tau_end, 2001))
    analytic = np.array([-local_ray_frequency(float(x), E) for x in sample[0]])
    E_samples = np.array([-(1.+geometry(float(x))[2])*p for x, p in sample.T])
    metric_null = np.array([
        1.0-((-1.0-geometry(float(x))[2])+geometry(float(x))[2])**2
        for x in sample[0]
    ])
    return {
        "rho_start": rho_start, "rho_end": float(solution.y[0, -1]),
        "pg_coordinate_arrival_time": tau_end,
        "final_local_frequency_over_E": float(-solution.y[1, -1]),
        "maximum_relative_momentum_vs_analytic": float(np.max(np.abs((sample[1]-analytic)/analytic))),
        "maximum_Killing_energy_defect": float(np.max(np.abs(E_samples-E))),
        "maximum_null_condition_roundoff": float(np.max(np.abs(metric_null))),
        "rhs_evaluations": solution.nfev, "rtol": rtol, "atol": atol, "max_step": max_step,
    }


def compute() -> dict:
    horizon = brentq(lambda x: geometry(x)[0], 1.0, 3.0, xtol=5e-15)
    assert abs(horizon-1.9006916054701435) < 1e-11
    stations = [("parent exterior", 8.), ("event horizon", horizon),
                ("minimum sphere", 0.), ("expanding child", -1.),
                ("expanding child", -4.), ("expanding child", -16.),
                ("expanding child", -64.)]
    nu_parent = local_ray_frequency(8.)
    ray_table = []
    for label, rho in stations:
        A, Ap, beta, radius = geometry(rho)
        frequency = local_ray_frequency(rho)
        ray_table.append({
            "region": label, "rho": rho, "radius": radius, "A": A,
            "beta": beta, "angular_frequency_kappa1": 1.0/radius,
            "angular_frequency_over_parent_angular_frequency": math.sqrt(65.0)/radius,
            "nu_over_Killing_frequency": frequency,
            "nu_over_parent_station_frequency": frequency/nu_parent,
            "ingoing_characteristic_speed": -1.-beta,
        })
    child_table = [child_data(x) for x in [0., -1., -4., -16., -64.]]
    coarse = integrate_ray(1e-9, 1e-11, .05)
    fine = integrate_ray(1e-12, 1e-14, .025)
    mode_rows = []
    for k in (0., 1., 4., 16.):
        lower = integrate_child_mode(k, 1, -16., 1e-9, 1e-11, .1)
        upper = integrate_child_mode(k, 1, -16., 1e-12, 1e-14, .05)
        mode_rows.append({"fine": upper, "coarse": lower,
            "absolute_projection_refinement_difference": abs(
                upper["negative_instantaneous_projection"]-lower["negative_instantaneous_projection"]),
            "endpoint_spinor_l2_refinement_difference": float(np.linalg.norm(
                np.array(upper["spinor_real"])+1j*np.array(upper["spinor_imag"])
                -np.array(lower["spinor_real"])-1j*np.array(lower["spinor_imag"])))})
    # Independent high-precision scalar evaluation, not a continuum certificate.
    mp.mp.dps = 70
    max_geometry_difference = 0.
    for _, x in stations:
        xx = mp.mpf(x)
        aa = 1+3*xx+3*(1+xx*xx)*(mp.atan(xx)-mp.pi/2)
        nn = 1/(1+mp.sqrt(1-aa))
        max_geometry_difference = max(max_geometry_difference, abs(float(nn)-local_ray_frequency(x)))
    # Hermiticity, instantaneous squared dispersion, and a no-mixing control.
    max_hermiticity = max_square_error = max_commutator_error = 0.
    for rho in [0., -1., -4., -16.]:
        A, _, _, r = geometry(rho)
        for k in [0., 1., 4., 16.]:
            H = child_hamiltonian(rho, k, 1)
            max_hermiticity = max(max_hermiticity, float(np.max(abs(H-H.conj().T))))
            max_square_error = max(max_square_error, float(np.max(abs(
                H@H-(k*k/(-A)+1/(r*r))*np.eye(2)))))
            Ap = geometry(rho)[1]
            a = math.sqrt(-A)
            hp, ht = Ap/(2*a), -a*rho/(r*r)
            Hdot = -(k/a)*hp*SIGMA3-(1/r)*ht*SIGMA1
            sigma2 = np.array([[0, -1j], [1j, 0]])
            predicted_commutator = 2j*k/(a*r)*(hp-ht)*sigma2
            max_commutator_error = max(max_commutator_error, float(np.max(abs(
                H@Hdot-Hdot@H-predicted_commutator))))
    assert fine["maximum_relative_momentum_vs_analytic"] < 1e-9
    assert max(row["fine"]["sampled_norm_defect"] for row in mode_rows) < 1e-9
    assert mode_rows[0]["fine"]["negative_instantaneous_projection"] < 1e-20
    return {
        "title": "Two-room frequency transport on the inherited NSC black-universe benchmark",
        "scope": "Fixed-background ray and free Dirac-mode calculation, generated separately from repository evidence",
        "reference_commit": COMMIT,
        "reference_files": [f"{REPO}/blob/{COMMIT}/{path}" for path in REFERENCE_PATHS],
        "primary_reference": "https://arxiv.org/html/gr-qc/0511109v1",
        "units": "c=hbar=L_throat=1; retained benchmark mass m=1",
        "model": {
            "A": "1+3rho+3(1+rho^2)(atan(rho)-pi/2)",
            "r": "sqrt(1+rho^2)", "beta": "sqrt(1-A)",
            "ray_observer": "u=partial_tau-beta*partial_rho; unit freely falling PG normal",
            "ray_transfer": "nu_child/nu_parent=(1+beta_parent)/(1+beta_child)",
            "child_clock": "dT=-drho/sqrt(-A)",
            "child_Hamiltonian": "H_(k,kappa)=(k/sqrt(-A))*sigma3+(kappa/r)*sigma1",
            "child_norm_rescaling": "canonical spinor = r*sqrt(a_parallel)*covariant spinor",
            "child_eigenfrequency_squared": "k^2/(-A)+kappa^2/r^2",
            "child_time_evolution": "i dpsi/dT=H psi; no reset after the initial neck state",
        },
        "horizon_rho": horizon, "ray_stations": ray_table,
        "child_stations": child_table,
        "ray_integration": {"coarse": coarse, "fine": fine,
            "arrival_time_refinement_difference": abs(coarse["pg_coordinate_arrival_time"]-fine["pg_coordinate_arrival_time"])},
        "child_mode_integrations": mode_rows,
        "checks": {
            "high_precision_ray_frequency_absolute_difference": max_geometry_difference,
            "Hamiltonian_hermiticity_maximum": max_hermiticity,
            "Hamiltonian_squared_identity_maximum": max_square_error,
            "geometric_commutator_identity_maximum": max_commutator_error,
            "k0_no_mixing_control": mode_rows[0]["fine"]["negative_instantaneous_projection"],
            "scope": "Finite arithmetic, solver and sampled-norm checks; not a propagated physical-model error bound",
        },
        "not_computed": [
            "self-consistent state--geometry backreaction",
            "full parent-to-child quantum scattering and transported state",
            "observational fitting or quantum-gravity completion",
            "particle creation for a specified in-vacuum and out-detector",
            "a global ancestry clock or conserved global cosmological energy",
            "necessity of nesting or a universal scale factor",
        ],
        "software": {"python": platform.python_version(), "numpy": np.__version__,
                     "scipy": scipy.__version__, "mpmath": mp.__version__},
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error(f"Refusing to replace existing output: {args.output}")
    result = compute()
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write("\n")
    print(f"Saved {args.output}")
    print("rho, ray frequency / parent frequency:")
    for row in result["ray_stations"]:
        print(f"{row['rho']: .9f}: {row['nu_over_parent_station_frequency']:.12g}")
    print("Child rho, proper T, instantaneous frequency ratio (k=1)/(k=0), kappa=1:")
    for row in result["child_stations"]:
        print(row["rho"], row["proper_time_from_neck"], row["frequency_ratio_k1_to_k0_kappa1"])
    print("Mode k, final negative instantaneous projection, refinement difference:")
    for row in result["child_mode_integrations"]:
        print(row["fine"]["k"], row["fine"]["negative_instantaneous_projection"], row["absolute_projection_refinement_difference"])
    print(json.dumps(result["checks"], indent=2))

if __name__ == "__main__":
    main()
