"""Compact evidence for cosmological matching and separate surrounding gravity."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from recursive_horizons.nsc_atlas_parent_observables import (
    DustParent, de_sitter_bulk_slice, population_balance,
    projected_brane_source, shadow_lensing,
    critical_tension_wall,
)

ROOT = Path(__file__).resolve().parents[1]


def calculation():
    started = time.process_time()
    counts = 0
    class CountedParent(DustParent):
        def state(self, tau):
            nonlocal counts
            counts += 1
            return super().state(tau)
    parent = CountedParent()
    samples = [parent.state(t) for t in (0., .01, 1., 4.)]
    arrival = parent.arrival()
    def rhs(t, y):
        return [(1-2*parent.mass/y[0])/parent.state(t)["f"]]
    def contact(t, y):
        return y[0]-parent.state(t)["radius"]
    contact.terminal = True
    contact.direction = 1.
    ode = solve_ivp(rhs, (0., 10.), [2.], events=contact, method="DOP853",
                    rtol=2e-12, atol=2e-13)
    if not ode.success or len(ode.t_events[0]) != 1:
        raise RuntimeError("independent null ODE did not reach the boundary")
    gauss_arrivals = {}
    for order in (32, 64):
        x, w = np.polynomial.legendre.leggauss(order)
        distance = parent.tortoise(4.)-parent.tortoise(2.)
        def phase(t):
            rates = [1/(1+parent.state(s)["radial_velocity"]) for s in (x+1)*t/2]
            return t/2*np.dot(w, rates)-distance
        gauss_arrivals[order] = brentq(phase, 0., 5., xtol=2e-13)
    T = np.array([[.7, .1, -.2, 0.], [.1, .2, .03, 0.],
                  [-.2, .03, .3, .01], [0., 0., .01, -.1]])
    projection = projected_brane_source(T, tension=6.6, kappa5_squared=1., L5=1.,
                                        weyl=np.diag([.6, .1, .2, .3]))
    warp = de_sitter_bulk_slice(np.array([0., .3, .7, 1.]))
    ode_gap = float(ode.t_events[0][0])-arrival["reception_time"]
    if abs(ode_gap) > 3e-10 or projection["projection_residual"] > 2e-14:
        raise RuntimeError("geodesic/projection comparison failed")
    files = (
        "src/recursive_horizons/nsc_atlas_parent_observables.py",
        "scripts/check_nsc_atlas_parent_observables.py",
        "scripts/plot_nsc_atlas_parent_observables.py",
        "tests/test_nsc_atlas_parent_observables.py",
        "src/recursive_horizons/nsc_atlas_nested_slice.py",
    )
    elapsed = time.process_time()-started
    if counts > 200000 or elapsed > 120:
        raise RuntimeError("finite calculation budget exhausted")
    return dict(
        schema="nsc-atlas-parent-observables-v1",
        attachment=dict(model="Four-dimensional flat dust Einstein-Straus exterior",
                        retained_Schwarzschild_mass=.4, initial_join_radius=4.,
                        reference_length="Original ell=2/k; G=c=1",
                        supplied="M, join radius, flat dust and zero parent Lambda",
                        branch="Original increasing-R right exterior beyond the bridge",
                        change="Alternative completion of AE end, not topology change of S3",
                        material_surface_at_join=False, new_vacuum_contribution=0.),
        cosmic_snapshots=samples, exterior_clock=arrival,
        independent_checks=dict(null_ODE_arrival_gap=ode_gap,
                                Gauss_32_64_arrival_gap=gauss_arrivals[32]-gauss_arrivals[64],
                                direct_Gauss_Israel_residual=projection["projection_residual"],
                                bulk_tangential_curvature_residual=float(max(abs(warp["tangential_bulk_curvature"]+1))),
                                bulk_normal_curvature_residual=float(max(abs(warp["normal_bulk_curvature"]+1)))),
        separate_bulk_candidate=dict(
            dimension=5, geometric_Lambda5_L5_squared=-6.,
            supplied_detuning=1.1, induced_Lambda4_L5_squared=float(warp["Lambda4_L5_squared"]),
            H4_L5=math.sqrt(float(warp["H4_squared"])), warp_coordinate_horizon_over_L5=float(warp["coordinate_horizon"]),
            scope="Exact supplied-tension AdS5/dS4 embedding. Not the original four-dimensional shell model or a scale-selection law."),
        shadow_matter=[shadow_lensing(f) for f in (0., .5, 1.)],
        shadow_domain="Positive observer brane, massless unstabilized radion, zero mode, static weak field. Same profile; f is dynamical-potential fraction.",
        open_population_counterexample=[population_balance(a, exponent=3.) for a in (1., 1.5, 2.)],
        original_wall_balance=dict(
            existing_mu=.8, maximum_static_mu=16/27,
            nested_outer_inner_mass=.01, nested_outer_delta_mu=.78,
            nested_outer_potential="Massless-cap V at mu=.78 minus .04/x; its upper bound is already negative for every u>=0.",
            existing_state_has_static_solution_for_any_nonnegative_H_squared=False,
            critical_example=critical_tension_wall(.25),
            proof="At every stationary potential point V_RR=-12A/R^6-6C<0. The critical mass decreases from mu=16/27 as u=H^2/k^2 grows.",
            inherited_invariants="kM and H/k remain free under (R,M,k,H)->(sR,sM,k/s,H/s)."),
        uncertainty="Fixed-rule and ODE differences are numerical indicators; no observed fit or physical error distribution is supplied.",
        charged_parent_state_evaluations=counts, process_cpu_seconds=elapsed,
        field_grids=0,
        source_sha256={f: hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in files},
        reproduce="python scripts/run_python.py scripts/check_nsc_atlas_parent_observables.py",
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    encoded = json.dumps(calculation(), indent=2, allow_nan=False)+"\n"
    if len(encoded.encode()) > 16000:
        raise RuntimeError("summary exceeds 16KB")
    if args.output:
        with args.output.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(encoded)
    print(encoded)


if __name__ == "__main__":
    main()
