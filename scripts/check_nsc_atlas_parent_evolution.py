"""Compute the finite common-clock parent and its evolving geometric trace."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time

from recursive_horizons.nsc_atlas_parent_evolution import NormalParent, ParentBudget
from recursive_horizons.nsc_atlas_geometric_average import KAPPA, averaged_gravity_balance


ROOT = Path(__file__).resolve().parents[1]


def calculation():
    budget = ParentBudget()
    models = {order: NormalParent(order, budget=budget) for order in (32, 64)}
    times = (0., .005, .01)
    cuts = {order: [model.snapshot(t) for t in times] for order, model in models.items()}
    comparisons = {
        key: max(abs(c["sources"][key]-f["sources"][key])
                 for c, f in zip(cuts[32], cuts[64]))
        for key in cuts[64][0]["sources"]}
    if max(comparisons.values()) > 2e-8:
        raise RuntimeError("normal-parent 32/64 comparison failed")
    balances = []
    mid = cuts[64][1]
    for h in (1e-4, 5e-5):
        lo, hi = (models[64].snapshot(.005+sign*h) for sign in (-1, 1))
        a = mid["sources"]
        rate = lambda k: (hi["sources"][k]-lo["sources"][k])/(2*h)
        balance = averaged_gravity_balance(H=a["H_normal"], Q=a["Q"],
            curvature=a["mean_intrinsic_curvature"], Q_rate=rate("Q"),
            curvature_rate=rate("mean_intrinsic_curvature"),
            mean_energy=a["kappa_mean_energy"]/KAPPA,
            mean_stress_trace=3*a["kappa_mean_pressure"]/KAPPA,
            mean_energy_rate=rate("kappa_mean_energy")/KAPPA)
        balance["time_difference_step"] = h
        balance["scope"] = "Evaluated on the constructed normal-clock parent; derivative differences are numerical indicators."
        balance["volume_rate_residual"] = (
            (hi["moments"][0]-lo["moments"][0])/(2*h)-mid["moments"][1])
        balance["H_rate_residual"] = rate("H_normal")-(
            a["three_volume_acceleration"]/3-a["H_normal"]**2)
        balances.append(balance)
    if abs(balances[-1]["integrability_residual"]) > 3e-9:
        raise RuntimeError("evolved source/geometry balance failed")
    condensed = []
    for r in cuts[64]:
        row = {k: r[k] for k in ("time", "moments", "sources", "walls", "controls")}
        row["central_flash"] = models[64].central_flash(r["time"])
        condensed.append(row)
    rb = models[64].rb
    files = (
        "src/recursive_horizons/nsc_atlas_parent_evolution.py",
        "src/recursive_horizons/nsc_atlas_nested_slice.py",
        "src/recursive_horizons/nsc_atlas_geometric_average.py",
        "src/recursive_horizons/nsc_atlas_nested_spacetime.py",
        "scripts/check_nsc_atlas_parent_evolution.py",
        "tests/test_nsc_atlas_parent_evolution.py")
    return dict(schema="nsc-atlas-parent-evolution-v1",
        inputs=dict(masses=[0., .01, .4], regional_H=[1., 1., 0.], common_k=2.,
                    initial_wall_radii=[.8, 1.5], original_outer_cut=4.,
                    field="phi=0", reference_length="L_ref=ell=2/k",
                    clock="Proper time of geodesics initially normal to the shared slice.",
                    outer_boundary="Original R=4 normal ray, moving after t=0."),
        snapshots=condensed, quadrature_comparisons=comparisons,
        midpoint_differential_controls=balances,
        finite_domain=dict(normal_time=[0., .01], original_wall_proper_limit=.04,
            host_minimum_initial_radius=rb,
            host_positive_radius_lower_bound=rb-2*.01*.01**2/rb**2,
            radial_Jacobi_argument="A(0)=1,B(0)>=0,A_tt=(2M/R^3+H^2)A>=0; each crossing adds k*A/abs(gamma*v) to B.",
            crossing_directions="Cap into host at inner wall; exterior into host at outer wall.",
            source_scale="The supplied regional H and membrane tension are unchanged."),
        interpretation=dict(
            computed="Lambda_g,tr=(Q-mean_R3)/4 on the transported domain.",
            endpoint_Lambda_g_trace_times_Lref2=condensed[-1]["sources"]["geometric_trace_projection"],
            total_trace_including_explicit_sources=condensed[-1]["sources"]["total_trace_projection"],
            generated_bulk_cosmological_constant=False,
            reason="Nonzero trace drift and vacuum-form defect; input cap H is not generated.",
            common_quiet_exterior_nesting_Lambda_difference=0.,
            observational_validation=False),
        controls_scope="32/64 and time-difference changes are numerical indicators, not rigorous rounding or model-error bounds.",
        scalar_evaluations=budget.count,
        process_cpu_seconds=time.process_time()-budget.started,
        new_field_grids=0,
        source_sha256={f: hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in files},
        reproduce="python scripts/run_python.py scripts/check_nsc_atlas_parent_evolution.py")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    encoded = json.dumps(calculation(), indent=2, allow_nan=False)+"\n"
    if len(encoded.encode()) > 24000:
        raise RuntimeError("compact parent record exceeded 24KB")
    if args.output:
        with args.output.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(encoded)
    print(encoded)


if __name__ == "__main__":
    main()
