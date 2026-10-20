"""Evaluate two fixed compact completions and their local-observer restriction."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time

from recursive_horizons.nsc_atlas_compact_parent import (
    CompactParent, reflected_initial_wall, quiet_throat_continuation,
)
from recursive_horizons.nsc_atlas_parent_evolution import ParentBudget
from recursive_horizons.nsc_atlas_geometric_average import KAPPA, averaged_gravity_balance


ROOT = Path(__file__).resolve().parents[1]


def calculation():
    budget = ParentBudget()
    times = (0., .005, .01)
    models = {q: CompactParent(q, budget=budget) for q in (32, 64)}
    runs = {q: {right: [m.snapshot(t, right=right) for t in times]
                for right in ("regular", "nested")} for q, m in models.items()}
    comparisons = {}
    for right in ("regular", "nested"):
        comparisons[right] = {
            key: max(abs(c["sources"][key]-f["sources"][key])
                     for c, f in zip(runs[32][right], runs[64][right]))
            for key in runs[64][right][0]["sources"]}
        if max(comparisons[right].values()) > 2e-8:
            raise RuntimeError("compact-parent 32/64 comparison failed")
    differences = []
    for a, b in zip(runs[64]["regular"], runs[64]["nested"]):
        differences.append(dict(time=a["time"],
            Delta_geometric_trace=b["sources"]["geometric_trace_projection"]-a["sources"]["geometric_trace_projection"],
            Delta_total_trace=b["sources"]["total_trace_projection"]-a["sources"]["total_trace_projection"],
            Delta_curved_volume_fraction=b["curved_region_volume_fraction"]-a["curved_region_volume_fraction"],
            Delta_left_flash_redshift=b["left_flash"]["one_plus_redshift"]-a["left_flash"]["one_plus_redshift"],
            Delta_left_flash_reception=b["left_flash"]["central_reception_time"]-a["left_flash"]["central_reception_time"]))
    balance_records = {}
    model = models[64]
    for right in ("regular", "nested"):
        mid = runs[64][right][1]
        a = mid["sources"]
        records = []
        for h in (1e-4, 5e-5):
            lo, hi = (model.snapshot(.005+sign*h, right=right) for sign in (-1, 1))
            rate = lambda k: (hi["sources"][k]-lo["sources"][k])/(2*h)
            b = averaged_gravity_balance(H=a["H_normal"], Q=a["Q"],
                curvature=a["mean_intrinsic_curvature"], Q_rate=rate("Q"),
                curvature_rate=rate("mean_intrinsic_curvature"),
                mean_energy=a["kappa_mean_energy"]/KAPPA,
                mean_stress_trace=3*a["kappa_mean_pressure"]/KAPPA,
                mean_energy_rate=rate("kappa_mean_energy")/KAPPA)
            records.append(dict(step=h, trace_rate=b["Lambda_trace_rate"],
                source_balance=b["explicit_source_balance"],
                integrability_residual=b["integrability_residual"],
                trace_balance_residual=b["trace_balance_residual"],
                volume_rate_residual=(hi["total_volume"]-lo["total_volume"])/(2*h)-mid["moments"][1],
                H_rate_residual=rate("H_normal")-(a["three_volume_acceleration"]/3-a["H_normal"]**2)))
        if abs(records[-1]["integrability_residual"]) > 2e-7:
            raise RuntimeError("compact evolving conservation failed")
        balance_records[right] = records
    causal = {}
    for right in ("regular", "nested"):
        initial = {side: model.connector_event(0., side=side, right=right)
                   for side in ("left", "right")}
        final = {side: model.connector_event(.01, side=side, right=right)
                 for side in ("left", "right")}
        product = initial["left"]["U"]*initial["right"]["V"]
        if product <= 1:
            raise RuntimeError("initial null-front common-future exclusion failed")
        causal[right] = dict(initial=initial, endpoint=final,
            initial_U_left_times_V_right=product,
            argument="Future causal paths have dU>=0,dV>=0. Opposite exterior quadrants cannot communicate. The two initial inward null fronts require UV>=this product>1, beyond the physical connector R>0.",
            local_signal_comparison="Exact retained left geometry and cap flash; no claimed signal from the remote extra nest.")
    condensed = {}
    for right in ("regular", "nested"):
        condensed[right] = [{k: r[k] for k in (
            "time", "total_volume", "curved_region_volume_fraction",
            "supplied_bulk_trace", "material_wall_trace", "sources", "wall_count",
            "controls", "left_flash")} for r in runs[64][right]]
    files = (
        "src/recursive_horizons/nsc_atlas_compact_parent.py",
        "src/recursive_horizons/nsc_atlas_parent_evolution.py",
        "src/recursive_horizons/nsc_atlas_nested_slice.py",
        "src/recursive_horizons/nsc_atlas_geometric_average.py",
        "scripts/check_nsc_atlas_compact_parent.py",
        "scripts/plot_nsc_atlas_compact_parent.py",
        "tests/test_nsc_atlas_compact_parent.py")
    throat = {q: quiet_throat_continuation(.01, 1., order=q, budget=budget) for q in (32, 64)}
    return dict(schema="nsc-atlas-compact-parent-v1",
        construction=dict(topology="S3", regular_centres=2, asymptotic_ends=0,
            added_bridge_surface_source=0, common_connector_mass=.4,
            common_tension=2., curved_region_H=1., vacuum_connector_H=0.,
            reference_length="L_ref=ell=2/k", scalar_field="zero",
            states="Nested left + regular right; nested left + nested right.",
            clock="Proper time of the shared original normal congruence.",
            supplied="Regional vacuum values, masses, membrane law, initial wall radii and slicing; no observed fraction or target Lambda."),
        right_initial_reflection=reflected_initial_wall(model.regular.initial_walls["outer"]),
        snapshots=condensed, matched_remote_nesting_differences=differences,
        comparisons_32_64=comparisons, midpoint_balance=balance_records,
        causal_comparison=causal,
        quiet_throat=throat[64],
        throat_lifetime_32_64_difference=throat[32]["normal_lifetime"]-throat[64]["normal_lifetime"],
        interpretation="Compactness removes isolated AE dilution; nesting changes global averages. Non-vacuum stress and the exact local causal screen prevent identifying those changes as a generated local cosmological Lambda.",
        uncertainty="Quadrature and derivative differences are numerical indicators, not certified bounds or observational/model uncertainties.",
        scalar_evaluations=budget.count, process_cpu_seconds=time.process_time()-budget.started,
        new_field_grids=0,
        source_sha256={f: hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in files},
        reproduce="python scripts/run_python.py scripts/check_nsc_atlas_compact_parent.py")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    encoded = json.dumps(calculation(), indent=2, allow_nan=False)+"\n"
    if len(encoded.encode()) > 24000:
        raise RuntimeError("compact-completion record exceeds 24KB")
    if args.output:
        with args.output.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(encoded)
    print(encoded)


if __name__ == "__main__":
    main()
